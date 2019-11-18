import sys
sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
import combine_results as cr
from parameters import Parameters
import rv_solver as rv

#end local imports
import os
from time import time
import wobble
import h5py
import numpy as np
import scipy as sp #for realigning time axis with a fit
import shutil
import yaml
#from compare_ws
from tqdm import tqdm
import barycorrpy as bary 
import matplotlib as mpl
import matplotlib.pyplot as plt #for plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.stats import sigma_clip
#from distutils.dir_util import copy_tree # not currently used
from scipy.constants import codata 

#if __name__ == "__main__": #NOTE If called via os.system this will trigger even if in queue
    #queue = True
    ##################### parameters
    
    #if queue:
        #parameter_filename = "yaml_temp/optimize_parameters.yaml"
        #parameters = Parameters(filename = parameter_filename)
        
        #starname = parameters.dictionary["starname"]
        #K_star = parameters.dictionary["K_star"]
        #K_t = parameters.dictionary["K_t"]
        #niter = parameters.dictionary["niter"]
        #start = parameters.dictionary["start"]
        #end = parameters.dictionary["end"]
        #chunk_size = parameters.dictionary["chunk_size"]
        #defaultQ = parameters.dictionary["defaultQ"]
        #if defaultQ:
            #default_str = "_def"
        #else:
            #default_str = ""
            
        #parameter_dict = parameters.dictionary    
        #if "reg_file_star" in parameter_dict and "reg_file_t" in parameter_dict:
            #star_reg_file = parameter_dict["reg_file_star"]
            #tellurics_reg_file = parameter_dict["reg_file_t"]
    #else:
        #starname = 'GJ876'
        #K_star = 0
        #K_t = 3
        #niter = 300
        #start = 11 #first order 
        #end = 14 # first order that is *not* included
        #chunk_size = 2 #use at least 2 so median is not out of bounds #NOTE: currrently needs to match reg file chunk_size
        ##wether to use default regularization
        #defaultQ = True
        #if defaultQ:
            #default_str = "_def"
        #else:
            #default_str = ""
            
#Pseudocode /Plan:
#-1- Create directory tree: all data products should be in one file
#0. Load base reg file
#1. Schedule optimizing the star with reg parametersshifted +-1 for all except L2_basis_weights: 6 * 3 optmizations
#   1.1. Generate test regulrization files 
#2. Compare each of the 6 sets for RV residual rms and chose the best one in  each case:
#   2.1. Generate new baseline regulrization file 
    
def update_reg_file(base_reg_file, new_reg_file, orders,  parameter_key, exponent):
    """
    Creates new regularization files with the selected parameter shifted by the exponent
    base_reg_file : `str`
        file path of the base regularization file
    new_reg_file : `str`
        file path of the new regularization file
    orders : list of `int`
        orders for which the regularization file should be adjusted (assumes files runing from order 0 to 60)
    parameter_key: `str`
        key of the parameter to be adjusted
    exponent: `int`
        exponent shift of the regularization parameter
    """
    factor = 10**exponent
    with h5py.File(base_reg_file,'r') as f:
        if not parameter_key in list(f.keys()):
            raise Exception("Invalid parameter_key in update_reg_file")
        with h5py.File(new_reg_file,'w') as g:
            for key in list(f.keys()):
                if key == parameter_key:
                    temp = f[key][()]
                    for o in  orders:
                        temp[o] = temp[o] * factor
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
                else: 
                    temp = f[key][()]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
                    
 #TODO  wite the compare function                   
def compare_results(file_list, parameter_change_list, bary_starname, orbital_parameters, objects, servaldir):
    """
    Compares result files to find the best rv scatter around literature fit returns change in parameter that yielded best results
    
    file_list: list of `str`
        list containing the file paths of the files to be compared
    parameter_change_list: list of `int`
        list containing the parameter exponent shifts used to create the files in file_list
    orbital_parameters : list of `float`
        orbital_parameters = [K, P, e, omega, T0]
        parameters of the keplerian fit to be used as "true" baseline
    """
    sigma_list = np.zeros(len(file_list)) + 100 # 100 is a fudge factor
    plots = True
    if plots:
        rec_loop_directory, key_name = os.path.split(os.path.split(file_list[0])[0]) # HACK go  up 2 directories to loop directory
        plot_directory = rec_loop_directory + "/compare_plots"
        os.makedirs(plot_directory, exist_ok = True)
        
        pp =PdfPages(plot_directory +"/"+ key_name  +".pdf")
        fig = plt.figure(figsize=(15, 9), dpi=200)
        mpl.rc('font', size=16)
        plt.clf()
        fig.clf()
        ax1=plt.gca()
    
    for f, fil in enumerate(file_list):
        #assumes order of file_listand parameter_change_list are matched. (maybe extract from file name?)
        wobble_res = h5py.File(fil,'r')
        w_dates = wobble_res['dates'][()] 
        w_dates_utc = wobble_res['dates_utc'][()]
        
        w_RVs = wobble_res['star_time_rvs'][()]
        w_RVs_original = w_RVs
        w_RVs_er = wobble_res['star_time_sigmas'][()]
        
        #barycorr for wobble_orig
        from scipy.constants import codata 
        lightvel = codata.value('speed of light in vacuum') #for barycorr
        # CAHA Coordinates for barycorr
        _lat = 37.2236
        _lon = -2.54625
        _elevation = 2168.
        
        w_RVs_original_barycorr = np.zeros(len(w_dates))
        for n in tqdm(range(len(w_RVs_original_barycorr))):
            w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel)[0]
            
        #Serval Correction 
        #read in SERVAL
        ser_rvc = np.loadtxt(servaldir+objects[1]+"/"+objects[1]+".rvc.dat")
        # remove entries with nan in drift
        ind_finitedrift = np.isfinite(ser_rvc[:,3])
        ser_rvc = ser_rvc[ind_finitedrift]
        ser_corr = - ser_rvc[:,8] - ser_rvc[:,3]
        #match wobble and serval
        indices_serval = [] 
        indices_wobble = []
        for n in range(len(w_dates)):
            ind_jd = np.where(np.abs(ser_rvc[:,0]-w_dates[n]) == np.nanmin(np.abs(ser_rvc[:,0]-w_dates[n])))[0][0]
            if (ser_rvc[ind_jd,0]-w_dates[n])*24*60<20.: #only takes matches closer than 20 minutes
                indices_serval.append(ind_jd)
                indices_wobble.append(n)
        print("#serval_ind:"+str(len(indices_serval)), "#wobble_ind:"+str(len(indices_wobble)))
        #now set up all the data according to the indices
        ser_rvc = ser_rvc[indices_serval]
        ser_corr = ser_corr[indices_serval]
        
        w_dates = w_dates[indices_wobble]
        w_dates_utc = w_dates_utc[indices_wobble]
        w_RVs_original_barycorr = w_RVs_original_barycorr[indices_wobble]       + ser_corr
        w_RVs_er = w_RVs_er[indices_wobble]
        
        def fit_func(t, T0_offset):
                return rv.radial_velocity(t , orbital_parameters[0], orbital_parameters[1], orbital_parameters[2],orbital_parameters[3], orbital_parameters[4] + T0_offset)
            
        #fit to Wobble
        xdata = w_dates
        ydata = w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)
        popt, pcov = sp.optimize.curve_fit(fit_func, xdata, ydata,  sigma = w_RVs_er, absolute_sigma = True)
        print("T0_offset Wobble = ", popt)
        T0_offset = popt[0]
        
        #make these weighted (maybe: thsi may not be a good idea if residuals are not strongly correlated to error (as with wobble results))
        sigma_wob = np.nanstd(sigma_clip(
        w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr) - fit_func(w_dates, T0_offset)
        ,sigma = 5))
        sigma_list[f] = sigma_wob
        
        sigma_wob_noclip = np.nanstd(
        w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr) - fit_func(w_dates, T0_offset)
        )
        if plots:
            #fit to serval:
            xdata = ser_rvc[:,0]
            ydata = ser_rvc[:,1] - np.nanmean(ser_rvc[:,1])
            popt_s, pcov_s = sp.optimize.curve_fit(fit_func, xdata, ydata, sigma = ser_rvc[:,2], absolute_sigma = True)
            print("T0_offset Serval = ", popt_s)
            T0_offset_s = popt_s[0]
            
            sigma_ser = np.nanstd(sigma_clip(
            ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) - fit_func(ser_rvc[:,0], T0_offset_s)
            ,sigma = 5) )
            
            sigma_ser_noclip = np.nanstd(
            ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) - fit_func(ser_rvc[:,0], T0_offset_s)
            )
            
            xlst = np.linspace(w_dates[0], w_dates[0] + orbital_parameters[1]*0.99999, num=100)
            ylst = [rv.radial_velocity(t , orbital_parameters[0], orbital_parameters[1], orbital_parameters[2],orbital_parameters[3], orbital_parameters[4] + T0_offset) for t in xlst]
            #sort by xlst
            pltlst = [[xlst[j],ylst[j]] for j in range(len(xlst))]
            def mod_sort(elem):
                return elem[0] % orbital_parameters[1]
            pltlst = sorted(pltlst, key = mod_sort)
            pltlst = np.asarray(pltlst)
            pltlst = [pltlst[:,0],pltlst[:,1]]
            
            
            ax1.plot(pltlst[0] % orbital_parameters[1], pltlst[1], "r-", label = "literature orbit (Wobble T0_offset)")
            ax1.errorbar((w_dates) % orbital_parameters[1], (w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)), yerr = w_RVs_er,fmt = "x", label="Wobble_Corr, clipped_sigma = {0:.3f}, noclip = {1:.3f} ".format(sigma_wob, sigma_wob_noclip))
            ax1.errorbar((ser_rvc[:,0]) % orbital_parameters[1], ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]),yerr = ser_rvc[:,2] ,fmt = "x", label= "SERVAL_Corr, clipped_sigma = {0:.3f}, noclip = {1:.3f}".format(sigma_ser, sigma_ser_noclip), color = "C2")
            ax1.set_ylabel("RVs [m/s]")
            ax1.set_xlabel('jd')
            # add the parameter change to the title
            title_pre = os.path.split(os.path.split(fil)[0])[1]
            plt.title(title_pre + ", Phased ("+str(orbital_parameters[1])+"d) RVs for "+objects[0]+" ("+objects[2]+") "+" - "+objects[1]+";")
            plt.grid(True)
            plt.tight_layout()
            plt.legend(shadow=True)
            plt.savefig(pp, format='pdf')
            plt.clf()
            fig.clf()
            ax1 = plt.gca()
        
        
    
    if plots:# include some nice progress plots. TODO make it not crudely placed inside this function?
        plt.close(fig)
        pp.close()
    
    best_index = np.argmin(sigma_list)
    return parameter_change_list[best_index]
        
        
        
    
if __name__ == "__main__":
    #def_star_filename = '../wobble/regularization/default_star.hdf5'
    #def_t_filename = '../wobble/regularization/default_t.hdf5'
    
    #new_test_file = '../wobble/regularization/update_test_star.hdf5'
    #orders = [o for o  in range(11,42)]
    
    #update_reg_file(def_star_filename, new_test_file, orders, "L1_template", 1) 
    
    ############################
    run_name = "test_run_0.1"#check where this  is used
    orbital_parameters = [17.38, 2.643859, 0.152, 325.8, 2450000 + 1551.72]#for compare 
    #orbital_parameters = [K, P, e, omega, T0]
    bary_starname = "GJ436" #for barycorr (simbad database name)
    objects = ["GJ436", "J11421+267", "vis", 2.643859] # for loading serval correction
    servaldir = "/data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"#used to correct data before compare
    
    #this only works with a previously created e.g. started from queue or first creating it here
    starting_star_reg = "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_1_star.hdf5"
    starting_t_reg = "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_1_t.hdf5"
    
    parameters = Parameters(starname = "GJ436", defaultQ = False, 
                        start = 11,
                        end = 53,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  starting_star_reg,
                        reg_file_t = starting_t_reg,
                        output_suffix = 'reg_search_test')
    
    parameter_filename = "yaml_temp/optimize_parameters.yaml"
    parameters.write(parameter_filename)
    parameters = Parameters(filename = parameter_filename)
    
    starname = parameters.dictionary["starname"] 
    K_star = parameters.dictionary["K_star"]
    K_t = parameters.dictionary["K_t"]
    niter = parameters.dictionary["niter"]
    output_suffix = parameters.dictionary["output_suffix"]
    start = parameters.dictionary["start"]
    end = parameters.dictionary["end"]

    #directory for lx39
    top_directory = "/data/cmatthe/wobble_reg_search/"
    if os.path.exists(top_directory):
        parameters.dictionary.update({"top_directory" : top_directory})
    else:
        raise Exception("top directory {0} does not exist".format(top_directory))
    
    # Currently created in top: will probably need change
    #results_directory = top_directory + "/results_{0}{1}/".format(starname, output_suffix)
    #parameters.dictionary.update({"results_directory" : results_directory})
    
    data_directory = "/data/cmatthe/wobble_data/data/"
    if os.path.exists(data_directory):
        parameters.dictionary.update({"data_directory" : data_directory})
    else:
        raise Exception("data directory {0} does not exist".format(data_directory))
    
    
    run_directory = top_directory + "{0}_{1}/".format(starname, run_name)
    os.makedirs(run_directory, exist_ok = True)
    parameters.dictionary.update({"run_directory" : run_directory})
    
    parameters.write(parameter_filename)
    
    #TODO :
    # Put reg fies at the right places for uptdate_reg to find them
    
    max_loop_iterations = 5
    parameter_changes = [-1,0,1]
    reg_types = ["star", "t"] #t for tellurics NOTE(DO NOT CHANGE THIS NOMENCLATURE IT IS HARD CODED LATER ON)HACK
    #############################
    
    base_star_reg = starting_star_reg
    base_t_reg = starting_t_reg
    
    #immplemented Jul 23: skip optmization if files already exist in case of rerun
    
    
    start_time = time()
    for i in range(max_loop_iterations):
        loop_directory = run_directory + "loop_{0}/".format(i)
        os.makedirs(loop_directory, exist_ok = True)
        
        for j,typ in enumerate(reg_types):
            if typ == "star":
                parameter_key_list = ["L2_template", "L1_template"]
                parameters.dictionary.update({"reg_file_t" : base_t_reg})
                base_reg_file = base_star_reg
            if typ == "t": # t= tellurics
                parameter_key_list = ["L2_template", "L1_template", "L2_basis_vectors", "L1_basis_vectors"]
                parameters.dictionary.update({"reg_file_star" : base_star_reg})
                base_reg_file = base_t_reg
            
            # HACK write current next_base_reg_file to a temporary reg file, for later update the next_base_reg_file
            temp_reg_file = loop_directory + "temp_{0}_reg".format(typ) + ".hdf5"
            update_reg_file(base_reg_file, temp_reg_file, orders = [0,1,2], parameter_key = "L2_template", exponent = 0)# just makes a copy
            
            for k,key in enumerate(parameter_key_list):
                file_list = ['' for change in parameter_changes]
                for c, change in enumerate(parameter_changes):
                    #create results_directory
                    results_directory = loop_directory + "results_{2}_{0}_{1}/".format(key, change, typ)# NOTE this name is hard coded below
                    os.makedirs(results_directory, exist_ok = True)
                    parameters.dictionary.update({"results_directory" : results_directory})
                    
                    #create regularization directory, updated reg files and pass them to the parameter file
                    regularization_directory = results_directory + "regularization/"
                    os.makedirs(regularization_directory, exist_ok = True)
                    new_reg_file = regularization_directory + "reg_{0}_{1}_{2}".format(typ, key, change) + ".hdf5"
                    orders = [o for o in range(start, end)] #HACK placeholder before this is changed to orderwise
                    update_reg_file(base_reg_file, new_reg_file, orders, parameter_key = key, exponent = change)
                    if typ == "star":
                        parameters.dictionary.update({"reg_file_star" : new_reg_file})
                    if typ == "t":
                        parameters.dictionary.update({"reg_file_t" : new_reg_file})
                    parameters.write(parameter_filename)
                    
                    #check if file already exists, and if so skip
                    #HACK FILE name should not be hardcoded as it is determined in opt_top
                    results_file_base = results_directory + 'results_{0}_Kstar{1}_Kt{2}_'.format(starname, K_star, K_t, niter) + parameters.dictionary["output_suffix"] 
                    results_file = results_file_base + '.hdf5'
                    if os.path.isfile(results_file):
                        latest_results_file = results_file
                    else:
                        #if it exists, copy existing results with identical base_reg. (for 0 change)
                        #name must be passed up from opt_top
                        if change == 0 and not (j == 0 and k == 0):
                            latest_results_file = parameters.dictionary["latest_results_file"]
                            file_basename = os.path.basename(latest_results_file)
                            file_to_copy = loop_directory + "results_{2}_{0}_{1}/".format(parameter_key_list[0], change, reg_types[0]) + file_basename
                            copy_at_dest = results_directory + file_basename
                            shutil.copy(file_to_copy, copy_at_dest) #if you just supplied the destination directory basename would be reused anyways
                        else:
                            #run wobble_optimization
                            parameters.dictionary.update({"script" : "optimize_top_3.2_reg.py"})
                            script_name = parameters.dictionary["script"]
                            #write parameters before switching script, so they can be loaded
                            parameters.write(parameter_filename)
                            os.system("python3 "  + script_name)
                            #read back parameter file, so that changes done in sub_script are updated
                            parameters = Parameters(filename = parameter_filename)
                        #insert correct filepath into file_list for compare_results    
                        latest_results_file = parameters.dictionary["latest_results_file"] #passed up from opt_top
                    file_list[c] = latest_results_file
                    print("finished loop {0}, typ {1}, key {2}, change {3}".format(i, typ, key, change))
                    print("time elapsed gradient_search: {0:.2f} h".format((time() - start_time)/3600.0))
                
                #Compare results, find "best", update reg
                # Comparison will initially result based on standart deviation from literature planet fit
                best_change = compare_results(file_list, parameter_changes, bary_starname, orbital_parameters, objects, servaldir)
                
                #update next_base_reg_file (for next loop)
                next_base_reg_file = loop_directory + "next_base_{0}_reg.hdf5".format(typ)
                update_reg_file(temp_reg_file, next_base_reg_file, orders, parameter_key = key, exponent = best_change)
                # HACK write current next_base_reg_file to a temporary reg file, for later update the next_base_reg_file
                temp_reg_file = loop_directory + "temp_{0}_reg".format(typ) + ".hdf5"
                update_reg_file(next_base_reg_file, temp_reg_file, orders = [0,1,2], parameter_key = "L2_template", exponent = 0)# just makes a copy
        
        #set new base_reg_file
        base_star_reg = loop_directory + "next_base_star_reg.hdf5"
        base_t_reg = loop_directory + "next_base_t_reg.hdf5"     
    

    
    
