#Optimize a set of regularization parameters for a given star using RV deviation from a "known" planet solution as the criterion
#basis is reg_gradient_search_old_template.py
#changes to implement:
#-only fit Serval (vis) T0_offset. (all deviations from this should be negligible or erraneous)
#-for nir match  serval order rvs to semi orders
#-adapt fit function for multiple planets

#local imports
import rv_solver as rv
from run_wobble import *
from process_results import Results_ws
#imports
import sys
import os
from time import time
import wobble
import h5py
import numpy as np
import scipy as sp #for realigning time axis with a fit
import shutil # TODO chekc if needed

#from compare_ws
from tqdm import tqdm
import barycorrpy as bary
# HACK https://github.com/astropy/astropy/issues/9427 workaround for broken iers mirror
from astropy.utils import iers
iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')
import matplotlib as mpl
import matplotlib.pyplot as plt #for plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.stats import sigma_clip
#from distutils.dir_util import copy_tree # not currently used
from scipy.constants import codata 

import warnings # to catch optimize warnings as errors
from scipy.optimize import OptimizeWarning

from time import sleep


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
                    
                    
def find_Serval_T0_offset(orbital_parameters_mult, objects, servaldir):
    #orbital parameters used as outside reference tend to have offset time axis, this function estimates this ofset with a fit
    
    ser_avcn = np.loadtxt(servaldir+objects[0][1]+"/"+objects[0][1]+".avcn.dat")
    # remove entries with nan in drift
    ind_finitedrift = np.isfinite(ser_avcn[:,3])
    ser_avcn = ser_avcn[ind_finitedrift]
    
    def keplarian_rv_mult(t):
        #Calculates keplarian rv for multiple supplied sets of planet parameters 
        total_rv = 0
        for parameters in orbital_parameters_mult:
            total_rv = total_rv + rv.radial_velocity_M0(t , parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
        return total_rv
    
    def fit_func(t, T0_offset):
        return keplarian_rv_mult(t + T0_offset)
    
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        try:
            #fit to Serval
            xdata = ser_avcn[:,0]
            ydata = ser_avcn[:,1] - np.nanmean(ser_avcn[:,1])
            popt, pcov = sp.optimize.curve_fit(fit_func, xdata, ydata,  sigma = ser_avcn[:,2], absolute_sigma = True)
            print("T0_offset Serval (wobble_fit_err) = ", popt)
            T0_offset = popt[0]
            return T0_offset
            #T0_source = "SERVAL (wob_fit_err)"
        except OptimizeWarning:
            print(OptimizeWarning)
            sys.exit("Serval T0 offset fit failed")
            
def compare_results(file_list, parameter_change_list, bary_starname, orbital_parameters_mult, objects, servaldir, serval_T0_offset):
    """
    Compares result files to find the best rv scatter around literature fit returns change in parameter that yielded best results
    
    file_list: list of `str`
        list containing the file paths of the files to be compared
    parameter_change_list: list of `int`
        list containing the parameter exponent shifts used to create the files in file_list
    orbital_parameters_mult : list of `float`
        orbital_parameters = [K, P, e, omega, T0]
        parameters of the keplerian fit to be used as "true" baseline
    """
    def keplarian_rv_mult(t):
        #Calculates keplarian rv for multiple supplied sets of planet parameters 
        total_rv = 0
        for parameters in orbital_parameters_mult:
            total_rv = total_rv + rv.radial_velocity_M0(t , parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
        return total_rv
    
    def fit_func(t, T0_offset):
        return keplarian_rv_mult(t + T0_offset)
    
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
        res = Results_ws(wobble_file = fil
                 , serval_dir = servaldir
                 , carmenes_object_ID = objects[0][1]
                 , bary_starname = bary_starname
                 , load_bary = False
                 , archive = False)
        res.apply_corrections()
        w_RVs_barycorr = res.w_RVs_barycorr
        w_dates = res.w_dates
        w_RVs_er = res.w_RVs_er
        #skip fitting test RV curve to wobble data. Just use SERVAL T0 offset
        T0_offset = serval_T0_offset
        sigma_wob = np.nanstd(sigma_clip(
        w_RVs_barycorr - np.nanmean(w_RVs_barycorr) - fit_func(w_dates, T0_offset)
        ,sigma = 5))
        sigma_list[f] = sigma_wob

        sigma_wob_noclip = np.nanstd(
        w_RVs_barycorr - np.nanmean(w_RVs_barycorr) - fit_func(w_dates, T0_offset)
        )
        
        if plots:
            T0_offset_s = T0_offset
            ser_avcn  = res.ser_avcn
            sigma_ser = np.nanstd(sigma_clip(
            ser_avcn[:,1] - np.nanmean(ser_avcn[:,1]) - fit_func(ser_avcn[:,0], T0_offset_s)
            ,sigma = 5) )
            
            sigma_ser_noclip = np.nanstd(
            ser_avcn[:,1] - np.nanmean(ser_avcn[:,1]) - fit_func(ser_avcn[:,0], T0_offset_s)
            )
            
            sigma_wob_Soffset = np.nanstd(sigma_clip(
            w_RVs_barycorr - np.nanmean(w_RVs_barycorr) - fit_func(w_dates, T0_offset_s)
            ,sigma = 5))
            sigma_list[f] = sigma_wob
            
            sigma_wob_noclip_Soffset = np.nanstd(
            w_RVs_barycorr - np.nanmean(w_RVs_barycorr) - fit_func(w_dates, T0_offset_s)
            )
            
            xlst = np.linspace(w_dates[0], w_dates[0] + orbital_parameters_mult[0][1]*0.99999, num=100)
            ylst = [fit_func(t, T0_offset) for t in xlst]
            #sort by xlst
            pltlst = [[xlst[j],ylst[j]] for j in range(len(xlst))]
            def mod_sort(elem):
                return elem[0] % orbital_parameters_mult[0][1]
            pltlst = sorted(pltlst, key = mod_sort)
            pltlst = np.asarray(pltlst)
            pltlst = [pltlst[:,0],pltlst[:,1]]
            
            T0_source = "SERVAL "
            
            ax1.plot(pltlst[0] % orbital_parameters_mult[0][1], pltlst[1], "r-", label = "literature orbit (" + T0_source + "T0_offset)")
            ax1.errorbar((w_dates) % orbital_parameters_mult[0][1], (w_RVs_barycorr-np.nanmean(w_RVs_barycorr)), yerr = w_RVs_er,fmt = "x", label="Wobble_Corr, clipped_sigma = {0:.3f}, noclip = {1:.3f} ".format(sigma_wob, sigma_wob_noclip))
            ax1.errorbar((ser_avcn[:,0]) % orbital_parameters_mult[0][1], ser_avcn[:,1] - np.nanmean(ser_avcn[:,1]),yerr = ser_avcn[:,2] ,fmt = "x", label= "SERVAL_Corr, clipped_sigma = {0:.3f}, noclip = {1:.3f}".format(sigma_ser, sigma_ser_noclip), color = "C2")
            ax1.plot([], [], ' ', label="Wobble_Corr_SERVAL_fit, clipped_sigma = {0:.3f}, noclip = {1:.3f} ".format(sigma_wob_Soffset, sigma_wob_noclip_Soffset))
            ax1.set_ylabel("RVs [m/s]")
            ax1.set_xlabel('jd')
            # add the parameter change to the title
            title_pre = os.path.split(os.path.split(fil)[0])[1]
            plt.title(title_pre + ", Phased ("+str(orbital_parameters_mult[0][1])+"d) RVs for "+objects[0][0]+" ("+objects[0][2]+") "+" - "+objects[0][1]+";")
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
    
                    
def regularization(parameters
                   ,run_name
                   ,objects
                   ,orbital_parameters_mult
                   ,bary_starname
                   ,servaldir
                   ,starting_star_reg
                   ,starting_t_reg
                   ,top_directory
                   ,data_directory
                   ,max_loop_iterations
                   ,parameter_changes
                   ):
    starname = objects[0][0]
    p = parameters
    
    reg_types = ["star", "t"]#t for tellurics NOTE(DO NOT CHANGE THIS NOMENCLATURE IT IS HARD CODED LATER ON)HACK
    serval_T0_offset = find_Serval_T0_offset(orbital_parameters_mult, objects, servaldir)
    base_star_reg = starting_star_reg
    base_t_reg = starting_t_reg
    
    #Generate global_epochs_list and attach to parameters object NOTE TODO must be attached to paameters object that is passed to run_wobble later
    data_file= data_directory + starname + p.data_suffix + '_e2ds.hdf5'
    data = wobble.Data(data_file, orders = np.arange(p.start, p.end), min_flux=10**-5, min_snr = p.min_snr,
                       parameters = p #passed to continuum_fit as **kwargs
                       )
    global_epochs_list = p.global_epochs_list = data.epochs.tolist()
    global_orders_list = p.global_orders_list = data.orders.tolist()
    orders = global_orders_list # Drop SNR limited orders #HACK bandaid fix to remove troublesome orders in  IR.
    
    #setup output directory for this run 
    run_directory = top_directory + "{0}_{1}/".format(starname, run_name)
    os.makedirs(run_directory, exist_ok = True)
    
    start_time = time()
    for i in range(max_loop_iterations):
        #initialize parameters for this loop from input parameters object
        parameters_loop = p
        #create loop_directory
        loop_directory = run_directory + "loop_{0}/".format(i)
        os.makedirs(loop_directory, exist_ok = True)
        for j,typ in enumerate(reg_types):
            if typ == "star":
                parameter_key_list = ["L2_template", "L1_template"]
                base_reg_file = base_star_reg
            if typ == "t": # t = tellurics
                parameter_key_list = ["L2_template", "L1_template", "L2_basis_vectors", "L1_basis_vectors"]
                base_reg_file = base_t_reg
            #HACK Save a copy of the current (from previous loop) base reg file to a temporary file where it can be updated to the next_base_reg_file 
            temp_reg_file = loop_directory + "temp_{0}_reg".format(typ) + ".hdf5"
            update_reg_file(base_reg_file, temp_reg_file, orders = [0,1,2], parameter_key = "L2_template", exponent = 0)# just makes a copy (inputs except for filenames are dummy placeholders)
            
            for o, order in enumerate(orders):
            #update loop parameters to only run 1 order
                parameters_loop.start = order
                parameters_loop.end = order + 1
                for k,key in enumerate(parameter_key_list):
                    file_list = ['' for change in parameter_changes]
                    for c, change in enumerate(parameter_changes):
                        #create results_directory
                        results_directory = loop_directory + "results_{0}_ord{1}_{2}_{3}/".format(typ, order, key, change)# NOTE this name format is hard coded below
                        os.makedirs(results_directory, exist_ok = True)
                        parameters_loop.results_dir = results_directory
                        #create regularization directory, updated reg files and pass them to the parameter object
                        regularization_directory = results_directory + "regularization/"
                        os.makedirs(regularization_directory, exist_ok = True)
                        new_reg_file = regularization_directory + "reg_{0}_ord{1}_{2}_{3}".format(typ,order,key, change) + ".hdf5"
                        update_reg_file(base_reg_file, new_reg_file, orders = [order], parameter_key = key, exponent = change)
                        if typ == "star":
                            parameters_loop.reg_file_star = new_reg_file
                        if typ == "t":
                            parameters_loop.reg_file_t = new_reg_file
                        
                        #check if file already exists, and if so skip(for time savings on reruns)
                        #NOTE HACK Assumes same hard coded file format as in chunk.py
                        results_file_base = results_directory + 'results_{0}_Kstar{1}_Kt{2}_'.format(starname, p.K_star, p.K_t, p.niter) + p.output_suffix
                        results_file = results_file_base + '.hdf5'
                        file_basename = os.path.basename(results_file)
                        if os.path.isfile(results_file):
                            latest_results_file = results_file
                            '''check if unreliable
                        else:
                            #if already present, copy existing results with identical regularization (for 0 change)
                            if change == 0 and not (j == 0 and k == 0):
                                file_to_copy = loop_directory + "results_{0}_ord{1}_{2}_{3}/".format(reg_types[0], order, parameter_key_list[0], change) + file_basename
                                copy_at_dest = results_directory + file_basename
                                shutil.copy(file_to_copy, copy_at_dest)
                                #change latest_results_file to copy_at_dest so it gets properly  put into file_list loaded in compare_results
                                latest_results_file = copy_at_dest
                                '''
                        else:
                            #run wobble_optimization
                            run_wobble(parameters_loop)
                            latest_results_file = results_directory + file_basename
                        file_list[c] = latest_results_file # TODO THis feature really needs an overhaul
                        print("finished loop {0}, typ {1}, key {2}, change {3}".format(i, typ, key, change))
                        print("time elapsed gradient_search: {0:.2f} h".format((time() - start_time)/3600.0))
                    
                    #Compare results, find "best", update reg
                    # Comparison will initially result based on standart deviation from literature planet fit
                    best_change = compare_results(file_list, parameter_changes, bary_starname, orbital_parameters_mult, objects, servaldir, serval_T0_offset)
                    
                    #update next_base_reg_file (for next loop)
                    next_base_reg_file = loop_directory + "next_base_{0}_reg.hdf5".format(typ)
                    update_reg_file(temp_reg_file, next_base_reg_file, orders = [order], parameter_key = key, exponent = best_change)
                    # HACK write current next_base_reg_file to a temporary reg file, for later update the next_base_reg_file
                    temp_reg_file = loop_directory + "temp_{0}_reg".format(typ) + ".hdf5"
                    update_reg_file(next_base_reg_file, temp_reg_file, orders = [0,1,2], parameter_key = "L2_template", exponent = 0)# just makes a copy
        
        #set new base_reg_file
        base_star_reg = loop_directory + "next_base_star_reg.hdf5"
        base_t_reg = loop_directory + "next_base_t_reg.hdf5"     
        #TODO Output some reg file graphs
            
        
                    
if __name__ == "__main__":
    #test_0
    '''
    run_name = "test_1_no_copy"
    objects = [["GJ1148", "J11417+427", "nir", 41.382]]
    orbital_parameters_mult = [[38.37, 41.380, 0.380, 258.1, 2450000 + 1581.046, 299.0],
                           [11.34, 532.58, 0.342, 210.4, 2450000 + 1581.046, 272.6]
                            ]# "known" planet solution to compare with
    bary_starname = "GJ1148"
    servaldir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/"#Always use vis data as comparison, since its more reliable for fitting the T0 offset (may not be true for very late star types)
    
    #chose which regularization files to start from
    starting_star_reg = 'regularization/flat_reg_star.hdf5'
    starting_t_reg = 'regularization/flat_reg_t.hdf5'
    
    #set up directory structure
    top_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/regularization/"
    os.makedirs(top_directory, exist_ok = True)
    #if not os.path.exists(top_directory):
        #raise Exception("top directory {0} does not exist".format(top_directory))
    data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    
    max_loop_iterations = 6
    parameter_changes = [-1,0,1]

    #TODO generate parameters file within regularization for individual order runs
    parameters = Parameters(starname = "GJ1148",
                        data_suffix = "_nir_drift_shift_split",
                        start = 20,
                        end = 21,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  starting_star_reg,
                        reg_file_t = starting_t_reg,
                        output_suffix = run_name,
                        results_dir = '../results/',
                        data_dir= '../data/'
                        )
    #This function can take a LONG time 13*max_loop_iterations*run_wobble_time
    regularization(parameters
                   ,run_name
                   ,objects
                   ,orbital_parameters_mult
                   ,bary_starname
                   ,servaldir
                   ,starting_star_reg
                   ,starting_t_reg
                   ,top_directory
                   ,data_directory
                   ,max_loop_iterations
                   ,parameter_changes
                   )
                   '''
    """
    #NIR Test
    run_name = "all_orders"
    objects = [["GJ1148", "J11417+427", "nir", 41.382]]
    orbital_parameters_mult = [[38.37, 41.380, 0.380, 258.1, 2450000 + 1581.046, 299.0],
                           [11.34, 532.58, 0.342, 210.4, 2450000 + 1581.046, 272.6]
                            ]# "known" planet solution to compare with
    bary_starname = "GJ1148"
    servaldir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/"#Always use vis data as comparison, since its more reliable for fitting the T0 offset (may not be true for very late star types)
    
    #chose which regularization files to start from
    starting_star_reg = 'regularization/flat_reg_star.hdf5'
    starting_t_reg = 'regularization/flat_reg_t.hdf5'
    
    #set up directory structure
    top_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/regularization/"
    os.makedirs(top_directory, exist_ok = True)
    #if not os.path.exists(top_directory):
        #raise Exception("top directory {0} does not exist".format(top_directory))
    data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    
    max_loop_iterations = 6
    parameter_changes = [-1,0,1]

    #TODO generate parameters file within regularization for individual order runs
    parameters = Parameters(starname = "GJ1148",
                        data_suffix = "_nir_drift_shift_split",
                        start = 0,
                        end = 56,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  starting_star_reg,
                        reg_file_t = starting_t_reg,
                        output_suffix = run_name,
                        results_dir = '../results/',
                        data_dir= '../data/'
                        )
    #This function can take a LONG time 13*max_loop_iterations*run_wobble_time
    regularization(parameters
                   ,run_name
                   ,objects
                   ,orbital_parameters_mult
                   ,bary_starname
                   ,servaldir
                   ,starting_star_reg
                   ,starting_t_reg
                   ,top_directory
                   ,data_directory
                   ,max_loop_iterations
                   ,parameter_changes
                   )
    """
    
    ##GJ1148 with log(reg)=0 as seed
    #run_name = "no_reg_seed"
    #objects = [["GJ1148", "J11417+427", "vis", 41.382]]
    #orbital_parameters_mult = [[38.37, 41.380, 0.380, 258.1, 2450000 + 1581.046, 299.0],
                           #[11.34, 532.58, 0.342, 210.4, 2450000 + 1581.046, 272.6]
                            #]# "known" planet solution to compare with
    #bary_starname = "GJ1148"
    #servaldir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/"#Always use vis data as comparison, since its more reliable for fitting the T0 offset (may not be true for very late star types)
    
    ##chose which regularization files to start from
    #starting_star_reg = 'regularization/dummy_star_K0_no_reg.hdf5'
    #starting_t_reg = 'regularization/dummy_t_K3_no_reg.hdf5'
    
    ##set up directory structure
    #top_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/regularization/"
    #os.makedirs(top_directory, exist_ok = True)
    ##if not os.path.exists(top_directory):
        ##raise Exception("top directory {0} does not exist".format(top_directory))
    #data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    
    #max_loop_iterations = 6
    #parameter_changes = [-1,0,1]

    ##TODO generate parameters file within regularization for individual order runs
    #parameters = Parameters(starname = "GJ1148",
                        #data_suffix = "_vis_drift_shift",
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  starting_star_reg,
                        #reg_file_t = starting_t_reg,
                        #output_suffix = run_name,
                        #results_dir = '../results/',
                        #data_dir= '../data/'
                        #)
    ##This function can take a LONG time 13*max_loop_iterations*run_wobble_time
    #regularization(parameters
                   #,run_name
                   #,objects
                   #,orbital_parameters_mult
                   #,bary_starname
                   #,servaldir
                   #,starting_star_reg
                   #,starting_t_reg
                   #,top_directory
                   #,data_directory
                   #,max_loop_iterations
                   #,parameter_changes
                   #)
#BEGIN
   #additional loops for GJ436 optimization
    run_name = "add_loops"
    objects = [["GJ436", "J11421+267", "vis", 2.644]]
    orbital_parameters_mult = [[17.38, 2.644, 0.152, 325.8, 2450000 + 1552.077, 78.3]]# "known" planet solution to compare with
    bary_starname = "GJ436"
    servaldir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/"#Always use vis data as comparison, since its more reliable for fitting the T0 offset (may not be true for very late star types)
    
    #chose which regularization files to start from
    starting_star_reg = '../../wobble_reg_search/GJ436_orderwise_snr+drift_shift_1/loop_5/next_base_star_reg.hdf5'
    starting_t_reg = '../../wobble_reg_search/GJ436_orderwise_snr+drift_shift_1/loop_5/next_base_t_reg.hdf5'
    
    #set up directory structure
    top_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/regularization/"
    os.makedirs(top_directory, exist_ok = True)
    #if not os.path.exists(top_directory):
        #raise Exception("top directory {0} does not exist".format(top_directory))
    data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    
    max_loop_iterations = 3
    parameter_changes = [-1,0,1]

    #TODO generate parameters file within regularization for individual order runs
    parameters = Parameters(starname = "GJ436",
                        min_snr = 60, # Reverted to 60.
                        data_suffix = "_vis_drift+nzp",
                        start = 11,
                        end = 53,
                        chunk_size = 1,
                        niter = 160,
                        reg_file_star =  starting_star_reg,
                        reg_file_t = starting_t_reg,
                        output_suffix = run_name,
                        results_dir = '../results/',
                        data_dir= '../data/'
                        )
    #This function can take a LONG time 13*max_loop_iterations*run_wobble_time
    regularization(parameters
                   ,run_name
                   ,objects
                   ,orbital_parameters_mult
                   ,bary_starname
                   ,servaldir
                   ,starting_star_reg
                   ,starting_t_reg
                   ,top_directory
                   ,data_directory
                   ,max_loop_iterations
                   ,parameter_changes
                   )

