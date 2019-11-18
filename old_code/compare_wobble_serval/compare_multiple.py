
"""

@author: cmatthe
"""
#%% imports
import sys
sys.path.append('/data/cmatthe/python/wobble_aux')# path to timported python fiels containing auxiliary functions, this  path only works on lx39
import rv_solver as rv

import matplotlib.pyplot as plt #for plt
from matplotlib.backends.backend_pdf import PdfPages # for PdfPages: suspiciously abscent in akaminskis code
#I suspect there is a hidden standart import for example also for numpy?
# perhaps in  sys.path.append('/home/akaminsk/pythontests')
import matplotlib as mpl #missing in akaminski
import h5py 
import numpy as np
import scipy as sp #for realigning time axis with a fit
from tqdm import tqdm
import barycorrpy as bary 
from scipy.constants import codata 
lightvel = codata.value('speed of light in vacuum') #for barycorr

import pickle # for writing barycorr results into a binary file

from astropy.stats import sigma_clip

from time import sleep
import os


#%%

# CAHA Coordinates for barycorr
_lat = 37.2236
_lon = -2.54625
_elevation = 2168.

objects = [["GJ436", "J11421+267", "vis", 2.644]]
orbital_parameters_mult = [
    [17.38, 2.644, 0.152, 325.8, 2450000 + 1552.077, 78.3]
    ]#(Try with M0 offset of 78.3 deg (see Carm paper))
#orbital_parameters = [17.38, 2.643859, 0.152, 325.8, 2450000 + 1551.72] #use these
#orbital_parameters = [K, P, e, omega, T0, MO]

#####
load_rvs = False #If true will load RVs from pickle saved files if they exist
use_avcn = True #Select whether to use NZP corretion if True or not if False
error_clipping = True
primary_planet_phased_plot = True

#select False if drift correction was already applied to spectrum
correct_w_for_drift = False

#####
#min_sleep = 30
#print("sleeping for {0} minutes".format(min_sleep))
#sleep(60*min_sleep)
#####

#####
if use_avcn == True:
    avcn_tag = "_avcn_"
else:
    avcn_tag = ""


#recalculate_baryQ = False #TODO reimplement pickle saving of barycorrected RVs
kt = "Kt3"
#kt = "Kt3__reg_def_chunk_roll1"
#kt = "Kt3_rs_orderwise_test_0_reg"
#kt = "Kt3_loop_2_reg"
#kt = 'Kt3_adrianutc_attached'
#kt='Kt3_order_39_test'
#kt='Kt3_opt3'
extra_identifier = "_snr+ds+cont1_regs"
pp =PdfPages("wobble_multi_comparison_" +objects[0][0] + avcn_tag + kt + extra_identifier+".pdf")
fig = plt.figure(figsize=(15, 9), dpi=200)
mpl.rc('font', size=16)
plt.clf()
fig.clf()
ax1=plt.gca()

#Generate list of result files to be included in comparison
dir_wobble = "/data/cmatthe/compare_wobble_serval/wobbledir/"

res_file_list = []
loop_nr_list = [i for i in range(0,6)]
#TODO put loop nr identifiers in the plot titles
for l in loop_nr_list:
    #next_file = dir_wobble + "results_GJ436_Kstar0_Kt3_loop_{0}_reg.hdf5".format(l)
    next_file = dir_wobble + "results_GJ436_Kstar0_Kt3_snr+ds_l{0}_reg.hdf5".format(l)
    if res_file_list == []:
        res_file_list = [next_file]
    else:
        res_file_list.append(next_file)



for i in objects:
    print(i[0], i[1], i[2])
    for f, fil in enumerate(res_file_list):
        #read in wobble result
        if i[2]=="vis":
            #servaldir = "/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
            #to run from lx39, perhas replace this with  an if statement
            servaldir_lx39="/data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
            servaldir=servaldir_lx39
            #wobble_res = h5py.File("/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_Kstar0_Kt3.hdf5", 'r')
            wobble_res_lx39 = h5py.File(fil,'r')
            wobble_res=wobble_res_lx39
            #N_ord=42 # note CARMENES has 61 orderes only 11-52 used here NOTE better to use N_ord supplied by file (note that this is actually R
            N_ord = len(wobble_res["orders"][()]) #should be equivalent to ["R"]
            k=np.arange(N_ord)
            keys_orderRVs=[]
            for j in k:
                keys_orderRVs.append('order'+str(j))
                
        else:
            #servaldir = "/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/servaldir/CARM_NIR/"
            #to run from lx39
            servaldir_lx39="/data/cmatthe/compare_wobble_serval/servaldir/CARM_NIR/"
            servaldir=servaldir_lx39
            #wobble_res = h5py.File("/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_nir_Kstar0_Kt3.hdf5")
            wobble_res_lx39 = h5py.File(fil,'r')
            wobble_res=wobble_res_lx39
            N_ord = len(wobble_res["orders"][()]) #should be equivalent to ["R"]
            k=np.arange(N_ord)
            keys_orderRVs=[]
            for j in k:
                keys_orderRVs.append('order'+str(j))


            
        #w_dates = wobble_res['dates'].value
        w_dates = wobble_res['dates'][()] 
        w_dates_utc = wobble_res['dates_utc'][()]
        w_orders = wobble_res['orders'][()]
        w_epochs = wobble_res['epochs'][()]
        #use selfcombined RVs a wobble RVs seem broken currently
        w_RVs = wobble_res['star_time_rvs'][()]
        w_RVs_original = w_RVs
        w_RVs_er = wobble_res['star_time_sigmas'][()]
        
        
        #barycorr for wobble_orig
        if 'bary_starname' in locals():
            w_RVs_original_barycorr = np.zeros(len(w_dates))
            for n in tqdm(range(len(w_RVs_original_barycorr))):
                w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel)[0]
        else:
            w_RVs_original_barycorr = np.zeros(len(w_dates))
            for n in tqdm(range(len(w_RVs_original_barycorr))):
                w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel)[0]
            
            
        #read in SERVAL
        ser_rvc = np.loadtxt(servaldir+i[1]+"/"+i[1]+".rvc.dat") # rvc
        if use_avcn == True:
            ser_avcn = np.loadtxt(servaldir+i[1]+"/"+i[1]+".avcn.dat")
            ser_rvc = ser_avcn#HACK
        # remove entries with nan in drift
        ind_finitedrift = np.isfinite(ser_rvc[:,3])
        ser_rvc = ser_rvc[ind_finitedrift]
        
        #correction list to apply serval correctios to Wobble
        #add to wobble to be comparable to ser_rvc (clearly not addapted to orderwise yet)
        # be  careful to subtract out mean correction for normalized data # NOTE this si still not a perfect method
        ser_corr = - ser_rvc[:,8] - ser_rvc[:,3]
        if use_avcn == True:
            #8 sa drift, 3 drift, 9 NZP
            ser_corr = - ser_rvc[:,8] - ser_rvc[:,3] - ser_rvc[:,9] # use with avcn
            ser_corr_wob = ser_corr
        
    #optionally remove drift correction if this has been performed during data file generation
        if correct_w_for_drift == False:
            ser_corr_wob = ser_corr + ser_rvc[:,3]
            
        #plot of RV differences
        #for  this we need to match the observation by the JDs -> start with wobble date and find the one with the lowest timediff from serval
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
        ser_corr = ser_corr[indices_serval] # corrcetion: ser_rvc[:,5] + ser_corr = ser_rvc[:,1] within 10^-10ms^-1
        ser_corr_wob = ser_corr_wob[indices_serval]
        
        w_dates = w_dates[indices_wobble]
        w_dates_utc = w_dates_utc[indices_wobble]
        w_epochs = w_epochs[indices_wobble]
        w_RVs_er = w_RVs_er[indices_wobble]
        w_RVs_original = w_RVs_original[indices_wobble]                         + ser_corr_wob
        w_RVs_original_barycorr = w_RVs_original_barycorr[indices_wobble]       + ser_corr_wob
        
        if len(i) >= 4:
            #NOTE Fixed:This should be in comparison to orbital solution vrad = K[cos(/theta + /omega) + e * cos(/omega)]
                            
            sigma_wob = np.nanstd(sigma_clip(
            w_RVs_original_barycorr -np.nanmean(w_RVs_original_barycorr)
            ,sigma = 5) )
            sigma_wob_no_clip = np.nanstd(
            w_RVs_original_barycorr -np.nanmean(w_RVs_original_barycorr))
            
            sigma_ser = np.nanstd(sigma_clip(
            ser_rvc[:,1] - np.nanmean(ser_rvc[:,1])
            ,sigma = 5) )
            sigma_ser_no_clip = np.nanstd(
            ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]))
            #Time is incorrectly aligned. probably a result of rounding period and periastron time. Idea: re-fit T0 for both wobble and Serval
                  
            def error_stats(residuals, errors, kappa = 5):
                residual_offset = np.average(residuals)
                # note standart deviation autopmatically corrects for the offset
                sigma = np.nanstd(residuals)
                sigma_clipped = np.nanstd(sigma_clip(residuals,sigma = 5) )
                ind_bad = np.where(np.absolute(residuals - residual_offset) > kappa * errors)[0]
                ind_good= np.where(np.absolute(residuals - residual_offset) <= kappa * errors)[0]
                sigma_error_clipped = np.nanstd(residuals[ind_good])
                return sigma, sigma_clipped,sigma_error_clipped, ind_bad, ind_good, residual_offset
           
                
            if 'orbital_parameters_mult' in locals():
                def keplarian_rv_mult(t):
                    total_rv = 0
                    for parameters in orbital_parameters_mult:
                        total_rv = total_rv + rv.radial_velocity_M0(t , parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
                    return total_rv
                
                def fit_func_mult(t, T0_offset):
                    return keplarian_rv_mult(t + T0_offset)
                
                # For single planet work e.g. phase sychronized data
                
                def keplarian_rv(t):
                    parameters = orbital_parameters_mult[0]
                    return rv.radial_velocity_M0(t , parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
                
                def fit_func(t, T0_offset):
                    return keplarian_rv(t + T0_offset)
                
                
            #with fitted model (x and (y offset removed))
                #fit to serval:
                xdata = ser_rvc[:,0]
                ydata = ser_rvc[:,1] - np.nanmean(ser_rvc[:,1])
                popt_s, pcov_s = sp.optimize.curve_fit(fit_func_mult, xdata, ydata, sigma = ser_rvc[:,2], absolute_sigma = True
                #,p0 = [0,residual_offset], bounds = ([-orbital_parameters_mult[0][1],residual_offset - 0.1*orbital_parameters_mult[0][0]],[orbital_parameters_mult[0][1], residual_offset + 0.1*orbital_parameters_mult[0][0]])
                ) # NOTE HACK p0 choice based on wobble can fail
                print("T0_offset Serval = ", popt_s)
                T0_offset_s = popt_s[0]
                
                #fit to wobble:
                xdata = w_dates
                ydata = w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)
                popt, pcov = sp.optimize.curve_fit(fit_func_mult, xdata, ydata, sigma = w_RVs_er, absolute_sigma = True
                #,p0 = [0, residual_offset], bounds = ([-orbital_parameters_mult[0][1],residual_offset - 0.1*orbital_parameters_mult[0][0]],[orbital_parameters_mult[0][1], residual_offset + 0.1*orbital_parameters_mult[0][0]])
                ) # NOTE HACK p0 choice based on wobble can fail
                print("T0_offset Wobble = ", popt)
                T0_offset = popt[0]
                
                wob_residuals_sign_mult = w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr) - fit_func_mult(w_dates, popt[0]
                            #, popt[1]
                            )
                wob_mult_stats = error_stats(wob_residuals_sign_mult, w_RVs_er)
                residual_offset = wob_mult_stats[5] # Extremely similar offsets for both wob and ser)
                
                #print(np.average(wob_residuals_sign_mult))
                #print( np.nanmean(w_RVs_original_barycorr))
                
                ser_residuals_sign_mult = ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) - fit_func_mult(ser_rvc[:,0], popt_s[0]
                #, popt_s[1]
                )
                ser_mult_stats = error_stats(ser_residuals_sign_mult, ser_rvc[:,2])
                #print(np.average(ser_residuals_sign_mult))
                
                w_dates_sorted = np.sort(w_dates)
                xlst = np.linspace(w_dates_sorted[0], w_dates_sorted[-1], num = 3000)
                ax1.plot(xlst , fit_func_mult(xlst, popt_s[0]
                                                #,popt_s[1]
                                                ), "k-", label = "lit orbit M0 (serval_fit)")
                ax1.errorbar(w_dates, (w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)) - residual_offset , yerr = w_RVs_er,fmt = "x", label="Wobble_Corr, no_clip ={0:.3f},clipped_sigma = {1:.3f}, err_clipped = {2:.3f} ".format(wob_mult_stats[0], wob_mult_stats[1], wob_mult_stats[2]))
                ax1.plot(ser_rvc[:,0], ser_rvc[:,5] - np.nanmean(ser_rvc[:,5]) - residual_offset, "x", label= "SERVAL")
                ax1.errorbar(ser_rvc[:,0], ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) - residual_offset,yerr = ser_rvc[:,2] ,fmt = "x", label= "SERVAL_Corr, no_clip ={0:.3f}, clipped_sigma = {1:.3f}, err_clipped = {2:.3f}".format(ser_mult_stats[0], ser_mult_stats[1], ser_mult_stats[2]))
                if error_clipping == True:
                    ind_bad_wob = wob_mult_stats[3]
                    ind_bad_ser = ser_mult_stats[3]
                    ax1.plot((w_dates[ind_bad_wob]), (w_RVs_original_barycorr[ind_bad_wob] - np.nanmean(w_RVs_original_barycorr)), "o",color = "C0", label = "error_clipped wobble points")
                    ax1.plot((ser_rvc[:,0][ind_bad_ser]), ser_rvc[:,1][ind_bad_ser] - np.nanmean(ser_rvc[:,1]), "o", color = "C2", label= "error_clipped serval points")
                #ax1.plot(w_dates % i[3], (w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)),"x", label="Wobble_Corr")
                ax1.set_ylabel("RVs [m/s]")
                ax1.set_xlabel('jd')
                plt.title("RVs and multiplanet fit for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
                plt.suptitle(os.path.basename(fil), y = 1.0)
                plt.grid(True)
                plt.tight_layout()
                plt.legend(shadow=True)
                plt.savefig(pp, format='pdf')
                plt.clf()
                fig.clf()
                ax1 = plt.gca()
                
                if primary_planet_phased_plot == True:
                    #NOTE uses same fit as multi planet model
                    orbital_parameters = orbital_parameters_mult[0]
                    
                    wob_residuals_sign = w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr) - fit_func(w_dates, popt[0]
                            #, popt[1]
                            )
                    wob_stats = error_stats(wob_residuals_sign, w_RVs_er)
                    residual_offset = wob_mult_stats[5] # Extremely similar offsets for both wob and ser)
                    
                    #print(np.average(wob_residuals_sign_mult))
                    #print( np.nanmean(w_RVs_original_barycorr))
                    
                    ser_residuals_sign = ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) - fit_func_mult(ser_rvc[:,0], popt_s[0]
                    #, popt_s[1]
                    )
                    ser_stats = error_stats(ser_residuals_sign, ser_rvc[:,2])
                    
                    
                    xlst = np.linspace(w_dates[0], w_dates[0] + orbital_parameters[1]*0.99999, num=100)
                    ylst = [keplarian_rv(t + T0_offset_s)  for t in xlst]
                    #sort by xlst
                    pltlst = [[xlst[j],ylst[j]] for j in range(len(xlst))]
                    def mod_sort(elem):
                        return elem[0] % orbital_parameters[1]
                    pltlst = sorted(pltlst, key = mod_sort)
                    pltlst = np.asarray(pltlst)
                    pltlst = [pltlst[:,0],pltlst[:,1]]
                    #to_print = pltlst[0] % i[3]
                    #print(to_print)
                    
                    ax1.plot(pltlst[0] % i[3], pltlst[1], "r-", label = "literature orbit (Serval T0_offset) single planet")
                    #ax1.plot(pltlst[0] % i[3], keplarian_rv(pltlst[0]), "k-", label = "lit orbit M0 (no fit) single planet")
                            
                    #ax1.plot(w_dates % i[3], (w_RVs_original+w_bervs-np.nanmean(w_RVs_original+w_bervs)),"x", label="Wobble_orig")
                    ax1.errorbar((w_dates) % i[3], (w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)), yerr = w_RVs_er,fmt = "x", label="Wobble_Corr, no_clip ={0:.3f},clipped_sigma = {1:.3f}, err_clipped = {2:.3f} ".format(wob_stats[0], wob_stats[1], wob_stats[2]))
                    ax1.plot((ser_rvc[:,0]) % i[3], ser_rvc[:,5] - np.nanmean(ser_rvc[:,5]), "x", label= "SERVAL")
                    ax1.errorbar((ser_rvc[:,0]) % i[3], ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]),yerr = ser_rvc[:,2] ,fmt = "x", label= "SERVAL_Corr, no_clip ={0:.3f}, clipped_sigma = {1:.3f}, err_clipped = {2:.3f}".format(ser_stats[0], ser_stats[1], ser_stats[2]))
                    if error_clipping == True:
                        ind_bad_wob = wob_stats[3]
                        ind_bad_ser = ser_stats[3]
                        ax1.plot((w_dates[ind_bad_wob]) % i[3], (w_RVs_original_barycorr[ind_bad_wob] - np.nanmean(w_RVs_original_barycorr)), "o",color = "C0", label = "error_clipped wobble points")
                        ax1.plot((ser_rvc[:,0][ind_bad_ser]) % i[3], ser_rvc[:,1][ind_bad_ser] - np.nanmean(ser_rvc[:,1]), "o", color = "C2", label= "error_clipped serval points")
                    #ax1.plot(w_dates % i[3], (w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)),"x", label="Wobble_Corr")
                    ax1.set_ylabel("RVs [m/s]")
                    ax1.set_xlabel('jd')
                    plt.title("Phased ("+str(i[3])+"d) RVs for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
                    plt.suptitle(os.path.basename(fil), y = 1)
                    plt.grid(True)
                    plt.tight_layout()
                    plt.legend(shadow=True)
                    plt.savefig(pp, format='pdf')
                    plt.clf()
                    fig.clf()
                    ax1 = plt.gca()
            
plt.close(fig)
pp.close()
                        
                    
