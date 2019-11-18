
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


#%%

# CAHA Coordinates for barycorr
_lat = 37.2236
_lon = -2.54625
_elevation = 2168.

# from carmenes data itself
#orbital_parameters = [37.02, 41.382, 0.392, 253.6, 2450000 + 1581.046, 303.9] #303.9
#objects = [["GJ1148", "J11417+427", "vis", 41.382]]
#alternate 2 planet dyn fit (as keplarian)
#orbital_parameters_mult = [[38.37, 41.380, 0.380, 258.1, 2450000 + 1581.046, 299.0],
                           #[11.34, 532.58, 0.342, 210.4, 2450000 + 1581.046, 272.6]
                            #]

#objects = [["GJ876", "J22532-142", "vis", 61.082]]
#orbital_parameters = [212.07, 61.082, 0.027, 35.1, 2450000 + 602.093, 341.1]
#orbital_parameters_mult = [[212.07, 61.082, 0.027, 35.1, 2450000 + 602.093, 341.1],
                            #[88.34, 30.126, 0.250, 51.6, 2450000 + 602.093 , 293.3],
                            #[6.14, 1.938, 0.082, 272.8, 2450000 + 602.093, 316.7],
                            #[3.39, 124.4, 0.040, 263.6, 2450000 + 602.093, 310.3]
                             #]

objects = [["GJ436", "J11421+267", "vis", 2.644]]
orbital_parameters = [17.38, 2.644, 0.152, 325.8, 2450000 + 1552.077, 78.3] #(Try with M0 offset of 78.3 deg (see Carm paper)) #currently usig these
orbital_parameters_mult = [[17.38, 2.644, 0.152, 325.8, 2450000 + 1552.077, 78.3]]
#orbital_parameters = [17.38, 2.643859, 0.152, 325.8, 2450000 + 1551.72] 
#orbital_parameters = [K, P, e, omega, T0, MO]


#objects = [["Wolf294", "J06548+332", "vis", 14.237]]
#orbital_parameters = [0.00, 14.237, 0.152, 325.8, 2450000 + 1551.72] #HACK I just set the RV amplitude to 0

#Teegarden catalogue name:
#bary_starname = "GAT 1370"
#orbital_parameters = [2.02, 4.91, 0.00, 77, 2458000 + 52.3]
#objects = [["Teegarden", "J02530+168", "vis", 4.91]]

#objects = [["GJ285", "J07446+035", "vis", 2.777]]


#####
load_rvs = True #If true will load RVs from pickle saved files if they exist
use_avcn = True #Select whether to use NZP corretion if True or not if False
error_clipping = True

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
#kt = "Kt3__reg_def_chunk_roll1"
#kt = "Kt3_rs_orderwise_test_0_reg"
kt = "Kt3_loop_7_reg_n1000"
#kt = "Kt3_loop_2_reg"
#kt = 'Kt3_adrianutc_attached'
#kt='Kt3_order_39_test'
#kt='Kt3_opt3'
extra_identifier = ""
pp =PdfPages("wobble_servalcomparison_" +objects[0][0] + avcn_tag + kt + extra_identifier+".pdf")
fig = plt.figure(figsize=(15, 9), dpi=200)
mpl.rc('font', size=16)
plt.clf()
fig.clf()
ax1=plt.gca()

#read in wobble and serval and compare RVs from Wobble (barycentric corrected) with SERVAL RVs, but also with RVCs and AVCs 

for i in objects:
    print(i[0], i[1], i[2])
    #read in wobble result
    if i[2]=="vis":
        #servaldir = "/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
        #to run from lx39, perhas replace this with  an if statement
        servaldir_lx39="/data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
        servaldir=servaldir_lx39
        #wobble_res = h5py.File("/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_Kstar0_Kt3.hdf5", 'r')
        wobble_res_lx39 = h5py.File("/data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_Kstar0_"+kt+".hdf5",'r')
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
        wobble_res_lx39 = h5py.File("/data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_nir_Kstar0_"+kt+".hdf5",'r')
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
    #Initialize arrays
    w_order_RVs = np.zeros((N_ord,len(w_dates)))
    w_order_RVs_ivars = np.zeros((N_ord,len(w_dates)))
    w_order_RV_scatter = np.zeros(len(w_dates))
    for j, k in enumerate(keys_orderRVs):
        temp = wobble_res[k] #load order dataset
        tempRVs = temp['star_rvs'][()] # load RVs for this order at all epochs
        tempRVs_ivars = temp['star_ivars_rvs'][()] # load error of RVs
        w_order_RVs[j,:] = tempRVs # fill values into Rows of array
        w_order_RVs_ivars[j,:] = tempRVs_ivars
        if np.any(tempRVs_ivars < 0.0)== True:
            print("Found negative ivars.")
        #we may want to drop  these orders as they may represent *maxima* in optimization
    for j in range (len(w_dates)):
        w_order_RV_scatter[j] = np.nanstd(w_order_RVs[:,j])
        
    #TEST 
    #w_order_RVs = w_order_RVs - np.nanmean(w_order_RVs)
    
    #initialize array for w_RVs_own which circumvent the broken combine in wobble
    w_RVs_own = w_RVs * 0.0
    w_RVs_unweighted = w_RVs * 0.0
    for ep in range (len(w_RVs_own)):
        # use the ivars; if they are negative then do not use the order. Otherwise use weighted mean with sigma= np.sqrt(1./ivars)
        t_RV = []
        t_RV_sigma = []
        for ord in range(len(w_order_RVs[:,ep])):
            if w_order_RVs_ivars[ord, ep] > 0.0:
                t_RV.append(w_order_RVs[ord, ep])
                t_RV_sigma.append(np.sqrt(1./w_order_RVs_ivars[ord, ep]))
        t_RV = np.asarray(t_RV)
        t_RV_sigma = np.asarray(t_RV_sigma)
        w_RVs_own[ep] = np.average(t_RV, weights=1./t_RV_sigma**2)
        w_RVs_unweighted[ep] = np.nanmedian(t_RV)  #alternate using unweighted median
        
    #w_RVs = w_RVs_own #overwrite the wobble combined RVs with self Combined
    w_RVs = w_RVs_unweighted #overwrite the wobble combined RVs with unweighted self Combined
    
    w_bervs = wobble_res['bervs'][()]
    #make own (proper) barycentic corrections
    w_RVs_barycorr = np.zeros(len(w_dates))
    
    #load barycorrected rvs from pickle archive to save time on recalculting them
    if load_rvs:
        try:
            f = open("rv_files/"+i[0] + "_" + i[2]+ avcn_tag +kt+ "_bervs_unweight.pickle", "rb")
            w_RVs_barycorr = pickle.load(f)
            f.close()
            
            f = open("rv_files/"+i[0]+"_"+i[2]+ avcn_tag+kt+"_bervs_ivar.pickle", "rb")
            w_RVs_barycorr_ivar = pickle.load(f)
            f.close()
            
            f = open("rv_files/"+i[0]+"_"+i[2]+ avcn_tag +kt+"_bervs_orig.pickle", "rb")
            w_RVs_original_barycorr = pickle.load(f)
            f.close()
            
            f = open("rv_files/"+i[0]+"_"+i[2]+ avcn_tag +kt+"_bervs_orders.pickle", "rb")
            w_order_RVs = pickle.load(f)
            w_order_RVs_barycorr = w_order_RVs
            f.close()
            
        except FileNotFoundError:
            print("RV file not found, recalculating barycentric corrections")
            load_rvs = False
            
    if not load_rvs:
        if 'bary_starname' in locals():
            if True:
                for n in tqdm(range(len(w_RVs_barycorr))):
                    w_RVs_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs[n]/lightvel)[0]
                
            w_RVs_barycorr_ivar = np.zeros(len(w_dates))
          
            if True:
                for n in tqdm(range(len(w_RVs_barycorr_ivar))):
                    w_RVs_barycorr_ivar[n]=bary.get_BC_vel(w_dates_utc[n], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_own[n]/lightvel)[0]
    
            
            #barycorr for wobble_orig
            w_RVs_original_barycorr = np.zeros(len(w_dates))
            if True:
                for n in tqdm(range(len(w_RVs_original_barycorr))):
                    w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel)[0]
    
                
            #barycorr for w_orders (Takes a while ca 4 mins)
            w_order_RVs_barycorr = w_order_RVs #just to initialize array shape
            if True:
                for order in tqdm(range(len(w_order_RVs))):
                    for ep in tqdm(range(len(w_order_RVs[order]))):
                        w_order_RVs_barycorr[order, ep] = bary.get_BC_vel(w_dates_utc[ep], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=w_order_RVs[order, ep]/lightvel)[0]
        else:
            if True:
                for n in tqdm(range(len(w_RVs_barycorr))):
                    w_RVs_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs[n]/lightvel)[0]
            
                
            w_RVs_barycorr_ivar = np.zeros(len(w_dates))
           
            if True:
                for n in tqdm(range(len(w_RVs_barycorr_ivar))):
                    w_RVs_barycorr_ivar[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_own[n]/lightvel)[0]
              
            #barycorr vor wobble_orig
            w_RVs_original_barycorr = np.zeros(len(w_dates))
            if True:
                for n in tqdm(range(len(w_RVs_original_barycorr))):
                    w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel)[0]
                
            #barycorr for w_orders (Takes a while ca 4 mins)
            w_order_RVs_barycorr = w_order_RVs
            if True:
                for order in tqdm(range(len(w_order_RVs))):
                    for ep in tqdm(range(len(w_order_RVs[order]))):
                        w_order_RVs_barycorr[order, ep] = bary.get_BC_vel(w_dates_utc[ep], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_order_RVs[order, ep]/lightvel)[0]
        
    #export to pickle
    
        f = open("rv_files/"+i[0]+"_"+i[2]+ avcn_tag +kt+"_bervs_unweight.pickle", "wb")
        pickle.dump(w_RVs_barycorr, f)
        f.close()
        
        f = open("rv_files/"+i[0]+"_"+i[2]+ avcn_tag +kt+"_bervs_ivar.pickle", "wb")
        pickle.dump(w_RVs_barycorr_ivar, f)
        f.close()
        
        f = open("rv_files/"+i[0]+"_"+i[2]+ avcn_tag +kt+"_bervs_orig.pickle", "wb")
        pickle.dump(w_RVs_original_barycorr, f)
        f.close()
        
        f = open("rv_files/"+i[0]+"_"+i[2]+ avcn_tag +kt+"_bervs_orders.pickle", "wb")
        pickle.dump(w_order_RVs , f)
        f.close()
        
    #HACK
    w_order_RVs = w_order_RVs_barycorr
    for j in range (len(w_dates)):
        w_order_RV_scatter[j] = np.nanstd(w_order_RVs[:,j])
        
    #read in SERVAL
    ser_rvc = np.loadtxt(servaldir+i[1]+"/"+i[1]+".rvc.dat") # rvc
    if use_avcn == True:
        ser_avcn = np.loadtxt(servaldir+i[1]+"/"+i[1]+".avcn.dat")
        ser_rvc = ser_avcn#HACK
    # remove entries with nan in drift
    ind_finitedrift = np.isfinite(ser_rvc[:,3])
    ser_rvc = ser_rvc[ind_finitedrift]
    #read in also info file, from which we get SNR and airmass
    ser_info = np.genfromtxt(servaldir+i[1]+"/"+i[1]+".info.cvs", delimiter=";")
    # now go throughall ser_rvc dates, find the lowest difference in jd, and store SNR and airmass
    ser_addinfo=np.zeros((len(ser_rvc), 2)) #initializes array
    for n in range(len(ser_rvc)):
            ind_jd = np.where(np.abs(ser_info[:,1]-ser_rvc[n,0]) == np.nanmin(np.abs(ser_info[:,1]-ser_rvc[n,0])))[0][0] #ser_info [:,1] is BJD and so is ser_rvc[n,0] this matches the ones closest to each other
            ser_addinfo[n,0] = ser_info[ind_jd, 3] # SNR
            ser_addinfo[n,1] = ser_info[ind_jd, 8] # Airmass
            
    #Import serval orderwise Rvs
    ser_rvo = np.loadtxt(servaldir+i[1]+"/"+i[1]+".rvo.dat")
    ser_rvo = ser_rvo[ind_finitedrift]
    
    ser_rvo_err = np.loadtxt(servaldir+i[1]+"/"+i[1]+".rvo.daterr")
    ser_rvo_err = ser_rvo_err[ind_finitedrift]
    
    #correction list to apply serval correctios to Wobble
    #add to wobble to be comparable to ser_rvc (clearly not addapted to orderwise yet)
    # be  careful to subtract out mean correction for normalized data # NOTE this si still not a perfect method
    ser_corr = - ser_rvc[:,8] - ser_rvc[:,3]
    if use_avcn == True:
        #8 sa drift, 3 drift, 9 NZP
        ser_corr = - ser_rvc[:,8] - ser_rvc[:,3] - ser_rvc[:,9] # use with avcn
    
    ### plots ###
    
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
    #apply corretion to all wobble results and serval orders
    ser_rvo = ser_rvo[indices_serval]     
    #HACK to correct ser_rvo
    for j in range(len(ser_rvo)):
        ser_rvo[j,5:] = ser_rvo[j,5:] + ser_corr[j]
        
    ser_rvo_err = ser_rvo_err[indices_serval]
    ser_addinfo = ser_addinfo[indices_serval]
    w_dates = w_dates[indices_wobble]
    w_dates_utc = w_dates_utc[indices_wobble]
    w_epochs = w_epochs[indices_wobble]
    w_RVs_barycorr = w_RVs_barycorr[indices_wobble]                         + ser_corr
    w_RVs_barycorr_ivar = w_RVs_barycorr_ivar[indices_wobble]               + ser_corr
    w_RVs = w_RVs[indices_wobble]                                           + ser_corr
    w_RVs_er = w_RVs_er[indices_wobble]
    w_RVs_own = w_RVs_own[indices_wobble]                                   + ser_corr
    w_RVs_unweighted = w_RVs_unweighted[indices_wobble]                     + ser_corr
    w_bervs = w_bervs[indices_wobble]
    w_order_RVs = w_order_RVs[:,indices_wobble]                             + ser_corr
    w_order_RV_scatter = w_order_RV_scatter[indices_wobble]
    w_RVs_original = w_RVs_original[indices_wobble]                         + ser_corr
    w_RVs_original_barycorr = w_RVs_original_barycorr[indices_wobble]       + ser_corr
    
    
    #simple plot of RVs (only matched)
    ax1.plot(w_dates, w_RVs_barycorr - np.nanmean(w_RVs_barycorr),"x", label="Wobble_bary")
    ax1.plot(ser_rvc[:,0], ser_rvc[:,5] - np.nanmean(ser_rvc[:,5]), "x", label= "SERVAL")
    ax1.plot(ser_rvc[:,0], ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]), "x", label= "SERVAL_Corr")
    ax1.set_ylabel("RVs [m/s]")
    ax1.set_xlabel('jd')
    plt.title("matched RVs for "+i[0]+" ("+i[2]+") "+" - "+i[1]+"; (RVs only BERV-corrected; no FP-drift, no SA, no NZP corrections)")
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    #Check Serval correction
    mean_corr = np.nanmean(ser_corr)
    #ax1.plot(w_dates, w_RVs_barycorr - np.nanmean(w_RVs_barycorr),"x", label="Wobble_bary")
    ax1.plot(ser_rvc[:,0], ser_rvc[:,5] + ser_corr - ser_rvc[:,1], "x", label= "SERVAL_manual -Scorr meancorr = "+ str(mean_corr))
    ax1.set_ylabel("RVs [m/s]")
    ax1.set_xlabel('jd')
    plt.title("check corrections")
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
# =============================================================================
#     #simple plot of matched RVs Serval corrected
#     ax1.plot(w_dates, w_RVs_barycorr - np.nanmean(w_RVs_barycorr),"x", label="Wobble")
#     ax1.plot(ser_rvc[:,0], ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]), "x", label= "SERVAL")
#     ax1.set_ylabel("RVs [m/s]")
#     ax1.set_xlabel('jd')
#     plt.title("matched RVs for "+i[0]+" ("+i[2]+") "+" - "+i[1]+"; Serval Corrected")
#     plt.grid(True)
#     plt.tight_layout()
#     plt.legend(shadow=True)
#     plt.savefig(pp, format='pdf')
#     plt.clf()
#     fig.clf()
#     ax1 = plt.gca()
# =============================================================================
    
    #simple barrycentric correction
    ax1.plot(w_dates, w_bervs ,"x", label="Wobble")
    ax1.plot(w_dates, w_RVs_barycorr - w_RVs_unweighted  ,"x", label="w_RVs_barycorr - w_RVs_unweighted")
    
    ax1.set_ylabel("w_bervs [m/s]")
    ax1.set_xlabel('jd')
    plt.title("w_bervs for "+i[0]+" ("+i[2]+") "+" - "+i[1]+"")
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    
    
    
     # plot also the scatter of all orders vs. JD
    #for o in range(N_ord):
    #    ax1.plot(w_dates, w_order_RVs[o,:]-w_RVs, "x", label="o"+str(o))
   
    for ep in range(len(w_dates)):
        ax1.plot(np.zeros(N_ord) + w_dates[ep], w_order_RVs[:,ep] - np.nanmean(w_RVs_original_barycorr), "x")
    
    
    
    wob_rv_scatter = np.median(w_order_RV_scatter)
    ax1.plot(w_dates, w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr), "ko",label = 'Wobble; scatter_median = ' + str(wob_rv_scatter))
    #ax1.plot(w_dates, w_RVs_unweighted, "ko")
    ax1.set_ylabel('RVs per order and combined [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV scatter of Wobble orders "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    orders = w_orders + 5
    for ep in range(len(w_dates)):
        ax1.plot(np.zeros(N_ord) + w_dates[ep], ser_rvo[ep, orders] - np.nanmean(ser_rvc[:,5]), "x")
    
    ser_order_scatter = [np.nanstd(ser_rvo[ep, orders]) for ep in range(len(w_dates))]
    ser_rv_scatter = np.median(ser_order_scatter)
    ax1.plot(w_dates, ser_rvc[:,5] - np.nanmean(ser_rvc[:,5]), "ko", label = 'Serval; scatter_median = ' + str(ser_rv_scatter) )
    #ax1.plot(w_dates, w_RVs_unweighted, "ko")
    ax1.set_ylabel('RVs per order and combined [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV scatter of SERVAL orders "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    #Plot comparing wob and serval measurements with errorbars showing order scatter
    ###repeats
    wob_rv_scatter = np.median(w_order_RV_scatter)
    ser_order_scatter = [np.nanstd(ser_rvo[ep, orders]) for ep in range(len(w_dates))]
    ser_rv_scatter = np.median(ser_order_scatter)
    ###
    ax1.errorbar(
        [i for i in range(len(w_dates))],
        w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr) ,
        yerr = w_order_RV_scatter,
        ls='none', marker = "x",label = 'Wobble; scatter_median = ' + str(wob_rv_scatter))
    
    ax1.errorbar(
        [i+ 0.5 for i in range(len(w_dates))],
        ser_rvc[:,5] - np.nanmean(ser_rvc[:,5]) ,
        yerr = ser_order_scatter,
        ls='none', marker = "x", label = 'Serval; scatter_median = ' + str(ser_rv_scatter))
    
    ax1.set_ylabel('RVs and scatter [m/s]')
    ax1.set_xlabel('epoch')
    plt.title("RV scatter of Wobble orders "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
     #Same with capped y axis
    ###repeats
    wob_rv_scatter = np.median(w_order_RV_scatter)
    ser_order_scatter = [np.nanstd(ser_rvo[ep, orders]) for ep in range(len(w_dates))]
    ser_rv_scatter = np.median(ser_order_scatter)
    ###
    ax1.set_ylim(-50,50)
    ax1.errorbar(
        [i for i in range(len(w_dates))],
        w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr),
        yerr = w_order_RV_scatter,
        ls='none', marker = "x",label = 'Wobble; scatter_median = ' + str(wob_rv_scatter))
    
    ax1.errorbar(
        [i+ 0.5 for i in range(len(w_dates))],
        ser_rvc[:,5] - np.nanmean(ser_rvc[:,5]) ,
        yerr = ser_order_scatter,
        ls='none', marker = "x", label = 'Serval; scatter_median = ' + str(ser_rv_scatter))
    
    ax1.set_ylabel('RVs and scatter [m/s]')
    ax1.set_xlabel('epoch')
    plt.title("RV scatter of Wobble orders "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    #RV differences serval-wobble
    # seems to be correlted with barycentric correction 
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_barycorr -np.nanmean(w_RVs_barycorr)), "x",label= "w_median_bary")
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_barycorr_ivar-np.nanmean(w_RVs_barycorr_ivar)), "x", label= "w_ivar_weighted_bary")
    ax1.set_ylabel('delta_RVs [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV differences (SERVAL-Wobble) for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    
    #using simplified additive barycorr -> differences ofup to ca 5m/s ~ 1%
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs+w_bervs-np.nanmean(w_RVs+w_bervs)), "x", label="w_median")
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_original+w_bervs-np.nanmean(w_RVs_original+w_bervs)), "x" , label="w_originals") #original wobble combine yields much worse scatter
    #ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_barycorr -np.nanmean(w_RVs_barycorr)), "x",label= "w_median_bary")
    #ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_barycorr_ivar-np.nanmean(w_RVs_barycorr_ivar)), "x", label= "w_ivar_weighted_bary")
    ax1.set_ylabel('delta_RVs [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV differences additive BERV (SERVAL-Wobble) for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
     #using simplified additive barycorr -> differences ofup to ca 5m/s ~ 1%
    ax1.plot(ser_rvc[:,0], w_RVs+w_bervs- w_RVs_barycorr , "x", label="w_median-w_median_bary")
    ax1.plot(ser_rvc[:,0], w_RVs_original+w_bervs - w_RVs_original_barycorr , "x", label="w_orig-w_orig_bary")
    #ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_barycorr -np.nanmean(w_RVs_barycorr)), "x",label= "w_median_bary")
    #ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_barycorr_ivar-np.nanmean(w_RVs_barycorr_ivar)), "x", label= "w_ivar_weighted_bary")
    ax1.set_ylabel('delta_RVs [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV differences Barycor vs. additive BERV (Wobble-Wobble) for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    


    #flat until 10km/s then anticorrelated (fairly significantly)
    ax1.plot(ser_rvc[:,7], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_barycorr-np.nanmean(w_RVs_barycorr)), "x")
    ax1.set_ylabel('delta_RVs [m/s]')
    ax1.set_xlabel('BERV [km/s]')
    plt.title("Correlation of RV diff (SERVAL-Wobble) with BERV for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    #plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
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
        if 'orbital_parameters' in locals():
            def fit_func(t, T0_offset):
                return rv.radial_velocity(t , orbital_parameters[0], orbital_parameters[1], orbital_parameters[2],orbital_parameters[3], orbital_parameters[4] + T0_offset)
            
            def keplarian_rv(t):
                return rv.radial_velocity_M0(t , orbital_parameters[0], orbital_parameters[1], orbital_parameters[2],orbital_parameters[3], orbital_parameters[4], orbital_parameters[5])
            
            if 'orbital_parameters_mult' in locals():
                def keplarian_rv_mult(t):
                    total_rv = 0
                    for parameters in orbital_parameters_mult:
                        total_rv = total_rv + rv.radial_velocity_M0(t , parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
                    return total_rv
                
                def fit_func_mult(t, T0_offset):
                    return keplarian_rv_mult(t + T0_offset)
                        
                #NOTE this is not a very sensible location for this definiton        
                def error_stats(residuals, errors, kappa = 5):
                    residual_offset = np.average(residuals)
                    # note standart deviation autopmatically corrects for the offset
                    sigma = np.nanstd(residuals)
                    sigma_clipped = np.nanstd(sigma_clip(residuals,sigma = 5) )
                    ind_bad = np.where(np.absolute(residuals - residual_offset) > kappa * errors)[0]
                    ind_good= np.where(np.absolute(residuals - residual_offset) <= kappa * errors)[0]
                    sigma_error_clipped = np.nanstd(residuals[ind_good])
                    return sigma, sigma_clipped,sigma_error_clipped, ind_bad, ind_good, residual_offset
                        
                #####################        
                #Plot all planet fit with data
                #NOTE There seems to be a systematic deviation from 0 mean in the residuals -> likely due to bad normalization with the mean of the data beforehand
                wob_residuals_sign_mult = w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr) - keplarian_rv_mult(w_dates)
                wob_mult_stats = error_stats(wob_residuals_sign_mult, w_RVs_er)
                residual_offset = wob_mult_stats[5] # Extremely similar offsets for both wob and ser)
                
                #print(np.average(wob_residuals_sign_mult))
                #print( np.nanmean(w_RVs_original_barycorr))
                
                ser_residuals_sign_mult = ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) - keplarian_rv_mult(ser_rvc[:,0])
                ser_mult_stats = error_stats(ser_residuals_sign_mult, ser_rvc[:,2])
                #print(np.average(ser_residuals_sign_mult))
                
                w_dates_sorted = np.sort(w_dates)
                xlst = np.linspace(w_dates_sorted[0], w_dates_sorted[-1], num = 500)
                ax1.plot(xlst , keplarian_rv_mult(xlst), "k-", label = "lit orbit M0 (no fit)")
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
                plt.grid(True)
                plt.tight_layout()
                plt.legend(shadow=True)
                plt.savefig(pp, format='pdf')
                plt.clf()
                fig.clf()
                ax1 = plt.gca()
                
                ax1.errorbar(w_dates, wob_residuals_sign_mult - residual_offset, yerr = w_RVs_er,fmt = "x", label = "wobble_corr")
                ax1.errorbar((ser_rvc[:,0]), ser_residuals_sign_mult - residual_offset, yerr = ser_rvc[:,2] ,fmt = "x",color = "C2", label = "SERVAL_corr")
                
                if error_clipping == True:
                    ax1.plot((w_dates[ind_bad_wob]), wob_residuals_sign_mult[ind_bad_wob] - residual_offset, "o",color = "C0", label = "error_clipped wobble points")
                
                ax1.set_ylabel("RVs [m/s]")
                ax1.set_xlabel('jd')
                plt.title("RVs Residuals multiplanet fit for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
                plt.grid(True)
                plt.tight_layout()
                plt.legend(shadow=True)
                plt.savefig(pp, format='pdf')
                plt.clf()
                fig.clf()
                ax1 = plt.gca()
                
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
                xlst = np.linspace(w_dates_sorted[0], w_dates_sorted[-1], num = 500)
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
                plt.grid(True)
                plt.tight_layout()
                plt.legend(shadow=True)
                plt.savefig(pp, format='pdf')
                plt.clf()
                fig.clf()
                ax1 = plt.gca()
                
                ax1.errorbar(w_dates, wob_residuals_sign_mult - residual_offset, yerr = w_RVs_er,fmt = "x", label = "wobble_corr")
                ax1.errorbar((ser_rvc[:,0]), ser_residuals_sign_mult - residual_offset, yerr = ser_rvc[:,2] ,fmt = "x",color = "C2", label = "SERVAL_corr")
                
                if error_clipping == True:
                    ax1.plot((w_dates[ind_bad_wob]), wob_residuals_sign_mult[ind_bad_wob] - residual_offset, "o",color = "C0", label = "error_clipped wobble points")
                
                ax1.set_ylabel("RVs [m/s]")
                ax1.set_xlabel('jd')
                plt.title("RVs Residuals multiplanet fit for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
                plt.grid(True)
                plt.tight_layout()
                plt.legend(shadow=True)
                plt.savefig(pp, format='pdf')
                plt.clf()
                fig.clf()
                ax1 = plt.gca()
            
                ######################
                
            ###########################     
            #fit to Wobble
            xdata = w_dates
            ydata = w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)
            popt, pcov = sp.optimize.curve_fit(fit_func, xdata, ydata,  sigma = w_RVs_er, absolute_sigma = True)
            print("T0_offset Wobble = ", popt)
            T0_offset = popt[0]
            
            #make these weighted (maybe: thsi may not be a good idea if residuals are not strongly correlated to error (as with wobble results))
            wob_residuals_sign = w_RVs_original_barycorr - np.nanmean(w_RVs_original_barycorr) - keplarian_rv(w_dates)
            wob_residuals = np.absolute(wob_residuals_sign)
            
            sigma_wob = np.nanstd(sigma_clip(wob_residuals_sign,sigma = 5) )
            sigma_wob_no_clip = np.nanstd(wob_residuals_sign)
            
            #residual vs error clipping
            error_clipping_kappa = 5
            ind_bad_wob = np.where(wob_residuals > error_clipping_kappa * w_RVs_er)[0]
            ind_good_wob = np.where(wob_residuals <= error_clipping_kappa * w_RVs_er)[0]
            
            sigma_wob_error_clipped = np.nanstd(wob_residuals_sign[ind_good_wob]) 
            
            #fit to serval:
            xdata = ser_rvc[:,0]
            ydata = ser_rvc[:,1] - np.nanmean(ser_rvc[:,1])
            popt_s, pcov_s = sp.optimize.curve_fit(fit_func, xdata, ydata, sigma = ser_rvc[:,2], absolute_sigma = True, p0 = T0_offset) # NOTE HACK p0 choice based on wobble can fail
            print("T0_offset Serval = ", popt_s)
            T0_offset_s = popt_s[0]
            
            ser_residuals_sign = ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) - keplarian_rv(ser_rvc[:,0])
            ser_residuals = np.absolute(ser_residuals_sign)
            
            sigma_ser = np.nanstd(sigma_clip(
            ser_residuals_sign
            ,sigma = 5) )
            sigma_ser_no_clip = np.nanstd(
            ser_residuals_sign)
            
            #residual vs error clipping
            ind_bad_ser = np.where(ser_residuals > error_clipping_kappa * ser_rvc[:,2])[0]
            ind_good_ser = np.where(ser_residuals <= error_clipping_kappa * ser_rvc[:,2])[0]
            
            sigma_ser_error_clipped = np.nanstd(ser_residuals_sign[ind_good_ser]) 
            
            xlst = np.linspace(w_dates[0], w_dates[0] + orbital_parameters[1]*0.99999, num=100)
            ylst = [rv.radial_velocity(t , orbital_parameters[0], orbital_parameters[1], orbital_parameters[2],orbital_parameters[3], orbital_parameters[4] + T0_offset) for t in xlst]
            #sort by xlst
            pltlst = [[xlst[j],ylst[j]] for j in range(len(xlst))]
            def mod_sort(elem):
                return elem[0] % orbital_parameters[1]
            pltlst = sorted(pltlst, key = mod_sort)
            pltlst = np.asarray(pltlst)
            pltlst = [pltlst[:,0],pltlst[:,1]]
            #to_print = pltlst[0] % i[3]
            #print(to_print)
            
            ax1.plot(pltlst[0] % i[3], pltlst[1], "r-", label = "literature orbit (Wobble T0_offset) single planet")
            ax1.plot(pltlst[0] % i[3], keplarian_rv(pltlst[0]), "k-", label = "lit orbit M0 (no fit) single planet")
                     
        #ax1.plot(w_dates % i[3], (w_RVs_original+w_bervs-np.nanmean(w_RVs_original+w_bervs)),"x", label="Wobble_orig")
        ax1.errorbar((w_dates) % i[3], (w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)), yerr = w_RVs_er,fmt = "x", label="Wobble_Corr, no_clip ={0:.3f},clipped_sigma = {1:.3f}, err_clipped = {2:.3f} ".format(sigma_wob_no_clip, sigma_wob, sigma_wob_error_clipped))
        ax1.plot((ser_rvc[:,0]) % i[3], ser_rvc[:,5] - np.nanmean(ser_rvc[:,5]), "x", label= "SERVAL")
        ax1.errorbar((ser_rvc[:,0]) % i[3], ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]),yerr = ser_rvc[:,2] ,fmt = "x", label= "SERVAL_Corr, no_clip ={0:.3f}, clipped_sigma = {1:.3f}, err_clipped = {2:.3f}".format(sigma_ser_no_clip, sigma_ser, sigma_ser_error_clipped))
        if error_clipping == True:
            ax1.plot((w_dates[ind_bad_wob]) % i[3], (w_RVs_original_barycorr[ind_bad_wob] - np.nanmean(w_RVs_original_barycorr)), "o",color = "C0", label = "error_clipped wobble points")
            ax1.plot((ser_rvc[:,0][ind_bad_ser]) % i[3], ser_rvc[:,1][ind_bad_ser] - np.nanmean(ser_rvc[:,1]), "o", color = "C2", label= "error_clipped serval points")
        #ax1.plot(w_dates % i[3], (w_RVs_original_barycorr-np.nanmean(w_RVs_original_barycorr)),"x", label="Wobble_Corr")
        ax1.set_ylabel("RVs [m/s]")
        ax1.set_xlabel('jd')
        plt.title("Phased ("+str(i[3])+"d) RVs for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        ##########################################
    
    #Plot residuals. Resuduals are relative to T0_offset fit to each dataset individually
    if 'orbital_parameters' in locals():
        ax1.errorbar((w_dates) % i[3], wob_residuals_sign, yerr = w_RVs_er,fmt = "x", label = "wobble_corr")
        ax1.errorbar((ser_rvc[:,0]) % i[3], ser_residuals_sign, yerr = ser_rvc[:,2] ,fmt = "x",color = "C2", label = "SERVAL_corr")
        
        if error_clipping == True:
            ax1.plot((w_dates[ind_bad_wob]) % i[3], wob_residuals_sign[ind_bad_wob], "o",color = "C0", label = "error_clipped wobble points")
        
        ax1.set_ylabel("RVs [m/s]")
        ax1.set_xlabel('jd')
        plt.title("RVs Residuals for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        
        ax1.errorbar((w_dates) % i[3], wob_residuals_sign / w_RVs_er, yerr = w_RVs_er/w_RVs_er ,fmt = "x", label = "wobble_corr")
        ax1.errorbar((ser_rvc[:,0]) % i[3], ser_residuals_sign/ser_rvc[:,2] ,yerr = ser_rvc[:,2]/ser_rvc[:,2]  ,fmt = "x",color = "C2", label = "SERVAL_corr")
        
        if error_clipping == True:
            ax1.plot((w_dates[ind_bad_wob]) % i[3], (wob_residuals_sign / w_RVs_er)[ind_bad_wob], "o",color = "C0", label = "error_clipped wobble points")
        
        ax1.set_ylabel("n*error residuals")
        ax1.set_xlabel('jd')
        plt.title("RVs Residuals as multiples of error for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        #Same Plot with cut x axis
        ax1.errorbar((w_dates) % i[3], wob_residuals_sign / w_RVs_er, yerr = 0*w_RVs_er/w_RVs_er ,fmt = "x", label = "wobble_corr")
        ax1.errorbar((ser_rvc[:,0]) % i[3], ser_residuals_sign/ser_rvc[:,2] ,yerr = 0*ser_rvc[:,2]/ser_rvc[:,2]  ,fmt = "x",color = "C2", label = "SERVAL_corr")
        if error_clipping == True:
            ax1.plot((w_dates[ind_bad_wob]) % i[3], (wob_residuals_sign / w_RVs_er)[ind_bad_wob], "o",color = "C0", label = "error_clipped wobble points")
        
        ax1.set_ylabel("n*error residuals")
        ax1.set_xlabel('jd')
        ax1.set_ylim(-10,+10)
        plt.title("RVs Residuals as multiples of error for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        
         #make list of bad epochs
        bad_epochs = w_epochs[np.where(wob_residuals > 10)[0]]
        bad_dates = w_dates[np.where(wob_residuals > 10)[0]]
        very_bad_dates = w_dates[np.where(wob_residuals > 20)[0]]
        print("bad epochs: ", bad_epochs)
        print("bad dates: ", bad_dates)
        print("very bad dates: ", very_bad_dates)
        
        ax1.plot(wob_residuals, w_RVs_er,"x", label = "wobble_corr")
        ax1.plot(ser_residuals, ser_rvc[:,2] , "x",color = "C2", label = "SERVAL_corr")
        if error_clipping == True:
            ax1.plot((wob_residuals[ind_bad_wob]), w_RVs_er[ind_bad_wob], "o",color = "C0", label = "error_clipped wobble points")
        
        ax1.set_ylabel("RV error [m/s]")
        ax1.set_xlabel('RV residuals [m/s]')
        plt.title("RVs error correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        ax1.plot(wob_residuals, ser_addinfo[:,0],"x", label = "wobble_corr")
        ax1.set_ylabel("SNR")
        ax1.set_xlabel('RV residuals [m/s]')
        plt.title("RV residuals vs SNR correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        ax1.plot(wob_residuals, ser_addinfo[:,1],"x", label = "wobble_corr")
        ax1.set_ylabel("Airmass")
        ax1.set_xlabel('RV residuals [m/s]')
        plt.title("RV residuals vs Airmass correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        ax1.plot(wob_residuals_sign, ser_corr,"x", label = "wobble_corr")
        if error_clipping == True:
            ax1.plot((wob_residuals_sign[ind_bad_wob]), ser_corr[ind_bad_wob], "o",color = "C0", label = "error_clipped wobble points")
        ax1.set_ylabel("Ser_corr [m/s]")
        ax1.set_xlabel('RV residuals [m/s]')
        plt.title("RV residuals vs Correction correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
    
        ax1.plot(wob_residuals_sign/ w_RVs_er, ser_corr,"x", label = "wobble_corr")
        ax1.set_ylabel("Ser_corr [m/s]")
        ax1.set_xlabel(r'RV residuals [$\sigma$]')
        plt.title("RV residuals vs Correction correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        ax1.plot(wob_residuals_sign/ w_RVs_er, ser_corr/w_RVs_er,"x", label = "wobble_corr")
        ax1.plot([1.1*float(x) for x in range(int(min(wob_residuals_sign/ w_RVs_er)),int(max(wob_residuals_sign/ w_RVs_er)))],[1.1*float(x) for x in range(int(min(wob_residuals_sign/ w_RVs_er)),int(max(wob_residuals_sign/ w_RVs_er)))], "-", label = "slope 1", color="r")
        if error_clipping == True:
            ax1.plot(((wob_residuals_sign/ w_RVs_er)[ind_bad_wob]), (ser_corr/w_RVs_er)[ind_bad_wob], "o",color = "C0", label = "error_clipped wobble points")
        ax1.set_ylabel(r"Ser_corr [$\sigma$]")
        ax1.set_xlabel(r'RV residuals [$\sigma$]')
        plt.title("RV residuals vs Correction correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        if use_avcn == True:
            #8 sa drift, 3 drift, 9 NZP
            #ser_corr = - ser_rvc[:,8] - ser_rvc[:,3] - ser_rvc[:,9] # use with avcn
            ax1.plot(wob_residuals_sign/ w_RVs_er, - ser_rvc[:,8], "x", label = "wobble_corr")
            ax1.set_ylabel("SA drift correction [m/s]")
            ax1.set_xlabel(r'RV residuals [$\sigma$]')
            plt.title("RV residuals vs SA drift correction correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
            plt.grid(True)
            plt.tight_layout()
            plt.legend(shadow=True)
            plt.savefig(pp, format='pdf')
            plt.clf()
            fig.clf()
            ax1 = plt.gca()
            
            ax1.plot(wob_residuals_sign/ w_RVs_er, - ser_rvc[:,3], "x", label = "wobble_corr")
            ax1.set_ylabel("drift correction [m/s]")
            ax1.set_xlabel(r'RV residuals [$\sigma$]')
            plt.title("RV residuals vs drift correction correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
            plt.grid(True)
            plt.tight_layout()
            plt.legend(shadow=True)
            plt.savefig(pp, format='pdf')
            plt.clf()
            fig.clf()
            ax1 = plt.gca() 
            
            ax1.plot(wob_residuals_sign/ w_RVs_er, - ser_rvc[:,9], "x", label = "wobble_corr")
            ax1.set_ylabel("NZP correction [m/s]")
            ax1.set_xlabel(r'RV residuals [$\sigma$]')
            plt.title("RV residuals vs NZP correction correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
            plt.grid(True)
            plt.tight_layout()
            plt.legend(shadow=True)
            plt.savefig(pp, format='pdf')
            plt.clf()
            fig.clf()
            ax1 = plt.gca() 
        
    # Just serval corrected - wobble_orig with barrycorr (tight view)
    #sigma clip before calculating std
    sigma_difference = np.nanstd(sigma_clip(
        (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_original_barycorr -np.nanmean(w_RVs_original_barycorr))
        ,sigma = 5) )
    sigma_difference_no_clip = np.nanstd(
        (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_original_barycorr -np.nanmean(w_RVs_original_barycorr)))
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_original_barycorr -np.nanmean(w_RVs_original_barycorr)) , "x", label="w_orig_bary, std of difference = {0:.3f}, clipped = {1:.3f} ".format(sigma_difference_no_clip, sigma_difference) )
    ax1.errorbar(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_original_barycorr -np.nanmean(w_RVs_original_barycorr)) , yerr = w_RVs_er, ls='none')
    ax1.set_ylim(bottom = -8 , top = 8)
    ax1.set_ylabel('delta_RVs [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV differences Serval - Wobble (y- windowed) for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca() 
    
    #print(w_orders)
    #last order (41, 55) makes issues here 
    good_orders = []
    for index, order in enumerate(w_orders[:-2]):
        #print(str(index), str(order))
        order_ind = index
        ind_rvo = order + 5
        sigma_difference = np.nanstd(
        (ser_rvo[:, ind_rvo] - np.nanmean(ser_rvo[:, ind_rvo])) - (w_order_RVs[order_ind ,:] -np.nanmean(w_order_RVs[order_ind ,:]))
         )
        sigma_difference_clipped = np.nanstd(sigma_clip(
        (ser_rvo[:, ind_rvo] - np.nanmean(ser_rvo[:, ind_rvo])) - (w_order_RVs[order_ind ,:]-np.nanmean(w_order_RVs[order_ind ,:]))
        ,sigma = 5) )
        serval_order_error = np.nanmedian(ser_rvo_err[:, ind_rvo])
        #Comparison with Serval rvc as "objective" standart
        sigma_difference_rvc = np.nanstd(
        (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_order_RVs[order_ind ,:] -np.nanmean(w_order_RVs[order_ind ,:]))
         )
        sigma_difference_clipped_rvc = np.nanstd(sigma_clip(
        (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_order_RVs[order_ind ,:]-np.nanmean(w_order_RVs[order_ind ,:]))
        ,sigma = 5) )
        #Compare Serval order vs Serval
        sigma_difference_rvc_S = np.nanstd(
        (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (ser_rvo[:, ind_rvo] - np.nanmean(ser_rvo[:, ind_rvo])) 
         )
        sigma_difference_clipped_rvc_S = np.nanstd(sigma_clip(
        (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (ser_rvo[:, ind_rvo] - np.nanmean(ser_rvo[:, ind_rvo])) 
        ,sigma = 5) )
        
        #Make list of "Good" orders with e.g ks < 10
        if sigma_difference_clipped < 5:
            good_orders.append(order_ind)
        
        ax1.plot(w_dates[:], (ser_rvo[:, ind_rvo] - np.nanmean(ser_rvo[:, ind_rvo])) - (w_order_RVs[order_ind ,:]  -np.nanmean(w_order_RVs[order_ind ,:]) ), "x", label= r"std_of_difference = {0:.3f}, $\kappa \sigma$ clipped = {1:.3f} ".format(sigma_difference, sigma_difference_clipped))
        ax1.plot(w_dates[:], (ser_rvo[:, ind_rvo] - np.nanmean(ser_rvo[:, ind_rvo])), "x", label= "Serval", alpha = 0.5)
        ax1.plot(w_dates[:], (w_order_RVs[order_ind ,:]  -np.nanmean(w_order_RVs[order_ind ,:])), "x", label= "wobble", alpha = 0.5)
        
        ax1.plot([], [], ' ', label="Serval order median error = {0:.3f}".format(serval_order_error))
        ax1.plot([], [], ' ', label= r"w_o vs S $\sigma ( \Delta)$ = {0:.3f}, $\kappa \sigma$ clipped = {1:.3f} ".format(sigma_difference_rvc, sigma_difference_clipped_rvc))
        ax1.plot([], [], ' ', label= r"S_o vs S $\sigma ( \Delta)$ = {0:.3f}, $\kappa \sigma$ clipped = {1:.3f} ".format(sigma_difference_rvc_S, sigma_difference_clipped_rvc_S))
        
        ax1.set_ylabel('RVs per order and combined [m/s]')
        ax1.set_xlabel('jd')
        plt.title("RV order " + str(order)+ "(i=" + str(order_ind) + ")" +" s-w (quick Bervs)" +i[0]+" ("+i[2]+") "+" - "+i[1])
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
    print("good_orders: ", w_orders[good_orders])
    
    if len(good_orders) > 0:
        #NOTE
        #RV scatter good orders with capped y axis
        #Use median of good order only
        w_order_RVs_good = np.zeros((len(good_orders),len(w_dates)))
        for ep in range (len(w_dates)):
            for index, order in enumerate(good_orders):
                w_order_RVs_good[index, ep] = w_order_RVs[order, ep]
        
        
        w_RVs_median_good = np.zeros(len(w_dates))  
        w_order_RVs_good_std = np.zeros(len(w_dates))
        for ep in range(len(w_dates)):
            w_RVs_median_good[ep] = np.nanmedian(w_order_RVs_good[:,ep])
            w_order_RVs_good_std[ep] = np.nanstd(w_order_RVs_good[:,ep])
        #NOTE
        
        ###repeats
        wob_rv_scatter = np.median(w_order_RVs_good_std)
        ser_order_scatter = [np.nanstd(ser_rvo[ep, w_orders[good_orders]]) for ep in range(len(w_dates))]
        ser_rv_scatter = np.median(ser_order_scatter)
        ###
        ax1.set_ylim(-50,50)
        ax1.errorbar(
            [i for i in range(len(w_dates))],
            w_RVs_median_good - np.nanmean(w_RVs_median_good),
            yerr = w_order_RVs_good_std,
            ls='none', marker = "x",label = 'Wobble; scatter_median = ' + str(wob_rv_scatter))
        
        ax1.errorbar(
            [i+ 0.5 for i in range(len(w_dates))],
            ser_rvc[:,1] - np.nanmean(ser_rvc[:,1]) ,
            yerr = ser_order_scatter,
            ls='none', marker = "x", label = 'Serval; scatter_median = ' + str(ser_rv_scatter))
        
        ax1.set_ylabel('RVs and scatter [m/s]')
        ax1.set_xlabel('epoch')
        plt.title("RV scatter good orders "+i[0]+" ("+i[2]+") "+" - "+i[1])
        plt.grid(True)
        plt.tight_layout()
        plt.legend(shadow=True)
        plt.savefig(pp, format='pdf')
        plt.clf()
        fig.clf()
        ax1 = plt.gca()
        
        # Just serval uncorrected - wobble_orig with barrycorr (tight view)
        #sigma clip before calculating std
        sigma_difference = np.nanstd(sigma_clip(
            (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_median_good -np.nanmean(w_RVs_median_good))
            ,sigma = 5) )
        ax1.plot(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_median_good -np.nanmean(w_RVs_median_good)) , "x", label="w_median_good, std of difference (clipped)= {0:.3f}".format(sigma_difference) )
        ax1.errorbar(ser_rvc[:,0], (ser_rvc[:,1] -np.nanmean(ser_rvc[:,1])) - (w_RVs_median_good -np.nanmean(w_RVs_median_good)) , yerr = w_order_RVs_good_std, ls='none')
        ax1.set_ylim(bottom = -8 , top = 8)
        ax1.set_ylabel('delta_RVs [m/s]')
        ax1.set_xlabel('jd')
        plt.title("RV differences Serval - Wobble (y- windowed) for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca() 
    

plt.close(fig)
pp.close()
    



