#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 11:36:17 2019

@author: cmatthe
"""
#%% imports
import matplotlib.pyplot as plt #for plt
from matplotlib.backends.backend_pdf import PdfPages # for PdfPages: suspiciously abscent in akaminskis code
#I suspect there is a hidden standart import for example also for numpy?
# perhaps in  sys.path.append('/home/akaminsk/pythontests')
import matplotlib as mpl #missing in akaminski
import h5py 
import numpy as np
from tqdm import tqdm
import barycorrpy as bary 
from scipy.constants import codata 
lightvel = codata.value('speed of light in vacuum') #for barycorr

import pickle # for writing barycorr results into a binary file



#%%

# CAHA Coordinates for barycorr
_lat = 37.2236
_lon = -2.54625
_elevation = 2168.


pp =PdfPages("wobble_servalcomparison_reconstruct.pdf")
fig = plt.figure(figsize=(15, 9), dpi=200)
mpl.rc('font', size=16)
plt.clf()
fig.clf()
ax1=plt.gca()

#read in wobble and serval and compare RVs from Wobble (barycentric corrected) with SERVAL RVs, but also with RVCs and AVCs 

objects = [["GJ876", "J22532-142", "vis"]]

for i in objects:
    print(i[0], i[1], i[2])
    #read in wobble result
    if i[2]=="vis":
        #servaldir = "/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
        #to run from lx39, perhas replace this with  an if statement
        servaldir_lx39="/data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
        servaldir=servaldir_lx39
        #wobble_res = h5py.File("/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_Kstar0_Kt3.hdf5", 'r')
        wobble_res_lx39 = h5py.File("/data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_Kstar0_Kt3.hdf5")
        wobble_res=wobble_res_lx39
        N_ord=42 # note CARMENES has 61 orderes only 11-52 used here
        k=np.arange(N_ord)
        keys_orderRVs=[]
        for j in k:
            keys_orderRVs.append('order'+str(j))
            
    else:
        #servaldir = "/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/servaldir/CARM_NIR/"
        #to run from lx39
        servaldir_lx39="/data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
        servaldir=servaldir_lx39
        #wobble_res = h5py.File("/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_nir_Kstar0_Kt3.hdf5")
        wobble_res_lx39 = h5py.File("/data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_nir_Kstar0_Kt3.hdf5")
        wobble_res=wobble_res_lx39
        N_ord=17
        k=np.arange(N_ord)
        keys_orderRVs=[]
        for j in k:
            keys_orderRVs.append('order'+str(j))
            
    #w_dates = wobble_res['dates'].value
    w_dates = wobble_res['dates'][()] #apparently new version of .value which horught depreciation warnings
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
        
    #initialize array for w_RVs_own which circumvent the broken combine in wobble
    w_RVs_own = w_RVs * 0.0
    w_RVs_unweighted = w_RVs * 0.0
    for ep in range (len(w_RVs_own)):
        w_RVs_own[ep] = np.nanmedian(w_order_RVs[:,ep])
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
        
    w_RVs = w_RVs_own #overwrite the wobble combined RVs with self Combined
    w_bervs = wobble_res['bervs'][()]
    #make own (proper) barycentic corrections
    w_RVs_barycorr = np.zeros(len(w_dates))
    try:
        f = open(i[0] + "_" + i[2] + "_bervs.pickle", "rb")
        w_RVs_barycorr = pickle.load(f)
        f.close()
    except FileNotFoundError:
        for n in tqdm(range(len(w_RVs_barycorr))):
            w_RVs_barycorr[n]=bary.get_BC_vel(w_dates[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs[n]/lightvel)[0]
        f = open(i[0]+"_"+i[2]+"_bervs.pickle", "wb")
        pickle.dump(w_RVs_barycorr, f)
        f.close()
        
    
    #read in SERVAL
    ser_rvc = np.loadtxt(servaldir+i[1]+"/"+i[1]+".rvc.dat")
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
    ser_addinfo = ser_addinfo[indices_serval]
    w_dates = w_dates[indices_wobble]
    w_RVs_barycorr = w_RVs_barycorr[indices_wobble]
    w_RVs = w_RVs[indices_wobble]
    w_RVs_own = w_RVs_own[indices_wobble]
    w_RVs_unweighted = w_RVs_unweighted[indices_wobble]
    w_bervs = w_bervs[indices_wobble]
    w_order_RVs = w_order_RVs[:,indices_wobble]
    w_order_RV_scatter = w_order_RV_scatter[indices_wobble]
    w_RVs_original = w_RVs_original[indices_wobble]
    
    
    #simple plot of RVs (only matched)
    ax1.plot(w_dates, w_RVs_barycorr - np.nanmean(w_RVs_barycorr),"x", label="Wobble")
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
    #for ep in range(len(w_dates)):
    #    ax1.plot(np.zeros(N_ord) + w_dates[ep], w_order_RVs[:,ep]-w_RVs[ep], "x")
    for ep in range(len(w_dates)):
        ax1.plot(np.zeros(N_ord) + w_dates[ep], w_order_RVs[:,ep] + w_bervs[ep], "x")
    ax1.plot(w_dates, w_RVs + w_bervs, "ko")
    ax1.plot(w_dates, w_RVs_unweighted + w_bervs, "wo")
    ax1.set_ylabel('RVs per order and combined [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV scatter of Wobble orders "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    #plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    #RV differences serval-wobble
    # seems to be correlted with barycentric correction 
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,5] -np.nanmean(ser_rvc[:,5])) - (w_RVs_barycorr-np.nanmean(w_RVs_barycorr)), "x")
    ax1.set_ylabel('delta_RVs [m/s]')
    ax1.set_xlabel('jd')
    plt.title("RV differences (SERVAL-Wobble) for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    #plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    
    #using simplified additive barycorr -> differences ofup to ca 5m/s ~ 1%
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,5] -np.nanmean(ser_rvc[:,5])) - (w_RVs+w_bervs-np.nanmean(w_RVs+w_bervs)), "x", label="w")
    ax1.plot(ser_rvc[:,0], (ser_rvc[:,5] -np.nanmean(ser_rvc[:,5])) - (w_RVs_original+w_bervs-np.nanmean(w_RVs_original+w_bervs)), "x" , label="w_originals") #original wobble combine yiels much worse scatter
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

    #flat until 10km/s then anticorrelated (fairly significantly)
    ax1.plot(ser_rvc[:,7], (ser_rvc[:,5] -np.nanmean(ser_rvc[:,5])) - (w_RVs_barycorr-np.nanmean(w_RVs_barycorr)), "x")
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
    
    ax1.plot(w_bervs, (ser_rvc[:,5] -np.nanmean(ser_rvc[:,5])) - (w_RVs_barycorr-np.nanmean(w_RVs_barycorr)), "x")
    ax1.set_ylabel('delta_RVs [m/s]')
    ax1.set_xlabel('w_BERV [m/s]')
    plt.title("Correlation of RV diff (SERVAL-Wobble) with BERV for "+i[0]+" ("+i[2]+") "+" - "+i[1])
    plt.grid(True)
    plt.tight_layout()
    #plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    
        
    
    
    

plt.close(fig)
pp.close()
    



