import numpy as np
import h5py
from tqdm import tqdm
import barycorrpy as bary 
from tabulate import tabulate
from scipy.constants import codata 
lightvel = codata.value('speed of light in vacuum') #for barycorr
# CAHA Coordinates for barycorr
_lat = 37.2236
_lon = -2.54625
_elevation = 2168.

vels_dir = "/data/cmatthe/vels_dir/"

#objects = [["GJ1148", "J11417+427", "vis", 41.382]]
#objects = [["GJ876", "J22532-142", "vis", 61.082]]
#objects = [["GJ436", "J11421+267", "vis", 2.644]]
objects = [["GJ3512", "J08413+594", "vis", 203.59]]

#bary_starname = "HD 199305"
#objects = [["GJ809A", "J20533+621", "vis"]]

corrections = True
correct_w_for_drift = True

for i in objects:
    #kt = "Kt3_snr+ds_l5_reg"
    #kt = "Kt3__l4reg_drift_shift"
    kt = "Kt3__loop4_reg_snr"
    #kt = "Kt3_rs_orderwise_test_0_reg_snr"
    #kt = "Kt3_loop_2_reg"
    #kt='Kt3_stitched_def'
    wobble_res_lx39 = h5py.File("/data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_Kstar0_"+kt+".hdf5",'r')
    wobble_res=wobble_res_lx39


    #w_dates = wobble_res['dates'].value
    w_dates = wobble_res['dates'][()]
    w_dates_utc = wobble_res['dates_utc'][()]
    #use selfcombined RVs a wobble RVs seem broken currently, actually maybe they are fine now (05.04.2019)
    w_RVs = wobble_res['star_time_rvs'][()]
    w_RVs_original = w_RVs
    w_RVs_er = wobble_res['star_time_sigmas'][()]
    w_bervs = wobble_res['bervs'][()]
    
     #barycorr vor wobble_orig
    w_RVs_original_barycorr = np.zeros(len(w_dates))
    if True:
        for n in tqdm(range(len(w_RVs_original_barycorr))):
            if "bary_starname" in locals():
                w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel)[0]
            else:    
                w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel)[0]
    
    
    if corrections == True:
        servaldir_lx39="/data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
        servaldir=servaldir_lx39
        #read in SERVAL
        #ser_rvc = np.loadtxt(servaldir+i[1]+"/"+i[1]+".rvc.dat") # rvc
        #if use_avcn == True:
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
                
        #8 sa drift, 3 drift, 9 NZP
        ser_corr = - ser_rvc[:,8] - ser_rvc[:,3] - ser_rvc[:,9] # use with avcn
        ser_corr_wob = ser_corr
        #optionally remove drift correction if this has been performed during data file generation
        if correct_w_for_drift == False:
            ser_corr_wob = ser_corr + ser_rvc[:,3]
        
        #for  this we need to match the observation by the JDs -> start with wobble date and find the one with the lowest timediff from serval
        indices_serval = [] 
        indices_wobble = []
        for n in range(len(w_dates)):
            ind_jd = np.where(np.abs(ser_rvc[:,0]-w_dates[n]) == np.nanmin(np.abs(ser_rvc[:,0]-w_dates[n])))[0][0]
            if (ser_rvc[ind_jd,0]-w_dates[n])*24*60<20.: #only takes matches closer than 20 minutes
                indices_serval.append(ind_jd)
                indices_wobble.append(n)
        #now set up all the data according to the indices
        ser_rvc = ser_rvc[indices_serval]
        ser_corr = ser_corr[indices_serval] # corrcetion: ser_rvc[:,5] + ser_corr = ser_rvc[:,1] within 10^-10ms^-1
        ser_corr_wob = ser_corr_wob[indices_serval]
        
        w_dates = w_dates[indices_wobble]
        w_RVs_original_barycorr = w_RVs_original_barycorr[indices_wobble]       + ser_corr_wob
        w_RVs_er = w_RVs_er[indices_wobble]
        
    #TODO rigorous error treatement (e.g. for additinal error due to correction)
    array = np.array([w_dates, w_RVs_original_barycorr, w_RVs_er])#NOTE errors assume additional error due to barycorr is small
    array = np.ndarray.transpose(array)
    
    np.savetxt(vels_dir + i[0] + kt + '.vels', array, fmt ='%.18f')
    
    if corrections == True:
        array = np.array([ser_rvc[:,0], ser_rvc[:,1],ser_rvc[:,2]])#NOTE errors assume additional error due to barycorr is small
        array = np.ndarray.transpose(array)
        np.savetxt(vels_dir + i[0] + "_serval_avcn" + '.vels', array, fmt ='%.18f')
    
    #f = open(vels_dir + i[0] + kt + '.vels', 'w')
    #f.write(tabulate(array))
    #f.close()
