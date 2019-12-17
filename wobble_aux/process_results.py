#reproduce functionality of make_vels and compare_ws
import numpy as np
import h5py
from tqdm import tqdm
import barycorrpy as bary 
from tabulate import tabulate
from scipy.constants import codata 
lightvel = codata.value('speed of light in vacuum') #for barycorr
import os
# CAHA Coordinates for barycorr
_lat = 37.2236
_lon = -2.54625
_elevation = 2168.


#1) make basic make_vels that just outputs a wobble vels file
def basic_make_vels(vels_dir, results_file):
    wobble_res = h5py.File(results_file, 'r')
    w_dates = wobble_res['dates'][()]
    w_dates_utc = wobble_res['dates_utc'][()]
    w_RVs = wobble_res['star_time_rvs'][()]
    w_RVs_er = wobble_res['star_time_sigmas'][()]
    
    w_RVs_barycorr = np.zeros(len(w_dates))
    for n in tqdm(range(len(w_RVs_barycorr))):
        try:
            w_RVs_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs[n]/lightvel,
                                                           leap_update = False #HACK barycorrpy issue 27
                                                           )[0]
        except:
            print("Barycentric correction during make_vels failed. Not outputting a .vels file")
            return
        
    array = np.array([w_dates, w_RVs_original_barycorr, w_RVs_er])#NOTE errors assume additional error due to barycorr is small
    array = np.ndarray.transpose(array)
    
    file_basename, = os.path.splitext(os.path.basename(results_file))
    np.savetxt(vels_dir + file_basename + '.vels', array, fmt ='%.18f')
    
#also include corrections
def make_vels(vels_dir, results_file, serval_dir, carmencita_ID, starname, correct_w_for_drift = False):
    #TODO find Carmencita ID from starname?
    #Import wobble results
    wobble_res = h5py.File(results_file, 'r')
    w_dates = wobble_res['dates'][()]
    w_dates_utc = wobble_res['dates_utc'][()]
    w_RVs = wobble_res['star_time_rvs'][()]
    w_RVs_er = wobble_res['star_time_sigmas'][()]
    
    w_RVs_barycorr = np.zeros(len(w_dates))
    for n in tqdm(range(len(w_RVs_barycorr))):
        try:
            w_RVs_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs[n]/lightvel,
                                                           leap_update = False #HACK barycorrpy issue 27
                                                           )[0]
        except:
            print("Barycentric correction during make_vels failed. Not outputting a .vels file")
            return
    #import serval results
    ser_avcn = np.loadtxt(servaldir+ carmencita_ID +"/"+ carmencita_ID +".avcn.dat")
    # remove entries with nan in drift
    ind_finitedrift = np.isfinite(ser_avcn[:,3])
    ser_avcn = ser_avcn[ind_finitedrift]
    
    #8 sa drift, 3 drift, 9 NZP
    ser_corr = - ser_avcn[:,8] - ser_avcn[:,3] - ser_avcn[:,9] # use with avcn
    ser_corr_wob = ser_corr
    #optionally remove drift correction if this has been performed during data file generation
    if correct_w_for_drift:
        ser_corr_wob = ser_corr + ser_rvc[:,3]
        
    #match the observation by the JDs -> start with wobble date and find the one with the lowest timediff from serval
    indices_serval = [] 
    indices_wobble = []
    for n in range(len(w_dates)):
        ind_jd = np.where(np.abs(ser_avcn[:,0]-w_dates[n]) == np.nanmin(np.abs(ser_avcn[:,0]-w_dates[n])))[0][0]
        if (ser_avcn[ind_jd,0]-w_dates[n])*24*60<20.: #only takes matches closer than 20 minutes
            indices_serval.append(ind_jd)
            indices_wobble.append(n)
            
    #now set up all the data according to the indices
    ser_avcn = ser_avcn[indices_serval]
    ser_corr = ser_corr[indices_serval]
    ser_corr_wob = ser_corr_wob[indices_serval]
    
    w_dates = w_dates[indices_wobble]
    w_RVs_barycorr = w_RVs_barycorr[indices_wobble] + ser_corr_wob
    w_RVs_er = w_RVs_er[indices_wobble]
    
    #TODO rigorous error treatement (e.g. for additinal error due to correction)
    array = np.array([w_dates, w_RVs_original_barycorr, w_RVs_er])#NOTE errors assume additional error due to barycorr is small
    array = np.ndarray.transpose(array)
    
    file_basename, = os.path.splitext(os.path.basename(results_file))
    np.savetxt(vels_dir + file_basename + '.vels', array, fmt ='%.18f')
    
    array = np.array([ser_rvc[:,0], ser_rvc[:,1],ser_rvc[:,2]])#NOTE errors assume additional error due to barycorr is small
    array = np.ndarray.transpose(array)
    np.savetxt(vels_dir + starname + "_serval_avcn" + '.vels', array, fmt ='%.18f')
    
