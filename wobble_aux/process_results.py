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
    w_dates = wobble_res['dates'][()]
    w_dates_utc = wobble_res['dates_utc'][()]
    w_RVs = wobble_res['star_time_rvs'][()]
    w_RVs_original = w_RVs
    w_RVs_er = wobble_res['star_time_sigmas'][()]
    
    w_RVs_original_barycorr = np.zeros(len(w_dates))
    for n in tqdm(range(len(w_RVs_original_barycorr))):
        try:
            w_RVs_original_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs_original[n]/lightvel,
                                                           leap_update = False #HACK barycorrpy issue 27
                                                           )[0]
        except:
            print("Barycentric correction during make_vels failed. Not outputting a .vels file")
            return
        
    array = np.array([w_dates, w_RVs_original_barycorr, w_RVs_er])#NOTE errors assume additional error due to barycorr is small
    array = np.ndarray.transpose(array)
    
    file_basename, = os.path.splitext(os.path.basename(results_file))
    np.savetxt(vels_dir + file_basename + '.vels', array, fmt ='%.18f')
            
    
