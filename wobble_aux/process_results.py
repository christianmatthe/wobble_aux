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

''' Doesn't work needs bary_starname
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
'''

#also include corrections
def make_vels(vels_dir, results_file, serval_dir, carmenes_object_ID, starname, correct_w_for_drift = False):
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
    ser_avcn = np.loadtxt(servaldir+ carmenes_object_ID +"/"+ carmenes_object_ID +".avcn.dat")
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
    
#plan: make standart results import function (wobble-Serval cross match)
class Results_ws():
    """
    Object contains results from both Wobble and SERVAL matched by date for correction of wobble with SERVAL corrections (NZP and SA (drifts already applied in data generation))) and for comparison.
    
    # TODO document all input parameters
    Parameters
    ----------
    wobble_file : `path`
        a wobble results filename
    serval_file : `path` (.avcn.dat)
        a serval avcn results filename
    """
    def __init__(self
                 , wobble_file
                 , serval_dir
                 , bary_starname
                 , load_bary = True
                 , archive = True
                 , bary_archive = os.path.dirname(__file__) + "/" + "../results/bary_archive/"
                 )
        self.bary_archive = bary_archive
        os.makedirs(bary_archive, exist_ok = True)
        self.bary_starname = bary_starname
        # Import wobble
        self.w_file_basename, = os.path.splitext(os.path.basename(wobble_file))
        self.wobble_res = h5py.(wobble_file, 'r')
        self.w_dates = wobble_res['dates'][()] 
        self.w_dates_utc = wobble_res['dates_utc'][()]
        self.w_orders = wobble_res['orders'][()]
        self.w_epochs = wobble_res['epochs'][()]
        self.w_RVs = wobble_res['star_time_rvs'][()]
        self.w_RVs_er = wobble_res['star_time_sigmas'][()]
        #orderwise RVs 
        N_ord = len(self.w_orders) #should be equivalent to ["R"]
        k=np.arange(N_ord)
        keys_orderRVs=[]
        for j in k:
            keys_orderRVs.append('order'+str(j))
        #initialize arrays
        self.w_order_RVs = np.zeros((N_ord,len(w_dates)))
        self.w_order_RVs_ivars = np.zeros((N_ord,len(w_dates)))
        self.w_order_RV_scatter = np.zeros(len(w_dates))
        for j, k in enumerate(keys_orderRVs):
            temp = wobble_res[k] #load order dataset
            tempRVs = temp['star_rvs'][()] # load RVs for this order at all epochs
            tempRVs_ivars = temp['star_ivars_rvs'][()] # load error of RVs
            self.w_order_RVs[j,:] = tempRVs # fill values into Rows of array
            self.w_order_RVs_ivars[j,:] = tempRVs_ivars
            if np.any(tempRVs_ivars < 0.0)== True:
                print("Found negative ivars.")
        #we may want to drop  these orders as they may represent *maxima* in optimization
        for j in range (len(w_dates)):
            self.w_order_RV_scatter[j] = np.nanstd(w_order_RVs[:,j])
        
        #barycentric correction
        #load barycorrected rvs from pickle archive to save time on recalculting them
        filename_RVs_barycorr = self.bary_archive + self.w_file_basename +"_barycorr.pkl"
        filename_order_RVs_barycorr = self.bary_archive + self.w_file_basename +"_barycorr_orders.pkl"
        if load_bary:
            try:
                with open(filename_RVs_barycorr, "rb") as f:
                    self.w_RVs_barycorr = pickle.load(f)
                with open(filename_order_RVs_barycorr, "rb") as f:
                    self.w_order_RVs_barycorr = pickle.load(f)
            except FileNotFoundError:
                print("RV file not found, recalculating barycentric corrections")
                load_bary = False
        if not load_rvs:
            self.w_RVs_barycorr = np.zeros(len(self.w_dates))
            for n in tqdm(range(len(w_RVs_barycorr))):
                self.w_RVs_barycorr[n]=bary.get_BC_vel(self.w_dates_utc[n], starname= self.bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas=self.w_RVs[n]/lightvel,
                                                        leap_update = False #HACK barycorrpy issue 27
                                                        )[0]
                
            self.w_order_RVs_barycorr = self.w_order_RVs
            for order in tqdm(range(len(self.w_order_RVs))):
                for ep in tqdm(range(len(self.w_order_RVs[order]))):
                    self.w_order_RVs_barycorr[order, ep] = bary.get_BC_vel(self.w_dates_utc[ep], starname=self.bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas= self.w_order_RVs[order, ep]/lightvel,
                                                        leap_update = False #HACK barycorrpy issue 27
                                                        )[0]
            if archive:
                with open(filename_RVs_barycorr) as f:
                    pickle.dump(self.w_RVs_barycorr, f)
                with open(filename_order_RVs_barycorr) as f:
                    pickle.dump(self.w_order_RVs_barycorr, f)
            
        # Import SERVAL
        ser_avcn = np.loadtxt(serval_file)
        # remove entries with nan in drift
        ind_finitedrift = np.isfinite(ser_avcn[:,3])
        ser_avcn = ser_avcn[ind_finitedrift]
        
        #match the observation by the JDs -> start with wobble date and find the one with the lowest timediff from serval
        indices_serval = [] 
        indices_wobble = []
        for n in range(len(w_dates)):
            ind_jd = np.where(np.abs(ser_avcn[:,0]-w_dates[n]) == np.nanmin(np.abs(ser_avcn[:,0]-w_dates[n])))[0][0]
            if (ser_avcn[ind_jd,0]-w_dates[n])*24*60<20.: #only takes matches closer than 20 minutes
                indices_serval.append(ind_jd)
                indices_wobble.append(n)
        #TODO finish this
        
        
        
        
     #TODO   
    def apply_corrections():
        #8 sa drift, 3 drift, 9 NZP
        ser_corr = - ser_avcn[:,8] - ser_avcn[:,3] - ser_avcn[:,9] # use with avcn
        ser_corr_wob = ser_corr
        #optionally remove drift correction if this has been performed during data file generation
        if correct_w_for_drift:
            ser_corr_wob = ser_corr + ser_rvc[:,3]
    
