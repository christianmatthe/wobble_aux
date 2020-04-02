#reproduce functionality of make_vels and compare_ws
import numpy as np
import h5py
from tqdm import tqdm
import barycorrpy as bary 
from tabulate import tabulate
from scipy.constants import codata 
lightvel = codata.value('speed of light in vacuum') #for barycorr
import os
import pickle

#https://github.com/astropy/astropy/issues/9427 workaround for broken iers mirror
from astropy.utils import iers
iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')

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
#def make_vels(vels_dir, results_file, serval_dir, carmenes_object_ID, starname, correct_w_for_drift = False):
    ##TODO find Carmencita ID from starname?
    ##Import wobble results
    #wobble_res = h5py.File(results_file, 'r')
    #w_dates = wobble_res['dates'][()]
    #w_dates_utc = wobble_res['dates_utc'][()]
    #w_RVs = wobble_res['star_time_rvs'][()]
    #w_RVs_er = wobble_res['star_time_sigmas'][()]
    
    #w_RVs_barycorr = np.zeros(len(w_dates))
    #for n in tqdm(range(len(w_RVs_barycorr))):
        #try:
            #w_RVs_barycorr[n]=bary.get_BC_vel(w_dates_utc[n], starname=i[0], lat=_lat, longi=_lon, alt=_elevation, zmeas=w_RVs[n]/lightvel,
                                                           #leap_update = False #HACK barycorrpy issue 27
                                                           #)[0]
        #except:
            #print("Barycentric correction during make_vels failed. Not outputting a .vels file")
            #return
    ##import serval results
    #ser_avcn = np.loadtxt(servaldir+ carmenes_object_ID +"/"+ carmenes_object_ID +".avcn.dat")
    ## remove entries with nan in drift
    #ind_finitedrift = np.isfinite(ser_avcn[:,3])
    #ser_avcn = ser_avcn[ind_finitedrift]
    
    ##8 sa drift, 3 drift, 9 NZP
    #ser_corr = - ser_avcn[:,8] - ser_avcn[:,3] - ser_avcn[:,9] # use with avcn
    #ser_corr_wob = ser_corr
    ##optionally remove drift correction if this has been performed during data file generation
    #if correct_w_for_drift:
        #ser_corr_wob = ser_corr + ser_rvc[:,3]
        
    ##match the observation by the JDs -> start with wobble date and find the one with the lowest timediff from serval
    #indices_serval = [] 
    #indices_wobble = []
    #for n in range(len(w_dates)):
        #ind_jd = np.where(np.abs(ser_avcn[:,0]-w_dates[n]) == np.nanmin(np.abs(ser_avcn[:,0]-w_dates[n])))[0][0]
        #if (ser_avcn[ind_jd,0]-w_dates[n])*24*60<20.: #only takes matches closer than 20 minutes
            #indices_serval.append(ind_jd)
            #indices_wobble.append(n)
            
    ##now set up all the data according to the indices
    #ser_avcn = ser_avcn[indices_serval]
    #ser_corr = ser_corr[indices_serval]
    #ser_corr_wob = ser_corr_wob[indices_serval]
    
    #w_dates = w_dates[indices_wobble]
    #w_RVs_barycorr = w_RVs_barycorr[indices_wobble] + ser_corr_wob
    #w_RVs_er = w_RVs_er[indices_wobble]
    
    ##TODO rigorous error treatement (e.g. for additinal error due to correction)
    #array = np.array([w_dates, w_RVs_original_barycorr, w_RVs_er])#NOTE errors assume additional error due to barycorr is small
    #array = np.ndarray.transpose(array)
    
    #file_basename, = os.path.splitext(os.path.basename(results_file))
    #np.savetxt(vels_dir + file_basename + '.vels', array, fmt ='%.18f')
    
    #array = np.array([ser_rvc[:,0], ser_rvc[:,1],ser_rvc[:,2]])#NOTE errors assume additional error due to barycorr is small
    #array = np.ndarray.transpose(array)
    #np.savetxt(vels_dir + starname + "_serval_avcn" + '.vels', array, fmt ='%.18f')


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
                 , carmenes_object_ID
                 , bary_starname
                 , load_bary = True
                 , archive = True
                 , bary_archive = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/bary_archive/"
                 , correct_w_for_drift = False
                 , correct_drift = True
                 , correct_NZP = True
                 , correct_SA = True
                 ):
        self.wobble_file = wobble_file
        self.bary_starname = bary_starname
        #if bary_archive == "/" + "../results/bary_archive/":
            #file_dir = os.path.dirname(os.path.abspath(__file__))
            #print(file_dir )
            #bary_archive = file_dir + bary_archive
        if archive == True:
            os.makedirs(bary_archive, exist_ok = True)
        # Import wobble
        wobble_res = h5py.File(wobble_file, 'r')
        w_dates = wobble_res['dates'][()] 
        w_dates_utc = wobble_res['dates_utc'][()]
        self.w_orders = w_orders = wobble_res['orders'][()]
        w_epochs = wobble_res['epochs'][()]
        w_RVs = wobble_res['star_time_rvs'][()]
        w_RVs_er = wobble_res['star_time_sigmas'][()]
        w_bervs = wobble_res['bervs'][()]
        #orderwise RVs 
        N_ord = len(w_orders) #should be equivalent to ["R"]
        k=np.arange(N_ord)
        keys_orderRVs=[]
        for j in k:
            keys_orderRVs.append('order'+str(j))
        #initialize arrays
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
        
        #barycentric correction
        #load barycorrected rvs from pickle archive to save time on recalculting them
        # TODO Isssue in next line
        w_file_basename = os.path.splitext(os.path.basename(wobble_file))[0]
        filename_RVs_barycorr = bary_archive + w_file_basename +"_barycorr.pkl"
        filename_order_RVs_barycorr = bary_archive + w_file_basename +"_barycorr_orders.pkl"
        if load_bary:
            try:
                with open(filename_RVs_barycorr, "rb") as f:
                    w_RVs_barycorr = pickle.load(f)
                with open(filename_order_RVs_barycorr, "rb") as f:
                    w_order_RVs_barycorr = pickle.load(f)
            except FileNotFoundError:
                print("RV file not found, recalculating barycentric corrections")
                load_bary = False
        if not load_bary:
            w_RVs_barycorr = np.zeros(len(w_dates))
            for n in tqdm(range(len(w_RVs_barycorr))):
                w_RVs_barycorr[n] = bary.get_BC_vel(w_dates_utc[n], starname = bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas = w_RVs[n]/lightvel,
                                                        leap_update = False #HACK barycorrpy issue 27
                                                        )[0]
                
            w_order_RVs_barycorr = w_order_RVs
            for order in tqdm(range(len(w_order_RVs))):
                for ep in tqdm(range(len(w_order_RVs[order]))):
                    w_order_RVs_barycorr[order, ep] = bary.get_BC_vel(w_dates_utc[ep], starname=bary_starname, lat=_lat, longi=_lon, alt=_elevation, zmeas= w_order_RVs[order, ep]/lightvel,
                                                        leap_update = False #HACK barycorrpy issue 27
                                                        )[0]
            if archive:
                with open(filename_RVs_barycorr, "wb") as f:
                    pickle.dump(w_RVs_barycorr, f)
                with open(filename_order_RVs_barycorr, "wb") as f:
                    pickle.dump(w_order_RVs_barycorr, f)
            
        # Import SERVAL TODO make some of these that aren't strictly necessary for corrections optional
        # TODO replace [i] with carmenes_object_ID
        ser_avcn = np.loadtxt(serval_dir+ carmenes_object_ID +"/"+ carmenes_object_ID +".avcn.dat")
        #read in also info file, from which we get SNR and airmass
        ser_info = np.genfromtxt(serval_dir+ carmenes_object_ID +"/"+ carmenes_object_ID +".info.cvs", delimiter=";")
        ser_addinfo=np.zeros((len(ser_avcn), 2)) #initializes array
        for n in range(len(ser_avcn)):
            ind_jd = np.where(np.abs(ser_info[:,1]-ser_avcn[n,0]) == np.nanmin(np.abs(ser_info[:,1]-ser_avcn[n,0])))[0][0] #ser_info [:,1] is BJD and so is ser_rvc[n,0] this matches the ones closest to each other
            ser_addinfo[n,0] = ser_info[ind_jd, 3] # SNR
            ser_addinfo[n,1] = ser_info[ind_jd, 8] # Airmass
        #Import serval orderwise Rvs
        ser_rvo = np.loadtxt(serval_dir+ carmenes_object_ID +"/"+ carmenes_object_ID +".rvo.dat")
        
        ser_rvo_err = np.loadtxt(serval_dir+ carmenes_object_ID +"/"+ carmenes_object_ID +".rvo.daterr")
        # remove entries with nan in drift
        ind_finitedrift = np.isfinite(ser_avcn[:,3])
        ser_avcn = ser_avcn[ind_finitedrift]
        ser_rvo = ser_rvo[ind_finitedrift]
        ser_rvo_err = ser_rvo_err[ind_finitedrift]
        
        #match the observation by the JDs -> start with wobble date and find the one with the lowest timediff from serval
        indices_serval = [] 
        indices_wobble = []
        for n in range(len(w_dates)):
            ind_jd = np.where(np.abs(ser_avcn[:,0]-w_dates[n]) == np.nanmin(np.abs(ser_avcn[:,0]-w_dates[n])))[0][0]
            if (ser_avcn[ind_jd,0]-w_dates[n])*24*60<20.: #only takes matches closer than 20 minutes
                indices_serval.append(ind_jd)
                indices_wobble.append(n)
                
        #apply to all arrays and write to object
        self.ser_avcn = ser_avcn[indices_serval]
        self.ser_rvo = ser_rvo[indices_serval]
        self.ser_rvo_err = ser_rvo_err[indices_serval]
        self.ser_addinfo = ser_addinfo[indices_serval]
        #wobble
        self.w_dates = w_dates[indices_wobble]
        self.w_dates_utc = w_dates_utc[indices_wobble]
        self.w_epochs = w_epochs[indices_wobble]
        self.w_RVs_barycorr = w_RVs_barycorr[indices_wobble]                         
        self.w_RVs = w_RVs[indices_wobble]                                           
        self.w_RVs_er = w_RVs_er[indices_wobble]                 
        self.w_bervs = w_bervs[indices_wobble]
        self.w_order_RVs = w_order_RVs[:,indices_wobble]
        self.w_order_RVs_barycorr = w_order_RVs_barycorr[:,indices_wobble]
        self.w_order_RV_scatter = w_order_RV_scatter[indices_wobble]
        
        '''
        # TODO finish this
            - test with rewritten make_vels
            
        '''
        
        
    def apply_corrections(self, 
                          correct_w_for_drift = False,
                          correct_w_for_SA = False,
                          correct_w_for_NZP = False,
                          
                          correct_drift = True, correct_NZP = True, correct_SA = True,
                          
                          inverse_NZP = False):
        ser_avcn = self.ser_avcn
            
        ser_corr = np.zeros(len(self.ser_avcn))
        if correct_drift: #3 drift
            ser_corr -=  ser_avcn[:,3]
        if correct_SA:    #8 SA 
            ser_corr -=  ser_avcn[:,8]
        if correct_NZP:   #9 NZP
            ser_corr -=  ser_avcn[:,9]
        ser_corr_wob = np.zeros(ser_corr.shape)
        #if (correct_drift and not correct_w_for_drift): # in case wobble was already drift corrected in make_data
            ##ser_corr_wob += ser_avcn[:,3] #FOR SOME reason this assignment also affects ser_corr
            #ser_corr_wob = ser_corr_wob + ser_avcn[:,3]
        if correct_w_for_drift: #3 drift
            ser_corr_wob -=  ser_avcn[:,3]
        if correct_w_for_SA:    #8 SA 
            ser_corr_wob -=  ser_avcn[:,8]
        if correct_w_for_NZP:   #9 NZP
            ser_corr_wob -=  ser_avcn[:,9]
            
        if inverse_NZP:
            ser_corr_wob = ser_corr_wob + 0.5 * ser_avcn[:,9]
        #apply to object
        self.ser_corr = ser_corr
        self.ser_corr_wob = ser_corr_wob
        
        #correct ser_rvo
        for j in range(len(self.ser_rvo)):
            self.ser_rvo[j,5:] = self.ser_rvo[j,5:] + ser_corr[j]
            
        # TODO error propagation
        self.w_RVs_barycorr += ser_corr_wob
        self.w_RVs += ser_corr_wob
        self.w_order_RVs += ser_corr_wob
        
    def make_vels(self, vels_dir , output_file_basename = None):
    #outputs vels file from results_ws object
        if output_file_basename is None:
            file_basename = os.path.splitext(os.path.basename(self.wobble_file))[0]
        else:
            file_basename = output_file_basename
        #output wobble results
        array = np.array([self.w_dates, self.w_RVs_barycorr, self.w_RVs_er])
        array = np.ndarray.transpose(array)
        np.savetxt(vels_dir + file_basename + '.vels', array, fmt ='%.18f')
        #output serval results
        array = np.array([self.ser_avcn[:,0], self.ser_avcn[:,1], self.ser_avcn[:,2]])
        array = np.ndarray.transpose(array)
        np.savetxt(vels_dir + self.bary_starname + "_serval_avcn" + '.vels', array, fmt ='%.18f')
        
    def eval_vels(self, vels_dir , output_file_basename = None):
    #outputs vels file from results_ws object
        if output_file_basename is None:
            file_basename = os.path.splitext(os.path.basename(self.wobble_file))[0]
        else:
            file_basename = output_file_basename
        #output wobble results
        array = np.array([self.w_dates, self.w_RVs_barycorr, self.w_RVs_er])
        array = np.ndarray.transpose(array)
        vels_filename = vels_dir + file_basename + '.vels'
        np.savetxt(vels_filename, array, fmt ='%.18f')
        return vels_filename


if __name__ == "__main__":
    #test
    #carmenes_object_ID = "J11417+427"
    #bary_starname = "GJ1148"
    #wobble_file = "/data/cmatthe/compare_wobble_serval/wobbledir/results_GJ1148_Kstar0_Kt3_recheck_all_orders.hdf5"
    #serval_dir = "/data/cmatthe/compare_wobble_serval/servaldir/CARM_VIS/"
    #res = Results_ws(wobble_file
                 #, serval_dir
                 #, carmenes_object_ID
                 #, bary_starname
                 #, load_bary = True
                 #, archive = True)
    #print(res.w_RVs_barycorr)
    
    #vels test on laptop
    carmenes_object_ID = "J11421+267"
    bary_starname = "GJ436"
    wobble_file = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/results_GJ436_Kstar0_Kt3_git_run_wobble_test0.hdf5"
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" +"../data/servaldir/CARM_VIS/"
    vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/"
    os.makedirs(vels_dir, exist_ok = True)
    
    res = Results_ws(wobble_file
                 , serval_dir
                 , carmenes_object_ID
                 , bary_starname
                 , load_bary = True
                 , archive = True)
    res.apply_corrections()
    res.make_vels(vels_dir)
    
