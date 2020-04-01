import numpy as np
from scipy.io.idl import readsav
from scipy.interpolate import interp1d
# from harps_hacks import read_harps
import h5py
import math
from astropy.io import fits
from astropy.time import Time
#from astropy.utils.iers import IERS_A_URL_MIRROR #IF main  source unavailable
import shutil
import glob
import os
import barycorrpy as bary
import pandas as pd
from time import time
from tqdm import tqdm

##https://github.com/astropy/astropy/issues/8981 alternate mirrors for iers (old)
#from astropy.utils import iers
#from astropy.utils.iers import conf as iers_conf
##iers_conf.iers_auto_url
#'https://maia.usno.navy.mil/ser7/finals2000A.all'
#iers_conf.iers_auto_url = 'https://astroconda.org/aux/astropy_mirror/iers_a_1/finals2000A.all'
#iers_conf.iers_auto_url_mirror = 'https://astroconda.org/aux/astropy_mirror/iers_a_2/finals2000A.all'

#works as of 05.March.2020
from astropy.utils import iers
iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')

from scipy.constants import codata 
lightvel = codata.value('speed of light in vacuum') #for barycorr

# CAHA Coordinates
_lat = 37.2236
_lon = -2.54625
_elevation = 2168.


def lambda_drift(v_drift, lambda_0):
    return v_drift * lambda_0 / lightvel
    

def dimensions(arm):
    if arm == "vis":
        M = 4096
        R = 61
    elif arm == "nir":
        # TODO:
        M = 4080
        R = 28
    else:
        print("{} not recognized. valid options are: \"vis\" or"
              " \"nir\"".format(arm))
        return
    return M, R

def read_data_from_fits(filelist, arm='vis', starname = None, serval_dir = None):
    names = pd.read_csv(os.path.dirname(os.path.abspath(__file__)) + '/carmenes_aux_files/name_conversion_list.csv')
    name_dict = dict(zip(names['#Karmn'], names['Name']))
    # input : a list of filenames
    N = len(filelist)  # number of epochs
    M, R = dimensions(arm)
    data = [np.zeros((N, M)) for r in range(R)]
    ivars = [np.zeros((N, M)) for r in range(R)]
    xs = [np.zeros((N, M)) for r in range(R)]
    empty = np.array([], dtype=int)
    pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, total_drifts = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N) ,np.zeros(N), np.zeros(N)
    
    #file headers sometimes list object not by Carmenes ID but instead with catalogue names (e.g. GJ436 obs 81 lists Ross905 istead of J11421+267. FIX:
    sp = fits.open(filelist[0])
    carmenes_object_ID_master = str(sp[0].header['OBJECT']).strip() #ID in header has extra space in front of it
    print("Object_ID: ", carmenes_object_ID_master) 
    
    for n, f in enumerate(tqdm(filelist)):
        sp = fits.open(f)
        
        if not serval_dir:
            print("no serval directory supplied. Not correcting for NZP")
            #include NZP by adding them to drifts before correction
        else:
            nzp_shift = True
            carmenes_object_ID = str(sp[0].header['OBJECT']).strip() #ID in header has extra space in front of it
            if carmenes_object_ID != carmenes_object_ID_master:
                print()
                print("mismatched object ID: " ,carmenes_object_ID)
                print("n, f:", n, f)
                print()
                #use master instead:
                carmenes_object_ID = carmenes_object_ID_master
            ser_avcn = np.loadtxt(serval_dir+ carmenes_object_ID +"/"+ carmenes_object_ID +".avcn.dat")
            nzp = ser_avcn[:,9]
        
        try:
            pipeline_rvs[n] = sp[0].header['HIERARCH CARACAL SERVAL RV'] * 1.e3 # m/s
            pipeline_sigmas[n] = sp[0].header['HIERARCH CARACAL SERVAL E_RV'] * 1.e3 # m/s
        except KeyError:
            pipeline_rvs[n] = 0
            pipeline_sigmas[n] = 0
        try:
            drifts[n] = sp[0].header['HIERARCH CARACAL DRIFT FP RV']
        except KeyError:
            print("WARNING: {0} Drift missing. Skipping this one.".format(f))
            empty = np.append(empty, n)
            continue
        if not starname:
            starname = name_dict[sp[0].header['OBJECT']]
        jd_start = Time(sp[0].header['DATE-OBS'])
        jd_mid = jd_start.jd + sp[0].header['HIERARCH CARACAL TMEAN'] * 1/(24*60*60)
        dates_utc[n] = jd_mid
        # for nir ignore all dates before 2016. recommended by Adrian
        date = bary.JDUTC_to_BJDTDB(jd_mid, starname,
                                                           leap_update = False #HACK barycorrpy issue 27
                                                           )[0]
        if date >=2457754.5:#1 JAN 2017
            dates[n] = date
        else:
            if arm == "vis":
                dates[n] = date
            elif arm == "nir":
                print("Date is before 2017 for NIR measurement. Skipping this one.")
                empty = np.append(empty, n)
                continue
            else:
                print("{} not recognized. valid options are: \"vis\" or"
                " \"nir\"".format(arm))
                return
        bervs[n] = bary.get_BC_vel(jd_mid, starname=starname, lat=_lat,
                                   longi=_lon, alt=_elevation,
                                                           leap_update = False #HACK barycorrpy issue 27
                                                           )[0]  # m/s
        airms[n] = sp[0].header['AIRMASS']
        try:
            wave = sp['WAVE'].data
            spec = sp['SPEC'].data
            sig = sp['SIG'].data
        except Exception as e:
            print('{} Skipping file {}.'.format(e, f))
            empty = np.append(empty, n)
            continue
        
    
        total_drifts[n] = drifts[n]
        if nzp_shift:
            #match only the nth date
            #match the observation by the JDs -> start with wobble date and find the one with the lowest timediff from serval
            indices_serval = [] 
            indices_wobble = []
            #for n in range(len(dates)):
            ind_jd = np.where(np.abs(ser_avcn[:,0]- dates[n]) == np.nanmin(np.abs(ser_avcn[:,0]- dates[n])))[0][0]
            if (ser_avcn[ind_jd,0]-dates[n])*24*60<20.: #only takes matches closer than 20 minutes
                indices_serval.append(ind_jd)
                indices_wobble.append(n)
                    
                # add NZP to drift corrections that match dates in SERVAL:
                total_drifts[n] = total_drifts[n] + nzp[indices_serval]
                #print("totals:", total_drifts)
        
        # save stuff
        for r in range(R):
            data[r][n, :] = spec[r, :]
            ivars[r][n, :] = 1 / sig[r, :]**2
            #xs[r][n, :] = wave[r, :] # replaced with drfit corrected version
    
            for l in range(len(data[r][n,:])):
                lambda_drifts = lambda_drift(total_drifts[n], wave[r, l])
                xs[r][n, l] = wave[r, l] - lambda_drifts

    # delete data with missing attributes:
    for r in range(R):
        data[r] = np.delete(data[r], empty, axis=0)
        ivars[r] = np.delete(ivars[r], empty, axis=0)
        xs[r] = np.delete(xs[r], empty, axis=0)
    pipeline_rvs = np.delete(pipeline_rvs, empty)
    pipeline_sigmas = np.delete(pipeline_sigmas, empty)
    dates = np.delete(dates, empty)
    bervs = np.delete(bervs, empty)
    airms = np.delete(airms, empty)
    drifts = np.delete(drifts, empty)
    dates_utc= np.delete(dates_utc, empty)

    # re-introduce BERVs to HARPS results:
    # pipeline_rvs -= bervs
    # pipeline_rvs -= np.mean(pipeline_rvs)

    return data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc

# Not used anywhere I think 03.12.2019
#def savfile_to_filelist(savfile, destination_dir='../data/'):
    ## copies CCF + E2DS files to destination_dir and returns a list of the CCFs
    #s = readsav(savfile)
    #filelist = []
    #files = [f.decode('utf8') for f in s.files]
    #for f in files:
        #shutil.copy2(f, destination_dir)
        #spec_file = str.replace(f, 'ccf_G2', 'e2ds')
        #shutil.copy2(spec_file, destination_dir)
        #basename = f[str.rfind(f, '/') + 1:]
        #filelist = np.append(filelist, destination_dir+basename)
    #return filelist


def write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs,
               airms, drifts, dates_utc, filenames, hdffile):
    h = h5py.File(hdffile, 'w')
    dset = h.create_dataset('data', data=data)
    dset = h.create_dataset('ivars', data=ivars)
    dset = h.create_dataset('xs', data=xs)
    dset = h.create_dataset('pipeline_rvs', data=pipeline_rvs)
    dset = h.create_dataset('pipeline_sigmas', data=pipeline_sigmas)
    dset = h.create_dataset('dates', data=dates)
    dset = h.create_dataset('bervs', data=bervs)
    dset = h.create_dataset('airms', data=airms)
    dset = h.create_dataset('drifts', data=drifts)
    dset = h.create_dataset('dates_utc', data=dates_utc)
    filenames = [a.encode('utf8') for a in filenames] # h5py workaround
    dset = h.create_dataset('filelist', data=filenames)
    h.close()
    
def split_orders(array):
    array_old = array
    shape_old = array.shape
    x_width_new = shape_old[2] // 2
    shape_new = (2 * shape_old[0], shape_old[1], x_width_new)
    array_new = np.zeros(shape_new)
    
    for i in range(2 * shape_old[0]):
        if i % 2 == 0:
            array_new[i] = array_old[i//2,:,:x_width_new]
        if i % 2 == 1:
            array_new[i] = array_old[i//2,:,x_width_new:]
    return array_new

def split_orders_file(filename):
    split_sets = ["data", "ivars", "xs"]
    with h5py.File(filename,'r') as f:
        with h5py.File(filename.split("e2ds")[0] + "split_e2ds.hdf5",'w') as g:
            for key in list(f.keys()):
                temp = f[key][()]
                if key in split_sets:
                    temp_split = split_orders(temp)
                    temp = temp_split
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)
                
def make_data(starname, arm, data_directory, simbad_name = None, serval_dir = None, nzp_shift = True):
    #if not simbad_name:
        #simbad_name = starname
    #print(starname)
    filelist = glob.glob(data_directory + 'CARM_raw_data/{0}/*sci-gtoc-{1}_A.fits'.format(starname, arm))
    if nzp_shift == True:
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm= arm, starname = simbad_name, serval_dir = serval_dir)
        hdffile = data_directory+'{0}_{1}_drift+nzp_e2ds.hdf5'.format(starname, arm)
    else:
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm= arm, starname= simbad_name)
        hdffile = data_directory+'{0}_{1}_drift_shift_e2ds.hdf5'.format(starname, arm)
    write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    if arm == "nir":
        split_orders_file(hdffile) # save  an aditional split copy



if __name__ == "__main__":
    data_directory="../data/"
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" +"../data/servaldir/CARM_VIS/" #read data already includes name dictionary
    
    
    #if True: # GJ1148 :vis
        #starname = "GJ1148"
        ##simbad_name = "Ross 1003"
        #arm = "vis"
        #make_data(starname, arm, data_directory, serval_dir = serval_dir)
        
######################### Laptop examples
if True: # GJ436 :vis
        starname = "GJ436"
        #simbad_name = "GJ436"
        arm = "vis"
        make_data(starname, arm, data_directory, serval_dir = serval_dir, nzp_shift = True)
        
if True: # GJ3473 :vis
        starname = "GJ3473"
        simbad_name = "G 50-16"
        arm = "vis"
        make_data(starname, arm, data_directory, simbad_name = simbad_name, serval_dir = serval_dir, nzp_shift = True)
        
######################### deprecated below    
    
    #if True: # Barnard vis
        #starname = "Barnard"
        #simbad_name = "GJ699"
        #arm = "vis"
        #make_data(starname, arm, data_directory, simbad_name)
    
    #if True: # Teegarden : nir
        #starname = "Teegarden"
        #simbad_name = "GAT 1370"
        #arm = "nir"
        #make_data(starname, arm, data_directory, simbad_name)
    
    #if True: # GJ1148 :vis
        #starname = "GJ1148"
        #arm = "vis"
        #make_data(starname, arm, data_directory)
    
    #if True: # GJ3473? / G050-16A / G 50-16 inn SIMBAD :vis
        #starname = "GJ3473"
        #simbad_name = "G 50-16"
        #arm = "vis"
        #make_data(starname, arm, data_directory, simbad_name)
    
    #if True: # GJ436 :vis
        #starname = "GJ436"
        #arm = "vis"
        #make_data(starname, arm, data_directory)
        
    #if True: # GJ436 :nir
        #starname = "GJ436"
        #arm = "nir"
        #make_data(starname, arm, data_directory)
