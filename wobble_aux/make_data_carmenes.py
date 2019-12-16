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

#https://github.com/astropy/astropy/issues/8981 alternate mirrors for iers
from astropy.utils import iers
from astropy.utils.iers import conf as iers_conf
#iers_conf.iers_auto_url
'https://maia.usno.navy.mil/ser7/finals2000A.all'
iers_conf.iers_auto_url = 'https://astroconda.org/aux/astropy_mirror/iers_a_1/finals2000A.all'
iers_conf.iers_auto_url_mirror = 'https://astroconda.org/aux/astropy_mirror/iers_a_2/finals2000A.all'

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

def read_data_from_fits(filelist, arm='vis', starname=None):
    names = pd.read_csv('carmenes_aux_files/name_conversion_list.csv')
    name_dict = dict(zip(names['#Karmn'], names['Name']))
    # input : a list of filenames
    N = len(filelist)  # number of epochs
    M, R = dimensions(arm)
    data = [np.zeros((N, M)) for r in range(R)]
    ivars = [np.zeros((N, M)) for r in range(R)]
    xs = [np.zeros((N, M)) for r in range(R)]
    empty = np.array([], dtype=int)
    pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N) ,np.zeros(N)
    for n, f in enumerate(filelist):
        sp = fits.open(f)
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
        date = bary.JDUTC_to_BJDTDB(jd_mid, starname)[0]
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
                                   longi=_lon, alt=_elevation)[0]  # m/s
        airms[n] = sp[0].header['AIRMASS']
        try:
            wave = sp['WAVE'].data
            spec = sp['SPEC'].data
            sig = sp['SIG'].data
        except Exception as e:
            print('{} Skipping file {}.'.format(e, f))
            empty = np.append(empty, n)
            continue
        # save stuff
        for r in range(R):
            data[r][n, :] = spec[r, :]
            ivars[r][n, :] = 1 / sig[r, :]**2
            #xs[r][n, :] = wave[r, :] # replaced with drfit corrected version
            for l in range(len(data[r][n,:])):
                lambda_drifts = lambda_drift(wave[r, l], drifts[n])
            xs[r][n, :] = wave[r, :] - lambda_drifts

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
                
def make_data(starname, arm, data_directory):
    filelist = glob.glob(data_directory + 'CARM_raw_data/{0}/*sci-gtoc-{1}_A.fits'.format(starname, arm))
    data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm= arm, starname= None)
    hdffile = data_directory+'{0}_{1}_drift_shift_e2ds.hdf5'.format(starname, arm)
    write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    if arm == "nir":
        split_orders_file(hdffile) # save  an aditional split copy



if __name__ == "__main__":
    data_directory="../data/"
    
    
    
    if True: # GJ876 :vis
        starname = "GJ436"
        arm = "vis"
        make_data(starname, arm, data_directory)
        
    if True: # GJ876 :nir
        starname = "GJ436"
        arm = "nir"
        make_data(starname, arm, data_directory)
