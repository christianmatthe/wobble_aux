import numpy as np
from scipy.io.idl import readsav
from scipy.interpolate import interp1d
# from harps_hacks import read_harps
import h5py
import math
from astropy.io import fits
from astropy.time import Time
import shutil
import glob
import os
import barycorrpy as bary
import pandas as pd

# CAHA Coordinates
_lat = 37.2236
_lon = -2.54625
_elevation = 2168.


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
    names = pd.read_csv('name_conversion_list.csv')
    name_dict = dict(zip(names['#Karmn'], names['Name']))
    # input : a list of filenames
    N = len(filelist)  # number of epochs
    M, R = dimensions(arm)
    data = [np.zeros((N, M)) for r in range(R)]
    ivars = [np.zeros((N, M)) for r in range(R)]
    xs = [np.zeros((N, M)) for r in range(R)]
    empty = np.array([], dtype=int)
    pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
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
        dates[n] = bary.JDUTC_to_BJDTDB(jd_mid, starname)[0]
        bervs[n] = bary.get_BC_vel(jd_mid, starname=starname, lat=_lat,
                                   longi=_lon, alt=_elevation)[0]  # m/s # why not use dates[n]?
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
            xs[r][n, :] = wave[r, :]

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

    # re-introduce BERVs to HARPS results:
    # pipeline_rvs -= bervs
    # pipeline_rvs -= np.mean(pipeline_rvs)

    return data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts


def savfile_to_filelist(savfile, destination_dir='../data/'):
    # copies CCF + E2DS files to destination_dir and returns a list of the CCFs
    s = readsav(savfile)
    filelist = []
    files = [f.decode('utf8') for f in s.files]
    for f in files:
        shutil.copy2(f, destination_dir)
        spec_file = str.replace(f, 'ccf_G2', 'e2ds')
        shutil.copy2(spec_file, destination_dir)
        basename = f[str.rfind(f, '/') + 1:]
        filelist = np.append(filelist, destination_dir+basename)
    return filelist


def write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs,
               airms, drifts, filenames, hdffile):
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
    filenames = [a.encode('utf8') for a in filenames] # h5py workaround
    dset = h.create_dataset('filelist', data=filenames)
    h.close()


if __name__ == "__main__":

    if False: #YZ Cet VIS
        filelist = glob.glob('/home/jkemmer/Documents/CARMENES/DATA/VIS/raw/YZ_Cet/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="vis", starname=None)
        hdffile = './YZ_Cet_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)

    if False: #YZ Cet NIR
        filelist = glob.glob('/home/jkemmer/Documents/CARMENES/DATA/NIR/YZ_Cet/*sci-gtoc-nir_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="nir", starname=None)
        hdffile = './YZ_Cet_nir_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)

    if True: #HD 119130 VIS
        filelist = glob.glob('/home/jkemmer/Documents/CARMENES/DATA/VIS/raw/HD119130/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="vis", starname='HD 119130')
        hdffile = './HD119130_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)
