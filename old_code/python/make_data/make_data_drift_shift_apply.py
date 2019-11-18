from make_data_carmenes_drift_shift import *
import sys
sys.path.append('/data/cmatthe/python/wobble_aux')# path to timported python fiels containing auxiliary functions, this  path only works on lx39
from data_file_split import split_orders_file

if __name__ == "__main__":
    data_directory="/data/cmatthe/wobble_data/data/"
    
    
    
    if True: # GJ876 :vis
        starname = "GJ876"
        arm = "vis"
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/{0}/*sci-gtoc-{1}_A.fits'.format(starname, arm))
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm= arm, starname= None)
        hdffile = data_directory+'{0}_{1}_drift_shift_e2ds.hdf5'.format(starname, arm)
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    
    if False: # GJ1148
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ1148/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'GJ1148_vis_drift_shift_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    
    
    if False: # GJ1148 /J11417+427 :NIR
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ1148/*sci-gtoc-nir_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="nir", starname= None)
        hdffile = data_directory+'GJ1148_nir_drift_shift_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        split_orders_file(hdffile) # save  an aditional split copy
        
        
    if False: ##J11421+267 / GJ436 / Ross 905
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ436/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= None) # try entering starname (catalogue identifier) if simbad can't find star
        hdffile = data_directory+'GJ436_vis_drift_shift_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        
    if False: ##J11421+267 / GJ436 / Ross 905: NIR
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ436/*sci-gtoc-nir_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="nir", starname= None) # try entering starname (catalogue identifier) if simbad can't find star
        hdffile = data_directory+'GJ436_nir_drift_shift_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        split_orders_file(hdffile) # save  an aditional split copy
