from make_data_carmenes import *
import sys
sys.path.append('/data/cmatthe/python/wobble_aux')# path to timported python fiels containing auxiliary functions, this  path only works on lx39
from data_file_split import split_orders_file

if __name__ == "__main__":
    data_directory="/data/cmatthe/wobble_data/data/"
    
    
    if False: # GJ3512 /J08413+594 :vis
        starname = "GJ3512"
        arm = "vis"
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/{0}/*sci-gtoc-{1}_A.fits'.format(starname, arm))
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm= arm, starname= None)
        hdffile = data_directory+'{0}_{1}_e2ds.hdf5'.format(starname, arm)
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    
    if False: # GJ3512 /J08413+594 :NIR
        starname = "GJ3512"
        arm = "nir"
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/{0}/*sci-gtoc-{1}_A.fits'.format(starname, arm))
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm= arm, starname= None)
        hdffile = data_directory+'{0}_{1}_e2ds.hdf5'.format(starname, arm)
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        split_orders_file(hdffile) # save  an aditional split copy
    
    if True: # GJ1148 /J11417+427 :NIR
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ1148/*sci-gtoc-nir_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="nir", starname= None)
        hdffile = data_directory+'GJ1148_nir_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        split_orders_file(hdffile) # save  an aditional split copy
    
    if False: ##J11421+267 / GJ436 / Ross 905: NIR
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ436/*sci-gtoc-nir_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="nir", starname= None) # try entering starname (catalogue identifier) if simbad can't find star
        hdffile = data_directory+'GJ436_nir_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    
    if False: ##J20566+621 / GJ809A (Adrians inkonsitentes Signal) /HD 199305 
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ809A/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= 'HD 199305') # try entering starname (catalogue identifier) if simbad can't find star
        hdffile = data_directory+'GJ809A_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)

    if True: # GJ1148
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ1148/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'GJ1148_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        
    if False: # GJ876
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ876/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'GJ876_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)

    if False: # Wolf294
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/Wolf294/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'Wolf294_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        
    if False: # Teegarden
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/Teegarden/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= 'GAT 1370') # try entering starname (catalogue identifier) if simbad can't find star
        hdffile = data_directory+'Teegarden_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    
    if False: ##J11421+267 / GJ436 / Ross 905
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ436/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= None) # try entering starname (catalogue identifier) if simbad can't find star
        hdffile = data_directory+'GJ436_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
        
    if False: ##J07446+035 / GJ285 / YZCMi
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ285/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc = read_data_from_fits(filelist, arm="vis", starname= 'GJ 285') # try entering starname (catalogue identifier) if simbad can't find star
        hdffile = data_directory+'GJ285_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, dates_utc, filelist, hdffile)
    

### deprecated below

    if False: # Wolf294
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/Wolf294/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'Wolf294_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)
        
    if False: # Barnards
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/Barnard/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'Barnard_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)
    
    if False: #J11421+267 / GJ436 / Ross 905
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ436/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'GJ436_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)
        
    if False: # GJ876
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/GJ876/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="vis", starname= None)
        hdffile = data_directory+'GJ876_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)
        
    if False: #YZ Cet VIS
        filelist = glob.glob('/data/cmatthe/CARM_raw_data/YZ_Cet/*sci-gtoc-vis_A.fits')
        data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts = read_data_from_fits(filelist, arm="vis", starname=None)
        hdffile = data_directory+'YZ_Cet_vis_e2ds.hdf5'
        write_data(data, ivars, xs, pipeline_rvs, pipeline_sigmas, dates, bervs, airms, drifts, filelist, hdffile)
