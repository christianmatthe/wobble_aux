import h5py
import numpy as np


#BROKEN AND DEPRECATED
#base_dir = "/home/christian/lx39_data/cmatthe/wobble_data/"

#starname = 'GJ1148'
#K_star = 0
#K_t = 3

#results_file_base = base_dir + 'results/results_{0}_Kstar{1}_Kt{2}_n{3}'.format(starname, K_star, K_t, niter)

#results_file = results_file_base + '_stitched' + defaultQ + '.hdf5'

#data_file= base_dir + 'data/' + starname+ '_vis'+'_e2ds.hdf5'

data_file = '/data/cmatthe/wobble_data/data/GJ1148_vis_e2ds.hdf5'
results_file = '/data/cmatthe/wobble_data/results/results_GJ1148_Kstar0_Kt3_stitched_def.hdf5' 

with h5py.File(results_file,'r+') as f:
    #workaround for wobble to include dates_utc, match as in compare_serval_wobble
    with h5py.File(data_file,'r') as g:
        indices_dates_utc = [] 
        indices_dates = []
        for n in range(len(f['dates'][()])):
            ind_jd = np.where(np.abs(g['dates_utc'][()]-f['dates'][n]) == np.nanmin(np.abs(g['dates_utc'][()]-f['dates'][n])))[0][0]
            if (g['dates_utc'][ind_jd]-f['dates'][n])*24*60<20.: #only takes matches closer than 20 minutes
                indices_dates_utc.append(ind_jd)
                indices_dates.append(n)
                
        w_dates_utc = g['dates_utc'][indices_dates_utc]        
        f['dates_utc'] = w_dates_utc
