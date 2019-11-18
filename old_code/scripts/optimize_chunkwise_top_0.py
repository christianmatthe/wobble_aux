import sys
sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
import combine_results as cr

import os
from time import time
import wobble
import h5py
import numpy as np
import shutil



#################### parameters
start = 11 #first order 
end = 14 # first order that is *not* included
n_orders = end - start
chunk_size = 2 #use at least 2 so median is not out of bounds #NOTE: currrently needs to match reg fille chunk_size

niter = 300

starname = 'GJ876'
K_star = 0
K_t = 3

#wether to use default regularization
defaultQ = True
if defaultQ:
    default_str = "_def"
else:
    default_str = ""
    
#TODO: add top_directory parameter 


results_file_base = '../results/results_{0}_Kstar{1}_Kt{2}_n{3}'.format(starname, K_star, K_t, niter)

results_file = results_file_base + '_importtest' + default_str + '.hdf5'

data_file= '../data/' + starname+ '_vis'+'_e2ds.hdf5'

wobble_dir = '../../compare_wobble_serval/wobbledir/'
wobble_dir_file_name = 'results_{0}_Kstar{1}_Kt{2}_n{3}'.format(starname, K_star, K_t, niter) + '_stitched' + default_str + '.hdf5'
compare_file = wobble_dir + wobble_dir_file_name

####################


start_time=time()

#make list of chunks

chunks = cr.chunk_list(start, end, chunk_size)    

# HACK global epoch cutting: use epoch (and order) cutting only once on all data, then pass epochs to chunks
data = wobble.Data(data_file, orders= np.arange(start, end), min_flux=10**-5)
epochs = np.asarray(data.epochs)
epochs_list_str = str(list(epochs))

for i in range(len(chunks)):
    start_order = chunks[i, 0]
    end_order = chunks[i, 1]
    os.system("python3 script_chunkwise_niter.py {0} {1} {2} {3} {4} \"{5}\" {6} {7}".format(start_order, end_order, starname, K_star, K_t, epochs_list_str, niter, int(defaultQ)))
    #TODO pass this asa parameter file, its getting unwieldy
    print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
    

print("all chunks optimized: writing combined file") 
cr.results_file_stitch(start, end, chunk_size, results_file, results_file_base, default_str)


results = wobble.Results(filename = results_file)
results.combine_orders('star')
    
print("final RVs calculated.")
print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))

results.write(results_file)
cr.append_dates_utc(results_file, data_file)# cannot be done before results.write, as this will remove dates_utc

#also write to wobble_dir
shutil.copy(results_file, compare_file)

    
print("all scripts executed")
print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
