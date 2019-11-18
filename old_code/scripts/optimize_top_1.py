import sys
sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
import combine_results as cr
from parameters import Parameters
#end local imports
import os
from time import time
import wobble
import h5py
import numpy as np
import shutil
import yaml

#TODO pass results_file_base to chunk
if __name__ == "__main__": #NOTE If called via os.system this will trigger even if in queue
    queue = True
    #################### parameters
    if queue:
        parameters = Parameters(filename = "yaml_temp/optimize_parameters.yaml")
        
        starname = parameters.dictionary["starname"]
        K_star = parameters.dictionary["K_star"]
        K_t = parameters.dictionary["K_t"]
        niter = parameters.dictionary["niter"]
        start = parameters.dictionary["start"]
        end = parameters.dictionary["end"]
        chunk_size = parameters.dictionary["chunk_size"]
        defaultQ = parameters.dictionary["defaultQ"]
        if defaultQ:
            default_str = "_def"
        else:
            default_str = ""
    else:
        starname = 'GJ876'
        K_star = 0
        K_t = 3
        niter = 300
        start = 11 #first order 
        end = 14 # first order that is *not* included
        chunk_size = 2 #use at least 2 so median is not out of bounds #NOTE: currrently needs to match reg file chunk_size
        #wether to use default regularization
        defaultQ = True
        if defaultQ:
            default_str = "_def"
        else:
            default_str = ""
        
        
    #TODO: add top_directory parameter 


    results_file_base = '../results/results_{0}_Kstar{1}_Kt{2}_'.format(starname, K_star, K_t, niter) + parameters.dictionary["output_suffix"] 

    results_file = results_file_base + parameters.dictionary["output_suffix"] + default_str + '.hdf5'

    data_file= '../data/' + starname+ '_vis'+'_e2ds.hdf5'

    wobble_dir = '../../compare_wobble_serval/wobbledir/'
    wobble_dir_file_name = 'results_{0}_Kstar{1}_Kt{2}_'.format(starname, K_star, K_t, niter) + parameters.dictionary["output_suffix"] +'.hdf5'
    compare_file = wobble_dir + wobble_dir_file_name

    ####################


start_time = time()

#make list of chunks

chunks = cr.chunk_list(start, end, chunk_size)    

# HACK global epoch cutting: use epoch (and order) cutting only once on all data, then pass epochs to chunks
data = wobble.Data(data_file, orders= np.arange(start, end), min_flux=10**-5)
epochs = np.asarray(data.epochs)
epochs_list = data.epochs.tolist()
epochs_list_str = str(list(epochs))

for i in range(len(chunks)):
    start_chunk = int(chunks[i, 0])
    end_chunk = int(chunks[i, 1])
    
    parameters.dictionary.update({
        "epochs_list" : epochs_list,
        "start_chunk" : start_chunk,
        "end_chunk" : end_chunk
        })
        
    #Write into YAML to pass to chunk script
    parameters.write("yaml_temp/optimize_parameters_chunk.yaml")
    
    
    os.system("python3 optimize_chunk_1.py")
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
