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
#NOTE WORKS only in Wob_env_2 (wobble_19_03_2019, deprecated in wobble_14_06_2019) 
#TODO pass results_file_base to chunk
if __name__ == "__main__": #NOTE If called via os.system this will trigger even if in queue
    #################### parameters
    
    parameter_filename = "yaml_temp/optimize_parameters.yaml"
    parameter_chunk_filename = "yaml_temp/optimize_parameters_chunk.yaml"
    parameters = Parameters(filename = parameter_filename)
    
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
        
    parameter_dict = parameters.dictionary    
    if "reg_file_star" in parameter_dict and "reg_file_t" in parameter_dict:
        star_reg_file = parameter_dict["reg_file_star"]
        tellurics_reg_file = parameter_dict["reg_file_t"]
        
    top_directory = parameter_dict["top_directory"]
    run_directory = parameter_dict["run_directory"]
    data_directory = parameter_dict["data_directory"]
    
    results_directory = parameter_dict["results_directory"]
        
    ####################

    results_file_base = results_directory + 'results_{0}_Kstar{1}_Kt{2}_'.format(starname, K_star, K_t, niter) + parameters.dictionary["output_suffix"] 

    results_file = results_file_base + default_str + '.hdf5'
    parameters.dictionary.update({"latest_results_file" : results_file})
    parameters.write(parameter_filename)

    data_file= data_directory + starname+ '_vis'+'_e2ds.hdf5'

    #wobble_dir = top_directory + '../compare_wobble_serval/wobbledir/'
    #wobble_dir_file_name = 'results_{0}_Kstar{1}_Kt{2}_'.format(starname, K_star, K_t, niter) + parameters.dictionary["output_suffix"] +'.hdf5'
    #compare_file = wobble_dir + wobble_dir_file_name

    ####################


start_time = time()

#make list of chunks

chunks = cr.chunk_list(start, end, chunk_size) 


# HACK global epoch cutting: use epoch (and order) cutting only once on all data, then pass epochs to chunks
data = wobble.Data(data_file, orders= np.arange(start, end), min_flux=10**-5)
epochs = np.asarray(data.epochs)
epochs_list = data.epochs.tolist()
epochs_list_str = str(list(epochs))
# 3.1: introduce custom epoch list from preselected good epochs
if "custom_epochs" in parameters.dictionary.keys():
    epochs_list = list(set(epochs_list).intersection(parameters.dictionary["custom_epochs"]))

##create file containing actually used regularization parameters forfuture reference
#reg_log_file = "/data/cmatthe/wobble_data/results/reg_parameters_record/reg_par_" + parameters.dictionary["output_suffix"] + ".txt"
#parameters.dictionary.update({
        #"reg_log_file" : reg_log_file,
        #})
#f = open(reg_log_file,"w+")
#f.write("order, L1_template, L2_template, L1_basis_vector, L2_basis_vector, L2_template\n")
#f.close()

for i in range(len(chunks)):
    start_chunk = int(chunks[i, 0])
    end_chunk = int(chunks[i, 1])
    
    #Create chunk reg files: #NOTE currently with customm reg file only, not default safe
    #NOTE assumes reg file starts at order 0 TODO implenent check that file is 61 orders long
    star_reg_file_chunk = results_directory + '/regularization/temp_star_chunk.hdf5'
    with h5py.File(star_reg_file,'r') as f:
        with h5py.File(star_reg_file_chunk,'w') as g:
            for key in list(f.keys()):
                    temp = f[key][()][start_chunk : end_chunk]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
    tellurics_reg_file_chunk = results_directory + '/regularization/temp_t_chunk.hdf5'
    with h5py.File(tellurics_reg_file,'r') as f:
        with h5py.File(tellurics_reg_file_chunk,'w') as g:
            for key in list(f.keys()):
                    temp = f[key][()][start_chunk : end_chunk]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
            
    
    parameters.dictionary.update({
        "epochs_list" : epochs_list,
        "start_chunk" : start_chunk,
        "end_chunk" : end_chunk,
        
        "reg_file_star_chunk" : star_reg_file_chunk,
        "reg_file_t_chunk" : tellurics_reg_file_chunk
        })
        
    #Write into YAML to pass to chunk script
    parameters.write(parameter_chunk_filename)
    
    
    
    
    os.system("python3 optimize_chunk_3.2.1_reg_mask.py")
    print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
    

print("all chunks optimized: writing combined file") 
cr.results_file_stitch(start, end, chunk_size, results_file, results_file_base, default_str)


results = wobble.Results(filename = results_file)
results.combine_orders('star')
    
print("final RVs calculated.")
print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))

results.write(results_file)
cr.append_dates_utc(results_file, data_file)# cannot be done before results.write, as this will remove dates_utc
cr.append_parameters(results_file, parameter_filename) # appends "parameter_dictionary" as string dataset to the .hdf5 for later use. Load with "cr.load_parameters_hdf5"

#also write to wobble_dir
#shutil.copy(results_file, compare_file)

    
print("all scripts executed")
print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
