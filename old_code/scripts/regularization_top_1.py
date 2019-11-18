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

if __name__ == "__main__": 
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
        output_suffix = parameters.dictionary["output_suffix"]
        if defaultQ:
            default_str = "_def"
        else:
            default_str = ""

    else:
        starname = 'GJ876'
        K_star = 0
        K_t = 3
        start = 0 #first order #reg seems to break if not startet at 0, TODO and at any rate wobble would not use the file propperly in that case
        end = 61 # first order that is *not* included
        chunk_size = 42 #use at least 2 so median is not out of bounds, make sure last chunk is not accidentaly 1 long

    data_file= '../data/' + starname+ '_vis'+'_e2ds.hdf5'

    reg_file_base = '../wobble/regularization/{0}_star_K{1}'.format(starname, K_star)# NOTE this is defined separately in top and chunk -> could cause issues
    reg_t_file_base = '../wobble/regularization/{0}_t_K{1}'.format(starname, K_t)
    
    star_filename_final = reg_file_base + '_orders[{0},{1})_'.format(start, end)+ output_suffix + '.hdf5'
    tellurics_filename_final = reg_t_file_base + '_orders[{0},{1})_'.format(start, end)+ output_suffix + '.hdf5'

    ####################


    start_time = time()

    #make list of chunks

    chunks = cr.chunk_list(start, end, chunk_size)    

    # HACK global epoch cutting: use epoch (and order) cutting only once on all data, then pass epochs to chunks
    orders = np.arange(start, end)
    data = wobble.Data(data_file, orders = orders, min_flux=10**-5)
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
        
        os.system("python3 regularization_chunk_1.py")
        print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
        

    print("all chunks optimized: writing combined file")

    ####Star combine
    # initialize with first chunk (save first chunk as stitched file)
    #using the first file as the base will mean assuming (hoping) they all have the same number of epochs
    start_order = chunks[0, 0]
    end_order = chunks[0, 1]
    filename_chunk = cr.reg_file_chunk(start_order, end_order, reg_file_base)
    shutil.copy(filename_chunk, star_filename_final)
    
    if not len(chunks) == 1:
        with h5py.File(star_filename_final,'r+') as f:
            
                        
            #append orders from other chunks
            for i in range(1, len(chunks)):
            # load files and combine them
                start_order = chunks[i, 0]
                end_order = chunks[i, 1]
                filename_chunk = cr.reg_file_chunk(start_order, end_order, reg_file_base)
                with h5py.File(filename_chunk,'r') as g:
                    for key in list(g.keys()):
                        # append new regularization parameters
                        temp = f[key][()]
                        del f[key] # HACK but works 
                        f[key] = np.append(temp, g[key][()])
                        
        ####tellurics combine
        # initialize with first chunk (save first chunk as stitched file)
        #using the first file as the base will mean assuming (hoping) they all have the same number of epochs
        start_order = chunks[0, 0]
        end_order = chunks[0, 1]
        filename_chunk = cr.reg_t_file_chunk(start_order, end_order, reg_t_file_base)
        shutil.copy(filename_chunk, tellurics_filename_final)

        with h5py.File(tellurics_filename_final,'r+') as f:
            
                        
            #append orders from other chunks
            for i in range(1, len(chunks)):
            # load files and combine them
                start_order = chunks[i, 0]
                end_order = chunks[i, 1]
                filename_chunk = cr.reg_t_file_chunk(start_order, end_order, reg_t_file_base) 
                with h5py.File(filename_chunk,'r') as g:
                    for key in list(g.keys()):
                        # append new regularization parameters
                        temp = f[key][()]
                        del f[key] # HACK but works 
                        f[key] = np.append(temp, g[key][()])

    print("regulization files finished")
    print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
