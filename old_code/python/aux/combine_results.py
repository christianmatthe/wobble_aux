import os
from time import time
import wobble
import h5py
import numpy as np
import shutil
import yaml

# WIP
#TODO reduce dependence of variables defined outside of function e.g. default_str, results_file

#def results_file_chunk(start_order, end_order, default_str = default_str):
def results_file_chunk(start_order, end_order, results_file_base, default_str):
    results_file_chunk_name = results_file_base + '_orders[{0},{1})'.format(start_order, end_order) + default_str + '.hdf5' #this name must be same as in chunkscript
    return results_file_chunk_name

def reg_file_chunk(start_order, end_order, reg_file_base):
    reg_file_chunk_name = reg_file_base + '_orders[{0},{1}).hdf5'.format(start_order, end_order) #this name must match name in regularization_chunks
    return reg_file_chunk_name

def reg_t_file_chunk(start_order, end_order, reg_t_file_base):
    reg_t_file_chunk_name = reg_t_file_base + '_orders[{0},{1}).hdf5'.format(start_order, end_order) #this name must match name in regularization_chunks
    return reg_t_file_chunk_name




def chunk_list(start, end, chunk_size):
    chunks=np.asarray([[]])
    #fewer than chunksize orders
    if start + chunk_size >= end:
        chunks = np.asarray([[start, end]])
    else:
        for i in range(start, end, chunk_size):
            if i == start:
                start_order = i
                end_order = i + chunk_size
                chunks = np.asarray([[start_order, end_order]])
            else:
                if i + chunk_size < end:
                    start_order = i
                    end_order = i + chunk_size
                else:
                    start_order = i
                    end_order = end
                chunks = np.append(chunks, [[start_order, end_order]], axis=0) 
    return chunks

def results_file_stitch(start, end, chunk_size, results_file, results_file_base, default_str):
    print("Results: writing to {0}".format(results_file))
    # initialize with first chunk (save first chunk as stitched file)
    #using the first file as the base will mean assuming (hoping) they all have the same number of epochs
    chunks = chunk_list(start, end, chunk_size)
    start_order = chunks[0, 0]
    end_order = chunks[0, 1]
    filename_chunk = results_file_chunk(start_order, end_order, results_file_base, default_str)
    shutil.copy(filename_chunk, results_file)

    # replace number of orders
    with h5py.File(results_file,'r+') as f:
        n_orders = end - start # note this is not the number of optimized orders if any orders  where dropped
        f['R'][()] = n_orders
        print('n_orders = R =',f['R'][()])
                    
        #append orders from other chunks
        for i in range(1, len(chunks)):
        # load files and combine them
            start_order = chunks[i, 0]
            end_order = chunks[i, 1]
            filename_chunk = results_file_chunk(start_order, end_order, results_file_base, default_str)
            with h5py.File(filename_chunk,'r') as g:
                # append new order numbers
                orders = f['orders'][()]
                del f['orders'] # HACK but works 
                f['orders'] = np.append(orders, g['orders'][()])
                
                #check which orders are already there because wobble might drop some
                # TODO clean way is probably to make some neat ordered list of orders in either file first 
                key_list = list(f.keys())
                for r in range(n_orders):
                    if 'order{0}'.format(r) not in key_list:  #find first order index not already present
                        for r_chunk in range(g['R'][()]):
                            r_tot=r+r_chunk
                            #f.create_group('order{0}'.format(r_tot))
                            g.copy('order{0}'.format(r_chunk), f,  name = 'order{0}'.format(r_tot))
                    
                        break
                    
def append_dates_utc(results_file, data_file):
  # epochwise matching of dates_utc
    with h5py.File(results_file,'r+') as f:
    
        #workaround for wobble to include dates_utc, match as in compare_serval_wobble
        with h5py.File(data_file,'r') as g:
            epoch_indices =list( f['epochs'][()])
            w_dates_utc = g['dates_utc'][epoch_indices]        
            f['dates_utc'] = w_dates_utc
            
def append_parameters(results_file, parameter_file):
    with h5py.File(results_file,'r+') as f:
    
        with open(parameter_file, 'r') as stream:
            data_loaded = yaml.safe_load(stream)
            data_dumped = yaml.dump(data_loaded)
            
        data_hdf5 = [data_dumped.encode('utf8')] # h5py workaround
        if 'parameter_dictionary' in list(f.keys()):
            del f['parameter_dictionary']
        f.create_dataset('parameter_dictionary', data = data_hdf5)
            
def load_parameters_hdf5(results_file):#returns parameter dictionary saved as string in hdf5 file via "append_parameters" function
        with h5py.File(results_file,'r+') as f:
            param_yaml = f['parameter_dictionary'][()][0]
            param_yaml = param_yaml.decode('utf8')
            param_yaml = yaml.safe_load(param_yaml)
        return param_yaml
            

if __name__=='__main__':
    
    #################### parameters
    start = 11 #first order 
    end = 53 # first order that is *not* included
    n_orders = end - start
    chunk_size = 5 #use at least 2 so median is not out of bounds #NOTE: currrently needs to match reg fille chunk_size

    niter = 160

    starname = 'GJ876'
    K_star = 0
    K_t = 3

    #wether to use default regularization
    defaultQ = False
    if defaultQ:
        default_str = "_def"
    else:
        default_str = ""

    base_dir = '/data/cmatthe/wobble_data/'

    #results_file_base = base_dir + 'results/results_{0}_Kstar{1}_Kt{2}_n{3}'.format(starname, K_star, K_t, niter)
    #results_file_base = base_dir + 'results/results_{0}_Kstar{1}_Kt{2}'.format(starname, K_star, K_t)

    #results_file = results_file_base + '_stitched' + default_str + '.hdf5'
    kt='Kt3_adrian'
    results_file = base_dir + "results/results_"+ starname +"_Kstar0_"+kt+".hdf5"

    data_file= base_dir + 'data/' + starname+ '_vis'+'_e2ds.hdf5'
    
    wobble_dir = '/data/cmatthe/compare_wobble_serval/wobbledir/'
    wobble_dir_file_name = "results_"+ starname +"_Kstar0_"+kt+ "utc_attached"+".hdf5"
    
    wobble_file = wobble_dir + wobble_dir_file_name

    ####################
    #results_file_stich(start, end, chunk_size, results_file, results_file_base, default_str)
    
    

    #results = wobble.Results(filename = results_file)
    #results.combine_orders('star')
    

    #print("final RVs calculated.")
    ##print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))

    #results.write(results_file) #This apperently deletes the appended dates_utc by overwriting the file (confirmed)
    

    append_dates_utc(results_file, data_file)
    
    #append_parameters(results_file, "/data/cmatthe/wobble_data/scripts/yaml_temp/optimize_parameters.yaml")
    #param_yaml = load_parameters_hdf5(results_file)
    #with open("/data/cmatthe/wobble_data/scripts/yaml_temp/attach_test.yaml", 'w') as outfile:
        #yaml.dump(param_yaml, outfile)
    
    #also write to wobble_dir
    shutil.copy(results_file, wobble_file)

        
    print("all scripts executed")
    #print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
