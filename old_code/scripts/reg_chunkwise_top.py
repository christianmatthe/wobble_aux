import os
from time import time
import wobble
import h5py
import numpy as np
import shutil

#################### parameters
start = 11 #first order #reg seems to break if not startet at 0, TODO and at any rate wobble would not use the file propperly in that case
end = 60 # first order that is *not* included
n_orders = end - start
chunk_size = 42 #use at least 2 so median is not out of bounds, make sure last chunk is not accidentaly 1 long

starname = 'GJ876'
K_star = 0
K_t = 3
reg_file_base = '../wobble/regularization/{0}_star_K{1}'.format(starname, K_star)
reg_t_file_base = '../wobble/regularization/{0}_t_K{1}'.format(starname, K_t)

####################

def reg_file_chunk(start_order, end_order):
    reg_file_chunk_name = reg_file_base + '_orders[{0},{1}).hdf5'.format(start_order, end_order) #this name must match name in regularization_chunks
    return reg_file_chunk_name

def reg_t_file_chunk(start_order, end_order):
    reg_t_file_chunk_name = reg_t_file_base + '_orders[{0},{1}).hdf5'.format(start_order, end_order) #this name must match name in regularization_chunks
    return reg_t_file_chunk_name


start_time=time()

#make list of chunks

chunks=np.asarray([[]])
for i in range(start, end, chunk_size):
    if i == start:
        start_order = i
        end_order = i + chunk_size
        chunks = [[start_order, end_order]]
    else:
        if i + chunk_size < end:
            start_order = i
            end_order = i + chunk_size
        else:
            start_order = i
            end_order = end
        chunks = np.append(chunks, [[start_order, end_order]], axis=0)    

'''
# HACK global epoch cutting: use epoch (and order) cutting only onnce on all data, then pass epochs to chunks
data = wobble.Data(starname+'_vis'+'_e2ds.hdf5', filepath='../data/', orders= np.arange(start, end), min_flux=10**-5)
epochs = np.asarray(data.epochs)
epochs_list_str = str(list(epochs))
'''

for i in range(len(chunks)):
    start_order = chunks[i, 0]
    end_order = chunks[i, 1]
    os.system("python3 regularization_carm_chunks.py {0} {1} {2} {3} {4}".format(start_order, end_order, starname, K_star, K_t))
    #os.system("python3 script_chunkwise.py {0} {1} {2} {3} {4} \"{5}\"".format(start_order, end_order, starname, K_star, K_t, epochs_list_str)) #version with epoch list (for optimize not for reg)
    #TODO pass this asa parameter file, its getting unwieldy
    print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
    

print("all chunks optimized: writing combined file") #TODO separate this off into another file that can be run separately 

####Star combine
# initialize with first chunk (save first chunk as stitched file)
#using the first file as the base will mean assuming (hoping) they all have the same number of epochs
star_filename_final = reg_file_base + '_orders[{0},{1})_stitched.hdf5'.format(start, end)
start_order = chunks[0, 0]
end_order = chunks[0, 1]
filename_chunk = reg_file_chunk(start_order, end_order)
shutil.copy(filename_chunk, star_filename_final)

with h5py.File(star_filename_final,'r+') as f:
    
                
    #append orders from other chunks
    for i in range(1, len(chunks)):
    # load files and combine them
        start_order = chunks[i, 0]
        end_order = chunks[i, 1]
        filename_chunk = reg_file_chunk(start_order, end_order) 
        with h5py.File(filename_chunk,'r') as g:
            for key in list(g.keys()):
                # append new regularization parameters
                temp = f[key][()]
                del f[key] # HACK but works 
                f[key] = np.append(temp, g[key][()])
                
####tellurics combine
# initialize with first chunk (save first chunk as stitched file)
#using the first file as the base will mean assuming (hoping) they all have the same number of epochs
tellurics_filename_final = reg_t_file_base + '_orders[{0},{1})_stitched.hdf5'.format(start, end)
start_order = chunks[0, 0]
end_order = chunks[0, 1]
filename_chunk = reg_t_file_chunk(start_order, end_order)
shutil.copy(filename_chunk, tellurics_filename_final)

with h5py.File(tellurics_filename_final,'r+') as f:
    
                
    #append orders from other chunks
    for i in range(1, len(chunks)):
    # load files and combine them
        start_order = chunks[i, 0]
        end_order = chunks[i, 1]
        filename_chunk = reg_t_file_chunk(start_order, end_order) 
        with h5py.File(filename_chunk,'r') as g:
            for key in list(g.keys()):
                # append new regularization parameters
                temp = f[key][()]
                del f[key] # HACK but works 
                f[key] = np.append(temp, g[key][()])

print("regulization files finished")
print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
