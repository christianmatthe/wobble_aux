import os
from time import time
import wobble
import h5py
import numpy as np
import shutil

#################### parameters
start = 0 #first order 
end = 60 # first order that is *not* included
n_orders = end - start
chunk_size = 5 #use at least 2 so median is not out of bounds #NOTE: currrently needs to match reg fille chunk_size

niter = 1000

starname = 'GJ876'
K_star = 0
K_t = 3

#wether to use default regularization
defaultQ = True
if defaultQ:
    default_str = "_def"
else:
    default_str = ""


results_file_base = '../results/results_{0}_Kstar{1}_Kt{2}_n{3}'.format(starname, K_star, K_t, niter)

results_file = results_file_base + '_stitched' + default_str + '.hdf5'

data_file= '../data/' + starname+ '_vis'+'_e2ds.hdf5'

####################

def results_file_chunk(start_order, end_order):
    results_file_chunk_name = results_file_base + '_orders[{0},{1})'.format(start_order, end_order) + default_str + '.hdf5' #this name must be same as in script_GJ876_Chunktest
    return results_file_chunk_name

def append_dates_utc(results_file, data_file):
  # epochwise matching of dates_utc
    with h5py.File(results_file,'r+') as f:
    
        #workaround for wobble to include dates_utc, match as in compare_serval_wobble
        with h5py.File(data_file,'r') as g:
            epoch_indices =list( f['epochs'][()])
            w_dates_utc = g['dates_utc'][epoch_indices]        
            f['dates_utc'] = w_dates_utc

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

# HACK global epoch cutting: use epoch (and order) cutting only onnce on all data, then pass epochs to chunks
data = wobble.Data(data_file, orders= np.arange(start, end), min_flux=10**-5)
epochs = np.asarray(data.epochs)
epochs_list_str = str(list(epochs))

for i in range(len(chunks)):
    start_order = chunks[i, 0]
    end_order = chunks[i, 1]
    os.system("python3 script_chunkwise_niter.py {0} {1} {2} {3} {4} \"{5}\" {6} {7}".format(start_order, end_order, starname, K_star, K_t, epochs_list_str, niter, int(defaultQ))) # TODO make this name context dependent/ work for any star params
    #TODO pass this asa parameter file, its getting unwieldy
    print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
    

print("all chunks optimized: writing combined file") #TODO separate this off into another file that can be run separately   

#results_file = results_file_base + '_stitched'+'.hdf5' moved to top
print("Results: writing to {0}".format(results_file))
# initialize with first chunk (save first chunk as stitched file)
#using the first file as the base will mean assuming (hoping) they all have the same number of epochs
start_order = chunks[0, 0]
end_order = chunks[0, 1]
filename_chunk = results_file_chunk(start_order, end_order)
shutil.copy(filename_chunk, results_file)

# replace number of orders
with h5py.File(results_file,'r+') as f:
    f['R'][()] = n_orders
    print('n_orders = R =',f['R'][()])
                
    #append orders from other chunks
    for i in range(1, len(chunks)):
    # load files and combine them
        start_order = chunks[i, 0]
        end_order = chunks[i, 1]
        filename_chunk = results_file_chunk(start_order, end_order)
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

    #workaround for wobble to include dates_utc, match as in compare_serval_wobble BROKEN
    #TODO replace with method from combine_results.py
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


results = wobble.Results(filename = results_file)
results.combine_orders('star')
    
print("final RVs calculated.")
print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))

results.write(results_file)
append_dates_utc(results_file, data_file)# cannot be done before resullts.write, as this will remove dates_utc
#also write to wobble_dir
wobble_dir = '../../compare_wobble_serval/wobbledir/'
wobble_dir_file_name = 'results_{0}_Kstar{1}_Kt{2}_n{3}'.format(starname, K_star, K_t, niter) + '_stitched' + default_str + '.hdf5'
compare_file = wobble_dir + wobble_dir_file_name)
shutil.copy(results_file, compare_file)

    
print("all scripts executed")
print("time elapsed total: {0:.2f} min".format((time() - start_time)/60.0))
