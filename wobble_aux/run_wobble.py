import numpy as np
import shutil
from time import time
import h5py
import wobble
import os
import matplotlib.pyplot as plt
import dill

#Reproduce combine results first: functions required for runnign wobble in small chunks of orders so as not to overflow RAM
def file_chunk_name(start_order, end_order, chunk_dir):
    results_file_chunk_name = chunk_dir + 'orders[{0},{1}]'.format(start_order, end_order) + '.hdf5' #this name must be same as in chunkscript
    return results_file_chunk_name

def chunk_list(start, end, chunk_size, order_list = None):
    #make chunks only from order list
    if order_list is not None:
        def chunks_from_list(lst, n):
            """Yield successive n-sized chunks from lst."""
            for i in range(0, len(lst), n):
                yield lst[i:i + n]
        chunks = list(chunks_from_list(order_list, chunk_size))
    else:
        raise Exception("order_list must be provided to chunk_list()")
        #legacy 
        #chunks=np.asarray([[]])
        ##fewer than chunksize orders
        #if start + chunk_size >= end:
            #chunks = np.asarray([[start, end]])
        #else:
            #for i in range(start, end, chunk_size):
                #if i == start:
                    #start_order = i
                    #end_order = i + chunk_size
                    #chunks = np.asarray([[start_order, end_order]])
                #else:
                    #if i + chunk_size < end:
                        #start_order = i
                        #end_order = i + chunk_size
                    #else:
                        #start_order = i
                        #end_order = end
                    #chunks = np.append(chunks, [[start_order, end_order]], axis=0) 
    return chunks

#def results_file_stitch(start, end, chunk_size, results_file, chunk_dir):
    #"""turns chunks into one continuous file, as if wobble had been run in one piece"""
    #print("Writing results to {0}".format(results_file))
    ##initialize results_file as copy of first chunk_list
    ## NOTE files must have the same epochs
    #chunks = chunk_list(start, end, chunk_size)
    #start_order = chunks[0, 0]
    #end_order = chunks[0, 1]
    #filename_chunk = file_chunk_name(start_order, end_order, chunk_dir)
    #shutil.copy(filename_chunk, results_file)
    
    ##replace number of orders since chunk will not have all of them
    #with h5py.File(results_file, 'r+') as f:
        ##n_orders = end-start  # note this is not the number of optimized orders if any orders  where dropped
        ##f['R'][()] = n_orders
        
        ##append orders from other chunks
        #if len(chunks) > 1:
            #for i in range(1, len(chunks)):
                ##load files and combine them
                #start_order = chunks[i, 0]
                #end_order = chunks[i, 1]
                #filename_chunk = file_chunk_name(start_order, end_order, chunk_dir)
                #with h5py.File(filename_chunk,'r') as g:
                    ##update total order number
                    #f['R'][()] = f['R'][()] + g['R'][()]
                    #n_orders = f['R'][()]
                    ## append new order numbers
                    #orders = f['orders'][()]
                    #del  f['orders'] #workaround to edit dimensions of entry in h5py file
                    #f['orders'] = np.append(orders, g['orders'][()])
                    #key_list = list(f.keys())
                    #for r in range(n_orders):
                        #if 'order{0}'.format(r) not in key_list: #find first intex not already present
                            #for r_chunk in range(g['R'][()]):
                                #r_tot  = r + r_chunk
                                #g.copy('order{0}'.format(r_chunk), f, name = 'order{0}'.format(r_tot))
                                
                            #break
                        
def results_file_stitch(chunks, results_file, chunk_dir):
    """turns chunks into one continuous file, as if wobble had been run in one piece"""
    print("Writing results to {0}".format(results_file))
    #initialize results_file as copy of first chunk_list
    # NOTE files must have the same epochs
    start_order = chunks[0][0]
    end_order = chunks[0][-1]
    filename_chunk = file_chunk_name(start_order, end_order, chunk_dir)
    shutil.copy(filename_chunk, results_file)
    
    #replace number of orders since chunk will not have all of them
    with h5py.File(results_file, 'r+') as f:
        #n_orders = end-start  # note this is not the number of optimized orders if any orders  where dropped
        #f['R'][()] = n_orders
        
        #append orders from other chunks
        if len(chunks) > 1:
            for i in range(1, len(chunks)):
                #load files and combine them
                start_order = chunks[i][0]
                end_order = chunks[i][-1]
                filename_chunk = file_chunk_name(start_order, end_order, chunk_dir)
                with h5py.File(filename_chunk,'r') as g:
                    #update total order number
                    f['R'][()] = f['R'][()] + g['R'][()]
                    n_orders = f['R'][()]
                    # append new order numbers
                    orders = f['orders'][()]
                    del  f['orders'] #workaround to edit dimensions of entry in h5py file
                    f['orders'] = np.append(orders, g['orders'][()])
                    key_list = list(f.keys())
                    for r in range(n_orders):
                        if 'order{0}'.format(r) not in key_list: #find first intex not already present
                            for r_chunk in range(g['R'][()]):
                                r_tot  = r + r_chunk
                                g.copy('order{0}'.format(r_chunk), f, name = 'order{0}'.format(r_tot))
                                
                            break

def append_dates_utc(results_file, data_file):
  # epochwise matching of dates_utc
  # workaround for wobble to include dates_utc, to match with SERVAL dates
    with h5py.File(results_file,'r+') as f:
        with h5py.File(data_file,'r') as g:
            epoch_indices =list( f['epochs'][()])
            w_dates_utc = g['dates_utc'][epoch_indices]        
            f['dates_utc'] = w_dates_utc
            
def append_parameters(parameters, results_file):
    #append pickled binary version of parameter object to results file
    #TODO also attach yaml for manual lookup
    pkl_string_parameters = dill.dumps(parameters)
    with h5py.File(results_file,'r+') as f:
        dset = f.create_dataset("parameters",(100,))
        dset.attrs["pkl"] = np.void(pkl_string_parameters)
        
def read_parameters_from_results(results_file):
    with h5py.File(results_file,'r') as f:
        out = f["parameters"].attrs["pkl"]
        pkl_string_parameters = out.tostring()
        parameters = dill.loads(pkl_string_parameters)
        return parameters
        
        
        
            
#TODO Parameter functions 

class Parameters: 
    """
    The Parameter object: contains the parameters specifying how to run wobble.
    
    Parameters
    ----------
    starname : `str`
        Name of the star to be worked on.
    K_star : `int` (default `0`)
        Number of components used to model the stellar spectrum.
    K_t : `int` (default `3`)
        Number of components used to model the telluric spectrum.
    niter : `int` (default `160`)
        number of iterations to be used when optimizing.
    start : `int` (default `11`)
        number of first order that will be optimized.
    end : `int` (default `53`)
        number of first order that will *not* be optimized.
    chunk_size : `int` (default `5`)
        Number of orders that will be optimized in one chunk.
    reg_file_star: `.hdf5 file` (default `None`)
        star regularization file to be used by wobble
    reg_file_t: `.hdf5 file` (default `None`)
        tellurics regularization file to be used by wobble
    output_suffix : `str` (default ``)
        string appended to output filenames. Serves as distinguishing name when rerunning the same dataset with different parameters
    data_suffix : `str` (default ``)
        string appended to data file to be imported. Allows for selecting different data files for one starname
    results_dir : `str` (default `../results/`)
        path to directory in which results will be saved
    data_dir : `str` (default `../data/`)
        path to directory from which data is loaded
    min_snr : `int` (default `60`)
        epochs with average SNR over range of optimized orders below this  number are dropped (5 is wobble default but this is too low in many cases)
    plots : `bool` (default `True`)
        whether or not plots of the generated synthetic spectra are saved
    continuum_order : `int` (default `1`)
        (wobble default is 6)
        order of the polynomial used for normalizing the spectra
    continuum_nsigma : [`float`, `float`] (default `[0.5,1]`)
        data points this many sigma deviations [down, up] from the continuum fit are cut every step during continuum normalization
    plot_continuum : `bool` (default `False`)
        whether or not to output plots of the continua during data import (Note: a LOT of plots (orders*epochs))
        
    
        
    
    """
    def __init__(self,
                 starname = None, #Required
                 K_star = 0,
                 K_t = 3,
                 niter = 160,
                 start = 11,
                 end = 53,
                 chunk_size = 5,
                 reg_file_star = None, #Required
                 reg_file_t = None,   #Required
                 output_suffix = "",
                 data_suffix = "",
                 results_dir = '../results/',
                 data_dir= '../data/',
                 min_snr = 60,
                 plots = True,
                 continuum_order = 1,
                 continuum_nsigma = [0.5,1],
                 plot_continuum  = False
                 ):
        self.starname = starname
        self.K_star = K_star
        self.K_t = K_t
        self.niter = niter
        self.start = start
        self.end = end
        self.chunk_size = chunk_size
        self.reg_file_star = reg_file_star
        self.reg_file_t = reg_file_t
        self.output_suffix = output_suffix
        self.data_suffix = data_suffix
        self.results_dir = results_dir
        self.data_dir= data_dir
        self.min_snr = min_snr
        self.plots = plots
        self.continuum_order = continuum_order
        self.continuum_nsigma = continuum_nsigma
        self.plot_continuum  = plot_continuum
        '''
        self.dictionary = {
            "starname" : starname,
            "K_star" : K_star,
            "K_t" : K_t,
            "niter" : niter,
            "start" : start,
            "end" : end,
            "chunk_size" : chunk_size,
            "reg_file_star" : reg_file_star
            "reg_file_t" : reg_file_t
            "output_suffix" : output_suffix,
            "data_suffix" : data_suffix
            }
            '''
        if starname is None:
            raise Exception("Must specify which star to run with the 'starname' keyword in parameters")
        if reg_file_star is None or reg_file_t is None:
            raise Exception("Must specify both reg_file_star and reg_file_t")


def reg_chunk(chunk, reg_file_star, reg_file_t):
    """Create subsets of the regularization files to be used for a chunk"""
    #start_chunk = int(chunk[0])
    #end_chunk = int(chunk[-1])
    #NOTE assumes reg file starts at order 0 TODO implenent check that file is 61 orders long i.e. that reg file is valid
    #TODO Remove hardcoded temp location?
    reg_file_star_chunk = os.path.dirname(os.path.abspath(__file__)) + '/regularization/temp_star_chunk.hdf5'
    with h5py.File(reg_file_star,'r') as f:
        with h5py.File(reg_file_star_chunk,'w') as g:
            for key in list(f.keys()):
                    #temp = f[key][()][start_chunk : end_chunk]
                    temp = f[key][()][chunk]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
    reg_file_t_chunk = os.path.dirname(os.path.abspath(__file__)) +'/regularization/temp_t_chunk.hdf5'
    with h5py.File(reg_file_t,'r') as f:
        with h5py.File(reg_file_t_chunk,'w') as g:
            for key in list(f.keys()):
                    #temp = f[key][()][start_chunk : end_chunk]
                    temp = f[key][()][chunk]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
    return reg_file_star_chunk, reg_file_t_chunk

def run_wobble(parameters):
    p = parameters
    
    results_name = 'results_{0}_Kstar{1}_Kt{2}_'.format(p.starname, p.K_star, p.K_t, p.niter) + p.output_suffix
    results_file_base = p.results_dir + results_name
    results_file = results_file_base + '.hdf5'
    data_file = p.data_file = p.data_dir + p.starname  + p.data_suffix + '_e2ds.hdf5'
    
    temp_dir = p.temp_dir = p.results_dir + '/temp_' + results_name + '/'
    plot_dir = p.plot_dir = p.results_dir + '/plots_' + results_name + '/'
    #make (output) directory
    #TODO these data and results dirs should be handlesd somewhere else
    os.makedirs(p.results_dir, exist_ok = True)
    os.makedirs(p.data_dir, exist_ok = True)
    
    os.makedirs(temp_dir, exist_ok = True)
    os.makedirs(plot_dir, exist_ok = True)
    
    start_time = p.start_time = time()
    #generate epoch list
    # if parameters has been passed a global epochs list use this going forward. This is primarily used by regularization.py
    try:
        epochs_list = p.epochs_list = p.global_epochs_list
    except AttributeError:
        print("Loading data. May take a few minutes")
        data = wobble.Data(data_file, orders = np.arange(p.start, p.end), min_flux=10**-5, min_snr = p.min_snr,
                            parameters = p
                            )
        epochs_list = p.epochs_list = data.epochs.tolist()
    #orders_list = p.orders_list = data.orders.tolist() #too agressive in visible # TODO implement alternate nir and vis handling 
    orders_list = p.orders_list = np.arange(p.start, p.end).tolist()
    
    
    chunks = p.chunks = chunk_list(p.start, p.end, p.chunk_size, p.orders_list)
    #Loop over chunks
    for i in range(len(chunks)):
        #pass parameters object to chunk script
        p.i = i
        with open(os.path.dirname(os.path.abspath(__file__)) + "/carmenes_aux_files/chunk_parameters.pkl", "wb") as f:
            dill.dump(p, f)
        #start chunk script
        os.system("python3  {0}/chunk.py".format(os.path.dirname(os.path.abspath(__file__))))
        
    print("all chunks optimized: writing combined file") 
    results_file_stitch(chunks, results_file, temp_dir)
    
    #Combine orders
    results = wobble.Results(filename = results_file)
    results.combine_orders('star')
    print("final RVs calculated.")
    print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
    results.write(results_file)
    append_dates_utc(results_file, data_file)# cannot be done before results.write, as .write will remove dates_utc
    append_parameters(p, results_file)
    '''
    # test
    loaded_parameters = read_parameters_from_results(results_file)
    attrs = dir(loaded_parameters)
    print(attrs)
    for attr in attrs:
        #print(attr, getattr(loaded_parameters, attr))
        q = getattr(p, attr) == getattr(loaded_parameters, attr)
        if q == True:
            print(q)
        else:
            print(attr, q)
    '''
    print("results saved as: {0}".format(results_file))
    print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
    
    #delete temp_dir which at this point only contains duplicates
    shutil.rmtree(temp_dir)
    print("deleted: {0}".format(temp_dir))
    
    
    
    


if __name__ == "__main__":
    #laptop test example
    parameters = Parameters(starname = "GJ436",
                            data_suffix = "_vis_drift+nzp",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            output_suffix = "laptop_example_0",
                            plot_continuum = False)
    run_wobble(parameters)
    
    parameters = Parameters(starname = "GJ3473",
                            data_suffix = "_vis_drift+nzp",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            output_suffix = "laptop_example_1",
                            plot_continuum = False)
    run_wobble(parameters)
    
    '''
    # Test append and read parameters to file
    parameters = Parameters(starname = "GJ1148",
                            data_suffix = "_vis_drift+nzp",
                            start = 30,
                            end = 31,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_0/next_base_star_reg.hdf5',
                            reg_file_t = '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_0/next_base_t_reg.hdf5',
                            output_suffix = "param_pkl_test",
                            plot_continuum = False)
    run_wobble(parameters)
    '''
    
    # run GJ1148 with different regs based on the no_reg seed on lx39 with new NZP corrected data
    #GJ1148
    #for i in range(0,6):
        #parameters = Parameters(starname = "GJ1148",
                                #data_suffix = "_vis_drift+nzp",
                                #start = 11,
                                #end = 53,
                                #chunk_size = 5,
                                #niter = 160,
                                #reg_file_star =  '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_{0}/next_base_star_reg.hdf5'.format(i),
                                #reg_file_t = '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_{0}/next_base_t_reg.hdf5'.format(i),
                                #output_suffix = "GJ1148_no_reg_seed_reg_loop_{0}_with nzp".format(i),
                                #plot_continuum = False)
        #run_wobble(parameters)
    
    ## run GJ436 with different regs based on the no_reg seed on lx39
    ##GJ436
    #for i in range(0,6):
        #parameters = Parameters(starname = "GJ436",
                                #data_suffix = "_vis_drift_shift",
                                #start = 11,
                                #end = 53,
                                #chunk_size = 5,
                                #niter = 160,
                                #reg_file_star =  '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_{0}/next_base_star_reg.hdf5'.format(i),
                                #reg_file_t = '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_{0}/next_base_t_reg.hdf5'.format(i),
                                #output_suffix = "GJ1148__no_reg_seed_reg_loop_{0}".format(i),
                                #plot_continuum = False)
        #run_wobble(parameters)
    
    ## run GJ1148 with different regs based on the no_reg seed on lx39
    ##GJ1148
    #for i in range(3,6):
        #parameters = Parameters(starname = "GJ1148",
                                #data_suffix = "_vis_drift_shift",
                                #start = 11,
                                #end = 53,
                                #chunk_size = 5,
                                #niter = 160,
                                #reg_file_star =  '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_{0}/next_base_star_reg.hdf5'.format(i),
                                #reg_file_t = '/data/cmatthe/wobble_aux/results//regularization/GJ1148_no_reg_seed/loop_{0}/next_base_t_reg.hdf5'.format(i),
                                #output_suffix = "GJ1148__no_reg_seed_reg_loop_{0}".format(i),
                                #plot_continuum = False)
        #run_wobble(parameters)
    
    
    #Barnard with dummy reg
    #parameters = Parameters(starname = "Barnard",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #output_suffix = "no_reg",
                            #plot_continuum = False)
    #run_wobble(parameters)
    
    #GJ436 with dummy reg
    #parameters = Parameters(starname = "GJ436",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #output_suffix = "no_reg",
                            #plot_continuum = True)
    #run_wobble(parameters)
    
    ########### check flat reg loop4 vs loop3
    #GJ1148
    #parameters = Parameters(starname = "GJ1148",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #output_suffix = "no_reg_all orders",
                            #plot_continuum = True)
    #run_wobble(parameters)
    
    # run Teegarden IR with different regs on lx39
    #Teegarden
    #for i in range(6):
        #parameters = Parameters(starname = "Teegarden",
                                #data_suffix = "_nir_drift_shift_split",
                                #start = 0,
                                #end = 56,
                                #chunk_size = 5,
                                #niter = 160,
                                #reg_file_star =  '/data/cmatthe/wobble_aux/results//regularization/GJ1148_all_orders/loop_{0}/next_base_star_reg.hdf5'.format(i),
                                #reg_file_t = '/data/cmatthe/wobble_aux/results//regularization/GJ1148_all_orders/loop_{0}/next_base_t_reg.hdf5'.format(i),
                                #output_suffix = "GJ1148_all_orders_reg_loop_{0}".format(i),
                                #plot_continuum = False)
        #run_wobble(parameters)
    
    # run GJ1148 IR with different regs on lx39
    #GJ1148
    #for i in range(6):
        #parameters = Parameters(starname = "GJ1148",
                                #data_suffix = "_nir_drift_shift_split",
                                #start = 0,
                                #end = 56,
                                #chunk_size = 5,
                                #niter = 160,
                                #reg_file_star =  '/data/cmatthe/wobble_aux/results//regularization/GJ1148_all_orders/loop_{0}/next_base_star_reg.hdf5'.format(i),
                                #reg_file_t = '/data/cmatthe/wobble_aux/results//regularization/GJ1148_all_orders/loop_{0}/next_base_t_reg.hdf5'.format(i),
                                #output_suffix = "GJ1148_all_orders_reg_loop_{0}".format(i),
                                #plot_continuum = False)
        #run_wobble(parameters)
    
    """
    ########### check flat reg loop4 vs loop3
    #GJ1148
    parameters = Parameters(starname = "GJ1148",
                            data_suffix = "_vis_drift_shift",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/flat_reg_star.hdf5',
                            reg_file_t = 'regularization/flat_reg_t.hdf5',
                            output_suffix = "flat_l4",
                            plot_continuum = True)
    run_wobble(parameters)
    
    parameters = Parameters(starname = "GJ1148",
                            data_suffix = "_vis_drift_shift",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/flat_reg_l3_star.hdf5',
                            reg_file_t = 'regularization/flat_reg_l3_t.hdf5',
                            output_suffix = "flat_l3",
                            plot_continuum = True)
    run_wobble(parameters)
    #GJ436
    parameters = Parameters(starname = "GJ436",
                            data_suffix = "_vis_drift_shift",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/flat_reg_star.hdf5',
                            reg_file_t = 'regularization/flat_reg_t.hdf5',
                            output_suffix = "flat_l4",
                            plot_continuum = True)
    run_wobble(parameters)
    
    parameters = Parameters(starname = "GJ436",
                            data_suffix = "_vis_drift_shift",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/flat_reg_l3_star.hdf5',
                            reg_file_t = 'regularization/flat_reg_l3_t.hdf5',
                            output_suffix = "flat_l3",
                            plot_continuum = True)
    run_wobble(parameters)
    """
    
    ###########
    #parameters = Parameters(starname = "GJ1148",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            #reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            #output_suffix = "recheck_all_orders",
                            #plot_continuum = True)
    
    #parameters = Parameters(starname = "GJ1148",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            #reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            #output_suffix = "continuum_recheck",
                            #plot_continuum = True)
    
    #parameters = Parameters(starname = "GJ1148",
                            #data_suffix = "_vis_drift_shift",
                            #start = 30,
                            #end = 34,
                            #chunk_size = 1,
                            #niter = 160,
                            #reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            #reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            #output_suffix = "quickcheck",
                            #plot_continuum = False)
    
    #parameters = Parameters(starname = "GJ3473",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 20,
                            #niter = 160,
                            #reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            #reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            #min_snr = 10,
                            #output_suffix = "l4_reg_snr_10",
                            #plot_continuum = False)
    
    #parameters = Parameters(starname = "GJ3473",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            #reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            #output_suffix = "l4_reg",
                            #plot_continuum = True)

    #parameters = Parameters(starname = "GJ436",
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            #reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            #output_suffix = "git_run_wobble_test0")
                            
    #quick example : use as standart test?
    #parameters = Parameters(starname = "GJ436",
                            #data_suffix = "_vis_drift_shift",
                            #start = 30,
                            #end = 34,
                            #chunk_size = 2,
                            #niter = 160,
                            #reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            #reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            #output_suffix = "test_2_continua",
                            #plot_continuum = False
                            #)
    
    #NIR test
    #parameters = Parameters(starname = "GJ1148",
                        #data_suffix = "_nir_drift_shift_split",
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  'regularization/flat_reg_star.hdf5',
                        #reg_file_t = 'regularization/flat_reg_t.hdf5',
                        #output_suffix = 'test_nir_continuum_0.5down_1up',
                        #plot_continuum = True
                        #)
    
    #parameters = Parameters(starname = "GJ1148",
                        #data_suffix = "_nir_drift_shift_split",
                        #start = 35,
                        #end = 40,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  'regularization/flat_reg_star.hdf5',
                        #reg_file_t = 'regularization/flat_reg_t.hdf5',
                        #output_suffix = 'order_cut_test',
                        #plot_continuum = True
                        #)
    
    
    #run_wobble(parameters)
