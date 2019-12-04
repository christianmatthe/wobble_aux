# Meant to recreate functionality of optimize_top_3.1.2._data_suffix but without calling subscipts via os.system
import numpy as np
import shutil
from time import time
import h5py
import wobble
import os
import matplotlib.pyplot as plt 

#Reproduce combine results first: functions required for runnign wobble in small chunks of orders so as not to overflow RAM
def file_chunk_name(start_order, end_order, chunk_dir):
    results_file_chunk_name = chunk_dir + 'orders[{0},{1})'.format(start_order, end_order) + '.hdf5' #this name must be same as in chunkscript
    return results_file_chunk_name

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

def results_file_stitch(start, end, chunk_size, results_file, chunk_dir):
    """turns chunks into one continuous file, as if wobble had been run in one piece"""
    print("Writing results to {0}".format(results_file))
    #initialize results_file as copy of first chunk_list
    # NOTE files must have the same epochs
    chunks = chunk_list(start, end, chunk_size)
    start_order = chunks[0, 0]
    end_order = chunks[0, 1]
    filename_chunk = file_chunk_name(start_order, end_order, chunk_dir)
    shutil.copy(filename_chunk, results_file)
    
    #replace number of orders since chunk will not have all of them
    with h5py.File(results_file, 'r+') as f:
        n_orders = end-start  # note this is not the number of optimized orders if any orders  where dropped
        f['R'][()] = n_orders
        
        #append orders from other chunks
        if len(chunks) > 1:
            for i in range(1, len(chunks)):
                #load files and combine them
                start_order = chunks[i, 0]
                end_order = chunks[i, 1]
                filename_chunk = file_chunk_name(start_order, end_order, chunk_dir)
                with h5py.File(filename_chunk,'r') as g:
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
            
#TODO Parameter functions 

class Parameters: 
    """
    The Parameter dictonary object: contains the parameters specifying how to run wobble.
    
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
                 min_snr = 60
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
    start_chunk = int(chunk[0])
    end_chunk = int(chunk[1])
    #NOTE assumes reg file starts at order 0 TODO implenent check that file is 61 orders long i.e. that reg file is valid
    reg_file_star_chunk = 'regularization/temp_star_chunk.hdf5'
    with h5py.File(reg_file_star,'r') as f:
        with h5py.File(reg_file_star_chunk,'w') as g:
            for key in list(f.keys()):
                    temp = f[key][()][start_chunk : end_chunk]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
    reg_file_t_chunk = 'regularization/temp_t_chunk.hdf5'
    with h5py.File(reg_file_t,'r') as f:
        with h5py.File(reg_file_t_chunk,'w') as g:
            for key in list(f.keys()):
                    temp = f[key][()][start_chunk : end_chunk]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
    return reg_file_star_chunk, reg_file_t_chunk

def run_wobble(parameters):
    #TODO make these parameters or eliminate
    plots = True
    movies = False
    
    p = parameters
    results_name = 'results_{0}_Kstar{1}_Kt{2}_'.format(p.starname, p.K_star, p.K_t, p.niter) + p.output_suffix
    results_file_base = p.results_dir + results_name
    results_file = results_file_base + '.hdf5'
    data_file = p.data_dir + p.starname  + p.data_suffix + '_e2ds.hdf5'
    
    temp_dir = p.results_dir + '/temp_' + results_name + '/'
    plot_dir = p.results_dir + '/plots_' + results_name + '/'
    #make (output) directory
    #TODO these data and results dirs should be handlesd somewhere else
    os.makedirs(p.results_dir, exist_ok = True)
    os.makedirs(p.data_dir, exist_ok = True)
    
    os.makedirs(temp_dir, exist_ok = True)
    os.makedirs(plot_dir, exist_ok = True)
    
    start_time = time()
    chunks = chunk_list(p.start, p.end, p.chunk_size)
    #generate epoch list
    data = wobble.Data(data_file, orders = np.arange(p.start, p.end), min_flux=10**-5, min_snr = p.min_snr)
    epochs_list = data.epochs.tolist()
    
    #Loop over chunks
    for i in range(len(chunks)):
        start_time_chunk = time()
        reg_file_star_chunk, reg_file_t_chunk = reg_chunk(chunks[i], p.reg_file_star, p.reg_file_t)
        start_order = chunks[i, 0]
        end_order = chunks[i, 1]
        
        print("running wobble on star {0} with K_star = {1}, K_t = {2}, orders[{3},{4})".format(p.starname, p.K_star, p.K_t, start_order, end_order))
        orders = np.arange(start_order, end_order)
        data = wobble.Data(data_file, orders=orders, epochs=epochs_list, min_flux=10**-5, min_snr=0)
        results = wobble.Results(data=data)
        
        print("data loaded")
        print("time elapsed: {0:.2f} min".format((time() - start_time)/60.0))
        elapsed_time = time() - start_time
        
        star_learning_rate = 0.1
        telluric_learning_rate = 0.01
        for r,o in enumerate(orders):
            model = wobble.Model(data, results, r)
            model.add_star('star', variable_bases = p.K_star, 
                            regularization_par_file = reg_file_star_chunk, 
                            learning_rate_template = star_learning_rate)
            model.add_telluric('tellurics', rvs_fixed = True, variable_bases = p.K_t, 
                                regularization_par_file = reg_file_t_chunk , 
                                learning_rate_template = telluric_learning_rate)
            print("--- ORDER {0} ---".format(o))
            
            if plots:
                epochs_to_plot = epochs_list[0:51:50]# alsways choses 0th and 50th epoch from actually used epochs
                wobble.optimize_order(model, niter = p.niter, save_history = True,
                                      basename = plot_dir + 'history', movies = movies,
                                      epochs_to_plot = epochs_to_plot)
                fig, ax = plt.subplots(1, 1, figsize=(8,5))
                ax.plot(data.dates, results.star_rvs[r] + data.bervs - np.mean(results.star_rvs[r] + data.bervs), 
                        'k.', alpha=0.8)
                ax.plot(data.dates, data.pipeline_rvs + data.bervs - np.mean(data.pipeline_rvs + data.bervs), 
                        'r.', alpha=0.5)   
                ax.set_ylabel('RV (m/s)', fontsize=14)     
                ax.set_xlabel('BJD', fontsize=14)   
                plt.savefig(plot_dir+'results_rvs_o{0}.png'.format(o))
                plt.close(fig)
                for ep in epochs_to_plot:
                    # TODO Check if line below is necessary
                    e = np.where(np.array(epochs_list) == ep)[0][0] # updates e to be the location in the epochlist without the cut epochs
                    
                    fig, (ax, ax2) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[4, 1]}, figsize=(12,5))
                    xs = np.exp(data.xs[r][e])
                    ax.scatter(xs, np.exp(data.ys[r][e]), marker=".", alpha=0.5, c='k', label='data', s=40)
                    mask = data.ivars[r][e] <= 1.e-8
                    ax.scatter(xs[mask], np.exp(data.ys[r][e][mask]), marker=".", alpha=1., c='white', s=20)
                    ax.plot(xs, np.exp(results.star_ys_predicted[r][e]), c='r', alpha=0.8)
                    ax.plot(xs, np.exp(results.tellurics_ys_predicted[r][e]), c='b', alpha=0.8)
                    ax2.scatter(xs, np.exp(data.ys[r][e]) - np.exp(results.star_ys_predicted[r][e]
                                                                + results.tellurics_ys_predicted[r][e]), 
                                marker=".", alpha=0.5, c='k', label='data', s=40)
                    ax2.scatter(xs[mask], np.exp(data.ys[r][e][mask]) - np.exp(results.star_ys_predicted[r][e]
                                                                + results.tellurics_ys_predicted[r][e])[mask], 
                                marker=".", alpha=1., c='white', s=20)
                    ax.set_ylim([0.0,1.3])
                    ax2.set_ylim([-0.08,0.08])
                    ax.set_xticklabels([])
                    fig.tight_layout()
                    fig.subplots_adjust(hspace=0.05)
                    plt.savefig(plot_dir+'results_synth_o{0}_e{1}.png'.format(o, epochs_list[e]))
                    plt.close(fig)
            else:
                wobble.optimize_order(model, niter= p.niter)
            del model # not sure if this does anything
            print("order {1} optimization finished. time elapsed: {0:.2f} min".format((time() - start_time)/60.0, o))
            print("this order took {0:.2f} min".format((time() - start_time - elapsed_time)/60.0))
            elapsed_time = time() - start_time
        print("all orders in chunk optimized.")
        chunk_dir = temp_dir
        results_chunk = file_chunk_name(start_order, end_order, chunk_dir)
        results.write(results_chunk)
        print("results saved as: {0}".format(results_chunk))
        print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
    print("all chunks optimized: writing combined file") 
    results_file_stitch(p.start, p.end, p.chunk_size, results_file, chunk_dir)
    
    #Combine orders
    results = wobble.Results(filename = results_file)
    results.combine_orders('star')
    print("final RVs calculated.")
    print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
    results.write(results_file)
    append_dates_utc(results_file, data_file)# cannot be done before results.write, as .write will remove dates_utc
    print("results saved as: {0}".format(results_file))
    print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
    
    
    
    
    
    


if __name__ == "__main__":
    

    parameters = Parameters(starname = "GJ436",
                            data_suffix = "_vis_drift_shift",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/GJ436_orderwise_avcn_l4_star.hdf5',
                            reg_file_t = 'regularization/GJ436_orderwise_avcn_l4_t.hdf5',
                            output_suffix = "git_run_wobble_test0")
    
    run_wobble(parameters)
