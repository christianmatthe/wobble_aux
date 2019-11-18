import sys
sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
import combine_results as cr
from parameters import Parameters

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import wobble
from time import time
import h5py
import os
import sys
import ast
import yaml
#NOTE WORKS only in Wob_env_2 (wobble_19_03_2019, deprecated in wobble_14_06_2019) 
#TODO pass results_file_base to chunk
#load parameters
parameters = Parameters(filename = "yaml_temp/optimize_parameters_chunk.yaml")
parameter_dict = parameters.dictionary

start_order= parameter_dict["start_chunk"]#TODO fix this name mismatch
end_order= parameter_dict["end_chunk"]
starname = parameter_dict["starname"]
K_star = parameter_dict["K_star"]
K_t = parameter_dict["K_t"]
epochs_list = parameter_dict["epochs_list"]
niter = parameter_dict["niter"]
defaultQ = parameter_dict["defaultQ"]

#reg_log_file = parameter_dict["reg_log_file"]


top_directory = parameter_dict["top_directory"]
run_directory = parameter_dict["run_directory"]
data_directory = parameter_dict["data_directory"]
results_directory = parameter_dict["results_directory"]
    

if True:
    plots = True
    epochs = [0,50] # to plot
    movies = False
    
    #NOTE: These reg file names will use the uncombined chunkwise regularisation files. the chunksize and range must therefore be matching
    if "reg_file_star_chunk" in parameter_dict and "reg_file_t_chunk" in parameter_dict:
        star_reg_file = parameter_dict["reg_file_star_chunk"]
        tellurics_reg_file = parameter_dict["reg_file_t_chunk"]
        
    #else:
        #star_reg_file = '../wobble/regularization/{0}_star_K{1}_orders[{2},{3}).hdf5'.format(starname, K_star, start_order, end_order)
        #tellurics_reg_file = '../wobble/regularization/{0}_t_K{1}_orders[{2},{3}).hdf5'.format(starname, K_t, start_order, end_order)
        
    plot_dir = results_directory + '/plots_{0}_Kstar{1}_Kt{2}_'.format(starname, K_star, K_t, niter)+ parameter_dict["output_suffix"] +'/'
    
    results_file = results_directory + 'results_{0}_Kstar{1}_Kt{2}_{5}_orders[{3},{4}).hdf5'.format(starname, K_star, K_t, start_order, end_order, parameters.dictionary["output_suffix"])
    
    

    print("running wobble on star {0} with K_star = {1}, K_t = {2}, orders[{3},{4})".format(starname, K_star, K_t, start_order, end_order))
    start_time = time()
    orders = np.arange(start_order, end_order)
    data = wobble.Data(starname+'_vis'+'_e2ds.hdf5', filepath= data_directory, orders=orders, epochs=epochs_list, min_flux=10**-5, min_snr=0)
    # min_snr=0 is part of a HACK that centralises epoch cutting to the top level script, will this interfere with order cutting?
    
    results = wobble.Results(data=data)
    
    ''' doesn't work becaus the enoumeration of orders will restart at 0
    #Load wobble results file should it already exist, so that only the current optimzed chunk gets replaced
    #note this may fail if the data is not the same as the data used to create the original results file
    try:
        results = wobble.Results(filename=results_file)
    except FileNotFoundError:
        results = wobble.Results(data=data)
    '''    
    
    print("data loaded")
    print("time elapsed: {0:.2f} min".format((time() - start_time)/60.0))
    elapsed_time = time() - start_time
    

    if plots:
        print("plots will be saved under directory: {0}".format(plot_dir))
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    star_learning_rate = 0.1
    telluric_learning_rate = 0.01
    for r,o in enumerate(orders):
        model = wobble.Model(data, results, r)
        model.add_star('star', variable_bases=K_star, 
                        regularization_par_file=star_reg_file, 
                        learning_rate_template=star_learning_rate)
        model.add_telluric('tellurics', rvs_fixed=True, variable_bases=K_t, 
                            regularization_par_file=tellurics_reg_file, 
                            learning_rate_template=telluric_learning_rate)
        print("--- ORDER {0} ---".format(o))
        '''
        with open(reg_log_file,"a") as f:
            lst=[]
            for name in model.components[1].regularization_par:
                lst.append(int(np.log10(getattr(model.components[1], name))))
            #print("{0},{1},{2},{3},{4}".format(lst[0],lst[1],lst[2],lst[3],lst[4]))            
            #f.write("{0},{1},{2},{3},{4}".format(lst[0],lst[1],lst[2],lst[3],lst[4]))
            #print("{0},{1},{2},{3},{4}".format(lst[0],lst[1],lst[2],lst[3],lst[4]))
            f.write("{5},{0},{1},{2},{3},{4}\n".format(lst[0],lst[1],lst[2],lst[3],lst[4],o))
        '''
        if plots:
            wobble.optimize_order(model, niter=niter, save_history=True, 
                                  basename=plot_dir+'history', movies=movies,
                                  epochs_to_plot=epochs, rv_uncertainties=True) 
            fig, ax = plt.subplots(1, 1, figsize=(8,5))
            ax.plot(data.dates, results.star_rvs[r] + data.bervs - np.mean(results.star_rvs[r] + data.bervs), 
                    'k.', alpha=0.8)
            ax.plot(data.dates, data.pipeline_rvs + data.bervs - np.mean(data.pipeline_rvs + data.bervs), 
                    'r.', alpha=0.5)   
            ax.set_ylabel('RV (m/s)', fontsize=14)     
            ax.set_xlabel('BJD', fontsize=14)   
            plt.savefig(plot_dir+'results_rvs_o{0}.png'.format(o))
            plt.close(fig)           
            for e in epochs:
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
                plt.savefig(plot_dir+'results_synth_o{0}_e{1}.png'.format(o, e))
                plt.close(fig)
        else:
            wobble.optimize_order(model, niter=niter)
        del model # not sure if this does anything
        print("order {1} optimization finished. time elapsed: {0:.2f} min".format((time() - start_time)/60.0, o))
        print("this order took {0:.2f} min".format((time() - start_time - elapsed_time)/60.0))
        elapsed_time = time() - start_time
    
    print("all orders optimized.")
    print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
    
    #makes issues in chunked version
    '''
    results.combine_orders('star')
    
    print("final RVs calculated.")
    print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
    '''
        
    results.write(results_file)
        
    print("results saved as: {0}".format(results_file))
    print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
