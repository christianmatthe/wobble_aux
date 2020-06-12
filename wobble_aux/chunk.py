#Runs a subset of orders in separate python script to fix a wobble memory leak that progressively fills RAM when too man orders are optimized consecutively
import numpy as np
import shutil
from time import time
import h5py
import wobble
import os
import matplotlib.pyplot as plt
import dill

#local import
from run_wobble import *

#import parameters
with open(os.path.dirname(os.path.abspath(__file__)) + "/" + "carmenes_aux_files/chunk_parameters.pkl", "rb") as f:
    p = dill.load(f)

#Variables generated in run_wobble efore chunk: 
#temp_dir, plot_dir, start_time, epochs_list, chunks, data_file
#transfer Variables to local names
i  = p.i
temp_dir  = p.temp_dir
plot_dir  = p.plot_dir
start_time = p.start_time
epochs_list = p.epochs_list
orders_list = p.orders_list
chunks  = p.chunks
data_file = p.data_file


###Start of wobble Chunk
start_time_chunk = time()
reg_file_star_chunk, reg_file_t_chunk = reg_chunk(chunks[i], p.reg_file_star, p.reg_file_t)
start_order = chunks[i][0]
end_order = chunks[i][-1]

print("running wobble on star {0} with K_star = {1}, K_t = {2}, orders[{3},{4})".format(p.starname, p.K_star, p.K_t, start_order, end_order))
#orders = np.arange(start_order, end_order)
#Take intersection of chunk orders and order list to drop orders previously cut
#orders = list(set([x for x in range(start_order, end_order)]) & set(orders_list))
orders = chunks[i]
print("orders: ", orders)
data = wobble.Data(data_file, orders=orders, epochs=epochs_list, min_flux=10**-5, 
                   #min_snr=0,
                   min_snr = p.min_snr,
                   chunkQ = True,
                   parameters = p
                   )
'''
 #################                  )
#TODO Apply continnuum norm failed arder cutting: delete data.drop_orders fromm chunk and order list
#HACK this ignores the underlying issue: WHY ARE there different results in the first place?
#DIDNN'T WORK Note self.drop_orders:  [15]
#data.drop_orders [15]
#p.drop_orders [[15]] IIs of different format

try:
    for order in data.drop_orders:
        p.drop_orders.append(order) #append to drop orders list
    print("data.drop_orders", data.drop_orders)
    print("p.drop_orders", p.drop_orders)
except AttributeError:
    print("data.drop_orders is not defined")
    
orders_list = p.orders_list = [x for x in p.orders_list if x not in p.drop_orders]
#delete order from chunk list
for j, chunk in enumerate(chunks):
    chunk = [x for x in chunk if x not in p.drop_orders]
    chunks[j] = chunk
p.chunks = chunks
#update orders for  this chunk
orders = chunks[i]
#initialize chunk based variables
reg_file_star_chunk, reg_file_t_chunk = reg_chunk(chunks[i], p.reg_file_star, p.reg_file_t)
start_order = chunks[i][0]
end_order = chunks[i][-1]
#NOTE Would require passing parammeters back up
with open(os.path.dirname(os.path.abspath(__file__)) + "/carmenes_aux_files/chunk_parameters.pkl", "wb") as f:
    dill.dump(p, f)
#Sometimes the chunkwise continuum normilization finds different drop orders DUE TO DIFFERENT MIN SNR SETTING This could be  relatively easily fixed, bbut is it actually a bigger proble since low snr points in actually used continuum?
#Solution: set chunk flag in Data class and disallow order cutting so mmin_snr can be set to p.min_snr : didn't work
##########
'''

print("data.R: ", data.R)
results = wobble.Results(data=data)
print("results.R: ", results.R)

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
    
    if p.plots:
        epochs_to_plot = epochs_list[0:51:50]# alsways chooses 0th and 50th epoch from actually used epochs
        wobble.optimize_order(model, niter = p.niter, save_history = True,
                                basename = plot_dir + 'history', movies = False,
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
results_chunk = file_chunk_name(start_order, end_order, temp_dir)
results.write(results_chunk)
print("results saved as: {0}".format(results_chunk))
print("time elapsed: {0:.2f} minutes".format((time() - start_time)/60.0))
###end of wobble chunk
