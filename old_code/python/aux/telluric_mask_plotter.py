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

data_directory = "/data/cmatthe/wobble_data/data/"
results_directory = "/data/cmatthe/wobble_data/results/"
plot_dir = "/data/cmatthe/python/wobble_aux/telluric_masks/"

#results_name = "results_GJ436_Kstar0_Kt3_loop_4_reg"
#results_name = "results_GJ1148_Kstar0_Kt3_rs_orderwise_test_0_reg" #broken
results_name = "results_GJ1148_Kstar0_Kt3_pres:GJ436_l4_reg"
#results_name = "results_GJ1148_Kstar0_Kt3_rs_orderwise_test_0_reg_snr"
#results_name = "results_GJ1148_Kstar0_Kt3__l4reg_drift_shift"

#results_name = "results_GJ809A_Kstar0_Kt3__loop4_reg_snr"
#results_name = "results_Wolf294_Kstar0_Kt3__loop4_reg_snr"
#results_name = "results_GJ436_Kstar0_Kt3__nir_split_flat_reg"
#results_name = "results_GJ1148_Kstar0_Kt3__nir_split_flat_reg"
#results_name = "results_GJ3512_Kstar0_Kt3__nir_split_flat_reg"
starname = "GJ1148"
vis = True

show_mask = False

#orders = np.arange(11,53)
#lowest_optimized_order = 11
#orders = [49,50]
if vis == True:
    orders = np.arange(11,53)
    lowest_optimized_order = 11
    data = wobble.Data(starname+'_vis'+'_e2ds.hdf5', filepath= data_directory, orders=orders, min_flux=10**-5, min_snr=0)
else:
    orders = np.arange(0,56)
    lowest_optimized_order = 0
    data = wobble.Data(starname+'_nir_split'+'_e2ds.hdf5', filepath= data_directory, orders=orders, min_flux=10**-5, min_snr=0)

######
#plot_dir = plot_dir + results_name + "/"#+ "_bad_ep/"
plot_dir = plot_dir + results_name + "_bad_ord(49)/"
#plot_dir = plot_dir + "GJ1148_o49_snr" + "/"
os.makedirs(plot_dir, exist_ok = True)

results = wobble.Results(filename = results_directory + results_name +".hdf5")

telluric_mask = np.genfromtxt("/data/cmatthe/python/wobble_aux/telluric_masks/" + 
                              "telluric_mask_carm_short.dat"
)


#epochs = [ 5, 24, 65, 67, 74, 88, 93] #bad GJ436 
#epochs = [10, 13, 32, 48, 66] # bad GJ1148
#epochs = [8, 20, 23,  24,  70,  97,107, 123, 141, 154] #bad Wolf294 (20 as "good standart")
epochs_results = results.epochs
#epochs = epochs_results[10:50:5]
#epochs = [0,9,11,10,12,60]

#epochs = list(set(epochs) & set(results.epochs)) #build intersectionn  to make sure epoch 
epochs = epochs_results

#orders_to_plot = orders
orders_to_plot = [49]
for o in orders_to_plot:
    r = o - lowest_optimized_order
    for e_ind ,e in enumerate(epochs):
        #print(e)
        #print(np.array(epochs_results))
        ep = np.where(np.array(epochs_results) == e)[0][0]
        #e for data ep for results. assumes data is not epoch cut
        fig, (ax, ax2) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[4, 1]}, figsize=(12,5))
        xs = np.exp(data.xs[r][e])
        ax.scatter(xs, np.exp(data.ys[r][e]), marker=".", alpha=0.5, c='k', label='data', s=40)
        mask = data.ivars[r][e] <= 1.e-8
        ax.scatter(xs[mask], np.exp(data.ys[r][e][mask]), marker=".", alpha=1., c='white', s=20)
        ax.plot(xs, np.exp(results.star_ys_predicted[r][ep]), c='r', alpha=0.8)
        ax.plot(xs, np.exp(results.tellurics_ys_predicted[r][ep]), c='b', alpha=0.8)
        ylims = ax.get_ylim()
        xlims = ax.get_xlim()
        #also plot telluric mask
        if show_mask == True:
            ax.plot(telluric_mask[:,0], telluric_mask[:,1], c = 'g', alpha = 0.8)
        
        
        
        ax2.scatter(xs, np.exp(data.ys[r][e]) - np.exp(results.star_ys_predicted[r][ep]
                                                    + results.tellurics_ys_predicted[r][ep]), 
                    marker=".", alpha=0.5, c='k', label='data', s=40)
        ax2.scatter(xs[mask], np.exp(data.ys[r][e][mask]) - np.exp(results.star_ys_predicted[r][ep]
                                                    + results.tellurics_ys_predicted[r][ep])[mask], 
                    marker=".", alpha=1., c='white', s=20)
        ax.set_ylim([0.0,1.3])
        #ax.set_xlim([xs[0],xs[-1]])
        ax.set_xlim(xlims)
        ax2.set_ylim([-0.08,0.08])
        ax.set_xticklabels([])
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.05)
        plt.savefig(plot_dir+'results_synth_o{0}_e{1}.png'.format(o, e))
        plt.close(fig)
