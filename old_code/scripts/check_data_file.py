#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:18:11 2019

@author: cmatthe
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
#import wobble

# =============================================================================
# import sys
# =============================================================================


modulename = 'wobble'
if modulename not in dir():
    print( 'You have not imported the {} module'.format(modulename))
    import wobble

# =============================================================================
# data = wobble.Data(filename='HD189733_e2ds.hdf5',filepath="../data/")
# print(data.component_names)
# =============================================================================
#%%
#f = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/data/J11421+267_vis_e2ds.hdf5','r')
#f = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/data/GJ1148_vis_e2ds.hdf5','r')
#f = h5py.File('/data/cmatthe/wobble_data/data/GJ876_vis_e2ds.hdf5','r')
f = h5py.File('/data/cmatthe/wobble_data/data/GJ436_nir_split_e2ds.hdf5','r')
#f = h5py.File('/data/cmatthe/wobble_data/data/GJ436_nir_drift_shift_split_e2ds.hdf5','r')
#g = h5py.File('/data/cmatthe/wobble_data/data/GJ436_nir_split_e2ds.hdf5','w+')
#g = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/data/Wolf294_vis_e2ds.hdf5','r')
#f = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/wobble/regularization/default_star.hdf5','r')
#h = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/wobble/regularization/default_t.hdf5','r')
#f = h5py.File('/data/cmatthe/wobble_data/wobble/regularization/temp_star_chunk.hdf5','r')
#h = h5py.File('/data/cmatthe/wobble_data/wobble/regularization/temp_t_chunk.hdf5','r')

#f = h5py.File('/data/cmatthe/wobble_data/results/results_Wolf294_Kstar0_Kt3_def_chunk_5_n60.hdf5','r')
#g = h5py.File('/data/cmatthe/wobble_data/results/results_Wolf294_Kstar0_Kt3_testing.hdf5','r+')
#g = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/results/results_GJ1148_Kstar0_Kt3_stitched_def.hdf5','r')
#g = h5py.File('/home/cmatthe/lx39_data/cmatthe/compare_wobble_serval/wobbledir/results_GJ876_Kstar0_Kt3_n1000_stitched_def.hdf5','r')
#g = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/results/results_GJ876_Kstar0_Kt3_n1000_stitched_def.hdf5','r')
#g = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/results/results_GJ876_Kstar0_Kt3_stitched.hdf5','r')
#g = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/results/results_Wolf294_Kstar0_Kt3_stitched_def.hdf5','r')
#g = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/wobble/regularization/GJ1148_star_K0_orders[11,53)_stitched.hdf5','r')
#i = h5py.File('/home/cmatthe/lx39_data/cmatthe/wobble_data/wobble/regularization/GJ1148_t_K3_orders[11,53)_stitched.hdf5','r')

print(list(f.keys()))
for key in list(f.keys()):
    print(key + " shape" , f[key][()].shape)
print(f["dates"][()])
#print(len(f["dates"][()]))
#print((f["data"][()].shape))
#print((f["xs"][()].shape))

#splitting function that splits each order in half
#orders are the 0th array axis. we want to split the wavelenght axis (2nd axis) and make 2 0th axis entries
#data types that need splitting: data, ivars, xs
#def split_orders(array):
    #array_old = array
    #shape_old = array.shape
    #x_width_new = shape_old[2] // 2
    #shape_new = (2 * shape_old[0], shape_old[1], x_width_new)
    #array_new = np.zeros(shape_new)
    
    #for i in range(2 * shape_old[0]):
        #if i % 2 == 0:
            #array_new[i] = array_old[i//2,:,:x_width_new]
        #if i % 2 == 1:
            #array_new[i] = array_old[i//2,:,x_width_new:]
    #return array_new
    
#arr = f["xs"][()]
#arr_split = split_orders(arr)#

#print(arr_split[0,0],arr_split[1,0])
#print(arr_split[0,0].shape,arr_split[1,0].shape)

#print(f['data'][()])
#print(f['dates'][()])
#xs = f['xs'][()]
#print(f['xs'][()])
#print(f['xs'][()][1][0][0])
#x_deltas = [[i,xs[0,0,i+1] -xs[0,0,i]]  for i in range(len(xs[0,0])-2)]
#print(x_deltas)
#print(f['data'][()].shape,f['xs'][()].shape)

#print(list(h.keys()))
#print(f['epochs'][()])
#print(f['epochs'][()]==g['epochs'][()])
#temp = f['epochs'][()]
#g.create_dataset('epochs', data = temp)
#print(f['L2_template'][()])
#print(len(f['L1_template'][()]))
#print(g['L2_template'][()])
#print(len(g['L1_template'][()]))
#print(list(h.keys()))
#print(list(i.keys()))
#print(h['L2_template'][()])
#print(h['L2_basis_vectors'][()])
#print(h['L2_basis_weights'][()])
#print(len(h['L1_template'][()]))
#print(i['L2_template'][()])
#print(i['L2_basis_vectors'][()])
#print(i['L2_basis_weights'][()])
#print(len(i['L2_template'][()]))
#list(g['order0'].keys()))
# =============================================================================
#print(list(g.keys()))
#print(g["orders"][()][[0,3,4]])

# =============================================================================
# dates_utc = f['dates_utc'][()]
# 
# print(g['orders'][()])
# g['dates_utc'] = dates_utc
# =============================================================================
# =============================================================================
# print(f['dates_utc'][()])
# dates_utc = f['dates_utc'][()]
# 
# g['dates_utc'] = dates_utc
# print(g['dates_utc'][()])
# =============================================================================

#print(g['N'][()])
#print(type(f['orders'][()]))
#print(list(f['order0'].keys()))
#print(f['order0']['star_rvs'][()])
#f['data'][()]
#f['data'][(8)]
#f['L1_template'][()]
#np.allclose(g['star_time_rvs'][()],f['star_time_rvs'][()] )
#f['orders'][()]

# =============================================================================
# #%%
# #works only with unstitched -> there must be a difference between these files
# results_file = '/home/cmatthe/lx39_data/cmatthe/wobble_data/results/results_GJ876_Kstar0_Kt3.hdf5'
# #results_file = '/home/cmatthe/lx39_data/cmatthe/wobble_data/results/results_GJ876_Kstar0_Kt3_orders[11,21).hdf5'
# #results_file = '/home/cmatthe/lx39_data/cmatthe/wobble_data/results/results_GJ876_Kstar0_Kt3_stitched.hdf5'
# results = wobble.Results(filename = results_file)
# 
# results.combine_orders('star')
# 
