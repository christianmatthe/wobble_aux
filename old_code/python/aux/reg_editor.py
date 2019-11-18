import sys
sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
import combine_results as cr
from parameters import Parameters
#end local imports

import numpy as np
import matplotlib.pyplot as plt
import wobble
import tensorflow as tf
from tqdm import tqdm
import h5py
import os


#def edit_reg_file(base_reg_file, new_reg_file, orders):
    #"""
    #Creates new regularization files with the selected parameter shifted by the exponent
    #base_reg_file : `str`
        #file path of the base regularization file
    #new_reg_file : `str`
        #file path of the new regularization file
    #orders : list of `int`
        #orders for which the regularization file should be adjusted (assumes files runing from order 0 to 60)
    #"""
    
    #with h5py.File(base_reg_file,'r') as f:
        #with h5py.File(new_reg_file,'w') as g:
            #for key in list(f.keys()):
                #if key == parameter_key:
                    #temp = f[key][()]
                    #if type(orders) == int:
                        #temp[orders] = temp[orders] * factor
                    #else:
                        #for o in  orders:
                            #temp[o] = temp[o] * factor
                    #if key in list(g.keys()):
                        #del g[key]
                    #g.create_dataset(key, data = temp)
                #else: 
                    #temp = f[key][()]
                    #if key in list(g.keys()):
                        #del g[key]
                    #g.create_dataset(key, data = temp)
''' # For GJ876 recovery
#Manual recovery of reg file from plot
#Star
Star_L2_template = [3,3,3,2,3,2,3,3,3,3,0,-2,1,0,-1,0,1,-1,2,-1,4,4,4,4,4,6,7,4,5,6,4,4,5,1,4,2,2,0,2,5,5,5]
Star_L1_template = [-4,-2,-2,1,-4,-2,-4,-1,-3,-4,-2,-2,2,-3,-1,-2,-2,-2,-3,-4,-3,-3,-2,-2,-3,-2,0,-3,-2,0,-3,-1,-3,-2,-1,-3,-1,-4,-4,0,-2,-3]

T_L1_basis_vectors = [4,2,3,2,0,2,3,0,-1,2,1,2,3,3,3,2,3,1,4,0,3,3,5,4,4,4,3,0,0,2,3,3,4,4,4,-1,3,-3,4,4,4,2]
T_L1_template = [4,4,4,4,2,3,4,5,4,4,4,3,5,2,4,4,2,4,17,14,4,4,17,5,4,5,27,36,5,3,3,3,5,4,5,4,4,3,4,3,16,19]
T_L2_basis_vectors = [12,8,7,6,8,7,10,7,6,5,7,7,7,6,6,7,6,7,6,5,7,7,7,7,6,7,7,8,7,6,7,7,7,6,6,7,6,7,6,6,6,2]
T_L2_basis_weights = [0 for i in range(11,53)]
T_L2_template = [12,11,12,11,13,10,11,11,11,11,11,11,11,10,12,11,10,12,11,8,9,10,11,22,10,11,12,20,9,10,9,11,11,13,11,11,12,11,11,11,10,11]

star_dictionary = {"L1_template" : Star_L1_template, "L2_template" : Star_L2_template}
t_dictionary = {"L1_template" : T_L1_template, "L2_template" : T_L2_template, "L1_basis_vectors" : T_L1_basis_vectors, "L2_basis_vectors" : T_L2_basis_vectors, "L2_basis_weights": T_L2_basis_weights}

if (len(Star_L1_template) == len(Star_L2_template) == len(T_L1_basis_vectors) == len(T_L1_template) == len(T_L2_basis_vectors) == len (T_L2_basis_weights) == len(T_L2_template)):
    print("All lenghts are equal")
else:
    print("Lengths unequal")
    quit()

                #Tellurics
#extend the [11,53) reg file for GJ876 to full length                    
orders = [i for i in range(11,53)]
#base reg_file is a dummy file with the correct format
base_reg_file = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5"
new_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ876_t_K3_orders[11,53)_stitched_reformatted.hdf5"
#old_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ876_t_K3_orders[11,53)_stitched.hdf5"
dic = t_dictionary

with h5py.File(base_reg_file,'r') as f:
    with h5py.File(new_reg_file,'w') as g:
        #with h5py.File(old_reg_file,'r') as h:
        for key in list(f.keys()):
            temp = f[key][()]
            for o in orders:
                temp[o] = 10** dic[key][o - orders[0]]
            if key in list(g.keys()):
                del g[key]
            g.create_dataset(key, data = temp)

        #Star
#extend the [11,53) reg file for GJ876 to full length                    
orders = [i for i in range(11,53)]
#base reg_file is a dummy file with the correct format
base_reg_file = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5"
new_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ876_star_K0_orders[11,53)_stitched_reformatted.hdf5"
#old_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ876_star_K0_orders[11,53)_stitched.hdf5"
dic = star_dictionary

with h5py.File(base_reg_file,'r') as f:
    with h5py.File(new_reg_file,'w') as g:
        #with h5py.File(old_reg_file,'r') as h:
        for key in list(f.keys()):
            temp = f[key][()]
            for o in orders:
                temp[o] = 10** dic[key][o - orders[0]]
            if key in list(g.keys()):
                del g[key]
            g.create_dataset(key, data = temp)                           
                        
'''
#GJ1148 reformat to 60 length

                #Tellurics
#extend the [11,53) reg file for GJ876 to full length                    
orders = [i for i in range(11,53)]
#base reg_file is a dummy file with the correct format
base_reg_file = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5"
#new_reg_file is the destination
new_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ1148_t_K3_orders[11,53)_stitched_reformatted.hdf5"
#old_reg_file contains the information to be inserted into the dummy file at the correct places
old_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ1148_t_K3_orders[11,53)_stitched.hdf5"


with h5py.File(base_reg_file,'r') as f:
    with h5py.File(new_reg_file,'w') as g:
        with h5py.File(old_reg_file,'r') as h:
            for key in list(f.keys()):
                temp = f[key][()]
                for o in orders:
                    r=o - orders[0]
                    temp[o] = h[key][()][r]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)

        #Star
#extend the [11,53) reg file for GJ876 to full length                    
orders = [i for i in range(11,53)]
#base reg_file is a dummy file with the correct format
base_reg_file = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5"
#new_reg_file is the destination
new_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ1148_star_K0_orders[11,53)_stitched_reformatted.hdf5"
#old_reg_file contains the information to be inserted into the dummy file at the correct places
old_reg_file = "/data/cmatthe/wobble_data/wobble/regularization/GJ1148_star_K0_orders[11,53)_stitched.hdf5"


with h5py.File(base_reg_file,'r') as f:
    with h5py.File(new_reg_file,'w') as g:
        with h5py.File(old_reg_file,'r') as h:
            for key in list(f.keys()):
                temp = f[key][()]
                for o in orders:
                    r=o - orders[0]
                    temp[o] = h[key][()][r]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)              
