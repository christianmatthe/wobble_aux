#import sys
#sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
#import combine_results as cr
#from parameters import Parameters
#end local imports

import numpy as np
import matplotlib.pyplot as plt
import wobble
import tensorflow as tf
from tqdm import tqdm
import h5py
import os

start = 0
end = 61

R = end - start


#star_filename = '../wobble/regularization/dummy_star_K0_no_reg.hdf5'
star_filename = '../../wobble_aux/regularization/dummy_star_K0_no_reg.hdf5'
#if not os.path.isfile(star_filename):
with h5py.File(star_filename,'w') as f:
    f.create_dataset('L1_template', data=np.zeros(R)+1.e0)
    f.create_dataset('L2_template', data=np.zeros(R)+1.e0)
    #if K_star > 0:
        #f.create_dataset('L1_basis_vectors', data=np.zeros(R)+1.e5)
        #f.create_dataset('L2_basis_vectors', data=np.zeros(R)+1.e6)
        #f.create_dataset('L2_basis_weights', data=np.ones(R)) # never tuned, just need to pass to wobble                

#tellurics_filename = '../wobble/regularization/dummy_t_K3_no_reg.hdf5'
tellurics_filename = '../../wobble_aux/regularization/dummy_t_K3_no_reg.hdf5'
print(os.path.abspath(tellurics_filename))
#if not os.path.isfile(tellurics_filename):                
with h5py.File(tellurics_filename,'w') as f:
    print("test")
    if True:
        f.create_dataset('L1_template', data=np.zeros(R)+1.e0)
        f.create_dataset('L2_template', data=np.zeros(R)+1.e0)
    if True:
        f.create_dataset('L1_basis_vectors', data=np.zeros(R)+1.e0)
        f.create_dataset('L2_basis_vectors', data=np.zeros(R)+1.e0)
        f.create_dataset('L2_basis_weights', data=np.ones(R)) # never tuned, just need to pass to wobble
            
#star_filename = '../wobble/regularization/dummy_star_K0_no_reg.hdf5'
#if not os.path.isfile(star_filename):
    #with h5py.File(star_filename,'w') as f:
        #f.create_dataset('L1_template', data=np.zeros(R))
        #f.create_dataset('L2_template', data=np.zeros(R))
        ##if K_star > 0:
            ##f.create_dataset('L1_basis_vectors', data=np.zeros(R)+1.e5)
            ##f.create_dataset('L2_basis_vectors', data=np.zeros(R)+1.e6)
            ##f.create_dataset('L2_basis_weights', data=np.ones(R)) # never tuned, just need to pass to wobble                

#tellurics_filename = '../wobble/regularization/dummy_t_K3_10**5.hdf5'
#if not os.path.isfile(tellurics_filename):                
    #with h5py.File(tellurics_filename,'w') as f:
        #if True:
            #f.create_dataset('L1_template', data=np.zeros(R)+1.e5)
            #f.create_dataset('L2_template', data=np.zeros(R)+1.e5)
        #if True:
            #f.create_dataset('L1_basis_vectors', data=np.zeros(R)+1.e5)
            #f.create_dataset('L2_basis_vectors', data=np.zeros(R)+1.e5)
            #f.create_dataset('L2_basis_weights', data=np.ones(R)) # never tuned, just need to pass to wobble


###########
#star_filename = '../wobble/regularization/dummy_star_K0.hdf5'
#if not os.path.isfile(star_filename):
    #with h5py.File(star_filename,'w') as f:
        #f.create_dataset('L1_template', data=np.zeros(R)+1.e-2)
        #f.create_dataset('L2_template', data=np.zeros(R)+1.e2)
        ##if K_star > 0:
            ##f.create_dataset('L1_basis_vectors', data=np.zeros(R)+1.e5)
            ##f.create_dataset('L2_basis_vectors', data=np.zeros(R)+1.e6)
            ##f.create_dataset('L2_basis_weights', data=np.ones(R)) # never tuned, just need to pass to wobble                

#tellurics_filename = '../wobble/regularization/dummy_t_K3.hdf5'
#if not os.path.isfile(tellurics_filename):                
    #with h5py.File(tellurics_filename,'w') as f:
        #if True:
            #f.create_dataset('L1_template', data=np.zeros(R)+1.e4)
            #f.create_dataset('L2_template', data=np.zeros(R)+1.e6)
        #if True:
            #f.create_dataset('L1_basis_vectors', data=np.zeros(R)+1.e3)
            #f.create_dataset('L2_basis_vectors', data=np.zeros(R)+1.e8)
            #f.create_dataset('L2_basis_weights', data=np.ones(R)) # never tuned, just need to pass to wobble
######################
'''
#recreate original chunked readout of first 5 entries of default regto see if it reproduces results
def_star_filename = '../wobble/regularization/default_star.hdf5'
def_t_filename = '../wobble/regularization/default_t.hdf5'
chunk_size = 5
roll = 1

reg_file_name = '../wobble/regularization/def_chunk_{0}_roll{1}_star.hdf5'.format(chunk_size, roll)
with h5py.File(def_star_filename,'r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
            temp1 = f[key][()][0:chunk_size]
            temp1 = np.roll(temp1,-roll)
            temp2 = temp1
            while len(temp2) < R:
                temp2 = np.append(temp2, temp1)
            temp2 = temp2[0:R]
            if key in list(g.keys()):
                del g[key]
            g.create_dataset(key, data = temp2)
            
reg_file_name = '../wobble/regularization/def_chunk_{0}_roll{1}_t.hdf5'.format(chunk_size, roll)
with h5py.File(def_t_filename,'r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
            temp1 = f[key][()][0:chunk_size]
            temp1 = np.roll(temp1,roll)
            temp2 = temp1
            while len(temp2) < R:
                temp2 = np.append(temp2, temp1)
            temp2 = temp2[0:R]
            if key in list(g.keys()):
                del g[key]
            g.create_dataset(key, data = temp2)
'''
'''
# Guesstimation of good parameters from default star reg and wolf294 telluric regs
star_filename = '../wobble/regularization/reg_guess_star_0.hdf5'
if not os.path.isfile(star_filename):
    with h5py.File(star_filename,'w') as f:
        f.create_dataset('L1_template', data=np.zeros(R)+1.e-2)
        f.create_dataset('L2_template', data=np.zeros(R)+0)
        
tellurics_filename = '../wobble/regularization/reg_guess_t_0.hdf5'
if not os.path.isfile(tellurics_filename):                
    with h5py.File(tellurics_filename,'w') as f:
        if True:
            f.create_dataset('L1_template', data=np.zeros(R)+1.e4)
            f.create_dataset('L2_template', data=np.zeros(R)+1.e11)
        if True:
            f.create_dataset('L1_basis_vectors', data=np.zeros(R)+1.e3)
            f.create_dataset('L2_basis_vectors', data=np.zeros(R)+1.e5)
            f.create_dataset('L2_basis_weights', data=np.ones(R)) # never tuned, just need to pass to wobble
'''

'''
#make each factor high by 10^5 and see what happens
reg_file_name = '../wobble/regularization/dummy_star_K0_high.hdf5'
factor = 10**5
with h5py.File('../wobble/regularization/dummy_star_K0.hdf5','r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
            if key == 'L2_basis_weights':
                temp = f[key][()]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)
            else: 
                temp = f[key][()]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp * factor)

reg_file_name = '../wobble/regularization/dummy_t_K3_high.hdf5'
factor = 10**5
with h5py.File('../wobble/regularization/dummy_t_K3.hdf5','r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
            if key == 'L2_basis_weights':
                temp = f[key][()]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)
            else: 
                temp = f[key][()]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp * factor)
                
                
                
#make each factor *low* by 10^5 and see what happens
reg_file_name = '../wobble/regularization/dummy_star_K0_low.hdf5'
factor = 10**-5
with h5py.File('../wobble/regularization/dummy_star_K0.hdf5','r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
                if key == 'L2_basis_weights':
                    temp = f[key][()]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
                else: 
                    temp = f[key][()]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp * factor)

reg_file_name = '../wobble/regularization/dummy_t_K3_low.hdf5'
factor = 10**-5
with h5py.File('../wobble/regularization/dummy_t_K3.hdf5','r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
            if key == 'L2_basis_weights':
                temp = f[key][()]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)
            else: 
                temp = f[key][()]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp * factor)

###################
reg_file_name = '../wobble/regularization/dummy_star_K0_step.hdf5'
factor = 10**-5
with h5py.File('../wobble/regularization/dummy_star_K0.hdf5','r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
                if key == 'L2_basis_weights':
                    temp = f[key][()]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
                else: 
                    temp = f[key][()]
                    temp[3] = temp[3] * factor ** -2
                    for i in range(25,61):
                        temp[i] = temp[i] * factor ** -2
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp * factor)




reg_file_name = '../wobble/regularization/dummy_t_K3_step.hdf5'
factor = 10**-5
with h5py.File('../wobble/regularization/dummy_t_K3.hdf5','r') as f:
    with h5py.File(reg_file_name,'w') as g:
        for key in list(f.keys()):
            if key == 'L2_basis_weights':
                temp = f[key][()]
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)
            else: 
                temp = f[key][()]
                temp[3] = temp[3] * factor ** -2
                for i in range(25,61):
                    temp[i] = temp[i] * factor ** -2
                
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp * factor)
'''
