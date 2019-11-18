import h5py
import numpy as np
import matplotlib.pyplot as plt
#import wobble

# =============================================================================
# import sys
# =============================================================================
def split_orders(array):
    array_old = array
    shape_old = array.shape
    x_width_new = shape_old[2] // 2
    shape_new = (2 * shape_old[0], shape_old[1], x_width_new)
    array_new = np.zeros(shape_new)
    
    for i in range(2 * shape_old[0]):
        if i % 2 == 0:
            array_new[i] = array_old[i//2,:,:x_width_new]
        if i % 2 == 1:
            array_new[i] = array_old[i//2,:,x_width_new:]
    return array_new

def split_orders_file(filename):
    split_sets = ["data", "ivars", "xs"]
    with h5py.File(filename,'r') as f:
        with h5py.File(filename.split("e2ds")[0] + "split_e2ds.hdf5",'w') as g:
            for key in list(f.keys()):
                temp = f[key][()]
                if key in split_sets:
                    temp_split = split_orders(temp)
                    temp = temp_split
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)

if __name__ == "__main__":


    #starname = "GJ1148"
    #starname = "GJ3512"
    starname = "GJ436"
    
    #data_suffix = "_nir"
    data_suffix = "_nir_drift_shift"

    modulename = 'wobble'
    if modulename not in dir():
        print( 'You have not imported the {} module'.format(modulename))
        import wobble
        
        
    #f = h5py.File('/data/cmatthe/wobble_data/data/GJ436_nir_e2ds.hdf5','r')
    #g = h5py.File('/data/cmatthe/wobble_data/data/GJ436_nir_split_e2ds.hdf5','w+')
    
    split_sets = ["data", "ivars", "xs"]
    
    with h5py.File('/data/cmatthe/wobble_data/data/{0}{1}_e2ds.hdf5'.format(starname, data_suffix),'r') as f:
        with h5py.File('/data/cmatthe/wobble_data/data/{0}{1}_split_e2ds.hdf5'.format(starname,data_suffix),'w') as g:
            for key in list(f.keys()):
                temp = f[key][()]
                if key in split_sets:
                    temp_split = split_orders(temp)
                    temp = temp_split
                if key in list(g.keys()):
                    del g[key]
                g.create_dataset(key, data = temp)
    
