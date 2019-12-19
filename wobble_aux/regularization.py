#Optimize a set of regularization parameters for a given star using RV deviation from a "known" planet solution as the criterion
#basis is reg_gradient_search.py
import h5py


def update_reg_file(base_reg_file, new_reg_file, orders,  parameter_key, exponent):
    """
    Creates new regularization files with the selected parameter shifted by the exponent
    base_reg_file : `str`
        file path of the base regularization file
    new_reg_file : `str`
        file path of the new regularization file
    orders : list of `int`
        orders for which the regularization file should be adjusted (assumes files runing from order 0 to 60)
    parameter_key: `str`
        key of the parameter to be adjusted
    exponent: `int`
        exponent shift of the regularization parameter
    """
    factor = 10**exponent
    with h5py.File(base_reg_file,'r') as f:
        if not parameter_key in list(f.keys()):
            raise Exception("Invalid parameter_key in update_reg_file")
        with h5py.File(new_reg_file,'w') as g:
            for key in list(f.keys()):
                if key == parameter_key:
                    temp = f[key][()]
                    for o in  orders:
                        temp[o] = temp[o] * factor
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
                else: 
                    temp = f[key][()]
                    if key in list(g.keys()):
                        del g[key]
                    g.create_dataset(key, data = temp)
                    
