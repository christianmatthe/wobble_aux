import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/wobble_aux" )
from make_data_carmenes import *
from run_wobble import *

#make data
#prepare wobble data sets from carmenes raw data
data_directory="data/"
#as default place carmenes raw data from GTO archive into folder: data_directory/CARM_raw_data/<starname>/
if True: # GJ1148 :vis
    starname = "GJ1148"
    arm = "vis"
    make_data(starname, arm, data_directory)
        
#run_wobble
    #GJ1148
    parameters = Parameters(starname = "GJ1148",
                            data_suffix = "_vis_drift_shift",#Default data suffix produced by make_data .indicates that the visual arm data will be used Doppler shifted by SERVAL drifts
                            start = 11, # first order optimized
                            end = 53, # first order *not* optimized
                            chunk_size = 5, #number of orders that will be optimized in one go (circumvents RAM overflow) will use as many chunks as necessary to optimize all requested orders
                            niter = 160, #number of optimization iterations
                            reg_file_star =  'wobble_aux/regularization/dummy_star_K0_no_reg.hdf5', #regularization files used. These dummy files have all regularization set to 1 (i.e. 0 in wobbles log space)
                            reg_file_t = 'wobble_aux/regularization/dummy_t_K3_no_reg.hdf5',
                            output_suffix = "no_reg", #appended to results filename
                            results_dir = 'results/', #where to write results
                            data_dir= data_directory, # where to read data from
                            
                            plot_continuum = False # If True outputs plots for all orders and epochs showing the continuum normalization NOTE thats 3450 plots for this example
                            )
                            #additional options available see Parameters class definition in run_wobble.py
    run_wobble(parameters)
