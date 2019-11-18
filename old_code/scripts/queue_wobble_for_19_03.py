import sys
sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
import combine_results as cr
from parameters import Parameters
#end local imports

import os
from time import time
from time import sleep
#TODO create parameter_dict class
# Have list of parameter dicts that are passed one by one to the scipts
# include script as parameter so different scripts can be executed

####
#min_sleep = 2*60
#print("sleeping for {0} minutes".format(min_sleep))
#sleep(60*min_sleep)
####


#scripts = ["optimize_top_1.py"]
#length = len(scripts)
#list_of_parameter_dictionaries
# #TODO add script 
#and regularization selection parameters to parameters class, optimize script needs to be adjusted to react

all_parameters = []


################################
# BEGIN For Presentation: All originally performed on order 6 continuum Normalization

######GJ876
parameters = Parameters(starname = "GJ876", defaultQ = False,
                        data_suffix = "_vis_drift_shift",
                        start = 11,
                        end = 53,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        output_suffix = "pres:GJ436_l4_reg+snr+drift")
                        
parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
all_parameters.append(parameters.dictionary)

######GJ1148

    
##reproduce initial "default" results
#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_5_roll1_star.hdf5",
                        #reg_file_t = "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_5_roll1_t.hdf5",
                        #output_suffix = 'pres:def_roll1')
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.3_data_suffix_no_snr.py"})
#all_parameters.append(parameters.dictionary)

##reproduce results with GJ436 optimized regularization
#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.3_data_suffix_no_snr.py"})
#all_parameters.append(parameters.dictionary)

##reproduce results with GJ436 optimized regularization + snr
#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg+snr")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

##reproduce results with GJ436 optimized regularization + drift
#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg+drift")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.3_data_suffix_no_snr.py"})
#all_parameters.append(parameters.dictionary)

#reproduce results with GJ436 optimized regularization + drift + snr
#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg+snr+drift")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)



############GJ436

##reproduce initial "default" results
#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_vis",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_5_roll1_star.hdf5",
                        #reg_file_t = "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_5_roll1_t.hdf5",
                        #output_suffix = 'pres:def_roll1')
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.3_data_suffix_no_snr.py"})
#all_parameters.append(parameters.dictionary)

##reproduce results with GJ436 optimized regularization
#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_vis",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.3_data_suffix_no_snr.py"})
#all_parameters.append(parameters.dictionary)

##reproduce results with GJ436 optimized regularization + snr
#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_vis",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg+snr")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

##reproduce results with GJ436 optimized regularization + drift
#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg+drift")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.3_data_suffix_no_snr.py"})
#all_parameters.append(parameters.dictionary)

##reproduce results with GJ436 optimized regularization + drift + snr
#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg+snr+drift")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#########
##Check GJ876
#parameters = Parameters(starname = "GJ876", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = "pres:GJ436_l4_reg+snr+drift_HighRV")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

##Dummy reg file tests

#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_data/wobble/regularization/dummy_star_K0_no_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_data/wobble/regularization/dummy_t_K3_no_reg.hdf5',
                        #output_suffix = "pres:dummy_reg_0")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_data/wobble/regularization/dummy_star_K0_no_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_data/wobble/regularization/dummy_t_K3_10**5.hdf5',
                        #output_suffix = "pres:dummy_reg_t_10**5")
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

# END
################################

#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_nir_drift_shift_split",
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_nir_split_flat_reg_drift_shift_post2017')
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#for i in range(5,-1,-1):
    #parameters = Parameters(starname = "GJ436", defaultQ = False, 
                            #data_suffix = "_vis_drift_shift",
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  "/data/cmatthe/wobble_reg_search/GJ436_orderwise_snr+drift_shift_1/loop_{0}/next_base_star_reg.hdf5".format(i),
                            #reg_file_t = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_snr+drift_shift_1/loop_{0}/next_base_t_reg.hdf5".format(i),
                            #output_suffix = 'snr+ds_l{0}_reg'.format(i))
    #parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
    #all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_l4reg_drift_shift_continuum1')
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_nir_drift_shift_split",
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_nir_split_flat_reg_drift_shift_continuum1')
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_l4reg_drift_shift')
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ1148", defaultQ = False,
                        #data_suffix = "_nir_drift_shift_split",
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_nir_split_flat_reg_drift_shift')
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_vis_drift_shift",
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_l4reg_drift_shift')
                        
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ436", defaultQ = False,
                        #data_suffix = "_nir_drift_shift_split",
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_nir_split_flat_reg_drift_shift_redata')
#parameters.dictionary.update({"script" : "optimize_top_3.1.2_data_suffix.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ3512", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_loop4_reg_snr')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ3512", defaultQ = False, 
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_nir_split_flat_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr_NIR.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ1148", defaultQ = False, 
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_nir_split_flat_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr_NIR.py"})
#all_parameters.append(parameters.dictionary)


#parameters = Parameters(starname = "GJ436", defaultQ = False, 
                        #start = 0,
                        #end = 56,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_nir_split_flat_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr_NIR.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ436", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_loop4_reg_snr')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_loop4_reg_snr')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ809A", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_loop4_reg_snr')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ809A", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = '_loop4_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)



#parameters = Parameters(starname = "GJ1148", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = 'rs_orderwise_test_0_reg_snr')
#parameters.dictionary.update({"script" : "optimize_top_3.1.1_snr.py"})
#all_parameters.append(parameters.dictionary)


#parameters = Parameters(starname = "GJ436", defaultQ = False, 
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_7/next_base_star_reg.hdf5",
                            #reg_file_t = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_7/next_base_t_reg.hdf5",
                            #output_suffix = 'servalmask_test')
#parameters.dictionary.update({"script" : "optimize_top_3.3_servalmask.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ436", defaultQ = False, 
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 1000,
                            #reg_file_star =  "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_7/next_base_star_reg.hdf5",
                            #reg_file_t = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_7/next_base_t_reg.hdf5",
                            #output_suffix = 'loop_7_reg_n1000')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ876", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  "/data/cmatthe/wobble_data/wobble/regularization/dummy_star_K0.hdf5",
                        #reg_file_t = "/data/cmatthe/wobble_data/wobble/regularization/dummy_t_K3.hdf5",
                        #output_suffix = '_reg_dummy')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ876", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_5_roll1_star.hdf5",
                        #reg_file_t = "/data/cmatthe/wobble_data/wobble/regularization/def_chunk_5_roll1_t.hdf5",
                        #output_suffix = '_reg_def_chunk_roll1')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#for i in range(0,5,1):
    #parameters = Parameters(starname = "GJ436", defaultQ = False, 
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  "/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_{0}/next_base_star_reg.hdf5".format(i),
                            #reg_file_t = "/data/cmatthe/wobble_reg_search/GJ436_test_run_0.1/loop_{0}/next_base_t_reg.hdf5".format(i),
                            #output_suffix = 'loop_{0}_reg'.format(i-5))
    #parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
    #all_parameters.append(parameters.dictionary)

#for i in range(3,15,2):
    #parameters = Parameters(starname = "GJ436", defaultQ = False, 
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_{0}/next_base_star_reg.hdf5".format(i),
                            #reg_file_t = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_{0}/next_base_t_reg.hdf5".format(i),
                            #output_suffix = 'loop_{0}_reg'.format(i))
    #parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
    #all_parameters.append(parameters.dictionary)

#for i in range(14,0,-2):
    #parameters = Parameters(starname = "GJ436", defaultQ = False, 
                            #start = 11,
                            #end = 53,
                            #chunk_size = 5,
                            #niter = 160,
                            #reg_file_star =  "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_{0}/next_base_star_reg.hdf5".format(i),
                            #reg_file_t = "/data/cmatthe/wobble_reg_search/GJ436_orderwise_avcn_0/loop_{0}/next_base_t_reg.hdf5".format(i),
                            #output_suffix = 'loop_{0}_reg'.format(i))
    #parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
    #all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ876", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  "/data/cmatthe/wobble_data/wobble/regularization/GJ876_star_K0_orders[11,53)_stitched_reformatted.hdf5",
                        #reg_file_t = "/data/cmatthe/wobble_data/wobble/regularization/GJ876_t_K3_orders[11,53)_stitched_reformatted.hdf5",
                        #output_suffix = 'wobble_reg_reformatted')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ285", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = 'orderwise_test_0_reg_YZ_CMi')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = 'rs_orderwise_test_0_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ1148", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = 'rs_orderwise_test_0_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ876", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = 'rs_orderwise_test_0_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Teegarden", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_star_reg.hdf5',
                        #reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_0/loop_4/next_base_t_reg.hdf5',
                        #output_suffix = 'rs_orderwise_test_0_reg')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)

'''
parameters = Parameters(starname = "GJ436", defaultQ = False, 
                        start = 11,
                        end = 53,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_2/loop_8/next_base_star_reg.hdf5',
                        reg_file_t = '/data/cmatthe/wobble_reg_search/GJ436_orderwise_test_2/loop_8/next_base_t_reg.hdf5',
                        output_suffix = 'orderwise_reg_search_test_2')
parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
all_parameters.append(parameters.dictionary)
'''
#parameters = Parameters(starname = "GJ436", defaultQ = False, 
                        #start = 39,
                        #end = 40,
                        #chunk_size = 1,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_5_roll1_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_5_roll1_t.hdf5',
                        #output_suffix = 'order_39_test')
#parameters.dictionary.update({"script" : "optimize_top_3.1.py"})
#all_parameters.append(parameters.dictionary)


#parameters = Parameters(starname = "GJ436", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_5_roll1_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_5_roll1_t.hdf5',
                        #output_suffix = 'def_chunk_5_roll1_n160')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/Wolf294_star_K0_orders[11,54)_regtest1406.hdf5',
                        #reg_file_t = '../wobble/regularization/Wolf294_t_K3_orders[11,54)_regtest1406.hdf5',
                        #output_suffix = 'regtest1406')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)



#default redo with new reg file matching (which should not really help for Harps reg, but lets try anyways)
#default redo with chunked readout of reg_file to simulate original implentation:


#parameters = Parameters(starname = "GJ1148", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_5_roll1_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_5_roll1_t.hdf5',
                        #output_suffix = 'def_chunk_5_roll1_n160')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "GJ436", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_5_roll1_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_5_roll1_t.hdf5',
                        #output_suffix = 'def_chunk_5_roll1_n160')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)


#parameters = Parameters(starname = "Teegarden", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_5_roll1_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_5_roll1_t.hdf5',
                        #output_suffix = 'opt3')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)

#NOTE all scripts below are only confirmed to be runnning under wob_env_2 (wobble_19_03_2019)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_5_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_5_t.hdf5',
                        #output_suffix = 'def_chunk_5_n160')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 8,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_8_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_8_t.hdf5',
                        #output_suffix = 'def_chunk_8_n160')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 16,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_16_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_16_t.hdf5',
                        #output_suffix = 'def_chunk_16_n160')
#parameters.dictionary.update({"script" : "optimize_top_3.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 60,
                        #reg_file_star =  '../wobble/regularization/def_chunk_8_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_8_t.hdf5',
                        #output_suffix = 'def_chunk_8_n60')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 300,
                        #reg_file_star =  '../wobble/regularization/def_chunk_8_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_8_t.hdf5',
                        #output_suffix = 'def_chunk_8_n300')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)
##16
#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_16_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_16_t.hdf5',
                        #output_suffix = 'def_chunk_16_n160')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 60,
                        #reg_file_star =  '../wobble/regularization/def_chunk_16_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_16_t.hdf5',
                        #output_suffix = 'def_chunk_16_n60')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 300,
                        #reg_file_star =  '../wobble/regularization/def_chunk_16_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_16_t.hdf5',
                        #output_suffix = 'def_chunk_16_n300')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)
##4
#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_4_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_4_t.hdf5',
                        #output_suffix = 'def_chunk_4_n160')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 60,
                        #reg_file_star =  '../wobble/regularization/def_chunk_4_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_4_t.hdf5',
                        #output_suffix = 'def_chunk_4_n60')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)

#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 5,
                        #niter = 300,
                        #reg_file_star =  '../wobble/regularization/def_chunk_4_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_4_t.hdf5',
                        #output_suffix = 'def_chunk_4_n300')
#parameters.dictionary.update({"script" : "optimize_top_2.py"})
#all_parameters.append(parameters.dictionary)
'''
parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  '../wobble/regularization/default_star.hdf5',
                        reg_file_t = '../wobble/regularization/default_t.hdf5',
                        output_suffix = 'def_opt2_n160')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

#default redo with chunked readout of reg_file to simulate original implentation:
parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  '../wobble/regularization/def_chunk_5_star.hdf5',
                        reg_file_t = '../wobble/regularization/def_chunk_5_t.hdf5',
                        output_suffix = 'def_chunk_5_n160')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 300,
                        reg_file_star =  '../wobble/regularization/def_chunk_5_star.hdf5',
                        reg_file_t = '../wobble/regularization/def_chunk_5_t.hdf5',
                        output_suffix = 'def_chunk_5_n300')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  '../wobble/regularization/def_chunk_10_star.hdf5',
                        reg_file_t = '../wobble/regularization/def_chunk_10_t.hdf5',
                        output_suffix = 'def_chunk_10_n160')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 300,
                        reg_file_star =  '../wobble/regularization/def_chunk_10_star.hdf5',
                        reg_file_t = '../wobble/regularization/def_chunk_10_t.hdf5',
                        output_suffix = 'def_chunk_10_n300')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  '../wobble/regularization/def_chunk_1_star.hdf5',
                        reg_file_t = '../wobble/regularization/def_chunk_1_t.hdf5',
                        output_suffix = 'def_chunk_1_n160')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 300,
                        reg_file_star =  '../wobble/regularization/def_chunk_1_star.hdf5',
                        reg_file_t = '../wobble/regularization/def_chunk_1_t.hdf5',
                        output_suffix = 'def_chunk_1_n300')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 5,
                        niter = 60,
                        reg_file_star =  '../wobble/regularization/def_chunk_5_star.hdf5',
                        reg_file_t = '../wobble/regularization/def_chunk_5_t.hdf5',
                        output_suffix = 'def_chunk_5_n60')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)
'''
'''
#dummy reg queue
parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        reg_file_star =  '../wobble/regularization/dummy_star_K0_high.hdf5',
                        reg_file_t = '../wobble/regularization/dummy_t_K3_high.hdf5',
                        output_suffix = 'dummy_reg_high')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)
'''
'''
parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        reg_file_star =  '../wobble/regularization/dummy_star_K0_low.hdf5',
                        reg_file_t = '../wobble/regularization/dummy_t_K3_low.hdf5',
                        output_suffix = 'dummy_reg_low')
parameters.dictionary.update({"script" : "optimize_top_1.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        reg_file_star =  '../wobble/regularization/dummy_star_K0_step.hdf5',
                        reg_file_t = '../wobble/regularization/dummy_t_K3_step.hdf5',
                        output_suffix = 'dummy_reg_step')
parameters.dictionary.update({"script" : "optimize_top_1.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        reg_file_star =  '../wobble/regularization/dummy_star_K0.hdf5',
                        reg_file_t = '../wobble/regularization/dummy_t_K3.hdf5',
                        output_suffix = 'dummy_reg')
parameters.dictionary.update({"script" : "optimize_top_1.py"})
all_parameters.append(parameters.dictionary)

# new (hopefully matched, and consistent) reg queue
parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 16,
                        reg_file_star =  '../wobble/regularization/Wolf294_star_K0_orders[11,54)_regtest.hdf5',
                        reg_file_t = '../wobble/regularization/Wolf294_t_K3_orders[11,54)_regtest.hdf5',
                        #HACK these are just placeholder right now because parameters.py wants to see them here, fix for case of regularization
                        output_suffix = 'regtest')
parameters.dictionary.update({"script" : "regularization_top_2.py"})
all_parameters.append(parameters.dictionary)

#TODO automate handoff of regularization file

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 16,
                        reg_file_star =  '../wobble/regularization/Wolf294_star_K0_orders[11,54)_regtest.hdf5',
                        reg_file_t = '../wobble/regularization/Wolf294_t_K3_orders[11,54)_regtest.hdf5',
                        output_suffix = 'regtest')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)
'''
''' manual reg change tests
parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 16,
                        reg_file_star =  '../wobble/regularization/reg_guess_star_0.hdf5',
                        reg_file_t = '../wobble/regularization/Wolf294_t_K3_orders[11,54)_regtest.hdf5',
                        output_suffix = 'reg_star_guess')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 16,
                        reg_file_star =  '../wobble/regularization/reg_guess_star_0.hdf5',
                        reg_file_t = '../wobble/regularization/reg_guess_t_0.hdf5',
                        output_suffix = 'reg_both_guess')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 54,
                        chunk_size = 16,
                        reg_file_star =  '../wobble/regularization/dummy_star_K0.hdf5',
                        reg_file_t = '../wobble/regularization/reg_guess_t_0.hdf5',
                        output_suffix = 'reg_t_guess')
parameters.dictionary.update({"script" : "optimize_top_2.py"})
all_parameters.append(parameters.dictionary)
'''
#print(all_parameters)


length = len(all_parameters)

start_time = time()

for i, params in enumerate(all_parameters):
    parameters = Parameters(starname = "placeholder")
    parameters.dictionary.update(params)
    parameters.write("yaml_temp/optimize_parameters.yaml")
    
    script_name = parameters.dictionary["script"]
    
    
    os.system("python3 "  + script_name)
    j = i+1
    print("script{0}/{1} finished".format(j, length))
    print("queue time elapsed: {0:.2f} min, {1:.2f} h".format((time() - start_time)/60.0, (time() - start_time)/3600.0))
