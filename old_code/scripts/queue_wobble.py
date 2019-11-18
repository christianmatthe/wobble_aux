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



#scripts = ["optimize_top_1.py"]
#length = len(scripts)
#list_of_parameter_dictionaries
# #TODO add script 
#and regularization selection parameters to parameters class, optimize script needs to be adjusted to react

all_parameters = []

parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        start = 11,
                        end = 53,
                        chunk_size = 5,
                        niter = 160,
                        reg_file_star =  '../wobble/regularization/Wolf294_star_K0_orders[11,54)_regtest1406.hdf5',
                        reg_file_t = '../wobble/regularization/Wolf294_t_K3_orders[11,54)_regtest1406.hdf5',
                        output_suffix = 'regtest1406_opt4')
parameters.dictionary.update({"script" : "optimize_top_4.py"})
all_parameters.append(parameters.dictionary)


#test whether new reg wll run on  old files
#parameters = Parameters(starname = "Wolf294", defaultQ = False, 
                        #start = 11,
                        #end = 54,
                        #chunk_size = 16,
                        #reg_file_star =  '../wobble/regularization/Wolf294_star_K0_orders[11,54)_regtest1406.hdf5',
                        #reg_file_t = '../wobble/regularization/Wolf294_t_K3_orders[11,54)_regtest1406.hdf5',
                        ##HACK these are just placeholder right now because parameters.py wants to see them here, fix for case of regularization
                        #output_suffix = 'regtest1406')
#parameters.dictionary.update({"script" : "regularization_top_2.py"})
#all_parameters.append(parameters.dictionary)


#default redo with new reg file matching (which should not really help for Harps reg, but lets try anyways)
#default redo with chunked readout of reg_file to simulate original implentation:

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

#parameters = Parameters(starname = "Teegarden", defaultQ = False, 
                        #start = 11,
                        #end = 53,
                        #chunk_size = 5,
                        #niter = 160,
                        #reg_file_star =  '../wobble/regularization/def_chunk_5_roll1_star.hdf5',
                        #reg_file_t = '../wobble/regularization/def_chunk_5_roll1_t.hdf5',
                        #output_suffix = 'opt4')
#parameters.dictionary.update({"script" : "optimize_top_4.py"})
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
