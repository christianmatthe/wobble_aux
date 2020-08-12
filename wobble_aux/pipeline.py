import make_data_carmenes as mdc 
import run_wobble as rw
import process_results as pr

import os
from tqdm import tqdm

if __name__ == "__main__":
    #setup dictionaries
    name_dict = name_dict_master = {
        
        "GJ436"     : "J11421+267",
        "GJ1148"    : "J11417+427",
        "GJ3473"    : "J08023+033",
        "YZ Cet"    : "J01125-169",
        "GJ15A"     : "J00183+440",
        "GJ176"     : "J04429+189",
        "GJ536"     : "J14010-026",
        
        "GJ581"     : "J15194-077", #seems to actually work fine didn't download data #TODO Retest this data seems present to me
        
        "GJ3512"    : "J08413+594", #issues due to low min_snr, drops all orders at snr 60
        "Wolf294"   : "J06548+332",
        "GJ876"     : "J22532-142",
        "Teegarden" : "J02530+168",
        "Barnard"   : "J17578+046"
        
        }
    
    name_dict_nir = name_dict_master_nir = {
        
        "GJ436"     : "J11421+267",
        "GJ1148"    : "J11417+427",
        
        "GJ3512"    : "J08413+594", #issues due to low min_snr, drops all orders at snr 60
        
        }

    
    simbad_dict = {
        "GJ436"     : "GJ436",
        "GJ1148"    : "GJ1148",
        "GJ3473"    : "G 50-16",
        "YZ Cet"    : "YZ Cet",
        "GJ15A"     : "GJ15A" ,
        "GJ176"     : "GJ176" ,
        "GJ536"     : "GJ536",
        "GJ581"     : "GJ581",
        "GJ3512"    : "GJ3512",
        "Wolf294"   : "Wolf294",
        "GJ876"     : "GJ876" ,
        "Teegarden" : "GAT 1370",
        "Barnard"   : "GJ699"
        
        }
    
    #Denote default niters
    #TODO implement this as default dictionary below
    niter_dict = {
    
    "GJ436"     : 160, #default
    "GJ1148"    : 160, #default
    "GJ3473"    : 160, #default
    "YZ Cet"    : 160, #default
    "GJ15A"     : 160, #default
    "GJ176"     : 160, #default
    "GJ536"     : 160, #default
    
    "GJ581"     : 160, #default
    
    "GJ3512"    : 300, #due to larger than usual RV amplitude TODO check exact value
    "Wolf294"   : 160, #default
    "GJ876"     : 700, #due to larger than usual RV amplitude NOTE above 600-700 fit history is sometimes erratic 
    "Teegarden" : 160, #default
    "Barnard"   : 160, #default
    
    }
    
    
    #BEGIN make_data_carmenes
    '''
    data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    #data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/test/"#NOTE testing issues only has GJ536 data
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #read data already includes name dictionary
    
    
    
    for star in tqdm(name_dict):
        starname = star
        simbad_name = simbad_dict[star]
        arm = "vis"
        mdc.make_data(starname, arm, data_directory, simbad_name = simbad_name, serval_dir = serval_dir, nzp_shift = True, drift_shift = True
    '''
    #END make_data_carmenes
        
    #run wobble below
    
    # results_dir_base = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/" 
    results_dir_base = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline_2/" #Change to v2 indicates change to default continuum nsigma. now [0.3,3] for vis
    
    '''
     #first test  
    #BEGIN baseline_0   
    pipeline_run_name = "pipeline_test_0/"
    results_dir = results_dir_base + pipeline_run_name
    os.makedirs(results_dir, exist_ok = True)
    
    
    for star in tqdm(name_dict):
        parameters = rw.Parameters(starname = star,
                                results_dir = results_dir,
                                
                            data_suffix = "_vis_drift+nzp",
                            start = 11,
                            end = 53,
                            chunk_size = 5,
                            niter = 160,
                            reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            output_suffix = "baseline_0",
                            plot_continuum = False)
        rw.run_wobble(parameters)
        #TODO automatically write file list fron  this?
    #END baseline_0
    '''
     
    
    # # BEGIN niter test
    # niter_list = [40,80,
    #     120,160,200] 
    # for niter in niter_list:
    #     pipeline_run_name = "pipeline_n{}/".format(niter)
    #     results_dir = results_dir_base + pipeline_run_name
    #     os.makedirs(results_dir, exist_ok = True)
    #     for star in tqdm(name_dict):
    #         parameters = rw.Parameters(starname = star,
    #                                 results_dir = results_dir,
                                    
    #                             data_suffix = "_vis_drift+nzp",
    #                             start = 11,
    #                             end = 53,
    #                             chunk_size = 5,
    #                             niter = niter,
    #                             reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
    #                             reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
    #                             output_suffix = "n_{}".format(niter),
    #                             plot_continuum = False)
    #         rw.run_wobble(parameters)
    ## END niter test
    
    
    
    ###BEGIN GJ876 more niter test
    ##name_dict = {"GJ876" : "J22532-142"}
    #name_dict = {"GJ3512"    : "J08413+594"}
    #niter_list = [250,300,350,400,450,500,600]
    ##niter_list = [600,700,800,900,1000] 
    #for niter in niter_list:
        #pipeline_run_name = "pipeline_n{}/".format(niter)
        #results_dir = results_dir_base + pipeline_run_name
        #os.makedirs(results_dir, exist_ok = True)
        #for star in tqdm(name_dict):
            #parameters = rw.Parameters(starname = star,
                                    #results_dir = results_dir,
                                    
                                #data_suffix = "_vis_drift+nzp",
                                #start = 11,
                                #end = 53,
                                #chunk_size = 5,
                                #niter = niter,
                                #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                                #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                                #output_suffix = "n_{}".format(niter),
                                #plot_continuum = False)
            #rw.run_wobble(parameters)
    ###END GJ876 more niter test
    
    ##BEGIN no snr cutting (min snr 5) test
    #name_dict = name_dict_master
    ##name_dict = {"Teegarden" : "J02530+168"}
    ##name_dict = {"GJ3512"    : "J08413+594"}
    
    
    #name_dict = name_dict_master
    #del name_dict["GJ3512"]
    
    #snr_list = [#5
                ##,30,15,
                #45,60
                #] 
    #for min_snr in snr_list:
        #pipeline_run_name = "pipeline_min_snr{}/".format(min_snr)
        #results_dir = results_dir_base + pipeline_run_name
        #os.makedirs(results_dir, exist_ok = True)
        #for star in tqdm(name_dict):
            #parameters = rw.Parameters(starname = star,
                                    #results_dir = results_dir,
                                    
                                    #min_snr = min_snr,
                                    
                                #data_suffix = "_vis_drift+nzp",
                                #start = 11,
                                #end = 53,  #CHNGE TIS BACK TO 53
                                #chunk_size = 5,
                                #niter = 160,
                                #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                                #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                                #output_suffix = "min_snr{}".format(min_snr),
                                #plot_continuum = False)
            #rw.run_wobble(parameters)
    #END no snr cutting (min snr 5) test
    
    ##BEGIN no snr cutting (min snr 5) test, with Niter list
    #name_dict = name_dict_master
    ##name_dict = {"Teegarden" : "J02530+168"}
    ##name_dict = {"GJ3512"    : "J08413+594"}
    #name_dict = {"GJ876" : name_dict_master["GJ876"]}
    
    
    #name_dict = name_dict_master
    #del name_dict["GJ3512"]
    
    #snr_list = [5
                #,30,15,
                #45,60
                #] 
    #for min_snr in snr_list:
        #pipeline_run_name = "pipeline_min_snr{}/".format(min_snr)
        #results_dir = results_dir_base + pipeline_run_name
        #os.makedirs(results_dir, exist_ok = True)
        #for star in tqdm(name_dict):
            #try:
                #niter = niter_dict[star]
            #except:
                #print(star, " not in  niter_dict, using niter = 160 as default")
                #niter = 160
        
            #parameters = rw.Parameters(starname = star,
                                    #results_dir = results_dir,
                                    
                                    #min_snr = min_snr,
                                    
                                #data_suffix = "_vis_drift+nzp",
                                #start = 11,
                                #end = 53,  #CHNGE TIS BACK TO 53
                                #chunk_size = 5,
                                #niter = 160,
                                #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                                #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                                #output_suffix = "min_snr{}".format(min_snr),
                                #plot_continuum = False)
            #rw.run_wobble(parameters)
    ##END no snr cutting (min snr 5) test


    ##BEGIN Prototype Niter_dict
    ##pipeline_run_name = "pipeline_niter_dict/"
    #pipeline_run_name = "pipeline_niter_dict_continuum_test_[0.5,1]/"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    #name_dict = name_dict_master
    ##name_dict = {"GJ876" : name_dict_master["GJ876"]}
    ##name_dict = {"GJ436" : name_dict_master["GJ436"]}
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_vis_drift+nzp",
                            #start = 11, #CHNGE TIS BACK TO 11
                            #end = 53,  #CHNGE TIS BACK TO 53
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #output_suffix = "n_{}".format(niter),
                            ##continuum_order = 1,  # NOTE
                            #continuum_nsigma = [0.5,1],# NOTE
                            #plot_continuum = True)
        #rw.run_wobble(parameters)
    ###END 
    
    ##BEGIN Prototype Niter_dict
    #pipeline_run_name = "pipeline_niter_dict/"
    ##pipeline_run_name = "pipeline_continuum_test_[0.3,1]/"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    #name_dict = name_dict_master
    ##name_dict = {"Teegarden" : name_dict_master["Teegarden"]}
    ##name_dict = {"GJ436" : name_dict_master["GJ436"]}
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_vis_drift+nzp",
                            #start = 11, #CHNGE TIS BACK TO 11
                            #end = 53,  #CHNGE TIS BACK TO 53
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #output_suffix = "n_{}".format(niter),
                            #continuum_order = 1,  # NOTE
                            ##continuum_nsigma = [0.3,1],# NOTE
                            #plot_continuum = True)
        #rw.run_wobble(parameters)
    ##END 
    
    #BEGIN Continuum tests
    c_order_list = [0,1,2,3,4,5,6]
    for continuum_order in c_order_list:
        pipeline_run_name = "pipeline_continuum_order_{}".format(continuum_order)
        results_dir = results_dir_base + pipeline_run_name
        os.makedirs(results_dir, exist_ok = True)
        name_dict = name_dict_master
        #name_dict = {"Teegarden" : name_dict_master["Teegarden"]}
        #name_dict = {"GJ436" : name_dict_master["GJ436"]}
        for star in tqdm(name_dict):
            try:
                niter = niter_dict[star]
            except:
                print(star, " not in  niter_dict, using niter = 160 as default")
                niter = 160
            parameters = rw.Parameters(starname = star,
                                    results_dir = results_dir,
                                    
                                    min_snr = 60,
                                    
                                data_suffix = "_vis_drift+nzp",
                                start = 11, #CHNGE TIS BACK TO 11
                                end = 53,  #CHNGE TIS BACK TO 53
                                chunk_size = 5,
                                niter = niter,
                                reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                                reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                                output_suffix = "",
                                continuum_order = continuum_order,  # NOTE
                                continuum_nsigma = [0.3,3],# NOTE
                                plot_continuum = True)
            rw.run_wobble(parameters)
    #END 
    
    
        
    #BEGIN No drift correction tests    
    ##BEGIN make_data_carmenes
    #data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    ##data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/test/"#NOTE testing issues only has GJ536 data
    #serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #read data already includes name dictionary
    
    
    #name_dict = name_dict_master
    #for star in tqdm(name_dict):
        #starname = star
        #simbad_name = simbad_dict[star]
        #arm = "vis"
        #mdc.make_data(starname, arm, data_directory, simbad_name = simbad_name, serval_dir = serval_dir, nzp_shift = False, drift_shift = False)
    
    ##END make_data_carmenes
    ##run regular opti on all stars without drift
    #pipeline_run_name = "pipeline_no_drift_corr/"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    #name_dict = name_dict_master
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_vis",
                            #start = 11,
                            #end = 53,  #CHNGE TIS BACK TO 53
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #output_suffix = "no_drift_corr",
                            #plot_continuum = False)
        #rw.run_wobble(parameters)
    ##END No drift correction tests  
    
    #######################################################################################################################################################################################
    #Reruns for powerpoint pres recreation section
    
    ##BEGIN GJ436 NIR with new and old continnuum
    
    ##BEGIN make_data_carmenes
    ##BEGIN No drift correction tests   
    #data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    #serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_NIR/" #read data already includes name dictionary
    
    
    #name_dict = {"GJ436" : name_dict_master_nir["GJ436"]}
    #for star in tqdm(name_dict):
        #starname = star
        #simbad_name = simbad_dict[star]
        #arm = "nir"
        #mdc.make_data(starname, arm, data_directory, simbad_name = simbad_name, serval_dir = serval_dir, nzp_shift = False, drift_shift = False)
    ##END make_data_carmenes
    
    ##NOTE "nir_nodrift_cont_5_[0.3,3]"
    ##May not work the same due to now implemented telluric dropping
    #continuum_order = 5
    #continuum_nsigma = [0.3,3]
    #pipeline_run_name = "pipeline_pp_reruns"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    #name_dict = {"GJ436" : name_dict_master_nir["GJ436"]}
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_nir_split",
                            #start = 0, 
                            #end = 56, 
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #continuum_order = continuum_order,  # NOTE
                            #continuum_nsigma = continuum_nsigma,# NOTE
                            #output_suffix = "nir_nodrift_cont_{0}_{1}".format(continuum_order, continuum_nsigma),
                            #plot_continuum = True)
        #rw.run_wobble(parameters)
        
    ##NOTE "nir_nodrift_cont_1_[0.5,1]" 
    ##May not work the same due to now fixed implementation issues
    #continuum_order = 1
    #continuum_nsigma = [0.5,1]
    #pipeline_run_name = "pipeline_pp_reruns"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    ##name_dict = name_dict_master
    #name_dict = {"GJ436" : name_dict_master_nir["GJ436"]}
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_nir_split",
                            #start = 0, 
                            #end = 56,  
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #continuum_order = continuum_order,  # NOTE
                            #continuum_nsigma = continuum_nsigma,# NOTE
                            #output_suffix = "nir_nodrift_cont_{0}_{1}".format(continuum_order, continuum_nsigma),
                            #plot_continuum = True)
        #rw.run_wobble(parameters)
        
        
    ##BEGIN make_data_carmenes
    #data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    #serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_NIR/" #read data already includes name dictionary
    
    
    #name_dict = {"GJ436" : name_dict_master_nir["GJ436"]}
    #for star in tqdm(name_dict):
        #starname = star
        #simbad_name = simbad_dict[star]
        #arm = "nir"
        #mdc.make_data(starname, arm, data_directory, simbad_name = simbad_name, serval_dir = serval_dir, nzp_shift = True, drift_shift = True)
    ##END make_data_carmenes
    
    
    ##NOTE "nir_cont_5_[0.3,3]"
    ##May not work the same due to now implemented telluric dropping
    #continuum_order = 5
    #continuum_nsigma = [0.3,3]
    #pipeline_run_name = "pipeline_pp_reruns"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    #name_dict = {"GJ436" : name_dict_master_nir["GJ436"]}
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_nir_drift+nzp_split",
                            #start = 0, 
                            #end = 56, 
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #continuum_order = continuum_order,  # NOTE
                            #continuum_nsigma = continuum_nsigma,# NOTE
                            #output_suffix = "nir_cont_{0}_{1}".format(continuum_order, continuum_nsigma),
                            #plot_continuum = True)
        #rw.run_wobble(parameters)
        
    ##NOTE "nir_nodrift_cont_1_[0.5,1]" 
    ##May not work the same due to now fixed implementation issues
    #continuum_order = 1
    #continuum_nsigma = [0.5,1]
    #pipeline_run_name = "pipeline_pp_reruns"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    ##name_dict = name_dict_master
    #name_dict = {"GJ436" : name_dict_master_nir["GJ436"]}
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_nir_drift+nzp_split",
                            #start = 0, 
                            #end = 56,  
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #continuum_order = continuum_order,  # NOTE
                            #continuum_nsigma = continuum_nsigma,# NOTE
                            #output_suffix = "nir_cont_{0}_{1}".format(continuum_order, continuum_nsigma),
                            #plot_continuum = True)
        #rw.run_wobble(parameters)
        
    ###END 
    
##BEGIN GJ436 and GJ1148 with no (dummy) reg, default reg, wobble_reg, GJ436 regs #TODO do this with loop
    #reg_dict = {
        #'no_reg'      : ['regularization/dummy_star_K0_no_reg.hdf5', 
                         #'regularization/dummy_t_K3_no_reg.hdf5'
                         #],
        #'default_reg' : [
                         #'../wobble_19_03_2019/wobble/wobble/regularization/default_star.hdf5',
                         #'../wobble_19_03_2019/wobble/wobble/regularization/default_t.hdf5'
                         #],
        #'wobble_reg'  : [
                        #'../../wobble_data/wobble/regularization/GJ1148_star_K0_orders[11,53)_stitched_reformatted.hdf5'
                        #,
                         #'../../wobble_data/wobble/regularization/GJ1148_t_K3_orders[11,53)_stitched_reformatted.hdf5'
                         #]
        #}
    #for i in range(0,6):
        #reg_dict['l{0}_reg'.format(i)] = [
        #'../../wobble_reg_search/GJ436_orderwise_snr+drift_shift_1/loop_{0}/next_base_star_reg.hdf5'.format(i),
        #'../../wobble_reg_search/GJ436_orderwise_snr+drift_shift_1/loop_{0}/next_base_t_reg.hdf5'.format(i)
        #]
    #for key, reg in reg_dict.items():
        ##May not work the same due to now fixed implementation issues
        #continuum_order = 1 # NOTE
        #continuum_nsigma = [0.3,3] # NOTE
        #pipeline_run_name = "pipeline_pp_reruns/reg"
        #results_dir = results_dir_base + pipeline_run_name
        #os.makedirs(results_dir, exist_ok = True)
        ##name_dict = name_dict_master
        #name_dict = {"GJ436" : name_dict_master["GJ436"],
                    #"GJ1148" : name_dict_master["GJ1148"]
                    #}
        #for star in tqdm(name_dict):
            #try:
                #niter = niter_dict[star]
            #except:
                #print(star, " not in  niter_dict, using niter = 160 as default")
                #niter = 160
            #parameters = rw.Parameters(starname = star,
                                    #results_dir = results_dir,
                                    
                                    #min_snr = 60,
                                    
                                #data_suffix = "_vis_drift+nzp",
                                #start = 11, 
                                #end = 53,  
                                #chunk_size = 5,
                                #niter = niter,
                                #reg_file_star =  reg[0],
                                #reg_file_t = reg[1],
                                #continuum_order = continuum_order,  
                                #continuum_nsigma = continuum_nsigma,
                                #output_suffix = "{0}".format(key),
                                #plot_continuum = False)
            #rw.run_wobble(parameters)
##END
#BEGIN Baseline+SNR: GJ436 without SNR and without drift
##May not work the same due to now fixed implementation issues
    #snr_list = [5,15,30,45,60]
    #for snr in snr_list:
        #continuum_order = 1 # NOTE
        #continuum_nsigma = [0.3,3] # NOTE
        #pipeline_run_name = "pipeline_pp_reruns/Baseline+SNR"
        #results_dir = results_dir_base + pipeline_run_name
        #os.makedirs(results_dir, exist_ok = True)
        ##name_dict = name_dict_master
        #name_dict = {"GJ436" : name_dict_master["GJ436"],
                    #"GJ1148" : name_dict_master["GJ1148"]
                    #}
        #for star in tqdm(name_dict):
            #try:
                #niter = niter_dict[star]
            #except:
                #print(star, " not in  niter_dict, using niter = 160 as default")
                #niter = 160
            #parameters = rw.Parameters(starname = star,
                                    #results_dir = results_dir,
                                    
                                    #min_snr = snr,
                                    
                                #data_suffix = "_vis",
                                #start = 11, 
                                #end = 53,  
                                #chunk_size = 5,
                                #niter = niter,
                                #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                                #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                                #continuum_order = continuum_order,  
                                #continuum_nsigma = continuum_nsigma,
                                #output_suffix = "min_snr_{}".format(snr),
                                #plot_continuum = True)
            #rw.run_wobble(parameters)
##END Baseline+SNR: GJ436 without SNR and without drift
##BEGIN GJ3512 as example for completely "bad" epoch

    #continuum_order = 1 # NOTE
    #continuum_nsigma = [0.3,3] # NOTE
    #pipeline_run_name = "pipeline_pp_reruns/GJ3512"
    #results_dir = results_dir_base + pipeline_run_name
    #os.makedirs(results_dir, exist_ok = True)
    ##name_dict = name_dict_master
    #name_dict = {"GJ3512" : name_dict_master["GJ3512"]
                #}
    #for star in tqdm(name_dict):
        #try:
            #niter = niter_dict[star]
        #except:
            #print(star, " not in  niter_dict, using niter = 160 as default")
            #niter = 160
        #parameters = rw.Parameters(starname = star,
                                #results_dir = results_dir,
                                
                                #min_snr = 60,
                                
                            #data_suffix = "_vis_drift+nzp",
                            #start = 11, 
                            #end = 53,  
                            #chunk_size = 5,
                            #niter = niter,
                            #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                            #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                            #continuum_order = continuum_order,  
                            #continuum_nsigma = continuum_nsigma,
                            #output_suffix = "n_{}".format(niter),
                            #plot_continuum = True)
        #rw.run_wobble(parameters)
        #Reruns for powerpoint pres recreation section
##END Baseline+SNR: GJ436 without SNR and without drift
        
##BEGIN TODO Baseline+SNR: GJ1148 o49 e31 continuum norm (order5) error 
###May not work the same due to now fixed implementation issues
    #snr_list = [5]
    #for snr in snr_list:
        #continuum_order = 5 # NOTE
        #continuum_nsigma = [0.3,3] # NOTE
        #pipeline_run_name = "pipeline_pp_reruns/Baseline+SNR"
        #results_dir = results_dir_base + pipeline_run_name
        #os.makedirs(results_dir, exist_ok = True)
        ##name_dict = name_dict_master
        #name_dict = {#"GJ436" : name_dict_master["GJ436"],
                    #"GJ1148" : name_dict_master["GJ1148"]
                    #}
        #for star in tqdm(name_dict):
            #try:
                #niter = niter_dict[star]
            #except:
                #print(star, " not in  niter_dict, using niter = 160 as default")
                #niter = 160
            #parameters = rw.Parameters(starname = star,
                                    #results_dir = results_dir,
                                    
                                    #min_snr = snr,
                                    
                                #data_suffix = "_vis",
                                #start = 49, 
                                #end = 50,  
                                #chunk_size = 5,
                                #niter = niter,
                                #reg_file_star =  'regularization/dummy_star_K0_no_reg.hdf5',
                                #reg_file_t = 'regularization/dummy_t_K3_no_reg.hdf5',
                                #continuum_order = continuum_order,  
                                #continuum_nsigma = continuum_nsigma,
                                #output_suffix = "cont_{}".format(continuum_order),
                                #plot_continuum = True)
            #rw.run_wobble(parameters)
##END TODO Baseline+SNR: GJ1148 o49 e31 continuum norm (order5) error 
    #######################################################################################################################################################################################
    #####################################
    
    
#use eval results after this
    
    
