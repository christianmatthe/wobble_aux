import make_data_carmenes as mdc 
import run_wobble as rw
import process_results as pr

import os
from tqdm import tqdm

if __name__ == "__main__":
    #setup dictionaries
    name_dict = {
        
        #"GJ436"     : "J11421+267",
        #"GJ1148"    : "J11417+427",
        #"GJ3473"    : "J08023+033",
        #"YZ Cet"    : "J01125-169",
        #"GJ15A"     : "J00183+440",
        #"GJ176"     : "J04429+189",
        #"GJ536"     : "J14010-026",
        
        #"GJ581"     : "J15194-077",
        
        "GJ3512"    : "J08413+594", #issues due to low min_snr, drops all orders at snr 60
        "Wolf294"   : "J06548+332",
        "GJ876"     : "J22532-142",
        "Teegarden" : "J02530+168",
        "Barnard"   : "J17578+046"
        
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
    '''
    #make_data_carmenes
    #BEGIN parameters mdc
    data_directory = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/"
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #read data already includes name dictionary
    #END parameters mdc
    
    for star in tqdm(name_dict):
        starname = star
        simbad_name = simbad_dict[star]
        arm = "vis"
        mdc.make_data(starname, arm, data_directory, simbad_name = simbad_name, serval_dir = serval_dir, nzp_shift = True)
    '''    
        
    #run wobble
    #BEGIN parameters rw
    results_dir_base = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/" 
    pipeline_run_name = "pipeline_test_0/"
    results_dir = results_dir_base + pipeline_run_name
    os.makedirs(results_dir, exist_ok = True)
    #END parameters rw
    
    
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
        
#use eval results after this
    
    
