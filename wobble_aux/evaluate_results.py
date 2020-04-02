#Goal is to us RVmod (from exostriker) as a auto fit library in order to compare multiple results file with respect to their rms to their auto fit.
#Heavily based on exostriker/Notebook_and_script_examples/RVmod_as_py_lib_auto_fit_example.py by Trifon Trifonov
'''
Goals:
1. Import from list of files
2. process_results
    a) load parameters object from results_file
3. auto fit 
    a) save fit parameters and rms
    b) plot commparison plots
        i) make this modular "e.g. if  plot_error_vs_corrections == True:"
        ii) Use the PROPER order numbers (for AQ)
 

'''
import run_wobble as rw
import process_results as pr

if __name__ == "__main__": 
    
    

    file_list = ["results_GJ436_Kstar0_Kt3_git_run_wobble_test1.hdf5","results_GJ436_Kstar0_Kt3_git_run_wobble_test0.hdf5" ] #laptop sample file
    
    ########### 
    
    results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/"
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" +"../data/servaldir/CARM_VIS/" #NOTE only for VIS
    
    
    ''' THis Dictionary does not always follow the namming connventions I used for "starname"
    #laptop test example:
    names = pd.read_csv(os.path.dirname(os.path.abspath(__file__)) + '/carmenes_aux_files/name_conversion_list.csv')
    name_dict_inverse = dict(zip(names['#Karmn'], names['Name'])) # yields Carm ID for catalogue name
    '''
    
    #HACK Either write Carm ID into results file or make a more permanent name dict
    # dictionary connecting results_file["parameters"].attrs["pkl"] -> parameters.starname to CARMENES ID for serval results matching
    name_dict = {
        "GJ436"  : "J11421+267",
        "GJ3473" : "J08023+033"
        
        }
    
    simbad_dict = {
        "GJ436"  : "GJ436",
        "GJ3473" : "G 50-16"
        
        }
    
    for f in file_list:
        wobble_file = results_dir + f
        
        parameters = rw.read_parameters_from_results(wobble_file)
        
        #basic functionality:
        #1. write functions to make basic plot with Serval and Wobble data and auto fit, and rms
        #2. plot all rms into one plot (for same star) 
        
        
    
    
