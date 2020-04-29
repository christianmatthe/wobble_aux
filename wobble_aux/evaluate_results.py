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
#temp
import h5py
#

import os
import run_wobble as rw
import process_results as pr

#imports fromm RVmod_as_py_lib_auto_fit_example
import sys 
sys.path.append('/home/christian/Documents/exostriker/lib/') #RV_mod directory must be in your path
sys.path.append('/home/christian/Documents/exostriker/lib/RV_mod/')
#on LX39
#sys.path.append('/data/cmatthe/python/exostriker/lib/')
#sys.path.append('/data/cmatthe/python/exostriker/lib/RV_mod/')
import RV_mod as rv
import gls as gls
import numpy as np
import dill

#for plots
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl

import numpy as np

def vels_to_kep_fit(dataset_name, vels_file):
    fit=rv.signal_fit()
    #option include e.g. name = starname to identify session

    #fit.cwd = '../' # NOTE it is also important that the ES current working directory (cwd) point to the "lib" directory. This will be fixed in future releases 
    fit.cwd = '/home/christian/Documents/exostriker/' # HACK on laptop
    #fit.cwd = '/data/cmatthe/python/exostriker/' # OnLX39
    
    
    #write RVmod auto fit example into a function
    #takes vels file and auto fits planets to it. outputs kep_fit object
    fit.add_dataset(dataset_name, vels_file, 0.0,0.0)  # the last two entries are initial offset and jitter
    
    # Lets not fit for jitters now, i.e. keep at the initial value of 0 m/s
    fit.use.use_jitters[0] = False
    fit.use.use_jitters[1] = False
    
    #  Run it once to find the RV offsets, no planets yet.
    fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=False, amoeba_starts=20, print_stat=False)
    # Run GLS once 
    rv.run_gls(fit)
    #TEST rv.run_gls("jo",123)
    rv.run_gls_o_c(fit)
        
    # now lets find the planets in our data!
    fit.auto_fit_max_pl = 2
    fit.auto_fit_FAP_level = 0.001 # this corresponds to FAP = 0.1%. GLS power with FAP below that level we take as planet candidate.
    fit.auto_fit_allow_ecc = True # otherwise will look only for circular planets

    fit = rv.find_planets(fit)
    
    # Lets fit one more time with RV jitters modeled
    fit.use.use_jitters[0] = True
    fit.use.use_jitters[1] = True

    fit.fitting(minimize_loglik=True) # This will run the Fortran Simplex, which optimizes the lnL and thus can fit for the jitters (but no errors so far)
    fit.fitting(minimize_loglik=False) # Now that the jiters are known we can get to the L-M method to gett Covar.Matrix par. uncertanties.
    
    #Lets print the best fit params  
    print("Loglik = %s"%fit.loglik)
    fit.print_info() #this is an obsolete function call, will be replaced!
    
    #lets copy the Keplerian object, we will need it later for plotting
    kep_fit = dill.copy(fit)
    
    return kep_fit

def plot_time_series(kep_fit, output_file):
    ###### For nice plotting ##############

    mpl.rcParams['axes.linewidth'] = 2.0 #set the value globally
    mpl.rcParams['xtick.major.pad']='8'
    mpl.rcParams['ytick.major.pad']='2'

    # set tick width
    mpl.rcParams['xtick.major.size'] = 8
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 5
    mpl.rcParams['xtick.minor.width'] = 2

    mpl.rcParams['ytick.major.size'] = 8
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['ytick.minor.size'] = 5
    mpl.rcParams['ytick.minor.width'] = 2
    
    # HACK usetex fails (probably because no latex is installed) mpl.rc('text',usetex=True)
    mpl.rc('text',usetex=False)
    font = {'family' : 'normal','weight' : 'bold','size'   : 18,'serif':['Helvetica']}
    mpl.rc('font', **font)
    ##### time series format ######
    f = plt.figure(0, figsize=(8,6.5))
    plt.subplots_adjust(hspace=0.005)
    format_im = 'pdf' #'pdf' or png
    dpi = 300

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    #gs.update(  wspace=0.05)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    
    color = ['b', 'r', 'g', 'r']
    symbol = ['o', 'D', 'o', 'o']
    markersize = [5, 5, 6, 6] 
    alpha = [1, 1, 1, 1]

    model_color = 'k'
    model_lw = '1.0'
    #### Get the time series (these below are self explanatory) ########     
    jd        = kep_fit.fit_results.rv_model.jd
    rvs       = kep_fit.fit_results.rv_model.rvs
    rv_err    = kep_fit.fit_results.rv_model.rv_err
    o_c       = kep_fit.fit_results.rv_model.o_c
    
    data_set  = kep_fit.filelist.idset

    # we can add the jitter
    add_jitter = True
    if add_jitter == True:
        rv_err = np.array([np.sqrt(rv_err[i]**2 + kep_fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])



    # Kep model time series #
    kep_model_x = kep_fit.fit_results.model_jd
    #HACK unknow issue with  the above makes it not equal to the below?
    #kep_model_x = np.linspace(min(jd), max(jd), 1000) #obj.fit_results.model_jd # 1000 is same dimensions but produces incorrect plots
    kep_model_y = kep_fit.fit_results.model

    

    ###################################################################

    offset_pre  = 250
    offset_post = 250

    zero_point_T = range((int(min(jd))-offset_pre),(int(max(jd))+offset_post),10)
    zero_point   = np.zeros(len(zero_point_T))


    ax1.plot(kep_model_x, kep_model_y,       '-', linewidth=model_lw, color=model_color)
    ax2.plot(zero_point_T,zero_point,'-', linewidth=model_lw, color=model_color)      

    

    for i in range(len(data_set)):

            ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)
            ax2.errorbar(jd[i],o_c[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])],color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)



    
    ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 
    ax1.set_xlim(min(jd)-offset_pre,max(jd)+offset_post)
    

    ax2.set_xlabel(r'JD [day]',fontsize=16)
    ax2.set_ylabel(r'o$-$c  [m/s]',fontsize=16, rotation = 'vertical') 
    ax2.set_xlim(min(jd)-offset_pre,max(jd)+offset_post)

    ax2.locator_params(axis="x", nbins=9)
    plt.setp( ax2.get_yticklabels(), fontsize=15,weight='bold')
    plt.setp( ax2.get_xticklabels(), fontsize=15,weight='bold')
    
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) 


    # HACK plt.savefig('RV_plot_example_time_series.%s'%(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    plt.savefig(output_file + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla() 
    ax2.cla()
    
    
    
    
        

if __name__ == "__main__": 
    
    
    #run_name = "test_1"
    #file_list = ["results_GJ436_Kstar0_Kt3_laptop_example_0.hdf5", "results_GJ3473_Kstar0_Kt3_laptop_example_1.hdf5" ] #laptop sample file
    
    #run_name = "test_2"
    #file_list = ["results_GJ1148_Kstar0_Kt3_eval_example_2.hdf5","results_GJ436_Kstar0_Kt3_laptop_example_0.hdf5", "results_GJ3473_Kstar0_Kt3_laptop_example_1.hdf5"]
    
    run_name ="baseline_0"
    file_list = ["results_GJ436_Kstar0_Kt3_baseline_0.hdf5",
        "results_GJ1148_Kstar0_Kt3_baseline_0.hdf5",
        "results_GJ3473_Kstar0_Kt3_baseline_0.hdf5",
        "results_YZ Cet_Kstar0_Kt3_baseline_0.hdf5",
        "results_GJ15A_Kstar0_Kt3_baseline_0.hdf5",
        "results_GJ176_Kstar0_Kt3_baseline_0.hdf5", "results_GJ536_Kstar0_Kt3_baseline_0.hdf5", "results_GJ3512_Kstar0_Kt3_baseline_0.hdf5", "results_Wolf294_Kstar0_Kt3_baseline_0.hdf5", "results_GJ876_Kstar0_Kt3_baseline_0.hdf5" , "results_Teegarden_Kstar0_Kt3_baseline_0.hdf5", "results_Barnard_Kstar0_Kt3_baseline_0.hdf5"
                 ]
    
    ########### 
    
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/"#
    results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    output_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/evaluate/{0}/".format(run_name)
    os.makedirs(output_dir, exist_ok = True)
    
    
    ''' THis Dictionary does not always follow the namming connventions I used for "starname"
    #laptop test example:
    names = pd.read_csv(os.path.dirname(os.path.abspath(__file__)) + '/carmenes_aux_files/name_conversion_list.csv')
    name_dict_inverse = dict(zip(names['#Karmn'], names['Name'])) # yields Carm ID for catalogue name
    '''
    
    #HACK Either write Carm ID into results file or make a more permanent name dict
    # dictionary connecting results_file["parameters"].attrs["pkl"] -> parameters.starname to CARMENES ID for serval results matching
    name_dict = {
        
        "GJ436"     : "J11421+267",
        "GJ1148"    : "J11417+427",
        "GJ3473"    : "J08023+033",
        "YZ Cet"    : "J01125-169",
        "GJ15A"     : "J00183+440",
        "GJ176"     : "J04429+189",
        "GJ536"     : "J14010-026",
        
        #"GJ581"     : "J15194-077", didn't download data
        
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
    
    for f in file_list:
        wobble_file = results_dir + f
        with h5py.File(wobble_file,'r') as fil:
            print("keys: ", fil.keys())
        
        parameters = rw.read_parameters_from_results(wobble_file)
        
        #TODO basic functionality:
        #1. write functions to make basic plot with Serval and Wobble data and auto fit, and rms
        #2. plot all rms into one plot (for same star) 
        
        #proof of concept: just plot time series auto fits for all stars individually
        #1. results file to vels
        #2. vels to kep_fit
        #3. plot
        ""
        #1. 
        #vels test on laptop
        carmenes_object_ID = name_dict[parameters.starname]
        bary_starname = simbad_dict[parameters.starname]
        #vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/" This does nothing right?
        os.makedirs(vels_dir, exist_ok = True)
        
        res = pr.Results_ws(wobble_file
                    , serval_dir
                    , carmenes_object_ID
                    , bary_starname
                    , load_bary = True
                    , archive = True)
        res.apply_corrections()
        vels_file = res.eval_vels(output_dir)

        
        #2.
        dataset_name = parameters.starname
        print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        kep_fit = vels_to_kep_fit(dataset_name, vels_file)
        
        #3.
        output_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] #TODO make into function?
        plot_time_series(kep_fit, output_file)
        
        #TODO Fix incorrect plot start for fit in time series, fixed with HACK for kep_model_x
        '''
        
        #HACK quck test onn laptop only
        
            #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/"#
    
    #base_dir = "/home/christian/Documents/wobble_aux/wobble_aux"
    #results_dir = base_dir + "/" + "../results/pipeline/pipeline_test_0/"
    #serval_dir = base_dir + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    #output_dir = base_dir + "/" + "../results/evaluate/{0}/".format(run_name)
    #os.makedirs(output_dir, exist_ok = True)
    
    
    results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    output_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/evaluate/{0}/".format(run_name)
    os.makedirs(output_dir, exist_ok = True)
    
    #vels_dir = output_dir
    #vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/"
    for f in file_list:
        wobble_file = results_dir + f
        vels_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] + ".vels"
        
        #2.
        dataset_name = os.path.splitext(f)[0]
        print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        kep_fit = vels_to_kep_fit(dataset_name, vels_file)
        
        #3.
        output_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] #TODO make into function?
        plot_time_series(kep_fit, output_file)
        
        
        
    
    
