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
import yaml

def vels_to_kep_fit(dataset_name, vels_file, n_planets_max = 2):
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
    ## Run GLS once 
    #rv.run_gls(fit)
    ##TEST rv.run_gls("jo",123)
    #rv.run_gls_o_c(fit)
    
    #add period limits
    rv.run_gls(fit,fend =2,fbeg=10000)
    rv.run_gls_o_c(fit,fend =2,fbeg=10000)
        
    # now lets find the planets in our data!
    fit.auto_fit_max_pl = n_planets_max
    fit.auto_fit_FAP_level = 0.001 # this corresponds to FAP = 0.1%. GLS power with FAP below that level we take as planet candidate.
    fit.auto_fit_allow_ecc = True # otherwise will look only for circular planets

    #fit = rv.find_planets(fit)
    fit = rv.find_planets_restricted(fit,fend=1.5)
    
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

def drop_points(kep_fit, vels_file, n = 4):
    #delete the n largest offsets     #NOTE Could use other criterion such as 5sigma
    o_c = kep_fit.fit_results.rv_model.o_c
    o_c_square = np.square(o_c)
    i_max_list = []
    for n in range(n):
        i_max = np.argmax(o_c_square)
        o_c_square[i_max] = 0 #eliminate counting it again
        i_max_list = np.append(i_max_list, i_max)
    print(i_max_list)
    vels = np.genfromtxt(vels_file)
    #delete entries fro vels
    vels  = np.delete(vels, i_max_list, axis = 0)
    new_vels_file = os.path.splitext(vels_file)[0] + "delete_{}".format(n) + ".vels"
    np.savetxt(new_vels_file, vels, fmt ='%.18f')
    return new_vels_file
    
    

def plot_time_series(kep_fit, output_file, table = True, width = 10, precision = 2):
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
    
    mpl.rc('text',usetex=True)
    #font = {'family' : 'normal','weight' : 'bold','size'   : 12,'serif':['Helvetica']}
    font = {'family' : 'normal','weight' : 'bold','size'   : 12,'serif':['Helvetica']}
    mpl.rc('font', **font)
    
    #print(mpl.rcParams.keys())
    mpl.rcParams['lines.markeredgecolor']='k'#black edges on symbols
    
    ##### time series format ######
    f = plt.figure(0, figsize=(8,6.5))
    plt.subplots_adjust(hspace=0.005)
    format_im = 'png' #'pdf' or png
    dpi = 300

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    #gs.update(  wspace=0.05)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    
    color = ['b', 'r', 'g', 'r']
    symbol = ['o', 'D', 'o', 'o']
    markersize = [5, 5, 6, 6] 
    alpha = [1, 1, 1, 1]
    
    #model_color = 'k'
    model_color = 'C7'
    model_lw = '1.0'
    jd_offset = 2450000 # for cleaner x-axis labels
    #### Get the time series (these below are self explanatory) ########  
    jd        = kep_fit.fit_results.rv_model.jd - jd_offset
    rvs       = kep_fit.fit_results.rv_model.rvs
    rv_err    = kep_fit.fit_results.rv_model.rv_err
    o_c       = kep_fit.fit_results.rv_model.o_c
    
    #print("jd: ",jd)
    print("type(jd): ", type(jd))
    
    data_set  = kep_fit.filelist.idset
    #print("data_set: ",data_set)

    # we can add the jitter
    add_jitter = True
    if add_jitter == True:
        rv_err = np.array([np.sqrt(rv_err[i]**2 + kep_fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])



    # Kep model time series #
    kep_model_x = kep_fit.fit_results.model_jd - jd_offset

    kep_model_y = kep_fit.fit_results.model

    

    ###################################################################

    offset_pre  = 250
    offset_post = 250

    zero_point_T = range((int(min(jd))-offset_pre),(int(max(jd))+offset_post),10)
    zero_point   = np.zeros(len(zero_point_T))


    ax1.plot(kep_model_x, kep_model_y,       '-', linewidth=model_lw, color=model_color)
    ax2.plot(zero_point_T,zero_point,'-', linewidth=model_lw, color=model_color)      

    
    #print("len(data_set): ",len(data_set)) #there seems to be an issue when just 
                                            #one dataset is included
    for i in range(len(data_set)):

            ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1
                         )
            ax2.errorbar(jd[i],o_c[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])],color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)



    
    ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 
    ax1.set_xlim(min(jd)-offset_pre,max(jd)+offset_post)
    
    #ax1.legend()

    ax2.set_xlabel(r'JD [d] - {}'.format(jd_offset),fontsize=16)
    ax2.set_ylabel(r'o$-$c  [m/s]',fontsize=16, rotation = 'vertical') 
    ax2.set_xlim(min(jd)-offset_pre,max(jd)+offset_post)

    ax2.locator_params(axis="x", nbins=9)
    
    plt.setp( ax2.get_yticklabels(), fontsize=15,weight='bold')
    plt.setp( ax2.get_xticklabels(), fontsize=15,weight='bold')
    
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
        
    
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) 
    
    ##TODO
    if table == True:
    #option 1 manual matplollib table
    #https://matplotlib.org/3.1.1/gallery/misc/table_demo.html#sphx-glr-gallery-misc-table-demo-py
        #adapt to number of planets
        npl = kep_fit.npl
        #HACK to prevent failure for no planets
        if npl == 0:
            npl = 1
        alpha_num_dict={1:'b', 2:'c', 3: 'd', 4: 'e'}
        print('npl: ', npl)
    
        
        #columns = ['Planet b', 'Planet c']
        columns = ['Planet {}'.format(alpha_num_dict[i]) for i in range(1,npl+1)]
        rows = ['K [m/s]', 'P [d]', 'r.m.s. [m/s]']
        
        cell_text_dict = {}
        for i in range(npl):
            try:
                parameter_ID =7*i
                cell_text_dict[(i,0)] = r'{0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f}'.format(kep_fit.params.planet_params[parameter_ID], max(np.abs(kep_fit.param_errors.planet_params_errors[parameter_ID])), width = width, precision = precision)
            except:
                cell_text_dict[(i,0)] = 'NA'
                
            try:
                parameter_ID =7*i +1
                cell_text_dict[(i,1)] = r'{0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f}'.format(kep_fit.params.planet_params[parameter_ID], max(np.abs(kep_fit.param_errors.planet_params_errors[parameter_ID])), width = width, precision = precision)
            except:
                cell_text_dict[(i,1)] = 'NA'
                
        try:
            rms = '{0:{width}.{precision}f}'.format(float(kep_fit.fit_results.rms), width = width, precision = precision)
        except:
            rms = 'NA'
        
        cell_text = [
        [cell_text_dict[(i,0)] for i in range(npl)], 
        [cell_text_dict[(i,1)] for i in range(npl)]
        ]
        rms_text = [rms]
        if npl !=1:
            rms_text.extend(['' for i in range(npl-1)])
        cell_text.append(rms_text)
        #[
        #[Kb, Kc],
        #[Pb, Pc],
        #[rms, '']
                    #]
        
        ax1.table(cellText=cell_text,
                      rowLabels=rows,
                      #rowColours=colors,
                      colLabels=columns,
                      loc='top')
                      #Adjust layout to make room for the table:
        plt.subplots_adjust(left=0.2,bottom=0.2) #TODO Check what this actually does
                      
                      
    #option 2 exostriker latex table probably possible but i need somehting dumb that works
        #table = latex_pl_param_tabular_string(kep_fit, width = 10, precision = 2)
        #dummy does nt work either without error? 
        #https://stackoverflow.com/questions/8524401/how-can-i-place-a-table-on-a-plot-in-matplotlib
        #table = latex_pl_param_tabular_string_dummy(kep_fit, width = 10, precision = 2)
        #plt.text(9,3.4,table,size=12)
                      


    # HACK plt.savefig('RV_plot_example_time_series.%s'%(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    plt.savefig(output_file + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla() 
    ax2.cla()
    
def latex_pl_param_tabular_string(obj, width = 10, precision = 2):
  
        #text =
        #'''       
    #\\begin{table}[ht]
    #% \\begin{adjustwidth}{-4.0cm}{} 
    #% \\resizebox{0.69\\textheight}{!}
    #% {\\begin{minipage}{1.1\\textwidth}
    
    #\centering   
    #\caption{{}}   
    #\label{table:}  
    
    text = r'''
    \\begin{tabular}{lrrrrrrrr}     % 2 columns 
    
    \hline\hline  \\noalign{\\vskip 0.7mm}      
    '''
    
     
    text = text + '''Parameter \hspace{0.0 mm}'''
    for i in range(obj.npl):     
        text = text + '''& Planet %s '''%chr(98+i)
    text = text + '''\\\\
    \hline \\noalign{\\vskip 0.7mm} 
    
    '''
    if obj.type_fit["RV"] == True:
        text = text + '''{0:{width}s}'''.format("$K$  [m\,s$^{-1}$]", width = 30)
        for i in range(obj.npl):     
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.params.planet_params[7*i], max(np.abs(obj.param_errors.planet_params_errors[7*i])), width = width, precision = precision)
        text = text + '''\\\\
        '''        
        
    text = text + '''{0:{width}s}'''.format("$P$  [day]", width = 30)
    for i in range(obj.npl):     
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.params.planet_params[7*i +1], max(np.abs(obj.param_errors.planet_params_errors[7*i +1])), width = width, precision = precision)
    text = text + '''\\\\
    '''  
    text = text + '''{0:{width}s}'''.format("$e$  ", width = 30)
    for i in range(obj.npl):     
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.params.planet_params[7*i +2], max(np.abs(obj.param_errors.planet_params_errors[7*i +2])), width = width, precision = precision)
    text = text + '''\\\\
    '''  
    text = text + '''{0:{width}s}'''.format("$\omega$  [deg]", width = 30)
    for i in range(obj.npl):     
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.params.planet_params[7*i +3], max(np.abs(obj.param_errors.planet_params_errors[7*i +3])), width = width, precision = precision)
    text = text + '''\\\\
    '''  
    text = text + '''{0:{width}s}'''.format("$M_{\\rm 0}$  [deg]", width = 30)
    for i in range(obj.npl):     
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.params.planet_params[7*i +4], max(np.abs(obj.param_errors.planet_params_errors[7*i +4])), width = width, precision = precision)
    text = text + '''\\\\
    '''         
    
    if obj.mod_dynamical == True:         
        text = text + '''{0:{width}s}'''.format("$i$  [deg]", width = 30)
        for i in range(obj.npl):     
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.params.planet_params[7*i +5], max(np.abs(obj.param_errors.planet_params_errors[7*i +5])), width = width, precision = precision)
        text = text + '''\\\\    
        '''      
        text = text + '''{0:{width}s}'''.format("$\Omega$  [deg]", width = 30)
        for i in range(obj.npl):     
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.params.planet_params[7*i +6], max(np.abs(obj.param_errors.planet_params_errors[7*i +6])), width = width, precision = precision)
        text = text + '''\\\\
        '''            
    
    if obj.type_fit["Transit"] == True:
        text = text + '''{0:{width}s}'''.format("$t_{\\rm 0}$  [day]", width = 30)
        for i in range(obj.npl):     
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.t0[i], max(abs(obj.t0_err[i])), width = width, precision = precision)
        text = text + '''\\\\
        '''            
        text = text + '''{0:{width}s}'''.format("Rad.  [$R_\oplus$]", width = 30)
        for i in range(obj.npl):     
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.pl_rad[i], max(abs(obj.pl_rad_err[i])), width = width, precision = precision)
        text = text + '''\\\\
        '''            
        text = text + '''{0:{width}s}'''.format("$a$  [$R_\odot$]", width = 30)
        for i in range(obj.npl):     
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.pl_a[i], max(abs(obj.pl_a_err[i])), width = width, precision = precision)
        text = text + '''\\\\
        '''   
        
    text = text + '''{0:{width}s}'''.format("$a$  [au]", width = 30)
    for i in range(obj.npl):     
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.fit_results.a[i], 0, width = width, precision = precision)
    text = text + '''\\\\
    '''      
    text = text + '''{0:{width}s}'''.format("$m \sin i$  [$M_{\\rm jup}$]", width = 30)
    for i in range(obj.npl):     
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.fit_results.mass[i], 0, width = width, precision = precision)
    text = text + '''\\\\
    '''          
    text = text + '''{0:{width}s}'''.format("$t_{\omega}$  [day]", width = 30)
    for i in range(obj.npl):     
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format((float(obj.epoch) - (np.radians(obj.params.planet_params[7*i + 4])/(2*np.pi))*obj.params.planet_params[7*i + 1] ), 0, width = width, precision = precision)
    text = text + '''\\\\ 
    '''          
    
    text = text + '''{0:{width}s}'''.format("RV lin. trend", width = 30)            
    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.params.linear_trend),float(max(np.abs(obj.param_errors.linear_trend_error))) , width = 30, precision = 6)
    text = text + '''\\\\
    '''   

    text = text + '''{0:{width}s}'''.format("RV quad. trend", width = 30)            
    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.rv_quadtr),float(max(np.abs(obj.rv_quadtr_err))) , width = 30, precision = 6)
    text = text + '''\\\\
    '''           
                    
    for i in range(obj.filelist.ndset):   
        text = text + '''{0:{width}s}'''.format("RV$_{\\rm off}$ %s"%(i+1), width = 30)            
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.params.offsets[i]), float(max(np.abs(obj.param_errors.offset_errors[i]))), width = width, precision = precision)
        text = text + '''\\\\
    '''   
    for i in range(obj.filelist.ndset):   
        text = text + '''{0:{width}s}'''.format("RV$_{\\rm jit}$ %s"%(i+1), width = 30)            
        text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.params.jitters[i]), float(max(np.abs(obj.param_errors.jitter_errors[i]))), width = width, precision = precision)
        text = text + '''\\\\
    '''   

    text = text + '''{0:{width}s}'''.format("$\chi^2$", width = 30)            
    text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.chi2), width = width, precision = precision)
    text = text + '''\\\\
    '''    
    text = text + '''{0:{width}s}'''.format("$\chi_{\\nu}^2$", width = 30)            
    text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.reduced_chi2), width = width, precision = precision)
    text = text + '''\\\\
    '''        
    text = text + '''{0:{width}s}'''.format("$r.m.s.$ [m\,s$^{-1}$]", width = 30)            
    text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.rms), width = width, precision = precision)
    text = text + '''\\\\
    '''            

    text = text + '''{0:{width}s}'''.format("$-\ln\mathcal{L}$", width = 30)            
    text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.loglik), width = width, precision = precision)
    text = text + '''\\\\
    '''        
    text = text + '''{0:{width}s}'''.format("N$_{\\rm RV}$ data", width = 30)            
    text = text + '''& {0:{width}.{precision}f} '''.format(len(obj.fit_results.jd), width = width, precision = 0)
    text = text + '''\\\\
    '''         
    
    text = text + '''{0:{width}s}'''.format("Epoch", width = 30)            
    text = text + '''& {0:{width}.{precision}f} '''.format(obj.epoch, width = width, precision = precision)
    text = text + '''\\\\
    '''           
    
    text = text + '''\\\\
    \hline \\noalign{\\vskip 0.7mm} 
    
    '''     
    
    text = text + '''        
    \end{tabular}  
    '''
    #% \end{minipage}}
    #% \end{adjustwidth}
    
    #%\\tablefoot{\small }
    
    #\end{table}
    #'''  
    

    return text


def latex_pl_param_tabular_string_dummy(obj, width = 10, precision = 2):
    
    test_string = r'''\begin{tabular}{ c | c | c | c } & col1 & col2 & col3 \\\hline row1 & 11 & 12 & 13 \\\hline row2 & 21 & 22 & 23 \\\hline  row3 & 31 & 32 & 33 \end{tabular}'''
    return test_string
  
#use as prototype for wobble equivalent:
def eval_serval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = {}, drop_pointsQ = False):
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    output_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/evaluate/{0}/".format(run_name)
    os.makedirs(output_dir, exist_ok = True)
    
    for index,f in enumerate(file_list):
        wobble_file = results_dir + f
        vels_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] + ".vels"

        #BEGIN EDIT
        #1 try to access results files on network drive from laptop: Works well, use as default
        parameters = rw.read_parameters_from_results(wobble_file)
        carmenes_object_ID = name_dict[parameters.starname]
        bary_starname = simbad_dict[parameters.starname]
        #vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/" #This does nothing right?
        #os.makedirs(vels_dir, exist_ok = True)

        res = pr.Results_ws(wobble_file
                    , serval_dir
                    , carmenes_object_ID
                    , bary_starname
                    , load_bary = True
                    , archive = True)
        res.apply_corrections()
        vels_file = res.eval_vels_serval(output_dir)#NOTE Only change fromm wobble implementation  so far
        #END EDIT

        #2.vels to kep_fit
        dataset_name = parameters.starname +" Serval_avcn"
        print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        if parameters.starname in n_planet_dict:
            n_planets_max = n_planet_dict[parameters.starname]
            print(parameters.starname, ": ", n_planets_max)
        else:
            n_planets_max = 2
        kep_fit = vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)
        
        #remove worst fitting poinnts and run again if option is True
        if drop_pointsQ:
            n = 4
            vels_file = drop_points(kep_fit, vels_file, n = n)
            dataset_name = dataset_name + "_dropped_{}".format(n)
            #fit again
            kep_fit = vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)

        #3.plot
        output_file = output_dir + dataset_name
        plot_time_series(kep_fit, output_file)

        #4 output latex_pl_param_table
        rv.latex_pl_param_table(kep_fit, width = 10, precision = 2, asymmetric = False, file_name= output_dir + dataset_name + '.tex', 
                                
                                #path= output_dir This apparently does nothing?
                                )
        #5 build rms list 
        starname = parameters.starname
        rms = float(kep_fit.fit_results.rms)
            
        if index == 0:
            rms_dict = {starname : rms}
        else:
            #rms_array = np.append(rms_array, [array_entry], axis = 0)
            rms_dict[starname] = rms 
            #print(rms_dict)
            
    rms_file = output_dir + "/rms_serval.yaml"
    with open(rms_file, 'w') as fp:
        yaml.dump(rms_dict, fp)

def eval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = {}, drop_pointsQ = False):
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    output_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/evaluate/{0}/".format(run_name)
    os.makedirs(output_dir, exist_ok = True)
    
    for index,f in enumerate(file_list):
        wobble_file = results_dir + f
        vels_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] + ".vels"

        #BEGIN EDIT
        #1 try to access results files on network drive from laptop: Works well, use as default
        parameters = rw.read_parameters_from_results(wobble_file)
        carmenes_object_ID = name_dict[parameters.starname]
        bary_starname = simbad_dict[parameters.starname]
        #vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/" #This does nothing right?
        #os.makedirs(vels_dir, exist_ok = True)

        res = pr.Results_ws(wobble_file
                    , serval_dir
                    , carmenes_object_ID
                    , bary_starname
                    , load_bary = True
                    , archive = True)
        res.apply_corrections()
        vels_file = res.eval_vels(output_dir)
        #END EDIT
        
        #2.vels to kep_fit
        dataset_name = os.path.splitext(f)[0]
        print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        if parameters.starname in n_planet_dict:
            n_planets_max = n_planet_dict[parameters.starname]
        else:
            n_planets_max = 2
        kep_fit = vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)

        #remove worst fitting poinnts and run again if option is True
        if drop_pointsQ:
            n = 4
            vels_file = drop_points(kep_fit, vels_file, n = n)
            dataset_name = dataset_name + "_dropped_{}".format(n)
            #fit again
            kep_fit = vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)    
    
        #3.plot
        output_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0]
        plot_time_series(kep_fit, output_file)
        
        #4 output latex_pl_param_table
        rv.latex_pl_param_table(kep_fit, width = 10, precision = 2, asymmetric = False, file_name= output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] + '.tex', 
                                
                                #path= output_dir This apparently does nothing?
                                )
        #5 build rms list 
        starname = parameters.starname
        rms = float(kep_fit.fit_results.rms)
        
        if index == 0:
            rms_dict = {starname : rms}
        else:
            rms_dict[starname] = rms 

            
    
    rms_file = output_dir + "/rms.yaml"
    with open(rms_file, 'w') as fp:
        yaml.dump(rms_dict, fp)
          
def plot_rv_comparison(kep_fit_list, output_file,
                       legend_labels = None,
                       phased = False,
                       format_im = 'png', #'pdf' or png
                       ip = 0, #declare index of primary data set
                       table = True, width = 10, precision = 2
                       ):
    if legend_labels == None:
        legend_labels = ['dataset {0}'.format(i) for i in range(len(kep_fit_list))]
        
# Edit of plot_time_series to include multiple fitted data sets
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
    
    mpl.rc('text',usetex=True)
    #font = {'family' : 'normal','weight' : 'bold','size'   : 12,'serif':['Helvetica']}
    font = {'family' : 'normal','weight' : 'bold','size'   : 12,'serif':['Helvetica']}
    mpl.rc('font', **font)
    
    #print(mpl.rcParams.keys())
    mpl.rcParams['lines.markeredgecolor']='k'#black edges on symbols
    
    ##### time series format ######
    f = plt.figure(0, figsize=(8,6.5))
    plt.subplots_adjust(hspace=0.005)
    format_im = format_im #'pdf' or png
    dpi = 300

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    #gs.update(  wspace=0.05)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    
    color = ['C0', 'C1', 'C2', 'C3']
    symbol = ['o', 'D', 'v', 's']
    markersize = [3, 3, 3, 3] 
    alpha = [1, 1, 1, 1]
    
    #model_color = 'k'
    model_color = 'C7'
    model_lw = '1.0'
###########################
    if phased == False:
    #Data reading
        jd_offset = 2450000 # for cleaner x-axis labels
        #### Get the time series (these below are self explanatory) ########  
        jd = np.array([kep_fit.fit_results.rv_model.jd - jd_offset 
                    for i,kep_fit in enumerate(kep_fit_list)])
        rvs = np.array([kep_fit.fit_results.rv_model.rvs 
                        for i,kep_fit in enumerate(kep_fit_list)])
        rv_err = np.array([kep_fit.fit_results.rv_model.rv_err 
                        for i,kep_fit in enumerate(kep_fit_list)])
        o_c = np.array([kep_fit.fit_results.rv_model.o_c 
                        for i,kep_fit in enumerate(kep_fit_list)])
        
        
        #print(jd)
        #data_set  = kep_fit.filelist.idset

        # we can add the jitter
        add_jitter = True
        if add_jitter == True:
            #if len(kep_fit_list) == 1:
                ##Fallback if only one kep_fit was supplied #TODO is this required?
                #kep_fit = kep_fit_list[0]
                #data_set  = kep_fit.filelist.idset
                #rv_err[0] = np.array([np.sqrt(rv_err[0][i]**2 + kep_fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])
            #else:
            for j in range(len(kep_fit_list)):
                kep_fit = kep_fit_list[j]
                data_set  = kep_fit.filelist.idset
                rv_err[j] = np.array([np.sqrt(rv_err[j][i]**2 + kep_fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])
                



        # Kep model time series #
        ip = 0 #primary fit index #TODO make input
        kep_fit_primary = kep_fit_list[ip]
        kep_model_x = kep_fit_primary.fit_results.model_jd - jd_offset

        kep_model_y = kep_fit_primary.fit_results.model

        

        ###################################################################
        
        #TODO make scaling, make fit plot run to edges?
        offset_pre  = 100
        offset_post = 100

        zero_point_T = range((int(min(jd[ip]))-offset_pre),(int(max(jd[ip]))+offset_post),10)
        zero_point   = np.zeros(len(zero_point_T))


        ax1.plot(kep_model_x, kep_model_y,       '-', linewidth=model_lw, color=model_color, label = "{0} fit".format(legend_labels[ip]))
        ax2.plot(zero_point_T,zero_point,'-', linewidth=model_lw, color=model_color)      

        
        #print("len(data_set): ",len(data_set)) #there seems to be an issue when just 
                                                #one dataset is included
        for i in range(len(kep_fit_list)):

                ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[i], fmt=symbol[i], linestyle='None', markersize = markersize[i], color=color[i], capsize = 0, elinewidth=1,mew=0.1,
                            label = legend_labels[i] + ', rms = ' + '{0:{width}.{precision}f}'.format(float(kep_fit_list[i].fit_results.rms), width = width, precision = precision)
                            )
                ax2.errorbar(jd[i],o_c[i], yerr=rv_err[i], alpha=alpha[i], fmt=symbol[i], linestyle='None', markersize = markersize[i],color=color[i], capsize = 0, elinewidth=1,mew=0.1)
                
        ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 
        ax1.set_xlim(min(jd[ip])-offset_pre,max(jd[ip])+offset_post)
        
        ax1.legend(loc = 'upper left')
        #add padding for legend
        ymin, ymax = ax1.get_ylim()
        ax1.set_ylim(ymin, ymax + (ymax-ymin)*0.15)#TODO make this scale with legend size

        ax2.set_xlabel(r'JD [d] - {}'.format(jd_offset),fontsize=16)
        ax2.set_ylabel(r'o$-$c  [m/s]',fontsize=16, rotation = 'vertical') 
        ax2.set_xlim(min(jd[ip])-offset_pre,max(jd[ip])+offset_post)

        ax2.locator_params(axis="x", nbins=9)
        
        plt.setp( ax2.get_yticklabels(), fontsize=15,weight='bold')
        plt.setp( ax2.get_xticklabels(), fontsize=15,weight='bold')
        
        # Fine-tune figure; make subplots close to each other and hide x ticks for
        # all but bottom plot.
            
        
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        
        ##TODO Input Table if desired
        plt.grid(True)
        plt.savefig(output_file + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
        ax1.cla() 
        ax2.cla()
        
    if phased == True:
        for j in range(kep_fit_list[ip].npl):

            f = plt.figure(1, figsize=(9,6))
        
            ax1 = plt.subplot(111)
            
            for i in range(len(kep_fit_list)):
                planet = j+1
                kep_fit = kep_fit_list[i]
                data, model = rv.phase_RV_planet_signal(kep_fit,planet)

                jd = data[0]
                rvs = data[1]
                rv_err = data[2]
                data_set = data[3]
                
                if i == ip:
                    ax1.plot(model[0],model[1], linestyle =  '-', linewidth=model_lw, color=model_color, label = "{0} fit".format(legend_labels[ip]))

                # we can add the jitter
                add_jitter = True
                if add_jitter == True:
                    rv_err = np.array([np.sqrt(rv_err[k]**2 + kep_fit.params.jitters[kk]**2)  for k,kk in enumerate(data_set)])

                
                ax1.errorbar(jd,rvs, yerr=rv_err, alpha=alpha[i], fmt=symbol[i], linestyle='None', markersize = markersize[i], color=color[i], capsize = 0, elinewidth=1,mew=0.1,
                label = legend_labels[i] + ', rms = ' + '{0:{width}.{precision}f}'.format(float(kep_fit_list[i].fit_results.rms), width = width, precision = precision)
                                )
                
            ax1.legend(loc = 'upper left')
            #add padding for legend
            ymin, ymax = ax1.get_ylim()
            ax1.set_ylim(ymin, ymax + (ymax-ymin)*0.15)#TODO make this scale with legend size
            
            ax1.set_xlabel(r'phase [d]',fontsize=16)
            ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 

            plt.grid(True)
            plt.savefig(output_file + '_planet_%s.%s'%(j+1,format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
            ax1.cla() 


def output_fit_parameters(kep_fit):
    #orbital_parameters = [K, P, e, omega, M0] NOTE T0 removed: requires change to function in compare_ws
    #mult [[planet_b],[planet_c],...]
    #also output rms for comparison
    #T_epoch (T_0) = obj.epoch
    obj = kep_fit
    orbital_parameters = np.zeros((obj.npl, 5))
    
    T_epoch = obj.epoch
    
    if obj.type_fit["RV"] == True:
        for i in range(obj.npl):
            orbital_parameters[i] = [obj.params.planet_params[7*i + j] for j in range(5)]
            
    return [T_epoch, orbital_parameters]
    

''' THis Dictionary does not always follow the naming conventions I used for "starname"
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
    
    "GJ581"     : "J15194-077",
    
    "GJ3512"    : "J08413+594", #issues due to low min_snr, drops all orders at snr 60
    "Wolf294"   : "J06548+332",
    "GJ876"     : "J22532-142",
    "Teegarden" : "J02530+168", # issues with differing min_snr recheck
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

n_planet_dict = {
    
    "GJ436"     : 1, #should be 1 ,2 just foor test
    "GJ1148"    : 2,
    "GJ3473"    : 2, #default
    "YZ Cet"    : 2, #default
    "GJ15A"     : 2, #default
    "GJ176"     : 2, #default
    "GJ536"     : 2, #default
    
    "GJ581"     : 2, #default
    
    "GJ3512"    : 2, #1 condirmed , 1 as long term trend 
    "Wolf294"   : 2, #default
    "GJ876"     : 4, #Carmenes Paper
    "Teegarden" : 2, #default
    "Barnard"   : 2, #default
    
    }

if __name__ == "__main__": 
    
    
    #run_name = "test_1"
    #file_list = ["results_GJ436_Kstar0_Kt3_laptop_example_0.hdf5", "results_GJ3473_Kstar0_Kt3_laptop_example_1.hdf5" ] #laptop sample file
    
    #run_name = "test_2"
    #file_list = ["results_GJ1148_Kstar0_Kt3_eval_example_2.hdf5","results_GJ436_Kstar0_Kt3_laptop_example_0.hdf5", "results_GJ3473_Kstar0_Kt3_laptop_example_1.hdf5"]
    

    #run_name ="baseline_0.3_complete_function_test"
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    #file_list = ["results_GJ436_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ1148_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ3473_Kstar0_Kt3_baseline_0.hdf5",
        #"results_YZ Cet_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ15A_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ176_Kstar0_Kt3_baseline_0.hdf5", "results_GJ536_Kstar0_Kt3_baseline_0.hdf5", "results_GJ3512_Kstar0_Kt3_baseline_0.hdf5", "results_Wolf294_Kstar0_Kt3_baseline_0.hdf5", "results_GJ876_Kstar0_Kt3_baseline_0.hdf5" , "results_Teegarden_Kstar0_Kt3_baseline_0.hdf5", "results_Barnard_Kstar0_Kt3_baseline_0.hdf5"
                 #]
    

    ##BEGIN n test
    ##loop over all iteration test pipeline results
    #for n in [40,80,120,160,200]:
    ##GJ876 with  higher n
    ##for n in [800,900,1000]:
       ## name_dict = {"GJ876"     : "J22532-142"}
        #name_dict = {"GJ436"     : "J11421+267"}
        #run_name = "n_{}".format(n)
        #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_n{}/".format(n)
        #file_list = ["results_{0}_Kstar0_Kt3_n_{1}.hdf5".format(starname,n) for starname in name_dict.keys()]
        
        #eval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = n_planet_dict)
        #eval_serval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = n_planet_dict)
    ##END n test
    
    ##BEGIN n test
    ##loop over all iteration test pipeline results
    ##for n in [40,80,120,160,200]:
    ##GJ876 with  higher n
    #for n in [#40,80,120,160,200,
              #250,300,350,400,450,500,600
              ##,700,800,900,1000
              #]:
        #name_dict = {"GJ876"     : "J22532-142"}
        ##name_dict = {"GJ3512"    : "J08413+594"}        
        #run_name = "n_{}".format(n)
        #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_n{}/".format(n)
        #file_list = ["results_{0}_Kstar0_Kt3_n_{1}.hdf5".format(starname,n) for starname in name_dict.keys()]
        
        #eval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = n_planet_dict)
        #eval_serval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = n_planet_dict)
    ###END n test
    
    ##BEGIN
    #run_name ="baseline_0.5_n_planets"
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    #file_list = ["results_GJ436_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ1148_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ3473_Kstar0_Kt3_baseline_0.hdf5",
        #"results_YZ Cet_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ15A_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ176_Kstar0_Kt3_baseline_0.hdf5", "results_GJ536_Kstar0_Kt3_baseline_0.hdf5", "results_GJ3512_Kstar0_Kt3_baseline_0.hdf5", "results_Wolf294_Kstar0_Kt3_baseline_0.hdf5", "results_GJ876_Kstar0_Kt3_baseline_0.hdf5" , "results_Teegarden_Kstar0_Kt3_baseline_0.hdf5", "results_Barnard_Kstar0_Kt3_baseline_0.hdf5"
                 #]
    #eval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = n_planet_dict)
    #eval_serval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = n_planet_dict)
    ##END
    
    #BEGIN
    run_name = "drop_points_test_TeX"
    results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    file_list = ["results_GJ436_Kstar0_Kt3_baseline_0.hdf5",
        "results_GJ1148_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ3473_Kstar0_Kt3_baseline_0.hdf5",
        #"results_YZ Cet_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ15A_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ176_Kstar0_Kt3_baseline_0.hdf5", "results_GJ536_Kstar0_Kt3_baseline_0.hdf5", "results_GJ3512_Kstar0_Kt3_baseline_0.hdf5", "results_Wolf294_Kstar0_Kt3_baseline_0.hdf5", "results_GJ876_Kstar0_Kt3_baseline_0.hdf5" , "results_Teegarden_Kstar0_Kt3_baseline_0.hdf5", "results_Barnard_Kstar0_Kt3_baseline_0.hdf5"
                 ]

    eval_serval_complete(run_name, file_list, name_dict, simbad_dict, results_dir, n_planet_dict = n_planet_dict, drop_pointsQ = True)
    #END
    
    
    ##BEGIN
    #run_name ="baseline_0.4_fend"
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    #file_list = ["results_GJ436_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ1148_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ3473_Kstar0_Kt3_baseline_0.hdf5",
        #"results_YZ Cet_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ15A_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ176_Kstar0_Kt3_baseline_0.hdf5", "results_GJ536_Kstar0_Kt3_baseline_0.hdf5", "results_GJ3512_Kstar0_Kt3_baseline_0.hdf5", "results_Wolf294_Kstar0_Kt3_baseline_0.hdf5", "results_GJ876_Kstar0_Kt3_baseline_0.hdf5" , "results_Teegarden_Kstar0_Kt3_baseline_0.hdf5", "results_Barnard_Kstar0_Kt3_baseline_0.hdf5"
                 #]
    ##eval_complete(run_name, file_list, name_dict, simbad_dict, results_dir)
    #eval_serval_complete(run_name, file_list, name_dict, simbad_dict, results_dir)
    ##END
    
    
    ############
    ##TEST drop_points
    #run_name  = "drop_points_test"
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    #serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    #output_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/evaluate/{0}/".format(run_name)
    #os.makedirs(output_dir, exist_ok = True)
    
    #file_list = ["results_GJ436_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ1148_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ3473_Kstar0_Kt3_baseline_0.hdf5",
        #"results_YZ Cet_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ15A_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ176_Kstar0_Kt3_baseline_0.hdf5", "results_GJ536_Kstar0_Kt3_baseline_0.hdf5", "results_GJ3512_Kstar0_Kt3_baseline_0.hdf5", "results_Wolf294_Kstar0_Kt3_baseline_0.hdf5", "results_GJ876_Kstar0_Kt3_baseline_0.hdf5" , "results_Teegarden_Kstar0_Kt3_baseline_0.hdf5", "results_Barnard_Kstar0_Kt3_baseline_0.hdf5"
                 #]
    
    #for index,f in enumerate(file_list):
        #wobble_file = results_dir + f
        #vels_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] + ".vels"

        ##BEGIN EDIT
        ##1 try to access results files on network drive from laptop: Works well, use as default
        #parameters = rw.read_parameters_from_results(wobble_file)
        #carmenes_object_ID = name_dict[parameters.starname]
        #bary_starname = simbad_dict[parameters.starname]
        ##vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/" #This does nothing right?
        ##os.makedirs(vels_dir, exist_ok = True)

        #res = pr.Results_ws(wobble_file
                    #, serval_dir
                    #, carmenes_object_ID
                    #, bary_starname
                    #, load_bary = True
                    #, archive = True)
        #res.apply_corrections()
        #vels_file = res.eval_vels(output_dir)
        ##END EDIT

        ##2.vels to kep_fit
        #dataset_name = os.path.splitext(f)[0]
        #print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        #kep_fit = vels_to_kep_fit(dataset_name, vels_file)
    
        ##3.plot
        #output_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0]
        #plot_time_series(kep_fit, output_file)
        
        #n = 4
        #vels_file = drop_points(kep_fit, vels_file, n = n)
        ##2.vels to kep_fit
        #dataset_name = os.path.splitext(f)[0] + "_dropped_{}".format(n)
        #print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        #kep_fit = vels_to_kep_fit(dataset_name, vels_file)
        
        ##3.plot
        #output_file = output_dir + dataset_name
        #plot_time_series(kep_fit, output_file)
    
    
    ########### 
    
    #basic functionality:
    #1. write functions to make basic plot with Serval and Wobble data and auto fit, and rms
    #2. plot all rms into one plot (for same star) 
    
    #proof of concept: just plot time series auto fits for all stars individually
    #1. results file to vels
    #2. vels to kep_fit
    #3. plot
    #4 output latex_pl_param_table
    #5 build rms list 
    
        
    
    #Serval eval only really needs to be  run  onnce #TODO write load(from baseline folder))
    #eval_serval_complete(run_name, file_list, name_dict, simbad_dict, results_dir)
    
    #eval_complete(run_name, file_list, name_dict, simbad_dict, results_dir)


    '''
    #vels_dir = output_dir
    #vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/"
    for index,f in enumerate(file_list):
        wobble_file = results_dir + f
        vels_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] + ".vels"
        
        #BEGIN EDIT
        #1 try to access results files on network drive from laptop: Works well, use as default
        parameters = rw.read_parameters_from_results(wobble_file)
        carmenes_object_ID = name_dict[parameters.starname]
        bary_starname = simbad_dict[parameters.starname]
        #vels_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/vels_dir/" #This does nothing right?
        #os.makedirs(vels_dir, exist_ok = True)
        
        res = pr.Results_ws(wobble_file
                    , serval_dir
                    , carmenes_object_ID
                    , bary_starname
                    , load_bary = True
                    , archive = True)
        res.apply_corrections()
        vels_file = res.eval_vels(output_dir)
        #END EDIT
        
        #2.vels to kep_fit
        dataset_name = os.path.splitext(f)[0]
        print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        kep_fit = vels_to_kep_fit(dataset_name, vels_file)
        
        #3.plot
        output_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] #TODO make into function?
        plot_time_series(kep_fit, output_file)
        
        #4 output latex_pl_param_table
        rv.latex_pl_param_table(kep_fit, width = 10, precision = 2, asymmetric = False, file_name= output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0] + '.tex', 
                                
                                #path= output_dir This apparently does nothing?
                                )
        
        
        #5 overview of rms and rms improvement 
        #a) write function that extracts rmms and writes it into a txt file 
        #b) thusly compile list of rms for both baseline
        #c) compile  list for each change
        #d) compare rms with baseline

        
        starname = parameters.starname
        rms = float(kep_fit.fit_results.rms)
        
        #TODO write these to a txt file
        #load existing file
        #if os.path.isfile(rms_file):
            #rms_array = np.genfromtxt(rms_file)
            #rms_array = np.append(rms_array, array_entry, axis = 0)
        #else:
            #rms_array = [array_entry] #this mmight  be a problem
        if index == 0:
            rms_dict = {starname : rms}
        else:
            #rms_array = np.append(rms_array, [array_entry], axis = 0)
            rms_dict[starname] = rms 
            print(rms_dict)
            
    
    rms_file = output_dir + "/rms.yaml"
    with open(rms_file, 'w') as fp:
        yaml.dump(rms_dict, fp)
    #np.savetxt(rms_file, rms_array , fmt ='%s , %s'
               #) # This would sav eth e number as a string which would require a backconversion when  reading
    
    with open(rms_file, 'r') as fp:
        rms_dict = yaml.load(fp, Loader=yaml.FullLoader)
    #rms_array = read_rms_file(output_dir + "/rms.txt")
    print(rms_dict)
    '''
