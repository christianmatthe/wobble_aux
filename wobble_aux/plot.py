
import os
import wobble
import numpy as np
import matplotlib.pyplot as plt

import run_wobble as rw
import process_results as pr
import cws_refactor as cws

#import exostriker functions from  evaluate results
import evaluate_results as er
import matplotlib as mpl

################
name_dict = er.name_dict #connects starname to carmenes ID
simbad_dict = er.simbad_dict #connects starname to simbad iidentifiers that work with barycorr
n_planet_dict = er.n_planet_dict #dictionary of max number of planets to be fit
###########


def op(order_index, arm = "vis"):
    #calculates physical CARMENES interference order ("order physical -> op" from order index
    if arm == "vis":
        order_max = 118
    elif arm == "nir":
        order_max = 63
    op = order_max - order_index
    return op

def plot_wobble_order(orders, epochs, results_file, plot_dir, show_mask=False):
    mpl.rcParams.update(mpl.rcParamsDefault) #reset to default plotting parameters
    
    p = rw.read_parameters_from_results(results_file)
    p.plot_continuum = False #prevent loaded parameters from replotting continuum
    #copy data file name convetion from run_wobble
    data_file = p.data_file
    
    data = wobble.Data(data_file, orders=orders, min_flux=10**-5, min_snr=0, parameters = p)
    
    results = wobble.Results(filename = results_file)
    
    epochs_results = results.epochs
    
    lowest_optimized_order = 11
    orders_to_plot = orders
    for o in orders_to_plot:
        #r = o - lowest_optimized_order
        r_r = np.where(np.array(results.orders) == o)[0][0]
        r_d = np.where(np.array(data.orders) == o)[0][0]
        for e_ind ,e in enumerate(epochs):
            if e not in epochs_results:
                print("epoch {} not in epochs_results. Skipping.".format(e))
                continue #Skip epochs that may have been dropped
            #print(e)
            #print(np.array(epochs_results))
            ep = np.where(np.array(epochs_results) == e)[0][0]
            #e for data ep for results. assumes data is not epoch cut
            fig, (ax, ax2) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[4, 1]}, figsize=(12,5))
            xs = np.exp(data.xs[r_d][e])
            ax.scatter(xs, np.exp(data.ys[r_d][e]), marker=".", alpha=0.5, c='k', label='data', s=40)
            mask = data.ivars[r_d][e] <= 1.e-8
            ax.scatter(xs[mask], np.exp(data.ys[r_d][e][mask]), marker=".", alpha=1., c='white', s=20)
            ax.plot(xs, np.exp(results.star_ys_predicted[r_r][ep]), c='r', alpha=0.8, label='Stellar model')
            ax.plot(xs, np.exp(results.tellurics_ys_predicted[r_r][ep]), c='b', alpha=0.8, label='Telluric model')
            ylims = ax.get_ylim()
            xlims = ax.get_xlim()
            #also plot telluric mask
            mask_str = ''
            if show_mask == True:
                mask_str = "_mask"
                file_dir = os.path.dirname(os.path.abspath(__file__))
                telluric_mask_file = file_dir + "/"+"carmenes_aux_files/" + "telluric_mask_carm_short.dat"
                telluric_mask = np.genfromtxt(telluric_mask_file)
                #ax.plot(telluric_mask[:,0], telluric_mask[:,1], c = 'g', alpha = 0.8, label = 'SERVAL telluric mask')
                ax.fill_between(telluric_mask[:,0], 1.3*telluric_mask[:,1], color = 'g', alpha = 0.6, label = 'SERVAL telluric mask')
                ax2.fill_between(telluric_mask[:,0], -0.08 + 0.16*telluric_mask[:,1],y2 = -0.08, color = 'g', alpha = 0.6, label = 'SERVAL telluric mask')
                
            ax.legend(loc = "lower left")
            
            ax2.scatter(xs, np.exp(data.ys[r_d][e]) - np.exp(results.star_ys_predicted[r_r][ep]
                                                        + results.tellurics_ys_predicted[r_r][ep]), 
                        marker=".", alpha=0.5, c='k', label='data', s=40)
            ax2.scatter(xs[mask], np.exp(data.ys[r_d][e][mask]) - np.exp(results.star_ys_predicted[r_r][ep]
                                                        + results.tellurics_ys_predicted[r_r][ep])[mask], 
                        marker=".", alpha=1., c='white', s=20)
            
            ax.set_ylabel('Normalized Flux')
            ax2.set_ylabel('Residuals')
            ax2.set_xlabel(r'Wavelength ($\AA$)')
            
            ax.set_ylim([0.0,1.3])
            #ax.set_xlim([xs[0],xs[-1]])
            ax.set_xlim(xlims)
            ax2.set_xlim(xlims)
            ax2.set_ylim([-0.08,0.08])
            ax.set_xticklabels([])
            fig.tight_layout()
            fig.subplots_adjust(hspace=0.05)
            
            #plt.title("JD {0}".format(data.dates[e])) #HACK
            
            plt.savefig(plot_dir + str(p.starname) + '_results_synth_op{0}_o{1}_e{2}{3}.pdf'.format(op(o), o, e, mask_str))
            plt.savefig(plot_dir + str(p.starname) + '_results_synth_op{0}_o{1}_e{2}{3}.png'.format(op(o), o, e, mask_str))
            plt.close(fig)
            
def plot_wobble_order_zoom(orders, epochs, results_file, plot_dir, show_mask=False, xlims=[], ylims= [0.0, 1.3] ):
    mpl.rcParams.update(mpl.rcParamsDefault) #reset to default parameters
    
    p = rw.read_parameters_from_results(results_file)
    p.plot_continuum = False #prevent loaded parameters from replotting continuum
    #copy data file name convetion from run_wobble
    data_file = p.data_file
    
    data = wobble.Data(data_file, orders=orders, min_flux=10**-5, min_snr=0, parameters = p)
    
    results = wobble.Results(filename = results_file)
    
    epochs_results = results.epochs
    
    lowest_optimized_order = 11
    orders_to_plot = orders
    for o in orders_to_plot:
        #r = o - lowest_optimized_order
        r_r = np.where(np.array(results.orders) == o)[0][0]
        r_d = np.where(np.array(data.orders) == o)[0][0]
        for e_ind ,e in enumerate(epochs):
            #print(e)
            #print(np.array(epochs_results))
            ep = np.where(np.array(epochs_results) == e)[0][0]
            #e for data ep for results. assumes data is not epoch cut
            fig, (ax, ax2) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[4, 1]}, figsize=(12,5))
            xs = np.exp(data.xs[r_d][e])
            ax.scatter(xs, np.exp(data.ys[r_d][e]), marker=".", alpha=0.5, c='k', label='data', s=40)
            mask = data.ivars[r_d][e] <= 1.e-8
            ax.scatter(xs[mask], np.exp(data.ys[r_d][e][mask]), marker=".", alpha=1., c='white', s=20)
            ax.plot(xs, np.exp(results.star_ys_predicted[r_r][ep]), c='r', alpha=0.8, label='Stellar model')
            ax.plot(xs, np.exp(results.tellurics_ys_predicted[r_r][ep]), c='b', alpha=0.8, label='Telluric model')
            ylims = ylims
            xlims = xlims
            #also plot telluric mask
            mask_str = ''
            if show_mask == True:
                mask_str = "_mask"
                file_dir = os.path.dirname(os.path.abspath(__file__))
                telluric_mask_file = file_dir + "/"+"carmenes_aux_files/" + "telluric_mask_carm_short.dat"
                telluric_mask = np.genfromtxt(telluric_mask_file)
                #ax.plot(telluric_mask[:,0], telluric_mask[:,1], c = 'g', alpha = 0.8, label = 'SERVAL telluric mask')
                ax.fill_between(telluric_mask[:,0], 1.3*telluric_mask[:,1], color = 'g', alpha = 0.6, label = 'SERVAL telluric mask')
                ax2.fill_between(telluric_mask[:,0], -0.08 + 0.16*telluric_mask[:,1],y2 = -0.08, color = 'g', alpha = 0.6, label = 'SERVAL telluric mask')
                
            ax.legend(loc = "lower left")
            
            ax2.scatter(xs, np.exp(data.ys[r_d][e]) - np.exp(results.star_ys_predicted[r_r][ep]
                                                        + results.tellurics_ys_predicted[r_r][ep]), 
                        marker=".", alpha=0.5, c='k', label='data', s=40)
            ax2.scatter(xs[mask], np.exp(data.ys[r_d][e][mask]) - np.exp(results.star_ys_predicted[r_r][ep]
                                                        + results.tellurics_ys_predicted[r_r][ep])[mask], 
                        marker=".", alpha=1., c='white', s=20)
            
            ax.set_ylabel('Normalized Flux')
            ax2.set_ylabel('Residuals')
            ax2.set_xlabel(r'Wavelength ($\AA$)')
            
            ax.set_ylim(ylims)
            #ax.set_xlim([xs[0],xs[-1]])
            ax.set_xlim(xlims)
            ax2.set_xlim(xlims)
            ax2.set_ylim([-0.08,0.08])
            ax.set_xticklabels([])
            fig.tight_layout()
            fig.subplots_adjust(hspace=0.05)
            
            #plt.title("JD {0}".format(data.dates[e])) #HACK
            
            plt.savefig(plot_dir + str(p.starname) + '_zoom_op{0}_o{1}_e{2}{3}.pdf'.format(op(o), o, e, mask_str))
            plt.savefig(plot_dir + str(p.starname) + '_zoom_op{0}_o{1}_e{2}{3}.png'.format(op(o), o, e, mask_str))
            plt.close(fig)
    
if __name__ == "__main__":
    #NOTE #####################
    output_dir_master = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/fin_plots/MT_git/figures/"
    #NOTE #####################
    
    
    #file_dir = os.path.dirname(__file__)
    #telluric_mask_file = file_dir + "/"+"carmenes_aux_files/" + "telluric_mask_carm_short.dat"
    
    """
    #dummy example res file in lx39 format
    results_file = "/home/christian/lx39_data/cmatthe/wobble_aux/results/pipeline/pipeline_min_snr60/results_GJ436_Kstar0_Kt3_min_snr60.hdf5"
    #results_file = "/data/cmatthe/wobble_aux/results/pipeline/pipeline_min_snr60/results_GJ436_Kstar0_Kt3_min_snr60.hdf5"
    
    run_name = "basic_function"
    plot_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/fin_plots/{0}/".format(run_name)
    os.makedirs(plot_dir, exist_ok = True)
    
    #orders = [i for i in range(11,52)]
    orders = [43]
    epochs = [10,40,70]
    
    plot_wobble_order(orders, epochs, results_file, plot_dir, show_mask = False)
    """

#Basic Wobble functionn presentation plot

    #results_file = "/home/christian/lx39_data/cmatthe/wobble_aux/results/pipeline_2/pipeline_niter_dict/results_GJ436_Kstar0_Kt3_n_160.hdf5"
    ##results_file = "/data/cmatthe/wobble_aux/results/pipeline_2/pipeline_niter_dict/results_GJ436_Kstar0_Kt3_n_160.hdf5"
    
    #run_name = "basic_function"
    #plot_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/fin_plots/{0}/".format(run_name)
    #os.makedirs(plot_dir, exist_ok = True)
    
    ##orders = [i for i in range(11,52)]
    #orders = [43]
    #epochs = [10,40,70]
    ##plot_wobble_order(orders, epochs, results_file, plot_dir, show_mask = False)
    #plot_wobble_order(orders, epochs, results_file, plot_dir, show_mask = True)
    
    #plot_wobble_order_zoom(orders, epochs, results_file, plot_dir, show_mask=True, xlims = [8090, 8118], ylims= [0.7, 1.1] )
 
# Check Teegarden emissionn in order 45 and 49 
    #results_file = "/data/cmatthe/wobble_aux/results/pipeline_2/pipeline_niter_dict/results_Teegarden_Kstar0_Kt3_n_160.hdf5"
    
    ##results_file = "/home/christian/lx39_data/cmatthe/wobble_aux/results/pipeline_2/pipeline_niter_dict_continuum_test_[0.5,1]/results_Teegarden_Kstar0_Kt3_n_160.hdf5"
    ##results_file = "/data/cmatthe/wobble_aux/results/pipeline_2/pipeline_niter_dict_continuum_test_[0.5,1]/results_Teegarden_Kstar0_Kt3_n_160.hdf5"
    
    #run_name = "Teegarden_emission"
    #plot_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/fin_plots/{0}/".format(run_name)
    #os.makedirs(plot_dir, exist_ok = True)
    
    ##orders = [i for i in range(11,52)]
    #orders = [45,49]
    #epochs = [i for i in range(72)]
    
    #plot_wobble_order(orders, epochs, results_file, plot_dir, show_mask = False)

#BEGIN GJ1148 o49 e31 bad SNR/Continuum example
    results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline_2/pipeline_pp_reruns/"
    file_list = [#"Baseline+SNRresults_GJ436_Kstar0_Kt3_min_snr_5.hdf5",
                 "Baseline+SNRresults_GJ1148_Kstar0_Kt3_min_snr_5.hdf5"
                 ]
    
    run_name = "bad_SNR/Cont"
    output_dir = output_dir_master + "{0}/".format(run_name)
    os.makedirs(output_dir, exist_ok = True)
    
    for f in file_list:
        results_file = results_dir + f
        
        #orders = [i for i in range(11,52)]
        orders = [49]
        #epochs = [10,40,70]
        #epochs = [i for i in range(25,35)]
        epochs = [26,31]
        
        plot_wobble_order(orders, epochs, results_file, output_dir, show_mask = True)
#END
########################################################################################
##BEGIN Initial RV Results
    #run_name = "Initial_RV"
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline_2/pipeline_pp_reruns/"
    #file_list = ["Baseline+SNRresults_GJ436_Kstar0_Kt3_min_snr_5.hdf5",
                 #"Baseline+SNRresults_GJ1148_Kstar0_Kt3_min_snr_5.hdf5"
                 #]

    
    ##serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    #output_dir = output_dir_master + "{0}/".format(run_name)
    #os.makedirs(output_dir, exist_ok = True)
    
    #for index,f in enumerate(file_list):
        #wobble_file = results_dir + f
        #dataset_name = os.path.splitext(os.path.split(wobble_file)[1])[0]
        #output_folder = output_dir + dataset_name + "/"
        #os.makedirs(output_folder, exist_ok = True)
        
        ##make res
        #res = cws.load_results(wobble_file, correct_w_for_drift = True, correct_w_for_NZP = True)
        
        ##make vels
        #vels_file = res.eval_vels_serval(output_dir)
        
        ##make kep_fit
        #parameters = rw.read_parameters_from_results(wobble_file)
        #n_planets_max = n_planet_dict[parameters.starname]
        #kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)
        ##print("dir(kep_fit): ", dir(kep_fit))
        
        ##make fit_parameters
        #fit_parameters = er.output_fit_parameters(kep_fit)
        
        
        ###plot_matched_RVs(res, output_folder)
        ###plot_fit(res, output_folder, fit_parameters)
        
        ###NOTE Plot Time Series
        ##cws.plot_time_series(res, output_folder, fit_parameters, format_im = "pdf"
                         ##)
        ##cws.plot_time_series(res, output_folder, fit_parameters, phased = True, format_im = "pdf"
                         ##)
                         
        ##NOTE Plot SNR correlation
        #cws.plot_residuals_vs_snr(res, output_folder, fit_parameters, format_im = "pdf")
        
###END
##BEGIN Initial RV Results : SNR
    #run_name = "SNR_fix"
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline_2/pipeline_pp_reruns/"
    #file_list = ["Baseline+SNRresults_GJ436_Kstar0_Kt3_min_snr_60.hdf5",
                 #"Baseline+SNRresults_GJ1148_Kstar0_Kt3_min_snr_60.hdf5"
                 #]

    #output_dir = output_dir_master + "{0}/".format(run_name)
    #os.makedirs(output_dir, exist_ok = True)
    
    #for index,f in enumerate(file_list):
        #wobble_file = results_dir + f
        #dataset_name = os.path.splitext(os.path.split(wobble_file)[1])[0]
        #output_folder = output_dir + dataset_name + "/"
        #os.makedirs(output_folder, exist_ok = True)
        
        ##make res
        #res = cws.load_results(wobble_file, correct_w_for_drift = True, correct_w_for_NZP = True)
        ##NOTE remove corrections for precorrected wobble
        
        ##make vels
        #vels_file = res.eval_vels_serval(output_dir)
        
        ##make kep_fit
        #parameters = rw.read_parameters_from_results(wobble_file)
        #n_planets_max = n_planet_dict[parameters.starname]
        #kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)
        
        ##make fit_parameters
        #fit_parameters = er.output_fit_parameters(kep_fit)
        
        ###NOTE Plot Time Series
        ##cws.plot_time_series(res, output_folder, fit_parameters, format_im = "pdf"
                         ##)
        ##cws.plot_time_series(res, output_folder, fit_parameters, phased = True, format_im = "pdf"
                         ##)
                         
        ##NOTE Plot SNR correlation
        #cws.plot_residuals_vs_snr(res, output_folder, fit_parameters, format_im = "pdf")
##END


