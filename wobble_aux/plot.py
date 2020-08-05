
import os
import wobble
import numpy as np
import matplotlib.pyplot as plt

import run_wobble as rw
import process_results as pr

#import exostriker functions from  evaluate results
import evaluate_results as er
import matplotlib as mpl

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
    
########################################################################################
#Test of  importing eval functions
    ##BEGIN
    #run_name = "drop_points_test_TeX_plot"
    #results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline/pipeline_test_0/"
    #file_list = ["results_GJ436_Kstar0_Kt3_baseline_0.hdf5",
        #"results_GJ1148_Kstar0_Kt3_baseline_0.hdf5",
        ##"results_GJ3473_Kstar0_Kt3_baseline_0.hdf5",
        ##"results_YZ Cet_Kstar0_Kt3_baseline_0.hdf5",
        ##"results_GJ15A_Kstar0_Kt3_baseline_0.hdf5",
        ##"results_GJ176_Kstar0_Kt3_baseline_0.hdf5", "results_GJ536_Kstar0_Kt3_baseline_0.hdf5", "results_GJ3512_Kstar0_Kt3_baseline_0.hdf5", "results_Wolf294_Kstar0_Kt3_baseline_0.hdf5", "results_GJ876_Kstar0_Kt3_baseline_0.hdf5" , "results_Teegarden_Kstar0_Kt3_baseline_0.hdf5", "results_Barnard_Kstar0_Kt3_baseline_0.hdf5"
                 #]

    #er.eval_complete(run_name, file_list, er.name_dict, er.simbad_dict, results_dir, n_planet_dict = er.n_planet_dict, drop_pointsQ = True)
    ##END
    
#Test of plot_rv_comparison: Basic impersonation of "INITIAL RV results" presentation plot
    #BEGIN
    run_name = "pres_rerun/baseline/"
    results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline_2/pipeline_pp_reruns/"
    file_list = ["Baseline+SNRresults_GJ436_Kstar0_Kt3_min_snr_5.hdf5",
        "Baseline+SNRresults_GJ436_Kstar0_Kt3_min_snr_5.hdf5"
                 ]

    #run_name, file_list, 
    name_dict = er.name_dict
    simbad_dict = er.simbad_dict
    n_planet_dict = er.n_planet_dict
    drop_pointsQ = False
    
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    output_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/fin_plots/{0}/".format(run_name)
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
        
        #HACK
        kep_fit_list = [[],[]]

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
        kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)

        #remove worst fitting poinnts and run again if option is True
        if drop_pointsQ:
            n = 4
            vels_file = er.drop_points(kep_fit, vels_file, n = n)
            dataset_name = dataset_name + "_dropped_{}".format(n)
            #fit again
            kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max) 
            
        kep_fit_list[0] = kep_fit
            
        #HACK SAME FOR SERVAL
        vels_file = res.eval_vels_serval(output_dir)#NOTE Only change from wobble implementation  so far
        #END EDIT

        #2.vels to kep_fit
        dataset_name = parameters.starname +" Serval_avcn"
        print("dataset_name: ", dataset_name ,"vels_vile: ", vels_file)
        if parameters.starname in n_planet_dict:
            n_planets_max = n_planet_dict[parameters.starname]
            print(parameters.starname, ": ", n_planets_max)
        else:
            n_planets_max = 2
        kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)
        
        #remove worst fitting poinnts and run again if option is True
        if drop_pointsQ:
            n = 4
            vels_file = er.drop_points(kep_fit, vels_file, n = n)
            dataset_name = dataset_name + "_dropped_{}".format(n)
            #fit again
            kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)
        
        kep_fit_list[1] = kep_fit
        
        #3.plot
        output_file = output_dir + os.path.splitext(os.path.split(wobble_file)[1])[0]
        #er.plot_rv_comparison(kep_fit_list, output_file)
        er.plot_rv_comparison(kep_fit_list, output_file, phased = True, legend_labels = [r'\texttt{wobble}', 'SERVAL' ], ip = 1)
    #END
