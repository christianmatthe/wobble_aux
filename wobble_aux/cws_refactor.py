"""

@author: cmatthe
"""
#%% imports
import sys
import rv_solver as rv
import process_results as pr
import evaluate_results as er
import run_wobble as rw

import matplotlib as mpl
import matplotlib.pyplot as plt

import os
import numpy as np
################
name_dict = er.name_dict #connects starname to carmenes ID
simbad_dict = er.simbad_dict #connects starname to simbad iidentifiers that work with barycorr
n_planet_dict = er.n_planet_dict #dictionary of max number of planets to be fit
###########
#simple plot of RVs (only matched)
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
#mpl.rcParams['lines.markeredgecolor']='k'#black edges on symbols
##########################

def load_results(wobble_file):
    serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    parameters = rw.read_parameters_from_results(wobble_file)
    carmenes_object_ID = name_dict[parameters.starname]
    bary_starname = simbad_dict[parameters.starname]
    
    res = pr.Results_ws(wobble_file
                , serval_dir
                , carmenes_object_ID
                , bary_starname
                , load_bary = True
                , archive = True)
    res.apply_corrections()
    return res

def keplerian_rv(fit_parameters, t):
        t_epoch = fit_parameters[0]
        orbital_parameters = fit_parameters[1]
        total_rv = 0
        for parameters in orbital_parameters:
            total_rv = total_rv + rv.radial_velocity_M0(t , parameters[0], parameters[1], parameters[2], parameters[3],t_epoch , parameters[4])
        return total_rv
    
def residuals(rvs, dates, fit_parameters):
    #TODO rewrite as straight finction of date and rv list in and residuals out
    fit_RVs = [keplerian_rv(fit_parameters, date) for date in dates]
    residuals = res.w_RVs_barycorr - np.nanmean(res.w_RVs_barycorr) - fit_RVs
    return residuals
    
###################### Plots
def plot_matched_RVs(res, output_folder):
    ##simple plot of RVs (only matched)
    
    #Plot definition
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    
    ax1.plot(res.w_dates, res.w_RVs_barycorr - np.nanmean(res.w_RVs_barycorr),"x", label="Wobble bary")
    ax1.plot(res.ser_avcn[:,0], res.ser_avcn[:,5] - np.nanmean(res.ser_avcn[:,5]), "x", label= "SERVAL")
    ax1.plot(res.ser_avcn[:,0], res.ser_avcn[:,1] - np.nanmean(res.ser_avcn[:,1]), "x", label= "SERVAL Corr")
    ax1.set_ylabel("RVs [m/s]")
    ax1.set_xlabel('jd')
    #plt.title("matched RVs for "+i[0]+" ("+i[2]+") "+" - "+i[1]+"; (RVs only BERV-corrected; no FP-drift, no SA, no NZP corrections)")
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    
    
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(output_folder + "basic_RV" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla()
    
def plot_fit(res, output_folder, fit_parameters):
    ##simple plot of RVs (only matched)
    
    #Plot definition
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    
    w_dates_sorted = np.sort(res.w_dates)
    xlst = np.linspace(w_dates_sorted[0], w_dates_sorted[-1], num = 3000)
    ylst = np.array([keplerian_rv(fit_parameters, date) for date in xlst])
    #date_sort = res.w_dates.argsort()
    #print("date_sort", date_sort)
    #print("res.w_dates[[1,2]]:", res.w_dates[[1,2]])
    #print("fit_y[[1,2]]:", fit_y[[1,2]]) 
    ax1.plot(xlst, ylst ,"-", color = "C7", label="fit")
    ax1.plot(res.w_dates, res.w_RVs_barycorr - np.nanmean(res.w_RVs_barycorr),"x", label="Wobble bary")
    ax1.plot(res.ser_avcn[:,0], res.ser_avcn[:,5] - np.nanmean(res.ser_avcn[:,5]), "x", label= "SERVAL")
    ax1.plot(res.ser_avcn[:,0], res.ser_avcn[:,1] - np.nanmean(res.ser_avcn[:,1]), "x", label= "SERVAL Corr")
    ax1.set_ylabel("RVs [m/s]")
    ax1.set_xlabel('jd')
    #plt.title("matched RVs for "+i[0]+" ("+i[2]+") "+" - "+i[1]+"; (RVs only BERV-corrected; no FP-drift, no SA, no NZP corrections)")
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    
    
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(output_folder + "fit" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla()
    

if __name__ == "__main__":
    run_name = "cws_test"
    results_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/pipeline_2/pipeline_pp_reruns/"
    file_list = ["Baseline+SNRresults_GJ436_Kstar0_Kt3_min_snr_5.hdf5",
                 "Baseline+SNRresults_GJ1148_Kstar0_Kt3_min_snr_5.hdf5"
                 ]

    
    #serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    output_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../results/fin_plots/{0}/".format(run_name)
    os.makedirs(output_dir, exist_ok = True)
    
    for index,f in enumerate(file_list):
        wobble_file = results_dir + f
        dataset_name = os.path.splitext(os.path.split(wobble_file)[1])[0]
        output_folder = output_dir + dataset_name + "/"
        os.makedirs(output_folder, exist_ok = True)
        
        #make res
        res = load_results(wobble_file)
        
        #make vels
        vels_file = res.eval_vels_serval(output_dir)
        
        #make kep_fit
        parameters = rw.read_parameters_from_results(wobble_file)
        n_planets_max = n_planet_dict[parameters.starname]
        kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)
        
        #make fit_parameters
        fit_parameters = er.output_fit_parameters(kep_fit)
        
        
        #plot_matched_RVs(res, output_folder)
        plot_fit(res, output_folder, fit_parameters)
        

