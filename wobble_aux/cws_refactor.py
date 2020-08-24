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
from astropy.stats import sigma_clip

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
mpl.rcParams['lines.markeredgecolor']='k'#black edges on symbols
##########################

def load_results(wobble_file, arm = "vis", **kwargs):
    if arm == "vis":
        serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_VIS/" #NOTE only for VIS
    elif arm == "nir":
        serval_dir = os.path.dirname(os.path.abspath(__file__)) + "/" + "../data/servaldir/CARM_NIR/" #NOTE only for NIR
    else:
        raise Exception("Invalid 'arm' specification. Must be 'vis' or 'nir'")
    parameters = rw.read_parameters_from_results(wobble_file)
    carmenes_object_ID = name_dict[parameters.starname]
    bary_starname = simbad_dict[parameters.starname]
    
    res = pr.Results_ws(wobble_file
                , serval_dir
                , carmenes_object_ID
                , bary_starname
                , load_bary = True
                , archive = True)
    res.apply_corrections(**kwargs)
    return res

def keplerian_rv(fit_parameters, t):
        t_epoch = fit_parameters[0]
        orbital_parameters = fit_parameters[1]
        total_rv = 0
        for parameters in orbital_parameters:
            total_rv = total_rv + rv.radial_velocity_M0(t , parameters[0], parameters[1], parameters[2], parameters[3],t_epoch , parameters[4])
        return total_rv
    
def fit_residuals(dates, rvs,  fit_parameters):
    if len(dates) == len(rvs):
        #residuals = np.array([rvs[i] - keplerian_rv(fit_parameters, dates[i]) 
                              #- np.nanmean(rvs) for i in range(len(rvs))])
        residuals = np.asarray([rvs[i] - keplerian_rv(fit_parameters, dates[i]) 
                               for i in range(len(rvs))])
    else:
        raise Exception("Missmatch: len(dates) =|= len(rvs)")
    return residuals

def error_stats(residuals, errors, kappa = 5):
    #correct for rv offset
    weights = 1 / errors**2
    residual_offset = np.average(residuals, weights = weights)
    # NOTE standart deviation autopmatically corrects for the offset -> maybe rms should be used or at least checked
    res_corr = residuals - residual_offset
    rms = np.sqrt(np.mean(res_corr**2))
    sigma = np.nanstd(res_corr)
    if len(sigma_clip(res_corr,sigma = kappa, masked = False)) < 2:
        sigma_clipped = np.nan
    else:
        sigma_clipped = np.nanstd(sigma_clip(res_corr,sigma = kappa, masked = False) )
    ind_bad = np.where(np.absolute(res_corr) > kappa * errors)[0]
    ind_good= np.where(np.absolute(res_corr) <= kappa * errors)[0]
    sigma_error_clipped = np.nanstd(res_corr[ind_good])
    return sigma, sigma_clipped,sigma_error_clipped, ind_bad, ind_good, residual_offset, rms
    
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
    fig.clf()
    
def plot_time_series(res, output_folder, fit_parameters,
                     format_im = 'png', #'pdf' or png
                     phased = False,
                     fit_label = "SERVAL fit"
                     ):
    # Plot wobble and Serval with fit and errors
    
    #Calculate residuals
    wob_residuals = fit_residuals(res.w_dates, res.w_RVs_barycorr, fit_parameters)
    wob_err_stats = error_stats(wob_residuals, res.w_RVs_er)
    ser_residuals = fit_residuals(res.ser_avcn[:,0], res.ser_avcn[:,1], fit_parameters)
    ser_err_stats = error_stats(ser_residuals, res.ser_avcn[:,2])
    
    #err_stats = [wob_err_stats, ser_err_stats]
    
    jd_offset = 2450000 # for cleaner x-axis labels
    
    
    #Plot definition
    color = ['C0', 'C1', 'C2', 'C3']
    symbol = ['o', 'D', 'v', 's']
    markersize = [3, 3, 3, 3] 
    alpha = [1, 1, 1, 1]
    
    model_color = 'C7'
    model_lw = '1.0'
    
    fig = plt.figure(0, figsize=(8,6.5))
    #ax1=plt.gca()
    
    #plt.subplots_adjust(hspace=0.005)
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    gs.update(#wspace=0.05
            hspace = 0.005
        )

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    w_dates_sorted = np.sort(res.w_dates)
    xlst = np.linspace(w_dates_sorted[0], w_dates_sorted[-1], num = 3000)
    ylst = np.array([keplerian_rv(fit_parameters, date) for date in xlst])
    
    if phased == True:
        # HACK
        jd_offset = w_dates_sorted[0]
        #set jd_offset to first entry in sorted dates list to lock phase even for slightly different period fits (so long as datasets contain same earliest entry (not guaranteed) #NOTE Does not quite solve the problem, as different peiod still de phases results 
        xlst = xlst - jd_offset
        period = fit_parameters[1][0][1]
        xlst_phase = xlst % period
        xlst = xlst_phase[np.argsort(xlst_phase)]
        ylst = ylst[np.argsort(xlst_phase)]
    else:
        xlst = xlst - jd_offset
        
    ax1.plot(xlst, ylst ,"-", color = model_color, linewidth=model_lw, label = fit_label)
    
    zero_point_T = range((int(min(res.w_dates - jd_offset))),
                         (int(max(res.w_dates - jd_offset))),10)
    if phased == True:
        zero_point_T = np.linspace(0, period, num = 10)
    zero_point   = np.zeros(len(zero_point_T))
    ax2.plot(zero_point_T ,zero_point,'-', linewidth=model_lw, color=model_color)
    
    #make plot lists
    plot_dates = np.array([res.w_dates, res.ser_avcn[:,0]]) - jd_offset
    if phased == True:
        plot_dates = plot_dates % period
    plot_rvs = [res.w_RVs_barycorr, res.ser_avcn[:,1]]
    plot_errs = [res.w_RVs_er, res.ser_avcn[:,2]]
    #plot_labels = [
                    #r"\texttt{wobble}" + ", rms = {0:.3f}, sigma\_clipped = {1:.3f}, err\_clipped = {2:.3f} ".format(wob_err_stats[0], wob_err_stats[1], wob_err_stats[2]),
                    #"SERVAL, rms = {0:.3f}, sigma\_clipped = {1:.3f}, err\_clipped = {2:.3f}".format(ser_err_stats[0], ser_err_stats[1], ser_err_stats[2])
                   #]
    plot_labels = [
                    r"\texttt{wobble}" + ", rms = {0:.3f}, 5$\cdot$error clipped = {2:.3f} ".format(wob_err_stats[0], wob_err_stats[1], wob_err_stats[2]),
                    "SERVAL, rms = {0:.3f}, 5$\cdot$error clipped = {2:.3f}".format(ser_err_stats[0], ser_err_stats[1], ser_err_stats[2])
                   ]
    plot_residuals = [wob_residuals, ser_residuals]
    plot_residual_offsets = [wob_err_stats[5], ser_err_stats[5]]
    ind_bad = [wob_err_stats[3], ser_err_stats[3]]
    cut_labels = [r"\texttt{wobble} 5$\cdot$error clipped points", "SERVAL 5$\cdot$error clipped points"]
    
    #Plot
    for i in range(len(plot_dates)):
        ax1.errorbar(plot_dates[i],
                    plot_rvs[i] 
                    - plot_residual_offsets[i],
                    yerr = plot_errs[i],
                    alpha=alpha[i], fmt=symbol[i], linestyle='None', markersize = markersize[i], color=color[i], capsize = 0, elinewidth=1,mew=0.1, 
                    label= plot_labels[i]
                    )
        ax2.errorbar(plot_dates[i],
                    plot_residuals[i] - plot_residual_offsets[i],
                    yerr = plot_errs[i],
                    alpha=alpha[i], fmt=symbol[i], linestyle='None', markersize = markersize[i], color=color[i], capsize = 0, elinewidth=1,mew=0.1
                    )
    
        if True: #error_clipping == True:
            ax1.plot(plot_dates[i][ind_bad[i]],
                        plot_rvs[i][ind_bad[i]]
                        - plot_residual_offsets[i],
                        "x",
                        #yerr = plot_errs[i],
                        alpha=alpha[i], linestyle='None', markersize = markersize[i]*2, color = "k", 
                        #capsize = 0, elinewidth=1,mew=0.1, 
                        label= cut_labels[i]
                        )
            
            ax2.plot(plot_dates[i][ind_bad[i]],
                        plot_residuals[i][ind_bad[i]]
                        - plot_residual_offsets[i],
                        "x",
                        #yerr = plot_errs[i],
                        alpha=alpha[i], linestyle='None', markersize = markersize[i]*2, color = "k", 
                        #capsize = 0, elinewidth=1,mew=0.1, 
                        label= cut_labels[i]
                        )

    ax1.set_ylabel(r'RV [m/s]',fontsize=15, rotation = 'vertical')
    #ax1.set_xlabel('jd')
    
    
    ax2.set_xlabel(r'JD [d] - {}'.format(jd_offset),fontsize=16)
    if phased == True:
        ax2.set_xlabel(r'Phase [d]',fontsize=16)
    ax2.set_ylabel(r'Residuals [m/s]',fontsize=15, rotation = 'vertical')
    
    ax2.locator_params(axis="x", nbins=9)
    ax2.locator_params(axis="y", min_n_ticks = 3
                       )
    #TODO make custom pruning of uppper tick (do not plot ticks in upper 10%)
    #so that ax2 tick does nto interfere with  ax1 tick
    y_loc = ax2.yaxis.get_majorticklocs()
    #print("y_loc: ", y_loc)
    #print("y_loc[1:-2]: ", y_loc[1:-2])
    #print("ylim: ", ax2.get_ylim())
    y2_min, y2_max = ax2.get_ylim()
    y_loc = [y for y in y_loc if y2_min < y < y2_max - (y2_max - y2_min)*0.1]
    #print("y_loc: ", y_loc)
    ax2.set_yticks(y_loc)
    
    
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    #ax1.grid(True)
    #ax2.grid(True)
    #plt.tight_layout()
    #add padding for legend
    ymin, ymax = ax1.get_ylim()
    ax1.set_ylim(ymin, ymax + (ymax-ymin)*0.14)#TODO make this scale with legend size
    
    h, l = ax1.get_legend_handles_labels()
    #print("labels:", l)
    l[1] = r"5$\cdot$error clipped points"
    select = [0,1,3,4]
    ax1.legend([h[i] for i in select], [l[i] for i in select],
               #shadow = True,
               framealpha = 0.5,
               loc = "upper left",
               fontsize = "small",
               
               ncol = 2
               )
    
    dpi = 300
    if phased == True:
        plt.savefig(output_folder + "time_series_phased_pb" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    else:
        plt.savefig(output_folder + "time_series" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla()
    fig.clf()
    plt.close()
    
def plot_time_series_parallel(res, output_folder, 
                              fit_parameters_wob,
                              fit_wob_label,
                              fit_parameters_ser,
                              fit_ser_label,
                     format_im = 'png', #'pdf' or png
                     phased = False
                     ):
    # Plot wobble and Serval with fit and errors
    
    #Plot definition
    color = ['C0', 'C1', 'C2', 'C3']
    symbol = ['o', 'D', 'v', 's']
    markersize = [3, 3, 3, 3] 
    alpha = [1, 1, 1, 1]
    
    model_color = 'C7'
    model_lw = '1.0'
    
    fig = plt.figure(0, figsize=(12,6.5))
    #ax1=plt.gca()
    
    #plt.subplots_adjust(hspace=0.005)
    gs = mpl.gridspec.GridSpec(2, 2, height_ratios=[3, 1]) 
    gs.update(wspace=0.005,
            hspace = 0.005
        )

    #ax1 = plt.subplot(gs[0,0])
    #ax2 = plt.subplot(gs[1,1])
    ############################
    
    #Calculate residuals
    fit_parameters_list = [fit_parameters_wob, fit_parameters_ser]
    fit_labels = [fit_wob_label, fit_ser_label]
    
    wob_residuals = fit_residuals(res.w_dates, res.w_RVs_barycorr, fit_parameters_wob)
    wob_err_stats = error_stats(wob_residuals, res.w_RVs_er)
    ser_residuals = fit_residuals(res.ser_avcn[:,0], res.ser_avcn[:,1], fit_parameters_ser)
    ser_err_stats = error_stats(ser_residuals, res.ser_avcn[:,2])
    
    jd_offset = 2450000 # for cleaner x-axis labels
    fit_periods = [fit_parameters[1][0][1] for fit_parameters in fit_parameters_list]
    
    #w_dates_sorted = np.sort(res.w_dates)
    #xlst = np.linspace(w_dates_sorted[0], w_dates_sorted[-1], num = 3000)
    #ylst = np.array([keplerian_rv(fit_parameters, date) for date in xlst])
    
    #if phased == True:
        ## HACK
        #jd_offset = w_dates_sorted[0]
        ##set jd_offset to first entry in sorted dates list to lock phase even for slightly different period fits (so long as datasets contain same earliest entry (not guaranteed) #NOTE Does not quite solve the problem, as different peiod still de phases results 
        #xlst = xlst - jd_offset
        #period = fit_parameters[1][0][1]
        #xlst_phase = xlst % period
        #xlst = xlst_phase[np.argsort(xlst_phase)]
        #ylst = ylst[np.argsort(xlst_phase)]
    #else:
        #xlst = xlst - jd_offset
    
    #zero_point_T = range((int(min(res.w_dates - jd_offset))),
                        #(int(max(res.w_dates - jd_offset))),10)
    #if phased == True:
        #zero_point_T = np.linspace(0, period, num = 10)
    #zero_point   = np.zeros(len(zero_point_T))
    
    #TODO Make this just work as loop
    #make plot lists
    plot_dates = np.array([res.w_dates, res.ser_avcn[:,0]]) - jd_offset
    w_dates_sorted = np.sort(res.w_dates)
    if phased == True:
        jd_offset = w_dates_sorted[0]
        #set jd_offset to first entry in sorted dates list to lock phase even for slightly different period fits (so long as datasets contain same earliest entry (not guaranteed) #NOTE Does not quite solve the problem, as different peiod still de phases results 
        plot_dates = np.array([res.w_dates, res.ser_avcn[:,0]]) - jd_offset
        plot_dates = [plot_dates[i] % fit_periods[i] for i in range(len(plot_dates))]
        
        
    plot_rvs = [res.w_RVs_barycorr, res.ser_avcn[:,1]]
    plot_errs = [res.w_RVs_er, res.ser_avcn[:,2]]
    #plot_labels = [
                    #r"\texttt{wobble}" + ", rms = {0:.3f}, sigma\_clipped = {1:.3f}, err\_clipped = {2:.3f} ".format(wob_err_stats[0], wob_err_stats[1], wob_err_stats[2]),
                    #"SERVAL, rms = {0:.3f}, sigma\_clipped = {1:.3f}, err\_clipped = {2:.3f}".format(ser_err_stats[0], ser_err_stats[1], ser_err_stats[2])
                   #]
    plot_labels = [
                    r"\texttt{wobble}" + ", rms = {0:.3f}, 5$\cdot$error clipped = {2:.3f} ".format(wob_err_stats[0], wob_err_stats[1], wob_err_stats[2]),
                    "SERVAL, rms = {0:.3f}, 5$\cdot$error clipped = {2:.3f}".format(ser_err_stats[0], ser_err_stats[1], ser_err_stats[2])
                   ]
    
    
    #plot_residuals = [wob_residuals, ser_residuals] = 
    #[fit_residuals(plot_dates[i], plot_rvs[i], fit_parameters[i])
                      #for i in range(len(plot_dates))]
    #plot_err_stats = [wob_err_stats, ser_err_stats] = 
    #[error_stats(plot_residuals[i], plot_errs[i])
                      #for i in range(len(plot_dates))]
    
    plot_residuals = [wob_residuals, ser_residuals]
    plot_residual_offsets = [wob_err_stats[5], ser_err_stats[5]]
    ind_bad = [wob_err_stats[3], ser_err_stats[3]]
    cut_labels = [r"\texttt{wobble} 5$\cdot$error clipped points", "SERVAL 5$\cdot$error clipped points"]
    
    #Plot
    for i in range(len(plot_dates)):
        ax1 = plt.subplot(gs[0,i])
        ax2 = plt.subplot(gs[1,i])
        
        #### Plot Fit
        xlst = np.linspace(w_dates_sorted[0], w_dates_sorted[-1], num = 3000)
        ylst = np.array([keplerian_rv(fit_parameters_list[i], date) for date in xlst])
        if phased == True:
            # HACK
            xlst = xlst - jd_offset
            period = fit_periods[i]
            xlst_phase = xlst % period
            xlst = xlst_phase[np.argsort(xlst_phase)]
            ylst = ylst[np.argsort(xlst_phase)]
        else:
            xlst = xlst - jd_offset
        
        zero_point_T = range((int(min(res.w_dates - jd_offset))),
                            (int(max(res.w_dates - jd_offset))),10)
        if phased == True:
            zero_point_T = np.linspace(0, period, num = 10)
        zero_point   = np.zeros(len(zero_point_T))
        
        
        ax1.plot(xlst, ylst ,"-", color = model_color, linewidth=model_lw,
                 label = fit_labels[i])
        ax2.plot(zero_point_T ,zero_point,'-', linewidth=model_lw, color=model_color)
        
        #### Plot data
        ax1.errorbar(plot_dates[i],
                    plot_rvs[i] 
                    - plot_residual_offsets[i],
                    yerr = plot_errs[i],
                    alpha=alpha[i], fmt=symbol[i], linestyle='None', markersize = markersize[i], color=color[i], capsize = 0, elinewidth=1,mew=0.1, 
                    label= plot_labels[i]
                    )
        ax2.errorbar(plot_dates[i],
                    plot_residuals[i] - plot_residual_offsets[i],
                    yerr = plot_errs[i],
                    alpha=alpha[i], fmt=symbol[i], linestyle='None', markersize = markersize[i], color=color[i], capsize = 0, elinewidth=1,mew=0.1
                    )
    
        if True: #error_clipping == True:
            ax1.plot(plot_dates[i][ind_bad[i]],
                        plot_rvs[i][ind_bad[i]]
                        - plot_residual_offsets[i],
                        "x",
                        #yerr = plot_errs[i],
                        alpha=alpha[i], linestyle='None', markersize = markersize[i]*2, color = "k", 
                        #capsize = 0, elinewidth=1,mew=0.1, 
                        label= cut_labels[i]
                        )
            
            ax2.plot(plot_dates[i][ind_bad[i]],
                        plot_residuals[i][ind_bad[i]]
                        - plot_residual_offsets[i],
                        "x",
                        #yerr = plot_errs[i],
                        alpha=alpha[i], linestyle='None', markersize = markersize[i]*2, color = "k", 
                        #capsize = 0, elinewidth=1,mew=0.1, 
                        label= cut_labels[i]
                        )  

        ax1.set_ylabel(r'RV [m/s]',fontsize=15, rotation = 'vertical')
        #ax1.set_xlabel('jd')
        
        
        ax2.set_xlabel(r'JD [d] - {}'.format(jd_offset),fontsize=16)
        if phased == True:
            ax2.set_xlabel(r'Phase [d]',fontsize=16)
        ax2.set_ylabel(r'Residuals [m/s]',fontsize=15, rotation = 'vertical')
        
        
        if i == 1: # Load y tick settings and put them onn the right side
            #ax1.set_yticks(yticks_ax1_left[0])
            #ax1.set_yticklabels(yticks_ax1_left[1])
            #ax2.set_yticks(yticks_ax2_left[0])
            #ax2.set_yticklabels(yticks_ax2_left[1])
            
            ax1.set_ylim(ax1_ylims_left)
            ax2.set_ylim(ax2_ylims_left)
            
            ax1.yaxis.tick_right()
            ax2.yaxis.tick_right()
            ax1.yaxis.set_label_position("right")
            ax2.yaxis.set_label_position("right")  
        
        ax2.locator_params(axis="x", nbins=9)
        ax2.locator_params(axis="y", min_n_ticks = 3
                        )
        #TODO make custom pruning of uppper tick (do not plot ticks in upper 10%)
        #so that ax2 tick does nto interfere with  ax1 tick
        y_loc = ax2.yaxis.get_majorticklocs()
        #print("y_loc: ", y_loc)
        #print("y_loc[1:-2]: ", y_loc[1:-2])
        #print("ylim: ", ax2.get_ylim())
        y2_min, y2_max = ax2.get_ylim()
        y_loc = [y for y in y_loc if y2_min < y < y2_max - (y2_max - y2_min)*0.1]
        #print("y_loc: ", y_loc)
        ax2.set_yticks(y_loc)
        
        
        #plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #ax1.grid(True)
        #ax2.grid(True)
        #plt.tight_layout()
        #add padding for legend
        if i == 0 : # Setup yaxis with first plot settings
            ymin, ymax = ax1.get_ylim()
            ax1.set_ylim(ymin, ymax + (ymax-ymin)*0.20)#TODO make this scale with legend size
            
            #HACK save Ytick locations to duplicate for Serval plot
            #print("yticks(): ", plt.yticks()) # Doesn't work for subplots
            #yticks_ax1_left = ax1.get_yticks() , ax1.get_yticklabels()
            #yticks_ax2_left = ax2.get_yticks() , ax2.get_yticklabels()
            ax1_ylims_left = ax1.get_ylim()
            ax2_ylims_left = ax2.get_ylim()
    
            
        #if i == 1: # Load y tick settings and put them onn the right side
            ##ax1.set_yticks(yticks_ax1_left[0])
            ##ax1.set_yticklabels(yticks_ax1_left[1])
            ##ax2.set_yticks(yticks_ax2_left[0])
            ##ax2.set_yticklabels(yticks_ax2_left[1])
            
            #ax1.set_ylim = ax1_ylims_left
            #ax2.set_ylim = ax2_ylims_left
            
            #ax1.yaxis.tick_right()
            #ax2.yaxis.tick_right()
            #ax1.yaxis.set_label_position("right")
            #ax2.yaxis.set_label_position("right")  
        
        #Make Legend
        ax1.legend(fontsize = "small")
        #h, l = ax1.get_legend_handles_labels()
        ##print("labels:", l)
        #l[1] = r"5$\cdot$error clipped points"
        #select = [0,1,3,4]
        #ax1.legend([h[i] for i in select], [l[i] for i in select],
                ##shadow = True,
                #framealpha = 0.5,
                #loc = "upper left",
                #fontsize = "small",
                
                #ncol = 2
                #)
                    
    dpi = 300
    if phased == True:
        plt.savefig(output_folder + "time_series_parallel_phased_pb" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    else:
        plt.savefig(output_folder + "time_series_parallel" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    fig.clf()
    plt.close()
    

def plot_residuals_vs_snr(res, output_folder, fit_parameters,
                     format_im = 'png', #'pdf' or png
                     ):
    fig = plt.figure(0, figsize=(8,6))
    ax1= plt.gca()
    
    #Calculate residuals
    wob_residuals = fit_residuals(res.w_dates, res.w_RVs_barycorr, fit_parameters)
    wob_err_stats = error_stats(wob_residuals, res.w_RVs_er)
    ser_residuals = fit_residuals(res.ser_avcn[:,0], res.ser_avcn[:,1], fit_parameters)
    ser_err_stats = error_stats(ser_residuals, res.ser_avcn[:,2])
    
    plot_residual_offsets = [wob_err_stats[5], ser_err_stats[5]]
    
    ax1.plot(res.ser_addinfo[:,0],
             np.absolute(wob_residuals - plot_residual_offsets[0]),
             "o", color = "C0",
             label = r"\texttt{wobble}")
    #ax1.plot(res.ser_addinfo[:,0],
             #np.absolute(ser_residuals - plot_residual_offsets[1]),
             #"D", color = "C1",
             #label = r"SERVAL")
    ax1.set_ylabel('Absolute RV residuals [m/s]')
    ax1.set_xlabel("SNR")
    #plt.title("RV residuals vs SNR correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
    plt.grid(True)
    #plt.tight_layout()
    plt.legend(shadow=True)
    
    dpi = 300
    plt.savefig(output_folder + "residuals_vs_SNR" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla()
    fig.clf()
    
def plot_residuals_vs_corr(res, output_folder, fit_parameters,
                     format_im = 'png', #'pdf' or png
                     ):
    fig = plt.figure(0, figsize=(8,6))
    ax1= plt.gca()
    
    #Calculate residuals
    wob_residuals = fit_residuals(res.w_dates, res.w_RVs_barycorr, fit_parameters)
    wob_err_stats = error_stats(wob_residuals, res.w_RVs_er)
    ser_residuals = fit_residuals(res.ser_avcn[:,0], res.ser_avcn[:,1], fit_parameters)
    ser_err_stats = error_stats(ser_residuals, res.ser_avcn[:,2])
    
    drift = res.ser_avcn[:,3]
    nzp = res.ser_avcn[:,9]
    corr = drift + nzp
    
    
    plot_residual_offsets = [wob_err_stats[5], ser_err_stats[5]]
    ind_bad = [wob_err_stats[3], ser_err_stats[3]]
    plot_errs = [res.w_RVs_er, res.ser_avcn[:,2]]
    
    markersize = 3
    model_color = 'C7'
    model_lw = '1.0'
    
    #plot 1 to - 1 correlation line
    xlst = np.linspace(np.amin(corr), np.amax(corr), num = 50)
    ax1.plot(xlst, -xlst ,"-", color = model_color, linewidth=model_lw, label="Slope -1")
    
    
    
    ax1.errorbar(corr,
             wob_residuals - plot_residual_offsets[0],
             yerr = plot_errs[0],
             fmt = "o", color = "C0", markersize = markersize,
             capsize = 0, elinewidth=1,mew=0.1, 
             label = r"\texttt{wobble}")
    #ax1.plot(res.ser_addinfo[:,0],
             #np.absolute(ser_residuals - plot_residual_offsets[1]),
             #"D", color = "C1",
             #label = r"SERVAL")
             
    ax1.plot(corr[ind_bad[0]],
            wob_residuals[ind_bad[0]] - plot_residual_offsets[0],
            "x",
            linestyle='None', markersize = markersize*2, color = "k", 
            label= r"5$\cdot$error clipped points"
            )
             
             
    ax1.set_ylabel('RV residuals [m/s]')
    ax1.set_xlabel("drift + NZP [m/s]")
    #plt.title("RV residuals vs SNR correlation for "+i[0]+" ("+i[2]+") "+" - "+i[1]+";")
    plt.grid(True)
    #plt.tight_layout()
    plt.legend()
    
    dpi = 300
    plt.savefig(output_folder + "residuals_vs_corr" + '.{}'.format(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla()
    fig.clf()

if __name__ == "__main__":
    
#BEGIN Initial RV Results
    run_name = "Initial_RV"
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
        res = load_results(wobble_file, correct_w_for_drift = True, correct_w_for_NZP = True)
        
        #make vels
        vels_file = res.eval_vels_serval(output_dir)
        
        #make kep_fit
        parameters = rw.read_parameters_from_results(wobble_file)
        n_planets_max = n_planet_dict[parameters.starname]
        kep_fit = er.vels_to_kep_fit(dataset_name, vels_file, n_planets_max = n_planets_max)
        #print("dir(kep_fit): ", dir(kep_fit))
        
        #make fit_parameters
        fit_parameters = er.output_fit_parameters(kep_fit)
        
        
        #plot_matched_RVs(res, output_folder)
        #plot_fit(res, output_folder, fit_parameters)
        plot_time_series(res, output_folder, fit_parameters, format_im = "pdf"
                         )
        plot_time_series(res, output_folder, fit_parameters, phased = True, format_im = "pdf"
                         )
##END
