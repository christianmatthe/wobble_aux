import numpy as np
import h5py
import pdb
import matplotlib.pyplot as plt

from .utils import fit_continuum

#edited imports
import os
import scipy.interpolate as interpolate
#import sys
#sys.path.append('../../../wobble_aux')
#from run_wobble import Parameters

class AllDataDropped(Exception):
    def __init__(self):
        pass

class Data(object):
    """
    The data object: contains the spectra and associated data.
    All objects in `data` are numpy arrays or lists of arrays.
    Includes all orders and epochs.
    
    Parameters
    ----------
    filename : `str`
        Name of HDF5 file storing the data.
    filepath : `str` (default `../data/`)
        Path to append to filename.
    orders : `list` (default `None`)
        Indices of spectral orders to read in. If `None` include all. 
        Even if it's only one index, it must be a list.
    epochs : `int` or `list` (default `None`)
        Indices of epochs to read in. If `None` include all.
    min_flux : `float` (default `1.`)
        Flux in counts/pixel below which a pixel is masked out.
    max_norm_flux : `float` (default `2.`)
        Flux in normalized counts/pixel above which a pixel is masked out.
    padding : `int` (default `2`)
        Number of pixels to additionally mask on either side of a bad high pixel.
    min_snr : `float` (default `5.`)
        Mean SNR below which which we discard sections of data.
    log_flux : `bool` (default `True`)
        Determines whether fitting will happen using logarithmic flux (default) 
        or linear flux.
    """
    def __init__(self, filename, filepath='', 
                    orders = None, 
                    epochs = None,
                    min_flux = 1.,
                    max_norm_flux = 2.,
                    padding = 2,
                    min_snr = 5.,
                    log_flux = True,
                    chunkQ = False, #Flag to set for different behavior during chunk inn run_wobble (main aim diasble order and epoch cutting while retaining data point cutting)
                    **kwargs):
        origin_file = filepath+filename
        self.read_data(origin_file, orders=orders, epochs=epochs)
        self.mask_low_pixels(min_flux=min_flux, padding=padding, min_snr=min_snr)
        
        orders = np.asarray(self.orders)
        orders_start = orders
        
        if chunkQ == False: #Centralizes order and epoch cutting to top level only
            #HACK NOTE #NOTE change 25/10/2019 snr+data_suffix , temporarily removed order cutting since it breaks orderswise optimization, and was ignored before anyways
            snrs_by_order = np.sqrt(np.nanmean(self.ivars, axis=(1,2)))
            #disable order cutting here as it is done in continuum normalize
            orders_to_cut = snrs_by_order < min_snr
            print("Low SNR orders {0}: average SNR < {1:.0f} Note: not dropped".format(orders[orders_to_cut], min_snr))
            '''
            #if np.sum(orders_to_cut) > 0:
                #print("Data: Dropping orders {0} because they have average SNR < {1:.0f}".format(orders[orders_to_cut], min_snr))
                #orders = orders[~orders_to_cut]
                #self.read_data(origin_file, orders=orders, epochs=epochs) # overwrite with new data
                #self.mask_low_pixels(min_flux=min_flux, padding=padding, min_snr=min_snr)
            #if len(orders) == 0:
                #print("All orders failed the quality cuts with min_snr={0:.0f}.".format(min_snr))
                #raise AllDataDropped
            '''
            
            
            epochs = np.asarray(self.epochs)
            snrs_by_epoch = np.sqrt(np.nanmean(self.ivars, axis=(0,2)))
            epochs_to_cut = snrs_by_epoch < min_snr
            if np.sum(epochs_to_cut) > 0:
                print("Data: Dropping epochs {0} because they have average SNR < {1:.0f}".format(epochs[epochs_to_cut], min_snr))
                epochs = epochs[~epochs_to_cut]
                self.read_data(origin_file, orders=orders, epochs=epochs) # overwrite with new data
                self.mask_low_pixels(min_flux=min_flux, padding=padding, min_snr=min_snr)
            if len(epochs) == 0:
                print("All epochs failed the quality cuts with min_snr={0:.0f}.".format(min_snr))
                raise AllDataDropped
                #return
                
        elif chunkQ == True:
            #just in case these ar eneeded anywhere
            snrs_by_order = np.sqrt(np.nanmean(self.ivars, axis=(1,2)))
            
            epochs = np.asarray(self.epochs)
            snrs_by_epoch = np.sqrt(np.nanmean(self.ivars, axis=(0,2)))
            
            #Pixel masking already done above
            
            
        # log and normalize:
        self.ys = np.log(self.fluxes) 
        self.drop_orders = [] # empty list to be used by continuum normalize
        print("self.drop_orders define: ", self.drop_orders)
        self.continuum_normalize(**kwargs) 
        
        # HACK - optionally un-log it:
        if not log_flux:
            self.ys = np.exp(self.ys)
            self.ivars = self.flux_ivars
                
        # mask out high pixels:
        for r in range(self.R):
            bad = self.ys[r] > max_norm_flux
            self.ys[r][bad] = 1.
            for pad in range(padding): # mask out neighbors of high pixels
                bad = np.logical_or(bad, np.roll(bad, pad+1))
                bad = np.logical_or(bad, np.roll(bad, -pad-1))
            self.ivars[r][bad] = 0.
   
   
    def read_data(self, origin_file, orders = None, epochs = None):
        """Read origin file and set up data attributes from it"""
        # TODO: add asserts to check data are finite, no NaNs, non-negative ivars, etc
        with h5py.File(origin_file, 'r') as f: #04.06.2020 added explicit mode 'r' due to h5py 3.0 deprecation
            self.origin_file = origin_file
            if orders is None:
                orders = np.arange(len(f['data']))
            if epochs is None:
                self.N = len(f['dates']) # all epochs
                self.epochs = np.arange(self.N)
            else:
                self.epochs = epochs
                self.N = len(epochs)
                for e in epochs:
                    assert (e >= 0) & (e < len(f['dates'])), \
                        "epoch #{0} is not in datafile {1}".format(e, self.origin_file)
            self.epoch_groups = [list(np.arange(self.N))]
            self.fluxes = [f['data'][i][self.epochs,:] for i in orders]
            self.xs = [np.log(f['xs'][i][self.epochs,:]) for i in orders]
            self.flux_ivars = [f['ivars'][i][self.epochs,:] for i in orders] # ivars for linear fluxes
            self.pipeline_rvs = np.copy(f['pipeline_rvs'])[self.epochs]
            self.pipeline_sigmas = np.copy(f['pipeline_sigmas'])[self.epochs]
            self.dates = np.copy(f['dates'])[self.epochs]
            self.bervs = np.copy(f['bervs'])[self.epochs]
            self.drifts = np.copy(f['drifts'])[self.epochs]
            self.airms = np.copy(f['airms'])[self.epochs]
            self.filelist = [a.decode('utf8') for a in np.copy(f['filelist'])[self.epochs]]
            self.R = len(orders) # number of orders
            self.orders = orders # indices of orders in origin_file
            self.ivars = [self.fluxes[i]**2 * self.flux_ivars[i] for i in range(self.R)] # ivars for log(fluxes)
    
    def mask_low_pixels(self, min_flux = 1., padding = 2, min_snr = 5.):
        """Set ivars to zero for pixels and edge regions that are bad."""
        # mask out low pixels:
        #for carmenes has to be adjusted for negative and nan values due to FOX (according to Jonas and Adrian)
        for r in range(self.R):
            #bad = self.fluxes[r] < min_flux
            bad = np.logical_or(self.fluxes[r] < min_flux, np.isnan(self.fluxes[r]))
            self.fluxes[r][bad] = min_flux
            for pad in range(padding): # mask out neighbors of low pixels
                bad = np.logical_or(bad, np.roll(bad, pad+1))
                bad = np.logical_or(bad, np.roll(bad, -pad-1))
            self.flux_ivars[r][bad] = 0.
            self.ivars[r][bad] = 0.
            
        # find bad regions in masked spectra:
        for r in range(self.R):
            self.trim_bad_edges(r, min_snr=min_snr) # HACK
            
    #NOTE not one of wobbles standart functions
    def mask_telluric_mask_pixels(self, **kwargs):
        """Set ivars to zero for masked regions that are bad."""
        #NOTE needs to be applied after continuum noorm if use for inverted selection (Masked regions only)
        
        try:
            p = kwargs["parameters"]
        except:
            raise Exception("must pass parameters object to data object via kwargs in .Data for telluric mask")
        
        #TODO make this not Hard Coded
        file_dir = os.path.dirname(__file__)
        telluric_mask_file = file_dir + "/"+"../../../wobble_aux/carmenes_aux_files/" + "telluric_mask_carm_short.dat"
        #mask tellurics
        if telluric_mask_file is not None:
            telluric_mask = np.genfromtxt(telluric_mask_file)
            #extend mask (with value 0) to include earliest CARMENES orders
            mask = telluric_mask
            mask = np.insert(mask, 0, [0.0,0],axis = 0)
            mask = np.append(mask, [[50000, 0]], axis = 0)
            #create decision function basedon interpolation of Mask
            mask_function = interpolate.interp1d(mask[:,0],mask[:,1])
            mask_array = np.zeros(np.array(self.xs).shape)
            mask_bool = np.ones(np.array(self.xs).shape, dtype = bool)
            #loop over all entries
            for o, xs_order in enumerate(self.xs):
                for e, xs_epoch in enumerate(xs_order):
                    for l, xs_lambda in enumerate(xs_epoch):
                        if mask_function(np.exp(xs_lambda)) == 1:
                            mask_array[o,e,l] = 1
                            mask_bool[o,e,l] = False
        #TODO
        #[not bool for bool in mask_bool] to invert selection
        if p.mask tellurics == "mask" or p.mask tellurics == "inverted_mask":
            if p.mask tellurics == "inverted_mask":
                [not bool for bool in mask_bool] to invert selection
            for r in range(self.R):
                for n in range(self.N):
                    #Could probablly be done without the loops
                    self.flux_ivars[r][n][mask_bool[r,n]] = 0.
                    self.ivars[r][n][mask_bool[r,n]] = 0.
                    #TODO Check where 0 Ivars are dropped. I think they are just not weighted?
        else:
            if p.mask_tellurics != "no_mask":
                raise Exception("invalid p.mask_tellurics")
            
        
        

    def trim_bad_edges(self, r, window_width = 128, min_snr = 5.):
        """
        Find edge regions that contain no information and trim them.
        
        Parameters
        ----------
        r : `int`
            order index
        window_width : `int`
            number of pixels to average over for local SNR            
        min_snr : `float`
            SNR threshold below which we discard the data
        """
        for n in range(self.N):
            n_pix = len(self.xs[0][n])
            for window_start in range(n_pix - window_width):
                mean_snr = np.sqrt(np.nanmean(self.ivars[r][n,window_start:window_start+window_width]))
                if mean_snr > min_snr:
                    self.ivars[r][n,:window_start] = 0. # trim everything to left of window
                    break
            for window_start in reversed(range(n_pix - window_width)):
                mean_snr = np.sqrt(np.nanmean(self.ivars[r][n,window_start:window_start+window_width]))
                if mean_snr > min_snr:
                    self.ivars[r][n,window_start+window_width:] = 0. # trim everything to right of window
                    break

    #### 
    '''
    def continuum_normalize(self, plot_continuum=False, plot_dir='../results/', **kwargs):
        """Continuum-normalize all spectra using a polynomial fit. Takes kwargs of utils.fit_continuum"""
        for r in range(self.R):
            for n in range(self.N):
                try:
                    fit = fit_continuum(self.xs[r][n], self.ys[r][n], self.ivars[r][n], **kwargs)
                    if plot_continuum:
                        fig, ax = plt.subplots(1, 1, figsize=(8,5))
                        ax.scatter(self.xs[r][n], self.ys[r][n], marker=".", alpha=0.5, c='k', s=40)
                        mask = self.ivars[r][n] <= 1.e-8
                        ax.scatter(self.xs[r][n][mask], self.ys[r][n][mask], marker=".", alpha=1., c='white', s=20)                        
                        ax.plot(self.xs[r][n], fit)
                        fig.savefig(plot_dir+'continuum_o{0}_e{1}.png'.format(r, n))
                        plt.close(fig)
                    self.ys[r][n] -= fit
                except:
                    print("ERROR: Data: order {0}, epoch {1} could not be continuum normalized!".format(r,n))
    '''
    #### Edited function below to change plot dir
    def continuum_normalize(self, **kwargs):
        # passinng parameters via kwargs do data class doew not require changes if adapting this to a new wobble version 
        try:
            p = kwargs["parameters"]
        except:
            raise Exception("must pass parameters object to data object via kwargs in .Data for continuum_normalize")
        plot_continuum = p.plot_continuum
        if plot_continuum:
            plot_dir_continuum = p.plot_dir + "/continuum/"
            os.makedirs(plot_dir_continuum, exist_ok = True)
        order = p.continuum_order
        nsigma = p.continuum_nsigma
        
        #TODO make this not Hard Coded
        file_dir = os.path.dirname(__file__)
        telluric_mask_file = file_dir + "/"+"../../../wobble_aux/carmenes_aux_files/" + "telluric_mask_carm_short.dat"
        #mask tellurics
        if telluric_mask_file is not None:
            telluric_mask = np.genfromtxt(telluric_mask_file)
            #extend mask (with value 0) to include earliest CARMENES orders
            mask = telluric_mask
            mask = np.insert(mask, 0, [0.0,0],axis = 0)
            mask = np.append(mask, [[50000, 0]], axis = 0)
            #create decision function basedon interpolation of Mask
            mask_function = interpolate.interp1d(mask[:,0],mask[:,1])
            mask_array = np.zeros(np.array(self.xs).shape)
            mask_bool = np.ones(np.array(self.xs).shape, dtype = bool)
            #loop over all entries
            for o, xs_order in enumerate(self.xs):
                for e, xs_epoch in enumerate(xs_order):
                    for l, xs_lambda in enumerate(xs_epoch):
                        if mask_function(np.exp(xs_lambda)) == 1:
                            mask_array[o,e,l] = 1
                            mask_bool[o,e,l] = False
            #mask data (all parts with this shape) (xs, ys, ((fluxes, flux_ivars,)) ivars)
            xs_masked = np.ma.masked_array(self.xs, mask = mask_array)
            ys_masked = np.ma.masked_array(self.ys, mask = mask_array)
            ivars_masked = np.ma.masked_array(self.ivars, mask = mask_array)
            
            #xs_deleted = np.delete(self.xs, mask_bool, axis = 0)
            #ys_deleted = np.delete(self.ys, mask_bool, axis = 0)
            #ivars_deleted = np.delete(self.ivars, mask_bool, axis = 0)
            #print(xs_deleted, xs_deleted.shape)
            
            
            #NOTE Move masking to loop below. I suspect the mask breaks fit_continuum
            
        """Continuum-normalize all spectra using a polynomial fit. Takes kwargs of utils.fit_continuum"""
        #self.drop_orders = []
        for r in range(self.R):
            for n in range(self.N):
                #try:
                    
                #fit = fit_continuum(xs_deleted[r][n], ys_deleted[r][n], ivars_deleted[r][n], order = order
                                ##, **kwargs
                                #)
                try:
                    #fit = fit_continuum(xs_masked[r][n].compressed(), ys_masked[r][n].compressed(), ivars_masked[r][n].compressed(), order = order, nsigma = nsigma
                                        #, **kwargs
                                        #) #pass compressed arrays, to get rid of masked sections
                    fit = fit_continuum(self.xs[r][n][mask_bool[r,n]], self.ys[r][n][mask_bool[r,n]], self.ivars[r][n][mask_bool[r,n]]
                    ,order = order
                    ,nsigma = nsigma
                    #, **kwargs
                    )
                    
                    #fit = fit_continuum(xs_masked[r][n], ys_masked[r][n], ivars_masked[r][n], order = order
                                        ##, **kwargs
                                        #)
                    #fit = fit_continuum(self.xs[r][n], self.ys[r][n], self.ivars[r][n], order = order
                                        ##, **kwargs
                                        #)
                except Exception as err:
                    print("Continuum normalization of order {0} epoch {1} failed. Dropping order {0}".format(self.orders[r], self.epochs[n]))
                    self.drop_orders.append(self.orders[r]) #append to drop orders list
                    print("self.drop_orders: ", self.drop_orders)
                    #continue
                    break #in contrast to continue this should break the entireloop and go to  the next order instead of next epoch
                    
                    
                
                                    
                if plot_continuum:
                    
                    fig, ax = plt.subplots(1, 1, figsize=(8,5))
                    ax.scatter(self.xs[r][n], self.ys[r][n], marker=".", alpha=0.1, c='b', s=40)
                    ax.scatter(xs_masked[r][n], ys_masked[r][n], marker=".", alpha=0.5, c='k', s=40)
                    mask_ivars = self.ivars[r][n] <= 1.e-8
                    ax.scatter(self.xs[r][n][mask_ivars], self.ys[r][n][mask_ivars], marker=".", alpha=1., c='white', s=20)                        
                    ax.plot(xs_masked[r][n].compressed(), fit)
                    fig.savefig(plot_dir_continuum +'continuum_o{0}_e{1}.png'.format(self.orders[r], self.epochs[n]))
                    plt.close(fig)
                #self.ys[r][n] -= fit
                #fit needs to be interpolated to match  grid before mask
                try: #just for testing where it breaks
                    fit_function = interpolate.interp1d(xs_masked[r][n].compressed(), fit, fill_value = "extrapolate")#requires numpy 1.17 -> creates deprecation warnings
                    #HACK If this fails for whatever reason drop the order
                except:                 
                    print("dropping order {0} due to epoch {1} continuum interpolation error. xs_masked[r][n].compressed():".format(self.orders[r], self.epochs[n]))
                    print(xs_masked[r][n].compressed())
                    
                    self.drop_orders.append(self.orders[r]) #append to drop orders list
                    print("self.drop_orders: ", self.drop_orders)
                    break # jump to next order
                    
                    ##still throw error
                    #fit_function = interpolate.interp1d(xs_masked[r][n].compressed(), fit, fill_value = "extrapolate")#requires numpy 1.17 -> creates deprecation warnings
                    
                    
                self.ys[r][n] -= fit_function(self.xs[r][n])
                #except:
                    #print("ERROR: Data: order {0}, epoch {1} could not be continuum normalized!".format(self.orders[r],self.epochs[n]))
                    
    def append(self, data2):
        """Append another dataset to the current one(s)."""
        assert self.R == data2.R, "ERROR: Number of orders must be the same."
        for attr in ['dates', 'bervs', 'pipeline_rvs', 'pipeline_sigmas', 
                        'airms', 'drifts', 'filelist', 'origin_file']:
            setattr(self, attr, np.append(getattr(self, attr), getattr(data2, attr)))
        for attr in ['fluxes', 'xs', 'flux_ivars', 'ivars', 'ys']:
            attr1 = getattr(self, attr)
            attr2 = getattr(data2, attr)
            full_attr = [np.append(attr1[i], attr2[i], axis=0) for i in range(self.R)]
            setattr(self, attr, full_attr)
        self.epochs = [self.epochs, data2.epochs] # this is a hack that needs to be fixed
        self.epoch_groups.append((self.N + data2.epochs))
        self.N = self.N + data2.N
                
