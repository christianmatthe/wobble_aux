import matplotlib.pyplot as plt #for plt
from matplotlib.backends.backend_pdf import PdfPages # for PdfPages: suspiciously abscent in akaminskis code
#I suspect there is a hidden standart import for example also for numpy?
# perhaps in  sys.path.append('/home/akaminsk/pythontests')
import matplotlib as mpl #missing in akaminski
import h5py 
import numpy as np

def plot_regularization(reg_file, pdf_name, start_order = 0, plot_title = ""):
    pp =PdfPages(pdf_name)
    fig = plt.figure(figsize=(15, 9), dpi=200)
    mpl.rc('font', size=16)
    plt.clf()
    fig.clf()
    ax1=plt.gca()
    
    f = h5py.File(reg_file, 'r')
    
    for key in list(f.keys()):
        ax1.plot(np.asarray(np.arange(start_order, len(f[key][()]) + start_order)), np.log10(f[key][()]), "-o", label = key)
        
    ax1.set_ylabel('regularization_amplitude')
    ax1.set_xlabel('order')
    plt.title(plot_title)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    #plt.close(fig)
    #pp.close()
        
if __name__ == "__main__":
    #reg_dir = '/data/cmatthe/wobble_data/wobble/regularization/'
    reg_dir = '/data/cmatthe/wobble_aux/wobble_aux/regularization/'
    pdf_dir = '/data/cmatthe/python/wobble_aux/reg_plots/'
    
    #name = 'GJ1148_t_K3_orders[11,53)_stitched'
    start_order = 11
    
    #name = 'def_chunk_5_roll1_star'
    #name= "Wolf294_star_K0_orders[11,54)_regtest1406"
    #name= "GJ876_star_K0_orders[11,53)_stitched_reformatted"
    #name= "GJ1148_star_K0_orders[11,53)_stitched_reformatted"
    name= "dummy_star_K0_no_reg"
    start_order = 0
    
    reg_file = reg_dir + name + '.hdf5'
    pdf_name = pdf_dir + name + '_reg_plots.pdf'
    
    plot_title = name + '_regularization parameters'
    
    #plot_regularization(reg_file, pdf_name, start_order = start_order, plot_title = plot_title)
    pp =PdfPages(pdf_name)
    fig = plt.figure(figsize=(15, 9), dpi=200)
    mpl.rc('font', size=16)
    plt.clf()
    fig.clf()
    ax1=plt.gca()
    
    f = h5py.File(reg_file, 'r')
    
    for key in list(f.keys()):
        ax1.plot(np.asarray(np.arange(start_order, len(f[key][()]) + start_order)), np.log10(f[key][()]), "-o", label = key)
        
    ax1.set_ylabel('regularization_amplitude')
    ax1.set_xlabel('order')
    plt.title(plot_title)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()
    
    
    #TODO clean this up (2 plots in one pdfpages
    #name = 'def_chunk_5_roll1_t'
    #name= "Wolf294_t_K3_orders[11,54)_regtest1406"
    #name= "GJ876_t_K3_orders[11,53)_stitched_reformatted"
    #name= "GJ1148_t_K3_orders[11,53)_stitched_reformatted"
    name= "dummy_t_K3_no_reg"
    
    reg_file = reg_dir + name + '.hdf5'
    plot_title = name + '_regularization parameters'
    
    f = h5py.File(reg_file, 'r')
    
    for key in list(f.keys()):
        ax1.plot(np.asarray(np.arange(start_order, len(f[key][()]) + start_order)), np.log10(f[key][()]), "-o", label = key)
        
    ax1.set_ylabel('regularization_amplitude')
    ax1.set_xlabel('order')
    plt.title(plot_title)
    plt.grid(True)
    plt.tight_layout()
    plt.legend(shadow=True)
    plt.savefig(pp, format='pdf')
    plt.clf()
    fig.clf()
    ax1 = plt.gca()

    plt.close(fig)
    pp.close()
