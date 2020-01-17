import matplotlib.pyplot as plt #for plt
from matplotlib.backends.backend_pdf import PdfPages # for PdfPages: suspiciously abscent in akaminskis code
#I suspect there is a hidden standart import for example also for numpy?
# perhaps in  sys.path.append('/home/akaminsk/pythontests')
import matplotlib as mpl #missing in akaminski
import h5py 
import numpy as np
import os

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
    loop_iterations = 6
    
    #name= "GJ436_orderwise_snr+drift_shift_1"
    #run_dir = '/data/cmatthe/wobble_reg_search/'+ name +'/'
    name= "GJ1148_all_orders"
    run_dir = '/data/cmatthe/wobble_aux/results/regularization/'+ name +'/'
    pdf_dir = '/data/cmatthe/python/wobble_aux/reg_plots/'
    
    start_order = 0
    
    pdf_name = pdf_dir + name + '_reg_plots.pdf'
    
    plot_title = name + '_regularization parameters'
    
    #plot_regularization(reg_file, pdf_name, start_order = start_order, plot_title = plot_title)
    pp =PdfPages(pdf_name)
    fig = plt.figure(figsize=(15, 9), dpi=200)
    mpl.rc('font', size=16)
    plt.clf()
    fig.clf()
    ax1=plt.gca()
    
    for i in range(loop_iterations):
        plot_title = "loop_{0}".format(i)+ '_star_regularization parameters'
        reg_file = run_dir + "loop_{0}/".format(i) + "next_base_star_reg" + '.hdf5'
        if os.path.isfile(reg_file):
        
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
    
    for i in range(loop_iterations):
        plot_title = "loop_{0}".format(i)+ '_t_regularization parameters'
        reg_file = run_dir + "loop_{0}/".format(i) + "next_base_t_reg" + '.hdf5'
        
        if os.path.isfile(reg_file):
        
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
