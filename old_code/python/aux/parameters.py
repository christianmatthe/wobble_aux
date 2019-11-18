import sys
sys.path.append('/data/cmatthe/python/wobble_aux')#this  path only works on lx39
import combine_results as cr

import yaml
#TODO should probably use separate classes for regularization and optimize parameter files
class Parameters:
    """
    The Parameter dictonary object: contains the parameters specifying how to run wobble scripts. It is intended to be written into a .yaml file to be passed between scripts
    
    Parameters
    ----------
    starname : `str`
        Name of the star to be worked on.
    K_star : `int` (default `0`)
        Number of components used to model the stellar spectrum.
    K_t : `int` (default `3`)
        Number of components used to model the telluric spectrum.
    defaultQ : `bool` (default `True`)
        Determines whether default regularization parameters will be used.
    niter : `int` (default `160`)
        number of iterations to be used when optimizing.
    start : `int` (default `11`)
        number of first order that will be optimized.
    end : `int` (default `53`)
        number of first order that will *not* be optimized.
    chunk_size : `int` (default `5`)
        Size of the chunks of orders that will be optimized in one chunk level script. 
    output_suffix : `str` (default ``)
        string appended to output filenames. Serves as distinguishing name when rerunning the same dataset with different parameters
    data_suffix : `str` (default ``)
        string appended to data file to be imported. Allows for selecting different data files for one starname NOTE only compatible with optimize with optimize_top_3.1.2 other versions my be hardcoded
    
        
    
    """
    def __init__(self, 
                 filename = None,
                 starname = None,
                 K_star = 0,
                 K_t = 3,
                 defaultQ = True,
                 niter = 160,
                 #needed for chunked execution
                 start = 11,
                 end = 53,
                 chunk_size = 5,
                 
                 reg_file_star = None,
                 reg_file_t = None,
                 
                 output_suffix = '',
                 
                 data_suffix = ""
                 ):
        if filename is None:
            self.dictionary =  {
            #loaded by optimize_chunk
            #both
            "starname" : starname,
            "K_star" : K_star,
            "K_t" : K_t,
            "defaultQ" : defaultQ,
            "niter" : niter,
            #needed by optimize_top
            "start" : start,
            "end" : end,
            "chunk_size" : chunk_size,
            
            "output_suffix" : output_suffix,
            
            "data_suffix" : data_suffix
            #not currently recalled
            #"default_str" : default_str,
            #"results_file_base" : results_file_base,
            #"results_file" : results_file,
            #"data_file" : data_file,
            #"wobble_dir" : wobble_dir,
            #"wobble_dir_file_name" : wobble_dir_file_name,
            #"compare_file" : compare_file
            }
            if defaultQ is False:
                if reg_file_star is not None and reg_file_t is not None:
                    self.dictionary.update({"reg_file_star" : reg_file_star})
                    self.dictionary.update({"reg_file_t" : reg_file_t})
                else:
                    raise Exception("When  using non-default regularization must supply both reg_file_star and reg_file_t keywords.")
                    
            return
        
        if starname is None:
            self.dictionary = self.read(filename)
            return
        print("Results: must supply either starname or filename keywords.")
        
    def write(self, filename):
        """ Write to .yaml file"""
        with open(filename, 'w') as outfile:
            yaml.dump(self.dictionary, outfile)
                 
    def read(self, filename):
        """ read from to .yaml file"""
        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)
            return data_loaded
        
    def alt_regularization(self, regularization_filename):
        self.dictionary.update({"alt_regularization" : regularization_filename})
                 
if __name__ == "__main__":
    parameters = Parameters(starname = "GJ876", start = 0, end = 61, chunk_size = 15, defaultQ = False)
    parameters.write("/data/cmatthe/wobble_data/scripts/yaml_temp/optimize_parameters.yaml")
                 
