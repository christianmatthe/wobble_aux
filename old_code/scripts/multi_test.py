#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 11:00:42 2019

@author: cmatthe
"""

import os    
import multiprocessing as mp                                                                   
from multiprocessing import Pool                                                
                                                                                
                                                                                
processes = ('script_51peg.py','script_barnards.py','script_HD189733.py')                                    
                                                  
                                                                                
def run_process(process):                                                             
    os.system('python {}'.format(process))                                       
                                                                                
                                                                                
pool = Pool(processes=3)                                                        
pool.map(run_process, processes) 

        
