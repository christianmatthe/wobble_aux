#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 11:13:04 2019

@author: cmatthe
"""
import numpy as np
from tabulate import tabulate

#f = open("/home/cmatthe/Documents/python/trifon-master/datafiles/hip5364_crires.vels","r+")

data = np.genfromtxt("/home/cmatthe/Documents/python/trifon-master/datafiles/hip5364_crires.vels")
print(data)


f = open('table.txt', 'w')
f.write(tabulate(data))
f.close()
     

data_re = np.genfromtxt("/home/cmatthe/Documents/python/trifon-master/datafiles/hip5364_crires.vels")
print(data_re)