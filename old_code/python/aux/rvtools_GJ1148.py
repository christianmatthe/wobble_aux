
# coding: utf-8

# ## CAVEMAN rvtools 
#   - ### rvstat - General Statistics and Linear Fit
#   - ### rvgls - (Generalized) Lomb-Scargle Periodogram with Floating Mean (based on [astroML](http://www.astroml.org/modules/generated/astroML.time_series.lomb_scargle.html))
#   - ### rvkepfit - Keplerian fit up to 4 planets

# In[ ]:

from pycarm import carmdata
targets = carmdata.get_targets('J11417+427')
print(targets)


# In[ ]:

from pycarm.target import Target

# create target object
my_target = Target(targets[0]) # or explicitely my_target=Target('J04153-076')
print("Target name: '{}'".format(my_target.name))
print("The default data is type '{}' and channel '{}'\n".format(my_target.default_datatype, my_target.default_channel))

# information about a data type can be retrieved as follows:
avc_info = carmdata.get_datatype(my_target.default_datatype, print_results=True)

# default data can be accessed as follows
my_target.data

#HACK in Wobble data for GJ1148

wobble_res_lx39 = h5py.File("/data/cmatthe/compare_wobble_serval/wobbledir/results_"+i[0]+"_Kstar0_"+kt+".hdf5",'r')
wobble_res=wobble_res_lx39
#w_dates = wobble_res['dates'].value
w_dates = wobble_res['dates'][()]
#use selfcombined RVs a wobble RVs seem broken currently, actually maybe they are fine now (05.04.2019)
w_RVs = wobble_res['star_time_rvs'][()]
w_RVs_original = w_RVs
w_RVs_er = wobble_res['star_time_sigmas'][()]
w_bervs = wobble_res['bervs'][()]

print(my_target.data[1])

# see carmdata.ipynb to learn more how to work with carmdata package
'''

# In[ ]:

import matplotlib
get_ipython().magic('matplotlib inline')
#enable 2x images (bette quality of plot) by just adding the line
get_ipython().magic("config InlineBackend.figure_format = 'retina'")

from rvtools import RvStat 
stat = RvStat(my_target)
stat.print()

# one can get each statistic value, e.g. mean and median, as follows:
# print(stat.mean, stat.median)

# plot avc data and linear fit
stat.plot()


# In[ ]:

# GLS periodogram 
from rvtools import RvGls

# create/initiate gls object for my target  
gls = RvGls(my_target) 
# create periodogram
gls.rungls()
# plot resutls
gls.plot()
print("The strongest peak: " + str(gls.pgm_peaks[0]) + " at " + str(gls.periods_peaks[0]))
#print(gls.Nf,gls.fmin,gls.fmax,gls.deltaf,gls.Tbase)


# In[ ]:

# Keplerian fit
from pycarm.carmencita import CarmencitaReport
from pycarm import cavelog
from rvtools import RvKepFit
import logging


# set log level
cavelog.setLog(loglevel = [logging.INFO, logging.DEBUG])

# create a rvfit object for my target 
rvfit = RvKepFit(my_target)
# set up a star mass value taken from carmencita
cita = CarmencitaReport(my_target.name)
rvfit.setStMass(cita.stmass.mass.value)
# set up slope parameter to true
rvfit.setSlopeUse(True)             
# set up P(days) and K(m/s) initial parameters for fitting
p = gls.periods_peaks[0] # take the strongest peak from the gls above
k = (stat.max - stat.min)/2.
# add planet with p and k
rvfit.add_planet(p,k)
# run keplerian fitting with one planet 
rvfit.runRVmod()
rvfit.print()
rvfit.plot()


# In[ ]:

# show phase folded plot
rvfit.plot_phase_folded()


# In[ ]:

# run gls for residuals
gls.setData(x=rvfit.results_jd, y=rvfit.results_o_c, yerr=rvfit.results_rv_err)
gls.rungls()
gls.plot()
print("The strongest peak: " + str(gls.pgm_peaks[0]) + " at " + str(gls.periods_peaks[0]))


# In[ ]:

# add a second planet to fit with the following p,k
p = gls.periods_peaks[0]
k = rvfit.rms   
# be careful, if you run this cells twise a new planet will be added
# to re-run it please launch notebook from the initialization of the RvGls object at least
rvfit.add_planet(p,k)
rvfit.runRVmod()
rvfit.plot()
rvfit.print()


# In[ ]:

# show phase folded plots
rvfit.plot_phase_folded()

'''

