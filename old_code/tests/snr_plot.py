
import os
import wobble
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#import run_wobble as rw


snr_file = "/data/cmatthe/wobble_aux/data/servaldir/CARM_VIS/J02530+168/J02530+168.snr.dat"

snr_array = np.genfromtxt(snr_file)

#print(snr_array[36:38])
#print("shape: ", np.shape(snr_array[37]))

ep_lst = [i for i in range(30,41)]
snr_avg = [snr_array[i][1] for i in ep_lst]

fig = plt.figure(figsize=(15, 9), dpi=200)
mpl.rc('font', size=16)
plt.clf()
fig.clf()
ax1=plt.gca()
ax1.plot(ep_lst, snr_avg,"x", label="Overall_snr")

ax1.set_ylabel('snr')
ax1.set_xlabel('epoch')
plt.title("")
plt.grid(True)
plt.tight_layout()
plt.legend(shadow=True)
plt.savefig("plot.pdf", format='pdf')
plt.clf()
fig.clf()
ax1 = plt.gca()

order = 45
snr_avg = [snr_array[i][order +1] for i in ep_lst]

fig = plt.figure(figsize=(15, 9), dpi=200)
mpl.rc('font', size=16)
plt.clf()
fig.clf()
ax1=plt.gca()
ax1.plot(ep_lst, snr_avg,"x", label="snr")

ax1.set_ylabel('snr')
ax1.set_xlabel('epoch')
plt.title("Order {0} snr".format(order))
plt.grid(True)
plt.tight_layout()
plt.legend(shadow=True)
plt.savefig("plot_order_{0}.pdf".format(order), format='pdf')
plt.clf()
fig.clf()
ax1 = plt.gca()
