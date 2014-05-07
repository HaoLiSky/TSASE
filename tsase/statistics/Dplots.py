#!/usr/bin/env python
import pylab as pl
from scipy import stats
import time
import sys
import numpy as np
params = {'legend.fontsize':10}
pl.rcParams.update(params)

tick_label_fontsize=12
axis_label_fontsize=12
pl.rc('xtick', labelsize=tick_label_fontsize)
pl.rc(('xtick.major','xtick.minor'),  pad=10)
pl.rc('ytick', labelsize=tick_label_fontsize)

fig = pl.figure()
ax1 = fig.add_subplot(111)

###  READ DATA ###
filename = sys.argv[1]
data = np.loadtxt( filename, skiprows=1)

steps = data[:, 0]
dx = data[:, 1]
dy = data[:, 2]
dz = data[:, 3]
msd = data[:, 4]

### CALC GRAD ###
time = steps*0.001
#msd_mets=msd*0.0000000001*0.0000000001
#time = steps*0.000000000000001
msd_m, msd_b, r_value, p_value, std_err = stats.linregress(time, msd)
msd_line=msd_m*time+msd_b
print msd_m
## PLOT DATA ###
#ax1.plot(time, dx, 'b-', label="dx")
#ax1.plot(time, dy, 'g-', label="dy")
#ax1.plot(time, dz, 'k-', label="dz")

#ax2 = ax1.twinx()

ax1.plot(time, msd, 'r-', label="MSD ")
ax1.plot(time, msd_line, 'k-', label="MSD ")

### AXIS AND LABELS ###
ax1.set_xlabel("Time, ps", fontsize=axis_label_fontsize, verticalalignment='top')
#ax1.set_ylabel(")",fontsize=axis_label_fontsize)
#ax2.set_ylabel("",fontsize=axis_label_fontsize)
ax1.axis([0, 5000, 0.0, 5.0])

lines1, labels1 = ax1.get_legend_handles_labels()
ax1.legend(lines1, labels1 , loc="upper left", fancybox=True, shadow=True)

#pl.legend()
fig.savefig("msdLiplot_time.pdf")

pl.show()

