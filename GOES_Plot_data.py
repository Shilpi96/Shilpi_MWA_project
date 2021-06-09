import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.gridspec import GridSpec


import sunpy
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net import hek
from sunpy.time import TimeRange, parse_time
import sunpy.timeseries as ts


####### Downloading the Data

tr = TimeRange(['2014-09-28 02:00', '2014-09-28 04:00'])
results = Fido.search(a.Time(tr), a.Instrument.xrs)
files = Fido.fetch(results, path = 'path of the directory', overwrite = False)

####### Creating the timeseries

goes_ts = ts.TimeSeries(files, source="XRS")
tr = TimeRange("2014/09/28T02:25:00", "2014/09/28T04:00:00") #change this to what time range you want to plot
goes_tr = goes_ts.truncate(tr). #### need to truncate if goes have two timeseries of GOES-13 and GOES-15
date_format_goes = dates.DateFormatter("%H:%M")

fig = plt.figure()
gs = GridSpec(1,1)

ax1 = fig.add_subplot(gs[0])

gax = goes_tr.plot(ax1, legend=False, fontsize=16, rot=0)

gax.xaxis.set_major_formatter(date_format_goes)

gax.set_xlabel("Time (UTC)",fontsize=16)
gax.set_ylabel(r"Watts m$^{-2}$",fontsize=16)
gax.text(0.05,0.9,'Solar X-ray Flux', fontdict={'size':16}, transform=gax.transAxes)
gax.set_yscale("log")
gax.autoscale(axis="x", tight=True)
gax2 = ax1.twinx()
gax2.set_yscale("log")
gax.set_ylim(1e-9, 1e-3)
gax2.set_ylim(1e-9, 1e-3)
gax2.set_yticks((1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))#((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
gax2.set_yticklabels(('A', 'B', 'C', 'M', 'X', ' '))#((' ', 'A', 'B', 'C', 'M', 'X', ' '))
handles, labels = gax.get_legend_handles_labels()
handles.reverse()

labels = [r"1-8 $\AA$", r"0.5-4 $\AA$"]
gax.yaxis.grid(True, 'major')
gax.legend(handles, labels, fontsize=16)
plt.savefig('/home/shilpi/Desktop/Goes_plot.png', dpi=200, format='png')
plt.show()

