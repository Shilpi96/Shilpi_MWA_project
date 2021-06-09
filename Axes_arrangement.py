import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure(figsize=(10, 9))
outer = gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.03)
outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.05)
inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:,:1], wspace=0.1, hspace=0.1)
inner1 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:1,1:], wspace=0.001, hspace=0.01)           ##### from 0th row and 1st column  

inner2 = gridspec.GridSpecFromSubplotSpec(6, 1,
        subplot_spec=outer[1:3,1:], wspace=0.001, hspace=0.01)
for i in range(1):
	ax = plt.Subplot(fig, inner0[i])
	t = ax.text(0.5,0.5, 'Spectrogram %i' % (i))
	t.set_ha('center')
	ax.set_xticks([])
	ax.set_yticks([])
	fig.add_subplot(ax)

for i in range(1):
	ax = plt.Subplot(fig, inner1[i])
	t = ax.text(0.5,0.5, 'image %i' % (i))
	t.set_ha('center')
	ax.set_xticks([])
	ax.set_yticks([])
	fig.add_subplot(ax)

for i in range(6):
	ax = plt.Subplot(fig, inner2[i])
	t = ax.text(0.5,0.5, 'image %i' % (i))
	t.set_ha('center')
	ax.set_xticks([])
	ax.set_yticks([])
	fig.add_subplot(ax)


plt.show()
