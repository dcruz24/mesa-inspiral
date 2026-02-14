from matplotlib import rcParams as rc
import scipy.integrate as spi
import numpy as np

lw1 = 3.6
lw2 = 2.4
lw3 = 2.4
ft1 = 30
ft2 = 30
lg = 26
tl1 = 12
tl2 = 8

rc.update({'figure.figsize': (12,9)})
rc.update({'axes.labelsize': ft1})
rc.update({'axes.titlesize': ft1})
rc.update({'axes.linewidth': lw2})
rc.update({'lines.linewidth': lw1})
rc.update({'axes.formatter.limits': (-4,4)})
rc.update({'lines.markersize': 20.})
rc.update({'lines.markeredgewidth': lw1})
rc.update({'markers.fillstyle': 'none'})
rc.update({'xtick.labelsize': ft2})
rc.update({'ytick.labelsize': ft2})
rc.update({'xtick.direction': 'out'})
rc.update({'ytick.direction': 'out'})
rc.update({'xtick.major.size': tl1})
rc.update({'ytick.major.size': tl1})
rc.update({'xtick.minor.size': tl2})
rc.update({'ytick.minor.size': tl2})
rc.update({'xtick.minor.width': lw3})
rc.update({'xtick.major.width': lw3})
rc.update({'xtick.minor.width': lw3})
rc.update({'ytick.major.width': lw3})
rc.update({'ytick.minor.width': lw3})
rc.update({'xtick.major.pad': '8'})
rc.update({'ytick.major.pad': '10'})
rc.update({'xtick.top': False})
rc.update({'ytick.right': False})
rc.update({'legend.fontsize': lg})
rc.update({'legend.numpoints': 1})
rc.update({'legend.frameon': False})
rc.update({'font.family':'STIXGeneral','mathtext.fontset':'stix'})

cols = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd',
        u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

eijk = np.zeros((3, 3, 3))
eijk[0,1,2] = eijk[1,2,0] = eijk[2,0,1] = 1
eijk[0,2,1] = eijk[2,1,0] = eijk[1,0,2] = -1
