import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
import codecs
sys.getdefaultencoding()
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-3,3))

acceptance='alice'

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}
plt.rc('font', **font)
plt.rc('text', usetex=False)
linestyles = [
  '_', '-', '--', ':'
]
mstyle = [u'*',u'o']
colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

plt.figure(figsize=(6,15))
fig = plt.figure(1)
Npanels=6
xx0=0.19
yy0=0.05
ww0=1.0-xx0-0.03
hh0=1.0-yy0-0.01
xmin=0
xmax=2.0

for ipanel in range (0,Npanels):
  ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
  plt.plot([0,5],[0,0],linestyle='dashed',color='grey')
  xlim(0,2.0) 
    
  if ipanel == 0:
    ptype='pipi'
    ymin=-0.04
    ymax=0.9
    ax.set_yticks(np.arange(-1,1,0.2), minor=False)
    ax.set_yticklabels(np.arange(-1,1,0.2), minor=False, family='serif', size=18)
    ax.set_yticks(np.arange(-1,1,0.1), minor=True)
    ax.annotate('(f) $\pi|\pi$',xy=(1.9,ymax-0.48*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')

  
  if ipanel == 1:
    ptype='piK'
    ymin=-0.004
    ymax=0.09
    ax.set_yticks(np.arange(-1,1,0.04), minor=False)
    ax.set_yticklabels(np.arange(-1,1,0.04), minor=False, family='serif', size=18)
    ax.set_yticks(np.arange(-1,1,0.02), minor=True)
    ax.annotate('(e) $\pi|K$',xy=(1.9,ymax-0.18*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')

  if ipanel == 2:
    ptype='pip'
    ymin=-0.002
    ymax=0.045
    ax.set_yticks(np.arange(0.0,.05,0.02), minor=False)
    ax.set_yticklabels(np.arange(0.0,.05,0.02), minor=False, family='serif', size=18)
    ax.set_yticks(np.arange(0.0,.05,0.01), minor=True)
    ax.annotate('(d) $\pi|p$',xy=(1.9,ymax-0.18*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
    
  if ipanel == 3:
    ptype='KK'
    ymin=-0.02
    ymax=0.45
    ax.set_yticks(np.arange(-1,1,0.2), minor=False)
    ax.set_yticklabels(np.arange(-1,1,0.2), minor=False, family='serif', size=18)
    ax.set_yticks(np.arange(-1,1,0.1), minor=True)
    ax.annotate('(c) $K|K$',xy=(1.9,ymax-0.18*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')

  if ipanel == 4:
    ptype='Kp'
    ymin=-0.005
    ymax=0.11
    ax.set_yticks(np.arange(-0.04,1,0.04), minor=False)
    ax.set_yticklabels(np.arange(-0.04,1,0.04), minor=False, family='serif', size=18)
    ax.set_yticks(np.arange(-0.04,1,0.02), minor=True)
    ax.annotate('(b) $K|p$',xy=(1.9,ymax-0.18*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')

  if ipanel == 5:
    ptype='pp'
    ymin=-0.025
    ymax=0.25
    ax.set_yticks(np.arange(-0.999999,1,0.1), minor=False)
    ax.set_yticklabels(np.arange(-1,1,0.1), minor=False, family='serif', size=18)
    ax.set_yticks(np.arange(-1,1,0.05), minor=True)
    ax.annotate('(a) $p|p$',xy=(1.9,ymax-0.18*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
    
  plt.ylim(ymin,ymax)
  plt.xlim(0,2.0)
  
  
  #ipanel = 0(pipi), 1(piK), 2(pip), 3(KK), 4(Kp), 5(pp)  
  for idata in range(0,2):
    #if idata == 0:
    #datadir='../exp_data_corrected/star_cent05/'
    # datadir='../exp_data/star_cent05/'
    if idata == 1:
      if ipanel == 0:
        centrality='alice_cent0_5'
      if ipanel == 1 or ipanel == 3:
        centrality='alice_cent0_10'
      if ipanel == 2 or ipanel == 4 or ipanel == 5:
        centrality='alice_cent0_20'
      filename='../D1_sigma1/model_output/default_sum/'+centrality+'/results_'+acceptance+'/'+ptype+'/bf_y.dat'
      data_fromcharges = np.loadtxt(filename,skiprows=0,unpack=True)
      filename='../fromcascade/model_output/default_sum/'+centrality+'/results_'+acceptance+'/'+ptype+'/bf_y.dat'
      data_fromcascade = np.loadtxt(filename,skiprows=0,unpack=True)
      dely=data_fromcharges[0]
      B=data_fromcharges[1]+data_fromcascade[1]
      plt.plot(dely,B,linestyle=linestyles[1],linewidth=3,
      color='k',markersize=6, marker=mstyle[idata],
      markerfacecolor='none',markeredgewidth=2,markeredgecolor='k')
      plt.plot(dely,data_fromcharges[1],linestyle=linestyles[1],linewidth=1,
      color='r',markersize=6, marker=mstyle[idata],
      markerfacecolor='none',markeredgewidth=2,markeredgecolor='r')
      plt.plot(dely,data_fromcascade[1],linestyle=linestyles[1],linewidth=1,
      color='g',markersize=6, marker=mstyle[idata],
      markerfacecolor='none',markeredgewidth=2,markeredgecolor='g')
    if idata == 0:
      if ipanel == 0:
        filename = 'alicedata/BF_dy_pipi_C0_5_ALICE.dat'
        expdata = np.loadtxt(filename,skiprows=0,unpack=True)
        expdely=expdata[0]
        expB=expdata[1]
        dB=expdata[2]
        plt.errorbar(expdely,2.0*expB,yerr=dB,linestyle='None',linewidth=1,color=colors[0],markersize=8, marker='*',markerfacecolor='none',markeredgewidth=2,markeredgecolor=colors[0],label='ALICE-Preliminary')
        plt.legend(loc='upper right',fontsize=24)
      if ipanel == 3:
        filename = 'alicedata/BF_dy_KK_C0_10_ALICE.dat'
        expdata = np.loadtxt(filename,skiprows=0,unpack=True)
        expdely=expdata[0]
        expB=expdata[1]
        dB=expdata[2]
        plt.errorbar(expdely,2.0*expB,yerr=dB,linestyle='None',linewidth=1,color=colors[0],markersize=8, marker='*',
        markerfacecolor='none',markeredgewidth=2,markeredgecolor=colors[0])

    if ipanel == 0 and idata==1:
      plt.xlabel('$\Delta y$',fontsize=24, family='serif')
      plt.ylabel('$B(\Delta y)$',fontsize=24,y=3.0, family='serif', labelpad=16)

    if ipanel!=0:
      ax.set_xticklabels([])
      
    if idata==1:
      ax.set_xticks(np.arange(0,4,0.5), minor=False)
      ax.set_xticks(np.arange(0,4,0.1), minor=True)
    
    if idata==1 and ipanel==0:
      ax.set_xticklabels(np.arange(0,4,0.5),minor=False,size=18, family='serif')
    else:
      ax.set_xticklabels([])
    
    xlim(0,xmax)
    
  ax.xaxis.set_major_formatter(sformatter)
  ax.yaxis.set_major_formatter(sformatter)
  if ipanel!=0:
    ax.set_xticklabels([])  

filename='general_y.pdf'
plt.savefig(filename,format='pdf')
os.system('open -a Preview '+filename)
#plt.show()

quit()
