import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(False)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(5,9))
fig = plt.figure(1)

npanels=6
for ipanel in range(0,npanels):
  if ipanel==0:
    cent='0_5'
    description='(l) 0-5%'
    efficiency=0.7
  if ipanel==1:
    cent='10_20'
    description='(k) 10-20%'
    efficiency=0.72
  if ipanel==2:
    cent='20_30'
    description='(j) 20-30%'
    efficiency=0.74
  if ipanel==3:
    cent='30_40'
    description='(i) 30-40%'
    efficiency=0.76
  if ipanel==4:
    cent='40_50'
    description='(h) 40-50%'
    efficiency=0.78
  if ipanel==5:
    cent='50_60'
    description='(g) 50-60%'
    efficiency=0.8
  
  bottom=0.06
  top=0.98
  left=0.14
  right=0.96
  
  ax = fig.add_axes([left,bottom+ipanel*(top-bottom)/npanels,right-left,(top-bottom)/npanels])
  ax.tick_params(axis='both', which='major', labelsize=12)

  chargesdata = np.loadtxt('fromcharges/rhic_cent'+cent+'/results_star/allcharges/bf_phi.dat',skiprows=0,unpack=True)
  x=chargesdata[0]
  y1=chargesdata[1]
  y2=y1
  i=0
  while(i<72):
    y2[i]=y1[71-i]
    i=i+1
  y=(y1+y2)
  y=efficiency*y*180.0/pi

  cascadedata = np.loadtxt('fromcascade/rhic_cent'+cent+'/results_star/allcharges/bf_phi.dat',skiprows=0,unpack=True)
  xc=cascadedata[0]
  yc1=cascadedata[1]
  yc2=yc1
  i=0
  while(i<72):
    yc2[i]=yc1[71-i]
    i=i+1
  yc=(yc1+yc2)
  yc=efficiency*yc*180.0/pi

  ysum=y+yc

  stardata=np.loadtxt('stardata/AuAuPhiCent'+cent+'.dat',skiprows=0,unpack=True)
  xstar=stardata[0]*180.0/pi
  ystar=stardata[1]

  normstar=0.0
  normmodel=0.0
  widthstar=0.0
  widthmodel=0.0
  for i in range(0,20):
    phimodel=(i+0.5)*(4.0*pi/180.0)
    normmodel+=(4.0*pi/180.0)*ysum[i]
    widthmodel+=(4.0*pi/180.0)*ysum[i]*phimodel

  for i in range(0,24):
    phistar=(i+0.5)*(7.5*pi/180.0)
    normstar+=(4.0*pi/180.0)*ystar[i]
    widthstar+=(4.0*pi/180.0)*ystar[i]*phistar
  
  widthstar=widthstar/normstar
  widthmodel=widthmodel/normmodel
  print('cent=',cent,' normstar=',normstar,' normmodel=',normmodel)
  print('----- widthstar=',widthstar,' widthmodel=',widthmodel)


  plt.plot(x,y,linestyle='--',linewidth=3,color='g')
  #plt.plot(-x,y,linestyle='-',linewidth=3,color='r')
  plt.plot(x,yc,linestyle=':',linewidth=4,color='g')
  #plt.plot(-x,yc,linestyle='-',linewidth=3,color='g')
  plt.plot(x,ysum,linestyle='-',linewidth=3,color='b')
  #plt.plot(-x,ysum,linestyle='-',linewidth=3,color='b')
  plt.plot(xstar,ystar,linestyle='',markersize=8,marker='*',markerfacecolor='r',markeredgecolor='r')


  #plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

  #plt.semilogy(x,y)

  ax.set_xticks(np.arange(0,200,45), minor=False)
  if ipanel==0:
    ax.set_xticklabels(np.arange(0,200,45), minor=False, family='serif')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
  else:
    ax.set_xticklabels([])
  ax.set_xticks(np.arange(0,200,15), minor=True)
  plt.xlim(0.0,180.0)

  ax.set_yticks(np.arange(0,1.0,0.1), minor=False)
  ax.set_yticklabels(np.arange(0.0,0.325,0.1), minor=False, family='serif')
  if ipanel==5:
    ax.set_yticklabels(np.arange(0,0.6,0.1), minor=False, family='serif')
  ax.set_yticks(np.arange(0,1.0,0.02), minor=True)
  plt.ylim(0.0,0.4)
  #ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
  #ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
  #ax.yaxis.set_major_formatter(sformatter)
  
  ax.text(172,0.325,description,fontsize=16,ha='right')

  if(ipanel==0):
    plt.xlabel('$\Delta\phi$ (deg)', fontsize=18, weight='normal')
  if(ipanel==3):
    plt.ylabel('$B(\Delta\phi)$',position=(0,-0.06),fontsize=20)

    
#end loop

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('bf_phi_allcents.pdf',format='pdf')
os.system('open -a Preview bf_phi_allcents.pdf')
#plt.show()
quit()
