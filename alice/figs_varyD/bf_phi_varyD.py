import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

efficiency=0.734
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(5,8))
fig=plt.figure(1)
height=0.306
width=0.79
x0=0.165
y0=0.072
icent=0

for ichargepair in range(1,4):
  #fig=plt.figure(ichargepair)
  if ichargepair==1:
    chargepair='allcharges'
  if ichargepair==2:
    chargepair='KK'
  if ichargepair==3:
    chargepair='pp'
    
  ax = fig.add_axes([x0+icent*width,y0+(ichargepair-1)*height,width,height])
    
  print('ichargepair=',ichargepair)

  centrality='rhic_cent0_5'
  #centrality='rhic_cent40_50'

  filename='../fromcascade_narrowphi/'+centrality+'/results_star/'+chargepair+'/bf_phi.dat'
  chargesdata = np.loadtxt(filename,skiprows=0,unpack=True)
  x=chargesdata[0]
  ycascade=chargesdata[1]
  
  filename='../fromcharges_narrowphi/'+centrality+'/results_star/'+chargepair+'/bf_phi.dat'
  chargesdata = np.loadtxt(filename,skiprows=0,unpack=True)
  x=chargesdata[0]
  ycharges=chargesdata[1]
  ysum=2.0*(ycascade+ycharges)
  
  filename='../fromcharges_D0.5/'+centrality+'/results_star/'+chargepair+'/bf_phi.dat'
  chargesdata = np.loadtxt(filename,skiprows=0,unpack=True)
  x_D0_5=chargesdata[0]
  ycharges=chargesdata[1]
  ysum_D0_5=2.0*(ycascade+ycharges)

  filename='../fromcharges_D2/'+centrality+'/results_star/'+chargepair+'/bf_phi.dat'
  chargesdata = np.loadtxt(filename,skiprows=0,unpack=True)
  x_2=chargesdata[0]
  ycharges=chargesdata[1]
  ysum_D2=2.0*(ycascade+ycharges)
    
  filename='../fromcharges_D4/'+centrality+'/results_star/'+chargepair+'/bf_phi.dat'
  chargesdata = np.loadtxt(filename,skiprows=0,unpack=True)
  x_4=chargesdata[0]
  ycharges=chargesdata[1]
  ysum_D4=2.0*(ycascade+ycharges)
    

  plt.plot(x,efficiency*ysum_D0_5*180.0/pi,linestyle='--',linewidth=2,color='r',label='$0.5D_{\\rm latt}$')
  plt.plot(x,efficiency*ysum*180.0/pi,linestyle='-',linewidth=2,color='k',label='$D_{\\rm latt}$')
  plt.plot(x,efficiency*ysum_D2*180.0/pi,linestyle='dotted',linewidth=3,color='g',label='$2D_{\\rm latt}$')
  plt.plot(x,efficiency*ysum_D4*180.0/pi,linestyle='-.',linewidth=3,color='b',label='$4D_{\\rm latt}$')
  

  #ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
  #ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
  ax.yaxis.set_major_formatter(sformatter)

  if(ichargepair==1):
    plt.xlabel('$\Delta\phi$ (degrees)', fontsize=18, weight='normal')
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xticks(np.arange(0,225,45), minor=False)
    ax.set_xticklabels(np.arange(0,225,45), minor=False, family='serif')
    ax.set_xticks(np.arange(0,225,15), minor=True)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
    plt.xlim(0,180)
    ax.set_yticks(np.arange(0,1.0,0.1), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.1), minor=False, family='serif')
    ax.set_yticks(np.arange(0,1.0,0.05), minor=True)
    plt.ylim(0.0,0.35)
    plt.text(5,0.02,'all charges',ha='left',fontsize=20)
  else:
    plt.xlabel('')
    ax.set_xticklabels([])
    plt.xlim(0.0,181)
  if(ichargepair==2):
    ax.set_xticks(np.arange(0,225,45), minor=False)
    ax.set_xticklabels([])
    ax.set_xticks(np.arange(0,225,15), minor=True)
    plt.xlim(0,180)
    plt.ylabel('$B(\Delta\phi)$',fontsize=20)
    ax.set_yticks(np.arange(0,0.2,0.01), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.01), minor=False, family='serif')
    ax.set_yticks(np.arange(0,1.0,0.005), minor=True)
    plt.ylim(0.0,0.0299)
    plt.text(5,0.002,'$K^+K^-$',ha='left',fontsize=20)
  else:
    plt.ylabel('')
  if(ichargepair==3):
    ax.set_xticks(np.arange(0,225,45), minor=False)
    ax.set_xticklabels([])
    ax.set_xticks(np.arange(0,225,15), minor=True)
    plt.xlim(0,180)
    ax.set_yticks(np.arange(0,1.0,0.01), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.01), minor=False, family='serif')
    ax.set_yticks(np.arange(0,1.0,0.005), minor=True)
    plt.ylim(0.0,0.055)
    plt.text(5,0.005,'$p\\bar{p}$',ha='left',fontsize=20)
    ax.legend()
  
  
  if ichargepair==1:
    stardata=np.loadtxt('stardata/AuAuPhiCent0_5.dat',skiprows=0,unpack=True)
    #stardata=np.loadtxt('stardata/AuAuPhiCent40_50.dat',skiprows=0,unpack=True)
    xstar=stardata[0]*180.0/pi
    ystar=stardata[1]
    #plt.plot(xstar,ystar,linestyle='-',linewidth=2,color='b')

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
    print('normstar=',normstar,' normmodel=',normmodel)
    print('widthstar=',widthstar,' widthmodel=',widthmodel)

                
    
    #plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
    #fontsize=12, color='gray')
    #plt.subplots_adjust(top=0.85)


plt.savefig('bf_phi_varyD.pdf',format='pdf')
os.system('open -a Preview bf_phi_varyD.pdf')
#plt.show()
quit()
