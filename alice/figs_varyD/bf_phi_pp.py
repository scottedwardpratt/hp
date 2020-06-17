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
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

centrality='rhic_cent0_5'

chargesdata = np.loadtxt('fromcharges/'+centrality+'/results_star/pp/bf_phi.dat',skiprows=0,unpack=True)
x=chargesdata[0]
y1=chargesdata[1]
y2=y1
i=0
while(i<72):
  y2[i]=y1[71-i]
  i=i+1
y=(y1+y2)
y=efficiency*y*180.0/pi

chargesdata_D0_5 = np.loadtxt('fromcharges_D0.5/'+centrality+'/results_star/pp/bf_phi.dat',skiprows=0,unpack=True)
x_D4=chargesdata_D0_5[0]
y1=chargesdata_D0_5[1]
y2=y1
i=0
while(i<72):
  y2[i]=y1[71-i]
  i=i+1
y_D0_5=(y1+y2)
y_D0_5=efficiency*y_D0_5*180.0/pi

filename='fromcharges_D2/'+centrality+'/results_star/pp/bf_phi.dat'
chargesdata_D2 = np.loadtxt(filename,skiprows=0,unpack=True)
x_D4=chargesdata_D2[0]
y1=chargesdata_D2[1]
y2=y1
i=0
while(i<72):
  y2[i]=y1[71-i]
  i=i+1
y_D2=(y1+y2)
y_D2=efficiency*y_D2*180.0/pi

filename='fromcharges_D4/'+centrality+'/results_star/pp/bf_phi.dat'
print 'filename=',filename
chargesdata_D4 = np.loadtxt(filename,skiprows=0,unpack=True)
x_D4=chargesdata_D4[0]
y1=chargesdata_D4[1]
y2=y1
i=0
while(i<72):
  y2[i]=y1[71-i]
  i=i+1
y_D4=(y1+y2)
y_D4=efficiency*y_D4*180.0/pi

cascadedata = np.loadtxt('fromcascade_narrowphi/'+centrality+'/results_star/pp/bf_phi.dat',skiprows=0,unpack=True)
xc=cascadedata[0]
yc1=cascadedata[1]
yc2=yc1
i=0
while(i<72):
  yc2[i]=yc1[71-i]
  i=i+1
yc=(yc1+yc2)
yc=efficiency*yc*180.0/pi

ysum=(y+yc)
ysum_D0_5=y_D0_5+yc
ysum_D2=y_D2+yc
ysum_D4=y_D4+yc

stardata=np.loadtxt('stardata/AuAuPhiCent0_5.dat',skiprows=0,unpack=True)
#stardata=np.loadtxt('stardata/AuAuPhiCent40_50.dat',skiprows=0,unpack=True)
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
print('normstar=',normstar,' normmodel=',normmodel)
print('widthstar=',widthstar,' widthmodel=',widthmodel)

#plt.plot(x,y,linestyle='-',linewidth=3,color='r')
#plt.plot(x_D4,y_D4,linestyle='-',linewidth=3,color='r')
#plt.plot(-x,y,linestyle='-',linewidth=3,color='r')
#plt.plot(x,yc,linestyle='-',linewidth=3,color='g')
#plt.plot(-x,yc,linestyle='-',linewidth=3,color='g')
plt.plot(x,ysum,linestyle='-',linewidth=2,color='k')
plt.plot(x,ysum_D0_5,linestyle='-',linewidth=2,color='r')
plt.plot(x,ysum_D2,linestyle='-',linewidth=2,color='g')
plt.plot(x,ysum_D4,linestyle='-',linewidth=2,color='b')
#plt.plot(-x,ysum,linestyle='-',linewidth=3,color='b')
#plt.plot(xstar,ystar,markersize=8,linestyle='None',marker='*',markerfacecolor='r',markeredgecolor='r')

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,225,45), minor=False)
ax.set_xticklabels(np.arange(0,225,45), minor=False, family='serif')
ax.set_xticks(np.arange(0,225,15), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.xlim(0,181)

ax.set_yticks(np.arange(0,1.0,0.02), minor=False)
ax.set_yticklabels(np.arange(0,1.0,0.01), minor=False, family='serif')
ax.set_yticks(np.arange(0,1.0,0.01), minor=True)
plt.ylim(-0.005,0.06)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\Delta\phi$ (degrees)', fontsize=18, weight='normal')
plt.ylabel('$B(\Delta\phi)$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('bf_phi_pp.pdf',format='pdf')
os.system('open -a Preview bf_phi_pp.pdf')
#plt.show()
quit()
