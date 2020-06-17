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

######
efficiency=0.75
cent='50_60'
######

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

chargesdata = np.loadtxt('fromcharges/rhic_cent'+cent+'/results_star/allcharges/bf_eta.dat',skiprows=0,unpack=True)
x=chargesdata[0]
y=chargesdata[1]*efficiency/(1.0-0.5*x)

cascadedata = np.loadtxt('fromcascade/rhic_cent'+cent+'/results_star/allcharges/bf_eta.dat',skiprows=0,unpack=True)
xc=cascadedata[0]
yc=cascadedata[1]*efficiency/(1.0-0.5*xc)

ysum=y+yc

stardata=np.loadtxt('stardata/AuAuCent'+cent+'.dat',skiprows=0,unpack=True)
xstar=stardata[0]
ystar=stardata[1]/(1.0-0.5*xstar)

normstar=0.0
normmodel=0.0
widthstar=0.0
widthmodel=0.0
for i in range(0,20):
  eta=(i+0.5)*0.1
  normstar+=0.1*ystar[i]
  widthstar+=0.1*ystar[i]*eta
  normmodel+=0.1*ysum[i]
  widthmodel+=0.1*ysum[i]*eta

widthstar=widthstar/normstar
widthmodel=widthmodel/normmodel
print('normstar=',normstar,' normmodel=',normmodel)
print('widthstar=',widthstar,' widthmodel=',widthmodel)


plt.plot(x,y,linestyle='-',linewidth=3,color='r')
#plt.plot(-x,y,linestyle='-',linewidth=3,color='r')
plt.plot(x,yc,linestyle='-',linewidth=3,color='g')
#plt.plot(-x,yc,linestyle='-',linewidth=3,color='g')
plt.plot(x,ysum,linestyle='-',linewidth=3,color='b')
#plt.plot(-x,ysum,linestyle='-',linewidth=3,color='b')
plt.plot(xstar,ystar,linestyle='None',markersize=8,marker='*',markerfacecolor='r',markeredgecolor='r')


#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,3,0.5), minor=False)
ax.set_xticklabels(np.arange(0,3,0.5), minor=False, family='serif')
ax.set_xticks(np.arange(0,3,0.3), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,2)

ax.set_yticks(np.arange(0,1.0,0.2), minor=False)
ax.set_yticklabels(np.arange(0,1.0,0.2), minor=False, family='serif')
ax.set_yticks(np.arange(0,1.0,0.05), minor=True)
plt.ylim(0.0,0.6)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\Delta\eta$', fontsize=18, weight='normal')
plt.ylabel('$B(\Delta\eta)$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('bf_eta.pdf',format='pdf')
os.system('open -a Preview bf_eta.pdf')
#plt.show()
quit()
