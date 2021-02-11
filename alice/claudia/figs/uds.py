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

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.17,0.12,0.8,0.8])

mydata = np.loadtxt('udsdens.dat',skiprows=1,unpack=True)
T=mydata[0]
udratio=mydata[1]
sratio=mydata[2]
#sratio_qgas=mydata[3]
#sratio_qgp=mydata[4]

plt.plot(T,0.5*udratio,linestyle='-',linewidth=3,color='r')
plt.plot(T,sratio,linestyle='-',linewidth=3,color='g')
Tqgp=[150.0,250]
sratio_qgas=[0.0816021,0.0816021]
sratio_qgp=[0.0590912,0.0590912]

xTc=[150.0,150.0]
yTc=[-0.2,0.2]
plt.plot(xTc,yTc,linestyle='--',linewidth=2,color='k')

#plt.plot(Tqgp,sratio_qgas,color='b',linestyle="-",linewidth=3)
plt.plot(Tqgp,sratio_qgp,color='b',linestyle="--",linewidth=3)

#plt.plot(T_l,chiuu_l-chiss_l+chiud_l,linestyle='-',color='k')
#plt.plot(T_h,chiuu_h-chiss_h+chiud_h,linestyle='-',color='k')

#plt.plot(T,chiud,linestyle='-',linewidth=2,color='y',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(T,chius,linestyle='-',linewidth=2,color='cyan',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(T,chiss,linestyle='-',linewidth=2,color='g',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,500,50), minor=False)
ax.set_xticklabels(np.arange(0,500,50), minor=False, family='serif')
ax.set_xticks(np.arange(0,500,25), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(50,250)

ax.set_yticks(np.arange(-1,1,0.05), minor=False)
ax.set_yticklabels(np.arange(-1,1,0.05), minor=False, family='serif')
ax.set_yticks(np.arange(-1,1,0.01), minor=True)
plt.ylim(0.0,0.2)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%2f'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$T$ (MeV)', fontsize=18, weight='normal')
plt.ylabel('$density/s$',fontsize=18)

text(100,0.035,"$n_{\\rm s}$",fontsize=20,color='g')
text(90,0.135,"$n_{\\rm u}=n_{\\rm d}$",fontsize=20,color='r')
text(180,0.065,"$n_{\\rm s}=n_{\\rm u}=n_{\\rm d}$",fontsize=20,color='b')
text(55,0.095,"hadron gas",fontsize=20,color='k')
text(185,0.045,"parton gas",fontsize=20,color='k')
text(152,0.005,"$T_{\\rm interface}$",fontsize=18,color='k')

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('uds.pdf',format='pdf')
os.system('open -a Preview uds.pdf')
#plt.show()
quit()
