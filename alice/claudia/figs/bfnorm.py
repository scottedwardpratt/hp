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
        'size'   : 18}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,4.8))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.145,0.83,0.82])

mydata = np.loadtxt('bfnorm.dat',skiprows=1,unpack=True)
T=mydata[0]
pipinorm=mydata[1];
Kpinorm=mydata[2];
KKnorm=mydata[3];
ppinorm=mydata[4];
pKnorm=mydata[5];
ppnorm=mydata[6];

plt.plot(T,pipinorm,linestyle='-',linewidth=3,color='r')
plt.plot(T,Kpinorm,linestyle='-',linewidth=3,color='purple')
plt.plot(T,KKnorm,linestyle='-',linewidth=3,color='g')
plt.plot(T,ppinorm,linestyle='-',linewidth=3,color='cyan')
plt.plot(T,pKnorm,linestyle='-',linewidth=3,color='orange')
plt.plot(T,ppnorm,linestyle='-',linewidth=3,color='b')


xTc=[155.0,155.0]
yTc=[-0.2,1.0]
plt.plot(xTc,yTc,linestyle='--',linewidth=2,color='k')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,200,10), minor=False)
ax.set_xticklabels(np.arange(0,200,10), minor=False, family='serif',size=18)
ax.set_xticks(np.arange(0,200,5), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(125,175)

ax.set_yticks(np.arange(-1,1.1,0.2), minor=False)
ax.set_yticklabels(np.arange(-1,1.1,0.2), minor=False, family='serif',size=18)
ax.set_yticks(np.arange(-1,1.1,0.05), minor=True)
plt.ylim(0.0,1.0)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%2f'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$T$ (MeV)', fontsize=22, weight='normal')
plt.ylabel('$Z_{h|h^\prime}$',fontsize=22,labelpad=-2)

text(172,0.9,"$\pi\pi$",fontsize=22,color='r',ha='right')
text(160,0.56,"$pp$",fontsize=22,color='b',ha='right')
text(174.5,0.462,"$KK$",fontsize=22,color='g',ha='right')
text(129,0.39,"$\pi K$",fontsize=22,color='purple',ha='left')
text(151,0.39,"$\pi p$",fontsize=22,color='cyan',ha='right')
text(172,0.04,"$Kp$",fontsize=22,color='orange',ha='right')
text(156,0.25,"$T_{\\rm interface}$",fontsize=24,color='k')

plt.arrow(128, 0.41, 0, 0.05, head_width=1, head_length=0.02, fc='purple', ec='purple', width=0.3,zorder=50)
#plt.arrow(173, 0.11, 0, -0.07, head_width=1, head_length=0.02, fc='orange', ec='orange', width=0.3, zorder=100)

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('bfnorm.pdf',format='pdf')
os.system('open -a Preview bfnorm.pdf')
#plt.show()
quit()
