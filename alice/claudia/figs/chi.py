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
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.19,0.14,0.78,0.82])

mydata = np.loadtxt('chi.dat',skiprows=1,unpack=True)
T=mydata[0]
chiuu=mydata[1]
chiud=mydata[2]
chius=mydata[3]
chiss=mydata[4]
mydatah = np.loadtxt('chihadron.dat',skiprows=1,unpack=True)
T_h=mydatah[0]
chiuu_h=mydatah[1]
chiud_h=mydatah[2]
chius_h=mydatah[3]
chiss_h=mydatah[9]
mydatal = np.loadtxt('chilattice.dat',skiprows=1,unpack=True)
T_l=mydatal[0]
chiuu_l=mydatal[1]
chiud_l=mydatal[2]
chius_l=mydatal[3]
chiss_l=mydatal[4]

chiBB=(2.0*chiuu+2.0*chiud+chiss+4.0*chius)/9.0;

plt.plot(T,chiuu,linestyle='-',linewidth=2,color='r')
plt.plot(T_h,chiuu_h,color='r',linestyle="None",markersize=6, marker='s',markevery=2)
plt.plot(T_l,chiuu_l,color='r',linestyle="None",markersize=6, marker='o',markevery=2)

plt.plot(T,chiud,linestyle='-',linewidth=2,color='cyan')
plt.plot(T_h,chiud_h,linestyle="None",color='cyan',markersize=6, marker='s',markevery=2)
plt.plot(T_l,chiud_l,linestyle="None",color='cyan',markersize=6, marker='o',markevery=2)

plt.plot(T,chius,linestyle='-',linewidth=2,color='b')
plt.plot(T_h,chius_h,linestyle="None",color='b',markersize=6, marker='s',markevery=2)
plt.plot(T_l,chius_l,linestyle="None",color='b',markersize=6, marker='o',markevery=2)

plt.plot(T,chiss,linestyle='-',linewidth=2,color='g')
plt.plot(T_h,chiss_h,linestyle='None',linewidth=2,color='g',markersize=6, marker='s',markevery=2)
plt.plot(T_l,chiss_l,linestyle='None',linewidth=2,color='g',markersize=6, marker='o',markevery=2)

#plt.plot(T,chiBB*3.0,linestyle='-',color='k')

xTc=[155.0,155.0]
yTc=[-0.2,0.2]
plt.plot(xTc,yTc,linestyle='--',linewidth=2,color='k')

Tqgp=[200.0,400]
#sratio_qgas=[0.0816021,0.0816021]
sratio_qgp=[0.0590912,0.0590912]
plt.plot(Tqgp,sratio_qgp,linestyle='--',linewidth=2,color='k')
zero_qgp=[0.0,0.0]
Tzero=[0.0,400]
plt.plot(Tzero,zero_qgp,linestyle='--',linewidth=1,color='k')

#plt.plot(T_l,chiuu_l-chiss_l+chiud_l,linestyle='-',color='k')
#plt.plot(T_h,chiuu_h-chiss_h+chiud_h,linestyle='-',color='k')

#plt.plot(T,chiud,linestyle='-',linewidth=2,color='y',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(T,chius,linestyle='-',linewidth=2,color='cyan',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(T,chiss,linestyle='-',linewidth=2,color='g',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,500,100), minor=False)
ax.set_xticklabels(np.arange(0,500,100), minor=False, family='serif', size='18')
ax.set_xticks(np.arange(0,500,50), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(50,405)

ax.set_yticks(np.arange(-1,1,0.05), minor=False)
ax.set_yticklabels(np.arange(-1,1,0.05), minor=False, family='serif',size='18')
ax.set_yticks(np.arange(-1,1,0.01), minor=True)
plt.ylim(-0.15,0.15)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%2f'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$T$ (MeV)', fontsize=22, weight='normal')
plt.ylabel('$\chi_{ab}/s$',fontsize=22,labelpad=-5)

text(55,-0.105,"$\chi_{ud}$",fontsize=22,color='cyan')
text(70,-0.03,"$\chi_{us}$",fontsize=22,color='b')
text(70,0.028,"$\chi_{ss}$",fontsize=22,color='green')
text(55,0.098,"$\chi_{uu}$",fontsize=22,color='red')
text(158,-0.125,"$T_{\\rm interface}$",fontsize=24,color='k')
text(275,0.04,"parton gas",fontsize=22,color='k')
text(370,0.13,"(a)",fontsize=22,color='k')

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('chi.pdf',format='pdf')
os.system('open -a Preview chi.pdf')
#plt.show()
quit()
