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
ax = fig.add_axes([0.17,0.14,0.8,0.82])

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
chiBB_h=(2.0*chiuu_h+2.0*chiud_h+chiss_h+4.0*chius_h)/9.0;
chiBB_l=(2.0*chiuu_l+2.0*chiud_l+chiss_l+4.0*chius_l)/9.0;

chiQQ=(5.0*chiuu-4.0*chiud+chiss-2.0*chius)/9.0;
chiQQ_h=(5.0*chiuu_h-4.0*chiud_h+chiss_h-2.0*chius_h)/9.0;
chiQQ_l=(5.0*chiuu_l-4.0*chiud_l+chiss_l-2.0*chius_l)/9.0;

plt.plot(T,chiBB,linestyle='-',color='g')
plt.plot(T_l,chiBB_l,linestyle="None",color='g',markersize=6, marker='o',markevery=2)
plt.plot(T_h,chiBB_h,linestyle="None",color='g',markersize=6, marker='s',markevery=2)

plt.plot(T,chiQQ,linestyle='-',color='r')
plt.plot(T_l,chiQQ_l,linestyle="None",color='r',markersize=6, marker='o',markevery=2)
plt.plot(T_h,chiQQ_h,linestyle="None",color='r',markersize=6, marker='s',markevery=2)

xTc=[155.0,155.0]
yTc=[-0.2,0.2]
plt.plot(xTc,yTc,linestyle='--',linewidth=2,color='k')

Tqgp=[200.0,400]
#sratio_qgas=[0.0816021,0.0816021]
sratio_qgp=[0.0590912/3.0,0.0590912/3.0]
plt.plot(Tqgp,sratio_qgp,linestyle='--',linewidth=2,color='k')
sratio_qgp=[2.0*0.0590912/3.0,2.0*0.0590912/3.0]
plt.plot(Tqgp,sratio_qgp,linestyle='--',linewidth=2,color='k')

#plt.plot(T_l,chiuu_l-chiss_l+chiud_l,linestyle='-',color='k')
#plt.plot(T_h,chiuu_h-chiss_h+chiud_h,linestyle='-',color='k')

#plt.plot(T,chiud,linestyle='-',linewidth=2,color='y',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(T,chius,linestyle='-',linewidth=2,color='cyan',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(T,chiss,linestyle='-',linewidth=2,color='g',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=18)

ax.set_xticks(np.arange(0,500,50), minor=False)
ax.set_xticklabels(np.arange(0,500,50), minor=False, family='serif',size=18)
ax.set_xticks(np.arange(0,500,25), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(50,405)

ax.set_yticks(np.arange(-1,1,0.03), minor=False)
ax.set_yticklabels(np.arange(-1,1,0.03), minor=False, family='serif',size=18)
ax.set_yticks(np.arange(-1,1,0.01), minor=True)
plt.ylim(0.0,0.125)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
#ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$T$ (MeV)', fontsize=22, weight='normal')
plt.ylabel('$\chi_{BB}/s$, $\chi_{QQ}/s$',fontsize=22,labelpad=1)

text(70,0.007,"$\chi_{BB}$",fontsize=22,color='g')
text(70,0.1,"$\chi_{QQ}$",fontsize=22,color='r')
text(158,0.005,"$T_{\\rm interface}$",fontsize=24,color='k')
text(280,0.03,"parton gas",fontsize=22,color='k')
text(370,0.115,"(b)",fontsize=22,color='k')

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('chiBBQQ.pdf',format='pdf')
os.system('open -a Preview chiBBQQ.pdf')
#plt.show()
quit()
