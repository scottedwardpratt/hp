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

efficiency=1.0

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

centrality='alice_cent0_10'
chargepair='KK'
alicedata_filename='alicedata/BF_dy_KK_C0_10_ALICE.dat'

chargesdata = np.loadtxt('../D1/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
x=chargesdata[0]
y1=chargesdata[1]
y_D1=efficiency*y1

cascadedata = np.loadtxt('../fromcascade/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
xc=cascadedata[0]
yc1=cascadedata[1]
yc=efficiency*yc1

ysum_D1=y_D1+yc

chargesdatatest = np.loadtxt('../Dtest/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
xtest=chargesdatatest[0]
y1test=chargesdatatest[1]
y_Dtest=efficiency*y1test

cascadedatatest = np.loadtxt('../fromcascadetest/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
xctest=cascadedatatest[0]
yc1test=cascadedatatest[1]
yctest=efficiency*yc1test

ysum_Dtest=y_Dtest+yctest


#cdata=np.loadtxt('../acc_correction/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
#cfactor=cdata[2]*4
#ysum_D1=ysum_D1/cfactor
#ysum_D0_5=ysum_D0_5/cfactor
#ysum_D2=ysum_D2/cfactor
#ysum_D4=ysum_D4/cfactor

Dy=0.1
D1norm=Dtestnorm=0.0
for i in range(0,20):
  D1norm+=ysum_D1[i]*Dy
  Dtestnorm+=ysum_Dtest[i]*Dy
print('norms are: ',D1norm,', ',Dtestnorm)

alicedata=np.loadtxt(alicedata_filename,skiprows=0,unpack=True)
alice_y=alicedata[0]
alice_bf=alicedata[1]
alice_errors=alicedata[2]

plt.plot(x,ysum_D1,linestyle='-',linewidth=2,marker='o',color='r',label='$D_{\\rm latt}$')
plt.plot(xtest,ysum_Dtest,linestyle='-',linewidth=2,marker='o',color='k',label='perfect acceptance')
plt.errorbar(alice_y,2.0*alice_bf,alice_errors,linestyle='None',markersize='12',marker='*',color='k',label='$\cdot$ALICE (prel)')
#plt.plot(aliceTOT_phi,aliceTOT_bf,linestyle='None',marker='*',color='k',label='ALICE (prel)')
ax.legend()

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,2,0.5), minor=False)
ax.set_xticklabels(np.arange(0,2,0.5), minor=False, family='serif')
ax.set_xticks(np.arange(0,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,1.2)

ax.set_yticks(np.arange(0,2.0,0.2), minor=False)
ax.set_yticklabels(np.arange(0,2.0,0.2), minor=False, family='serif')
ax.set_yticks(np.arange(0,2.0,0.05), minor=True)
plt.ylim(0.0,0.95)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\Delta y$', fontsize=18, weight='normal')
plt.ylabel('$B(\Delta y)$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('bf_y_KK_test.pdf',format='pdf')
#os.system('open -a Preview bf_y_KK.pdf')
os.system('okular bf_y_KK_test.pdf &')
#plt.show()
quit()
