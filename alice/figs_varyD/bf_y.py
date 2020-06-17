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


chargepair='KK'
#ALICE_dNdy_pi=ALICE_dNdy_K=ALICE_dNdy_p=ALICE_dNdy=0

if chargepair == 'pipi':
  centrality='alice_cent0_5'
  alicedata_filename='alicedata/BF_dy_pipi_C0_5_ALICE.dat'
  YMAX=0.45
if chargepair == 'KK' or chargepair=='piK':
  centrality='alice_cent0_10'
  alicedata_filename='alicedata/BF_dy_KK_C0_10_ALICE.dat'
  YMAX=0.225
#if chargepair == 'pip' or chargepair=='Kp' or chargepair=='pp':

if centrality == 'alice_cent0_5':
  ALICE_dNdy_pi=732.5
  ALICE_dNdy_K=109.0
  ALICE_dNdy_p=33.5
if centrality == 'alice_cent0_10':
  ALICE_dNdy_pi=668.75
  ALICE_dNdy_K=99.75
  ALICE_dNdy_p=30.75

if chargepair == 'pipi':
  ALICE_dNdy=ALICE_dNdy_pi
if chargepair == 'KK' or chargepair=='piK':
  ALICE_dNdy=ALICE_dNdy_K
if chargepair == 'pip' or chargepair=='Kp' or chargepair=='pp':
   ALICE_dNdy=ALICE_dNdy_p

print('ALICE_dNdy=',ALICE_dNdy)

chargesdata = np.loadtxt('../D1/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
x=chargesdata[0]
y1=chargesdata[1]
y_D1=efficiency*y1

chargesdata_D0_5 = np.loadtxt('../D0.5/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
x_D4=chargesdata_D0_5[0]
y1=chargesdata_D0_5[1]
y_D0_5=efficiency*y1

filename='../D2/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat'
chargesdata_D2 = np.loadtxt(filename,skiprows=0,unpack=True)
x_D4=chargesdata_D2[0]
y1=chargesdata_D2[1]
y_D2=efficiency*y1

filename='../D4/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat'
print ('filename=',filename)
chargesdata_D4 = np.loadtxt(filename,skiprows=0,unpack=True)
x_D4=chargesdata_D4[0]
y1=chargesdata_D4[1]
y_D4=efficiency*y1

cascadedata = np.loadtxt('../fromcascade/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)
xc=cascadedata[0]
yc1=cascadedata[1]
yc=efficiency*yc1

ysum_D1=y_D1+yc
ysum_D0_5=y_D0_5+yc
ysum_D2=y_D2+yc
ysum_D4=y_D4+yc

cdata=np.loadtxt('../acc_correction/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_y.dat',skiprows=0,unpack=True)

print('ALICE_dNdy=',ALICE_dNdy)

ysum_D0_5=2.7*ALICE_dNdy*ysum_D0_5/cdata[2]
ysum_D1=2.7*ALICE_dNdy*ysum_D1/cdata[2]
ysum_D2=2.7*ALICE_dNdy*ysum_D2/cdata[2]
ysum_D4=2.7*ALICE_dNdy*ysum_D4/cdata[2]

Dy=0.1
D1norm=D2norm=D0_5norm=D4norm=0.0
for i in range(0,20):
  D0_5norm+=ysum_D0_5[i]*Dy
  D1norm+=ysum_D1[i]*Dy
  D2norm+=ysum_D2[i]*Dy
  D4norm+=ysum_D4[i]*Dy
print('norms are: ',D0_5norm,', ',D1norm,', ',D2norm,', ',D4norm)

alicedata=np.loadtxt(alicedata_filename,skiprows=0,unpack=True)
alice_y=alicedata[0]
alice_bf=alicedata[1]
alice_errors=alicedata[2]

plt.plot(x,0.5*ysum_D0_5,linestyle='-',linewidth=2,marker='o',color='r',label='$0.5D_{\\rm latt}$')
plt.plot(x,0.5*ysum_D1,linestyle='-',linewidth=2,marker='o',color='k',label='$D_{\\rm latt}$')
plt.plot(x,0.5*ysum_D2,linestyle='-',linewidth=2,marker='o',color='g',label='$2D_{\\rm latt}$')
plt.plot(x,0.5*ysum_D4,linestyle='-',linewidth=2,marker='o',color='b',label='$4D_{\\rm latt}$')
plt.errorbar(alice_y,alice_bf,alice_errors,linestyle='None',marker='*',color='k',label='ALICE (prel)')
#plt.plot(aliceTOT_phi,aliceTOT_bf,linestyle='None',marker='*',color='k',label='ALICE (prel)')
ax.legend()

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,2,0.5), minor=False)
ax.set_xticklabels(np.arange(0,2,0.5), minor=False, family='serif')
ax.set_xticks(np.arange(0,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.xlim(0.0,1.8)

ax.set_yticks(np.arange(0,1.0,0.1), minor=False)
ax.set_yticklabels(np.arange(0,1.0,0.1), minor=False, family='serif')
ax.set_yticks(np.arange(0,1.0,0.05), minor=True)
plt.ylim(0.0,YMAX)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\Delta y$', fontsize=18, weight='normal')
plt.ylabel('$B(\Delta y)$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
filename='bf_y_'+chargepair
plt.savefig(filename,format='pdf')
command='open -a Preview '+filename
os.system(command)
#plt.show()
quit()
