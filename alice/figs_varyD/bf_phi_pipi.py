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

centrality='alice_cent0_5'
chargepair='pipi'

chargesdata = np.loadtxt('../D1/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_phi.dat',skiprows=0,unpack=True)
x=chargesdata[0]
y1=chargesdata[1]
y2=y1
i=0
while(i<28):
  y2[i]=y1[27-i]
  i=i+1
y_D1=(y1+y2)*0.5
y_D1=efficiency*y_D1*180.0/pi

chargesdata_D0_5 = np.loadtxt('../D0.5/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_phi.dat',skiprows=0,unpack=True)
x_D4=chargesdata_D0_5[0]
y1=chargesdata_D0_5[1]
y2=y1
i=0
while(i<28):
  y2[i]=y1[27-i]
  i=i+1
y_D0_5=(y1+y2)*0.5
y_D0_5=efficiency*y_D0_5*180.0/pi

filename='../D2/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_phi.dat'
chargesdata_D2 = np.loadtxt(filename,skiprows=0,unpack=True)
x_D4=chargesdata_D2[0]
y1=chargesdata_D2[1]
y2=y1
i=0
while(i<28):
  y2[i]=y1[27-i]
  i=i+1
y_D2=(y1+y2)*0.5
y_D2=efficiency*y_D2*180.0/pi

filename='../D4/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_phi.dat'
print ('filename=',filename)
chargesdata_D4 = np.loadtxt(filename,skiprows=0,unpack=True)
x_D4=chargesdata_D4[0]
y1=chargesdata_D4[1]
y2=y1
i=0
while(i<28):
  y2[i]=y1[27-i]
  i=i+1
y_D4=(y1+y2)*0.5
y_D4=efficiency*y_D4*180.0/pi

cascadedata = np.loadtxt('../fromcascade/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_phi.dat',skiprows=0,unpack=True)
xc=cascadedata[0]
yc1=cascadedata[1]
yc2=yc1
i=0
while(i<28):
  yc2[i]=yc1[27-i]
  i=i+1
yc=(yc1+yc2)*0.5
yc=efficiency*yc*180.0/np.pi

ysum_D1=y_D1+yc
ysum_D0_5=y_D0_5+yc
ysum_D2=y_D2+yc
ysum_D4=y_D4+yc

Dphi=pi/14.0
D1norm=D2norm=D0_5norm=D4norm=0.0
for i in range(0,18):
  D0_5norm+=ysum_D0_5[i]*Dphi
  D1norm+=ysum_D1[i]*Dphi
  D2norm+=ysum_D2[i]*Dphi
  D4norm+=ysum_D4[i]*Dphi
print('norms are: ',D0_5norm,', ',D1norm,', ',D2norm,', ',D4norm)

alicedata=np.loadtxt('alicedata/alice_pipidata.txt',skiprows=0,unpack=True)
alice_phi=alicedata[0]
alice_phi=alice_phi*180.0/np.pi
alice_bf=alicedata[1]
alice_errors=alicedata[2]
aliceTOT_bf=zeros(15)
aliceTOT_phi=zeros(15)
aliceTOT_errors=zeros(15)

aliceTOT_bf[0]=alice_bf[8]
aliceTOT_errors[0]=alice_errors[8]
aliceTOT_phi[0]=alice_phi[8]
aliceTOT_bf[14]=alice_bf[22]
aliceTOT_errors[14]=alice_errors[22]
aliceTOT_phi[14]=alice_phi[22]
for i in range(1,8):
  aliceTOT_bf[i]=0.5*(alice_bf[8-i]+alice_bf[8+i])
  aliceTOT_errors[i]=0.5*(alice_errors[8-i]+alice_errors[8+i])
  print('i=',i,' angles are ',alice_phi[8+i],' ',alice_phi[8-i])
  aliceTOT_phi[i]=alice_phi[8+i]
for i in range(9,14):
  aliceTOT_bf[i]=0.5*(alice_bf[8+i]+alice_bf[28-i])
  aliceTOT_errors[i]=0.5*(alice_errors[8+i]+alice_errors[28-i])*sqrt(2.0)
  aliceTOT_phi[i]=alice_phi[8+i]
  print('i=',i,' angles are ',alice_phi[8+i],' ',alice_phi[28-i]-360.0)

plt.plot(x,ysum_D0_5,linestyle='-',linewidth=2,marker='o',color='r',label='$0.5D_{\\rm latt}$')
plt.plot(x,ysum_D1,linestyle='-',linewidth=2,marker='o',color='k',label='$D_{\\rm latt}$')
plt.plot(x,ysum_D2,linestyle='-',linewidth=2,marker='o',color='g',label='$2D_{\\rm latt}$')
plt.plot(x,ysum_D4,linestyle='-',linewidth=2,marker='o',color='b',label='$4D_{\\rm latt}$')
plt.errorbar(aliceTOT_phi,aliceTOT_bf,aliceTOT_errors,markersize='10',linestyle='None',marker='*',color='k',label='ALICE (prel)')
#plt.plot(aliceTOT_phi,aliceTOT_bf,linestyle='None',marker='*',color='k',label='ALICE (prel)')
ax.legend()

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,225,45), minor=False)
ax.set_xticklabels(np.arange(0,225,45), minor=False, family='serif')
ax.set_xticks(np.arange(0,225,15), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.xlim(0.0,181)

ax.set_yticks(np.arange(0,1.0,0.05), minor=False)
ax.set_yticklabels(np.arange(0,1.0,0.05), minor=False, family='serif')
ax.set_yticks(np.arange(0,1.0,0.025), minor=True)
plt.ylim(0.0,0.4)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\Delta\phi$ (degrees)', fontsize=18, weight='normal')
plt.ylabel('$B(\Delta\phi)$',fontsize=18,weight='normal')
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('bf_phi_pipi.pdf',format='pdf')
os.system('open -a Preview bf_phi_pipi.pdf')
#plt.show()
quit()
