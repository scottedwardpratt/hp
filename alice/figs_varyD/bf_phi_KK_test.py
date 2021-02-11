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

chargesdata_Dtest = np.loadtxt('../D0.5/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_phi.dat',skiprows=0,unpack=True)
x_Dtest=chargesdata_Dtest[0]
y1test=chargesdata_Dtest[1]
y2test=y1test
i=0
while(i<28):
  y2test[i]=y1test[27-i]
  i=i+1
y_Dtest=(y1test+y2test)*0.5
y_Dtest=efficiency*y_Dtest*180.0/pi

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

cascadetestdata = np.loadtxt('../fromcascade/model_output/default_sum/'+centrality+'/results_alice/'+chargepair+'/bf_phi.dat',skiprows=0,unpack=True)
xctest=cascadetestdata[0]
yc1test=cascadetestdata[1]
yc2test=yc1test
i=0
while(i<28):
  yc2test[i]=yc1test[27-i]
  i=i+1
yctest=(yc1test+yc2test)*0.5
yctest=efficiency*yctest*180.0/np.pi

ysum_D1=y_D1+yc
ysum_Dtest=y_Dtest+yctest


#ysum_D0_5=ysum_D1=ysum_D2=ysum_D4=yc

Dphi=pi/14.0
D1norm=D1testnorm=0.0
for i in range(0,18):
  D1norm+=ysum_D1[i]*Dphi
  D1testnorm+=ysum_D1[i]*Dphi
print('norms are: ',D1norm,', ',D1testnorm)

alicedata=np.loadtxt('alicedata/alice_KKdata.txt',skiprows=0,unpack=True)
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
alice_norm=0.0
for i in range(1,8):
  aliceTOT_bf[i]=0.5*(alice_bf[8-i]+alice_bf[8+i])
  aliceTOT_errors[i]=0.5*(alice_errors[8-i]+alice_errors[8+i])
  print('i=',i,' angles are ',alice_phi[8+i],' ',alice_phi[8-i])
  aliceTOT_phi[i]=alice_phi[8+i]
  alice_norm=alice_norm+aliceTOT_bf[i]
for i in range(9,14):
  aliceTOT_bf[i]=0.5*(alice_bf[8+i]+alice_bf[28-i])
  aliceTOT_errors[i]=0.5*(alice_errors[8+i]+alice_errors[28-i])*sqrt(2.0)
  aliceTOT_phi[i]=alice_phi[8+i]
  print('i=',i,' angles are ',alice_phi[8+i],' ',alice_phi[28-i]-360.0)
  alice_norm=alice_norm+aliceTOT_bf[i]

alice_norm=alice_norm+aliceTOT_bf[8]+aliceTOT_bf[0]
alice_norm=alice_norm*pi/14
print('alice_norm=',alice_norm)






#stardata=np.loadtxt('stardata/AuAuPhicent0_10.dat',skiprows=0,unpack=True)
#stardata=np.loadtxt('stardata/AuAuPhiCent40_50.dat',skiprows=0,unpack=True)
#xstar=stardata[0]*180.0/pi
#ystar=stardata[1]

#normstar=0.0
#normmodel=0.0
#widthstar=0.0
#widthmodel=0.0
#for i in range(0,20):
#  phimodel=(i+0.5)*(4.0*pi/180.0)
#  normmodel+=(4.0*pi/180.0)*ysum[i]
#  widthmodel+=(4.0*pi/180.0)*ysum[i]*phimodel

#for i in range(0,24):
#  phistar=(i+0.5)*(7.5*pi/180.0)
#  normstar+=(4.0*pi/180.0)*ystar[i]
#  widthstar+=(4.0*pi/180.0)*ystar[i]*phistar
  
#widthstar=widthstar/normstar
#widthmodel=widthmodel/normmodel
#print('normstar=',normstar,' normmodel=',normmodel)
#print('widthstar=',widthstar,' widthmodel=',widthmodel)

#plt.plot(x,ysum_D0_5,linestyle='-',linewidth=2,marker='o',color='r',label='$0.5D_{\\rm latt}$')
plt.plot(x,ysum_D1,linestyle='-',linewidth=2,marker='o',color='r',label='$D_{\\rm latt}$')
plt.plot(x,ysum_Dtest,linestyle='-',linewidth=2,marker='o',color='g',label='(perfect acceptance)')
#plt.plot(x,ysum_D2,linestyle='-',linewidth=2,marker='o',color='g',label='$2D_{\\rm latt}$')
#plt.plot(x,ysum_D4,linestyle='-',linewidth=2,marker='o',color='b',label='$4D_{\\rm latt}$')
plt.errorbar(aliceTOT_phi,aliceTOT_bf*1.5,aliceTOT_errors*1.5,markersize='10',linestyle='None',marker='*',color='k',label='1.5$\cdot$ALICE (prel)')
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
plt.ylim(0.0,0.25)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\Delta\phi$ (degrees)', fontsize=18, weight='normal')
plt.ylabel('$B(\Delta\phi)$',fontsize=18,weight='normal')
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('bf_phi_KK_test.pdf',format='pdf')
#os.system('open -a Preview bf_phi_KK.pdf')
os.system('okular bf_phi_KK_test.pdf &')
#plt.show()
quit()