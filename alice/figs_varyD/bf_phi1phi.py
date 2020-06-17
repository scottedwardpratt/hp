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
plt.figure(figsize=(5,10))
fig = plt.figure(1)

############################
centa=40
centb=50
efficiency=0.78
#############################
##### First do from cascade
ax = fig.add_axes([0.17,0.06,0.77,0.31])

if centa == 0 and centb == 5:
  centrality='rhic_cent0_5'
if centa == 10 and centb == 20:
  centrality='rhic_cent10_20'
if centa == 20 and centb == 30:
  centrality='rhic_cent20_30'
if centa == 30 and centb == 40:
  centrality='rhic_cent30_40'
if centa == 40 and centb == 50:
  centrality='rhic_cent40_50'
if centa == 50 and centb == 60:
  centrality='rhic_cent50_60'
  
datafilename='fromcascade/'+centrality+'/results_star/allcharges/bf_phi.dat'
print 'datafilename=',datafilename
mydata = np.loadtxt(datafilename,skiprows=0,unpack=True)
x=mydata[0]
y=mydata[1]*(180/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+y[i]
    phi=-178+4*i
    phibar=phibar+y[i]*phi
    cosdphibar+=y[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CASCADE: All angles: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm
  
datafilename='fromcascade/'+centrality+'/results_star/allcharges_phi0/bf_phi.dat'
mydata0 = np.loadtxt(datafilename,skiprows=0,unpack=True)
x0=mydata0[0]
y0=mydata0[1]*(180.0/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+y0[i]
    phi=-178+4*i
    phibar=phibar+y0[i]*phi
    cosdphibar+=y0[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CASCADE: phia~0: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm

datafilename='fromcascade/'+centrality+'/results_star/allcharges_phi45/bf_phi.dat'
mydata45 = np.loadtxt(datafilename,skiprows=0,unpack=True)
x45=mydata45[0]
y45=mydata45[1]*(180.0/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+y45[i]
    phi=-178+4*i
    phibar=phibar+y45[i]*phi
    cosdphibar+=y45[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CASCADE: phia~45: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm
  
datafilename='fromcascade/'+centrality+'/results_star/allcharges_phi90/bf_phi.dat'
mydata90 = np.loadtxt(datafilename,skiprows=0,unpack=True)
x90=mydata90[0]
y90=mydata90[1]*(180.0/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+y90[i]
    phi=-178+4*i
    phibar=phibar+y90[i]*phi
    cosdphibar+=y90[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CASCADE: phia~90: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm

plt.xlabel('$\Delta\phi$ (degrees)', fontsize=18, weight='normal')
ax.set_xticks(np.arange(-225,225,45), minor=False)
ax.set_xticklabels(np.arange(-225,225,45), minor=False, family='serif')
ax.set_xticks(np.arange(-225,225,15), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.xlim(-180,180)

ax.set_yticks(np.arange(0,2.0,0.04), minor=False)
ax.set_yticklabels(np.arange(0,2.0,0.04), minor=False, family='serif')
ax.set_yticks(np.arange(0,2.0,0.01), minor=True)
plt.ylim(-0.01,0.10)

plt.text(170,0.07,"From\n Cascade",color='k',fontsize='20',ha='right')

#x = np.arange(0,40.1,1)
#y =100.0* x**2
#z=0.5*y
#Use linestyle='None' for no line...
#plt.plot(x,y,marker='None',linestyle='-',linewidth=4,color='k')
plt.plot(x0,y0,marker='o',linestyle='None',color='r')
plt.plot(x45,y45,marker='^',linestyle='None',color='g')
plt.plot(x90,y90,marker='s',linestyle='None',color='b')

#############################
##### No do fromcharges
ax = fig.add_axes([0.17,0.37,0.77,0.31])

datafilename='fromcharges/'+centrality+'/results_star/allcharges/bf_phi.dat'
mydata = np.loadtxt(datafilename,skiprows=0,unpack=True)
x=mydata[0]
yy=efficiency*mydata[1]*(180/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+yy[i]
    phi=-178+4*i
    phibar=phibar+yy[i]*phi
    cosdphibar+=yy[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CHARGES: All angles: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm

datafilename='fromcharges/'+centrality+'/results_star/allcharges_phi0/bf_phi.dat'
mydata0 = np.loadtxt(datafilename,skiprows=0,unpack=True)
x0=mydata0[0]
yy0=efficiency*mydata0[1]*(180.0/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+yy0[i]
    phi=-178+4*i
    phibar=phibar+yy0[i]*phi
    cosdphibar+=yy0[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CHARGES: phia~0: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm

datafilename='fromcharges/'+centrality+'/results_star/allcharges_phi45/bf_phi.dat'
mydata45 = np.loadtxt(datafilename,skiprows=0,unpack=True)
x45=mydata45[0]
yy45=efficiency*mydata45[1]*(180.0/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+yy45[i]
    phi=-178+4*i
    phibar=phibar+yy45[i]*phi
    cosdphibar+=yy45[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CHARGES: phia~45: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm

  
datafilename='fromcharges/'+centrality+'/results_star/allcharges_phi90/bf_phi.dat'
mydata90 = np.loadtxt(datafilename,skiprows=0,unpack=True)
x90=mydata90[0]
yy90=efficiency*mydata90[1]*(180.0/pi)
phibar=0.0
norm=0.0
cosdphibar=0.0
for i in range(0,90):
    norm=norm+yy90[i]
    phi=-178+4*i
    phibar=phibar+yy90[i]*phi
    cosdphibar+=yy90[i]*cos(phi*pi/180.0)
phibar=phibar/norm
cosdphibar=cosdphibar/norm
print 'FROM CHARGES: phia~90: <phi>=',phibar,' <cosDphi>=',cosdphibar/norm

plt.ylabel('$B(\Delta\phi|\phi_1)$',fontsize=22,labelpad=1)
#plt.xlabel()
ax.set_xticks(np.arange(-225,225,45), minor=False)
ax.set_xticks(np.arange(-225,225,15), minor=True)
ax.set_xticklabels([])
plt.xlim(-180,180)

ax.set_yticks(np.arange(0,2.0,0.04), minor=False)
ax.set_yticklabels(np.arange(0,2.0,0.04), minor=False, family='serif')
ax.set_yticks(np.arange(0,2.0,0.01), minor=True)
plt.ylim(-0.01,0.10)

plt.text(-86,0.05,"$\phi_1\\approx 0$",color='r',fontsize='22',ha='right')
plt.text(-85,0.06,"$\phi_1\\approx 45$",color='g',fontsize='22',ha='right')
plt.text(-85,0.07,"$\phi_1\\approx 90$",color='b',fontsize='22',ha='right')
plt.text(170,0.07,"From\n Hydro",color='k',fontsize='20',ha='right')

#dist.s from charges were written in terms of phi2-phi1 instead of phi1-phi2

#x = np.arange(0,40.1,1)
#y =100.0* x**2
#z=0.5*y
#Use linestyle='None' for no line...
#plt.plot(x,yy,marker='None',linestyle='-',linewidth=4,color='k')
plt.plot(x0,yy0,marker='o',linestyle='None',color='r')
plt.plot(x45,yy45,marker='^',linestyle='None',color='g')
plt.plot(x90,yy90,marker='s',linestyle='None',color='b')

################################
####### Now do sum
ax = fig.add_axes([0.17,0.68,0.77,0.31])

yyy=y+yy
yyy0=y0+yy0
yyy45=y45+yy45
yyy90=y90+yy90
#plt.plot(x,yyy,marker='None',linestyle='-',linewidth=4,color='k')
plt.plot(x0,yyy0,marker='o',linestyle='None',color='r')
plt.plot(x45,yyy45,marker='^',linestyle='None',color='g')
plt.plot(x90,yyy90,marker='s',linestyle='None',color='b')

norm=0.0
for i in arange(0,90):
    norm+=4.0*(pi/180.0)*yyy[i]

print('norm=',norm)

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)


ax.set_yticks(np.arange(0,2.0,0.04), minor=False)
ax.set_yticklabels(np.arange(0,2.0,0.04), minor=False, family='serif')
ax.set_yticks(np.arange(0,2.0,0.01), minor=True)
plt.ylim(-0.01,0.18)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

ax.set_xticks(np.arange(-225,225,45), minor=False)
ax.set_xticks(np.arange(-225,225,15), minor=True)
ax.set_xticklabels([])
plt.xlim(-180,180)

plt.text(170,0.145,"Sum",color='k',fontsize='20',ha='right')

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)


plt.savefig('bf_phi1phi.pdf',format='pdf')
os.system('open -a Preview bf_phi1phi.pdf')
#plt.show()
quit()
