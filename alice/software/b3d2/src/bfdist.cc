#include "balancearrays.h"
#include "b3d.h"
#include "acceptance.h"

using namespace std;

void CBFRapDist::Increment(double y,int charge,double eff){
	int iy;
	iy=floorl((y+BF_YMAX)/dely);
	if(iy<0 || iy>=nmax){
		printf("iy=%d, nmax=%d\n",iy,nmax);
		exit(1);
	}
	if(charge>0){
		Nplus[iy]+=eff;
		eff2Nplus[iy]+=eff*eff;
	}
	else{
		Nminus[iy]+=eff;
		eff2Nminus[iy]+=eff*eff;
	}
}

void CBFPhiDist::Increment(double phi,int charge,double eff){
	int iphi;
	while(phi>PI)
		phi-=2.0*PI;
	while(phi<-PI)
		phi+=2.0*PI;
	phi*=180.0/PI;
	iphi=floorl((phi+180.0)/delphi);
	if(iphi<0 || iphi>=nmax){
		printf("iphi=%d, nmax=%d\n",iphi,nmax);
		exit(1);
	}
	if(charge>0){
		Nplus[iphi]+=eff;
		eff2Nplus[iphi]+=eff*eff;
	}
	else{
		Nminus[iphi]+=eff;
		eff2Nminus[iphi]+=eff*eff;
	}	
}

void CBFRapDist::Print(){
	for(int iy=0;iy<nmax;iy++){
		printf("%2d %6.3f %8.5f %8.5f\n",iy,(0.5+iy)*dely,Nplus[iy],Nminus[iy]);
	}
}

void CBFPhiDist::Print(){
	for(int iphi=0;iphi<nmax;iphi++){
		printf("%2d %6.3f %8.5f %8.5f\n",iphi,(0.5+iphi)*delphi*180.0/PI,Nplus[iphi],Nminus[iphi]);
	}
}
