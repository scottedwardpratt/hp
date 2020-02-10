#include "balance.h"
#include "hydro2uds.h"
#include "constants.h"

using namespace std;
using namespace boost::math;

bool CHydroBalance::FakeReadOSCAR(CHydroMesh *hydromesh){
	double deltau=0.02,tau0=0.6,xmin=-13.0,xmax=13.0,ymin=-13.0,ymax=13.0;
	double udotr,x,y,rmax=0.0,highestT,tau,tauf,R,R0=10.0,uR0=1.0;
	double xbar=0.0,ybar=0.0,norm=0.0,uR=1.0;
	double Stot=0.0;
	bool keepgoing=true;
	int olditau,ix,iy,iline,flag,alpha,beta,nplus=0,nminus=0;
	char dummy[300];
	double **pitilde;
	FourVector u;
	double e,p,t,vx,vy,r;
	double etaovers,taupi,tauPI,u0;
	pitilde=new double*[4];
	for(alpha=0;alpha<4;alpha++){
		pitilde[alpha]=new double[4];
		for(beta=0;beta<4;beta++){
			pitilde[alpha][beta]=0.0;
		}
	}
	tau=TAU0+itauread*DELTAU;
	tauf=10.0+TAU0;
	hydromesh->tau=tau;
	R=R0-R0*(tau-TAU0)/(tauf-TAU0);
	highestT=0.0;
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			hydromesh->pitildexx[ix][iy]=pitilde[1][1];
			hydromesh->pitildexy[ix][iy]=pitilde[1][2];
			hydromesh->pitildeyy[ix][iy]=pitilde[2][2];
			hydromesh->GetXY(ix,iy,x,y);
			r=sqrt(x*x+y*y);
			t=Tf-0.10*(r-R)/R0;
			if(t>highestT)
				highestT=t;
			hydromesh->T[ix][iy]=t;
			u[1]=uR0*(x/R0);
			u[2]=uR0*(y/R0);
			u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
			u[3]=0.0;
			hydromesh->UX[ix][iy]=u[1];
			hydromesh->UY[ix][iy]=u[2];
		}
	}
	//printf("tau=%g, highestT=%g\n",tau,highestT);
	for(alpha=0;alpha<4;alpha++){
		delete pitilde[alpha];
	}
	delete pitilde;
	itauread+=1;
	if(highestT<Tf)
		keepgoing=false;
	tau0readcheck=false;
	tau0check=false;
	return keepgoing;
}

double CHydroBalance::SpectraFromHyper(double mass,double px,double py){
	double spectra,E,k0,k1,Omega0sum=0.0,DTAU=0.02;
	double ux,uy,u0,tau,dOmegaX,dOmegaY,dOmega0,udotp,pdotOmega,x,y,r,Z;
	int ix,iy;
	CHyperElement *hyper;
	double factor=pow(2.0*PI*HBARC/1000.0,-3);
	spectra=0.0;
	E=sqrt(mass*mass+px*px+py*py);
	list<CHyperElement *>::iterator it;
	Omega0sum=0.0;
	for(it=hyperlist.begin();it!=hyperlist.end();++it){
		hyper=*it;
		tau=hyper->tau;
		dOmegaX=hyper->dOmegaX; dOmegaY=hyper->dOmegaY; dOmega0=hyper->dOmega0;
		Omega0sum+=dOmega0/tau;
		//Omega0sum+=dOmega0/(tau-0.5*DTAU);
		ux=hyper->ux; uy=hyper->uy;
		x=hyper->x; y=hyper->y;
		mesh->GetIXIY_lower(x,y,ix,iy);
		r=sqrt(hyper->x*hyper->x+hyper->y*hyper->y);
			
		u0=sqrt(1.0+ux*ux+uy*uy);
		udotp=u0*E-ux*px-uy*py;
		Z=u0*E/Tf;
		k0 = cyl_bessel_k(0,Z);
		k1 = cyl_bessel_k(1,Z);
		/*
		double k0test=0.0,k1test=0.0,deta=0.001,eta;
		for(eta=0.5*deta;eta<100;eta+=deta){
			k0test+=deta*exp(-Z*cosh(eta));
			k1test+=deta*cosh(eta)*exp(-Z*cosh(eta));
		}
		printf("Z=%g, k0=%g=?%g,  k1=%g=?%g\n",Z,k0,k0test,k1,k1test);
		*/
		
		pdotOmega=E*dOmega0+px*dOmegaX+py*dOmegaY;
		//pdotOmega=E*dOmega0;
		if(pdotOmega>0.0){
			spectra+=2.0*k0*factor*(px*dOmegaX+py*dOmegaY)
				*exp((ux*px+uy*py)/Tf);
			spectra+=2.0*k1*factor*(E*dOmega0)*exp((ux*px+uy*py)/Tf);
		}
	}
	//printf("Omega0sum=%g, pi*R^2=%g\n",Omega0sum,100.0*PI);
	
	return spectra;
}
