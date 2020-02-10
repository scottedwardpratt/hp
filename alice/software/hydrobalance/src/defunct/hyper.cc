#include "hyper.h"

void CHyperElement::Copy(CHyperElement *oldhyper){
	tau=oldhyper->tau;
	x=oldhyper->x; y=oldhyper->y;
	dOmegaX=oldhyper->dOmegaX;
	dOmegaY=oldhyper->dOmegaY;
	dOmega0=oldhyper->dOmega0;
	ux=oldhyper->ux;
	uy=oldhyper->uy;
	pitildexx=oldhyper->pitildexx;
	pitildexy=oldhyper->pitildexy;
	pitildeyy=oldhyper->pitildeyy;
	T=oldhyper->T;
}

void CHyperElement::CalcDOmegaMax(){
	double u0=sqrt(1.0+ux*ux+uy*uy);
	double dOmega2=dOmega0*dOmega0-dOmegaX*dOmegaX-dOmegaY*dOmegaY;
	udotdOmega=u0*dOmega0+ux*dOmegaX+uy*dOmegaY;
	//if(fabs(ux)>0.1)
	//	printf("ux=%g, dOmegaX=%g, dOmega0=%g\n",ux,dOmegaX,dOmega0);
	if(udotdOmega<0.0){
		printf("udotdOMega<0!!! = %g\n",udotdOmega);
		exit(1);
	}
	dOmegaVec=sqrt(-dOmega2+udotdOmega*udotdOmega);
	dOmegaMax=fabs(udotdOmega)+dOmegaVec;
	vOmega=-udotdOmega/dOmegaVec;
}

void CHyperElement::Print(){
	printf("HyperElement Info:\n");
	printf("tau=%g, T=%g, x=%g, y=%g, ux=%g, uy=%g, dOmega0=%g, dOmegaX=%g, dOmegaY=%g\n",tau,T,x,y,ux,uy,dOmega0,dOmegaX,dOmegaY);
	printf("udotdOmega=%g, dOmegaMax=%g, vOmega=%g, dOmegaVec=%g\n",udotdOmega,dOmegaMax,vOmega,dOmegaVec);
	printf("pitildexx=%g, pitildeyy=%g, pitildexy=%g\n",
	pitildexx,pitildeyy,pitildexy);
	printf("---------------------------------------------------------------\n");
}
