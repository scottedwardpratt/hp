#include "balance.h"
#include "hbdefs.h"
#include "constants.h"

using namespace std;

void CResInfo::MakePart(CHyperElement *hyper,FourVector &p,double maxweight){
	CBalance *cb=reslist->cb;
	double weight,P=cb->Pf,epsilon=cb->epsilonf,T=cb->Tf,lambda=cb->lambdaf;
	bool VISCOUSCORRECTIONS=true;
	bool success=false;
	int nparts=0;
	CHBPart *part;
	CRandy *randy=(reslist->cb)->randy;
	double h=P+epsilon;
	double pdotdOmega,udotOmega,Omega2,eta,y,OmegaMax,pmag,udotp;
	double dOmega0,dOmegaX,dOmegaY;
	double pitilde[4][4];
	int alpha,beta,intweight,n,nsample;
	FourVector pnoviscous,u;
	double m,delN,r[3],w[3];
	CResInfoMap *resmap=&(reslist->resmap);
	CResInfoMap::iterator rpos;
	pitilde[1][1]=hyper->pitildexx;
	pitilde[1][2]=pitilde[2][1]=hyper->pitildexy;
	pitilde[1][3]=pitilde[3][1]=0.0;
	pitilde[2][2]=hyper->pitildeyy;
	pitilde[2][3]=pitilde[3][2]=0.0;
	pitilde[3][3]=-hyper->pitildexx-hyper->pitildeyy;
	dOmega0=hyper->dOmega0;
	dOmegaX=hyper->dOmegaX;
	dOmegaY=hyper->dOmegaY;
	double vOmega=hyper->vOmega;
	//printf("YYYYYYYYYYYYYYYYYYYYYY, lambda=%g, h=%g, P=%g, epsilon=%g\n",
	//lambda,h,P,epsilon);
	//hyper->Print();
	do{
		if(maxweight<0.0){
			m=mass;
			randy->generate_boltzmann(m,T,pnoviscous);
		}
		else{
			m=GenerateThermalMass(maxweight,T);
			randy->generate_boltzmann(m,T,pnoviscous);
		}
		if(VISCOUSCORRECTIONS){
			/*
			pitilde[1][1]=pitilde[2][2]=0.25*h*lambda;
			pitilde[1][2]=pitilde[2][1]=0.0;
			pitilde[3][3]=-0.5*h*lambda;
			*/
			
			for(alpha=1;alpha<4;alpha++){
				p[alpha]=pnoviscous[alpha];
				for(beta=1;beta<4;beta++){
					p[alpha]+=pitilde[alpha][beta]*pnoviscous[beta]/(h*lambda);
				}
			}
			p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+m*m);
		}
		else{
			for(alpha=0;alpha<4;alpha++)
				p[alpha]=pnoviscous[alpha];
		}

		u[1]=hyper->ux; u[2]=hyper->uy;	u[3]=0.0;
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
		udotp=p[0];
		Misc::Boost(u,p);

		pdotdOmega=p[0]*dOmega0+p[1]*dOmegaX+p[2]*dOmegaY;
		weight=pdotdOmega/(udotp*hyper->dOmegaMax);
		if(fabs(weight)>1.0){
			printf("weight out of line, =%g\n",weight);
			Print();
			printf("pdotdOmega=%g, p=(%g,%g,%g,%g)\n",pdotdOmega,p[0],p[1],p[2],p[3]);
			printf("udotOmega=%g, u=(%g,%g,%g,%g), dOmega=(%g,%g,%g)\n",
			hyper->udotdOmega,u[0],u[1],u[2],u[3],dOmega0,dOmegaX,dOmegaY);
			printf("dOmegaMax=%g, udotdOmega=%g\n",hyper->dOmegaMax,hyper->udotdOmega);
			exit(1);
		}
	}while(randy->ran()>weight);
}
