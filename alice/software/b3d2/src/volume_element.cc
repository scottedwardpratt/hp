#include "sampler.h"
#include "resonances.h"
#include "part.h"
#include "randy.h"
#include "misc.h"

using namespace std;

void CvolumeElement2D::Initialize(){
	int alpha,beta;
	T=0.0;
	(*density).clear();
	density=NULL;
	for(alpha=0;alpha<4;alpha++){
		Omega[alpha]=0.0;
		for(beta=0;beta<4;beta++){
			pitilde[alpha][beta]=0.0;
		}
	}
}

CvolumeElement2D::CvolumeElement2D(){
	int alpha,beta;
	pitilde=new double *[4];
	for(alpha=0;alpha<4;alpha++){
		Omega[alpha]=0.0;
		pitilde[alpha]=new double[4];
		for(beta=0;beta<4;beta++){
			pitilde[alpha][beta]=0.0;
		}
	}
}

void CvolumeElement2D::Print(){
	int alpha,beta;
	printf("tau=%g, xy=(%g,%g), u=(%g,%g)\n",tau,x,y,ux,uy);
	printf("Omega=(%g,%g,%g,%g), OmegaMax=%g, udotdOmega=%g\n",Omega[0],Omega[1],Omega[2],Omega[3],OmegaMax,udotdOmega);
	printf("T=%g, epsilon=%g, P=%g, nhadrons=%g, lambda=%g, Xscale=%g\n",T,epsilon,P,nhadrons,lambda,Xscale);
	printf(" --------  pi_ij  --------\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			printf("%10.4f ",pitilde[alpha][beta]);
		}
		printf("\n");
	}
}

// danger below: if density is a pointer to some a density array used by some other element
// the density will be deleted everywhere

void CvolumeElement2D::CopyEquilibriumQuantities(CvolumeElement2D *element){
	T=element->T;
	P=element->P;
	epsilon=element->epsilon;
	lambda=element->lambda;
	nhadrons=element->nhadrons;
	if(density!=NULL && density!=element->density){
		(*density).clear();
		density=NULL;
	}
	density=element->density;
}

int CvolumeElement2D::MakeParts(){
	int nparts=0;
	CPart *part;
	double h=P+epsilon;
	FourVector g={1,-1,-1,-1};
	double weight,pdotdOmega,mass,et,eta,y,ran1,ran2;
	int ires,alpha,beta,intweight,n,nsample=sampler->b3d->NSAMPLE; 
	bool idecay=0;
	bool create_negparts=sampler->b3d->COOPERFRYE_CREATENEGPARTS;
	FourVector p,plab,pnoviscous,u;
	double delN,r[3],w[3];
	double udotdOmega;
	double delNtot=nhadrons*OmegaMax*Xscale*nsample;
	//printf("------ delNtot=%g,OmegaMax=%g,Xscale=%g,nhadrons=%g,nsample=%d\n",delNtot,OmegaMax,Xscale,nhadrons,nsample);
	CResInfoMap *resmap=&sampler->reslist->resmap;
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	CRandy *randy=sampler->randy;
	//printf("In CvolumeElement2D::MakeParts, OmegaMax=%g,nhadrons=%g,Xscale=%g,delNtot=%g\n",OmegaMax,nhadrons,Xscale,delNtot);
	//printf("Omega=(%g,%g,%g,%g)\n",Omega[0],Omega[1],Omega[2],Omega[3]);
	//Misc::Pause();

	if(sampler->cummulative_N+delNtot > sampler->cummulative_random){
		for(rpos=resmap->begin();rpos!=resmap->end();rpos++){
			resinfo=rpos->second;
			if(resinfo->code==22){
				rpos++;
				resinfo=rpos->second;
			}
			ires=resinfo->ires;
			delN=(*density)[ires]*OmegaMax*Xscale*nsample;
			sampler->cummulative_N+=delN;
			while(sampler->cummulative_N>sampler->cummulative_random){
				ires=resinfo->ires;
				double width=resinfo->width;
				do{
					if(width>1.0)
						mass=resinfo->GenerateThermalMass(sampler->maxweightf[ires], sampler->Tf);
					else
						mass=resinfo->mass;
					n=1;
					randy->generate_boltzmann(resinfo->mass,T/double(n),pnoviscous);
					if(sampler->VISCOUSCORRECTIONS){
						for(alpha=1;alpha<4;alpha++){
							p[alpha]=pnoviscous[alpha];
							for(beta=1;beta<4;beta++){
								p[alpha]-=pitilde[alpha][beta]*g[beta]*pnoviscous[beta]/((epsilon+P)*lambda);
							}
						}
						p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+mass*mass);
					}
					else{
						for(alpha=0;alpha<4;alpha++)
							p[alpha]=pnoviscous[alpha];
					}
					ran1=randy->ran();
					ran2=sqrt(randy->ran());
					w[0]=(1.0-ran1)*ran2;
					w[1]=ran1*ran2;
					w[2]=1.0-ran2;
					if(sampler->TRIANGLE_FORMAT){
						u[1]=w[0]*vertex[0]->ux+w[1]*vertex[1]->ux+w[2]*vertex[2]->ux;
						u[2]=w[0]*vertex[0]->uy+w[1]*vertex[1]->uy+w[2]*vertex[2]->uy;
					}
					else{
						u[1]=ux; u[2]=uy;
					}
					u[3]=0.0;
					u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
					Misc::Boost(u,p,plab);
					pdotdOmega=plab[0]*Omega[0]-plab[1]*Omega[1]-plab[2]*Omega[2];
					weight=pdotdOmega/(p[0]*OmegaMax);
				}while(weight<0.0);
				if(abs(weight)>1.0){
					printf("weight out of line, =%g\n",weight);
					Print();
					printf("p=(%g,%g,%g,%g), mass=%g\n",p[0],p[1],p[2],p[3],mass);
					exit(1);
				}
				if(randy->ran()<fabs(weight)){
					part=sampler->b3d->GetDeadPart();
					if(sampler->TRIANGLE_FORMAT){
						for(alpha=0;alpha<3;alpha++)
							r[alpha]=w[0]*vertex[0]->r[alpha]+w[1]*vertex[1]->r[alpha]+w[2]*vertex[2]->r[alpha];
					
						if(randy->ran()<0.5){
							r[1]=-r[1];
							plab[1]=-plab[1];
						}
						if(randy->ran()<0.5){
							r[2]=-r[2];
							plab[2]=-plab[2];
						}
					}
					else{
						r[0]=tau; r[1]=x; r[2]=y;
					}
					intweight=1;
					y=atanh(plab[3]/plab[0]);
					eta=(1.0-2.0*randy->ran())*sampler->ETAMAX;
					y+=eta;
#ifdef __SAMPLER_WRITE_XY__
					fprintf(sampler->xyfptr,"%g %g %g %g\n",r[0],r[1],r[2],eta);
#endif
					if(create_negparts || intweight>0){
						//part->Init(resinfo->code,r[1],r[2],r[0],eta,plab[1],plab[2],mass,y,intweight);
						part->InitBalance(resinfo->code,r[1],r[2],r[0],eta,plab[1],plab[2],mass,y,intweight,-1);
					}
					nparts+=1;
				}
				sampler->cummulative_random+=-log(randy->ran());
			}
		}
	}
	else{
		sampler->cummulative_N+=delNtot;
	}
	return nparts;
}

void CvolumeElement2D::CalcOmegaMax(){
	double omax,Omega2,udotdOmega,u0;
	Omega2=Omega[0]*Omega[0]-Omega[1]*Omega[1]-Omega[2]*Omega[2];
	if(sampler->TRIANGLE_FORMAT){
		OmegaMax=0.0;
		for(int i=0;i<3;i++){
			u0=sqrt(1.0+vertex[i]->ux*vertex[i]->ux+vertex[i]->uy*vertex[i]->uy);
			udotdOmega=u0*Omega[0]-vertex[i]->ux*Omega[1]-vertex[i]->uy*Omega[2];
			omax=fabs(udotdOmega+sqrt(-Omega2+udotdOmega*udotdOmega));
			if(omax>OmegaMax)
				OmegaMax=omax;
		}
	}
	else{
		u0=sqrt(1.0+ux*ux+uy*uy);
		udotdOmega=u0*Omega[0]-ux*Omega[1]-uy*Omega[2];
		OmegaMax=fabs(udotdOmega+sqrt(-Omega2+udotdOmega*udotdOmega));
	}
}

void CvolumeElement2D::CalcOmegaMax3D(){
	double Omega2,udotdOmega,u0;
	Omega2=Omega[0]*Omega[0]-Omega[1]*Omega[1]-Omega[2]*Omega[2]-Omega[3]*Omega[3];
	u0=sqrt(1.0+ux*ux+uy*uy+uz*uz);
	udotdOmega=u0*Omega[0]-ux*Omega[1]-uy*Omega[2]-uz*Omega[3];
	OmegaMax=fabs(udotdOmega+sqrt(-Omega2+udotdOmega*udotdOmega));
}

/*
// This only used to fill out shear tensor in lab frame
void CvolumeElement::FillOutShearTensor(double &pixx,double &pixy,double &pixz,double &piyy,double &piyz,double &pizz){
//Beginning with pixx,pixy,pixz,piyy,piyz,pizz, fills out remaining tensor and makes it traceless
int alpha,beta;
double trpi;
const int g[4]={1,-1,-1,-1};
pi[1][1]=pixx; pi[1][2]=pixy; pi[1][3]=pixz;
pi[2][2]=piyy; pi[2][3]=piyz;
pi[3][3]=pizz;
pi[2][1]=pi[1][2];
pi[3][1]=pi[1][3];
pi[3][2]=pi[2][3];
// choose other elements so that u_mu pi^{mu nu}=0
for(alpha=1;alpha<4;alpha++){
pi[0][alpha]=0.0;
for(beta=1;beta<4;beta++)
pi[0][alpha]+=pi[beta][alpha]*u[beta]/u[0];
pi[alpha][0]=pi[0][alpha];
}
pi[0][0]=0.0;
for(alpha=1;alpha<4;alpha++)
pi[0][0]+=pi[0][alpha]*u[alpha]/u[0];
// Make tensor traceless -- no bulk correction
trpi=pi[0][0]-pi[1][1]-pi[2][2]-pi[3][3];
for(alpha=0;alpha<4;alpha++){
pi[alpha][alpha]-=(trpi/3.0)*g[alpha];
for(beta=0;beta<4;beta++){
pi[alpha][beta]+=(trpi/3.0)*u[alpha]*u[beta];
}
}
}
*/
