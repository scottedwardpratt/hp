#include "hydro2uds.h"
#include "hyper.h"
#include "misc.h"

using namespace std;

void CHydroBalance::HyperFind(){
	int ix,iy,a,b;
	double u0,omega,gradT2,udotgradT;
	double w,dTdx,dTdy,dTdt,dV;
	CHyperElement hyper;
	CHyperElement *newhyper;
	bool GGTt,GGTx,GGTy;
	if(!tau0check){
		for(ix=1;ix<mesh->NX-1;ix++){
			for(iy=1;iy<mesh->NY-1;iy++){
				if(ix<0 || iy<0 || ix==CHydroMesh::NX || iy==CHydroMesh::NY){
					printf("CHydroBalance::HyperFind() -- ix=%d, iy=%d\n",ix,iy);
					exit(1);
				}
				GetGradT(ix,iy,dTdt,dTdx,dTdy,GGTt,GGTx,GGTy);
				if(GGTt || GGTx || GGTy){
					hyper.tau=0.5*(newmesh->tau+mesh->tau);
					GetXYBar(ix,iy,hyper.x,hyper.y);
					GetUxyBar(ix,iy,hyper.ux,hyper.uy);
					GetPiTildeBar(ix,iy,hyper.pitildexx,hyper.pitildexy,
					hyper.pitildeyy);
					GetTBar(ix,iy,hyper.T);
					if(GetDOmega(dTdt,dTdx,dTdy,
					hyper.dOmega0,hyper.dOmegaX,hyper.dOmegaY,GGTt,GGTx,GGTy)){
						for(a=0;a<3;a++){
							for(b=0;b<3;b++){
								chitothyper(a,b)+=hyper.udotdOmega*chif(a,b);
							}
						}
						newhyper=new CHyperElement;
						newhyper->Copy(&hyper);
						hyperlist.push_back(newhyper);
					}
				}
			}
		}
	}
}

bool CHydroBalance::GetGradT(int ix,int iy,
double &dTdt,double &dTdx,double &dTdy,bool &GGTt,bool &GGTx,bool &GGTy){
	bool hypercheck=false;
	double Txplus,Txminus,Typlus,Tyminus,Ttplus,Ttminus;
	double DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX;
	int NX,NY;
	mesh->GetDimensions(NX,NY,DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX);
	
	Txminus=0.25*(mesh->T[ix][iy]+newmesh->T[ix][iy]
		+mesh->T[ix][iy+1]+newmesh->T[ix][iy+1]);
	Txplus=0.25*(mesh->T[ix+1][iy]+newmesh->T[ix+1][iy]
		+mesh->T[ix+1][iy+1]+newmesh->T[ix+1][iy+1]);
	
	Tyminus=0.25*(mesh->T[ix][iy]+newmesh->T[ix][iy]
		+mesh->T[ix+1][iy]+newmesh->T[ix+1][iy]);
	Typlus=0.25*(mesh->T[ix][iy+1]+newmesh->T[ix][iy+1]
		+mesh->T[ix+1][iy+1]+newmesh->T[ix+1][iy+1]);
	
	Ttminus=0.25*(mesh->T[ix][iy]+mesh->T[ix][iy+1]
		+mesh->T[ix+1][iy]+mesh->T[ix+1][iy+1]);
	Ttplus=0.25*(newmesh->T[ix][iy]+newmesh->T[ix][iy+1]
		+newmesh->T[ix+1][iy]+newmesh->T[ix+1][iy+1]);
	
	dTdx=-(Txplus-Txminus)/DX;
	dTdy=-(Typlus-Tyminus)/DY;
	dTdt=(Ttplus-Ttminus)/DELTAU;
	
	GGTt=GGTx=GGTy=false;
	if((Ttplus-Tf)*(Ttminus-Tf)<0.0)
		GGTt=true;
	if((Txplus-Tf)*(Txminus-Tf)<0.0)
		GGTx=true;
	if((Typlus-Tf)*(Tyminus-Tf)<0.0)
		GGTy=true;
	if(GGTt || GGTx || GGTy)
		hypercheck=true;
	return hypercheck;
}

bool CHydroBalance::GetDOmega(double dTdt,double dTdx,double dTdy,
double &dOmega0,double &dOmegaX,double &dOmegaY,bool GGTt,bool GGTx,
bool GGTy){
	double dV,tau=0.5*(newmesh->tau+mesh->tau);
	double dTdr=sqrt(dTdx*dTdx+dTdy*dTdy);
	bool success=false;
	dOmega0=dOmegaX=dOmegaY=0.0;
	if(fabs(dTdt)>(DX/DELTAU)*dTdr){
		if(GGTt){
			dV=tau*DX*DY;
			dOmega0=-dV*dTdt/fabs(dTdt);
			dOmegaX=-dV*dTdx/fabs(dTdt);
			dOmegaY=-dV*dTdy/fabs(dTdt);
			success=true;
		}
	}
	else if(fabs(dTdx)>fabs(dTdy)){
		if(GGTx){
			dV=tau*DELTAU*DY;
			dOmega0=-dV*dTdt/fabs(dTdx);
			dOmegaX=-dV*dTdx/fabs(dTdx);
			dOmegaY=-dV*dTdy/fabs(dTdx);
			success=true;
		}
	}
	else{
		if(GGTy){
			dV=tau*DELTAU*DX;
			dOmega0=-dV*dTdt/fabs(dTdy);
			dOmegaX=-dV*dTdx/fabs(dTdy);
			dOmegaY=-dV*dTdy/fabs(dTdy);
			success=true;
		}
	}
	return success;
}

