#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
//#include "b3d.h"

using namespace std;
int main(){
	const int NID=3;
	int iT=0,id1,id2,pid[NID]={211,321,2212};
	double strangecontent,udcontent,sq,sg,nq;  //ratios to entropy
	FILE *fptr=fopen("udsdens.dat","w");
	CparameterMap parmap;
	parmap.ReadParsFromFile("respars.dat");
	CResList *reslist=new CResList(&parmap);
	double rz3=1.20205693406684822643647,rz4=pow(PI,4/90.0);
	double Tf[40],bfnorm[40][NID][NID]={0.0};
	double taumax=1.0E20;
	
	CResInfoMap::iterator rpos;
	Eigen::Vector3d netnorm;
	int a;
	for(a=0;a<3;a++)
		netnorm[a]=0.0;
	reslist->FindFinalProducts(taumax);
	
	for(iT=2;iT<36;iT++){
		Tf[iT]=5.0*iT;
		nq=(3.0/4.0)*12.0*(4.0*PI/pow(2.0*PI/HBARC,3))*2*pow(Tf[iT],3.0)*rz3;
		sq=4.0*(7.0/8.0)*12.0*(4.0*PI/pow(2.0*PI/HBARC,3))*2*pow(Tf[iT],3.)*rz4;
		sg=4.0*12.0*(4.0*PI/pow(2.0*PI/HBARC,3))*2*pow(Tf[iT],3)*rz4;
		reslist->Tf=Tf[iT];
		reslist->CalcEoSandChiandQdens(reslist->Tf,reslist->Pf,reslist->epsilonf,reslist->nf,reslist->densityf,
		reslist->maxweightf,reslist->chif,strangecontent,udcontent);
		fprintf(fptr,"%5.1f %9.7f %9.7f %9.7f %9.7f\n",Tf[iT],udcontent,strangecontent,nq/(3.0*sq),nq/(3.0*sq+sg));
		reslist->chiinvf=(reslist->chif).inverse();
		//printf("Tf[iT]=%g: norm=%g\n",Tf[iT],norm);
		printf("---------------------------------------------------------------\n");
		for(id1=0;id1<NID;id1++){
			for(id2=id1;id2<NID;id2++){
				bfnorm[iT][id1][id2]=reslist->CalcBalanceNorm(pid[id1],pid[id2],taumax);
				if(id1!=id2)
					bfnorm[iT][id2][id1]=bfnorm[iT][id1][id2];
			}
		}
	}
	fclose(fptr);
	
	fptr=fopen("bfnorm.dat","w");
	fprintf(fptr,"  T     pipi       Kpi       KK         pp        ppi       pK\n");
	for(iT=2;iT<36;iT++){
		fprintf(fptr,"%5.1f ",Tf[iT]);
		for(id1=0;id1<NID;id1++){
			for(id2=0;id2<=id1;id2++){
				fprintf(fptr,"%9.6f ",bfnorm[iT][id1][id2]);
			}
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);	
	return 0;
}