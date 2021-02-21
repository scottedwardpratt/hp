#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
//#include "b3d.h"

using namespace std;
int main(){
	//const int NID=27;
	//int pid[NID]={-3112,3112,-3122,3122,-3212,3212,-3222,3222,-3312,3312,-3322,3322,-3334,3334,111,-211,211,-311,311,-321,321,-2112,2112,-2212,2212,221,331};
	const int NID=12;
	int pid[NID]={3112,3122,3212,3222,3312,3322,3334,211,311,321,2112,2212};
	//const int NID=4;
	//int pid[NID]={211,321,2112,2212};
	
	int iT=0,id1,id2;
	double strangecontent,udcontent,sq,sg,nq;  //ratios to entropy
	FILE *fptr=fopen("udsdens.dat","w");
	CparameterMap parmap;
	parmap.ReadParsFromFile("respars.dat");
	CResList *reslist=new CResList(&parmap);
	double rz3=1.20205693406684822643647,rz4=pow(PI,4/90.0);
	double Tf[40],bfnorm[40][NID][NID]={0.0};
	double taumax=100.0;
	
	CResInfoMap::iterator rpos;
	CResInfo *resinfo1,*resinfo2;
	double netcharge,netbaryon,netstrange;
	reslist->FindFinalProducts(taumax);
	
	//for(iT=2;iT<36;iT++){
	for(iT=30;iT<31;iT++){
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
		for(id1=0;id1<NID;id1++){
			for(id2=0;id2<NID;id2++){
				bfnorm[iT][id1][id2]=reslist->CalcBalanceNorm(pid[id1],pid[id2],taumax);
			}
		}
	}
	fclose(fptr);
	
	fptr=fopen("bfnorm.dat","w");
	fprintf(fptr,"  T     pipi       Kpi       KK         pp        ppi       pK\n");
	//for(iT=2;iT<36;iT++){
	for(iT=30;iT<31;iT++){
		fprintf(fptr,"%5.1f ",Tf[iT]);
		for(id2=0;id2<NID;id2++){
			resinfo2=reslist->GetResInfoPtr(pid[id2]);
			printf("Given %d:\n",resinfo2->code);
			netcharge=netbaryon=netstrange=0.0;
			for(id1=0;id1<NID;id1++){
				resinfo1=reslist->GetResInfoPtr(pid[id1]);
				fprintf(fptr,"%9.6f ",bfnorm[iT][id1][id2]);
				netcharge+=bfnorm[iT][id1][id2]*resinfo1->charge;
				netbaryon+=bfnorm[iT][id1][id2]*resinfo1->baryon;
				netstrange+=bfnorm[iT][id1][id2]*resinfo1->strange;
			}
			printf("%5d: netQ=%12.8f, netB=%12.8f, netS=%12.8f\n",resinfo2->code,netcharge,netbaryon,netstrange);
			printf("-------------------\n");
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);	
	return 0;
}