#include "parametermap.h"
#include "resonances.h"
//#include "b3d.h"

using namespace std;
int main(){
	int pid1,pid2;
	
	CparameterMap parmap;
	parmap.ReadParsFromFile("respars.dat");
	CResList *reslist=new CResList(&parmap);
	
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	Eigen::Vector3d netnorm;
	int a;
	for(a=0;a<3;a++)
		netnorm[a]=0.0;
	reslist->FindFinalProducts();
	
	printf("Enter pid1 & pid2: ");
	scanf("%d %d",&pid1,&pid2);
	reslist->Tf=155.0;
	reslist->CalcEoSandChi(reslist->Tf,reslist->Pf,reslist->epsilonf,reslist->nf,reslist->densityf,reslist->maxweightf,reslist->chif);
	reslist->chiinvf=(reslist->chif).inverse();
	double norm=reslist->CalcBalanceNorm(pid1,pid2);
	printf("norm=%g\n",norm);
	printf("---------------------------------------------------------------\n");
	/*
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	Eigen::Vector3d netnorm;
	int a;
	for(a=0;a<3;a++)
		netnorm[a]=0.0;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		if(!resinfo->decay){
			pid1=resinfo->code;
			norm=reslist->CalcBalanceNorm(pid1,pid2);
			printf("PID1=%d, norm=%g\n",pid1,norm);
			for(a=0;a<3;a++){
				netnorm[a]+=norm*resinfo->q[a];
			}
		}
	}
	for(a=0;a<3;a++)
		printf("netnorm[%d]=%g\n",a,netnorm[a]);
	*/
	return 0;
}