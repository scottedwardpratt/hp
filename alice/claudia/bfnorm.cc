#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
//#include "b3d.h"

using namespace std;
int main(){
	int pid1,pid2;
	double strangecontent,udcontent,sq,sg,nq;  //ratios to entropy
	FILE *fptr=fopen("udsdens.dat","w");
	CparameterMap parmap;
	parmap.ReadParsFromFile("respars.dat");
	CResList *reslist=new CResList(&parmap);
	double rz3=1.20205693406684822643647,rz4=pow(PI,4/90.0);
	
	CResInfoMap::iterator rpos;
	Eigen::Vector3d netnorm;
	int a;
	for(a=0;a<3;a++)
		netnorm[a]=0.0;
	reslist->FindFinalProducts();
	
	pid1=pid2=321;
	//printf("Enter pid1 & pid2: ");
	//scanf("%d %d",&pid1,&pid2);
	for(double Tf=10.0;Tf<176;Tf+=5){
		nq=(3.0/4.0)*12.0*(4.0*PI/pow(2.0*PI/HBARC,3))*2*pow(Tf,3.0)*rz3;
		sq=4.0*(7.0/8.0)*12.0*(4.0*PI/pow(2.0*PI/HBARC,3))*2*pow(Tf,3.)*rz4;
		//ng=16.0*(4.0*PI/pow(2.0*PI/HBARC,3))*2*pow(Tf,3)*rz3;
		sg=4.0*12.0*(4.0*PI/pow(2.0*PI/HBARC,3))*2*pow(Tf,3)*rz4;
		reslist->Tf=Tf;
		reslist->CalcEoSandChiandQdens(reslist->Tf,reslist->Pf,reslist->epsilonf,reslist->nf,reslist->densityf,reslist->maxweightf,reslist->chif,strangecontent,udcontent);
		printf("Tf=%g, strangecontent=%g, udcontent=%g\n",Tf,strangecontent,udcontent);
		fprintf(fptr,"%5.1f %9.7f %9.7f %9.7f %9.7f\n",Tf,udcontent,strangecontent,nq/(3.0*sq),nq/(3.0*sq+sg));
		reslist->chiinvf=(reslist->chif).inverse();
		double norm=reslist->CalcBalanceNorm(pid1,pid2);
		printf("Tf=%g: norm=%g\n",Tf,norm);
		printf("---------------------------------------------------------------\n");
	}
	fclose(fptr);
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