#include "balancearrays.h"
#include "part.h"
#include "resonances.h"
#include "parametermap.h"
#include "misc.h"

using namespace std;

CBFNumer::CBFNumer(CparameterMap *parmapset){
	int ieta;
	npairs=0.0;
	parmap=parmapset;
	Netabins=parmap->getD("BF_NETABINS",50);
	Nybins=parmap->getD("BF_NYBINS",50);
	Nqbins=parmap->getD("BF_NQINVBINS",100);
	Nphibins=parmap->getD("BF_NPHIBINS",180);
	Deta=parmap->getD("BF_DETA",0.1);
	Dy=parmap->getD("BF_DELY",0.1);
	Dqinv=parmap->getD("BF_DQINV",10);
	Dphi=360.0/double(Nphibins);
	
	Bqinv.resize(Nqbins,0);
	Bqout.resize(Nqbins,0);
	Bqside.resize(Nqbins,0);
	Bqlong.resize(Nqbins,0);
	Beta.resize(Netabins,0);
	By.resize(Nybins,0);
	Betas.resize(Netabins,0);
	Bphi.resize(Nphibins,0);
	
	Cqinv.resize(Nqbins,0);
	Cqout.resize(Nqbins,0);
	Cqside.resize(Nqbins,0);
	Cqlong.resize(Nqbins,0);
	Ceta.resize(Netabins,0);
	Cy.resize(Nybins,0);
	Cetas.resize(Netabins,0);
	Cphi.resize(Nphibins,0);
	
	Nyphi.resize(Netabins);
	for(ieta=0;ieta<Netabins;ieta++){
		Nyphi[ieta].resize(Nphibins,0);
	}
}

void CBFNumer::Reset(){
	Bqinv.assign(Nqbins,0);
	Bqout.assign(Nqbins,0);
	Bqside.assign(Nqbins,0);
	Bqlong.assign(Nqbins,0);
	Beta.assign(Netabins,0);
	By.assign(Nybins,0);
	Bphi.assign(Nphibins,0);
	Betas.assign(Netabins,0);
	Cqinv.assign(Nqbins,0);
	Cqout.assign(Nqbins,0);
	Cqside.assign(Nqbins,0);
	Cqlong.assign(Nqbins,0);
	Ceta.assign(Netabins,0);
	Cy.assign(Nybins,0);
	Cphi.assign(Nphibins,0);
	Cetas.assign(Netabins,0);
	for(int ieta=0;ieta<Netabins;ieta++){
		Nyphi[ieta].assign(Nphibins,0);
	}
}

void CBFNumer::Increment(CPart *parta,CPart *partb,double effa,double effb){
	int ibin,iphi,iy;
	double qinv,qout,qside,qlong,deleta,dely,delphi,deletas;
	double QaQb,CaCb;
	Misc::outsidelong(parta->p,partb->p,qinv,qout,qside,qlong,deleta,dely,delphi);
	qout=fabs(qout); qside=fabs(qside); qlong=fabs(qlong);
	QaQb=(parta->resinfo->charge*parta->bweight)*(partb->resinfo->charge*partb->bweight);
	QaQb*=effa*effb;
	CaCb=parta->bweight*partb->bweight*effa*effb;
	//printf("QaQb=%g, CaCb=%g\n",QaQb,CaCb);
	deletas=fabs(parta->eta-partb->eta);
	npairs+=QaQb;
	
	if(dely<0.0 || deleta<0.0 || qout<0.0 || qside<0.0 || qlong <0.0 || qinv<0.0 || deletas<0.0){
		printf("bad sign: dely=%g, deleta=%g, deletas=%g, q=(%g,%g,%g,%g)\n",dely,deleta,deletas,qout,qside,qlong,qinv);
		exit(1);
	}
	
	
	ibin=floorl(qinv/Dqinv);	
	if(ibin>=0 && ibin<Nqbins){
		Bqinv[ibin]-=QaQb;
		Cqinv[ibin]+=CaCb;
	}
	
	ibin=floorl(qout/Dqinv);
	if(ibin>=0 && ibin<Nqbins){
		Bqout[ibin]-=QaQb;
		Cqinv[ibin]+=CaCb;
	}
	
	ibin=floorl(qside/Dqinv);
	if(ibin>=0 && ibin<Nqbins){
		Bqside[ibin]-=QaQb;
		Cqside[ibin]+=CaCb;
	}
	
	ibin=floorl(qlong/Dqinv);
	if(ibin>=0 && ibin<Nqbins){
		Bqlong[ibin]-=QaQb;
		Cqlong[ibin]+=CaCb;
	}
	
	ibin=floorl(deleta/Deta);
	if(ibin>=0 && ibin<Netabins){
		Beta[ibin]-=QaQb;
		Ceta[ibin]+=CaCb;
	}
	
	ibin=floorl(dely/Dy);
	if(ibin>=0 && ibin<Nybins){
		By[ibin]-=QaQb;
		Cy[ibin]+=CaCb;
	}
	iy=ibin;
	
	double phia=atan2(parta->p[2],parta->p[1]);
	if(sin(2.0*phia)<0.0){
		delphi=-delphi;
	}
	ibin=floorl((180.0+delphi)/Dphi);
	if(ibin>=0 && ibin<Nphibins){
		Bphi[ibin]-=QaQb;
		Cphi[ibin]+=CaCb;
	}
	iphi=ibin;
	
	ibin=floorl(deletas/Dy);
	if(ibin>=0 && ibin<Nybins){
		Betas[ibin]-=QaQb;
		Cetas[ibin]+=CaCb;
	}
	
	if(iphi<Nphibins && iy<Netabins)
		Nyphi[iy][iphi]-=QaQb;
	
}

CBFDenom::CBFDenom(CparameterMap *parmapset){
	parmap=parmapset;
	Nplus=Nminus=dNdy=0.0;
}

void CBFDenom::Reset(){
	Nplus=Nminus=dNdy=0.0;
}

void CBFDenom::Increment(CPart *part,double eff){
	int charge=part->resinfo->charge;
	if(charge==1)
		Nplus+=eff;
	else if(charge==-1)
		Nminus+=eff;
	else if (fabs(charge)>1){
		printf("charge in CBFDenom::Increment >1!! = %d", charge);
	}
}

void CBFNumer::WriteNumer(string dirname,string numertype,bool NoQ){
	string filename;
	FILE *fptr;
	int ibin,jbin;

	string command="mkdir -p "+dirname+"/"+name;
	system(command.c_str());
	
	if(!NoQ){
		filename=dirname+"/"+name+"/"+numertype+"_qinv.dat";
		fptr=fopen(filename.c_str(),"w");
		for(ibin=0;ibin<Nqbins;ibin++){
			fprintf(fptr,"%7.2f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
			(0.5+ibin)*Dqinv,Bqinv[ibin],Bqout[ibin],Bqside[ibin],Bqlong[ibin],
			Cqinv[ibin],Cqout[ibin],Cqside[ibin],Cqlong[ibin]);
		}
		fclose(fptr);
	}
	
	filename=dirname+"/"+name+"/"+numertype+"_y.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Nybins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Dy,By[ibin],Cy[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_eta.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Netabins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Deta,Beta[ibin],Ceta[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_etas.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Netabins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",(0.5+ibin)*Deta,Betas[ibin],Betas[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_phi.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Nphibins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e\n",-180.0+(0.5+ibin)*Dphi,Bphi[ibin],Cphi[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/"+numertype+"_yphi.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#Nybins= %d Nphibins= %d\n",Netabins,Nphibins);
	for(ibin=0;ibin<Netabins;ibin++){
		for(jbin=0;jbin<Nphibins;jbin++)
			fprintf(fptr,"%12.4e ",Nyphi[ibin][jbin]);
		fprintf(fptr,"\n");
	}
	fclose(fptr);
	
}

void CBFNumer::Print(){
	int ibin;
	/*
	for(ibin=0;ibin<Nqbins;ibin++){
		printf("%6.1f %10.3e %10.3e %10.3e %10.3e\n",
		(0.5+ibin)*Dqinv,Bqinv[ibin],Bqout[ibin],Bqside[ibin],
		Bqlong[ibin]);
	}
	*/
	printf("----- By ------\n");
	for(ibin=0;ibin<Nybins;ibin++){
		printf("%6.1f %10.3e\n",(0.5+ibin)*Dy,By[ibin]);
	}
	printf("----- Beta -----\n");
	for(ibin=0;ibin<Netabins;ibin++){
		printf("%6.1f %10.3e\n",(0.5+ibin)*Deta,Beta[ibin]);
	}
	printf("----- Bphi -----\n");
	for(ibin=0;ibin<Nphibins;ibin++){
		printf("%3d: %6.1f %10.3e\n",ibin,-180.0+(0.5+ibin)*Dphi,Bphi[ibin]);
	}
	
}
