#include "balance.h"
using namespace std;

CBFNumer::CBFNumer(CparameterMap *parmapset){
	npairs=0.0;
	parmap=parmapset;
	Netabins=parmap->getD("BF_NETABINS",50);
	Nybins=parmap->getD("BF_NYBINS",50);
	Nqbins=parmap->getD("BF_NQINVBINS",100);
	Nphibins=parmap->getD("BF_NPHIBINS",36);
	Deta=parmap->getD("BF_DETA",0.1);
	Dy=parmap->getD("BF_DY",0.1);
	Dqinv=parmap->getD("BF_DQINV",10);
	Dphi=180.0/double(Nphibins);
	Bqinv.resize(Nqbins);
	Bqout.resize(Nqbins);
	Bqside.resize(Nqbins);
	Bqlong.resize(Nqbins);
	Beta.resize(Netabins);
	By.resize(Nybins);
	Bphi.resize(2*Nphibins);
	Betas.resize(Netabins);
}

void CBFNumer::Increment(CHBPart *parta,CHBPart *partb,double w){
	int ibin;
	double qinv,qout,qside,qlong,deleta,dely,delphi,deletas;
	double QaQb;
	Misc::outsidelong(parta->p,partb->p,qinv,qout,qside,qlong,deleta,dely,delphi);
	QaQb=(parta->resinfo->charge*parta->bweight)*(partb->resinfo->charge*partb->bweight)*w;
	deletas=fabs(parta->eta-partb->eta);
	npairs+=QaQb;
		
	ibin=floorl(qinv/Dqinv);	
	if(ibin<Nqbins)
		Bqinv[ibin]-=QaQb;
	
	ibin=floorl(qout/Dqinv);
	if(ibin<Nqbins)
		Bqout[ibin]-=QaQb;
	
	ibin=floorl(qside/Dqinv);
	if(ibin<Nqbins)
		Bqside[ibin]-=QaQb;
	
	ibin=floorl(qlong/Dqinv);
	if(ibin<Nqbins)
		Bqlong[ibin]-=QaQb;
	
	ibin=floorl(deleta/Deta);
	if(ibin<Netabins)
		Beta[ibin]-=QaQb;
	
	ibin=floorl(dely/Dy);
	if(ibin<Nybins)
		By[ibin]-=QaQb;
	
	double phia=atan2(parta->p[2],parta->p[1]);
	if(sin(2.0*phia)<0.0)
		delphi=-delphi;

	ibin=Nphibins+floorl(delphi/Dphi);
	if(ibin>=2*Nphibins || ibin<0){
		printf("ibin=%d, >Nphbins, delphi=%g\n",ibin,delphi);
		exit(1);
	}
	if(ibin<2*Nphibins)
		Bphi[ibin]-=QaQb;
	
	ibin=floorl(dely/Dy);
	if(ibin<Nybins){
		Betas[ibin]-=QaQb;
	}
}

CBFDenom::CBFDenom(CparameterMap *parmapset){
	parmap=parmapset;
	Nplus=Nminus=0;
}

void CBFDenom::Increment(CHBPart *part,double w){
	int charge=part->resinfo->charge;
	if(charge==1)
		Nplus+=w;
	else if(charge==-1)
		Nminus+=w;
	else if (fabs(charge)>1){
		printf("charge in CBFDenom::Increment >1!! = %d", charge);
	}
}

void CBalance::IncrementGammaP(CHBPart *parta,CHBPart *partb,double w){
	double phia,phib;
	double QaQb;
	QaQb=(parta->resinfo->charge*parta->bweight)*(partb->resinfo->charge*partb->bweight);
	phia=atan2(parta->p[2],parta->p[1]);
	phib=atan2(partb->p[2],partb->p[1]);
	gammap-=cos(phia+phib)*QaQb*w;
	normtest-=QaQb;
}

void CBalance::PrintBFNumer(){
	
}

void CBalance::PrintBFDenom(){
	printf("pi+,pi-: nplus=%g, minus=%g\n",denom_pi->Nplus,denom_pi->Nminus);
	printf("K+,K-  : nplus=%g, minus=%g\n",denom_K->Nplus,denom_K->Nminus);
	printf("p,pbar : nplus=%g, minus=%g\n",denom_p->Nplus,denom_p->Nminus);
	printf("N+,N-  : nplus=%g, minus=%g\n",denom_charged->Nplus,denom_charged->Nminus);
}

void CBalance::CreateBFArrays(){
	numer_pipi=new CBFNumer(&parmap);
	numer_piK=new CBFNumer(&parmap);
	numer_pip=new CBFNumer(&parmap);
	numer_KK=new CBFNumer(&parmap);
	numer_Kp=new CBFNumer(&parmap);
	numer_pp=new CBFNumer(&parmap);
	numer_charged=new CBFNumer(&parmap);
	numer_charged_phi0=new CBFNumer(&parmap);
	numer_charged_phi45=new CBFNumer(&parmap);
	numer_charged_phi90=new CBFNumer(&parmap);
	
	numer_pipi->name="pipi";
	numer_piK->name="piK";
	numer_pip->name="pip";
	numer_KK->name="KK";
	numer_Kp->name="Kp";
	numer_pp->name="pp";
	numer_charged->name="charged";
	
	bf_pipi=new CBFNumer(&parmap);
	bf_piK=new CBFNumer(&parmap);
	bf_pip=new CBFNumer(&parmap);
	bf_KK=new CBFNumer(&parmap);
	bf_Kp=new CBFNumer(&parmap);
	bf_pp=new CBFNumer(&parmap);
	bf_charged=new CBFNumer(&parmap);
	bf_charged_phi0=new CBFNumer(&parmap);
	bf_charged_phi45=new CBFNumer(&parmap);
	bf_charged_phi90=new CBFNumer(&parmap);
	
	bf_pipi->name="pipi";
	bf_piK->name="piK";
	bf_pip->name="pip";
	bf_KK->name="KK";
	bf_Kp->name="Kp";
	bf_pp->name="pp";
	bf_charged->name="charged";
	bf_charged_phi0->name="charged_phi0";
	bf_charged_phi45->name="charged_phi45";
	bf_charged_phi90->name="charged_phi90";
	
	denom_pi=new CBFDenom(&parmap);
	denom_K=new CBFDenom(&parmap);
	denom_p=new CBFDenom(&parmap);
	denom_charged=new CBFDenom(&parmap);
	denom_charged_phi0=new CBFDenom(&parmap);
	denom_charged_phi45=new CBFDenom(&parmap);
	denom_charged_phi90=new CBFDenom(&parmap);
}

void CBalance::ConstructBFs(){
	double Ntot;
	ConstructBF(numer_pipi,denom_pi,bf_pipi,1.0);
	ConstructBF(numer_piK,denom_pi,bf_piK,0.5);
	ConstructBF(numer_pip,denom_pi,bf_pip,0.5);
	ConstructBF(numer_KK,denom_K,bf_KK,1.0);
	ConstructBF(numer_Kp,denom_K,bf_Kp,0.5);
	ConstructBF(numer_pp,denom_p,bf_pp,1.0);
	//printf("DenomCharged=%g\n",(denom_charged->Nplus+denom_charged->Nminus)
		// /(NSAMPLE_HYDRO2UDS*NSAMPLE_UDS2BAL*TotalVolume));
	ConstructBF(numer_charged,denom_charged,bf_charged,1.0);
	ConstructBF(numer_charged_phi0,denom_charged_phi0,bf_charged_phi0,1.0);
	ConstructBF(numer_charged_phi45,denom_charged_phi45,bf_charged_phi45,1.0);
	ConstructBF(numer_charged_phi90,denom_charged_phi90,bf_charged_phi90,1.0);
	Ntot=denom_charged->Nplus+denom_charged->Nminus;
	//printf("normtest=%g, Ntot=%g, gammap/normtest=%g\n",normtest,Ntot,gammap/normtest);
	gammap=gammap/(Ntot*Ntot/(NSAMPLE_HYDRO2UDS*NSAMPLE_UDS2BAL));
	v2=v2/Ntot;
	printf("gammap=%g, v2=%g\n",gammap,v2);
	FILE *fptr=fopen("gammap.dat","a");
	fprintf(fptr,"%s: gammap=%g\n",qualifier.c_str(),gammap);
	fclose(fptr);
}

void CBalance::ConstructBF(CBFNumer *numer,CBFDenom *denom,CBFNumer *bf,double doublecount){
	bf->Bqinv=numer->Bqinv;
	bf->Bqout=numer->Bqout;
	bf->Bqside=numer->Bqside;
	bf->Bqlong=numer->Bqlong;
	bf->Beta=numer->Beta;
	bf->By=numer->By;
	bf->Bphi=numer->Bphi;
	bf->npairs=numer->npairs;
	int ibin;
	double norm;
	double N=denom->Nplus+denom->Nminus;
	printf("Denom=%g\n",(denom->Nplus+denom->Nminus)/(NSAMPLE_HYDRO2UDS*NSAMPLE_UDS2BAL*TotalVolume));
	for(ibin=0;ibin<numer->Nqbins;ibin++){
		bf->Bqinv[ibin]=doublecount*numer->Bqinv[ibin]/(N*numer->Dqinv);
		bf->Bqout[ibin]=doublecount*numer->Bqout[ibin]/(N*numer->Dqinv);
		bf->Bqside[ibin]=doublecount*numer->Bqside[ibin]/(N*numer->Dqinv);
		bf->Bqlong[ibin]=doublecount*numer->Bqlong[ibin]/(N*numer->Dqinv);
	}
	norm=0.0;
	for(ibin=0;ibin<numer->Netabins;ibin++){
		bf->Beta[ibin]=doublecount*numer->Beta[ibin]/(N*numer->Deta);
		bf->Betas[ibin]=doublecount*numer->Betas[ibin]/(N*numer->Deta);
	}
	for(ibin=0;ibin<numer->Nybins;ibin++){
		bf->By[ibin]=doublecount*numer->By[ibin]/(N*numer->Dy);
	}
	for(ibin=0;ibin<2*numer->Nphibins;ibin++){
		bf->Bphi[ibin]=doublecount*numer->Bphi[ibin]/(N*numer->Dphi);
		norm+=numer->Dphi*bf->Bphi[ibin];
	}
	printf("%7s: normalization=%g, npairs=%g\n",bf->name.c_str(),norm,bf->npairs);
}

void CBalance::WriteBFs(){
	string dirname=bf_results_dirname+"/"+qualifier;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	bf_pipi->Write(dirname);
	bf_piK->Write(dirname);
	bf_pip->Write(dirname);
	bf_KK->Write(dirname);
	bf_Kp->Write(dirname);
	bf_pp->Write(dirname);
	bf_charged->Write(dirname);
	bf_charged_phi0->Write(dirname);
	bf_charged_phi45->Write(dirname);
	bf_charged_phi90->Write(dirname);
}

void CBFNumer::Write(string dirname){
	string filename;
	FILE *fptr;
	int ibin;
	
	string command="mkdir -p "+dirname+"/"+name;
	system(command.c_str());
	
	filename=dirname+"/"+name+"/bf_qinv.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Nqbins;ibin++){
		fprintf(fptr,"%7.2f %10.3e %10.3e %10.3e %10.3e\n",
		(0.5+ibin)*Dqinv,Bqinv[ibin],Bqout[ibin],Bqside[ibin],
		Bqlong[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/bf_y.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Nybins;ibin++){
		fprintf(fptr,"%7.2f %10.3e\n",(0.5+ibin)*Dy,By[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/bf_eta.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Netabins;ibin++){
		fprintf(fptr,"%7.2f %10.3e\n",(0.5+ibin)*Deta,Beta[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/bf_etas.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<Netabins;ibin++){
		fprintf(fptr,"%7.2f %10.3e\n",(0.5+ibin)*Deta,Betas[ibin]);
	}
	fclose(fptr);
	
	filename=dirname+"/"+name+"/bf_phi.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ibin=0;ibin<2*Nphibins;ibin++){
		fprintf(fptr,"%7.2f %10.3e\n",-180.0+(0.5+ibin)*Dphi,Bphi[ibin]);
	}
	fclose(fptr);
		
}

void CBFNumer::Print(){
	int ibin;

	for(ibin=0;ibin<Nqbins;ibin++){
		printf("%6.1f %10.3e %10.3e %10.3e %10.3e\n",
		(0.5+ibin)*Dqinv,Bqinv[ibin],Bqout[ibin],Bqside[ibin],
		Bqlong[ibin]);
	}

	for(ibin=0;ibin<Nybins;ibin++){
		printf("%6.1f %10.3e\n",(0.5+ibin)*Dy,By[ibin]);
	}

	for(ibin=0;ibin<Netabins;ibin++){
		printf("%6.1f %10.3e\n",(0.5+ibin)*Deta,Beta[ibin]);
	}

	for(ibin=0;ibin<Nphibins;ibin++){
		printf("%6.1f %10.3e\n",-180.0+(0.5+ibin)*Dphi,Bphi[ibin]);
	}
	
}
