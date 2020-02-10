#ifndef __PART_CC__
#define __PART_CC__

#include "part.h"
#include "constants.h"

CBalance *CHBPart::cb=NULL;

CHBPart::CHBPart(){
	tau0=0.0;
	for(int alpha=0;alpha<4;alpha++){
		r[alpha]=p[alpha]=0.0;
	}
	msquared=0.0;
	y=eta=0.0;
	bweight=1.0;
	balanceID=-1;
	resinfo=NULL;
	partmap=NULL;
}

CHBPart::~CHBPart(){
}

void CHBPart::Copy(CHBPart *part){  //copies all info 
	tau0=part->tau0;
	for(int alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
		r[alpha]=part->r[alpha];
	}
	msquared=part->msquared;
	y=part->y;
	eta=part->eta;
	bweight=part->bweight;
	balanceID=part->balanceID;
	resinfo=part->resinfo;
	//partmap=part->partmap;
}

void CHBPart::CopyPositionInfo(CHBPart *part){  //copies all info except actionmap
	tau0=part->tau0;
	eta=part->eta;
	for(int alpha=0;alpha<4;alpha++){
		r[alpha]=part->r[alpha];
	}
}

void CHBPart::CopyMomentumInfo(CHBPart *part){  //copies all info except actionmap
	y=part->y;
	msquared=part->msquared;
	for(int alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
	}
}

void CHBPart::Init(int IDset,double rxset,double ryset,double tauset,
double etaset,double pxset,double pyset,double mset,double rapidityset,
int bweightset,int balanceIDset){
	badmother=0;
	double et;
	CResInfo *resinfoptr;
	int ID;
	resinfo=cb->reslist->GetResInfoPtr(IDset);
	ID=resinfo->code;
	if(ID!=IDset){
		printf("ID mismatch, ID=%d, resinfo->codeID=%d\n",IDset,ID);
	}
	p[1]=pxset; p[2]=pyset; msquared=mset*mset; y=rapidityset;
	r[1]=rxset; r[2]=ryset; tau0=tauset; eta=etaset;
	bweight=bweightset;
	resinfoptr=cb->reslist->GetResInfoPtr(ID);
	if(resinfoptr->decay==false){
		msquared=resinfoptr->mass*resinfoptr->mass;
	}
	r[3]=tau0*sinh(eta);
	r[0]=tau0*cosh(eta);
	et=sqrt(p[1]*p[1]+p[2]*p[2]+msquared);
	p[3]=et*sinh(y);
	Setp0();
	if(tau0<0.0){
		printf("FATAL: tau0<0, tau0^2=%g\n",tau0);
		Print();
		exit(1);
	}
	balanceID=balanceIDset;
}

void CHBPart::Print(){
	printf("________________ PART INFO: balanceID=%d  _____________________________\n",balanceID);
	printf("ID=%d, m_onshell=%g, M=%g, tau0=%g=?%g\n r=(%g,%g,%g,%g) eta=%g=?%g, bweight=%g\n", 
	resinfo->code,resinfo->mass,sqrt(msquared),double(tau0),sqrt(r[0]*r[0]-r[3]*r[3]),r[0],r[1],r[2],r[3],eta,GetEta(tau0),bweight);
	printf("p=(%15.9e,%15.9e,%15.9e,%15.9e), y=%g =? %g\n",p[0],p[1],p[2],p[3],y,atanh(p[3]/p[0]));
	printf("________________________________________________________________________\n");
}

double CHBPart::GetMass(){
	if(resinfo->code==22)
		return 0.0;
	else
		return sqrt(msquared);
}

/* in sampler there is an array of densities and temperature. also make array of densities of minmass*/

void CHBPart::Setp0(){
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+msquared);
}

void CHBPart::SetY(){
	y=asinh(p[3]/GetMT());
}

void CHBPart::SetMass(){
	msquared=p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
}

double CHBPart::GetEta(double tau){
	double dy,deta,dtau0,dtau;
	dy=y;
	deta=eta;
	dtau0=tau0;
	dtau=tau;
	deta=dy-asinh((dtau0/dtau)*sinh(dy-deta));
	return deta;
}

double CHBPart::GetMT(){
	return sqrt(p[0]*p[0]-p[3]*p[3]);
}

void CHBPart::FindDecay(){
	double t,gamma,vz,newt,newz,taudecay;
	t=HBARC/resinfo->width;
	gamma=p[0]/sqrt(msquared);
	t=-t*gamma*log(cb->randy->ran());
	vz=p[3]/p[0];
	newt=r[0]+t;
	newz=r[3]+vz*t;
	taudecay=sqrt(newt*newt-newz*newz);
}

double CHBPart::GetPseudoRapidity(){
	double pmag,eta_ps;
	pmag=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	eta_ps=atanh(p[3]/pmag);
	return eta_ps;
}

double CHBPart::GetRapidity(){
	return 0.5*log((p[0]+p[3])/(p[0]-p[3]));
}

void CHBPart::Boost(FourVector &u){
	BoostP(u);
	BoostR(u);
}

void CHBPart::BoostP(FourVector &u){
	Misc::Boost(u,p);
	int mu;
	double ndotp,udotn,udotp;
	
	ndotp=p[0];
	udotn=u[0];
	udotp=u[0]*p[0]-u[1]*p[1]-u[2]*p[2]-u[3]*p[3];
	p[0]=p[0]-((udotp+ndotp)/(1.0+udotn));
	for(mu=0;mu<4;mu++){
		p[mu]=-((udotp+ndotp)/(1.0+udotn))*u[mu]+2*ndotp*u[mu]+p[mu];
	}
}

void CHBPart::BoostR(FourVector &u){
	int mu;
	double ndotr,udotn,udotr;
	ndotr=r[0];
	udotn=u[0];
	udotr=u[0]*r[0]-u[1]*r[1]-u[2]*r[2]-u[3]*r[3];
	r[0]=r[0]-((udotr+ndotr)/(1.0+udotn));
	for(mu=0;mu<4;mu++){
		r[mu]=-((udotr+ndotr)/(1.0+udotn))*u[mu]+2*ndotr*u[mu]+r[mu];
	}
}

CHBPart* CBalance::GetNewPart(mapip *pmap,int balanceid){
	int ipart;
	CHBPart *part;
	if (deadpartlist.size()==0){
		for(ipart=0;ipart<500;ipart++){
			part=new CHBPart();
			part->partmap=NULL;
			deadpartlist.push_back(part);
		}
	}
	part=deadpartlist.front();
	part->balanceID=balanceid;
	part->partmap=pmap;
	pmap->insert(pairip(balanceid,part));
	deadpartlist.pop_front();
	
	return part;
}

void CHBPart::Kill(){
	pair<mapip::iterator,mapip::iterator> itpair;
	mapip::iterator it,it1,it2;
	itpair=partmap->equal_range(balanceID);
	it1=itpair.first;
	it2=itpair.second;
	it=it1;
	while(it!=it2 && it->second!=this){
		++it;
	}
	if(it->second!=this){
		printf("CHBPart::Kill() -- Can't find particle in map to kill!!!!!\n");
		Print();
		exit(1);
	}
	partmap->erase(it);
	partmap=NULL;
	cb->deadpartlist.push_back(this);
}

void CHBPart::BoostRap(double dely){
	double gamma=cosh(dely),gammav=sinh(dely);
	double r0=r[0],p0=p[0];
	y=y+dely;
	eta=eta+dely;
	r[0]=gamma*r0+gammav*r[3];
	r[3]=gamma*r[3]+gammav*r0;
	p[0]=gamma*p0+gammav*p[3];
	p[3]=gamma*p[3]+gammav*p0;
}

void CHBPart::CalcDCA(double *dca){
	int alpha;
	char nantestc[20];
	string nantests;
	double pdotr,p2;
	p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
	pdotr=(p[1]*r[1]+p[2]*r[2]+p[3]*r[3])/p2;
	for(alpha=1;alpha<4;alpha++)
		dca[alpha]=(r[alpha]-pdotr*p[alpha])/1.0E13;
	dca[0]=sqrt(dca[1]*dca[1]+dca[2]*dca[2]+dca[3]*dca[3]);
	sprintf(nantestc,"%g",r[0]);
	nantests=nantestc;
	if(nantests=="NaN" || nantests=="nan" || nantests=="inf" || nantests=="INF"){
		printf("::: dca=(%g,%g,%g,%g)\n",dca[0],dca[1],dca[2],dca[3]);
		printf("::: r=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
		printf("::: p=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
		exit(1);
	}
}

#endif
