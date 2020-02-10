//
//  acceptance_CHEAP.cc
//

#include "acceptance.h"
using namespace std;

CAcceptance_CHEAP::CAcceptance_CHEAP(CparameterMap *parmapin) : CAcceptance(){
	ETAMIN=-6;
	ETAMAX=6;
	ptmin=0;
	ptmax=20000000.0;
}

void CAcceptance_CHEAP::CalcAcceptance(bool &accept,double &efficiency,CHBPart *part,int centrality){
	double pt,y,*p=part->p;
	int pid=part->resinfo->code;
	double gammav,pmag,m;
	double eta;
	double ctau_kaon=3.7,ctau_pion=7.8,lmin=1.0;
	double A0;
	
	/* accept=false;
	if(dca[0]<1.5){
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	//printf("pt=%g\n",pt);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	//y=atanh(p[3]/p[0]);
	m=part->resinfo->mass;
	 
	if(eta>ETAMIN && eta<ETAMAX && pt>ptmin && pt<ptmax){
	if(pid==211 || pid==-211){
	A0=0.759;
	accept=true;
	m=sqrt(p[0]*p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_pion));
	}
	if(pid==2212 || pid==-2212){
	if(pid==2212) A0=0.221;
	else A0=0.246;
	accept=true;
	efficiency=A0;
	}
	if(pid==321 || pid==-321){
	A0=0.45;
	accept=true;
	pmag=sqrt(pt*pt+p[3]*p[3]);
	m=sqrt(p[0]*p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_kaon));
	}
	if(pt>600.0) efficiency*=0.65;
	}
	}
	else{
	printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	}
	if(m!=m){
	part->Print();
	exit(1);
	}
	 
	*/
	
	/*
	accept=false;
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	//y=atanh(p[3]/p[0]);
	//if(dca[0]>0.000001) printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	if(pt>ptmin && pt<ptmax && eta>ETAMIN && eta<ETAMAX){
	accept=true;
	efficiency=0.8;
	}
	**/
	
	accept=true;
	efficiency=1.0;
}

void CAcceptance_CHEAP::CalcAcceptanceNoID(bool &accept,double &efficiency,CHBPart *part,int centrality){
	double pt,y,*p=part->p;
	int pid=part->resinfo->code;
	double gammav,pmag,m;
	double eta;
	double ctau_kaon=3.7,ctau_pion=7.8,lmin=1.0;
	double A0;
	
	/* accept=false;
	if(dca[0]<1.5){
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	//printf("pt=%g\n",pt);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	//y=atanh(p[3]/p[0]);
	m=part->resinfo->mass;
	 
	if(eta>ETAMIN && eta<ETAMAX && pt>ptmin && pt<ptmax){
	if(pid==211 || pid==-211){
	A0=0.759;
	accept=true;
	m=sqrt(p[0]*p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_pion));
	}
	if(pid==2212 || pid==-2212){
	if(pid==2212) A0=0.221;
	else A0=0.246;
	accept=true;
	efficiency=A0;
	}
	if(pid==321 || pid==-321){
	A0=0.45;
	accept=true;
	pmag=sqrt(pt*pt+p[3]*p[3]);
	m=sqrt(p[0]*p[0]-pmag*pmag);
	gammav=pmag/m;
	efficiency=A0*exp(-lmin/(gammav*ctau_kaon));
	}
	if(pt>600.0) efficiency*=0.65;
	}
	}
	else{
	printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	}
	if(m!=m){
	part->Print();
	exit(1);
	}
	 
	*/
	
	/*
	accept=false;
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	//y=atanh(p[3]/p[0]);
	//if(dca[0]>0.000001) printf("dca=%g,%g,%g,%g\n",dca[0],dca[1],dca[2],dca[3]);
	if(pt>ptmin && pt<ptmax && eta>ETAMIN && eta<ETAMAX){
	accept=true;
	efficiency=0.8;
	}
	**/
	
	accept=true;
	efficiency=1.0;
	
	
}
