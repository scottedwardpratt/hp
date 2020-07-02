#include "acceptance.h"
#include "resonances.h"
#include "part.h"
#include "parametermap.h"
#include "misc.h"

using namespace std;

CAcceptance_ALICE::CAcceptance_ALICE(CparameterMap *parmapin) : CAcceptance(){
	ETAMIN=-0.9; // Don't bother calling Acceptance Routine if outside these boundaries.
	ETAMAX=0.9;
	PTMIN=200.0; 
	PTMAX=2500.0;
	parmap=parmapin;    
	CENTRALITY=parmap->getI("ALICE_CENTRALITY",0);   // CENTRALITY=0 is most central
}

void CAcceptance_ALICE::CalcAcceptance(bool &accept,double &efficiency,CPart *part){
	double pt,*p=part->p;
	double dca[4],dcaxy;
	int pid=part->resinfo->code;
	
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	accept=false;
	efficiency=0.0;
	if(pt>PTMIN && pt<PTMAX && part->resinfo->charge!=0){
		part->CalcDCA(dca);
		if(fabs(dca[3])<2.0){
			dcaxy=sqrt(dca[1]*dca[1]+dca[2]*dca[2]);
			if(abs(pid)==211){
				if(pt>200.0 && pt<2000.0 && dcaxy<0.04){ 
					accept=true; efficiency=1.0;
				}
			}
			else if(abs(pid)==321){
				if(pt>200.0 && pt<2000.0 && dcaxy<2.0){
					accept=true; efficiency=1.0;
				}
			}
			else if(abs(pid)==2212){
				if(pt>500.0 && pt<2500 && dcaxy<0.04){
					accept=true; efficiency=1.0;
				}
			}
		}
	}
}

void CAcceptance_ALICE::CalcAcceptance_Realistic(bool &accept,double &efficiency,CPart *part){
	double pt,*p=part->p,y=part->y;
	double dca[4],dcaxy;
	int pid=part->resinfo->code;
	printf("this is dead code in acceptance_ALICE.cc, shouldn't be here\n");
	exit(1);
	
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	//pmag=sqrt(pt*pt+p[3]*p[3]);
	//eta=atanh(p[3]/pmag);
	accept=false;
	efficiency=0.0;
	if(pt>PTMIN && pt<PTMAX && part->resinfo->charge!=0){
		part->CalcDCA(dca);
		if(fabs(dca[3])<2.0){
			dcaxy=sqrt(dca[1]*dca[1]+dca[2]*dca[2]);
			if(abs(pid)==211){
				if(pt<2000.0 && dcaxy<0.04 && fabs(y)<0.8){ // for cross-species fabs(y)<0.7
					accept=true; efficiency=1.0;
				}
			}
			else if(abs(pid)==321){
				if(pt<2000.0 && dcaxy<2.0 && fabs(y)<0.7){
					accept=true; efficiency=1.0;
				}
			}
			else if(abs(pid)==2212){
				if(pt>500.0 && pt<2500 && dcaxy<0.04 && fabs(y)<0.7){  // for pp BF, it was required that fabs(y)<0.7
					accept=true; efficiency=1.0;
				}
			}
		}
	}
}

void CAcceptance_ALICE::CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part){
	double eta,pt,pmag,*p=part->p;
	//double dca[4],dcaxy;
	//int pid=part->resinfo->code;
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	accept=false;
	efficiency=0.0;
	if(pt>PTMIN && pt<PTMAX && part->resinfo->charge!=0){
		accept=true; efficiency=1.0;
	}
}

double CAcceptance_ALICE::GetDelYMax(int pida,int pidb){
	if(abs(pida)==211 && abs(pidb)==211){
		return 1.6;
	}
	else if(abs(pida)==2212 && abs(pidb)==2212){
		return 1.2;
	}
	else
		return 1.4;
}
