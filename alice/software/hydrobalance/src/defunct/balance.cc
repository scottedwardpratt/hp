#include "commondefs.h"
#include "part.h"
#include "balance.h"
#include "constants.h"

using namespace std;

CBalance *CPart::cb=NULL;

CBalance::CBalance(CHydroBalance *hb){
	parmap=(hb->parmap);
	qualifier=hb->qualifier;
	Init();
}

CBalance::CBalance(string parfilename,string qualifier_set){
	qualifier=qualifier_set;
	parmap.ReadParsFromFile(parfilename);
	Init();
}

void CBalance::Init(){
	CResInfoMap::iterator itr;
	CResInfo *resinfo;
	randy=new CRandy(time(NULL));
	CPart::cb=this;
	CResList::cb=this;
	ReadResonanceList();
	Tf=parmap.getD("FREEZEOUT_TEMP",140);
	NSAMPLE_HYDRO2UDS=parmap.getD("NSAMPLE_HYDRO2UDS",4);
	NSAMPLE_UDS2BAL=parmap.getD("NSAMPLE_UDS2BAL",4);
	densityf.resize(int(reslist->resmap.size()));
	maxweightf.resize(int(reslist->resmap.size()));
	chif.setZero();
	chiinvf.setZero();
	reslist->CalcEoSandChi(Tf,Pf,epsilonf,nf,densityf,maxweightf,chif);
	printf("chif\n");
	cout << chif << endl;
	double Q2=(4*chif(0,0)+chif(1,1)+chif(2,2)-2*chif(0,1)
		-2*chif(1,0)-2*chif(0,2)-2*chif(2,0)+chif(1,2)+chif(2,1))/9.0;
	printf("chiQQ=%g\n",Q2);
	lambdaf=reslist->GetLambda(Tf,epsilonf,Pf);
	ransum=0.0;
	ranthresh=randy->ran_exp();
	Nerror=0.0;
	gammap=v2=normtest=0.0;
	Nchi=parmap.getI("NCHI",200);
	CalcChiInv(-1.0,chiinvf);
	CreateBFArrays();
	badmothercount=0;
	TotalVolume=0.0;
	acceptance_description=parmap.getS("BF_ACCEPTANCE","PERFECT");
	if(acceptance_description=="PERFECT"){
		acceptance=new CAcceptance_CHEAP(&parmap);
		bf_results_dirname="results_perfect";
	}
	else if(acceptance_description=="STAR"){
		acceptance=new CAcceptance_STAR(&parmap);
		bf_results_dirname="results_star";
	}
	else if(acceptance_description=="ALICE"){
		acceptance=new CAcceptance_ALICE(&parmap);
		bf_results_dirname="results_alice";
	}
	else{
		printf("Define BF_ACCEPTANCE in parameters.dat\n");
		exit(1);
	}
}

void CBalance::ReadResonanceList(){
	reslist=new CResList(&parmap);
}

void CBalance::GenDecayProcessHadronsFromCharges(){
	int maxbid,bid,newcount=0,dummy;
	CCharge *chargea,*chargeb;
	pair<mapic::iterator,mapic::iterator> icpair_even,icpair_odd;
	// for testing
	Eigen::Matrix3d chitesth; chitesth.setZero();
	CPart *parta,*partb;
	CResInfo *resinfoa,*resinfob;
	mapip::iterator it,ita0,itaf,itb0,itbf,ita,itb;
	pair<mapip::iterator,mapip::iterator> itpair_even,itpair_odd;
	int a,b;
	//
	mapic::iterator itc;
	mapip::iterator itp;
	CPart *part;
	itc=emap.end(); itc--;
	maxbid=itc->first;
	printf("maxbid=%d\n",maxbid);
	double meanpt;
	int npt=0;
	for(bid=0;bid<maxbid;bid+=2){
		partmap.clear();
		//printf("-----------\n");
		newcount+=2;
		icpair_even=emap.equal_range(bid);
		itc=icpair_even.first;
		chargea=itc->second;
		icpair_odd=emap.equal_range(bid+1);
		itc=icpair_odd.first;
		chargeb=itc->second;
		GenHadronsFromCharge(bid,chargea);
		GenHadronsFromCharge(bid+1,chargeb);

		// For testing
		int na=0,nb=0;
		itpair_even=partmap.equal_range(bid);
		itpair_odd=partmap.equal_range(bid+1);
		ita0=itpair_even.first;
		itaf=itpair_even.second;
		itb0=itpair_odd.first;
		itbf=itpair_odd.second;
		if(ita0!=itaf && itb0!=itbf){
			for(ita=ita0;ita!=itaf;++ita){
				parta=ita->second;
				resinfoa=parta->resinfo;
				na+=1;
				nb=0;
				for(itb=itb0;itb!=itbf;++itb){
					partb=itb->second;
					resinfob=partb->resinfo;
					nb+=1;
					for(a=0;a<3;a++){
						for(b=0;b<3;b++)
							chitesth(a,b)-=0.5*(resinfoa->q[a]*resinfob->q[b]+resinfoa->q[b]*resinfob->q[a])*parta->bweight*partb->bweight;
					}
				}
			}
		}
		//printf("na=%d, nb=%d, sum =? %d\n",na,nb,int(partmap.size()));
		
		DecayParts(partmap);
		for(itp=partmap.begin();itp!=partmap.end();++itp){
			part=itp->second;
			if(abs(part->resinfo->code)==211){
				meanpt+=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				npt+=1;
			}
		}
		ProcessPartMap();
		DeletePartMapParts();
		if(newcount*10>maxbid){
			printf("finished %ld percent of PartMap Processing, bid=%d\n",
			lrint(bid*100.0/double(maxbid)),bid);
			newcount=0;
		}
	}
	printf("<pt>=%g\n",meanpt/double(npt));
	// for testing
	printf("chitesth from hadrons (before decays)\n");
	chitesth=chitesth/(TotalVolume*double(NSAMPLE_HYDRO2UDS*NSAMPLE_UDS2BAL*NSAMPLE_UDS2BAL));
	printf("-----chitesth -----\n");
	cout << chitesth << endl;
	double Q2=(4*chitesth(0,0)+chitesth(1,1)+chitesth(2,2)-2*chitesth(0,1)
		-2*chitesth(1,0)-2*chitesth(0,2)-2*chitesth(2,0)+chitesth(1,2)+chitesth(2,1))/9.0;
	printf("chiQQ=%g\n",Q2);
	//
	DeleteEMapCharges();
}

void CBalance::GenHadronsFromCharge(int balanceID,CCharge *charge){
	int a,b,ires,Nres,ntries,ntrymax=100000,iv;
	double ux,uy,u0,delN,bweight,dOmega0,dOmegaX,dOmegaY;
	double rapidity,delv,vOmega,w,dNbar;
	Eigen::Vector3d Qtest;
	double ntest,mass;
	int Nsample;
	CHyperElement *hyper=charge->hyper;
	hyper->CalcDOmegaMax();
	vOmega=charge->hyper->vOmega;
	//if(vOmega<-1.0)
	//hyper->Print(); // for testing
	if(vOmega>0.0){
		printf("vOmega>0.0, =%g\n",vOmega);
		exit(1);
	}
	CPart *part;
	FourVector p;
	Nres=reslist->resmap.size();
	ux=hyper->ux;
	uy=hyper->uy;
	u0=sqrt(1.0+ux*ux+uy*uy);
	mapic::iterator it;
	Eigen::Vector3d Qprime,Q,q,Qcheck;
	CResInfoMap::iterator itr;
	CResInfo *resinfo;
	
	Q(0)=charge->q[0];
	Q(1)=charge->q[1];
	Q(2)=charge->q[2];
	//Qtest(0)=Qtest(1)=Qtest(2)=0.0;
	Qprime=chiinvf*Q;
	for(itr=reslist->resmap.begin();itr!=reslist->resmap.end();++itr){
		resinfo=itr->second;
		ires=resinfo->ires;
		if(resinfo->baryon!=0 || resinfo->charge!=0 || resinfo->strange!=0){
			q(0)=resinfo->q[0]; q(1)=resinfo->q[1]; q(2)=resinfo->q[2];
			delN=densityf[ires]*q.dot(Qprime); // number of hadrons to create
			bweight=delN/fabs(delN);
			ransum+=fabs(delN*NSAMPLE_UDS2BAL);
			//Qtest+=delN*q;
			while(ransum>ranthresh){
				ranthresh+=randy->ran_exp();
				hyper->GetP(resinfo,p,mass,maxweightf[ires]);
				part=GetNewPart(&partmap,balanceID);
				rapidity=charge->eta+asinh(p[3]/sqrt(mass*mass
					+p[1]*p[1]+p[2]*p[2]));
				part->InitBalance(resinfo->code,charge->x,charge->y,charge->tau,charge->eta,
				p[1],p[2],mass,rapidity,bweight,balanceID);
			}
		}
	}
	//printf("Q=(%g,%g,%g) =? (%g,%g,%g)\n",Q(0),Q(1),Q(2),Qtest(0),Qtest(1),Qtest(2));
}

void CBalance::GenHadronsFromHyperElements(){
	list<CHyperElement *>::iterator it;
	int nhyper=0;
	for(it=hyperlist.begin();it!=hyperlist.end();++it){
		GenHadronsFromHyperElement(*it);
		nhyper+=1;
	}
	hyperlist.clear();
	printf("N Hyper-elements=%d\n",nhyper);
	CalcChiTotFromH();
}

void CBalance::GenHadronsFromHyperElement(CHyperElement *hyper){
	int a,b,ires,Nres,balanceID;
	double ux,uy,u0,delN,eta,rapidity,mass;
	double nsample=NSAMPLE_UDS2BAL*NSAMPLE_HYDRO2UDS;
	CPart *part;
	FourVector p;
	Nres=reslist->resmap.size();
	ux=hyper->ux;
	uy=hyper->uy;
	u0=sqrt(1.0+ux*ux+uy*uy);
	mapic::iterator it;
	CResInfoMap::iterator itr;
	CResInfo *resinfo;
	if(ransum+nf*hyper->udotdOmega*nsample>ranthresh){
		for(itr=reslist->resmap.begin();itr!=reslist->resmap.end();++itr){
			resinfo=itr->second;
			ires=resinfo->ires;
			if(resinfo->code!=22){
				delN=densityf[ires]*hyper->udotdOmega;
				ransum+=nsample*fabs(delN);
				while(ransum>ranthresh){
					ranthresh+=randy->ran_exp();
					hyper->GetP(resinfo,p,mass,maxweightf[ires]);
					balanceID=ncpartmap.size();
					part=GetNewPart(&ncpartmap,balanceID);
					//eta=YMIN+(YMAX-YMIN)*randy->ran();
					eta=0.0;
					rapidity=eta+asinh(p[3]/sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]));
					if(hyper->tau<0.0){
						printf("hyper->tau<0, =%g\n",hyper->tau);
						hyper->Print();
						exit(1);
					}
					part->InitBalance(resinfo->code,hyper->x,hyper->y,hyper->tau,eta,p[1],p[2],mass,rapidity,1.0,balanceID);
				}
			}
		}
	}
	else
		ransum+=nf*hyper->udotdOmega*nsample;
	delete hyper;
}

void CBalance::ProcessPartMap(){
	int balanceID,maxbid,pida,pidb;
	pair<mapip::iterator,mapip::iterator> itpair_even,itpair_odd;
	mapip::iterator it,ita0,itaf,itb0,itbf,ita,itb;
	it=partmap.end(); it--;
	maxbid=it->first;
	CPart *parta,*partb;
	double w=1.0/NSAMPLE_UDS2BAL;
	double meanpt=0.0;
	int npt=0;

	for(balanceID=0;balanceID<maxbid;balanceID+=2){
		itpair_even=partmap.equal_range(balanceID);
		itpair_odd=partmap.equal_range(balanceID+1);
		ita0=itpair_even.first;
		itaf=itpair_even.second;
		itb0=itpair_odd.first;
		itbf=itpair_odd.second;
		if(ita0!=itaf && itb0!=itbf){
			for(ita=ita0;ita!=itaf;++ita){
				parta=ita->second;
				pida=parta->resinfo->code;
				if(parta->resinfo->charge!=0){
					if(abs(pida)!=211 && abs(pida)!=321 && abs(pida)!=2212){
						printf("mystery pida=%d\n",pida);
						exit(1);
					}
					for(itb=itb0;itb!=itbf;++itb){
						partb=itb->second;
						pidb=partb->resinfo->code;
						if(partb->resinfo->charge!=0){
							if(abs(pidb)!=211 && abs(pidb)!=321 && abs(pidb)!=2212){
								printf("mystery pidb=%d\n",pidb);
								exit(1);
							}
							IncrementNumer(parta,partb,w);
						}
					}
				}
			}
		}
	}
}

void CBalance::ProcessNCPartMap(){
	int balanceID,maxbid;
	double w;
	mapip::iterator it,itb;
	CPart *part,*partb;
	pair<mapip::iterator,mapip::iterator> itpair;
	it=ncpartmap.end(); it--;
	maxbid=it->first;
	long long int count=0,newcount=0;
	for(balanceID=0;balanceID<=maxbid;balanceID++){
		count+=1;
		newcount+=1;
		itpair=ncpartmap.equal_range(balanceID);
		for(it=itpair.first;it!=itpair.second;++it){
			part=it->second;
			part->bweight=1.0;
			if(abs(part->resinfo->charge)==1){
				w=1.0;
				IncrementDenom(part,w);
				for(itb=itpair.first;itb!=it;++itb){
					partb=itb->second;
					if(abs(partb->resinfo->charge)==1){
						partb->bweight=1.0;
						w=1.0;
						IncrementNumer(part,partb,w);
						IncrementNumer(partb,part,w);
					}
				}
			}
		}
		if(newcount*10>maxbid || count==maxbid-1){
			printf("finished %ld percent of NCPartMap Processing\n",
			lrint(count*100.0/double(maxbid+1)));
			newcount=0;
		}
	}
}

void CBalance::IncrementDenom(CPart *part,double w){
	int pid;
	bool accepta;
	double effa,dely,ya,phia;
	CPart parta;
	pid=part->resinfo->code;
	if(abs(pid)!=211 && abs(pid)!=321 && abs(pid)!=2212){
		printf("Hmmm, pid=%d, %s\n",pid,part->resinfo->name.c_str());
		part->resinfo->Print();
		exit(1);
	}
	parta.Copy(part);
	ya=acceptance->ETAMIN+randy->ran()*(acceptance->ETAMAX-acceptance->ETAMIN);
	dely=ya-parta.y;
	parta.BoostRap(dely);
	acceptance->CalcAcceptance(accepta,effa,&parta);
	phia=atan2(parta.p[2],parta.p[1]);
	v2+=w*cos(2.0*phia);
	if(accepta && randy->ran()<effa){	
		if(abs(pid)==211)
			denom_pi->fw=
				&parta,w);
		if(abs(pid)==321)
			denom_K->Increment(&parta,w);
		if(abs(pid)==2212)
			denom_p->Increment(&parta,w);
	}
	acceptance->CalcAcceptanceNoID(accepta,effa,&parta);
	if(accepta && randy->ran()<effa){	
		denom_charged->Increment(&parta,w);
		phia=atan2(parta.p[2],parta.p[1]);
		phia=phia*180.0/PI;
		if(phia>180.0)
			phia=phia-180.0;
		if(phia<0.0)
			phia=-phia;
		if(phia>90.0)
			phia=180-phia;
		if(phia<15.0)
			denom_charged_phi0->Increment(&parta,w);
		if(phia>30.0 && phia<60.0)
			denom_charged_phi45->Increment(&parta,w);
		if(phia>75.0)
			denom_charged_phi90->Increment(&parta,w);
	}
}

void CBalance::IncrementNumer(CPart *parta,CPart *partb,double w){
	double effa,effb,ya,dely,phia;
	bool accepta,acceptb;
	int pida,pidb;
	CPart partaa,partbb;
	if(parta->badmother==0 || partb->badmother==0 ||
	(parta->badmother!=partb->badmother)){
		partaa.Copy(parta);
		ya=acceptance->ETAMIN+randy->ran()*(acceptance->ETAMAX-acceptance->ETAMIN);
		dely=ya-partaa.y;
		partaa.BoostRap(dely);
		acceptance->CalcAcceptance(accepta,effa,&partaa);
		if(accepta && randy->ran()<effa){
			partbb.Copy(partb);
			partbb.BoostRap(dely);
			acceptance->CalcAcceptance(acceptb,effb,&partbb);
			if(acceptb && randy->ran()<effb){
				pida=parta->resinfo->code;
				pidb=partb->resinfo->code;
				if(abs(pida)==211 && abs(pidb)==211)
					numer_pipi->Increment(parta,partb,w);
				if( (abs(pida)==211 && abs(pidb)==321) 
					|| (abs(pidb)==211 && abs(pida)==321) )
						numer_piK->Increment(parta,partb,w);
				if( (abs(pida)==211 && abs(pidb)==2212)
					|| (abs(pidb)==211 && abs(pida)==2212) )
						numer_pip->Increment(parta,partb,w);
				if(abs(pida)==321 && abs(pidb)==321)
					numer_KK->Increment(parta,partb,w);
				if( (abs(pida)==321 && abs(pidb)==2212)
					|| (abs(pidb)==321 && abs(pida)==2212) )
						numer_Kp->Increment(parta,partb,w);
				if(abs(pida)==2212 && abs(pidb)==2212)
					numer_pp->Increment(parta,partb,w);
			}
		}
		
		if(abs(parta->resinfo->charge)==1 && abs(partb->resinfo->charge)==1){
			partaa.Copy(parta);
			ya=acceptance->ETAMIN+randy->ran()*(acceptance->ETAMAX-acceptance->ETAMIN);
			dely=ya-partaa.y;
			partaa.BoostRap(dely);
			acceptance->CalcAcceptanceNoID(accepta,effa,&partaa);
			if(accepta && randy->ran()<effa){
				partbb.Copy(partb);
				partbb.BoostRap(dely);
				acceptance->CalcAcceptanceNoID(acceptb,effb,&partbb);
				if(acceptb && randy->ran()<effb){
					pida=parta->resinfo->code;
					pidb=partb->resinfo->code;
					numer_charged->Increment(parta,partb,w);
					phia=atan2(parta->p[2],parta->p[1]);
					phia=phia*180.0/PI;
					if(phia>180.0)
						phia=phia-180.0;
					if(phia<0.0){
						phia=-phia;
					}
					if(phia>90.0){
						phia=180-phia;
					}
					if(phia<15.0)
						numer_charged_phi0->Increment(parta,partb,w);
					if(phia>30.0 && phia<60.0)
						numer_charged_phi45->Increment(parta,partb,w);
					if(phia>75.0)
						numer_charged_phi90->Increment(parta,partb,w);
					IncrementGammaP(parta,partb,w);
				}
			}
		}
	}
}

void CBalance::CheckPartMap(mapip *partmap){
	mapip::iterator it;
	int bid;
	CPart *part;
	for(it=partmap->begin();it!=partmap->end();++it){
		part=it->second;
		bid=it->first;
		if(part->balanceID!=bid){
			printf("balanceID=%d, key=%d\n",part->balanceID,bid);
			part->Print();
			exit(1);
		}
	}
}

void CBalance::CalcChiTotFromQ(){
	pair<mapic::iterator,mapic::iterator> icpair_even,icpair_odd;
	Eigen::Matrix3d chitot=;
	mapic::iterator itc;
	int a,b,bid,maxbid;
	int NSAMPLE_HYDRO2UDS=parmap.getI("NSAMPLE_HYDRO2UDS",1);
	Eigen::Vector3d qa;
	Eigen::Vector3d qb;
	itc=emap.end(); itc--;
	maxbid=itc->first;
	CCharge *chargea,*chargeb;
	for(bid=0;bid<maxbid;bid+=2){
		icpair_even=emap.equal_range(bid);
		itc=icpair_even.first;
		chargea=itc->second;
		icpair_odd=emap.equal_range(bid+1);
		itc=icpair_odd.first;
		chargeb=itc->second;
		for(a=0;a<3;a++){
			for(b=0;b<3;b++){
				chitot(a,b)-=0.5*(chargea->q[a]*chargeb->q[b]+chargea->q[b]*chargeb->q[a]);
			}
		}
	}

	chitot=chitot/(NSAMPLE_HYDRO2UDS*TotalVolume);
	printf("Chitot From Charges\n");
	cout << chitot << endl;
	double Q2=(4*chitot(0,0)+chitot(1,1)+chitot(2,2)-2*chitot(0,1)
		-2*chitot(1,0)-2*chitot(0,2)-2*chitot(2,0)+chitot(1,2)+chitot(2,1))/9.0;
	printf("chiQQ=%g\n",Q2);
}

void CBalance::CalcChiTotFromH(){
	mapip::iterator it;
	Eigen::Matrix3d chiratio,chitot;
	chiratio.setZero(); chitot.setZero()
	CResInfo *resinfo;
	int a,b;
	double chiQQ=0.0;
	CPart *part;
	for(it=ncpartmap.begin();it!=ncpartmap.end();++it){
		part=it->second;
		resinfo=part->resinfo;
		chiQQ+=resinfo->charge*resinfo->charge;
		for(a=0;a<3;a++){
			for(b=0;b<3;b++){
				chitot(a,b)+=resinfo->q[a]*resinfo->q[b];
			}
		}
	}
	chiQQ=chiQQ/(TotalVolume*NSAMPLE_HYDRO2UDS*NSAMPLE_UDS2BAL);
	chitot=chitot/(TotalVolume*NSAMPLE_HYDRO2UDS*NSAMPLE_UDS2BAL);
	printf("chif from NC hadrons\n");
	cout << chitot << endl;
	double Q2=(4*chitot(0,0)+chitot(1,1)+chitot(2,2)-2*chitot(0,1)
		-2*chitot(1,0)-2*chitot(0,2)-2*chitot(2,0)+chitot(1,2)+chitot(2,1))/9.0;
	printf("chiQQ=%g =? %g\n",Q2,chiQQ);
}

void CBalance::DeletePartMapParts(){
	mapip::iterator it,itnext;
	CPart *part;
	it=partmap.begin();
	while(it!=partmap.end()){
		itnext=it;
		++itnext;
		part=it->second;
		delete part;
		it=itnext;
	}
	partmap.clear();
}

void CBalance::DeleteEMapCharges(){
	mapic::iterator it,itnext;
	CCharge *charge;
	it=emap.begin();
	while(it!=emap.end()){
		itnext=it;
		++itnext;
		charge=it->second;
		if(charge!=NULL)
			delete charge;
		it=itnext;
	}
	emap.clear();
}

void CBalance::DeleteNCPartMapParts(){
	mapip::iterator it,itnext;
	CPart *part;
	it=ncpartmap.begin();
	while(it!=ncpartmap.end()){
		itnext=it;
		++itnext;
		part=it->second;
		delete part;
		it=itnext;
	}
	ncpartmap.clear();
}

void CBalance::DeleteHyperElements(){
	list<CHyperElement *>::iterator it,itnext;
	CHyperElement *hyper;
	
	it=hyperlist.begin();
	while(it!=hyperlist.end()){
		itnext=it;
		++itnext;
		hyper=*it;
		delete hyper;
		it=itnext;
	}
	printf("readyto clear\n");
	hyperlist.clear();
}

CPart* CBalance::GetNewPart(mapip *pmap,int balanceid){
	int ipart;
	CPart *part;
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

void CBalance::WriteChiTrans(){
	bool WDECAY=true;
	int a,b,itau,ires,dires,dpid,itb;
	double branching;
	double densityp=0.0,densityK=0.0,densitypi=0.0;
	FILE *fptr;
	CResInfoMap::iterator it;
	CResInfo *resinfo,*dresinfo;
	CBranchInfo *binfo;
	Eigen::Vector3d q;
	q.setZero();
	Eigen::Vector3d qpi={1,-1,0};
	Eigen::Vector3d qK={1,0,-1};
	Eigen::Vector3d qp={2,1,0};
	Eigen::Vector3d dqtilde; qtilde.setZero();
	Eigen::Vector3d qtildepi; qtildepi.setZero();
	Eigen::Vector3d qtildeK; qtildeK.setZero();
	Eigen::Vector3d qtildep; qtildep.setZero();
	
	for(it=reslist->resmap.begin();it!=reslist->resmap.end();++it){
		resinfo=it->second;
		ires=resinfo->ires;
		q(0)=resinfo->q[0];
		q(1)=resinfo->q[1];
		q(2)=resinfo->q[2];
		dqtilde=chiinvf*q*densityf[ires];
		if(!resinfo->decay){
			if(resinfo->code==211){
				qtildepi+=dqtilde;
				densitypi+=densityf[ires];
			}
			if(resinfo->code==321){
				qtildeK+=dqtilde;
				densityK+=densityf[ires];
			}
			if(resinfo->code==2212){
				qtildep+=dqtilde;
				densityp+=densityf[ires];
			}
		}
		else{
			if(WDECAY){
				for(itb=0;itb<resinfo->branchlist.size();itb++){
					branching=resinfo->branchlist[itb]->branching;
				
					for(dires=0;dires<resinfo->branchlist[itb]->resinfoptr.size();dires++){
						dresinfo=resinfo->branchlist[itb]->resinfoptr[dires];
						if(dresinfo->code==211){
							qtildepi+=dqtilde*branching;
							densitypi+=densityf[ires]*branching;
						}
						if(dresinfo->code==321){
							qtildeK+=dqtilde*branching;
							densityK+=densityf[ires]*branching;
						}
						if(dresinfo->code==2212){
							qtildep+=dqtilde*branching;
							densityp+=densityf[ires]*branching;
						}
					}	
				}
			}
		}
	}
	
	bool SCALINGBYDENS=false;
	if(SCALINGBYDENS){
		qtildepi=qtildepi/densitypi;
		qtildeK=qtildeK/densityK;
		qtildep=qtildep/densityp;
	}
	
	
	string dirname="data";
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename;
	printf("----- qtildepi -----\n");
	cout << qtildepi << endl;
	printf("----- qtildeK -----\n");
	cout << qtildeK << endl;
	printf("----- qtildep -----\n");
	cout << qtildep << endl;
	
	qtildepi*=2.0; qtildeK*=2.0; qtildep*=2.0;
	
	
	if(WDECAY)
		filename=dirname+"/chitrans_wdecays.dat";
	else
		filename=dirname+"/chitrans.dat";
	fptr=fopen(filename.c_str(),"w");
	
	fprintf(fptr,"#         uu          dd         ud          du             us           su         ds         sd         ss\n");
	
	fprintf(fptr,"0.5 %11.4e %11.4e %11.4e %11.4e %11.4e  %11.4e %11.4e %11.4e %11.4e\n",
	qtildepi(0)*qtildepi(0),qtildepi(1)*qtildepi(1),
	qtildepi(0)*qtildepi(1),qtildepi(1)*qtildepi(0),
	qtildepi(0)*qtildepi(2),qtildepi(2)*qtildepi(0),
	qtildepi(1)*qtildepi(2),qtildepi(2)*qtildepi(1),
	qtildepi(2)*qtildepi(2));
	
	fprintf(fptr,"1.5 %11.4e %11.4e %11.4e %11.4e %11.4e  %11.4e %11.4e %11.4e %11.4e\n",
	qtildepi(0)*qtildeK(0),qtildepi(1)*qtildeK(1),
	qtildepi(0)*qtildeK(1),qtildepi(1)*qtildeK(0),
	qtildepi(0)*qtildeK(2),qtildepi(2)*qtildeK(0),
	qtildepi(1)*qtildeK(2),qtildepi(2)*qtildeK(1),
	qtildepi(2)*qtildeK(2));
	
	fprintf(fptr,"2.5 %11.4e %11.4e %11.4e %11.4e %11.4e  %11.4e %11.4e %11.4e %11.4e\n",
	qtildepi(0)*qtildep(0),qtildepi(1)*qtildep(1),
	qtildepi(0)*qtildep(1),qtildepi(1)*qtildep(0),
	qtildepi(0)*qtildep(2),qtildepi(2)*qtildep(0),
	qtildepi(1)*qtildep(2),qtildepi(2)*qtildep(1),
	qtildepi(2)*qtildep(2));
	
	fprintf(fptr,"3.5 %11.4e %11.4e %11.4e %11.4e %11.4e  %11.4e %11.4e %11.4e %11.4e\n",
	qtildeK(0)*qtildeK(0),qtildeK(1)*qtildeK(1),
	qtildeK(0)*qtildeK(1),qtildeK(1)*qtildeK(0),
	qtildeK(0)*qtildeK(2),qtildeK(2)*qtildeK(0),
	qtildeK(1)*qtildeK(2),qtildeK(2)*qtildeK(1),
	qtildeK(2)*qtildeK(2));
	
	fprintf(fptr,"4.5 %11.4e %11.4e %11.4e %11.4e %11.4e  %11.4e %11.4e %11.4e %11.4e\n",
	qtildeK(0)*qtildep(0),qtildeK(1)*qtildep(1),
	qtildeK(0)*qtildep(1),qtildeK(1)*qtildep(0),
	qtildeK(0)*qtildep(2),qtildeK(2)*qtildep(0),
	qtildeK(1)*qtildep(2),qtildeK(2)*qtildep(1),
	qtildeK(2)*qtildep(2));
	
	fprintf(fptr,"5.5 %11.4e %11.4e %11.4e %11.4e %11.4e  %11.4e %11.4e %11.4e %11.4e\n",
	qtildep(0)*qtildep(0),qtildep(1)*qtildep(1),
	qtildep(0)*qtildep(1),qtildep(1)*qtildep(0),
	qtildep(0)*qtildep(2),qtildep(2)*qtildep(0),
	qtildep(1)*qtildep(2),qtildep(2)*qtildep(1),
	qtildep(2)*qtildep(2));
	
	fclose(fptr);
	
	if(WDECAY)
		filename=dirname+"/chitrans_short_wdecays.dat";
	else
		filename=dirname+"/chitrans_short.dat";
	fptr=fopen(filename.c_str(),"w");
	
	fprintf(fptr,"#       uu+dd      ud+du   us+su+ds+sd    ss\n");
	
	fprintf(fptr,"0.5 %11.4e %11.4e %11.4e %11.4e\n",
	qtildepi(0)*qtildepi(0)+qtildepi(1)*qtildepi(1),
	qtildepi(0)*qtildepi(1)+qtildepi(1)*qtildepi(0),
	qtildepi(0)*qtildepi(2)+qtildepi(2)*qtildepi(0)+qtildepi(1)*qtildepi(2)+qtildepi(2)*qtildepi(1),
	qtildepi(2)*qtildepi(2));
	
	fprintf(fptr,"1.5 %11.4e %11.4e %11.4e %11.4e\n",
	qtildepi(0)*qtildeK(0)+qtildepi(1)*qtildeK(1),
	qtildepi(0)*qtildeK(1)+qtildepi(1)*qtildeK(0),
	qtildepi(0)*qtildeK(2)+qtildepi(2)*qtildeK(0)+qtildepi(1)*qtildeK(2)+qtildepi(2)*qtildeK(1),
	qtildepi(2)*qtildeK(2));
	
	fprintf(fptr,"2.5 %11.4e %11.4e %11.4e %11.4e\n",
	qtildepi(0)*qtildep(0)+qtildepi(1)*qtildep(1),
	qtildepi(0)*qtildep(1)+qtildepi(1)*qtildep(0),
	qtildepi(0)*qtildep(2)+qtildepi(2)*qtildep(0)+qtildepi(1)*qtildep(2)+qtildepi(2)*qtildep(1),
	qtildepi(2)*qtildep(2));
	
	fprintf(fptr,"3.5 %11.4e %11.4e %11.4e %11.4e\n",
	qtildeK(0)*qtildeK(0)+qtildeK(1)*qtildeK(1),
	qtildeK(0)*qtildeK(1)+qtildeK(1)*qtildeK(0),
	qtildeK(0)*qtildeK(2)+qtildeK(2)*qtildeK(0)+qtildeK(1)*qtildeK(2)+qtildeK(2)*qtildeK(1),
	qtildeK(2)*qtildeK(2));
	
	fprintf(fptr,"4.5 %11.4e %11.4e %11.4e %11.4e\n",
	qtildeK(0)*qtildep(0)+qtildeK(1)*qtildep(1),
	qtildeK(0)*qtildep(1)+qtildeK(1)*qtildep(0),
	qtildeK(0)*qtildep(2)+qtildeK(2)*qtildep(0)+qtildeK(1)*qtildep(2)+qtildeK(2)*qtildep(1),
	qtildeK(2)*qtildep(2));
	
	fprintf(fptr,"5.5 %11.4e %11.4e %11.4e %11.4e\n",
	qtildep(0)*qtildep(0)+qtildep(1)*qtildep(1),
	qtildep(0)*qtildep(1)+qtildep(1)*qtildep(0),
	qtildep(0)*qtildep(2)+qtildep(2)*qtildep(0)+qtildep(1)*qtildep(2)+qtildep(2)*qtildep(1),
	qtildep(2)*qtildep(2));
	
	fclose(fptr);

}
