#ifndef __CHARGE_CC__
#define __CHARGE_CC__

#include "b3d.h"
#include "charge.h"
#include "part.h"
#include "sampler.h"
#include "constants.h"
#include "randy.h"

using namespace std;

void CB3D::GenHadronsFromCharges(){
	int maxbid,bid,newcount=0;
	CCharge *chargea,*chargeb;
	pair<CChargeMap::iterator,CChargeMap::iterator> icpair_even,icpair_odd;
	CChargeMap::iterator itc;
	CPartMap::iterator itp;
	itc=chargemap.end(); itc--;
	maxbid=itc->first;

	sampler->cummulative_N=0.0;
	sampler->cummulative_random=randy->ran_exp();
	
	printf("maxbid=%d, number of charges=%d\n",maxbid,int(chargemap.size()));
	for(bid=0;bid<maxbid;bid+=2){
		//printf("-----------\n");
		newcount+=2;
		icpair_even=chargemap.equal_range(bid);
		itc=icpair_even.first;
		chargea=itc->second;
		icpair_odd=chargemap.equal_range(bid+1);
		itc=icpair_odd.first;
		chargeb=itc->second;
		GenHadronsFromCharge(bid,chargea);
		GenHadronsFromCharge(bid+1,chargeb);
	}
}

void CB3D::GenHadronsFromCharge(int balanceID,CCharge *charge){
	int ires;
	double delN,bweight;
	double rapidity,mass;
	CHyperElement *hyper=&(charge->hyper);
	CPart *part;
	FourVector p;
	CChargeMap::iterator it;
	Eigen::Vector3d Qprime,Q,q;
	CResInfoMap::iterator itr;
	CResInfo *resinfo;
	
	Q(0)=charge->q[0];
	Q(1)=charge->q[1];
	Q(2)=charge->q[2];
	Qprime=sampler->chiinvf*Q;
	for(itr=reslist->resmap.begin();itr!=reslist->resmap.end();++itr){
		resinfo=itr->second;
		ires=resinfo->ires;
		if(resinfo->baryon!=0 || resinfo->charge!=0 || resinfo->strange!=0){ 
			q[0]=resinfo->q[0]; q[1]=resinfo->q[1]; q[2]=resinfo->q[2];
			delN=sampler->densityf[ires]*q.dot(Qprime); // number of hadrons to create
			bweight=charge->weight*delN/fabs(delN);
			sampler->cummulative_N+=fabs(delN*NSAMPLE_UDS2BAL);
			while(sampler->cummulative_N>sampler->cummulative_random){
				sampler->cummulative_random+=randy->ran_exp();
				hyper->GetP(resinfo,p,mass,sampler->maxweightf[ires]);
				part=GetDeadPart();
				rapidity=charge->eta+asinh(p[3]/sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]));
				part->InitBalance(resinfo->code,charge->x,charge->y,charge->tau,charge->eta,p[1],p[2],mass,rapidity,bweight,balanceID);
			}
		}
	}
}

void CB3D::AnalyzeCharges(){
	double phi1,phi2,delphi;
	int nphi=18,iphi;
	if(int(phicount.size())!=nphi)
		phicount.resize(nphi);
	CChargeMap::iterator it;
	CCharge *charge1,*charge2;
	for(it=chargemap.begin();it!=chargemap.end();it++){
		charge1=it->second;
		++it;
		charge2=it->second;
		if(charge1->q[2]!=0 && charge2->q[2]!=0){
			phi1=atan2(charge1->y,charge1->x);
			phi2=atan2(charge2->y,charge2->x);
			delphi=fabs(phi1-phi2);
			while (delphi>PI){
				delphi-=2.0*PI;
				delphi=fabs(delphi);
			}
			//printf("delphi=%g\n",delphi);
			iphi=lrint(floor((delphi/PI)*nphi));
			phicount[iphi]-=charge1->q[2]*charge2->q[2];
		}
	}
}

void CB3D::ReadCharges(int ichargefile){
	string dirname="udsdata/"+qualifier;
	char chargefile[10];
	sprintf(chargefile,"%d",ichargefile);
	//string filename=dirname+"/"+parmap.getS("CHARGESINFO_FILENAME","uds.dat");
	string filename=dirname+"/"+"uds"+chargefile+".dat";
	printf("reading charges from %s\n",filename.c_str());
	CHyperElement *hyper;
	double etaspread=0.0;
	int neta=0,maxbid=0;
	char dummy[120];
	vector<double> etaboost;
	CChargeMap::iterator it;
	CCharge *charge;
	int balanceID,qu,qd,qs,bidcharge;
	double u0,ux,uy,x,y,tau_read,eta,w,dOmega0,dOmegaX,dOmegaY,pitildexx;
	double pitildexy,pitildeyy;
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,120,fptr);
	chargemap.clear();
	double norm=0.0,ux2bar=0.0,uy2bar=0.0,udotrbar=0.0,x2bar=0.0,y2bar=0.0;
	do{
		fscanf(fptr,"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&balanceID,&qu,&qd,&qs,&w,&tau_read,&eta,&x,&y,&ux,&uy,&dOmega0,&dOmegaX,&dOmegaY,&pitildexx,&pitildeyy,&pitildexy);
		if(!feof(fptr)){
			charge=new CCharge();
			hyper=&(charge->hyper);
			charge->q[0]=qu;
			charge->q[1]=qd;
			charge->q[2]=qs;
			charge->weight=w;
			charge->tau=tau_read;
			charge->eta=eta;
			if(fabs(charge->q[2])==1){
				etaspread+=eta*eta;
				neta+=1;
			}
			charge->x=x;
			charge->y=y;
			hyper->tau=tau_read;
			hyper->x=x;
			hyper->y=y;
			hyper->dOmega0=dOmega0;
			hyper->dOmegaX=dOmegaX;
			hyper->dOmegaY=dOmegaY;
			u0=sqrt(1.0+ux*ux+uy*uy);
			hyper->ux=ux;
			hyper->uy=uy;
			hyper->udotdOmega=u0*dOmega0-ux*dOmegaX-uy*dOmegaY;
			hyper->pitildexx=1000.0*pitildexx;
			hyper->pitildeyy=1000.0*pitildeyy;
			hyper->pitildexy=1000.0*pitildexy;
			hyper->T=sampler->Tf;
			hyper->P=sampler->Pf;
			hyper->epsilon=sampler->epsilonf;
			hyper->h=hyper->P+hyper->epsilon;
			hyper->lambda=sampler->lambdaf;
			chargemap.insert(CChargePair(balanceID,charge));
			norm+=1.0;
			ux2bar+=ux*ux;
			uy2bar+=uy*uy;
			x2bar+=x*x;
			y2bar+=y*y;
			udotrbar+=ux*x+uy*y;
			if(balanceID>maxbid)
				maxbid=balanceID;
		}
	}while(!feof(fptr));
	printf("Read %d Charges\n",int(chargemap.size()));
	printf("etaspread=%g\n",sqrt(etaspread/double(neta)));
	printf("<x^2>=%g, <y^2>=%g, <ux^2>=%g, <uy^2>=%g, <u.r>=%g\n",
	x2bar/norm,y2bar/norm,ux2bar/norm,uy2bar/norm,udotrbar/norm);
	fclose(fptr);
	
	etaboost.resize((maxbid+1)/2);
	for(bidcharge=0;bidcharge<(maxbid+1)/2;bidcharge+=1){
		etaboost[bidcharge]=ETAMAX*(1.0-2.0*randy->ran());
	}
	for(it=chargemap.begin();it!=chargemap.end();++it){
		charge=it->second;
		balanceID=it->first;
		bidcharge=floorl(balanceID/2);
		charge->eta+=etaboost[bidcharge];
		while(charge->eta>ETAMAX){
			charge->eta-=2.0*ETAMAX;
		}
		while(charge->eta<-ETAMAX){
			charge->eta+=2.0*ETAMAX;
		}
	}
	etaboost.clear();
	//CalcChiTotFromQ();
}

void CB3D::DeleteCharges(){
	CChargeMap::iterator it,itnext;
	CCharge *charge;
	it=chargemap.begin();
	while(it!=chargemap.end()){
		itnext=it;
		++itnext;
		charge=it->second;
		if(charge!=NULL)
			delete charge;
		it=itnext;
	}
	chargemap.clear();
}

void CB3D::IncrementChiTotFromCharges(){
	pair<CChargeMap::iterator,CChargeMap::iterator> icpair_even,icpair_odd;
	CChargeMap::iterator itc;
	int a,b,bid,maxbid;
	int NSAMPLE_HYDRO2UDS=parmap.getI("NSAMPLE_HYDRO2UDS",1);
	Eigen::Vector3d qa;
	Eigen::Vector3d qb;
	itc=chargemap.end(); itc--;
	maxbid=itc->first;
	CCharge *chargea,*chargeb;
	for(bid=0;bid<maxbid;bid+=2){
		icpair_even=chargemap.equal_range(bid);
		itc=icpair_even.first;
		chargea=itc->second;
		icpair_odd=chargemap.equal_range(bid+1);
		itc=icpair_odd.first;
		chargeb=itc->second;
		for(a=0;a<3;a++){
			for(b=0;b<3;b++){
				chitotQ(a,b)-=0.5*(chargea->q[a]*chargeb->q[b]+chargea->q[b]*chargeb->q[a]);
			}
		}
	}
	printf("-- ChiTot From Charges\n");
	cout << chitotQ/(NSAMPLE_HYDRO2UDS*TotalVolume) <<endl;
}

void CB3D::IncrementChiTotFromHadrons(){
	CPartMap bfpartmap;
	int maxbid=-1;
	pair<CPartMap::iterator,CPartMap::iterator> itpair_even,itpair_odd;
	CPartMap::iterator it,ita0,itaf,itb0,itbf,ita,itb;
	CPart *parta,*partb,*part;
	for(it=PartMap.begin();it!=PartMap.end();++it){
		part=it->second;
		if(part->balanceID>maxbid)
			maxbid=part->balanceID;
		if(part->balanceID>=0)
			bfpartmap.insert(CPartPair(part->balanceID,part));
	}
	printf("maxbid=%d, bfpartmap.size=%d\n",maxbid,int(bfpartmap.size()));
	
	int bid,count=0;
	int NSAMPLE_HYDRO2UDS=parmap.getI("NSAMPLE_HYDRO2UDS",1);
	double dchi;
	
	for(bid=0;bid<=maxbid;bid+=2){
		itpair_even=bfpartmap.equal_range(bid);
		ita0=itpair_even.first;
		itaf=itpair_even.second;
		itpair_odd=bfpartmap.equal_range(bid+1);
		itb0=itpair_odd.first;
		itbf=itpair_odd.second;
		if(ita0!=itaf && itb0!=itbf){
			for(ita=ita0;ita!=itaf;++ita){
				parta=ita->second;
				if(parta->balanceID!=bid){
					printf("parta bids don't match %d !=%d\n",parta->balanceID,bid);
					exit(1);
				}
				for(itb=itb0;itb!=itbf;++itb){
					partb=itb->second;
					if(partb->balanceID!=bid+1){
						printf("parta bids don't match %d !=%d\n",partb->balanceID,bid+1);
						exit(1);
					}
					for(int a=0;a<3;a++){
						for(int b=0;b<3;b++){
							dchi=0.5*parta->resinfo->q[a]*partb->resinfo->q[b]*parta->bweight*partb->bweight
								+0.5*parta->resinfo->q[b]*partb->resinfo->q[a]*parta->bweight*partb->bweight;
							chitotH(a,b)-=0.5*dchi;
							chitotH(b,a)-=0.5*dchi;
							//printf("bids=%d, %d\n",parta->balanceID,partb->balanceID);
							if(parta->resinfo->q[a]*partb->resinfo->q[b]!=0)
								count+=1;
						}
					}
				}
			}
		}
	}
	bfpartmap.clear();
	printf("-- ChiTot from Hadrons --  count=%d\n",count);
	cout << chitotH/(NSAMPLE_HYDRO2UDS*TotalVolume*NSAMPLE_UDS2BAL*NSAMPLE_UDS2BAL) << endl;
}

#endif