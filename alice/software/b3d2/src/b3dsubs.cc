#ifndef __B3DSUBS_CC__
#define __B3DSUBS_CC__

#include "b3d.h"
#include "part.h"
#include "misc.h"
#include "resonances.h"
#include "cell.h"
#include "constants.h"
#include "randy.h"

void CB3D::WriteAnnihilationData(){
	if(BARYON_ANNIHILATION){
		double total=0.0;
		int itau,imax=lrint(TAUCOLLMAX);
		for(int itau=0;itau<imax;itau++){
			printf("%6.2f %g\n",(itau+0.5),annihilation_array[itau]);
			total+=annihilation_array[itau];
		}
		printf("%g total annihilations, nbaryons=%d, annihilation fraction=%g\n",total,nbaryons,2.0*total/double(nbaryons));
	}
}

void CB3D::PerformAllActions(){
	if(DENSWRITE){
		for(int itau=0;itau<DENSWRITE_NTAU;itau++){
			AddAction_DensCalc((itau+1.0)*DENSWRITE_DELTAU);
		}
	}
	CAction *action;
	nscatter=nbscatter=ndecay=npass=nmerge=nswallow=npass=nexit=nactivate=ninelastic=ncheck=nactionkills=nbaryons=0;
	ncollisions=nannihilate=nregenerate=0;
	tau=0.0;
	nactions=0;	
	CActionMap::iterator epos=ActionMap.begin();
	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Perform();
		epos=ActionMap.begin();
	}
}

void CB3D::KillAllActions(){
	CAction *action;
	CActionMap::iterator epos=ActionMap.begin();
	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Kill();
		epos=ActionMap.begin();
	}
	nactions=0;
}

void CB3D::KillAllParts(){
	CPart *part;
	CPartMap::iterator ppos;
	
	ppos=PartMap.begin();
	while(ppos!=PartMap.end()){
		part=ppos->second;
		if(part->currentmap!=&PartMap){
			printf("KillAllParts:  currentpart not listed as PartMap\n");
			part->Print();
			printf("PartMap.size=%d, DeadPartMap.size=%d\n",int(PartMap.size()),int(DeadPartMap.size()));
			part->currentmap=&PartMap;
			//exit(1);
			Misc::Pause();
		}
		part->Kill();
		ppos=PartMap.begin();
	}
	
	// Recheck 
	for(ppos=DeadPartMap.begin();ppos!=DeadPartMap.end();++ppos){
		part=ppos->second;
		if(part->currentmap!=&DeadPartMap){
			printf("particle in dead part map has wrong current map\n");
			part->Print();
			printf("PartMap.size=%d, DeadPartMap.size=%d\n",int(PartMap.size()),int(DeadPartMap.size()));
			exit(1);
		}
	}

	int ix,iy,ieta;
	if(COLLISIONS){
		CPartMap *partmap;
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				for(ieta=0;ieta<2*NETA;ieta++){
					partmap=&(cell[ix][iy][ieta]->partmap);
					ppos=partmap->begin();
					while(ppos!=partmap->end()){
						part=ppos->second;
						part->currentmap=&DeadPartMap;
						partmap->erase(ppos);
						//part->Kill();
						ppos=partmap->begin();
					}
				}
			}
		}
	}
}

void CB3D::PrintActionMap(CActionMap *actionmap){
	CActionMap::iterator epos;
	CAction *action;
	int iaction=0;
	printf("_________________ ACTIONMAP %d actions _________________________\n",int(actionmap->size()));
	for(epos=actionmap->begin();epos!=actionmap->end();++epos){
		iaction+=1;
		action=epos->second;
		printf("iaction=%d : ",iaction);
		action->Print();
	}
}

void CB3D::FindAllCollisions(){
	double taucoll;
	CPartMap::iterator ppos1,ppos2;
	CPart *part1,*part2;
	CActionMap::iterator epos;
	CAction *action;
	for(ppos1=PartMap.begin();ppos1!=PartMap.end();++ppos1){
		part1=ppos1->second;
		part1->KillActions();
	}
	ppos1=PartMap.begin();
	part1=ppos1->second;
	ppos2=ppos1; ++ppos2;
	while(ppos2!=PartMap.end()){
		part1=ppos1->second; part2=ppos2->second;
		if(part1->balanceID<0 || part2->balanceID<0)
			FindCollision(part1,part2,taucoll);
		ppos1=ppos2;
		++ppos2;
	}
}

void CB3D::PrintPartList(){
	CPartMap::iterator ppos2,ppos1=PartMap.begin();
	while(ppos1!=PartMap.end()){
		printf("%d ",ppos1->second->listid);
		ppos2=ppos1; ++ppos2;
		if(ppos2!=PartMap.end()){
			if(ppos1->second->actionmother!=ppos2->second->actionmother) printf("| ");
		}
		++ppos1;
	}
	printf("\n");
}

void CB3D::PrintMuTInfo(){
	int ix,iy,itau,ntau=lrint(TAUCOLLMAX/CMuTInfo::DELTAU);
	double Tpi,mupi,uxpi,uypi,TK,muK,uxK,uyK;
	double tau,estarpi,estarK;
	char filename[40];
	FILE *fptr;
	CMuTInfo *mti;
	iy=NXY;
	for(itau=0;itau<ntau;itau++){
		tau=(itau+1)*MUTCALC_DELTAU;
		sprintf(filename,"mucalc_results/mutinfo_tau%g.dat",tau);
		fptr=fopen(filename,"w");
		fprintf(fptr,"#   x     Npi     E*pi     Tpi     mupi    uxpi    NK      EK*      TK     muK     uxK\n");
		for(ix=0;ix<2*NXY;ix++){
			mti=muTinfo[ix][iy];
			mti->FindMuTInfo_pi(itau);
			mti->FindMuTInfo_K(itau);
			estarpi=sqrt(mti->Epi[itau]*mti->Epi[itau]-mti->Pxpi[itau]*mti->Pxpi[itau]-mti->Pypi[itau]*mti->Pypi[itau])/mti->Npi[itau];
			estarK=sqrt(mti->EK[itau]*mti->EK[itau]-mti->PxK[itau]*mti->PxK[itau]-mti->PyK[itau]*mti->PyK[itau])/mti->NK[itau];
			fprintf(fptr,"%6.2f %7d %7.2f %7.2f %7.2f %7.4f %7d %7.2f %7.2f %7.2f %7.4f\n",
			DXY*(-NXY+ix+0.5),mti->Npi[itau],estarpi,mti->Tpi,mti->mupi,mti->uxpi,
			mti->NK[itau],estarK,mti->TK,mti->muK,mti->uxK);
		}
		fclose(fptr);
	}
}

void CB3D::ListFutureCollisions(){
	CActionMap::iterator epos=ActionMap.begin();
	CAction *action;
	CPartMap::iterator p1,p2;
	printf("------------------- LIST OF FUTURE COLLISIONS ---------------------\n");
	while(epos!=ActionMap.end()){
		action=epos->second;
		if(action->type==2){
			p1=action->partmap.begin();
			p2=p1; ++p2;
			printf("%d  %d  will collide at %g\n",p1->second->listid,p2->second->listid,double(action->tau));
		}
		epos++;
	}
}

double CB3D::CalcSigma(CPart *part1,CPart *part2){
	double sigma=SIGMADEFAULT;
	CPartMap::iterator ppos;
	CPart *part3,*part4;
	const double g[4]={1,-1,-1,-1};
	double sigma_annihilation,sigma_inel,Gamma,G,G2,MR,M,b,q2,q3,q4,qR2;
	double tan2delta;
	double mt,P[4],q[4],P2,Pdotq;
	const int NWMAX=5000;
	double inel_weight[NWMAX]={0.0};
	double inel_d=0.0,q_prime;
	int ir1,ir2,irflip,alpha,G_Value,L_merge, pmq=0, pmb=0, pms=0;
	CMerge *merge;
	list<CInelasticInfo>::iterator inel;
	list<CInelasticInfo> inel_list;
	bool G_Parity = false;
	int netb=0,netq=0,nets=0,ndaughters;
	int qpions,iK,ipair,npaircheck;
	double Plab,p1dotp2;
	int itau;
	double jR,j1,j2,j1_i,j2_i,rstrange;

	ir1=part1->resinfo->ires; ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1; ir1=ir2; ir2=irflip;
	}
	merge=reslist->MergeArray[ir1][ir2];
	if(merge!=NULL || (BARYON_ANNIHILATION && (part1->resinfo->baryon*part2->resinfo->baryon)<0) || INELASTIC){
		p1dotp2=0.0;
		for(alpha=0;alpha<4;alpha++){
			p1dotp2+=part1->p[alpha]*part2->p[alpha]*g[alpha];
		}
		P2=part1->msquared+part2->msquared+2.0*p1dotp2;
		q2=Misc::triangle2(P2,part1->msquared,part2->msquared);
		M=sqrt(P2);
	}
	else{
		return sigma;
	}

	// Annihilation
	if(BARYON_ANNIHILATION && (part1->resinfo->baryon*part2->resinfo->baryon)<0){
		Plab=sqrt((p1dotp2*p1dotp2/(part2->msquared))-part1->msquared);
		sigma_annihilation=6.7*pow(Plab/1000.0,-0.7)/double(NSAMPLE);
		rstrange=0.5*sqrt(sigma_annihilation);
		rstrange*=pow(ANNIHILATION_SREDUCTION,abs(part1->resinfo->strange))+pow(ANNIHILATION_SREDUCTION,abs(part2->resinfo->strange));
		sigma+=rstrange*rstrange;
	}

	//Calculate quantities used for both inel scattering and merging
	if(merge!=NULL || inel!=inel_list.end()){
		j1=part1->resinfo->spin;
		j2=part2->resinfo->spin;
	}

	//Check for merging
	while(merge!=NULL){
		Gamma=merge->resinfo->width;
		b=merge->branching;
		jR=merge->resinfo->spin;
		MR=merge->resinfo->mass;
		L_merge = merge->L;
		qR2=Misc::triangle2(MR*MR,part1->msquared,part2->msquared);
		q3=pow(q2/qR2,(2*L_merge + 1)/2);
		q4=pow(q2/qR2,(2*L_merge)/2);
		G=Gamma*(MR/M)*q3*1.2/(1.0+0.2*q4);
		tan2delta=pow(0.5*G/(M-MR),2);

		sigma+=b*((4.0*PI*HBARC*HBARC/q2)*(tan2delta/(1.0+tan2delta))
			*((2.0*jR+1.0)/((2.0*j1+1.0)*(2.0*j2+1.0))));
		merge=merge->next;
	}
	
	//Check for Inelastic Scatering
	//inel_d = (2.0*j1+1.0)*(2.0*j2+1.0)*q2;
	if(INELASTIC){
		if(part1->resinfo->G_Parity && part2->resinfo->G_Parity){
			G_Parity = true;
			G_Value = part1->resinfo->G_Parity * part2->resinfo->G_Parity;
		}
		// First calculate denominator
		int iw;
		if(inelasticlist->UseInelasticArray){
			inel_list = inelasticlist->InelasticArray[ir1][ir2];
		}else{
			netq = part1->resinfo->charge+part2->resinfo->charge;
			netb = part1->resinfo->baryon+part2->resinfo->baryon;
			nets = part1->resinfo->strange+part2->resinfo->strange;
			inel_list = inelasticlist->ThermalArray[abs(netb)][abs(netq)][abs(nets)][Misc::Sign(netb)][Misc::Sign(netq)][Misc::Sign(nets)];
		}
		inel_d = Q0*Q0;
		inel = inel_list.begin();
		iw=0;
		while(inel!=inel_list.end() && (M>inel->min_mass)){
			if((G_Parity && (inel->resinfo_1->G_Parity * inel->resinfo_2->G_Parity == G_Value)) || (!G_Parity)){
				j1_i=inel->resinfo_1->spin;
				j2_i=inel->resinfo_2->spin;
				if(inel->resinfo_1->mass+inel->resinfo_2->mass<M){
					//if(inel->resinfo_1->baryon==0 || inel->resinfo_2->baryon==0){
					q_prime = Misc::triangle(M,inel->resinfo_1->mass, inel->resinfo_2->mass);
					inel_weight[iw] = (2.0*j1_i+1.0)*(2.0*j2_i+1.0)*q_prime;
					//}
				}
				inel_d += inel_weight[iw];
			}
			iw+=1;
			if(iw==NWMAX){
				printf("MUST INCREASE NWMAX in int CB3D::Collide\n");
				exit(1);
			}
			inel++;

			// now thumb through
			inel = inel_list.begin();
			iw=0;
			while(inel!=inel_list.end() && (M>inel->min_mass)){
				if((G_Parity && (inel->resinfo_1->G_Parity * inel->resinfo_2->G_Parity == G_Value)) || !G_Parity){
					sigma+=SIGMAINELASTIC*inel_weight[iw]/(inel_d*double(NSAMPLE));
				}
				iw+=1;
				inel++;
			}
		}
	}
	return sigma;
}

// Note part2 is fake and part1 is real
void CB3D::SplitPart(CPart *part1,CPart *part2){
	double oldeta,mt,g1,g2;
	CB3DCell *cell;
	part2->Copy(part1); // does not change reality or weights
	if(BJORKEN){
		oldeta=part1->eta;
		part1->eta=-ETAMAX+2.0*ETAMAX*randy->ran();
		part1->y+=(part1->eta-oldeta);
		part1->r[3]=part1->tau0*sinh(part1->eta);
		part1->r[0]=sqrt(part1->tau0*part1->tau0+part1->r[3]*part1->r[3]);
		mt=sqrt(part1->msquared+part1->p[1]*part1->p[1]+part1->p[2]*part1->p[2]);
		part1->p[3]=mt*sinh(part1->y);
		part1->Setp0();
	}
	else{
		randy->ran_gauss2(g1,g2);
		part1->r[1]+=0.5*g1;
		part1->r[2]+=0.5*g2;
		oldeta=part1->eta;
		part1->eta+=0.5*randy->ran_gauss()/tau;
		part1->y+=(part1->eta-oldeta);
		part1->r[3]=part1->tau0*sinh(part1->eta);
		part1->r[0]=sqrt(part1->tau0*part1->tau0+part1->r[3]*part1->r[3]);
		mt=sqrt(part1->msquared+part1->p[1]*part1->p[1]+part1->p[2]*part1->p[2]);
		part1->p[3]=mt*sinh(part1->y);
		part1->Setp0();
	}
	
	cell=part1->FindCell();
	part1->ChangeCell(cell);
	if(part1->currentmap!=&PartMap)
		part1->ChangeMap(&PartMap);
	
	cell=part2->FindCell();
	part2->ChangeCell(cell);
	if(part2->currentmap!=&PartMap)
		part2->ChangeMap(&PartMap);
}

CPart* CB3D::GetDeadPart(){
	CPart *newpart;
	if(DeadPartMap.size()==0){
		for(int ipart=0;ipart<DELNPARTSTOT*NSAMPLE;ipart++)
			newpart=new CPart(npartstot);
		printf("made new parts, npartstot=%d, tau=%g\n",npartstot,tau);
	}
	return DeadPartMap.begin()->second;
}

void CB3D::GetDeadParts(CPart *&part1,CPart *&part2){
	int ipart;
	CPart *newpart;
	while(DeadPartMap.size()<2){
		for(ipart=0;ipart<DELNPARTSTOT*NSAMPLE;ipart++)
			newpart=new CPart(npartstot);
		printf("made new parts, npartstot=%d, tau=%g\n",npartstot,tau);
	}
	CPartMap::iterator ppos=DeadPartMap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
}

void CB3D::GetDeadParts(array<CPart*,5> &product){
	int ipart;
	CPart *newpart;
	while(DeadPartMap.size()<5){
		for(ipart=0;ipart<DELNPARTSTOT*NSAMPLE;ipart++)
			newpart=new CPart(npartstot);
		printf("made new parts, npartstot=%d, tau=%g\n",npartstot,tau);
	}
	CPartMap::iterator ppos=DeadPartMap.begin();
	for(ipart=0;ipart<5;ipart++){
		product[ipart]=ppos->second;
		++ppos;
	}
}

CAction* CB3D::GetDeadAction(){
	CAction *action;
	if(DeadActionMap.size()==0){
		CPartMap::iterator ppos=PartMap.begin();
		int nreal=0;
		CPart *part;
		while(ppos!=PartMap.end()){
			part=ppos->second;
			++ppos;
		}
		for(int iaction=0;iaction<DELNACTIONSTOT*NSAMPLE;iaction++)
			action=new CAction(nactionstot);
		printf("created new actions, nactionstot=%d\n",nactionstot);
	}
	return DeadActionMap.begin()->second;
}

void CB3D::CheckPartMap(){
	CPartMap::iterator iter;
	CPart *part;
	for(iter=PartMap.begin();iter!=PartMap.end();iter++){
		part=iter->second;
		if(part->currentmap!=&PartMap){
			printf("----- FAILED CheckPartMap-----\n");
			part->Print();
			exit(1);
		}
	}
}


#endif