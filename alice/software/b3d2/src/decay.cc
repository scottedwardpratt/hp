#ifndef __DECAY_CC__
#define __DECAY_CC__

#include "b3d.h"
#include "randy.h"
#include "resonances.h"
#include "part.h"
#include "misc.h"
#include "constants.h"

void CB3D::Decay(CPart *mother,int &nbodies,array<CPart *,5> &daughter){
	const double HBARC=197.326;
	int ibody,jbody,alpha;
	double mass[6],mtot,mprime,wmaxmass,wmass,mguess,kmass,kmaxmass,kguess;
	CPart *dptr;

	FourVector *p[6],kprime,qprime,ptot,pprime,u12,pp,u;
	double q,weight,wmax,sthet,cthet,phi;
	double p3mag,kprimemax,p3max,ppmax,kprimemax2,kprimemag2,qprimemax,qprimemax2,qprimemag2,ppmag;
	double e1prime,e2prime,e3prime,e4prime,e1max,e2max,e3max,e4max,e12;

	mass[0]=mother->GetMass();
	p[0]=&mother->p;
	
	/* Create daughter objects */
	mtot=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		mass[ibody+1]=daughter[ibody]->resinfo->mass;
	}
	for(ibody=0;ibody<nbodies;ibody++){
		if(daughter[ibody]->resinfo->decay){
			//generate mass according to density of states, ~ rho(m)*k*E1*E2
			mprime=0;
			for(jbody=0;jbody<nbodies;jbody++){
				if(jbody!=ibody)
					mprime+=mass[jbody+1];
			}
			kmaxmass=pow(mass[0],4)+pow(mprime,4)-2.0*mass[0]*mass[0]*mprime*mprime;
			kmaxmass=0.5*sqrt(kmaxmass)/mass[0];
			wmaxmass=kmaxmass*kmaxmass*sqrt(kmaxmass*kmaxmass+mprime*mprime);
			do{
				mguess=daughter[ibody]->resinfo->GenerateMass();
				if(mass[0]>mguess+mprime){
					kguess=pow(mass[0],4)+pow(mprime,4)+pow(mguess,4)-2.0*mass[0]*mass[0]*mprime*mprime
						-2.0*mass[0]*mass[0]*mguess*mguess-2.0*mguess*mguess*mprime*mprime;
					kguess=0.5*sqrt(kguess)/mass[0];
					wmass=kguess*sqrt(mprime*mprime+kguess*kguess)*sqrt(mguess*mguess+kguess*kguess);
				}
				else
					wmass=-1.0;
			}while(wmass<0.0 || randy->ran()>wmass/wmaxmass);
			if(wmass>wmaxmass){
				printf("In  CB3D::Decay, wmass=%g > wmaxmass=%g\n",wmass,wmaxmass);
				printf("kguess=%g, mguess=%g, mprime=%g, E1=%g, E2=%g, E1+E2=%g=?%g\n",kguess,mguess,mprime,
				sqrt(kguess*kguess+mguess*mguess),sqrt(kguess*kguess+mprime*mprime),
				sqrt(kguess*kguess+mguess*mguess)+sqrt(kguess*kguess+mprime*mprime),mass[0]);
				exit(1);
			}
			mass[ibody+1]=mguess;
		}
		else{
			mass[ibody+1]=daughter[ibody]->resinfo->mass;
		}
		mtot+=mass[ibody+1];
		p[ibody+1]=&daughter[ibody]->p;
	}

	/* TWO-BODY DECAYS */
	if(nbodies==2){
		cthet=1.0-2.0*randy->ran();
		sthet=sqrt(1.0-cthet*cthet);
		phi=2.0*PI*randy->ran();
		q=sqrt(Misc::triangle(mass[0],mass[1],mass[2]));
		(*p[1])[3]=q*cthet;
		(*p[1])[1]=q*sthet*cos(phi);
		(*p[1])[2]=q*sthet*sin(phi);
		(*p[2])[3]=-(*p[1])[3];
		(*p[2])[2]=-(*p[1])[2];
		(*p[2])[1]=-(*p[1])[1];
		(*p[1])[0]=sqrt(mass[1]*mass[1]+(*p[1])[1]*(*p[1])[1]+(*p[1])[2]*(*p[1])[2]+(*p[1])[3]*(*p[1])[3]);
		(*p[2])[0]=sqrt(mass[2]*mass[2]+(*p[2])[1]*(*p[2])[1]+(*p[2])[2]*(*p[2])[2]+(*p[2])[3]*(*p[2])[3]);
	}
	/* THREE-BODY DECAYS */
	else if(nbodies==3){
		kprimemax2=Misc::triangle(mass[0]-mass[3],mass[1],mass[2]);
		kprimemax=sqrt(kprimemax2);
		p3max=sqrt(Misc::triangle(mass[0],mass[1]+mass[2],mass[3]));
		e1max=sqrt(pow(mass[1],2)+p3max*p3max);
		e2max=sqrt(pow(mass[2],2)+p3max*p3max);
		e3max=sqrt(pow(mass[3],2)+p3max*p3max);
		//wmax=p3max*(e1max*e2max/(mass[1]*mass[2]))*(mass[1]+mass[2])/(e1max+e2max);
		wmax=p3max*pow(e1max+e2max,2)*e3max/(mass[1]+mass[2]);
		do{
			TRY_AGAIN:
			do{
				kprime[1]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[2]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[3]=kprimemax*(2.0*randy->ran()-1.0);
				kprimemag2=kprime[1]*kprime[1]+
					kprime[2]*kprime[2]+kprime[3]*kprime[3];
			} while(kprimemag2>kprimemax2);
			e1prime=sqrt(kprimemag2+mass[1]*mass[1]);
			e2prime=sqrt(kprimemag2+mass[2]*mass[2]);
			if(e1prime+e2prime+mass[3]>mass[0]) goto TRY_AGAIN;
			p3mag=sqrt(Misc::triangle(mass[0],e1prime+e2prime,mass[3]));
			cthet=1.0-2.0*randy->ran();
			sthet=sqrt(1.0-cthet*cthet);
			phi=2.0*PI*randy->ran();
			(*p[3])[3]=p3mag*cthet;
			(*p[3])[1]=p3mag*sthet*cos(phi);
			(*p[3])[2]=p3mag*sthet*sin(phi);
			(*p[3])[0]=sqrt(p3mag*p3mag+mass[3]*mass[3]);
			e12=sqrt(pow(e1prime+e2prime,2)+p3mag*p3mag);
			for(alpha=1;alpha<4;alpha++)
				u12[alpha]=-(*p[3])[alpha]/(e1prime+e2prime);
			u12[0]=sqrt(1.0+u12[1]*u12[1]+u12[2]*u12[2]+u12[3]*u12[3]);
			kprime[0]=e1prime;
			Misc::lorentz(u12,kprime,*p[1]);
			kprime[0]=e2prime;
			for(alpha=1;alpha<=3;alpha++) kprime[alpha]=-kprime[alpha];
			Misc::lorentz(u12,kprime,(*p[2]));
			weight=p3mag*pow((*p[1])[0]+(*p[2])[0],2)*(*p[3])[0]/(e1prime+e2prime);
		} while(randy->ran()>weight/wmax);
	}
	/* FOUR-BODY DECAYS */
	else if(nbodies==4){
		kprimemax2=Misc::triangle(mass[0]-mass[3]-mass[4],mass[1],mass[2]);
		kprimemax=sqrt(kprimemax2);
		qprimemax2=Misc::triangle(mass[0]-mass[1]-mass[2],mass[3],mass[4]);
		qprimemax=sqrt(qprimemax2);
		
		ppmax=sqrt(Misc::triangle(mass[0],mass[1]+mass[2],mass[3]+mass[4]));
		e1max=sqrt(pow(mass[1],2)+ppmax*ppmax);
		e2max=sqrt(pow(mass[2],2)+ppmax*ppmax);
		e3max=sqrt(pow(mass[3],2)+ppmax*ppmax);
		e4max=sqrt(pow(mass[4],2)+ppmax*ppmax);
		wmax=ppmax*pow(e1max+e2max,2)*pow(e3max+e4max,2)/((mass[1]+mass[2])*(mass[3]+mass[4]));
		
		do{
			TRY_AGAIN_4:
			do{
				kprime[1]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[2]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[3]=kprimemax*(2.0*randy->ran()-1.0);
				kprimemag2=kprime[1]*kprime[1]+kprime[2]*kprime[2]+kprime[3]*kprime[3];
			} while(kprimemag2>kprimemax2);
			e1prime=sqrt(kprimemag2+mass[1]*mass[1]);
			e2prime=sqrt(kprimemag2+mass[2]*mass[2]);
			do{
				qprime[1]=qprimemax*(2.0*randy->ran()-1.0);
				qprime[2]=qprimemax*(2.0*randy->ran()-1.0);
				qprime[3]=qprimemax*(2.0*randy->ran()-1.0);
				qprimemag2=qprime[1]*qprime[1]+qprime[2]*qprime[2]+qprime[3]*qprime[3];
			} while(qprimemag2>qprimemax2);
			e3prime=sqrt(qprimemag2+mass[3]*mass[3]);
			e4prime=sqrt(qprimemag2+mass[4]*mass[4]);
			
			if(e1prime+e2prime+e3prime+e4prime>mass[0]) goto TRY_AGAIN_4;

			ppmag=Misc::triangle(mass[0],e1prime+e2prime,e3prime+e4prime);
			if(ppmag>0){
				ppmag=sqrt(ppmag);
				cthet=1.0-2.0*randy->ran();
				sthet=sqrt(1.0-cthet*cthet);
				phi=2.0*PI*randy->ran();
				pp[3]=ppmag*cthet;
				pp[1]=ppmag*sthet*cos(phi);
				pp[2]=ppmag*sthet*sin(phi);

				pp[0]=sqrt(ppmag*ppmag+pow(e1prime+e2prime,2));
				for(alpha=0;alpha<4;alpha++)
					u[alpha]=pp[alpha]/(e1prime+e2prime);
				kprime[0]=sqrt(mass[1]*mass[1]+kprimemag2);
				Misc::lorentz(u,kprime,*p[1]);
				kprime[0]=sqrt(mass[2]*mass[2]+kprimemag2);
				for(alpha=1;alpha<4;alpha++)
					kprime[alpha]=-kprime[alpha];
				Misc::lorentz(u,kprime,*p[2]);

				for(alpha=1;alpha<4;alpha++)
					pp[alpha]=-pp[alpha];
				pp[0]=sqrt(ppmag*ppmag+pow(e3prime+e4prime,2));
				for(alpha=0;alpha<4;alpha++)
					u[alpha]=pp[alpha]/(e3prime+e4prime);
				qprime[0]=sqrt(mass[3]*mass[3]+qprimemag2);
				Misc::lorentz(u,qprime,(*p[3]));
				qprime[0]=sqrt(mass[4]*mass[4]+qprimemag2);
				for(alpha=1;alpha<4;alpha++)
					qprime[alpha]=-qprime[alpha];
				Misc::lorentz(u,qprime,*p[4]);

				weight=ppmag*pow((*p[1])[0]+(*p[2])[0],2)*pow((*p[3])[0]+(*p[4])[0],2)/((e1prime+e2prime)*(e3prime+e4prime));
			}
			else weight=0.0;
		} while(randy->ran()>weight/wmax);
	}

	/* Boost the new particles */
	for(alpha=0;alpha<4;alpha++) u[alpha]=mother->p[alpha]/mass[0];
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		Misc::lorentz(u,*p[ibody+1],pprime);
		dptr->active=true;
		for(alpha=0;alpha<4;alpha++)
			dptr->p[alpha]=pprime[alpha];
		dptr->CopyPositionInfo(mother);
		dptr->msquared=pow(dptr->resinfo->mass,2);
		dptr->Setp0();
		dptr->SetY();
		dptr->eta0=mother->eta0;
		dptr->phi0=mother->phi0;
	}
}

#endif
