#include "b3d.h"
#include "part.h"
#include "resonances.h"
#include "commondefs.h"
#include "sampler.h"
#include "randy.h"
#include "misc.h"
#include "constants.h"
int ReadParts(string filename,CPartMap &partmap,CResList *reslist);

using namespace std;

int ReadParts(string filename,CPartMap &partmap,CResList *reslist){
	FILE *oscarfile;
	CB3DBinaryPartInfo bpart;
	CPart *part;
	int nparts_read,ipart=0,ievent_read;
	double et;
	CResInfo *resinfo;
	oscarfile=fopen(filename.c_str(),"rb");
	fread(&ievent_read,sizeof(int),1,oscarfile);
	fread(&nparts_read,sizeof(int),1,oscarfile);
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		fread(&bpart,sizeof(bpart),1,oscarfile);
		part=new CPart();

		resinfo=reslist->GetResInfoPtr(bpart.ID);
		part->resinfo=resinfo;
		part->p[1]=bpart.px;
		part->p[2]=bpart.py;
		part->y=bpart.rapidity;
		part->r[1]=bpart.x;
		part->r[2]=bpart.y;
		part->tau0=bpart.tau;
		part->eta0=bpart.eta;
		part->msquared=resinfo->mass*resinfo->mass;
		et=sqrt(part->msquared+part->p[1]*part->p[1]+part->p[2]*part->p[2]);
		part->p[3]=et*sinh(part->y);
		part->p[0]=sqrt(et*et+part->p[3]*part->p[3]);
		part->phi0=atan2(part->p[2],part->p[1]);
		part->phi0*=180.0/PI;
		part->weight=part->bweight=1.0;
		partmap.insert(CPartPair(part->y,part));
	}
	fclose(oscarfile);
	return nparts_read;
}

int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("model_output/fixed_parameters.dat");
	int nparts0,nparts1;
	CResList *reslist=new CResList(&parmap);
	CPartMap partmap0,partmap1;
	string filename0="model_output/default_0/alice_cent0_5/oscar.dat";
	string filename1="model_output/default_1/alice_cent0_5/oscar.dat";
	nparts0=ReadParts(filename0,partmap0,reslist);
	printf("nparts0=%d\n",nparts0);
	nparts1=ReadParts(filename1,partmap1,reslist);
	printf("nparts1=%d\n",nparts1);
	
	CPartMap::iterator iter0,iter1;
	CPart *part0,*part1;
	double qinv,qout,qout_lcms,qlong,qside,delphi,dely,deleta;
	double DELQ=2.5,DELY=0.1,DELPHI;
	const int NQMAX=200,NY=50,NPHI=18;
	long long int NPIONS=0;
	int ibin;
	DELPHI=180.0/NPHI;
	double mixed_inv[NQMAX]={0.0},mixed_out[NQMAX]={0.0},mixed_side[NQMAX]={0.0},mixed_long[NQMAX]={0.0};
	double mixed_y[NY]={0.0},mixed_phi[NPHI]={0.0};
	
	for(iter0=partmap0.begin();iter0!=partmap0.end();iter0++){
		part0=iter0->second;
		if(abs(part0->resinfo->code)==211){
			NPIONS+=1;
			for(iter1=partmap1.begin();iter1!=partmap1.end();iter1++){
				part1=iter1->second;
				if(abs(part1->resinfo->code)==211){
					Misc::outsidelong_lcms(part0->p,part1->p,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi);
					ibin=lrint(floor(dely/DELY));
					if(ibin<NY){
						mixed_y[ibin]+=1.0;
					}
					ibin=lrint(floor(delphi/DELPHI));
					if(ibin<NPHI){
						mixed_phi[ibin]+=1.0;
					}
					ibin=lrint(floor(qinv/DELQ));
					if(ibin<NQMAX)
						mixed_inv[ibin]+=1.0;
					ibin=lrint(floor(fabs(qout/DELQ)));
					if(ibin<NQMAX)
						mixed_out[ibin]+=1.0;
					ibin=lrint(floor(fabs(qside/DELQ)));
					if(ibin<NQMAX)
						mixed_side[ibin]+=1.0;
					ibin=lrint(floor(fabs(qlong/DELQ)));
					if(ibin<NQMAX)
						mixed_long[ibin]+=1.0;
				}
			}
		}
	}
	printf("NPIONS=%lld\n",NPIONS);
	printf("--- vs rapidity ---\n");
	for(ibin=0;ibin<NY;ibin++){
		printf("%5.2f %g\n",(ibin+0.5)*DELY,mixed_y[ibin]/double(NPIONS));
	}
	printf("--- vs phi ---\n");
	for(ibin=0;ibin<NPHI;ibin++){
		printf("%5.2f %g\n",(ibin+0.5)*DELPHI,mixed_phi[ibin]/double(NPIONS));
	}
	printf("--- vs qout ---\n");
	for(ibin=0;ibin<NQMAX;ibin++){
		printf("%5.2f %g\n",(ibin+0.5)*DELQ,mixed_out[ibin]/double(NPIONS));
	}
	printf("--- vs qside ---\n");
	for(ibin=0;ibin<NQMAX;ibin++){
		printf("%5.2f %g\n",(ibin+0.5)*DELQ,mixed_side[ibin]/double(NPIONS));
	}
	printf("--- vs qlong ---\n");
	for(ibin=0;ibin<NQMAX;ibin++){
		printf("%5.2f %g\n",(ibin+0.5)*DELQ,mixed_long[ibin]/double(NPIONS));
	}
	
	return 0;
}
