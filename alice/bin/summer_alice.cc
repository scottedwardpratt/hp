#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include "misc.h"

const double PI=2.0*atan(1.0);

using namespace std;

int main(int argc,char *argv[]){
	const int NPHIBINS=28,NYBINS=20;
	const int NTYPES=10,NJDIR=24;
	const double YMAX=1.0;
	string qualifier;
	string acceptance="results_alice";
	string dirprefix="default";
	//string dirprefix="default";
	char qs[200];
	const int NEVENTS=1000;
	char dummy[200];
	
	FILE *fptr_read;
	FILE *fptr_write;
	char oldfilename[100],newfilename[100];
	vector<double> Bysum(NYBINS),Bphisum(NPHIBINS);
	vector<double> Cysum(NYBINS),Cphisum(NPHIBINS);
	double DY=2.0*YMAX/NYBINS,DPHI=360.0/double(NPHIBINS);
	double bf,normy=double(NYBINS)/(NEVENTS*NJDIR*YMAX),normphi=double(NPHIBINS)/(NEVENTS*NJDIR*2.0*PI);
	string olddatafilename;
	int iphi,iy;
	int jdir,itype;
	//printf("Enter qualifier: ");
	//scanf("%s",qs);
	qualifier="alice_cent0_5";
	//qualifier=qs;
	string dirname[NTYPES]={"KK","Kp","piK","pip","pipi","pp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90"};
	char command[120];
	for(itype=0;itype<NTYPES;itype++){
		sprintf(command,"mkdir -p %s_sum/%s/%s/%s",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str(),dirname[itype].c_str());
		system(command);
	}
	
	// Read/sum BF denominators
	olddatafilename="bf_Cyphi.dat";
	for(itype=0;itype<NTYPES;itype++){
		std::fill(Cysum.begin(), Cysum.end(), 0.0);
		std::fill(Cphisum.begin(), Cphisum.end(), 0.0);
	
		for(jdir=0;jdir<NJDIR;jdir++){
			sprintf(oldfilename,"%s_%d/%s/%s/%s/%s",dirprefix.c_str(),jdir,qualifier.c_str(),acceptance.c_str(),
			dirname[itype].c_str(),olddatafilename.c_str());
			//printf("oldfilename=%s\n",oldfilename);
			fptr_read=fopen(oldfilename,"r");
			fgets(dummy,200,fptr_read);
			for(iy=0;iy<NYBINS;iy++){
				for(iphi=0;iphi<NPHIBINS;iphi++){
					fscanf(fptr_read,"%lf",&bf);
					//printf("%6.3f ",bf);
					Cysum[iy]+=bf*normy;
					Cphisum[iphi]+=bf*normphi;
				}
				//printf("\n");
			}
			fclose(fptr_read);
		}
		//for(iy=0;iy<NYBINS;iy++)
		//	printf("Cysum[%d]=%g\n",iy,Cysum[iy]);
	}
	// Read/sum BF numerators
	olddatafilename="bf_Byphi.dat";
	for(itype=0;itype<NTYPES;itype++){
		std::fill(Bysum.begin(), Bysum.end(), 0.0);
		std::fill(Bphisum.begin(), Bphisum.end(), 0.0);
		
		for(jdir=0;jdir<NJDIR;jdir++){
			sprintf(oldfilename,"%s_%d/%s/%s/%s/%s",dirprefix.c_str(),jdir,qualifier.c_str(),acceptance.c_str(),
			dirname[itype].c_str(),olddatafilename.c_str());
			//printf("oldfilename=%s\n",oldfilename);
			fptr_read=fopen(oldfilename,"r");
			fgets(dummy,200,fptr_read);
			for(iy=0;iy<NYBINS;iy++){
				for(iphi=0;iphi<NPHIBINS;iphi++){
					fscanf(fptr_read,"%lf",&bf);
					//printf("%6.3f ",bf);
					if(Cysum[iy]>100){
						Bysum[iy]+=bf*normy/Cysum[iy];
						Bphisum[iphi]+=bf*normphi/Cysum[iy];
					}
				}
				//printf("\n");
			}
			fclose(fptr_read);
		}
		
		
		
		sprintf(newfilename,"%s_sum/%s/%s/%s/bfsum_y.dat",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str(),dirname[itype].c_str());
		fptr_write=fopen(newfilename,"w");
		printf("newfilename=%s\n",newfilename);
		for(iy=0;iy<NYBINS;iy++){
			fprintf(fptr_write,"%4.2f %g\n",(iy+0.5)*DY,Bysum[iy]);
			printf("%4.2f %g\n",(iy+0.5)*DY,Bysum[iy]);
		}
		fclose(fptr_write);
		
		sprintf(newfilename,"%s_sum/%s/%s/%s/bfsum_phi.dat",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str(),dirname[itype].c_str());
		fptr_write=fopen(newfilename,"w");
		printf("newfilename=%s\n",newfilename);
		for(iphi=0;iphi<NPHIBINS;iphi++){
			fprintf(fptr_write,"%4.2f %g\n",-180.0+(iphi+0.5)*DPHI,Bphisum[iphi]);
			printf("%4.2f %g\n",-180.0+(iphi+0.5)*DPHI,Bphisum[iphi]);
		}
		fclose(fptr_write);
		
	}
	return 0;
}


