#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	string qualifier="alice_cent0_10";
	string acceptance="results_alice";
	string dirprefix="default";
	char qs[200];
	const int NTYPES=10;
	int NJDIR;
	bool CalcGammaP=false;
	FILE *fptr_read;
	FILE *fptr_write;
	char oldfilename[100],newfilename[100];
	vector<double> ysum[4],xsum;
	double x,y[4],bfnorm,error,dx;
	int jdir,itype,ifn,ix;
	
	printf("Enter qualifier: ");
	scanf("%s",qs);
	printf("Enter NJDIR: (the number of default_xxx dirs)");
	scanf("%d",&NJDIR);
	
	qualifier=qs;
	string olddatafilename[4]={"bf_phi.dat","bf_eta.dat","bf_y.dat","bf_qinv.dat"};
	string newdatafilename[4]={"bf_phi.dat","bf_eta.dat","bf_y.dat","bf_qinv.dat"};
	int Nqcut=4;
	string dirname[NTYPES]={"KK","Kp","piK","pip","pipi","pp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90"};
	char command[120];
	for(itype=0;itype<NTYPES;itype++){
		sprintf(command,"mkdir -p %s_sum/%s/%s/%s",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str(),dirname[itype].c_str());
		system(command);
	}
	
	// Read denom sum
	double denom_pi=0,denom_K=0,denom_p=0,denom_allcharges=0;
	double denom_allcharges_phi0=0,denom_allcharges_phi45=0,denom_allcharges_phi90=0;
	char dummy[100];
	string gpstring,dumbo;
	double nplus,nminus;
	printf("NJDIR=%d\n",NJDIR);
	for(jdir=0;jdir<NJDIR;jdir++){
		sprintf(oldfilename,"%s_%d/%s/%s/denom.dat",dirprefix.c_str(),jdir,qualifier.c_str(),acceptance.c_str());
		fptr_read=fopen(oldfilename,"r");
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_pi+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_K+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_p+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges_phi0+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges_phi45+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges_phi90+=nplus+nminus;
		fclose(fptr_read);		
		double normtest,mult;
	}
	sprintf(newfilename,"%s_sum/%s/%s/denom.dat",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str());
	fptr_write=fopen(newfilename,"w");
	fprintf(fptr_write,"denom_pi: %g\ndenom_K: %g\ndenom_p: %g\ndenom_allcharges: %g\ndenom_allcharges_phi0: %g\ndenom_allcharges_phi45: %g\ndenom_allcharges_phi90 %g\n",
	denom_pi,denom_K,denom_p,denom_allcharges,denom_allcharges_phi0,denom_allcharges_phi45,denom_allcharges_phi90);
	fclose(fptr_write);
	
	
	// Now average balance funtions
	for(itype=0;itype<NTYPES;itype++){
		if(itype<6)
			Nqcut=4;
		else
			Nqcut=3;
		for(ifn=0;ifn<Nqcut;ifn++){
			dx=0.0;
			if(ifn==0)
				dx=360.0/28.0;
			if(ifn==1)
				dx=0.1;
			if(ifn==2)
				dx=0.1;
			if(ifn==3)
				dx=10.0;
			sprintf(newfilename,"%s_sum/%s/%s/%s/%s",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str(),dirname[itype].c_str(),
			newdatafilename[ifn].c_str());
			fptr_write=fopen(newfilename,"w");
			xsum.clear();
			ysum[0].clear();
			ysum[1].clear();
			ysum[2].clear();
			ysum[3].clear();
			for(jdir=0;jdir<NJDIR;jdir++){
				sprintf(oldfilename,"%s_%d/%s/%s/%s/%s",dirprefix.c_str(),jdir,qualifier.c_str(),acceptance.c_str(),
				dirname[itype].c_str(),olddatafilename[ifn].c_str());
				fptr_read=fopen(oldfilename,"r");
				ix=0;
				fscanf(fptr_read,"%lf %lf",&x,&y[0]);
				if(ifn>2)
					fscanf(fptr_read,"%lf %lf %lf",&y[1],&y[2],&y[3]);
				fgets(dummy,100,fptr_read);

				while(!feof(fptr_read)){
					if(jdir==0){
						xsum.push_back(x);
						ysum[0].push_back(y[0]);
						if(ifn>2){
							ysum[1].push_back(y[1]);
							ysum[2].push_back(y[2]);
							ysum[3].push_back(y[3]);
						}
					}
					else{
						ysum[0][ix]+=y[0];
						if(ifn>2){
							ysum[1][ix]+=y[1];
							ysum[2][ix]+=y[2];
							ysum[3][ix]+=y[3];
						}
					}
					ix+=1;
					fscanf(fptr_read,"%lf %lf",&x,&y[0]);
					if(ifn>2)
						fscanf(fptr_read,"%lf %lf %lf",&y[1],&y[2],&y[3]);
					fgets(dummy,100,fptr_read);
				}
				fclose(fptr_read);
			}
			int nx=xsum.size();
			bfnorm=0.0;
			for(ix=0;ix<nx;ix++){
				y[0]=ysum[0][ix];
				if(ifn==0)
					y[0]=0.5*(ysum[0][ix]+ysum[0][nx-ix-1]);
				y[0]=y[0]/double(NJDIR);
				if(ifn>2){
					y[1]=ysum[1][ix]/double(NJDIR);
					y[2]=ysum[2][ix]/double(NJDIR);
					y[3]=ysum[3][ix]/double(NJDIR);
				}
				fprintf(fptr_write,"%5.2f %10.3e",xsum[ix],y[0]);
				if(ifn>2)
					fprintf(fptr_write," %10.3e %10.3e %10.3e",xsum[ix],y[0]);
				fprintf(fptr_write,"\n");
				bfnorm+=y[0]*dx;
				if(itype==4 && ifn==0){
					printf("%2d %g\n",ix,y[0]);
				}
			}
			fclose(fptr_write);
			printf("=============== bfnorm=%g ===============\n",bfnorm);
		}
	}
	
	
	
	return 0;
}


