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
	//string dirprefix="default";
	char qs[200];
	const int NTYPES=10,NJDIR=720,J0DIR=0;
	bool CalcGammaP=false;
	FILE *fptr_read;
	FILE *fptr_write;
	char oldfilename[100],newfilename[100];
	vector<double> ysum,xsum;
	double x,y,xbar,xbarnorm,x2bar,fudge,bfnorm,error,dx;
	int jdir,itype,ifn,ix;
	printf("Enter qualifier: ");
	scanf("%s",qs);
	qualifier=qs;
	string olddatafilename[3]={"bf_phi.dat","bf_eta.dat","bf_y.dat"};
	string newdatafilename[3]={"bf_phi.dat","bf_eta.dat","bf_y.dat"};
	string dirname[NTYPES]={"KK","Kp","piK","pip","pipi","pp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90"};
	char command[120];
	for(itype=0;itype<NTYPES;itype++){
		sprintf(command,"mkdir -p %s_sum/%s/%s/%s",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str(),dirname[itype].c_str());
		system(command);
	}
	
	// Write denom sum
	double denom_pi=0,denom_K=0,denom_p=0,denom_allcharges=0;
	double denom_allcharges_phi0=0,denom_allcharges_phi45=0,denom_allcharges_phi90=0;
	double gammap=0.0,gammap2=0.0,gammap1;
	char dummy[100];
	string gpstring,dumbo;
	double nplus,nminus;
	for(jdir=J0DIR;jdir<NJDIR;jdir++){
		sprintf(oldfilename,"%s_%d/%s/%s/denom.dat",dirprefix.c_str(),jdir,qualifier.c_str(),acceptance.c_str());
		//printf("oldfilename=%s\n",oldfilename);
		fptr_read=fopen(oldfilename,"r");
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_pi+=nplus+nminus;
		//printf("nplus=%lf nminus=%lf, denom_pi=%lf\n",nplus,nminus,denom_pi);
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
		if(CalcGammaP){
			sprintf(oldfilename,"%s_%d/%s/%s/gammap.dat",dirprefix.c_str(),jdir,qualifier.c_str(),acceptance.c_str());
			fptr_read=fopen(oldfilename,"r");
			fscanf(fptr_read,"%lf",&gammap1);
			//fscanf(fptr_read,"%s",dummy);
			//dumbo=dummy;
			//gpstring=dumbo.substr(7,dumbo.length());
			//gammap1=atof(gpstring.c_str());
			//printf("%4d: gammap for this event=%g\n",jdir,gammap1);
			fgets(dummy,100,fptr_read);
			gammap+=gammap1;
			gammap2+=gammap1*gammap1;
			fclose(fptr_read);
		}
	}
	//
	sprintf(newfilename,"%s_sum/%s/%s/denom.dat",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str());
	fptr_write=fopen(newfilename,"w");
	fprintf(fptr_write,"denom_pi: %g\ndenom_K: %g\ndenom_p: %g\ndenom_allcharges: %g\ndenom_allcharges_phi0: %g\ndenom_allcharges_phi45: %g\ndenom_allcharges_phi90 %g\n",
	denom_pi,denom_K,denom_p,denom_allcharges,denom_allcharges_phi0,denom_allcharges_phi45,denom_allcharges_phi90);
	fclose(fptr_write);
	
	if(CalcGammaP){
		sprintf(newfilename,"%s_sum/%s/%s/gammap.dat",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str());
		fptr_write=fopen(newfilename,"w");
		gammap=gammap/double(NJDIR-J0DIR);
		gammap2=gammap2/double(NJDIR-J0DIR);
		gammap2=sqrt((gammap2-gammap*gammap)/double(NJDIR-J0DIR));
		fprintf(fptr_write,"%g  %g\n",gammap,gammap2);
		printf("%g\n",gammap/double(NJDIR-J0DIR));
		fclose(fptr_write);
	}
	
	double denom;
	for(itype=0;itype<NTYPES;itype++){
		for(ifn=0;ifn<3;ifn++){
			dx=0.0;
			if(ifn==0)
				dx=360.0/28.0;
			if(ifn==1)
				dx=0.1;
			if(ifn==2)
				dx=0.1;
			sprintf(newfilename,"%s_sum/%s/%s/%s/%s",dirprefix.c_str(),qualifier.c_str(),acceptance.c_str(),dirname[itype].c_str(),
			newdatafilename[ifn].c_str());
			printf("newfilename=%s\n",newfilename);
			fptr_write=fopen(newfilename,"w");
			xsum.clear();
			ysum.clear();
			xbar=xbarnorm=x2bar=0.0;
			for(jdir=J0DIR;jdir<NJDIR;jdir++){
				denom=double(NJDIR-J0DIR);
				sprintf(oldfilename,"%s_%d/%s/%s/%s/%s",dirprefix.c_str(),jdir,qualifier.c_str(),acceptance.c_str(),
				dirname[itype].c_str(),olddatafilename[ifn].c_str());
				//printf("oldfilename=%s\n",oldfilename);
				fptr_read=fopen(oldfilename,"r");
				ysum.clear();
				ix=0;
				fscanf(fptr_read,"%lf %lf",&x,&y);
				fgets(dummy,100,fptr_read);
				fudge=1.0;
				//if(ifn>0 && jdir<200)
					//fudge=0.25;
				y=y*fudge;
				do{
					if(jdir==J0DIR){
						xsum.push_back(x);
						ysum.push_back(y);
					}
					else{
						ysum[ix]+=y;
					}
					xbar+=y*x;
					x2bar+=x*x*y;
					xbarnorm+=y;
					ix+=1;
					fscanf(fptr_read,"%lf %lf",&x,&y);
					fgets(dummy,100,fptr_read);
					fudge=1.0;
					//if(ifn>0 && jdir<200)
					//	fudge=0.25;
					y=y*fudge;
				}while(!feof(fptr_read));
				fclose(fptr_read);
			}
			int nx=xsum.size();
			bfnorm=0.0;
			for(ix=0;ix<nx;ix++){
				double ybar=ysum[ix];
				if(ifn==0 && (itype!=8))
					ybar=0.5*(ysum[ix]+ysum[nx-ix-1]);
				ybar=ybar/double(NJDIR-J0DIR);
				fprintf(fptr_write,"%5.2f %10.3e\n",xsum[ix],ybar);
				bfnorm+=ybar*dx;
				if(itype==4 && ifn==0){
					printf("%2d %g\n",ix,ybar);
				}
			}
			fclose(fptr_write);
			//printf("ifn=%d, <x>=%g,  <x^2>=%g\n",ifn,xbar/xbarnorm,x2bar/xbarnorm);
			//if(ifn==2 && itype==6)
			printf("bfnorm=%g\n------------------------\n",bfnorm);
		}
	}
	
	
	
	return 0;
}


