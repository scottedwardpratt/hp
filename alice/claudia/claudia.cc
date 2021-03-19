#include "b3d.h"
#include "resonances.h"
#include "constants.h"

using namespace std;
void GetLatticeEoS(vector<double> &Parray,vector<double> &sarray,vector<double> &epsilonarray);

int main(){
	vector<double> Parray,sarray,epsilonarray;
	vector<double> density,maxweight;
	vector< vector< vector<double> > > chiharray,chiarray;
	Eigen::Matrix3d chih;
	char dummy[250];

	Parray.resize(401); sarray.resize(401); epsilonarray.resize(401);
	chiarray.resize(401); chiharray.resize(401);
	double T,chi,error,Ph,epsilonh,sh,nh,TH=155,TQGP=165.0,delT,w;
	delT=TQGP-TH;
	int iT,a,b,iupper,ilower;
	for(iT=0;iT<401;iT+=5){
		chiharray[iT].resize(3);
		chiarray[iT].resize(3);
		for(a=0;a<3;a++){
			chiharray[iT][a].resize(3);
			chiarray[iT][a].resize(3);
			for(b=0;b<3;b++){
				chiarray[iT][a][b]=-1000000.0;
			}
		}
	}

	CparameterMap parmap;
	parmap.set("RESONANCE_INFO_FILE",string("progdata/resinfo/resonances_pdg_weak.dat"));
	parmap.set("RESONANCE_DECAYS_FILE",string("progdata/resinfo/decays_pdg_weak.dat"));
	CResList *reslist=new CResList(&parmap);
	string command="rm -f chi_BQS.dat; rm -f chi_uds.dat";
	system(command.c_str());
	//  This prints chi/s 3x3 matrices for these temperatures to chi_uds.dat

	FILE *fptrh=fopen("chihadron.dat","w");
	for(T=20;T<TQGP+20.1;T+=5.0){
		reslist->CalcEoSandChi(T,Ph,epsilonh,nh,density,maxweight,chih);
		sh=(Ph+epsilonh)/T;
		iT=lrint(T);
		fprintf(fptrh,"%5.1f ",T);
		for(a=0;a<3;a++){
			for(b=0;b<3;b++){
				chiharray[iT][a][b]=chih(a,b)/sh;
				fprintf(fptrh,"%8.6f ",chiharray[iT][a][b]);
			}
		}
		fprintf(fptrh,"\n");
	}
	fclose(fptrh);
	
	//command="mv chi_uds.dat chihadron.dat";
	//system(command.c_str());
	
	GetLatticeEoS(Parray,sarray,epsilonarray);
	
	FILE *fptr=fopen("latticedata/chi_uu.dat","r");
	printf("----- chi_uu/s --------\n");
	do{
		fscanf(fptr,"%lf %lf %lf",&T,&chi,&error);
		fgets(dummy,250,fptr);
		iT=lrint(T);
		chiarray[iT][0][0]=chiarray[iT][1][1]=chi*T*T*T/sarray[iT];
		printf("%5.1f %g %g\n",
		T,T*T*T*chi/sarray[iT],chiarray[iT][0][0]);
	}while(!feof(fptr) && T<399.9);
	fclose(fptr);
	
	fptr=fopen("latticedata/c2ud.cont","r");
	printf("----- chi_ud/s --------\n");
	do{
		fscanf(fptr,"%lf %lf %lf",&T,&chi,&error);
		fgets(dummy,250,fptr);
		iT=lrint(T);
		chiarray[iT][0][1]=chiarray[iT][1][0]=chi*T*T*T/sarray[iT];
		printf("%5.1f %g %g\n",T,T*T*T*chi/sarray[iT],chiarray[iT][0][1]);
	}while(!feof(fptr) && T<399.9);
	fclose(fptr);
	
	fptr=fopen("latticedata/c2us.cont","r");
	printf("----- chi_us/s --------\n");
	do{
		fscanf(fptr,"%lf %lf %lf",&T,&chi,&error);
		fgets(dummy,250,fptr);
		iT=lrint(T);
		chiarray[iT][0][2]=chiarray[iT][1][2]=chiarray[iT][2][0]
			=chiarray[iT][2][1]=chi*T*T*T/sarray[iT];
		printf("%5.1f %g %g\n",T,T*T*T*chi,chiarray[iT][0][2]);
	}while(!feof(fptr) && T<399.9);
	fclose(fptr);
	
	fptr=fopen("latticedata/c2S.cont","r");
	printf("----- chi_ss/s --------\n");
	do{
		fscanf(fptr,"%lf %lf %lf",&T,&chi,&error);
		fgets(dummy,250,fptr);
		iT=lrint(T);
		chiarray[iT][2][2]=chi*T*T*T/sarray[iT];
		printf("%5.1f %g %g\n",T,T*T*T*chi/sarray[iT],chiarray[iT][2][2]);
	}while(!feof(fptr) && T<399.9);
	fclose(fptr);
	
	// NOW DO OUTPUT

	FILE *output=fopen("chi_overs.dat","w");
	fprintf(output,"#  T  chi_uu/s chi_du/s chi_su/s chi_ss/s\n");
	for(T=20.0;T<400.5;T+=5.0){
		fprintf(output,"%5.1f",T);
		iT=lrint(T);
		for(a=0;a<3;a++){
			for(b=0;b<=a;b++){
				if((a==0&&b==0) || (a==1&&b==0) || (a==2&&b==0) || (a==2&&b==2)){
					if(T<TH){
						chi=chiharray[iT][a][b];
					}
					else{
						//chi=chiarray[iT][a][b];	
						if(chiarray[iT][a][b]<-10000){
							ilower=iT;
							do{
								ilower=ilower-5;
							}while(chiarray[ilower][a][b]<-10000);
							iupper=iT;
							do{
								iupper=iupper+5;
							}while(chiarray[iupper][a][b]<-10000);
							w=double(iupper-iT)/double(iupper-ilower);
							chi=w*chiarray[ilower][a][b]
								+(1.0-w)*chiarray[iupper][a][b];
						}
						else{
							chi=chiarray[iT][a][b];
						}
						
						if(T<TQGP)
							chi=chiharray[iT][a][b]*(TQGP-T)/delT+chi*(T-TH)/delT;
					}
					fprintf(output," %8.5f",chi);
				}
			}
		}
		fprintf(output,"\n");
	}
	fclose(output);
	
	return 0;
}

void GetLatticeEoS(vector<double> &Parray,vector<double> &sarray,vector<double> &epsilonarray){
	const double h0=0.1396,h1=-0.1800,h2=0.0350,f0=2.76,f1=6.79,f2=-5.29,g1=-0.47,g2=1.04;
	double dT=1.0,poverT4,T,t,I,cs2,hbarc3=pow(HBARC,-3);
	FILE *fptr;
	int iT;
	poverT4=0.0;
	Parray[0]=epsilonarray[0]=sarray[0]=0.0;
	printf("  T   s/T^3\n");
	fptr=fopen("eos.dat","w");
	fprintf(fptr,"# T P(MeV/fm^3) epsilon  s(fm^-3) c_s^2\n");
	for(iT=0;iT<400;iT++){
		T=(iT+0.5)*dT;
		t=T/200.0;
		I=pow(T,4.0)*exp(-(h1/t)-h2/(t*t))*( h0+f0*(tanh(f1*t+f2)+1.0)/(1.0+g1*t+g2*t*t) );
		poverT4+=(dT/T)*I/pow(T,4.0);
		
		T=(iT+1.0)*dT;
		Parray[iT+1]=pow(T,4)*poverT4;
		I=pow(T,4.0)*exp(-(h1/t)-h2/(t*t))*( h0+f0*(tanh(f1*t+f2)+1.0)/(1.0+g1*t+g2*t*t) );
		epsilonarray[iT+1]=I+3.0*Parray[iT+1];
		sarray[iT+1]=(Parray[iT+1]+epsilonarray[iT+1])/T;
		printf("%5.1f %g\n",T,sarray[iT]/pow(T,3));
		if(iT>0){
			cs2=(Parray[iT+1]-Parray[iT-1])/(epsilonarray[iT+1]-epsilonarray[iT-1]);
			fprintf(fptr,"%5.1f %14.6e %14.6e %14.6e %14.6e\n",iT*dT,Parray[iT]*hbarc3,epsilonarray[iT]*hbarc3,sarray[iT]*hbarc3,cs2);
		}
	}
	fclose(fptr);
	printf("-------------------------\n");
}