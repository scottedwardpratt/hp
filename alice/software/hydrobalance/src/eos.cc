#include "eos.h"
#include "resonances.h"
#include "parametermap.h"
#include "constants.h"

using namespace std;

vector<double> CEoS::epsilon_PST,CEoS::P_PST,CEoS::s_PST,CEoS::T_PST;
vector<double> CEoS::epsilon_claudia,CEoS::P_claudia,CEoS::s_claudia,CEoS::T_claudia;
vector<double> CEoS::twopiTD,CEoS::Tdiff;
vector<double> CEoS::chill_claudia,CEoS::chiud_claudia,CEoS::chils_claudia,CEoS::chiss_claudia;
vector<double> CEoS::chill_HSC,CEoS::chiud_HSC,CEoS::chils_HSC,CEoS::chiss_HSC;
vector<double> CEoS::dDdT;
mapdi CEoS::etmap;

CResList *CEoS::reslist=NULL;

CEoS::CEoS(){	
}

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;
	reslist=new CResList(parmap);
	ReadDiffusionData();
	ReadChiData_Claudia();
	FillOutdDdT();
};

void CEoS::ReadDiffusionData(){
	string dirname=parmap->getS("LATTICEDATA_DIRNAME","../latticedata");
	string filename=dirname+"/diffusion.dat";
	char dummy[100];
	char voldummy;
	int ntaudummy;
	int ndata=7;
	double errsysdummy,errsysstatdummy,t,td;
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,100,fptr);
	//printf("%s\n",dummy);
	fscanf(fptr,"%s",&voldummy);
	while(!feof(fptr)){
		fscanf(fptr,"%lf %d %lf %lf %lf",&t,&ntaudummy,&td,&errsysdummy,&errsysstatdummy);
		Tdiff.push_back(t*0.001);
		twopiTD.push_back(td);
		fscanf(fptr,"%s",&voldummy);
	}	
	fclose(fptr);
}

void CEoS::FillOutdDdT(){
	int ndata=epsilon_PST.size();
	dDdT.resize(ndata);
	Eigen::Matrix3d chi0,chi1,chi2,chiinv0,chiinv1,chiinv2,Dmat;
	chi0.setZero(); chi1.setZero(); chi2.setZero(); chiinv0.setZero(); chiinv1.setZero(); chiinv2.setZero();
	Dmat.setZero();
	int ie;
	double T0,T1,T2;
	for(ie=1;ie<ndata-1;ie++){
		T0=T_PST[ie-1];
		T1=T_PST[ie];
		T2=T_PST[ie+1];
		T=T0;
		GetChiOverS_Claudia();  // check whether this gives Chi or ChiOverS
		chi0(0,0)=chi0(1,1)=chill;
		chi0(0,1)=chi0(1,0)=chiud;
		chi0(2,2)=chiss;
		chi0(0,2)=chi0(1,2)=chi0(2,0)=chi0(2,1)=chils;
		chiinv0=chi0.inverse();
		T=T1;
		GetChiOverS_Claudia();
		chi1(0,0)=chi1(1,1)=chill;
		chi1(0,1)=chi1(1,0)=chiud;
		chi1(2,2)=chiss;
		chi1(0,2)=chi1(1,2)=chi1(2,0)=chi1(2,1)=chils;
		chiinv1=chi1.inverse();
		T=T2;
		GetChiOverS_Claudia();
		chi2(0,0)=chi2(1,1)=chill;
		chi2(0,1)=chi2(1,0)=chiud;
		chi2(2,2)=chiss;
		chi2(0,2)=chi2(1,2)=chi2(2,0)=chi2(2,1)=chils;
		chiinv2=chi2.inverse();
		Dmat=(chiinv2-chiinv0)/(T2-T0);
		Dmat=GetD(T1)*chi1*Dmat;
		/*
		printf("---- T=(%g,%g,%g) ----\n",T0,T1,T2);
		chi0.print("chi0");
		chiinv0.print("chiinv0");
		chi1.print("chi1");
		chiinv1.print("chiinv1");
		chi2.print("chi2:");
		chiinv2.print("chiinv2");
		Dmat.print("Dmat\n");
		*/
	}
}

double CEoS::GetD(double T){
	int i0;
	double result,delT;
	if(T>=Tdiff[0] && T<=Tdiff[Tdiff.size()-1]){
		i0=0;
		while(T>Tdiff[i0+1]){
			i0+=1;
		}
		delT=Tdiff[i0+1]-Tdiff[i0];
		result=((T-Tdiff[i0])*twopiTD[i0+1]/delT)+(Tdiff[i0+1]-T)*twopiTD[i0]/delT;
	}
	else if(T<Tdiff[0]){
		delT=Tdiff[1]-Tdiff[0];
		result=twopiTD[0]-(Tdiff[0]-T)*(twopiTD[1]-twopiTD[0])/delT;
	}
	else{
		i0=Tdiff.size()-1;
		delT=Tdiff[i0]-Tdiff[i0-1];
 		result=twopiTD[i0]+(T-Tdiff[i0])*(twopiTD[i0]-twopiTD[i0-1])/delT;
	}
	result=result*HBARC*0.001/(2.0*PI*T);
	//result*=0.001;
	return result;
}

void CEoS::ReadChiData_HSC(){
	string dirname=parmap->getS("LATTICEDATA_DIRNAME","../latticedata");
	string filename;
	double error,chi0;
	int idata;
	const int ndata=7;
	double chiB[ndata],chiI[ndata],chiQ[ndata],chiSS[ndata],chiLL[ndata],T[ndata];
	FILE *fptr;
	filename=dirname+"/chi-B.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&T[idata],&chiB[idata],&error);
	}
	fclose(fptr);
	filename=dirname+"/chi-I.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&T[idata],&chiI[idata],&error);
	}
	fclose(fptr);
	fclose(fptr);
	filename=dirname+"/chi-LL.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&T[idata],&chiLL[idata],&error);
	}
	fclose(fptr);
	fclose(fptr);
	filename=dirname+"/chi-Q.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&T[idata],&chiQ[idata],&error);
	}
	fclose(fptr);
	fclose(fptr);
	filename=dirname+"/chi-SS.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&T[idata],&chiSS[idata],&error);
	}
	fclose(fptr);
	
	chill_HSC.resize(ndata);
	chils_HSC.resize(ndata);
	chiss_HSC.resize(ndata);
	chiud_HSC.resize(ndata);
	for(idata=0;idata<ndata;idata++){
		chi0=pow(T[idata]/HBARC,3)*(pow(PI,4)/90.0)*(4.0/(2.0*PI*PI))*2.0;
		chiLL[idata]*=chi0;
		chiB[idata]*=chi0/3.0;
		chiSS[idata]*=chi0;
		chiI[idata]*=0.5*chi0;
		chiQ[idata]*=2.0*chi0/3.0;
		chill_HSC[idata]=chiLL[idata];
		chiss_HSC[idata]=chiSS[idata];
		chiud_HSC[idata]=chiLL[idata]-2.0*chiI[idata];
		chils_HSC[idata]=2.25*(chiB[idata]-(2.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]-(2.0/9.0)*chiud_HSC[idata]);
		//chils[idata]=-4.5*(chiQ[idata]-(5.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]+(4.0/9.0)*chiud[idata]);
		printf("-------  T=%g --------\n",T[idata]);
		printf("%8.5f %8.5f %8.5f\n",chill_HSC[idata],chiud_HSC[idata],chils_HSC[idata]);
		printf("%8.5f %8.5f %8.5f\n",chiud_HSC[idata],chill_HSC[idata],chils_HSC[idata]);
		printf("%8.5f %8.5f %8.5f\n",chils_HSC[idata],chils_HSC[idata],chiss_HSC[idata]);
		/*
		chils_HSC[idata]=2.25*(chiB[idata]-(2.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]-(2.0/9.0)*chiud_HSC[idata]);
		printf("T=%g: chils_HSC[%d]=%g\n",T[idata],idata,chils_HSC[idata]);
		chils_HSC[idata]=-4.5*(chiQ[idata]-(5.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]+(4.0/9.0)*chiud[idata]);
		printf("T=%g: chils_HSC[%d]=%g\n",T[idata],idata,chils_HSC[idata]);
		printf("--------\n");
		*/
		
	}
}

void CEoS::ReadChiData_Claudia(){
	string dirname=parmap->getS("LATTICEDATA_DIRNAME","latticedata");
	string filename;
	double error,chi0,delT=5.0;
	int idata;
	const int ndata=81;
	char dummy[100];
	FILE *fptr;
	// You will read in chi/s, not chi
	filename=dirname+"/chi_claudia.dat";
	fptr=fopen(filename.c_str(),"r");
	T_claudia.resize(ndata);
	chill_claudia.resize(ndata);
	chils_claudia.resize(ndata);
	chiss_claudia.resize(ndata);
	chiud_claudia.resize(ndata);
	fgets(dummy,100,fptr);
	for(idata=4;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf %lf %lf",
		&T_claudia[idata],&chill_claudia[idata],&chiud_claudia[idata],&chils_claudia[idata],&chiss_claudia[idata]);
		//printf("T=%g, chiuu=%g, chiud=%g, chiss=%g\n",T_claudia[idata],chill_claudia[idata],chiud_claudia[idata],chiss_claudia[idata]);
	}
	fclose(fptr);
}

void CEoS::GetChiOverS_Claudia(){
	double delT=5.0,Tmax=400.0,w0,w1,sc,Tm=1000.0*T;
	const int ndata=81;
	double answer;
	int iT0;
	if(Tm<20){
		chill=chill_claudia[5]*Tm/20.0;
		chiud=chiud_claudia[5]*Tm/20.0;
		chils=chils_claudia[5]*Tm/20.0;
		chiss=chiss_claudia[5]*Tm/20.0;
	}
	else if(Tm>=Tmax){
		chill=chill_claudia[ndata-1];
		chiud=chiud_claudia[ndata-1];
		chils=chils_claudia[ndata-1];
		chiss=chiss_claudia[ndata-1];
	}
	else{
		iT0=lrint(floor(Tm/delT));
		if(iT0>=ndata-1){
			printf("iT0 out of range in CEoS::GetChiOverS_Claudia()\n");
			exit(1);
		}
		w1=(Tm-delT*iT0)/delT;
		w0=1.0-w1;
		chill=chill_claudia[iT0]*w0+chill_claudia[iT0+1]*w1;
		chiud=chiud_claudia[iT0]*w0+chiud_claudia[iT0+1]*w1;
		chils=chils_claudia[iT0]*w0+chils_claudia[iT0+1]*w1;
		chiss=chiss_claudia[iT0]*w0+chiss_claudia[iT0+1]*w1;
		//printf("iT0=%d, T=%g, chiuu=%g, chiud=%g, chiss=%g\n",iT0,T_claudia[iT0],chill_claudia[iT0],chiud_claudia[iT0],chiss_claudia[iT0]);
	}
	
}

void CEoS::CalcEoS_PST(){
	string filename
		=parmap->getS("EOS_PSTDATA_FILENAME","../eos/EOS_tables/EOS_PST.dat");
	FILE *fptr=fopen(filename.c_str(),"r");
	FILE *fptrh=fopen("hadroneos.dat","w");
	double epsilon;
	const int ndata=155500;
	epsilon_PST.resize(ndata);
	P_PST.resize(ndata);
	s_PST.resize(ndata);
	T_PST.resize(ndata);
	
	int ie=0,ne;
	fscanf(fptr,"%lf ",&epsilon);
	while(!feof(fptr)){
		epsilon_PST[ie]=epsilon;
		fscanf(fptr,"%lf %lf %lf",&P_PST[ie],&s_PST[ie],&T_PST[ie]);
		//printf("epsilon=%g, P=%g, s=%g,T=%g\n",epsilon_PST[ie],P_PST[ie],s_PST[ie],T_PST[ie]);
		ie+=1;
		fscanf(fptr,"%lf ",&epsilon);
	}
	
	double TH=0.155, TQGP=0.175;
	double T,Ph,eh,nh,sh,w;
	int a,nres=reslist->resmap.size();
	vector<double> density;
	vector<double> maxweight;
	Eigen::Matrix3d chi;
	density.resize(nres);
	maxweight.resize(nres);
	ne=ie;
	for(ie=0;ie<ne;ie++){
		if(T_PST[ie]<TQGP){
			T=T_PST[ie];
			reslist->CalcEoSandChi(T*1000.0,Ph,eh,nh,density,maxweight,chi);
			Ph=Ph/1000.0;
			eh=eh/1000.0;
			sh=(Ph+eh)/T;
			if(T<TH){
				P_PST[ie]=Ph;
				s_PST[ie]=sh;
				epsilon_PST[ie]=eh;
			}
			else{
				w=(TQGP-T)/(TQGP-TH);
				P_PST[ie]=w*Ph+(1.0-w)*P_PST[ie];
				s_PST[ie]=w*sh+(1.0-w)*s_PST[ie];
				epsilon_PST[ie]=w*eh+(1.0-w)*epsilon_PST[ie];
				//printf("Tpst=%g, Ppst=%g, spst=%g, epst=%g\n",
				//T_PST[ie],P_PST[ie],s_PST[ie],epsilon_PST[ie]);
				//printf("T=%g, Ph=%g, sh=%g, eh=%g\n",T,Ph,sh,eh);
			}
		}
		fprintf(fptrh,"%13.7e %13.7e %13.7e %13.7e\n",epsilon_PST[ie],P_PST[ie],s_PST[ie],T_PST[ie]);
	}
	//printf("ie=%d =? %d\n",ie,ndata);
	//printf("epsilon=%g, P=%g, s=%g,T=%g\n",epsilon_PST[ie-1],P_PST[ie-1],s_PST[ie-1],T_PST[ie-1]);
	fclose(fptr);
	fclose(fptrh);
	
}

void CEoS::ReadEoS_PST(){
	string filename
		=parmap->getS("EOS_PSTDATA_FILENAME","../eos/EOS_tables/EOS_PST.dat");
	FILE *fptr=fopen(filename.c_str(),"r");
	double epsilon;
	const int ndata=155500;
	epsilon_PST.resize(ndata);
	P_PST.resize(ndata);
	s_PST.resize(ndata);
	T_PST.resize(ndata);
	
	int ie=0,ne;
	fscanf(fptr,"%lf ",&epsilon);
	while(!feof(fptr)){
		epsilon_PST[ie]=epsilon;
		fscanf(fptr,"%lf %lf %lf",&P_PST[ie],&s_PST[ie],&T_PST[ie]);
		//printf("epsilon=%g, P=%g, s=%g,T=%g\n",epsilon_PST[ie],P_PST[ie],s_PST[ie],T_PST[ie]);
		ie+=1;
		fscanf(fptr,"%lf ",&epsilon);
	}
	
	int a,nres=reslist->resmap.size();
	vector<double> density;
	vector<double> maxweight;
	Eigen::Matrix3d chi;
	density.resize(nres);
	maxweight.resize(nres);
	fclose(fptr);
	
}

void CEoS::BuildMap(){
	int ie;
	//printf("check BuildMap, epsilon_PST.size()=%d\n",int(epsilon_PST.size()));
	for(ie=0;ie<epsilon_PST.size();ie++){
		etmap.insert(pairdi(T_PST[ie],ie));
		//printf("ie=%d, etmap.size()=%d\n",ie,int(etmap.size()));
	}
}

void CEoS::GetEoSFromEpsilon_PST(double epsilonset){
	double depsilon=0.02,w0;
	int ie0,ndata;
	epsilon=epsilonset;
	ie0=lrint(floor(-0.5+epsilon/depsilon));
	if(ie0<0)
		ie0=0;
	if(ie0>=epsilon_PST.size()-1){
		ndata=epsilon_PST.size();
		P=P_PST[ndata]+0.33*(epsilon_PST[ndata]-epsilon);
		T=T_PST[ndata]*pow(epsilon/epsilon_PST[ndata],0.25);
		s=(P+epsilon)/T;
	}
	else{
		w0=(epsilon-epsilon_PST[ie0+1])/depsilon;
		P=w0*P_PST[ie0]+(1.0-w0)*P_PST[ie0+1];
		T=w0*T_PST[ie0]+(1.0-w0)*T_PST[ie0+1];
		s=w0*s_PST[ie0]+(1.0-w0)*s_PST[ie0+1];
	}
	if(fabs(T-0.140)<0.01){
		printf("TESTING, s=%g\n",s);
		double T,Ph,eh,nh,sh;
		int a,nres=reslist->resmap.size();
		vector<double> density,maxweight;;
		Eigen::Matrix3d chi;
		density.resize(nres);
		maxweight.resize(nres);
		reslist->CalcEoSandChi(T*1000.0,Ph,eh,nh,density,maxweight,chi);
		Ph=Ph/1000.0;
		eh=eh/1000.0;
		sh=(Ph+eh)/T;
	}
}

void CEoS::Print(){
	printf("---- T=%g, P=%g, epsilon=%g, s=%g  ----\n",T,P,epsilon,s);
	printf("chill=%g, chiud=%g, chils=%g, chiss=%g\n",chill,chiud,chils,chiss);
	printf("chill/s=%g, chiud/s=%g, chilsovers=%g, chissovers=%g\n",
	chillovers,chiudovers,chilsovers,chissovers);
}

void CEoS::PrintChi(){
	Eigen::Matrix3d chimat;
	int a,b;
	chimat(0,0)=chimat(1,1)=chill;
	chimat(1,0)=chimat(0,1)=chiud;
	chimat(0,2)=chimat(1,2)=chimat(2,1)=chimat(2,0)=chils;
	chimat(2,2)=chiss;
	printf("--------- Chi(Tf)-------------\n");
	cout<< chimat << endl;
	chimat=chimat/chimat(0,0);
	cout << chimat << endl;
}

void CEoS::GetEoSFromT_PST(double Tset){
	mapdi::iterator iter;
	int ie;
	double w0,delT;
	T=Tset;	
	iter=etmap.lower_bound(T);
	ie=iter->second;
	//printf("first=%g, second=%d\n",iter->first,iter->second);
	//printf("ie=%d\n",ie);
	if(ie!=0)
		ie-=1;
	delT=T_PST[ie+1]-T_PST[ie];
	w0=(T_PST[ie+1]-T)/delT;
	epsilon=w0*epsilon_PST[ie]+(1.0-w0)*epsilon_PST[ie+1];
	P=w0*P_PST[ie]+(1.0-w0)*P_PST[ie+1];
	s=w0*s_PST[ie]+(1.0-w0)*s_PST[ie+1];
	/*
	if(fabs(T-0.140)<0.000001){
		printf("TESTING, T=%g, P=%g, s=%g,nres=%d\n",
		T,P,s,int(reslist->resmap.size()));
		double Ph=0,eh=0,nh,sh;
		int a,nres=reslist->resmap.size();
		vector<double> density;
		vector< vector<double> > chi;
		density.resize(nres);
		chi.resize(3);
		for(a=0;a<3;a++)
			chi[a].resize(3);
		reslist->CalcEoSandChi(T*1000.0,Ph,eh,nh,density,chi);
		sh=(Ph+eh)/(1000.0*T);
		printf("s=%g, sh=%g\n",s,sh);
	}*/	
}
