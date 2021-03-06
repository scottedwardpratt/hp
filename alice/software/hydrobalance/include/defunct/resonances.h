#ifndef __RESONANCES_H__
#define __RESONANCES_H__

#include <map>
#include <cmath>
#include <array>
#include <gsl/gsl_sf.h>
#include <armadillo>
#include "parametermap.h"
#include "randy.h"
#include "hbdefs.h"

using namespace std;

class CBranchInfo{
public:
	vector<CResInfo *> resinfoptr; //pointers for resinfo
	double branching; 
	CBranchInfo();
};

class CMerge{
public:
	CMerge(CResInfo *resinfo,double branching, int L);
	CResInfo *resinfo;
	int L;
	double branching;
	CMerge *next;
};

class CResInfo{
public:
	int ires;
	double mass;
	double spin;
	double width; 
	double minmass;
	double netchi;
	string name;
	int code;
	int charge;
	int strange;
	int baryon;
	int q[3];
	int up,down;
	int G_Parity;
	bool decay; //false if stable, true if can decay. check if true
	CBranchList branchlist; 
	CBranchInfo	*bptr_minmass;
	void Print();
	void DecayGetResInfoPtr(int &nbodies,array<CResInfo *,5> &daughterresinfo);
	void DecayGetResInfoPtr_minmass(int &nbodies,array<CResInfo *,5> &daughterresinfo);
	bool CheckForDaughters(int code);
	bool CheckForNeutral();
	double GenerateMass();
	double GenerateThermalMass(double maxweight, double T);
	double ChiOmega(double vOmega,double T); //Calculates \int d^3p/(hbar^3 E) p.dOmega exp(-u.p/T)
	void MakePart(CHyperElement *hyper,FourVector &p,double maxweight);
	double ChiInt(double T,double vmax); // Integral used by ChiOmega
	double ChiTilde(double T,double vmax); // Integral used by ChiOmega
	CResInfo();
	static CRandy *randy;
	static CResList *reslist;
	static double **ChiA; // stored array used by ChiOmegaInt
};

class CResList{
public:
	CResList();
	~CResList();
	CResList(CparameterMap* parmap_in);
	CResInfoMap resmap;
	CResInfo *GetResInfoPtr(int ID);
	void ReadResInfo();
	void CalcEoSandChi(double T,double &P,double &epsilon,double &nh,vector<double> &density,vector<double> &maxweight,arma::mat &chi);
	void freegascalc_onespecies(double m,double T,double &e,double &p,double &dens,double &sigma2,double &dedt);
	void freegascalc_onespecies_finitewidth(double m, double m1, double m2, double T,double width,double minmass,double &epsilon,double &P,double &dens,double &sigma2,double &dedt, double &maxweight);
	double GetLambda(double T,double P,double epsilon);
	CparameterMap *parmap;
	CMerge ***MergeArray;
	double **SigmaMaxArray;
	double RESWIDTH_ALPHA;
	bool RESONANCE_DECAYS;
	bool USEPOLEMASS;
	//void freegascalc_onespecies_offshell(CResInfo *resinfo,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt);
	double triangle(double m0,double m1,double m2);
	static CBalance *cb;
};

#endif
