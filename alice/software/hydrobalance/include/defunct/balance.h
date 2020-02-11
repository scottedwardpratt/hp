#ifndef __BALANCE_H__
#define __BALANCE_H__

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <random>
#include <list>
#include <map>

#include "commondefs.h"
#include "part.h"
#include "randy.h"
#include "misc.h"
#include "acceptance.h"

using namespace std;

class CBalance{
public:
	CBalance(CHydroBalance *hbset);
	CBalance(string parfilename,string qualifierset);
	void Init();
	CparameterMap parmap;
	CResList *reslist;
	CAcceptance *acceptance;
	string acceptance_description;
	string bf_results_dirname;
	string qualifier;
	double NSAMPLE_HYDRO2UDS,NSAMPLE_UDS2BAL;
	CRandy *randy;
	double ransum,ranthresh,Nerror;  //Nerror = number of particles not emitted due to giving up in MC
	void ReadResonanceList();
	void ReadCharges();
	void ReadHyper();
	void CalcChiTotFromQ();
	void CalcChiTotFromH();
	
	// Generates partmap from uds map to build BF numerator
	void GenDecayProcessHadronsFromCharges();
	void GenHadronsFromCharge(int balanceID,CCharge *charge);
	// Generates ncpartmap from dOmega list to build BF denominator
	void GenHadronsFromHyperElements();
	void GenHadronsFromHyperElement(CHyperElement *hyper);
	void CalcChiInv(double vOmega,Eigen::Matrix3d &chiinv);
	int Nchi;
	
	void DecayParts(mapip &partmap);
	bool GetDecayProducts(CPart *part0,int &ni,array<CPart,10> &parti);
	void BoostPart(CPart *part,FourVector &u,double &ytarget,bool jdecay);
	void DecayPart(CPart *mother,int &nbodies,array<CPart,7> &daughter);
	void CalcDCA(bool decay,CPart *part,double *dca);
	void IncrementNumer(CPart *parta,CPart *partb,double w);
	void IncrementDenom(CPart *part,double w);
	void IncrementGammaP(CPart *parta,CPart *partb,double w);
	void DeletePartMapParts();
	void DeleteEMapCharges();
	void DeleteNCPartMapParts();
	void DeleteHyperElements();
	
	void ProcessPartMap();
	void ProcessNCPartMap();
	
	CPart* GetNewPart(mapip *partmap,int balanceid);
	list<CPart*> deadpartlist;
	mapip partmap; // map of correlated particles (generated from correlated charges)
	mapip ncpartmap; // map of uncorrelated particles (generated from dOmega)
	
	mapic emap; // emitted charges
	list<CHyperElement *> hyperlist;   // list of hyper-elements
	
	CBFNumer *numer_pipi,*numer_piK,*numer_pip,*numer_KK;
	CBFNumer *numer_Kp,*numer_pp,*numer_charged;
	CBFDenom *denom_pi,*denom_K,*denom_p,*denom_charged;
	CBFNumer *bf_pipi,*bf_piK,*bf_pip,*bf_KK,*bf_Kp,*bf_pp,*bf_charged;
	CBFNumer *numer_charged_phi0,*numer_charged_phi45,*numer_charged_phi90;
	CBFNumer *bf_charged_phi0,*bf_charged_phi45,*bf_charged_phi90;
	CBFDenom *denom_charged_phi0,*denom_charged_phi45,*denom_charged_phi90;
	double gammap,normtest,v2;

	double Tf,epsilonf,Pf,lambdaf,nf; // freezeout props
	double TotalVolume;
	vector<double> densityf,maxweightf;
	Eigen::Matrix3d chif,chiinvf;
	
	void PrintBFNumer();
	void PrintBFDenom();
	void CreateBFArrays();
	void ConstructBFs();
	void ConstructBF(CBFNumer *numer,CBFDenom *denom,CBFNumer *bf,double doublecount);
	void WriteBFs();
	void CheckPartMap(mapip *partmap);
	void CalcChiTot();
	void WriteChiTrans();
	
	double YMIN,YMAX;  // Dependent particle boosted randomly to within YMIN/YMAX
private:
	int badmothercount;
	
};

#endif
