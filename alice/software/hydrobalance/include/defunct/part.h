#ifndef __PART_H__
#define __PART_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <sys/stat.h>
#include "b3d.h"
#include "hbdefs.h"
#include "hydro2uds.h"
#include "balance.h"

using namespace std;

class CHBPart{
public:
	CPart();
	~CPart();
	double tau0;
	double y,eta,msquared,bweight;
	FourVector p,r;
	int balanceID;  // key for balance maps
	CResInfo *resinfo;
	mapip *partmap;
	int badmother;
	
	void Init(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity,int weight,int balanceid);
	void Propagate(double tau);
	void FindDecay();
	void CalcDCA(double *dca);
	void Setp0();
	void SetMass();
	void Copy(CPart *part);
	void CopyPositionInfo(CPart *part);
	void CopyMomentumInfo(CPart *part);
	double GetMass();
	void Kill();
	void Print();
	void BoostRap(double dely);
	void Boost(FourVector &u);
	void BoostP(FourVector &u);
	void BoostR(FourVector &u);
	//~CPart();

	double GetEta(double tau);
	double GetPseudoRapidity();
	double GetRapidity();
	double GetMT();
	void SetY();
	void SetEta(double neweta);
	static CBalance *cb;
};

class CB3DBinaryBalancePartInfo{
public:
	int ID,balanceID;
	double px,py,rapidity;
	void Print(){
		printf("ID=%d, balanceID=%d\n",ID,balanceID);
		printf("px=%g,py=%g, rapidity=%g\n",px,py,rapidity);
	}
};

#endif
