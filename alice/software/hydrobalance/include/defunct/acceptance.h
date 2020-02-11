#ifndef __ACCEPTANCE_H__
#define __ACCEPTANCE_H__

#include "hbdefs.h"
#include "part.h"
#include <iostream>
#include <fstream>

class CAcceptance{
public:
	CAcceptance();
	CAcceptance(CparameterMap *parmapin);
	double ETAMIN,ETAMAX,PTMIN,PTMAX;
	int CENTRALITY;
	CparameterMap *parmap;
	virtual void CalcAcceptance(bool &accept,double &efficiency,CPart *part);
	virtual void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
};

class CAcceptance_CHEAP : public CAcceptance{
public:
	CAcceptance_CHEAP(CparameterMap *parmapin);
	double ptmin,ptmax;
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part);
	void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
};

class CAcceptance_ALICE  : public CAcceptance{
public:
	CAcceptance_ALICE(CparameterMap *parmapin);
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part);
	void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
};

class CAcceptance_STAR : public CAcceptance{
public:
	CAcceptance_STAR(CparameterMap *parmapin);
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part);
	void CalcAcceptanceNoID(bool &accept,double &efficiency,CPart *part);
	void star_acc_eff(int pid,double pt,double eta,double phi,int cen,bool &accept,double &eff);
	
	// For efficiencies of non-Identified parts
  //double ManuelEff(double eta, double pt);
	double ScottEff(double eta, double pt);
	
private:
	double ManuelData[20][30];
};


#endif
