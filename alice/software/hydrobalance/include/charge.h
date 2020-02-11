#ifndef __CHARGE_H__
#define __CHARGE_H__

#include "commondefs.h"
#include "hyper.h"

using namespace std;

class CCharge{
public:
	CCharge();
	~CCharge(){};
	bool active;
	int q[3];
	double x,y,eta,tau,weight,vx,vy,rapidity;
	CHyperElement hyper;
	void Propagate(double newtau);
	void SetV(double ux,double uy);
	void Print();
	static CHydroBalance *hb;
	static CB3D *b3d;
};

#endif