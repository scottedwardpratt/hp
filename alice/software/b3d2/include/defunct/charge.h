#ifndef __CHARGE_H__
#define __CHARGE_H__

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <random>
#include "b3ddefs.h"
#include "hyper.h"
#include "sampler.h"
#include "b3d.h"
#include "misc.h"

using namespace std;

class CHydroBalance;

class CCharge{
public:
	CCharge();
	~CCharge(){};
	int q[3];
	double x,y,eta,tau,weight,vx,vy,rapidity;
	CHyperElement *hyper;
	void Print();
	static CB3D *b3d;
};

#endif