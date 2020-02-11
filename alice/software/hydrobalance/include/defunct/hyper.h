#ifndef __HYPER_H__
#define __HYPER_H__

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include "commondefs.h"

// Info for HyperSurface
class CHyperElement{
public:
	~CHyperElement(){};
	double tau,x,y,dOmegaX,dOmegaY,dOmega0,ux,uy,T;
	double pitildexx,pitildexy,pitildeyy;
	double udotdOmega,dOmegaMax,vOmega,dOmegaVec;
	void Copy(CHyperElement *oldhyper);
	void Print();
	void CalcDOmegaMax();
	// Note dOmegaX and dOmegaY = d\Omega_\mu
	// (subscript mu, so p.Omega = p0*dOmega0+px*dOmegaX+py*dOmegaY)
	static CHydroBalance *hb;
};

#endif