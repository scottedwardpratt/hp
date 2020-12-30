#ifndef __COMMON_DEFS_H__
#define __COMMON_DEFS_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <complex>
#include <sys/stat.h>
#include <ctime>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <cmath>
#include <random>
#include <vector>
#include <map>
#include <list>
#include <unordered_map>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <Eigen/Dense>

using namespace std;

class CCHCalc;
class CCHArray;
class CHydroBalance;
class CHydroMesh;
class CHyperMesh;
class CEoS;
class CCharge;
class CPart;
class CResList;
class CResInfo;
class CBranchInfo;
class CHyperElement;
class CBalance;
class CAcceptance;
class CAction;
class CB3D;
class CB3DCell;
class CRandy;
class CparameterMap;
class CBalanceArrays;
class CSampler;
class CLocalInfo;
class CAction;
class CRegenerate;
class CSEInfo;
class CHYDROtoB3D;
class CMCList;
class CKernel;
class C3DArray;
class CWaveFunction;
class CSourceCalc;
class CKernelWF;

typedef unordered_map<long int,CResInfo *> CResInfoMap;
typedef pair<long int, CResInfo*> CResInfoPair;
typedef vector<CBranchInfo *> CBranchList; //gives branchlist name
typedef multimap<int,CCharge* > CChargeMap;
typedef pair<int,CCharge* > CChargePair;
typedef multimap<int,CPart* > CPartMap;
typedef pair<int,CPart* > CPartPair;
typedef multimap<int,CCharge* > mapic;
typedef multimap<int,CPart* > mapip;
//typedef array<double,4> FourVector;
typedef double FourVector[4];
typedef multimap<double,CAction *> CActionMap;
typedef pair<double,CAction*> CActionPair;
typedef pair<int,CCharge* > pairic;
typedef pair<int,CCharge* > pairip;

#endif
