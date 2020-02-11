#ifndef __BALANCE_DEFS_H__
#define __BALANCE_DEFS_H__

#include <map>
#include <cmath>
#include <vector>

using namespace std;

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

typedef map<long int,CResInfo *> CResInfoMap;
typedef pair<long int, CResInfo*> CResInfoPair;
typedef vector<CBranchInfo *> CBranchList; //gives branchlist name
typedef multimap<int,CCharge* > mapic;
typedef pair<int,CCharge* > pairic;
typedef multimap<int,CPart* > mapip;
typedef pair<int,CPart* > pairip;
//typedef array<double,4> FourVector;
typedef double FourVector[4];
typedef pair<double,CHyperElement*> pairdh;


#endif