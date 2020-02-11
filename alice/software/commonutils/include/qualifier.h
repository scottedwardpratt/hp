#ifndef __QUALIFIER_H__
#define __QUALIFIER_H__
#include "commondefs.h"

using namespace std;

class CQualifiers{
public:
	int nqualifiers;
	string qualifier[100];
	int npars[100];
	string type[100][10],parname[100][10],value[100][10];
	void Read(string qfilename);
	void SetPars(CparameterMap *pmap,int iqualifier);
};

#endif
