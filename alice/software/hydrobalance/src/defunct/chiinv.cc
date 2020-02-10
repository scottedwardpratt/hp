#include "hydro2uds.h"

using namespace std;

void CEoS::CalcChiInv(double vOmega,Eigen::Matrix3d &chiinv){
	int a,b;
	CResInfo *resinfo;
	double chiomega;
	Eigen::Matrix3d chi;
	Eigen::Vector3d q;
	chi.setZero();
	CResInfoMap::iterator itr;
	for(itr=reslist->resmap.begin();itr!=reslist->resmap.end();++itr){
		resinfo=itr->second;
		q(0)=resinfo->q[0]; q(1)=resinfo->q[1]; q(2)=resinfo->q[2];
		if(resinfo->baryon!=0 || resinfo->charge!=0 || resinfo->strange!=0){		
			//chiomega=resinfo->ChiOmega(vOmega,Tf);
			chiomega=densityf[resinfo->ires];
			for(a=0;a<3;a++)
				for(b=0;b<=a;b++)
					chi(a,b)+=chiomega*q(a)*q(b);
		}
	}
	for(a=0;a<2;a++)
		for(b=a+1;b<3;b++)
			chi(a,b)=chi(b,a);
	if(vOmega<0.999)
		chiinv=chi.inverse();
	else
		chiinv.setZero();
}
