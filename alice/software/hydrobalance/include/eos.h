#ifndef __EoS_H__
#define __EoS_H__

#include "commondefs.h"

using namespace std;
typedef multimap<double,int> mapdi;
typedef pair<double,int> pairdi;

class CEoS{
public:
	CparameterMap *parmap;
	double T,P,epsilon,s;
	double chill,chiud,chils,chiss;
	double chillovers,chiudovers,chilsovers,chissovers;
	CEoS(CparameterMap *parmapset);
	CEoS();
	void ReadDiffusionData();
	// get D in fm
	static double GetD(double T);
		
	void ReadEoS_PST();
	void CalcEoS_PST();
	void GetEoSFromT_PST(double Tset);
	void GetEoSFromEpsilon_PST(double epsilonset);
	
	void ReadChiData_HSC();
	void ReadChiData_Claudia();
	void GetChiOverS_Claudia();
	void Print();
	void PrintChi();
	void FillOutdDdT();
	
	void BuildMap();
	
private:
	static vector<double> epsilon_PST,P_PST,s_PST,T_PST;
	static vector<double> epsilon_claudia,P_claudia,s_claudia,T_claudia;
	static vector<double> twopiTD,Tdiff;
	static vector<double> chill_claudia,chiud_claudia,chils_claudia,chiss_claudia;
	static vector<double> chill_HSC,chiud_HSC,chils_HSC,chiss_HSC;
	static vector<double> dDdT;
	static mapdi etmap;
	static CResList *reslist;

};


#endif
