#include "charge.h"
#include "misc.h"
#include "hyper.h"
#include "hydro2uds.h"
#include "randy.h"
#include "constants.h"

using namespace std;

CTrajInfo::CTrajInfo(int balanceIDset){
	balanceID=balanceIDset;
	string filename="trajectories/traj"+to_string(balanceID)+".txt";
	fptr=fopen(filename.c_str(),"w");
}

void CTrajInfo::add(double x1,double y1,double eta1,double tau1){
	x.push_back(x1);
	y.push_back(y1);
	eta.push_back(eta1);
	tau.push_back(tau1);
}

