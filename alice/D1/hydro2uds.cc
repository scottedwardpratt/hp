#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "hydro2uds.h"
#include "hyper.h"
#include "qualifier.h"

using namespace std;

int main(int argc,char *argv[]){
	bool oscarfile=true;
	int run_number=atoi(argv[1]);
	string udsfilename="uds"+string(argv[1])+".dat";
	CHydroBalance hb("udsdata/udsparameters.dat",run_number);
	hb.parmap.set("CHARGESINFO_FILENAME",udsfilename);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	printf("nqualifiers=%d\n",qualifiers.nqualifiers);
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		hb.qualifier=qualifiers.qualifier[iqual];
		printf("--------- BEGIN CALC FOR %s ---------\n",hb.qualifier.c_str());
		oscarfile=hb.ReadOSCAR(hb.mesh);
		printf("tau=%g, cmap.size=%d, emap.size=%d\n",
		hb.mesh->tau,int(hb.cmap.size()),int(hb.emap.size()));
		hb.HyperFind();
		oscarfile=hb.ReadOSCAR(hb.newmesh);
		hb.HyperFind();
		hb.MakeCharges();
		printf("tau=%g, cmap.size=%d, emap.size=%d\n",
		hb.newmesh->tau,int(hb.cmap.size()),int(hb.emap.size()));
		hb.PropagateCharges();
		do{
			hb.SwapMeshes();
			oscarfile=hb.ReadOSCAR(hb.newmesh);
			hb.HyperFind();
			hb.MakeCharges();
			hb.PropagateCharges();
			hb.ScatterCharges();
			if(fabs(lrint(hb.mesh->tau)-hb.mesh->tau)<0.001)
				printf("tau=%g, cmap.size=%d, emap.size=%d\n",hb.mesh->tau,int(hb.cmap.size()),int(hb.emap.size()));
		}while(oscarfile);
		hb.WriteCharges();
		if(run_number==0)
			hb.WriteHyper();
		printf("---- ChiTot(ij)/ChiTotHyper(ij)\n");
		Eigen::Matrix3d ctratio(3,3);
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				ctratio(i,j)=hb.chitot(i,j)/hb.chitothyper(i,j);
		cout << ctratio << endl;
	
		printf("------ ChiTot from Volume -------\n");
		cout << hb.chitot << endl;
		hb.chitot=hb.chitot/hb.chitot(0,0);
		cout << hb.chitot << endl;
		printf("------ ChiTot from Hyper-Surf -------\n");
		cout << hb.chitothyper << endl;
		hb.chitothyper=hb.chitothyper/hb.chitothyper(0,0);
		cout << hb.chitothyper << endl;
	
		hb.ClearCharges();
	}
	
  return 0;
}


