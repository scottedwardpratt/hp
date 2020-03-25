#include "b3d.h"
#include "sampler.h"
#include "balancearrays.h"
#include "qualifier.h"
#include "randy.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 4) {
		printf("Usage: b3d_fromcascade run_name ievent0 ieventf\n");
		exit(-1);
  }
	CBalanceArrays *barray;
	long long int npartstot;
	int ievent,iqual,nevents;
	double norm;
	string run_name=argv[1];
	int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	nevents=1+ieventf-ievent0;
	printf("ievent0=%d, ieventf=%d\n",ievent0,ieventf);
	CB3D *b3d=new CB3D(run_name);
	b3d->InitCascade();
	barray=b3d->balancearrays;
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		npartstot=0;
		b3d->SetQualifier(qualifiers.qualifier[iqual]);
		qualifiers.SetPars(&(b3d->parmap),iqual);
		printf("_________________ iqual=%d, nevents=%d ________________\n",iqual,nevents);
		b3d->sampler->ReadHyperElements2D_OSU();
		for(ievent=ievent0;ievent<ieventf;ievent++){
			printf("------ beginning, ievent=%d --------\n",ievent);
			b3d->Reset();
			b3d->randy->reset(-ievent);
			b3d->ReadOSCAR(ievent);
			//b3d->PerformAllActions();
			printf("nparts=%d\n",int(b3d->PartMap.size()));
			npartstot+=b3d->PartMap.size();
			barray->ProcessPartMap();
		}
		norm=nevents*b3d->NSAMPLE;
		printf("npartstot=%g\n",double(npartstot)/norm);
		barray->ConstructBFs();
		barray->WriteBFs();
		barray->WriteDenoms();
		barray->WriteGammaP();
	}
	return 0;
}
