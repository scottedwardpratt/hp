#include "b3d.h"
#include "qualifier.h"

using namespace std;

int main(int argc, char *argv[]){
	bool FROMCHARGES=true;
	if (argc != 4) {
		printf("Usage: phitest run_name ievent0 ieventf\n");
		exit(-1);
  }
	CBalanceArrays *barray;
	double dnchdy=0.0,tau;
	int nparts,nchargesample=1,isample;
	long long int npartstot;
	long long int ncolls=0,nannihilate=0,nregen=0,nbaryons=0,norm;
	int ievent,iqual,nevents;
	string run_name=argv[1];
	int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	nevents=1+ieventf-ievent0;
	printf("ievent0=%d, ieventf=%d\n",ievent0,ieventf);
	CB3D *b3d=new CB3D(run_name);
	b3d->InitCascade();
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		dnchdy=0;
		ncolls=0;
		npartstot=0;
		b3d->SetQualifier(qualifiers.qualifier[iqual]);
		qualifiers.SetPars(&(b3d->parmap),iqual);
		printf("_________________ iqual=%d, nevents=%d ________________\n",iqual,nevents);
		b3d->sampler->ReadHyperElements2D_OSU();
		for(ievent=ievent0;ievent<=ieventf;ievent++){
			printf("------ beginning, ievent=%d --------\n",ievent);
			for(isample=0;isample<nchargesample;isample++){
				b3d->Reset();
				//b3d->randy->reset(-isample-ievent*nchargesample-10000000);
				//nparts=b3d->sampler->GenHadronsFromHyperSurface(); // Generates particles from hypersurface, each has bid=-1
				if(FROMCHARGES){
					b3d->ReadCharges(ievent);
					b3d->AnalyzeCharges();
					//b3d->GenHadronsFromCharges(); // Generates inter-correlated parts, with bids = (0,1),(2,3)....
					b3d->DeleteCharges();
				}							
			}
		}
	}
	for(int iphi=0;iphi<18;iphi++){
		printf("%2d %10d\n",iphi,int(b3d->phicount[iphi]));
	}
	return 0;
}
