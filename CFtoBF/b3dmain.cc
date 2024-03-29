#include "b3d.h"
#include "sampler.h"
#include "balancearrays.h"
#include "qualifier.h"
#include "randy.h"
#include "misc.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 4) {
		printf("Usage: b3d run_name ievent0 ieventf\n");
		exit(-1);
  }
	long long int npartstot;
	long long int ncolls=0,nannihilate=0,nregen=0,nbaryons=0,norm;
	int ievent,iqual,nevents,nparts;
	string run_name=argv[1];
	int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	nevents=1+ieventf-ievent0;
	printf("ievent0=%d, ieventf=%d\n",ievent0,ieventf);
	CB3D *b3d=new CB3D(run_name);
	b3d->InitCascade();
	
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		ncolls=0;
		npartstot=0;
		b3d->SetQualifier(qualifiers.qualifier[iqual]->qualname);
		qualifiers.SetPars(&(b3d->parmap),iqual);
		printf("_________________ iqual=%d, nevents=%d ________________\n",iqual,nevents);
		b3d->sampler->ReadHyperElements2D_OSU();
		for(ievent=ievent0;ievent<=ieventf;ievent++){
			printf("------ beginning, ievent=%d --------\n",ievent);
			b3d->Reset();
			b3d->randy->reset(ievent);
			nparts=b3d->sampler->GenHadronsFromHyperSurface(); // Generates particles from hypersurface, each has bid=-1
			b3d->PerformAllActions();
			printf("nparts=%d=?%d\n",int(b3d->PartMap.size()),nparts);
			printf("nscatter=%lld, nbscatter=%lld, nmerges=%lld, ndecays=%lld,  ncellexits=%lld, nregenerate=%lld\n",
			b3d->nscatter,b3d->nbscatter,b3d->nmerge,b3d->ndecay,b3d->nexit,b3d->nregenerate);
			ncolls+=b3d->nscatter+b3d->nmerge;
			nbaryons+=b3d->nbaryons;
			nannihilate+=b3d->nannihilate;
			nregen+=b3d->nregenerate;
			npartstot+=b3d->PartMap.size();
			nparts=b3d->WriteOSCAR(ievent);
			printf("final nparts=%d\n",nparts);
		}
		norm=nevents*b3d->NSAMPLE;
		printf("nparts=%g, <# collisions>=%g \n",double(npartstot)/norm,double(ncolls)/norm);

	}
	return 0;
}
