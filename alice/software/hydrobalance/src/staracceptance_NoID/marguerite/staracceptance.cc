class Tstaracceptance{
public:
  int NVERTEXBINS;
  double zvertexmin,zvertexmax,delzvertex;
  double *zvertexweight;
  int NETABINS;
  double etamin,etamax,deleta;
  int icentrality; // 0,1,2,3 for the four centrality bins: 0 is central
  double **acceptprob; // will become acceptprob[NVERTEXBINS][NETABINS]
  bool acceptance(double zvertex,double *p); // returns 1 if accepted, 0 otherwise
  void setup(); // Run once at beginning
  double getzvertex();
};

void Tstaracceptance::setup(){
  char filename[100];
  FILE *infile;
  int ivertex,ieta;
  double eta,dndeta; // these are not stored
  double zvertexmaxweight;

  icentrality=0; // most central bin only
  NVERTEXBINS=10;
  zvertexweight=new double[NVERTEXBINS];
  zvertexweight[0]=0.081848;
  zvertexweight[1]=0.096890;
  zvertexweight[2]=0.095670;
  zvertexweight[3]=0.104070;
  zvertexweight[4]=0.107120;
  zvertexweight[5]=0.111120;
  zvertexweight[6]=0.108880;
  zvertexweight[7]=0.103670;
  zvertexweight[8]=0.096755;
  zvertexweight[9]=0.093977;
  zvertexmaxweight=zvertexweight[5];
  for(ivertex=0;ivertex<NVERTEXBINS;ivertex++){
    zvertexweight[ivertex]=zvertexweight[ivertex]/zvertexmaxweight;
  }

  zvertexmin=-75.0;
  zvertexmax=75.0;
  NETABINS=100;
  etamin=-1.5;
  etamax=1.5; 

  acceptprob=new double*[NETABINS];
  for(ivertex=0;ivertex<NVERTEXBINS;ivertex++)
    acceptprob[ivertex]=new double[NETABINS];

  delzvertex=(zvertexmax-zvertexmin)/double(NVERTEXBINS);
  deleta=(etamax-etamin)/double(NETABINS);

  for(ivertex=0;ivertex<NVERTEXBINS;ivertex++){
    printf("Initializing, ivertex=%d\n",ivertex);
    //strcpy(filename,form("../staracceptance/hist/etahist_c%d_v%d.txt\0",
    //		 icentrality,ivertex));
    sprintf(filename,"../staracceptance/marguerite/hist/etahist_c%d_v%d.txt\0",
	    icentrality,ivertex);
    infile=fopen(filename,"r");
    printf("Opening %s\n",filename);
    for(ieta=0;ieta<NETABINS;ieta++){
      fscanf(infile,"%lf %lf %lf",&eta,&dndeta,&acceptprob[ivertex][ieta]);
      //printf("%3.2f ",acceptprob[ivertex][ieta]);
    }
    //printf("\n");
    fclose(infile);
  }
}

bool Tstaracceptance::acceptance(double zvertex,double *p){
  int ieta,ivertex;
  bool answer;
  double eta,pt,pmag;

  pt=p[1]*p[1]+p[2]*p[2];
  pmag=sqrt(p[3]*p[3]+pt);
  pt=sqrt(pt);
  eta=0.5*log((pmag+p[3])/(pmag-p[3]));

  ivertex=int(double(NVERTEXBINS)*(zvertex-zvertexmin)
	      /(zvertexmax-zvertexmin));
  ieta=int((eta-etamin)/deleta);
  answer=0;
  if(ieta>=0 && ivertex < NVERTEXBINS && ieta < NETABINS){
    if(pmag<700.0 && fabs(eta)<1.3 && pt>100.0){
      if(ran2()<acceptprob[ivertex][ieta]) answer=1;
    }
  }
  return answer;
}

double Tstaracceptance::getzvertex(){
  int ivertex;
  double zvertex;
 TRY_AGAIN:
  zvertex=zvertexmin+(zvertexmax-zvertexmin)*ran2();
  //zvertex=-25.0+50.0*ran2();
  ivertex=int(double(NVERTEXBINS)*(zvertex-zvertexmin)
	      /(zvertexmax-zvertexmin));
  if(zvertexweight[ivertex]<ran2()) goto TRY_AGAIN;
  return zvertex;
}
