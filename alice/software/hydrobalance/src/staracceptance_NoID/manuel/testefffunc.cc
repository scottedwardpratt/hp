//
// testefffunc.cc
// file to test StEffAccFunc
//
// MCBS
//
#include <iostream>
#include <cstring>
#include <cstdio>
using namespace std;
#include "StEffAccFunc.h"
int main() {
  double pt,eta;
  char filename[80],identstring[80];

  // Create the StEffAccFunc class with the file name as argument.
  // The file names follow the convention are explained in StEffFunc.cc

  sprintf(identstring,"PiMinusHi\0");
  sprintf(filename,"full_field/efficiency%s.txt\0",identstring);
  cout << "filename " << filename << endl;
  StEffAccFunc efficiency(filename);

  printf("     ");
  for(eta=-1.25;eta<0;eta+=0.1)
    printf("%4.2f ",eta);
  printf("\n");
  for (double pt=0.05; pt<2; pt+=0.01) {
    printf("%4.2f: ",pt);
    for(eta=-1.25;eta<0;eta+=0.1)
      printf("%4.3f ",efficiency(eta,pt));
    printf("\n");
  }
    
  return 0;
}

#include "StEffAccFunc.cc"
