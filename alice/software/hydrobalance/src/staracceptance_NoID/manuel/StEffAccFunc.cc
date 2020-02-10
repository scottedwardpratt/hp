//
// StEffAccFunc.cc
//
// MCBS
//
// Class to return the values of the Star Efficiency and Acceptance
// for a given pt and a given eta.
//
// The constructor takes as input the name of the file where it should read
// the data.
//
// The file names follow the convention:
// - efficiency or acceptance
// - Particle Name, options are: PiPlus, PiMinus, KPlus, KMinus, Proton, Pbar
// - Multiplicity region, options are: Hi, Me, Lo.
// - end with the extension ".txt"

// This gives e.g. efficiencyPiPlusHi.txt
// 
// The files are simple ascii files where each line is
// eta_index   pt_index   value    error
//
// In general, the eta indices go from 0-19 and the pt indices range from 0-29
// the actual values are -1 < eta < 1 and 0 < pt < 3 GeV/c
//
//
// The Multiplicity regions correspond to the following ranges of the Au+Au hadronic cross section:
// These correspond to the 130 GeV data, with a magnetic field in STAR of ~0.25 Tesla.
// Hi : 0 - 20%   (i.e. central events)
// Me : 20 - 50%
// Lo : 50 - 80%  (i.e. peripheral events)
 
#include "StEffAccFunc.h"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <cmath>

StEffAccFunc::StEffAccFunc(const char* file) {

    // initialize array to zeros
    for (size_t i = 0; i<20; ++i)
	for (size_t j = 0; j<30; ++j) {
	    mData[i][j] = 0;
	}

    // read actual data from file
    ifstream ifs(file);
    assert(ifs);
    int eta_index, pt_index;
    double value, error;
    while (!ifs.eof()) {
	ifs >> eta_index >> pt_index >> value >> error;
	mData[eta_index][pt_index] = value;

	// I'm not using an array for the error, but if needed, it can be
	// treated just like the mData array.
    }
    ifs.close();
    
}

double StEffAccFunc::operator()(double eta, double pt) {
    if (fabs(eta)>1.0) {
      //cout << "Eta out of range " << eta << endl;
	return 0;
    }
    if (pt>3.0 || pt<0) {
      //cout << "pT out of range " << pt << endl;
	return 0;
    }
    eta+=1.0;
    size_t eta_index = static_cast<int>(eta/0.1);
    size_t pt_index  = static_cast<int>(pt/0.1);
    if (eta_index>=20) {
	// this should not happen if the eta_index
	// is properly done! But let's be safe...
	cout << "Eta index out of range! " << eta_index << " " << eta << endl;
	return 0;
    }
    if (pt_index>=30) {
	// ditto
	cout << "pT index out of range! " << pt_index << " " << pt << endl;
	return 0;
    }
    return mData[eta_index][pt_index];
}
