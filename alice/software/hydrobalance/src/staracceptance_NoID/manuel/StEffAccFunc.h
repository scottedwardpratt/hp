//
// StEffAccFunc.h
//
// MCBS

#ifndef StEffAccFunc_h
#define StEffAccFunc_h
class StEffAccFunc {
public:
    
    StEffAccFunc(const char*);
    virtual ~StEffAccFunc() {}

    double operator()(double, double);

private:
    double mData[20][30];
};

#endif
