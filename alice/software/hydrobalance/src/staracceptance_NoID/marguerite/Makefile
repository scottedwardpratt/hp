CXX 		= g++
CPPFLAGS        = -DGNU_GCC
CXXFLAGS	= -g -O
LDFLAGS		= 

PROGRAM         = testefffunc

$(PROGRAM):	$(PROGRAM).cc StEffAccFunc.o
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROGRAM).cc -o $(PROGRAM) StEffAccFunc.o $(LDFLAGS)

StEffAccFunc.o: StEffAccFunc.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) StEffAccFunc.cc -c

clean:
	rm -f $(PROGRAM) core *.o
