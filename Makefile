EXEC     = runHZZanalysis\
	#Add other executables here if needed...
CXX      = g++
CXXFLAGS = -g -O2 -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -Wno-deprecated-declarations -m64 -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/include
SRC      = $(filter-out $(wildcard Loopers_Sources/*Looper.cc), $(wildcard Loopers_Sources/*.cc)) $(wildcard Common/*.cc) $(wildcard Tools/*.cc) #All the analysis looper should end with Looper.cc
OBJ      = $(SRC:.cc=.o)

all: $(EXEC)

runHZZanalysis: runHZZ2l2nu.cc $(OBJ)
	$(CXX) $< $(CXXFLAGS) $(OBJ) -o $@

%.o: %.cc %.h
	$(CXX) -c $< $(CXXFLAGS) -o $@

.PHONY: clean mrproper

clean:
	rm -f *.o *.d Common/*.o Loopers_Sources/*.o Tools/*.o

mrproper: clean
	rm -f $(EXEC)
