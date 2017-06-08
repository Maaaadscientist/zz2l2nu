runHZZanalysis: runHZZ2l2nu.cc HZZ2l2nuLooper.o SmartSelectionMonitor.o
	g++ runHZZ2l2nu.cc -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -Wno-deprecated-declarations -m64 -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/include HZZ2l2nuLooper.o SmartSelectionMonitor.o -o runHZZanalysis

HZZ2l2nuLooper.o: Sources/HZZ2l2nuLooper.cc Includes/HZZ2l2nuLooper.h
	g++ -c Sources/HZZ2l2nuLooper.cc -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -Wno-deprecated-declarations -m64 -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/include -o HZZ2l2nuLooper.o

SmartSelectionMonitor.o: Sources/SmartSelectionMonitor.cc Includes/SmartSelectionMonitor.h
	g++ -c Sources/SmartSelectionMonitor.cc -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -Wno-deprecated-declarations -m64 -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed6/include -o SmartSelectionMonitor.o

clean:
	rm -f *.o *.d
