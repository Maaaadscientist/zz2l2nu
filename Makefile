EXEC     = runHZZanalysis\
	#Add other executables here if needed...
CXX      = g++
CXXFLAGS = $(shell root-config --cflags) $(shell root-config --libs) -O2
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
