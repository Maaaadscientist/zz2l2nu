EXEC     = runHZZanalysis\
	#Add other executables here if needed...
CXX      = g++
CXXFLAGS = -Iinclude $(shell root-config --cflags) $(shell root-config --libs) -O2
SRC      = $(filter-out $(wildcard src/*Looper.cc), $(wildcard src/*.cc)) #All the analysis looper should end with Looper.cc
OBJ      = $(SRC:.cc=.o)

all: $(EXEC)

runHZZanalysis: runHZZ2l2nu.cc $(OBJ)
	$(CXX) $< $(CXXFLAGS) $(OBJ) -o $@

%.o: %.cc %.h
	$(CXX) -c $< $(CXXFLAGS) -o $@

.PHONY: clean mrproper

clean:
	rm -f *.o *.d src/*.o

mrproper: clean
	rm -f $(EXEC)
