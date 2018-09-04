
CXX = g++
CXX_FLAGS = -std=c++17 -O3
OUTPUTNAME = popdyn

all: popdyn.cpp Makefile
	$(CXX) $(CXX_FLAGS) -o $(OUTPUTNAME) popdyn.cpp

clean:
	rm -f $(OUTPUTNAME)

