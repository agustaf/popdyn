
CXX = g++
CXX_FLAGS = -O3
OUTPUTNAME = popdyn

all:
	$(CXX) $(CXX_FLAGS) -o $(OUTPUTNAME) popdyn.cpp

clean:
	rm -f $(OUTPUTNAME)

