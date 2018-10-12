
CXX = g++
CXX_FLAGS = -std=c++11 -O3 -Wall -Wextra -Wpedantic -Wno-unused-function \
            -Wno-unused-const-variable -fvisibility=hidden

CXX_GCC = g++
CXX_FLAGS_GCC = -std=c++11 -O3 -Wall -Wextra -Wpedantic -Waggregate-return \
                -Wcast-align -Wcast-qual -Wconversion -Wconversion-null \
                -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond \
                -Wextra-semi -Wfloat-conversion -Wfloat-equal -Wformat=2 \
                -Wformat-nonliteral -Wformat-security -Wformat-signedness \
                -Winline -Winit-self -Winvalid-pch -Wlogical-op \
                -Wmissing-format-attribute -Wmissing-include-dirs \
                -Wnull-dereference -Wold-style-cast -Wopenmp-simd \
                -Woverloaded-virtual -Wpadded -Wpacked -Wpointer-arith \
                -Wredundant-decls -Wrestrict -Wshadow -Wshift-overflow=2 \
                -Wsign-compare -Wsign-conversion -Wswitch-default \
                -Wswitch-enum -Wsynth -Wtrampolines -Wunused-local-typedefs \
                -Wundef -Wunreachable-code -Wunsafe-loop-optimizations \
                -Wuseless-cast -Wunused -Wvarargs -Wvla \
                -Wvector-operation-performance -Wzero-as-null-pointer-constant

CXX_CLANG = clang++
CXX_FLAGS_CLANG = -std=c++11 -O3 -Wall -Wextra -Wpedantic -Wold-style-cast \
                  -Wshorten-64-to-32 -Wsign-conversion -Wpadded \
                  -Wzero-as-null-pointer-constant -Wnull-dereference \
                  -Wdouble-promotion -Wshadow -Wformat=2 -Wfloat-equal \
                  -Wshift-overflow -Wconversion -Wunused-const-variable

OUTPUTNAME = popdyn

all: popdyn.cpp Makefile
	$(CXX) $(CXX_FLAGS) -o $(OUTPUTNAME).exe popdyn.cpp

all-gcc: popdyn.cpp Makefile
	$(CXX_GCC) $(CXX_FLAGS_GCC) -o $(OUTPUTNAME)-gcc.exe popdyn.cpp

all-clang: popdyn.cpp Makefile
	$(CXX_CLANG) $(CXX_FLAGS_CLANG) -o $(OUTPUTNAME)-clang.exe popdyn.cpp

clean:
	rm -f $(OUTPUTNAME).exe $(OUTPUTNAME)-gcc.exe $(OUTPUTNAME)-clang.exe
