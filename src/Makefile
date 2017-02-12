CC=gcc
CXX=g++ -std=c++11
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DVERBOSE
BIN=../bin

all: GAL

GAL: GAL.cpp
	$(CXX) $(CFLAGS) common_functions.cpp -lm chromosome.cpp population.cpp manager.cpp crossover_functions.cpp GAL.cpp -lpthread -o $(BIN)/GAL

clean:
	rm -rf $(BIN)/GAL
