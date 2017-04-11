CC=gcc
CXX=g++ -std=c++11
CFLAGS=-g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -Wno-literal-suffix #-DVERBOSE
BIN=../bin

all: GAL

GAL: GAL.cpp
	$(CXX) $(CFLAGS) memory_pool.cpp chromosome.cpp local_search_functions.cpp mutation_functions.cpp readstream.cpp -lm common_functions.cpp -lm population.cpp manager.cpp crossover_functions.cpp GAL.cpp -lpthread -o $(BIN)/GAL

GAL_PX: GAL_PX.cpp
	$(CXX) $(CFLAGS) memory_pool.cpp chromosome.cpp local_search_functions.cpp mutation_functions.cpp readstream.cpp -lm common_functions.cpp -lm population.cpp manager.cpp crossover_functions.cpp GAL_PX.cpp -lpthread -o $(BIN)/GAL_PX

clean:
	rm -rf $(BIN)/GAL $(BIN)/GAL_PX
