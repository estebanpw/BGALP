CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DVERBOSE
BIN=../bin

all: GAL

GAL: GAL.cpp
	$(CXX) $(CFLAGS) common_functions.cpp -lm chromosome.cpp GAL.cpp -lpthread -o $(BIN)/GAL

clean:
	rm -rf $(BIN)/GAL
