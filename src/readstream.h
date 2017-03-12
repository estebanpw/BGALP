#pragma once
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <inttypes.h>
#include "defs.h"

#define TEXT "rt"
#define BINARY "rb"

class Readstream{
    private:
    FILE * input = NULL;
    void (*reading_function)(FILE * input, void * type_structure);
    void * type_structure;

    public:
    Readstream(const char * file, void (*reading_function)(FILE * input, void * type_structure), void * type_structure);
    void read();
    ~Readstream();
};

void reading_function_TSP(FILE * input, void * type_structure);