#pragma once
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <inttypes.h>
#include "common_functions.h"
#include "defs.h"

#define TEXT "rt"
#define BINARY "rb"
#define __STDC_FORMAT_MACROS
#define INIT_SEQS 10000

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

void reading_function_TSP(FILE * input, void * type_structure); // For symmetric TSP 
void reading_function_LB_reads(FILE * input, void * type_structure); // Load balancing of metagenomic reads