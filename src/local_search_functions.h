#pragma once
#include <cfloat>
#include "chromosome.h"
#define __STDC_FORMAT_MACROS
//Struct for PTHREADs

struct two_opt_args{
    Chromosome<uint64_t> * a;
    Chromosome<uint64_t> * b;
    void * solution;
};

void _2optSwap(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, uint64_t i, uint64_t k);
void run_2opt(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, void * solution_info);
void * run_pthreads_two_opt(void * a);