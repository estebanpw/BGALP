#pragma once
#include <cfloat>
#include "chromosome.h"

void _2optSwap(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, uint64_t i, uint64_t k);
void run_2opt(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, void * solution_info);