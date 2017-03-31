#pragma once


#include "manager.h"
#include "chromosome.h"
#include "population.h"
#define __STDC_FORMAT_MACROS
template <class T>
void mutation_function_TSP(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p);