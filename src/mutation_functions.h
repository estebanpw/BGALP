#pragma once


#include "manager.h"
#include "chromosome.h"
#include "population.h"
#define __STDC_FORMAT_MACROS
template <class T>
void mutation_function_TSP(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p);

template <class T>
void mutation_function_LB(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p);

template <class T>
void simple_mutation(Chromo_VRP<T> * p1, uint64_t n_nodes, std::default_random_engine * g, std::uniform_real_distribution<double> * u_r);