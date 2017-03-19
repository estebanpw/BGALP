#pragma once


#include "manager.h"
#include "chromosome.h"
#include "memory_pool.h"
#include "defs.h"





template <class T>
void single_point_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m);

template <class T>
void ordered_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m);

template <class T>
void fill_edge_table(Chromosome<T> * a, Edge_T<T> ** e_table, memory_pool * mp);