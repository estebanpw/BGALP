#pragma once

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <random>
#include "common_functions.h"
#include "defs.h"
#include "population.h"
#include "chromosome.h"

#define MAXIMIZE true 
#define MINIMIZE false

template <class T> class Manager;

template <class T>
class Manager{

private:
    Pair<Chromosome<T>> * pair;
    bool maximize;
    Chromosome<T> * replacement;
    //TODO parametrize the selection function as well
    Population<T> ** population;
    uint64_t n_populations;
    uint64_t mix_every;
    void * solution_info;
    void (*crossover_function)(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m);
    void (*mutation_function)(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p);
        
public:

    // Uniform generator for population
    std::default_random_engine uniform_generator;
    std::uniform_real_distribution<long double> u_d;
    // For ordered crossover
    unsigned char * marks; 

    Manager(uint64_t n_populations, uint64_t mix_every, void (*cf)(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m), void (*mut)(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p), void * solution_info, bool maximize);
    void set_populations(Population<T> * p, uint64_t index);
    void select_tournament(Population<T> * p1, Population<T> * p2, Pair<Chromosome<T>> * c_pair);
    void run(uint64_t n_itera);
    void generate_marks_for_ordered_crossover(uint64_t n_alleles);
    Chromosome<T> * get_best_individual();
    void update_fitnesses(uint64_t curr_population, uint64_t curr_indv);
    ~Manager();

};

