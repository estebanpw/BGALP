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


template <class T> class Manager;

template <class T>
class Manager{

private:
    Pair<T> * pair;
    bool maximize;
    //TODO parametrize the selection function as well

protected:
    Population<T> ** population;
    uint64_t n_populations;
    void * solution_info;
    void (*crossover_function)(T * a, T * b, T * replacement);
    // Uniform generator for population
    std::default_random_engine uniform_generator;
    std::uniform_int_distribution<uint64_t> u_d_p;
    // Uniform generator for chromosome length
    std::uniform_int_distribution<uint64_t> u_d_c;



public:
    Manager(uint64_t n_populations, void (*cf)(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager * m), void * solution_info, bool maximize);
    void set_populations(Population<T> * p, uint64_t index);
    void select_tournament(Population<T> * p1, Population<T> * p2, Pair<T> * c_pair);
    void run(uint64_t n_itera);
    uint64_t get_uniform_random_individual(){ return u_d_p(uniform_generator); }
    uint64_t get_uniform_random_allele(){ return u_d_c(uniform_generator); }
    void update_fitnesses(uint64_t curr_population, uint64_t curr_indv);
    ~Manager();

};

