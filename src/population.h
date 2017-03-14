#pragma once

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include "common_functions.h"
#include "defs.h"
#include "chromosome.h"


template <class T> class Population;

template <class T>
class Population{

protected:
    Chromosome<T> * individuals;
    Chromosome<T> ** ptr_individuals;
    uint64_t n_individuals;
    bool (*distance_function)(Position * p1, Position * p2);

    uint64_t index_worst;
    uint64_t index_best;


public:
    Population(uint64_t n_individuals, Chromosome<T> * individuals);
    void set_panmictic();
    uint64_t get_size(){ return this->n_individuals; }
    void set_individual_position(uint64_t ith, Position p);
    void set_neighborhood_function(bool (*dst)(Position * p1, Position * p2));
    bool is_in_neighborhood(uint64_t i1, uint64_t i2);
    Chromosome<T> * get_individual_at(uint64_t index){ return this->ptr_individuals[index]; }
    void set_individual_pointer_to(uint64_t index, Chromosome<T> * new_ptr){ this->ptr_individuals[index] = new_ptr; }
    void set_best(uint64_t index){ this->index_best = index; }
    uint64_t get_best(){ return this->index_best; }
    Chromosome<T> * get_best_individual(){ return this->ptr_individuals[index_best]; }
    Chromosome<T> * get_worst_individual(){ return this->ptr_individuals[index_worst]; }
    void replace_worst(Chromosome<T> * replacement);
    void replace_best(Chromosome<T> * replacement);
    void set_worst(uint64_t index){ this->index_worst = index; }
    uint64_t get_worst(){ return this->index_worst; }
    void print_all_fitness();
    ~Population();

};



