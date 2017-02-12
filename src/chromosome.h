#pragma once

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include "common_functions.h"
#include "defs.h"


template <class T> class Chromosome;
template <class T> class Chromo_rucksack;
template <class T> class Chromo_subsetsum;

template <class T>
class Chromosome{

protected:
    T * chromosome;
    uint64_t length;
    long double fitness;
    Position position;

public:
    T * get_chromosome(){ return this->chromosome; }
    long double * get_fitness(){ return &this->fitness; }
    void set_fitness(long double f){ this->fitness = f; }
    virtual void compute_fitness(void * solution_info) = 0;
    void set_allele(uint64_t index, T * value); 
    T * get_allele(uint64_t index); //Returns pointer to be modified
    Position * get_position(){ return &this->position;}
    void set_position(Position pos){ this->position = pos;}
    uint64_t get_length() { return this->length; }
    ~Chromosome();
};


template <class T>
class Chromo_rucksack : public Chromosome<T> {
public:
    Chromo_rucksack(uint64_t alleles, Position p);
    void compute_fitness(void * solution_info);
};

template <class T>
class Chromo_subsetsum : public Chromosome<T> {
public:
    Chromo_subsetsum(uint64_t alleles, Position p);
    void compute_fitness(void * solution_info);
};

