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

template <class T>
class Chromosome{

protected:
    T * chromosome;
    uint64_t length;
    long double fitness;

public:
    T * get_chromosome(){ return this->chromosome; }
    long double * get_fitness(){ return &this->fitness; }
    void set_fitness(long double f){ this->fitness = f; }
    virtual void compute_fitness() = 0;
    void set_allele(uint64_t index, T * value); //Makes a hard copy of value
    T * get_allele(uint64_t index); //Returns pointer to be modified
    ~Chromosome();
};


template <class T>
class Chromo_rucksack : public Chromosome<T> {
public:
    Chromo_rucksack(uint64_t alleles);
    void compute_fitness();
};


