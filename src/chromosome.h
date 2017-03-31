#pragma once

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <iostream>
#include <cfloat>
#include <random>
#include "common_functions.h"
#include "defs.h"
#define __STDC_FORMAT_MACROS

template <class T> class Chromosome;
template <class T> class Chromo_rucksack;
template <class T> class Chromo_subsetsum;
template <class T> class Chromo_TSP;

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
    void random_bit_fill();
    void random_bit_fill(uint64_t max_l);
    void set_allele(uint64_t index, T * value); 
    T * get_allele(uint64_t index);
    Position * get_position(){ return &this->position;}
    void set_position(Position pos){ this->position = pos;}
    uint64_t get_length() { return this->length; }
    void print_chromosome();
    void write_chromosome();
    ~Chromosome();
};


template <class T>
class Chromo_rucksack : public Chromosome<T> {
public:
    Chromo_rucksack(uint64_t alleles, Position p, INITIALIZER init_type);
    void compute_fitness(void * solution_info);
};

template <class T>
class Chromo_subsetsum : public Chromosome<T> {
public:
    Chromo_subsetsum(uint64_t alleles, Position p, INITIALIZER init_type);
    void compute_fitness(void * solution_info);
};

template <class T>
class Chromo_TSP : public Chromosome<T> {
public:
    Chromo_TSP(uint64_t alleles, Position p, INITIALIZER init_type, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d);
    void compute_fitness(void * solution_info);
    void verify_chromosome(char * step);
};

