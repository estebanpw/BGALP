#pragma once


#include "chromosome.h"
#include "manager.h"

template <class T>
void single_point_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager * m){
    uint64_t midpoint = m->get_uniform_random_allele();
    uint64_t i;
    for(i=0;i<midpoint;i++){
        replacement->set_allele(i, a->get_allele(i));
    }
    for(i=midpoint;i<replacement->get_length();i++){
        replacement->set_allele(i, a->get_allele(i));
    }
}