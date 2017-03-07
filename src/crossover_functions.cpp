#include "crossover_functions.h"

template <class T>
void single_point_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m){
    uint64_t midpoint = a->get_length()*m->u_d(m->uniform_generator);
    uint64_t i;
    for(i=0;i<midpoint;i++){
        replacement->set_allele(i, a->get_allele(i));
    }
    for(i=midpoint;i<replacement->get_length();i++){
        replacement->set_allele(i, b->get_allele(i));
    }
}

template void single_point_crossover<unsigned char>(Chromosome<unsigned char> * a, Chromosome<unsigned char> * b, Chromosome<unsigned char> * replacement, Manager<unsigned char> * m);