#include "mutation_functions.h"

template <class T>
void mutation_function_TSP(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p){
    
    
    for(uint64_t i=0;i<a->get_length();i++){
        if(m->u_d(m->uniform_generator) <= p){
            uint64_t swap_pos = i;
            while(swap_pos == i) swap_pos = (a->get_length())*m->u_d(m->uniform_generator);
            pop->swap_individuals(i, swap_pos);
        }
    }
}

template void mutation_function_TSP<uint64_t>(Chromosome<uint64_t> * a, Population<uint64_t> * pop, Manager<uint64_t> * m, long double p);
