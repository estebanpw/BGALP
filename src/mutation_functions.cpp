#include "mutation_functions.h"
#define __STDC_FORMAT_MACROS
template <class T>
void mutation_function_TSP(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p){
    
    
    for(uint64_t i=0;i<a->get_length();i++){
        if(m->u_d(m->uniform_generator) <= p){
            uint64_t swap_pos = i;
            while(swap_pos == i) swap_pos = (uint64_t)(a->get_length())*m->u_d(m->uniform_generator);
            // Warning: only for structures with non pointers
            T aux = *a->get_allele(i);
            *a->get_allele(i) = *a->get_allele(swap_pos);
            *a->get_allele(swap_pos) = aux;
            
        }
    }
}

template void mutation_function_TSP<uint64_t>(Chromosome<uint64_t> * a, Population<uint64_t> * pop, Manager<uint64_t> * m, long double p);
