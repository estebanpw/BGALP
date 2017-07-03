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

template <class T>
void mutation_function_LB(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p){
    return;
}

template <class T>
void simple_mutation(Chromo_VRP<T> * p1, uint64_t n_nodes, std::default_random_engine * g, std::uniform_real_distribution<double> * u_r){
    
    uint64_t n1 = (uint64_t)((double) n_nodes * ((*u_r)(*g)));
    uint64_t n2 = (uint64_t)((double) n_nodes * ((*u_r)(*g)));

    while(n1 == n2) n2 = (uint64_t) ((double)n_nodes * (*u_r)(*g));

    T aux = *p1->get_allele(n1);
    p1->set_allele(n1, p1->get_allele(n2));
    p1->set_allele(n2, &aux);

    /*

    if(((*u_r)(*g)) < 0.5){
        // Swap inside nodes
        uint64_t n1 = (uint64_t)((double) n_nodes * ((*u_r)(*g)));

        while(*p1->get_allele(n1) == 0) n1 = (uint64_t) ((double)n_nodes * (*u_r)(*g));
        

        uint64_t nleft = n1, nright = n1;
        while(nleft > 0 && *p1->get_allele(nleft) != 0){ // until left depot
            if(nleft == 0) nleft = n_nodes - 1; else nleft = (nleft - 1) % n_nodes;
        }
        while(nright < n_nodes && *p1->get_allele(nright) != 0){ // until right depot
            nright = (nright+1) % n_nodes;
        }

        uint64_t n2 = nleft + (uint64_t)((double) (nright-nleft) * ((*u_r)(*g)));
        while(n2 == n1) n2 = nleft + (uint64_t)((double) (nright-nleft) * ((*u_r)(*g)));

        // swap them now

        T aux = *p1->get_allele(n1);
        p1->set_allele(n1, p1->get_allele(n2));
        p1->set_allele(n2, &aux);

    }else{
        // Swap a 0 to left or right
        
        uint64_t n1 = (uint64_t)((double) n_nodes * ((*u_r)(*g)));
        

        uint64_t nright = n1;
        while(*p1->get_allele(nright) != 0){ 
            nright = (nright+1) % n_nodes;
        }

        if(((*u_r)(*g)) < 0.5){
            
            //p1->print_chromosome();
            
            T aux = *p1->get_allele(nright - 1);
            p1->set_allele(nright-1, p1->get_allele(nright));
            p1->set_allele(nright, &aux);            


            //p1->print_chromosome();
            //std::cout<<"change left at " << nright << "\n";
            
            //getchar();


        }else{
            if(nright + 1 < n_nodes){
                //p1->print_chromosome();
                T aux = *p1->get_allele(nright + 1);
                p1->set_allele(nright+1, p1->get_allele(nright));
                p1->set_allele(nright, &aux);            
                //p1->print_chromosome();
                //std::cout<<"change right at " << nright << "\n";
                //getchar();
            }
            
        }
        
        
        

        

    }
    */

    

}


template void mutation_function_TSP<uint64_t>(Chromosome<uint64_t> * a, Population<uint64_t> * pop, Manager<uint64_t> * m, long double p);
template void mutation_function_TSP<unsigned char>(Chromosome<unsigned char> * a, Population<unsigned char> * pop, Manager<unsigned char> * m, long double p);
template void mutation_function_LB(Chromosome<unsigned char> * a, Population<unsigned char> * pop, Manager<unsigned char> * m, long double p);
template void simple_mutation(Chromo_VRP<uint64_t> * p1, uint64_t n_nodes, std::default_random_engine * g, std::uniform_real_distribution<double> * u_r);