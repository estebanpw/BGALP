#include "crossover_functions.h"

template <class T>
void single_point_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m){
    uint64_t midpoint1 = a->get_length()*m->u_d(m->uniform_generator);
    uint64_t midpoint2 = (a->get_length()-midpoint1)*m->u_d(m->uniform_generator);
    midpoint2 += midpoint1;
    uint64_t i;

    for(i=0;i<midpoint1;i++){
        replacement->set_allele(i, a->get_allele(i));
    }
    for(i=midpoint1;i<midpoint2;i++){
        replacement->set_allele(i, b->get_allele(i));
    }
    for(i=midpoint2;i<replacement->get_length();i++){
        replacement->set_allele(i, a->get_allele(i));
    }
}

template <class T>
void ordered_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m){
    
    uint64_t low = (a->get_length())*m->u_d(m->uniform_generator);
    uint64_t high = low + (a->get_length()+1-low)*m->u_d(m->uniform_generator);
    
    for(uint64_t i=0;i<a->get_length();i++){
        m->marks[i] = 0;
    }
    
    /*
    if(high == a->get_length()){ std::cout << "it does reach high "<< high << std::endl;  getchar();}
    if(low == 0){ std::cout << "it does low " << low << std::endl;  getchar();}
    */

    for(uint64_t i=0;i<low;i++){
        replacement->set_allele(i, a->get_allele(i));
        m->marks[*a->get_allele(i)] = 1;
    }
    
    for(uint64_t i=high;i<a->get_length();i++){
        replacement->set_allele(i, a->get_allele(i));
        m->marks[*a->get_allele(i)] = 1;
    }

    
    
    uint64_t pos = low;
    //Copy the rest from the other
    
    for(uint64_t i=0;i<a->get_length();i++){
        if(m->marks[*b->get_allele(i)] == 0){
            //printf("pos: %" PRIu64"\n", pos);
            replacement->set_allele(pos, b->get_allele(i));
            m->marks[*b->get_allele(i)] = 1;
            pos++;
        }
    }
}

template <class T>
void fill_edge_table(Chromosome<T> * a, Edge_T<T> ** e_table, memory_pool * mp){
    
    uint64_t i;
    Edge_T<T> * previous;
    uint64_t current_allele, next_allele, previous_allele;
    // Add edges from circuit
    for(i=0;i<a->get_length();i++){
        
        current_allele = *a->get_allele(i);
        
        if(e_table[current_allele] == NULL){
            // Add node
            e_table[current_allele] = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
        }

        if(i > 0){ // The previous adjacent

            next_allele = *a->get_allele(i-1);
            Edge_T<T> * et_ptr = e_table[current_allele]->next;
            previous = NULL;
            while(et_ptr != NULL && et_ptr->node != next_allele){

                previous = et_ptr;
                et_ptr = et_ptr->next;
            }

            // Either we have a null node (vertex not added yet) or we found it from the other genome
            if(et_ptr == NULL){
                // Its a NULL one (vertex not added)
                et_ptr = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                et_ptr->node = next_allele;
                et_ptr->common = UNCOMMON;
                et_ptr->next = NULL;
                if(previous != e_table[current_allele] && previous != NULL) previous->next = et_ptr; else e_table[current_allele]->next = et_ptr;

            }else{
                // It exists, so its a common one
                et_ptr->common = COMMON;
            }


        }
        if(i < a->get_length()-1){ // The next adjacent

            previous_allele = *a->get_allele(i+1);
            Edge_T<T> * et_ptr = e_table[current_allele]->next;
            previous = NULL;
            while(et_ptr != NULL && et_ptr->node != previous_allele){

                previous = et_ptr;
                et_ptr = et_ptr->next;
            }

            // Either we have a null node (vertex not added yet) or we found it from the other genome
            if(et_ptr == NULL){
                // Its a NULL one (vertex not added)
                et_ptr = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                et_ptr->node = previous_allele;
                et_ptr->common = UNCOMMON;
                et_ptr->next = NULL;
                if(previous != e_table[current_allele] && previous != NULL) previous->next = et_ptr; else e_table[current_allele]->next = et_ptr;

            }else{
                // It exists, so its a common one
                et_ptr->common = COMMON;
            }

        }
    }
    
}

template void single_point_crossover<unsigned char>(Chromosome<unsigned char> * a, Chromosome<unsigned char> * b, Chromosome<unsigned char> * replacement, Manager<unsigned char> * m);
template void ordered_crossover<uint64_t>(Chromosome<uint64_t> * a, Chromosome<uint64_t> * b, Chromosome<uint64_t> * replacement, Manager<uint64_t> * m);
template void fill_edge_table(Chromosome<uint64_t> * a, Edge_T<uint64_t> ** e_table, memory_pool * mp);