#include "common_functions.h"


void terror(const char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

template <class T>
void random_shuffle_templated(uint64_t n_elements, T * vector, uint64_t seed, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d){
    
    T aux;
    uint64_t curr_num;
    for(uint64_t i=0;i<n_elements;i++){
        
        aux = vector[i];
        curr_num = u_d->operator()(*g);
        vector[i] = vector[curr_num];
        vector[curr_num] = aux;
    }

}
template <class T>
void restart_edge_tables(uint64_t n_nodes, Edge_T<T> ** e_table, memory_pool * mp){
    for(uint64_t i=0;i<n_nodes;i++){
        e_table[i] = NULL;
    }
    mp->full_reset();
}

template <class T>
void print_edge_tables(uint64_t n_nodes, Edge_T<T> ** e_table){
    Edge_T<T> * et_ptr;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL){

            std::cout << e_table[i]->partition << " -> (" << e_table[i]->degree << ") @" << i << ": ";
            et_ptr = e_table[i]->next;
            while(et_ptr != NULL){
                std::cout << et_ptr->node;
                if(et_ptr->common == COMMON) std::cout << "*";
                std::cout << ", ";
                et_ptr = et_ptr->next;
            }
            std::cout << std::endl;
        }
    }
    std::cout << "# ---------------" << std::endl;
}

template <class T>
bool find_in_vector(std::vector<T> * v, T key){
    for(typename std::vector<T>::iterator it = v->begin() ; it != v->end(); ++it){
        if(*it == key) return true;
    }
    return false;
}

template void random_shuffle_templated<uint64_t>(uint64_t n_elements, uint64_t * vector, uint64_t seed, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d);
template void restart_edge_tables(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, memory_pool * mp);
template void print_edge_tables(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template bool find_in_vector(std::vector<uint64_t> * v, uint64_t key);