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

template void random_shuffle_templated<uint64_t>(uint64_t n_elements, uint64_t * vector, uint64_t seed, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d);