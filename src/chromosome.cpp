#include "chromosome.h"


template <class T>
void Chromosome<T>::set_allele(uint64_t index, T * value){
    if(index < this->length){
        memcpy(&this->chromosome[index], value, sizeof(T));
    }else{
        throw "Index out of bounds";
    }
}

template <class T>
T * Chromosome<T>::get_allele(uint64_t index){
    if(index < this->length) return &this->chromosome[index]; else throw "Index out of bounds";
}

template <class T>
Chromosome<T>::~Chromosome(){
    std::free(this->chromosome);
}

template <class T>
Chromo_rucksack<T>::Chromo_rucksack(uint64_t alleles){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    this->length = alleles;
    this->fitness = 0;
}


template <class T>
void Chromo_rucksack<T>::compute_fitness(){

}

template class Chromosome<double *>;
template class Chromosome<uint64_t *>;
template class Chromo_rucksack<double *>;
template class Chromo_rucksack<uint64_t *>;