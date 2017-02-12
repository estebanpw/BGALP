#include "chromosome.h"


// Basic chromosome methods
template <class T>
void Chromosome<T>::set_allele(uint64_t index, T * value){
    
    this->chromosome[index] = *value;
    
}

template <class T>
T * Chromosome<T>::get_allele(uint64_t index){
    if(index < this->length) return &this->chromosome[index]; else throw "Index out of bounds";
}

template <class T>
Chromosome<T>::~Chromosome(){
    if(this->chromosome != NULL) std::free(this->chromosome);
}


// Rucksack chromosome (not done)
template <class T>
Chromo_rucksack<T>::Chromo_rucksack(uint64_t alleles, Position p){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    this->length = alleles;
    this->fitness = 0;
    this->position = p;
}


template <class T>
void Chromo_rucksack<T>::compute_fitness(void * solution_info){

}


// Subset sum chromosome
template <class T>
Chromo_subsetsum<T>::Chromo_subsetsum(uint64_t alleles, Position p){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    this->length = alleles;
    this->fitness = 0;
    this->position = p;
}


template <class T>
void Chromo_subsetsum<T>::compute_fitness(void * solution_info){
    Sol_subsetsum * ss = (Sol_subsetsum *) solution_info;
    int64_t t_sum = 0;
    for(uint64_t i=0; i<this->length; i++){
        t_sum +=  (int64_t) this->chromosome[i] * ss->values[i];
    }
    this->fitness = (ss->c - t_sum) * (ss->c - t_sum);
}



template class Chromosome<double>;
template class Chromo_rucksack<double>;
template class Chromo_subsetsum<unsigned char>;