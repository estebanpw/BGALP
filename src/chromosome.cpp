#include "chromosome.h"


// Basic chromosome methods
template <class T>
void Chromosome<T>::set_allele(uint64_t index, T * value){
    
    this->chromosome[index] = *value;
    
}

template <class T>
T * Chromosome<T>::get_allele(uint64_t index){
    return &this->chromosome[index];
}

template <class T>
void Chromosome<T>::random_bit_fill(){
    std::random_device engine;
    for(uint64_t i=0;i<this->length;i++){
        this->chromosome[i] = engine(); // sizeof(unsigned) * CHAR_BIT random bits
        this->chromosome[i] = (unsigned char) (this->chromosome[i] > 127);
    }
    
}

template <class T>
Chromosome<T>::~Chromosome(){
    if(this->chromosome != NULL) std::free(this->chromosome);
}


// Rucksack chromosome (not done)
template <class T>
Chromo_rucksack<T>::Chromo_rucksack(uint64_t alleles, Position p, INITIALIZER init_type){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    this->length = alleles;
    this->fitness = 0;
    this->position = p;
    if(init_type == RANDOM) this->random_bit_fill();
}


template <class T>
void Chromo_rucksack<T>::compute_fitness(void * solution_info){

}

template <class T>
void Chromo_rucksack<T>::print_chromosome(){

}

// Subset sum chromosome
template <class T>
Chromo_subsetsum<T>::Chromo_subsetsum(uint64_t alleles, Position p, INITIALIZER init_type){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    this->length = alleles;
    this->fitness = 0;
    this->position = p;
    if(init_type == RANDOM) this->random_bit_fill();
}


template <class T>
void Chromo_subsetsum<T>::compute_fitness(void * solution_info){
    Sol_subsetsum * ss = (Sol_subsetsum *) solution_info;
    int64_t t_sum = 0;
    for(uint64_t i=0; i<this->length; i++){
        t_sum +=  ((int64_t) this->chromosome[i]) * ss->values[i];
    }
    //printf("Hyped: %" PRId64"\n", ss->c - t_sum);
    //printf("Before: %" PRId64"\n", ((ss->c - t_sum) * (ss->c - t_sum)));
    this->fitness = (long double) ((ss->c - t_sum) * (ss->c - t_sum));
    //this->fitness = (long double) ((ss->c - t_sum) * (ss->c - t_sum));
    //printf("Update: %.3Le\n", this->fitness);
    
}

template <class T>
void Chromo_subsetsum<T>::print_chromosome(){
    fprintf(stdout, "\t@(%" PRId64", %" PRId64", %" PRId64") L: %" PRIu64"\n", this->position.x, this->position.y, this->position.z, this->length);
    fprintf(stdout, "\tF: %.3Le\n\t", this->fitness);
    for(uint64_t i=0;i<this->length;i++){
        fprintf(stdout, "%d,", this->chromosome[i]);
    }
    fprintf(stdout, "\n");
}


template class Chromosome<double>;
template class Chromosome<unsigned char>;
template class Chromo_rucksack<double>;
template class Chromo_subsetsum<unsigned char>;