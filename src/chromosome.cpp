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
    }
    
}

template <class T>
Chromosome<T>::~Chromosome(){
    if(this->chromosome != NULL) std::free(this->chromosome);
}

template <class T>
void Chromosome<T>::print_chromosome(){
    std::cout << "\t@(" << this->position.x <<  ", " << this->position.y << ", " << this->position.z << ") L: " << this->length << std::endl;
    std::cout << "\tF: " << this->fitness << " at " << this << std::endl;
    for(uint64_t i=0;i<this->length;i++){
        std::cout << this->chromosome[i] << ", ";
    }
    std::cout << std::endl;
}

template <class T>
void Chromosome<T>::write_chromosome(){
    std::cout << "\t@(" << this->position.x <<  ", " << this->position.y << ", " << this->position.z << ") L: " << this->length;
    std::cout << "\tF: " << this->fitness << " at " << this << std::endl;
    std::cout << std::endl;
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

// Subset sum chromosome
template <class T>
Chromo_subsetsum<T>::Chromo_subsetsum(uint64_t alleles, Position p, INITIALIZER init_type){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    this->length = alleles;
    this->fitness = 0;
    this->position = p;
    if(init_type == RANDOM) this->random_bit_fill();
    for(uint64_t i=0;i<this->length;i++){
        this->chromosome[i] = (unsigned char) (this->chromosome[i] > 127);
    }
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

// Chromosome for TSP

// Subset sum chromosome
template <class T>
Chromo_TSP<T>::Chromo_TSP(uint64_t alleles, Position p, INITIALIZER init_type, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    this->length = alleles;
    this->fitness = LDBL_MAX;
    this->position = p;
    uint64_t seed = 0; // Do not fix
    for(uint64_t i=0;i<alleles;i++){ this->chromosome[i] = i; }
    random_shuffle_templated(this->length, this->chromosome, seed, g, u_d);
    #ifdef VERBOSE
    this->print_chromosome();
    #endif
    
}


template <class T>
void Chromo_TSP<T>::compute_fitness(void * solution_info){
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix *) solution_info;
    long double path_sum = 0;
    for(uint64_t i=1; i<this->length; i++){
        path_sum +=  tsp->dist[this->chromosome[i-1]][this->chromosome[i]]; //Distance between node i and node j
    }
    path_sum += tsp->dist[this->chromosome[0]][this->chromosome[this->length-1]]; //Return to initial node
    this->fitness = path_sum;    
}

template <class T>
void Chromo_TSP<T>::verify_chromosome(char * step){
    uint64_t verification[this->length];
    for(uint64_t i=0;i<this->length;i++){
        verification[i] = 0;
    }

    for(uint64_t i=0;i<this->length;i++){
        verification[this->chromosome[i]]++;
    }
    for(uint64_t i=0;i<this->length;i++){
        if(verification[i] != 1){
            std::cout << "Found error at " << i << " at " << step << " having " << verification[i] << std::endl;
            throw "Aborting";
        }
    }
}


template class Chromosome<double>;
template class Chromosome<unsigned char>;
template class Chromosome<uint64_t>;
template class Chromo_rucksack<double>;
template class Chromo_subsetsum<unsigned char>;
template class Chromo_TSP<uint64_t>;