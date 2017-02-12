#include "population.h"


template <class T>
Population<T>::Population(uint64_t n_individuals, Chromosome<T> * individuals){
    
    this->individuals = individuals;
    this->n_individuals = n_individuals;
    this->index_best = 0;
    this->index_worst = 0;
}

template <class T>
void Population<T>::set_panmictic(){
    for(uint64_t i=0;i<this->n_individuals;i++){
        this->individuals[i].set_position(Position());
    }
}

template <class T>
Population<T>::~Population(){
    if(this->individuals != NULL) std::free(this->individuals);
}

template <class T>
void Population<T>::set_individual_position(uint64_t ith, Position p){
    this->individuals[ith].set_position(p);
}

template <class T>
void Population<T>::set_neighborhood_function(bool (*dst)(Position * p1, Position * p2)){
    this->distance_function = dst;
}

template <class T>
bool Population<T>::is_in_neighborhood(uint64_t i1, uint64_t i2){
    return distance_function(this->individuals[i1].get_position(), this->individuals[i2].get_position());
}


template class Population<double>;
template class Population<unsigned char>;