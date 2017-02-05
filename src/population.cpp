#include "population.h"


template <class T>
Population<T>::Population(uint64_t n_individuals){
    this->individuals = (Chromosome<T> *) std::malloc(n_individuals * sizeof(Chromosome<T>));
    if(this->individuals == NULL) throw "Could not allocate individuals";
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




template class Population<double>;