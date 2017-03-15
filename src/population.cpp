#include "population.h"


template <class T>
Population<T>::Population(uint64_t n_individuals, Chromosome<T> * individuals){
    
    this->individuals = individuals;
    this->n_individuals = n_individuals;
    this->index_best = 0;
    this->index_worst = 0;
    //Copy pointers for fast swapping
    this->ptr_individuals = (Chromosome<T> **) std::malloc(n_individuals*sizeof(Chromosome<T> *));
    for(uint64_t i=0;i<n_individuals;i++){
        this->ptr_individuals[i] = &individuals[i];
    }
}

template <class T>
void Population<T>::set_panmictic(){
    for(uint64_t i=0;i<this->n_individuals;i++){
        this->ptr_individuals[i]->set_position(Position());
    }
}

template <class T>
void Population<T>::replace_worst(Chromosome<T> * replacement){
    Chromosome<T> * aux;
    aux = this->ptr_individuals[this->get_worst()];
    this->ptr_individuals[this->get_worst()] = replacement;
    replacement = aux;
}

template <class T>
void Population<T>::replace_best(Chromosome<T> * replacement){
    Chromosome<T> * aux;
    aux = this->ptr_individuals[this->get_best()];
    this->ptr_individuals[this->get_best()] = replacement;
    replacement = aux;
}

template <class T>
Population<T>::~Population(){
    
    if(this->individuals != NULL) std::free(this->ptr_individuals);

}

template <class T>
void Population<T>::set_individual_position(uint64_t ith, Position p){
    this->ptr_individuals[ith]->set_position(p);
}

template <class T>
void Population<T>::set_neighborhood_function(bool (*dst)(Position * p1, Position * p2)){
    this->distance_function = dst;
}

template <class T>
bool Population<T>::is_in_neighborhood(uint64_t i1, uint64_t i2){
    return distance_function(this->ptr_individuals[i1]->get_position(), this->ptr_individuals[i2]->get_position());
}

template <class T>
void Population<T>::print_all_fitness(){
    for(uint64_t i=0; i<this->n_individuals; i++){
        fprintf(stdout, "@[%" PRIu64"]%Le\n", i, *this->get_individual_at(i)->get_fitness());
    }
}

template class Population<double>;
template class Population<unsigned char>;
template class Population<uint64_t>;