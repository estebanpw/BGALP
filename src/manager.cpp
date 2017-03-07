#include "manager.h"


template <class T>
Manager<T>::Manager(uint64_t n_populations, void (*cf)(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m), void * solution_info, bool maximize){
    
    this->maximize = maximize;
    this->n_populations = n_populations;
    this->population = (Population<T> **) std::malloc(n_populations*sizeof(Population<T> *));
    if(this->population == NULL) throw "Could not allocate population pointers";
    this->crossover_function = cf;
    this->u_d = std::uniform_real_distribution<long double>(0.0, 1.0); 
    this->pair = (Pair<Chromosome<T>> *) std::malloc(n_populations*sizeof(Pair<Chromosome<T>>));
    this->solution_info = solution_info;
}

template <class T>
void Manager<T>::set_populations(Population<T> * pop, uint64_t index){
    this->population[index] = pop;
}

template <class T>
Chromosome<T> * Manager<T>::get_best_individual(){
    Chromosome<T> * best_c = this->population[0]->get_individual_at(this->population[0]->get_best());
    for(uint64_t i=0;i<this->n_populations;i++){
        if(this->population[i]->get_individual_at(this->population[i]->get_best())->get_fitness() > best_c->get_fitness()){
            best_c = this->population[i]->get_individual_at(this->population[i]->get_best());
        }
    }
    return best_c;
}

template <class T>
void Manager<T>::update_fitnesses(uint64_t curr_population, uint64_t curr_indv){

    if(this->maximize){
        if(this->population[curr_population]->get_individual_at(this->population[curr_population]->get_best())->get_fitness() < this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_best(curr_indv);
        }
        if(this->population[curr_population]->get_individual_at(this->population[curr_population]->get_worst())->get_fitness() > this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_worst(curr_indv);
        }
    }else{
        if(this->population[curr_population]->get_individual_at(this->population[curr_population]->get_best())->get_fitness() > this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_best(curr_indv);
        }
        if(this->population[curr_population]->get_individual_at(this->population[curr_population]->get_worst())->get_fitness() < this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_worst(curr_indv);
        }
    }
}

template <class T>
void Manager<T>::run(uint64_t n_itera){
    uint64_t i, j;

    // Initialize all populations
    for(i=0;i<n_populations;i++){
        for(j=0;j<this->population[i]->get_size();j++){
            this->population[i]->get_individual_at(j)->compute_fitness(this->solution_info);
            update_fitnesses(i,j);
        }
    }

    Chromosome<T> * replacement;

    // Run generations
    for(i=1;i<n_itera;i++){
        for(j=0;j<n_populations;j++){

            // Selection method here
            this->select_tournament(this->population[j], this->population[j], &this->pair[j]);
            replacement = this->population[j]->get_individual_at(this->population[j]->get_worst());
            // Crossover function here
            this->crossover_function(this->pair[j]._e1, this->pair[j]._e2, replacement, this);

            replacement->compute_fitness(&this->solution_info);
            update_fitnesses(j, this->population[j]->get_worst());
            printf("%Le\n", *replacement->get_fitness()); getchar();
        }
    }
}

template <class T>
void Manager<T>::select_tournament(Population<T> * p1, Population<T> * p2, Pair<Chromosome<T>> * c_pair){
    
    c_pair->_e1 = p1->get_individual_at(p1->get_size()*u_d(uniform_generator));
    c_pair->_e2 = c_pair->_e1; //In case p1 and p2 are the same pop
    
    while(p2 != p1) c_pair->_e2 = p2->get_individual_at(p2->get_size()*u_d(uniform_generator));
}

template <class T>
Manager<T>::~Manager(){
    for(uint64_t i=0;i<n_populations;i++){
        delete population[i];
    }
    std::free(population);
    std::free(pair);

}

template class Manager<unsigned char>;