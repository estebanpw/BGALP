#include "manager.h"
#define __STDC_FORMAT_MACROS

template <class T>
Manager<T>::Manager(uint64_t n_populations, uint64_t mix_every, void (*cf)(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m), void (*mut)(Chromosome<T> * a, Population<T> * pop, Manager<T> * m, long double p), void * solution_info, bool maximize){
    
    this->maximize = maximize;
    this->mix_every = mix_every;
    this->n_populations = n_populations;
    this->population = (Population<T> **) std::malloc(n_populations*sizeof(Population<T> *));
    if(this->population == NULL) throw "Could not allocate population pointers";
    this->crossover_function = cf;
    this->mutation_function = mut;
    this->u_d = std::uniform_real_distribution<long double>(0.0, 1.0); 
    this->pair = (Pair<Chromosome<T>> *) std::malloc(n_populations*sizeof(Pair<Chromosome<T>>));
    this->solution_info = solution_info;
    this->marks = NULL;
}

template <class T>
void Manager<T>::set_populations(Population<T> * pop, uint64_t index){
    this->population[index] = pop;
}

template <class T>
void Manager<T>::generate_marks_for_ordered_crossover(uint64_t n_alleles){
    this->marks = (unsigned char *) std::malloc(n_alleles*sizeof(unsigned char));
    if(this->marks == NULL) throw "Could not generate marks for ordered crossover";
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
        if(*this->population[curr_population]->get_individual_at(this->population[curr_population]->get_best())->get_fitness() <= *this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_best(curr_indv);
        }
        if(*this->population[curr_population]->get_individual_at(this->population[curr_population]->get_worst())->get_fitness() > *this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_worst(curr_indv);
        }
        if(*this->population[curr_population]->get_individual_at(this->population[curr_population]->get_worst())->get_fitness() > *this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_worst(curr_indv);
        }
    }else{
        if(*this->population[curr_population]->get_individual_at(this->population[curr_population]->get_best())->get_fitness() >= *this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
            this->population[curr_population]->set_best(curr_indv);
        }
        if(*this->population[curr_population]->get_individual_at(this->population[curr_population]->get_worst())->get_fitness() < *this->population[curr_population]->get_individual_at(curr_indv)->get_fitness()){
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
        // In case they coincide
        if(this->population[i]->get_best() == this->population[i]->get_worst()){
            this->population[i]->set_worst((this->population[i]->get_worst() + 1) % this->population[i]->get_size());
        }

        /*
        printf("INIT CONFIG:\n");
        for(j=0;j<this->population[i]->get_size();j++){
            this->population[i]->get_individual_at(j)->print_chromosome();
        }
        printf("Current best f: %.3Le, worst : %.3Le\n", *this->population[i]->get_best_individual()->get_fitness(), *this->population[i]->get_worst_individual()->get_fitness());
        getchar();
        */
    }


    uint64_t replace_pos, k, second_worst;
    uint64_t print_time = (n_itera/20 != 0) ? (n_itera/20) : (1);


    // Run generations
    for(i=1;i<n_itera;i++){

        if(this->n_populations > 1 && i % this->mix_every == 0){
            //printf("mix is %" PRIu64"\n", this->mix_every);
            // Mix two populations randomly
            uint64_t pop1 = (uint64_t) ((long double)this->n_populations)*this->u_d(this->uniform_generator);
            uint64_t pop2 = pop1;
            do{
                pop2 = (uint64_t) ((long double)this->n_populations)*this->u_d(this->uniform_generator);
            }while(pop1 == pop2);

            // Swap them pointers
            uint64_t ptr1 = (uint64_t) (this->population[pop1]->get_size()-1)*this->u_d(this->uniform_generator);
            uint64_t ptr2 = (uint64_t) (this->population[pop2]->get_size()-1)*this->u_d(this->uniform_generator);
            
            Chromosome<T> * an_auxiliary_pointer_1 = this->population[pop1]->get_individual_at(ptr1);
            Chromosome<T> * an_auxiliary_pointer_2 = this->population[pop2]->get_individual_at(ptr2);

            
            this->population[pop1]->set_individual_pointer_to(ptr1, &(*an_auxiliary_pointer_2));
            this->population[pop2]->set_individual_pointer_to(ptr2, &(*an_auxiliary_pointer_1));
            
            
            //std::cout << "\tMix took place between pool " << pop1 << " and " << pop2 << "\n";

        }

        for(j=0;j<n_populations;j++){

            // Selection method here
            replace_pos = this->population[j]->get_worst();
            replacement = this->population[j]->get_individual_at(replace_pos);
            do{
                this->select_tournament(this->population[j], this->population[j], &this->pair[j]);
            }while(this->pair[j]._e1 == replacement || this->pair[j]._e2 == replacement);

            /*
            #ifdef VERBOSE
            // Some info
            fprintf(stdout, "Mated:\n"); this->pair[j]._e1->print_chromosome();
            fprintf(stdout, "with:\n"); this->pair[j]._e2->print_chromosome();
            #endif
            */

            // Crossover function here
            this->crossover_function(this->pair[j]._e1, this->pair[j]._e2, replacement, this);

            // Mutation here
            this->mutation_function(replacement, this->population[j], this, (long double)1/replacement->get_length());

            replacement->compute_fitness(this->solution_info);
            
            /*
            #ifdef VERBOSE
            //printf("looking for %" PRId64"\n", ((Sol_subsetsum *)this->solution_info)->c);
            fprintf(stdout, "Results in:\n"); replacement->print_chromosome();
            //getchar();
            #endif
            */

            // The worst individual is always replaced
            this->population[j]->replace_worst(replacement);

                        

            // Find new worst 
            second_worst = 0;
            for(k=0;k<this->population[j]->get_size();k++){
                if(k != replace_pos){
                    if(this->maximize && *this->population[j]->get_individual_at(second_worst)->get_fitness() >= *this->population[j]->get_individual_at(k)->get_fitness()) second_worst = k;
                    if(!this->maximize && *this->population[j]->get_individual_at(second_worst)->get_fitness() <= *this->population[j]->get_individual_at(k)->get_fitness()) second_worst = k;
                }
            }
            this->population[j]->set_worst(second_worst);
            
            
            // Update best index
            if(this->maximize && *replacement->get_fitness() >= *this->population[j]->get_best_individual()->get_fitness()) this->population[j]->set_best(replace_pos);
            if(!this->maximize && *replacement->get_fitness() <= *this->population[j]->get_best_individual()->get_fitness()) this->population[j]->set_best(replace_pos);
            
            //printf("Current best f: %.3Le, worst : %.3Le\n", *this->population[j]->get_best_individual()->get_fitness(), *this->population[j]->get_worst_individual()->get_fitness());
            
            if(i % print_time == 0){
                //  this->population[j]->print_all_fitness();
                //fprintf(stdout, "I(%" PRIu64") :: %.3Le (%" PRIu64") @%" PRIu64"\n", i, *this->population[j]->get_best_individual()->get_fitness(), (uint64_t)*this->population[j]->get_best_individual()->get_fitness() , this->population[j]->get_best());
                std::cout << "I(" << i<< ") :: " <<  *this->population[j]->get_best_individual()->get_fitness() <<" @" << this->population[j]->get_best()<<"\n"; 
                //getchar();
            }     
        }
    }
}

template <class T>
void Manager<T>::select_tournament(Population<T> * p1, Population<T> * p2, Pair<Chromosome<T>> * c_pair){
    
    c_pair->_e1 = p1->get_individual_at(p1->get_size()*u_d(uniform_generator));
    c_pair->_e2 = c_pair->_e1; //In case p1 and p2 are the same pop
    
    while(c_pair->_e1 == c_pair->_e2) c_pair->_e2 = p2->get_individual_at(p2->get_size()*u_d(uniform_generator));
}

template <class T>
Chromosome<T> ** Manager<T>::retrieve_k_best_solutions(uint64_t k){
    
    Chromosome<T> ** v = (Chromosome<T> **) std::calloc(k, sizeof(Chromosome<T> *));
    if(v == NULL) throw "Could not allocate best solutions";

    uint64_t j, z;
    long double curr_fitness;
    for(uint64_t i=0;i<this->n_populations;i++){
        for(j=0;j<this->population[i]->get_size();j++){

            curr_fitness = *this->population[i]->get_individual_at(j)->get_fitness();

            for(z = 0; z < k; z++){

                if(v[z] != NULL){
                    if(this->maximize && curr_fitness > *v[z]->get_fitness()){
                        v[z] = this->population[i]->get_individual_at(j);
                        break;
                    }
                    if(!this->maximize && curr_fitness < *v[z]->get_fitness()){
                        v[z] = this->population[i]->get_individual_at(j);
                        break;
                    }
                }else{
                    v[z] = this->population[i]->get_individual_at(j);
                    break;
                }
            }
        }
    }
    return v;
}


template <class T>
Manager<T>::~Manager(){
    for(uint64_t i=0;i<n_populations;i++){
        delete population[i];
    }
    if(this->marks != NULL) std::free(this->marks);
    std::free(population);
    std::free(pair);
}

template class Manager<unsigned char>;
template class Manager<uint64_t>;