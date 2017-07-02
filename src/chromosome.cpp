#include "chromosome.h"
#define __STDC_FORMAT_MACROS


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
        std::cout << (uint64_t) this->chromosome[i] << ", ";
    }
    std::cout << std::endl;
}

template <class T>
void Chromosome<T>::write_chromosome(){
    std::cout << "\t@(" << this->position.x <<  ", " << this->position.y << ", " << this->position.z << ") L: " << this->length;
    std::cout << "\tF: " << this->fitness << " at " << this << std::endl;
    std::cout << std::endl;
}

template <class T>
void Chromosome<T>::hard_copy_no_pointers(Chromosome<T> * source){
    for(uint64_t i=0;i<this->get_length();i++){
        this->chromosome[i] = source->chromosome[i];
        //destination->set_allele(source->get_allele());
    }
    
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
    if(this->chromosome == NULL) throw "Could not allocate chromosome";
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
    if(this->chromosome == NULL) throw "Could not allocate chromosome";
    this->length = alleles;
    this->fitness = LDBL_MAX;
    this->position = p;
    uint64_t seed; // Do not fix
    for(uint64_t i=0;i<alleles;i++){ this->chromosome[i] = i; }
    random_shuffle_templated(this->length, this->chromosome, seed, g, u_d);
        
}


template <class T>
void Chromo_TSP<T>::compute_fitness(void * solution_info){
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix *) solution_info;
    long double path_sum = 0;
    for(uint64_t i=1; i<this->length; i++){
        //std::cout << "Computing from " << this->chromosome[i-1] << " to " << this->chromosome[i] << " adds " << tsp->dist[this->chromosome[i-1]][this->chromosome[i]] << std::endl;
        path_sum +=  tsp->dist[this->chromosome[i-1]][this->chromosome[i]]; //Distance between node i and node j
    }
    //std::cout << "Computing from " << this->chromosome[0] << " to " << this->chromosome[this->length-1] << " adds " << tsp->dist[this->chromosome[0]][this->chromosome[this->length-1]] << std::endl;
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

// CHromosome for VRP

// Subset sum chromosome
template <class T>
Chromo_VRP<T>::Chromo_VRP(uint64_t alleles, uint64_t n_trucks, long double capacity, T depot, Position p, INITIALIZER init_type, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d, void * sol_VRP, uint64_t node_shift){
    this->chromosome = (T *) std::calloc(alleles, sizeof(T));
    this->lookup = (uint64_t *) std::calloc(alleles, sizeof(uint64_t));
    if(this->chromosome == NULL) throw "Could not allocate chromosome";
    this->length = alleles;
    this->fitness = LDBL_MAX;
    this->n_trucks = n_trucks;
    this->capacity = capacity;
    this->position = p;
    this->depot = depot;
    uint64_t seed; // Do not fix
    if(init_type == RANDOM){
        for(uint64_t i=1;i<alleles - (n_trucks - 1);i++){ this->chromosome[i] = i-1; }
        for(uint64_t i=0 ;i<(n_trucks - 1); i++){ this->chromosome[(alleles -(n_trucks - 1)) + i] = this->depot; } // Add depots
        random_shuffle_templated(this->length, this->chromosome, seed, g, u_d); // shuffle them 
    }
    if(init_type == PETALS){
        generate_petals_from_points(this->chromosome, sol_VRP, node_shift);
    }
    
    
        
}


template <class T>
void Chromo_VRP<T>::add_lookup(){
    // Add lookup
    for(uint64_t i=0; i<this->length; i++){
        
        //std::cout << "I have value " << this->chromosome[i] << " going to " <<  i <<  "\n"; 
        //getchar();
        this->lookup[(uint64_t) this->chromosome[i]] = i;


    }
    /*
    for(uint64_t i=0; i<this->length; i++){
        std::cout << " " << (uint64_t) this->chromosome[i] << "@" << i;
    }
    std::cout<<"\n";
    getchar();
    */
}


template <class T>
void Chromo_VRP<T>::compute_fitness(void * solution_info){
    Sol_VRP_matrix * vrp = (Sol_VRP_matrix *) solution_info;
    long double path_sum = 0;
    //long double capacity_sum = 0;

    /*
struct Sol_VRP_matrix{
    long double ** dist;
    uint64_t n;
    uint64_t * demands; // Customer demands
    uint64_t depot; // Node depot
};
    */

    for(uint64_t i=1; i<this->length; i++){
        //std::cout << "Computing from " << thvoid Chromo_VRP<T>::verify_chromosome(char * step)is->chromosome[i-1] << " to " << this->chromosome[i] << " adds " << tsp->dist[this->chromosome[i-1]][this->chromosome[i]] << std::endl;
        path_sum +=  vrp->dist[this->chromosome[i-1]][this->chromosome[i]]; //Distance between node i and node j
        //capacity_sum += vrp->demands[this->chromosome[i]];
    }
    //std::cout << "Computing from " << this->chromosome[0] << " to " << this->chromosome[this->length-1] << " adds " << tsp->dist[this->chromosome[0]][this->chromosome[this->length-1]] << std::endl;
    path_sum += vrp->dist[this->depot][this->chromosome[0]]; //First vehicle depot 
    //capacity_sum += vrp->demands[this->chromosome[0]]; // First customer
    path_sum += vrp->dist[this->chromosome[this->length-1]][this->depot];
    //std::cout << "Sum of capacity : " << capacity_sum << " and distance travelled " << path_sum << std::endl;
    this->fitness = path_sum;    
}

template <class T>
void Chromo_VRP<T>::verify_chromosome(char * step){
    uint64_t verification[this->length];
    for(uint64_t i=0;i<this->length;i++){
        verification[i] = 0;
    }

    for(uint64_t i=0;i<this->length;i++){
        verification[this->chromosome[i]]++;
    }
    for(uint64_t i=0;i<this->length;i++){
        if(verification[i] != 1 && i != 0 && i < 150){
            this->print_chromosome();
            std::cout << "Found error at " << i << " at " << step << " having " << verification[i] << std::endl;
            throw "Aborting";
        }
    }
}


// Chromosome for Load balancing

template <class T>
Chromo_LB<T>::Chromo_LB(uint64_t alleles, Position p, uint64_t n_threads, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d){
    this->chromosome = (T *) std::malloc(alleles * sizeof(T));
    if(this->chromosome == NULL) throw "Could not allocate chromosome";
    this->length = alleles;
    this->fitness = LDBL_MAX;
    this->position = p;
    uint64_t i;
    for(i=0;i<alleles;i++){
        this->chromosome[i] = (unsigned char) (i % n_threads);
    }
    
    srand(time(NULL));
    random_shuffle_templated(this->length, this->chromosome, rand()%100, g, u_d);
    
}

template <class T>
void Chromo_LB<T>::compute_fitness(void * solution_info){
    Sol_LB_reads * lb = (Sol_LB_reads *) solution_info;
    long double fitness = 0;
    uint64_t trails[lb->threads];
    memset(trails, 0, lb->threads*sizeof(uint64_t));
    // Compute trails 
    for(uint64_t i=0;i<lb->n;i++){
        trails[(uint64_t)this->chromosome[i]] += lb->lengths[i];
    }
    uint64_t max_len = 0;
    // Find max trail
    for(uint64_t i=0;i<lb->threads;i++){
        if(max_len < trails[i]) max_len = trails[i];
    }
    for(uint64_t i=0;i<lb->threads;i++){
        fitness += max_len - trails[i];
    }
    // Minimize distance with the longest execution trace    
    this->fitness = fitness;
    
}


template class Chromosome<double>;
template class Chromosome<unsigned char>;
template class Chromosome<uint64_t>;
template class Chromo_rucksack<double>;
template class Chromo_subsetsum<unsigned char>;
template class Chromo_TSP<uint64_t>;
template class Chromo_VRP<uint64_t>;
template class Chromo_LB<unsigned char>;