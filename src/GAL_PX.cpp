#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <cfloat>
#include <pthread.h>
#include <queue>
#include <time.h>
#include <random>
#include <iostream>
#include <float.h>
#include "manager.h"
#include "chromosome.h"
#include "distance_functions.h"
#include "crossover_functions.h"
#include "population.h"
#include "common_functions.h"
#include "mutation_functions.h"
#include "readstream.h"
#include "local_search_functions.h"

#define THE_MAX 1000
#define __STDC_FORMAT_MACROS
int DEBUG_ACTIVE = 0;


/*
    TODO
    replace all random generators by interval [0,1] and just multiply by lenghts
*/


void init_args(int argc, char ** av, char * data, uint64_t * n_itera, uint64_t * n_individuals, uint64_t * part, uint64_t * mix_every, uint64_t * n_trucks, uint64_t * random_sols, long double * capacity, uint64_t * node_shift);

int main(int argc, char **av) {

    // Number of iterations
    uint64_t n_itera = 10000;
    // Size of individuals
    uint64_t n_individuals = 100;
    // Partitions of population
    uint64_t part = 1;
    // Mix an individual between populations every k iterations
    uint64_t mix_every = n_itera/20;
    // Number of random solutions to recombine
    uint64_t random_sols = 10;

    // Number of trucks for the VRP
    uint64_t n_trucks = 30; // 1 is a TSP

    // Shifts in the petal algorithm 
    uint64_t node_shift = 0;

    // TSP structure
    Sol_VRP_matrix vrp;
    Sol_TSP_matrix tsp;
    // Path to tsp lib
    char tsp_lib[MAX_PATH];
    tsp_lib[0] = '\0';

    // Init arguments
    vrp.capacity = 0;
    init_args(argc, av, tsp_lib, &n_itera, &n_individuals, &part, &mix_every, &n_trucks, &random_sols, &vrp.capacity, &node_shift);
    long double temp_capacity = vrp.capacity;
    //mix_every = (n_itera/mix_every != 0) ? (n_itera/mix_every) : (100);

    // Readstream to load data 
    Readstream * rs = new Readstream(tsp_lib, &reading_function_VRP, (void *) &vrp);
    rs->read();
    tsp.dist = vrp.dist;
    tsp.n = vrp.n;
    // Override to adjust capacity
    if(temp_capacity != 0) vrp.capacity = temp_capacity;
    //vrp.capacity = 1000;



    // Number of alleles per individual
    uint64_t n_alleles_vrp = tsp.n + n_trucks - 1; // Considering *depot -> path -> depot -> path -> depot -> path -> *depot
    uint64_t n_alleles = tsp.n;
    long double best_fitness = DBL_MAX, average_fitness = DBL_MAX;
    uint64_t n_improved_local = 0, n_2opt_local = 0, strictly_improved = 0;
    
    // A generic position (0,0,0) for the chromosomes implies no geometry
    Position p = Position();


    
    // Random numbers for the manager
    uint64_t random_seed;
    std::default_random_engine generator = std::default_random_engine((uint64_t) &random_seed);
    std::uniform_int_distribution<uint64_t> u_d = std::uniform_int_distribution<uint64_t>(0, n_alleles_vrp-1);
    std::uniform_int_distribution<uint64_t> u_d_i = std::uniform_int_distribution<uint64_t>(0, n_individuals-1);
    std::uniform_real_distribution<double> u_r = std::uniform_real_distribution<double>(0, 1);

    // Allocate chromosomes
    
    Chromo_VRP<uint64_t> * ind = (Chromo_VRP<uint64_t> *) std::malloc(n_individuals*sizeof(Chromo_VRP<uint64_t>));
    Chromo_VRP<uint64_t> aux_vrp(n_alleles_vrp, n_trucks, vrp.capacity, vrp.depot, p, PETALS, &generator, &u_d, (void *) &vrp, 0);
    Chromo_VRP<uint64_t> aux_vrp_replacementA(n_alleles_vrp, n_trucks, vrp.capacity, vrp.depot, p, PETALS, &generator, &u_d, (void *) &vrp, 0);
    Chromo_VRP<uint64_t> aux_vrp_replacementB(n_alleles_vrp, n_trucks, vrp.capacity, vrp.depot, p, PETALS, &generator, &u_d, (void *) &vrp, 0);
    if(node_shift == 0) node_shift = max((uint64_t)1, n_alleles/n_individuals);
    uint64_t curr_shift = 0;
    if(ind == NULL) throw "Could not allocate individuals";
    for(uint64_t i=0;i<n_individuals;i++){
        //new (&ind_vrp[0]) Chromo_VRP<uint64_t>(n_alleles_vrp, n_trucks, vrp.capacity, vrp.depot, p, PETALS, &generator, &u_d, (void *) &vrp);
        new (&ind[i]) Chromo_VRP<uint64_t>(n_alleles_vrp, n_trucks, vrp.capacity, vrp.depot, p, PETALS, &generator, &u_d, (void *) &vrp, curr_shift);
        curr_shift = (curr_shift + node_shift) % (tsp.n - 1);
        ind[i].compute_fitness((void *) &vrp);

        local_swap_search(&ind[i], n_alleles_vrp, &generator, &u_r, (void *) &vrp);

        #ifdef VERBOSE
        ind[i].print_chromosome();
        //ind[i].verify_chromosome(" at init!!");
        #endif
        run_2opt_vrp(&ind[i], &aux_vrp, (void *) &vrp, n_trucks);
        #ifdef VERBOSE
        std::cout << "2-OPT improved ";
        ind[i].print_chromosome();
        #endif
        //getchar();
    }
        

    
    //exit(-1);

        

    // Get k random solutions
    
    // Allocate structure to hold best optimal paths generated in TSP
    
    // Pthreads 
    //two_opt_args args[random_sols];
    //pthread_t threads[random_sols];
    
    
    
    
    // Allocate edge tables, FIFO queue, table of partitions, and memory pool
    Edge_T<uint64_t> ** e_table = (Edge_T<uint64_t> **) std::calloc(2*n_alleles, sizeof(Edge_T<uint64_t> *));
    if(e_table == NULL) throw "Could not allocate edges table";
    // To have it sorted 
    Edge_T<uint64_t> ** e_table_sorted = (Edge_T<uint64_t> **) std::calloc(2*n_alleles, sizeof(Edge_T<uint64_t> *));
    if(e_table_sorted == NULL) throw "Could not allocate edges table";
    
    memory_pool * mp = new memory_pool(POOL_SIZE);
    
    

    std::queue<uint64_t> * FIFO_queue = new std::queue<uint64_t>();

    // Queues for multipartitioning
    std::queue<Edge_T<uint64_t> *> * entries_A, * entries_B, * exits_A, * exits_B;
    entries_A = new std::queue<Edge_T<uint64_t> *>();
    entries_B = new std::queue<Edge_T<uint64_t> *>();
    exits_A = new std::queue<Edge_T<uint64_t> *>();
    exits_B = new std::queue<Edge_T<uint64_t> *>();

    // For all possible combinations 
    
    
    //std::cout << best_fitness << "\t" << vrp.capacity << "\t" << node_shift << "\t" << n_improved_local << "\t" << feasible_local << "\n";
    
    unsigned char * chosen_entries = (unsigned char *) std::malloc(2*n_alleles*sizeof(unsigned char));
    long double * scores_A = (long double *) std::malloc(n_alleles * sizeof(long double));
    long double * scores_B = (long double *) std::malloc(n_alleles * sizeof(long double));
    if(chosen_entries == 0x0 || scores_A == 0x0 || scores_B == 0x0) throw "Could not allocate scores for A and B subtours";

    Chromo_VRP<uint64_t> * p1, * p2;
    Chromo_VRP<uint64_t> * best = &ind[0];
    Chromo_VRP<uint64_t> * worst = &ind[0];
    best_fitness = *ind[0].get_fitness();
    average_fitness = 0;
    
    uint64_t curr_pool = 0;
    uint64_t indiv_per_pool = n_individuals / part;
    double elitism_policy = 0.2;
    uint64_t local_optima_when_stagnant = 0;
    
    for(uint64_t i=0;i<n_itera;i++){

        // next population
        curr_pool = (curr_pool + 1) % part;


        // Get worst
        for(uint64_t t=indiv_per_pool*curr_pool;t<indiv_per_pool*(curr_pool+1);t++){
            if(*worst->get_fitness() < *ind[t].get_fitness()){
                worst = &ind[t];
            }
            average_fitness += *ind[t].get_fitness();
            if(best_fitness > *ind[t].get_fitness()){
                best_fitness = *ind[t].get_fitness();
                best = &ind[t];
            }
        }
        average_fitness /= (long double) indiv_per_pool;

        
        
        // Select tournament

        if(i % mix_every == 0){
            p1 = &ind[u_d_i(generator)];
            p2 = &ind[u_d_i(generator)];
            while(p1 == p2) p2 = &ind[u_d_i(generator)];
        }else{
            p1 = &ind[(uint64_t) (indiv_per_pool*curr_pool + (indiv_per_pool * u_r(generator)))];
            p2 = &ind[(uint64_t) (indiv_per_pool*curr_pool + (indiv_per_pool * u_r(generator)))];
            while(p1 == p2) p2 = &ind[(uint64_t) (indiv_per_pool*curr_pool + (indiv_per_pool * u_r(generator)))];
        }

        // Modify best policy
        if(u_r(generator) < elitism_policy){
            p1 = best;
            p2 = &ind[u_d_i(generator)];
            while(p1 == p2) p2 = &ind[u_d_i(generator)];
        }

        

        aux_vrp_replacementA.copy(p1); 
        aux_vrp_replacementB.copy(p2); 
    
        simple_mutation(&aux_vrp_replacementA, n_alleles_vrp, &generator, &u_r); // Mutate

        //local_swap_search(&aux_vrp_replacementA, n_alleles_vrp, &generator, &u_r, (void *) &vrp);
        //local_swap_search(&aux_vrp_replacementB, n_alleles_vrp, &generator, &u_r, (void *) &vrp);

        aux_vrp_replacementA.compute_fitness((void *) &vrp);
        aux_vrp_replacementB.compute_fitness((void *) &vrp);
    
        if(i > n_itera*0.99 && elitism_policy != 0){
            elitism_policy = 0;
            std::cout<<"#ELITISM OUT\n";
        }
            

        // Recombination part
        Chromo_VRP<uint64_t> * offspring;
        if(i > n_itera*0.99 || u_r(generator) < 0.0005){
            

            // Restart memory pool, FIFO queue and Edge table
            mp->full_reset();
            while(!FIFO_queue->empty()) FIFO_queue->pop();
            memset(e_table, 0x0, 2*n_alleles*sizeof(Edge_T<uint64_t> *)); // Contents are handled by the pool
            memset(&chosen_entries[0], (unsigned char) 2, 2*n_alleles); // Just a value to indicate that it has not been selected
            memset(&scores_A[0], (long double) 0, n_alleles);
            memset(&scores_B[0], (long double) 0, n_alleles);

            // Fill edge table for two random solutions
            aux_vrp_replacementA.add_lookup();
            aux_vrp_replacementB.add_lookup();

            fill_edge_table_vrp(&aux_vrp_replacementA, e_table, mp, CIRCUIT_A, vrp.depot);
            fill_edge_table_vrp(&aux_vrp_replacementB, e_table, mp, CIRCUIT_B, vrp.depot);
            
            
            // Calculate degree
            generate_degree(n_alleles, e_table);
            

            
            #ifdef VERBOSE
            print_edge_tables(n_alleles, e_table);
            #endif

            

            // Insert ghost vertices
            add_ghost_vertices_vrp(n_alleles, e_table, mp, vrp.depot);

            
            #ifdef VERBOSE
            //print_edge_tables_ghosted(n_alleles, e_table);
            //getchar();
            #endif

            


            // Locate partitions (attempt to find only two parts.)
            uint64_t node_id = 0, current_label = 0;
            bool keep_partitioning = true;
            e_table[vrp.depot]->already_tried_to_partitionate = true;
            //sort_edges_table_lookup(n_alleles, e_table, e_table_sorted);
            do{
                keep_partitioning = get_highest_node_unpartitioned_ghosted(n_alleles, e_table, &node_id);
                //keep_partitioning = get_highest_node_unpartitioned(n_alleles, e_table, &node_id);
                if(keep_partitioning){
                    if(true == find_connected_components_vrp(node_id, current_label, e_table, FIFO_queue, vrp.depot)) current_label++;
                    node_id++;
                }
            }while(keep_partitioning);

            
            
            // Get number of partitions 
            //uint64_t n_parts = get_number_of_partitions(n_alleles, e_table)+1;
            uint64_t n_parts = get_number_of_partitions_ghosted(n_alleles, e_table)+1;

            #ifdef VERBOSE
            std::cout << "Number of partitions: " << n_parts << std::endl;
            #endif
            // For good print
            //std::cout << recombinations++ << "\t" << n_parts << "\t";

            
            #ifdef VERBOSE
            print_edge_tables_ghosted(n_alleles, e_table);
            #endif
            
            

            // Generate surrogate edges for partitions 
            
            shorten_common_tours_ghosted(e_table, n_alleles); // TODO does this one need to be changed?

            #ifdef VERBOSE
            std::cout << "Shoretned tours \n";
            #endif
            

            // Mark entries and exists to tell if a multipartitioning approach is feasible
            // Empty queues first 
            while(!entries_A->empty()) entries_A->pop(); while(!entries_B->empty()) entries_B->pop();
            while(!exits_A->empty()) exits_A->pop(); while(!exits_B->empty()) exits_B->pop();

            

            mark_entries_and_exists_ghosted_vrp(n_alleles, e_table, entries_A, entries_B, exits_A, exits_B, vrp.depot);


            
            offspring = (*aux_vrp_replacementA.get_fitness() <= *aux_vrp_replacementB.get_fitness()) ? (&aux_vrp_replacementA) : (&aux_vrp_replacementB);
            offspring->compute_fitness( (void *) &vrp);
            long double prev_fit = *offspring->get_fitness();
            long double next_fit;
            uint64_t t_recomb = 0;
            Feasible<uint64_t> feasibility_partitioning = verify_entries_and_exits_vrp(n_parts, entries_A, entries_B, exits_A, exits_B, mp, e_table, (void *) &tsp, n_alleles, &aux_vrp_replacementA, &aux_vrp_replacementB, vrp.depot, offspring, (void *) &vrp, &t_recomb);

        
            offspring->compute_fitness( (void *) &vrp);
            next_fit = *offspring->get_fitness();
            if(t_recomb > 0 && elitism_policy == 0) n_improved_local++;
            if(t_recomb > 0 && prev_fit > *offspring->get_fitness()) strictly_improved++;
            if(elitism_policy == 0 && t_recomb > 0 && prev_fit > *offspring->get_fitness()) local_optima_when_stagnant++;
            run_2opt_vrp(offspring, &aux_vrp, (void *) &vrp, n_trucks);

        }else{
            offspring = (*aux_vrp_replacementA.get_fitness() <= *aux_vrp_replacementB.get_fitness()) ? (&aux_vrp_replacementA) : (&aux_vrp_replacementB);
            run_2opt_vrp(offspring, &aux_vrp, (void *) &vrp, n_trucks);
        }
        
        

        

        
        /*
        if(t_recomb > 0 && prev_fit > *offspring->get_fitness()) n_2opt_local++;
        */

        // Always replace worst
        if(*offspring->get_fitness() < *worst->get_fitness() || elitism_policy == 0){
            worst->copy(offspring);
            worst->compute_fitness((void *) &vrp);

        }
        


        if(i % 100 == 0){
            std::cout << "@" << i << " B: " << best_fitness << " AVG: " << average_fitness << std::endl;
            if(i % 1000 == 0) best->print_vrp_chromosome((void *) &vrp);
            //getchar();
        }

            
    }

    std::cout<< "@N local optima when stagnant "<< local_optima_when_stagnant << "\n";

    best->print_vrp_chromosome((void *) &vrp);

    std::cout << vrp.capacity << "\t" << node_shift << "\t" << n_improved_local<< "\t" << strictly_improved << "\t"  << n_2opt_local << "\n";

    
    

    return 0;
}


void init_args(int argc, char ** av, char * data, uint64_t * n_itera, uint64_t * n_individuals, uint64_t * part, uint64_t * mix_every, uint64_t * n_trucks, uint64_t * random_sols, long double * capacity, uint64_t * node_shift){
    
    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--debug") == 0) DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           GAL\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -data       [Path to data]\n");
            fprintf(stdout, "           -itera      [Integer > 0] def: 10000\n");
            fprintf(stdout, "           -indiv      [Integer > 0] def: 100\n");
            fprintf(stdout, "           -part       [Integer > 0] def: 1\n");
            fprintf(stdout, "           -mix        [Integer > 0] def: 10000/20\n");
            fprintf(stdout, "           -trucks     [Integer > 0] def: 1\n");
            fprintf(stdout, "           -capacity   [Double > 0]  def: specified in file\n");
            fprintf(stdout, "           -shift      [Integer > 0] def: dynamic\n");
            fprintf(stdout, "           -rsols      [Integer > 0] def: 10\n");
            fprintf(stdout, "           --debug     Turns debug on\n");
            fprintf(stdout, "           --help      Shows the help for program usage\n");
            exit(1);
        }
        
        if(strcmp(av[pNum], "-data") == 0){
            strncpy(data, av[pNum+1], strlen(av[pNum+1]));
            data[strlen(av[pNum+1])] = '\0';
        }
        if(strcmp(av[pNum], "-itera") == 0){
            *n_itera = (uint64_t) atoi(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-indiv") == 0){
            *n_individuals = (uint64_t) atoi(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-part") == 0){
            *part = (uint64_t) atoi(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-shift") == 0){
            *node_shift = (uint64_t) atoi(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-capacity") == 0){
            *capacity = (long double) atof(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-rsols") == 0){
            *random_sols = (uint64_t) atoi(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-mix") == 0){
            *mix_every = (uint64_t) atoi(av[pNum+1]);
        }
        if(strcmp(av[pNum], "-trucks") == 0){
            *n_trucks = (uint64_t) atoi(av[pNum+1]);
        }

        pNum++;
    }
    /*
    if(*multifrags==NULL || *out_file==NULL || path_files[0] == '\0'){
        terror("A frags file, a path to the fasta files and an output file must be specified");
    }
    */
    if(data[0] == '\0') throw "No input data selected";
}

