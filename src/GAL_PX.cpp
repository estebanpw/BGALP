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


void init_args(int argc, char ** av, char * data, uint64_t * n_itera, uint64_t * n_individuals, uint64_t * part, uint64_t * mix_every, uint64_t * n_trucks, uint64_t * random_sols);

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
    uint64_t n_trucks = 10; // 1 is a TSP

    // TSP structure
    Sol_VRP_matrix vrp;
    Sol_TSP_matrix tsp;
    // Path to tsp lib
    char tsp_lib[MAX_PATH];
    tsp_lib[0] = '\0';

    // Init arguments
    init_args(argc, av, tsp_lib, &n_itera, &n_individuals, &part, &mix_every, &n_trucks, &random_sols);
    mix_every = (n_itera/mix_every != 0) ? (n_itera/mix_every) : (100);

    // Readstream to load data 
    Readstream * rs = new Readstream(tsp_lib, &reading_function_VRP, (void *) &vrp);
    rs->read();
    tsp.dist = vrp.dist;
    tsp.n = vrp.n;
    // Disable this 
    vrp.capacity = 300;


    // Number of alleles per individual
    uint64_t n_alleles_vrp = tsp.n + n_trucks - 1; // Considering *depot -> path -> depot -> path -> depot -> path -> *depot
    uint64_t n_alleles = tsp.n;
    
    
    // A generic position (0,0,0) for the chromosomes implies no geometry
    Position p = Position();


    
    // Random numbers for the manager
    std::default_random_engine generator = std::default_random_engine(time(NULL));
    std::uniform_int_distribution<uint64_t> u_d = std::uniform_int_distribution<uint64_t>(0, n_alleles-1);

    // Allocate chromosomes
    Chromo_TSP<uint64_t> * ind = (Chromo_TSP<uint64_t> *) std::malloc(n_individuals*sizeof(Chromo_TSP<uint64_t>));
    if(ind == NULL) throw "Could not allocate individuals";
    for(uint64_t i=0;i<n_individuals;i++){
        new (&ind[i]) Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    }
    
    // VRP chromosome
    Chromo_VRP<uint64_t> * ind_vrp = (Chromo_VRP<uint64_t> *) std::malloc(2*sizeof(Chromo_VRP<uint64_t>));
    if(ind_vrp == NULL) throw "Could not allocate individual vrp";
    new (&ind_vrp[0]) Chromo_VRP<uint64_t>(n_alleles_vrp, n_trucks, vrp.capacity, vrp.depot, p, PETALS, &generator, &u_d, (void *) &vrp);
    new (&ind_vrp[1]) Chromo_VRP<uint64_t>(n_alleles_vrp, n_trucks, vrp.capacity, vrp.depot, p, PETALS, &generator, &u_d, (void *) &vrp);

    ind_vrp[0].compute_fitness((void *) &vrp);
    ind_vrp[0].print_chromosome();
    //exit(-1);


    // Add manager
    Manager<uint64_t> * manager = new Manager<uint64_t>(part, mix_every, &ordered_crossover, &mutation_function_TSP, (void *) &tsp, MINIMIZE);
    manager->generate_marks_for_ordered_crossover(n_alleles);



    // Partitionate population
    Population<uint64_t> ** population = (Population<uint64_t> **) std::malloc(part*sizeof(Population<uint64_t> *));
    if(population == NULL) throw "Could not allocate populations";
    uint64_t rate = n_individuals/part;
    for(uint64_t i=0;i<part;i++){
        // Assign chromosomes to population
        population[i] = new Population<uint64_t>(rate, &ind[i*rate]);
        population[i]->set_neighborhood_function(&all_together);

        // Add the populations to the manager
        manager->set_populations(population[i], i);
    }

    

    // Put the manager to run
    manager->run(n_itera);

    #ifdef VERBOSE
    std::cout << "Best individual fitness:" << *manager->get_best_individual()->get_fitness() << std::endl;
    
    manager->get_best_individual()->print_chromosome();
    #endif
    //getchar();

    // Get k best solutions
    /*
    uint64_t n_best_sols = 5; 
    Chromo_TSP<uint64_t> ** best_chromos = (Chromo_TSP<uint64_t> **) manager->retrieve_k_best_solutions(n_best_sols);
    for(uint64_t i=0;i<n_best_sols;i++){
        best_chromos[i]->print_chromosome();
    }
    getchar();
    */
    /*
        template <typename T>
    struct swath{
        Chromosome<T> * origin;
        uint64_t pos;
        uint64_t length;
    };

    template <typename T>
    struct optimal_path{
        swath<T> * nodes;
    };
    */
        

    // Get k random solutions
    
    uint64_t n_combs = (random_sols*(random_sols-1))/2;
    uint64_t n_neighbours = 20;

    // Allocate structure to hold best optimal paths generated in TSP
    optimal_path<Chromo_TSP<uint64_t>> best_paths;
    best_paths.nodes = (swath<Chromo_TSP<uint64_t>> **) std::malloc(n_alleles * sizeof(swath<Chromo_TSP<uint64_t>> *));
    best_paths.indexes = (uint64_t * ) std::calloc(n_alleles, sizeof(uint64_t));
    if(best_paths.indexes == NULL) throw("Could not allocate indexes for table of best paths");
    for(uint64_t i=0; i<n_alleles; i++){

        best_paths.nodes[i] = (swath<Chromo_TSP<uint64_t>> *) std::calloc(n_combs, sizeof(swath<Chromo_TSP<uint64_t>>));
        if(best_paths.nodes[i] == NULL) throw("Could not allocate table of best paths");
    }
    

    //build_neighbours_matrix_and_DLB(n_alleles, (void *) &tsp);
    //two_opt_DLB_NL(n_alleles, (void *) &tsp, &ind[0], n_neighbours);


    // Pthreads 
    //two_opt_args args[random_sols];
    //pthread_t threads[random_sols];
    
    #ifdef VERBOSE
    std::cout << "Random solutions: " << std::endl;
    #endif
    for(uint64_t i=0;i<random_sols;i++){
        #ifdef VERBOSE
        std::cout<< "F(1): " << *ind[i].get_fitness();
        #endif
        build_neighbours_matrix_and_DLB(n_alleles, (void *) &tsp, n_neighbours);
        two_opt_DLB_NL(n_alleles, (void *) &tsp, &ind[i], n_neighbours);
        for(uint64_t j=0;j<n_alleles;j++){
            std::free(tsp.neighbours[j]);
        }
        std::free(tsp.neighbours);
        #ifdef VERBOSE
        std::cout<<"F(2): " << *ind[i].get_fitness() << "\n";
        #endif
    }
    
    

    // Allocate edge tables, FIFO queue, table of partitions, and memory pool
    Edge_T<uint64_t> ** e_table = (Edge_T<uint64_t> **) std::calloc(2*n_alleles, sizeof(Edge_T<uint64_t> *));
    if(e_table == NULL) throw "Could not allocate edges table";
    PXTable<uint64_t> * part_table = (PXTable<uint64_t> *) std::malloc(n_alleles * sizeof(PXTable<uint64_t>));
    if(part_table == NULL) throw "Could not allocate partition table";
    memory_pool * mp = new memory_pool(POOL_SIZE);
    
    std::queue<uint64_t> * FIFO_queue = new std::queue<uint64_t>();

    // Queues for multipartitioning
    std::queue<Edge_T<uint64_t> *> * entries_A, * entries_B, * exits_A, * exits_B;
    entries_A = new std::queue<Edge_T<uint64_t> *>();
    entries_B = new std::queue<Edge_T<uint64_t> *>();
    exits_A = new std::queue<Edge_T<uint64_t> *>();
    exits_B = new std::queue<Edge_T<uint64_t> *>();

    // For all possible combinations 
    
    std::cout << "N\tPart(s)\tP1\tP2\tO\n";
    uint64_t recombinations = 0;
    unsigned char * chosen_entries = (unsigned char *) std::malloc(2*n_alleles*sizeof(unsigned char));
    long double * scores_A = (long double *) std::malloc(n_alleles * sizeof(long double));
    long double * scores_B = (long double *) std::malloc(n_alleles * sizeof(long double));
    if(chosen_entries == 0x0 || scores_A == 0x0 || scores_B == 0x0) throw "Could not allocate scores for A and B subtours";

    for(uint64_t i=0;i<random_sols;i++){
        for(uint64_t j=i+1;j<random_sols;j++){
            // Restart memory pool, FIFO queue and Edge table
            mp->full_reset();
            while(!FIFO_queue->empty()) FIFO_queue->pop();
            memset(e_table, 0x0, 2*n_alleles*sizeof(Edge_T<uint64_t> *)); // Contents are handled by the pool
            memset(&chosen_entries[0], (unsigned char) 2, 2*n_alleles); // Just a value to indicate that it has not been selected
            memset(&scores_A[0], (long double) 0, n_alleles);
            memset(&scores_B[0], (long double) 0, n_alleles);

            // Fill edge table for two random solutions
            fill_edge_table(&ind[i], e_table, mp, CIRCUIT_A);
            fill_edge_table(&ind[j], e_table, mp, CIRCUIT_B);
            
            // Calculate degree
            generate_degree(n_alleles, e_table);

            #ifdef VERBOSE
            //print_edge_tables(n_alleles, e_table);
            #endif

            // Insert ghost vertices
            add_ghost_vertices(n_alleles, e_table, mp);


            #ifdef VERBOSE
            //print_edge_tables_ghosted(n_alleles, e_table);
            //getchar();
            #endif

            


            // Locate partitions (attempt to find only two parts.)
            uint64_t node_id, current_label = 0;
            bool keep_partitioning = true;
            do{
                
                keep_partitioning = get_highest_node_unpartitioned_ghosted(n_alleles, e_table, &node_id);
                //keep_partitioning = get_highest_node_unpartitioned(n_alleles, e_table, &node_id);
                if(keep_partitioning){
                    find_connected_components(node_id, current_label, e_table, FIFO_queue);
                    current_label++;
                }
            }while(keep_partitioning);
            
            // Get number of partitions 
            //uint64_t n_parts = get_number_of_partitions(n_alleles, e_table)+1;
            uint64_t n_parts = get_number_of_partitions_ghosted(n_alleles, e_table)+1;

            //std::cout << "Number of partitions: " << n_parts << std::endl;
            std::cout << recombinations++ << "\t" << n_parts << "\t";

            #ifdef VERBOSE
            print_edge_tables_ghosted(n_alleles, e_table);
            #endif
            

            // Reset partition table 
            //memset(part_table, 0, n_alleles*sizeof(PXTable<uint64_t>));

            // Generate surrogate edges for partitions 
            
            shorten_common_tours_ghosted(e_table, n_alleles);

            

            // Mark entries and exists to tell if a multipartitioning approach is feasible
            // Empty queues first 
            while(!entries_A->empty()) entries_A->pop(); while(!entries_B->empty()) entries_B->pop();
            while(!exits_A->empty()) exits_A->pop(); while(!exits_B->empty()) exits_B->pop();

            

            mark_entries_and_exists_ghosted(n_alleles, e_table, entries_A, entries_B, exits_A, exits_B);

            
            #ifdef VERBOSE
            std::cout << "================================================================" << std::endl;
            std::cout << "Combining " << std::endl;
            ind[i].print_chromosome();
            std::cout << " with: " << std::endl;
            ind[j].print_chromosome();
            #endif

            

            // Verify that all entries conduct to the same exit 
            Feasible<uint64_t> feasibility_partitioning = verify_entries_and_exits(n_parts, entries_A, entries_B, exits_A, exits_B, mp, e_table, (void *) &tsp, n_alleles, &best_paths, &ind[i], &ind[j]);

            

            long double best_offspring;
            if(*ind[i].get_fitness() < *ind[j].get_fitness()){
                best_offspring = *ind[i].get_fitness();
            }else{
                best_offspring = *ind[j].get_fitness();
            }
            
            #ifdef VERBOSE
            std::cout << "Summary of partitioning\n";
            #endif

            Edge_T<uint64_t> * first_start = NULL;
            uint64_t feas_index = 0;
            long double t_score_A = *ind[i].get_fitness(), t_score_B = *ind[j].get_fitness();
            for(uint64_t w=0;w<n_parts;w++){
                if(feasibility_partitioning.feasible._e1[w] != 0x0){
                    long double score_part_A = 0.0;
                    long double score_part_B = 0.0;
                    #ifdef VERBOSE
                    std::cout << "Partition " << w << " has " << feasibility_partitioning.n_entries[w] << " entries and exits: \n";
                    #endif
                    for(uint64_t k=0;k<feasibility_partitioning.n_entries[w];k++){

                        if(first_start == NULL) first_start = feasibility_partitioning.feasible._e1[w][k].entry;
                        

                        #ifdef VERBOSE
                        std::cout << "(-> "<< feasibility_partitioning.feasible._e1[w][k].entry->node << " -> " << feasibility_partitioning.feasible._e1[w][k].exit->node << " )\n";
                        #endif
                        //score_part_A += evaluate_partition_subtours_multiple(feasibility_partitioning.feasible._e1[w][k].entry, feasibility_partitioning.feasible._e1[w][k].exit, feasibility_partitioning.feasible._e1[w][k].reverse, &ind[i], (void *) &tsp, e_table);
                        //score_part_A = (long double) evaluate_partition_subtours_multiple_ghosted(feasibility_partitioning.feasible._e1[w][k].entry, feasibility_partitioning.feasible._e1[w][k].exit, feasibility_partitioning.feasible._e1[w][k].reverse, &ind[i], (void *) &tsp, e_table);
                        score_part_A = feasibility_partitioning.feasible._e1[w][k].score;

                        //best_paths.nodes[i]
                        
                        #ifdef VERBOSE
                        std::cout << "(-> "<< feasibility_partitioning.feasible._e1[w][k].entry->node << " -> " << feasibility_partitioning.feasible._e1[w][k].exit->node << " )\n";
                        #endif
                        //score_part_B += evaluate_partition_subtours_multiple(feasibility_partitioning.feasible._e1[w][k].entry, feasibility_partitioning.feasible._e1[w][k].exit, feasibility_partitioning.feasible._e2[w][k].reverse, &ind[j], (void *) &tsp, e_table);
                        //score_part_B = (long double) evaluate_partition_subtours_multiple_ghosted(feasibility_partitioning.feasible._e1[w][k].entry, feasibility_partitioning.feasible._e1[w][k].exit, feasibility_partitioning.feasible._e2[w][k].reverse, &ind[j], (void *) &tsp, e_table);
                        score_part_B = feasibility_partitioning.feasible._e2[w][k].score;

                        #ifdef VERBOSE
                        std::cout << "Scores: " << score_part_A << " :: " << score_part_B << "\n";
                        #endif

                        if(score_part_A <= score_part_B){
                            chosen_entries[feasibility_partitioning.feasible._e1[w][k].entry->node] = CIRCUIT_A;
                        }else{
                            chosen_entries[feasibility_partitioning.feasible._e2[w][k].entry->node] = CIRCUIT_B;
                        }
                        
                        scores_A[feas_index] = score_part_A;
                        scores_B[feas_index] = score_part_B;
                        feas_index++;
                        

                        t_score_A -= score_part_A; // To calculate the residual graph
                        t_score_B -= score_part_B;
                    }
                    /*
                    for(uint64_t k=0;k<feasibility_partitioning.n_entries[w];k++){
                        std::cout << "(-> "<< feasibility_partitioning.feasible._e1[w][k].entry->node << " -> " << feasibility_partitioning.feasible._e1[w][k].exit->node << " )\n";
                        //score_part_B += evaluate_partition_subtours_multiple(feasibility_partitioning.feasible._e1[w][k].entry, feasibility_partitioning.feasible._e1[w][k].exit, feasibility_partitioning.feasible._e2[w][k].reverse, &ind[j], (void *) &tsp, e_table);
                        score_part_B += evaluate_partition_subtours_multiple_ghosted(feasibility_partitioning.feasible._e1[w][k].entry, feasibility_partitioning.feasible._e1[w][k].exit, feasibility_partitioning.feasible._e2[w][k].reverse, &ind[j], (void *) &tsp, e_table);
                    }

                    std::cout << "Scores:\n(A): " << score_part_A << "\n(B): " << score_part_B << std::endl;
                    */
                }
            }

            
            if(t_score_A <= t_score_B){
                best_offspring = t_score_A;
            }else{
                best_offspring = t_score_B;
            }    
            for(uint64_t w=0; w<feas_index; w++){
                if(scores_A[w] <= scores_B[w]) best_offspring += scores_A[w]; else best_offspring += scores_B[w];
            }
            
            
            if(feas_index > 0){
                std::cout << *ind[i].get_fitness() << "\t" << *ind[j].get_fitness() << "\t" << best_offspring << "\t" << feas_index;
                if(best_offspring < *ind[i].get_fitness() && best_offspring < *ind[j].get_fitness()) std::cout << "\t*"; 


                // Re-check that they are local optima under 2-opt
                //build_neighbours_matrix_and_DLB(n_alleles, (void *) &tsp);
                //two_opt_DLB_NL(n_alleles, (void *) &tsp, &ind[i], n_neighbours);

            }
            std::cout << "\n";
            //std::cout << "Fitnesses:\n(A): " << *ind[i].get_fitness() << "\n(B): " << *ind[j].get_fitness() << "\n(O): " << best_offspring << std::endl;
            
            #ifdef VERBOSE
            getchar();
            #endif
            
            continue;
            



            
            // To hold pairs of surrogates
            //Quartet<Edge_T<uint64_t>> current_px;


            
            
        }
    }
    
    /*
    T * origin;
    uint64_t pos;
    uint64_t length;
    long double score;
    */

    // Sort suboptimal tours 
    for(uint64_t i=0; i<n_alleles; i++){
        qsort(&best_paths.nodes[i][0], best_paths.indexes[i], sizeof(swath<uint64_t>), compare_swaths);
    }
    // Print suboptimal tours
    for(uint64_t i=0; i<n_alleles; i++){
        std::cout << "At allele " << i << "---------" << std::endl;
        for(uint64_t j=0; j<best_paths.indexes[i]; j++){
            std::cout << "@" << best_paths.nodes[i][j].pos << " -> " << *best_paths.nodes[i][j].origin->get_allele(best_paths.nodes[i][j].pos) << " has " << best_paths.nodes[i][j].score << "::veri:" << best_paths.nodes[i][j].verifier << std::endl << "\t->";
            for(uint64_t k=0; k<best_paths.nodes[i][j].length; k++){
                std::cout << *best_paths.nodes[i][j].origin->get_allele(best_paths.nodes[i][j].pos+k) << ", ";
            }
            std::cout << "\n";
        }
    }

    // Generate vrp with suboptimal paths inside
    ind_vrp[0].compute_fitness((void *) &vrp);
    ind_vrp[0].print_chromosome();
    generate_petals_from_points_and_suboptimal( ind_vrp[0].get_chromosome(), (void *) &vrp, &best_paths);
    ind_vrp[0].compute_fitness((void *) &vrp);
    ind_vrp[0].print_chromosome();
    ind_vrp[0].verify_chromosome(" after petals ");
    //generate_petals_from_points(ind_vrp[0].get_chromosome(), (void *) &vrp);
    run_2opt_vrp(&ind_vrp[0], &ind_vrp[1], (void *) &vrp, n_trucks);
    ind_vrp[0].compute_fitness((void *) &vrp);
    ind_vrp[0].print_chromosome();



    
    
    for(uint64_t i=0;i<tsp.n;i++){
        std::free(tsp.dist[i]);
    }
    std::free(tsp.dist);


    delete manager;
    delete mp;
    delete rs;
    //delete offspring_1;
    //delete offspring_2;
    std::free(ind);
    std::free(population);
    std::free(e_table);
    std::free(part_table);
    std::free(scores_A);
    std::free(scores_B);
    std::free(chosen_entries);

    return 0;
}


void init_args(int argc, char ** av, char * data, uint64_t * n_itera, uint64_t * n_individuals, uint64_t * part, uint64_t * mix_every, uint64_t * n_trucks, uint64_t * random_sols){
    
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

