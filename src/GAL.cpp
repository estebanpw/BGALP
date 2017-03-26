#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
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

int DEBUG_ACTIVE = 0;


/*
    TODO
    replace all random generators by interval [0,1] and just multiply by lenghts
*/


void init_args(int argc, char ** av, char * data, uint64_t * n_itera, uint64_t * n_individuals, uint64_t * part, uint64_t * mix_every);

int main(int argc, char **av) {

    // Number of iterations
    uint64_t n_itera = 10000;
    // Size of individuals
    uint64_t n_individuals = 100;
    // Partitions of population
    uint64_t part = 1;
    // Mix an individual between populations every k iterations
    uint64_t mix_every = n_itera/20;

    // TSP structure
    Sol_TSP_matrix tsp;
    // Path to tsp lib
    char tsp_lib[MAX_PATH];
    tsp_lib[0] = '\0';

    // Init arguments
    init_args(argc, av, tsp_lib, &n_itera, &n_individuals, &part, &mix_every);
    mix_every = (n_itera/mix_every != 0) ? (n_itera/mix_every) : (100);

    // Readstream to load data 
    Readstream * rs = new Readstream(tsp_lib, &reading_function_TSP, (void *) &tsp);
    rs->read();

    // Number of alleles per individual
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

    fprintf(stdout, "Best individual fitness: %.3Le\n", *manager->get_best_individual()->get_fitness());
    manager->get_best_individual()->print_chromosome();
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

    // Get k random solutions
    uint64_t random_sols = 10;
    
    std::cout << "Random solutions: " << std::endl;
    for(uint64_t i=0;i<random_sols;i++){
        ind[i].print_chromosome();
    }
    //getchar();
    
    //uint64_t combinations = ((random_sols)*(random_sols-1))/2;
    // Local search 
    std::cout << "Creating offsprings: " << std::endl;
    Chromo_TSP<uint64_t> * locally_optimals = (Chromo_TSP<uint64_t> *) std::malloc(random_sols*sizeof(Chromo_TSP<uint64_t>));
    if(locally_optimals == NULL) throw "Could not allocate locally optimal individuals";
    for(uint64_t i=0;i<random_sols;i++){
        new (&locally_optimals[i]) Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    }
    
    /*
    for(uint64_t i=0;i<random_sols;i++){
        for(uint64_t j=0;j<n_alleles;j++){
            locally_optimals->set_allele(j, ind[i].get_allele(j));
        }
    }
    */
    
    //Chromo_TSP<uint64_t> * aux1 = new Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    //Chromo_TSP<uint64_t> * aux2 = new Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    
    // For example 1 Whitley et al
    /*
    uint64_t values[12] = {0,1,2,3,4,5,6,7,8,9,10,11};
    aux1->set_allele(0, &values[0]);
    aux1->set_allele(1, &values[6]);
    aux1->set_allele(2, &values[8]);
    aux1->set_allele(3, &values[10]);
    aux1->set_allele(4, &values[11]);
    aux1->set_allele(5, &values[9]);
    aux1->set_allele(6, &values[7]);
    aux1->set_allele(7, &values[3]);
    aux1->set_allele(8, &values[5]);
    aux1->set_allele(9, &values[2]);
    aux1->set_allele(10, &values[1]);
    aux1->set_allele(11, &values[4]);

    aux2->set_allele(0, &values[0]);
    aux2->set_allele(1, &values[6]);
    aux2->set_allele(2, &values[10]);
    aux2->set_allele(3, &values[8]);
    aux2->set_allele(4, &values[9]);
    aux2->set_allele(5, &values[11]);
    aux2->set_allele(6, &values[7]);
    aux2->set_allele(7, &values[3]);
    aux2->set_allele(8, &values[2]);
    aux2->set_allele(9, &values[5]);
    aux2->set_allele(10, &values[4]);
    aux2->set_allele(11, &values[1]);

    aux1->print_chromosome();
    aux2->print_chromosome();
    */

    
    // Run 2-opt 
    std::cout << "Running 2-opt: \n\t";
    for(uint64_t i=0;i<random_sols;i++){
        run_2opt(&ind[i], &locally_optimals[i], (void *) &tsp);
        std::cout << i << std::endl;
        fflush(stdout);
    }
    std::cout << std::endl;

    std::cout << "Locally optimal: " << std::endl;
    for(uint64_t i=0;i<random_sols;i++){
        locally_optimals[i].print_chromosome();
    }
    
    // Allocate edge tables, FIFO queue and memory pool
    Edge_T<uint64_t> ** e_table = (Edge_T<uint64_t> **) std::calloc(n_alleles, sizeof(Edge_T<uint64_t> *));
    if(e_table == NULL) throw "Could not allocate edges table";
    memory_pool * mp = new memory_pool(POOL_SIZE);
    Chromo_TSP<uint64_t> * offspring_1 = new Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    Chromo_TSP<uint64_t> * offspring_2 = new Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    Chromo_TSP<uint64_t> * after_px_2opt = new Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    std::queue<uint64_t> FIFO_queue;

    // For all possible combinations 
    
    for(uint64_t i=0;i<random_sols;i++){
        for(uint64_t j=i+1;j<random_sols;j++){
            // Restart memory pool, FIFO queue and Edge table
            mp->full_reset();
            while(!FIFO_queue.empty()) FIFO_queue.pop();
            memset(e_table, 0x0, n_alleles*sizeof(Edge_T<uint64_t> *)); // Contents are handled by the pool

            // Fill edge table for two random solutions
            fill_edge_table(&locally_optimals[i], e_table, mp);
            fill_edge_table(&locally_optimals[j], e_table, mp);
            
            // Calculate degree
            generate_degree(n_alleles, e_table);

            // Locate partitions (attempt to find only two parts.)
            uint64_t node_id = get_highest_node_unpartitioned(n_alleles, e_table);
            find_connected_components(node_id, 0, e_table, &FIFO_queue);
            node_id = get_highest_node_unpartitioned(n_alleles, e_table);
            find_connected_components(node_id, 1, e_table, &FIFO_queue);
            
            // Find if partition crossover is feasible
            Quartet<Edge_T<uint64_t>> px;
            find_surrogate_edge_that_partitionates(n_alleles, e_table, &px);
            if(px._p1._e1 != NULL){
                apply_PX_chromosomes(n_alleles, e_table, &px, &locally_optimals[i], &locally_optimals[j], offspring_1, offspring_2);
                // Recompute fitness 
                offspring_1->compute_fitness((void *) &tsp);
                offspring_2->compute_fitness((void *) &tsp);

                std::cout << "---------------" << std::endl;
                locally_optimals[i].print_chromosome();
                locally_optimals[j].print_chromosome();
                std::cout << "\tgenerates" << std::endl;
                offspring_1->print_chromosome();
                run_2opt(offspring_1, after_px_2opt, (void *) &tsp);
                //std::cout << "\tRespect to 2opt:\t" << *after_px_2opt->get_fitness() << std::endl;
                offspring_2->print_chromosome();

                std::cout << "#\t" << *locally_optimals[i].get_fitness() << "\t";
                std::cout << *locally_optimals[j].get_fitness() << "\t";
                std::cout << *offspring_1->get_fitness() << "\t" << *offspring_2->get_fitness() << "\t";
                std::cout << *after_px_2opt->get_fitness() << "\t";

                run_2opt(offspring_2, after_px_2opt, (void *) &tsp);

                std::cout << *after_px_2opt->get_fitness() << std::endl;
                //std::cout << "\tRespect to 2opt:\t" << *after_px_2opt->get_fitness() << std::endl;

            }
        }
    }
    
    /*
    std::cout << "Building edges table: " << std::endl;
    Edge_T<uint64_t> ** e_table = (Edge_T<uint64_t> **) std::calloc(n_alleles, sizeof(Edge_T<uint64_t> *));
    if(e_table == NULL) throw "Could not allocate edges table";
    memory_pool * mp = new memory_pool(POOL_SIZE);

    // Fill edge table for two random solutions
    fill_edge_table(aux1, e_table, mp);
    fill_edge_table(aux2, e_table, mp);
    // Calculate degree
    generate_degree(n_alleles, e_table);
    //print_edge_tables(n_alleles, e_table);

    // Locate partitions
    std::queue<uint64_t> FIFO_queue;
    uint64_t node_id = get_highest_node_unpartitioned(n_alleles, e_table);
    find_connected_components(node_id, 0, e_table, &FIFO_queue);
    //std::cout << "First partition " << std::endl;
    //print_edge_tables(n_alleles, e_table);

    node_id = get_highest_node_unpartitioned(n_alleles, e_table);
    find_connected_components(node_id, 1, e_table, &FIFO_queue);
    //std::cout << "Second partition " << std::endl;
    //print_edge_tables(n_alleles, e_table);

    // Find if partition crossover is feasible
    Quartet<Edge_T<uint64_t>> px;
    find_surrogate_edge_that_partitionates(n_alleles, e_table, &px);

    Chromo_TSP<uint64_t> * offspring_1 = new Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    Chromo_TSP<uint64_t> * offspring_2 = new Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);

    if(px._p1._e1 != NULL) apply_PX_chromosomes(n_alleles, e_table, &px, aux1, aux2, offspring_1, offspring_2);

    */
    // Deallocating
    /*
    std::cout << "Results from PX: " << std::endl;
    offspring_1->print_chromosome();
    offspring_2->print_chromosome();
    */ 
    
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

    return 0;
}


void init_args(int argc, char ** av, char * data, uint64_t * n_itera, uint64_t * n_individuals, uint64_t * part, uint64_t * mix_every){
    
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
        if(strcmp(av[pNum], "-mix") == 0){
            *mix_every = (uint64_t) atoi(av[pNum+1]);
        }
        /*
        if(strcmp(av[pNum], "-pathfiles") == 0){
            strncpy(path_files, av[pNum+1], strlen(av[pNum+1]));
            path_files[strlen(av[pNum+1])] = '\0';
        }
        if(strcmp(av[pNum], "-annotations") == 0){
            strncpy(path_annotations, av[pNum+1], strlen(av[pNum+1]));
            path_annotations[strlen(av[pNum+1])] = '\0';
        }
        if(strcmp(av[pNum], "-write_blocks_bps") == 0){
            char templine[READLINE]; templine[0] = '\0';
            strcpy(templine, av[pNum+1]);
            strcat(templine, ".blocks");
            *out_blocks = fopen64(templine, "wt");

            templine[0] = '\0';
            strcpy(templine, av[pNum+1]);
            strcat(templine, ".breakpoints");
            *out_breakpoints = fopen64(templine, "wt");
            if(*out_blocks==NULL || *out_breakpoints==NULL) terror("Could not open blocks/breakpoints file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *out_file = fopen64(av[pNum+1], "wt");
            if(out_file==NULL) terror("Could not open output file");
        }
        if(strcmp(av[pNum], "-min_len_trimming") == 0){
            *min_len_trimming = (uint64_t) atoi(av[pNum+1]);
            if(*min_len_trimming < 0) terror("Minimum trimming length must be zero or more");
        }
        if(strcmp(av[pNum], "-min_trim_itera") == 0){
            *min_trim_itera = (uint64_t) atoi(av[pNum+1]);
            if(*min_trim_itera < 0) terror("Minimum number of trimming iterations must be zero or more");
        }
        if(strcmp(av[pNum], "-hash_table_divisor") == 0){
            *ht_size = (uint64_t) atoi(av[pNum+1]);
            if(*ht_size < 1) terror("The hash table divisor must be one at least");
        }
        */
        pNum++;
    }
    /*
    if(*multifrags==NULL || *out_file==NULL || path_files[0] == '\0'){
        terror("A frags file, a path to the fasta files and an output file must be specified");
    }
    */
    if(data[0] == '\0') throw "No input data selected";
}

// For Subset sum 
/*
// Size of individuals and number of alleles per individual
    uint64_t n_individuals = 100;
    uint64_t n_alleles = 1000;
    // A generic position (0,0,0) for the chromosomes
    Position p = Position();

    // Solution to the problem
    Sol_subsetsum sss;
    sss.values = (int64_t *) std::malloc(n_alleles * sizeof(int64_t));
    if(sss.values == NULL) throw "Could not allocate solution";
    sss.c = 0;
    std::default_random_engine uniform_generator;
    std::uniform_int_distribution<uint64_t> u_d (0, THE_MAX);
    for(uint64_t i=0;i<n_alleles;i++){
        sss.values[i] = u_d(uniform_generator);
        if(sss.values[i] % 3 == 0) sss.c += sss.values[i];
    }
    for(uint64_t i=0;i<n_alleles;i++){
        fprintf(stdout, "%" PRId64", ", sss.values[i]);
    }
    fprintf(stdout, "\nC=%" PRId64"\n", sss.c);

    //Allocate chromosomes
    Chromo_subsetsum<unsigned char> * ind = (Chromo_subsetsum<unsigned char> *) std::malloc(n_individuals*sizeof(Chromo_subsetsum<unsigned char>));
    if(ind == NULL) throw "Could not allocate individuals";
    for(uint64_t i=0;i<n_individuals;i++){
        new (&ind[i]) Chromo_subsetsum<unsigned char>(n_alleles, p, RANDOM);

    }
    
    //Assign chromosomes to population
    Population<unsigned char> * population = new Population<unsigned char>(n_individuals, ind);
    population->set_neighborhood_function(&all_together);
    
    //Add manager
    Manager<unsigned char> * manager = new Manager<unsigned char>(1, &single_point_crossover, (void *) &sss, MINIMIZE);

    //Set population and put the manager to run
    manager->set_populations(population, 0);
    manager->run(10000000);

    fprintf(stdout, "Best individual fitness: %Le\n", *manager->get_best_individual()->get_fitness());



    return 0;
*/