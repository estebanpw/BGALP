#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
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

#define THE_MAX 1000

int DEBUG_ACTIVE = 0;


/*
    TODO
    replace all random generators by interval [0,1] and just multiply by lenghts
*/


void init_args(int argc, char ** av, char * data, uint64_t * n_itera);

int main(int argc, char **av) {

    // Number of iterations
    uint64_t n_itera = 10000;

    // TSP structure
    Sol_TSP_matrix tsp;
    // Path to tsp lib
    char tsp_lib[MAX_PATH];
    tsp_lib[0] = '\0';

    // Init arguments
    init_args(argc, av, tsp_lib, &n_itera);

    // Readstream to load data 
    Readstream * rs = new Readstream(tsp_lib, &reading_function_TSP, (void *) &tsp);
    rs->read();

    // Size of individuals and number of alleles per individual
    uint64_t n_individuals = 100;
    uint64_t n_alleles = tsp.n;
    
    
    // A generic position (0,0,0) for the chromosomes implies no geometry
    Position p = Position();


    // Allocate chromosomes

    std::default_random_engine generator = std::default_random_engine(time(NULL));
    std::uniform_int_distribution<uint64_t> u_d = std::uniform_int_distribution<uint64_t>(0, n_alleles-1);

    Chromo_TSP<uint64_t> * ind = (Chromo_TSP<uint64_t> *) std::malloc(n_individuals*sizeof(Chromo_TSP<uint64_t>));
    if(ind == NULL) throw "Could not allocate individuals";
    for(uint64_t i=0;i<n_individuals;i++){
        new (&ind[i]) Chromo_TSP<uint64_t>(n_alleles, p, RANDOM, &generator, &u_d);
    }


    // Assign chromosomes to population
    Population<uint64_t> * population = new Population<uint64_t>(n_individuals, ind);
    population->set_neighborhood_function(&all_together);

    // Add manager
    Manager<uint64_t> * manager = new Manager<uint64_t>(1, &ordered_crossover, &mutation_function_TSP, (void *) &tsp, MINIMIZE);
    manager->generate_marks_for_ordered_crossover(n_alleles);


    //Set population and put the manager to run
    manager->set_populations(population, 0);
    manager->run(n_itera);

    fprintf(stdout, "Best individual fitness: %.3Le\n", *manager->get_best_individual()->get_fitness());

    
    delete manager;

    return 0;
}


void init_args(int argc, char ** av, char * data, uint64_t * n_itera){
    
    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--debug") == 0) DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           GAL\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -data       [Path to data]\n");
            fprintf(stdout, "           -itera      [Integer > 0] def: 10000\n");
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