#include "local_search_functions.h"
#define __STDC_FORMAT_MACROS


void _2optSwap(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, uint64_t i, uint64_t k){
        
    // Take route 1 to i and add in order
    uint64_t t, backward_pos;
    for(t=0;t<i;t++){
        two_opt_chrom->set_allele(t, route->get_allele(t));
    }
    // Take route i to k and add reversed
    for(t=i;t<k+1;t++){
        backward_pos = k-(t-i);
        two_opt_chrom->set_allele(t, route->get_allele(backward_pos));
    }
    // Take route k to end and add in order
    for(t=k+1;t<route->get_length();t++){
        two_opt_chrom->set_allele(t, route->get_allele(t));
    }
}

void run_2opt(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, void * solution_info){

    uint64_t i, k;
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix *) solution_info;
    

    
    //uint64_t improve = 0;
    long double best_distance = LDBL_MAX;
    while(*route->get_fitness() < best_distance){
    //while(improve < 20){
        
        start_again:

        best_distance = *route->get_fitness();
        for(i = 0; i < route->get_length() - 1; i++){
            for(k = i + 1; k < route->get_length(); k++){
                _2optSwap(route, two_opt_chrom, i, k);
                two_opt_chrom->compute_fitness(tsp);
                if (*two_opt_chrom->get_fitness() < best_distance){
                    for(i=0;i<two_opt_chrom->get_length();i++){
                        route->set_allele(i, two_opt_chrom->get_allele(i));
                    }
                    route->set_fitness(*two_opt_chrom->get_fitness());
                    /*
                    fprintf(stdout, "Improved at %" PRIu64", %" PRIu64"\n", i, k);
                    route->print_chromosome();
                    getchar();
                    */
                    //improve = 0;
                    goto start_again;
                }
            }
        }
        //improve++;
    }

    //printf("last: %Le %p\n", *route->get_fitness(), route);
    //printf("last: %Le %p\n", *two_opt_chrom->get_fitness(), two_opt_chrom);
}

void * run_pthreads_two_opt(void * a){
    two_opt_args * args = (two_opt_args *) a;
    run_2opt(args->a, args->b, args->solution);
    return NULL;
}

