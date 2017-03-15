#include "local_search_functions.h"



void 2optSwap(Chromosome<uint64_t> * route, Chromosome<uint64_t> * 2opt_chrom, uint64_t i, uint64_t k){
        
    // Take route 1 to i and add in order
    uint64_t t, backward_pos;
    for(t=0;t<i;t++){
        2opt_chrom->set_allele(t, route->get_allele(t));
    }
    // Take route i to k and add reversed
    for(t=i;t<k;t++){
        backward_pos = k-t;
        2opt_chrom->set_allele(t, route->get_allele(backward_pos));
    }
    // Take route k to end and add in order
    for(t=k;t<route->get_length();t++){
        2opt_chrom->set_allele(t, route->get_allele(t));
    }
    2opt_chrom->compute_fitness();

}

void run_2opt(Chromosome<uint64_t> * route, Chromosome<uint64_t> * 2opt_chrom, Chromosome<uint64_t> * aux){

    uint64_t i, k;
    long double prev_fitness = LDBL_MAX, best_distance;
    while(2opt_chrom->get_fitness() < prev_fitness){
        prev_fitness = 2opt_chrom->get_fitness();

        start_again:

        best_distance = 2opt_chrom->get_fitness();
        for (i = 0; i < route->get_length() - 1; i++) {
            for (k = i + 1; k < route->get_length(); k++) {
                new_route = 2optSwap(existing_route, i, k)
                new_distance = calculateTotalDistance(new_route)
                if (new_distance < best_distance) {
                    existing_route = new_route
                    goto start_again
                }
            }
        }
    }
}