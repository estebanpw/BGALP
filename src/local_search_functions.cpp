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

void _2optSwap_VRP(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, uint64_t i, uint64_t k, uint64_t stop){
        
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
    for(t=k+1;t<stop;t++){
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

                    // deactivate this for less complexity
                    goto start_again;
                }
            }
        }
        //improve++;
    }

    //printf("last: %Le %p\n", *route->get_fitness(), route);
    //printf("last: %Le %p\n", *two_opt_chrom->get_fitness(), two_opt_chrom);
}


void run_2opt_vrp(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, void * solution_info, uint64_t n_trucks){

    uint64_t i, k;
    Sol_VRP_matrix * vrp = (Sol_VRP_matrix *) solution_info;
    

    
    //uint64_t improve = 0;
    long double best_distance = LDBL_MAX;
    // To separate 2-OPT per chromosome 
    uint64_t stops[n_trucks];
    uint64_t index = 1;
    stops[0] = 1;
    for(i = 0; i< route->get_length(); i++){
        if(*route->get_allele(i) == 0){
            stops[index++] = i;
            //std::cout << "stopping at " << i << std::endl; getchar();
        }
    }


    for(uint64_t stop_idx = 0; stop_idx < index-1; stop_idx++){

        best_distance = LDBL_MAX;

        while(*route->get_fitness() < best_distance){
        //while(improve < 20){
            
            start_again:

            best_distance = *route->get_fitness();
            for(i = stops[stop_idx]; i < stops[stop_idx+1] - 2; i++){
                if(*route->get_allele(i) == 0) continue;

                for(k = i + 1; k < stops[stop_idx+1]-1; k++){
                    _2optSwap_VRP(route, two_opt_chrom, i, k, stops[stop_idx+1]-1);
                    two_opt_chrom->compute_fitness(vrp);
                    if (*two_opt_chrom->get_fitness() < best_distance){

                        //fprintf(stdout, "Improved at %" PRIu64", %" PRIu64" which is (%"PRIu64", %"PRIu64"\n", i, k, *route->get_allele(i), *route->get_allele(k));
                        

                        for(i=stops[stop_idx];i<stops[stop_idx+1]-1;i++){
                            route->set_allele(i, two_opt_chrom->get_allele(i));
                        }
                        route->set_fitness(*two_opt_chrom->get_fitness());
                        
                        route->verify_chromosome(" 2opt ");

                        //route->print_chromosome();
                        //getchar();
                        
                        //improve = 0;

                        // deactivate this for less complexity
                        goto start_again;
                    }
                }
            }
            //improve++;
        }
    }
    

    //printf("last: %Le %p\n", *route->get_fitness(), route);
    //printf("last: %Le %p\n", *two_opt_chrom->get_fitness(), two_opt_chrom);
}

int compare_neighbours(const void * p1, const void * p2){
    neighbour_type * n1 = (neighbour_type *) p1, * n2 = (neighbour_type *) p2;
    if(n1->distance < n2->distance) return -1;
    if(n1->distance == n2->distance) return 0;
    return 1;
}

void build_neighbours_matrix_and_DLB(uint64_t n_neighbours, void * sol_tsp, uint64_t to_keep_n){ // Build with n_nodes and use with less 
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix *) sol_tsp;
    tsp->neighbours = (neighbour_type **) std::malloc(n_neighbours * sizeof(neighbour_type *));
    tsp->DLB = (bool *) std::malloc(n_neighbours*sizeof(bool));
    if(tsp->neighbours == NULL || tsp->DLB == NULL) throw "Could not built neighbours matrix";
    uint64_t j;
    for(uint64_t i=0; i<n_neighbours; i++){
        tsp->DLB[i] = false;
        tsp->neighbours[i] = (neighbour_type *) std::malloc((n_neighbours) * sizeof(neighbour_type));
        if(tsp->neighbours[i] == NULL) throw "Could not allocate subloop of neighbours matrix";
        uint64_t index = 0;
        for(j=0; j<n_neighbours; j++){
            if(i != j){
                tsp->neighbours[i][index].id = j;
                tsp->neighbours[i][index].distance = tsp->dist[i][j];
                index++;
            } 
        }

        // Sort neighbours according to distance
        qsort(tsp->neighbours[i], n_neighbours-1, sizeof(neighbour_type), compare_neighbours);
        // Realloc 
        tsp->neighbours[i] = (neighbour_type *) std::realloc(tsp->neighbours[i], to_keep_n * sizeof(neighbour_type));
    }

    
    
}

/*
// Structs for TSP
struct Sol_TSP_matrix{
    long double ** dist;
    uint64_t n;
    neighbour_type ** neighbours; // For 2opt
    bool * DLB; // Dont Look Bits array for 2opt
};
*/

bool improve_city_2_opt(Chromosome<uint64_t> * tour, uint64_t base_pos, void * sol_tsp, uint64_t n_nodes, uint64_t n_neighbours, uint64_t * tour_positions){
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix *) sol_tsp;
    bool improved = false;
    uint64_t i, j, X1, X2, Y1, Y2; // Indexes and city numbers 
    for(uint64_t dir=0; dir<2; dir++){ //0 is forward, 1 is backward
        if(dir == 0){
            i = base_pos;
            X1 = *tour->get_allele(i);
            X2 = *tour->get_allele((i+1)%n_nodes);
        }else{
            i = (n_nodes + base_pos - 1) % n_nodes;
            X2 = *tour->get_allele(i);
            X1 = *tour->get_allele((i+1)%n_nodes);
        }

        for(uint64_t w=0;w<n_neighbours;w++){
            Y1 = tsp->neighbours[X1][w].id;
            //std::cout << " w is " << w << "\n";
            if(dir == 0){
                j = tour_positions[Y1];
                Y2 = *tour->get_allele((j+1)%n_nodes);
            }else{
                j = (n_nodes + tour_positions[Y1] - 1) % n_nodes;
                Y2 = *tour->get_allele(j);
            }
            if(X2 == Y1 || Y2 == X1) continue;
            if(X1 == Y1 && X2 == Y2) continue;

            // If there is gain
            if( (tsp->dist[X1][X2] + tsp->dist[Y1][Y2]) - (tsp->dist[X1][Y1] + tsp->dist[X2][Y2]) > 0){
                
                // i: startIndex, j: endIndex
                // Reverse segment 
                uint64_t left = i, right = j;
                uint64_t inversion_size = ((right - left + 1 + n_nodes) % n_nodes) / 2;
                
                left = (left + 1) % n_nodes;
                //right = (n_nodes + right - 1) % n_nodes;

                //std::cout << "Swapping (" << X1 << ", " << X2 << ") with ("<< Y1 << ", "<< Y2 << ") because gain is " << (tsp->dist[X1][X2] + tsp->dist[Y1][Y2]) - (tsp->dist[X1][Y1] + tsp->dist[X2][Y2]) << "\n";
                //std::cout << "in positions " << i << ", " << j << "\n";

                for(uint64_t counter=0; counter<inversion_size; counter++){
                    // Swap

                    uint64_t aux = *tour->get_allele(left);
                    tour->set_allele(left, tour->get_allele(right));

                    //std::cout << "swap is " << aux << " with " << *tour->get_allele(right) << "\n";

                    tour->set_allele(right, &aux);

                    tour_positions[*tour->get_allele(left)] = left;
                    tour_positions[*tour->get_allele(right)] = right;

                    left = (left + 1) % n_nodes;
                    right = (n_nodes + right - 1) % n_nodes;

                }

                //tour->print_chromosome();
                //getchar();
                // Set neighbours lists off 
                tsp->DLB[X1] = false; tsp->DLB[X2] = false; tsp->DLB[Y1] = false; tsp->DLB[Y2] = false;
                improved = true;
                goto out_2opt_loops;

            }
        }
    }
    goto otherwise_out;
    

    out_2opt_loops:
    //tour->compute_fitness((void *) tsp);
    //std::cout << "Improved!" << *tour->get_fitness() << std::endl;
    //getchar();
    otherwise_out:
    return improved;
}


void two_opt_DLB_NL(uint64_t n_nodes, void * sol_tsp, Chromosome<uint64_t> * route, uint64_t n_neighbours){
    bool locally_optimal = false;
    bool improved = false;
    bool any_change = false;

    uint64_t tour_positions[n_nodes];
    for(uint64_t i=0;i<n_nodes;i++){
        tour_positions[*route->get_allele(i)] = i;
    }
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix *) sol_tsp;
    while(!locally_optimal){
        any_change = false;
        for(uint64_t i=0; i<n_nodes; i++){
            
            if(tsp->DLB[i]) continue;
            improved = improve_city_2_opt(route, i, sol_tsp, n_nodes, n_neighbours, tour_positions);
            any_change = true;
            //std::cout << "f(" << i<<"): " << *route->get_fitness()<<"\n";
            if(!improved) tsp->DLB[i] = true; else locally_optimal = false;
        }
        if(any_change == false) break;
    }
    //std::cout << "------------------\n";
    route->compute_fitness(sol_tsp);
}


void * run_pthreads_two_opt(void * a){
    two_opt_args * args = (two_opt_args *) a;
    run_2opt(args->a, args->b, args->solution);
    return NULL;
}

void generate_petals_from_points_and_suboptimal_circular(uint64_t * c, void * sol_VRP, optimal_path<Chromo_TSP<uint64_t>> * best_paths){
    /*
    

    template <typename T, typename K>
    struct DPair{
        T _e1;
        K _e2;
    };

    struct Sol_VRP_matrix{
        CPair<long double> * points;
        long double ** dist;
        uint64_t n;
        uint64_t * demands; // Customer demands
        uint64_t depot; // Node depot
        long double capacity;
    };
    */
    
    
    Sol_VRP_matrix * vrp = (Sol_VRP_matrix *) sol_VRP;
    uint64_t * circular_buffer = (uint64_t *) std::calloc(vrp->n, sizeof(uint64_t));
    if(circular_buffer == NULL) throw "Could not allocate circular buffer";

    // node, angle
    unsigned char * used_alphas = (unsigned char *) std::calloc(vrp->n, sizeof(unsigned char)); // To tell which ones were used 
    DPair<uint64_t, long double> * alpha_sorted_table = (DPair<uint64_t, long double> *) std::malloc(vrp->n * sizeof(DPair<uint64_t, long double>));
    if(alpha_sorted_table == NULL) throw "Could not allocate table of sorted angles for nodes";

    // First, convert coordinates to make depot be the reference. Transform to polar coordinates. 
    long double x_prime, y_prime, alpha;
    for(uint64_t i=1; i<vrp->n; i++){
        // Transform respect to 0,0
        x_prime = vrp->points[i]._e1 - vrp->points[0]._e1;
        y_prime = vrp->points[i]._e2 - vrp->points[0]._e2;

        // Convert to polar (ignore radius)
        alpha = atan2(y_prime, x_prime);

        // Insert into node table
        alpha_sorted_table[i-1]._e1 = i;
        alpha_sorted_table[i-1]._e2 = alpha;
        #ifdef VERBOSE
        std::cout << i << ": " << alpha << std::endl;
        #endif
    }

    // Sort by angle and then just generate the solution from a random point in order while capacity constrains hold
    qsort((void *) &alpha_sorted_table[0], vrp->n - 1, sizeof(DPair<uint64_t, long double>), &compare_alpha_petals);


    /*
    struct Sol_VRP_matrix{
        CPair<long double> * points;
        long double ** dist;
        uint64_t n;
        uint64_t * demands; // Customer demands
        uint64_t depot; // Node depot
        long double capacity;
    };
    */

    // Fill the petals
    uint64_t allele_tracker = 0;
    uint64_t i_in_chromo = 0;
    long double current_cap = 0;
    uint64_t n_trucks = 1;
    
    
    while(allele_tracker < vrp->n - 1){
        
        /*
        optimal_path<Chromo_TSP<uint64_t>> * best_paths
template <typename T>
struct swath{
    T * origin; // points to the chromosome 
    uint64_t pos;
    uint64_t length;
    long double score;
};

template <typename T>
struct optimal_path{
    swath<T> ** nodes;
    uint64_t * indexes;
};
        */
        // conditions
        // 1 -> there is an optimal path
        if(best_paths->indexes[alpha_sorted_table[allele_tracker]._e1] > 0){
            uint64_t c_len = 0;
            uint64_t chromo_track = best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].pos;

            std::cout << *best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].origin->get_allele(chromo_track) << " <-> " << alpha_sorted_table[allele_tracker]._e1 << std::endl;

            while(c_len < best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].length && current_cap + vrp->demands[*best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].origin->get_allele(chromo_track)] < vrp->capacity){
                c[i_in_chromo] = *best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].origin->get_allele(chromo_track); // Add node                
                current_cap += vrp->demands[*best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].origin->get_allele(chromo_track)];
                std::cout << "Partially inserted " << c[i_in_chromo] << " and capacity is " << current_cap << " is it? " << best_paths->indexes[*best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].origin->get_allele(chromo_track)] << std::endl;
                used_alphas[*best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].origin->get_allele(chromo_track)] = 1; // used 
                i_in_chromo++;
                allele_tracker++;
                chromo_track++;
                c_len++;
            }

        } // deadlock with allele_tracker

        if(current_cap + vrp->demands[alpha_sorted_table[allele_tracker]._e1] > vrp->capacity){
            std::cout << "Reached cap: " << current_cap << " vs total " << vrp->capacity << std::endl;
            c[i_in_chromo] = vrp->depot;
            n_trucks++;
            current_cap = 0;
        }else if(used_alphas[allele_tracker] == 0){
            current_cap += vrp->demands[alpha_sorted_table[allele_tracker]._e1];
            c[i_in_chromo] = alpha_sorted_table[allele_tracker]._e1;
            std::cout << "Copied      " << c[i_in_chromo] << " and capacity is " << current_cap << std::endl;
            allele_tracker++;
        }else{
            allele_tracker++;
        }
        i_in_chromo++;
    }
    
    std::free(circular_buffer);
}

void generate_petals_from_points_and_suboptimal(uint64_t * c, void * sol_VRP, optimal_path<Chromo_TSP<uint64_t>> * best_paths){
    /*
    


    template <typename T, typename K>
    struct DPair{
        T _e1;
        K _e2;
    };
    */
    
    
    Sol_VRP_matrix * vrp = (Sol_VRP_matrix *) sol_VRP;
    
    // node, angle
    unsigned char * used_alphas = (unsigned char *) std::calloc(vrp->n, sizeof(unsigned char)); // To tell which ones were used 
    DPair<uint64_t, long double> * alpha_sorted_table = (DPair<uint64_t, long double> *) std::malloc(vrp->n * sizeof(DPair<uint64_t, long double>));
    if(alpha_sorted_table == NULL) throw "Could not allocate table of sorted angles for nodes";

    // First, convert coordinates to make depot be the reference. Transform to polar coordinates. 
    long double x_prime, y_prime, alpha;
    for(uint64_t i=1; i<vrp->n; i++){
        // Transform respect to 0,0
        x_prime = vrp->points[i]._e1 - vrp->points[0]._e1;
        y_prime = vrp->points[i]._e2 - vrp->points[0]._e2;

        // Convert to polar (ignore radius)
        alpha = atan2(y_prime, x_prime);

        // Insert into node table
        alpha_sorted_table[i-1]._e1 = i;
        alpha_sorted_table[i-1]._e2 = alpha;
        #ifdef VERBOSE
        std::cout << i << ": " << alpha << std::endl;
        #endif
    }

    // Sort by angle and then just generate the solution from a random point in order while capacity constrains hold
    qsort((void *) &alpha_sorted_table[0], vrp->n - 1, sizeof(DPair<uint64_t, long double>), &compare_alpha_petals);


    /*
    struct Sol_VRP_matrix{
        CPair<long double> * points;
        long double ** dist;
        uint64_t n;
        uint64_t * demands; // Customer demands
        uint64_t depot; // Node depot
        long double capacity;
    };
    */

    // Fill the petals
    uint64_t allele_tracker = 0;
    uint64_t i_in_chromo = 0;
    long double current_cap = 0;
    uint64_t n_trucks = 1;
    
    
    while(allele_tracker < vrp->n - 1){
        
        /*
        optimal_path<Chromo_TSP<uint64_t>> * best_paths
template <typename T>
struct swath{
    T * origin; // points to the chromosome 
    uint64_t pos;
    uint64_t length;
    long double score;
};

template <typename T>
struct optimal_path{
    swath<T> ** nodes;
    uint64_t * indexes;
};
        */
        // conditions
        // 1 -> there is an optimal path
        if(best_paths->indexes[alpha_sorted_table[allele_tracker]._e1] > 0){
            uint64_t c_len = 0;
            uint64_t chromo_track = best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].pos;
            uint64_t t_len = best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].length;
            Chromo_TSP<uint64_t> * cmt = best_paths->nodes[alpha_sorted_table[allele_tracker]._e1][0].origin;

            std::cout << *cmt->get_allele(chromo_track) << " <-> " << alpha_sorted_table[allele_tracker]._e1 << std::endl;

            while(c_len < t_len && (current_cap + vrp->demands[*cmt->get_allele(chromo_track)]) <= vrp->capacity && *cmt->get_allele(chromo_track) != vrp->depot && used_alphas[*cmt->get_allele(chromo_track)] == 0){
                c[i_in_chromo] = *cmt->get_allele(chromo_track); // Add node                
                current_cap += vrp->demands[*cmt->get_allele(chromo_track)];
                std::cout << "Partially inserted " << c[i_in_chromo] << " and capacity is " << current_cap << " is it? " << best_paths->indexes[*cmt->get_allele(chromo_track)] << std::endl;
                used_alphas[*cmt->get_allele(chromo_track)] = 1; // used 
                i_in_chromo++;
                //allele_tracker++;
                chromo_track++;
                c_len++;
            }
            if(c_len > 0) allele_tracker++;

        } // deadlock with allele_tracker

        if(current_cap + vrp->demands[alpha_sorted_table[allele_tracker]._e1] > vrp->capacity){
            std::cout << "Reached cap: " << current_cap << " vs total " << vrp->capacity << std::endl;
            c[i_in_chromo] = vrp->depot;
            n_trucks++;
            current_cap = 0;
            i_in_chromo++;
        }else if(used_alphas[alpha_sorted_table[allele_tracker]._e1] == 0){
            current_cap += vrp->demands[alpha_sorted_table[allele_tracker]._e1];
            used_alphas[alpha_sorted_table[allele_tracker]._e1] = 1;
            c[i_in_chromo] = alpha_sorted_table[allele_tracker]._e1;
            std::cout << "Copied      " << c[i_in_chromo] << " and capacity is " << current_cap << std::endl;
            allele_tracker++;
            i_in_chromo++;
        }else{
            allele_tracker++;
        }
        
    }
    
}

