#pragma once
#include <cfloat>
#include <cstdlib>
#include <cmath>
#include "chromosome.h"
#define __STDC_FORMAT_MACROS
//Struct for PTHREADs

struct two_opt_args{
    Chromosome<uint64_t> * a;
    Chromosome<uint64_t> * b;
    void * solution;
};

void _2optSwap(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, uint64_t i, uint64_t k);
void _2optSwap_VRP(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, uint64_t i, uint64_t k, uint64_t stop);
void run_2opt(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, void * solution_info);
void run_2opt_vrp(Chromosome<uint64_t> * route, Chromosome<uint64_t> * two_opt_chrom, void * solution_info, uint64_t n_trucks);
void build_neighbours_matrix_and_DLB(uint64_t n_neighbours, void * sol_tsp, uint64_t to_keep_n);
bool improve_city_2_opt(Chromosome<uint64_t> * tour, uint64_t base_pos, void * sol_tsp, uint64_t n_nodes, uint64_t n_neighbours, uint64_t * tour_positions);
void two_opt_DLB_NL(uint64_t n_nodes, void * sol_tsp, Chromosome<uint64_t> * route, uint64_t n_neighbours);
void * run_pthreads_two_opt(void * a);
void generate_petals_from_points_and_suboptimal(uint64_t * c, void * sol_VRP, optimal_path<Chromo_TSP<uint64_t>> * best_paths);
