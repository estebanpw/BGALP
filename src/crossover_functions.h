#pragma once

#include <queue>
#include <cstdlib>
#include "manager.h"
#include "chromosome.h"
#include "memory_pool.h"
#include "defs.h"
#define __STDC_FORMAT_MACROS




template <class T>
void single_point_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m);

template <class T>
void ordered_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m);

template <class T>
void fill_edge_table(Chromosome<T> * a, Edge_T<T> ** e_table, memory_pool * mp, uint64_t cycle_id);

template <class T>
void generate_degree(uint64_t n_nodes, Edge_T<T> ** e_table);

template <class T>
void mark_entries_and_exists(uint64_t n_nodes, Edge_T<T> ** e_table, std::queue<Edge_T<T> *> * entries_A, std::queue<Edge_T<T> *> * entries_B, std::queue<Edge_T<T> *> * exits_A, std::queue<Edge_T<T> *> * exits_B);

template <class T>
void mark_entries_and_exists_ghosted(uint64_t n_nodes, Edge_T<T> ** e_table, std::queue<Edge_T<T> *> * entries_A, std::queue<Edge_T<T> *> * entries_B, std::queue<Edge_T<T> *> * exits_A, std::queue<Edge_T<T> *> * exits_B);

template <class T>
Pair<Edge_T<T>> exit_from_entry(Edge_T<T> ** e_table, Edge_T<T> * entry, unsigned char CIRCUIT);

template <class T>
bool get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<T> ** e_table, uint64_t * node_id);

template <class T>
bool get_highest_node_unpartitioned_ghosted(uint64_t n_nodes, Edge_T<T> ** e_table, uint64_t * node_id);

template <class T>
void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<T> ** e_table, std::queue<T> * FIFO_queue);

template <class T>
Pair<Edge_T<T>> replace_surrogate_by_one(Edge_T<T> ** e_table, uint64_t i);

template <class T>
Pair<Edge_T<T>> abstract_replace_surrogate_by_one(Edge_T<T> ** e_table, uint64_t i);

template <class T>
Pair<Edge_T<T>> abstract_replace_surrogate_by_one_circuited(Edge_T<T> ** e_table, uint64_t i, uint64_t CIRCUIT, uint64_t * length);

template <class T>
Feasible<T> verify_entries_and_exits(uint64_t n_partitions, std::queue<Edge_T<T> *> * entries_A, std::queue<Edge_T<T> *> * entries_B, std::queue<Edge_T<T> *> * exits_A, std::queue<Edge_T<T> *> * exits_B, memory_pool * mp, Edge_T<T> ** e_table);

template <class T>
void add_ghost_vertices(uint64_t n_nodes, Edge_T<T> ** e_table, memory_pool * mp);

template <class T>
bool is_connected_to(Edge_T<T> ** e_table, uint64_t node_1, uint64_t node_2);

template <class T>
void find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * surrogates);

template <class T>
void generate_partitions(PXTable<T> * px_table, Edge_T<T> ** e_table, uint64_t n_nodes, memory_pool * mp);

template <class T>
void shorten_common_tours_ghosted(Edge_T<T> ** e_table, uint64_t n_nodes);

template <class T>
T evaluate_partition_subtours_multiple(Edge_T<T> * start, Edge_T<T> * end, bool reverse, Chromosome<T> * c, void * solution_info, Edge_T<T> ** e_table);

template <class T>
long double evaluate_partition_subtours_multiple_ghosted(Edge_T<T> * start, Edge_T<T> * end, bool reverse, Chromosome<T> * c, void * solution_info, Edge_T<T> ** e_table);

template <class T>
T evaluate_partition_subtours(Surrogate_Edge_T<T> * start, Surrogate_Edge_T<T> * end, Chromosome<T> * c, void * solution_info, int64_t partition1, int64_t partition2, Edge_T<T> ** e_table);

template <class T>
void apply_PX_chromosomes(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * px, Chromosome<T> * P1, Chromosome<T> * P2, Chromosome<T> * offspring_1, Chromosome<T> * offspring_2);

template <class T>
void apply_PX_chromosomes_best(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * px, Chromosome<T> * P1, Chromosome<T> * P2, Chromosome<T> * offspring_1);