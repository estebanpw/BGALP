#pragma once

#include <queue>
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
void fill_edge_table(Chromosome<T> * a, Edge_T<T> ** e_table, memory_pool * mp);

template <class T>
void generate_degree(uint64_t n_nodes, Edge_T<T> ** e_table);

template <class T>
bool get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<T> ** e_table, uint64_t * node_id);

template <class T>
void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<T> ** e_table, std::queue<T> * FIFO_queue);

template <class T>
Pair<Edge_T<T>> replace_surrogate_by_one(Edge_T<T> ** e_table, uint64_t i);

template <class T>
bool is_connected_to(Edge_T<T> ** e_table, uint64_t node_1, uint64_t node_2);

template <class T>
void find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * surrogates);

template <class T>
void generate_partitions(PXTable<T> * px_table, Edge_T<T> ** e_table, uint64_t n_nodes, memory_pool * mp);

template <class T>
void apply_PX_chromosomes(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * px, Chromosome<T> * P1, Chromosome<T> * P2, Chromosome<T> * offspring_1, Chromosome<T> * offspring_2);
