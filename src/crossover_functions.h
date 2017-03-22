#pragma once

#include <queue>
#include "manager.h"
#include "chromosome.h"
#include "memory_pool.h"
#include "defs.h"





template <class T>
void single_point_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m);

template <class T>
void ordered_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m);

template <class T>
void fill_edge_table(Chromosome<T> * a, Edge_T<T> ** e_table, memory_pool * mp);

template <class T>
void generate_degree(uint64_t n_nodes, Edge_T<T> ** e_table);

template <class T>
uint64_t get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<T> ** e_table);

template <class T>
void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<T> ** e_table, std::queue<T> * FIFO_queue);

/*
template <class T>
Part_list<T> * generate_lists_from_G(uint64_t n_nodes, Edge_T<T> ** e_table, memory_pool * mp);
*/