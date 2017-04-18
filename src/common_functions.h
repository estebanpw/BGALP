#pragma once
#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <random>
#include <iostream>
#include <math.h>
#include "chromosome.h"
#include "defs.h"
#include "memory_pool.h"
#define __STDC_FORMAT_MACROS
/**
 * Print the error message 's' and exit(-1)
 */
void terror(const char *s);

/**
 * Function to read char by char buffered from a FILE
 */
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

template <class T>
void random_shuffle_templated(uint64_t n_elements, T * vector, uint64_t seed, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d);

template <class T>
void restart_edge_tables(uint64_t n_nodes, Edge_T<T> ** e_table, memory_pool * mp);

template <class T>
void print_edge_tables(uint64_t n_nodes, Edge_T<T> ** e_table);

template <class T>
void print_edge_tables_ghosted(uint64_t n_nodes, Edge_T<T> ** e_table);

template <class T>
bool find_in_vector(std::vector<T> * v, T key);

template <class T>
uint64_t get_number_of_partitions(uint64_t n_nodes, Edge_T<T> ** e_table);

