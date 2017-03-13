#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <random>
#include <math.h>
#include "defs.h"

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