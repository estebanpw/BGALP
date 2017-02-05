#pragma once

#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include "common_functions.h"
#include "defs.h"
#include "chromosome.h"


template <class T> class Population;

template <class T>
class Population{

protected:
    Chromosome<T> * individuals;
    uint64_t n_individuals;


public:
    Population(uint64_t n_individuals);
    void set_panmictic();
    ~Population();
    
};



