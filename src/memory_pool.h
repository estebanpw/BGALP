#pragma once
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <inttypes.h>
#include <vector>
#include <string.h>
#include <ctype.h>
#define __STDC_FORMAT_MACROS
class memory_pool{

private:
    std::vector<char *> * mem_pool;
    std::vector<uint64_t> * base;
    uint64_t current_pool;
    uint64_t max_pools;
    uint64_t pool_size;

public:
    memory_pool(uint64_t pool_size);
    void * request_bytes(uint64_t bytes);
    void reset_n_bytes(uint64_t bytes);
    void reset_to(uint64_t pool, uint64_t position){ this->current_pool = 0; this->base->at(this->current_pool) = 0;}
    void full_reset();
    ~memory_pool();
};