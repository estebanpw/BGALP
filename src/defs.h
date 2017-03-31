#pragma once


#define READBUF 1000
#define POOL_SIZE 1024*1024*1024 //1 MB 
#define LIST_POOL_SIZE 128
#define MAX_PATH 2048
#define MAX_LINE 2048

#define UNCOMMON 0
#define COMMON 1
#define __STDC_FORMAT_MACROS
enum INITIALIZER { RANDOM, CLEAR };

// Double-evaluation free 
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

template <typename T>
struct Unit{
    T elem;
    struct Unit * next;
};

template <typename T>
struct Pair{
    T * _e1;
    T * _e2;
};

template <typename T>
struct Quartet{
    Pair<T> _p1;
    Pair<T> _p2;
};

struct Position{
    int64_t x;
    int64_t y;
    int64_t z;
    Position(int64_t x = 0, int64_t y = 0, int64_t z = 0): x(x), y(y), z(z) {}
};

// Structs for Subset sum
struct Sol_subsetsum{
    int64_t * values; //The array of values
    int64_t c; //The value to reach
};

// Structs for TSP
struct Sol_TSP_matrix{
    long double ** dist;
    uint64_t n;
};

// Struct for edge tables 
template <typename T>
struct Edge_T{
    uint64_t node;
    unsigned char common;
    struct Edge_T * next;
    int64_t partition;
    uint64_t degree;
    uint64_t n_commons;
    bool already_tried_to_partitionate;
};

template <typename T>
struct Surrogate_Edge_T{
    Edge_T<T> * left;
    Edge_T<T> * right;
};

template <typename T>
struct List{
    struct Edge_T<T> * v;
    struct List * next;
    struct List * prev;
};

template <typename T>
struct Part_list{
    struct List<T> * last_ptr;
    struct List<T> * LIST;
    struct List<T> ** const_access;
};

