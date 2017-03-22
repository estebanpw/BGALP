#pragma once


#define READBUF 1000
#define POOL_SIZE 1024*1024*1024 //1 MB 
#define LIST_POOL_SIZE 128
#define MAX_PATH 2048
#define MAX_LINE 2048

#define UNCOMMON 0
#define COMMON 1

enum INITIALIZER { RANDOM, CLEAR };

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