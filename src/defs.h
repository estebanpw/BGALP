#pragma once


#define READBUF 1024*1024*1024 // 1 MB
#define POOL_SIZE 1024*1024*1024 //1 MB 
#define LIST_POOL_SIZE 128
#define MAX_PATH 2048
#define MAX_LINE 2048

#define UNCOMMON 0
#define COMMON 1
#define CIRCUIT_A 0
#define CIRCUIT_B 1
#define __STDC_FORMAT_MACROS
enum INITIALIZER { RANDOM, CLEAR, PETALS };

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
struct CPair{
    T _e1;
    T _e2;
};

template <typename T>
struct Quartet{
    Pair<T> _p1;
    Pair<T> _p2;
};

template <typename T, typename K>
struct DPair{
    T _e1;
    K _e2;
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

struct neighbour_type{
    uint64_t id;
    long double distance;
};

// Structs for TSP
struct Sol_TSP_matrix{
    long double ** dist;
    uint64_t n;
    neighbour_type ** neighbours; // For 2opt
    bool * DLB; // Dont Look Bits array for 2opt
};


struct Sol_VRP_matrix{
    CPair<long double> * points;
    long double ** dist;
    uint64_t n;
    uint64_t * demands; // Customer demands
    uint64_t depot; // Node depot
    long double capacity;
    neighbour_type ** neighbours; // For 2opt
    bool * DLB; // Dont Look Bits array for 2opt
};

// Struct for load balancing metagenomic reads 
struct Sol_LB_reads{
    uint64_t * lengths;
    uint64_t n;
    long double current_max;
    unsigned char threads;
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
    uint64_t n_breakpoints;
    bool already_tried_to_partitionate; // To mark if it was used to generate connected components
    bool already_surrogate;             // To mark if the edge was used in a surrogate 
    bool is_entry_cycle_A;              // To mark whether the edge is an entry of the connected component
    bool is_entry_cycle_B;
    bool is_exit_cycle_A;               // Same but with exit 
    bool is_exit_cycle_B;
    bool incoming_A;                      // Tells if the edge is incoming or outcoming (sequential order in TSP)
    bool incoming_B;
    uint64_t belongs_to_cycle;          // Holds the ID from the hamiltonian cycle that generated it (0 or 1 currently using only two cycles)
    uint64_t subsolution_in_A;          // The subsolution (aka truck) that owns this node in solution A 
    uint64_t subsolution_in_B;          // The subsolution (aka truck) that owns this node in solution B 
    int64_t connects_partition;         // If its a common edge, to which partition does it connect?
    struct Edge_T * out_node;
    uint64_t orig_pos_A;                // Original position in the chromosome (use for constant access to allele in chromsome)
    uint64_t orig_pos_B;                // Original position in the chromosome (use for constant access to allele in chromsome)
    bool breakpoint;                    // As if it was a common one 
};

template <typename T>
struct swath{
    T * origin; // points to the chromosome 
    uint64_t pos;
    uint64_t length;
    long double score;
    uint64_t verifier;
};

template <typename T>
struct optimal_path{
    swath<T> ** nodes;
    uint64_t * indexes;
};

template <typename T>
struct Surrogate_Edge_T{
    Edge_T<T> * left;
    Edge_T<T> * right;
};

template <typename T>
struct List{
    T v;
    struct List * next;
    struct List * prev;
};

template <typename T>
struct PXTable{
    uint64_t n_surrogate_edges;
    List<Surrogate_Edge_T<T>> * su_gates; 
};

template <typename T>
struct feasible_partition{
    Edge_T<T> * entry;
    Edge_T<T> * exit;
    bool reverse;
    long double score;
};

template <typename T>
struct Feasible{
    CPair<feasible_partition<T> **> feasible;
    uint64_t * n_entries;
};


