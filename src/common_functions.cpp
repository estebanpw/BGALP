#include "common_functions.h"
#define __STDC_FORMAT_MACROS

void terror(const char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

template <class T>
void random_shuffle_templated(uint64_t n_elements, T * vector, uint64_t seed, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d){
    
    T aux;
    uint64_t curr_num;
    for(uint64_t i=0;i<n_elements;i++){
        
        aux = vector[i];
        curr_num = u_d->operator()(*g);
        vector[i] = vector[curr_num];
        vector[curr_num] = aux;
    }
}



template <class T>
void restart_edge_tables(uint64_t n_nodes, Edge_T<T> ** e_table, memory_pool * mp){
    for(uint64_t i=0;i<n_nodes;i++){
        e_table[i] = NULL;
    }
    mp->full_reset();
}

template <class T>
void print_edge_tables(uint64_t n_nodes, Edge_T<T> ** e_table){
    Edge_T<T> * et_ptr;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL){

            std::cout << e_table[i]->partition << " -> (" << e_table[i]->degree << ") @" << i << ": ";
            et_ptr = e_table[i]->next;
            while(et_ptr != NULL){
                std::cout << et_ptr->node;
                if(et_ptr->common == COMMON) std::cout << "*";
                if(et_ptr->incoming_A) std::cout<<"+a";
                if(et_ptr->incoming_B) std::cout<<"+b";
                if((et_ptr->common == COMMON || et_ptr->belongs_to_cycle == CIRCUIT_A) && et_ptr->incoming_A == false) std::cout<<"-a";
                if((et_ptr->common == COMMON || et_ptr->belongs_to_cycle == CIRCUIT_B) && et_ptr->incoming_B == false) std::cout<<"-b";
                std::cout << ", ";
                et_ptr = et_ptr->next;
            }
            std::cout << std::endl;
        }
    }
    std::cout << "# ---------------" << std::endl;
}

template <class T>
void print_edge_tables_ghosted(uint64_t n_nodes, Edge_T<T> ** e_table){
    Edge_T<T> * et_ptr;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL){

            std::cout << e_table[i]->partition << " -> (" << e_table[i]->degree << ") @" << i << ": ";
            et_ptr = e_table[i]->next;
            while(et_ptr != NULL){
                if(et_ptr->node > n_nodes) std::cout << et_ptr->node - n_nodes << "'"; else std::cout << et_ptr->node;
                if(et_ptr->common == COMMON) std::cout << "*";
                if(et_ptr->incoming_A) std::cout<<"+a";
                if(et_ptr->incoming_B) std::cout<<"+b";
                if((et_ptr->common == COMMON || et_ptr->belongs_to_cycle == CIRCUIT_A) && et_ptr->incoming_A == false) std::cout<<"-a";
                if((et_ptr->common == COMMON || et_ptr->belongs_to_cycle == CIRCUIT_B) && et_ptr->incoming_B == false) std::cout<<"-b";
                std::cout << ", ";
                et_ptr = et_ptr->next;
            }
            std::cout << std::endl;
        }
    }

    for(uint64_t i=n_nodes;i<2*n_nodes;i++){
        if(e_table[i] != NULL){

            std::cout << e_table[i]->partition << "$$ -> (" << e_table[i]->degree << ") @" << i-n_nodes << ": ";
            et_ptr = e_table[i]->next;
            while(et_ptr != NULL){
                if(et_ptr->node > n_nodes) std::cout << et_ptr->node - n_nodes << "'"; else std::cout << et_ptr->node;
                if(et_ptr->common == COMMON) std::cout << "*";
                if(et_ptr->incoming_A) std::cout<<"+a";
                if(et_ptr->incoming_B) std::cout<<"+b";
                if((et_ptr->common == COMMON || et_ptr->belongs_to_cycle == CIRCUIT_A) && et_ptr->incoming_A == false) std::cout<<"-a";
                if((et_ptr->common == COMMON || et_ptr->belongs_to_cycle == CIRCUIT_B) && et_ptr->incoming_B == false) std::cout<<"-b";
                std::cout << ", ";
                et_ptr = et_ptr->next;
            }
            std::cout << std::endl;
        }
    }

    std::cout << "# ---------------" << std::endl;
}

template <class T>
uint64_t get_number_of_partitions(uint64_t n_nodes, Edge_T<T> ** e_table){
    int64_t high = -1;
    for(uint64_t i=0;i<n_nodes;i++){
        if(high < e_table[i]->partition) high = e_table[i]->partition;
    }
    if(high == -1) return 0;
    return (uint64_t) high;
}

template <class T>
uint64_t get_number_of_partitions_ghosted(uint64_t n_nodes, Edge_T<T> ** e_table){
    int64_t high = -1;
    for(uint64_t i=0;i<2*n_nodes;i++){
        if(e_table[i] != NULL){
            if(high < e_table[i]->partition) high = e_table[i]->partition;
        }
        
    }
    if(high == -1) return 0;
    return (uint64_t) high;
}

template <class T>
bool find_in_vector(std::vector<T> * v, T key){
    for(typename std::vector<T>::iterator it = v->begin() ; it != v->end(); ++it){
        if(*it == key) return true;
    }
    return false;
}

int compare_alpha_petals(const void * p1, const void * p2){
    if( ((DPair<uint64_t, long double> *) p1)->_e2 <= ((DPair<uint64_t, long double> *) p2)->_e2 ) return -1; else return 1;
}

int compare_edges_degree(const void * p1, const void * p2){
    if(((Edge_T<uint64_t> *)p1) == NULL) return 1;
    if(((Edge_T<uint64_t> *)p2) == NULL) return -1;
    if(((Edge_T<uint64_t> *)p1)->degree > ((Edge_T<uint64_t> *)p2)->degree) return -1;
    if(((Edge_T<uint64_t> *)p1)->degree == ((Edge_T<uint64_t> *)p2)->degree) return 0;
    return 1;
}

template <class T>
void sort_edges_table_lookup(uint64_t n_nodes, Edge_T<T> ** e_table, Edge_T<T> ** e_table_sorted){
    // Copy all references 
    memcpy(e_table_sorted, e_table, 2*n_nodes*sizeof(Edge_T<T> *));
    // Sort them 
    qsort((void *) &e_table_sorted[0], 2*n_nodes, sizeof(Edge_T<T> *), &compare_edges_degree);

    #ifdef VERBOSE

    for(uint64_t i=0; i<2*n_nodes; i++){
        if(e_table_sorted[i] != NULL) std::cout << "Node " << i << " -> " << e_table_sorted[i]->degree << "\n";
    }
    //getchar();

    #endif 

}

template <class T>
void generate_petals_from_points(T * c, void * sol_VRP, uint64_t node_shift){
    /*
    

    template <typename T, typename K>
    struct DPair{
        T _e1;
        K _e2;
    };
    */
    bool sarcastic_jump = (node_shift % 2 == 0) ? (true) : (false);
    bool jump = true;
    
    Sol_VRP_matrix * vrp = (Sol_VRP_matrix *) sol_VRP;
    
    // node, angle
    DPair<uint64_t, long double> * alpha_sorted_table = (DPair<uint64_t, long double> *) std::malloc(vrp->n * sizeof(DPair<uint64_t, long double>));
    if(alpha_sorted_table == NULL) throw "Could not allocate table of sorted angles for nodes";

    // First, convert coordinates to make depot be the reference. Transform to polar coordinates. 
    long double x_prime, y_prime, alpha;
    for(uint64_t i=1; i<vrp->n; i++){
        // Transform respect to 0,0
        x_prime = vrp->points[i]._e1 - vrp->points[0]._e1;
        y_prime = vrp->points[i]._e2 - vrp->points[0]._e2;

        // Convert to polar (ignore radius)
        alpha = atan2(y_prime, x_prime);

        // Insert into node table
        alpha_sorted_table[i-1]._e1 = i;
        alpha_sorted_table[i-1]._e2 = alpha;
        #ifdef VERBOSE
        //std::cout << i << ": " << alpha << std::endl;
        #endif
    }

    // Sort by angle and then just generate the solution from a random point in order while capacity constrains hold
    qsort((void *) &alpha_sorted_table[0], vrp->n - 1, sizeof(DPair<uint64_t, long double>), &compare_alpha_petals);

    
    /*
    struct Sol_VRP_matrix{
        CPair<long double> * points;
        long double ** dist;
        uint64_t n;
        uint64_t * demands; // Customer demands
        uint64_t depot; // Node depot
        long double capacity;
    };
    */

    // Fill the petals
    uint64_t allele_tracker = 0;
    uint64_t take_shifted_node = node_shift;
    uint64_t i_in_chromo = 0;
    long double current_cap = 0;
    uint64_t n_trucks = 1;
    uint64_t times_jumped = 1;
    uint64_t initial = node_shift;
    uint64_t last_added;
    
    while(allele_tracker < vrp->n - 1){
        
        if(current_cap + vrp->demands[alpha_sorted_table[take_shifted_node]._e1] > vrp->capacity){
            #ifdef VERBOSE
            //std::cout << "Reached cap: " << current_cap << " vs total " << vrp->capacity << std::endl;
            #endif
            c[i_in_chromo] = vrp->depot;
            n_trucks++;
            current_cap = 0;
        }else{

            current_cap += vrp->demands[alpha_sorted_table[take_shifted_node]._e1];
            c[i_in_chromo] = alpha_sorted_table[take_shifted_node]._e1;
            last_added = take_shifted_node;
            //std::cout << ", " << take_shifted_node;
            
            if(sarcastic_jump){
                if(jump){
                    take_shifted_node = (take_shifted_node + 2) % (vrp->n - 1);
                    times_jumped++;
                    
                    
                    if(initial == take_shifted_node){
                        if(take_shifted_node == 0) take_shifted_node = vrp->n - 2; else take_shifted_node = (take_shifted_node - 1) % (vrp->n - 1);
                    }
                    
                    
                }else{
                    if(take_shifted_node == 0) take_shifted_node = vrp->n - 2; else take_shifted_node = (take_shifted_node - 1) % (vrp->n - 1);
                    jump = true;
                    times_jumped = 0;
                }            
                
            }else{
                take_shifted_node = (take_shifted_node + 1) % (vrp->n - 1);
            }

            

            
            
            
            if(times_jumped == 2) jump = false;
            allele_tracker++;
            
        }
        i_in_chromo++;
    }
    /*
    if(last_added != take_shifted_node){
        // Have to add last
        
        current_cap += vrp->demands[alpha_sorted_table[take_shifted_node]._e1];
        c[i_in_chromo] = alpha_sorted_table[take_shifted_node]._e1;
        
    } 
    */
    
    
}







template void random_shuffle_templated<uint64_t>(uint64_t n_elements, uint64_t * vector, uint64_t seed, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d);
template void random_shuffle_templated<unsigned char>(uint64_t n_elements, unsigned char * vector, uint64_t seed, std::default_random_engine * g, std::uniform_int_distribution<uint64_t> * u_d);
template void restart_edge_tables(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, memory_pool * mp);
template void print_edge_tables(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template void print_edge_tables_ghosted(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template bool find_in_vector(std::vector<uint64_t> * v, uint64_t key);
template uint64_t get_number_of_partitions(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template uint64_t get_number_of_partitions_ghosted(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template void sort_edges_table_lookup(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, Edge_T<uint64_t> ** e_table_sorted);
template void generate_petals_from_points(uint64_t * c, void * sol_VRP, uint64_t node_shift);