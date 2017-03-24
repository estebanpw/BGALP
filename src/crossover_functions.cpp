#include "crossover_functions.h"

template <class T>
void single_point_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m){
    uint64_t midpoint1 = a->get_length()*m->u_d(m->uniform_generator);
    uint64_t midpoint2 = (a->get_length()-midpoint1)*m->u_d(m->uniform_generator);
    midpoint2 += midpoint1;
    uint64_t i;

    for(i=0;i<midpoint1;i++){
        replacement->set_allele(i, a->get_allele(i));
    }
    for(i=midpoint1;i<midpoint2;i++){
        replacement->set_allele(i, b->get_allele(i));
    }
    for(i=midpoint2;i<replacement->get_length();i++){
        replacement->set_allele(i, a->get_allele(i));
    }
}

template <class T>
void ordered_crossover(Chromosome<T> * a, Chromosome<T> * b, Chromosome<T> * replacement, Manager<T> * m){
    
    uint64_t low = (a->get_length())*m->u_d(m->uniform_generator);
    uint64_t high = low + (a->get_length()+1-low)*m->u_d(m->uniform_generator);
    
    for(uint64_t i=0;i<a->get_length();i++){
        m->marks[i] = 0;
    }
    
    /*
    if(high == a->get_length()){ std::cout << "it does reach high "<< high << std::endl;  getchar();}
    if(low == 0){ std::cout << "it does low " << low << std::endl;  getchar();}
    */

    for(uint64_t i=0;i<low;i++){
        replacement->set_allele(i, a->get_allele(i));
        m->marks[*a->get_allele(i)] = 1;
    }
    
    for(uint64_t i=high;i<a->get_length();i++){
        replacement->set_allele(i, a->get_allele(i));
        m->marks[*a->get_allele(i)] = 1;
    }

    
    
    uint64_t pos = low;
    //Copy the rest from the other
    
    for(uint64_t i=0;i<a->get_length();i++){
        if(m->marks[*b->get_allele(i)] == 0){
            //printf("pos: %" PRIu64"\n", pos);
            replacement->set_allele(pos, b->get_allele(i));
            m->marks[*b->get_allele(i)] = 1;
            pos++;
        }
    }
}

template <class T>
void fill_edge_table(Chromosome<T> * a, Edge_T<T> ** e_table, memory_pool * mp){
    
    uint64_t i;
    Edge_T<T> * previous;
    uint64_t current_allele, next_allele, previous_allele;
    // Add edges from circuit
    for(i=0;i<a->get_length();i++){
        
        current_allele = *a->get_allele(i);
        
        if(e_table[current_allele] == NULL){
            // Add node
            e_table[current_allele] = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
            e_table[current_allele]->partition = -1; // No partition
        }

        if(i > 0){ // The previous adjacent

            next_allele = *a->get_allele(i-1);
            Edge_T<T> * et_ptr = e_table[current_allele]->next;
            previous = NULL;
            while(et_ptr != NULL && et_ptr->node != next_allele){

                previous = et_ptr;
                et_ptr = et_ptr->next;
            }

            // Either we have a null node (vertex not added yet) or we found it from the other genome
            if(et_ptr == NULL){
                // Its a NULL one (vertex not added)
                et_ptr = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                et_ptr->node = next_allele;
                et_ptr->common = UNCOMMON;
                et_ptr->next = NULL;
                if(previous != e_table[current_allele] && previous != NULL) previous->next = et_ptr; else e_table[current_allele]->next = et_ptr;

            }else{
                // It exists, so its a common one
                et_ptr->common = COMMON;
            }


        }
        if(i < a->get_length()-1){ // The next adjacent

            previous_allele = *a->get_allele(i+1);
            Edge_T<T> * et_ptr = e_table[current_allele]->next;
            previous = NULL;
            while(et_ptr != NULL && et_ptr->node != previous_allele){

                previous = et_ptr;
                et_ptr = et_ptr->next;
            }

            // Either we have a null node (vertex not added yet) or we found it from the other genome
            if(et_ptr == NULL){
                // Its a NULL one (vertex not added)
                et_ptr = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                et_ptr->node = previous_allele;
                et_ptr->common = UNCOMMON;
                et_ptr->next = NULL;
                if(previous != e_table[current_allele] && previous != NULL) previous->next = et_ptr; else e_table[current_allele]->next = et_ptr;

            }else{
                // It exists, so its a common one
                et_ptr->common = COMMON;
            }
        }
    }
    // We have to add an edge between first and last 
    
    bool was_found = false;
    Edge_T<T> * et_ptr;
    et_ptr = e_table[*a->get_allele(0)]->next;
    while(et_ptr != NULL){
        if(et_ptr->node == *a->get_allele(a->get_length()-1)){ et_ptr->common = COMMON; was_found = true; break; }
        if(et_ptr->next == NULL) break;
        et_ptr = et_ptr->next;
    }
    // Add ending edge
    if(!was_found){
        et_ptr->next = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
        et_ptr->next->node = *a->get_allele(a->get_length()-1);
        et_ptr->next->common = UNCOMMON;
        et_ptr->next->next = NULL;
        printf("added (1) %" PRIu64"\n", et_ptr->next->node);
    }
    

    was_found = false;
    et_ptr = e_table[*a->get_allele(a->get_length()-1)]->next;
    while(et_ptr != NULL){
        if(et_ptr->node == *a->get_allele(0)){ et_ptr->common = COMMON; was_found = true; break; }
        if(et_ptr->next == NULL) break;
        et_ptr = et_ptr->next;
    }
    // Add ending edge
    if(!was_found){
        et_ptr->next = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
        et_ptr->next->node = *a->get_allele(0);
        et_ptr->next->common = UNCOMMON;
        et_ptr->next->next = NULL;
        printf("added (2) %" PRIu64"\n", et_ptr->next->node);
    }
    
    
}

template <class T>
void generate_degree(uint64_t n_nodes, Edge_T<T> ** e_table){
    Edge_T<T> * ptr = NULL;
    uint64_t degree;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL){
            degree = 0;
            ptr = e_table[i]->next;
            while(ptr != NULL){
                degree++;
                ptr = ptr->next;
            }
            e_table[i]->degree = degree;
        }
    }
}

template <class T>
uint64_t get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<T> ** e_table){
    
    uint64_t node_id = 0;
    uint64_t degree = 0;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL && e_table[i]->partition == -1){
            if(degree < e_table[i]->degree){
                degree = e_table[i]->degree;
                node_id = i;
            }
        }
    }
    return node_id;
}

template <class T>
void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<T> ** e_table, std::queue<T> * FIFO_queue){

    // Add first node
    Edge_T<T> * ptr = NULL;
    if(e_table[init_node] == NULL) return;


    ptr = e_table[init_node]->next;
    e_table[init_node]->partition = partition_label;
    while(ptr != NULL){
        if(e_table[ptr->node] != NULL && e_table[ptr->node]->partition == -1 && ptr->common == UNCOMMON){
            FIFO_queue->push(ptr->node);
            e_table[ptr->node]->partition = partition_label;
        }
        ptr = ptr->next;
    }

    while(!FIFO_queue->empty()){
        uint64_t target_node = FIFO_queue->front();
        FIFO_queue->pop();
        if(e_table[target_node] != NULL && e_table[target_node]->partition == -1){
            ptr = e_table[target_node]->next;
            while(ptr != NULL){
                if(e_table[ptr->node] != NULL && e_table[ptr->node]->partition == -1 && ptr->common == UNCOMMON){
                    FIFO_queue->push(ptr->node);
                    e_table[ptr->node]->partition = partition_label;
                }
                ptr = ptr->next;
            }
        }
    }
}

template <class T>
Pair<Edge_T<T>> replace_surrogate_by_one(Edge_T<T> ** e_table, uint64_t i){

    Edge_T<T> * ptr, * master_ptr, * previous;
    Pair<Edge_T<T>> surrogate; surrogate._e1 = NULL; surrogate._e2 = NULL;

    // Extend in both directions until one node with degree > 2 is found
    master_ptr = e_table[i];
    previous = master_ptr;
    ptr = master_ptr->next;
    while(master_ptr->degree == 2){
        printf("Can i get stuck here? %" PRIu64"\n", ptr->node);
        previous = master_ptr;
        // Traverse in one direction
        if(ptr->node != previous->node){
            master_ptr = e_table[ptr->node];
            ptr = master_ptr->next;
        }else{
            // In case the first one is the previous one
            master_ptr = e_table[ptr->next->node];
            ptr = master_ptr->next;
        }
    }
    
    // The surrogate should be in previous (since master_ptr is no longer of degree 2)
    surrogate._e1 = previous;
    // Same path but backwards
    master_ptr = e_table[i];
    previous = master_ptr;
    ptr = master_ptr->next->next; // This is the only difference!!
    while(master_ptr->degree == 2){
        
        previous = master_ptr;
        // Traverse in one direction
        if(ptr->node != previous->node){
            master_ptr = e_table[ptr->node];
            ptr = master_ptr->next;
        }else{
            // In case the first one is the previous one
            master_ptr = e_table[ptr->next->node];
            ptr = master_ptr->next;
        }
    }
    // The other end of the surrogate is in previous
    surrogate._e2 = previous;

    return surrogate;
}

template <class T>
Quartet<T> find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<T> ** e_table){

    Quartet<T> surrogates;
    Pair<Edge_T<T>> p1, p2;
    
    int64_t up_surrogate_1, down_surrogate_1, up_surrogate_2, down_surrogate_2;
    uint64_t j;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL){
            if(e_table[i]->degree == 2 && e_table[i]->partition == -1){
                // Its a surrogate edge 

                p1 = replace_surrogate_by_one(e_table, i);

                up_surrogate_1 = down_surrogate_1 = -1;

                if(e_table[p1._e1->next->node]->degree != 2 && e_table[p1._e1->next->node]->partition != - 1){
                    // Connects to a partition
                    up_surrogate_1 = p1._e1->next->node;
                }
                if(e_table[p1._e1->next->next->node]->degree != 2 && e_table[p1._e1->next->next->node]->partition != - 1){
                    // Connects to a partition 
                    up_surrogate_1 = p1._e1->next->node;
                }

                if(e_table[p1._e2->next->node]->degree != 2 && e_table[p1._e2->next->node]->partition != - 1){
                    // Connects to a partition
                    down_surrogate_1 = p1._e2->next->node;
                }
                if(e_table[p1._e2->next->next->node]->degree != 2 && e_table[p1._e2->next->next->node]->partition != - 1){
                    // Connects to a partition 
                    down_surrogate_1 = p1._e2->next->node;
                }
                // We have in up_surrogate_1 and down_surrogate_1 the node IDs that connect different partitions

                if(up_surrogate_1 != -1 && down_surrogate_1 != -1){
                    for(j=i+1;j<n_nodes;j++){
                        if(e_table[j] != NULL){
                            if(e_table[j]->degree == 2 && e_table[j]->partition == -1){
                                p2 = replace_surrogate_by_one(e_table, j);

                                // At this point we are building over p1, p2 all combinations of surrogated edges
                                // We should halt at first two surrogates that connect two partitions

                                up_surrogate_2 = down_surrogate_2 = -1;

                                if(e_table[p2._e1->next->node]->degree != 2 && e_table[p2._e1->next->node]->partition != - 1){
                                    // Connects to a partition
                                    up_surrogate_2 = p2._e1->next->node;
                                }
                                if(e_table[p2._e1->next->next->node]->degree != 2 && e_table[p2._e1->next->next->node]->partition != - 1){
                                    // Connects to a partition 
                                    up_surrogate_2 = p2._e1->next->node;
                                }

                                if(e_table[p2._e2->next->node]->degree != 2 && e_table[p2._e2->next->node]->partition != - 1){
                                    // Connects to a partition
                                    down_surrogate_2 = p2._e2->next->node;
                                }
                                if(e_table[p2._e2->next->next->node]->degree != 2 && e_table[p2._e2->next->next->node]->partition != - 1){
                                    // Connects to a partition 
                                    down_surrogate_2 = p2._e2->next->node;
                                }

                                if(up_surrogate_2 != -1 && down_surrogate_2 != -1){
                                    // Final check, are they connecting the same partitions?
                                    std::cout << "PX might be feasible: " << std::endl;
                                    std::cout << "SG(1) = " << up_surrogate_1 << ", " << down_surrogate_1 << std::endl;
                                    std::cout << "SG(2) = " << up_surrogate_2 << ", " << down_surrogate_2 << std::endl;
                                }
                                
                            }
                        }
                    }
                }
                                 
            }
        }
    }

    return surrogates;
}



template void single_point_crossover<unsigned char>(Chromosome<unsigned char> * a, Chromosome<unsigned char> * b, Chromosome<unsigned char> * replacement, Manager<unsigned char> * m);
template void ordered_crossover<uint64_t>(Chromosome<uint64_t> * a, Chromosome<uint64_t> * b, Chromosome<uint64_t> * replacement, Manager<uint64_t> * m);
template void fill_edge_table(Chromosome<uint64_t> * a, Edge_T<uint64_t> ** e_table, memory_pool * mp);
template void generate_degree(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template uint64_t get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<uint64_t> ** e_table, std::queue<uint64_t> * FIFO_queue);
template Quartet<uint64_t> find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template Pair<Edge_T<uint64_t>> replace_surrogate_by_one(Edge_T<uint64_t> ** e_table, uint64_t i);
//template Part_list<uint64_t> * generate_lists_from_G(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, memory_pool * mp);