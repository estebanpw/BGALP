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
Pair<Edge_T<T>> find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<T> ** e_table){

    Edge_T<T> * surrogate_ptr_1 = NULL, * surrogate_ptr_2 = NULL;
    Pair<Edge_T<T>> surrogates; surrogates._e1 = NULL; surrogates._e2 = NULL;
    int64_t surrogate_1_belongs_to, surrogate_2_belongs_to;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL){
            if(e_table[i]->degree == 2 && e_table[i]->partition == -1){
                // Its a surrogate edge 

                surrogate_1_belongs_to = e_table[e_table[i]->next->node]->partition;
                surrogate_2_belongs_to = e_table[e_table[i]->next->next->node]->partition;
                
                if(surrogate_1_belongs_to >= 0 && surrogate_2_belongs_to >= 0){
                    // If both are assigned to a partition 
                    if(surrogate_1_belongs_to != surrogate_2_belongs_to){
                        // They connect different partitions 

                        // Check if we already had one surrogate 
                        if(surrogate_ptr_1 == NULL){
                            surrogate_ptr_1 = e_table[i];
                        }
                        else if(surrogate_ptr_2 == NULL){
                            surrogate_ptr_2 = e_table[i];
                            // Now we already have two surrogate edges
                            surrogates._e1 = surrogate_ptr_1;
                            surrogates._e2 = surrogate_ptr_2;
                            return surrogates;
                        }
                    }
                }
            }
        }
    }

    #ifdef VERBOSE
    std::cout << "Found partitioning surrogate edges at " << surrogate_ptr_1->node << " and " << surrogate_ptr_2->node << std::endl;
    #endif
    return surrogates;
}


/*
template <class T>
Part_list<T> * generate_lists_from_G(uint64_t n_nodes, Edge_T<T> ** e_table, memory_pool * mp){

    // Traverse the edge table 
    // Add those edges with degree 2 at the head, and those with degree > 2 at the end
    Edge_T<T> * et_ptr;
    List<T> * LIST1 = NULL;
    List<T> * head_ptr = NULL, * last_ptr = NULL;
    Part_list<T> * part_list = (Part_list<T> *) mp->request_bytes(sizeof(Part_list<T>));
    part_list->const_access = (List<T> **) mp->request_bytes(n_nodes * sizeof(List<T> *)); // To access in constant time

    uint64_t degree;
    for(uint64_t i=0;i<n_nodes;i++){

        if(e_table[i] != NULL){
            degree = 0;
            et_ptr = e_table[i]->next;
            while(et_ptr != NULL){
                if(et_ptr->common != COMMON) degree++;
                et_ptr = et_ptr->next;
            }


            if(degree > 2){
                // Add at the end of the list
                if(last_ptr == NULL){
                    last_ptr = (List<T> *) mp->request_bytes(sizeof(List<T>));
                    last_ptr->v = e_table[i];
                    last_ptr->next = NULL;
                    last_ptr->prev = NULL;
                    LIST1 = last_ptr;
                    
                }else{
                    last_ptr->next = (List<T> *) mp->request_bytes(sizeof(List<T>));
                    last_ptr->next->v = e_table[i];
                    last_ptr->next->next = NULL;
                    last_ptr->next->prev = last_ptr;
                    last_ptr = last_ptr->next;
                    
                }
                part_list->last_ptr = last_ptr;
                last_ptr->v->degree = degree;
                part_list->const_access[i] = last_ptr;
            }else{
                // Add at the head
                head_ptr = (List<T> *) mp->request_bytes(sizeof(List<T>));
                head_ptr->v = e_table[i];
                head_ptr->next = LIST1;
                head_ptr->prev = NULL;
                if(LIST1 != NULL) LIST1->prev = head_ptr;
                LIST1 = head_ptr;
                // To have the first last reference
                if(last_ptr == NULL) last_ptr = head_ptr;
                head_ptr->v->degree = degree;
                part_list->const_access[i] = head_ptr;
            }
        }    
    }
    part_list->LIST = LIST1;
    return part_list;
}
*/

/*
template <class T>
Part_list<T> * locate_partitions_in_G(uint64_t n_nodes, Edge_T<T> ** e_table, Part_list<T> * part_list, memory_pool * mp){
    List<T> * C = NULL;
    List<T> * I = part_list->last_ptr; // Initialize i to last 
    uint64_t partitions = 0;

    // Find first with degree > 2
    C = part_list->LIST;
    while(C != NULL){
        if(C->v->degree > 2) break;
        C = C->next;
    }
    // Assign to first partition
    C->v->partition = partitions++;

    

    // Swap vertex C with LIST1(i)
    List<T> aux; List<T> * aux_ptr = C;
    aux.v = C->v;
    aux.next = C->next;
    aux.prev = C->prev;

    // Update pointers in const_access
    part_list[C->v->node] = I;
    part_list[I->v->node] = aux_ptr;

    C->v = I->v; C->next = I->next; C->prev = I->prev;
    I->v = aux.v; I->next = aux.next; I->prev = aux.prev;


    // Decrement i 
    I = I->prev;

    // Put connected vertices from C in FIFO
    List<T> * FIFO = NULL;
    List<T> * last_in_FIFO = FIFO;
    Edge_T<T> * add_connected = C->v->next;
    while(add_connected != NULL){
        
        if(add_connected->common == UNCOMMON && e_table[add_connected->node]->partition == -1){
            if(FIFO == NULL){
                // Its the fist one
                last_in_FIFO = (List<T> *) mp->request_bytes(sizeof(List<T>));
                last_in_FIFO->v = add_connected;
                last_in_FIFO->prev = NULL;
                last_in_FIFO->next = NULL;
            }else{
                last_in_FIFO->next = (List<T> *) mp->request_bytes(sizeof(List<T>));
                last_in_FIFO->next->v = add_connected;
                last_in_FIFO->next->prev = last_in_FIFO;
                last_in_FIFO->next->next = NULL;
                last_in_FIFO = last_in_FIFO->next;
            }
        }
        
        add_connected = add_connected->next;
    }

}
*/

template void single_point_crossover<unsigned char>(Chromosome<unsigned char> * a, Chromosome<unsigned char> * b, Chromosome<unsigned char> * replacement, Manager<unsigned char> * m);
template void ordered_crossover<uint64_t>(Chromosome<uint64_t> * a, Chromosome<uint64_t> * b, Chromosome<uint64_t> * replacement, Manager<uint64_t> * m);
template void fill_edge_table(Chromosome<uint64_t> * a, Edge_T<uint64_t> ** e_table, memory_pool * mp);
template void generate_degree(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template uint64_t get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<uint64_t> ** e_table, std::queue<uint64_t> * FIFO_queue);
template Pair<Edge_T<uint64_t>> find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
//template Part_list<uint64_t> * generate_lists_from_G(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, memory_pool * mp);