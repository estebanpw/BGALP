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
            e_table[current_allele]->node = current_allele;
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
        //printf("added (1) %" PRIu64"\n", et_ptr->next->node);
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
        //printf("added (2) %" PRIu64"\n", et_ptr->next->node);
    }
    
    
}

template <class T>
void generate_degree(uint64_t n_nodes, Edge_T<T> ** e_table){
    Edge_T<T> * ptr = NULL;
    uint64_t degree;
    uint64_t n_commons;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL){
            degree = 0;
            n_commons = 0;
            ptr = e_table[i]->next;
            e_table[i]->already_tried_to_partitionate = false;
            while(ptr != NULL){
                degree++;
                if(ptr->common == COMMON) n_commons++;
                ptr = ptr->next;
            }
            e_table[i]->degree = degree;
            e_table[i]->n_commons = n_commons;
        }
    }
}

template <class T>
bool get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<T> ** e_table, uint64_t * node_id){
    
    *node_id = 0;
    bool found = false;
    uint64_t degree = 0;
    for(uint64_t i=0;i<n_nodes;i++){
        if(e_table[i] != NULL && e_table[i]->partition == -1){
            if(degree < e_table[i]->degree && e_table[i]->already_tried_to_partitionate == false){
                degree = e_table[i]->degree;
                *node_id = i;
                found = true;
                e_table[i]->already_tried_to_partitionate = true;
            }
        }
    }
    return found;
}

template <class T>
void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<T> ** e_table, std::queue<T> * FIFO_queue){

    // Add first node
    Edge_T<T> * ptr = NULL;
    if(e_table[init_node] == NULL) return;

    // Check if at least it has uncommon edges 
    ptr = e_table[init_node]->next;
    uint64_t count_uncommon = 0;
    while(ptr != NULL){
        if(ptr->common == UNCOMMON) count_uncommon++;
        ptr = ptr->next;
    }
    if(count_uncommon == 0) return;

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
        if(e_table[target_node] != NULL){
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

    Edge_T<T> * ptr, * last_replaced, * route_start, * route_end;
    Pair<Edge_T<T>> surrogate; surrogate._e1 = NULL; surrogate._e2 = NULL;
    
    // Extend in both directions until one node with degree > 2 is found
    uint64_t master_node = i;
    route_start = e_table[i];
    route_end = NULL;
    last_replaced = e_table[i];
    //printf("%" PRId64", ", e_table[last_replaced->node]->partition);
    ptr = e_table[i]->next;
    

    while(ptr != NULL){ // Until we find no more common edges
        
        // Loop until next common edge
        while(ptr != NULL && ptr->common != COMMON){
            ptr = ptr->next;
        }
        if(ptr != NULL && ptr->node != last_replaced->node){
            // Update previous so we wont traverse it again 
            
            //printf("%" PRId64", ", e_table[ptr->node]->partition);
            last_replaced = e_table[master_node];
            master_node = ptr->node;
            route_end = e_table[ptr->node];

            if(route_end->partition != -1 && route_start->partition != route_end->partition) goto finish;

            // Not null implies we found common edge 
            // Save this node as last and look for more 
            ptr = e_table[ptr->node]->next;
            
        }else if(ptr != NULL){
            // It is the same common edge as before, try to get next 
            ptr = ptr->next;
        }
    }

    finish:
    //printf("\n");
    surrogate._e1 = route_start;
    surrogate._e2 = route_end;

    
    return surrogate;
}

template <class T>
bool is_connected_to(Edge_T<T> ** e_table, uint64_t node_1, uint64_t node_2){
    Edge_T<T> * ptr = e_table[node_1]->next;
    while(ptr != NULL){
        
        if(ptr->node == node_2) return true;
        ptr = ptr->next;
    }
    return false;
}

template <class T>
void find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * surrogates){

    Pair<Edge_T<T>> p1, p2, p3;
    uint64_t j,k;
    for(uint64_t i=0;i<n_nodes;i++){

        // Only find surrogate edges from one vertex that is connected to "something"
        if(e_table[i]->n_commons == 1){
            p1 = replace_surrogate_by_one(e_table, i);
            //printf("This previous\n");

            if(p1._e1 != NULL && p1._e2 != NULL){

                if(p1._e1->partition != -1 && p1._e2->partition != -1 && p1._e1->partition != p1._e2->partition){
                    for(j=i+1;j<n_nodes;j++){
                        if(i != j){

                            if(e_table[j]->n_commons == 1){
                                p2 = replace_surrogate_by_one(e_table, j);
                                //printf("With this one\n");

                                if(p2._e1 != NULL && p2._e2 != NULL){
                                    if(p2._e1->partition != -1 && p2._e2->partition != -1 && p2._e1->partition != p2._e2->partition){
                                        if(p1._e1->node != p2._e1->node && p1._e1->node != p2._e2->node){
                                            if((p1._e1->partition == p2._e1->partition && p1._e2->partition == p2._e2->partition) || (p1._e1->partition == p2._e2->partition && p1._e2->partition == p2._e1->partition)){
                                                
                                                bool adjacency_condition_1 = is_connected_to(e_table, p1._e1->node, p2._e1->node);
                                                bool adjacency_condition_2 = is_connected_to(e_table, p1._e1->node, p2._e2->node);
                                                bool adjacency_condition_3 = is_connected_to(e_table, p1._e2->node, p2._e1->node);
                                                bool adjacency_condition_4 = is_connected_to(e_table, p1._e2->node, p2._e2->node);
                                                
                                                if(!adjacency_condition_1 && !adjacency_condition_2 && !adjacency_condition_3 && !adjacency_condition_4){

                                                    // Now check that is not another surrogate edge that escapes connecting the two partitions 
                                                    for(k=j+1;k<n_nodes;k++){
                                                        if(e_table[k]->n_commons == 1){
                                                            p3 = replace_surrogate_by_one(e_table, k);
                                                            if(p3._e1 != NULL && p3._e2 != NULL){
                                                                if(p3._e1->partition != -1 && p3._e2->partition != -1 && p3._e1->partition != p3._e2->partition){
                                                                    if(p3._e1->node != p1._e1->node && p3._e1->node != p1._e2->node && p3._e1->node != p2._e1->node && p3._e1->node != p2._e2->node){
                                                                        if(p3._e2->node != p1._e1->node && p3._e2->node != p1._e2->node && p3._e2->node != p2._e1->node && p3._e2->node != p2._e2->node){
                                                                            // Only two partitions so it must connect the two (0 and 1)
                                                                            std::cout << "Cancelled: " << std::endl;
                                                                            std::cout << "SG(1) = " << p1._e1->node << "(" << p1._e1->partition << " ), " << p1._e2->node << "(" << p1._e2->partition << ")" << std::endl;
                                                                            std::cout << "SG(2) = " << p2._e1->node << "(" << p2._e1->partition << " ), " << p2._e2->node << "(" << p2._e2->partition << ")" << std::endl;
                                                                            std::cout << "Because of: " << std::endl;
                                                                            std::cout << "SG(3) = " << p3._e1->node << "(" << p3._e1->partition << " ), " << p3._e2->node << "(" << p3._e2->partition << ")" << std::endl;
                                                                            surrogates->_p1._e1 = NULL;
                                                                            return;
                                                                        }
                                                                    }
                                                                    
                                                                }
                                                            }
                                                        }
                                                    }

                                                    // Partition is feasible if they connect the same partitions 
                                                    std::cout << "PX is feasible: " << std::endl;
                                                    std::cout << "SG(1) = " << p1._e1->node << "(" << p1._e1->partition << " ), " << p1._e2->node << "(" << p1._e2->partition << ")" << std::endl;
                                                    std::cout << "SG(2) = " << p2._e1->node << "(" << p2._e1->partition << " ), " << p2._e2->node << "(" << p2._e2->partition << ")" << std::endl;
                                                    surrogates->_p1 = p1;
                                                    surrogates->_p2 = p2;
                                                    return;
                                                    //getchar();
                                                }
                                            }
                                        }
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    surrogates->_p1._e1 = NULL;
}

template <class T>
void apply_PX_chromosomes(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * px, Chromosome<T> * P1, Chromosome<T> * P2, Chromosome<T> * offspring_1, Chromosome<T> * offspring_2){

    uint64_t swath_start_P1 = 0, swath_start_P2 = 0, swath_end_P1 = 0, swath_end_P2;
    uint64_t i, j;
    

    // pair1(e1,e2), pair2(e1,e2)
    // One swath goes from pair1(e2) to pair2(e1)
    // And the other goe sfrom pair2(e2) to pair1(e1)

    for(i=0;i<n_nodes;i++){
        // Find first swath in P1 
        if(*P1->get_allele(i) == px->_p1._e1->node) swath_start_P1 = i;
        if(*P1->get_allele(i) == px->_p1._e2->node) swath_end_P1 = i;
        // Find the first swath in P2 
        if(*P2->get_allele(i) == px->_p1._e1->node) swath_start_P2 = i;
        if(*P2->get_allele(i) == px->_p1._e2->node) swath_end_P2 = i;
        
        // Copy P1 in both offsprings
        offspring_1->set_allele(i, P1->get_allele(i));
        offspring_2->set_allele(i, P2->get_allele(i));
    }

    // Check that between swath_start and swath_end there are no other nodes in between in P2
    bool go_backwards_P2 = false;
    for(i=1;i<n_nodes;i++){
        if(*P2->get_allele((swath_start_P2 + i) % n_nodes) == px->_p2._e1->node){ go_backwards_P2 = true; goto escape_2_O1; }
        if(*P2->get_allele((swath_start_P2 + i) % n_nodes) == px->_p2._e2->node){ go_backwards_P2 = true; goto escape_2_O1; }
        if(*P2->get_allele((swath_start_P2 + i) % n_nodes) == px->_p1._e2->node){ goto escape_2_O1; }
    }

    escape_2_O1:
    // Check that between swath_start and swath_end there are no other nodes in between in P2
    bool go_backwards_P1 = false;
    for(i=1;i<n_nodes;i++){
        if(*P1->get_allele((swath_start_P1 + i) % n_nodes) == px->_p2._e1->node){ go_backwards_P1 = true; goto escape_1_O1; }
        if(*P1->get_allele((swath_start_P1 + i) % n_nodes) == px->_p2._e2->node){ go_backwards_P1 = true; goto escape_1_O1; }
        if(*P1->get_allele((swath_start_P1 + i) % n_nodes) == px->_p1._e2->node){ goto escape_1_O1; }
    }

    escape_1_O1:

    // At this point we have the swath position starting and ending in both P1 and P2 
    // and we know if the swaths are to be copied from left to right or right to left 

    // Proceed to copy 
    i = swath_end_P1;
    j = swath_end_P2;

    do{

        // Advance
        if(go_backwards_P1){
            if(i>0) i--; else i=n_nodes-1;
        }else{
            i = (i + 1) % n_nodes;
        }
        if(go_backwards_P2){
            if(j>0) j--; else j=n_nodes-1;
        }else{
            j = (j + 1) % n_nodes;
        }

        // Copy swath from P2 in offspring 1
        offspring_1->set_allele(i, P2->get_allele(j));

    }while(*P2->get_allele(j) != px->_p2._e1->node && *P2->get_allele(j) != px->_p2._e2->node);
    
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REPEAT 
    for(i=0;i<n_nodes;i++){
        // Find first second in P1 
        if(*P1->get_allele(i) == px->_p2._e1->node) swath_start_P1 = i;
        if(*P1->get_allele(i) == px->_p2._e2->node) swath_end_P1 = i;
        // Find the second swath in P2 
        if(*P2->get_allele(i) == px->_p2._e1->node) swath_start_P2 = i;
        if(*P2->get_allele(i) == px->_p2._e2->node) swath_end_P2 = i;
        
    }

    // Check that between swath_start and swath_end there are no other nodes in between in P2
    go_backwards_P2 = false;
    for(i=1;i<n_nodes;i++){
        if(*P2->get_allele((swath_start_P2 + i) % n_nodes) == px->_p1._e1->node){ go_backwards_P2 = true; goto escape_2_O2; }
        if(*P2->get_allele((swath_start_P2 + i) % n_nodes) == px->_p1._e2->node){ go_backwards_P2 = true; goto escape_2_O2; }
        if(*P2->get_allele((swath_start_P2 + i) % n_nodes) == px->_p2._e2->node){ goto escape_2_O2; }
    }

    escape_2_O2:
    // Check that between swath_start and swath_end there are no other nodes in between in P2
    go_backwards_P1 = false;
    for(i=1;i<n_nodes;i++){
        if(*P1->get_allele((swath_start_P1 + i) % n_nodes) == px->_p1._e1->node){ go_backwards_P1 = true; goto escape_1_O2; }
        if(*P1->get_allele((swath_start_P1 + i) % n_nodes) == px->_p1._e2->node){ go_backwards_P1 = true; goto escape_1_O2; }
        if(*P1->get_allele((swath_start_P1 + i) % n_nodes) == px->_p2._e2->node){ goto escape_1_O2; }
    }

    escape_1_O2:

    // At this point we have the swath position starting and ending in both P1 and P2 
    // and we know if the swaths are to be copied from left to right or right to left 

    // Proceed to copy 
    i = swath_end_P1;
    j = swath_end_P2;

    do{

        // Advance
        if(go_backwards_P1){
            if(i>0) i--; else i=n_nodes-1;
        }else{
            i = (i + 1) % n_nodes;
        }
        if(go_backwards_P2){
            if(j>0) j--; else j=n_nodes-1;
        }else{
            j = (j + 1) % n_nodes;
        }

        // Copy swath from P2 in offspring 1
        offspring_2->set_allele(j, P1->get_allele(i));

    }while(*P1->get_allele(i) != px->_p1._e1->node && *P1->get_allele(i) != px->_p1._e2->node);
    
}



template void single_point_crossover<unsigned char>(Chromosome<unsigned char> * a, Chromosome<unsigned char> * b, Chromosome<unsigned char> * replacement, Manager<unsigned char> * m);
template void ordered_crossover<uint64_t>(Chromosome<uint64_t> * a, Chromosome<uint64_t> * b, Chromosome<uint64_t> * replacement, Manager<uint64_t> * m);
template void fill_edge_table(Chromosome<uint64_t> * a, Edge_T<uint64_t> ** e_table, memory_pool * mp);
template void generate_degree(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template bool get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, uint64_t * node_id);
template void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<uint64_t> ** e_table, std::queue<uint64_t> * FIFO_queue);
template bool is_connected_to(Edge_T<uint64_t> ** e_table, uint64_t node_1, uint64_t node_2);
template void find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, Quartet<Edge_T<uint64_t>> * surrogates);
template Pair<Edge_T<uint64_t>> replace_surrogate_by_one(Edge_T<uint64_t> ** e_table, uint64_t i);
template void apply_PX_chromosomes(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, Quartet<Edge_T<uint64_t>> * px, Chromosome<uint64_t> * P1, Chromosome<uint64_t> * P2, Chromosome<uint64_t> * offspring_1, Chromosome<uint64_t> * offspring_2);
//template Part_list<uint64_t> * generate_lists_from_G(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, memory_pool * mp);