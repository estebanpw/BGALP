#include "crossover_functions.h"
#define __STDC_FORMAT_MACROS

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
void fill_edge_table(Chromosome<T> * a, Edge_T<T> ** e_table, memory_pool * mp, uint64_t cycle_id){
    
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
                if(cycle_id == CIRCUIT_A) et_ptr->incoming_A = true; else et_ptr->incoming_B = true;
                et_ptr->belongs_to_cycle = cycle_id;
                et_ptr->next = NULL;
                if(previous != e_table[current_allele] && previous != NULL) previous->next = et_ptr; else e_table[current_allele]->next = et_ptr;

            }else{
                // It exists, so its a common one
                et_ptr->common = COMMON;
                if(cycle_id == CIRCUIT_A) et_ptr->incoming_A = true; else et_ptr->incoming_B = true;
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
                if(cycle_id == CIRCUIT_A) et_ptr->incoming_A = false; else et_ptr->incoming_B = false; // It goes out
                et_ptr->belongs_to_cycle = cycle_id;
                et_ptr->next = NULL;
                if(previous != e_table[current_allele] && previous != NULL) previous->next = et_ptr; else e_table[current_allele]->next = et_ptr;

            }else{
                // It exists, so its a common one
                et_ptr->common = COMMON;
                if(cycle_id == CIRCUIT_A) et_ptr->incoming_A = false; else et_ptr->incoming_B = false;
            }
        }
    }
    // We have to add an edge between first and last 
    
    bool was_found = false;
    Edge_T<T> * et_ptr;
    et_ptr = e_table[*a->get_allele(0)]->next;
    while(et_ptr != NULL){
        if(et_ptr->node == *a->get_allele(a->get_length()-1)){ 
            et_ptr->common = COMMON;
            if(cycle_id == CIRCUIT_A) et_ptr->incoming_A = true; else et_ptr->incoming_B = true;
            was_found = true; 
            break; 
        }
        if(et_ptr->next == NULL) break;
        et_ptr = et_ptr->next;
    }
    // Add ending edge
    if(!was_found){
        et_ptr->next = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
        et_ptr->next->node = *a->get_allele(a->get_length()-1);
        et_ptr->next->common = UNCOMMON;
        if(cycle_id == CIRCUIT_A) et_ptr->next->incoming_A = true; else et_ptr->next->incoming_B = true;
        et_ptr->next->belongs_to_cycle = cycle_id;
        et_ptr->next->next = NULL;
        //printf("added (1) %" PRIu64"\n", et_ptr->next->node);
    }
    

    was_found = false;
    et_ptr = e_table[*a->get_allele(a->get_length()-1)]->next;
    while(et_ptr != NULL){
        if(et_ptr->node == *a->get_allele(0)){ 
            et_ptr->common = COMMON; 
            if(cycle_id == CIRCUIT_A) et_ptr->incoming_A = false; else et_ptr->incoming_B = false; 
            was_found = true; 
            break; 
        }
        if(et_ptr->next == NULL) break;
        et_ptr = et_ptr->next;
    }
    // Add ending edge
    if(!was_found){
        et_ptr->next = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
        et_ptr->next->node = *a->get_allele(0);
        et_ptr->next->common = UNCOMMON;
        if(cycle_id == CIRCUIT_A) et_ptr->next->incoming_A = false; else et_ptr->next->incoming_B = false;
        et_ptr->next->belongs_to_cycle = cycle_id;
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
            e_table[i]->already_surrogate = false;
            e_table[i]->is_entry_cycle_A = false;
            e_table[i]->is_entry_cycle_B = false;
            e_table[i]->is_exit_cycle_A = false;
            e_table[i]->is_exit_cycle_B = false;
            e_table[i]->out_node = NULL;
            e_table[i]->connects_partition = -1;
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
void mark_entries_and_exists(uint64_t n_nodes, Edge_T<T> ** e_table, std::queue<Edge_T<T> *> * entries_A, std::queue<Edge_T<T> *> * entries_B, std::queue<Edge_T<T> *> * exits_A, std::queue<Edge_T<T> *> * exits_B){
    uint64_t i;
    uint64_t n_edges_circuit_A, n_edges_circuit_B;
    Edge_T<T> * satisfies_common;
    //bool is_entry_A = false, is_entry_B = false;
    Edge_T<T> * ptr = NULL;
    for(i=0;i<n_nodes;i++){
        if(e_table[i]->degree == 3 && e_table[i]->partition != -1 && e_table[i]->n_commons == 1){ // It has to have 1 common and 1 for each circuit, and be in a partition
            ptr = e_table[i]->next;
            n_edges_circuit_A = 0;
            n_edges_circuit_B = 0;
            satisfies_common = NULL;
            Pair<Edge_T<uint64_t>> my_pair = abstract_replace_surrogate_by_one(e_table, i);
            if(my_pair._e1 == NULL || my_pair._e2 == NULL) continue;
            while(ptr != NULL){
                
                if(ptr->common == COMMON) satisfies_common = ptr;
                if(ptr->common == UNCOMMON && ptr->belongs_to_cycle == CIRCUIT_A) n_edges_circuit_A++;
                if(ptr->common == UNCOMMON && ptr->belongs_to_cycle == CIRCUIT_B) n_edges_circuit_B++;
                               

                
                ptr = ptr->next;
            }
            if(satisfies_common != NULL && n_edges_circuit_A == 1){
                if(satisfies_common->incoming_A){ // entry
                    e_table[i]->is_entry_cycle_A = true; e_table[i]->is_exit_cycle_A = false;
                    entries_A->push(e_table[i]);    
                }else{ // exit
                    e_table[i]->is_entry_cycle_A = false; e_table[i]->is_exit_cycle_A = true; 
                    exits_A->push(e_table[i]);
                }                    
            }
            if(satisfies_common != NULL && n_edges_circuit_B == 1){
                if(satisfies_common->incoming_B){ // entry
                    e_table[i]->is_entry_cycle_B = true; e_table[i]->is_exit_cycle_B = false;
                    entries_B->push(e_table[i]);

                }else{ // exit
                    e_table[i]->is_entry_cycle_B = false; e_table[i]->is_exit_cycle_B = true; 
                    exits_B->push(e_table[i]);
                }                    
            }
            
        }
    }
}

template <class T>
void mark_entries_and_exists_ghosted(uint64_t n_nodes, Edge_T<T> ** e_table, std::queue<Edge_T<T> *> * entries_A, std::queue<Edge_T<T> *> * entries_B, std::queue<Edge_T<T> *> * exits_A, std::queue<Edge_T<T> *> * exits_B){
    uint64_t i;
    uint64_t n_edges_circuit_A, n_edges_circuit_B;
    Edge_T<T> * satisfies_common;
    //bool is_entry_A = false, is_entry_B = false;
    Edge_T<T> * ptr = NULL;
    for(i=0;i<2*n_nodes;i++){
        if(e_table[i] != NULL && e_table[i]->degree == 3 && e_table[i]->partition != -1 && e_table[i]->n_commons == 1){ // It has to have 1 common and 1 for each circuit, and be in a partition
            ptr = e_table[i]->next;
            n_edges_circuit_A = 0;
            n_edges_circuit_B = 0;
            satisfies_common = NULL;
            Pair<Edge_T<uint64_t>> my_pair = abstract_replace_surrogate_by_one(e_table, i);
            if(my_pair._e1 == NULL || my_pair._e2 == NULL) continue;
            while(ptr != NULL){
                
                if(ptr->common == COMMON) satisfies_common = ptr;
                if(ptr->common == UNCOMMON && ptr->belongs_to_cycle == CIRCUIT_A) n_edges_circuit_A++;
                if(ptr->common == UNCOMMON && ptr->belongs_to_cycle == CIRCUIT_B) n_edges_circuit_B++;
                               

                
                ptr = ptr->next;
            }
            if(satisfies_common != NULL && n_edges_circuit_A == 1){
                if(satisfies_common->incoming_A){ // entry
                    e_table[i]->is_entry_cycle_A = true; e_table[i]->is_exit_cycle_A = false;
                    entries_A->push(e_table[i]);    
                }else{ // exit
                    e_table[i]->is_entry_cycle_A = false; e_table[i]->is_exit_cycle_A = true; 
                    exits_A->push(e_table[i]);
                }                    
            }
            if(satisfies_common != NULL && n_edges_circuit_B == 1){
                if(satisfies_common->incoming_B){ // entry
                    e_table[i]->is_entry_cycle_B = true; e_table[i]->is_exit_cycle_B = false;
                    entries_B->push(e_table[i]);

                }else{ // exit
                    e_table[i]->is_entry_cycle_B = false; e_table[i]->is_exit_cycle_B = true; 
                    exits_B->push(e_table[i]);
                }                    
            }
            
        }
    }
}



template <class T>
Pair<Edge_T<T>> exit_from_entry(Edge_T<T> ** e_table, Edge_T<T> * entry, unsigned char CIRCUIT){
    Edge_T<T> * ptr = NULL;
    Edge_T<T> * ptr_traveler = entry;
    Edge_T<T> * last = entry;
    ptr = entry->next;

    if(CIRCUIT == CIRCUIT_A){
        while(ptr_traveler->is_exit_cycle_A == false){ // Until finding the exit
            while(ptr != NULL && ptr->belongs_to_cycle != CIRCUIT && last->node == ptr->node) ptr = ptr->next;
            if(ptr == NULL) terror("Reached null?");
            last = ptr;
            ptr_traveler = e_table[ptr->node];
        }
    }else{
        while(ptr_traveler->is_exit_cycle_B == false){ // Until finding the exit
            while(ptr != NULL && ptr->belongs_to_cycle != CIRCUIT && last->node == ptr->node) ptr = ptr->next;
            if(ptr == NULL) terror("Reached null?");
            last = ptr;
            ptr_traveler = e_table[ptr->node];
        }
    }
    
    Pair<Edge_T<T>> entry_and_exit;
    entry_and_exit._e1 = entry;
    entry_and_exit._e2 = ptr_traveler;
    return entry_and_exit;
}

int compare_Edges_T(const void * a, const void * b){
    Edge_T<uint64_t> * A = (Edge_T<uint64_t> *) a;
    Edge_T<uint64_t> * B = (Edge_T<uint64_t> *) b;
    if(A->node < B->node) return -1;
    if(A->node == B->node) return 0;
    return 1;
}

template <class T>
Feasible<T> verify_entries_and_exits(uint64_t n_partitions, std::queue<Edge_T<T> *> * entries_A, std::queue<Edge_T<T> *> * entries_B, std::queue<Edge_T<T> *> * exits_A, std::queue<Edge_T<T> *> * exits_B, memory_pool * mp, Edge_T<T> ** e_table){
    
    std::queue<Edge_T<T> *> aux_queue;
    Feasible<T> feasible; 
    feasible.n_entries = (uint64_t *) mp->request_bytes(n_partitions*sizeof(uint64_t));
    for(uint64_t i=0;i<n_partitions;i++){ feasible.n_entries[i] = 0; }
    uint64_t n_nodes_part_B[n_partitions]; for(uint64_t i=0;i<n_partitions;i++){ n_nodes_part_B[i] = 0; }
    //Pair<feasible_partition<T> **> feasibility;
    feasible_partition<T> ** parts_A = (feasible_partition<T> **) mp->request_bytes(n_partitions*sizeof(feasible_partition<T> *));
    feasible_partition<T> ** parts_B = (feasible_partition<T> **) mp->request_bytes(n_partitions*sizeof(feasible_partition<T> *));
    
    
    feasible.feasible._e1 = parts_A;
    feasible.feasible._e2 = parts_B;
    uint64_t part_indexes[n_partitions]; for(uint64_t i=0;i<n_partitions;i++){ part_indexes[i] = 0; }

    // Check how many for each partition 
    while(!entries_A->empty()){
        Edge_T<T> * ptr = entries_A->front();
        while(!entries_A->empty() && ptr->partition == entries_A->front()->partition){
            ptr = entries_A->front();
            #ifdef VERBOSE
            std::cout << "$ " << entries_A->front()->node << ", " << entries_A->front()->partition << std::endl;
            #endif
            aux_queue.push(ptr);
            entries_A->pop();
            feasible.n_entries[ptr->partition]++;
        }    
    }
    // Reinsert
    while(!aux_queue.empty()){ entries_A->push(aux_queue.front()); aux_queue.pop(); }
    // Generate arrays for partitioning
    for(uint64_t i=0;i<n_partitions;i++){
        parts_A[i] = (feasible_partition<T> *) mp->request_bytes(feasible.n_entries[i] * sizeof(feasible_partition<T>));
    }
    #ifdef VERBOSE
    std::cout << " -----------\n";
    #endif
    // Copy references
    while(!entries_A->empty()){
        Edge_T<T> * ptr = entries_A->front();
        while(!entries_A->empty() && ptr->partition == entries_A->front()->partition){
            ptr = entries_A->front();
            
            parts_A[ptr->partition][part_indexes[ptr->partition]].entry = &(*entries_A->front());
            part_indexes[ptr->partition]++;
            entries_A->pop();
        }    
    }

    #ifdef VERBOSE
    getchar();
    #endif

    for(uint64_t i=0;i<n_partitions;i++){ part_indexes[i] = 0; }
    // Check how many for each partition 
    while(!entries_B->empty()){
        Edge_T<T> * ptr = entries_B->front();
        while(!entries_B->empty() && ptr->partition == entries_B->front()->partition){
            ptr = entries_B->front();
            #ifdef VERBOSE
            std::cout << "& " << entries_B->front()->node << ", " << entries_B->front()->partition << std::endl;
            #endif
            aux_queue.push(ptr);
            entries_B->pop();
            n_nodes_part_B[ptr->partition]++;
        }    
    }
    // Reinsert
    while(!aux_queue.empty()){ entries_B->push(aux_queue.front()); aux_queue.pop(); }
    // Generate arrays for partitioning
    for(uint64_t i=0;i<n_partitions;i++){
        parts_B[i] = (feasible_partition<T> *) mp->request_bytes(n_nodes_part_B[i] * sizeof(feasible_partition<T>));
    }
    #ifdef VERBOSE
    std::cout << " -----------\n";
    #endif
    // Copy references
    while(!entries_B->empty()){
        Edge_T<T> * ptr = entries_B->front();
        while(!entries_B->empty() && ptr->partition == entries_B->front()->partition){
            ptr = entries_B->front();
            parts_B[ptr->partition][part_indexes[ptr->partition]].entry = &(*entries_B->front());
            part_indexes[ptr->partition]++;
            entries_B->pop();
        } 
    }
    

    // Go through all 
    uint64_t j;
    for(uint64_t i=0;i<n_partitions;i++){
        // If different number of entries, deactivate partition
        if(feasible.n_entries[i] != n_nodes_part_B[i]){ parts_A[i] = NULL; parts_B[i] = NULL; continue; }
        #ifdef VERBOSE
        std::cout << "at this level: " << feasible.n_entries[i]  << "\n";
        #endif
        for(j = 0; j<feasible.n_entries[i]; j++){

            Pair<Edge_T<T>> pA, pB;

            #ifdef VERBOSE
            std::cout << i << " - " << j << std::endl;
            std::cout << "(A) " << parts_A[i][j].entry->node << ", " << parts_A[i][j].entry->partition << " ends in ";
            #endif
            
            uint64_t length_A, length_B;

            pA = abstract_replace_surrogate_by_one_circuited(e_table, parts_A[i][j].entry->node, CIRCUIT_A, &length_A);
            
            #ifdef VERBOSE
            std::cout << "\n";
            #endif

            if(pA._e2 == NULL){ 
                #ifdef VERBOSE
                std::cout << "Make null (1)\n"; 
                #endif
                //parts_A[i] = NULL; parts_B[i] = NULL; goto out_of_part; 
                parts_A[i][j].active = false; parts_B[i][j].active = false; continue;
            }
            
            #ifdef VERBOSE
            std::cout << "(B) " << parts_B[i][j].entry->node << ", " << parts_B[i][j].entry->partition << " ends in ";
            #endif
            
            pB = abstract_replace_surrogate_by_one_circuited(e_table, parts_B[i][j].entry->node, CIRCUIT_B, &length_B);

            #ifdef VERBOSE
            std::cout << "\n";
            #endif
            
            if(pB._e2 == NULL){ 
                #ifdef VERBOSE
                std::cout << "Make null (2)\n";
                #endif
                //parts_A[i] = NULL; parts_B[i] = NULL; goto out_of_part; 
                parts_A[i][j].active = false; parts_B[i][j].active = false; continue;
            }
            
            if(length_A != length_B){ 
                #ifdef VERBOSE
                std::cout << "the lengths differ " << length_A <<" "  <<length_B << "\n"; 
                #endif
                //parts_A[i] = NULL; parts_B[i] = NULL; goto out_of_part; 
                parts_A[i][j].active = false; parts_B[i][j].active = false; continue;
            }

            // If exiting vertices are different this partition is not feasible
            //if(pA._e2->node != pB._e2->node && pA._e2->node != pB._e1->node || (pA._e1->node == pB._e1 || pA._e1->node == pB._e2) ){ 
            if((pA._e1 != pB._e1 && pA._e1 != pB._e2) || (pA._e2 != pB._e1 && pA._e2 != pB._e2)){
                #ifdef VERBOSE
                std::cout << "Destroying" << std::endl;
                #endif
                //parts_A[i] = NULL; parts_B[i] = NULL; goto out_of_part; 
                parts_A[i][j].active = false; parts_B[i][j].active = false; continue;
            }else{
                
                parts_A[i][j].exit = pA._e2;
                parts_B[i][j].exit = pB._e2;
                parts_B[i][j].reverse = false;

                if(pA._e2->node == pB._e1->node){
                    // They are reversed 

                    // USE ONLY FOR SYMMETRIC TSP 
                    /*
                    parts_B[i][j].reverse = true;
                    std::cout << "reversed to: " << pB._e1->node << "; " << pB._e2->node << "\n"; 
                    */
                    // FOR ASYMMETRIC
                    // Since they are reversed, destroy 
                    #ifdef VERBOSE
                    std::cout << "Destroying because of asymmetry" << std::endl;
                    #endif
                    //parts_A[i] = NULL; parts_B[i] = NULL; goto out_of_part; 
                    parts_A[i][j].active = false; parts_B[i][j].active = false; continue;

                }
                
                
            }
        }
        #ifdef VERBOSE
        std::cout << "Partition " << i << " is feasible!" << std::endl;
        getchar();
        #endif

        out_of_part:
        ;
    }
    /*    
    Multipartitioning<T> * multiparts = (Multipartitioning<T> *) mp->request_bytes(sizeof(Multipartitioning<T>));
    multiparts->n_parts = 0;

    Edge_T<T> * e_ptr_A = NULL;
    Edge_T<T> * e_ptr_B = NULL;
    */

    


    /*
    std::cout << "Entries A: "; while(!entries_A->empty()){ std::cout << entries_A->front()->node << ", "; entries_A->pop(); }
    std::cout << "Exits A: "; while(!exits_A->empty()){ std::cout << exits_A->front()->node << ", "; exits_A->pop(); }
    std::cout << "Entries B: "; while(!entries_B->empty()){ std::cout << entries_B->front()->node << ", "; entries_B->pop(); }
    std::cout << "Exits B: "; while(!exits_B->empty()){ std::cout << exits_B->front()->node << ", "; exits_B->pop(); }
    */
    
    //feasible.feasible = feasibility;
    return feasible;
}

template <class T>
void add_ghost_vertices(uint64_t n_nodes, Edge_T<T> ** e_table, memory_pool * mp){

    uint64_t i;
    Edge_T<T> * b, * c, * d, * e;
    
    for(i=0;i<n_nodes;i++){
        if(e_table[i]->degree > 3){
            // Found a vertex with degree 4, split it and add a ghost vertex

            // Make the connections i.e. if we have: 
            // a | b c d e 
            // then
            // a | b c 
            // and 
            // a' | d e 


            // Lower degree 
            e_table[i]->degree = 3;
            e_table[i]->n_commons = 1;


            // First, insert the new node at position n_nodes+i 
            e_table[n_nodes+i] = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
            e_table[n_nodes+i]->node = n_nodes+i;
            e_table[n_nodes+i]->degree = 3;
            e_table[n_nodes+i]->n_commons = 1; // The weighted zero
            e_table[n_nodes+i]->partition = -1; // Yet unpartitioned
            
            
            // Get references 
            b = e_table[i]->next;
            c = e_table[i]->next->next;
            d = e_table[i]->next->next->next;
            e = e_table[i]->next->next->next->next;

            // Disconnect all
            e_table[n_nodes+i]->next = NULL;
            b->next = NULL;
            c->next = NULL;
            d->next = NULL;

            Edge_T<T> * last_added = NULL;
            bool incoming[4] = { false, false, false, false };

            // Reconnect to a those that are incoming 
            if((b->incoming_A && b->belongs_to_cycle == CIRCUIT_A) || (b->incoming_B && b->belongs_to_cycle == CIRCUIT_B)){
                e_table[i]->next = b;
                last_added = b;
                incoming[0] = true;
            }
            if((c->incoming_A && c->belongs_to_cycle == CIRCUIT_A) || (c->incoming_B && c->belongs_to_cycle == CIRCUIT_B)){
                if(last_added == NULL){
                    e_table[i]->next = c;
                }else{
                    last_added->next = c;
                }
                last_added = c;
                incoming[1] = true;
            }
            if((d->incoming_A && d->belongs_to_cycle == CIRCUIT_A) || (d->incoming_B && d->belongs_to_cycle == CIRCUIT_B)){
                if(last_added == NULL){
                    e_table[i]->next = d;
                }else{
                    last_added->next = d;
                }
                last_added = d;
                incoming[2] = true;
            }
            if((e->incoming_A && e->belongs_to_cycle == CIRCUIT_A) || (e->incoming_B && e->belongs_to_cycle == CIRCUIT_B)){
                if(last_added == NULL){
                    e_table[i]->next = e;
                }else{
                    last_added->next = e;
                }
                last_added = e;
                incoming[3] = true;
            }

            // Add connections between a and a'
            if(last_added == NULL){ std::cout<< "With " << i << " impossible\n"; terror("Impossible"); }
            last_added->next = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
            last_added->next->node = n_nodes+i;
            last_added->next->common = COMMON;
            last_added->next->next = NULL;
            last_added->next->incoming_A = false; // Its going out in both cases
            last_added->next->incoming_B = false;
            last_added->next->connects_partition = -1;
            last_added->next->out_node = NULL;
            
            
            // Those that are false in "incoming" should be added to a'
            // And their links to a should be removed
            last_added = NULL;
            if(!incoming[0]){
                
                e_table[n_nodes+i]->next = b;
                last_added = b;
                Edge_T<T> * edge_new = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                edge_new->node = n_nodes+i;
                edge_new->common = UNCOMMON;
                edge_new->belongs_to_cycle = b->belongs_to_cycle;
                edge_new->incoming_A = false; edge_new->incoming_B = false;
                if(b->belongs_to_cycle == CIRCUIT_A) edge_new->incoming_A = true; else edge_new->incoming_B = true;
                edge_new->next = e_table[b->node]->next;
                e_table[b->node]->next = edge_new;
                
                Edge_T<T> * previous_node = e_table[b->node];
                Edge_T<T> * killer_node = e_table[b->node]->next;
                while(killer_node != NULL && killer_node->node != i){
                    previous_node = killer_node;
                    killer_node = killer_node->next;
                }
                if(killer_node != NULL) previous_node->next = killer_node->next; else previous_node->next = NULL;
                
            }
            if(!incoming[1]){
                if(last_added == NULL){
                    e_table[n_nodes+i]->next = c;
                }else{
                    last_added->next = c;
                }
                last_added = c;
                Edge_T<T> * edge_new = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                edge_new->node = n_nodes+i;
                edge_new->common = UNCOMMON;
                edge_new->belongs_to_cycle = c->belongs_to_cycle;
                edge_new->incoming_A = false; edge_new->incoming_B = false;
                if(c->belongs_to_cycle == CIRCUIT_A) edge_new->incoming_A = true; else edge_new->incoming_B = true;
                edge_new->next = e_table[c->node]->next;
                e_table[c->node]->next = edge_new;
                
                Edge_T<T> * previous_node = e_table[c->node];
                Edge_T<T> * killer_node = e_table[c->node]->next;
                while(killer_node != NULL && killer_node->node != i){
                    previous_node = killer_node;
                    killer_node = killer_node->next;
                }
                if(killer_node != NULL) previous_node->next = killer_node->next; else previous_node->next = NULL;
                
            }
            if(!incoming[2]){
                if(last_added == NULL){
                    e_table[n_nodes+i]->next = d;
                }else{
                    last_added->next = d;
                }
                last_added = d;
                Edge_T<T> * edge_new = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                edge_new->node = n_nodes+i;
                edge_new->common = UNCOMMON;
                edge_new->belongs_to_cycle = d->belongs_to_cycle;
                edge_new->incoming_A = false; edge_new->incoming_B = false;
                if(d->belongs_to_cycle == CIRCUIT_A) edge_new->incoming_A = true; else edge_new->incoming_B = true;
                edge_new->next = e_table[d->node]->next;
                e_table[d->node]->next = edge_new;
                
                Edge_T<T> * previous_node = e_table[d->node];
                Edge_T<T> * killer_node = e_table[d->node]->next;
                while(killer_node != NULL && killer_node->node != i){
                    previous_node = killer_node;
                    killer_node = killer_node->next;
                }
                if(killer_node != NULL) previous_node->next = killer_node->next; else previous_node->next = NULL;
                
            }
            if(!incoming[3]){
                if(last_added == NULL){
                    e_table[n_nodes+i]->next = e;
                }else{
                    last_added->next = e;
                }
                last_added = e;
                Edge_T<T> * edge_new = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
                edge_new->node = n_nodes+i;
                edge_new->common = UNCOMMON;
                edge_new->belongs_to_cycle = e->belongs_to_cycle;
                edge_new->incoming_A = false; edge_new->incoming_B = false;
                if(e->belongs_to_cycle == CIRCUIT_A) edge_new->incoming_A = true; else edge_new->incoming_B = true;
                edge_new->next = e_table[e->node]->next;
                e_table[e->node]->next = edge_new;
                
                Edge_T<T> * previous_node = e_table[e->node];
                Edge_T<T> * killer_node = e_table[e->node]->next;
                while(killer_node != NULL && killer_node->node != i){
                    previous_node = killer_node;
                    killer_node = killer_node->next;
                }
                if(killer_node != NULL) previous_node->next = killer_node->next; else previous_node->next = NULL;
                
            }
            
            // Finally, add connections between a' and a
            last_added->next = (Edge_T<T> *) mp->request_bytes(sizeof(Edge_T<T>));
            last_added->next->node = i;
            last_added->next->common = COMMON;
            last_added->next->next = NULL;
            last_added->next->incoming_A = true; // Its going in for both cases
            last_added->next->incoming_B = true;
            last_added->next->connects_partition = -1;
            last_added->next->out_node = NULL;

            /*
            std::cout<< "In " << i << "\n";
            print_edge_tables_ghosted(n_nodes, e_table);
            */


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
bool get_highest_node_unpartitioned_ghosted(uint64_t n_nodes, Edge_T<T> ** e_table, uint64_t * node_id){
    
    *node_id = 0;
    bool found = false;
    uint64_t degree = 0;
    for(uint64_t i=0;i<2*n_nodes;i++){
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
    std::queue<Edge_T<T> *> FIFO_queue;
    FIFO_queue.push(route_start);

    while(ptr != NULL){ // Until we find no more common edges
        
        // Loop until next common edge
        while(ptr != NULL && ptr->common != COMMON){
            ptr = ptr->next;
        }
        
        if(ptr != NULL && ptr->node != last_replaced->node){
            
            // If it has already been used, abort 
            if(e_table[master_node]->already_surrogate){
                surrogate._e1 = NULL;
                surrogate._e2 = NULL;
                return surrogate;
            }
        
            // Update previous so we wont traverse it again 
            
            //printf("%" PRId64", ", e_table[ptr->node]->partition);
            last_replaced = e_table[master_node];
            master_node = ptr->node;
            route_end = e_table[ptr->node];
            
            FIFO_queue.push(e_table[ptr->node]);
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

    // Mark surrogate as used so that it does not generate more fake surrogate edges
    if(route_start != NULL && route_end != NULL){
        while(!FIFO_queue.empty()){

            FIFO_queue.front()->already_surrogate = true;
            FIFO_queue.pop();
        }
    }

    
    return surrogate;
}

template <class T>
Pair<Edge_T<T>> abstract_replace_surrogate_by_one_circuited(Edge_T<T> ** e_table, uint64_t i, uint64_t CIRCUIT, uint64_t * length){

    Edge_T<T> * ptr, * last_replaced, * route_start, * route_end;
    Pair<Edge_T<T>> surrogate; surrogate._e1 = NULL; surrogate._e2 = NULL;
    
    // Extend in both directions until one node with degree > 2 is found
    uint64_t master_node = i;
    int64_t start_part = e_table[i]->partition;
    route_start = e_table[i];
    *length = 0;
    route_end = NULL;
    last_replaced = e_table[i];
    //printf("%" PRId64", ", e_table[last_replaced->node]->partition);
    ptr = e_table[i]->next;

    while(ptr != NULL){ // Until we find no more common edges
        
        // Loop until next uncommon and from the circuit
        /*
        while(ptr->belongs_to_cycle != CIRCUIT && ptr->common == COMMON){
            if(ptr == NULL) break;
            ptr = ptr->next;
        }
        */
        // It doesnt matter if its common as long as it connects the same partition
        //if(ptr != NULL && ptr->node != last_replaced->node && (ptr->belongs_to_cycle == CIRCUIT || ptr->common == COMMON)){
        if(ptr != NULL && ptr->node != last_replaced->node && (ptr->belongs_to_cycle == CIRCUIT || ptr->common == COMMON)){
            if((CIRCUIT == CIRCUIT_A && ptr->incoming_A == false) || (CIRCUIT == CIRCUIT_B && ptr->incoming_B == false)){
                // Update previous so we wont traverse it again 
            
                #ifdef VERBOSE
                std::cout<< ptr->node << ", ";
                #endif
                last_replaced = e_table[master_node];
                master_node = ptr->node;
                route_end = e_table[ptr->node];

                if(CIRCUIT == CIRCUIT_A && route_end->is_exit_cycle_A) goto finish; 
                if(CIRCUIT == CIRCUIT_B && route_end->is_exit_cycle_B) goto finish; 

                // Not null implies we found common edge 
                // Save this node as last and look for more 
                // Jump common tour
                
                (*length)++;
                ptr = e_table[ptr->node]->next;
                
            }else{
                if(ptr != NULL) ptr = ptr->next;
            }
        

            
        }else{
            if(ptr != NULL) ptr = ptr->next;
        }
        
        
    }
    std::cout << std::endl;
    finish:
    //printf("\n");
    if(route_end->is_exit_cycle_A || route_end->is_exit_cycle_B){
        surrogate._e1 = route_start;
        surrogate._e2 = route_end;
        #ifdef VERBOSE
        std::cout << "The good way\n";
        #endif
    }else{
        #ifdef VERBOSE
        std::cout << "The bad way\n";
        #endif
        surrogate._e1 = NULL;
        surrogate._e2 = NULL;
    }
    return surrogate;
}

template <class T>
Pair<Edge_T<T>> abstract_replace_surrogate_by_one(Edge_T<T> ** e_table, uint64_t i){

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
    if(route_start->partition != route_end->partition){
        surrogate._e1 = route_start;
        surrogate._e2 = route_end;
    }else{
        surrogate._e1 = NULL;
        surrogate._e2 = NULL;
    }
    

    
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
void generate_partitions(PXTable<T> * px_table, Edge_T<T> ** e_table, uint64_t n_nodes, memory_pool * mp){
    uint64_t i;
    Pair<Edge_T<T>> pair;
    List<Surrogate_Edge_T<T>> * list_su_p1 = NULL, * list_su_p2 = NULL;

    for(i=0; i<n_nodes; i++){

        // Only find surrogate edges from one vertex that is connected to "something"
        if(e_table[i]->n_commons == 1){
            pair = replace_surrogate_by_one(e_table, i);
            
            

            if(pair._e1 == NULL || pair._e2 == NULL) {
                continue; 
            }else{
                e_table[pair._e1->node]->connects_partition = pair._e2->partition;
                e_table[pair._e1->node]->out_node = e_table[pair._e2->node];
                e_table[pair._e2->node]->connects_partition = pair._e1->partition;
                e_table[pair._e2->node]->out_node = e_table[pair._e1->node];
            }
            if(pair._e1->partition == pair._e2->partition) continue;
            

            
            std::cout << "SG " << i << " connects (" << pair._e1->partition << ", " << pair._e2->partition << ") " << std::endl;
            //std::cout << " \t " << pair._e1->node << "(" << pair._e1->partition << " ), " << pair._e2->node << "(" << pair._e2->partition << ")" << std::endl;

            // Get memory 
            list_su_p1 = (List<Surrogate_Edge_T<T>> *) mp->request_bytes(sizeof(List<Surrogate_Edge_T<T>>));
            list_su_p1->v.left = pair._e1;
            list_su_p1->v.right = pair._e2;
            list_su_p2 = (List<Surrogate_Edge_T<T>> *) mp->request_bytes(sizeof(List<Surrogate_Edge_T<T>>));
            list_su_p2->v.left = pair._e1;
            list_su_p2->v.right = pair._e2;

            
            // Link next for p1 
            if(px_table[pair._e1->partition].su_gates != NULL){
                list_su_p1->next = px_table[pair._e1->partition].su_gates;
                px_table[pair._e1->partition].su_gates->prev = list_su_p1;
                px_table[pair._e1->partition].su_gates = list_su_p1;
                list_su_p1->prev = NULL;
            }else{
                px_table[pair._e1->partition].su_gates = list_su_p1;
                list_su_p1->prev = NULL;
                list_su_p1->next = NULL;
            }

            
            // Link next for p2 
            if(px_table[pair._e2->partition].su_gates != NULL){
                list_su_p2->next = px_table[pair._e2->partition].su_gates;
                px_table[pair._e2->partition].su_gates->prev = list_su_p2;
                px_table[pair._e2->partition].su_gates = list_su_p2;
                list_su_p2->prev = NULL;
            }else{
                px_table[pair._e2->partition].su_gates = list_su_p2;
                list_su_p2->prev = NULL;
                list_su_p2->next = NULL;
            }

            px_table[pair._e1->partition].n_surrogate_edges++;
            px_table[pair._e2->partition].n_surrogate_edges++;
            
        }
    }
}

template <class T>
void shorten_common_tours_ghosted(Edge_T<T> ** e_table, uint64_t n_nodes){
    uint64_t i;
    Pair<Edge_T<T>> pair;

    for(i=0; i<2*n_nodes; i++){

        // Only find surrogate edges from one vertex that is connected to "something"
        if(e_table[i] != NULL && e_table[i]->n_commons == 1){
            pair = replace_surrogate_by_one(e_table, i);
            
            if(pair._e1 == NULL || pair._e2 == NULL) {
                continue; 
            }else{
                e_table[pair._e1->node]->connects_partition = pair._e2->partition;
                e_table[pair._e1->node]->out_node = e_table[pair._e2->node];
                e_table[pair._e2->node]->connects_partition = pair._e1->partition;
                e_table[pair._e2->node]->out_node = e_table[pair._e1->node];
            }

        }
    }
}

template <class T>
T evaluate_partition_subtours_multiple(Edge_T<T> * start, Edge_T<T> * end, bool reverse, Chromosome<T> * c, void * solution_info, Edge_T<T> ** e_table){
    T score = 0;
    uint64_t swath_start_P1 = 0;
    uint64_t i;
    uint64_t n_nodes = c->get_length();
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix * ) solution_info;

    for(i=0;i<n_nodes;i++){
        // Find first swath in P1 
        if(reverse == false && *c->get_allele(i) == start->node) swath_start_P1 = i;
        if(reverse == true && *c->get_allele(i) == end->node) swath_start_P1 = i;
    }

    i = swath_start_P1;

    std::cout << " adding: ";

    if(reverse == false){
        std::cout << " direct ";
        do{

            
            i = (i + 1) % n_nodes;
            
            //std::cout << *c->get_allele(i) << "(" << e_table[*c->get_allele(i)]->partition << "), ";

            // Compute here
            uint64_t prescore = score;
            if(i == 0) score += tsp->dist[*c->get_allele(n_nodes-1)][*c->get_allele(i)];
            else score += tsp->dist[*c->get_allele(i-1)][*c->get_allele(i)];
            std::cout << *c->get_allele(i) << ",";

        //}while(e_table[*c->get_allele(i)]->partition == partition1 ||  e_table[*c->get_allele(i)]->partition == -1);
        //}while(e_table[*c->get_allele(i)]->partition != partition2);
        //}while(*c->get_allele(i) != end->left->node && *c->get_allele(i) != end->right->node);
        }while(*c->get_allele(i) != end->node);
    }else{
        std::cout << " reverse ";
        do{

            //if(i==0) i = n_nodes-1; else i--;
            
            //std::cout << *c->get_allele(i) << "(" << e_table[*c->get_allele(i)]->partition << "), ";

            // Compute here
            /*
            uint64_t prescore = score;
            if(i == n_nodes-1) score += tsp->dist[*c->get_allele(0)][*c->get_allele(n_nodes-1)];
            else score += tsp->dist[*c->get_allele(i+1)][*c->get_allele(i)];
            std::cout << *c->get_allele(i) << ",";
            */
            i = (i + 1) % n_nodes;
            
            //std::cout << *c->get_allele(i) << "(" << e_table[*c->get_allele(i)]->partition << "), ";

            // Compute here
            uint64_t prescore = score;
            if(i == 0) score += tsp->dist[*c->get_allele(n_nodes-1)][*c->get_allele(i)];
            else score += tsp->dist[*c->get_allele(i-1)][*c->get_allele(i)];
            std::cout << *c->get_allele(i) << ",";

        //}while(e_table[*c->get_allele(i)]->partition == partition1 ||  e_table[*c->get_allele(i)]->partition == -1);
        //}while(e_table[*c->get_allele(i)]->partition != partition2);
        //}while(*c->get_allele(i) != end->left->node && *c->get_allele(i) != end->right->node);
        }while(*c->get_allele(i) != start->node);
    }

    
    
    std::cout << std::endl;
    return score;
} 


template <class T>
long double evaluate_partition_subtours_multiple_ghosted(Edge_T<T> * start, Edge_T<T> * end, bool reverse, Chromosome<T> * c, void * solution_info, Edge_T<T> ** e_table){
    long double score = 0;
    uint64_t swath_start_P1 = 0;
    uint64_t i;
    uint64_t n_nodes_unghosted = c->get_length();
    uint64_t n_nodes = n_nodes_unghosted; //2*c->get_length();
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix * ) solution_info;

    for(i=0;i<n_nodes_unghosted;i++){
        // Find first swath in P1 
        if(reverse == false && (*c->get_allele(i) == start->node || n_nodes_unghosted+(*c->get_allele(i)) == start->node)) swath_start_P1 = i; 
        if(reverse == true && (*c->get_allele(i) == start->node || n_nodes_unghosted+(*c->get_allele(i)) == start->node)) swath_start_P1 = i;
    }

    i = swath_start_P1;

    #ifdef VERBOSE
    std::cout << " adding: (" << i << ")";
    #endif

    if(reverse == false){
        #ifdef VERBOSE
        std::cout << " direct ";
        #endif
        do{

            
            i = (i + 1) % n_nodes;
            
            //std::cout << *c->get_allele(i) << "(" << e_table[*c->get_allele(i)]->partition << "), ";

            // Compute here
            
            if(i == 0) score += (long double) tsp->dist[*c->get_allele(n_nodes-1)][*c->get_allele(i)];
            else score += (long double) tsp->dist[*c->get_allele(i-1)][*c->get_allele(i)];
            #ifdef VERBOSE
            std::cout << *c->get_allele(i) << ",";
            #endif

        //}while(e_table[*c->get_allele(i)]->partition == partition1 ||  e_table[*c->get_allele(i)]->partition == -1);
        //}while(e_table[*c->get_allele(i)]->partition != partition2);
        //}while(*c->get_allele(i) != end->left->node && *c->get_allele(i) != end->right->node);
        }while(*c->get_allele(i) != end->node && n_nodes_unghosted+(*c->get_allele(i)) != end->node);
    }else{
        #ifdef VERBOSE
        std::cout << " reverse ";
        #endif
        do{

            //if(i==0) i = n_nodes-1; else i--;
            
            //std::cout << *c->get_allele(i) << "(" << e_table[*c->get_allele(i)]->partition << "), ";

            // Compute here
            /*
            uint64_t prescore = score;
            if(i == n_nodes-1) score += tsp->dist[*c->get_allele(0)][*c->get_allele(n_nodes-1)];
            else score += tsp->dist[*c->get_allele(i+1)][*c->get_allele(i)];
            std::cout << *c->get_allele(i) << ",";
            */
            i = (i + 1) % n_nodes;
            
            //std::cout << *c->get_allele(i) << "(" << e_table[*c->get_allele(i)]->partition << "), ";

            // Compute here
            
            if(i == 0) score += (long double) tsp->dist[*c->get_allele(n_nodes-1)][*c->get_allele(i)];
            else score += (long double) tsp->dist[*c->get_allele(i-1)][*c->get_allele(i)];
            #ifdef VERBOSE
            std::cout << *c->get_allele(i) << ",";
            #endif

        //}while(e_table[*c->get_allele(i)]->partition == partition1 ||  e_table[*c->get_allele(i)]->partition == -1);
        //}while(e_table[*c->get_allele(i)]->partition != partition2);
        //}while(*c->get_allele(i) != end->left->node && *c->get_allele(i) != end->right->node);
        }while(*c->get_allele(i) != start->node && n_nodes_unghosted+(*c->get_allele(i)) != end->node);
    }

    
    #ifdef VERBOSE
    std::cout << std::endl;
    #endif
    return score;
} 

template <class T>
T evaluate_partition_subtours(Surrogate_Edge_T<T> * start, Surrogate_Edge_T<T> * end, Chromosome<T> * c, void * solution_info, int64_t partition1, int64_t partition2, Edge_T<T> ** e_table){
    T score = 0;
    uint64_t swath_start_P1 = 0, swath_end_P1 = 0;
    uint64_t i, aux;
    uint64_t n_nodes = c->get_length();
    Sol_TSP_matrix * tsp = (Sol_TSP_matrix * ) solution_info;

    for(i=0;i<n_nodes;i++){
        // Find first swath in P1 
        if(*c->get_allele(i) == start->left->node) swath_start_P1 = i;
        if(*c->get_allele(i) == start->right->node) swath_end_P1 = i;
    }

    // Check that between swath_start and swath_end there are no other nodes in between in P2
    uint64_t minus_pos = (swath_start_P1 == 0) ? (n_nodes-1) : (swath_start_P1 - 1);
    bool go_backwards_P1 = false;
    std::cout << ":: " << e_table[*c->get_allele((swath_start_P1 + 1) % n_nodes)]->partition << std::endl;
    std::cout << ":: " << e_table[*c->get_allele(minus_pos)]->partition << std::endl;

    if(e_table[*c->get_allele((swath_start_P1 + 1) % n_nodes)]->partition == partition1){  goto escape_1_O1; }
    else if(e_table[*c->get_allele(minus_pos)]->partition == partition1){go_backwards_P1 = true; goto escape_1_O1; }
    // If it gets here, it means that we have to reverse start with end 
    else{
        aux = swath_start_P1;
        swath_start_P1 = swath_end_P1;
        swath_end_P1 = aux;
    }
    // And re-check orders
    minus_pos = (swath_start_P1 == 0) ? (n_nodes-1) : (swath_start_P1 - 1);
    if(e_table[*c->get_allele((swath_start_P1 + 1) % n_nodes)]->partition == partition1){  goto escape_1_O1; }
    else if(e_table[*c->get_allele(minus_pos)]->partition == partition1){ go_backwards_P1 = true;goto escape_1_O1; }
    
    std::cout << ":: " << e_table[*c->get_allele((swath_start_P1 + 1) % n_nodes)]->partition << std::endl;
    std::cout << ":: " << e_table[*c->get_allele(minus_pos)]->partition << std::endl;


    escape_1_O1:

    // At this point we have the swath position starting and ending in both P1 and P2 
    // and we know if the swaths are to be copied from left to right or right to left 

    // Proceed to copy 
    i = swath_start_P1;

    std::cout << "Arriving partitions: " << partition1 << " - " << partition2 << std::endl;

    do{

        // Advance
        if(go_backwards_P1){
            if(i>0) i--; else i=n_nodes-1;
        }else{
            i = (i + 1) % n_nodes;
        }
        std::cout << *c->get_allele(i) << "(" << e_table[*c->get_allele(i)]->partition << "), ";

        // Compute here
        if(i == 0) score += tsp->dist[*c->get_allele(n_nodes-1)][*c->get_allele(i)];
        else score += tsp->dist[*c->get_allele(i-1)][*c->get_allele(i)];

    //}while(e_table[*c->get_allele(i)]->partition == partition1 ||  e_table[*c->get_allele(i)]->partition == -1);
    //}while(e_table[*c->get_allele(i)]->partition != partition2);
    //}while(*c->get_allele(i) != end->left->node && *c->get_allele(i) != end->right->node);
    }while(e_table[*c->get_allele(i)]->partition != partition2 && *c->get_allele(i) != end->left->node && *c->get_allele(i) != end->right->node);
    
    std::cout << std::endl;
    return score;
} 




template <class T>
void apply_PX_chromosomes(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * px, Chromosome<T> * P1, Chromosome<T> * P2, Chromosome<T> * offspring_1, Chromosome<T> * offspring_2){

    uint64_t swath_start_P1 = 0, swath_start_P2 = 0, swath_end_P1 = 0, swath_end_P2 = 0;
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

template <class T>
void apply_PX_chromosomes_best(uint64_t n_nodes, Edge_T<T> ** e_table, Quartet<Edge_T<T>> * px, Chromosome<T> * P1, Chromosome<T> * P2, Chromosome<T> * offspring_1){

    uint64_t swath_start_P1 = 0, swath_start_P2 = 0, swath_end_P1 = 0, swath_end_P2 = 0;
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
    
}



template void single_point_crossover<unsigned char>(Chromosome<unsigned char> * a, Chromosome<unsigned char> * b, Chromosome<unsigned char> * replacement, Manager<unsigned char> * m);
template void ordered_crossover<uint64_t>(Chromosome<uint64_t> * a, Chromosome<uint64_t> * b, Chromosome<uint64_t> * replacement, Manager<uint64_t> * m);
template void ordered_crossover<unsigned char>(Chromosome<unsigned char> * a, Chromosome<unsigned char> * b, Chromosome<unsigned char> * replacement, Manager<unsigned char> * m);
template void fill_edge_table(Chromosome<uint64_t> * a, Edge_T<uint64_t> ** e_table, memory_pool * mp, uint64_t cycle_id);
template void generate_degree(uint64_t n_nodes, Edge_T<uint64_t> ** e_table);
template void add_ghost_vertices(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, memory_pool * mp);
template void mark_entries_and_exists(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, std::queue<Edge_T<uint64_t> *> * entries_A, std::queue<Edge_T<uint64_t> *> * entries_B, std::queue<Edge_T<uint64_t> *> * exits_A, std::queue<Edge_T<uint64_t> *> * exits_B);
template void mark_entries_and_exists_ghosted(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, std::queue<Edge_T<uint64_t> *> * entries_A, std::queue<Edge_T<uint64_t> *> * entries_B, std::queue<Edge_T<uint64_t> *> * exits_A, std::queue<Edge_T<uint64_t> *> * exits_B);
template Pair<Edge_T<uint64_t>> exit_from_entry(Edge_T<uint64_t> ** e_table, Edge_T<uint64_t> * entry, unsigned char CIRCUIT);
template bool get_highest_node_unpartitioned(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, uint64_t * node_id);
template bool get_highest_node_unpartitioned_ghosted(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, uint64_t * node_id);
template void find_connected_components(uint64_t init_node, int64_t partition_label, Edge_T<uint64_t> ** e_table, std::queue<uint64_t> * FIFO_queue);
template bool is_connected_to(Edge_T<uint64_t> ** e_table, uint64_t node_1, uint64_t node_2);
template void find_surrogate_edge_that_partitionates(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, Quartet<Edge_T<uint64_t>> * surrogates);
template Feasible<uint64_t> verify_entries_and_exits(uint64_t n_partitions, std::queue<Edge_T<uint64_t> *> * entries_A, std::queue<Edge_T<uint64_t> *> * entries_B, std::queue<Edge_T<uint64_t> *> * exits_A, std::queue<Edge_T<uint64_t> *> * exits_B, memory_pool * mp, Edge_T<uint64_t> ** e_table);
template Pair<Edge_T<uint64_t>> abstract_replace_surrogate_by_one(Edge_T<uint64_t> ** e_table, uint64_t i);
template Pair<Edge_T<uint64_t>> abstract_replace_surrogate_by_one_circuited(Edge_T<uint64_t> ** e_table, uint64_t i, uint64_t CIRCUIT, uint64_t * length);
template Pair<Edge_T<uint64_t>> replace_surrogate_by_one(Edge_T<uint64_t> ** e_table, uint64_t i);
template void generate_partitions(PXTable<uint64_t> * px_table, Edge_T<uint64_t> ** e_table, uint64_t n_nodes, memory_pool * mp);
template void shorten_common_tours_ghosted(Edge_T<uint64_t> ** e_table, uint64_t n_nodes);
template uint64_t evaluate_partition_subtours_multiple(Edge_T<uint64_t> * start, Edge_T<uint64_t> * end, bool reverse, Chromosome<uint64_t> * c, void * solution_info, Edge_T<uint64_t> ** e_table);
template long double evaluate_partition_subtours_multiple_ghosted(Edge_T<uint64_t> * start, Edge_T<uint64_t> * end, bool reverse, Chromosome<uint64_t> * c, void * solution_info, Edge_T<uint64_t> ** e_table);
template uint64_t evaluate_partition_subtours(Surrogate_Edge_T<uint64_t> * start, Surrogate_Edge_T<uint64_t> * end, Chromosome<uint64_t> * c, void * solution_info, int64_t partition1, int64_t partition2, Edge_T<uint64_t> ** e_table);
template void apply_PX_chromosomes(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, Quartet<Edge_T<uint64_t>> * px, Chromosome<uint64_t> * P1, Chromosome<uint64_t> * P2, Chromosome<uint64_t> * offspring_1, Chromosome<uint64_t> * offspring_2);
template void apply_PX_chromosomes_best(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, Quartet<Edge_T<uint64_t>> * px, Chromosome<uint64_t> * P1, Chromosome<uint64_t> * P2, Chromosome<uint64_t> * offspring_1);
//template Part_list<uint64_t> * generate_lists_from_G(uint64_t n_nodes, Edge_T<uint64_t> ** e_table, memory_pool * mp);