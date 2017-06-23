#include "readstream.h"
#include <inttypes.h>
#define __STDC_FORMAT_MACROS
Readstream::Readstream(const char * file, void (*reading_function)(FILE * input, void * type_structure), void * type_structure){
    this->input = fopen64(file, "rt");
    if(this->input == NULL) throw "Could not open input file";
    this->reading_function = reading_function;
    this->type_structure = type_structure;
}

void Readstream::read(){
    this->reading_function(this->input, this->type_structure);
    fclose(this->input);
}


void reading_function_TSP(FILE * input, void * type_structure){
    Sol_TSP_matrix * tsp_mat = (Sol_TSP_matrix *) type_structure;

    uint64_t total = 0;
    uint64_t id;
    long double x, y;
    
    // Count lines
    char buffer[MAX_LINE];
    while(!feof(input) && fgets(buffer, MAX_LINE, input) != 0){

        //printf("Look pussy u read: %s\n", buffer);
        if(3 == sscanf(buffer, "%lu %Le %Le", &id, &x, &y)){
            total++;
            #ifdef VERBOSE
            fprintf(stdout, "%lu, %.3Le %.3Le\n", id, x, y);
            #endif
        }
    }
    // Go to start
    rewind(input);

    // Allocate space
    long double * aux = (long double *) std::malloc(total * sizeof(long double));
    long double * aux2 = (long double *) std::malloc(total * sizeof(long double));
    if(aux == NULL || aux2 == NULL) throw "Could not allocate temporary vector for weight matrix";

    tsp_mat->n = total;
    
    tsp_mat->dist = (long double **) std::malloc(total * sizeof(long double *));
    if(tsp_mat->dist == NULL) throw "Could not allocate TSP lib nodes (1)";
    for(uint64_t i=0;i<total;i++){
        tsp_mat->dist[i] = (long double *) malloc(total * sizeof(long double));
        if(tsp_mat->dist[i] == NULL) throw "Could not allocate TSP lib nodes (2)";
    }

    // Add nodes
    while(!feof(input) && fgets(buffer, MAX_LINE, input) != 0){

        if(3 == sscanf(buffer, "%lu %Le %Le", &id, &x, &y)){
            aux[id-1] = x;
            aux2[id-1] = y;
        }
    }
    /*
    int64_t xd, yd;
    long double rij;
    long double tij;
    */
    // Calculate distances
    for(uint64_t i=0;i<total;i++){
        for(uint64_t j=0;j<total;j++){

            // ATTENTION TODO
            // Currently using pseudo euc distance 
            // TODO parameterize this so that you dont have to re-compile

            // EUC-2D
            tsp_mat->dist[i][j] = sqrtl(powl(aux[i] - aux[j], 2.0) + powl(aux2[i] - aux2[j], 2.0));

            // ATT 
            
            /*
            xd = (int64_t)aux[i] - (int64_t)aux[j];
            yd = (int64_t)aux2[i] - (int64_t)aux2[j];
            rij = sqrtl( (xd*xd + yd*yd) / 10.0);
            tij = (long double) (round(rij));
            if(tij < rij) tsp_mat->dist[i][j] = tij+1; else tsp_mat->dist[i][j] = tij;
            */
            
            //std::cout << "Between (" << aux[i] << ", " << aux[j]  << ") and (" << aux2[i] << ", " << aux2[j]  << ") yields " << tsp_mat->dist[i][j] << std::endl;
            //std::cout << tsp_mat->dist[i][j] << " ";
            //getchar();
        }
        //std::cout << std::endl;
    }

    std::free(aux);
    std::free(aux2);
}

void reading_function_VRP(FILE * input, void * type_structure){
    /*
    struct Sol_VRP_matrix{
        long double ** dist;
        uint64_t n;
        uint64_t * demands; // Customer demands
        uint64_t depot; // Node depot
    };
    */
    
    Sol_VRP_matrix * vrp_mat = (Sol_VRP_matrix *) type_structure;

    uint64_t total = 0;
    uint64_t id;
    long double x, y;
    
    // Count lines
    char buffer[MAX_LINE];
    bool node_coord_section = false;
    while(!feof(input) && fgets(buffer, MAX_LINE, input) != 0){

        //printf("Look pussy u read: %s\n", buffer);
        if(strncmp(buffer, "CAPACITY", 8) == 0){
            sscanf(buffer, "CAPACITY : %Le", &x);
            vrp_mat->capacity = x;
        }
        if(strncmp(buffer, "DEMAND_SECTION", 14) == 0) break;

        if(node_coord_section && 3 == sscanf(buffer, "%lu %Le %Le", &id, &x, &y)){
            total++;
            #ifdef VERBOSE
            fprintf(stdout, "%lu, %.3Le %.3Le\n", id, x, y);
            #endif
        }
        if(strncmp(buffer, "NODE_COORD_SECTION", 18) == 0) node_coord_section = true;
    }
    // Go to start
    rewind(input);

    // Allocate space
    vrp_mat->points = (CPair<long double> *) std::malloc(total * sizeof(CPair<long double>));
    if(vrp_mat->points == NULL) throw "Could not allocate vector for 2D points";
    //long double * aux = (long double *) std::malloc(total * sizeof(long double));
    //long double * aux2 = (long double *) std::malloc(total * sizeof(long double));
    //if(aux == NULL || aux2 == NULL) throw "Could not allocate temporary vector for weight matrix";

    vrp_mat->n = total;
    vrp_mat->demands = (uint64_t *) std::malloc(total * sizeof(uint64_t));
    if(vrp_mat->demands == NULL) throw "Could not allocate demands";
    
    vrp_mat->dist = (long double **) std::malloc(total * sizeof(long double *));
    if(vrp_mat->dist == NULL) throw "Could not allocate TSP lib nodes (1)";
    for(uint64_t i=0;i<total;i++){
        vrp_mat->dist[i] = (long double *) malloc(total * sizeof(long double));
        if(vrp_mat->dist[i] == NULL) throw "Could not allocate TSP lib nodes (2)";
    }

    // Add nodes
    while(!feof(input) && fgets(buffer, MAX_LINE, input) != 0){
        if(strncmp(buffer, "DEMAND_SECTION", 14) == 0) break;
        if(3 == sscanf(buffer, "%lu %Le %Le", &id, &x, &y)){
            vrp_mat->points[id-1]._e1 = x;
            vrp_mat->points[id-1]._e2 = y;
        }
    }


    
    // Add demands
    while(!feof(input) && fgets(buffer, MAX_LINE, input) != 0){
        if(strncmp(buffer, "DEPOT_SECTION", 13) == 0) break;
        if(2 == sscanf(buffer, "%lu %Le", &id, &x)){
            vrp_mat->demands[id-1] = x;
            /*
            #ifdef VERBOSE
            std::cout << vrp_mat->demands[id-1] << std::endl;
            getchar();
            #endif
            */
        }
    }

    if(fgets(buffer, MAX_LINE, input) == 0) terror("No depot found");
    if(1 == sscanf(buffer, "%lu", &id)){
        vrp_mat->depot = id-1;
    }

    

    /*
    int64_t xd, yd;
    long double rij;
    long double tij;
    */
    // Calculate distances
    for(uint64_t i=0;i<total;i++){
        for(uint64_t j=0;j<total;j++){

            // ATTENTION TODO
            // Currently using pseudo euc distance 
            // TODO parameterize this so that you dont have to re-compile

            // EUC-2D
            vrp_mat->dist[i][j] = sqrtl(powl(vrp_mat->points[i]._e1 - vrp_mat->points[j]._e1, 2.0) + powl(vrp_mat->points[i]._e2 - vrp_mat->points[j]._e2, 2.0));

            // ATT 
            
            /*
            xd = (int64_t)aux[i] - (int64_t)aux[j];
            yd = (int64_t)aux2[i] - (int64_t)aux2[j];
            rij = sqrtl( (xd*xd + yd*yd) / 10.0);
            tij = (long double) (round(rij));
            if(tij < rij) tsp_mat->dist[i][j] = tij+1; else tsp_mat->dist[i][j] = tij;
            */
            
            //std::cout << "Between (" << aux[i] << ", " << aux[j]  << ") and (" << aux2[i] << ", " << aux2[j]  << ") yields " << tsp_mat->dist[i][j] << std::endl;
            //std::cout << vrp_mat->dist[i][j] << " ";
            //getchar();
        }
        //std::cout << std::endl;
    }

    /*
    std::cout << "Demands" << std::endl;
    for(uint64_t i=0;i<total;i++){
        std::cout << vrp_mat->demands[i] << " ";
    }

    std::cout << "Depot" << vrp_mat->depot << std::endl;
    */
}


void reading_function_LB_reads(FILE * input, void * type_structure){
    Sol_LB_reads * LB_mat = (Sol_LB_reads *) type_structure;


    //Current length of array and variables for the buffer
    uint64_t idx = 0, r = 0;
    //To force reading from the buffer
    idx = READBUF + 1;
    
    // Vector to read in batches
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = (char *) std::calloc(READBUF, sizeof(char))) == NULL) {
        throw "Could not allocate memory for read buffer";
    }
    // Starting size of lengths of sequences
    LB_mat->lengths = (uint64_t *) std::malloc(INIT_SEQS * sizeof(uint64_t));
    if(LB_mat->lengths == NULL) throw "Could not allocate lengths of reads";
    LB_mat->current_max = 0;
    // Number of reallocs
    uint64_t reallocs = 1;



    char c = '\0'; //Initialize
    uint64_t curr_len, n_sequences = 0;


    c = buffered_fgetc(temp_seq_buffer, &idx, &r, input);
    while((!feof(input) || (feof(input) && idx < r))){

            if(c == '>'){
                    n_sequences++;
                    curr_len = 0;
                    while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, input);

                    while((!feof(input) || (feof(input) && idx < r)) && c != '>'){
                        c = buffered_fgetc(temp_seq_buffer, &idx, &r, input);
                        curr_len++;
                    }

                    if(n_sequences == INIT_SEQS*reallocs){
                        reallocs = reallocs + 1;
                        LB_mat->lengths = (uint64_t *) std::realloc(LB_mat->lengths, reallocs*INIT_SEQS*sizeof(uint64_t));
                        if(LB_mat->lengths == NULL) throw "Could not reallocate sequence lengths";
                    }
                    LB_mat->lengths[n_sequences-1] = curr_len;
                    // The initial fitness maximum is assigning all tasks to just one thread
                    LB_mat->current_max += (long double) curr_len; 

            }else{
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, input);
            }
    }
    LB_mat->n = n_sequences;

    std::free(temp_seq_buffer);
}