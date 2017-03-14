#include "readstream.h"

Readstream::Readstream(const char * file, void (*reading_function)(FILE * input, void * type_structure), void * type_structure){
    this->input = fopen64(file, "rt");
    if(this->input == NULL) throw "Could not open input file";
    this->reading_function = reading_function;
    this->type_structure = type_structure;
}

void Readstream::read(){
    this->reading_function(this->input, this->type_structure);
}

Readstream::~Readstream(){
    if(this->input != NULL) fclose(this->input);
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
        if(3 == sscanf(buffer, "%" PRIu64" %Le %Le", &id, &x, &y)){
            total++;
            #ifdef VERBOSE
            fprintf(stdout, "%" PRIu64", %.3Le %.3Le\n", id, x, y);
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

        if(3 == sscanf(buffer, "%" PRIu64" %Le %Le", &id, &x, &y)){
            aux[id-1] = x;
            aux2[id-1] = y;
        }
    }

    // Calculate distances
    for(uint64_t i=0;i<total;i++){
        for(uint64_t j=0;j<total;j++){
            tsp_mat->dist[i][j] = sqrtl(powl(aux[i] - aux[j], 2.0) + powl(aux2[i] - aux2[j], 2.0));
        }
    }

    std::free(aux);
    std::free(aux2);
}