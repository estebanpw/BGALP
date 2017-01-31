#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "chromosome.h"
#include "common_functions.h"

int DEBUG_ACTIVE = 0;

void init_args(int argc, char ** av);

int main(int ac, char **av) {

    uint64_t n_alleles = 100;

    Chromo_rucksack<uint64_t *> * ind = new Chromo_rucksack<uint64_t *>(n_alleles);
    printf("%Le <- \n", *ind->get_fitness());


    return 0;
}


void init_args(int argc, char ** av){
    
    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--debug") == 0) DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           GAL\n");
            fprintf(stdout, "OPTIONAL:\n");
            //fprintf(stdout, "           -min_len_trimming   [Integer:   0<=X] (default 100)\n");
            fprintf(stdout, "           --debug     Turns debug on\n");
            fprintf(stdout, "           --help      Shows the help for program usage\n");
            exit(1);
        }
        /*
        if(strcmp(av[pNum], "-multifrags") == 0){
            *multifrags = fopen64(av[pNum+1], "rb");
            strncpy(path_frags, av[pNum+1], strlen(av[pNum+1]));
            path_frags[strlen(av[pNum+1])] = '\0';
            if(multifrags==NULL) terror("Could not open multifrags file");
        }
        if(strcmp(av[pNum], "-pathfiles") == 0){
            strncpy(path_files, av[pNum+1], strlen(av[pNum+1]));
            path_files[strlen(av[pNum+1])] = '\0';
        }
        if(strcmp(av[pNum], "-annotations") == 0){
            strncpy(path_annotations, av[pNum+1], strlen(av[pNum+1]));
            path_annotations[strlen(av[pNum+1])] = '\0';
        }
        if(strcmp(av[pNum], "-write_blocks_bps") == 0){
            char templine[READLINE]; templine[0] = '\0';
            strcpy(templine, av[pNum+1]);
            strcat(templine, ".blocks");
            *out_blocks = fopen64(templine, "wt");

            templine[0] = '\0';
            strcpy(templine, av[pNum+1]);
            strcat(templine, ".breakpoints");
            *out_breakpoints = fopen64(templine, "wt");
            if(*out_blocks==NULL || *out_breakpoints==NULL) terror("Could not open blocks/breakpoints file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *out_file = fopen64(av[pNum+1], "wt");
            if(out_file==NULL) terror("Could not open output file");
        }
        if(strcmp(av[pNum], "-min_len_trimming") == 0){
            *min_len_trimming = (uint64_t) atoi(av[pNum+1]);
            if(*min_len_trimming < 0) terror("Minimum trimming length must be zero or more");
        }
        if(strcmp(av[pNum], "-min_trim_itera") == 0){
            *min_trim_itera = (uint64_t) atoi(av[pNum+1]);
            if(*min_trim_itera < 0) terror("Minimum number of trimming iterations must be zero or more");
        }
        if(strcmp(av[pNum], "-hash_table_divisor") == 0){
            *ht_size = (uint64_t) atoi(av[pNum+1]);
            if(*ht_size < 1) terror("The hash table divisor must be one at least");
        }
        */
        pNum++;
    }
    /*
    if(*multifrags==NULL || *out_file==NULL || path_files[0] == '\0'){
        terror("A frags file, a path to the fasta files and an output file must be specified");
    }
    */
}

