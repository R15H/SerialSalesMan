#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

#define ERROR(message) printf(stderr, "message");


struct city {
    int id;
    double x;
    double y;
};

void parse_inputs(int argc, char ** argv){
    if(argc < 3) {
        ERROR("Not enough arguments passed. What are you doing?")
    };
    char *cities_file = argv[1];
    int lower_bound = atoi(argv[2]);
    // read file line by line








}

int main(int argc, char *argv[]) {
    double exec_time;
    parse_inputs(argc, argv);
    exec_time = -omp_get_wtime();

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(); // to the stdout!
}
