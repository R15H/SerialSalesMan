#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#define ERROR(message) printf(stderr, "message");

// print macro for debugging when DEBUG is defined
#ifdef DEBUG
#define MESSAGE(message, ...) printf(stderr, message, __VA_ARGS__)
#else
#define MESSAGE(message, ...)
#endif

#define NR_OF_SOLUTIONS 100
struct transversal {
    int *path; // array of cities <-- can be as big as 2x the size of the number of cities
    double cost;
};

struct transversal solution; // we only store on solution, the best, which is a complete transversal



struct city {
    int id;

    struct city *cities; // index i -> city i
    double *cost;       // index i -> city i
};



void parse_inputs(int argc, char ** argv){
    if(argc < 3) {
        ERROR("Not enough arguments passed. What are you doing?")
    };
    char *cities_file = argv[1];
    int lower_bound = atoi(argv[2]);

    FILE * cities_fp = fopen(cities_file, "r");
    if (cities_fp == NULL) {
        ERROR("File not found")
        exit(-1);
    }
    char * buffer[1024];
    fgets(buffer, 1024, cities_fp) == NULL;
    int number_of_cities = atoi(strtok(buffer, " "));
    int number_of_roads = atoi(strtok(NULL, " "));
    while(fgets(buffer, 1024, cities_fp) == NULL) {
        MESSAGE("Read line: %s", buffer);

        int city_number = atoi(strtok(buffer, " "));
        int city_destination = atoi(strtok(NULL, " "));
        int city_cost = atoi(strtok(NULL, " "));
    }
    fclose(cities_fp);
}

void print_result() {

}

int main(int argc, char *argv[]) {
    double exec_time;
    parse_inputs(argc, argv);
    exec_time = -omp_get_wtime();

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(); // to the stdout!
}
