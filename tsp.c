#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <stdbool.h>
#include "tscp.h"
#include "queue.h"

#define ERROR(message) printf(stderr, "message");

// print macro for debugging when DEBUG is defined
#ifdef DEBUG
#define MESSAGE(message, ...) printf(stderr, message, __VA_ARGS__)
#else
#define MESSAGE(message, ...)
#endif



void tscp(struct AlgorithmState *algo_state){


    struct transversal_step *step  = calloc(1,sizeof(struct transversal_step)); // calloc to initialize all values to 0
    step->current_city = algo_state->cities; // start from the first city, same as step->current_city = &AlgorithmState->cities[0];
    queue_push(algo_state->queue, step);

    struct transversal_step *solution;

    struct transversal_step *current_step;
    while(current_step = queue_pop(algo_state->queue)){

        int step_exceeds_max_cost = current_step->cost > algo_state->max_lower_bound;
        if(step_exceeds_max_cost) {
            free(current_step);
            continue;
        }
        int step_worst_than_found_solution  = current_step->cost > solution->cost;
        if(step_worst_than_found_solution) {
            // TODO free all memory except the solution
            continue;
        }
        int visited_all_cities = algo_state->number_of_cities == (current_step->nr_cities_visited+1);
        if(visited_all_cities) {
            if(current_step->cost < solution->cost) {
                solution = current_step;
                // TODO free solution
            }
            continue;
        }
        // for each neighbor of the current city that is not in the current path
        for(int i = 0; i < current_step->current_city->nr_cities; i++){

            struct city* found_city =  bsearch(current_step->current_city->cities[i], current_step->path, current_step->nr_cities_visited, sizeof(struct city *), compare_cities);
            int city_is_in_path = found_city != NULL;
            if(city_is_in_path) continue;
            int new_cost = found_city->cost[i] + current_step->cost;


            int worst_than_found_solution = new_cost > solution->cost;
            int exceeds_max_cost = new_cost > algo_state->max_lower_bound;
            if(worst_than_found_solution || exceeds_max_cost) continue; // we let the compiler short circuit the expression, for we trust in the compiler!

            // create a new step
            struct transversal_step *step  = calloc(1,sizeof(struct transversal_step));
            step->current_city = current_step->current_city->cities[i];
            current_step->current_city->cities;
            queue_push(algo_state->queue, step);
        }
    }






}





int number_of_cities;
int number_of_roads;
struct city *cities;


// Initializes the AlgorithmState variable
void parse_inputs(int argc, char ** argv, struct AlgorithmState *algo_state){
    if(argc < 3) {
        ERROR("Not enough arguments passed. What are you doing?")
    };
    char *cities_file = argv[1];
    algo_state->max_lower_bound = atoi(argv[2]);

    char * buffer[1024];
    FILE * cities_fp = fopen(cities_file, "r");
    if (cities_fp == NULL) {
        ERROR("File not found")
        exit(-1);
    }
    fgets(buffer, 1024, cities_fp);
    algo_state->number_of_cities = atoi(strtok(buffer, " "));
    algo_state->cities = malloc(sizeof(struct city) * number_of_cities);
    algo_state->number_of_roads = atoi(strtok(NULL, " "));

    while(fgets(buffer, 1024, cities_fp) == NULL) {
        MESSAGE("Read line: %s", buffer);

        int city_number = atoi(strtok(buffer, " "));
        int city_destination = atoi(strtok(NULL, " "));
        int city_cost = atoi(strtok(NULL, " "));
    }
    fclose(cities_fp);
}

void print_result(struct AlgorithmState *algo_state) {
    if(algo_state->solution == NULL) {
        printf("NO SOLUTION");
        return;
    }
    printf("%.1f\n", algo_state->solution->cost);

    struct transversal_step *step = algo_state->solution;
    struct city *city;
    while( (city = step->current_city) != NULL) {
        printf("%d ", city->id);
        step = step->previous_step;
    }
}

void dealloc_data() {
    free(cities);
}

int compare_paths(struct transversal_step *a, struct transversal_step *b) {
    // resolves which path is better (and therefore will have the higher priority)
    // Paths with smaller costs are better
    return a->cost - b->cost;
}

// continuo só de paths --> uma pagina na memoria só vai ter paths e nao vai ter outros lixos (quando se tiver a ver um path)
// steps do mesmo path vao estar juntos na memoria a principio (pk só troca quando há um outro de maior prioridade)
// trocar entre paths é menos dispendioso pk uma pagina tem os paths de varios caminhos
// pode ser preciso ter bues paginas para correr o algo uma x ....

int main(int argc, char *argv[]) {
    double exec_time;
    struct AlgorithmState algo_state;
    parse_inputs(argc, argv, &algo_state);
    exec_time = -omp_get_wtime();


    algo_state.queue = queue_create((char (*)(void *, void *)) compare_paths);
    tscp(&algo_state);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(&algo_state); // to the stdout!

    queue_delete(algo_state.queue);
    dealloc_data();
}
