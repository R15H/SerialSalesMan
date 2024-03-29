#include <stdio.h>
#include <omp.h>
#include <string.h>
#include "tscp.h"
#include "queue.h"

#define ERROR(message) fprintf(stderr, #message);


// print macro for debugging when DEBUG is defined
#ifdef DEBUG
#define MESSAGE(message, ...) printf(stderr, message, __VA_ARGS__)
#else
#define MESSAGE(message, ...)
#endif


void free_step(struct step_middle *step) {
    if (step == NULL) return;
    step->ref_counter--;
    if (step->ref_counter == 0) {
        struct step_middle *prev = step->previous_step;
        free(step);
        free_step(prev);
    }
}

void free_tour(struct Tour *tour) {
    if (tour == NULL) return;
    struct step_middle *step = tour->previous_step;
    free(tour);
    free_step(step);
}

/* Function to calculate x raised to the power y in O(logn)
    Time Complexity of optimized solution: O(logn)
*/
int power2(int x, unsigned int y) {
    int temp;
    if (y == 0)
        return 1;

    temp = power2(x, y / 2);
    if ((y % 2) == 0)
        return temp * temp;
    else
        return x * temp * temp;
}

double get_global_lower_bound(int number_of_cities, struct city *cities) {
    double lower_bound = 0;
    for (int i = 0; i < number_of_cities; i++) {
        lower_bound = cities[i].min_cost + cities[i].min_cost2;
    }
    return lower_bound/2;
}

double compute_updated_lower_bound(double lower_bound, int source_city,int destination_city){
    double jump_cost = cities[source_city].cost[destination_city];
    double ct = jump_cost >= cities[source_city].cost[destination_city] ? cities[source_city].cost[destination_city] : cities[source_city].cost[source_city];
    double cf = jump_cost >= cities[source_city].min_cost2 ? cities[source_city].min_cost2 :  cities[source_city].min_cost;
    return lower_bound + jump_cost - (ct+cf)/2;
}

char inline get_visited_all_cities(struct Tour *tour, struct AlgorithmState *algo_state) {
    return (tour->cities_visited & algo_state->all_cities_visited_mask) == algo_state->all_cities_visited_mask;
}

char inline get_was_visited(struct Tour *tour, int city_id) {
    return tour->cities_visited & binary_masks[city_id];
}


double get_cost_from_city_to_city(int from, int to) {
    return cities[from].cost[to]; // possible bug here...
}

struct Tour *go_to_city(struct Tour *tour, int city_id, struct AlgorithmState *algo_state) {
    // create a new tour and convert the current tour to a step middle
    struct Tour *new_tour = malloc(sizeof(union step));
    new_tour->current_city = city_id;
    new_tour->cities_visited = tour->cities_visited | binary_masks[city_id];
    new_tour->cost = compute_updated_lower_bound(tour->cost, tour->current_city, city_id);
    new_tour->previous_step = tour;
    return new_tour;
}

#define DOUBLE_MAX 1.7976931348623157e+308

void tscp(struct AlgorithmState *algo_state) {


    algo_state->solution = malloc(sizeof(struct Tour));
    algo_state->solution->cost = DOUBLE_MAX;
    algo_state->solution->cities_visited = algo_state->all_cities_visited_mask;
    algo_state->solution->current_city = 0;
    algo_state->solution->previous_step = NULL;

    struct Tour *first_step = malloc(sizeof(union step)); // calloc to initialize all values to 0
    first_step->current_city = 0;
    first_step->cities_visited = 1;
    first_step->cost = get_global_lower_bound(algo_state->number_of_cities, cities);
    first_step->previous_step = NULL;

    queue_push(algo_state->queue, first_step);

    struct Tour *current_tour;
    while ((current_tour = queue_pop(algo_state->queue))) {

        char step_exceeds_max_cost = current_tour->cost > algo_state->max_lower_bound;
        char step_worst_than_found_solution = current_tour->cost > algo_state->solution->cost;
        if (step_exceeds_max_cost ||
            step_worst_than_found_solution) {// we let the compiler short circuit the expression, for we trust in the compiler!
            free_tour(current_tour);
            continue;
        }
        char returned_to_start = current_tour->current_city == 0;
        if (get_visited_all_cities(current_tour, algo_state) && returned_to_start) {
            char current_tour_is_better = current_tour->cost < algo_state->solution->cost;
            if (current_tour_is_better) {
                free_tour(algo_state->solution);
                algo_state->solution = current_tour;
            } else {
                free_tour(current_tour);
            }

            continue;
        }

        // for each neighbor of the current city that is not in the current path


        //struct Tour *new_tours[algo_state->number_of_cities];
        struct Tour *newTour;
        if (get_visited_all_cities(current_tour,
                                   algo_state)) { // if we visited all cities, we must return to the start city (0)
            newTour = go_to_city(current_tour, 0, algo_state);
            queue_push(algo_state->queue, newTour);
        } else
            for (int i = 0; i < algo_state->number_of_cities; i++) {
                if (current_tour->current_city == i || get_was_visited(current_tour, i)) continue;

                newTour = go_to_city(current_tour, i, algo_state);
                queue_push(algo_state->queue, newTour);


                //new_tours[i] = new_tour;
            }
        // convert current tour to a step middle
        struct step_middle *new_step_middle = (struct step_middle *) current_tour;
        new_step_middle->ref_counter = 1000000000000000; // -> its equal to the number of tours added!! algo_state->number_of_cities;  // iterations of i
        // create mutex <-- for part2

        // push all new tours to the queue
        // This must be done only after the step_middle is created, otherwise another thread may grab the one of the new tours, try to free it and end up accessing the current tour as a step_middle!
        //for (int i = 0; i < algo_state->number_of_cities; i++) {
        //    queue_push(algo_state->queue, new_tours[i]);
        // }

    }

}


#define CITY_POINTER_SIZE sizeof(int)
// define a city pointer type


void place_cost_in_city(int city_source, int city_destination, double cost, int number_of_cities) {
    struct city *current_city = &cities[city_source];
    if (current_city->cost == NULL) {
        current_city->cost = calloc(number_of_cities, sizeof(double));
    }
    (current_city->cost)[city_destination] = cost;
}


// Initializes the AlgorithmState variable
void parse_inputs(int argc, char **argv, struct AlgorithmState *algo_state) {
    printf("Parsing inputs...\n");
    if (argc < 3) {
        ERROR("Not enough arguments passed. What are you doing?\n")
        exit(-1);
    };
    char *cities_file = argv[1];
    algo_state->max_lower_bound = atoi(argv[2]);

    char buffer[1024];
    FILE *cities_fp = fopen(cities_file, "r");

    if (cities_fp == NULL) {
        ERROR(File not found)
        exit(-1);
    }
    fgets((char *) &buffer, 1024, cities_fp);
    algo_state->number_of_cities = atoi(strtok(buffer, " "));
    cities = calloc(algo_state->number_of_cities, sizeof(struct city));
    algo_state->number_of_roads = atoi(strtok(NULL, " "));

    all_cities_visited_mask = 0;
    for (int i = 0; i < algo_state->number_of_cities; ++i) {
        binary_masks[i] = power2(2, i);
        all_cities_visited_mask += binary_masks[i];
    }
    algo_state->all_cities_visited_mask = all_cities_visited_mask;

    while (fgets(buffer, 1024, cities_fp) != NULL) {
        MESSAGE("Read line: %s", buffer);

        int city_number = atoi(strtok(buffer, " "));
        int city_destination = atoi(strtok(NULL, " "));
        double city_cost = strtod(strtok(NULL, " "), NULL);
        struct city current_city = cities[city_number];
        current_city.id = city_number;
        place_cost_in_city(city_number, city_destination, city_cost, algo_state->number_of_cities);
        place_cost_in_city(city_destination, city_number, city_cost, algo_state->number_of_cities);
    }
    fclose(cities_fp);
}

void print_result(struct AlgorithmState *algo_state) {
    if (algo_state->solution == NULL) {
        printf("NO SOLUTION");
        return;
    }
    printf("%.1f\n", algo_state->solution->cost);

    struct step_middle *step = (struct step_middle *) algo_state->solution; // this is actually a Tour but its okay
    int cities_visited = algo_state->number_of_cities+1;
    char pathTaken[cities_visited*2+1];
    pathTaken[cities_visited*2] = '\0';

    for(int i = cities_visited-1; i >= 0; i--){

        pathTaken[i*2+1] = ' ';
        pathTaken[i*2] = (char)(step->current_city + '0');
        step = step->previous_step;
    }



    printf("%s\n", pathTaken);
}

void dealloc_data() {
    free(cities);
}

char compare_paths(struct Tour *a, struct Tour *b) {
    // resolves which path is better (and therefore will have the higher priority)
    // Paths with smaller costs are better
    if (a->cost > b->cost) return '\0';
    return '\1';
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
