#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <assert.h>
#include "tscp.h"
#include "queue.h"
#include <stdlib.h>

#define ERROR(message) fprintf(stderr, #message);

#ifdef DEBUG
#define MESSAGE(message, ...) printf(stderr, message, __VA_ARGS__)
#else

#define MESSAGE(message, ...)
#endif

#define DOUBLE_MAX 1.7976931348623155e+308

#pragma runtime_checks("", off)

union step *get_clean_step() {
    return malloc(sizeof(union step));
}


int get_visited_all_cities(struct Tour *tour, struct AlgorithmState *algo_state) {
    return (tour->cities_visited & algo_state->all_cities_visited_mask) == algo_state->all_cities_visited_mask;
}

unsigned int get_was_visited(struct Tour *tour, int city_id) {
    return tour->cities_visited & binary_masks[city_id];
}

double get_cost_from_city_to_city(unsigned int from, unsigned int to) {
    return cities[from].cost[to];
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

void print_result(struct AlgorithmState *algo_state) {
    if (algo_state->solution == NULL) {
        printf("NO SOLUTION");
        return;
    }

    struct step_middle *step = (struct step_middle *) algo_state->solution; // this is actually a Tour but its okay
    int cities_visited = algo_state->number_of_cities + 1;
    int values[50]; // max 50 cities

    printf("%.1f\n", algo_state->solution->cost);

    for (int i = cities_visited - 1; i >= 0; i--) {
        values[i] = step->current_city;
        step = step->previous_step;
        if (step == NULL) break;
    }
    for (int i = 0; i < cities_visited - 1; i++) printf("%d ", values[i]);
    printf("0 \n");
}

void free_step(struct step_middle *step) {
    if (step == NULL)return;
    step->ref_counter--;
    if (step->ref_counter <= 0) {
        free_step(step->previous_step);
        free(step);
    }
}

void free_tour(struct Tour *tour) {
    if (tour == NULL) return;
    struct step_middle *step = tour->previous_step;
    free_step(step);
    free(tour);
}

double get_global_lower_bound(int number_of_cities, struct city *cities) {
    double lower_bound = 0;
    for (int i = 0; i < number_of_cities; i++) lower_bound = cities[i].min_cost + cities[i].min_cost2;
    return lower_bound / 2;
}

inline double compute_updated_lower_bound(double lower_bound, unsigned int source_city, unsigned int destination_city) {
    double jump_cost = get_cost_from_city_to_city(source_city, destination_city);
    int comp_1 = jump_cost >= cities[source_city].min_cost2;
    double ct = (comp_1) * cities[source_city].min_cost2
                + (comp_1 == 0) * cities[source_city].min_cost;
    int comp_2 = jump_cost >= cities[source_city].min_cost2;
    double cf =
            (comp_2) * cities[destination_city].min_cost2 +
            (comp_2 == 0) * cities[destination_city].min_cost;
    return lower_bound + jump_cost - (ct + cf) / 2;
}

inline struct Tour *go_to_city(struct Tour *tour, short city_id, struct AlgorithmState *algo_state, double cost) {
    // create a new tour and convert the current tour to a step middle
    struct Tour *new_tour = (struct Tour *) get_clean_step();
    new_tour->current_city = city_id;
    new_tour->cities_visited = tour->cities_visited | binary_masks[city_id];
    new_tour->cost = cost;
    new_tour->previous_step = (struct step_middle *) tour;
    return new_tour;
}

int shouldnt_create_tour(double cost, int current_city, struct AlgorithmState *algo_state) {
    int step_exceeds_max_cost = cost > algo_state->max_lower_bound;
    int step_worst_than_found_solution = cost > algo_state->solution->cost;
    return step_exceeds_max_cost ||
           step_worst_than_found_solution;
}

void visit_city(struct Tour *tour,int destination, struct AlgorithmState *algo_state, &tours_created){
    if (!get_was_visited(tour, i)) {
        double new_cost = compute_updated_lower_bound(tour->cost, tour->current_city, i);
        int discard_tour = shouldnt_create_tour(new_cost, i, algo_state);

        if (!discard_tour) {
            struct Tour *new_tour = go_to_city(tour, i, algo_state, new_cost);
            int finished = get_visited_all_cities(new_tour, algo_state);
            if (!finished) {
                tours_created++;
                queue_push(algo_state->queue, new_tour);
            } else {
                ((struct step_middle *) new_tour)->ref_counter = 1;
                double final_cost = compute_updated_lower_bound(new_tour->cost, i, 0);
                int current_tour_is_better = final_cost < algo_state->solution->cost;
                if (current_tour_is_better) {
                    free_tour(algo_state->solution);
                    algo_state->solution = go_to_city(new_tour, 0, algo_state, final_cost);
                    queue_trim(algo_state->queue, final_cost);
                }
            }
        }
    }
}

int analyseTour(struct Tour *tour, struct AlgorithmState *algo_state) {
    int tours_created = 0;
    int loops = algo_state->number_of_cities - 1;
    struct Tour tours_to_free[64];
    int free_tours = 0;
    int i = 0;
    for (; i < loops; i += 2) {


        if (get_was_visited(tour, i + 1)) continue;
        double new_cost = compute_updated_lower_bound(tour->cost, tour->current_city, i + 1);

        int discard_tour = shouldnt_create_tour(new_cost, i + 1, algo_state);
        if (!discard_tour) {
            struct Tour *new_tour = go_to_city(tour, i + 1, algo_state, new_cost);
            int finished = get_visited_all_cities(new_tour, algo_state);
            if (finished) {
                ((struct step_middle *) new_tour)->ref_counter = 1;
                double final_cost = compute_updated_lower_bound(new_tour->cost, i + 1, 0);
                int current_tour_is_better = final_cost < algo_state->solution->cost;
                if (current_tour_is_better) {
                    free_tour(algo_state->solution);
                    algo_state->solution = go_to_city(new_tour, 0, algo_state, final_cost);
                    queue_trim(algo_state->queue, final_cost);
                }
                continue;
            }
            tours_created++;
            queue_push(algo_state->queue, new_tour);
        }
    }
    if (algo_state->number_of_cities % 2 != 0) {
        i = algo_state->number_of_cities - 1;
        struct Tour *new_tour;
        if (!get_was_visited(tour, i)) {
            double new_cost = compute_updated_lower_bound(tour->cost, tour->current_city, i);

            int discard_tour = shouldnt_create_tour(new_cost, i, algo_state);
            if (discard_tour) return tours_created;
            new_tour = go_to_city(tour, i, algo_state, new_cost);
            int finished = get_visited_all_cities(new_tour, algo_state);
            if (finished) {
                ((struct step_middle *) new_tour)->ref_counter = 1;
                double final_cost = compute_updated_lower_bound(new_tour->cost, i, 0);

                int current_tour_is_better = final_cost < algo_state->solution->cost;
                if (current_tour_is_better) {
                    free_tour(algo_state->solution);
                    algo_state->solution = go_to_city(new_tour, 0, algo_state, final_cost);
                    queue_trim(algo_state->queue, final_cost);
                }
                return tours_created;
            }
            tours_created++;
            queue_push(algo_state->queue, new_tour);
        }
    }

    return tours_created;
}


void tscp(struct AlgorithmState *algo_state) {
    algo_state->solution = (struct Tour *) get_clean_step();
    algo_state->solution->cost = 100000000;
    algo_state->solution->cities_visited = algo_state->all_cities_visited_mask;
    algo_state->solution->current_city = 0;
    algo_state->solution->previous_step = NULL;

    struct Tour *first_step = (struct Tour *) get_clean_step();
    first_step->current_city = 0;
    first_step->cities_visited = 1;
    first_step->cost = get_global_lower_bound(algo_state->number_of_cities, cities);
    first_step->previous_step = NULL;

    queue_push(algo_state->queue, first_step);

    struct Tour *current_tour;
    while ((current_tour = queue_pop(algo_state->queue))) {
        int newToursCreated = analyseTour(current_tour, algo_state);
        if (newToursCreated == 0) {
            free_tour(current_tour);
            continue;
        }
        ((struct step_middle *) current_tour)->ref_counter = newToursCreated;
    }
}


void place_cost_in_city(int city_source, int city_destination, double cost, int number_of_cities) {
    struct city *current_city = &cities[city_source];
    if (current_city->cost == NULL) {
        current_city->cost = calloc(number_of_cities, sizeof(double));
        // initialize every entry of the array to DOUBLE_MAX/2
        for (int i = 0; i < number_of_cities; i++) {
            (current_city->cost)[i] = DOUBLE_MAX / 4;
        }
    }
    if (current_city->min_cost > cost)
        current_city->min_cost = cost;
    else if (current_city->min_cost2 > cost)
        current_city->min_cost2 = cost;
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
    algo_state->max_lower_bound = atof(argv[2]);

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
        int city_number = atoi(strtok(buffer, " "));
        int city_destination = atoi(strtok(NULL, " "));
        double city_cost = strtod(strtok(NULL, " "), NULL);
        place_cost_in_city(city_number, city_destination, city_cost, algo_state->number_of_cities);
        place_cost_in_city(city_destination, city_number, city_cost, algo_state->number_of_cities);
    }
    fclose(cities_fp);
}

void dealloc_data() {
    free(cities);
}


int main(int argc, char *argv[]) {
    double exec_time;
    struct AlgorithmState algo_state;
    parse_inputs(argc, argv, &algo_state);
    exec_time = -omp_get_wtime();
    algo_state.queue = queue_create((char (*)(void *, void *)) NULL);
    tscp(&algo_state);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(&algo_state);
    queue_delete(algo_state.queue);
    dealloc_data();
}
