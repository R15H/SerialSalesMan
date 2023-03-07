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

void free_step(struct step_middle *step)
{
    if (step == NULL)
        return;
    step->ref_counter--;
    if (step->ref_counter == 0)
    {
        struct step_middle *prev = step->previous_step;
        free(step);
        free_step(prev);
    }
}

void free_tour(struct Tour *tour)
{
    if (tour == NULL)
        return;
    struct step_middle *step = tour->previous_step;
    free(tour);
    free_step(step);
}

/* Function to calculate x raised to the power y in O(logn)
    Time Complexity of optimized solution: O(logn)
*/
int power2(int x, unsigned int y)
{
    int temp;
    if (y == 0)
        return 1;

    temp = power2(x, y / 2);
    if ((y % 2) == 0)
        return temp * temp;
    else
        return x * temp * temp;
}

char inline get_visited_all_cities(struct Tour *tour, struct AlgorithmState *algo_state)
{
    return tour->cities_visited & algo_state->all_cities_visited_mask;
}

char inline get_was_visited(struct Tour *tour, int city_id)
{
    return tour->cities_visited & binary_masks[city_id];
}

double get_cost_from_city_to_city(int from, int to)
{
    return (*cities[from].cost)[to]; // possible bug here...
}

void tscp(struct AlgorithmState *algo_state)
{

    struct Tour *first_step = calloc(1, sizeof(union step)); // calloc to initialize all values to 0
    queue_push(algo_state->queue, first_step);

    algo_state->solution = calloc(1, sizeof(struct Tour));
    struct Tour *solution = algo_state->solution;

    struct Tour *current_tour;
    while ((current_tour = queue_pop(algo_state->queue)))
    {

        char step_exceeds_max_cost = current_tour->cost > algo_state->max_lower_bound;
        char step_worst_than_found_solution = current_tour->cost > solution->cost;
        if (step_exceeds_max_cost ||
            step_worst_than_found_solution)
        { // we let the compiler short circuit the expression, for we trust in the compiler!
            free_tour(current_tour);
            continue;
        }
        if (get_visited_all_cities(current_tour, algo_state))
        {
            char current_tour_is_better = current_tour->cost < solution->cost;
            if (current_tour_is_better)
            {
                free_tour(solution);
                solution = current_tour;
            }
            else
            {
                free_tour(current_tour);
            }

            continue;
        }

        // for each neighbor of the current city that is not in the current path

        struct Tour *new_tours[algo_state->number_of_cities];
        for (int i = 0; i < algo_state->number_of_cities; i++)
        {
            if (current_tour->current_city == i)
                continue; // skip the current city
            if (get_was_visited(current_tour, i))
                continue;

            // create a new tour and convert the current tour to a step middle
            struct Tour *new_tour = malloc(sizeof(union step)); // malloc because we initialize the values
            new_tour->previous_step = current_tour->previous_step;
            new_tour->current_city = i;
            new_tour->cities_visited = current_tour->cities_visited | binary_masks[i];
            new_tour->cost = current_tour->cost + get_cost_from_city_to_city(current_tour->current_city, i);

            new_tours[i] = new_tour;
        }
        // convert current tour to a step middle
        struct step_middle *new_step_middle = (struct step_middle *)current_tour;
        new_step_middle->ref_counter = algo_state->number_of_cities; // iterations of i
        // create mutex <-- for part2

        // push all new tours to the queue
        // This must be done only after the step_middle is created, otherwise another thread may grab the one of the new tours, try to free it and end up accessing the current tour as a step_middle!
        for (int i = 0; i < algo_state->number_of_cities; i++)
        {
            queue_push(algo_state->queue, new_tours[i]);
        }
    }
}

#define CITY_POINTER_SIZE sizeof(int)
// define a city pointer type

void place_cost_in_city(int city_source, int city_destination, double cost, int number_of_cities)
{
    struct city current_city = cities[city_source];
    printf("place_cost_debug0...\n");
    if (current_city.cities == NULL)
    {
        current_city.cities = calloc(number_of_cities, CITY_POINTER_SIZE);
        current_city.cost = calloc(number_of_cities, sizeof(double));
    }
    printf("%d %d %d \n", current_city.cities, *current_city.cities, 20);
    printf("%d \n", current_city.cities[0]);
    current_city.cities[0] = 10;
    printf("%d %d %d \n", current_city.cities, *current_city.cities, 10);
    printf("%d \n", current_city.cities[0]);
    printf("place_cost_debug1...\n");
    current_city.cities[city_destination] = city_destination;
    printf("place_cost_debug2...\n");
    (*current_city.cost)[city_destination] = cost;
    printf("place_cost_debug3...\n");
    (*(current_city.cost))[city_destination] = cost;
}

// Initializes the AlgorithmState variable
void parse_inputs(int argc, char **argv, struct AlgorithmState *algo_state)
{
    printf("Parsing inputs...\n");
    if (argc < 3)
    {
        ERROR("Not enough arguments passed. What are you doing?\n")
        exit(-1);
    };
    char *cities_file = argv[1];
    algo_state->max_lower_bound = atoi(argv[2]);

    char buffer[1024];
    FILE *cities_fp = fopen(cities_file, "r");
    printf("debug1...\n");
    if (cities_fp == NULL)
    {
        ERROR(File not found)
        exit(-1);
    }
    printf("debug2...\n");
    fgets((char *)&buffer, 1024, cities_fp);
    algo_state->number_of_cities = atoi(strtok(buffer, " "));
    cities = calloc(algo_state->number_of_cities, sizeof(struct city));
    algo_state->number_of_roads = atoi(strtok(NULL, " "));
    printf("debug3...\n");
    all_cities_visited_mask = 0;
    for (int i = 0; i < algo_state->number_of_cities; ++i)
    {
        binary_masks[i] = power2(2, i);
        all_cities_visited_mask += binary_masks[i];
    }
    algo_state->all_cities_visited_mask = all_cities_visited_mask;
    printf("debug4...\n");

    while (fgets(buffer, 1024, cities_fp) != NULL)
    {
        printf("debug5...\n");
        printf("debug5...\n");
        MESSAGE("Read line: %s", buffer);
        printf("debug4...\n");

        int city_number = atoi(strtok(buffer, " "));
        int city_destination = atoi(strtok(NULL, " "));
        double city_cost = strtod(strtok(NULL, " "), NULL);
        struct city current_city = cities[city_number];
        current_city.id = city_number;
        printf("debug4...\n");
        place_cost_in_city(city_number, city_destination, city_cost, algo_state->number_of_cities);
        printf("debug4...\n");
        place_cost_in_city(city_destination, city_number, city_cost, algo_state->number_of_cities);
        printf("debug4...\n");
        printf("debug6...\n");
    }
    fclose(cities_fp);
    printf("debug7\n");
}

void print_result(struct AlgorithmState *algo_state)
{
    if (algo_state->solution == NULL)
    {
        printf("NO SOLUTION");
        return;
    }
    printf("%.1f\n", algo_state->solution->cost);

    struct step_middle *step = (struct step_middle *)algo_state->solution; // this is actually a Tour but its okay
    do
    {
        printf("%d ", step->current_city);
    } while ((step = step->previous_step) != NULL);
}

void dealloc_data()
{
    free(cities);
}

int compare_paths(struct Tour *a, struct Tour *b)
{
    // resolves which path is better (and therefore will have the higher priority)
    // Paths with smaller costs are better
    return a->cost - b->cost;
}

// continuo só de paths --> uma pagina na memoria só vai ter paths e nao vai ter outros lixos (quando se tiver a ver um path)
// steps do mesmo path vao estar juntos na memoria a principio (pk só troca quando há um outro de maior prioridade)
// trocar entre paths é menos dispendioso pk uma pagina tem os paths de varios caminhos
// pode ser preciso ter bues paginas para correr o algo uma x ....

int main(int argc, char *argv[])
{
    double exec_time;
    struct AlgorithmState algo_state;
    parse_inputs(argc, argv, &algo_state);
    exec_time = -omp_get_wtime();

    algo_state.queue = queue_create((char (*)(void *, void *))compare_paths);
    tscp(&algo_state);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(&algo_state); // to the stdout!

    queue_delete(algo_state.queue);
    dealloc_data();
}
