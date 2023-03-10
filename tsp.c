#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <assert.h>
#include "tscp.h"
#include "queue.h"

#define ERROR(message) fprintf(stderr, #message);

// print macro for debugging when DEBUG is defined
#ifdef DEBUG
#define MESSAGE(message, ...) printf(stderr, message, __VA_ARGS__)
#else
#define MESSAGE(message, ...)
#endif

#define DOUBLE_MAX 1.7976931348623155e+308
#define INT_MAX 2147483647

#pragma runtime_checks("", off)
union step *steps;
int *free_steps;
int current_free_step = -1;

union step *get_clean_step()
{
    assert(current_free_step > 0);
    return &steps[free_steps[current_free_step--]];
}

void free_union_step(int stepID)
{
    current_free_step++;
    free_steps[current_free_step] = stepID;
}

// Tour interface
int get_visited_all_cities(struct Tour *tour, struct AlgorithmState *algo_state)
{
    // printf("%d\n",algo_state->all_cities_visited_mask);
    // printf("%d\n",tour->cities_visited);
    return (tour->cities_visited & algo_state->all_cities_visited_mask) == algo_state->all_cities_visited_mask;
}

// Cities interface

unsigned int get_was_visited(struct Tour *tour, int city_id)
{
    // if(cities[tour->current_city].cost[city_id] == DOUBLE_MAX/4) return 1; // if the cost is infinity then the city was visited
    return tour->cities_visited & binary_masks[city_id];
}

double get_cost_from_city_to_city(unsigned int from, unsigned int to)
{
    return cities[from].cost[to];
}

void print_result(struct AlgorithmState *algo_state)
{
    if (algo_state->solution == NULL)
    {
        printf("NO SOLUTION");
        return;
    }

    struct step_middle *step = (struct step_middle *)algo_state->solution; // this is actually a Tour but its okay
    int cities_visited = algo_state->number_of_cities + 1;
    char pathTaken[30 * 2 + 1];
    printf("%.1f\n",algo_state->solution->cost);

    for (int i = cities_visited - 1; i >= 0; i--)
    {
        // printf("Printing inside\n");
        if (step->previous_step == NULL)
            break;
        printf("%d ", step->current_city);
        step = step->previous_step;
    }
    // printf("Printing\n");

    for (int i = cities_visited - 1; i >= 0; i--)
    {

        pathTaken[i * 2 + 1] = ' ';
        pathTaken[i * 2] = (char)(step->current_city + '0');
        if (step->previous_step == NULL)
        {
            // printf("Unexpected null step!\n");
            break;
        }
        step = step->previous_step;
    }
    pathTaken[cities_visited * 2] = '\0';
    // printf("%s\n", pathTaken);
    printf("0 \n");
}

void free_step(struct step_middle *step)
{
    if (step == NULL)
        return;
    // omp_set_lock(&step->decrease_counter_lock);
    step->ref_counter--;
    if (step->ref_counter == 0)
    { // WARNING: HAZARDOUS CODE WHICH RELIES ON SEMANTICS
        // printf("Freeing step_middle...");
        free_step(step->previous_step);
        free_union_step(step->stepID);

        // omp_unset_lock(&step->decrease_counter_lock); // por causa da nossa boa consciencia
        // omp_destroy_lock(&step->decrease_counter_lock);
        // free(step);
    }
    // omp_unset_lock(&step->decrease_counter_lock);
}

void free_tour(struct Tour *tour)
{
    // printf("Freeing tour...\n");
    if (tour == NULL)
        return;
    struct step_middle *step = tour->previous_step;
    free_step(step);
    free_union_step(tour->stepID);
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

double get_global_lower_bound(int number_of_cities, struct city *cities)
{
    double lower_bound = 0;
    for (int i = 0; i < number_of_cities; i++)
    {
        lower_bound = cities[i].min_cost + cities[i].min_cost2;
    }
    return lower_bound / 2;
}

double compute_updated_lower_bound(double lower_bound, unsigned int source_city, unsigned int destination_city)
{
    double jump_cost = get_cost_from_city_to_city(source_city, destination_city);
    double ct = jump_cost >= cities[source_city].min_cost2 ? cities[source_city].min_cost2
                                                           : cities[source_city].min_cost;
    double cf =
        jump_cost >= cities[destination_city].min_cost2 ? cities[destination_city].min_cost2 : cities[destination_city].min_cost;
    return lower_bound + jump_cost - (ct + cf) / 2;
}

struct Tour *go_to_city(struct Tour *tour, int city_id, struct AlgorithmState *algo_state)
{
    // create a new tour and convert the current tour to a step middle
    struct Tour *new_tour = (struct Tour *)get_clean_step();
    new_tour->current_city = city_id;
    new_tour->cities_visited = tour->cities_visited | binary_masks[city_id];
    new_tour->cost = compute_updated_lower_bound(tour->cost, tour->current_city, city_id);
    new_tour->previous_step = (struct step_middle *)tour;
    return new_tour;
}

void print_tour(struct Tour *tour)
{
    static int count = 0;
    if (tour == NULL)
    {
        printf("NULL");
        return;
    }
    printf("Tour nr %d: current_city: %d, cities_visited: %d, cost: %f, previous_step: %p\n", count++, tour->current_city,
           tour->cities_visited, tour->cost, tour->previous_step);
}

void print_step_middle(struct step_middle *step)
{
    if (step == NULL)
    {
        printf("NULL");
        return;
    }
    printf("Step: current_city: %d, previous_step: %p, ref_counter: %d\n",
           step->current_city,
           step->previous_step, step->ref_counter);
}

void print_algo_state(struct AlgorithmState *algo)
{
    fprintf(stderr, "Algorithm state: \n");
    fprintf(stderr, "number_of_cities: %d\n", algo->number_of_cities);
    fprintf(stderr, "all_cities_visited_mask: %d\n", algo->all_cities_visited_mask);
    fprintf(stderr, "max_lower_bound: %f\n", algo->max_lower_bound);
    fprintf(stderr, "Current best solution: ");
    print_tour(algo->solution);
    print_result(algo);
    fprintf(stderr, "\n Queue:");
    queue_print(algo->queue, stderr, (void *)print_step_middle);
}

void tscp(struct AlgorithmState *algo_state)
{

    algo_state->solution = (struct Tour *)get_clean_step();
    algo_state->solution->cost = 100000000; // DOUBLE_MAX;
    algo_state->solution->cities_visited = algo_state->all_cities_visited_mask;
    algo_state->solution->current_city = 0;
    algo_state->solution->previous_step = NULL;

    struct Tour *first_step = (struct Tour *)get_clean_step();
    first_step->current_city = 0;
    first_step->cities_visited = 1;
    first_step->cost = get_global_lower_bound(algo_state->number_of_cities, cities);
    first_step->previous_step = NULL;

    queue_push(algo_state->queue, first_step);

    struct Tour *current_tour;
    while ((current_tour = queue_pop(algo_state->queue)))
    {
        // print_algo_state(algo_state);
        // print_tour(current_tour);
        int step_exceeds_max_cost = current_tour->cost > algo_state->max_lower_bound;
        int step_worst_than_found_solution = current_tour->cost > algo_state->solution->cost;
        if (step_exceeds_max_cost ||
            step_worst_than_found_solution)
        { // we let the compiler short circuit the expression, for we trust in the compiler!
            free_tour(current_tour);
            continue;
        }
        int returned_to_start = current_tour->current_city == 0;
        if (get_visited_all_cities(current_tour, algo_state) && returned_to_start)
        {
            // print_algo_state(algo_state);
            int current_tour_is_better = current_tour->cost < algo_state->solution->cost;
            if (current_tour_is_better)
            {
                free_tour(algo_state->solution);
                algo_state->solution = current_tour;
                // print_tour(current_tour);
                // print_algo_state(algo_state);
            }
            else
                free_tour(current_tour);

            continue;
        }

        struct Tour *newTour;
        int newToursCreated = 0;
        int loops = algo_state->number_of_cities;

        if (get_visited_all_cities(current_tour, algo_state) && !returned_to_start)
        {
            newTour = go_to_city(current_tour, 0, algo_state);
            queue_push(algo_state->queue, newTour);
            newToursCreated = 1;
            // printf("returned to start\n");
        }
        else
        {
            int i = 0;
            for (; i < loops-1; i += 2)
            {
                int should_go_to_city = !(current_tour->current_city == i || get_was_visited(current_tour, i));
                if (should_go_to_city)
                {
                    newTour = go_to_city(current_tour, i, algo_state);
                    queue_push(algo_state->queue, newTour);
                    newToursCreated++;
                }
                // printf("Evaling city %d and %d\n", i, i+1);

                int should_go_to_city2 = !(current_tour->current_city == i + 1 || get_was_visited(current_tour, i + 1));
                if (should_go_to_city2)
                {
                    newTour = go_to_city(current_tour, i + 1, algo_state);
                    queue_push(algo_state->queue, newTour);
                    newToursCreated++;
                }
            }
            // print_tour(current_tour);
            //  loop one by one with the remainder of the cities
            if (0 != algo_state->number_of_cities % 2)
            {
                int should_go_to_city = !(current_tour->current_city == algo_state->number_of_cities - 1 || get_was_visited(current_tour, algo_state->number_of_cities - 1));
                if (should_go_to_city)
                {
                    newTour = go_to_city(current_tour, algo_state->number_of_cities - 1, algo_state);
                    queue_push(algo_state->queue, newTour);
                    newToursCreated++;
                }
            }
        }

        // convert this tour to a step middle
        ((struct step_middle *)current_tour)->ref_counter = newToursCreated;
        //    omp_init_lock(&((struct step_middle *) current_tour)->decrease_counter_lock);

        // push all new tours to the queue that this thread found
    }
}

#define CITY_POINTER_SIZE sizeof(int)
// define a city pointer type

void place_cost_in_city(int city_source, int city_destination, double cost, int number_of_cities)
{
    struct city *current_city = &cities[city_source];
    if (current_city->cost == NULL)
    {
        current_city->cost = calloc(number_of_cities, sizeof(double));
        // initialize every entry of the array to DOUBLE_MAX/2
        for (int i = 0; i < number_of_cities; i++)
        {
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
void parse_inputs(int argc, char **argv, struct AlgorithmState *algo_state)
{
    printf("Parsing inputs...\n");
    if (argc < 3)
    {
        ERROR("Not enough arguments passed. What are you doing?\n")
        exit(-1);
    };
    char *cities_file = argv[1];
    algo_state->max_lower_bound = atof(argv[2]);

    char buffer[1024];
    FILE *cities_fp = fopen(cities_file, "r");

    if (cities_fp == NULL)
    {
        ERROR(File not found)
        exit(-1);
    }
    fgets((char *)&buffer, 1024, cities_fp);
    algo_state->number_of_cities = atoi(strtok(buffer, " "));
    cities = calloc(algo_state->number_of_cities, sizeof(struct city));
    algo_state->number_of_roads = atoi(strtok(NULL, " "));

    all_cities_visited_mask = 0;
    for (int i = 0; i < algo_state->number_of_cities; ++i)
    {
        binary_masks[i] = power2(2, i);
        all_cities_visited_mask += binary_masks[i];
    }

    algo_state->all_cities_visited_mask = all_cities_visited_mask;

    while (fgets(buffer, 1024, cities_fp) != NULL)
    {
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

void dealloc_data()
{
    free(cities);
}

char compare_paths(struct Tour *a, struct Tour *b)
{
    // resolves which path is better (and therefore will have the higher priority)
    // Paths with smaller costs are better
    if (a->cost > b->cost)
        return '\0';
    return '\1';
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
    // init with max value
    steps = malloc(INT_MAX);
    free_steps = malloc(INT_MAX);

    printf("Step pointer %p\n", steps);
    printf("Free steps pointer %p\n", free_steps);
    fflush(stdout);
    // inicialize the step ids
    for (int i = 0; i <= INT_MAX / sizeof(union step); i++)
    {
        steps[i].Tour.stepID = i;
        free_union_step(i);
    }

    algo_state.queue = queue_create((char (*)(void *, void *))compare_paths);
    tscp(&algo_state);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(&algo_state); // to the stdout!

    queue_delete(algo_state.queue);
    dealloc_data();
}
