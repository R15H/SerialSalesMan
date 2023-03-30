#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "tscp.h"
#include <stdlib.h>

#define ERROR(message) fprintf(stderr, #message);

#ifdef DEBUG
#define MESSAGE(message, ...) printf(stderr, message, __VA_ARGS__)
#else

#define MESSAGE(message, ...)
#endif

#define DOUBLE_MAX 1.7976931348623155e+308

#pragma runtime_checks("", off)



#define REALLOC_SIZE 1024
#define SWAP(x, y) void* tmp = x; x = y; y = tmp;

double deviation = 0;
// Return the index of the parent node
static size_t parent_of(size_t i)
{
    return (i - 1) / 2;
}



inline int compare_paths(struct Tour * path1, struct Tour* path2){

    return path1->cost < path2->cost;
    //return (path1)->cost + (-path1->nr_visited + path2->nr_visited)* deviation < ((struct Tour*) path2)->cost;
}


// Bubble-down the element to the correct position
// (i.e., compare it to its child and then swap them if necessary).
// Assume that all the elements in the subtree is already sorted.
void bubble_down(priority_queue_t *queue, size_t node)
{
    size_t left_child = 2 * node + 1;
    size_t right_child = 2 * node + 2;
    size_t i = node;

    // Compare with the left node
    if (left_child < queue->size && compare_paths(queue->buffer[node], queue->buffer[left_child])
        //compare_paths(queue->buffer[node], queue->buffer[left_child])
            )
    {
        i = left_child;
    }

    // Compare with the right node
    if (right_child < queue->size && compare_paths(queue->buffer[i], queue->buffer[right_child])
       // (((struct Tour*) queue->buffer[i])->cost < ((struct Tour*)(queue->buffer[right_child]))->cost)

        //compare_paths(queue->buffer[i], queue->buffer[right_child])
            )
    {
        i = right_child;
    }

    // If node is not in the correct position, swap and then sort the subtree
    if (i != node)
    {
        SWAP(queue->buffer[i], queue->buffer[node])
        bubble_down(queue, i);
    }
}

// Create a new priority queue
priority_queue_t *queue_create(char (*cmp)(void *, void *))
{
    priority_queue_t *queue;

    queue = malloc(sizeof(priority_queue_t));

    queue->buffer = malloc(REALLOC_SIZE * sizeof(void*));
    queue->max_size = REALLOC_SIZE;
    queue->size = 0;
    queue->cmpfn = cmp;

    return queue;
}

// Delete the priority queue
void queue_delete(priority_queue_t *queue)
{
    queue->size = -1;
    queue->max_size = -1;
    free(queue->buffer);
}

// Insert a new element in the queue and then sort its contents.
void queue_push(priority_queue_t *queue, void* new_element)
{
    // Reallocate buffer if necessary
    if (queue->size + 1 > queue->max_size)
    {
        queue->max_size += REALLOC_SIZE;
        queue->buffer = realloc(queue->buffer, queue->max_size * sizeof(void*));
    }

    // Insert the new_element at the end of the buffer
    size_t node = queue->size;
    queue->buffer[queue->size++] = new_element;

    // Bubble-up the new element to the correct position
    // (i.e., compare it to the parent and then swap them if necessary)
    while (node > 0 &&
    compare_paths(queue->buffer[node], queue->buffer[parent_of(node)])
            )

    {
        size_t parent = parent_of(node);
        SWAP(queue->buffer[node], queue->buffer[parent])
        node = parent;
    }
}

// Return the element with the lowest value in the queue, after removing it.
void* queue_pop(priority_queue_t *queue)
{
    if(queue->size == 0)
        return NULL;

    // Stores the lowest element in a temporary
    void* top_val = queue->buffer[0];

    // Put the last element in the queue in the front.
    queue->buffer[0] = queue->buffer[queue->size - 1];

    // Remove the duplicated element in the back.
    --queue->size;

    // Sort the queue based on the value of the nodes.
    bubble_down(queue, 0);

    return top_val;
}

// Duplicate queue
priority_queue_t *queue_duplicate(priority_queue_t* queue)
{
    priority_queue_t *other;
    other = malloc(sizeof(priority_queue_t));
    other->max_size = queue->max_size;
    other->size = queue->size;
    other->cmpfn = queue->cmpfn;
    other->buffer = malloc(queue->max_size * sizeof(void*));
    memcpy(other->buffer, queue->buffer, queue->max_size * sizeof(void*));

    return other;
}

// Print the contents of the priority queue
void queue_print(priority_queue_t* queue, FILE *fp,
                 void (*print_node)(FILE *, void*))
{
    priority_queue_t *queue_copy = queue_duplicate(queue);

    while (queue_copy->size > 0)
    {
        void* node = queue_pop(queue_copy);
        print_node(fp, node);
    }

    queue_delete(queue_copy);
    free(queue_copy);
}
void* remove_element(priority_queue_t *queue, size_t node)
{
    if (node >= queue->size) {
        return NULL;
    }

    // Swap the node we want to remove with the root node
    SWAP(queue->buffer[0], queue->buffer[node]);

    // Remove the root node (which is now the node we want to remove)
    void* removed_val = queue_pop(queue);

    return removed_val;
}


void queue_trim(priority_queue_t *queue, double maxCost){
    if(queue->size < 2000) return;
    printf("Solution found... %f\n", maxCost);
    for(int j = 0; j < queue->size-1; j++){
        if(((struct Tour*)queue->buffer[j])->cost >= maxCost){
            void *toFree = queue->buffer[j];
            remove_element(queue, j);
            free_tour(toFree);
        }
    }
    /*
    for(int j = 0; j< queue->size-1; j += 2){
        if(((struct Tour*)queue->buffer[j+1])->cost >= maxCost){
            free_tour(queue->buffer[j+1]);
            remove_element(queue, j+1);
        }
    }
     */
}
































inline int get_visited_all_cities(struct Tour *tour, struct AlgorithmState *algo_state) {
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

    int cities_visited = algo_state->number_of_cities + 1;
    int values[50]; // max 50 cities

    printf("%.1f\n", algo_state->solution->cost/2);

    for(int i = 0; i < cities_visited; i++)
        printf("%d ", algo_state->solution->visited_list[i]);
}

inline void free_tour(struct Tour *tour) {
    free(tour);
   // printf("Tour freed...");
}

double get_global_lower_bound(int number_of_cities, struct city *cities) {
    double lower_bound = 0;
    for (int i = 0; i < number_of_cities; i++) lower_bound = cities[i].min_cost + cities[i].min_cost2;
    return lower_bound ;
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
    return lower_bound + jump_cost - (ct + cf) ;
}

inline struct Tour *go_to_city(struct Tour *tour, short city_id, struct AlgorithmState *algo_state, double cost) {
    struct Tour *new_tour = malloc(sizeof(struct Tour) + sizeof(short) * (algo_state->number_of_cities+1));
    memcpy(new_tour, tour, sizeof(struct Tour) + sizeof(short) * (tour->nr_visited+1));
    //new_tour->nr_visited = tour->nr_visited + 1; // this plus can be done only once

    new_tour->visited_list[++new_tour->nr_visited] = city_id;
    new_tour->cities_visited = tour->cities_visited | binary_masks[city_id];
    new_tour->cost = cost;
    return new_tour;
}

void reused_go_to_city(struct Tour *tour, short city_id, struct AlgorithmState *algo_state, double cost){
    tour->visited_list[++tour->nr_visited] = city_id;
    tour->cost = cost;
    tour->cities_visited = tour->cities_visited | binary_masks[city_id];
}

/*
int shouldnt_create_tour(double cost, int current_city, struct AlgorithmState *algo_state) {
    int step_exceeds_max_cost = cost > algo_state->max_lower_bound;
    int step_worst_than_found_solution = cost > algo_state->solution->cost;
    return step_exceeds_max_cost ||
           step_worst_than_found_solution;
}
 */

void visit_city(struct Tour *tour,int destination, struct AlgorithmState *algo_state, int *tours_created){
    int i = destination;
    if (!get_was_visited(tour, i)) {
        double new_cost = compute_updated_lower_bound(tour->cost, tour->visited_list[tour->nr_visited], i);
        int discard_tour = new_cost > algo_state->solution->cost;

        if (!discard_tour) {
            struct Tour *new_tour = tour;
            if(tours_created == 0) reused_go_to_city(tour, i, algo_state, new_cost);
            else new_tour = go_to_city(tour, i, algo_state, new_cost);
            (*tours_created)++;
            int finished = get_visited_all_cities(new_tour, algo_state);
            if (!finished) {
                queue_push(algo_state->queue, new_tour);
            } else {
                double final_cost = compute_updated_lower_bound(new_tour->cost, i, 0);
                int current_tour_is_better = final_cost < algo_state->solution->cost;
                if (current_tour_is_better) {
                    free_tour(algo_state->solution);
                    reused_go_to_city(new_tour, 0, algo_state, final_cost);
                    algo_state->solution = new_tour;//go_to_city(new_tour, 0, algo_state, final_cost);
                    queue_trim(algo_state->queue, final_cost);
                } else {
                    free(new_tour); // not free_tour because we only want to delete this piece, and do not wnat to look to prev step
                }
            }
        }
    }
}

void reused_visit_city(struct Tour *tour,int destination, struct AlgorithmState *algo_state, int *tours_created){
    int i = destination;
    if (!get_was_visited(tour, i)) {
        double new_cost = compute_updated_lower_bound(tour->cost, tour->visited_list[tour->nr_visited], i);
        int discard_tour = new_cost > algo_state->solution->cost;

        if (!discard_tour) {
            reused_go_to_city(tour, i, algo_state, new_cost);
            struct Tour *new_tour = tour;
            int finished = get_visited_all_cities(new_tour, algo_state);
            if (!finished) {
                (*tours_created)++;
                queue_push(algo_state->queue, new_tour);
            } else {
                double final_cost = compute_updated_lower_bound(new_tour->cost, i, 0);
                int current_tour_is_better = final_cost < algo_state->solution->cost;
                if (current_tour_is_better) {
                    (*tours_created)++;
                    free_tour(algo_state->solution);
                    reused_go_to_city(new_tour, 0, algo_state, final_cost);
                    algo_state->solution = new_tour;
                    queue_trim(algo_state->queue, final_cost);
                } else {
                    free(new_tour); // not free_tour because we only want to delete this piece, and do not wnat to look to prev step
                }
            }
        }
    }
}
inline int  analyseTour(struct Tour *tour, struct AlgorithmState *algo_state) {
    int tours_created = 0;
    int loops = algo_state->number_of_cities - 1;
    struct Tour this_tour = *tour;
    int i; // ignore the city 0, which we will only return in the end
    for(i =1; i < algo_state->number_of_cities; i++) {
        visit_city(tour, i, algo_state, &tours_created);
        //visit_city(tour, i+1, algo_state, &tours_created);
        if(tours_created) break;
    }
    for(; i < algo_state->number_of_cities; i++) {
        visit_city(&this_tour, i, algo_state, &tours_created);
    }

    queue_trim(algo_state->queue, algo_state->solution->cost);


    //if(loops % 2 == 0) loops--;
    //if (algo_state->number_of_cities % 2 == 0) {
        //reused_visit_city(tour,algo_state->number_of_cities - 1, algo_state, &tours_created); // loops is ajusted so this always runs
    //}

    return tours_created;
}


void tscp(struct AlgorithmState *algo_state) {
    algo_state->solution = malloc(sizeof (struct Tour)
            + sizeof(short) * (algo_state->number_of_cities+1)
            );
    algo_state->solution->cost = algo_state->max_lower_bound;
    algo_state->solution->cities_visited = algo_state->all_cities_visited_mask;
    algo_state->solution->nr_visited = 63;

    struct Tour *first_step = malloc(sizeof (struct Tour)
                                     + sizeof(short) * (algo_state->number_of_cities+1)
            );
    first_step->visited_list[0] = 0;
    first_step->cities_visited = 1;
    first_step->nr_visited = 0;
    first_step->cost = get_global_lower_bound(algo_state->number_of_cities, cities);

    queue_push(algo_state->queue, first_step);

    struct Tour *current_tour;
    while ((current_tour = queue_pop(algo_state->queue))) {
        int newToursCreated = analyseTour(current_tour, algo_state);
        if (newToursCreated == 0) {
            free_tour(current_tour);
        }
    }
}

double average_cost = 0;

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
    (current_city->cost)[city_destination] = cost*2;
}

// Initializes the AlgorithmState variable
void parse_inputs(int argc, char **argv, struct AlgorithmState *algo_state) {
    printf("Parsing inputs...\n");
    if (argc < 3) {
        ERROR("Not enough arguments passed. What are you doing?\n")
        exit(-1);
    };
    char *cities_file = argv[1];
    algo_state->max_lower_bound = atof(argv[2])*2;

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

    double *road_cost = malloc(sizeof(double)* algo_state->number_of_roads); //algo_state->number_of_roads];
    int number_of_roads_read = 0;
    while (fgets(buffer, 1024, cities_fp) != NULL) {
        int city_number = atoi(strtok(buffer, " "));
        int city_destination = atoi(strtok(NULL, " "));
        double city_cost = strtod(strtok(NULL, " "), NULL);
        place_cost_in_city(city_number, city_destination, city_cost, algo_state->number_of_cities);
        place_cost_in_city(city_destination, city_number, city_cost, algo_state->number_of_cities);
        average_cost += city_cost;
        road_cost[number_of_roads_read++] = city_cost;
    }
    // compute the standard variance of the costs
    average_cost = average_cost/number_of_roads_read;
    for (int i = 0; i < number_of_roads_read; ++i) {
        deviation += (road_cost[i] - average_cost) * (road_cost[i] - average_cost);
    }
    deviation = sqrt(deviation/number_of_roads_read);

    fclose(cities_fp);
    free(road_cost);
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
