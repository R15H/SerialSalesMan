#ifndef UNTITLED_TSCP_H
#define UNTITLED_TSCP_H

#include "queue.h"
#include <omp.h>


struct city {
    double min_cost;
    double min_cost2;
    double *cost;
};

struct city *cities;

#define MAX_CITIES 32
struct city CitiesTable[MAX_CITIES];
unsigned int binary_masks[MAX_CITIES];
unsigned int all_cities_visited_mask;


struct Tour {
    short nr_visited;
    double cost;
    unsigned int cities_visited; // bit map of the cities visited --> needs to be 64 bits unsigned long long
    short visited_list[15];
};





struct AlgorithmState {
    double max_lower_bound;
    int number_of_cities;
    int number_of_roads;
    //struct city *cities;
    priority_queue_t *queue;
    struct Tour *solution;
    unsigned int all_cities_visited_mask;
};

void free_tour(struct Tour *tour);
#endif // UNTITLED_TSCP_H
