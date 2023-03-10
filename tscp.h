#ifndef UNTITLED_TSCP_H
#define UNTITLED_TSCP_H

#include "queue.h"


struct city {
    short id;

    // n-1 connections --> all nodes are connected to each other
//    struct road **roads;

    //int **cities; // se tiverem organizadas por cost ent precisamos desta linha caso contrario podemos ter simplesment os costs
    double min_cost;
    double min_cost2;
    double *cost; // --> test if this is better -> entry with index "id" contains the second lowest edge cost
};

struct city *cities;

#define MAX_CITIES 32
struct city CitiesTable[MAX_CITIES];
unsigned int binary_masks[MAX_CITIES];
unsigned int all_cities_visited_mask;

struct step_middle {
    short current_city;
    int stepID;
    struct step_middle *previous_step;
    omp_lock_t decrease_counter_lock; // 8 bytes, replaces cost
    unsigned int ref_counter; // nr of paths it belongs to
};

struct Tour {
    short current_city;
    int stepID;
    struct step_middle *previous_step;
    double cost;
    unsigned int cities_visited; // bit map of the cities visited
};



union step {
    struct step_middle step_middle; // we know a step middle is middle because it is accessed through a pointer to a step middle
    struct Tour Tour;        // we know a step head is head because it is retrieved from the queue
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

#endif // UNTITLED_TSCP_H
