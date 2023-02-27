#ifndef UNTITLED_TSCP_H
#define UNTITLED_TSCP_H

#include "queue.h"

struct city {
    int id;

    // n-1 connections --> all nodes are connected to each other
    short **cities;
    double **cost;
};

#define MAX_CITIES 32
struct city CitiesTable[MAX_CITIES];
unsigned int binary_mask[MAX_CITIES];


struct Tour {
    short current_city;
    struct step_middle *previous_step;
    double cost;
    unsigned int cities_visited; // bit map of the cities visited
};

struct step_middle {
    short current_city;
    struct step_middle *previous_step;
    omp_lock_t decrease_counter_lock; // 8 bytes, replaces cost
    unsigned int ref_counter; // nr of paths it belongs to
};

union step {
    struct Tour;        // we know a step head is head because it is retrieved from the queue
    struct step_middle; // we know a step middle is middle because it is accessed through a pointer to a step middle
};


struct AlgorithmState {
    int max_lower_bound;
    int number_of_cities;
    int number_of_roads;
    struct city *cities;
    priority_queue_t *queue;
    struct Tour *solution;
};

#endif UNTITLED_TSCP_H
