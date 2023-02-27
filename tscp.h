//
// Created by Tester on 2/23/2023.
//

#ifndef UNTITLED_TSCP_H
#define UNTITLED_TSCP_H

#include "queue.h"
struct city {
    int id;

    // stores the indexes of the cities that are connected to this city
    // We use indexes because ints are smaller than pointers
    // They sorted in ascending order, so we can perform binary search
    int **cities;

    // if we just want to access linearly -> stop flag -> we don't need to store the size, nor increase "i"

    // we will iterate through cities, if we find a city with this flag, we will stop
    // the index we calculate from here will be a valid index for the cost array
    int size;
    double *cost;       // index i -> city i
    int nr_cities;      // number of cities
};

//struct transversal {
    //int *path; // array of cities <-- can be as big as 2x the size of the number of cities      <--- we need to check
    //double cost; // total cost until now
//};


struct transversal_step { // union de step_head com step_middle  --> step_head precisa de ter o custo do caminho mas nao precisa de ter o integrates_paths nem o bit lock
    struct city *current_city;
    struct transversal_step *previous_step;
    double cost;
    unsigned int nr_cities_visited; // this may or may not be the length of the path
    unsigned short integrates_paths; // nr of paths that depend (directly) on this step  --> one thread adds all paths --> subtracing by multiple threads
    // if 1 --> we can free this step (last thread does not need to aquire the lock)


    // bit that locks
    omp_lock_t decrease_counter_lock; // 1 if locked, 0 if not locked   <--- 8 bytes
    // MAX nr of paths unsigned short --> 65535
};


// Space used:

// com pooling --> step addr can become relative to the pool --> fewer bytes used


// 8 bytes = pointer
// 2 bytes -> unsigned short poolID
// 2 bytes -> unsigned short stepID  6 bytes total


struct stand_alone_path { // multiple steps in one datastructure
    int **cities; // <--- relative pointers -> ocupy fewer bytes -> with SMID instructions we can perform  operations in 8 addresses at the same time
    double cost;
    short int path_length;
    //unsigned int nr_cities_visited; // this may or may not be the length of the path
};

struct step_head {
    int *current_city;
    struct transversal_step *previous_step;
    double cost;
    unsigned int nr_cities_visited; // this may or may not be the length of the path
};
// cache line -> 64 bytes
struct step_middle {
    struct city *current_city; // set of cities to support SIMD
    struct transversal_step *previous_step;
    omp_lock_t decrease_counter_lock; // 8 bytes, replaces cost
    unsigned int nr_belongs_to; // nr of paths it belongs to  , replaces nr_of_cities_visited
};

union step {
    struct step_head;   // we know a step head is head because it is retreived from the queue
    struct step_middle; // we know a step middle is middle because it is accessed through a pointer to a step middle
};

// each thread will have a pool of spaces for transversal steps
// this way we do not have to call free all the time & memory is not fragmented

struct AlgorithmState {
    int max_lower_bound;
    int number_of_cities;
    int number_of_roads;
    struct city *cities;
    priority_queue_t *queue;
    struct transversal_step *solution;
};

#endif //UNTITLED_TSCP_H
