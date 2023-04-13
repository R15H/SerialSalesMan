#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "tscp.h"
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>

#define ERROR(message) fprintf(stderr, #message);

#ifdef DEBUG
#define MESSAGE(message, ...) printf(stderr, message, __VA_ARGS__)
#else

#define MESSAGE(message, ...)
#endif

#define DOUBLE_MAX 1.7976931348623155e+308

#pragma runtime_checks("", off)


struct SerialTour {
    unsigned short cities_visited;
    int nr_visited;
    int current_city;
    double cost;
    double lb;
    short cities[1];
};
int nr_of_cities;

MPI_Datatype MPI_SerialTour;
MPI_Datatype MPI_Order;

int nr_processes = 0;
// requests for sending the solution
MPI_Request sol_requests[64] = {0};
double sol_values[64] = {999999999999999};

// variables for receiving soluitons
double solution_cost;
MPI_Request solution_cost_request_send;
int p,id;

#define REALLOC_SIZE 1024
#define SWAP(x, y) void* tmp = x; x = y; y = tmp;

double deviation = 0;
// Return the index of the parent node
static size_t parent_of(size_t i)
{
    return (i - 1) / 2;
}

void print_tour(struct Tour* tour){
    printf("Cost %f\nCities visited: %d\nPrev step %p", tour->cost, tour->cities_visited, tour->previous_step);
}


int compare_paths(struct Tour * path1, struct Tour* path2){

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

void pq(FILE *file, struct Tour * node){
    print_tour(node);
}

void queue_trim(priority_queue_t *queue, double maxCost){
    //printf("Solution found... %d", rand());
    if(queue->size < 1000){
        return;
    }
    for(int j = queue->size/2; j< queue->size-1; j += 2){
        if(((struct Tour*)queue->buffer[j])->cost >= maxCost){
            free_tour(queue->buffer[j]);
            remove_element(queue, j);
        }
        if(((struct Tour*)queue->buffer[j+1])->cost >= maxCost){
            free_tour(queue->buffer[j+1]);
            remove_element(queue, j+1);
        }
    }
}






























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

void serializeTour(struct Tour *tour, struct SerialTour *serial_tour){
    serial_tour->cities_visited =  tour->cities_visited;
    serial_tour->nr_visited = tour->nr_visited;
    serial_tour->cities_visited = tour->cities_visited;
    serial_tour->current_city = tour->current_city;
    serial_tour->cost = tour->cost;
    serial_tour->lb = tour->lb;
    struct step_middle *step = tour->previous_step;
    int i = tour->nr_visited;
    while(step != NULL){
        assert(i >=0 );
        serial_tour->cities[i] = step->current_city;
        i--;
        step = step->previous_step;
    }
}


void print_result(struct AlgorithmState *algo_state, struct SerialTour *tour) {
    if (algo_state->solution->cost == algo_state->max_lower_bound) {
        printf("NO SOLUTION");
        return;
    }

    struct step_middle *step = (struct step_middle *) algo_state->solution; // this is actually a Tour but its okay
    int cities_visited = algo_state->number_of_cities + 1;
    int values[50]; // max 50 cities

    printf("%.1f\n", tour->cost/2);

    //print the path
    for(int i = 0; i < cities_visited; i++){
        printf("%d ", tour->cities[i]);
    }

    /*for (int i = cities_visited - 1; i >= 0; i--) {
        values[i] = step->current_city;
        step = step->previous_step;
        if (step == NULL) break;
    } */
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
    return lower_bound ;
}

double compute_updated_lower_bound(double lower_bound, double jump_cost,  unsigned int source_city, unsigned int destination_city) {
    int comp_1 = jump_cost >= cities[source_city].min_cost2;
    double ct = (comp_1) * cities[source_city].min_cost2
                + (comp_1 == 0) * cities[source_city].min_cost;
    int comp_2 = jump_cost >= cities[source_city].min_cost2;
    double cf =
            (comp_2) * cities[destination_city].min_cost2 +
            (comp_2 == 0) * cities[destination_city].min_cost;
    return lower_bound + jump_cost - (ct + cf) ;
}

struct Tour *go_to_city(struct Tour *tour, short city_id, struct AlgorithmState *algo_state, double lb, double cost) {
    // create a new tour and convert the current tour to a step middle
    struct Tour *new_tour = (struct Tour *) get_clean_step();
    new_tour->current_city = city_id;
    new_tour->nr_visited = tour->nr_visited + 1;
    new_tour->cities_visited = tour->cities_visited | binary_masks[city_id];
    new_tour->lb = lb;
    new_tour->cost = cost;
    new_tour->previous_step = (struct step_middle *) tour;
    return new_tour;
}

//int shouldnt_create_tour(double cost, int current_city, struct AlgorithmState *algo_state) {
    //int step_exceeds_max_cost = cost > algo_state->max_lower_bound;
    //int step_worst_than_found_solution = cost > algo_state->solution->cost;
    //return step_exceeds_max_cost ||
           //step_worst_than_found_solution;
//}

#define SOLUTION_FOUND_TO_MASTER 13456

void scatter_solution(double solution_cost, double sol_lb){
    static MPI_Request request = NULL;
    struct solution_data current_solution;
    MPI_Status status;
    int flag;
    sol_values[id] = solution_cost;
    printf("Solution found by %d with cost: %f\n", id, solution_cost);
    fflush(stdout);
    current_solution.sol_cost = solution_cost;
    current_solution.sol_lb = sol_lb;
    if(request != NULL){
        MPI_Test(&request, &flag, &status);
        // can only broadcast if previous broadcast finished
        if(flag){
            MPI_Isend(&current_solution , 2, MPI_DOUBLE, 0, SOLUTION_FOUND_TO_MASTER ,MPI_COMM_WORLD, &request);//, &sol_requests[id]); // the variable being broadcasted is a pointer, as so it cannot be in the stack!
            printf("id: %d scattering a solution with cost: %f\n", id, solution_cost);
        }
        // TODO either wait for the previous broadcast to finish, or sequechule send... (OMP task?)
    } else MPI_Isend(&current_solution , 2, MPI_DOUBLE, 0, SOLUTION_FOUND_TO_MASTER ,MPI_COMM_WORLD, &request);//, &sol_requests[id]); // the variable being broadcasted is a pointer, as so it cannot be in the stack!
}


double queue_health[64];
double queue_health_requests[64];
double health_deviation[64];


#define QUEUE_HEALTH_STATUS 10

#define MAX_QUEUE_PACKETS 10
#define QUEUE_PACKETS_SIZE 1
struct remote_proc {
    double q_health[MAX_QUEUE_PACKETS][QUEUE_PACKETS_SIZE];
    double q_health_deviation[QUEUE_PACKETS_SIZE];
};

#define MAX_PROCESSES 64

int orders[MAX_PROCESSES];

#define QUEUE_DISTRIBUTE_GIVE 123048921
#define QUEUE_DISTRIBUTE_TAKE 123048922
#define QUEUE_DISTRIBUTION_ORDERS_SEND 21938824
#define QUEUE_DISTRIBUTION_ORDERS_RECEIVE 21933884
struct remote_proc remotes[64];


void distribute_load(){
    if(id == 0){

    }
}



int compare(const double *a, const double *b){
    return *a < *b;
}


struct Tour* deserializeTour(struct SerialTour *serialTour){
    // Another option would be to have the tour store a pointer to an array, instead of building the linked list, since the cities visited are only usefull when presenting the solution
    // Drawbacks:
    // more complex free (one more if), (but deserialization is would be more imidiate)
    struct Tour *tour = malloc(sizeof(struct Tour));
    tour->cities_visited = serialTour->cities_visited;
    tour->nr_visited = serialTour->nr_visited;
    tour->current_city = serialTour->current_city;
    tour->cost = serialTour->cost;
    tour->lb = serialTour->lb;
    struct step_middle *this_step = (struct step_middle *)  tour;
    int i = serialTour->nr_visited;
    for(; i < tour->nr_visited; i--){
        struct step_middle *step = malloc(sizeof(struct step_middle));
        step->current_city = serialTour->cities[i];
        this_step->previous_step = step;
        this_step = step;
    }
    //free(serialTour->cities);
    return tour;
}



struct order {
    int from;
    int to;
};



MPI_Datatype MPI_Health_Packet;
struct queue_health_packet {
   int from;
   double sum;
};

struct SerialTour * initialize_serial_tours(struct SerialTour** tours, int nr_of_tours){
    printf("Initializing serial tours vector\n");
    for(int t=0; t < nr_of_tours;t++ ){
        tours[t] = malloc((sizeof(struct SerialTour) +
                ((nr_of_cities-1) * sizeof(unsigned short) )
               )
                       );
        printf("%p\n", tours[t]);
    }

}

struct queue_health_packet  *packets_sums = NULL;
#define MAX_ORDERS_PER_PROCESS 10


int order_recv[MAX_PROCESSES][MAX_ORDERS_PER_PROCESS] =    {0,0};
int order_recv_num[MAX_PROCESSES]  = {0,0};
int order_send[MAX_PROCESSES][MAX_ORDERS_PER_PROCESS] = {0,0};
int order_send_num[MAX_PROCESSES] = {0,0};
void load_balance(priority_queue_t *queue){ // reports the queue healthy and executes orders from the master if there are any
    // iterate through the first 100 items of the queue and average their LB


    struct queue_health_packet sum[MAX_QUEUE_PACKETS];
    // Access the health of each queue packet
    int j = 0;
    int i = 0;
    printf("Here %d %d %ld %d\n", MAX_QUEUE_PACKETS,j, queue->size, id);
    while(j < MAX_QUEUE_PACKETS) {
        sum[j].from = id;
        sum[j].sum = 0;
        int this_it =0;
        for(; this_it < QUEUE_PACKETS_SIZE && queue->size > i ; i++ ){
            printf("LB %f  i: %d j: %d\n",(queue->buffer[i])->lb,i,j);
            sum[j].sum += (queue->buffer[i])->lb;
            this_it++;
        }
        if(this_it < QUEUE_PACKETS_SIZE) // did not finish
            sum[j].sum = 0;
        j++;
        printf("[%d] %d Sum: %f\n", id,j, sum[j].sum);
    }


    // Get all queue packets from all processes
    //MPI_Health_Packet
    MPI_Gather( &sum, (MAX_QUEUE_PACKETS)*sizeof(struct queue_health_packet), MPI_BYTE, packets_sums,  (MAX_QUEUE_PACKETS) *(sizeof(struct queue_health_packet)), MPI_BYTE, 0, MPI_COMM_WORLD); // TODO check optimization diferent values for counts
    // print the contents of packets_sums

    if(id == 0) {
        for(int k = 0; k < nr_processes * MAX_QUEUE_PACKETS; k++){
            printf("%d from: %d sum: %f\n",k, packets_sums[k].from, packets_sums[k].sum);
        }

        // Match queue packets with each other
        qsort(packets_sums, MAX_QUEUE_PACKETS * nr_processes, sizeof(struct queue_health_packet),
              (int (*)(const void *, const void *)) compare);
        int start = 0;
        int end = MAX_QUEUE_PACKETS*nr_processes-1;
        for (; start >= end; start++, end--) {
            int from = packets_sums[start].from;
            int receiver = packets_sums[end].from;
            if(packets_sums[start].sum == 0 ){
                order_send[from][order_send_num[from]++] = -1;
                order_recv[receiver][order_recv_num[receiver]++] = -1;
                break; // do not send empty packets! or half empty (this way is simpler, no need to do the average 4 example)
            }

            order_send[from][order_send_num[from]++] = packets_sums[end].from;
            order_recv[receiver][order_recv_num[receiver]++] = packets_sums[start].from;
        }
        // This is a scatter but since order is not guaranteed we need to send the data individually to each process
        static MPI_Request orders_req_send[MAX_PROCESSES];
        static MPI_Request orders_req_recv[MAX_PROCESSES];

        for(int i = 0; i < nr_processes; i++){
            MPI_Isend(&order_send[i], MAX_ORDERS_PER_PROCESS, MPI_INT, i, QUEUE_DISTRIBUTE_GIVE, MPI_COMM_WORLD, &orders_req_send[i]); //
            MPI_Isend(&order_recv[i], MAX_ORDERS_PER_PROCESS, MPI_INT, i, QUEUE_DISTRIBUTE_TAKE, MPI_COMM_WORLD, &orders_req_recv[i]);
        }
    }

    MPI_Recv(&order_send[id], MAX_ORDERS_PER_PROCESS, MPI_INT, 0, QUEUE_DISTRIBUTE_GIVE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&order_recv[id], MAX_ORDERS_PER_PROCESS, MPI_INT, 0, QUEUE_DISTRIBUTE_TAKE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Iterate through all send orders and asynchronously send the queue packets to the other process
    for(int i = 0; i < MAX_ORDERS_PER_PROCESS; i++){
        int to = order_send[id][i];
        if(to == -1) break;
        // Serialize the queue packets
        static struct SerialTour *tours_to_send[MAX_QUEUE_PACKETS * QUEUE_PACKETS_SIZE] = {NULL}; // we guarantee this is not being reused since this loadbalancing is called very sporadically
        if(tours_to_send[0] == NULL) initialize_serial_tours(tours_to_send, MAX_QUEUE_PACKETS*QUEUE_PACKETS_SIZE);

        for(int j = 0; j < MAX_QUEUE_PACKETS * QUEUE_PACKETS_SIZE; j++){
            printf("j %p\n", tours_to_send[j]);
            struct Tour *tour = queue_pop(queue);
            serializeTour(tour, tours_to_send[j]);
            free(tour);
        }
        // Send the queue packets to the other process
        //MPI_Send(&tours_to_send, MAX_QUEUE_PACKETS * QUEUE_PACKETS_SIZE, MPI_Tour, to, QUEUE_DISTRIBUTION_ORDERS_SEND, MPI_COMM_WORLD);
        static MPI_Request send_tours_req[MAX_ORDERS_PER_PROCESS];
        MPI_Isend(&tours_to_send, MAX_QUEUE_PACKETS * QUEUE_PACKETS_SIZE, MPI_SerialTour, to, QUEUE_DISTRIBUTION_ORDERS_SEND, MPI_COMM_WORLD, &send_tours_req[i]);
    }

    // Iterate through each possible receive order
    for(int i = 0; i < MAX_ORDERS_PER_PROCESS; i++){
        int from = order_recv[id][i];
        if(from == -1) break;
        // Receive the queue packets from the other process
        static struct SerialTour *tours_to_receive[MAX_QUEUE_PACKETS*QUEUE_PACKETS_SIZE] = {NULL}; // TODO AJUST FOR VARIABLE SIZE
        if(tours_to_receive[0] == NULL) initialize_serial_tours(&tours_to_receive[0], MAX_QUEUE_PACKETS);
        MPI_Recv(&tours_to_receive, MAX_QUEUE_PACKETS * QUEUE_PACKETS_SIZE, MPI_SerialTour, from, QUEUE_DISTRIBUTION_ORDERS_RECEIVE, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // TODO make this async?
        // Deserialize the queue packets and push them to the queue
        for(int j = 0; j < MAX_QUEUE_PACKETS * QUEUE_PACKETS_SIZE; j++){
            struct Tour *tour = deserializeTour(tours_to_receive[j]);
            queue_push(queue, tour);
        }
    }

}


void visit_city(struct Tour *tour,int destination, struct AlgorithmState *algo_state, int *tours_created){
    int i = destination;
    if (!get_was_visited(tour, i)) {
        double jump_cost = get_cost_from_city_to_city(tour->current_city, i);
        double new_lb = compute_updated_lower_bound(tour->lb, jump_cost,tour->current_city, i);

        if (new_lb <= algo_state->sol_lb) {
            struct Tour *new_tour = go_to_city(tour, i, algo_state, new_lb, tour->cost + jump_cost);
            int finished = get_visited_all_cities(new_tour, algo_state);
            if (!finished) {
                (*tours_created)++;
                queue_push(algo_state->queue, new_tour);
            } else {
                jump_cost = get_cost_from_city_to_city(i, 0);
                double final_lb = compute_updated_lower_bound(new_tour->lb,jump_cost ,i, 0);
                double final_cost = new_tour->cost + jump_cost;
                int current_tour_is_better = final_cost < algo_state->sol_cost;
                if (current_tour_is_better) {
                    (*tours_created)++;
                    free_tour(algo_state->solution);
                    algo_state->solution = go_to_city(new_tour, 0, algo_state, final_lb,final_cost);
                    algo_state->sol_cost = final_cost;
                    algo_state->sol_lb = final_lb;
                    scatter_solution(final_cost, final_lb);
                    queue_trim(algo_state->queue, final_cost);
                } else {
                    free(new_tour); // not free_tour because we only want to delete this piece, and do not wnat to look to prev step
                }
            }
        }
    }
}




int  analyseTour(struct Tour *tour, struct AlgorithmState *algo_state) {
    int tours_created = 0;
    int loops = algo_state->number_of_cities - 1;
    int i = 1;
    for (; i < loops; i += 2) {
        visit_city(tour, i, algo_state, &tours_created);
        visit_city(tour, i+1, algo_state, &tours_created);
    }
    if (algo_state->number_of_cities % 2 == 0) {
        visit_city(tour,algo_state->number_of_cities - 1, algo_state, &tours_created);
    }

    return tours_created;
}

int process_with_solution;

void gather_best_solution(struct AlgorithmState * algo_state){
    static double vol_cost[128];
    static struct solution_data solutionData = {999999999999999,9999999999999};
    static MPI_Request gather_solution[64] = {NULL};
    // check if any of the processes has a better solution, iterate over all processes
    // we use this if to periodically check, a new if would consume more resources...
    static MPI_Request solution_broadcast = NULL;
    static int first = 0;
    if(id == 0){

        // receive solution from all processes using MPI_Irecv
        // if a process has a better solution, update the solution
        for (int i = 0; i < nr_processes; ++i) {
            //if(i == id) continue; // do not send to self (we already have the value)
            MPI_Status status;
            int flag;
            int res = 0;
            if(gather_solution[i] == NULL) MPI_Irecv(&vol_cost[i*2], 2, MPI_DOUBLE, i,  SOLUTION_FOUND_TO_MASTER,MPI_COMM_WORLD,  &gather_solution[i]); // get next receive going
            else res = MPI_Test(&gather_solution[i], &flag, &status);

            if (flag) {
                printf("---Received solution from slave %d with cost %.2f--\n", i, vol_cost[i*2]);

                if(algo_state->sol_cost > vol_cost[i*2] && vol_cost[i*2] != 0){
                    solutionData.sol_cost = vol_cost[i*2];
                    solutionData.sol_lb = vol_cost[i*2+1];
                    process_with_solution = i;
                    printf("UPDATED best solution\n");
                }
                MPI_Irecv(&vol_cost[i*2], 2, MPI_DOUBLE, i,  SOLUTION_FOUND_TO_MASTER,MPI_COMM_WORLD,  &gather_solution[i]); // get next receive going
                //MPI_Irecv(&vol_cost[i*2], 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,  &gather_solution[i]); // get next receive going
            }
        }

        if(solutionData.sol_cost < algo_state->sol_cost && solutionData.sol_cost != 0) {
            algo_state->sol_cost = solutionData.sol_cost;
            algo_state->sol_lb = solutionData.sol_lb;
            printf("MASTER new global solution: %f", algo_state->sol_cost);
            MPI_Status status;
            if(solution_broadcast) MPI_Wait(&solution_broadcast, &status);
            MPI_Ibcast(&solutionData , 2,MPI_DOUBLE, 0, MPI_COMM_WORLD, &solution_broadcast);
        } // else do not broadcast, no solution has been found
    }
    else { // a normal process
        MPI_Status status;
        if(solution_broadcast) MPI_Wait(&solution_broadcast, &status); // MPI_Test(&solution_broadcast, &flag, &status);
        MPI_Ibcast(&solutionData , 2,MPI_DOUBLE, 0, MPI_COMM_WORLD, &solution_broadcast);
        if(algo_state->sol_cost > solutionData.sol_cost){
            algo_state->sol_cost = solutionData.sol_cost;
            algo_state->sol_lb = solutionData.sol_lb;
            printf("Process %d received a solution from master: %f\n", id, algo_state->sol_cost);
            queue_trim(algo_state->queue, solutionData.sol_cost);
        }
    }




    /*
    for(int i=0; i < nr_processes; i++){
        MPI_Status status;
        int flag;
        if(i == id) {
            printf("algo_state->solution->cost: %f   sol_values[%d] %f\n", algo_state->solution->cost,i, sol_values[id]);
            if (algo_state->sol_cost < sol_values[id]){
                sol_values[id] = algo_state->sol_cost;
            } else continue;
        }
        //sol_values[id] >= algo_state->solution->cost) continue; // do not send to self (we already have the value)
        if(sol_requests[i] != NULL)
        {
            MPI_Test(&sol_requests[i], &flag, &status);
            if(flag == 1){
                printf("Process %d received solution cost from process %d: %f and has %f local cost\n", id, i, sol_values[i], algo_state->solution->cost);
                if(sol_values[i] < algo_state->solution->cost){
                    algo_state->solution->cost = sol_values[i];
                    process_with_solution = i;
                    //solution_cost = vol_cost[i];
                    printf("Updated solution cost on process %d to %f\n", id, sol_values[i]);
                }
                MPI_Ibcast(&sol_values[i] , 1, MPI_DOUBLE, i, MPI_COMM_WORLD, &sol_requests[i]);
            }
        } else  MPI_Ibcast(&sol_values[i] , 1, MPI_DOUBLE, i, MPI_COMM_WORLD, &sol_requests[i]);
    }
     */
}


void tscp(struct AlgorithmState *algo_state) {
    algo_state->sol_cost = algo_state->max_lower_bound*2;;
    algo_state->sol_lb = algo_state->max_lower_bound;
    algo_state->solution = (struct Tour *) get_clean_step();
    algo_state->solution->cost = algo_state->max_lower_bound*2;
    algo_state->solution->cities_visited = algo_state->all_cities_visited_mask;
    algo_state->solution->current_city = 0;
    algo_state->solution->previous_step = NULL;
    algo_state->solution->lb = algo_state->max_lower_bound;;

    struct Tour *first_step = (struct Tour *) get_clean_step();
    first_step->current_city = 0;
    first_step->nr_visited = 0;
    first_step->cities_visited = 1;
    first_step->lb = get_global_lower_bound(algo_state->number_of_cities, cities);
    first_step->cost = 0;
    first_step->previous_step = NULL;

    queue_push(algo_state->queue, first_step);

    struct Tour *current_tour;
    int i = 0; // nr of tours to analyse before spliting the load
    while ((current_tour = queue_pop(algo_state->queue))) {
        int newToursCreated = analyseTour(current_tour, algo_state);
        i++;
        if (newToursCreated == 0) {
            free_tour(current_tour);
            continue;
        }
        ((struct step_middle *) current_tour)->ref_counter = newToursCreated;
        if(i > 1000) break;
    }


    priority_queue_t *final_queue = queue_create(NULL);
    if(id == 0) printf("There are %ld tours to distribute\n", algo_state->queue->size);
    int current_proc = 0;
    int direction = -1;
    while ((current_tour = queue_pop(algo_state->queue))) {
        //printf("Current proc %d, current_tour %p\n", current_proc, current_tour);
        //fflush(stdout);
        if(current_proc == id) queue_push(final_queue,current_tour);
        else free_tour(current_tour);
        if(current_proc == nr_processes-1 || current_proc == 0) direction = -direction;
        current_proc += direction;
    }
    printf("Process %d has %ld tours to analyse\n", id, final_queue->size);
    queue_delete(algo_state->queue);
    algo_state->queue = final_queue;

    i =0;
    while ((current_tour = queue_pop(algo_state->queue))) {
        int newToursCreated = analyseTour(current_tour, algo_state);
        i++;
        if (newToursCreated == 0) {
            free_tour(current_tour);
            continue;
        }
        ((struct step_middle *) current_tour)->ref_counter = newToursCreated;
        if(i > 5000000) {
            printf("Syncronizing & loadbalancing...\n");
            gather_best_solution(algo_state);
            load_balance(algo_state->queue);
            i = 0;
        }
    }
    // wrap up
    printf("Syncronizing & loadbalancing...\n");
    gather_best_solution(algo_state);
    load_balance(algo_state->queue);

    // TODO check if more work needs to be done before terminating

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
        current_city->min_cost = DOUBLE_MAX;
        current_city->min_cost2 = DOUBLE_MAX;
    }
    if (current_city->min_cost > cost)
        current_city->min_cost = cost;
    else if (current_city->min_cost2 > cost)
        current_city->min_cost2 = cost;
    (current_city->cost)[city_destination] = cost*2;

    current_city = &cities[city_destination];
    if (current_city->min_cost > cost)
        current_city->min_cost = cost;
    else if (current_city->min_cost2 > cost)
        current_city->min_cost2 = cost;
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

    double road_cost[algo_state->number_of_roads];
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
}

void dealloc_data() {
    free(cities);
}

struct SerialTour *send_solutions;
struct SerialTour *receive_solutions;

void create_mpi_struct_tour(int nr_of_cities){
    int block_lengths[6] = {1, 1, 1, 1, 1, nr_of_cities};
    MPI_Aint displacements[6];
    displacements[0] = 0;
    displacements[1] = sizeof(short);
    displacements[2] = displacements[1] + sizeof(int);
    displacements[3] = displacements[2] + sizeof(int);
    displacements[4] = displacements[3] + sizeof(double);
    displacements[5] = displacements[4] + sizeof(double);

    //= {0, sizeof(unsigned short), sizeofsizeof(int), sizeof(int) + sizeof(unsigned int), sizeof(double) + sizeof(int) + sizeof(unsigned int), sizeof(double) + sizeof(int) + sizeof(unsigned int) + sizeof(double)};
    MPI_Datatype types[6] = {MPI_UNSIGNED_SHORT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_SHORT};
    MPI_Type_create_struct(6, block_lengths, displacements, types, &MPI_SerialTour);
    MPI_Type_commit(&MPI_SerialTour);
}

void create_mpi_struct_order(){
    // Create a type for the struct order
    MPI_Aint offsets[2];
    offsets[0] = 0;
    offsets[1] = sizeof(int);
    int block_lengths[2] = {1, 1};
    MPI_Datatype types2[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, block_lengths, offsets, types2, &MPI_Order);
    MPI_Type_commit(&MPI_Order);
}

void create_mpi_struct_packet(){
    // Create a type for the struct order
    MPI_Aint offsets[2];
    offsets[0] = 0;
    offsets[1] = sizeof(int);
    int block_lengths[2] = {1, 1};
    MPI_Datatype types2[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Type_create_struct(2, block_lengths, offsets, types2, &MPI_Health_Packet);
    MPI_Type_commit(&MPI_Health_Packet);
}

int main(int argc, char *argv[]) {
    double exec_time;
    struct AlgorithmState algo_state;
    parse_inputs(argc, argv, &algo_state);


    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nr_processes);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    if(id == 0){
        packets_sums = malloc(sizeof (struct queue_health_packet)* MAX_QUEUE_PACKETS*nr_processes);
    }


    printf("Nr of proccesses %d\n", nr_processes);

    send_solutions =
            malloc((sizeof(struct SerialTour) +
            ((algo_state.number_of_cities-1) * sizeof(unsigned short) )
            ) * nr_processes);
    receive_solutions = malloc((sizeof(struct SerialTour) +
            ((algo_state.number_of_cities-1) * sizeof(unsigned short))
            ) * nr_processes);



    create_mpi_struct_order();
    create_mpi_struct_tour(algo_state.number_of_cities);
    create_mpi_struct_packet();


    nr_of_cities = algo_state.number_of_cities;





    MPI_Barrier(MPI_COMM_WORLD);

    printf("Process %d ready, %d processes in total\n", id, p);
    exec_time = -MPI_Wtime();

//exec_time = -omp_get_wtime();
    algo_state.queue = queue_create((char (*)(void *, void *)) NULL);
    tscp(&algo_state);
    // gather algo_state solution from all processes
    serializeTour(algo_state.solution, &send_solutions[id]);
    printf("GATHER ALL SOLUTIONS Sending solution from process %d, cost: %.2f, lb: %.2f\n", id, send_solutions[id].cost, send_solutions[id].lb);
    exec_time += MPI_Wtime();
    fprintf(stderr, "%.1fs\n", exec_time);
    fflush(stderr);
    MPI_Gather( &send_solutions[id], 1, MPI_SerialTour, receive_solutions, 1, MPI_SerialTour, 0, MPI_COMM_WORLD);
    int best_solution_id = 0;
    printf("ALL GOOD\n");


    fflush(stdout);
    if(id == 0){
    for(int i =0; i < nr_processes; i++){
        struct SerialTour* serial_tour = &receive_solutions[i];
        printf("GATHER ALL SOLUTIONS Received on root: process: %d, cost: %.2f, lb: %.2f.\n", i, serial_tour->cost, serial_tour->lb);
        serial_tour->cost < receive_solutions[best_solution_id].cost ? best_solution_id = i : 0;
    }
    }

//exec_time += omp_get_wtime();
    //if(id == 0) print_result(&algo_state, &receive_solutions[best_solution_id]);
    queue_delete(algo_state.queue);
    dealloc_data();
    //MPI_Type_free(&MPI_SerialTour);
    MPI_Finalize();
}
