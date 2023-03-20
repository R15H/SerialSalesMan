#include "queue.h"
#include "tscp.h"

#define REALLOC_SIZE 1024
#define SWAP(x, y) void* tmp = x; x = y; y = tmp;

// Return the index of the parent node
static size_t parent_of(size_t i)
{
	return (i - 1) / 2;
}

void print_tour(struct Tour* tour){
    printf("Cost %f\nCities visited: %d\nPrev step %p", tour->cost, tour->cities_visited, tour->previous_step);
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
	if (left_child < queue->size &&
            (((struct Tour*) queue->buffer[node])->cost < ((struct Tour* ) queue->buffer[left_child])->cost)
        //compare_paths(queue->buffer[node], queue->buffer[left_child])
            )
	{
		i = left_child;
	}
	
	// Compare with the right node
	if (right_child < queue->size &&
            (((struct Tour*) queue->buffer[i])->cost < ((struct Tour*)(queue->buffer[right_child]))->cost)

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
            (((struct Tour*)(queue->buffer[parent_of(node)]))->cost < ((struct Tour*) queue->buffer[node] )->cost))
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
    //
    //printf("Solution found, trimming the tree....\nQUEUE SIZE:%zu\n\n",queue->size);
    int i = 1;
    print_tour((struct Tour*)queue->buffer[i]);

    //queue_print(queue, stdout, pq);
    for(int j = 0; j< queue->size; j++){
        if(((struct Tour*)queue->buffer[j])->cost >= maxCost){
            free_tour(queue->buffer[j]);
            remove_element(queue, j);
        }
    }
}


