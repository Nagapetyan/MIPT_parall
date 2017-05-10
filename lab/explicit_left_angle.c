#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


static inline int assert (int cond, const char *message);
static inline void Free(void *buf);
static inline void free_2d(int **buf, const int size);
static inline void init_2d(int ***buf, const int size);


void matrixGenerator(int  **firstMatrix, int d1, int d2);


int const root = 0;

int main(int argc, char **argv){
	int d1 = strtol(argv[1], NULL, 10);
	int d2 = strtol(argv[2], NULL, 10);
	int tau = strtol(argv[3], NULL, 10);
	int h = strtol(argv[4], NULL, 10);

	MPI_Init(NULL, NULL);
	

	return 0;
}


static inline int assert (int cond, const char *message){ 
	if (!cond){ 
	 	fprintf(stderr, "Error in %s in  %d\n", message, __LINE__); 
	 	abort(); 
	}
	return 0;
}


static inline void Free(void *buf){
	free(buf);
	buf = NULL;
}

static inline void free_2d(int **buf, const int size){
	for(int i = 0; i < size; ++i)
		Free(buf[i]);
	Free(buf);
}

static inline void init_2d(int ***buf, const int size){
	int **data = (int**)calloc(size, sizeof(int*));
	assert(data != NULL, "Bad allocation");
	for(int i = 0; i < size; i++){
		data[i] = (int*)calloc(size, sizeof(int));
		assert(data[i] != NULL, "Bad allocation");
	}
	*buf = data;
}


void matricesGenerator(int **matrix, int d1, int 2d){
	srand(time(NULL));
	int min = INT_MAX;
	for(int i = 0; i < d1; ++i){
		for (int j = 0; j < d2; ++j){
			matrix[i][j] = rand();
			if (matrix[i][j]  < min)
				min = matrix[i][j];
		}
	}
	for (int i = 0; i < d1; ++i){
		for(int j = 0; j < d2; ++j){
			matrix[i][j] = (int)(matrix[i][j]/min);
		}
	}
}