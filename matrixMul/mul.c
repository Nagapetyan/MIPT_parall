#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <stddef.h> 

static inline int ASSERT (int cond, const char *message){ 
	if (!cond){ 
	 	fprintf(stderr, "Error in %s in  %d\n", message, __LINE__); 
	 	abort(); 
	}
	return 0;
}

typedef struct Cartesian_grind_info{
	MPI_Comm grid_comm;
	MPI_Comm rowComm;
	MPI_Comm colComm;
	int num_proc;
	int rank;
	int col, row;
	int grid_size;
} grid_info_t;

static inline void Free(void *buf);
static inline void free_2d(int **buf, const int size);
static inline void init_2d(int ***buf, const int size);

void matricesGenerator(int  **firstMatrix, int **secondMatrix, const int size);
void createCartesianTopology(grid_info_t *grid_info);
void foxMultiply(grid_info_t *grid_info, int **local_a, int **local_b, int **local_c, int size);
void packageMatrix(int **matrix, int *package, int const dimension_size);
void unpackageMaxtix(int *package, int **matrix, int const dimension_size);
void localMultiply(int **a, int **b, int  **c, int const size);

int CHECK(int **true, int *res, int size);

int const root = 0;

int main(int argc, char **argv){
	int size = strtol(argv[1], NULL, 10);
	ASSERT(size > 0, "Wrong size");

	int **firstMatrix, **secondMatrix;
	init_2d(&firstMatrix, size);
	init_2d(&secondMatrix, size);
	matricesGenerator(firstMatrix, secondMatrix, size);

	MPI_Init(NULL, NULL);

	grid_info_t grid_info; 
	createCartesianTopology(&grid_info);

	int block_size = size/grid_info.grid_size;
	int base_row = grid_info.row * block_size;
	int base_col = grid_info.col * block_size;

	int **local_a, **local_b, **local_c;
	init_2d(&local_a, block_size);
	init_2d(&local_b, block_size);
	init_2d(&local_c, block_size);

	for (int i = base_row; i < base_row + block_size; i++){
		for (int j = base_col; j < base_col + block_size; j++) {
			local_a[i - (base_row)][j - (base_col)] = firstMatrix[i][j];
			local_b[i - (base_row)][j - (base_col)] = secondMatrix[i][j];
		}
	}

	double start_time;
	MPI_Barrier(MPI_COMM_WORLD);
	if (grid_info.rank == root)
		start_time = MPI_Wtime();

	foxMultiply(&grid_info, local_a, local_b, local_c, size);

	free_2d(local_a, block_size);
	free_2d(local_b, block_size);

	int **true;
	init_2d(&true, size);
	localMultiply(firstMatrix, secondMatrix, true, size);

	free_2d(firstMatrix, size);
	free_2d(secondMatrix, size);

	MPI_Barrier(MPI_COMM_WORLD);
	if(grid_info.rank == root)
		printf("Execution time: %lf\n", MPI_Wtime() - start_time);

	int *result = (int*)calloc(block_size*block_size, sizeof(int));
	ASSERT(result != NULL, "Bad allocation of local result package");
	packageMatrix(local_c, result, block_size);
	free_2d(local_c, block_size);

	int *totalresult = (int*)calloc(size*size, sizeof(int));
	ASSERT(totalresult != NULL, "Bad allocation of total result package");

	MPI_Gather(result, block_size*block_size, MPI_INT, totalresult, block_size*block_size, MPI_INT, root, grid_info.grid_comm);

	if (grid_info.rank == root){
		int *data = (int*)calloc(size*size, sizeof(int));
		ASSERT(data != NULL, " ");
		        
		int k = 0;
		for (int bi = 0; bi < grid_info.grid_size; bi++)
			for (int bj = 0; bj < grid_info.grid_size; bj++)
				for (int i = bi * block_size; i < bi * block_size + block_size; i++)
					for (int j = bj * block_size; j < bj * block_size + block_size; j++) {
						data[i * size + j] = totalresult[k];
						k++;
					}
		int check;
		check = CHECK(true, data, size);
		if(check == 0)
			printf("Smth wrong with program\n");
		else{
			printf("ALL OKAY\n");
		}
	}


	MPI_Finalize();
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
	ASSERT(data != NULL, "Bad allocation");
	for(int i = 0; i < size; i++){
		data[i] = (int*)calloc(size, sizeof(int));
		ASSERT(data[i] != NULL, "Bad allocation");
	}
	*buf = data;
}

void matricesGenerator(int **A, int **B, const int size){
	srand(time(NULL));
	int min = INT_MAX;
	for(int i = 0; i < size; ++i){
		for (int j = 0; j < size; ++j){
			A[i][j] = rand();
			if (A[i][j]  < min)
				min = A[i][j];
			B[i][j] = rand();
			if (B[i][j] < min)
				min = B[i][j];
		}
	}
	for (int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			A[i][j] = (int)(A[i][j]/min);
			B[i][j] = (int)(B[i][j]/min);
		}
	}
}

void createCartesianTopology(grid_info_t *grid_info){
	MPI_Comm_rank(MPI_COMM_WORLD, &(grid_info->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(grid_info->num_proc));
	ASSERT(grid_info->num_proc %2 == 0, "Wrong num of proc");

	grid_info->grid_size = (int) sqrt(grid_info->num_proc);
	ASSERT((grid_info->grid_size * (grid_info->grid_size)) == (grid_info->num_proc), "Wrong num of proc");


	int dims[2], periods[2], reorder = 0;
	dims[0] =  dims[1] = grid_info->grid_size;
	periods[0] = periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &(grid_info->grid_comm));

	int coords[2];
	MPI_Cart_coords(grid_info->grid_comm, grid_info->rank, 2, coords);
	grid_info->row = coords[0]; 
	grid_info->col = coords[1];

	int subdims[2];
	subdims[0] = 0;  
	subdims[1] = 1; 
	MPI_Cart_sub(grid_info->grid_comm, subdims, &(grid_info->rowComm));
	subdims[0] = 1;
	subdims[1] = 0;
	MPI_Cart_sub(grid_info->grid_comm, subdims, &(grid_info->colComm));
}


void foxMultiply(grid_info_t *grid_info, int **local_a, int **local_b, int **local_c, int size){
	int **tmp_a;
	int stage;
	int bcast_root;
	int tmp_size;
	int source;
	int dest;

	tmp_size = size / grid_info->grid_size;
	source = (grid_info->row + 1) % (grid_info->grid_size);
	dest = (grid_info->row + grid_info->grid_size  - 1) % (grid_info->grid_size);

	init_2d(&tmp_a, tmp_size);

	int *package = (int*)calloc(tmp_size*tmp_size, sizeof(int));
	ASSERT(package != NULL, "Bad package allocation in foxMultiply");


	for(stage = 0; stage < grid_info->grid_size; ++stage){
		bcast_root = (grid_info->row + stage) % grid_info->grid_size;
		if (bcast_root == grid_info->col){
			packageMatrix(local_a, package, tmp_size);
			MPI_Bcast(package, tmp_size*tmp_size, MPI_INT, bcast_root, grid_info->rowComm);
			unpackageMaxtix(package, local_a, tmp_size);   
			localMultiply(local_a, local_b, local_c, tmp_size);
		} else{
			packageMatrix(tmp_a, package, tmp_size);
			MPI_Bcast(package, tmp_size*tmp_size, MPI_INT, bcast_root, grid_info->rowComm);
			unpackageMaxtix(package, tmp_a, tmp_size);
			localMultiply(tmp_a, local_b, local_c, tmp_size);
		}

		packageMatrix(local_b, package, tmp_size);
		MPI_Sendrecv_replace(package, tmp_size*tmp_size, MPI_INT, dest, 0, source, 0, grid_info->colComm, MPI_STATUS_IGNORE);
		unpackageMaxtix(package, local_b, tmp_size);
	}
}

void packageMatrix(int **matrix, int *package, int const dimension_size){
	for(int i = 0; i < dimension_size; ++i){
		for(int j = 0; j < dimension_size; ++j)
			package[j + dimension_size*i] = matrix[i][j];
	}
}

void unpackageMaxtix(int *package, int **matrix, const int dimension_size){
	for(int i = 0; i < dimension_size; ++i){
		for(int j = 0; j < dimension_size; ++j)
			matrix[i][j] = package[i * dimension_size + j];
	}
}

void localMultiply(int **a, int **b, int  **c, int const size){
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			for(int k = 0; k < size; ++k)
				c[i][j]  += a[i][k] * b[k][j];
		}
	}
}


int CHECK(int **true, int *res, int size){
	int check = 1;
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++){
			if(true[i][j]!=res[i*size+j]){
			check = 0;
		}
	}
	return check;
}
