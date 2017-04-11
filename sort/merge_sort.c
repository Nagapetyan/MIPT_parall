#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int const root = 0;

int *mergeSort(int *data, int *buf, int left, int right);
void merge(int *data, int *left, int left_size, int *right, int right_size);
void parallelMerge(int *data, int size, int height, int rank, int num_proc);
void fillData(int *data, int size);

static inline void assert(int cond, const char *message);
static inline void Free(void *buf);

int log_2(int x);

int CHECK(int *check, int *sorted_data, int size);

int main(int argc, char **argv){
	assert( (argc == 2), "Wrong input number");
	int size = strtol(argv[1], NULL, 10);
	assert(size > 0, "Wrong size");

	int rank, num_proc;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	MPI_Barrier(MPI_COMM_WORLD);
	double start_time;
	if(rank == root)
		start_time = MPI_Wtime();

//	assert(size % (num_proc) == 0, "Meet the requirement: size of data mod num_proc = 0");

	if(rank == root){
		int *unsorted_data = (int*)calloc(size, sizeof(int));
		assert(unsorted_data != NULL, "Bad allocation of initial data");
		fillData(unsorted_data, size);
		
		int *check = (int*)calloc(size, sizeof(int));
		memmove(check, unsorted_data, size*sizeof(int));

		int root_Ht = log_2(num_proc);
		parallelMerge(unsorted_data, size, root_Ht, rank, num_proc);


		if(rank == root)
			fprintf(stdout, " Execution time: %lf\n",  MPI_Wtime() - start_time);

		int control_flag = CHECK(check, unsorted_data, size);
		if(control_flag == 1)
			fprintf(stdout, "OK\n");
		else
			fprintf(stdout, "FEELSBADMAN\n");


		Free(unsorted_data);
		Free(check);

		MPI_Finalize();
	}else{
		int info[2];
		MPI_Recv(info, 2, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		int *receive_data = (int*)calloc(info[0], sizeof(int));
		assert(receive_data != NULL, "Bad allocation of array for receive data in main");

		MPI_Recv(receive_data, info[0], MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
		parallelMerge(receive_data, info[0], info[1], rank, num_proc);

		Free(receive_data);

		MPI_Finalize();
	}
	return 0;
}

void parallelMerge(int *data, int size, int height, int rank, int num_proc){
	int parent = rank & ~(1 << height);
	int next = height - 1;
	int child = rank | (1 << next);

	if(height > 0){
		if(child >= num_proc){
			parallelMerge(data, size, next, rank, num_proc);
		}
		else{
			int left_size = size/2;
			int right_size = size - left_size;

			int *left_array = (int*) calloc(left_size, sizeof(int));
			assert(left_array != NULL, "Bad alloc left arr");
			
			int *right_array = (int*)calloc(right_size, sizeof(int));
			assert(right_array != NULL, "Bad alloc right arr");

			memmove(left_array, data, left_size*sizeof(int));
			memmove(right_array, data + left_size, right_size*sizeof(int));

			int info[2];
			info[0] = right_size;
			info[1] = next;

			MPI_Send(info, 2, MPI_INT, child, 0, MPI_COMM_WORLD);
			MPI_Send(right_array, right_size, MPI_INT, child, 0, MPI_COMM_WORLD);

			parallelMerge(left_array, left_size, next, rank, num_proc);

			MPI_Recv(right_array, right_size, MPI_INT, child, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			merge(data, left_array, left_size, right_array, right_size);
		}	
	}else{

		int *tmp = (int*)calloc(size, sizeof(int));
		tmp = mergeSort(data, tmp, 0, size-1);
		memmove(data, tmp, size*sizeof(int));
	}


	if(parent != rank){
		MPI_Send(data, size,  MPI_INT, parent, 0, MPI_COMM_WORLD);
	}
}

void merge(int *data, int *left, int left_size, int *right, int right_size){
	int rc = 0, lc = 0, count = 0;
	while(rc < right_size && lc < left_size){

		data[count++] = right[rc] < left[lc] ? right[rc++] : left[lc++];
	}
	while(rc < right_size)
		data[count++] = right[rc++];
	while(lc < left_size)
		data[count++] = left[lc++];
}


int log_2(int x){
	int count = 1, ret = 0;
	while(count < x){
		count+= count;
		ret++;
	}
	return ret;
}

int *mergeSort(int *data, int *buf, int left, int right){
	

	if(left == right){
		buf[left] = data[left];
		return buf;
	}

	int  middle = (left + right) / 2;
	int *l_buf = mergeSort(data, buf, left, middle);
	int *r_buf = mergeSort(data, buf, middle+1, right);

	int *target = l_buf == data ? buf : data;
	int l_cur = left, r_cur = middle + 1;

	for (int i = left; i <= right; i++){
		if (l_cur <= middle && r_cur <= right){
			if (l_buf[l_cur] < r_buf[r_cur]){
				target[i] = l_buf[l_cur];
				l_cur++;
			}else{
				target[i] = r_buf[r_cur];
				r_cur++;
			}
		}
		else if (l_cur <= middle){
			target[i] = l_buf[l_cur];
			l_cur++;
		}
		else{
			target[i] = r_buf[r_cur];
			r_cur++;
		}
	}
	return target;
}

void fillData(int *data, int size){
	srand(time(NULL));
	for(int i = 0; i < size; ++i)
		data[i] = rand();
}


static inline void assert(int cond, const char *message){
	if(!cond){
		fprintf(stderr, "Error: %s \n", message); 
	 	abort(); 
	}
}

static inline void Free(void *buf){
	free(buf);
	buf = NULL;
}


int CHECK(int *check, int *sorted_data, int size){
	int control_flag = 1;
	int *tmp = (int*)calloc(size, sizeof(int));
	check = mergeSort(check, tmp, 0, size-1);
	for(int i = 0; i < size; i++)
		if (sorted_data[i] != check[i])
			control_flag = 0;
	return control_flag;
}


/*******************************************************************************************************************************************************************/
/* 						Process Ranks in a Four-Level Sorting Tree									*    /


/*							       NODE 0							
								
				NODE0							NODE4

		NODE0			NODE2					NODE4			NODE6

	NODE0		NODE1	    NODE2	NODE3			NODE4		NODE5	    NODE6	NODE7

																					*/