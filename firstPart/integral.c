#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

static inline double func(double x){return x*x;}
static inline void assert(int cond, const char *message);

int const up_limit = 1;
int const down_limit = 0;
int const root = 0;

double integral(double step, int rank, int num_proc);

int main(int argc, char **argv){
	assert(argc == 1 || argc == 2, "Wrong number of arguments");
	double step = 0.0001;
	if(argc == 2)
		step = strtod(argv[1], NULL);

	MPI_Init(NULL,NULL);
	
	int rank, num_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Barrier(MPI_COMM_WORLD);
	double start_time;
	if(rank == root)
		start_time = MPI_Wtime();
	

	double res = integral(step, rank, num_proc);
	double total_res = 0;
	MPI_Reduce(&res, &total_res, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == root)
		printf("Execution time: %lf\n", MPI_Wtime() - start_time);

	MPI_Finalize();
	return 0;
}

static inline void assert(int cond, const char *message){
	if (!cond){ 
	 	fprintf(stderr, "Error in %s in  %d\n", message, __LINE__); 
	 	abort(); 
	}
}


double integral(double step, int rank, int num_proc){
	long unsigned num_of_steps = (up_limit - down_limit)/step;
	double down_edge = rank * up_limit * 1./num_proc;
	double up_edge = (rank +1)  * 1. * up_limit / num_proc;
	num_of_steps /= num_proc;

	double res = 0;
	for(long unsigned i = 0; i < num_of_steps; ++i)
		res += func(down_edge + i * step);
	res += (up_edge + down_edge)/2;
	res *= step;

	return res;
}