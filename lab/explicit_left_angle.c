#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static inline void assert (int cond, const char *message);
static inline void Free(void *buf);
static inline double next_step(double r, double u_prevT_currX, double u_prevT_prevX);

int const root = 0;


int main(int argc, char **argv){
	assert(argc == 3 || argc == 5, "Wrong number of arguments. Please enter number of steps by x and t");
	double X = 1, T = 1;
	if(argc == 5){
		X = strtol(argv[3], NULL, 10);
		T = strtol(argv[4], NULL, 10);
	}

	int x_num_steps = strtol(argv[1], NULL, 10);
	int t_num_steps = strtol(argv[2], NULL, 10);
	double h = X/x_num_steps;
	double tau = T/t_num_steps;

	double r = tau/h;

	MPI_Init(NULL, NULL);
	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	MPI_Barrier(MPI_COMM_WORLD);
	double start_time;
	if(rank == root)
		start_time = MPI_Wtime();

	int width_by_x = x_num_steps / num_proc;
	int width_by_t = t_num_steps / num_proc;

	double **u_local;

	if(rank == root){
		u_local =  (double**)calloc(x_num_steps, sizeof(double*));
		for(int i = 0; i < x_num_steps; ++i)
			u_local[i] = (double*)calloc(width_by_t, sizeof(double));

		for(int i = 0; i < x_num_steps; ++i)
			u_local[i][0] = i*h;
		for(int i = 0; i < width_by_t; ++i)
			u_local[0][i] = i*tau;

		for(int iter = 0; iter < num_proc; iter++){
			for(int i = 1*(iter+1); i < width_by_x*(iter+1); ++i)
				for(int j = 1; j < width_by_t; ++j)
					u_local[i][j] = next_step(r, u_local[i][j-1], u_local[i-1][j-1]);


			if(rank != num_proc-1){
				double *send_data = (double*)calloc(width_by_x, sizeof(double));
				assert(send_data != NULL, "Couldn't allocate memory for send data buffer in root process");
				for(int i = iter*width_by_x; i < width_by_x*(iter+1); ++i)
					send_data[i-iter*width_by_x] = u_local[i][width_by_t-1];
				MPI_Send(send_data, width_by_x, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
				Free(send_data);
			}
		}

	}else{
		u_local =  (double**)calloc(x_num_steps, sizeof(double*));
		for(int i = 0; i < x_num_steps; ++i)
			u_local[i] = (double*)calloc(width_by_t+1, sizeof(double));

		u_local[0][1] = tau*width_by_t*rank;
		for(int i = 2; i < width_by_t+1; ++i)
			u_local[0][i] =u_local[0][i-1] + tau;

		for(int iter = 0; iter < num_proc; ++iter){
			double *receive_data = (double*)calloc(width_by_x, sizeof(double));
			assert(receive_data != NULL, "Couldn't allocate memory for receive data buffer");
			MPI_Recv(receive_data, width_by_x, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int i = 0; i < width_by_x; ++i)
				u_local[i+ iter*width_by_x][0] = receive_data[i];

			for(int i = width_by_x*iter; i < width_by_x*(iter+1); ++i)
				for(int j = 1; j < width_by_t+1; ++j)
					if(i != 0)
						u_local[i][j] = next_step(r, u_local[i][j-1], u_local[i-1][j-1]);

			if(rank != num_proc-1){
				double *send_data = (double*)calloc(width_by_x, sizeof(double));
				assert(send_data != NULL, "Couldn't allocate memory for send data buffer in non root process");
				for(int i = iter*width_by_x; i < width_by_x*(iter+1); ++i)
					send_data[i-iter*width_by_x] = u_local[i][width_by_t];
				MPI_Send(send_data, width_by_x, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
				Free(send_data);
			}
			Free(receive_data);
		}
	}


	int barrier=0;
	FILE *f_write=fopen("output.txt","a");
	assert(f_write != NULL, "Couldn't open file");
	if(rank == root){
		for(int i = 0; i < width_by_t; ++i){
			for(int j = 0; j < x_num_steps; ++j)
				fprintf(f_write, "%lf ", u_local[j][i]);
			fprintf(f_write, "\n");
		}
		fflush(NULL);
		if(rank != num_proc-1)
			MPI_Send(&barrier, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	}else{
		MPI_Recv(&barrier, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for(int i = 1; i < width_by_t+1; ++i){
				for(int j = 0; j < x_num_steps; ++j)
					fprintf(f_write, "%lf ", u_local[j][i]);
				fprintf(f_write, "\n");
			}
		fflush(NULL);
		if(rank != num_proc-1)
			MPI_Send(&barrier, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	}


	if(rank == root)
		printf("Execution time: %lf\n", MPI_Wtime() - start_time);


	if(rank == root){
		for(int i = 0; i < width_by_t; ++i)
			Free(u_local[i]);
		Free(u_local);
	}else{
		for(int i = 0; i < width_by_t+1; ++i)
			Free(u_local[i]);
		Free(u_local);
	}

	MPI_Finalize();
	return 0;
}


static inline void assert (int cond, const char *message){ 
	if (!cond){ 
	 	fprintf(stderr, "Error in %s \n", message); 
	 	abort(); 
	}
}


static inline void Free(void *buf){
	free(buf);
	buf = NULL;
}

static inline double next_step(double r, double u_prevT_currX, double u_prevT_prevX){
	return (1-r)*u_prevT_currX + r* u_prevT_prevX;
}

