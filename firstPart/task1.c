#include <mpi.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stddef.h>


const int root = 0;

typedef struct Bounds{
	int down;
	int up;
}bound_t;

static inline int nDigits(int number)
{
	return (number > 0) ? (floor(log10(abs(number))) + 1)  :  (number < 0) ? (floor(log10(abs(number))) + 2) : 1;
}

static inline void Free(void *buf)
{
	free(buf);
	buf = NULL;
}

void determineBounds(int rank, int num_proc, bound_t *bounds, int  n)
{
	if (rank == 0)
   		bounds->down = 1;
   	else 
   		bounds->down = rank * n / num_proc + 1;

   	if (rank == num_proc -1)
   		bounds->up = n;
   	else
   		bounds->up = (rank+1)* n/num_proc ;
}

void rankResult(const bound_t *bounds, mpf_t result, mpf_t current_number)
{
	mpf_t one;
	mpf_init_set_str(one, "1", 10);


   	for (int  i = bounds->down; i < bounds->up+1; i++)
   	{
   		mpf_t tmp;

   		char *tmpString = (char*)calloc(nDigits(i), sizeof(char));   
   		assert(tmpString != NULL);		
   		sprintf(tmpString, "%d", i);             
   		
   		mpf_init_set_str(tmp,  tmpString, 10);
   		
   		mpf_mul(current_number, current_number, tmp);
   		mpf_div(tmp, one, current_number);
   		mpf_add(result, result, tmp);
   		
   		mpf_clear(tmp);
   		Free(tmpString);
   	}

   	mpf_clear(one);
}



void Receive(int rank, mpf_t receivedData, int key)
{
	MPI_Status status;

	int len;
	if(key == 1)
		MPI_Recv(&len, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, &status);
	if(key == 0)
		MPI_Recv(&len, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);

	char *recvMessage = (char*) calloc(len, sizeof(char));
	if(key == 1)
		MPI_Recv(recvMessage, len, MPI_CHAR, rank, 0, MPI_COMM_WORLD, &status);
	if(key == 0)
		MPI_Recv(recvMessage, len, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, &status);

	mpf_init_set_str(receivedData, recvMessage, 10);

	Free(recvMessage);

}




void Send(mpf_t data, int rank, int key)
{
	char *data_stringType;
   	mp_exp_t exp;
   	data_stringType = mpf_get_str(NULL, &exp, 10, 0, data);
   	exp  = exp - (int)strlen(data_stringType);
   	char *exp_stringType = (char*)calloc(nDigits(exp), sizeof(char));
   	assert(exp_stringType != NULL);
   	sprintf(exp_stringType, "%d", exp);

   	char *message = (char*)calloc(strlen(data_stringType) + nDigits(exp)+1, sizeof(char));
   	assert(message != NULL);
   	memmove(message, data_stringType, strlen(data_stringType));
	message[strlen(data_stringType)] = 'e';
	memmove(message+strlen(data_stringType)+1, exp_stringType, strlen(exp_stringType));

	int len = (int)strlen(message);

	if (key == 0){
		MPI_Send(&len, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
		MPI_Send(message, strlen(message), MPI_CHAR, rank+1, 0, MPI_COMM_WORLD);
	}
	else{
		MPI_Send(&len, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
		MPI_Send(message, strlen(message), MPI_CHAR, root, 0, MPI_COMM_WORLD);
	}

	Free(message);
	Free(data_stringType);
	Free(exp_stringType);
}


int main(int argc, char** argv)
{
	assert(argc == 2);
	int n = strtol(argv[1], NULL, 10);
	assert(n>=0);

	MPI_Init(NULL, NULL);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int num_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

   	MPI_Bcast(&n, 1, MPI_INT, 0,  MPI_COMM_WORLD);

   	bound_t bounds;
   	determineBounds(rank, num_proc, &bounds, n);

   	mpf_t result, totalResult ;
   	mpf_init(result);
   	mpf_init(totalResult);
   	mpf_t currentNumber;
   	mpf_init_set_str(currentNumber, "1", 10);


   	rankResult(&bounds, result, currentNumber);


   	if(rank < num_proc - 1)
   		Send(currentNumber, rank,0);
   	if(rank > 0)
   	{
   		mpf_t ratio;
   		Receive(rank, ratio, 0);
   		mpf_mul(currentNumber, currentNumber, ratio);
   		mpf_div(result, result, ratio);
   		mpf_clear(ratio);
   	}




   	if (rank != 0){
   		Send(result, rank,1);
   	}
   	else{
   		for (int i = 1; i < num_proc; ++i)
   		{
   			mpf_t result_iRank;
   			mpf_init(result_iRank);
   			Receive(i, result_iRank,1);
   			mpf_add(totalResult, totalResult, result_iRank);
   			mpf_clear(result_iRank);

   		}

   		mpf_add(totalResult, totalResult, result);
   		gmp_printf("%Ff\n", totalResult);
   		fflush(NULL);
   	}

   	mpf_clear(result);
   	mpf_clear(totalResult);

   	MPI_Finalize();
	return 0;
}