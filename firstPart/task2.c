#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define  n 2000

typedef struct Bounds
{
	int up;
	int down;
} bound_t;

static inline void Free(void *buf);
static inline int toDigit(char c) {return c - '0';}
static inline char toChar(int c) {return c + '0';} 


void determineBounds(int rank, int num_proc, int count, bound_t *bounds, const int *index);

int digitAddition(const bound_t bounds, const char *fisrstNumber, const char *secondNumber, char *rankSum);

int countingOverflow(const char *fisrstNumber, const char *secondNumber, int *index);

int main(int argc, char **argv)
{
	char *fisrstNumber = (char*) calloc(n, sizeof(char));
	assert(fisrstNumber != NULL);
	char *secondNumber = (char*) calloc(n, sizeof(char));
	assert(secondNumber != NULL);

	FILE *inputData = fopen("dataForTask2.txt", "r");
	assert(inputData != NULL);

	fread(fisrstNumber, n,1, inputData);

	char tmp;
	fread(&tmp, 1,1,inputData);
	fread(secondNumber, n, 1, inputData);

	fclose(inputData);


	MPI_Init(NULL, NULL);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int num_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);


	bound_t bounds;
   	int  *index = (int*) calloc(n, sizeof(int));
   	assert(index != NULL);
   	int counts = countingOverflow(fisrstNumber, secondNumber, index);
 
   	determineBounds(rank, num_proc, counts, &bounds, index);	
   	
   	char *result = (char*)calloc(counts, sizeof(char));
   	assert(result != NULL);

   	int carry;
   	carry = digitAddition(bounds, fisrstNumber, secondNumber, result);

   	FILE *answer = fopen("answer.txt", "a");

   	if (rank == num_proc - 1) {
   		fprintf(answer, "%d", carry);
   		fprintf(answer, "%s", result);	
   	}

	if (rank != num_proc - 1)
   	{
   		int control;
   		MPI_Recv(&control, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   		fprintf(answer, "%s", result);
   	}
   	if (rank != 0)
   	{
   		int control = 1;
   		MPI_Send(&control, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
   	}
   	

	MPI_Finalize();


	Free(fisrstNumber);
	Free(secondNumber);
	Free(result);
	return 0;
}

int countingOverflow(const char *fisrstNumber, const char *secondNumber, int *index){
	int count = 0;
	for (int i = 0; i < n ; i++)
	{
		if (i != n-1){
			if ( ((toDigit(fisrstNumber[i]) * toDigit(secondNumber[i])) == 1) && ((toDigit(fisrstNumber[i+1]) * toDigit(secondNumber[i+1])) != 1)){
					index[count] = i;
					count++;
				}
		}
		else if ((toDigit(fisrstNumber[i]) * toDigit(secondNumber[i])) == 1){
			index[count] = i;
			count++;
		}
	}	
	return count;
}

void determineBounds(int rank, int num_proc, int count, bound_t *bounds, const int *index){
	if (count - num_proc >= -1){
		if(rank == 0)
			bounds->down = n-1;
		else 
		{
			bounds->down = index[count*(num_proc - rank) / num_proc];
		}
		if (rank == num_proc-1)
			bounds->up = 0;
		else
			bounds->up = index[count*(num_proc - rank - 1) / num_proc] + 1;
	}
	else{
		if (num_proc - rank <= count){
			if (rank == num_proc - 1)
				bounds->up = 0;
			else
				bounds->up = index[num_proc - rank - 2] ;
			bounds->down = index[num_proc - rank-1] ;
		}
		else{
			if (rank == 0)
				bounds->down = n-1;
			else
				bounds->down = n -1 -  (n - index[count-1]) *rank / (num_proc - count);
			if(num_proc - rank == count+1)
				bounds->up = index[count-1] + 1;
			else
				bounds->up = n  -  (n - index[count-1]) *(rank+1) / (num_proc - count) ;
		}
	}
}



int digitAddition(const bound_t bounds, const char *fisrstNumber, const char *secondNumber, char *result){
	int iCarry = 0;
	for(int i = bounds.down; i > bounds.up-1; i--){
		result[i - bounds.up] = toChar(toDigit(fisrstNumber[i]) ^ toDigit(secondNumber[i]));
		result[i  - bounds.up] = toChar(toDigit(result[i  - bounds.up]) ^ iCarry);
		iCarry = (toDigit(fisrstNumber[i]) * toDigit(secondNumber[i])) || (toDigit(fisrstNumber[i]) * iCarry) || (toDigit(secondNumber[i] )* iCarry);
	}
	return iCarry;
}



static inline void Free(void *buf){
	free(buf);
	buf = NULL;
}

