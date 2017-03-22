#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define  n 10

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

	int check;
	check = fread(fisrstNumber, n,1, inputData);
//	if(!check) {
//		fprintf(stderr, "Unable to read first \n" );
//		exit(0);
//	}
	char tmp;
	fread(&tmp, 1,1,inputData);
	check = fread(secondNumber, n, 1, inputData);
//	if(!check) {
//		fprintf(stderr, "Unable to read second\n" );
//		exit(0);
//	}

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
 
 //	printf("%d\n", counts);
  // 	for (int i = 0; i < counts; i++)
   //		printf("index[%d] =%d\n",i, index[i]);
   	determineBounds(rank, num_proc, counts, &bounds, index);	
   	
   	char *result = (char*)calloc(counts, sizeof(char));
   	assert(result != NULL);

   //	printf("rank %d  up %d   down   %d\n", rank,  bounds.up, bounds.down);

   	int carry;
   	carry = digitAddition(bounds, fisrstNumber, secondNumber, result);

  //	printf("rank %d  %d  %s\n",rank,  (int)strlen(result), result);
   //	if(rank == num_proc-1)
   //		result[n] = carry;

   
   
   //	printf("%s\n+\n%s\n=\n%d%s\n",fisrstNumber, secondNumber,  carry, result);

   
   //	printf("rank = %d  result=%d%s\n",rank, carry,result );

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
   	//	printf("rank=%d result:%s\n", rank, result);
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
				//	printf("%d and i = %d\n", count,i);

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
	/*
	else{
		if( rank < count){
			if (rank == 0)
				bounds->down = n-1;
			else
				bounds->down = index[count - rank];
			bounds->up = index[count - rank -1] + 1;
		}	
		else{
			if(rank == num_proc - 1)
				bounds->up = 0;
			else
				bounds->up = index[0]  * (num_proc - 1 - rank)/ (num_proc - count) + 1;
			bounds->down = index[0] * (num_proc -  rank)/(num_proc - count);

		}
	}
	*/
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
//		printf("fisrstNumber[%d]=%c secondNumber[%d]=%c\n", i,fisrstNumber[i], i, secondNumber[i]);
		result[i - bounds.up] = toChar(toDigit(fisrstNumber[i]) ^ toDigit(secondNumber[i]));
	//	printf("first   %c\n", result[i]);
		result[i  - bounds.up] = toChar(toDigit(result[i  - bounds.up]) ^ iCarry);
	//	printf("second   %c\n", result[i- bounds.down]);
		iCarry = (toDigit(fisrstNumber[i]) * toDigit(secondNumber[i])) || (toDigit(fisrstNumber[i]) * iCarry) || (toDigit(secondNumber[i] )* iCarry);
	//	printf("%d\n", iCarry);
	}
//	printf("%s\n", result);
//	printf("\n");
	return iCarry;
}



static inline void Free(void *buf)
{
	free(buf);
	buf = NULL;
}

