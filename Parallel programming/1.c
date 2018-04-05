#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

#define FUNC(x) (1.0 / (1.0 + x *  x))

typedef struct
{
	int begin, end, num;
	double result;
}block;

double calc(int begin, int end, int num);

int main_proc(int size)
{
	double start_time = 0.0, end_time = 0.0;
	int num = 0, i = 0;
	printf("num = \n");
	scanf("%d", &num);

	start_time = MPI_Wtime();
	double result = calc(0, num, num) * 4.0;
	end_time = MPI_Wtime();
	printf("main_proc result = %lf; elapsed time = %lf\n", result, end_time - start_time);


	start_time = MPI_Wtime();
	block * b_list = calloc(size, sizeof(block));
	
	int block_sum = 0;
	for(i = 0; i < size; ++i)
	{
		b_list[i].num = num / size;
		if(i < num % size)
			b_list[i].num += 1;
		b_list[i].begin = block_sum;
		block_sum += b_list[i].num;
		b_list[i].end = block_sum;
		b_list[i].num = num;
	}

	for(i = 1; i < size; ++i)
	{
		MPI_Send(&b_list[i], sizeof(block), MPI_BYTE, i, 0, MPI_COMM_WORLD);
	}

	b_list[0].result = calc(b_list[0].begin, b_list[0].end, b_list[0].num);
	MPI_Status stat = {};
	for(i = 1; i < size; ++i)
	{
		MPI_Recv(&b_list[i], sizeof(block), MPI_BYTE, i, 0, MPI_COMM_WORLD, &stat);
		b_list[0].result += b_list[i].result;
	}
	b_list[0].result *= 4.0;	
	end_time = MPI_Wtime();
	printf("multiple proc result = %lf; elapsed time = %lf\n", b_list[0].result, end_time - start_time);

	free(b_list);
	return 0;
}

int proc(int rank, int size)
{
	block blk = {};
	MPI_Status stat = {};
	MPI_Recv(&blk, sizeof(block), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &stat);

	blk.result = calc(blk.begin, blk.end, blk.num);

	MPI_Send(&blk, sizeof(block), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
	
	printf("process %d sent result = %lf\n", rank, blk.result);
	return 0;
}


int main(int argc, char * argv[])
{
	MPI_Init(&argc, &argv);
	int rank = 0, size = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int ret = 0;
	if(rank == 0)
		ret = main_proc(size);
	else
		ret = proc(rank, size);
	
	MPI_Finalize();
	return ret;
}

double calc(int begin, int end, int num)
{
	if(begin >= end)
		return 0.0;
	int i = 0;
	double x = 0.0;
	double step = 1.0 / (double)num;
	double result = FUNC((double)begin * step) + FUNC((double)end * step);
	
	for(i = begin + 1; i < end; ++i)
	{
		x = (double)i * step;
		result += FUNC(x) * 2.0;
	}
	result *= step / 2.0;

	return result;
}


