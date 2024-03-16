#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void getArgs(int argc, char **argv, char **in, char **out, int *steps)
{
	if (argc < 4) 
    {
		printf("Not enough paramters:  ./homework IN_FILENAME OUT_FILENAME NUM_STEPS\n");
		exit(1);
	}
	(*in) = (char *)malloc(strlen(argv[1])* sizeof(char));
    memcpy((*in), argv[1], strlen(argv[1]));

    (*out) = (char *)malloc(strlen(argv[2])* sizeof(char));
    memcpy((*out), argv[2], strlen(argv[2]));

    (*steps) = atoi(argv[3]);
}

void allocateContinuousMatrix(int ***matrix, int rows, int cols) 
{
    // Allocate memory for array of pointers (rows)
    *matrix = (int **)malloc(rows * sizeof(int *));
    if (*matrix == NULL) {
        printf("Error: Memory allocation failed for rows.\n");
        exit(1);
    }

    // Allocate memory for data (continuous block)
    (*matrix)[0] = (int *)malloc(rows * cols * sizeof(int));
    if ((*matrix)[0] == NULL) 
    {
        printf("Error: Memory allocation failed for matrix data.\n");
        exit(1);
    }

    // Set pointers in array of pointers to point to rows in continuous block
    for (int i = 1; i < rows; ++i) {
        (*matrix)[i] = (*matrix)[i - 1] + cols;
    }
}

void readMatrixFromFile(char *in, int ***matrix, int *rows, int *cols) 
{
    // This function allocates and reads the matrix from the in file.
    FILE *file = fopen(in, "r");

    if (file == NULL) {
        printf("Error: Unable to open file %s for reading.\n", in);
        exit(1);
    }

    if (fscanf(file, "%d %d", rows, cols) != 2) 
    {
        printf("Error: Failed to read matrix dimensions from file.\n");
        exit(1);
    }

    allocateContinuousMatrix(matrix, *rows, *cols);

    for (int i = 0; i < *rows; ++i) {
        for (int j = 0; j < *cols; ++j) {
            if (fscanf(file, "%d", &(*matrix)[i][j]) != 1) 
            {
                printf("Error: Failed to read element at position (%d, %d).\n", i, j);
                exit(1);
            }
        }
    }

    fclose(file);
}

int writeMatrixToFile(char* out, int **matrix, int rows, int cols)
{
    // writes to the output file the number of lines and columns, and the matrix
    FILE *file = fopen(out, "w");
    if (file == NULL) {
        printf("Error: Unable to open file %s for writing.\n", out);
        exit(1);
    }

    if (fprintf(file, "%d ", rows) < 0)
    {
        printf("Ops! Something went wrong1!\n");
    }

    if (fprintf(file, "%d\n", cols) < 0)
    {
        printf("Ops! Something went wrong2!\n");
    }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (fprintf(file, "%d ", matrix[i][j]) < 2)
            {
                printf("Ops! Something went wrong!\n");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
    return 0;
}

void divide_matrix(int rows, int cols, int nProcesses, int** sendcounts, int** displs, int** local_rows, int** starts)
{
    // This function calculates the number of rows for each process
    int chunk_size = rows / nProcesses;
    int remaining_size = rows % nProcesses;

    for (int i = 0; i < nProcesses; i++) 
    {
        (*local_rows)[i] = chunk_size;
        if (i < remaining_size)
        {
            (*local_rows)[i]++;
        }

        (*sendcounts)[i] = (*local_rows)[i] * cols;

        if (i == 0)
        {
            (*starts)[i] = 0;
            (*displs)[i] = 0;
        }
        else
        {
            (*starts)[i] = (*starts)[i - 1] + (*local_rows)[i - 1];
            (*displs)[i] = (*displs)[i - 1] + (*sendcounts)[i - 1];
        }
    }
}

int valid(int i, int j, int rows, int cols) 
{
    // Checks if the index is in range
    return i >= 0 && i < rows && j >= 0 && j < cols;
}

int main(int argc, char **argv)
{
	int rank;
	int nProcesses;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

    char *in, *out;
    int steps;
    getArgs(argc, argv, &in, &out, &steps);

    int **matrix;
    int rows, cols;
    MPI_Status status;
    MPI_Request request;
	int* sendcounts = (int*)malloc(nProcesses * sizeof(int));
	int* displs = (int*)malloc(nProcesses * sizeof(int));

    if (rank == 0)
    {
        // Process 0 reads the matrix and calculates the displacements and counts for every process
        readMatrixFromFile(in, &matrix, &rows, &cols);

        int *local_rows = (int *)calloc(nProcesses, sizeof(int));
        int *starts = (int *)calloc(nProcesses, sizeof(int));

        divide_matrix(rows,cols,nProcesses,&sendcounts,&displs,&local_rows,&starts);

        for (int i = 0; i < nProcesses; i++)
        {
            MPI_Isend(&(local_rows[i]), 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            MPI_Isend(&(starts[i]), 1, MPI_INT, i, 1, MPI_COMM_WORLD, &request);
        }
    }

    int local_row;
    int start;

    // recives the individual data
    MPI_Recv(&local_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&start, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

    // The needed info is broadcasted to all processes
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(sendcounts, nProcesses, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, nProcesses, MPI_INT, 0, MPI_COMM_WORLD);

    // Every process has to allocate the matrix before receiving it from process 0
    if (rank != 0)
    {
        allocateContinuousMatrix(&matrix, rows, cols);
    }

    MPI_Bcast(*matrix, rows * cols, MPI_INT, 0, MPI_COMM_WORLD);

    // Every processes creates a new matrix in which the results of each step will be calculated
    int** new_matrix;
    allocateContinuousMatrix(&new_matrix, rows, cols);

    int k = 0;
    while(k < steps)
    {
        //printf("Process %d \n",rank);
        for(int i = start; i < start + local_row; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                int count = 0;
                int dir[8][2] = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}}; // Sus, jos, dreapta, stanga, diagonale

                for (int l = 0; l < 8; ++l) 
                {
                    int new_i = i + dir[l][0];
                    int new_j = j + dir[l][1];

                    if (valid(new_i, new_j, rows, cols) && matrix[new_i][new_j] == 1) 
                    {
                        count++;
                    }
                }

                if (((count == 3) || (count == 6)) && (matrix[i][j] == 0))
                {
                    new_matrix[i][j] = 1;
                }
                else if (((count == 2) || (count == 3)) && (matrix[i][j] == 1))
                {
                    new_matrix[i][j] = 1;
                }
                else
                {
                    new_matrix[i][j] = 0;
                }
            }
        }

        for(int i = start; i < start + local_row; i++)
            for(int j = 0; j < cols; j++)
                matrix[i][j] = new_matrix[i][j];

        
        if (nProcesses > 1)
        {
            // if there are more than one process then each process will update the marginal lines of its section with the lines calculated by its neighbours
            if (rank > 0)
            {
                MPI_Isend(new_matrix[start], cols, MPI_INT, rank - 1, 9, MPI_COMM_WORLD, &request);
            }
            if (rank < nProcesses - 1)
            {
                printf("%d   %d\n",rank,start + local_row -1);
                MPI_Isend(new_matrix[start + local_row - 1], cols, MPI_INT, rank + 1, 10, MPI_COMM_WORLD, &request);
            }

            if (rank > 0)
            {
                MPI_Recv(matrix[start - 1], cols, MPI_INT, rank - 1, 10, MPI_COMM_WORLD, &status);
            }
            if (rank < nProcesses - 1)
            {
                MPI_Recv(matrix[start + local_row], cols, MPI_INT, rank + 1, 9, MPI_COMM_WORLD, &status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        k++;
        
    }

    //printf("%d    %d\n",rank,local_row);
    
    MPI_Gatherv(*matrix + start*cols, sendcounts[rank], MPI_INT, *matrix, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank==0)
    {
        writeMatrixToFile(out,matrix,rows,cols);
    }
	MPI_Finalize();
	return 0;
}