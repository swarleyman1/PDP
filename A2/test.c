#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv){
  int rank, size;
	MPI_Status status;
	MPI_Request send_request, recv_request;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
 
   int n;
 
  if (rank == 0){
    if (argc != 3)
  	{
  		printf("Usage: stencil input_file output_file\n");
  		return 1;
  	}
    char *input_name = argv[1];
	  char *output_name = argv[2];
   
   // Read input file
		if (0 > (num_values = read_input(input_name, &input)))
		{
			return 1;
		}

    n = num_values[0];
 }
 
 MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
 
 if (rank == 0){
   float A[n*n], B[n*n], C[n*n];
   for (int i = 0; i < n*n; ++i){
     A[i] = num_values[i+1];
     B[i] = num_values[i+1+n*n];
   }
 
 // Define the datatype for a column
    MPI_Datatype column_type;
    MPI_Type_vector(
        n, // count
        1, // blocklength
        n, // stride
        MPI_FLOAT, // oldtype
        &column_type
    );
    MPI_Type_commit(&column_type);
    
    //Send columns of b
    float B_own[n*(n/size)];
    MPI_Scatter(B, n/size, column_type, B_own, n*(n/size), MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    //Send out rows for A
    int rows = 1;
    float A_own[n*rows];
    if (rank == 0){
    for (int row = 0; row < n; row += rows){
          for (int col = 0; col < n; ++col){
            A_own[col] = A[row*n+col];
          }
    }
    }
    //MPI_Bcast(A_own, n*rows, 
    //MPI_Scatter(A, n*(n/size), MPI_FLOAT, A_own, n*(n/size), MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    //Do matrix multiplication
    //for
    
 }
 
 free(B);
  return 0;
}

int read_input(const char *file_name, float **values)
{
	FILE *file;
	if (NULL == (file = fopen(file_name, "r")))
	{
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values))
	{
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(float))))
	{
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i = 0; i < num_values; i++)
	{
		if (EOF == fscanf(file, "%lf", &((*values)[i])))
		{
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file))
	{
		perror("Warning: couldn't close input file");
	}
	return num_values;
}