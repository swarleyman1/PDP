#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>

int main(int argc, char** argv) {

	if (argc != 3) {
        	printf("Expected: Input filename and Output filename \n");
        	return 0;
    	}

    char* input = argv[1];
    char* output = argv[2];
    float **A, **B, **C, **A_local, **B_local;

	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    ////////Read matrix inputs 
    if(rank == 0){

        FILE* file;
        file = fopen(input, "r");
            int n;
            fscanf(file, "%d", &n);

            A = malloc(n * sizeof(float *));
            B = malloc(n * sizeof(float *));
            C = malloc(n * sizeof(float *));

            for (int i = 0; i < n; i++) {
                    A[i] = malloc(n * sizeof(float));
                    B[i] = malloc(n * sizeof(float));
            }

        for (int i=0; i<(n); i++){
            for (int j=0; j<(n); j++){
                        fscanf(file, "%f", &A[i][j]);		}
            }
        for (int i=0; i<(n); i++){
                    for (int j=0; j<(n); j++){
                            fscanf(file, "%f", &B[i][j]);
                    }
            }

        for (int i=0; i<(n); i++){
                for (int j=0; j<(n); j++){
                        printf("%f ", A[i][j]);
                }
            printf("\n");
            }
        printf("\n");
        for (int i=0; i<(n); i++){
                for (int j=0; j<(n); j++){
                        printf("%f ", B[i][j]);
                }
            printf("\n");
            }

    MPI_Comm new_communicator;
    //MPI_Cart_create(MPI_COMM_WORLD, 2, (int[]){n, n}, (int[]){true, true}, true, &new_communicator);

    A_local = malloc((n/size) * sizeof(float *));
    B_local = malloc((n/size) * sizeof(float *));
    //C_local = malloc((n/size) * sizeof(float *));
    
    for (int i = 0; i < n; i++) {
        A_local[i] = malloc(n * sizeof(float));
        B_local[i] = malloc(n * sizeof(float));
        }

    // Scatter matrix A to all processes
    MPI_Scatter(A, n / size * n, MPI_INT, A_local[0], n / size * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Scatter matrix B to all processes
    MPI_Scatter(B, n / size * n, MPI_INT, B_local[0], n / size * n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (int i=0; i<(n/size); i++){
                for (int j=0; j<(n/size); j++){
                        printf("%f ", A_local[i][j]);
                }
            printf("\n");
            }
        printf("\n");
        for (int i=0; i<(n/size); i++){
                for (int j=0; j<(n/size); j++){
                        printf("%f ", B_local[i][j]);
                }
            printf("\n");
            }

    }

    


    if(rank == 0){
        free(A); 
        free(B);
    }
   
    //free(C);
    MPI_Finalize();
}


/*
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                C[i][j] = A[i][j]*B[i][j];
                    }
            }

        for (int i=0; i<n; i++){
                for (int j=0; j<n; j++){
                        printf("%f ", C[i][j]);
                }
            printf("\n");
            }
            */