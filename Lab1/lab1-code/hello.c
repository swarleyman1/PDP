/**********************************************************************
 * A simple "hello world" program for MPI/C
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

  int size, rank;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */

  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processes  */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get the rank of this process */


  if (rank == 0) {                      /* Only the first process prints this  */
    printf("There are %d processes\n", size);
  }

  printf("Hello from %d!\n", rank);     /* Print a message              */

  MPI_Finalize();                       /* Shut down and clean up MPI   */

  return 0;
}
