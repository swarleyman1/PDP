#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

// Parallel quicksort using MPI

/*
Algorithm:
1 Divide the data into p equal parts, one per process
2 Sort the data locally for each process
3 Perform global sort
    3.1 Select pivot element within each process set
    3.2 Locally in each process, divide the data into two sets according to the
    pivot (smaller or larger)
    3.3 Split the processes into two groups and exchange data pairwise between
    them so that all processes in one group get data less than the pivot and
    the others get data larger than the pivot.
    3.4 Merge the two sets of numbers in each process into one sorted list
4 Repeat 3.1 - 3.4 recursively for each half until each
group consists of one single process.
*/

// Read the input file and return the number of elements (from A1)
int read_input(const char *file_name, double **values)
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
    if (NULL == (*values = malloc(num_values * sizeof(double))))
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

// Write the output file (from A1)
int write_output(char *file_name, const double *output, int num_values)
{
    FILE *file;
    if (NULL == (file = fopen(file_name, "w")))
    {
        perror("Couldn't open output file");
        return -1;
    }
    for (int i = 0; i < num_values; i++)
    {
        if (0 > fprintf(file, "%.4f ", output[i]))
        {
            perror("Couldn't write to output file");
        }
    }
    if (0 > fprintf(file, "\n"))
    {
        perror("Couldn't write to output file");
    }
    if (0 != fclose(file))
    {
        perror("Warning: couldn't close output file");
    }
    return 0;
}

// Compare function for qsort
int compare(const void *a, const void *b)
{
    return (*(double *)a - *(double *)b);
}

// Main
int main(int argc, char **argv)
{
    if (4 != argc)
    {
        printf("Usage: stencil input_file output_file number_of_applications\n");
        return 1;
    }
    char *input_name = argv[1];
    char *output_name = argv[2];
    int pivot_method = atoi(argv[3]);

    int rank, size;
    MPI_Status status;
    MPI_Request send_request, recv_request;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

    int num_values;
    double *input;
    double *output;
    // Read input file
    if (rank == 0)
    {
        if (0 > (num_values = read_input(input_name, &input)))
        {
            return 2;
        }
    }

    // Broadcast the number of values to all processes
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate memory for output
    if (NULL == (output = malloc(num_values * sizeof(double))))
    {
        perror("Couldn't allocate memory for output");
        return -1;
    }

    // initialize variables
    int chunk_size = num_values / size;
    int remainder = num_values % size; // insertion sort the remainder after the parallel sort??

    // Scatter the input data to all processes
    double *local_list = malloc(chunk_size * sizeof(double));
    MPI_Scatter(input, chunk_size, MPI_DOUBLE, local_list, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*
    if(rank == 0 && remainder != 0){
        realloc(local_list, (chunk_size + remainder) * sizeof(double));
        for(int i = 0; i < remainder; i++){
            local_list[chunk_size + i] = input[num_values - remainder + i];
        }
    }
    */

    // Sort the local list
    qsort(local_list, chunk_size, sizeof(double), compare);

    int depth = 0;
    int num_groups = 2;
    int group_size = size / num_groups;

    // Perform global sort
    while (group_size > 1)
    {
        // Select pivot element within each process set
        double pivot;
        switch (pivot_method)
        {
        case 1:
            // Select the median in one processor (first) in each group of processors
            if (rank % group_size == 0)
            {
                if (chunk_size % 2 == 0)
                {
                    pivot = (local_list[chunk_size / 2] + local_list[chunk_size / 2 - 1]) / 2;
                }
                else
                {
                    pivot = local_list[chunk_size / 2];
                }
            }
            break;

        case 2:
            // Select the median of all medians in each processor group
            int median;
            double *medians = malloc(size * sizeof(double));
            if (chunk_size % 2 == 0)
            {
                median = (local_list[chunk_size / 2] + local_list[chunk_size / 2 - 1]) / 2;
            }
            else
            {
                median = local_list[chunk_size / 2];
            }

            MPI_Allgather(&median, 1, MPI_DOUBLE, medians, 1, MPI_DOUBLE, MPI_COMM_WORLD);
            if (rank % group_size == 0)
            {
                double temp[group_size];
                for (int i = rank; i < rank + group_size; i++)
                {
                    temp[i] = medians[i];
                }
                qsort(temp, group_size, sizeof(int), compare); // Can we use qsort here? or do we need to implement our own?
                if (group_size % 2 == 0)
                {
                    pivot = (temp[group_size / 2] + temp[group_size / 2 - 1]) / 2;
                }
                else
                {
                    pivot = temp[group_size / 2];
                }
            }
            break;

        case 3:
            // Select the mean value of all medians in each processor group
            int median;
            double *medians = malloc(size * sizeof(double));
            if (chunk_size % 2 == 0)
            {
                median = (local_list[chunk_size / 2] + local_list[chunk_size / 2 - 1]) / 2;
            }
            else
            {
                median = local_list[chunk_size / 2];
            }

            MPI_Allgather(&median, 1, MPI_DOUBLE, medians, 1, MPI_DOUBLE, MPI_COMM_WORLD);
            if (rank % group_size == 0)
            {
                double temp = 0;
                for (int i = rank; i < rank + group_size; i++)
                {
                    temp += medians[i];
                }
                pivot = temp / group_size;
            }
            break;

        default:
            printf("Invalid pivot method\n");
            return 1;
        }

        // Broadcast the pivot to all processes in the group
        // MPI_Bcast syntax: MPI_Bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator)
    }



    return 0;
}
