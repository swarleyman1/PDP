#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

// Parallel quicksort using MPI (unfinished)

/*
Run with: mpirun -np 4 ./quicksort input_file output_file pivot_method
(number of processes must be a perfect square)

pivot_method:
1 = median of first process in each group
2 = median of medians
3 = mean of medians


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

// Swap function
void swap(double *a, double *b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

// Select pivot element
double pivot_selection(int *data, int length, int pivot_method, MPI_Comm comm)
{
    double pivot;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    switch (pivot_method)
    {
    case 1:
        // Select the median in one processor (first) in each group of processors
        if (rank == 0)
        {
            if (length % 2 == 0)
            {
                pivot = (local_list[length / 2] + local_list[length / 2 - 1]) / 2;
            }
            else
            {
                pivot = local_list[length / 2];
            }

            MPI_Bcast(&pivot, 1, MPI_DOUBLE, rank, comm);
        }
        break;

    case 2:
        // Select the median of all medians in each processor group
        double median, medians[size];

        if (length % 2 == 0)
        {
            median = (local_list[length / 2] + local_list[length / 2 - 1]) / 2;
        }
        else
        {
            median = local_list[length / 2];
        }

        // MPI_Allgather(&median, 1, MPI_DOUBLE, medians, 1, MPI_DOUBLE, comm);
        MPI_Gather(&median, 1, MPI_DOUBLE, medians, 1, MPI_DOUBLE, 0, comm);

        if (rank == 0)
        {
            qsort(medians, size, sizeof(double), compare);
            if (size % 2 == 0) // Will always be even since size is a power of 2
            {
                pivot = (medians[size / 2] + medians[size / 2 - 1]) / 2;
            }
            else
            {
                pivot = medians[size / 2];
            }

            MPI_Bcast(&pivot, 1, MPI_DOUBLE, rank, comm);
        }
        break;

    case 3:
        // Select the mean value of all medians in each processor group
        double median, medians[size];

        if (length % 2 == 0)
        {
            median = (local_list[length / 2] + local_list[length / 2 - 1]) / 2;
        }
        else
        {
            median = local_list[length / 2];
        }

        // MPI_Allgather(&median, 1, MPI_DOUBLE, medians, 1, MPI_DOUBLE, comm);
        MPI_Gather(&median, 1, MPI_DOUBLE, medians, 1, MPI_DOUBLE, 0, comm);

        if (rank == 0)
        {
            double temp = 0;
            for (int i = 0; i < size; i++)
            {
                temp += medians[i];
            }
            pivot = temp / size;

            MPI_Bcast(&pivot, 1, MPI_DOUBLE, rank, comm);
        }
        break;

    default:
        printf("ERROR: Invalid pivot method\n");
        break;
    }

    return pivot;
}

// recursive parallel quicksort
void QuicksortInner(int *data, int length, MPI_Comm comm, int pivot_method, int depth)
{
    int rank, size;
    double pivot;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Base case
    if (size == 1 || depth == 0)
    {
        // TODO
        return;
    }

    // Select pivot
    pivot = pivot_selection(data, length, pivot_method, comm);

    // Partition data
    int i = 0, j = length - 1;
    while (i < j)
    {
        while (data[i] < pivot)
        {
            i++;
        }
        while (data[j] > pivot)
        {
            j--;
        }
        if (i < j)
        {
            swap(&data[i], &data[j]);
        }
    }

    // Swap data between processors
    // (left half of processors send data larger than pivot
    // to right half of processors and vice versa)

    // ------ TODO ------

    // Split communicator into two groups
    int color = (rank < size / 2) ? 0 : 1;
    MPI_Comm new_comm;
    MPI_Comm_split(comm, color, rank, &new_comm);

    // Recursively call quicksort on each group
    QuicksortInner(data, length, new_comm, pivot_method, depth - 1);

    // Free communicator
    MPI_Comm_free(&new_comm);

    return;
}

void Quicksort(int *data, int length, int pivot_method, int max_depth)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunk_size = length / size;
    double *local_list = malloc(chunk_size * sizeof(double));
    MPI_Scatter(data, length, MPI_INT, local_list, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    qsort(local_list, chunk_size, sizeof(double), compare);

    QuicksortInner(local_list, chunk_size, MPI_COMM_WORLD, pivot_method, max_depth);

    MPI_Gather(local_list, chunk_size, MPI_INT, data, length, MPI_INT, 0, MPI_COMM_WORLD);

    free(local_list);

    return;
}

// Main
int main(int argc, char **argv)
{
    if (4 != argc)
    {
        printf("Usage: quicksort input_file output_file pivot_method\n");
        return 1;
    }
    char *input_name = argv[1];
    char *output_name = argv[2];
    int pivot_method = atoi(argv[3]);

    int rank, size;
    MPI_Status status;
    MPI_Request send_request, recv_request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int num_values;
    double *input;
    // double *output; // Not used
    //  Read input file
    if (rank == 0)
    {
        if (0 > (num_values = read_input(input_name, &input)))
        {
            printf("ERROR: Could not read input file\n");
            return 2;
        }
        // output = malloc(num_values * sizeof(double));
    }

    // Broadcast number of values to all processors
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Start timer
    double start_time = MPI_Wtime();

    // Go into parallel quicksort
    Quicksort(input, num_values, pivot_method, 2 * log2(size));

    // Stop timer
    double end_time = MPI_Wtime();

    // Print time
    if (rank == 0)
    {
        printf("%f\n", end_time - start_time);
    }

    // Write output file and free memory
    if (rank == 0)
    {
        write_output(output_name, input, num_values);

        free(input);
        // free(output);
    }

    MPI_Finalize();

    return 0;
}