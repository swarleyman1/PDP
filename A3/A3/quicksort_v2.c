#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <string.h>

/**************************************************
Parallel quicksort using MPI for A3 in PDP course

Authors:
    - Jacob Malmenstedt
    - Viveka Olsson
    - Emy Engström
**************************************************/



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
int read_input(const char *file_name, int **values)
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
    if (NULL == (*values = malloc(num_values * sizeof(int))))
    {
        perror("Couldn't allocate memory for input");
        return -1;
    }
    for (int i = 0; i < num_values; i++)
    {
        if (EOF == fscanf(file, "%d", &((*values)[i])))
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
int write_output(char *file_name, const int *output, int num_values)
{
    FILE *file;
    if (NULL == (file = fopen(file_name, "w")))
    {
        perror("Couldn't open output file");
        return -1;
    }
    for (int i = 0; i < num_values; i++)
    {
        if (0 > fprintf(file, "%d ", output[i]))
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
    return (*(int *)a - *(int *)b);
}

// Swap function
void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Insertion sort
void insertion_sort(int *arr, int n)
{ // loop over all elements
    for (int i = 1; i < n; i++)
    {
        int key = arr[i]; // element to be inserted
        int j = i - 1;    // element before key
        // loop over elements before key moving them to the right
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key; // insert key
    }
}

void insertion_sort_double(double *arr, double n)
{ // loop over all elements
    for (int i = 1; i < n; i++)
    {
        double key = arr[i]; // element to be inserted
        int j = i - 1;       // element before key
        // loop over elements before key moving them to the right
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key; // insert key
    }
}

// Merge two sorted arrays
void merge(int *arr1, int *arr2, int n1, int n2, int *arr3)
{
    int i = 0, j = 0, k = 0;
    // Traverse both array
    while (i < n1 && j < n2)
    {
        // Check if current element of first
        // array is smaller than current element
        // of second array. 
        if (arr1[i] < arr2[j]) // if yes, store first array element and increment first array index
        {
            arr3[k] = arr1[i];
            k++;
            i++;
        }
        else
        {
            arr3[k] = arr2[j]; // else store second array element and increment second array index
            k++;
            j++;
        }
    }

    // Store remaining elements of first array
    while (i < n1)
    {
        arr3[k] = arr1[i];
        k++;
        i++;
    }

    // Store remaining elements of second array
    while (j < n2)
    {
        arr3[k] = arr2[j];
        k++;
        j++;
    }
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
                pivot = (data[length / 2] + data[length / 2 - 1]) / 2;
            }
            else
            {
                pivot = data[length / 2];
            }
        }

        MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);

        break;

    case 2:
        // Select the median of all medians in each processor group
        double median1;
        double *medians1 = malloc(size * sizeof(double));

        if (length % 2 == 0)
        {
            median1 = (data[length / 2] + data[length / 2 - 1]) / 2;
        }
        else
        {
            median1 = data[length / 2];
        }

        // MPI_Allgather(&median, 1, MPI_DOUBLE, medians, 1, MPI_DOUBLE, comm);
        MPI_Gather(&median1, 1, MPI_DOUBLE, medians1, 1, MPI_DOUBLE, 0, comm);

        if (rank == 0)
        {
            insertion_sort_double(medians1, size);

            if (size % 2 == 0) // Will always be even since size is a power of 2
            {
                pivot = (medians1[size / 2] + medians1[size / 2 - 1]) / 2;
            }
            else
            {
                pivot = medians1[size / 2];
            }
        }

        MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);
        free(medians1);

        break;

    case 3:
        // Select the mean value of all medians in each processor group
        double median2;
        double *medians2 = malloc(size * sizeof(double));

        if (length % 2 == 0)
        {
            median2 = (data[length / 2] + data[length / 2 - 1]) / 2;
        }
        else
        {
            median2 = data[length / 2];
        }

        MPI_Gather(&median2, 1, MPI_DOUBLE, medians2, 1, MPI_DOUBLE, 0, comm);

        if (rank == 0)
        {
            double temp = 0;
            for (int i = 0; i < size; i++)
            {
                temp += medians2[i];
            }
            pivot = temp / size;
        }

        MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);
        free(medians2);

        break;

    default:
        printf("ERROR: Invalid pivot method\n");
        break;
    }

    return pivot;
}

// recursive parallel quicksort
int QuicksortInner(int *data, int length, MPI_Comm comm, int pivot_method, int depth)
{
    int rank, size;
    double pivot;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Base case
    if (size == 1 || depth == 0)
    {
        return length;
    }

    // Select pivot
    pivot = pivot_selection(data, length, pivot_method, comm);

    // Partition data around pivot
    int i = 0;
    int j = length - 1;
    while (i < j)
    {
        while (data[i] < pivot && i < length)
        {
            i++;
        }
        while (data[j] >= pivot && j > 0)
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
    int send_count, recv_count;
    if (rank < size / 2)
    {
        send_count = length - i;
    }
    else
    {
        send_count = i;
    }

    // find partner rank
    int partner_rank;
    if (rank < size / 2)
    {
        partner_rank = rank + size / 2;
    }
    else
    {
        partner_rank = rank - size / 2;
    }

    // send and receive data
    if (rank < size / 2)
    {
        MPI_Ssend(&send_count, 1, MPI_INT, partner_rank, 0, comm);
        MPI_Recv(&recv_count, 1, MPI_INT, partner_rank, 0, comm, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Recv(&recv_count, 1, MPI_INT, partner_rank, 0, comm, MPI_STATUS_IGNORE);
        MPI_Ssend(&send_count, 1, MPI_INT, partner_rank, 0, comm);
    }

    // Create send, receive and temp buffers
    int *send_buffer = malloc(send_count * sizeof(int));
    int *recv_buffer = malloc(recv_count * sizeof(int));
    int *temp_buffer = malloc(length * sizeof(int));

    // Copy data to send buffer
    if (rank < size / 2)
    {
        memcpy(send_buffer, &data[i], send_count * sizeof(int)); // Send data larger than pivot
    }
    else
    {
        memcpy(send_buffer, &data[0], send_count * sizeof(int)); // Send data smaller than pivot
    }

    if (rank < size / 2)
    {
        MPI_Send(send_buffer, send_count, MPI_INT, partner_rank, depth, comm);
        MPI_Recv(recv_buffer, recv_count, MPI_INT, partner_rank, depth, comm, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Recv(recv_buffer, recv_count, MPI_INT, partner_rank, depth, comm, MPI_STATUS_IGNORE);
        MPI_Send(send_buffer, send_count, MPI_INT, partner_rank, depth, comm);
    }

    // Place remaining data in temp buffer
    if (rank < size / 2)
    {
        memcpy(temp_buffer, &data[0], i * sizeof(int));
    }
    else
    {
        memcpy(temp_buffer, &data[i], (length - i) * sizeof(int));
    }

    // Clear data
    memset(data, 0, length * sizeof(int));

    // Merge data
    // merge(int *arr1, int *arr2, int n1, int n2, int *arr3)
    if (rank < size / 2)
    {
        merge(temp_buffer, recv_buffer, i, recv_count, data);
    }
    else
    {
        merge(temp_buffer, recv_buffer, length - i, recv_count, data);
    }

    // Free buffers
    free(send_buffer);
    free(recv_buffer);
    free(temp_buffer);

    // Update length
    length = length - send_count + recv_count;

    // Merge data in each processor
    // (actually uses insertion sort to sort data instead of merging.
    // Should be fast enough since data is already partially sorted)
    //insertion_sort(data, length);

    //qsort(data, length, sizeof(int), compare);

    // Split communicator into two groups
    int color = (rank < size / 2) ? 0 : 1;
    MPI_Comm new_comm;
    MPI_Comm_split(comm, color, rank, &new_comm);

    // Recursively call quicksort on each group
    length = QuicksortInner(data, length, new_comm, pivot_method, depth - 1);

    // Free communicator
    MPI_Comm_free(&new_comm);

    return length;
}

void Quicksort(int *data, int length, int pivot_method, int max_depth)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Special case where length < size and Scatterv will fail
    // (for all reasonable values of size and length, this should never happen
    // but in case it does insertion sort will be faster than quicksort)
    if (length < size)
    {
        if (rank == 0)
        {
            insertion_sort(data, length);
        }
        return;
    }

    int chunk_size = length / size;
    int remainder = length % size;

    int *local_list = malloc(length * sizeof(int));
    int *send_counts = malloc(size * sizeof(int));
    int *send_displacements = malloc(size * sizeof(int));

    // Calculate send counts and displacements
    for (int i = 0; i < size; i++)
    {
        send_counts[i] = chunk_size;
        if (i < remainder)
        {
            send_counts[i]++;
        }
        send_displacements[i] = i * chunk_size;

        if (i < remainder)
        {
            send_displacements[i] += i;
        }
        else
        {
            send_displacements[i] += remainder;
        }
    }

    // Scatter data
    MPI_Scatterv(data, send_counts, send_displacements, MPI_INT, local_list, length, MPI_INT, 0, MPI_COMM_WORLD);

    // Sort local list
    qsort(local_list, send_counts[rank], sizeof(int), compare);

    // Go into recursion and collect lengths of local lists
    int local_length = QuicksortInner(local_list, send_counts[rank], MPI_COMM_WORLD, pivot_method, max_depth);

    // gather local lengths
    int *local_lengths = malloc(size * sizeof(int));
    MPI_Gather(&local_length, 1, MPI_INT, local_lengths, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Create displacement array and receive buffer
    int *displacements = malloc(size * sizeof(int));
    int *recv_buffer = malloc(length * size * sizeof(int));

    if (rank == 0)
    {
        displacements[0] = 0;
        for (int i = 1; i < size; i++)
        {
            displacements[i] = displacements[i - 1] + local_lengths[i - 1];
        }
    }

    // Gather using gatherv does not work
    // Gatherv syntax: MPI_Gatherv(send_buffer, send_count, send_type, recv_buffer, recv_counts, displacements, recv_type, root, comm)
    //MPI_Gatherv(local_list, length, MPI_INT, data, local_lengths, displacements, MPI_INT, 0, MPI_COMM_WORLD);

    // Gather data into recv_buffer instead of gatherv
    MPI_Gather(local_list, length, MPI_INT, recv_buffer, length, MPI_INT, 0, MPI_COMM_WORLD);

    // Pick data from recv_buffer
    if (rank == 0)
    {
        for (int i = 0; i < size; i++)
        {
            int k = i * length;
            for (int j = displacements[i]; j < displacements[i] + local_lengths[i]; j++)
            {
                data[j] = recv_buffer[k];
                k++;
            }
        }
    }

    free(local_list);
    free(send_counts);
    free(send_displacements);
    free(local_lengths);
    free(displacements);
    free(recv_buffer);

    return;
}

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

    // Stop error messages in UPPMAX
    char* env_var = "OMPI_MCA_btl=^openib";
    int ret = putenv(env_var);
    if (ret != 0) {
        perror("putenv");
        return 1;
    }

    int rank, size;
    MPI_Status status;
    MPI_Request send_request, recv_request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int num_values;
    int *list_to_sort;

    //  Read input file
    if (rank == 0)
    {
        if (0 > (num_values = read_input(input_name, &list_to_sort)))
        {
            printf("ERROR: Could not read input file\n");
            return 2;
        }
    }

    // Broadcast number of values to all processors
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Start timer
    double time = MPI_Wtime();

    // Go into parallel quicksort
    Quicksort(list_to_sort, num_values, pivot_method, 2 * log2(size));

    // Stop timer
    time = MPI_Wtime() - time;

    // Get max time
    double max_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    // Write output and free memory
    if (rank == 0)
    {
        printf("%f\n", max_time);
        write_output(output_name, list_to_sort, num_values);
        free(list_to_sort);
    }    

    MPI_Finalize();

    return 0;
}