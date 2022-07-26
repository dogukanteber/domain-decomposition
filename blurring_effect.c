#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>

#define NUMBER_OF_DIMS 2
#define GLOBAL_GRID_SIZE 8
#define LOCAL_GRID_SIZE 4
#define P_LOCAL_GRID_SIZE LOCAL_GRID_SIZE+2 // PADDED_LOCAL_GRID_SIZE
#define PADDING 1

/*
    A structure that holds neighbors of a process
*/
typedef struct neighbors {
    int north;
    int east;
    int west;
    int south;
    int north_east;
    int north_west;
    int south_east;
    int south_west;
} neighbors;

/*
    To understand what these arrays are, go to line 148.
*/
int init_x_pos[2] = { 0.4 * GLOBAL_GRID_SIZE, 0.5 * GLOBAL_GRID_SIZE };
int init_y_pos[2] = { 0.4 * GLOBAL_GRID_SIZE, 0.5 * GLOBAL_GRID_SIZE };

/*
    Initializing subroutine for neighbors structure. Note that
    we are assuming the length of the coords is always 2
*/
void get_neighbors(int* coords, neighbors* n, MPI_Comm comm) {
    int west_coords[2] = { coords[0], coords[1] - 1 };
    MPI_Cart_rank(comm, west_coords, &(n->west));

    int east_coords[2] = { coords[0], coords[1] + 1 };
    MPI_Cart_rank(comm, east_coords, &(n->east));

    int north_coords[2] = { coords[0] + 1, coords[1] };
    MPI_Cart_rank(comm, north_coords, &(n->north));

    int south_coords[2] = { coords[0] - 1, coords[1] };
    MPI_Cart_rank(comm, south_coords, &(n->south));

    int north_east_coords[2] = { coords[0] + 1, coords[1] + 1};
    MPI_Cart_rank(comm, north_east_coords, &(n->north_east));

    int north_west_coords[2] = { coords[0] + 1, coords[1] - 1};
    MPI_Cart_rank(comm, north_west_coords, &(n->north_west));

    int south_east_coords[2] = { coords[0] - 1, coords[1] + 1};
    MPI_Cart_rank(comm, south_east_coords, &(n->south_east));

    int south_west_coords[2] = { coords[0] - 1, coords[1] - 1};
    MPI_Cart_rank(comm, south_west_coords, &(n->south_west));
}

// start -- debugging functions
void print_neighbors(neighbors* n, int rank) {
    printf("Mine, process %d's, west neighbor's rank is %d\n", rank, n->west);
    printf("Mine, process %d's, east neighbor's rank is %d\n", rank, n->east);
    printf("Mine, process %d's, north neighbor's rank is %d\n", rank, n->north);
    printf("Mine, process %d's, south neighbor's rank is %d\n", rank, n->south);
    printf("Mine, process %d's, north_west neighbor's rank is %d\n", rank, n->north_west);
    printf("Mine, process %d's, north_east neighbor's rank is %d\n", rank, n->north_east);
    printf("Mine, process %d's, south_west neighbor's rank is %d\n", rank, n->south_west);
    printf("Mine, process %d's, south_east neighbor's rank is %d\n", rank, n->south_east);
}

void print_array(int n, int m, int (*arr)[m]) {
    for ( int i=0; i<n; ++i ) {
        for ( int j=0; j<m; ++j ) {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_array_double(int n, int m, double (*arr)[m]) {
    for ( int i=0; i<n; ++i ) {
        for ( int j=0; j<m; ++j ) {
            printf("%lf ", arr[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
// end -- debugging functions

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // start -- setup for cartesian topology
    int dims[NUMBER_OF_DIMS] = { 0, 0 };
    MPI_Dims_create(np, NUMBER_OF_DIMS, dims);
    
    int periods[NUMBER_OF_DIMS] = { true, true };
    int reorder = true;

    MPI_Comm new_comm;
    MPI_Cart_create(MPI_COMM_WORLD, NUMBER_OF_DIMS, dims, periods, reorder, &new_comm);

    int rank;
    MPI_Comm_rank(new_comm, &rank);

    int my_coords[NUMBER_OF_DIMS];
    MPI_Cart_coords(new_comm, rank, NUMBER_OF_DIMS, my_coords);
    // end -- setup for cartesian topology

    // start -- setup for our custom datatype

    /*
        In this block we are defining our custom datatype to ease writing to the
        file or reading from the file. Also, this datatype is going to make our lives
        much easier when doing the calculation for blurring effect.

        Below you can see the representation of our custom datatype:
        g g g g
        g x x g
        g x x g
        g g g g

        g represents a guard cell, d represents a data
    */
    MPI_Datatype local_grid;
    int start[2] = {PADDING, PADDING};
    // we are adding 2 because we wrap guard cells around our grid
    int arrsize[2] = {P_LOCAL_GRID_SIZE * PADDING, P_LOCAL_GRID_SIZE * PADDING};
    int gridsize[2] = {LOCAL_GRID_SIZE, LOCAL_GRID_SIZE};
    MPI_Type_create_subarray(2, arrsize, gridsize,
                                start, MPI_ORDER_C, MPI_DOUBLE, &local_grid);
    MPI_Type_commit(&local_grid);
    // end -- setup for our custom datatype

    // start -- initializing local grid in each process
    
    /*
        Beware that in both of the loops we are not starting from 0 and not ending at
        length of the row/col. This is because we wrapped our grid with guard cells.
        Hence, we do not touch those cells.

        The if-else condition and i_global, j_global variables are used for initializing
        the matrix. We want our initial matrix to contain 1's in the middle as the shape
        of a square. For instance, if the global matrix is 4x4, we want it to look like this:

            0 0 0 0
            0 1 1 0
            0 1 1 0
            0 0 0 0
        
        Therefore, these variables and the if statement is used for this purpose.
    */

    double start_time, end_time;

    start_time = MPI_Wtime();
    double local_arr[P_LOCAL_GRID_SIZE * PADDING][P_LOCAL_GRID_SIZE * PADDING] = { 0 };
    for ( int i = PADDING; i <= LOCAL_GRID_SIZE; ++i ) {
        int i_global = (i - 1) + LOCAL_GRID_SIZE * my_coords[0];
        for ( int j = PADDING; j <= LOCAL_GRID_SIZE; ++j ) {
            int j_global = (j - 1) + LOCAL_GRID_SIZE * my_coords[1];
            if ( (i_global == init_x_pos[0] || i_global == init_x_pos[1]) && (j_global == init_y_pos[0] || j_global == init_y_pos[1]) ) {
                local_arr[i][j] = 1;
            }
            else {
                local_arr[i][j] = 0;
            }
        }
    }
    // end -- initializing local grid in each process

    // start -- setup our another custom datatype for file operations

    /*
        Since we are not writing contiguous block of data (we are using 2d arrays),
        we cannot use the existing file operations without creating our own datatype.
        This datatype lets us write 2d arrays into the file and read from file to 2d
        arrays.
    */
    MPI_Datatype view_local_arr;
    int startV[2] = { my_coords[0]*LOCAL_GRID_SIZE, my_coords[1]*LOCAL_GRID_SIZE };
    int arrsizeV[2] = { dims[0]*LOCAL_GRID_SIZE, dims[1]*LOCAL_GRID_SIZE };
    int gridsizeV[2] = { LOCAL_GRID_SIZE, LOCAL_GRID_SIZE };

    MPI_Type_create_subarray(2, arrsizeV, gridsizeV,
                                startV, MPI_ORDER_C, MPI_DOUBLE, &view_local_arr);
    MPI_Type_commit(&view_local_arr);
    // end -- setup our another custom datatype for file operations

    // start -- putting initial values to our file
    MPI_File fh;
    MPI_File_open(new_comm, "result.dat", MPI_MODE_CREATE | MPI_MODE_WRONLY,
                    MPI_INFO_NULL, &fh);

    MPI_File_set_view(fh, 0, MPI_DOUBLE, view_local_arr, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, &local_arr[0][0], 1, local_grid, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    // end -- putting initial values to our file

    // start -- setup datatypes for communicating with other processes in one call

    /*
        row_type will let us send/receive a block of row in one call.
    */
    MPI_Datatype row_type;
    MPI_Type_vector(LOCAL_GRID_SIZE, 1, 1, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    /*
        col_type will let us send/receive a block of col in one call.
    */
    MPI_Datatype col_type;
    // we are adding 2 because of guard cells
    MPI_Type_vector(LOCAL_GRID_SIZE, 1, LOCAL_GRID_SIZE + 2, MPI_DOUBLE, &col_type);
    MPI_Type_commit(&col_type);
    // end -- setup datatypes for communicating with other processes in one call

    int steps = 3;
    for ( int step=0; step<steps; ++step) {
        // start -- printing the results

        /*
            Below is very inefficient in terms of performance. But this is just for
            illustration purposes. It is used to print the result to the screen
            since we write and read binary data when working with files. We can
            extract this block of code and the program will still give us the
            same result as it was given earlier.
        */
        // double res[GLOBAL_GRID_SIZE][GLOBAL_GRID_SIZE];
        // MPI_File_open(new_comm, "result.dat", MPI_MODE_RDONLY,
        //             MPI_INFO_NULL, &fh);
        //     MPI_File_read(fh, &res[0][0], GLOBAL_GRID_SIZE*GLOBAL_GRID_SIZE, MPI_DOUBLE, MPI_STATUS_IGNORE);
        //     MPI_File_close(&fh);

        // if ( rank == 0) {
        //     print_array_double(GLOBAL_GRID_SIZE, GLOBAL_GRID_SIZE, res);
        // }
        // end -- printing the results

        // start -- getting data from the file
        MPI_File_open(new_comm, "result.dat", MPI_MODE_RDONLY,
                    MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, MPI_DOUBLE, view_local_arr, "native", MPI_INFO_NULL);

        double arr[P_LOCAL_GRID_SIZE][P_LOCAL_GRID_SIZE] = { 0 };
        MPI_File_read_all(fh, &arr[0][0], 1, local_grid, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        // end -- getting data from the file

        // start -- communication between processees

        /*
            For each process, we are making a call to all
            neighbor processes to give us what we want.
        */
        neighbors n;
        get_neighbors(my_coords, &n, new_comm);

        MPI_Request requests[16];
        MPI_Status statuses[16];

        MPI_Isend(&arr[1][1], 1, col_type, n.west, 0, new_comm, &requests[0]);
        MPI_Irecv(&arr[1][LOCAL_GRID_SIZE+1], 1, col_type, n.east, 0, new_comm, &requests[1]);

        MPI_Isend(&arr[1][LOCAL_GRID_SIZE], 1, col_type, n.east, 0, new_comm, &requests[2]);
        MPI_Irecv(&arr[1][0], 1, col_type, n.west, 0, new_comm, &requests[3]);

        MPI_Isend(&arr[LOCAL_GRID_SIZE][1], 1, row_type, n.south, 0, new_comm, &requests[4]);
        MPI_Irecv(&arr[0][1], 1, row_type, n.north, 0, new_comm, &requests[5]);

        MPI_Isend(&arr[1][1], 1, row_type, n.north, 0, new_comm, &requests[6]);
        MPI_Irecv(&arr[LOCAL_GRID_SIZE+1][1], 1, row_type, n.south, 0, new_comm, &requests[7]);

        MPI_Isend(&arr[LOCAL_GRID_SIZE][LOCAL_GRID_SIZE], 1, MPI_DOUBLE, n.south_east, 0, new_comm, &requests[8]);
        MPI_Irecv(&arr[0][0], 1, MPI_DOUBLE, n.north_west, 0, new_comm, &requests[9]);

        MPI_Isend(&arr[LOCAL_GRID_SIZE][1], 1, MPI_DOUBLE, n.south_west, 0, new_comm, &requests[10]);
        MPI_Irecv(&arr[0][LOCAL_GRID_SIZE+1], 1, MPI_DOUBLE, n.north_east, 0, new_comm, &requests[11]);

        MPI_Isend(&arr[1][LOCAL_GRID_SIZE], 1, MPI_DOUBLE, n.north_east, 0, new_comm, &requests[12]);
        MPI_Irecv(&arr[LOCAL_GRID_SIZE+1][0], 1, MPI_DOUBLE, n.south_west, 0, new_comm, &requests[13]);

        MPI_Isend(&arr[1][1], 1, MPI_DOUBLE, n.north_west, 0, new_comm, &requests[14]);
        MPI_Irecv(&arr[LOCAL_GRID_SIZE+1][LOCAL_GRID_SIZE+1], 1, MPI_DOUBLE, n.south_east, 0, new_comm, &requests[15]);

        MPI_Waitall(16, requests, statuses);
        // end -- communication between processees

        // start -- blurring effect calculation

        /*
            This is where the actual calculation happens. Note that we are ignoring
            guard cells when looping. In the inner loop, we are calculating a value and 
            are putting it in new_array we have created earlier.
        */
        double new_arr[P_LOCAL_GRID_SIZE][P_LOCAL_GRID_SIZE] = { 0 };
        for ( int i=1; i<=LOCAL_GRID_SIZE; ++i ) {
            for ( int j=1; j<=LOCAL_GRID_SIZE; ++j ) {
                new_arr[i][j] = ( arr[i-1][j-1] + arr[i-1][j] + arr[i-1][j+1] + \
                                  arr[i][j-1]   + arr[i][j]   + arr[i][j+1] + \
                                  arr[i+1][j-1] + arr[i+1][j] + arr[i+1][j+1] ) / 9;
            }
        }
        // end -- blurring effect calculation

        // start -- writing our results to a file
        MPI_File_open(new_comm, "result.dat", MPI_MODE_CREATE | MPI_MODE_WRONLY,
                 MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, MPI_DOUBLE, view_local_arr, "native", MPI_INFO_NULL);
        MPI_File_write_all(fh, &new_arr[0][0], 1, local_grid, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        // end -- writing our results to a file
    }
    end_time = MPI_Wtime();

    if ( rank == 0 ) {
        printf("Elapsed time: %lf - Number of processes: %d\n", end_time - start_time, np);
    }

    MPI_Type_free(&local_grid);
    MPI_Type_free(&view_local_arr);
    MPI_Type_free(&row_type);
    MPI_Type_free(&col_type);


    MPI_Finalize();

    return EXIT_SUCCESS;
}
