#include <iostream>
#include <mpi.h>
#include <omp.h>

const int Nx = 100;
const int Ny = 100;
const int num_steps = 100;
const double dt = 0.01;
const double dx = 1.0;
const double dy = 1.0;
const double nu = 0.01;
const double U_top = 1.0;

void set_boundary(double* grid, int start_i, int local_nx);
void exchange_halo(double* grid, int local_nx, int rank, int size);

void set_boundary(double* grid, int start_i, int local_nx) {
    #pragma omp parallel for
    for (int i = 0; i < local_nx + 2; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;
            if (j == 0 || j == Ny - 1) {
                grid[idx] = (j == Ny - 1) ? U_top : 0.0;
            }
        }
    }
}

void exchange_halo(double* grid, int local_nx, int rank, int size) {
    if (rank > 0) {
        MPI_Sendrecv(&grid[Ny], Ny, MPI_DOUBLE, rank - 1, 0,
                     &grid[0], Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Sendrecv(&grid[local_nx * Ny], Ny, MPI_DOUBLE, rank + 1, 0,
                     &grid[(local_nx + 1) * Ny], Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int num_threads = omp_get_max_threads(); // Fixed thread count retrieval

    if (rank == 0) {
        std::cout << "Number of MPI processes: " << size << std::endl;
        std::cout << "Number of OpenMP threads per process: " << num_threads << std::endl;
        std::cout << "Total Processing Units Used: " << size * num_threads << std::endl;
    }

    int local_nx = Nx / size;
    int start_i = rank * local_nx;

    double* u = new double[(local_nx + 2) * Ny]();
    set_boundary(u, start_i, local_nx);

    double start_time = MPI_Wtime();

    for (int t = 0; t < num_steps; t++) {
        exchange_halo(u, local_nx, rank, size);

        #pragma omp parallel for collapse(2)
        for (int i = 1; i <= local_nx; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                int idx = i * Ny + j;
                int idx_left = (i - 1) * Ny + j;
                int idx_right = (i + 1) * Ny + j;
                int idx_up = i * Ny + (j + 1);
                int idx_down = i * Ny + (j - 1);

                u[idx] = u[idx] - dt * (
                    0.5 * (u[idx] + u[idx_left]) * (u[idx_right] - u[idx_left]) / (2.0 * dx) +
                    0.5 * (u[idx] + u[idx_down]) * (u[idx_up] - u[idx_down]) / (2.0 * dy)
                ) + nu * dt * (
                    (u[idx_right] - 2.0 * u[idx] + u[idx_left]) / (dx * dx) +
                    (u[idx_up] - 2.0 * u[idx] + u[idx_down]) / (dy * dy)
                );
            }
        }
    }

    double end_time = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Total execution time: " << end_time - start_time << " seconds" << std::endl;
    }

    delete[] u;
    MPI_Finalize();
    return 0;
}


// #include <iostream>
// #include <mpi.h>
// #include <omp.h>

// const int Nx = 100;
// const int Ny = 100;
// const int num_steps = 100;
// const double dt = 0.01;
// const double dx = 1.0;
// const double dy = 1.0;
// const double nu = 0.01;
// const double U_top = 1.0;

// void set_boundary(double* grid, int start_i, int local_nx);
// void exchange_halo(double* grid, int local_nx, int rank, int size);

// void set_boundary(double* grid, int start_i, int local_nx) {
//     #pragma omp parallel for
//     for (int i = 0; i < local_nx + 2; i++) {
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             if (j == 0 || j == Ny - 1) {
//                 grid[idx] = (j == Ny - 1) ? U_top : 0.0;
//             }
//         }
//     }
// }

// void exchange_halo(double* grid, int local_nx, int rank, int size) {
//     if (rank > 0) {
//         MPI_Sendrecv(&grid[Ny], Ny, MPI_DOUBLE, rank - 1, 0,
//                      &grid[0], Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }
//     if (rank < size - 1) {
//         MPI_Sendrecv(&grid[local_nx * Ny], Ny, MPI_DOUBLE, rank + 1, 0,
//                      &grid[(local_nx + 1) * Ny], Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }
// }

// int main(int argc, char* argv[]) {
//     MPI_Init(&argc, &argv);
//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     int local_nx = Nx / size;
//     int start_i = rank * local_nx;

//     double* u = new double[(local_nx + 2) * Ny]();
//     set_boundary(u, start_i, local_nx);

//     double start_time = MPI_Wtime();

//     for (int t = 0; t < num_steps; t++) {
//         exchange_halo(u, local_nx, rank, size);
        
//         #pragma omp parallel for collapse(2)
//         for (int i = 1; i <= local_nx; i++) {
//             for (int j = 1; j < Ny - 1; j++) {
//                 int idx = i * Ny + j;
//                 int idx_left = (i - 1) * Ny + j;
//                 int idx_right = (i + 1) * Ny + j;
//                 int idx_up = i * Ny + (j + 1);
//                 int idx_down = i * Ny + (j - 1);

//                 u[idx] = u[idx] - dt * (
//                     0.5 * (u[idx] + u[idx_left]) * (u[idx_right] - u[idx_left]) / (2.0 * dx) +
//                     0.5 * (u[idx] + u[idx_down]) * (u[idx_up] - u[idx_down]) / (2.0 * dy)
//                 ) + nu * dt * (
//                     (u[idx_right] - 2.0 * u[idx] + u[idx_left]) / (dx * dx) +
//                     (u[idx_up] - 2.0 * u[idx] + u[idx_down]) / (dy * dy)
//                 );
//             }
//         }
//     }

//     double end_time = MPI_Wtime();
//     if (rank == 0) {
//         std::cout << "Total time: " << end_time - start_time << " seconds" << std::endl;
//     }

//     delete[] u;
//     MPI_Finalize();
//     return 0;
// }
