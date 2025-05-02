#include <iostream>
#include <mpi.h>
#include <fstream>
#include <string>

const int Nx = 100;
const int Ny = 100;
const int num_steps = 100;
const double dt = 0.01;
const double dx = 1.0;
const double dy = 1.0;
const double nu = 0.01;
const double U_top = 1.0;

void set_boundary(double* grid, int start_i, int local_nx, int rank);
void exchange_halo(double* grid, int local_nx, int rank, int size);
void save_velocity_to_csv(double* u, int start_i, int local_nx, int rank, int t);

void set_boundary(double* grid, int start_i, int local_nx, int rank) {
    for (int i = 0; i < local_nx + 2; i++) {
        int global_i = start_i + i - 1;
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;
            if (global_i < 0 || global_i >= Nx || j == 0 || j == Ny - 1) {
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

void save_velocity_to_csv(double* u, int start_i, int local_nx, int rank, int t) {
    std::string filename = "velocity_step_" + std::to_string(t) + "_rank_" + std::to_string(rank) + ".csv";
    std::ofstream outfile(filename);
    outfile << "x,y,u\n";
    for (int i = 1; i <= local_nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;
            outfile << start_i + i - 1 << "," << j << "," << u[idx] << "\n";
        }
    }
    outfile.close();
}

int main() {
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_nx = Nx / size;
    int start_i = rank * local_nx;

    double* u = new double[(local_nx + 2) * Ny]();
    set_boundary(u, start_i, local_nx, rank);

    // Start timer
    double start_time = MPI_Wtime();

    for (int t = 0; t < num_steps; t++) {
        exchange_halo(u, local_nx, rank, size);
        save_velocity_to_csv(u, start_i, local_nx, rank, t);
    }

    // End timer
    double end_time = MPI_Wtime();

    if (rank == 0) {
        std::cout << "Total execution time: " << end_time - start_time << " seconds" << std::endl;
        std::cout << "Number of MPI processes used: " << size << std::endl;
        std::cout << "Total Processing Units Used: " << size << std::endl;
    }

    delete[] u;
    MPI_Finalize();
    return 0;
}


// #include <iostream>
// #include <mpi.h>
// #include <fstream>
// #include <string>

// const int Nx = 100;
// const int Ny = 100;
// const int num_steps = 100;
// const double dt = 0.01;
// const double dx = 1.0;
// const double dy = 1.0;
// const double nu = 0.01;
// const double U_top = 1.0;

// void set_boundary(double* grid, int start_i, int local_nx, int rank);
// void exchange_halo(double* grid, int local_nx, int rank, int size);
// void save_velocity_to_csv(double* u, int start_i, int local_nx, int rank, int t);

// void set_boundary(double* grid, int start_i, int local_nx, int rank) {
//     for (int i = 0; i < local_nx + 2; i++) {
//         int global_i = start_i + i - 1;
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             if (global_i < 0 || global_i >= Nx || j == 0 || j == Ny - 1) {
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

// void save_velocity_to_csv(double* u, int start_i, int local_nx, int rank, int t) {
//     std::string filename = "velocity_step_" + std::to_string(t) + "_rank_" + std::to_string(rank) + ".csv";
//     std::ofstream outfile(filename);
//     outfile << "x,y,u\n";
//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             outfile << start_i + i - 1 << "," << j << "," << u[idx] << "\n";
//         }
//     }
//     outfile.close();
// }

// int main() {
//     MPI_Init(NULL, NULL);
//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     int local_nx = Nx / size;
//     int start_i = rank * local_nx;

//     double* u = new double[(local_nx + 2) * Ny]();
//     set_boundary(u, start_i, local_nx, rank);

//     double start_time = MPI_Wtime();

//     for (int t = 0; t < num_steps; t++) {
//         exchange_halo(u, local_nx, rank, size);
//         save_velocity_to_csv(u, start_i, local_nx, rank, t);
//     }

//     double end_time = MPI_Wtime();
    
//     if (rank == 0) {
//         std::cout << "Total execution time: " << end_time - start_time << " seconds" << std::endl;
//         std::cout << "Number of MPI processes used: " << size << std::endl;
//     }

//     int rank, size;
//         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//         MPI_Comm_size(MPI_COMM_WORLD, &size);
//         if (rank == 0) {
//             std::cout << "Total Processing Units Used: " << size << std::endl;
//     }

//     delete[] u;
//     MPI_Finalize();
//     return 0;
// }



// #include <iostream>
// #include <mpi.h>
// #include <fstream>
// #include <string>

// const int Nx = 100;
// const int Ny = 100;
// const int num_steps = 100;
// const double dt = 0.01;
// const double dx = 1.0;
// const double dy = 1.0;
// const double nu = 0.01;
// const double U_top = 1.0;

// void set_boundary(double* grid, int start_i, int local_nx, int rank);
// void exchange_halo(double* grid, int local_nx, int rank, int size);
// void save_velocity_to_csv(double* u, int start_i, int local_nx, int rank, int t);

// void set_boundary(double* grid, int start_i, int local_nx, int rank) {
//     for (int i = 0; i < local_nx + 2; i++) {
//         int global_i = start_i + i - 1;
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             if (global_i < 0 || global_i >= Nx || j == 0 || j == Ny - 1) {
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

// void save_velocity_to_csv(double* u, int start_i, int local_nx, int rank, int t) {
//     std::string filename = "velocity_step_" + std::to_string(t) + "_rank_" + std::to_string(rank) + ".csv";
//     std::ofstream outfile(filename);
//     outfile << "x,y,u\n";
//     for (int i = 1; i <= local_nx; i++) {
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             outfile << start_i + i - 1 << "," << j << "," << u[idx] << "\n";
//         }
//     }
//     outfile.close();
// }

// int main() {
//     MPI_Init(NULL, NULL);
//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     int local_nx = Nx / size;
//     int start_i = rank * local_nx;

//     double* u = new double[(local_nx + 2) * Ny]();
//     set_boundary(u, start_i, local_nx, rank);

//     double start_time = MPI_Wtime();

//     for (int t = 0; t < num_steps; t++) {
//         exchange_halo(u, local_nx, rank, size);
//         save_velocity_to_csv(u, start_i, local_nx, rank, t);
//     }

//     double end_time = MPI_Wtime();
//     if (rank == 0) {
//         std::cout << "Total time: " << end_time - start_time << " seconds" << std::endl;
//     }

//     delete[] u;
//     MPI_Finalize();
//     return 0;
// }
