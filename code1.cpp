#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <string> 


const int Nx = 100;
const int Ny = 100;
const int num_steps = 100;
const double dt = 0.01;
const double dx = 1.0;
const double dy = 1.0;
const double nu = 0.01; // kinematic viscosity
const double U_top = 1.0; // Top wall velocity

void set_boundary(double* grid, int start_i, int local_nx, int rank, int size, int Ny, const char* type); 
void exchange_halo(double* grid, int local_nx, int Ny, int rank, int size);

void set_boundary(double* grid, int start_i, int local_nx, int rank, int size, int ny, const char* var_type) {
    // Set boundary conditions based on var_type (u, v, or p)
    int local_grid_size_x = local_nx + 2;
    for (int i = 0; i < local_grid_size_x; i++) {
        int global_i = start_i + i - 1; // Map local to global
        for (int j = 0; j < ny; j++) {
            int idx = i * ny + j;
            if (global_i < 0 || global_i >= Nx || j == 0 || j == ny - 1) {
                if (strcmp(var_type, "u") == 0) {
                    if (j == ny - 1) grid[idx] = U_top; // Top wall
                    else grid[idx] = 0.0; // No-slip elsewhere
                } else if (strcmp(var_type, "v") == 0) {
                    grid[idx] = 0.0; // No-slip for v
                } else if (strcmp(var_type, "p") == 0) {
                    // Zero gradient for pressure at boundaries
                    if (j == 0 && i > 0 && i < local_grid_size_x - 1) grid[idx] = grid[idx + ny];
                    else if (j == ny - 1 && i > 0 && i < local_grid_size_x - 1) grid[idx] = grid[idx - ny];
                    // Handle x boundaries via halos
                }
            }
        }
    }
}

void exchange_halo(double* grid, int local_nx, int ny, int rank, int size) {
    MPI_Request request[4];
    int num_requests = 0;

    if (rank > 0) {
        MPI_Irecv(&grid[0], ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request[num_requests++]);
        MPI_Send(&grid[ny], ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    }
    if (rank < size - 1) {
        MPI_Irecv(&grid[(local_nx + 1) * ny], ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request[num_requests++]);
        MPI_Send(&grid[local_nx * ny], ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    }

    if (num_requests > 0) MPI_Waitall(num_requests, request, MPI_STATUSES_IGNORE);
}

// Function to save velocity field to CSV
void save_velocity_to_csv(double* u, double* v, int start_i, int local_nx, int ny, int rank, int t) {
    std::string filename = "velocity_step_" + std::to_string(t) + "_rank_" + std::to_string(rank) + ".csv";
    std::ofstream outfile(filename);
    outfile << "x,y,u,v\n";
    for (int i = 1; i <= local_nx; i++) {
        int global_x = start_i + i - 1;
        for (int j = 0; j < ny; j++) {
            int idx = i * ny + j;
            outfile << global_x << "," << j << "," << u[idx] << "," << v[idx] << "\n";
        }
    }
    outfile.close();
}

int main() {
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int base_size = Nx / size;
    int remainder = Nx % size;
    int start_i = rank * base_size + std::min(rank, remainder);
    int end_i = start_i + base_size + (rank < remainder ? 1 : 0) - 1;
    int local_nx = end_i - start_i + 1;

    int local_grid_size_x = local_nx + 2;
    double* u_current = new double[local_grid_size_x * Ny]();
    double* v_current = new double[local_grid_size_x * Ny]();
    double* p_current = new double[local_grid_size_x * Ny]();
    double* u_next = new double[local_grid_size_x * Ny]();
    double* v_next = new double[local_grid_size_x * Ny]();
    double* p_next = new double[local_grid_size_x * Ny]();
    double* u_star = new double[local_grid_size_x * Ny]();
    double* v_star = new double[local_grid_size_x * Ny]();
    double* div_u_star = new double[local_grid_size_x * Ny]();

    // Initialize grids
    #pragma omp parallel for
    for (int i = 0; i < local_grid_size_x * Ny; i++) {
        u_current[i] = 0.0;
        v_current[i] = 0.0;
        p_current[i] = 0.0;
    }

    set_boundary(u_current, start_i, local_nx, rank, size, Ny, "u");
    set_boundary(v_current, start_i, local_nx, rank, size, Ny, "v");
    set_boundary(p_current, start_i, local_nx, rank, size, Ny, "p");

    double time_start = MPI_Wtime();

    for (int t = 0; t < num_steps; t++) {
        // Step 1: Exchange halo data for u_current and v_current
        exchange_halo(u_current, local_nx, Ny, rank, size);
        exchange_halo(v_current, local_nx, Ny, rank, size);

        // Set boundary conditions
        set_boundary(u_current, start_i, local_nx, rank, size, Ny, "u");
        set_boundary(v_current, start_i, local_nx, rank, size, Ny, "v");

        // Step 2: Compute u_star and v_star
        #pragma omp parallel for
        for (int i = 1; i <= local_nx; i++) {
            for (int j = 0; j < Ny; j++) {
                int idx = i * Ny + j;
                int idx_left = (i-1) * Ny + j;
                int idx_right = (i+1) * Ny + j;
                int idx_up = i * Ny + (j+1);
                int idx_down = i * Ny + (j-1);
                u_star[idx] = u_current[idx] - dt * (
                    0.5 * (u_current[idx] + u_current[idx_left]) * (u_current[idx_right] - u_current[idx_left]) / (2.0 * dx) +
                    0.5 * (v_current[idx] + v_current[idx_down]) * (u_current[idx_up] - u_current[idx_down]) / (2.0 * dy)
                ) + nu * dt * (
                    (u_current[idx_right] - 2.0 * u_current[idx] + u_current[idx_left]) / (dx * dx) +
                    (u_current[idx_up] - 2.0 * u_current[idx] + u_current[idx_down]) / (dy * dy)
                );
                v_star[idx] = v_current[idx] - dt * (
                    0.5 * (u_current[idx] + u_current[idx_left]) * (v_current[idx_right] - v_current[idx_left]) / (2.0 * dx) +
                    0.5 * (v_current[idx] + v_current[idx_down]) * (v_current[idx_up] - v_current[idx_down]) / (2.0 * dy)
                ) + nu * dt * (
                    (v_current[idx_right] - 2.0 * v_current[idx] + v_current[idx_left]) / (dx * dx) +
                    (v_current[idx_up] - 2.0 * v_current[idx] + v_current[idx_down]) / (dy * dy)
                );
            }
        }

        // Step 3: Compute div_u_star
        #pragma omp parallel for
        for (int i = 1; i <= local_nx; i++) {
            for (int j = 0; j < Ny; j++) {
                int idx = i * Ny + j;
                int idx_right = (i+1) * Ny + j;
                int idx_up = i * Ny + (j+1);
                div_u_star[idx] = (u_star[idx_right] - u_star[idx]) / dx + (v_star[idx_up] - v_star[idx]) / dy;
            }
        }

        // Step 4: Solve pressure Poisson equation iteratively (Jacobi)
        for (int iter = 0; iter < 10; iter++) {
            exchange_halo(p_current, local_nx, Ny, rank, size);
            #pragma omp parallel for
            for (int i = 1; i <= local_nx; i++) {
                for (int j = 1; j < Ny - 1; j++) {
                    int idx = i * Ny + j;
                    int idx_left = (i-1) * Ny + j;
                    int idx_right = (i+1) * Ny + j;
                    int idx_down = i * Ny + (j-1);
                    int idx_up = i * Ny + (j+1);
                    p_next[idx] = 0.25 * (
                        p_current[idx_left] + p_current[idx_right] +
                        p_current[idx_down] + p_current[idx_up]
                    ) - (dx * dx / (4.0 * dt)) * div_u_star[idx];
                }
            }
            double* temp = p_current;
            p_current = p_next;
            p_next = temp;
        }

        // Set boundary conditions for p_current
        set_boundary(p_current, start_i, local_nx, rank, size, Ny, "p");

        // Step 5: Correct velocities
        #pragma omp parallel for
        for (int i = 1; i <= local_nx; i++) {
            for (int j = 0; j < Ny; j++) {
                int idx = i * Ny + j;
                int idx_right = (i+1) * Ny + j;
                int idx_up = i * Ny + (j+1);
                u_next[idx] = u_star[idx] - dt / dx * (p_current[idx_right] - p_current[idx]);
                v_next[idx] = v_star[idx] - dt / dy * (p_current[idx_up] - p_current[idx]);
            }
        }

        // Swap current and next grids
        double* temp_u = u_current;
        u_current = u_next;
        u_next = temp_u;
        double* temp_v = v_current;
        v_current = v_next;
        v_next = temp_v;

        // Save velocity data at each time step
        save_velocity_to_csv(u_current, v_current, start_i, local_nx, Ny, rank, t);
    }

    double time_end = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Total time: " << time_end - time_start << " seconds" << std::endl;
    }

    delete[] u_current;
    delete[] v_current;
    delete[] p_current;
    delete[] u_next;
    delete[] v_next;
    delete[] p_next;
    delete[] u_star;
    delete[] v_star;
    delete[] div_u_star;

    MPI_Finalize();
    return 0;
}
