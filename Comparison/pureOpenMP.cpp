#include <iostream>
#include <omp.h>
#include <fstream>
#include <string>
#include <cstring>

const int Nx = 100;
const int Ny = 100;
const int num_steps = 100;
const double dt = 0.01;
const double dx = 1.0;
const double dy = 1.0;
const double nu = 0.01; // Kinematic viscosity
const double U_top = 1.0; // Top wall velocity

void set_boundary(double* grid, const char* var_type);
void save_velocity_to_csv(double* u, double* v, int t);

void set_boundary(double* grid, const char* var_type) {
    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;
            if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {
                grid[idx] = (strcmp(var_type, "u") == 0 && j == Ny - 1) ? U_top : 0.0;
            }
        }
    }
}

void save_velocity_to_csv(double* u, double* v, int t) {
    std::string filename = "velocity_step_" + std::to_string(t) + ".csv";
    std::ofstream outfile(filename);
    outfile << "x,y,u,v\n";
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            int idx = i * Ny + j;
            outfile << i << "," << j << "," << u[idx] << "," << v[idx] << "\n";
        }
    }
    outfile.close();
}

int main() {
    double* u = new double[Nx * Ny]();
    double* v = new double[Nx * Ny]();
    double* u_star = new double[Nx * Ny]();
    double* v_star = new double[Nx * Ny]();

    set_boundary(u, "u");
    set_boundary(v, "v");

    // Get the number of OpenMP threads
    int num_threads = 1;
    #pragma omp parallel
    {
        #pragma omp master
        {
            num_threads = omp_get_max_threads();
        }
    }

    double start_time = omp_get_wtime();

    for (int t = 0; t < num_steps; t++) {
        #pragma omp parallel for
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                int idx = i * Ny + j;
                int idx_left = (i - 1) * Ny + j;
                int idx_right = (i + 1) * Ny + j;
                int idx_up = i * Ny + (j + 1);
                int idx_down = i * Ny + (j - 1);

                u_star[idx] = u[idx] - dt * (
                    0.5 * (u[idx] + u[idx_left]) * (u[idx_right] - u[idx_left]) / (2.0 * dx) +
                    0.5 * (v[idx] + v[idx_down]) * (u[idx_up] - u[idx_down]) / (2.0 * dy)
                ) + nu * dt * (
                    (u[idx_right] - 2.0 * u[idx] + u[idx_left]) / (dx * dx) +
                    (u[idx_up] - 2.0 * u[idx] + u[idx_down]) / (dy * dy)
                );
            }
        }

        std::swap(u, u_star);
        save_velocity_to_csv(u, v, t);
    }

    double end_time = omp_get_wtime();
    
    std::cout << "Total execution time: " << end_time - start_time << " seconds" << std::endl;
    std::cout << "Number of OpenMP threads used: " << num_threads << std::endl;
    std::cout << "Total Processing Units Used: " << num_threads << std::endl;

    delete[] u;
    delete[] v;
    delete[] u_star;
    delete[] v_star;

    return 0;
}


// #include <iostream>
// #include <omp.h>
// #include <fstream>
// #include <string>
// #include <cstring>

// const int Nx = 100;
// const int Ny = 100;
// const int num_steps = 100;
// const double dt = 0.01;
// const double dx = 1.0;
// const double dy = 1.0;
// const double nu = 0.01; // kinematic viscosity
// const double U_top = 1.0; // Top wall velocity

// void set_boundary(double* grid, const char* type);
// void save_velocity_to_csv(double* u, double* v, int t);

// void set_boundary(double* grid, const char* var_type) {
//     #pragma omp parallel for
//     for (int i = 0; i < Nx; i++) {
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {
//                 if (strcmp(var_type, "u") == 0) {
//                     grid[idx] = (j == Ny - 1) ? U_top : 0.0;
//                 } else {
//                     grid[idx] = 0.0;
//                 }
//             }
//         }
//     }
// }

// void save_velocity_to_csv(double* u, double* v, int t) {
//     std::string filename = "velocity_step_" + std::to_string(t) + ".csv";
//     std::ofstream outfile(filename);
//     outfile << "x,y,u,v\n";
//     for (int i = 0; i < Nx; i++) {
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             outfile << i << "," << j << "," << u[idx] << "," << v[idx] << "\n";
//         }
//     }
//     outfile.close();
// }

// int main() {
//     double* u = new double[Nx * Ny]();
//     double* v = new double[Nx * Ny]();
//     double* u_star = new double[Nx * Ny]();
//     double* v_star = new double[Nx * Ny]();

//     set_boundary(u, "u");
//     set_boundary(v, "v");

//     int num_threads = 1;
//     #pragma omp parallel
//     {
//         #pragma omp master
//         {
//             num_threads = omp_get_max_threads();
//         }
//     }

//     double start_time = omp_get_wtime();

//     for (int t = 0; t < num_steps; t++) {
//         #pragma omp parallel for
//         for (int i = 1; i < Nx - 1; i++) {
//             for (int j = 1; j < Ny - 1; j++) {
//                 int idx = i * Ny + j;
//                 int idx_left = (i - 1) * Ny + j;
//                 int idx_right = (i + 1) * Ny + j;
//                 int idx_up = i * Ny + (j + 1);
//                 int idx_down = i * Ny + (j - 1);

//                 u_star[idx] = u[idx] - dt * (
//                     0.5 * (u[idx] + u[idx_left]) * (u[idx_right] - u[idx_left]) / (2.0 * dx) +
//                     0.5 * (v[idx] + v[idx_down]) * (u[idx_up] - u[idx_down]) / (2.0 * dy)
//                 ) + nu * dt * (
//                     (u[idx_right] - 2.0 * u[idx] + u[idx_left]) / (dx * dx) +
//                     (u[idx_up] - 2.0 * u[idx] + u[idx_down]) / (dy * dy)
//                 );
//             }
//         }

//         std::swap(u, u_star);
//         save_velocity_to_csv(u, v, t);
//     }

//     double end_time = omp_get_wtime();
//     std::cout << "Total execution time: " << end_time - start_time << " seconds" << std::endl;
//     std::cout << "Number of OpenMP threads used: " << num_threads << std::endl;

//     int num_threads = omp_get_max_threads();
//     std::cout << "Total Processing Units Used: " << num_threads << std::endl;
    
//     delete[] u;
//     delete[] v;
//     delete[] u_star;
//     delete[] v_star;

//     return 0;
// }



// #include <iostream>
// #include <omp.h>
// #include <fstream>
// #include <string>
// #include <cstring>

// const int Nx = 100;
// const int Ny = 100;
// const int num_steps = 100;
// const double dt = 0.01;
// const double dx = 1.0;
// const double dy = 1.0;
// const double nu = 0.01; // kinematic viscosity
// const double U_top = 1.0; // Top wall velocity

// void set_boundary(double* grid, const char* type);
// void save_velocity_to_csv(double* u, double* v, int t);

// void set_boundary(double* grid, const char* var_type) {
//     #pragma omp parallel for
//     for (int i = 0; i < Nx; i++) {
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {
//                 if (strcmp(var_type, "u") == 0) {
//                     grid[idx] = (j == Ny - 1) ? U_top : 0.0;
//                 } else {
//                     grid[idx] = 0.0;
//                 }
//             }
//         }
//     }
// }

// void save_velocity_to_csv(double* u, double* v, int t) {
//     std::string filename = "velocity_step_" + std::to_string(t) + ".csv";
//     std::ofstream outfile(filename);
//     outfile << "x,y,u,v\n";
//     for (int i = 0; i < Nx; i++) {
//         for (int j = 0; j < Ny; j++) {
//             int idx = i * Ny + j;
//             outfile << i << "," << j << "," << u[idx] << "," << v[idx] << "\n";
//         }
//     }
//     outfile.close();
// }

// int main() {
//     double* u = new double[Nx * Ny]();
//     double* v = new double[Nx * Ny]();
//     double* u_star = new double[Nx * Ny]();
//     double* v_star = new double[Nx * Ny]();

//     set_boundary(u, "u");
//     set_boundary(v, "v");

//     double start_time = omp_get_wtime();

//     for (int t = 0; t < num_steps; t++) {
//         #pragma omp parallel for
//         for (int i = 1; i < Nx - 1; i++) {
//             for (int j = 1; j < Ny - 1; j++) {
//                 int idx = i * Ny + j;
//                 int idx_left = (i - 1) * Ny + j;
//                 int idx_right = (i + 1) * Ny + j;
//                 int idx_up = i * Ny + (j + 1);
//                 int idx_down = i * Ny + (j - 1);

//                 u_star[idx] = u[idx] - dt * (
//                     0.5 * (u[idx] + u[idx_left]) * (u[idx_right] - u[idx_left]) / (2.0 * dx) +
//                     0.5 * (v[idx] + v[idx_down]) * (u[idx_up] - u[idx_down]) / (2.0 * dy)
//                 ) + nu * dt * (
//                     (u[idx_right] - 2.0 * u[idx] + u[idx_left]) / (dx * dx) +
//                     (u[idx_up] - 2.0 * u[idx] + u[idx_down]) / (dy * dy)
//                 );
//             }
//         }

//         std::swap(u, u_star);
//         save_velocity_to_csv(u, v, t);
//     }

//     double end_time = omp_get_wtime();
//     std::cout << "Total time: " << end_time - start_time << " seconds" << std::endl;

//     delete[] u;
//     delete[] v;
//     delete[] u_star;
//     delete[] v_star;

//     return 0;
// }
