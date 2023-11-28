// Assessment 1 - Investigating the QU Algorithm for Eigenvalue Computation in Square Matrices, with Application to Coupled Harmonic Oscillators
// This program consists of two distinct functionalities:

// Program 1: Eigenvalue Computation
// This program uses the QU iterative algorithm to find eigenvalues of an nxn square matrix provided by the user. The computed solutions are displayed on the screen.

// Program 2: Eigenvalue Analysis for Coupled Harmonic Oscillators
// This program applies the QU algorithm to solve the eigenvalue problem of a pair of coupled harmonic oscillators. 
// It derives a 2x2 matrix from the system's equations of motion and computes both iterated and true eigenvalues. 
// Additionally, the oscillation modes of the system are calculated from the eigenvalues. The results are displayed on the screen, and an output file is generated.
//  A Python script can plot this file, enabling the exploration of the mass-frequency relationship of the system.


// Author: Tom Batchelor
// Date: 03/11/2023

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

void program1(int matrix_size, double** input_matrix, double** output_matrix);
void program2(int matrix_size, double** input_matrix, double** output_matrix);
void AppendDataToFile(FILE *file, double m, double angular_freq_1, double angular_freq_2);
void get_input_matrix(int matrix_size, double** input_matrix);
void print_matrix(int matrix_size, double** matrix);
void calc_F_values(int matrix_size, double** input_matrix, double** F_matrix);
void calc_Q_matrix(int matrix_size, double** F_matrix, double** Q_matrix);
void calc_U_matrix(int matrix_size, double** input_matrix, double** F_matrix, double** Q_matrix, double** U_matrix);
void matrix_multiplier(int matrix_size, double** matrix_1, double** matrix_2, double** output_matrix);
void perform_QU_decomposition(int matrix_size, double** input_matrix, double** output_matrix);
void free_matrix(int matrix_size, double **matrix);
void by_hand_eigen(double k, double m, double* true_eigenvalues);
void matrix_copy(int matrix_size, double **source, double **destination);
double** allocate_matrix(int matrix_size);

int main() {
    int user_selection;
    int matrix_size;
    bool menuRunning = true;

    while (menuRunning) {
        printf("\nPROGRAM MENU\n");
        printf("---------------\n");
        printf("Program 1: Find eigenvalues of an nxn square matrix\n");
        printf("Program 2: Select spring constant (k) and mass (m) to find solutions to the equations of motion of a pair of coupled harmonic oscillators\n");
        printf("Enter 1 or 2 to select which program or 0 to exit: ");
        scanf("%d", &user_selection);

        // initialise input and output matrices as double pointers
        double** input_matrix = NULL;
        double** output_matrix = NULL;

        switch (user_selection) {
            case 1:
                printf("\nPROGRAM 1 SELECTED\n---------------------\n");
                printf("Enter the dimension (n) for an nxn square matrix: ");
                scanf("%i", &matrix_size);

                if (matrix_size > 0){
                    // allocate memory for input and output matrices of size nxn as user specified
                    input_matrix = allocate_matrix(matrix_size);
                    output_matrix = allocate_matrix(matrix_size);

                    // run QU decomposition to find eigenvalues of input matrix
                    program1(matrix_size, input_matrix, output_matrix);

                    // free the memory used by input and output matrices
                    free_matrix(matrix_size, input_matrix);
                    free_matrix(matrix_size, output_matrix);
                    break;
                }
                else{
                    printf("\nPlease enter a valid dimension. Returning to menu.\n");
                    break;
                }
                
            case 2:
                printf("\nPROGRAM 2 SELECTED\n---------------------\n");
                matrix_size = 2;

                // allocate memory for input and output matrices of size 2x2 as this is size of matrix in 2-body harmonic oscillator eigenvalue problem
                input_matrix = allocate_matrix(matrix_size);
                output_matrix = allocate_matrix(matrix_size);

                // run QU decomposition to calculate angular frequency of pair of harmonic oscillators 
                program2(matrix_size, input_matrix, output_matrix);

                // free the memory used by input and output matrices
                free_matrix(matrix_size, input_matrix);
                free_matrix(matrix_size, output_matrix);
                break;
            case 0:
                menuRunning = false;
                printf("Exiting the program. Goodbye!\n");
                break;
            default:
                printf("Invalid selection. Please enter 1, 2, or 0 to exit.\n");
                break;
        }
    }
    return 0;
}

// run QU decomposition to find eigenvalues of an input matrix
void program1(int matrix_size, double **input_matrix, double **output_matrix) {

    // retrieve matrix entries from user
    get_input_matrix(matrix_size, input_matrix);
    printf("Input Matrix:\n");
    print_matrix(matrix_size, input_matrix);

    // run QU decomposition algorithm
    perform_QU_decomposition(matrix_size, input_matrix, output_matrix);

    // print output matrix and eigenvalues
    printf("Final Output Matrix:\n");
    print_matrix(matrix_size, output_matrix);
    printf("---------------------\n");
    for (int i = 0; i < matrix_size; i++) {
        printf("Iterated Eigenvalue %i = %lf \n", i + 1, output_matrix[i][i]);
    }
}

// run QU decomposition to calculate angular frequency of a pair of harmonic oscillators 
void program2(int matrix_size, double** input_matrix, double** output_matrix){

    // retrieve physical values from user
    // k: Spring constant in N/m
    // lower_mass: Lower bound mass in kg
    // upper_mass: Upper bound mass in kg
    // steps: Number of steps in between the masses
    // true_eigenvalues: Array to store the true eigenvalues
    double k,lower_mass, upper_mass, steps, true_eigenvalues[2];
    printf("Enter a value for the spring constant k in N/m: ");
    scanf("%lf", &k);
    printf("Enter an integer value for the lower bound mass m in kg: ");
    scanf("%lf", &lower_mass);
    printf("Enter an integer value for the upper bound mass m in kg: ");
    scanf("%lf", &upper_mass);
    printf("Enter the number of steps in between the masses: ");
    scanf("%lf", &steps);
    
    // Check if input values are valid for the calculation
    if (k > 0 && lower_mass > 0 && upper_mass > 0 && steps > 1 && upper_mass > lower_mass){

        // initialise arrays for outputed angular frequency data 
        double* angular_freq_1 = (double*)malloc(((upper_mass - lower_mass)/steps) * sizeof(double));
        double* angular_freq_2 = (double*)malloc(((upper_mass - lower_mass)/steps) * sizeof(double));

        // open file for writing mass and angular frequency to
        const char *file_path = "data.txt";
        FILE *file = fopen(file_path, "w");
        if (file == NULL) {
            printf("Error opening the file.\n");
            return;
        }
        
        // loop through masses in number of steps specified by user
        for (double m = lower_mass; m <= upper_mass; m += (upper_mass - lower_mass)/steps) {
            int counter = 0;

            // calculate true eigenvalues of each matrix using result from characterstic equation
            by_hand_eigen(k, m, true_eigenvalues);

            // build 'input matrix' for QU decomposition by applying spring constant (k) and mass values (m) to the physical eigenvalue problem
            for (int i = 0; i < matrix_size; i++) {
                for (int j = 0; j < matrix_size; j++) {
                    if (i == j) {
                        input_matrix[j][i] = (-2) * (k/ m);
                    } else {
                        input_matrix[j][i] = (k / m);
                    }
                }
            }

            // run QU decomposition algorithm
            perform_QU_decomposition(matrix_size, input_matrix, output_matrix);
            printf("---------------------\n");
            printf("mass = %fkg; spring constant = %fN/m\n", m, k);
            printf("---------------------\n");

            // print eigenvalues:
            for (int i = 0; i < matrix_size; i++) {
                printf("True Eigenvalue %i: %lf\n", i + 1, true_eigenvalues[i]);
                printf("Iterated Eigenvalue %i = %lf \n", i + 1, output_matrix[i][i]);

                // calculate angular frequency modes from eigenvalues
                if (i == 0){
                    angular_freq_1[counter] = sqrt(-output_matrix[i][i]);
                }
                else{
                    angular_freq_2[counter] = sqrt(-output_matrix[i][i]);
                }
            }
            printf("angular frequency 1 = %lf : angular frequency 2 = %lf \n", angular_freq_1[counter], angular_freq_2[counter]);

            // append mass and frequency modes to file
            AppendDataToFile(file, m, angular_freq_1[counter], angular_freq_2[counter]);
            counter ++;
        }
        fclose(file);

        // call python plotting script that plots mass against both frequency modes
        system("python3 plotting.py");

        // free memory used by angular frequency arrays
        free(angular_freq_1);
        free(angular_freq_2);
    }
    else{
        printf("invalid input");
    }

}

// function to get the input matrix values of an nxn square matrix from user
void get_input_matrix(int matrix_size, double** input_matrix) {
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            printf("Enter element at row %d, column %d: ", i + 1, j + 1);
            scanf("%lf", &input_matrix[i][j]);
        }
    }
}

// function to print a nxn square matrix
void print_matrix(int matrix_size, double** matrix) {
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            printf("%lf\t", matrix[i][j]);
        }
        printf("\n");
    }
}

// function to append mass and frequency data to a file
void AppendDataToFile(FILE *file, double m, double angular_freq_1, double angular_freq_2) {
    // check file opened successfully
    if (file == NULL) {
        printf("Error opening the file.\n");
        return;
    }
    // write mass and angular frequency to file
    fprintf(file, "%lf %lf %lf\n", m, angular_freq_1, angular_freq_2);
}

// function to calculate the F matrix values from the input matrix in the first stage of QU decomposition
void calc_F_values(int matrix_size, double** input_matrix, double** F_matrix) {

    // initialize the temporary column array with zeros
    double* temporary_column = (double*)calloc(matrix_size, sizeof(double));

    // iterate through each column of the input matrix
    for (int input_col = 0; input_col < matrix_size; input_col++) {

        // iterate through each column of the F matrix
        for (int F_col = 0; F_col < matrix_size; F_col++) {

            // perform calculations up to F matrix column = input matrix column - 1
            if (input_col > F_col) {
                double dot_product = 0.0;
                double mod_squared = 0.0;

                // calculate dot product and modulus squared of the F matrix column
                for (int i = 0; i < matrix_size; i++) {
                    dot_product += input_matrix[i][input_col] * F_matrix[i][F_col];
                    mod_squared += pow(F_matrix[i][F_col],2);
                }

                // create the temporary column with the prefactor * relative F matrix column
                for (int i = 0; i < matrix_size; i++) {
                    temporary_column[i] += (dot_product / mod_squared) * F_matrix[i][F_col];
                }
            }
        }
        // compute the final F matrix column by subtracting the temporary column from the input matrix column
        for (int i = 0; i < matrix_size; i++) {
            F_matrix[i][input_col] = input_matrix[i][input_col] - temporary_column[i];
        }
    }
    free(temporary_column);
}

// function to calculate the Q matrix using the F matrix
void calc_Q_matrix(int matrix_size, double** F_matrix, double** Q_matrix) {

    for (int i = 0; i < matrix_size; i++) {
        double mod_squared = 0.0;

        // calculate the modulus squared for the current column of F
        for (int k = 0; k < matrix_size; k++) {
            mod_squared += pow(F_matrix[k][i], 2);
        }

        // calculate Q values by dividing each F value by the modulus squared of its column
        for (int j = 0; j < matrix_size; j++) {
            Q_matrix[j][i] = (F_matrix[j][i] / sqrt(mod_squared));
        }
    }
}

// function to construct the upper triangular U matrix using the F, Q and input matrices
void calc_U_matrix(int matrix_size, double** input_matrix, double** F_matrix, double** Q_matrix, double** U_matrix) {

    // iterate through each column of the input matrix
    for (int input_col = 0; input_col < matrix_size; input_col++) {

        // iterate through each column of the F matrix
        for (int F_col = 0; F_col < matrix_size; F_col++) {

            // diagonal values = modulus of relative column F column
            if (input_col == F_col) {
                double mod_squared = 0.0;
                for (int i = 0; i < matrix_size; i++) {
                    mod_squared += pow(F_matrix[i][input_col], 2);
                }
                U_matrix[F_col][input_col] = sqrt(mod_squared);
            }

            // values below the diagonal (F_col > input_col) are set to 0
            else if (F_col > input_col) {
                U_matrix[F_col][input_col] = 0.0;
            }

            // values above the diagonal computed by taking the dot product of the relative column of the input matrix with the relative row of the Q matrix
            else {
                double dot_product = 0.0;
                for (int i = 0; i < matrix_size; i++) {
                    dot_product += input_matrix[i][input_col] * Q_matrix[i][F_col];
                }
                U_matrix[F_col][input_col] = dot_product;
            }
        }
    }
}

// function to multiply two nxn matrices and store the result in output_matrix.
void matrix_multiplier(int matrix_size, double** matrix_1, double** matrix_2, double** output_matrix) {
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            for (int k = 0; k < matrix_size; k++) {
                output_matrix[i][j] += matrix_2[i][k] * matrix_1[k][j];
            }
        }
    }
}

// function to copy a matrix of size nxn from source to destination
void matrix_copy(int matrix_size, double **source, double **destination) {
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            destination[i][j] = source[i][j];
        }
    }
}

// function to perform QU decomposition until tolerence on diagonals is reached
void perform_QU_decomposition(int matrix_size, double** input_matrix, double** output_matrix) {

    // intitialise variables for iteration tolerance
    double diagonal_difference = 0.0;
    double tolerance = 0.000001;

    // allocate memory for matrices used in stages of QU decomposition
    double** F_matrix = allocate_matrix(matrix_size);
    double** Q_matrix = allocate_matrix(matrix_size);
    double** U_matrix = allocate_matrix(matrix_size);
    double** temp_output_matrix = allocate_matrix(matrix_size);

    // perform decomposition steps
    calc_F_values(matrix_size, input_matrix, F_matrix);
    calc_Q_matrix(matrix_size, F_matrix, Q_matrix);
    calc_U_matrix(matrix_size, input_matrix, F_matrix, Q_matrix, U_matrix);
    matrix_multiplier(matrix_size, Q_matrix, U_matrix, temp_output_matrix);

    // copy temporary output matrix to output matrix
    matrix_copy(matrix_size, temp_output_matrix, output_matrix);

    // loop through the diagonal elements of output matrix and input matrix and calculate the absolute difference for each diagonal and sum them
    for (int i = 0; i < matrix_size; i++) {
        diagonal_difference += sqrt(pow((output_matrix[i][i] - input_matrix[i][i]),2));
    }

    // check if the diagonal difference is less than the tolerance
    if (diagonal_difference <= tolerance) {
        // free the matrices used in this iteration
        free_matrix(matrix_size, F_matrix);
        free_matrix(matrix_size, Q_matrix);
        free_matrix(matrix_size, U_matrix);
        free_matrix(matrix_size, temp_output_matrix);
        return;
    }
    
    // update the input_matrix for the next iteration
    matrix_copy(matrix_size, output_matrix, input_matrix);

    // recursively call the function update output_matrix
    perform_QU_decomposition(matrix_size, input_matrix, output_matrix);

    // free the matrices used in this iteration
    free_matrix(matrix_size, F_matrix);
    free_matrix(matrix_size, Q_matrix);
    free_matrix(matrix_size, U_matrix);
}

// function to calculate the eigenvalues by hand, derived from characteristic equation solutions to the eigenvalue problem of a pair of harmonic oscillators
void by_hand_eigen(double k, double m, double* true_eigenvalues){
    double value = k/m;
    true_eigenvalues[0] = -3*value;
    true_eigenvalues[1] = -1*value;
}

// function to allocate memory for an nxn square matrix
double** allocate_matrix(int matrix_size) {
    double** matrix = (double**)malloc(matrix_size * sizeof(double*));
    for (int i = 0; i < matrix_size; i++) {
        matrix[i] = (double*)malloc(matrix_size * sizeof(double));
    }
    return matrix;
}

// function to free memory allocated for a 2D matrix.
void free_matrix(int matrix_size, double **matrix) {

    // check that the matrix hasnt been freed already to avoid double freeing
    if (matrix) {
        for (int i = 0; i < matrix_size; i++) {
            if (matrix[i]) {
                free(matrix[i]);
            }
        }
        free(matrix);
    }
}

