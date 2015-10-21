#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <time.h>
#include <jacobi.h>

using namespace std;
using namespace arma;


void WriteToFile(const vec& eigenvalues, int p_max, const int n_step, const int number_of_iterations, const double time, string FileName, const int converge_test)
{
    ofstream myfile;
        string filename = "EigvalSolver_" + FileName + "_pMax" + to_string(p_max) + "_nStep" + to_string(n_step) + ".txt";
        myfile.open (filename);
        myfile << "Solution of the eigenvalueproblem for the " << FileName << " algorithm." << endl;
        myfile << "Dimention of matrix, one less than n_step: " << n_step << endl;
        myfile << "Value of p_max: " << p_max << endl;
        myfile << "Execution time: " << time << endl;
        if(number_of_iterations != -1){
            myfile << "Number of iterations for jacobi algoritm: " << number_of_iterations << endl;
        }
        else
        {
            myfile << endl;
        }
        //Writing eigenvalues to file if the jacobi method did converge.
        if(converge_test == true)
        {
            myfile << "---------------------" << endl;
            myfile << "Eigenvalues (sorted)" << "  "<< "" << "     " << "Eigenvectors (soon)" << endl;
            myfile << "---------------------" << endl;
            int number_of_eigenvalues_printed = 10;
            for (int i=0; i < number_of_eigenvalues_printed; i++)
            {
                if(i == size(eigenvalues,0))
                {
                    cout << "Length of eigenval vec is shorter than number_of_eigenvalues_printed for " << FileName << " solver, exiting loop." << endl;
                    break;
                }
                myfile << eigenvalues[i] << "    " << "--" << endl;
            }
            myfile << "(Maximum writing " << number_of_eigenvalues_printed << " eigenvalues to file.)" << endl;
            myfile.close();

            cout << "Datafile done for n_step=" << n_step << " with the " << FileName <<" solver." << endl;
            cout << endl;
        }
        //Writing error message til file if jacobi method did not converge.
        if(converge_test == false)
        {
            myfile << endl; myfile << endl;
            myfile << "Jacobi method did not converge after " << number_of_iterations << " iterations!" << endl;
            myfile << endl; myfile << endl;

            cout << "Datafile done for n_step=" << n_step << " with the " << FileName <<" solver which did not converge." << endl;
            cout << endl;
        }
}




//jacobi myJacobiSolver;

//myJacobiSolver.


int main(int argc, char *argv[])
{
    if(argc == 1)
    {
        cout << "No arguments. Give n_step, maximum number of iterations and p_max on command line "
                "Eks: ./main n_step no_of_iterations p_max or: ./main 10 10000 5" << endl;
        exit(1);
    }
    else{
        //variables from command line:
        const int n_step = atof(argv[1]);
        int maxNumberOfIterations = atof(argv[2]);
        const double p_max = atof(argv[3]); //writing p instead of rho

    //-------------------------------------------------------------
        const double p_min = 0;
        const double h = (p_max - p_min)/n_step;
        vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h
        vec V = p%p;

        // Constructing matrix B:
        double e = -1/(h*h); // all elements of the e vec is the same
        vec d = 2/(h*h) + V;
        mat B = zeros<mat>(n_step-1,n_step-1);
        for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
        {
            B(i,j) = e;
            B(i,i) = d[i+1];
            B(j,i) = e;
        }
        B(n_step-2,n_step-2) = d[n_step-1];

    //-------------------------------------------------------------
        //clocking the operations (only solve, not making file):
        clock_t start_arma, finish_arma; //declaring start and finish time
        start_arma = clock();

        //solving equations with armadillo lib:
        vec eigval_arma = eig_sym(B);

        //stopping timer:
        finish_arma = clock();
        double time_arma = ( (finish_arma - start_arma)/((double)CLOCKS_PER_SEC ) );
        //cout << "Armadillo lib. eigenvalue solver: Time for n_step="
        //     << n_step << ":  " << time_arma << " seconds" << endl;

        int converge_test = 1; //initializing to true (false if it does not converge) for jacobi method
        WriteToFile(eigval_arma, p_max, n_step, -1, time_arma, "arma", converge_test);
    //-------------------------------------------------------------
        //clocking the operations (only solve, not making file):
        clock_t start_jacobi, finish_jacobi; //declaring start and finish time
        start_jacobi = clock();

        //solving equations with jacobi rotation:
        int numberOfIterations = solve_eq_jacobi_rotation(B, n_step, maxNumberOfIterations, converge_test);

        //stopping timer:
        finish_jacobi = clock();
        double time_jacobi = ( (finish_jacobi - start_jacobi)/((double)CLOCKS_PER_SEC ) );
    //    cout << "Jacobi rotation eigenvalue solver: Time for n_step="
    //         << n_step << ":  " << time_jacobi << " seconds" << endl;
        cout << "Total number of iterations with Jacobi rotation: " << numberOfIterations << endl;

        //retriving eigenvalues, interested in the three first:
        vec eigval_jacobi_rot = B.diag();
        eigval_jacobi_rot = sort(eigval_jacobi_rot);

        WriteToFile(eigval_jacobi_rot, p_max, n_step, numberOfIterations, time_jacobi, "jacobi", converge_test);
    //-------------------------------------------------------------

        /// fix write to file. write all data? or just the first ones?
        ///
        /// how big n_step need to get three lowest eigvals to 4 leading digits?
        /// dependency og p_max? save/write to file?
        /// how many transformations (before 0) as function of n_step
        /// also write time to file as function of n_step for both solvers
        /// if did not converge write that to file!!


    }

    return 0;
}

