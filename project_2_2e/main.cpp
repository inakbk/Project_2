#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <time.h>
#include <jacobi.h>
#include <jacobisolver.h>
#include <writetofile.h>

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
    if(argc < 4)
    {
        cout << "Not enough command line arguments given. "
                "Give 3, in the following order: n_step, maximum number of iterations and p_max on command line." << endl;
        cout << "Eks: >> ./main 10 10000 5" << endl;
        exit(1);
    }
    else{
        //variables from command line:
        const int n_step = atof(argv[1]);
        int maxNumberOfIterations = atof(argv[2]);
        const double p_max = atof(argv[3]); //writing p instead of rho

//-------------------------------------------------------------
        //constructing matrix B:
        const double p_min = 0;
        const double h = (p_max - p_min)/n_step;
        vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h
        double w_r = 0;
        vec V = w_r*w_r*p%p + 1.0/p; //new potential

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
        //clocking the operations:
        clock_t start_arma, finish_arma;
        start_arma = clock();

        //solving equations with armadillo lib:
        vec eigval_arma = zeros<vec>(n_step-1);
        mat eigvec_arma = zeros<mat>(n_step-1,n_step-1);
        eig_sym(eigval_arma, eigvec_arma, B);

        finish_arma = clock();
        double time_arma = ( (finish_arma - start_arma)/((double)CLOCKS_PER_SEC ) );

        writetofile fileArma(eigval_arma, eigvec_arma.col(0), p_max, n_step, time_arma, "arma");

//-------------------------------------------------------------
        vec eigval_jacobi = zeros<vec>(n_step-1);
        mat eigvec_jacobi = zeros<mat>(n_step-1,n_step-1);
        int numberOfIterations = 0;
        int converge_test = 1; //initializing to true, false if jacobi did not converge
        double time_jacobi = 0;

        //solving equations with jacobi rotation:
        jacobisolver solveTestMatrix(B, n_step, maxNumberOfIterations);
        solveTestMatrix.solve_w_jacobi_rotation(eigval_jacobi, eigvec_jacobi, numberOfIterations, converge_test, time_jacobi);

        cout << "Total number of iterations with Jacobi rotation: " << numberOfIterations << endl;

        cout << "Time Jacobi: " << time_jacobi << endl;
        writetofile fileJacobi(eigval_jacobi, eigvec_jacobi.col(0), p_max, n_step, time_jacobi, "jacobi", numberOfIterations, converge_test);

//-------------------------------------------------------------
    }

    return 0;
}
