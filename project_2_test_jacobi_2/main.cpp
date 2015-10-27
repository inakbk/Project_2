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
        // Constructing test matrix B:
        mat B = randu<mat>(n_step-1,n_step-1);
        B = B.t()*B;  //making shure it is symetric

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

        writetofile fileJacobi(eigval_jacobi, eigvec_jacobi.col(0), p_max, n_step, time_jacobi, "jacobi", numberOfIterations, converge_test);

//-------------------------------------------------------------

//        eigvec_arma.print();
//        cout << "---" << endl;
//        eigvec_jacobi.print();
//        cout << "eigenvaules:" << endl;
//        eigval_arma.print();
//        cout << "----" << endl;
//        eigval_jacobi.print();

        cout << "diff. eigenvec:" << endl;
        //cout << size(eigvec_arma) << size(eigvec_jacobi) << endl;
        vec a = eigvec_arma.col(0) - eigvec_jacobi.col(0);
        a.print();
        cout << "----" << endl;
        cout << "diff. eigval:" << endl;
        vec b = eigval_arma - eigval_jacobi;
        b.print();
    }

    return 0;
}

