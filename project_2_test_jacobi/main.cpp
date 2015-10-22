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

        mat B = ones<mat>(n_step-1,n_step-1);
        for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
        {
            B(i,i) = 0;
        }
        B(n_step-2,n_step-2) = 0;
        B(1,2) = 9;
        //B(2,1) = 11;
        //B.print();

//-------------------------------------------------------------
        //clocking the operations:
        clock_t start_arma, finish_arma;
        start_arma = clock();

        //solving equations with armadillo lib:
        vec eigval_arma = eig_sym(B);

        finish_arma = clock();
        double time_arma = ( (finish_arma - start_arma)/((double)CLOCKS_PER_SEC ) );

        int converge_test = 1; //this is only false if jacobi did not converge (always true for arma)
        writetofile fileArma;
        fileArma.funcWriteToFile(eigval_arma, p_max, n_step, -1, time_arma, "arma", converge_test);

//-------------------------------------------------------------
        vec eigval_jacobi = zeros<vec>(n_step-1);
        int numberOfIterations = 0;
        double time_jacobi = 0;
        //solving equations with jacobi rotation:
        jacobisolver solveTestMatrix(B, n_step, maxNumberOfIterations);
        solveTestMatrix.solve_w_jacobi_rotation(eigval_jacobi, numberOfIterations, converge_test, time_jacobi);

        cout << "Total number of iterations with Jacobi rotation: " << numberOfIterations << endl;

        writetofile fileJacobi;
        fileJacobi.funcWriteToFile(eigval_jacobi, p_max, n_step, numberOfIterations, time_jacobi, "jacobi", converge_test);

//-------------------------------------------------------------
        eigval_arma.print();
        cout << "----" << endl;
        eigval_jacobi.print();
    }

    return 0;
}

