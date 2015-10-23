#ifndef JACOBISOLVER
#define JACOBISOLVER

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <time.h>
#include <jacobi.h>

class jacobisolver
{
public:
    int n_step;
    int maxNumberOfIterations;
    mat B;
    mat R;

    // constructor, initiating input values and sets default values if not any given
    jacobisolver(mat firstMatrix = zeros<mat>(2,2), int numberOfStep = 2, int maxIterations=1)
    {
        if(maxIterations == 1)
        {
            cout << "Not enough arguments given to jacobisolver. Aborting." << endl;
            exit(1);
        }
        B = firstMatrix;
        n_step = numberOfStep;
        maxNumberOfIterations = maxIterations;
        //Eigenvector initiated as identity
        R = eye<mat>(n_step-1, n_step-1);
    }

    //solving and timing the whole jacobi algorithm:
    void solve_w_jacobi_rotation(vec& eigval_jacobi, mat& eigvec_jacobi, mat& R, int& numberOfIterations, int& converge_test, double time_jacobi)
    {
        clock_t start_jacobi, finish_jacobi;
        start_jacobi = clock();

        iterations(eigval_jacobi, R, numberOfIterations, converge_test);
        retriveSortEigenValVec(eigval_jacobi, eigvec_jacobi, B, R);

        finish_jacobi = clock();
        time_jacobi = ( (finish_jacobi - start_jacobi)/((double)CLOCKS_PER_SEC ) );
    }

    void iterations(vec& eigval_jacobi, mat& R, int& numberOfIterations, int& converge_test)
    {
        double tolerance = 1.0e-14;
        double max_off_diagonal = 0;
        int k = 0;
        int l = 0;

        //instantiating an object of class jacobi:
        jacobi partsOfJacobiMethod;

        partsOfJacobiMethod.find_max_elem_index(k, l, max_off_diagonal, B, n_step); //initial value for max_off_diagonal to enter loop

        numberOfIterations = 0;
        while(tolerance < max_off_diagonal)
        {
            if(++numberOfIterations > maxNumberOfIterations) {
                cout << "Jacobi algorithm did not converge after " << maxNumberOfIterations << " iterations for n_step= " << n_step << ". Exiting jacobi rotation solver!" << endl;
                converge_test = 0;
                B = zeros<mat>(n_step-1,n_step-1);
                break;
            }
            //finding the value and index(k,l) of the maximum element in B:
            partsOfJacobiMethod.find_max_elem_index(k, l, max_off_diagonal, B, n_step);

            //finding the values of c ans s (the S transformation matrix):
            double c = 0;
            double s = 0;
            partsOfJacobiMethod.transformation_matrix(c, s, B, k, l);
            //doing one transformation:
            partsOfJacobiMethod.jacobi_rotation(B, R, c, s, k, l, n_step);

        }
    }

    void retriveSortEigenValVec(vec& eigval_jacobi, mat& eigvec_jacobi, mat B, mat& R)
    {
        eigval_jacobi = B.diag();
        eigval_jacobi = sort(eigval_jacobi);
        uvec indexes_sorted_eigenval = sort_index(B.diag());
        for(int i=0; i<=n_step-2; ++i)
        {
            eigvec_jacobi.col(i) = R.col(indexes_sorted_eigenval(i));
        }
        cout << "eigvec jacobi" << endl;
    }


};

#endif // JACOBISOLVER

