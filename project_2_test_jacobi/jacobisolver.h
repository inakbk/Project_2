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
    //variables:
    int k;
    int l;
    double c;
    double s;
    double max_off_diagonal;
    int n_step;
    mat B;
    int maxNumberOfIterations;
    int converge_test;

    //double time;
    vec eigval_jacobi;
    int numberOfIterations;

    // constructor, initiate values and sets default values if not any given
    jacobisolver(mat firstMatrix = zeros<mat>(2,2), int numberOfStep = 2, int maxIterations=1)
    {
        B = firstMatrix;
        n_step = numberOfStep;
        maxNumberOfIterations = maxIterations;
    }

    //methods:
    void solve_w_jacobi_rotation(vec& eigval_jacobi, int& numberOfIterations, int& converge_test)
    {
        double tolerance = 1.0e-08;
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
            partsOfJacobiMethod.jacobi_rotation(B, c, s, k, l, n_step);

        }

        retriveSortEigenvalues(eigval_jacobi, B);
    }

    void retriveSortEigenvalues(vec& eigval_jacobi, mat B)
    {
        eigval_jacobi = B.diag();
        eigval_jacobi = sort(eigval_jacobi);
    }


};

#endif // JACOBISOLVER

