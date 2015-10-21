#ifndef JACOBI
#define JACOBI

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <time.h>

using namespace std;
using namespace arma;

class jacobi
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

    // constructor goes here

    //functions:
    void find_max_elem_index(int& index_k, int& index_l, double& value_maximum_off_diagonal, const mat &first_matrix, const int number_of_steps)
    {
        k = index_k;
        l = index_l;
        max_off_diagonal = value_maximum_off_diagonal;
        B = first_matrix;
        n_step = number_of_steps;

        max_off_diagonal = -1.0;
        //Checking all off-diagonal elements:
        for(int i=0; i<n_step-1; i++)  {
            for(int j=i+1; j<n_step-1; j++) {
                //storing value and index of max off-diag-element
                if(fabs(B(i,j)) > max_off_diagonal)
                {
                    max_off_diagonal = fabs(B(i,j));
                    k = i;
                    l = j;
                }
            }
        }
        if(max_off_diagonal < 0) {
            cout << "Warning, find_max_elem_index did not find any offdiagonal element with absolute value larger than zero. Aborting!" << endl;
            exit(1);
        }
    }

    void transformation_matrix(double& cos, double& sin, const mat &first_matrix, const int index_k, const int index_l)
    {
        c = cos;
        s = sin;
        B = first_matrix;
        k = index_k;
        l = index_l;

        if(B(k,l) != 0)
        {
            double t = 0;
            double tau = (B(l,l) - B(k,k))/(2.0*B(k,l));
            // ensuring that theta<=pi/4:
            if(tau>0)
            {
                t = -tau + sqrt(1 + (tau*tau));
            }
            else
            {
                t = -tau - sqrt(1 + (tau*tau));
            }
            c = 1.0/sqrt(1.0 + (t*t));
            s = t*c;
        }
        else
        {
            c = 1.0;
            s = 0.0;
            cout << "Element chosen is zero: B(" << k << "," << l << ") = 0 ==> c = 1 and s = 0" << endl;
        }
    }

    void jacobi_rotation(mat& first_matrix, const double cos, const double sin, const int index_k, const int index_l, const int number_of_steps)
    {
        c = cos;
        s = sin;
        k = index_k;
        l = index_l;
        B = first_matrix;
        n_step = number_of_steps;

        double B_kk = B(k,k);
        double B_ll = B(l,l);
        B(k,k) = B_kk*c*c - 2*B(k,l)*c*s + B_ll*s*s;
        B(l,l) = B_ll*c*c + 2*B(k,l)*c*s + B_kk*s*s;
        B(k,l) = 0.0;//(B_kk - B_ll)*c*s + B(k,l)*(c*c - s*s); //or just set to 0
        B(l,k) = B(k,l); //symetrix matrix

        //changing the remaining elements:
        for(int i=0; i<n_step-1; ++i)
        {
            if( (i!=k) && (i!=l) )
            {
                double B_ik = B(i,k);
                double B_il = B(i,l);
                B(i,k) = B_ik*c - B_il*s;
                B(k,i) = B(i,k);
                B(i,l) = B_il*c + B_ik*s;
                B(l,i) = B(i,l);
            }
        }
    }

    int solve_eq_jacobi_rotation(mat& first_matrix, const int number_of_steps, const int max_nr_it, int& test_conv)
    {
        B = first_matrix;
        n_step = number_of_steps;
        maxNumberOfIterations = max_nr_it;
        converge_test = test_conv;

        double tolerance = 1.0e-08;
        double max_off_diagonal = 0;
        int k = 0;
        int l = 0;

        //instantiating an object of class jacobi method

        find_max_elem_index(k, l, max_off_diagonal, B, n_step); //initial value for max_off_diagonal to enter loop

        int numberOfIterations = 0;
        while(tolerance < max_off_diagonal)
        {
            if(++numberOfIterations > maxNumberOfIterations) {
                cout << "Jacobi algorithm did not converge after " << maxNumberOfIterations << " iterations for n_step= " << n_step << ". Exiting jacobi rotation solver!" << endl;
                converge_test = 0;
                B = zeros<mat>(n_step-1,n_step-1);
                break;
            }
            //finding the value and index(k,l) of the maximum element in B:
            find_max_elem_index(k, l, max_off_diagonal, B, n_step);

            //finding the values of c ans s (the S transformation matrix):
            double c = 0;
            double s = 0;
            transformation_matrix(c, s, B, k, l);
            //transformation of B:
            jacobi_rotation(B, c, s, k, l, n_step);

        }
        cout << "numb: " << numberOfIterations << endl;
        return numberOfIterations;
    }

};

#endif // JACOBI

