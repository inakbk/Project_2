//Program to find the magnitude of p_max

#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <time.h>

using namespace std;
using namespace arma;

int main()
{
    int n = 10;
    vec rho_max = linspace(4,6,n);
    mat eigen_values = zeros<mat>(3,n);
    mat A = zeros<mat>(3,n);

    for(int i=0; i<=n-1; i++)
    {
        const int n_step = 1000; //?
        const double p_max = rho_max[i];
        const double p_min = 0;

        const double h = (p_max - p_min)/n_step;
        vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h
        vec V = p%p;

        // initializing B
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

        vec eigval = eig_sym(B);
        eigen_values(0,i) = eigval[0];
        A(0,i) = 3;
        eigen_values(1,i) = eigval[1];
        A(1,i) = 7;
        eigen_values(2,i) = eigval[2];
        A(2,i) = 11;
        //eigval.print();
    }

    //A.print();
    mat C = A - eigen_values;
    //eigen_values.print();

    C.print();

}
