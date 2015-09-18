#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <time.h>


using namespace std;
using namespace arma;

int main()
{
    const int n_step = 11; //?
    const double p_max = 10.0; //writing p instead of rho
    const double p_min = 0;

    const double h = (p_max - p_min)/n_step;
    vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h
    vec V = p%p;

    // initializing B
    double e = -1/(h*h); // all elements of the e vec is the same
    vec d = 2/(h*h) + V;
    d.print();

    mat B = zeros<mat>(n_step-1,n_step-1);
    for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
    {
        B(i,j) = e;
        B(i,i) = d[i+1];
        B(j,i) = e;
    }
    B(n_step-2,n_step-2) = d[n_step-1];
    B.print();
/*

    ////
    double t = 0;
    double tau = (a[l][l] - a[k][k])/(2*a[k][l]); //=cot(2*theta)

    // ensuring that theta<=pi/4:
    if(tau>0)
    {
        t = -tau + sqrt(1 + (tau*tau)); //set sqrt outside to compute only once?
    }
    if(tau<0)
    {
        t = -tau - sqrt(1 + (tau*tau));
    }
    double c = 1/sqrt(1 + (t*t));
    double s = t*c;

    ///make a loop over i,l,k:

    b[i][i]*/

    return 0;
}

