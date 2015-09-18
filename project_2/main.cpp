#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <time.h>


using namespace std;
using namespace arma;

int main()
{
    int n_step = 11; //?
    double p_max = 10.0; //writing p instead of rho
    double p_min = 0;

    double h = (p_max - p_min)/n_step;
    vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h


    ////
    double t = 0;
    double tau = (a[l][l] - a[k][k])/(2*a[k][l]); //=cot(2*theta)
    if(tau>0)
    {
        t = -tau + sqrt(1 + (tau*tau)); //set sqrt outside to compute only once?

    }
    if(tau<0)
    {
        t = -tau - sqrt(1 + (tau*tau));
    }


    theta = //..
    s = sin(theta);
    c = cos(theta);

    return 0;
}

