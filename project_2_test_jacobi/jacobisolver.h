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

    // constructor goes here
    //jacobi()

    //methods:
};






#endif // JACOBISOLVER

