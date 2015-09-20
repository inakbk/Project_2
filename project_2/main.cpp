#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <time.h>

using namespace std;
using namespace arma;

void find_max_elem_index(int& k, int& l, double& max_off_diagonal, mat B, int n_step)
{
    //Checking all off-diagonal elements:
    for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
    {
        //storing value and index of max off-diag-element
        if( ((B(i,j)*B(i,j)) > max_off_diagonal) && (i!=j) )
        {
            max_off_diagonal = B(i,j)*B(i,j);
            k = i;
            l = j;
        }
    }
}

void transformation_matrix(double& c, double& s, const mat B, const int k, const int l)
{
    double t = 0;
    if(B(k,l) != 0)
    {
        double tau = (B(l,l) - B(k,k))/(2*B(k,l));
        // ensuring that theta<=pi/4:
        if(tau>0)
        {
            t = -tau + sqrt(1 + (tau*tau));
        }
        else
        {
            t = -tau - sqrt(1 + (tau*tau)); //corect choice of +- tau?
        }
        c = 1/sqrt(1 + (t*t));
        s = t*c;
    }
    else
    {
        c = 1.0;
        s = 0.0;
    }
}

void jacobi_rotation(mat& B, const double c, const double s, const int k, const int l, const int n_step)
{
    double B_kk = B(k,k);
    double B_ll = B(l,l);
    double B_ik, B_il;
    B(k,k) = B_kk*c*c - 2*B(k,l)*c*s + B_ll*s*s;
    B(l,l) = B_ll*c*c + 2*B(k,l)*c*s + B_kk*s*s;
    B(k,l) = 0.0;//(B_kk - B_ll)*c*s + B(k,l)*(c*c - s*s); //or just set to 0
    B(l,k) = B(k,l); //symetrix matrix

    //changing the remaining elements:
    for(int i=0; i<=n_step-2; ++i)
    {
        if((i!=k) && (i!=l))
        {
            B_ik = B(i,k);
            B_il = B(i,l);
            B(i,k) = B_ik*c - B_il*s;
            B(k,i) = B(i,k);
            B(i,l) = B_il*c - B_ik*s;
            B(l,i) = B(i,l);
        }
    }
}

int main()
{
    const int n_step = 10; //?
    const double p_max = 5; //writing p instead of rho
    const double p_min = 0;

    const double h = (p_max - p_min)/n_step;
    vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h
    cout << p[1]-p[0] << h << endl;
    vec V = p%p;

//-------------------------------------------------------------
    //Constructing test matrix
    mat B = ones<mat>(n_step-1, n_step-1);
    for(int i=0; (i<=n_step-2); ++i)
    {
        B(i,i) = 0.0;
    }
    //B.print();


//-------------------------------------------------------------
/*
    // Constructing B:
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
    //B.print();
    */
//-------------------------------------------------------------
    //solving with armadillo lib:
    vec eigval = eig_sym(B);
    //eigval.print();

//-------------------------------------------------------------
    //algorithm with jacobi rotation:
    double tolerance = 1.0e-08;
    double max_off_diagonal = tolerance + 0.1; //initial val to enter loop
    int k = 0;
    int l = 0;

    while(tolerance < max_off_diagonal)
    {
        //finding the value and index(k,l) of the maximum element in B:
        max_off_diagonal = B(0,1)*B(0,1)-1; //why do I need this one? initializing?

        find_max_elem_index(k, l, max_off_diagonal, B, n_step);

        //finding the values of c ans s (the S transformation matrix):
        double c = 0;
        double s = 0;
        transformation_matrix(c, s, B, k, l);

        //transformation of B:
        jacobi_rotation(B, c, s, k, l, n_step);
    }
    //retriving eigenvalues:
    cout << "------" << endl;
    vec a = B.diag();
    a = sort(a);
    //eigval = sort(eigval);
    a.print();

//-------------------------------------------------------------
    //comparing with old matrix:
    cout << "------" << endl;
    //B.print();
    cout << "------" << endl;

//-------------------------------------------------------------
    //comparing with arma lib eigval
    cout << "------" << endl;
    //cout << a[0] << endl;
    //cout << eigval[0] << endl;
    eigval.print();

//-------------------------------------------------------------
    //comparing with start out diagonal
    //d.shed_row(0);
    //d.shed_row(n_step-1);
    //cout << "------" << endl;
    //(a-d).print();


    return 0;
}

