#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <time.h>

using namespace std;
using namespace arma;

void find_max_elem_index(int& k, int& l, double& max_off_diagonal, const mat &B, const int n_step)
{
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

void transformation_matrix(double& c, double& s, const mat &B, const int k, const int l)
{
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
        cout << "Element chosen is zero: B(" << k << "," << l << ") = 0 ==> c = 1 and s = 0"
    }
}

void jacobi_rotation(mat& B, const double c, const double s, const int k, const int l, const int n_step)
{
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

int main()
{
    const int n_step = 100;
    const double p_max = 5; //writing p instead of rho
    const double p_min = 0;

    const double h = (p_max - p_min)/n_step;
    vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h
    vec V = p%p;

//-------------------------------------------------------------
    //Constructing test matrix
    mat B = ones<mat>(n_step-1, n_step-1);
    for(int i=0; (i<n_step-1); ++i)
    {
        B(i,i) = 0.0;
    }
    //B.print();


//-------------------------------------------------------------

    // Constructing B:
  /*  double e = -1/(h*h); // all elements of the e vec is the same
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
    //solving equations with armadillo lib:
    vec eigval = eig_sym(B);

//-------------------------------------------------------------
    //solving equations with jacobi rotation:
    double tolerance = 1.0e-08;
    double max_off_diagonal = 1.0;
    int k = 0;
    int l = 0;
    find_max_elem_index(k, l, max_off_diagonal, B, n_step); //initial val to enter loop
    int maxNumberOfIterations = 10000;

    int numberOfIterations = 0;
    while(tolerance < max_off_diagonal)
    {
        // cout << endl << endl << "Iteration " << numberOfIterations << " with matrix:" << endl;
        // B.print();

        if(++numberOfIterations > maxNumberOfIterations) {
            cout << "Jacobi algorithm did not converge after " << maxNumberOfIterations << " iterations. Aborting!" << endl;
            exit(1);
        }

        //finding the value and index(k,l) of the maximum element in B:
        // cout << "Finding maximum offdiagonal element..." << endl;
        find_max_elem_index(k, l, max_off_diagonal, B, n_step);
        // if( (numberOfIterations % 1000) == 0 ) {
        // cout << "Found element " << k << ", " << l << " with value " << max_off_diagonal << endl;
        // }

        //finding the values of c ans s (the S transformation matrix):
        double c = 0;
        double s = 0;
        // cout << "Finding c and s" << endl;
        transformation_matrix(c, s, B, k, l);
        // cout << "c=" << c << "  s  = " << s << endl;

        //transformation of B:
        // cout << "Now rotating..." << endl;
        jacobi_rotation(B, c, s, k, l, n_step);
        // cout << "After rotating, matrix looks like: " << endl;
        // B.print();
        //cout << "----" << endl;
    }
    //retriving eigenvalues:
    cout << "------" << endl;
    vec a = B.diag();
    a = sort(a);
    //eigval = sort(eigval);

    //a.print();

//-------------------------------------------------------------
    //comparing with old matrix:
    cout << "------" << endl;
    //B.print();
    cout << "------" << endl;

//-------------------------------------------------------------
    //comparing with arma lib eigval
    cout << a[0] << endl;
    cout << a[1] << endl;
    cout << a[2] << endl;
    cout << "------" << endl;
    cout << eigval[0] << endl;
    cout << eigval[1] << endl;
    cout << eigval[2] << endl;
    //B.print();

    //eigval.print();

//-------------------------------------------------------------
    //comparing with start out diagonal
    //d.shed_row(0);
    //d.shed_row(n_step-1);
    //cout << "------" << endl;
    //(a-d).print();


    return 0;
}

