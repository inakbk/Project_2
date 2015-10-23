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

    //methods:
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
                t = 1.0/(tau + sqrt(1.0 + tau*tau));
            }
            else
            {
                t = -1.0/(-tau + sqrt(1.0 + tau*tau));
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

    void jacobi_rotation(mat& B, mat& R, const double c, const double s, const int k, const int l, const int n_step)
    {
        double B_kk = B(k,k);
        double B_ll = B(l,l);
        B(k,k) = B_kk*c*c - 2*B(k,l)*c*s + B_ll*s*s;
        B(l,l) = B_ll*c*c + 2*B(k,l)*c*s + B_kk*s*s;
        B(k,l) = 0.0;//(B_kk - B_ll)*c*s + B(k,l)*(c*c - s*s) just set to 0 exact
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

        double r_ik = R(i,k);
        double r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
        }
    }
};

#endif // JACOBI

