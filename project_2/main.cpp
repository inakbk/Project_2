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
    //d.print();

    mat B = zeros<mat>(n_step-1,n_step-1);
    for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
    {
        B(i,j) = e;
        B(i,i) = d[i+1];
        B(j,i) = e;
    }
    B(n_step-2,n_step-2) = d[n_step-1];
    B.print();

    double tolerance = 0.00000001;
    double max_off_diagonal = B(0,1)*B(0,1)-1;
    int k = 0;
    int l = 0;

    //this test is not working properly, or the change to the matrix is not working
    //while(tolerance < max_off_diagonal)
    //{
        //finding the value and index(k,l) of the maximum element in B:
        for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
        {
            //Checking all off-diagonal elements:
            if( ((B(i,j)*B(i,j)) > max_off_diagonal) && (i!=j) )
            {
                max_off_diagonal = B(i,j)*B(i,j);
                k = i;
                l = j;
                cout << "here2" << endl;
            }
        }
        //cout << k << endl; //this value does not change?
        //cout << l << endl;
        //cout << max_off_diagonal << endl;
        cout << "----------" << endl;

        //finding the values of c ans s (the S matrix):
        double t = 0;
        double tau = (B(l,l) - B(k,k))/(2*B(k,l)); //=cot(2*theta)

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

        cout << "This should be zero:" << endl;
        cout << (B(k,k) - B(l,l) )*c*s + B(k,l)*(c*c - s*s) << endl;


        //think this is not working, B is unchanged after the loop...
        for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
        {
            B(i,i) = B(i,i); //dont need to write?
            double temp_B = B(i,k);
            B(i,k) = B(i,k)*c - B(i,l)*s;
            B(i,l) = B(i,l)*c - temp_B*s;
            temp_B = B(k,k);
        }
        B(k,k) = B(k,k)*c*c - 2*B(k,l)*c*s + B(l,l)*s*s;
        B(l,l) = B(l,l)*c*c - 2*B(k,l)*c*s + temp_B*s*s;
        B(k,l) = (B(k,k) - B(l,l) )*c*s + B(k,l)*(c*c - s*s);
        B(l,k) = (B(k,k) - B(l,l) )*c*s + B(k,l)*(c*c - s*s);
    //}
    B.print();

    return 0;
}

