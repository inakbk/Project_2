#ifndef WRITETOFILE
#define WRITETOFILE

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>

using namespace std;
using namespace arma;

class writetofile
{
public:
    //variables
    vec eigenvalues;
    int p_max;
    int n_step;
    int number_of_iterations;
    double time;
    string FileName;
    int converge_test;

    // constructor goes here
    //writetofile()

    void funcWriteToFile(const vec& eigenvalues, int p_max, const int n_step, const int number_of_iterations, const double time, string FileName, const int converge_test)
    {
        ofstream myfile;
            string filename = "EigvalSolver_" + FileName + "_pMax" + to_string(p_max) + "_nStep" + to_string(n_step) + ".txt";
            myfile.open (filename);
            myfile << "Solution of the eigenvalueproblem for the " << FileName << " algorithm." << endl;
            myfile << "Dimention of matrix, one less than n_step: " << n_step << endl;
            myfile << "Value of p_max: " << p_max << endl;
            myfile << "Execution time: " << time << endl;
            if(number_of_iterations != -1){
                myfile << "Number of iterations for jacobi algoritm: " << number_of_iterations << endl;
            }
            else
            {
                myfile << endl;
            }
            //Writing eigenvalues to file if the jacobi method did converge.
            if(converge_test == true)
            {
                myfile << "---------------------" << endl;
                myfile << "Eigenvalues (sorted)" << "  "<< "" << "     " << "Eigenvectors (soon)" << endl;
                myfile << "---------------------" << endl;
                int number_of_eigenvalues_printed = 10;
                for (int i=0; i < number_of_eigenvalues_printed; i++)
                {
                    if(i == size(eigenvalues,0))
                    {
                        cout << "Length of eigenval vec is shorter than number_of_eigenvalues_printed for " << FileName << " solver, exiting loop." << endl;
                        break;
                    }
                    myfile << eigenvalues[i] << "    " << "--" << endl;
                }
                myfile << "(Maximum writing " << number_of_eigenvalues_printed << " eigenvalues to file.)" << endl;
                myfile.close();

                cout << "Datafile done for n_step=" << n_step << " with the " << FileName <<" solver." << endl;
                cout << endl;
            }
            //Writing error message til file if jacobi method did not converge.
            if(converge_test == false)
            {
                myfile << endl; myfile << endl;
                myfile << "Jacobi method did not converge after " << number_of_iterations << " iterations!" << endl;
                myfile << endl; myfile << endl;

                cout << "Datafile done for n_step=" << n_step << " with the " << FileName <<" solver which did not converge." << endl;
                cout << endl;
            }
    }

};



#endif // WRITETOFILE

