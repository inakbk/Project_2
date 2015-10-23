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
    //overloading constructor, one for jacobi and one for arma, this for jacobi:
    writetofile(const vec eigenvalues, const vec eigenvector_gs, const int p_max, const int n_step, const double time, string FileName, const int number_of_iterations, const int converge_test)
    {
        ofstream myfile;
            string filename = "EigenValVecSolver_" + FileName + "_pMax" + to_string(p_max) + "_nStep" + to_string(n_step) + ".txt";
            myfile.open (filename);
            myfile << "Equations solved with the " << FileName << " algorithm." << endl;
            myfile << "Dimention of matrix + 1, n_step = " << n_step << endl;
            myfile << "Value of p_max: " << p_max << endl;
            myfile << "Execution time: " << time << endl;

            myfile << "Number of iterations for jacobi algoritm: " << number_of_iterations << endl;

            //Writing eigenvalues/vectors to file if the jacobi method did converge.
            if(converge_test == true)
            {
                myfile << "---------------------" << endl;
                myfile << "Eigenvalues (sorted, maximum writing 10 eigenvalues to file)" << endl;
                //myfile << "---------------------" << endl;
                int number_of_eigenvalues_printed = 10;
                for (int i=0; i < number_of_eigenvalues_printed; i++)
                {
                    if(i == size(eigenvalues,0))
                    {
                        cout << "There is less than 10 eigenvalues for the " << FileName << " solver, printed only the first to file. Exiting loop." << endl;
                        break;
                    }
                    myfile << eigenvalues[i] << endl;
                }
                myfile << "---------------------" << endl;
                myfile << "Eigenvector ground state" << endl;
                for (int i=0; i < n_step-2; i++)
                {
                    myfile << eigenvector_gs[i] << endl;
                }
                myfile.close();
            }
            //Writing error message to file if jacobi method did not converge.
            if(converge_test == false)
            {
                myfile << endl; myfile << endl;
                myfile << "Jacobi method did not converge after " << number_of_iterations << " iterations!" << endl;
                myfile << endl; myfile << endl;

                cout << "Warning! Jacobi method did not converge after " << number_of_iterations << " iterations!" << endl;
                cout << endl;
            }
            cout << "Datafile done for n_step=" << n_step << " with the " << FileName <<" solver." << endl;
            cout << endl;
    }

    //this for arma
    writetofile(const vec eigenvalues, const vec eigenvector_gs, const int p_max, const int n_step, const double time, string FileName)
    {
        ofstream myfile;
            string filename = "EigenValVecSolver_" + FileName + "_pMax" + to_string(p_max) + "_nStep" + to_string(n_step) + ".txt";
            myfile.open (filename);
            myfile << "Equations solved with the " << FileName << " algorithm." << endl;
            myfile << "Dimention of matrix + 1, n_step = " << n_step << endl;
            myfile << "Value of p_max: " << p_max << endl;
            myfile << "Execution time: " << time << endl;

            myfile << endl; //number of iterations for jacobi goes here

            myfile << "---------------------" << endl;
            myfile << "Eigenvalues (sorted, maximum writing 10 eigenvalues to file)" << endl;
            //myfile << "---------------------" << endl;
            int number_of_eigenvalues_printed = 10;
            for (int i=0; i < number_of_eigenvalues_printed; i++)
            {
                if(i == size(eigenvalues,0))
                {
                    cout << "There is less than 10 eigenvalues for the " << FileName << " solver, printed only the first to file. Exiting loop." << endl;
                    break;
                }
                myfile << eigenvalues[i] << endl;
            }
            myfile << "---------------------" << endl;
            myfile << "Eigenvector ground state" << endl;
            for (int i=0; i < n_step-2; i++)
            {
                myfile << eigenvector_gs[i] << endl;
            }
            myfile.close();

            cout << "Datafile done for n_step= " << n_step << " with the " << FileName <<" solver." << endl;
            cout << endl;
    }
};

#endif // WRITETOFILE

