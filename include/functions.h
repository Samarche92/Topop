#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Eigen/Dense"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

using namespace Eigen;

class functions
{
    public:
        functions();
        virtual ~functions();

        void setnelx();
        int getnelx();

        void setnely();
        int getnely();

        void setvolfrac();
        double getvolfrac();

        void setpenal();
        double getpenal();

        void setrmin();
        double getrmin();

        void intialize();
        // initializing user-defined parameters

        void defbeam();
        // define loads and support

        void lk(const double& E=1.0,const double& nu=0.3);
        /* building element stiffness matrix
        input (optional) : Young modulus E and Poisson ratio nu */

        MatrixXd& getKE();

        ArrayXXd& OC(const ArrayXXd &x,const ArrayXXd &dc);
        /* optimality criteria update
        input : original array of design variables and filtered sensitivities
        output : updated array of design variables */

        VectorXd& FE_dense(const ArrayXXd &x);
        /* FE resolution using dense matrices
        input : array of design variables
        output : global displacement vector */

        ArrayXXd& check(const ArrayXXd &x,const ArrayXXd &dc);
        /* mesh-independecy filter
        input : array of design variables and original sensitivities
        output : modified sensitivites */

    protected:

    private:

        int _nelx,_nely; // ratio of horizontal and vertical elements
        double _volfrac,_penal,_rmin; // volume fraction, penalization power and filter size

        VectorXi freedofs; // free degrees of freedom
        MatrixXd KE; // stiffness matrix
};

#endif // FUNCTIONS_H
