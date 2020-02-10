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

        void defbeam();

        void lk();
        MatrixXd getKE();

        ArrayXXd OC(const ArrayXXd &x,const ArrayXXd &dc);

        VectorXd FE_dense(const ArrayXXd &x);

        ArrayXXd check(const ArrayXXd &x,const ArrayXXd &dc);

    protected:

    private:

        int _nelx,_nely;
        double _volfrac,_penal,_rmin;

        VectorXi freedofs; // free degrees of freedom
        MatrixXd KE; // stiffness matrix
};

#endif // FUNCTIONS_H
