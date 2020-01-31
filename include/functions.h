#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "Eigen/Dense"

#include <time.h>
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <thread>
#include <string>
#include <math.h>
#include <cmath>
#include <complex>

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

        ArrayXXd OC(ArrayXXd &x,const double &dc);

        VectorXd FE();

    protected:

    private:

        int _nelx,_nely;
        double _volfrac,_penal,_rmin;
};

MatrixXd lk();

#endif // FUNCTIONS_H