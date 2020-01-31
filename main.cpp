#include "functions.h"


#define Pi 3.14159265358979323846

using namespace std;

using namespace Eigen;

/*typedef Matrix<std::complex<double>, 3, 3> MatriX3c;
typedef Matrix<std::complex<double>, 3, 1> Vectorx3c;
typedef Matrix<std::complex<double>,Dynamic,Dynamic> MatrixXdc;
typedef Matrix<std::complex<double>, Dynamic, 1> VectorXdc;*/


int main()
{
    cout<<"hello hello"<<endl;

    functions* MyFunc = new functions();

    ///Initializing parameters
    MyFunc->intialize();

    int nelx=MyFunc->getnelx(),nely=MyFunc->getnely();
    double volfrac=MyFunc->getvolfrac();

    MatrixXd xold,x=MatrixXd::Constant(nely,nelx,volfrac);
    int loop=0;
    double change=1.0,c=0.0;

    while (change>0.01)
    {
        loop++;
        xold=x;

        change=(x-xold).maxCoeff();
        cout<<"It.: "<<loop<<" Obj.: "<<c<<" Vol.: "<<x.mean()<<" ch.: "<<change<<endl;
    }

    return 0;
}
