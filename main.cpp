#include "functions.h"


#define Pi 3.14159265358979323846

using namespace std;

using namespace Eigen;

int main()
{
    cout<<"hello hello"<<endl;

    functions* MyFunc = new functions();

    ///Initializing parameters
    MyFunc->intialize();

    int nelx=MyFunc->getnelx(),nely=MyFunc->getnely();
    int n1,n2;
    double volfrac=MyFunc->getvolfrac(),penal=MyFunc->getpenal();

    ArrayXXd xold,x=MatrixXd::Constant(nely,nelx,volfrac);
    VectorXd Ue(8);
    MatrixXd KE,dc;
    int loop=0;
    double change=1.0,c=0.0;
    SpVec U;

    /// Start iteration
    while (change>0.01)
    {
        loop++;
        xold=x;

        U=MyFunc->FE(x);

        /// Objective function and sensitivity analysis
        KE=lk();

        for (int ely=0; ely<nely; ++ely)
        {
            for (int elx=0; elx<nelx; ++elx)
            {
                n1=(nely+1)*elx+ely+1;
                n2=(nely+1)*(elx+1)+ely+1;

                Ue(0)=U.coeff(2*n1-2);
                Ue(1)=U.coeff(2*n1-1);
                Ue(2)=U.coeff(2*n2-2);
                Ue(3)=U.coeff(2*n2-1);
                Ue(4)=U.coeff(2*n2);
                Ue(5)=U.coeff(2*n2+1);
                Ue(6)=U.coeff(2*n1);
                Ue(7)=U.coeff(2*n1+1);

                c+=pow(x(ely,elx),penal)*Ue.dot(KE*Ue);
                dc(ely,elx)=-penal*pow(x(ely,elx),penal-1.0)*Ue.dot(KE*Ue);
            };
        };


        change=(x-xold).maxCoeff();
        cout<<"It.: "<<loop<<" Obj.: "<<c<<" Vol.: "<<x.mean()<<" ch.: "<<change<<endl;
    }

    return 0;
}
