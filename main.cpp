#include "functions.h"

using namespace std;

using namespace Eigen;

int main()
{
    functions* MyFunc = new functions();

    ///Initializing parameters
    MyFunc->intialize();
    /* the beam and stiffness matrix are constant
    they can be initialized outside the loop */
    MyFunc->defbeam();
    MyFunc->lk();

    int nelx=MyFunc->getnelx(),nely=MyFunc->getnely();
    int n1,n2;
    double volfrac=MyFunc->getvolfrac(),penal=MyFunc->getpenal();

    ArrayXXd xold,x=ArrayXXd::Constant(nely,nelx,volfrac),dc(nely,nelx);
    VectorXd Ue(8),U_dense;
    MatrixXd KE=MyFunc->getKE();
    int loop=0;
    double change=1.0,c;

    std::clock_t start;
    double duration;

    /// Start iteration
    while (change>0.01)
    {
        loop++;
        xold=x;
        c=0.0;

        /// FE resolution (timed to assess performance)
        start = std::clock();
        U_dense=MyFunc->FE_dense(x);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        cout<<"FE time "<<duration<<endl;

        /// Objective function and sensitivity analysis

        for (int ely=0; ely<nely; ++ely)
        {
            for (int elx=0; elx<nelx; ++elx)
            {
                n1=(nely+1)*elx+ely+1;
                n2=(nely+1)*(elx+1)+ely+1;

                Ue(0)=U_dense(2*n1-2);
                Ue(1)=U_dense(2*n1-1);
                Ue(2)=U_dense(2*n2-2);
                Ue(3)=U_dense(2*n2-1);
                Ue(4)=U_dense(2*n2);
                Ue(5)=U_dense(2*n2+1);
                Ue(6)=U_dense(2*n1);
                Ue(7)=U_dense(2*n1+1);

                c+=pow(x(ely,elx),penal)*Ue.dot(KE*Ue);
                dc(ely,elx)=-penal*pow(x(ely,elx),penal-1.0)*Ue.dot(KE*Ue);
            };
        };

        /// filtering
        dc=MyFunc->check(x,dc);

        ///design update by the OC method
        x=MyFunc->OC(x,dc);

        /// print results
        change=(x-xold).abs().maxCoeff();

        cout<<"It.: "<<loop<<" Obj.: "<<c<<" Vol.: "<<x.mean()<<" ch.: "<<change<<endl;
    }

    delete MyFunc;

    return 0;
}
