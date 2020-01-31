#include "functions.h"

functions::functions()
{
    //ctor
}

functions::~functions()
{
    //dtor
}

void functions::setnelx()
{
    std::cout<<"enter ratio of elements in the horizontal direction"<<endl;
    cin>>_nelx;
}

int functions::getnelx()
{
    return _nelx;
}

void functions::setnely()
{
    cout<<"enter ratio of elements in the vertical direction"<<endl;
    cin>>_nely;
}

int functions::getnely()
{
    return _nely;
}

void functions::setvolfrac()
{
    cout<<"enter volume fraction"<<endl;
    cin>>_volfrac;
}

double functions::getvolfrac()
{
    return _volfrac;
}

void functions::setpenal()
{
    cout<<"enter penalizing power"<<endl;
    cin>>_penal;
}

double functions::getpenal()
{
    return _penal;
}

void functions::setrmin()
{
    cout<<"enter filter size"<<endl;
    cin>>_rmin;
}

double functions::getrmin()
{
    return _rmin;
}

void functions::intialize()
{
    setnelx();
    setnely();
    setvolfrac();
    setpenal();
    setrmin();
}

ArrayXXd functions::OC(ArrayXXd &x,const double &dc)
{
    ArrayXXd xnew(_nely,_nelx);
    ArrayXXd xmin;
    double moov=0.2;
    double l1=0.0,l2=100000,lmid;

    while (l2-l1>1.e-4)
    {
        lmid=0.5*(l2+l1);
        xnew=(x+moov).min(x*sqrt(-dc/lmid));
        xnew=xnew.max(x-moov);
        xnew=xnew.max(0.001);
        if (xnew.sum()-_volfrac*_nelx*_nely>0)
        {
            l1=lmid;
        }
        else
        {
            l2=lmid;
        }

    }

    return xnew;
}

VectorXd functions::FE()
{

}

MatrixXd lk()
{
    double E=1.0, nu=0.3;
    double k[8];
    k[0]=0.5-nu/6.0;
    k[1]=1.0/8.0+nu/8.0;
    k[2]=-0.25-nu/12.0;
    k[3]=-1.0/8.0+3.0*nu/8.0;
    k[4]=-0.25+nu/12.0;
    k[5]=-k[1];
    k[6]=nu/6.0;
    k[7]=-k[3];

    MatrixXd KE=k[0]*MatrixXd::Identity(8,8);

    KE(0,1)=2.0*k[1];
    KE(0,2)=2.0*k[2];
    KE(0,3)=2.0*k[3];
    KE(0,4)=2.0*k[4];
    KE(0,5)=2.0*k[5];
    KE(0,6)=2.0*k[6];
    KE(0,7)=2.0*k[7];

    KE(1,2)=2.0*k[8];
    KE(1,3)=2.0*k[7];
    KE(1,4)=2.0*k[6];
    KE(1,5)=2.0*k[5];
    KE(1,6)=2.0*k[4];
    KE(1,7)=2.0*k[3];

    KE(2,3)=2.0*k[6];
    KE(2,4)=2.0*k[7];
    KE(2,5)=2.0*k[4];
    KE(2,6)=2.0*k[5];
    KE(2,7)=2.0*k[2];

    KE(3,4)=2.0*k[8];
    KE(3,5)=2.0*k[3];
    KE(3,6)=2.0*k[2];
    KE(3,7)=2.0*k[5];

    KE(4,5)=2.0*k[2];
    KE(4,6)=2.0*k[3];
    KE(4,7)=2.0*k[4];

    KE(5,6)=2.0*k[8];
    KE(5,7)=2.0*k[7];

    KE(6,7)=2.0*k[6];

    //matrix is symmetric
    KE+=KE.transpose();
    KE*=0.5;
    KE*=E/(1.0-nu*nu);
    return KE;
}
