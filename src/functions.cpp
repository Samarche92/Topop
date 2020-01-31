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
