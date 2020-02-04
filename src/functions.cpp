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

ArrayXXd functions::OC(const ArrayXXd &x,const MatrixXd &dc)
{
    ArrayXXd xnew(_nely,_nelx);
    ArrayXXd dcn=dc.array(),xmin;
    double moov=0.2;
    double l1=0.0,l2=100000,lmid;

    while (l2-l1>1.e-4)
    {
        lmid=0.5*(l2+l1);
        xnew=(x+moov).min(x*(-dcn/lmid).sqrt());
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

SpVec functions::FE(const ArrayXXd &x) ///FE solver using sparse matrices
{
    MatrixXd KE=lk();
    /* note : block writing operations for sparse matrices
    are not available with Eigen (unless the columns are contiguous) */
    //MatrixXd K=MatrixXd::Zero(2*(_nelx+1)*(_nely+1),2*(_nelx+1)*(_nely+1));
    SpMat K(2*(_nelx+1)*(_nely+1),2*(_nelx+1)*(_nely+1));
    SpVec F(2*(_nelx+1)*(_nely+1)), U(2*(_nelx+1)*(_nely+1));
    int n1,n2, ind1,ind2;
    int edof[8];
    SparseLU<SpMat, COLAMDOrdering<int> > solver;
    //std::vector<int> edof2;
    //Eigen::ArrayXi edof3(8);
    std::vector<T> tripletList;
    tripletList.reserve(64);

    for (int ely=0; ely<_nely; ++ely)
    {
        for (int elx=0; elx<_nelx; ++elx)
        {
            n1=(_nely+1)*elx+ely+1;
            n2=(_nely+1)*(elx+1)+ely+1;
            //cout<<n1<<'\t'<<n2<<endl;
            /// defining indexes to write in
            edof[0]=2*n1-2;
            edof[1]=2*n1-1;
            edof[2]=2*n2-2;
            edof[3]=2*n2-1;
            edof[4]=2*n2;
            edof[5]=2*n2+1;
            edof[6]=2*n1;
            edof[7]=2*n1+1;

            //K(edof,edof)=K(edof,edof)+pow(x(ely,elx),_penal)*KE;

            ///K matrix (slow) assembly

            for (int i=0; i<8; ++i)
            {
                for (int j=0; j<8; ++j)
                {
                    ind1=edof[i];
                    ind2=edof[j];
                    //cout<<ind1<<'\t'<<ind2<<endl;
                    //K.coeffRef(ind1,ind2) += pow(x(ely,elx),_penal)*KE(i,j);
                    tripletList.push_back(T(ind1,ind2,K.coeff(ind1,ind2)+pow(x(ely,elx),_penal)*KE(i,j)));
                };
            };

            K.setFromTriplets(tripletList.begin(), tripletList.end());

        };
    };

    //cout<<K.nonZeros()<<endl;
    //cout<<K.coeff(0,0)<<'\t'<<K.coeff(1,0)<<'\t'<<K.coeff(2,0)<<'\t'<<K.coeff(3,0)<<endl;
    //cout<<K.coeff(20,64)<<'\t'<<K.coeff(21,64)<<'\t'<<K.coeff(22,64)<<'\t'<<K.coeff(23,64)<<endl;

    F.insert(1,0)=-1.0;
    /// determining indices of dofs which need solving
    VectorXi fixeddofs(_nely+2);
    fixeddofs.head(_nely+1)=VectorXi::LinSpaced(_nely+1,0,2*_nely+1);
    VectorXi alldofs=VectorXi::LinSpaced(2*(_nely+1)*(_nelx+1),0,2*(_nely+1)*(_nelx+1)-1);
    VectorXi freedofs(2*(_nely+1)*(_nelx+1));

    auto it = std::set_difference(alldofs.data(), alldofs.data() + alldofs.size(),
                fixeddofs.data(), fixeddofs.data() + fixeddofs.size(),freedofs.data());

    freedofs.conservativeResize(std::distance(freedofs.data(), it)); // resize the result
    //cout<<freedofs(0)<<endl;
    //cout<<freedofs(103)<<endl;


    /// creating smaller matrices for solving system
    int Nfree=freedofs.size();
    //cout<<Nfree<<endl;
    SpMat K_sub(Nfree,Nfree);
    SpVec U_sub(Nfree), F_sub(Nfree);
    std::vector<T> tripletU,tripletF,tripletK;
    tripletK.reserve(Nfree*Nfree);

    for (int i=0; i<Nfree; ++i)
    {
        F_sub.insert(i)=F.coeff(freedofs(i));

        for (int j=0; j<Nfree; ++j)
        {
            tripletK.push_back(T(i,j,K.coeff(freedofs(i),freedofs(j))));
        };
    };

    K_sub.setFromTriplets(tripletK.begin(), tripletK.end());
    K_sub.makeCompressed(); //necessary for proper factorization

    /// Solving system using LU decomposition
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(K_sub);
    // Compute the numerical factorization
    solver.factorize(K_sub);
    //Use the factors to solve the linear system
    U_sub = solver.solve(F_sub);

    /// Casting solution into U vector
    for (int i=0; i<Nfree; ++i)
    {
        U.insert(freedofs(i))=U_sub.coeff(i);
    }

    return U;
}

VectorXd functions::FE_dense(const ArrayXXd &x) ///FE solver using dense matrices
{
    MatrixXd KE=lk();
    /* note : block writing operations for sparse matrices
    are not available with Eigen (unless the columns are contiguous) */
    MatrixXd K=MatrixXd::Zero(2*(_nelx+1)*(_nely+1),2*(_nelx+1)*(_nely+1));
    VectorXd F(2*(_nelx+1)*(_nely+1)), U(2*(_nelx+1)*(_nely+1));
    int n1,n2, ind1,ind2;
    int edof[8];
    //std::vector<int> edof2;
    //Eigen::ArrayXi edof3(8);

    for (int ely=0; ely<_nely; ++ely)
    {
        for (int elx=0; elx<_nelx; ++elx)
        {
            n1=(_nely+1)*elx+ely+1;
            n2=(_nely+1)*(elx+1)+ely+1;
            //cout<<n1<<'\t'<<n2<<endl;
            /// defining indexes to write in
            edof[0]=2*n1-2;
            edof[1]=2*n1-1;
            edof[2]=2*n2-2;
            edof[3]=2*n2-1;
            edof[4]=2*n2;
            edof[5]=2*n2+1;
            edof[6]=2*n1;
            edof[7]=2*n1+1;

            //K(edof,edof)+=pow(x(ely,elx),_penal)*KE;

            ///K matrix (slow) assembly

            for (int i=0; i<8; ++i)
            {
                for (int j=0; j<8; ++j)
                {
                    ind1=edof[i];
                    ind2=edof[j];
                    //cout<<ind1<<'\t'<<ind2<<endl;
                    K(ind1,ind2) += pow(x(ely,elx),_penal)*KE(i,j);
                };
            };
        };
    };

    //cout<<K.nonZeros()<<endl;
    //cout<<K(0,0)<<'\t'<<K(1,0)<<'\t'<<K(2,0)<<'\t'<<K(3,0)<<endl;
    //cout<<K(20,64)<<'\t'<<K(21,64)<<'\t'<<K(22,64)<<'\t'<<K(23,64)<<endl;

    F(1,0)=-1.0;
    /// determining indices of dofs which need solving
    VectorXi fixeddofs(_nely+2);
    fixeddofs.head(_nely+1)=VectorXi::LinSpaced(_nely+1,0,2*_nely+1);
    VectorXi alldofs=VectorXi::LinSpaced(2*(_nely+1)*(_nelx+1),0,2*(_nely+1)*(_nelx+1)-1);
    VectorXi freedofs(2*(_nely+1)*(_nelx+1));

    auto it = std::set_difference(alldofs.data(), alldofs.data() + alldofs.size(),
                fixeddofs.data(), fixeddofs.data() + fixeddofs.size(),freedofs.data());

    freedofs.conservativeResize(std::distance(freedofs.data(), it)); // resize the result
    //cout<<freedofs(0)<<endl;
    //cout<<freedofs(103)<<endl;


    /// creating smaller matrices for solving system
    int Nfree=freedofs.size();
    //cout<<Nfree<<endl;
    MatrixXd K_sub(Nfree,Nfree);
    VectorXd U_sub(Nfree), F_sub(Nfree);

    for (int i=0; i<Nfree; ++i)
    {
        F_sub(i)=F(freedofs(i));

        for (int j=0; j<Nfree; ++j)
        {
            K_sub(i,j)=K(freedofs(i),freedofs(j));
        };
    };

    /// Solving system using LU decomposition
    // Initializing solver
    Eigen::FullPivLU<MatrixXd> solver(K_sub);
    // Compute the numerical factorization
    K_sub=solver.matrixLU();
    //Use the factors to solve the linear system
    U_sub = solver.solve(F_sub);

    /// Casting solution into U vector
    for (int i=0; i<Nfree; ++i)
    {
        U(freedofs(i))=U_sub(i);
    }

    return U;
}

MatrixXd functions::check(const ArrayXXd &x,const MatrixXd &dc)
{
    MatrixXd dcn=MatrixXd::Zero(_nely,_nelx);
    double sum=0.0,fac;
    int mi,mj,Mi,Mj;

    for (int i=0; i<_nelx; ++i)
    {
        Mi=max(i-int(round(_rmin)),0);
        mi=min(i+int(round(_rmin)),_nelx);
        for (int j=0; j<_nely; ++j)
        {
            Mj=max(j-int(round(_rmin)),0);
            mj=min(j+int(round(_rmin)),_nely);
            for (int k=Mi; k<mi; ++k)
            {
                for (int l=Mj; l<mj; ++l)
                {
                    fac=_rmin-sqrt(double((i-k)*(i-k)+(j-l)*(j-l)));
                    sum+=max(0.0,fac);
                    dcn(j,i)+=max(0.0,fac)*x(l,k)*dc(l,k);
                };
            };

            dcn(j,i)/=(x(j,i)*sum);

        };
    };

    return dcn;
}

MatrixXd lk()
{
    double E=1.0, nu=0.3;
    std::array<double,8> k;
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

    KE(1,2)=2.0*k[7];
    KE(1,3)=2.0*k[6];
    KE(1,4)=2.0*k[5];
    KE(1,5)=2.0*k[4];
    KE(1,6)=2.0*k[3];
    KE(1,7)=2.0*k[2];

    KE(2,3)=2.0*k[5];
    KE(2,4)=2.0*k[6];
    KE(2,5)=2.0*k[3];
    KE(2,6)=2.0*k[4];
    KE(2,7)=2.0*k[1];

    KE(3,4)=2.0*k[7];
    KE(3,5)=2.0*k[2];
    KE(3,6)=2.0*k[1];
    KE(3,7)=2.0*k[4];

    KE(4,5)=2.0*k[1];
    KE(4,6)=2.0*k[2];
    KE(4,7)=2.0*k[3];

    KE(5,6)=2.0*k[7];
    KE(5,7)=2.0*k[6];

    KE(6,7)=2.0*k[5];

    //matrix is symmetric
    MatrixXd KEt=KE.transpose();
    KE+=KEt;
    KE*=0.5;
    KE*=E/(1.0-nu*nu);

    return KE;
}
