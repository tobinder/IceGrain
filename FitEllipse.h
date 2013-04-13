/*! \file FitEllipse.h
 *  \brief Class representing an ellipse fit.
 */ 
//**************************************************************************
// This java code is  an interactive demo of the first ellipse-specific
// direct fitting method presented in the papers:

//    M. Pilu, A. Fitzgibbon, R.Fisher ``Ellipse-specific Direct
//    least-square Fitting '' , IEEE International Conference on Image
//    Processing, Lausanne, September 1996. (poscript) (HTML) 

//    A. Fitzgibbon, M. Pilu , R.Fisher ``Direct least-square fitting of
//    Ellipses '' , International Conference on Pattern Recognition, Vienna,
//    August 1996. (poscript) - Extended version available as DAI Research
//    Paper #794

// The demo can be tried out at   
//   http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/PILU1/demo.html

// The code was written by Maurizio Pilu , University of Edinburgh.

// Some math routines are from the "Numerical Recipes in C" by
// Press/Teukolsky/Vettering/Flannery, Cambridge Uiniversity Press,
// Second Edition (1988). PLEASE READ COPYRIGHT ISSUES ON THE NUMERICAL
// RECIPES BOOK.

// NOTE: Some parts of the program are rather scruffy. The author
//       will tidy it up whan he has some spare time.

// DISCLAIMER: The authors and the department assume no responsabilities
//             whatsoever for any wrong use of this code. 

// COPYRIGHT: Any commercial use of the code and the method is forbidden
//            without written authorization from the authors.
//
// NOTE: Translated to C++ by Tobias Binder, University of Heidelberg
//**************************************************************************

class FitEllipse
{
    public:

    void fit(std::vector<point>, std::vector<float> &);
    void drawConic(std::vector<float>, int, std::vector<point> &);

    private:

    void rotate(double **, int, int, int, int, double, double);
    void jacobi(double **, int, double *, double **, int);
    void choldc(double **, int, double **);
    int inverse(double **, double **, int);
    void AperB(double **, double **, double **, int, int, int, int);
    void A_TperB(double **, double **, double **, int, int, int, int);
    void AperB_T(double **, double **, double **, int, int, int, int);
    void pv(double *, std::string);
    void pm(double **, std::string);
};

//fits ellpise with a0 * x^2 + a1 * xy + a2 * y^2 + a3 * x + a4 * y + a5 = 0
void FitEllipse::fit(std::vector<point> points, std::vector<float> & pvec)
{
    int np = points.size();// number of points

    double ** D = new double * [np+1];
    for (int dim=0; dim<np+1; dim++)
    {
        D[dim] = new double[7];
    }

    double ** S = new double * [7];
    double ** Const = new double * [7];
    double ** temp = new double * [7];
    double ** L = new double * [7]; 
    double ** C = new double * [7]; 
    double ** invL = new double * [7];
    double ** V = new double * [7]; 
    double ** sol = new double * [7];
    for (int dim=0; dim<7; dim++)
    {
        S[dim] = new double[7];
        Const[dim] = new double[7];
        for (int dim2=0; dim2<7; dim2++)
            Const[dim][dim2]=0.0;

        temp[dim] = new double[7];
        L[dim] = new double[7]; 
        C[dim] = new double[7]; 
        invL[dim] = new double[7];
        V[dim] = new double[7]; 
        sol[dim] = new double[7];
    }

    double * d = new double[7];
    double tx,ty;
    int nrot=0;

    Const[1][3]=-2;
    Const[2][2]=1;
    Const[3][1]=-2;	

    	if (np<6)
    {
        for (int dim=0; dim<np+1; dim++)
        {
            delete D[dim];
        }
        delete D;

        for (int dim=0; dim<7; dim++)
        {
            delete S[dim];
            delete Const[dim];
            delete temp[dim];
            delete L[dim];
            delete C[dim]; 
            delete invL[dim];
            delete V[dim]; 
            delete sol[dim];
        }
        delete S;
        delete Const;
        delete temp;
        delete L; 
        delete C; 
        delete invL;
        delete V; 
        delete sol;

        delete d;
 
        return;
    }

    // Now first fill design matrix
    for (int i=1; i <= np; i++)
    { 
        tx = points[i-1].x;
        ty = points[i-1].y;
        D[i][1] = tx*tx;
        D[i][2] = tx*ty;
        D[i][3] = ty*ty;
        D[i][4] = tx;
        D[i][5] = ty;
        D[i][6] = 1.0;
    }

    //pm(Const,"Constraint");
    // Now compute scatter matrix  S
    A_TperB(D,D,S,np,6,np,6);
    //pm(S,"Scatter");

    choldc(S,6,L);    
    //pm(L,"Cholesky");

    inverse(L,invL,6);
    	//pm(invL,"Inverse");

    AperB_T(Const,invL,temp,6,6,6,6);
    AperB(invL,temp,C,6,6,6,6);
    //pm(C,"The C matrix");

    jacobi(C,6,d,V,nrot);
    //pm(V,"The eigenvectors");  //OK
    //pv(d,"The eigenvalues");

    A_TperB(invL,V,sol,6,6,6,6);
    //pm(sol,"The GEV solution unnormalized");  //SOl

    // Now normalize them 
    for (int j=1;j<=6;j++)  //Scan columns
    {
        double mod = 0.0;
        for (int i=1;i<=6;i++)
            mod += sol[i][j]*sol[i][j];
            for (int i=1;i<=6;i++)
                sol[i][j] /=  sqrt(mod); 
    }

    //pm(sol,"The GEV solution");  //SOl

    double zero=10e-20;
    double minev=10e+20;
    int  solind=0;

    for (int i=1; i<=6; i++)
        if (d[i]<0 && fabs(d[i])>zero)	
            solind = i;

    // Now fetch the right solution
    for (int j=1;j<=6;j++)
        pvec.push_back(sol[j][solind]);

    //std::cout<<"------------The solution--------------"<<std::endl;
    //std::cout<<" "<<pvec[0]<<" "<<pvec[1]<<" "<<pvec[2]<<" "<<pvec[3]<<" "<<pvec[4]<<" "<<pvec[5]<<std::endl;

    for (int dim=0; dim<np+1; dim++)
    {
        delete D[dim];
    }
    delete D;

    for (int dim=0; dim<7; dim++)
    {
        delete S[dim];
        delete Const[dim];
        delete temp[dim];
        delete L[dim];
        delete C[dim]; 
        delete invL[dim];
        delete V[dim]; 
        delete sol[dim];
    }
    delete S;
    delete Const;
    delete temp;
    delete L; 
    delete C; 
    delete invL;
    delete V; 
    delete sol;

    delete d;
}

//calculates points on the fitted ellipse (number of points given)
void FitEllipse::drawConic(std::vector<float> pvec, int nptsk, std::vector<point> & XY)
{
    XY.resize(nptsk);
    int npts=nptsk/2;	

    double ** u = new double * [3];
    double ** Aiu = new double * [3];
    double ** L = new double * [3];
    double ** B = new double * [3];
    double ** Xpos = new double * [3];
    double ** Xneg = new double * [3];
    double ** ss1 = new double * [3];
    double ** ss2 = new double * [3];
    double ** uAiu = new double * [3];
    double ** A = new double * [3];
    double ** Ai = new double * [3];
    double ** Aib = new double * [3];
    double ** b = new double * [3];
    double ** r1 = new double * [2];
    for (int dim=0; dim<3; dim++)
    {
        u[dim] = new double[npts+1];
        Aiu[dim] = new double[npts+1];
        L[dim] = new double[npts+1];
        B[dim] = new double[npts+1];
        Xpos[dim] = new double[npts+1];
        Xneg[dim] = new double[npts+1];
        ss1[dim] = new double[npts+1];
        ss2[dim] = new double[npts+1];
        uAiu[dim] = new double[npts+1];
        A[dim] = new double[3];
        Ai[dim] = new double[3];
        Aib[dim] = new double[2];
        b[dim] = new double[2];
        if (dim<2) r1[dim] = new double[2];
    }

    double * lambda = new double[npts+1];
    float Ao, Ax, Ay, Axx, Ayy, Axy;
         
    double pi = 3.141592654f;      
    double theta;
    int i;
    int j;
    double kk;
  
    Ao = pvec[5];
    Ax = pvec[3];
    Ay = pvec[4];
    Axx = pvec[0];
    Ayy = pvec[2];
    Axy = pvec[1];

    A[1][1] = Axx;    A[1][2] = Axy/2;
    A[2][1] = Axy/2;  A[2][2] = Ayy;
    b[1][1] = Ax; b[2][1] = Ay;  

    // Generate normals linspace
    for (i=1, theta=0.0; i<=npts; i++, theta+=(pi/npts))
    {
        u[1][i] = cos(theta);
        u[2][i] = sin(theta);
    }
  
    inverse(A,Ai,2);
  
    AperB(Ai,b,Aib,2,2,2,1);
    A_TperB(b,Aib,r1,2,1,2,1);      
    r1[1][1] = r1[1][1] - 4*Ao;

    AperB(Ai,u,Aiu,2,2,2,npts);
    for (i=1; i<=2; i++)
        for (j=1; j<=npts; j++)
            uAiu[i][j] = u[i][j] * Aiu[i][j];

    for (j=1; j<=npts; j++)
    {
        if ( (kk=(r1[1][1] / (uAiu[1][j]+uAiu[2][j]))) >= 0.0)
            lambda[j] = sqrt(kk);
        else
            lambda[j] = -1.0;
    }

    // Builds up B and L
    for (j=1; j<=npts; j++)
        L[1][j] = L[2][j] = lambda[j];      
    for (j=1; j<=npts; j++)
    {
        B[1][j] = b[1][1];
        B[2][j] = b[2][1];
    }
  
    for (j=1; j<=npts; j++)
    {
        ss1[1][j] = 0.5 * (  L[1][j] * u[1][j] - B[1][j]);
        ss1[2][j] = 0.5 * (  L[2][j] * u[2][j] - B[2][j]);
        ss2[1][j] = 0.5 * ( -L[1][j] * u[1][j] - B[1][j]);
        ss2[2][j] = 0.5 * ( -L[2][j] * u[2][j] - B[2][j]);
    }

    AperB(Ai,ss1,Xpos,2,2,2,npts);
    AperB(Ai,ss2,Xneg,2,2,2,npts);

    for (j=1; j<=npts; j++)
    {
        if (lambda[j]==-1.0)
        {
            XY[j-1].x = -1;
            XY[j-1].y = -1;
            XY[j-1+npts].x = -1;
            XY[j-1+npts].y = -1;
        }
        else
        {
            XY[j-1].x = Xpos[1][j];
            XY[j-1].y = Xpos[2][j];
            XY[j-1+npts].x = Xneg[1][j];
            XY[j-1+npts].y = Xneg[2][j];
        }	  	                 
    }

    for (int dim=0; dim<3; dim++)
    {
        delete u[dim];
        delete Aiu[dim];
        delete L[dim];
        delete B[dim];
        delete Xpos[dim];
        delete Xneg[dim];
        delete ss1[dim];
        delete ss2[dim];
        delete uAiu[dim];
        delete A[dim];
        delete Ai[dim];
        delete Aib[dim];
        delete b[dim];
        if (dim<2) delete r1[dim];
    }
    delete u;
    delete Aiu;
    delete L;
    delete B;
    delete Xpos;
    delete Xneg;
    delete ss1;
    delete ss2;
    delete uAiu;
    delete A;
    delete Ai;
    delete Aib;
    delete b;
    delete r1;

    delete lambda;
}

void FitEllipse::rotate(double ** a, int i, int j, int k, int l, double tau, double s) 
{
    	double g,h;
    	g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);
    	a[k][l]=h+s*(g-h*tau);
}
    
void FitEllipse::jacobi(double ** a, int n, double * d , double ** v, int nrot)      
{
    	int j,iq,ip,i;
    	double tresh,theta,tau,t,sm,s,h,g,c;

    	double * b = new double[n+1];
    	double * z = new double[n+1];

    	for (ip=1;ip<=n;ip++)
    {
        for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
    	}

    	for (ip=1;ip<=n;ip++)
    {
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
    	}

    	nrot=0;
    	for (i=1;i<=50;i++)
    {
        sm=0.0;
        for (ip=1;ip<=n-1;ip++)
        {
            for (iq=ip+1;iq<=n;iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0)
        {
        	    /* free_vector(z,1,n);
            free_vector(b,1,n);  */

            delete b;
            	delete z;
            return;
        }
        if (i < 4)
            tresh=0.2*sm/(n*n);
        else
        	    tresh=0.0;
        for (ip=1;ip<=n-1;ip++)
        {
            for (iq=ip+1;iq<=n;iq++)
            {
                g=100.0*fabs(a[ip][iq]);
                if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
	                && fabs(d[iq])+g == fabs(d[iq]))
	                a[ip][iq]=0.0;
                else if (fabs(a[ip][iq]) > tresh)
                {
	                h=d[iq]-d[ip];
	                if (fabs(h)+g == fabs(h))
	                    t=(a[ip][iq])/h;
	                else
                    {
	                    theta=0.5*h/(a[ip][iq]);
	                    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	                    if (theta < 0.0) t = -t;
	                }
	                c=1.0/sqrt(1+t*t);
	                s=t*c;
	                tau=s/(1.0+c);
	                h=t*a[ip][iq];
	                z[ip] -= h;
                    z[iq] += h;
	                d[ip] -= h;
	                d[iq] += h;
	                a[ip][iq]=0.0;
                    for (j=1;j<=ip-1;j++)
                    {
	                    rotate(a,j,ip,j,iq,tau,s);
	                }
	                for (j=ip+1;j<=iq-1;j++)
                    {
	                    rotate(a,ip,j,j,iq,tau,s);
	                }
	                for (j=iq+1;j<=n;j++)
                    {
	                    rotate(a,ip,j,iq,j,tau,s);
	                }
	                for (j=1;j<=n;j++)
                    {
	                    rotate(v,j,ip,j,iq,tau,s);
	                }
	                ++nrot;
                }
            }
        }
        for (ip=1;ip<=n;ip++)
        {
            b[ip] += z[ip];
            d[ip]=b[ip];
            z[ip]=0.0;
        }
    }
    //std::cout<<"Too many iterations in routine JACOBI"<<std::endl;

    delete b;
    	delete z;
}  

//  Perform the Cholesky decomposition    
// Return the lower triangular L  such that L*L'=A  
void FitEllipse::choldc(double ** a, int n, double ** l)
{
    int i,j,k;
    double sum;
    double * p = new double[n+1];

    for (i=1; i<=n; i++)
    {
        for (j=i; j<=n; j++)
        {
            for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
            if (i == j)
            {
                if (sum<=0.0)  
	            //std::cout<<"A is not poitive definite!"<<std::endl;
	            {}
                else 
	                p[i]=sqrt(sum);
            }
            else 
            {
	            a[j][i]=sum/p[i];
            }
        }
    }       
    for (i=1; i<=n; i++)  
        for (j=i; j<=n; j++)  
            if (i==j)
                l[i][i] = p[i];
            else
            {
	            l[j][i]=a[j][i];  
	            l[i][j]=0.0;
            }    

    delete p;
}

//Calcola la inversa della matrice  B mettendo il risultato
//in InvB . Il metodo usato per l'inversione e' quello di
//Gauss-Jordan.   N e' l'ordine della matrice .
//ritorna 0 se l'inversione  corretta altrimenti ritorna
//SINGULAR .
int FitEllipse::inverse(double ** TB, double ** InvB, int N)
{  
    int k,i,j,p,q;
    double mult;
    double D,temp;
    double maxpivot;
    int npivot;
    double ** B = new double * [N+1];
    double ** A = new double * [N+1];
    double ** C = new double * [N+1];
    for (int dim=0; dim<N+1; dim++)
    {
        B[dim] = new double[N+2];
        A[dim] = new double[2*N+2];
        C[dim] = new double[N+1];
    }
    double eps = 10e-20;
  
    for(k=1;k<=N;k++)
        for(j=1;j<=N;j++)
            B[k][j]=TB[k][j];
  
    for (k=1;k<=N;k++)
    {
        for (j=1;j<=N+1;j++)
            A[k][j]=B[k][j];
        for (j=N+2;j<=2*N+1;j++)
            A[k][j]=(float)0;
        A[k][k-1+N+2]=(float)1;
    }
    for (k=1;k<=N;k++)
    {
        maxpivot=fabs((double)A[k][k]);
        npivot=k;
        for (i=k;i<=N;i++)
        if (maxpivot<fabs((double)A[i][k]))
        {
	        maxpivot=fabs((double)A[i][k]);
	        npivot=i;
        }
        if (maxpivot>=eps)
        {
            if (npivot!=k)
	            for (j=k;j<=2*N+1;j++)
	            {
		            temp=A[npivot][j];
		            A[npivot][j]=A[k][j];
		            A[k][j]=temp;
	            } ;
	        D=A[k][k];
	        for (j=2*N+1;j>=k;j--)
	            A[k][j]=A[k][j]/D;
	        for (i=1;i<=N;i++)
	        {
	            if (i!=k)
		        {
		            mult=A[i][k];
		            for (j=2*N+1;j>=k;j--)
		                A[i][j]=A[i][j]-mult*A[k][j] ;
		        }
	        }
	    }
        else
        {  //std::cout<<"The matrix may be singular !!"<<std::endl;
            for (int dim=0; dim<N+1; dim++)
            {
                delete B[dim];
                delete A[dim];
                delete C[dim];
            }
            delete B;
            delete A;
            delete C;

            return(-1);
        };
    }
    /**   Copia il risultato nella matrice InvB  ***/
    for (k=1,p=1;k<=N;k++,p++)
        for (j=N+2,q=1;j<=2*N+1;j++,q++)
            InvB[p][q]=A[k][j];

    for (int dim=0; dim<N+1; dim++)
    {
        delete B[dim];
        delete A[dim];
        delete C[dim];
    }
    delete B;
    delete A;
    delete C;

    return(0);
}
    
void FitEllipse::AperB(double ** _A, double ** _B, double ** _res, int _righA, int _colA, int _righB, int _colB)
{
    int p,q,l;                                      
    for (p=1;p<=_righA;p++)                        
        for (q=1;q<=_colB;q++)                        
        {
            _res[p][q]=0.0;                            
            for (l=1;l<=_colA;l++)                     
                _res[p][q]=_res[p][q]+_A[p][l]*_B[l][q];  
        }                                            
}                                                 

void FitEllipse::A_TperB(double ** _A, double ** _B, double ** _res, int _righA, int _colA, int _righB, int _colB)
{
    int p,q,l;                                      
    for (p=1;p<=_colA;p++)                        
        for (q=1;q<=_colB;q++)                        
        {
            _res[p][q]=0.0;                            
            for (l=1;l<=_righA;l++)                    
                _res[p][q]=_res[p][q]+_A[l][p]*_B[l][q];  
        }                                            
}

void FitEllipse::AperB_T(double ** _A, double ** _B, double ** _res, int _righA, int _colA, int _righB, int _colB)
{
    int p,q,l;                                      
    for (p=1;p<=_colA;p++)                         
        for (q=1;q<=_colB;q++)                        
        {
            _res[p][q]=0.0;                            
            for (l=1;l<=_righA;l++)                    
                _res[p][q]=_res[p][q]+_A[p][l]*_B[q][l];  
        }                                            
}

void FitEllipse::pv(double * v, std::string str)
{
    std::cout<<"------------"<<str<<"--------------"<<std::endl;
    std::cout<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<" "<<v[4]<<" "<<v[5]<<" "<<v[6]<<std::endl;
}
    
void FitEllipse::pm(double ** S, std::string str)
{
    std::cout<<"------------"<<str<<"--------------"<<std::endl;
    std::cout<<" "<<S[1][1]<<" "<<S[1][2]<<" "<<S[1][3]<<" "<<S[1][4]<<" "<<S[1][5]<<" "<<S[1][6]<<std::endl;
    std::cout<<" "<<S[2][1]<<" "<<S[2][2]<<" "<<S[2][3]<<" "<<S[2][4]<<" "<<S[2][5]<<" "<<S[2][6]<<std::endl;
    std::cout<<" "<<S[3][1]<<" "<<S[3][2]<<" "<<S[3][3]<<" "<<S[3][4]<<" "<<S[3][5]<<" "<<S[3][6]<<std::endl;
    std::cout<<" "<<S[4][1]<<" "<<S[4][2]<<" "<<S[4][3]<<" "<<S[4][4]<<" "<<S[4][5]<<" "<<S[4][6]<<std::endl;
    std::cout<<" "<<S[5][1]<<" "<<S[5][2]<<" "<<S[5][3]<<" "<<S[5][4]<<" "<<S[5][5]<<" "<<S[5][6]<<std::endl;
    std::cout<<" "<<S[6][1]<<" "<<S[6][2]<<" "<<S[6][3]<<" "<<S[6][4]<<" "<<S[6][5]<<" "<<S[6][6]<<std::endl;
}
