#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "blas.h"
#include "lapack.h"
#include "mex.h"

const size_t ione = 1;
const size_t itwo = 2;
const size_t ithree = 3;
const size_t iseven = 7;
const double fone = 1;
const double ftwo = 2;
const double fzero = 0;
const double fmone = -1;
size_t quiet = 0;

void printmat(double *A, size_t m, size_t n);
void fmpcsolve(double *A, double *B, double *At, double *Bt, double *eyen, 
        double *eyem, double *Q, double *R, double *Qf, double *g, double *w, double *zmax, double *zmin,
        double *x, double *z, size_t T, size_t n, size_t m, size_t nz, size_t niters, double kappa);

void gfgphp(double *Q, double *R, double *Qf, double *g, double *zmax, double *zmin, double *z,
        size_t T, size_t n, size_t m, size_t nz, double *gf, double *gp, double *hp);

void rdrp(double *A, double *B, double *Q, double *R, double *Qf, double *z, double *nu, 
        double *gf, double *gp, double *b, size_t T, size_t n, size_t m, size_t nz, 
        double kappa, double *rd, double *rp, double *Ctnu);

void resdresp(double *rd, double *rp, size_t T, size_t n, size_t nz, double *resd, 
        double *resp, double *res);
        
void dnudz(double *A, double *B, double *At, double *Bt, double *eyen, 
        double *eyem, double *Q, double *R, double *Qf, double *hp, double *rd, 
        double *rp, size_t T, size_t n, size_t m, size_t nz, double kappa, double *dnu, double *dz);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* problem setup */
    size_t i, j, m, n, nz, T, niters;
    double kappa;
    double *dptr, *dptr1, *dptr2;
    double *A, *B, *At, *Bt, *Q, *R, *Qf, *xmax, *xmin, *umax, *umin, *x;
    double *q = NULL, *qf = NULL, *r = NULL, *w = NULL;
    double *zmax, *zmin, *zmaxp, *zminp, *X, *U, *z, *eyen, *eyem, *x0, *g;
    double *X0, *U0;
    double *telapsed;
    clock_t t1, t2;

    A = mxGetPr(mxGetField(prhs[0],0,"A"));
    B = mxGetPr(mxGetField(prhs[0],0,"B"));
    Q = mxGetPr(mxGetField(prhs[0],0,"Q"));
    R = mxGetPr(mxGetField(prhs[0],0,"R"));
    Qf = mxGetPr(mxGetField(prhs[1],0,"Qf"));
    xmax = mxGetPr(mxGetField(prhs[0],0,"xmax"));
    xmin = mxGetPr(mxGetField(prhs[0],0,"xmin"));
    umax = mxGetPr(mxGetField(prhs[0],0,"umax"));
    umin = mxGetPr(mxGetField(prhs[0],0,"umin"));
    n = (int)mxGetScalar(mxGetField(prhs[0],0,"n"));
    m = (int)mxGetScalar(mxGetField(prhs[0],0,"m"));
    q = mxGetPr(mxGetField(prhs[0],0,"q"));
    qf = mxGetPr(mxGetField(prhs[0],0,"qf"));
    r = mxGetPr(mxGetField(prhs[0],0,"r"));
    w = mxGetPr(mxGetField(prhs[0],0,"w"));

    T = (int)mxGetScalar(mxGetField(prhs[1],0,"T"));
    kappa = (double)mxGetScalar(mxGetField(prhs[1],0,"kappa"));
    niters = (int)mxGetScalar(mxGetField(prhs[1],0,"niters"));
    quiet = (int)mxGetScalar(mxGetField(prhs[1],0,"quiet"));

    X0 = mxGetPr(prhs[2]);
    U0 = mxGetPr(prhs[3]);
    x0 = mxGetPr(prhs[4]);

    nz = T*(n+m);

    /* outputs */
    plhs[0] = mxCreateDoubleMatrix(n,T,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m,T,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    X = mxGetPr(plhs[0]);
    U = mxGetPr(plhs[1]);
    telapsed = mxGetPr(plhs[2]);

    At = malloc(sizeof(double)*n*n);
    Bt = malloc(sizeof(double)*n*m);
    eyen = malloc(sizeof(double)*n*n);
    eyem = malloc(sizeof(double)*m*m);
    z = malloc(sizeof(double)*nz);
    x = malloc(sizeof(double)*n);
    zmax = malloc(sizeof(double)*nz);
    zmin = malloc(sizeof(double)*nz);
    zmaxp = malloc(sizeof(double)*nz);
    zminp = malloc(sizeof(double)*nz);
    g = malloc(sizeof(double)*nz);
    
    /* eyen, eyem */
    dptr = eyen;
    for (i = 0; i < n*n; i++)
    {
        *dptr = 0;
        dptr++;
    }
    dptr = dptr-n*n;
    for (i = 0; i < n; i++)
    {
        *dptr = 1;
        dptr = dptr+n+1;
    }

    dptr = eyem;
    for (i = 0; i < m*m; i++)
    {
        *dptr = 0;
        dptr++;
    }
    dptr = dptr-m*m;
    for (i = 0; i < m; i++)
    {
        *(dptr+i*m+i) = 1;
    }
    dptr = x; dptr1 = x0;
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }
    dptr = z;
    for (i = 0; i < T; i++)
    {
        for (j = 0; j < m; j++)
        {
            *dptr = *(U0+i*m+j);
            dptr++;
        }
        for (j = 0; j < n; j++)
        {
            *dptr = *(X0+i*n+j);
            dptr++; 
        }
    }  

    /* At, Bt */
    F77_CALL(dgemm)("t","n",&n,&n,&n,&fone,A,&n,eyen,&n,&fzero,At,&n);
    F77_CALL(dgemm)("n","t",&m,&n,&m,&fone,eyem,&m,B,&n,&fzero,Bt,&m);

    /* zmax, zmin */
    dptr1 = zmax;
    dptr2 = zmin;
    for (i = 0; i < T; i++)
    {
        for (j = 0; j < m; j++)
        {
            *dptr1 = *(umax+j);
            *dptr2 = *(umin+j);
            dptr1++; dptr2++;
        }
        for (j = 0; j < n; j++)
        {
            *dptr1 = *(xmax+j);
            *dptr2 = *(xmin+j);
            dptr1++; dptr2++;
        }
    }  

    /* zmaxp, zminp */
    for (i = 0; i < nz; i++) zminp[i] = zmin[i] + 0.01*(zmax[i]-zmin[i]);
    for (i = 0; i < nz; i++) zmaxp[i] = zmax[i] - 0.01*(zmax[i]-zmin[i]);

    /* project z */
    for (i = 0; i < nz; i++) z[i] = z[i] > zmaxp[i] ? zmaxp[i] : z[i];
    for (i = 0; i < nz; i++) z[i] = z[i] < zminp[i] ? zminp[i] : z[i];

    /* construct g = [r q ... r qf]*/
    if(r != NULL && q != NULL && qf != NULL)
    {
    	size_t idx = 0;
		for (i = 0; i < T; ++i)
		{
			for(j = 0 ; j < m; ++j)
				g[idx++] = r[j];
			for(j = 0 ; j < n; ++j)
				g[idx++] = q[j];
		}
		idx--;
		for(i = 0; i < n; ++i)
		{
			g[idx--] = qf[n-i-1];
		}
    }
    else
    {
    	memset(g,0,sizeof(double)*nz);
    }


    t1 = clock();
    fmpcsolve(A,B,At,Bt,eyen,eyem,Q,R,Qf,g,w,zmax,zmin,x,z,T,n,m,nz,niters,kappa);
    t2 = clock();
    *telapsed = (double)(t2-t1)/(CLOCKS_PER_SEC);

    dptr = z;
    for (i = 0; i < T; i++)
    {
        for (j = 0; j < m; j++)
        {
            *(U+i*m+j) = *dptr;
            dptr++;
        }
        for (j = 0; j < n; j++)
        {
            *(X+i*n+j) = *dptr;
            dptr++; 
        }
    }  

    free(At); free(Bt); free(eyen); free(eyem);
    free(z); free(x); free(zmax); free(zmin); free(g);
    return;
}

void printmat(double *A, size_t m, size_t n)
{
    double *dptr;
    size_t j, i;
    dptr = A;
    for (j = 0; j < m; j++)
    {
        for (i = 0; i < n; i++)
        {
            printf("%5.4f\t", *(dptr+m*i+j));
        }
        printf("\n");
    }
    return;
}
    
void fmpcsolve(double *A, double *B, double *At, double *Bt, double *eyen,
         double *eyem, double *Q, double *R, double *Qf, double * g, double * w, double *zmax, double *zmin,
         double *x, double *z0, size_t T, size_t n, size_t m, size_t nz, size_t niters, double kappa) 
{
    size_t maxiter = niters;
    size_t iiter, i, cont;
    double alpha = 0.01;
    double beta = 0.9;
    double tol = 0.1;
    double s;
    double resd, resp, res, newresd, newresp, newres;
    double *b, *z, *nu, *Ctnu;
    double *dnu, *dz;
    double *gf, *gp, *hp, *newgf, *newgp, *newhp;
    double *rd, *rp, *newrd, *newrp;
    double *dptr, *dptr1, *dptr2, *dptr3;
    double *newnu, *newz, *newCtnu;

    /* memory allocation */
    b = malloc(sizeof(double)*T*n);
    dnu = malloc(sizeof(double)*T*n);
    dz = malloc(sizeof(double)*nz);
    nu = malloc(sizeof(double)*T*n);
    Ctnu = malloc(sizeof(double)*nz);
    z = malloc(sizeof(double)*nz);
    gp = malloc(sizeof(double)*nz);
    hp = malloc(sizeof(double)*nz);
    gf = malloc(sizeof(double)*nz);
    rp = malloc(sizeof(double)*T*n);
    rd = malloc(sizeof(double)*nz);
    newnu = malloc(sizeof(double)*T*n);
    newz = malloc(sizeof(double)*nz);
    newCtnu = malloc(sizeof(double)*nz);
    newgf = malloc(sizeof(double)*nz);
    newhp = malloc(sizeof(double)*nz);
    newgp = malloc(sizeof(double)*nz);
    newrp = malloc(sizeof(double)*T*n);
    newrd = malloc(sizeof(double)*nz);


    /* initialize dual variable to 0
     * Shouldn't we be able to warm start nu as well?
     */
    for (i = 0; i < n*T; i++) nu[i] = 0;

    if(w == NULL)
    {
        for (i = 0; i < n*T; i++)
        	b[i] = 0;
    }
    else
    {
        for (i = 0; i < n*T; i++)
        	b[i] = w[i%n];
    }

    F77_CALL(dgemv)("n",&n,&n,&fone,A,&n,x,&ione,&fone,b,&ione);
    dptr = z; dptr1 = z0;
    for (i = 0; i < nz; i++) 
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }
    if (quiet == 0)
    {   
        printf("\n iteration \t step \t\t rd \t\t\t rp\n");
    }
    for (iiter = 0; iiter < maxiter; iiter++)
    {
        gfgphp(Q,R,Qf,g,zmax,zmin,z,T,n,m,nz,gf,gp,hp);
        rdrp(A,B,Q,R,Qf,z,nu,gf,gp,b,T,n,m,nz,kappa,rd,rp,Ctnu);
        resdresp(rd,rp,T,n,nz,&resd,&resp,&res);

//        printf("gf_mpc = [\n");
//        printmat(gf,nz,1);
//        printf("];\n");

//        printf("gp_mpc = [\n");
//        printmat(gp,nz,1);
//        printf("];\n");

//        printf("hp_mpc = [\n");
//        printmat(hp, nz, 1);
//        printf("];\n");

        if (res < tol) break;
        dnudz(A,B,At,Bt,eyen,eyem,Q,R,Qf,hp,rd,rp,T,n,m,nz,kappa,dnu,dz);

//        printf("kappa_mpc = %g;\n",kappa);
//        printf("A_mpc = [\n");
//        printmat(A,n,n);
//        printf("];\n");
//
//        printf("At_mpc = [\n");
//        printmat(At,n,n);
//        printf("];\n");
//
//        printf("B_mpc = [\n");
//        printmat(B,n,m);
//        printf("];\n");
//
//        printf("Bt_mpc = [\n");
//        printmat(Bt,m,n);
//        printf("];\n");
//
//        printf("eyen_mpc = [\n");
//        printmat(eyen,n,n);
//        printf("];\n");
//
//        printf("eyem_mpc = [\n");
//        printmat(eyem,m,m);
//        printf("];\n");
//
//        printf("Q_mpc = [\n");
//        printmat(Q,n,n);
//        printf("];\n");
//
//        printf("Qf_mpc = [\n");
//        printmat(Qf,n,n);
//        printf("];\n");
//
//        printf("R_mpc = [\n");
//        printmat(R,m,m);
//        printf("];\n");
//
//
//
//
//        printf("rd_mpc = [\n");
//        printmat(rd, nz, 1);
//        printf("];\n");
//
//        printf("rp_mpc = [\n");
//        printmat(rp, T*n, 1);
//        printf("];\n");
//
//
//                printf("dz_mpc = [\n");
//                printmat(dz, nz, 1);
//                printf("];\n");
//
//                printf("dnu_mpc = [\n");
//                printmat(dnu, 500, 1);
//                printf("];\n");



        s = 1; 
        /* feasibility search */
        while (1)
        {
            cont = 0;
            dptr = z; dptr1 = dz; dptr2 = zmax; dptr3 = zmin;
            for (i = 0; i < nz && cont == 0; i++)
            {
                if (*dptr+s*(*dptr1) >= *dptr2) cont = 1;
                if (*dptr+s*(*dptr1) <= *dptr3) cont = 1;
                dptr++; dptr1++; dptr2++; dptr3++;
            }
            if (cont == 1)
            {
            	if(s<0.01)
            	{
//            		printf("s = %g for z[%i]=%g is infeasible !\n",s,i-1,dz[i-1]);

            		s = 0;
            		break;
            	}
                s = beta*s;
                continue;
            }
            else
                break;
        }

        dptr = newnu; dptr1 = nu; dptr2 = dnu;
        for (i = 0; i < T*n; i++)
        {
            *dptr = *dptr1+s*(*dptr2);
            dptr++; dptr1++; dptr2++;
        }
        dptr = newz; dptr1 = z; dptr2 = dz;
        for (i = 0; i < nz; i++)
        {
            *dptr = *dptr1+s*(*dptr2);
            dptr++; dptr1++; dptr2++;
        }

        /* insert backtracking line search */
        while (1)
        {
            gfgphp(Q,R,Qf,g,zmax,zmin,newz,T,n,m,nz,newgf,newgp,newhp);
            rdrp(A,B,Q,R,Qf,newz,newnu,newgf,newgp,b,T,n,m,nz,kappa,newrd,newrp,newCtnu);
            resdresp(newrd,newrp,T,n,nz,&newresd,&newresp,&newres);
            if (newres <= (1-alpha*s)*res || (s == 0))
            	break;
            else if(s<1e-3)
            {
            	s = 0;
//            	printf("s is too small!!\n");
            }
            else
            {
//            	printf("s = %g, newres = %g , oldres = %g\n",s, newres, res);
            }
            s = beta*s;
            dptr = newnu; dptr1 = nu; dptr2 = dnu;
            for (i = 0; i < T*n; i++)
            {
                *dptr = *dptr1+s*(*dptr2);
                dptr++; dptr1++; dptr2++;
            }
            dptr = newz; dptr1 = z; dptr2 = dz;
            for (i = 0; i < nz; i++)
            {
                *dptr = *dptr1+s*(*dptr2);
                dptr++; dptr1++; dptr2++;
            }
        }
        
        dptr = nu; dptr1 = newnu; 
        for (i = 0; i < T*n; i++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = z; dptr1 = newz;
        for (i = 0; i < nz; i++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        if (quiet == 0)
        {
            printf("    %d \t\t %5.4f \t %0.5e \t\t %0.5e\n",iiter,s,newresd,newresp);
        }
    }
    dptr = z0; dptr1 = z;
    for (i = 0; i < nz; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }

    free(b); free(dnu); free(dz); free(nu); free(Ctnu);
    free(z); free(gp); free(hp); free(gf); free(rp); free(rd);
    free(newnu); free(newz); free(newCtnu); free(newgf); free(newhp);
    free(newgp); free(newrp); free(newrd);
    return;
}

/* computes the search directions dz and dnu */
void dnudz(double *A, double *B, double *At, double *Bt, double *eyen,
        double *eyem, double *Q, double *R, double *Qf, double *hp, double *rd, 
        double *rp, size_t T, size_t n, size_t m, size_t nz, double kappa, double *dnu, double *dz)
{
    size_t i,j,info,nT;
    double *dptr, *dptr1, *dptr2, *dptr3, *temp, *tempmatn, *tempmatm;
    double *PhiQ, *PhiR, *Yd, *Yud, *Ld, *Lld, *Ctdnu, *gam, *v, *be, *rdmCtdnu;
    double *PhiinvQAt, *PhiinvRBt, *PhiinvQeye, *PhiinvReye, *CPhiinvrd;
    nT = n*T;

    /* allocate memory */
    PhiQ = malloc(sizeof(double)*n*n*T);
    PhiR = malloc(sizeof(double)*m*m*T);
    PhiinvQAt = malloc(sizeof(double)*n*n*T);
    PhiinvRBt = malloc(sizeof(double)*m*n*T);
    PhiinvQeye = malloc(sizeof(double)*n*n*T);
    PhiinvReye = malloc(sizeof(double)*m*m*T);
    CPhiinvrd = malloc(sizeof(double)*n*T);
    Yd = malloc(sizeof(double)*n*n*T);
    Yud = malloc(sizeof(double)*n*n*(T-1));
    Ld = malloc(sizeof(double)*n*n*T);
    Lld = malloc(sizeof(double)*n*n*(T-1));
    gam = malloc(sizeof(double)*n*T);
    v = malloc(sizeof(double)*n*T);
    be = malloc(sizeof(double)*n*T);
    temp = malloc(sizeof(double)*n);
    tempmatn = malloc(sizeof(double)*n*n);
    tempmatm = malloc(sizeof(double)*m*m);
    Ctdnu = malloc(sizeof(double)*nz);
    rdmCtdnu = malloc(sizeof(double)*nz);

    /* form PhiQ and PhiR */
    for (i = 0; i < T-1; i++)
    {
        dptr = PhiQ+n*n*i; dptr1 = Q;
        for (j = 0; j < n*n; j++)
        {
            *dptr = 2*(*dptr1);
            dptr++; dptr1++;
        }
        dptr = PhiQ+n*n*i; dptr1 = hp+m*(i+1)+n*i;
        for (j = 0; j < n; j++)
        {
            *dptr = *dptr+kappa*(*dptr1);
            dptr = dptr+n+1; dptr1++;
        }
        dptr = PhiR+m*m*i; dptr1 = R;
        for (j = 0; j < m*m; j++)
        {
            *dptr = 2*(*dptr1);
            dptr++; dptr1++;
        }
        dptr = PhiR+m*m*i; dptr1 = hp+i*(n+m);
        for (j = 0; j < m; j++)
        {
            *dptr = *dptr+kappa*(*dptr1);
            dptr = dptr+m+1; dptr1++;
        }
    }
    
    dptr = PhiR+m*m*(T-1); dptr1 = R;
    for (j = 0; j < m*m; j++)
    {
        *dptr = 2*(*dptr1);
        dptr++; dptr1++;
    }
    dptr = PhiR+m*m*(T-1); dptr1 = hp+(T-1)*(n+m);
    for (j = 0; j < m; j++)
    {
        *dptr = *dptr+kappa*(*dptr1);
        dptr = dptr+m+1; dptr1++;
    }
    dptr = PhiQ+n*n*(T-1); dptr1 = Qf;
    for (j = 0; j < n*n; j++)
    {
        *dptr = 2*(*dptr1);
        dptr++; dptr1++;
    }
    dptr = PhiQ+n*n*(T-1); dptr1 = hp+m*T+n*(T-1);
    for (j = 0; j < n; j++)
    {
        *dptr = *dptr+kappa*(*dptr1);
        dptr = dptr+n+1; dptr1++;
    }

    /* compute PhiinvQAt, PhiinvRBt, PhiinvQeye, PhiinvReye */
    for (i = 0; i < T; i++)
    {
        dptr = PhiinvQAt+n*n*i; dptr1 = At;
        for (j = 0; j < n*n; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = dptr-n*n; dptr1 = PhiQ+n*n*i;
        F77_CALL(dposv)("l",&n,&n,dptr1,&n,dptr,&n,&info);
        dptr = PhiinvQeye+n*n*i; dptr1 = eyen;
        for (j = 0; j < n*n; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = dptr-n*n; dptr1 = PhiQ+n*n*i;
        F77_CALL(dtrtrs)("l","n","n",&n,&n,dptr1,&n,dptr,&n,&info);
        F77_CALL(dtrtrs)("l","t","n",&n,&n,dptr1,&n,dptr,&n,&info);
    }
    for (i = 0; i < T; i++)
    {
        dptr = PhiinvRBt+m*n*i; dptr1 = Bt;
        for (j = 0; j < n*m; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = dptr-m*n; dptr1 = PhiR+m*m*i;
        F77_CALL(dposv)("l",&m,&n,dptr1,&m,dptr,&m,&info);
        dptr = PhiinvReye+m*m*i; dptr1 = eyem;
        for (j = 0; j < m*m; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = dptr-m*m; dptr1 = PhiR+m*m*i;
        F77_CALL(dtrtrs)("l","n","n",&m,&m,dptr1,&m,dptr,&m,&info);
        F77_CALL(dtrtrs)("l","t","n",&m,&m,dptr1,&m,dptr,&m,&info);
    }
    
    /* form Yd and Yud */
    dptr = Yud; dptr1 = PhiinvQAt; 
    for (i = 0; i < n*n*(T-1); i++)
    {
        *dptr = -(*dptr1);
        dptr++; dptr1++;
    }
    dptr2 = Yd; dptr3 = PhiinvQeye;
    for (i = 0; i < n*n; i++)
    {
        *dptr2 = *dptr3;
        dptr2++; dptr3++;
    }
    dptr2 = dptr2-n*n;
    F77_CALL(dgemm)("n","n",&n,&n,&m,&fone,B,&n,PhiinvRBt,&m,&fone,dptr2,&n);
    for (i = 1; i < T; i++)
    {
        dptr = Yd+n*n*i; dptr1 = PhiinvQeye+n*n*i;
        for (j = 0; j < n*n; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr1 = PhiinvRBt+m*n*i; dptr = dptr-n*n;
        F77_CALL(dgemm)("n","n",&n,&n,&m,&fone,B,&n,dptr1,&m,&fone,dptr,&n); 
        dptr1 = PhiinvQAt+n*n*(i-1);
        F77_CALL(dgemm)("n","n",&n,&n,&n,&fone,A,&n,dptr1,&n,&fone,dptr,&n);
    }

    /* compute Lii */
    dptr = Ld; dptr1 = Yd; 
    for (i = 0; i < n*n; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++; 
    }
    dptr = dptr-n*n; 
    F77_CALL(dposv)("l",&n,&ione,dptr,&n,temp,&n,&info);
    for (i = 1; i < T; i++)
    {
        dptr = Ld+n*n*(i-1); dptr1 = Yud+n*n*(i-1); dptr2 = Lld+n*n*(i-1);
        for (j = 0; j < n*n; j++)
        {
            *dptr2 = *dptr1;
            dptr2++; dptr1++;
        }
        dptr2 = dptr2-n*n;
        F77_CALL(dtrtrs)("l","n","n",&n,&n,dptr,&n,dptr2,&n,&info);
        dptr1 = tempmatn;
        for (j = 0; j < n*n; j++)
        {
            *dptr1 = *dptr2;
            dptr1++; dptr2++;
        }
        dptr1 = dptr1-n*n; dptr2 = dptr2-n*n;
        F77_CALL(dgemm)("t","n",&n,&n,&n,&fone,dptr1,&n,eyen,&n,&fzero,dptr2,&n);
        dptr = Ld+n*n*i; dptr1 = Yd+n*n*i;
        for (j = 0; j < n*n; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = dptr-n*n;
        F77_CALL(dgemm)("n","t",&n,&n,&n,&fmone,dptr2,&n,dptr2,&n,&fone,dptr,&n);
        F77_CALL(dposv)("l",&n,&ione,dptr,&n,temp,&n,&info);
    }

    /* compute CPhiinvrd */
    dptr = CPhiinvrd; dptr1 = rd+m;
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }
    dptr = dptr-n;
    F77_CALL(dtrsv)("l","n","n",&n,PhiQ,&n,dptr,&ione);
    F77_CALL(dtrsv)("l","t","n",&n,PhiQ,&n,dptr,&ione);
    dptr2 = temp; dptr1 = rd;
    for (i = 0; i < m; i++)
    {
        *dptr2 = *dptr1;
        dptr2++; dptr1++;
    }
    dptr2 = dptr2-m;
    F77_CALL(dtrsv)("l","n","n",&m,PhiR,&m,dptr2,&ione);
    F77_CALL(dtrsv)("l","t","n",&m,PhiR,&m,dptr2,&ione);
    F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,temp,&ione,&fone,dptr,&ione);
    
    for (i = 1; i < T; i++)
    {
        dptr = CPhiinvrd+n*i; dptr1 = rd+m+i*(n+m);
        for (j = 0; j < n; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = dptr-n; dptr3 = PhiQ+n*n*i;
        F77_CALL(dtrsv)("l","n","n",&n,dptr3,&n,dptr,&ione);
        F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr,&ione);
        dptr2 = temp; dptr1 = rd+i*(m+n);
        for (j = 0; j < m; j++)
        {
            *dptr2 = *dptr1;
            dptr2++; dptr1++;
        }
        dptr3 = PhiR+m*m*i; dptr2 = dptr2-m;
        F77_CALL(dtrsv)("l","n","n",&m,dptr3,&m,dptr2,&ione);
        F77_CALL(dtrsv)("l","t","n",&m,dptr3,&m,dptr2,&ione);
        F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,temp,&ione,&fone,dptr,&ione);
        dptr2 = temp; dptr1 = rd+(i-1)*(n+m)+m;
        for (j = 0; j < n; j++)
        {
            *dptr2 = *dptr1;
            dptr2++; dptr1++;
        }
        dptr3 = PhiQ+n*n*(i-1); dptr2 = dptr2-n;
        F77_CALL(dtrsv)("l","n","n",&n,dptr3,&n,dptr2,&ione);
        F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr2,&ione);
        F77_CALL(dgemv)("n",&n,&n,&fmone,A,&n,temp,&ione,&fone,dptr,&ione);
    }

    /* form be */
    dptr = be; dptr1 = rp; dptr2 = CPhiinvrd;
    for (i = 0; i < n*T; i++)
    {
        *dptr = (*dptr2)-(*dptr1);
        dptr++; dptr1++; dptr2++;
    }

    /* solve for dnu */
    dptr = v; dptr1 = be;
    for (i = 0; i < n; i++)
    {
        *dptr = -(*dptr1);
        dptr++; dptr1++;
    }
    dptr = dptr-n;
    F77_CALL(dtrsv)("l","n","n",&n,Ld,&n,dptr,&ione);
    for (i = 1; i < T; i++)
    {
        dptr = v+i*n; dptr1 = v+(i-1)*n; dptr2 = be+i*n; 
        for (j = 0; j < n; j++)
        {
            *dptr = *dptr2;
            dptr++; dptr2++;
        }
        dptr = dptr-n; dptr3 = Lld+n*n*(i-1);
        F77_CALL(dgemv)("n",&n,&n,&fmone,dptr3,&n,dptr1,&ione,&fmone,dptr,&ione);
        dptr3 = Ld+n*n*i;
        F77_CALL(dtrsv)("l","n","n",&n,dptr3,&n,dptr,&ione);
    }
    dptr = dnu+n*(T-1); dptr1 = v+n*(T-1);
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }
    dptr = dptr-n; dptr3 = Ld+n*n*(T-1);
    F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr,&ione);
    for (i = T-1; i > 0; i--)
    {
        dptr = dnu+n*(i-1); dptr1 = dnu+n*i; dptr2 = v+n*(i-1); 
        for (j = 0; j < n; j++)
        {
            *dptr = *dptr2;
            dptr++; dptr2++;
        }
        dptr = dptr-n; dptr3 = Lld+n*n*(i-1);
        F77_CALL(dgemv)("t",&n,&n,&fmone,dptr3,&n,dptr1,&ione,&fone,dptr,&ione);
        dptr3 = Ld+n*n*(i-1);
        F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr,&ione);
    }

    /* form Ctdnu */
    for (i = 0; i < T-1; i++)
    {
        dptr = Ctdnu+i*(n+m); dptr1 = dnu+i*n;
        F77_CALL(dgemv)("n",&m,&n,&fmone,Bt,&m,dptr1,&ione,&fzero,dptr,&ione);
        dptr = Ctdnu+i*(n+m)+m; dptr2 = dnu+(i+1)*n;
        for (j = 0; j < n; j++)
        {
            *dptr = *dptr1;
            dptr++; dptr1++;
        }
        dptr = dptr-n;
        F77_CALL(dgemv)("n",&n,&n,&fmone,At,&n,dptr2,&ione,&fone,dptr,&ione);
    }
    
    dptr = Ctdnu+(T-1)*(n+m); dptr1 = dnu+(T-1)*n;
    F77_CALL(dgemv)("n",&m,&n,&fmone,Bt,&m,dptr1,&ione,&fzero,dptr,&ione);
    dptr = dptr+m; 
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }
    dptr = rdmCtdnu; dptr1 = Ctdnu; dptr2 = rd;
    for (i = 0; i < nz; i++)
    {
        *dptr = -(*dptr1)-(*dptr2);
        dptr++; dptr1++; dptr2++;
    }

    /* solve for dz */
    for (i = 0; i < T; i++)
    {
        dptr = dz+(i+1)*m+i*n; dptr1 = rdmCtdnu+(i+1)*m+i*n;
        dptr2 = PhiinvQeye+n*n*i;
        F77_CALL(dgemv)("n",&n,&n,&fone,dptr2,&n,dptr1,&ione,&fzero,dptr,&ione);
    }
    for (i = 0; i < T; i++)
    {
        dptr = dz+i*(m+n); dptr1 = rdmCtdnu+i*(m+n);
        dptr2 = PhiinvReye+m*m*i;
        F77_CALL(dgemv)("n",&m,&m,&fone,dptr2,&m,dptr1,&ione,&fzero,dptr,&ione);
    }
    free(PhiQ); free(PhiR); free(PhiinvQAt); free(PhiinvRBt); free(PhiinvQeye);
    free(PhiinvReye); free(CPhiinvrd); free(Yd); free(Yud); free(Ld); free(Lld);
    free(gam); free(v); free(be); free(temp); free(tempmatn); free(tempmatm); 
    free(Ctdnu); free(rdmCtdnu);
    return;
}

/* computes rd and rp */
void rdrp(double *A, double *B, double *Q, double *R, double *Qf, double *z, double *nu, 
        double *gf, double *gp, double *b, size_t T, size_t n, size_t m, size_t nz, 
        double kappa, double *rd, double *rp, double *Ctnu)
{
    size_t i, j;
    double *Cz;
    double *dptr, *dptr1, *dptr2;

    Cz = malloc(sizeof(double)*T*n);
    
    /* compute Cz */
    dptr = Cz; dptr1 = z+m;
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }
    F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,z,&ione,&fone,Cz,&ione);
    for (i = 2; i <= T; i++)
    {
        dptr = Cz+(i-1)*n; dptr1 = z+m+(i-2)*(n+m); 
        dptr2 = z+m+(i-1)*(m+n);
        for (j = 0; j < n; j++)
        {
            *dptr = *dptr2;
            dptr++; dptr2++;
        }
        dptr = dptr-n; 
        F77_CALL(dgemv)("n",&n,&n,&fmone,A,&n,dptr1,&ione,&fone,dptr,&ione);
        dptr1 = dptr1+n;
        F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,dptr1,&ione,&fone,dptr,&ione);
    }
    /*
    dptr = Cz+(T-1)*n; dptr1 = z+m+(T-2)*(n+m);
    F77_CALL(dgemv)("n",&n,&n,&fmone,A,&n,dptr1,&ione,&fzero,dptr,&ione);
    dptr1 = dptr1+n;
    F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,dptr1,&ione,&fone,dptr,&ione);
    dptr1 = z+nz-n;
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr+*dptr1;
        dptr++; dptr1++;
    }
    */

    /* compute Ctnu */
    dptr = Ctnu; dptr1 = Ctnu+m; dptr2 = nu;
    for (i = 1; i <= T-1; i++)
    {
        F77_CALL(dgemv)("t",&n,&m,&fmone,B,&n,dptr2,&ione,&fzero,dptr,&ione);
        dptr = dptr+n+m;
        for (j = 0; j < n; j++)
        {
            *dptr1 = *dptr2;
            dptr1++; dptr2++;
        }
        dptr1 = Ctnu+m+(i-1)*(n+m);
        F77_CALL(dgemv)("t",&n,&n,&fmone,A,&n,dptr2,&ione,&fone,dptr1,&ione);
        dptr1 = dptr1+n+m;
    }
    F77_CALL(dgemv)("t",&n,&m,&fmone,B,&n,dptr2,&ione,&fzero,dptr,&ione);
    dptr = Ctnu+nz-n; dptr1 = nu+(T-1)*n;
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr1;
        dptr++; dptr1++;
    }

    dptr = rp; dptr1 = Cz; dptr2 = b;
    for (i = 0; i < n*T; i++)
    {
        *dptr = *dptr1-*dptr2;
        dptr++; dptr1++; dptr2++;
    }
    dptr = rd; dptr1 = Ctnu; dptr2 = gf;
    for (i = 0; i < nz; i++)
    {
        *dptr = *dptr1+*dptr2;
        dptr++; dptr1++; dptr2++;
    }
    F77_CALL(daxpy)(&nz,&kappa,gp,&ione,rd,&ione);
    free(Cz);
    return;
}

/* computes gf, gp and hp */
void gfgphp(double *Q, double *R, double *Qf, double * g, double *zmax, double *zmin, double *z,
        size_t T, size_t n, size_t m, size_t nz, double *gf, double *gp, double *hp)
{
    size_t i;
    double *dptr, *dptr1, *dptr2;
    double *gp1, *gp2;

    gp1 = malloc(sizeof(double)*nz);
    gp2 = malloc(sizeof(double)*nz);

    dptr = gp1; dptr1 = zmax; dptr2 = z;
    for (i = 0; i < nz; i++)
    {
        *dptr = 1.0/(*dptr1-*dptr2);
        dptr++; dptr1++; dptr2++;
    }
    dptr = gp2; dptr1 = zmin; dptr2 = z;
    for (i = 0; i < nz; i++)
    {
        *dptr = 1.0/(*dptr2-*dptr1);
        dptr++; dptr1++; dptr2++;
    }
    dptr = hp; dptr1 = gp1; dptr2 = gp2;
    for (i = 0; i < nz; i++)
    {
        *dptr = (*dptr1)*(*dptr1) + (*dptr2)*(*dptr2);
        dptr++; dptr1++; dptr2++;
    }
    dptr = gp; dptr1 = gp1; dptr2 = gp2;
    for (i = 0; i < nz; i++)
    {
        *dptr = *dptr1-*dptr2;
        dptr++; dptr1++; dptr2++;
    }
    
    for(i = 0 ; i < nz; ++i)
    {
    	gf[i] = g[i];
    }

    dptr = gf; dptr1 = z; 
    for (i = 0; i < T-1; i++)
    {
        F77_CALL(dgemv)("n",&m,&m,&ftwo,R,&m,dptr1,&ione,&fone,dptr,&ione);
        dptr = dptr+m; dptr1 = dptr1+m;
        F77_CALL(dgemv)("n",&n,&n,&ftwo,Q,&n,dptr1,&ione,&fone,dptr,&ione);
        dptr = dptr+n; dptr1 = dptr1+n;
    }
    F77_CALL(dgemv)("n",&m,&m,&ftwo,R,&m,dptr1,&ione,&fone,dptr,&ione);
    dptr = dptr+m; dptr1 = dptr1+m;
    F77_CALL(dgemv)("n",&n,&n,&ftwo,Qf,&n,dptr1,&ione,&fone,dptr,&ione);

    free(gp1); free(gp2);
    return;
}

/* computes resd, resp, and res */
void resdresp(double *rd, double *rp, size_t T, size_t n, size_t nz, double *resd, 
        double *resp, double *res)
{
    size_t nnu = T*n;
    *resp = F77_CALL(dnrm2)(&nnu,rp,&ione);
    *resd = F77_CALL(dnrm2)(&nz,rd,&ione);
    *res = sqrt((*resp)*(*resp)+(*resd)*(*resd));
    return;
}
