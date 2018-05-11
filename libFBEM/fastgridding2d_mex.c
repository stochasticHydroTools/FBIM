#include "mex.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>


void fastgridding2d(double x, double y, double alpha, double h, mwSize P, 
                    double *a, double *ix, double *iy)
{
    int i0,j0,lstart, l2, i,j,l;
    double dx, dy, e1x, e1y, e2x, e2y, e0, Vy;
    double *E0, *Ex, *Ey;
    E0 =  (double *)malloc( sizeof(double) * P );   
    Ex =  (double *)malloc( sizeof(double) * P );
    Ey =  (double *)malloc( sizeof(double) * P );
    
    i0 = (int) (x/h); 
    j0 = (int) (y/h);

    dx = x - i0*h;
    dy = y - j0*h;

    e0 = exp( -alpha*h*h );    

    e1x = exp( 2*alpha*dx*h );
    e1y = exp( 2*alpha*dy*h );

    e2x = exp( -alpha*dx*dx );
    e2y = exp( -alpha*dy*dy );

    lstart = -(P-1)/2; 
    for (l = 0; l < P ; l++) 
    {
        l2 = (lstart+l)*(lstart+l);
        E0[l] = pow(e0, l2); 
        Ex[l] = pow(e1x, lstart + l) * e2x ;
        Ey[l] = pow(e1y, lstart + l) * e2y ; 
    }
 
    for ( i =0; i  < P; i++ )
    {   
        Vy = E0[i] * Ey[i];
        for ( j=0; j<P; j++ ){
            a[j+i*P] = Vy * E0[j] * Ex[j];
        }

    }
    
    for (i = 0; i < P ; i++)
    {
        ix[i] = i0 + i + lstart ;
        iy[i] = j0 + i + lstart ;
    }
 
    free(E0); free(Ex); free(Ey);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *outMatrix, *ix, *iy;
    double x, y, alpha, h ; 
    mwSize P; 
    
    x = mxGetScalar(prhs[0]);
    y = mxGetScalar(prhs[1]);
    alpha = mxGetScalar(prhs[2]);
    h = mxGetScalar(prhs[3]);
    
    P = (mwSize) mxGetScalar(prhs[4]);
    
    plhs[0] = mxCreateDoubleMatrix(P, P, mxREAL);
    outMatrix = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(1, P, mxREAL);
    ix = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(1, P, mxREAL);
    iy = mxGetPr(plhs[2]);
    
    
    fastgridding2d(x, y, alpha,  h, P, outMatrix, ix, iy);
    
}