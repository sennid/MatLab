#include "mex.h"
#include "math.h"

void funcCheckInputParams(int nrhs,const mxArray **ar);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    funcCheckInputParams(nrhs,prhs);

    double *A, *B, *C;
    
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    int r = mxGetM(prhs[0]);
    int c = mxGetN(prhs[0]);

    const mwSize dim[2] = {c, r};
    mxArray *D;
    D = mxCreateNumericArray(2, dim, mxDOUBLE_CLASS, mxREAL);
    double *ptr_d;
    ptr_d = mxGetPr(D);
    int i, j;
    for(i = 0; i < c; i++)
    {
        for(j = 0; j < r; j++)
        {
            ptr_d[i*r+j] = B[i*r+j]*B[i*r+j] - 4*A[i*r+j]*C[i*r+j];
        }
    }


     mxArray *S;
    S = mxCreateNumericArray(2, dim, mxDOUBLE_CLASS, mxCOMPLEX);
    double *ptr_ans;
    ptr_ans = mxGetPr(S);
    for(i = 0; i < c; i++)
    {
        for(j = 0; j < r; j++)
        {   
             if (fabs(A[i*r+j]) > pow(10, -6))
            ptr_ans[i*r+j] = (-B[i*r+j] + sqrt(ptr_d[i*r+j]))/(2*A[i*r+j]);
            else if (fabs(B[i*r+j]) > pow(10, -6))
                ptr_ans[i*r+j] = -C[i*r+j] / B[i*r+j];
            else if (fabs(C[i*r+j]) > pow(10, -6))
                ptr_ans[i*r+j] = NAN;
            else {
                mexWarnMsgTxt("Hyper continuum of roots!!");
                ptr_ans[i*r+j] = 0;
            }
         }
    }

    plhs[0] = S;

     mxArray *S1;
    S1 = mxCreateNumericArray(2, dim, mxDOUBLE_CLASS, mxCOMPLEX);
    double *ptr_ans1;
    ptr_ans1 = mxGetPr(S1);
    for(i = 0; i < c; i++)
    {
        for(j = 0; j < r; j++)
        {
            if (fabs(A[i*r+j]) > pow(10, -6))
            ptr_ans1[i*r+j] = (-B[i*r+j] - sqrt(ptr_d[i*r+j]))/(2*A[i*r+j]);
            else if (fabs(B[i*r+j]) > pow(10, -6))
                ptr_ans1[i*r+j] = -C[i*r+j] / B[i*r+j];
            else if (fabs(C[i*r+j]) > pow(10, -6))
                ptr_ans1[i*r+j] = NAN;
            else
                ptr_ans1[i*r+j] = 0;
        }
    }

    plhs[1] = S1;

    plhs[2] = D;
}   

void funcCheckInputParams(int nrhs,const mxArray **ar)
{
    if(nrhs != 3)
        mexErrMsgTxt("Wrong amount of input data!");

    int rows0, cols0, rows1, cols1, rows2, cols2;
    rows0 = mxGetM(ar[0]);
    cols0 = mxGetN(ar[0]);
    rows1 = mxGetM(ar[1]);
    cols1 = mxGetN(ar[1]);
    rows2 = mxGetM(ar[2]);
    cols2 = mxGetN(ar[2]);
    if((rows0 != rows1) || (cols0 != cols1) || (rows1 != rows2) || (cols1 != cols2)) 
        mexErrMsgTxt("matrix has incorrect dimension!");

}