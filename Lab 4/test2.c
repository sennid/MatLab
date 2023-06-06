#include "mex.h"
#include "math.h"

void funcCheckInputParams(int nrhs,const mxArray **ar);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    funcCheckInputParams(nrhs,prhs);

    double *ptrA, *ptrM, *ptrN;
    unsigned int M, N;
    ptrA = mxGetPr(prhs[0]);
    ptrM = mxGetPr(prhs[1]);
    ptrN = mxGetPr(prhs[2]);

    M = (unsigned int)ptrM[0];
    N = (unsigned int)ptrN[0];
    
    mexPrintf("A[0, 0] = %f\n", ptrA[0]);


    //output
    const mwSize dim[2] = {M, N};
    mxArray *ansMat;
    ansMat = mxCreateNumericArray(2, dim, mxDOUBLE_CLASS, mxREAL);

    double *ptr_ans;
    ptr_ans = mxGetPr(ansMat);
    int i, j;
    for(i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            ptr_ans[i*N + j] = -3;
        }
    }
    plhs[0] = ansMat;
}

void funcCheckInputParams(int nrhs,const mxArray **ar)
{
    if(nrhs != 3)
        mexErrMsgTxt("Wrong amount of input data!");

    int rows, cols;
    rows = mxGetM(ar[1]);
    cols = mxGetN(ar[1]);
    if((rows != 1) || (cols != 1))
        mexErrMsgTxt("second input parameter is not scalar!");

    rows = mxGetM(ar[2]);
    cols = mxGetN(ar[2]);
    if((rows != 1) || (cols != 1))
        mexErrMsgTxt("third input parameter is not scalar!");
}