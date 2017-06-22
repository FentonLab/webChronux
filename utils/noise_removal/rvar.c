/* rvar.c */
#include "mex.h"
/* Input Arguments */
#define	X_IN	prhs[0] /*input vector*/
#define	NWIN	prhs[1] /*length of window length*/
#define NLEN	prhs[2] /*input vector length*/
/* Output Arguments */
#define	Y_OUT	plhs[0] /*output vector*/

/*Get variance from xp[st:ed]*/
static double varWin(double xp[], int st, int ed)
{
    int i;
    double mean, len, vari, x;

	len = ed-st+1; /*length of the vector part*/

	/*get mean first*/
	mean = 0;
    for (i=st; i <= ed; i++)
    {
		mean = mean + xp[i];
	}
	mean = mean / len;
	
	/*get variance*/
	vari = 0;
    for (i=st; i <= ed; i++)
    {
        x = xp[i];
		vari = vari + (x - mean)*(x - mean);
	}
	vari = vari / (len-1);
    return vari;
}

/*taking a sliding variance of x[], using window length (half)*/
/*length is the total vector length*/
static void variance(double	yp[], double x[], int winLenHalf, int length)
{
    double	xnew, xmax;
    int i, st, ed, winST, winED;

    /*starting index*/
	st = winLenHalf;
	/*ending index*/
	ed = length - winLenHalf - 1;
	
	for (i=st; i <= ed; i++)
	{
		winST = i - winLenHalf;
		winED = i + winLenHalf;
		
		yp[i] = varWin(x, winST, winED);
	}
}

/*nrhs = input, nlhs = output*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    double *yp, *x,*y, *pnwin, *pnlen;
    int m,n, nwin, nlen;

    if (nrhs != 3)
		mexErrMsgTxt("Three input arguments required.");
    if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
    if (!mxIsDouble(X_IN))
		mexErrMsgTxt("RMAX requires X to be a double vector.");
    
	m = mxGetM(X_IN); /*get rows*/
    n = mxGetN(X_IN); /*get columns*/
    
	/*get NWIN*/
	pnwin = mxGetPr(NWIN);
	nwin = (int)*pnwin;

	/*get NLEN*/
	pnlen = mxGetPr(NLEN);
	nlen = (int)*pnlen;

    /*get X_in matrix*/
	x = mxGetPr(X_IN);

    /* Create a matrix for the return argument */
    Y_OUT = mxCreateDoubleMatrix(1, n, mxREAL);
    yp = mxGetPr(Y_OUT); 
    variance(yp, x, nwin, nlen); 
}