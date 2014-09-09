//
// Copyright 1999 Molecular Structure Corporation
//                9009 New Trails Drive
//                The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved
//
// dtsvd.h      Initial author: Thom Hendrixson        September 1999
/*
 *
 * Copyright (C) 2014 Rigaku Americas Corporation
 *                    9009 New Trails Drive
 *                    The Woodlands, TX, USA  77381
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice(s), this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice(s), this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of the Rigaku Americas Corporation nor the 
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL RIGAKU AMERICAS CORPORATION BE LIABLE 
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA OR PROFITS; OR BUSINESS INTERUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 */


/****************************************************************************
 *                          Function Prototypes                             *
 ****************************************************************************/

#ifndef DTSVD_DEFS
#define DTSVD_DEFS

#include "dtrekdefs.h"

#ifndef LINUX
#ifndef OSF1
double   DTREK_EXPORT copysignf(double a, double b);
#endif
#endif

#include "dtarray.h"

DTREK_EXPORT int     gaussj(double **a, int n, double **b, int m);
DTREK_EXPORT double   pythag(double a, double b);
DTREK_EXPORT void    svdfit(double *x, double *y, double *sig, int ndata, double *a, int na,
               double **u, double **v, double *w, double *chisq,
               void (*function)(double, double *, int));
DTREK_EXPORT void    svbksb(double **u, double *w, double **v, int m, int n, double *b,
               double *x);
DTREK_EXPORT int     svdcmp(double **a, int m, int n, double *w, double **v);

DTREK_EXPORT double **matrix(int nrl, int nrh, int ncl, int nch);
DTREK_EXPORT int    *ivector(int low_index, int high_index);
DTREK_EXPORT double  *vectorDT(int low_index, int high_index);
DTREK_EXPORT void    free_ivector(int *v, int nl, int nh);
DTREK_EXPORT void    free_vector(double *v, int nl, int nh);
DTREK_EXPORT void    free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
DTREK_EXPORT bool    bCheck_matrix(double **m, int nrl, int nrh, int ncl, int nch);

DTREK_EXPORT int nInvMatND_svd(int nN,float* fANxN,float* fANxNInv);
DTREK_EXPORT int nSolveAxB_svd(int nN,float* fANxN,float* fXN,float* fBN);
DTREK_EXPORT int nSolveAxB_svd(int nN,double* fANxN,double* fXN,double* fBN);
DTREK_EXPORT int nInvMatND_svd(int nN,double* fANxN,double* fANxNInv,double* fANxNEigenVecs = NULL,double* fANEigenValues = NULL);


DTREK_EXPORT int nSolveAxB_gauss(int nN,float* fANxN,float* fXN,float* fBN);
DTREK_EXPORT float fScaleMatMN(int nN,int nM,float* fANxN,bool bApply);
DTREK_EXPORT double fScaleMatMN(int nN,int nM,double* fANxN,bool bApply);

// Routine only works for normal matrices (i.e., those that are 
// unit. diagonalizable)
// Correct eigenvalues are returned for non-normal matrices, 
// but the eigenvectors will be bogus.

DTREK_EXPORT void   vEigen(int nN,double* anxnfMat,double* anfEigenVals,double* anxnfEigenVecs);

// Solve a least squares system.  The ppfFuncs can point to one of the following "standard" functions.
// or can point to an array of user defined values.


const int    ms_nNumDummyFuncs = 20;

static double* ms_pfSqrtX = (double*) 0x1;
static double* ms_pfConst = (double*) 0x2;
static double* ms_pfX = (double*) 0x3;
static double* ms_pfX2 = (double*) 0x4;
static double* ms_pfXN = (double*) 0x2;

static char*   ms_pcLeastSquaresOutputFile = NULL;
static int     ms_nLeastSquaresOutputMode = 0;
static double  ms_fLeastSquaresOutlierRejectIoverSig = 3.0;
DTREK_EXPORT   int nSetLeastSquaresOutput(char* pcFile,int nMode);


DTREK_EXPORT  int nConvolve(itr<double>& afData,itr<double>& afResp,bool bConvolve,itr<double>& afResult);
DTREK_EXPORT  double fStandardFunction(double* pfFunction,int nPoint,int nNumPoints);

DTREK_EXPORT int nSolveLeastSquares(int nNumFuncs, int nNumPoints, double* pfVals, 
		       double* pfValSigmas,double** ppfFuncs, double* pfArgs, 
		       double* pfSigmas,int* pnUse,double* pfTotalVariance = NULL,itr<int>* panOutlierReject = NULL);

// Quadratic tools.
DTREK_EXPORT int nInitQuadratic(double a9x9fCoeffA[9][9],double a9fCoeffb[9]);
DTREK_EXPORT int nAddQuadratic(double fVal,double fSigma,double a3fVec[3],double a9x9fCoeffA[9][9],double a9fCoeffb[9]);
DTREK_EXPORT int nSolveQuadratic(double a3x3fQuad[3][3],double a9x9fCoeffA[9][9],double a9fCoeffb[9]);

DTREK_EXPORT int nInitQuadratic(int nN,bool bAllocateMemory,double** pfCoeffA /* [nN*(nN+1)+1][nN*(nN+1)+1] */ ,double** pfCoeffb /* [nN*(nN+1)+1] */);
DTREK_EXPORT int nAddQuadratic(int nN,double fVal,double fSigma,double* pfVec /* [nN] */,double* pfCoeffA /* [nN*(nN+1)][nN*(nN+1)+1] */,double* pfCoeffb  /* [nN*(nN+1)+1] */);
DTREK_EXPORT int nSolveQuadratic(int nN,double* pfGMat,double* pfhVec,double* pfConst,double* pfCoeffA,double* pfCoeffb);
DTREK_EXPORT double fEvalQuadratic(int nN,double* pfCoeffA,double* pfCoeffb,double fCoeffC,double* pfVector);

DTREK_EXPORT int nGridInterpClear(int nEntries);
DTREK_EXPORT int nGridInterpAdd(double fPos0,double fPos1,int nEntries,double* pfData);
DTREK_EXPORT int nGridInterpSolve(double fPos0,double fPos1,int nEntries,double* pfData);

DTREK_EXPORT int nBisection(int& nState,double fMin,double fMax,double& fTest,double fTestResult,double fResultTarget,double fTestTol);

DTREK_EXPORT bool bLinearFit(float x[], float y[], int ndata, 
                             float sig[], bool mwt, float *a, float *b, 
                             float *siga, float *sigb, float *chi2, float *q);

bool gammq(float a, float x, float& fOut);
bool gser(float *gamser, float a, float x, float *gln);
float gammln(float xx);
bool gcf(float *gammcf, float a, float x, float *gln);
#endif
