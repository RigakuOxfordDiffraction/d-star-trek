#ifndef DT_DTREKVEC_H
#define DT_DTREKVEC_H
//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTREKVEC.h        Initial author: J.W. Pflugrath           03-May-1995
//    This file is the header file for prototyping the d*TREK vector functions
//    The vector functions may become a class later.
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

//+Include files

#include <math.h>
#include "Dtrek.h"
#include "Cstring.h"
#include "dtarray.h"

//+Definitions and constants

#define SIGNI(a, b)  ((b)< 0 ? -abs(a) : abs(a))
#define SIGNF(a, b)  ((b)< 0 ? -fabs((double)a) : fabs((double)a))
#define SIGND(a, b)  ((b)< 0 ? -fabs(a) : fabs(a))

static float fRADIANS_PER_DEGREE = (3.14159265f/180.0f);
static float fPI                 = (3.14159265f);

static float  Gs_fRADIANS_PER_DEGREE = (3.14159265f/180.0f);
static float  Gs_fPI                 = (3.14159265f);
static double Gs_dRADIANS_PER_DEGREE = (3.1415926535897932384626433832795/180.0);
static double Gs_dPI                 = (3.1415926535897932384626433832795);

#define DTREK_DET_RESO_MAX          999.99f

//+Code begin


DTREK_EXPORT float fDot3D (const float fVec1[3], const float fVec2[3]);
DTREK_EXPORT float fDotND (const int nDim, const float fVec1[], const float fVec2[]);

DTREK_EXPORT void vCross3D (const float fVec1[3], const float fVec2[3], float fResult[3]);

//
// removed the const from the float mat[3][3] parameters in the function
// prototypes.  The DEC compiler is very unhappy if a non-const is passed.
// Could use the following workaround:
//      typedef float MATRIX[3][3];
//      vfunction(const MATRIX fMat1, ... ){...}
// as the DEC compiler seems to be happy if a non-const is passed in this
// manner.  Jim thinks that this is too weird, plus it would probably require
// changing a bunch of code, so this may be best left for a decision at a
// later date.
// tlh - 08 Feb 1996
//

//void vMulMat3DMat3D(const float fMat1[3][3], const float fMat2[3][3],
DTREK_EXPORT void vMulMat3DMat3D(float fMat1[3][3], float fMat2[3][3],
                    float fMat3[3][3]);
DTREK_EXPORT void vMulMatNDMatND(const int nDim, const float *pfMat1, const float *pfMat2,
                    float *pfMat3);

//void vMulMat3DVec3D(const float fMat[3][3], const float fVec[3],
DTREK_EXPORT void vMulMat3DVec3D(float fMat[3][3], const float *pfVec,
                    float *pfVec2);
DTREK_EXPORT void vMulMatNDVecND(const int nDim, const float *pfMat, const float *pfVec1,
                    float *pfVec2);
DTREK_EXPORT void vMulMatMDNDMatMDKD(const int nDimM, const int nDimN, const int nDimK, const float* pfMat1,
                        const float* pfMat2, float* pfMat3);


DTREK_EXPORT void vMulMat3DScalar(const float fMat1[9], const float fScalar, float fMat2[9]);
DTREK_EXPORT void vMulVec3DScalar(const float fVec1[3], const float fScalar, float fVec2[3]);
DTREK_EXPORT void vMulVecNDScalar(const int nDim, const float *pfVec1, const float fScalar,
		     float *pfVec2);

DTREK_EXPORT float fInvMat3D(const float *pfMat1, float *pfMat2);
DTREK_EXPORT float fInvMat2D(const float *pfMat1, float *pfMat2);
inline float fDetMat3D(const float *pfMat1) { return fInvMat3D(pfMat1, (float*)NULL); };
inline float fDetMat2D(const float *pfMat1) { return fInvMat2D(pfMat1, (float*)NULL); };
DTREK_EXPORT float fInvMatND(const int nDim, const float *pfMat1, float *pfMat2);

DTREK_EXPORT void vAddVec3DVec3D(const float fVec1[3], const float fVec2[3], float fVec3[3]);
DTREK_EXPORT void vAddVecNDVecND(const int nDim, const float fVec1[], const float fVec2[],
		    float fVec3[]);

DTREK_EXPORT void vSubVec3DVec3D(const float fVec1[3], const float fVec2[3], float fVec3[3]);
DTREK_EXPORT void vSubVecNDVecND(const int nDim, const float fVec1[], const float fVec2[],
		    float fVec3[]);



DTREK_EXPORT void vTranMat3D(float fMat[3][3]);
DTREK_EXPORT void vTranMatND(const int nDim, float *pfMat);
DTREK_EXPORT void vTranMatMDND(const int nDimM, const int nDimN, const float* pfMat1,float* pfMat2);

DTREK_EXPORT void vCopyMat3D(const float *pfMat1, float *pfMat2);
DTREK_EXPORT void vCopyVec3D(const float *pfVec1, float *pfVec2);
DTREK_EXPORT void vCopyVecND(const int nDim, const float *pfVec1, float *pfVec2);

DTREK_EXPORT float fNormVec3D(float fVec[3]);
DTREK_EXPORT float fLenVec3D(const float fVec[3]);
DTREK_EXPORT float fLenVecND(const int nDim, const float *pfVec);
DTREK_EXPORT void vNormMat3D(float a3x3fMat[3][3]);     // Orthogonalizes a matrix.

DTREK_EXPORT void vListVec3D(float fVec[3]);

DTREK_EXPORT void vListMat3D(float *pfMat);

DTREK_EXPORT void vListMatMN(const int nM, const int nN, float *pfMat);

// Convert a rotation around a vector to a matrix and its derivatives
DTREK_EXPORT void vConvRotVec3DMat3D(const float fRot, const float *pfVec,
			float *pfMat0,
			float *pfMat1=(float*)NULL, float *pfMat2=(float*)NULL,
			float *pfMat3=(float*)NULL, float *pfMat4=(float*)NULL
            );


DTREK_EXPORT void vConv3Rot3Vec3DMat3D(const float fRot0,
			  const float fRot1,
			  const float fRot2,
			  const float fVec0[3],
			  const float fVec1[3],
			  const float fVec2[3],
			  float *fMat0,                 // The resultant 3x3 matrix
              float *pfMatDeriv0 = (float*)NULL,    // And it's derivatives wrt fRot0,fRot1 and fRot2
              float *pfMatDeriv1 = (float*)NULL,
              float *pfMatDeriv2 = (float*)NULL);
        

DTREK_EXPORT int nCalcGetGonConstants(const float *pfVec1,  // Input rotation vector 1
                         const float *pfVec2,  // Input rotation vector 2
                         const float *pfVec3,  // Input rotation vector 3
                         float *pfSubConst,
                         float *pfDivConst,
                         float *pfZeroConst);

//int nDeconvMat3DVec3DRot3(const float fMatrix[3][3],
DTREK_EXPORT int nDeconvMat3DVec3DRot3(float fMatrix[3][3],
			  const float *pfVec0,
			  const float *pfVec1,
			  const float *pfVec2,
			  float *pfRot0, float *pfRot1, float *pfRot2,
			  const float fSign=1.0);

//void vDeconvMat3D3XYZ(const float fMatrix[3][3],
DTREK_EXPORT void vDeconvMat3D3XYZ(float fMatrix[3][3],
		      float *pfRotX, float *pfRotY, float *pfRotZ);

DTREK_EXPORT void vZeroMat3D(float fMatrix[9]);
DTREK_EXPORT void vZeroMat(const int nDim0, const int nDim1, float fMatrix[]);

DTREK_EXPORT void vIdentMat3D(float fMatrix[3][3]);
DTREK_EXPORT void vBuildUpperRightND(const int nDim, float *pfMat);
DTREK_EXPORT void vBuildBasis3D(float a3fVec[3],float a3fVecBasis[3][3]);

DTREK_EXPORT float fSign(const float fVar);

DTREK_EXPORT void
vAffineTransform(const int n0Start,
                 const int n0End,
                 const int n1Start,
                 const int n1End,
                 const int n2Start,
                 const int n2End,
                 const float *pfMatrix,
                 const float *pfTrans,
                 float *pfResult);

DTREK_EXPORT void vEigen(const int nDim, const int nMode, float *pfA, float *pfR);

DTREK_EXPORT int nSolveQuadratic(const double fA, const double fB, const double fC,double *pfV);


DTREK_EXPORT void vTRED2(const int nDim, float **ppfMat, float *pfDiag, float *pfOffDiag);
DTREK_EXPORT void vTQLI(const int nDim, float *pfDiag, float *pfOffDiag,
	   float **ppfEigenVecs);
DTREK_EXPORT void vEigsrt(const int nDim, float *pfEigval, float **ppfEigvec);


DTREK_EXPORT int
nMatrixCompress(const int nNumRowsCols, const int *pnDelFlags,
		float *pfMatrix, float *pfVector);

DTREK_EXPORT int
nint (const float fTemp);

DTREK_EXPORT int nFFT1D(const int nSize, const int nSign, float *pfData);



// Double versions of routines.

DTREK_EXPORT double fDot3D (const double fVec1[3], const double fVec2[3]);
DTREK_EXPORT double fDotND (const int nDim, const double fVec1[], const double fVec2[]);

DTREK_EXPORT void vCross3D (const double fVec1[3], const double fVec2[3], double fResult[3]);

//
// removed the const from the float mat[3][3] parameters in the function
// prototypes.  The DEC compiler is very unhappy if a non-const is passed.
// Could use the following workaround:
//      typedef float MATRIX[3][3];
//      vfunction(const MATRIX fMat1, ... ){...}
// as the DEC compiler seems to be happy if a non-const is passed in this
// manner.  Jim thinks that this is too weird, plus it would probably require
// changing a bunch of code, so this may be best left for a decision at a
// later date.
// tlh - 08 Feb 1996
//

DTREK_EXPORT void vMulMat3DMat3D(double fMat1[3][3], double fMat2[3][3],
                    double fMat3[3][3]);
DTREK_EXPORT void vMulMatNDMatND(const int nDim, const double *pfMat1, const double *pfMat2,
                    double *pfMat3);

DTREK_EXPORT void vMulMat3DVec3D(double fMat[3][3], const double *pfVec,
                    double *pfVec2);
DTREK_EXPORT void vMulMatNDVecND(const int nDim, const double *pfMat, const double *pfVec1,
                    double *pfVec2);
DTREK_EXPORT void vMulMatMDNDMatMDKD(const int nDimM, const int nDimN, const int nDimK, const double* pfMat1,
                        const double* pfMat2, double* pfMat3);


DTREK_EXPORT void vMulMat3DScalar(const double fMat1[9], const double fScalar, double fMat2[9]);
DTREK_EXPORT void vMulVec3DScalar(const double fVec1[3], const double fScalar, double fVec2[3]);
DTREK_EXPORT void vMulVecNDScalar(const int nDim, const double *pfVec1, const double fScalar,
		     double *pfVec2);

DTREK_EXPORT double fInvMat3D(const double *pfMat1, double *pfMat2);
DTREK_EXPORT double fInvMat2D(const double *pfMat1, double *pfMat2);
inline double fDetMat3D(const double *pfMat1) { return fInvMat3D(pfMat1, (double*)NULL); };
inline double fDetMat2D(const double *pfMat1) { return fInvMat2D(pfMat1, (double*)NULL); };
DTREK_EXPORT double fInvMatND(const int nDim, const double *pfMat1, double *pfMat2);

DTREK_EXPORT void vAddVec3DVec3D(const double fVec1[3], const double fVec2[3], double fVec3[3]);
DTREK_EXPORT void vAddVecNDVecND(const int nDim, const double fVec1[], const double fVec2[],
		    double fVec3[]);

DTREK_EXPORT void vSubVec3DVec3D(const double fVec1[3], const double fVec2[3], double fVec3[3]);
DTREK_EXPORT void vSubVecNDVecND(const int nDim, const double fVec1[], const double fVec2[],
		    double fVec3[]);



DTREK_EXPORT void vTranMat3D(double fMat[3][3]);
DTREK_EXPORT void vTranMatND(const int nDim, double *pfMat);
DTREK_EXPORT void vTranMatMDND(const int nDimM, const int nDimN, const double* pfMat1,double* pfMat2);

DTREK_EXPORT void vCopyMat3D(const double *pfMat1, double *pfMat2);
DTREK_EXPORT void vCopyVec3D(const double *pfVec1, double *pfVec2);
DTREK_EXPORT void vCopyVecND(const int nDim, const double *pfVec1, double *pfVec2);

DTREK_EXPORT double fNormVec3D(double fVec[3]);
DTREK_EXPORT double fLenVec3D(const double fVec[3]);
DTREK_EXPORT double fLenVecND(const int nDim, const double *pfVec);
DTREK_EXPORT void vNormMat3D(double a3x3fMat[3][3]);     // Orthogonalizes a matrix.

DTREK_EXPORT int nInvMat3D(int* pnMat1,int* pnMat2);
DTREK_EXPORT void vMulMat3DVec3D(int* pnMat,int* pnVec1,int* pnVec2);
DTREK_EXPORT bool bDivides(int nDim,int* pnVec,int nFactor);
DTREK_EXPORT bool bDivides(int nDim,int* pnVec);
DTREK_EXPORT void vCopyMat3D(int* pnMat1,int* pnMat2);
DTREK_EXPORT void vCopyMat3D(int* pnMat1,double* pfMat2);


DTREK_EXPORT void vListVec3D(double fVec[3]);

DTREK_EXPORT void vListMat3D(double *pfMat);

DTREK_EXPORT void vListMatMN(const int nM, const int nN, double *pfMat);

// Convert a rotation around a vector to a matrix and its derivatives
DTREK_EXPORT void vConvRotVec3DMat3D(const double fRot, const double *pfVec,
			double *pfMat0,
			double *pfMat1=(double*)NULL, double *pfMat2=(double*)NULL,
			double *pfMat3=(double*)NULL, double *pfMat4=(double*)NULL);

DTREK_EXPORT void vConv3Rot3Vec3DMat3D(const double fRot0,
			  const double fRot1,
			  const double fRot2,
			  const double fVec0[3],
			  const double fVec1[3],
			  const double fVec2[3],
   			  double *fMat0,                 // The resultant 3x3 matrix
              double *pfMatDeriv0 = (double*)NULL,    // And it's derivatives wrt fRot0,fRot1 and fRot2
              double *pfMatDeriv1 = (double*)NULL,
              double *pfMatDeriv2 = (double*)NULL);

DTREK_EXPORT int nCalcGetGonConstants(const double *pfVec1,  // Input rotation vector 1
                         const double *pfVec2,  // Input rotation vector 2
                         const double *pfVec3,  // Input rotation vector 3
                         double *pfSubConst,
                         double *pfDivConst,
                         double *pfZeroConst);

DTREK_EXPORT int nDeconvMat3DMat3D(double a3fRots[3],
                                   double a3x3fRotVecs[3][3],
                                   double a3x3fDeltaRot[3][3],
                                   bool bPreMultiply,
                                   double* pfResultMat = NULL);

DTREK_EXPORT int nDeconvMat3DVec3DRot1(double fMatrix[3][3],double* pfVector,double* pfRot);

//int nDeconvMat3DVec3DRot3(const double fMatrix[3][3],
DTREK_EXPORT int nDeconvMat3DVec3DRot3(double fMatrix[3][3],
			  const double *pfVec0,
			  const double *pfVec1,
			  const double *pfVec2,
			  double *pfRot0, double *pfRot1, double *pfRot2,
			  const double fSign=1.0);

//void vDeconvMat3D3XYZ(const double fMatrix[3][3],
DTREK_EXPORT void vDeconvMat3D3XYZ(double fMatrix[3][3],
		      double *pfRotX, double *pfRotY, double *pfRotZ);

DTREK_EXPORT void vZeroMat3D(double fMatrix[9]);
DTREK_EXPORT void vZeroMat(const int nDim0, const int nDim1, double fMatrix[]);

DTREK_EXPORT void vIdentMat3D(double fMatrix[3][3]);
DTREK_EXPORT void vBuildUpperRightND(const int nDim, double *pfMat);
DTREK_EXPORT void vBuildBasis3D(double a3fVec[3],double a3fVecBasis[3][3]);
DTREK_EXPORT void vConvertEllipsoidRelative(float* a2x2fEllipsoidAIn_,float* a2fEllipsoidb_,float* pfEllipsoidc_,float fCentIn0,float fCentIn1,float fCentOut0,float fCentOut1);
DTREK_EXPORT void vConvertEllipsoidRelative(double* a2x2fEllipsoidAIn,double* a2fEllipsoidb,double* pfEllipsoidc,double fCentIn0,double fCentIn1,double fCentOut0,double fCentOut1);
DTREK_EXPORT double fCalcRotationOffset(double a3fE3[3],double a3fS0[3],double a3fXRC[3],double* pfXRCRot = NULL);
DTREK_EXPORT int nGetEllipsoidTiltAngleCenter(double* a2x2fEllipsoidA,double* a2fEllispoidb,double fEllipsoidc,double* a2fRelOffset,double* a2fMajorMinor,double* pfRotOffset);
DTREK_EXPORT int nMergeEllipsoids(double* pfOffset1,double* a2x2fEllipsoidA1,double* a2fEllipsoidb1,double fEllipsoidc1,double* pfOffset2,double* a2x2fEllipsoidA2,double* a2fEllipsoidb2,double fEllipsoidc2,double* a2x2fEllipsoidAOut,double* a2fEllipsoidbOut,double* pfEllipsoidcOut);
DTREK_EXPORT int nBuildEllipse(double fRadius,double* a2x2fEllipsoidA,double* a2fEllispoidb,double* pfEllipsoidc);
DTREK_EXPORT int nBuildEllipse(double a2fMajorMinor[2],double a2fRelOffset[2],double fRotOffset,double* a2x2fEllipsoidAOut,double* a2fEllipsoidbOut,double* pfEllipsoidcOut);
DTREK_EXPORT double fEvalEllipse(double fPoint0,double fPoint1,double* a2x2fEllipsoidA,double* a2fEllispoidb,double* pfEllipsoidc);
DTREK_EXPORT double fSign(const double fVar);
DTREK_EXPORT int nUpperLowerEnvelope(int nMode,itr<double>& afDataYIn,itr<double>& afDataYOut,int nRangeLower = 0,int nRangeUpper = -1,itr<double>* pafDataX = NULL);
const int g_nFindApplyUpperEnvelope = 0;
const int g_nFindApplyLowerEnvelope = 1;


DTREK_EXPORT void
vAffineTransform(const int n0Start,
                 const int n0End,
                 const int n1Start,
                 const int n1End,
                 const int n2Start,
                 const int n2End,
                 const double *pfMatrix,
                 const double *pfTrans,
                 double *pfResult);

DTREK_EXPORT void vEigen(const int nDim, const int nMode, double *pfA, double *pfR);
DTREK_EXPORT int  nEigen2D(double a2x2fMat[2][2],double a2x2fEigenVec[2][2],double a2fEignVals[2]);

DTREK_EXPORT void vTRED2(const int nDim, double **ppfMat, double *pfDiag, double *pfOffDiag);
DTREK_EXPORT void vTQLI(const int nDim, double *pfDiag, double *pfOffDiag,
	   double **ppfEigenVecs);
DTREK_EXPORT void vEigsrt(const int nDim, double *pfEigval, double **ppfEigvec);

DTREK_EXPORT int nCompVecND(const int nDim,double *pfComp1,double *pfComp2);

DTREK_EXPORT int
nMatrixCompress(const int nNumRowsCols, const int *pnDelFlags,
		double *pfMatrix, double *pfVector);

DTREK_EXPORT int nFFT1D(const int nSize, const int nSign, double *pfData);
DTREK_EXPORT int nRealFFT1D(const int nSize, const int nSign, double *pfData);

DTREK_EXPORT int nSplineCompute(double* pfX,double* pfY,int nPoints,double fYp0,double fYpN,double* pfYpp);
DTREK_EXPORT double fSplineCalc(double fX,double fX1,double fY1,double fYpp1,double fX2,double fY2,double fYpp2);
DTREK_EXPORT int nConvexHull(int nNumPoints,double* aa2fPoints,int* anPointsInHull);
DTREK_EXPORT int nInvFunc(double fStartX,double fStepX,itr<double>& afFunctionY,double fInvStartY,double fInvEndY,double fInvStepY,itr<double>& afInvFunctionX,bool bForceMonotonic = false);
DTREK_EXPORT int nPlot(const char* pcName,itr<double>& afY,itr<double>* pafX = NULL);
DTREK_EXPORT int nPlot(const char* pcName,itr<double>& afY,double fStartX,double fStepX);


DTREK_EXPORT void vCopyMat3D(const float *pfMat1, double *pfMat2);

DTREK_EXPORT void vCopyMat3D(const double *pfMat1, float *pfMat2);

DTREK_EXPORT void vCopyVec3D(const float *pfVec1, double *pfVec2);

DTREK_EXPORT void vCopyVec3D(const double *pfVec1, float *pfVec2);

DTREK_EXPORT void vCopyVecND(const int nDim, const float *pfVec1, double *pfVec2);

DTREK_EXPORT void vCopyVecND(const int nDim, const double *pfVec1, float *pfVec2);

DTREK_EXPORT float fDot3D (const double fVec1[3], const float fVec2[3]);
DTREK_EXPORT float fDot3D (const float fVec1[3], const double fVec2[3]);

DTREK_EXPORT float* pfCast3D(const double*,int nWhich = 0);            // Use these routines when you need to send an argument to a function that requires a float as input
DTREK_EXPORT float* pfCastND(int nDim,const double*,int nWhich = 0);   // Multidimensional versino of above.  Maximum dimension on this, so be carefull.  It's only a helper!
DTREK_EXPORT double* pfCast3D(const float*,int nWhich = 0);            
DTREK_EXPORT double* pfCastND(int nDim,const float*,int nWhich = 0);   

// Sort routines for qsort().

#ifndef WIN32
extern "C"
{
#endif // !WIN32


DTREK_EXPORT int float_cmp(const void* a,const void* b);
DTREK_EXPORT int double_cmp(const void* a,const void* b);
DTREK_EXPORT int int_cmp(const void* a,const void* b);
DTREK_EXPORT int float_cmp_rel(const void* a,const void* b);
DTREK_EXPORT int double_cmp_rel(const void* a,const void* b);

DTREK_EXPORT int int_cmp_rel(const void* a,const void* b);
DTREK_EXPORT int int_cmp_rel_2(const void* a,const void* b);


DTREK_EXPORT int unsigned_short_int_cmp(const void* a,const void* b);
DTREK_EXPORT int Cstring_cmp_rel(const void* a,const void* b);
DTREK_EXPORT int float_cmp_tol(const void* a,const void* b);
DTREK_EXPORT int double_cmp_tol(const void* a,const void* b);

#ifndef WIN32
}
#endif // !WIN32

DTREK_EXPORT extern float g_fCmpFloatsTol;
DTREK_EXPORT extern double g_fCmpDoublesTol;
DTREK_EXPORT extern int* g_pnCmpInts;
DTREK_EXPORT extern float* g_pfCmpFloats;
DTREK_EXPORT extern double* g_pfCmpDoubles;
DTREK_EXPORT extern Cstring* g_psCmpCstrings;

DTREK_EXPORT extern itr<int> g_anSwapSizes;
DTREK_EXPORT extern itr<void*> g_apvSwapPointers;
DTREK_EXPORT int qsort_swap_arrays(const void  *elem1, const void *elem2,size_t width);
DTREK_EXPORT void qsortswap(void* base,size_t num,size_t width,int (*compare )(const void *elem1, const void *elem2 ),int(*myswap)(const void  *elem1, const void *elem2,size_t width));

// Mapping functions.
int nCreateMapping(itr<int>& anMap,int nNumEntries);
int nAddToMapping(itr<int>& anMap,int nFrom,int nTo,bool bPruneDuplicates = false);
int nGetMapping(itr<int>& anMap,int nFrom,itr<int>& anTo,bool bFIFO = true);


// Simplex function.
int nSimplexRun(double fValue,itr<double>& afDelta,itr<double>* pafStep = NULL);

DTREK_EXPORT extern int m_nMultiOptimizeTrialCount;
DTREK_EXPORT extern bool m_bMultiOptimizePrint;
DTREK_EXPORT extern double m_fMultiOptimizeMaxDist;
DTREK_EXPORT int nMultiOptimize(double& fValue,itr<double>& afParams,bool bMinimize = true,int nBatchSize = -1,int nBestToPickInBatch =  -1,int nNumTrials = -1,int nNumInitialSolutions = -1);

#define DTREK_VEC_ROUNDOFF_CEILING      0x0001
#define DTREK_VEC_ROUNDOFF_FLOOR        0x0002
DTREK_EXPORT double dRoundOff(double dIn, int nDecDigits, DTREK_WORD wCtrl=0U);
DTREK_EXPORT double dGetNearestMultipleInteger(double dIn, int nMultInteger, DTREK_WORD wCtrl=0U);
DTREK_EXPORT int nGetNearestIntFromVector(double dIn, std::vector<int>& anCheckValues, DTREK_WORD wCtrl=0U);

DTREK_EXPORT bool bIsSkipCommandLineArgument(const char* sArg);
#endif   // DT_DTREKVEC_H
