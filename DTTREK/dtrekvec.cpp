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
// dtrekvec.cc     Initial author: J.W. Pflugrath           03-May-1995
//    This file contains some vector and matrix functions
//    that do not belong to any class.
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
//+Description

//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "dtrekvec.h"           // Function prototypes
#include "dtarray.h"
#include "dtreksys.h"
#include "dtsvd.h"

#ifdef SUNOS
// This gets the SunOS swap routine
#include <algorithm>
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


//+Definitions, constants, and initialization of static member variables

//+Code begin

//+Public functions

float fDot3D (const float fVec1[3], const float fVec2[3])
{
  return (fVec1[0] * fVec2[0]  +  fVec1[1] * fVec2[1]  + fVec1[2] * fVec2[2]);
}

float fDotND (const int nDim, const float fVec1[], const float fVec2[])
{
  float fSum = 0.0;
  for (int i = 0; i < nDim; i++) fSum = fSum + (fVec1[i] * fVec2[i]);
  return (fSum);
}

void vCross3D (const float fVec1[3], const float fVec2[3], float fResult[3])
{
  fResult[0] = fVec1[1] * fVec2[2]  - fVec1[2] * fVec2[1];
  fResult[1] = fVec1[2] * fVec2[0]  - fVec1[0] * fVec2[2];
  fResult[2] = fVec1[0] * fVec2[1]  - fVec1[1] * fVec2[0];
}

//void vMulMat3DMat3D(const float fMat1[3][3], const float fMat2[3][3],
void vMulMat3DMat3D(float fMat1[3][3], float fMat2[3][3],
                    float fMat3[3][3])
{
  int i, j, k;
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  fMat3[j][i] = 0.0;
	  for (k = 0; k < 3; k++)
	    {
	      fMat3[j][i] = fMat3[j][i]  + fMat1[k][i] * fMat2[j][k];
	    }
	}
    }
}

void vMulMatNDMatND(const int nDim, const float *pfMat1, const float *pfMat2,
                    float *pfMat3)
{
  // Multiply an nDim x nDim matrix by another matrix
  // The matrix is stored as in Fortran, so Mat[column][row]
  // pfMat3 must be different from pfMat1 and pfMat2;

  int i, j, k;
  const float *pfTempMat1, *pfTempMat2;
  float *pfTempMat3;

  pfTempMat3 = pfMat3;

  for (i = 0; i < nDim; i++)
    {
      pfTempMat3 = pfMat3 + i;
      for (j = 0; j < nDim; j++)
	{
	  *pfTempMat3 = 0.0;
	  pfTempMat2 = pfMat2 + (j * nDim);
	  pfTempMat1 = pfMat1 + i;
	  for (k = 0; k < nDim; k++)
	    {
	      *pfTempMat3 = *pfTempMat3  + *pfTempMat1 * (*pfTempMat2++);
	      pfTempMat1  = pfTempMat1 + nDim;
	    }
	  pfTempMat3 = pfTempMat3 + nDim;
	}
    }
}


void vMulMatMDNDMatMDKD(const int nDimM, const int nDimN, const int nDimK, const float* pfMat1,
                        const float* pfMat2, float* pfMat3) {
  // Multiply an nDimM x nDimN matrix by another matrix with dimensions nDimN x nDimK
  // The matrix is stored as in Fortran, so Mat[column][row]
  // pfMat3 must be different from pfMat1 and pfMat2;

  int i, j, k;
  const float *pfTempMat1, *pfTempMat2;
  float *pfTempMat3;

  pfTempMat3 = pfMat3;

  for (i = 0; i < nDimM; i++)
    {
      pfTempMat3 = pfMat3 + i;
      for (j = 0; j < nDimK; j++)
	{
	  *pfTempMat3 = 0.0;
	  pfTempMat2 = pfMat2 + (j * nDimN);
	  pfTempMat1 = pfMat1 + i;
	  for (k = 0; k < nDimN; k++)
	    {
	      *pfTempMat3 = *pfTempMat3  + *pfTempMat1 * (*pfTempMat2++);
	      pfTempMat1  = pfTempMat1 + nDimM;
	    }
	  pfTempMat3 = pfTempMat3 + nDimM;
	}
    }
};



//void vMulMat3DVec3D(const float fMat[3][3], const float fVec1[3],
void vMulMat3DVec3D(float fMat[3][3], const float *pfVec1,
                    float *pfVec2)
{
  // pfVec1 and pfVec2 cannot be the same!
  int i, j;
  for (i = 0; i < 3; i++)
    {
      pfVec2[i] = 0.0;
      for (j = 0; j < 3; j++)
	{
	  pfVec2[i] = pfVec2[i] + (fMat[j][i] * pfVec1[j]);
	}
    }
}

void vMulMatNDVecND(const int nDim, const float *pfMat, const float *pfVec1,
                    float *pfVec2)
{
  // Multiply an nDim x nDim matrix by nDim vector.
  // The matrix is stored as in Fortran, so Mat[column][row]
  // pfVec1 and pfVec2 cannot be the same!

  const float *pfTempMat;
  const float *pfTempVec1;
  float       *pfTempVec2;
  int i, j;

  pfTempVec2 = pfVec2;
  for (i = 0; i < nDim; i++)
    {
      *pfTempVec2 = 0.0;
      pfTempVec1  = pfVec1;
      pfTempMat   = pfMat + i;
      for (j = 0; j < nDim; j++)
	{
	  *pfTempVec2 = *pfTempVec2 + (*pfTempMat * *pfTempVec1++);
	  pfTempMat  = pfTempMat + nDim;
	}
      pfTempVec2++;
    }
}

void vMulMat3DScalar(const float fMat1[9], const float fScalar, float fMat2[9])
{
  for (int i = 0; i < 9; i++)
    fMat2[i] = fMat1[i] * fScalar;
}

void vMulVec3DScalar(const float fVec1[3], const float fScalar, float fVec2[3])
{
  for (int i = 0; i < 3; i++)
    fVec2[i] = fVec1[i] * fScalar;
}

void vMulVecNDScalar(const int nDim, const float *pfVec1, const float fScalar,
		     float *pfVec2)
{
  for (int i = 0; i < nDim; i++)
    pfVec2[i] = pfVec1[i] * fScalar;
}

float fInvMat3D(const float *pfMat1, float *pfMat2)
{
  // Invert a 3x3 matrix.  Return the determinant.
  // Mat1 and Mat2 must be different matrices

  float fDet;

  fDet  =   pfMat1[0] * (pfMat1[4] * pfMat1[8]  -  pfMat1[7] * pfMat1[5])
          - pfMat1[1] * (pfMat1[3] * pfMat1[8]  -  pfMat1[5] * pfMat1[6])
          + pfMat1[2] * (pfMat1[3] * pfMat1[7]  -  pfMat1[4] * pfMat1[6]);
  if ( (0.0 != fDet) && (NULL != pfMat2) )
    {
      pfMat2[0] =  (pfMat1[4] * pfMat1[8]  -  (pfMat1[5] * pfMat1[7])) / fDet;
      pfMat2[1] = -(pfMat1[1] * pfMat1[8]  -  (pfMat1[2] * pfMat1[7])) / fDet;
      pfMat2[2] =  (pfMat1[1] * pfMat1[5]  -  (pfMat1[2] * pfMat1[4])) / fDet;
      pfMat2[3] = -(pfMat1[3] * pfMat1[8]  -  (pfMat1[5] * pfMat1[6])) / fDet;
      pfMat2[4] =  (pfMat1[0] * pfMat1[8]  -  (pfMat1[2] * pfMat1[6])) / fDet;
      pfMat2[5] = -(pfMat1[0] * pfMat1[5]  -  (pfMat1[2] * pfMat1[3])) / fDet;
      pfMat2[6] =  (pfMat1[3] * pfMat1[7]  -  (pfMat1[4] * pfMat1[6])) / fDet;
      pfMat2[7] = -(pfMat1[0] * pfMat1[7]  -  (pfMat1[1] * pfMat1[6])) / fDet;
      pfMat2[8] =  (pfMat1[0] * pfMat1[4]  -  (pfMat1[1] * pfMat1[3])) / fDet;
    }
  return (fDet);
}

float fInvMat2D(const float *pfMat1, float *pfMat2)
{
    float fDet;

    fDet = pfMat1[0]*pfMat1[3] - pfMat1[1]*pfMat1[2];

    if ( (0.0 != fDet) && (NULL != pfMat2))
    {
        pfMat2[0] = pfMat1[3] / fDet;
        pfMat2[3] = pfMat1[0] / fDet;
        pfMat2[1] = - pfMat1[1] / fDet;
        pfMat2[2] = - pfMat1[2] / fDet;
    };
    return (fDet);
};

float fInvMatND(const int nDim, const float *pfMat1, float *pfMat2)
{
  // Invert a nDim x nDim matrix by Gauss-Jordan elimination.
  // We assume the matrix is small and we don't care about using
  // extra storage.  See "Numerical Recipes" ....
  // Matrices are stored as in Fortran  A[column][row] so be careful.
  // pfMat1 and pfMat2 may be the same storage locations matrices
  // Returns (1.) if all Ok, (0.0) if not OK.

  int i, j, k, l, ll;
  float fTemp;
//  float fMat[nDim][nDim];
//  int   nPiv[nDim], nIndexRow[nDim], nIndexCol[nDim];  // Bookkeeping
  int   *nPiv, *nIndexRow, *nIndexCol;  // Bookkeeping
  int   nRow, nCol;
  float fBig;

  float *fMatLinear,*fptr1,*fptr2;

  fMatLinear = new float[nDim*nDim];
  nPiv       = new int[nDim];
  nIndexRow  = new int[nDim];
  nIndexCol  = new int[nDim];

  // Make a copy of the input matrix

//  vCopyVecND(nDim * nDim, pfMat1, (float *) fMat);
  vCopyVecND(nDim * nDim, pfMat1, fMatLinear);

  // "Zero" the bookkeeping pivot array

  for (i = 0; i < nDim; i++)
    nPiv[i] = 0;

  for (i = 0; i < nDim; i++)
    {
      fBig = 0.0;
      for (j = 0; j < nDim; j++)
	{
	  if (nPiv[j] != 1)
	    {
	      fptr1 = fMatLinear+j;
	      for (k = 0; k < nDim; k++)
		{
		  if (nPiv[k] == 0)
		    {
		      // fTemp = fMat[k][j];
		      fTemp = *(fptr1+k*nDim);
		      if (fTemp < 0.0) fTemp = -fTemp; // abs function
		      if (fTemp >= fBig)
			{
			  fBig = fTemp;
			  nRow = j;
			  nCol = k;
			}
		    }
		  else if (nPiv[k] > 1)
		    {
		      cout << "Singular matrix" << endl;
		      return (0.0);
		    }
		} // end k loop
	    }
	} // end j loop

      nPiv[nCol]++;

      // We now have the pivot element ...

      if (nRow != nCol)
	{
	  fptr1 = fMatLinear+nRow;
	  fptr2 = fMatLinear+nCol;
	  for (l = 0; l < nDim; l++, fptr1 += nDim, fptr2 += nDim)
	    {
	      //	fTemp         = fMat[l][nRow];
	      //	fMat[l][nRow] = fMat[l][nCol];
	      //	fMat[l][nCol] = fTemp;
	      fTemp  = *fptr1;
	      *fptr1 = *fptr2;
	      *fptr2 = fTemp;
	    } // end l loop
	}

      nIndexRow[i] = nRow;
      nIndexCol[i] = nCol;

      //    if (fMat[nCol][nCol] == 0.0) {
      fptr1 = fMatLinear+nCol*(nDim+1);
      if(0.0 == *fptr1)
	{
	  cout << "Singular matrix" << endl;
	  return (0.0);
	}
      //    fTemp = 1.0 / fMat[nCol][nCol];
      //    fMat[nCol][nCol] = 1.0;
      fTemp = 1.0f/(*fptr1);
      *fptr1 = 1.0f;
      fptr2 = fMatLinear+nCol;
      for (l = 0; l < nDim; l++, fptr2 += nDim)
	{
	  //      fMat[l][nCol] = fMat[l][nCol] * fTemp;
	  *fptr2 *= fTemp;
	}
      for (ll = 0; ll < nDim; ll++)
	{
	  fptr1 = fMatLinear+ll;
	  fptr2 = fptr1+nCol*nDim;
	  if (ll != nCol)
	    {
	      //	fTemp          = fMat[nCol][ll];
	      //	fMat[nCol][ll] = 0.0;
	      fTemp = *fptr2;
	      *fptr2 = 0.0;
	      fptr2 = fMatLinear+nCol;
	      for (l = 0; l < nDim; l++, fptr1 += nDim, fptr2 += nDim)
		{
		  //	  fMat[l][ll] = fMat[l][ll]  -  fMat[l][nCol] * fTemp;
		  *fptr1 -= (*fptr2)*fTemp;
		}
	    }
	} // end ll loop
    } // end i loop (main loop over columns of the reduction)
  for (l = nDim-1; l >= 0; l--)
    {
      if (nIndexRow[l] != nIndexCol[l])
	{
	  fptr1 = fMatLinear+nDim*nIndexRow[l];
	  fptr2 = fMatLinear+nDim*nIndexCol[i];
	  for (k = 0; k < nDim; k++, fptr1++, fptr2++)
	    {
	      //	fTemp                 = fMat[nIndexRow[l]][k];
	      //	fMat[nIndexRow[l]][k] = fMat[nIndexCol[l]][k];
	      //	fMat[nIndexCol[l]][k] = fTemp;
	      fTemp = *fptr1;
	      *fptr1 = *fptr2;
	      *fptr2 = fTemp;
	    }
	}
    }
//  vCopyVecND(nDim * nDim, (float *) fMat, pfMat2);
  vCopyVecND(nDim * nDim, fMatLinear, pfMat2);

  delete [] fMatLinear;

  //  delete [] fMat;
  delete [] nPiv;
  delete [] nIndexRow;
  delete [] nIndexCol;

  return (1.0);
}

void vAddVec3DVec3D(const float fVec1[3], const float fVec2[3], float fVec3[3])
{
  for (int i = 0; i < 3; i++)	
     fVec3[i] = fVec1[i] + fVec2[i];
}

void vAddVecNDVecND(const int nDim, const float fVec1[], const float fVec2[],
		    float fVec3[])
{
  for (int i = 0; i < nDim; i++)	
     fVec3[i] = fVec1[i] + fVec2[i];
}

void vSubVec3DVec3D(const float fVec1[3], const float fVec2[3], float fVec3[3])
{
  for (int i = 0; i < 3; i++)	
     fVec3[i] = fVec1[i] - fVec2[i];
}

void vSubVecNDVecND(const int nDim, const float fVec1[], const float fVec2[],
		    float fVec3[])
{
  for (int i = 0; i < nDim; i++)	
     fVec3[i] = fVec1[i] - fVec2[i];
}

void vTranMat3D(float fMat[3][3])
{
  // Transpose a 3x3 matrix in place

  float fTemp;
  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = i+1; j < 3; j++) {
      fTemp      = fMat[i][j];
      fMat[i][j] = fMat[j][i];
      fMat[j][i] = fTemp;
    }
  }
}

void vTranMatND(const int nDim, float *pfMat)
{
  // Transpose a nDim by nDim matrix in place

  float fTemp;
  int i, j, k, m;
  for (i = 0; i < nDim; i++)
    {
      for (j = i+1; j < nDim; j++)
	{
	  m = i * nDim + j;             // m is [i][j]
	  k = j * nDim + i;             // k is [j][i]
	  fTemp = pfMat[k];
	  pfMat[k] = pfMat[m];
	  pfMat[m] = fTemp;
	}
    }
}

void vTranMatMDND(const int nDimM, const int nDimN, const float* pfMat1,float* pfMat2) 
{
  // Transpose an nDimM x nDimN matrix in place
  int i, j;

  for (i = 0; i < nDimM; i++)
      for (j = 0; j< nDimN; j++)
          pfMat2[ j + nDimN * i ] = pfMat1[ i + nDimM * j ];
} 


void vCopyMat3D(const float *pfMat1, float *pfMat2)
{
  for (int i = 0; i < 9; i++) pfMat2[i] = pfMat1[i];
}

void vCopyVec3D(const float *pfVec1, float *pfVec2)
{
  for (int i = 0; i < 3; i++) pfVec2[i] = pfVec1[i];
}


void vCopyVecND(const int nDim, const float *pfVec1, float *pfVec2)
{
  for (int i = 0; i < nDim; i++) pfVec2[i] = pfVec1[i];
}


float fNormVec3D(float fVec[3])
{
  float fTemp;
  fTemp = fDot3D(fVec, fVec);
  if (fTemp > 0.0) {
    fTemp = 1.0 / sqrt((double)fTemp);
    fVec[0] = fVec[0] * fTemp;
    fVec[1] = fVec[1] * fTemp;
    fVec[2] = fVec[2] * fTemp;
  }
  return (fTemp);
}

float fLenVec3D(const float fVec[3])
{
  float fTemp;
  fTemp = fDot3D(fVec, fVec);
  fTemp = sqrt((double)fTemp);
  return (fTemp);
}

float fLenVecND(const int nDim, const float *pfVec)
{
  float fTemp;
  fTemp = fDotND(nDim, pfVec, pfVec);
  fTemp = sqrt((double)fTemp);
  return (fTemp);
}

void vNormMat3D(float a3x3fMat[3][3])
{
    float a3fSub[3];

    // gram-schmid (spell?) the matrix.
    fNormVec3D(a3x3fMat[0]);
    vMulVec3DScalar(a3x3fMat[0],fDot3D(a3x3fMat[0],a3x3fMat[1]),a3fSub);
    vSubVec3DVec3D(a3x3fMat[1],a3fSub,a3x3fMat[1]);
    fNormVec3D(a3x3fMat[1]);
    vMulVec3DScalar(a3x3fMat[0],fDot3D(a3x3fMat[0],a3x3fMat[2]),a3fSub);
    vSubVec3DVec3D(a3x3fMat[2],a3fSub,a3x3fMat[2]);
    vMulVec3DScalar(a3x3fMat[1],fDot3D(a3x3fMat[1],a3x3fMat[2]),a3fSub);
    vSubVec3DVec3D(a3x3fMat[2],a3fSub,a3x3fMat[2]);
    fNormVec3D(a3x3fMat[2]);
    return;
};

void vListVec3D(float fVec[3])
{
 cout << "Vector is: " << fVec[0] << ", " << fVec[1] << ", " << fVec[2] << endl;
}

void vListMat3D(float *pfMat)
{
  float *pfTemp;
  pfTemp = pfMat;
  cout << "Matrix is:" << endl;
  for (int i = 0; i < 3; i++)
    {
      cout << "[" << i << "][*]: " << pfTemp[0] << ", "
           << pfTemp[1] << ", " << pfTemp[2] << endl;
      pfTemp += 3;
    }
}

void vListMatMN(const int nDim0, const int nDim1, float *pfMat)
{
  // List an MxN matrix, where M is fastest changing address

  int i, j;
  float *pfTemp;
  pfTemp = pfMat;
  cout << "Matrix is:" << endl;
  for (i = 0; i < nDim1; i++)
    {
      cout << "[" << i << "][*]: ";
      for (j = 0; j < nDim0-1; j++)
	{
	  cout << *pfTemp++ << ", ";
	}
      cout << *pfTemp++ << '\n';
    }
  cout << endl;
}

// Convert a rotation (in degrees) around a vector to a matrix and
// possibly its derivatives

void vConvRotVec3DMat3D(const float fRot, const float *pfVec,
			float *pfMat0,float *pfMat1,float *pfMat2,float *pfMat3,float *pfMat4)
{


  double a5x3x3fMats[5][3][3];
  double a3fVec[3];

  vCopyVec3D(pfVec,&a3fVec[0]);
  vConvRotVec3DMat3D(fRot,&a3fVec[0],
      (pfMat0)?(&a5x3x3fMats[0][0][0]):NULL,
      (pfMat1)?(&a5x3x3fMats[1][0][0]):NULL,
      (pfMat2)?(&a5x3x3fMats[2][0][0]):NULL,
      (pfMat3)?(&a5x3x3fMats[3][0][0]):NULL,
      (pfMat4)?(&a5x3x3fMats[4][0][0]):NULL);

  if (pfMat0) 
      vCopyMat3D(&a5x3x3fMats[0][0][0],pfMat0);
  if (pfMat1) 
      vCopyMat3D(&a5x3x3fMats[1][0][0],pfMat1);
  if (pfMat2) 
      vCopyMat3D(&a5x3x3fMats[2][0][0],pfMat2);
  if (pfMat3) 
      vCopyMat3D(&a5x3x3fMats[3][0][0],pfMat3);
  if (pfMat4) 
      vCopyMat3D(&a5x3x3fMats[4][0][0],pfMat4);

}

// Returns a matrix composed of rotating 3 angles (in degrees) around 3 vectors


void vConv3Rot3Vec3DMat3D(const float fRot0,
			  const float fRot1,
			  const float fRot2,
			  const float fVec0[3],
			  const float fVec1[3],
			  const float fVec2[3],
			  float *pfMat0,
              float *pfMatDeriv0,   // 1st derivative wrt fRot0
              float *pfMatDeriv1,   // 1st derivative wrt fRot1
              float *pfMatDeriv2    // 1st derivative wrt fRot2              
              )   // The resultant 3x3 matrix
{
  float fTmat0[3][3], fTmat1[3][3],  fTmat2[3][3], fTmat3[3][3];
  float fTDmat0[3][3], fTDmat1[3][3],  fTDmat2[3][3];
  float *pfTD0 = (float*)NULL ,*pfTD1 = (float*)NULL,*pfTD2 = (float*)NULL;
  float fTmatX[3][3],fTmatXX[3][3];

  if (pfMatDeriv0) {
      // Assume that if we are asking for one derivative, we want the other ones as well.
      pfTD0 = &fTDmat0[0][0];
      pfTD1 = &fTDmat1[0][0];
      pfTD2 = &fTDmat2[0][0];
  };

  vConvRotVec3DMat3D(fRot0, fVec0, (float*)fTmat0,pfTD0);
  vConvRotVec3DMat3D(fRot1, fVec1, (float*)fTmat1,pfTD1);
  vConvRotVec3DMat3D(fRot2, fVec2, (float*)fTmat2,pfTD2);

  if (pfMatDeriv0) {
      vMulMat3DMat3D(fTmat2,fTmat1,fTmatX);
      vMulMat3DMat3D(fTmatX,fTDmat0,fTmatXX);
      vCopyMat3D(&fTmatXX[0][0],pfMatDeriv0);
      vMulMat3DMat3D(fTmat2,fTDmat1,fTmatX);
      vMulMat3DMat3D(fTmatX,fTmat0,fTmatXX);
      vCopyMat3D(&fTmatXX[0][0],pfMatDeriv1);
      vMulMat3DMat3D(fTDmat2,fTmat1,fTmatX);
      vMulMat3DMat3D(fTmatX,fTmat0,fTmatXX);
      vCopyMat3D(&fTmatXX[0][0],pfMatDeriv2);
  };

  vMulMat3DMat3D(fTmat2, fTmat1, fTmat3);   // Is this correct order?
  vMulMat3DMat3D(fTmat3, fTmat0, fTmatX);

  float *pfTemp1, *pfTemp2;
  pfTemp1 = pfMat0;
  pfTemp2 = &fTmatX[0][0];

  for (int i = 0; i < 9; i ++)
    *pfTemp1++ = *pfTemp2++;
}


// This function tries to adjust the rotation parameters (a3fRots) so that they accomodate a product of rotations.
// The data pair (a3x3fRotVecs,a3fRots) define an initial rotation matrix, which shall be pre-multiplied, or post-multiplied by
// the new rotation matrix a3x3fDeltaRot.  We want to obtain the new values for the parameters (a3fRots) that will effect this and produce the correct product
// rotation matrix.  This also can be done using the nDeconvMat3DVec3DRot3() function, but there is not garuantee that
// the resulting rotation values will be near the intial ones.  This routine claims to obtain nearby values by using derivatives.
// Normally, only 3 or 4 passes are needed to reach convergence.

int nDeconvMat3DMat3D(double a3fRots[3],double a3x3fRotVecs[3][3],double a3x3fDeltaRot[3][3],bool bPreMultiply,double* pfResultMat) {
    double a3x3fPrevRot[3][3];
    double a3x3x3fPrevRotD[3][3][3];
    double a3x3fNormMat[3][3];
    double a3x3fNormMatInv[3][3];
    double a3fNormVec[3];
    double a3x3fTargetRot[3][3];
    double a3fDelta[3];

    int nPass;
    int nRow,nCol,nDeriv,nVar;
    bool bConverged;


    nPass = -1;
    do {
        nPass++;

        // Build the rotation matrix from the last iteration.
        // Also, discover the derivatives.
        vConv3Rot3Vec3DMat3D(a3fRots[0],a3fRots[1],a3fRots[2],
            a3x3fRotVecs[0],a3x3fRotVecs[1],a3x3fRotVecs[2],
            &a3x3fPrevRot[0][0],
            &a3x3x3fPrevRotD[0][0][0],
            &a3x3x3fPrevRotD[1][0][0],
            &a3x3x3fPrevRotD[2][0][0]);

        // If this is the first pass, then build the target matrix.
        if (nPass==0) {
            if (bPreMultiply)
                vMulMat3DMat3D(a3x3fDeltaRot,a3x3fPrevRot,a3x3fTargetRot);
            else
                vMulMat3DMat3D(a3x3fPrevRot,a3x3fDeltaRot,a3x3fTargetRot);
        };

        vZeroMat(3,3,&a3x3fNormMat[0][0]);
        vZeroMat(3,1,&a3fNormVec[0]);
        for (nRow=0;nRow<3;nRow++) {
            for (nCol=0;nCol<3;nCol++) {
                for (nDeriv=0;nDeriv<3;nDeriv++) {
                    for (nVar=0;nVar<3;nVar++) {
                        a3x3fNormMat[nVar][nDeriv] += a3x3x3fPrevRotD[nDeriv][nCol][nRow]*a3x3x3fPrevRotD[nVar][nCol][nRow];
                    };
                    a3fNormVec[nDeriv] += 
                        (a3x3fTargetRot[nCol][nRow] - a3x3fPrevRot[nCol][nRow])*
                        a3x3x3fPrevRotD[nDeriv][nCol][nRow];
                };

            };
        };
        fInvMat3D(&a3x3fNormMat[0][0],&a3x3fNormMatInv[0][0]);
        vMulMat3DVec3D(a3x3fNormMatInv,a3fNormVec,a3fDelta);
        vAddVec3DVec3D(a3fDelta,a3fRots,a3fRots);
        bConverged = (((fabs(a3fDelta[0])<0.0001) && (fabs(a3fDelta[0])<0.0001) && (fabs(a3fDelta[0])<0.0001)));

    } while ((!bConverged) && (nPass<50));

    if ((bConverged) && (pfResultMat)) {
        vCopyMat3D(&a3x3fTargetRot[0][0],pfResultMat);
    };

    return (int) (!bConverged);
};




//int nDeconvMat3DVec3DRot3(const float fMatrix[3][3],
int nDeconvMat3DVec3DRot3(float fMatrix[3][3],
			   const float* pfVec0,
			   const float* pfVec1,
			   const float* pfVec2,
			   float *pfRot0, float *pfRot1, float *pfRot2,
			   const float fSignRot1)
{
  // Take a 3x3 rotation matrix and decompose it into rotations around the
  // fVec0, fVec1, fVec2 axes.  This is the inverse of vConv3Rot3...
  // This code is derived from work of David Thomas, while he was at
  // MRC, Cambridge.
  // NEEDS REPAIR: when when vecs are XYZ and Y=+-90, then X and Z are ||
  //               and this algorithm fails.

  int   i, j;
  int   nStat;
  int   nTemp;
  float fSubConst, fDivConst, fZeroConst;
  float fTempMat0[3][3], fTempMat1[3][3], fTempMat2[3][3];
  float fInvTempMat01[3][3];
  float fTempVec0[3];
  float fB[2];
  float fTemp;

  // Setup bogus results in case of early return

  *pfRot0 = -999.0;
  *pfRot1 = -999.0;
  *pfRot2 = -999.0;

  nStat = nCalcGetGonConstants(pfVec0, pfVec1, pfVec2,
			       &fSubConst, &fDivConst, &fZeroConst);
  if (0 == nStat)
    {
      // Constants OK, now try to do the math
      // First try figure out what the second axis angle would be.

      fTemp = 0.0;
      for (j = 0; j < 3; j++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      fTemp = fTemp + pfVec2[i] * fMatrix[j][i] * pfVec0[j];
	    }
	}
      if (0.0 != fDivConst)
	{
	  fTemp = (fTemp - fSubConst) / fDivConst;
	}
      else
	{
	  return (-2);  // Failure
	}
      if (1.0001 < fTemp)
	{
	  return (-2);
	}
      else if (-1.0001 > fTemp)
	return (-2);
      else
	{
	  if (1.0  < fTemp) fTemp =  1.0;
	  if (-1.0 > fTemp) fTemp = -1.0;
	  fTemp = acos((double)fTemp);
	}

      // fTemp contains the angle difference, now put the requested sign on it
      // and add in the offset

      fTemp = fZeroConst + fTemp * fSignRot1;
      nTemp =  (int) (fTemp / (2.0 * Gs_dPI)); // These two lines: fTemp mod 2pi
      fTemp = fTemp - (float) nTemp * 2.0 * Gs_dPI;

      *pfRot1  = fTemp / Gs_dRADIANS_PER_DEGREE;

      // Get rotation matrix for the second angle

      vConvRotVec3DMat3D(*pfRot1, pfVec1, (float*)fTempMat1);

      // Rotate 3rd axis by this rotation matrix

      vMulMat3DVec3D(fTempMat1, pfVec0, fTempVec0);

      // Evaluate matrix which when operated on (1,cos(Rot2),sin(Rot2)) will
      // be equal to Rot0

      vCross3D(pfVec2, fTempVec0, &fTempMat2[0][0]);       // 1st column
      fTemp = fDot3D(pfVec2, fTempVec0);
      for (i = 0; i < 3; i++)
	{
	  fTempMat2[2][i] = pfVec2[i] * fTemp;               // 3rd column
	  fTempMat2[1][i] = fTempVec0[i] - fTempMat2[2][i];  // 2nd column
	}
      // Rotate pfVec0 through required rotation

      vMulMat3DVec3D(fMatrix, pfVec0, fTempVec0);

      // Now work out the normal vector

      fB[0] = fDot3D(fTempVec0, &fTempMat2[0][0]);
      fB[1] = fDot3D(fTempVec0, &fTempMat2[1][0]);

      if (fB[0]*fB[0] + fB[1]*fB[1] < 1.0E-20)
	return (-2);

      *pfRot2 = atan2((double)fB[0], (double)fB[1]) / Gs_dRADIANS_PER_DEGREE;

      // Get rotation matrix for the third angle

      vConvRotVec3DMat3D(*pfRot2, pfVec2, (float*)fTempMat2);

      // Combine with rotation matrix for second angle

      vMulMat3DMat3D(fTempMat2, fTempMat1, fTempMat0);
      fTemp = fInvMat3D((float*) fTempMat0, (float*)fInvTempMat01);
      if (0.0 == fTemp) return (-2);
      vMulMat3DMat3D(fInvTempMat01, fMatrix, fTempMat2);

      // Work out Rot2 injector matrix (in fTempMat1)
      //                   T
      //  = pfVec2 * pfVec2

      for (i = 0; i < 3; i++)
	{
	  vMulVec3DScalar(pfVec0, pfVec0[i], &fTempMat1[i][0]);
	  // Convert Rot2 injector matrix to minus Rot2 projector matrix
	  fTempMat1[i][i] = fTempMat1[i][i] - 1.0f;
	}

      // Dot product of the minus Rot2 projector matrix with the altered matrix
      fB[1] = -fDotND(9, &fTempMat1[0][0], &fTempMat2[0][0]);

      // Dot product of the Rot2-cross matrix with the altered matrix

      fB[0] =  pfVec0[0] * (fTempMat2[1][2] - fTempMat2[2][1])
             + pfVec0[1] * (fTempMat2[2][0] - fTempMat2[0][2])
             + pfVec0[2] * (fTempMat2[0][1] - fTempMat2[1][0]);

      if (fB[0]*fB[0] + fB[1]*fB[1] < 1.0E-20)
	return (-3);

      // Finally solve for the first axis angle
      *pfRot0 = atan2((double)fB[0], (double)fB[1]) / Gs_dRADIANS_PER_DEGREE;
    }

  return (nStat);
}

//void vDeconvMat3D3XYZ(const float fMatrix[3][3],
void vDeconvMat3D3XYZ(float fMatrix[3][3],
		      float *pfRotX, float *pfRotY, float *pfRotZ)
{
  //     Extract the missetting angles from fMatrix
  //     If X, Y, Z are rotation angles around the x,y,z axes, then fMatrix =
  //
  //    cosYcosZ        sinXsinYcosZ - cosXsinZ     cosXsinYcosZ - sinXsinZ
  //    cosYsinZ        sinXsinYsinZ + cosXcosZ     cosXsinYsinZ - sinXcosZ
  //     -sinY           sinXcosY                    cosXcosY
  //
  //     Watch out for special cases when any angle is near 0, 90, 180, 270!

  float fSinY, fCosY;
  fSinY = -fMatrix[0][2];

  //  vListMat3D(&fMatrix[0][0]);

  if (1.0 < fSinY)
    fSinY = 1.0;
  else if (-1.0 > fSinY)
    fSinY = -1.0;

  fCosY = 1.0f  -  fSinY * fSinY;

  if (fCosY < 1.0E-6)
    {
      // Special case, rotation in Y is either -90 or +90 degrees.
      // This means X and Z are coincident (parallel or anti-parallel).
      // Let rotation in X be 0.0, and the missetting all in Z.  For this
      // case:  fMatrix[1][1] = cosZ, and fMatrix[0][1] = -sinZ
      // jwp: actually we should let input rotX stay the same and put
      //      all the rotation in Z!

      *pfRotZ = atan2((double)-fMatrix[1][0], (double)fMatrix[1][1]) / Gs_dRADIANS_PER_DEGREE;
      *pfRotY = 90.0;
      if (0.0 > fSinY) *pfRotY = -90.0;

      if ( (*pfRotX >= -360.0) && (*pfRotX <= +360.0) )
	{
	  // input rotX is believable, so leave it alone and adjust rotZ
	  *pfRotZ = *pfRotZ + (*pfRotX * fSinY);
	}
      else
	{
	  *pfRotX = 0.0;
	}
    }
  else
    {
      //  Let CY be positive so that -90 < Y < +90 degrees.

      float fRotX, fRotY, fRotZ;
      fCosY  = sqrt((double)fCosY);

      // Try to keep results close to inputs, (i.e. in same octant) if
      // you think inputs are legit.

      fRotX = atan2((double)fMatrix[1][2], (double)fMatrix[2][2]) / Gs_dRADIANS_PER_DEGREE;
      fRotY = atan2((double)fSinY, (double)fCosY)                 / Gs_dRADIANS_PER_DEGREE;
      fRotZ = atan2((double)fMatrix[0][1], (double)fMatrix[0][0]) / Gs_dRADIANS_PER_DEGREE;

      if ( (*pfRotX >= -360.0) && (*pfRotX <= +360.0) )
	{
	  *pfRotX = fRotX;
	}
      else
	{
	  *pfRotX = fRotX;
	}
      if ( (*pfRotY >= -360.0) && (*pfRotY <= +360.0) )
	{
	  *pfRotY = fRotY;
	}
      else
	{
	  *pfRotY = fRotY;
	}
      if ( (*pfRotZ >= -360.0) && (*pfRotZ <= +360.0) )
	{
	  *pfRotZ = fRotZ;
	}
      else
	{
	  *pfRotZ = fRotZ;
	}
    }
}

void vZeroMat3D(float fMatrix[9])
{
  // Set a 3x3 matrix to all zeroes

  int i;
  for (i = 0; i < 9; i++) fMatrix[i] = 0.0;
}

void vZeroMat(const int nDim0, const int nDim1, float fMatrix[])
{
  int i;
  for (i = 0; i < nDim0*nDim1; i++) fMatrix[i] = 0.0;
}


void vIdentMat3D(float fMatrix[3][3])
{
  // Set a 3x3 matrix to the identity matrix

  int i;
  vZeroMat3D(&fMatrix[0][0]);
  for (i = 0; i < 3; i++) fMatrix[i][i] = 1.0;

}

int nCalcGetGonConstants(const float *pfVec1,  // Input rotation vector 1
                         const float *pfVec2,  // Input rotation vector 2
                         const float *pfVec3,  // Input rotation vector 3
                         float *pfSubConst,
                         float *pfDivConst,
                         float *pfZeroConst)
{
  // Compute 3 constants needed by nDeconvMat3DVec3DRot3.

  float fVec1dotVec2;
  float fVec2dotVec3;
  float fVec1dotVec3;
  float fVec2crossVec3[3];
  float fDet;

  // Compute the triple scalar product Vec1 . (Vec2 x Vec3)

  vCross3D(pfVec2, pfVec3, fVec2crossVec3);
  fDet = fDot3D(pfVec1, fVec2crossVec3);

  // Compute some dot products

  fVec1dotVec2 = fDot3D(pfVec1, pfVec2);
  fVec2dotVec3 = fDot3D(pfVec2, pfVec3);
  fVec1dotVec3 = fDot3D(pfVec1, pfVec3);

  *pfSubConst  = fVec1dotVec3 * fVec2dotVec3;
  fVec1dotVec3 = fVec1dotVec3 - *pfSubConst;  // lhs is temporary variable now
  *pfDivConst = (float)sqrt((double)(fDet * fDet + fVec1dotVec3 * fVec1dotVec3));
  if (*pfDivConst > 1e-15) {
    *pfZeroConst = atan2((double)fDet, (double)fVec1dotVec3);
    return (0);
  }
  else
    return (1);
}

float fSign(const float fVar)
{
  // Return the sign of fVar as +1.0 or -1.0

  if (fVar >= 0.0)
    return (1.0);
  else
    return (-1.0);
}

void vBuildUpperRightND(const int nDim, float *pfMat)
{
  // Copy the lower left triangle of a nDim by nDim matrix to the upper right
  // This makes the matrix symmetric and equal to its tranpose

  int   i, j, k, l;
  for (i = 0; i < nDim; i++)
    {
      for (j = i+1; j < nDim; j++)
	{
	  k = i * nDim + j;
	  l = j * nDim + i;
	  pfMat[l] = pfMat[k];
	}
    }
}

void vBuildBasis3D(float a3fVec[3],float a3fVecBasis[3][3])\
{
    int nx;
    int nLargest,nSmallest,nOther;

    // Find the largest and smallest values.  These will get swapped.
    nLargest=0;
    nSmallest=1;
    for (nx=0;nx<3;nx++) {
        if (a3fVec[nx]>a3fVec[nLargest])
            nLargest = nx;
        if (a3fVec[nx]<a3fVec[nSmallest])
            nSmallest = nx;
    };
    nOther = (-nLargest -nSmallest + 6) % 3;
    if ((nLargest == nSmallest) ||
        (nLargest == nOther) ||
        (nSmallest == nOther)) {
        printf("Error in vBuildBasis3D()\n");
        return;
    };

    vCopyVec3D(a3fVec,a3fVecBasis[0]);
    a3fVecBasis[1][nSmallest] = a3fVecBasis[0][nLargest];
    a3fVecBasis[1][nLargest] = - a3fVecBasis[0][nSmallest];    
    a3fVecBasis[1][nOther] = 0.0;
    vCross3D(a3fVecBasis[0],a3fVecBasis[1],a3fVecBasis[2]);
    return;
};

void
vAffineTransform(const int n0Start,
                 const int n0End,
                 const int n1Start,
                 const int n1End,
                 const int n2Start,
                 const int n2End,
                 const float *pfMatrix,
                 const float *pfTrans,
                 float *pfResult)
{
  //  Transform one 3D grid to another 3D grid by an affine transformation.
  //  That is:
  //   G2  =   M I1 + T    where M is a 3x3 matrix and T is a 3D vector
  //   --      = --   -          =                     -
  //  The integer indices of the elements of grid 1 are found in I1 and the
  //  result is grid 2 in G2.
  //  Indices are numbered from 0 as in the C language (not from 1 as in
  //  Fortran).
  //
  //  Do it with minimal multiplication for speed.
  // This routine is derived from work of Gerard Bricogne.

  // Does matrix come in row-major or column-major?
  //
  // Set up constants

  float fM00     = pfMatrix[0];
  float fM01     = pfMatrix[1];
  float fM02     = pfMatrix[2];
  float fM10     = pfMatrix[3];
  float fM11     = pfMatrix[4];
  float fM12     = pfMatrix[5];
  float fM20     = pfMatrix[6];
  float fM21     = pfMatrix[7];
  float fM22     = pfMatrix[8];
  float fT0      = pfTrans[0];
  float fT1      = pfTrans[1];
  float fT2      = pfTrans[2];

  float f00Start = fM00 * (float) n0Start;
  float f01Start = fM01 * (float) n0Start;
  float f02Start = fM02 * (float) n0Start;

  float f10Start = fM10 * (float) n1Start;
  float f11Start = fM11 * (float) n1Start;
  float f12Start = fM12 * (float) n1Start;

  float f20Start = fM20 * (float) n2Start;
  float f21Start = fM21 * (float) n2Start;
  float f22Start = fM22 * (float) n2Start;

  float fT00, fT01, fT02, fT10, fT11, fT12, fT20, fT21, fT22;
  float f0Sum, f1Sum, f2Sum;

  fT20           = f20Start + fT0;
  fT21           = f21Start + fT1;
  fT22           = f22Start + fT2;

  float *pfTemp;

  pfTemp = pfResult;

  int i, j, k;   // Loop counters

  for (k = n2Start; k <= n2End; k++)
    {
      fT20 = fT20 + fM20;
      fT21 = fT21 + fM21;
      fT22 = fT22 + fM22;

      fT10 = f10Start;
      fT11 = f11Start;
      fT12 = f12Start;

      for (j = n1Start; j <= n1End; j++)
        {

          fT10 = fT10 + fM10;
          fT11 = fT11 + fM11;
          fT12 = fT12 + fM12;

          f0Sum = fT20 + fT10;
          f1Sum = fT21 + fT11;
          f2Sum = fT22 + fT12;

          fT00 = f00Start;
          fT01 = f01Start;
          fT02 = f02Start;

          for (i = n0Start; i <= n0End; i++)
            {
              fT00      = fT00 + fM00;
              fT01      = fT01 + fM01;
              fT02      = fT02 + fM02;
              *pfTemp++ = f0Sum + fT00;
              *pfTemp++ = f1Sum + fT01;
              *pfTemp++ = f2Sum + fT02;
            } // end i
        }  // end j
    }  // end k
};

                 

// Converted from old Fortran.  In this C++ version, please remember
// that the first element of an array has an index of 0.


void
vEigen(const int nDim, const int nMode, float *pfA, float *pfR)

// Routine vEIGEN
//
//        Purpose
//           Compute eigenvalues and eigenvectors of a real symmetric
//           matrix
//
//        Usage
//           int nDim, nMode;
//           float *pfA, *pfR;
//           vEigen(nDim, nMode, pfA, pfR)
//
//        Description of parameters
//         pfA - original matrix (symmetric), destroyed in computation.
//               resultant eigenvalues are developed in diagonal of
//               matrix pfA in descending order.
//  P.S. The equation for the indices of pfA is A[j][i] -> A[i+(j*j+j)/2])
//         pfRr - resultant matrix of eigenvectors (stored columnwise,
//                in same sequence as eigenvalues)
//         nDim - order of matrices pfA and pfR
//         nMode- input code
//                   0   Compute eigenvalues and eigenvectors
//                   1   Compute eigenvalues only (pfR need not be
//                       dimensioned but must still appear in calling
//                       sequence)
//
//        Remarks
//           original matrix pfA must be real symmetric (storage mode=1)
//           matrix pfA cannot be in the same location as matrix pfR
//
//        Subroutines and function subprograms required
//           sqrtf, abs
//
//        Method
//           Diagonalization method originated by Jacobi and adapted
//           by von Neumann for large computers as found in 'Mathematical
//           Methods for Digital Computers', edited by A. Ralston and
//           H.S. Wilf, John Wiley and Sons, New York, 1962, Chapter 7
//          See also Numerical Recipes, Chapter 11, subroutine JACOBI.
//        If a double precision version of this routine is desired, the
//        following variables should be declared double:
//
//     double precision a,r,anorm,anrmx,thr,x,y,sinx,sinx2,cosx,
//    1                 cosx2,sincs,range
//
//        The double precision version of this subroutine must also
//        contain double precision functions.  sqrtf in statements
//        below must be changed to sqrt.  abs below
//        must be changed to absd. The fRANGE constant should
//        be changed to 1.0d-12.

{

  int   i, j, k, l, m;
  int   nIQ, nJQ, nLL, nMM, nLM, nILR, nIMR, nILQ, nIMQ, nIL, nIM, nLQ, nMQ;
  int   nInd;
  float fX, fY, fSinX, fCosX, fSinX2, fCosX2, fSinCS, fTHR, fANORM, fANRMX;
  float *pfTemp;

  float fRANGE = (float)1.0E-6;

  if (pfA == pfR)
    cout << "WARNING in vEigen, pfA == pfR!\n" << flush;

  // Generate identity matrix in pfR

  pfTemp = pfR;
  if (0 == nMode)
    {
      for (j = 0; j < nDim; j++)
	{
	  for (i = 0; i < nDim; i++)
	    {
	      if (i == j)
		*pfTemp++ = 1.0;
	      else
		*pfTemp++ = 0.0;
	    }
	}
    }

  // Compute initial and final norms (anorm and anormx)

  fANORM = 0.0;
  for (i = 0; i < nDim; i++)
    {
      for (j = i; j < nDim; j++)
	{
	  if (i != j)
	    {
	      nIQ = i + (j * j  + j) / 2;
	      fANORM = fANORM + pfA[nIQ] * pfA[nIQ];
	    }
	}
    }
  if (fANORM > 0.0)
    {
      fANORM = 1.41421356f * sqrt((double)fANORM);
      fANRMX = fANORM * fRANGE / (float) nDim;

      // Initialize indicators and compute threshold, fThr

      fTHR = fANORM;

      do
	{                                   // until (fTHR <= fANRMX)
	  fTHR = fTHR / (float) nDim;
	  do
	    {                                 // until (nInd == 0){
	      nInd = 0;
	      for (l = 0; l < nDim - 1; l++)
		{   // from start to next to last col
		  for (m = l + 1; m < nDim; m++)
		    { // from l+1 to last col

		      // Compute sin and cos

		      nMQ = (m * m  + m) / 2;        //     [m][m]
		      nLQ = (l * l  + l) / 2;        //     [l][l]
		      nLM = l + nMQ;
		      if ( fabs(pfA[nLM]) >= fTHR)
			{
			  nInd = 1;         // Set flag that we did some calcs
			  nLL  = l + nLQ;
			  nMM  = m + nMQ;
			  fX = 0.5f * (pfA[nLL] - pfA[nMM]);
			  fY = -pfA[nLM] / sqrt((double)(pfA[nLM] * pfA[nLM]
						 +  fX * fX));
			  if (fX < 0.0) fY = -fY;
			  fSinX  = fY / sqrt(2.0 * (1.0 +
						     ( sqrt((double)(1.0 - fY * fY)))));
			  fSinX2 = fSinX * fSinX;
			  fCosX  = sqrt(1.0 - (double)fSinX2);
			  fCosX2 = fCosX * fCosX;
			  fSinCS = fSinX * fCosX;

			  // Rotate l and m columns
			
			  nILQ = nDim * l;
			  nIMQ = nDim * m;
			  for (i = 0; i < nDim; i++)
			    {
			      nIQ = (i * i  +  i) / 2;
			      if (i != l)
				{
				  if (i != m)
				    {
				      if (i < m)
					nIM = i + nMQ;
				      else
					nIM = m + nIQ;
				      if (i < l)
					nIL = i + nLQ;
				      else
					nIL = l + nIQ;
				      fX       = pfA[nIL] * fCosX
					         -  pfA[nIM] * fSinX;
				      pfA[nIM] = pfA[nIL] * fSinX
					         +  pfA[nIM] * fCosX;
				      pfA[nIL] = fX;
				    }
				}
			      if (0 == nMode)
				{
				  nILR      = nILQ + i;
				  nIMR      = nIMQ + i;
				  fX        = pfR[nILR] * fCosX
				              -  pfR[nIMR] * fSinX;
				  pfR[nIMR] = pfR[nILR] * fSinX
				              +  pfR[nIMR] * fCosX;
				  pfR[nILR] = fX;
				}
			    } // end of i loop
			
			  fX       = 2.0f * pfA[nLM] * fSinCS;
			  fY       =  pfA[nLL] * fCosX2  +  pfA[nMM] * fSinX2
			              -  fX;
			  fX       =  pfA[nLL] * fSinX2  +  pfA[nMM] * fCosX2
			              +  fX;
			  pfA[nLM] = (pfA[nLL] - pfA[nMM]) * fSinCS  +
			              pfA[nLM] * (fCosX2 - fSinX2);
			  pfA[nLL] =  fY;
			  pfA[nMM] =  fX;
			} // end if
		    } // end m loop
		} // end l loop
	    } while (1 == nInd);     // Repeat until inner loops not entered

	  //   Compare threshold with final norm

	} while (fTHR > fANRMX);

      //  Sort eigenvalues and eigenvectors

      nIQ = -nDim;
      for (i = 0; i < nDim; i++)
	{
	  nIQ = nIQ + nDim;
	  nLL  = i + (i * i + i) / 2;          // element [i][i]
	  nJQ  = nDim * (i - 1);
	  for (j = i; j < nDim; j++)
	    {
	      nJQ = nJQ + nDim;
	      nMM = j + (j * j + j) / 2;         // element [j][j]
	      if ( pfA[nLL] < pfA[nMM])
		{
		  // Swap diagonal elements A[nLL} and A[nMM]

		  fX      = pfA[nLL];
		  pfA[nLL] = pfA[nMM];
		  pfA[nMM] = fX;
		  if (0 == nMode)
		    {

		      // Swap off-diagonal elements

		      for (k = 0; k < nDim; k++)
			{
			  nILR      = nIQ + k;
			  nIMR      = nJQ + k;
			  fX        = pfR[nILR];
			  pfR[nILR] = pfR[nIMR];
			  pfR[nIMR] = fX;
			}
		    }
		}
	    } // end j
	} // end i
    }
}

void
vTRED2(const int nDim, float **ppfMat,
		   float *pfDiag, float *pfOffDiag)
{
  //  Householder reduction of a real, symmetric matrix in
  //	**pfMat[1..nDim][1..nDim] (NOTE FORTRAN-style indexing!)
  //  Matrix **pfMat is destroyed in the computation and on return
  //  holds the orthogonal Q matrix effecting the transformation.
  //  pfDiag[1..nDim] holds the diagonal elements
  //  pfOffDiag[1] is set to 0
  //  pfOffDiag[2..n] holds the off-diagonal elements
  //
  // This code is derived from algorithms published in
  // Numerical Recipes in C: The Art of Scientific Programming
  // (c) 1988, Cambridge University Press.
  // See pp. 373ff.

  int    i, j, k, L;          // Loop counters
                              // (L uppercase so not confused with 1
  float  fScale, fHH, fH, fG, fF;

  for (i = nDim; i >=2; i--)
    {
      L = i - 1;
      fH = fScale = 0.0;
      if (L > 1)
        {
          for (k = 1; k <= L; k++)
            {
              fScale += fabs(ppfMat[i][k]);
            }
          if (fScale == 0.0)
            {
              // Skip transformation
              pfOffDiag[i] = ppfMat[i][L];
            }
          else
            {
              for (k = 1; k <= L; k++)
                {
                  ppfMat[i][k] /= fScale;               // Use scale a's f tra
                  fH += ppfMat[i][k] * ppfMat[i][k];  // Form sigma in h
                }
              fF             = ppfMat[i][L];
              fG             = fF > 0.0 ? -sqrt((double)fH) : sqrt((double)fH);
              pfOffDiag[i]   = fScale * fG;
              fH            -= fF * fG;
              ppfMat[i][L] = fF - fG;
              fF             = 0.0;
              for (j = 1; j <= L; j++)
                {
                  ppfMat[j][i] = ppfMat[i][j] / fH;  // Store u/H in ith col
                  fG             = 0.0;
                  for (k = 1; k <= j; k++)
                    {
                      fG        += ppfMat[j][k] * ppfMat[i][k];
                    }
                  for (k = j+1; k <= L; k++)
                    {
                      fG        += ppfMat[k][j] * ppfMat[i][k];
                    }
                  pfOffDiag[j]  = fG / fH;// Form element of p in temp unused el
                  fF      += pfOffDiag[j] * ppfMat[i][j];

                } // end j loop
              fHH = fF / (fH + fH);
              for (j = 1; j <= L; j++)
                {
                  fF           = ppfMat[i][j];
                  pfOffDiag[j] = fG = pfOffDiag[j] - fHH * fF;
                  for (k = 1; k <= j; k++)
                    {
                      ppfMat[j][k] -= (fF * pfOffDiag[k] +  fG * ppfMat[i][k]);
                    }
                } // end j loop
            }
        }
      else
        {
          pfOffDiag[i] = ppfMat[i][L];
        }
      pfDiag[i] = fH;
    }
  pfDiag[1]    = 0.0;
  pfOffDiag[1] = 0.0;
  for (i = 1; i <= nDim; i++)
    {
      L = i - 1;
      if (pfDiag[i])
        {               // This skipped when i=0 jwp: I do not like this! :(
          for (j = 1; j <= L; j++)
            {
              fG = 0.0;
              for (k = 1; k <= L; k++)
                {
                  fG += ppfMat[i][k] * ppfMat[k][j];
                }
              for (k = 1; k <= L; k++)
                {
                  ppfMat[k][j] -= fG * ppfMat[k][i];
                }
            }
        }
      pfDiag[i]    = ppfMat[i][i];
      ppfMat[i][i] = 1.0;
      for (j = 1; j <= L; j++)
        {
          ppfMat[j][i] = ppfMat[i][j] = 0.0;
        }
    }
}

void
vTQLI(const int nDim, float *pfDiag, float *pfOffDiag,
      float **ppfEigenVecs)

{
  // QL algorithm with implicit shifts, to determine the eigenvalues and
  // eigenvectors of a real, symmetric matrix previously reduced by vTRED2.
  // On input pfDiag[1..nDim] contains the diagonal elements of the tridiagonal
  // matrix.  On output, it returns the eigenvalues.  The vector
  // pfOffDiag[1..nDim] inputs the subdiagonal elements of the tridiagonal
  // matrix with pfOffDiag[1] arbitrary.  On output pfOffDiag is destroyed.  The
  // matrix ppfEigenVecs is input as the identity matrix OR if the eigenvectors
  // of a matrix that has been reduced by vTRED2 are required, then
  // ppfEigenVecs is the matrix output by vTRED2.  In either case, the kth
  // column of ppfEigenVecs returns the normalized eigenvector corresponding
  // to ppfDiag[k].
  //
  // This code is derived from algorithms published in
  // Numerical Recipes in C: The Art of Scientific Programming
  // (c) 1988, Cambridge University Press.
  // See pp. 380ff.

  int m, L, iter, i, k;   // Loop counters  lowercase l is too easily confused
                          // with digit 1
  float fS, fR, fP, fG, fF, fDD, fC, fB;

  for (i = 2; i <= nDim; i++)
    {
      pfOffDiag[i-1] = pfOffDiag[i]; // Convenient to renumber array offdiag
    }
  pfOffDiag[nDim] = 0.0;

  for (L = 1; L <= nDim; L++)
    {
      iter = 0;
      do
	{
	  for (m = L; m <= nDim-1; m++)
	    {
	      fDD = fabs(pfDiag[m]) + fabs(pfDiag[m+1]);

	      // Look for a single small subdiagonal element to split the matrix

	      if ( (float) (fabs(pfOffDiag[m]) + fDD) == fDD) break;
	    }

	  // If there was no break, then m = n hopefully

	  if (m != L)
	    {
	      if (iter++ == 30)
		{
		  cout << "Too many iterations in vTQLI!\n";
		  return;
		}
	      fG = (pfDiag[L+1] - pfDiag[L]) / (2.0f * pfOffDiag[L]);
	      fR = sqrt((double) (fG * fG) + 1.0);
	      fG = pfDiag[m] - pfDiag[L] + pfOffDiag[L] / (fG + SIGNF(fR, fG));
	      fS = fC = 1.0;
	      fP = 0.0;
	      for (i = m-1; i >= L; i--)
		{
		  fF = fS * pfOffDiag[i];
		  fB = fC * pfOffDiag[i];
		  if (fabs(fF) >= fabs(fG))
		    {
		      fC              = fG / fF;
		      fR              = sqrt((double) (fC * fC) + 1.0);
		      pfOffDiag[i+1] = fF * fR;
		      fC             *= (fS = 1.0f / fR);
		    }
		  else
		    {
		      fS              = fF / fG;
		      fR              = sqrt( (double)(fS * fS) + 1.0);
		      pfOffDiag[i+1]  = fG * fR;
		      fS             *= (fC = 1.0f / fR);
		    }
		  fG           =  pfDiag[i+1] - fP;
		  fR           = (pfDiag[i]   - fG) * fS  +  2.0f * fC * fB;
		  fP           = fS * fR;
		  pfDiag[i+1]  = fG + fP;
		  fG           = fC * fR  - fB;
		  for (k = 1; k <= nDim; k++)
		    {
		      fF                   = ppfEigenVecs[k][i+1];
		      ppfEigenVecs[k][i+1] = fS * ppfEigenVecs[k][i] + fC * fF;
		      ppfEigenVecs[k][i]   = fC * ppfEigenVecs[k][i] - fS * fF;
		    }
		}
	      pfDiag[L]    = pfDiag[L] - fP;
	      pfOffDiag[L] = fG;
	      pfOffDiag[m] = 0.0;
	    } // endif  m != L
	} while (m != L);
    } // end l loop
}

void
vEigsrt(const int nDim, float *pfEigval, float **ppfEigvec)
{
  // Given the eigenvalues pfEigval[1...nDim] and eigenvectors
  // ppfEigvec[1...nDim][1...nDim] as output from vTQLI, this routine
  // sorts the eigenvalues into ASCENDING order, and rearranges the
  // columns of ppfEigvec accordingly.
  // From Numerical Recipes in C Section 11.1 p. 366.

  int i, j, k;
  float fTemp;

  for (i = 1; i <= nDim; i++)
    {
      fTemp = pfEigval[k=i];
      for (j = i + 1; j <= nDim; j++)
	{
	  if (pfEigval[j] >= fTemp) fTemp = pfEigval[k=j];     // Save largest
	}
      if (k != i)
	{
	  pfEigval[k] = pfEigval[i];                  // Need to swap
	  pfEigval[i] = fTemp;
	  for (j = 1; j <= nDim; j++)
	    {
	      fTemp           = ppfEigvec[j][i];
	      ppfEigvec[j][i] = ppfEigvec[j][k];
	      ppfEigvec[j][k] = fTemp;
	    }
	}
    }
}

int
nMatrixCompress(const int nNumRowsCols, const int *pnDelFlags,
		float *pfMatrix, float *pfVector)
{

  // Delete rows and columns of a square matrix
  //       pfMatrix[nNumRowsCols][nNumRowsCols]
  // according to values of pnDelFlags.  If pnDelFlags[i] is non-zero,
  // then delete row and column i.
  // Return the number of rows and columns
  // (nNum) left in pfMatrix which is now pfMatrix[nNum][nNum];
  // Delete elements pfVector, so that it matches ...

  int i, j;             // Loop counters
  int nNum;
  float *pfTempMat;
  float *pfTempVec;

  nNum = 0;
  pfTempMat = pfMatrix;
  pfTempVec = pfVector;

  for (i = 0; i < nNumRowsCols; i++)
    {
      if (0 == pnDelFlags[i])
	{
	  nNum++;
	  *pfTempVec++ = pfVector[i];   // This element stays in the vector
	}
      for (j = 0; j < nNumRowsCols; j++)
	{
	  if ( (0 == pnDelFlags[i]) && (0 == pnDelFlags[j]) )
	    {
	      // This element stays in matrix

	      *pfTempMat++ = pfMatrix[i * nNumRowsCols + j];
	    }
	}
    }
  return (nNum);
}

int nint(const float fTemp)
{
  if (0.0 < fTemp)
    return ( (int) (fTemp + 0.5));
  else
    return ( (int) (fTemp - 0.5));
}


#define SWAPFFT(a, b) tempr=(a); (a)=(b); (b)=tempr;

int
nFFT1D(const int nSize, const int nSign, float *pfData)
{
    // 1D Fast Fourier Transform routine from "Numerical Recipes in C", pp. 411-
    // nSize MUST be an integer power of 2
    // No checks are made on the input arguments.
    
    
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;
    int nn;
    int isign;

    // The array convention here is Fortran not C.  Decrement the array pointer to compensate.
    pfData--;

    
    isign = nSign;
    nn = nSize;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i+=2)
    {
        if (j > i)
        {
            SWAPFFT(pfData[j],   pfData[i]);
            SWAPFFT(pfData[j+1], pfData[i+1]);
        }
        m = n >> 1;
        while (m >=2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j+= m;
    }
    mmax = 2;
    while (n > mmax)
    {
        istep = 2 * mmax;
        theta = 6.28318530717959 / (double)(isign * mmax);
        wtemp = sin(0.5 * theta);
        wpr   = -2.0 * wtemp * wtemp;
        wpi   = sin(theta);
        wr    = 1.0;
        wi    = 0.0;
        for (m = 1; m < mmax; m+=2)
        {
            for (i = m; i <=n; i+=istep)
            {
                j          = i + mmax;
                tempr      = (float)(wr * pfData[j]    -  wi * pfData[j+1]);
                tempi      = (float)(wr * pfData[j+1]  +  wi * pfData[j]);
                pfData[j]    = pfData[i]   - tempr;
                pfData[j+1]  = pfData[i+1] - tempi;
                pfData[i]   += tempr;
                pfData[i+1] += tempi;
            }
            wr = (wtemp=wr) * wpr  -  wi * wpi  + wr;
            wi = wi * wpr  + wtemp * wpi  + wi;
        }
        mmax = istep;
    }
    return (0);
}






double fDot3D (const double fVec1[3], const double fVec2[3])
{
  return (fVec1[0] * fVec2[0]  +  fVec1[1] * fVec2[1]  + fVec1[2] * fVec2[2]);
}

double fDotND (const int nDim, const double fVec1[], const double fVec2[])
{
  double fSum = 0.0;
  for (int i = 0; i < nDim; i++) fSum = fSum + (fVec1[i] * fVec2[i]);
  return (fSum);
}

void vCross3D (const double fVec1[3], const double fVec2[3], double fResult[3])
{
  fResult[0] = fVec1[1] * fVec2[2]  - fVec1[2] * fVec2[1];
  fResult[1] = fVec1[2] * fVec2[0]  - fVec1[0] * fVec2[2];
  fResult[2] = fVec1[0] * fVec2[1]  - fVec1[1] * fVec2[0];
}

//void vMulMat3DMat3D(const double fMat1[3][3], const double fMat2[3][3],
void vMulMat3DMat3D(double fMat1[3][3], double fMat2[3][3],
                    double fMat3[3][3])
{
  int i, j, k;
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  fMat3[j][i] = 0.0;
	  for (k = 0; k < 3; k++)
	    {
	      fMat3[j][i] = fMat3[j][i]  + fMat1[k][i] * fMat2[j][k];
	    }
	}
    }
}

void vMulMatNDMatND(const int nDim, const double *pfMat1, const double *pfMat2,
                    double *pfMat3)
{
  // Multiply an nDim x nDim matrix by another matrix
  // The matrix is stored as in Fortran, so Mat[column][row]
  // pfMat3 must be different from pfMat1 and pfMat2;

  int i, j, k;
  const double *pfTempMat1, *pfTempMat2;
  double *pfTempMat3;

  pfTempMat3 = pfMat3;

  for (i = 0; i < nDim; i++)
    {
      pfTempMat3 = pfMat3 + i;
      for (j = 0; j < nDim; j++)
	{
	  *pfTempMat3 = 0.0;
	  pfTempMat2 = pfMat2 + (j * nDim);
	  pfTempMat1 = pfMat1 + i;
	  for (k = 0; k < nDim; k++)
	    {
	      *pfTempMat3 = *pfTempMat3  + *pfTempMat1 * (*pfTempMat2++);
	      pfTempMat1  = pfTempMat1 + nDim;
	    }
	  pfTempMat3 = pfTempMat3 + nDim;
	}
    }
}


void vMulMatMDNDMatMDKD(const int nDimM, const int nDimN, const int nDimK, const double* pfMat1,
                        const double* pfMat2, double* pfMat3) {
  // Multiply an nDimM x nDimN matrix by another matrix with dimensions nDimN x nDimK
  // The matrix is stored as in Fortran, so Mat[column][row]
  // pfMat3 must be different from pfMat1 and pfMat2;

  int i, j, k;
  const double *pfTempMat1, *pfTempMat2;
  double *pfTempMat3;

  pfTempMat3 = pfMat3;

  for (i = 0; i < nDimM; i++)
    {
      pfTempMat3 = pfMat3 + i;
      for (j = 0; j < nDimK; j++)
	{
	  *pfTempMat3 = 0.0;
	  pfTempMat2 = pfMat2 + (j * nDimN);
	  pfTempMat1 = pfMat1 + i;
	  for (k = 0; k < nDimN; k++)
	    {
	      *pfTempMat3 = *pfTempMat3  + *pfTempMat1 * (*pfTempMat2++);
	      pfTempMat1  = pfTempMat1 + nDimM;
	    }
	  pfTempMat3 = pfTempMat3 + nDimM;
	}
    }
};



//void vMulMat3DVec3D(const double fMat[3][3], const double fVec1[3],
void vMulMat3DVec3D(double fMat[3][3], const double *pfVec1,
                    double *pfVec2)
{
  // pfVec1 and pfVec2 cannot be the same!
  int i, j;
  if (pfVec1 == pfVec2)
    cout << "WARNING in vMulMat3DVec3D, pfVec1 == pfVec2!\n" << flush;    
  for (i = 0; i < 3; i++)
    {
      pfVec2[i] = 0.0;
      for (j = 0; j < 3; j++)
	{
	  pfVec2[i] = pfVec2[i] + (fMat[j][i] * pfVec1[j]);
	}
    }
}

void vMulMatNDVecND(const int nDim, const double *pfMat, const double *pfVec1,
                    double *pfVec2)
{
  // Multiply an nDim x nDim matrix by nDim vector.
  // The matrix is stored as in Fortran, so Mat[column][row]
  // pfVec1 and pfVec2 cannot be the same!

  const double *pfTempMat;
  const double *pfTempVec1;
  double       *pfTempVec2;
  int i, j;

  if (pfVec1 == pfVec2)
    cout << "WARNING in vMulMatNDVecND, pfVec1 == pfVec2!\n" << flush;    

  pfTempVec2 = pfVec2;
  for (i = 0; i < nDim; i++)
    {
      *pfTempVec2 = 0.0;
      pfTempVec1  = pfVec1;
      pfTempMat   = pfMat + i;
      for (j = 0; j < nDim; j++)
	{
	  *pfTempVec2 = *pfTempVec2 + (*pfTempMat * *pfTempVec1++);
	  pfTempMat  = pfTempMat + nDim;
	}
      pfTempVec2++;
    }
}

void vMulMat3DScalar(const double fMat1[9], const double fScalar, double fMat2[9])
{
  for (int i = 0; i < 9; i++)
    fMat2[i] = fMat1[i] * fScalar;
}

void vMulVec3DScalar(const double fVec1[3], const double fScalar, double fVec2[3])
{
  for (int i = 0; i < 3; i++)
    fVec2[i] = fVec1[i] * fScalar;
}

void vMulVecNDScalar(const int nDim, const double *pfVec1, const double fScalar,
		     double *pfVec2)
{
  for (int i = 0; i < nDim; i++)
    pfVec2[i] = pfVec1[i] * fScalar;
}

double fInvMat4D(const double *pfMat1, double *pfMat2) {

    return 0.0;
};

double fInvMat3D(const double *pfMat1, double *pfMat2)
{
  // Invert a 3x3 matrix.  Return the determinant.
  // Mat1 and Mat2 must be different matrices

  double fDet;

  fDet  =   pfMat1[0] * (pfMat1[4] * pfMat1[8]  -  pfMat1[7] * pfMat1[5])
          - pfMat1[1] * (pfMat1[3] * pfMat1[8]  -  pfMat1[5] * pfMat1[6])
          + pfMat1[2] * (pfMat1[3] * pfMat1[7]  -  pfMat1[4] * pfMat1[6]);
  if ( (0.0 != fDet) && (NULL != pfMat2) )
    {
      pfMat2[0] =  (pfMat1[4] * pfMat1[8]  -  (pfMat1[5] * pfMat1[7])) / fDet;
      pfMat2[1] = -(pfMat1[1] * pfMat1[8]  -  (pfMat1[2] * pfMat1[7])) / fDet;
      pfMat2[2] =  (pfMat1[1] * pfMat1[5]  -  (pfMat1[2] * pfMat1[4])) / fDet;
      pfMat2[3] = -(pfMat1[3] * pfMat1[8]  -  (pfMat1[5] * pfMat1[6])) / fDet;
      pfMat2[4] =  (pfMat1[0] * pfMat1[8]  -  (pfMat1[2] * pfMat1[6])) / fDet;
      pfMat2[5] = -(pfMat1[0] * pfMat1[5]  -  (pfMat1[2] * pfMat1[3])) / fDet;
      pfMat2[6] =  (pfMat1[3] * pfMat1[7]  -  (pfMat1[4] * pfMat1[6])) / fDet;
      pfMat2[7] = -(pfMat1[0] * pfMat1[7]  -  (pfMat1[1] * pfMat1[6])) / fDet;
      pfMat2[8] =  (pfMat1[0] * pfMat1[4]  -  (pfMat1[1] * pfMat1[3])) / fDet;
    }
  return (fDet);
}

double fInvMat2D(const double *pfMat1, double *pfMat2)
{
    double fDet;

    fDet = pfMat1[0]*pfMat1[3] - pfMat1[1]*pfMat1[2];

    if ( (0.0 != fDet) && (NULL != pfMat2))
    {
        pfMat2[0] = pfMat1[3] / fDet;
        pfMat2[3] = pfMat1[0] / fDet;
        pfMat2[1] = - pfMat1[1] / fDet;
        pfMat2[2] = - pfMat1[2] / fDet;
    };
    return (fDet);
};


double fInvMatND(const int nDim, const double *pfMat1, double *pfMat2)
{
  // Invert a nDim x nDim matrix by Gauss-Jordan elimination.
  // We assume the matrix is small and we don't care about using
  // extra storage.  See "Numerical Recipes" ....
  // Matrices are stored as in Fortran  A[column][row] so be careful.
  // pfMat1 and pfMat2 may be the same storage locations matrices
  // Returns (1.) if all Ok, (0.0) if not OK.

  int i, j, k, l, ll;
  double fTemp;
//  double fMat[nDim][nDim];
//  int   nPiv[nDim], nIndexRow[nDim], nIndexCol[nDim];  // Bookkeeping
  int   *nPiv, *nIndexRow, *nIndexCol;  // Bookkeeping
  int   nRow, nCol;
  double fBig;

  double *fMatLinear,*fptr1,*fptr2;

  fMatLinear = new double[nDim*nDim];
  nPiv       = new int[nDim];
  nIndexRow  = new int[nDim];
  nIndexCol  = new int[nDim];

  // Make a copy of the input matrix

//  vCopyVecND(nDim * nDim, pfMat1, (double *) fMat);
  vCopyVecND(nDim * nDim, pfMat1, fMatLinear);

  // "Zero" the bookkeeping pivot array

  for (i = 0; i < nDim; i++)
    nPiv[i] = 0;

  for (i = 0; i < nDim; i++)
    {
      fBig = 0.0;
      for (j = 0; j < nDim; j++)
	{
	  if (nPiv[j] != 1)
	    {
	      fptr1 = fMatLinear+j;
	      for (k = 0; k < nDim; k++)
		{
		  if (nPiv[k] == 0)
		    {
		      // fTemp = fMat[k][j];
		      fTemp = *(fptr1+k*nDim);
		      if (fTemp < 0.0) fTemp = -fTemp; // abs function
		      if (fTemp >= fBig)
			{
			  fBig = fTemp;
			  nRow = j;
			  nCol = k;
			}
		    }
		  else if (nPiv[k] > 1)
		    {
		      cout << "Singular matrix" << endl;
		      return (0.0);
		    }
		} // end k loop
	    }
	} // end j loop

      nPiv[nCol]++;

      // We now have the pivot element ...

      if (nRow != nCol)
	{
	  fptr1 = fMatLinear+nRow;
	  fptr2 = fMatLinear+nCol;
	  for (l = 0; l < nDim; l++, fptr1 += nDim, fptr2 += nDim)
	    {
	      //	fTemp         = fMat[l][nRow];
	      //	fMat[l][nRow] = fMat[l][nCol];
	      //	fMat[l][nCol] = fTemp;
	      fTemp  = *fptr1;
	      *fptr1 = *fptr2;
	      *fptr2 = fTemp;
	    } // end l loop
	}

      nIndexRow[i] = nRow;
      nIndexCol[i] = nCol;

      //    if (fMat[nCol][nCol] == 0.0) {
      fptr1 = fMatLinear+nCol*(nDim+1);
      if(0.0 == *fptr1)
	{
	  cout << "Singular matrix" << endl;
	  return (0.0);
	}
      //    fTemp = 1.0 / fMat[nCol][nCol];
      //    fMat[nCol][nCol] = 1.0;
      fTemp = 1.0/(*fptr1);
      *fptr1 = 1.0;
      fptr2 = fMatLinear+nCol;
      for (l = 0; l < nDim; l++, fptr2 += nDim)
	{
	  //      fMat[l][nCol] = fMat[l][nCol] * fTemp;
	  *fptr2 *= fTemp;
	}
      for (ll = 0; ll < nDim; ll++)
	{
	  fptr1 = fMatLinear+ll;
	  fptr2 = fptr1+nCol*nDim;
	  if (ll != nCol)
	    {
	      //	fTemp          = fMat[nCol][ll];
	      //	fMat[nCol][ll] = 0.0;
	      fTemp = *fptr2;
	      *fptr2 = 0.0;
	      fptr2 = fMatLinear+nCol;
	      for (l = 0; l < nDim; l++, fptr1 += nDim, fptr2 += nDim)
		{
		  //	  fMat[l][ll] = fMat[l][ll]  -  fMat[l][nCol] * fTemp;
		  *fptr1 -= (*fptr2)*fTemp;
		}
	    }
	} // end ll loop
    } // end i loop (main loop over columns of the reduction)
  for (l = nDim-1; l >= 0; l--)
    {
      if (nIndexRow[l] != nIndexCol[l])
	{
	  fptr1 = fMatLinear+nDim*nIndexRow[l];
	  fptr2 = fMatLinear+nDim*nIndexCol[i];
	  for (k = 0; k < nDim; k++, fptr1++, fptr2++)
	    {
	      //	fTemp                 = fMat[nIndexRow[l]][k];
	      //	fMat[nIndexRow[l]][k] = fMat[nIndexCol[l]][k];
	      //	fMat[nIndexCol[l]][k] = fTemp;
	      fTemp = *fptr1;
	      *fptr1 = *fptr2;
	      *fptr2 = fTemp;
	    }
	}
    }
//  vCopyVecND(nDim * nDim, (double *) fMat, pfMat2);
  vCopyVecND(nDim * nDim, fMatLinear, pfMat2);

  delete [] fMatLinear;

  //  delete [] fMat;
  delete [] nPiv;
  delete [] nIndexRow;
  delete [] nIndexCol;

  return (1.0);
}

void vAddVec3DVec3D(const double fVec1[3], const double fVec2[3], double fVec3[3])
{
  for (int i = 0; i < 3; i++)	
     fVec3[i] = fVec1[i] + fVec2[i];
}

void vAddVecNDVecND(const int nDim, const double fVec1[], const double fVec2[],
		    double fVec3[])
{
  for (int i = 0; i < nDim; i++)	
     fVec3[i] = fVec1[i] + fVec2[i];
}

void vSubVec3DVec3D(const double fVec1[3], const double fVec2[3], double fVec3[3])
{
  for (int i = 0; i < 3; i++)	
     fVec3[i] = fVec1[i] - fVec2[i];
}

void vSubVecNDVecND(const int nDim, const double fVec1[], const double fVec2[],
		    double fVec3[])
{
  for (int i = 0; i < nDim; i++)	
     fVec3[i] = fVec1[i] - fVec2[i];
}

void vTranMat3D(double fMat[3][3])
{
  // Transpose a 3x3 matrix in place

  double fTemp;
  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = i+1; j < 3; j++) {
      fTemp      = fMat[i][j];
      fMat[i][j] = fMat[j][i];
      fMat[j][i] = fTemp;
    }
  }
}

void vTranMatND(const int nDim, double *pfMat)
{
  // Transpose a nDim by nDim matrix in place

  double fTemp;
  int i, j, k, m;
  for (i = 0; i < nDim; i++)
    {
      for (j = i+1; j < nDim; j++)
	{
	  m = i * nDim + j;             // m is [i][j]
	  k = j * nDim + i;             // k is [j][i]
	  fTemp = pfMat[k];
	  pfMat[k] = pfMat[m];
	  pfMat[m] = fTemp;
	}
    }
}

void vTranMatMDND(const int nDimM, 
                  const int nDimN, 
                  const double* pfMat1,
                  double* pfMat2) 
{
  // Transpose an nDimM x nDimN matrix in place
  int i, j;

  for (i = 0; i < nDimM; i++)
      for (j = 0; j< nDimN; j++)
          pfMat2[ j + nDimN * i ] = pfMat1[ i + nDimM * j ];
} 


void vCopyMat3D(int* pnMat1,int* pnMat2) 
{
    for (int i = 0; i < 9; i++) pnMat2[i] = pnMat1[i];
}
void vCopyMat3D(int* pnMat1,double* pfMat2) 
{
    for (int i = 0; i < 9; i++) pfMat2[i] = pnMat1[i];
}

void vCopyMat3D(const double *pfMat1, double *pfMat2)
{
  for (int i = 0; i < 9; i++) pfMat2[i] = pfMat1[i];
}

void vCopyVec3D(const double *pfVec1, double *pfVec2)
{
  for (int i = 0; i < 3; i++) pfVec2[i] = pfVec1[i];
}


void vCopyVecND(const int nDim, const double *pfVec1, double *pfVec2)
{
  for (int i = 0; i < nDim; i++) pfVec2[i] = pfVec1[i];
}


double fNormVec3D(double fVec[3])
{
  double fTemp;
  fTemp = fDot3D(fVec, fVec);
  if (fTemp > 0.0) {
    //fTemp = sqrtf(fTemp);  // BUG! since sqrtf is for float
    fTemp = sqrt(fTemp);
    vMulVec3DScalar(fVec, 1.0/fTemp, fVec);
  }
  return (fTemp);
}

double fLenVec3D(const double fVec[3])
{
  double fTemp;
  fTemp = fDot3D(fVec, fVec);
  fTemp = sqrt(fTemp);
  return (fTemp);
}

double fLenVecND(const int nDim, const double *pfVec)
{
  double fTemp;
  fTemp = fDotND(nDim, pfVec, pfVec);
  fTemp = sqrt(fTemp);
  return (fTemp);
}

void vNormMat3D(double a3x3fMat[3][3])
{
    double a3fSub[3];

    // gram-schmid (spell?) the matrix.
    fNormVec3D(a3x3fMat[0]);
    vMulVec3DScalar(a3x3fMat[0],fDot3D(a3x3fMat[0],a3x3fMat[1]),a3fSub);
    vSubVec3DVec3D(a3x3fMat[1],a3fSub,a3x3fMat[1]);
    fNormVec3D(a3x3fMat[1]);
    vMulVec3DScalar(a3x3fMat[0],fDot3D(a3x3fMat[0],a3x3fMat[2]),a3fSub);
    vSubVec3DVec3D(a3x3fMat[2],a3fSub,a3x3fMat[2]);
    vMulVec3DScalar(a3x3fMat[1],fDot3D(a3x3fMat[1],a3x3fMat[2]),a3fSub);
    vSubVec3DVec3D(a3x3fMat[2],a3fSub,a3x3fMat[2]);
    fNormVec3D(a3x3fMat[2]);
    return;
};

void vListVec3D(double fVec[3])
{
 cout << "Vector is: " << fVec[0] << ", " << fVec[1] << ", " << fVec[2] << endl;
}

void vListMat3D(double *pfMat)
{
  double *pfTemp;
  pfTemp = pfMat;
  cout << "Matrix is:" << endl;
  for (int i = 0; i < 3; i++)
    {
      cout << "[" << i << "][*]: " << pfTemp[0] << ", "
           << pfTemp[1] << ", " << pfTemp[2] << endl;
      pfTemp += 3;
    }
}

void vListMatMN(const int nDim0, const int nDim1, double *pfMat)
{
  // List an MxN matrix, where M is fastest changing address

  int i, j;
  double *pfTemp;
  pfTemp = pfMat;
  cout << "Matrix is:" << endl;
  for (i = 0; i < nDim1; i++)
    {
      cout << "[" << i << "][*]: ";
      for (j = 0; j < nDim0-1; j++)
	{
	  cout << *pfTemp++ << ", ";
	}
      cout << *pfTemp++ << '\n';
    }
  cout << endl;
}

// Convert a rotation (in degrees) around a vector to a matrix and
// possibly its derivatives

void vConvRotVec3DMat3D(const double fRot, const double *pfVec,
                        double *pfMat0,   // The resultant 3x3 matrix 
                        double *pfMat1,   // 1st derivative of *pfMat0 (w.r.t. fRot)
                        double *pfMat2,   // 2nd derivative of *pfMat0 (w.r.t. pfVec[0])
                        double *pfMat3,   // 3rd derivative of *pfMat0 (w.r.t. pfVec[1])
                        double *pfMat4)   // 4th derivative of *pfMat0 (w.r.t. pfVec[2])
{
    
    
    double fCos, fSin, fOneMinusCos;
    double fTemp;
    
    fTemp        = fRot * Gs_dRADIANS_PER_DEGREE;
    fCos         = cos(fTemp);
    fSin         = sin(fTemp);
    fOneMinusCos = 1.0 - fCos;
    
    // Careful -- row, column not checked out!
    
    pfMat0[0] = fOneMinusCos * pfVec[0] * pfVec[0]  + fCos;
    pfMat0[1] = fOneMinusCos * pfVec[1] * pfVec[0]  + fSin * pfVec[2];
    pfMat0[2] = fOneMinusCos * pfVec[2] * pfVec[0]  - fSin * pfVec[1];
    pfMat0[3] = fOneMinusCos * pfVec[0] * pfVec[1]  - fSin * pfVec[2];
    pfMat0[4] = fOneMinusCos * pfVec[1] * pfVec[1]  + fCos;
    pfMat0[5] = fOneMinusCos * pfVec[2] * pfVec[1]  + fSin * pfVec[0];
    pfMat0[6] = fOneMinusCos * pfVec[0] * pfVec[2]  + fSin * pfVec[1];
    pfMat0[7] = fOneMinusCos * pfVec[1] * pfVec[2]  - fSin * pfVec[0];
    pfMat0[8] = fOneMinusCos * pfVec[2] * pfVec[2]  + fCos;
    
    //  cout << "Rot ang: " << fRot << " ";
    //  vListMat3D(pfMat0);
    if (NULL != pfMat1)
    {
        // First derivative requested
        
        pfMat1[0] = fSin * pfVec[0] * pfVec[0]  - fSin;
        pfMat1[1] = fSin * pfVec[1] * pfVec[0]  + fCos * pfVec[2];
        pfMat1[2] = fSin * pfVec[2] * pfVec[0]  - fCos * pfVec[1];
        pfMat1[3] = fSin * pfVec[0] * pfVec[1]  - fCos * pfVec[2];
        pfMat1[4] = fSin * pfVec[1] * pfVec[1]  - fSin;
        pfMat1[5] = fSin * pfVec[2] * pfVec[1]  + fCos * pfVec[0];
        pfMat1[6] = fSin * pfVec[0] * pfVec[2]  + fCos * pfVec[1];
        pfMat1[7] = fSin * pfVec[1] * pfVec[2]  - fCos * pfVec[0];
        pfMat1[8] = fSin * pfVec[2] * pfVec[2]  - fSin;
        vMulVecNDScalar(9, pfMat1, Gs_dRADIANS_PER_DEGREE, pfMat1);
    };
    if (NULL != pfMat2)
    {       
        // First derivative requested (w.r.t. );
        pfMat2[0] = 2.0 * fOneMinusCos * pfVec[0];
        pfMat2[1] = fOneMinusCos * pfVec[1];
        pfMat2[2] = fOneMinusCos * pfVec[2];
        pfMat2[3] = fOneMinusCos * pfVec[1];
        pfMat2[4] = 0.0;
        pfMat2[5] = fSin;
        pfMat2[6] = fOneMinusCos * pfVec[2];
        pfMat2[7] = - fSin;
        pfMat2[8] = 0.0;
    };
    if (NULL != pfMat3)
    {    
        // First derivative requested (w.r.t. );
        pfMat3[0] = 0.0;
        pfMat3[1] = fOneMinusCos * pfVec[0];
        pfMat3[2] = - fSin;
        pfMat3[3] = fOneMinusCos * pfVec[0];
        pfMat3[4] = 2.0 * fOneMinusCos * pfVec[1];
        pfMat3[5] = fOneMinusCos * pfVec[2];
        pfMat3[6] = fSin;
        pfMat3[7] = fOneMinusCos * pfVec[2];
        pfMat3[8] = 0.0;
    };
    if (NULL != pfMat4)
    {
        // First derivative requested (w.r.t. );
        pfMat4[0] = 0.0;
        pfMat4[1] = fSin;
        pfMat4[2] = fOneMinusCos * pfVec[0];
        pfMat4[3] = - fSin;
        pfMat4[4] = 0.0;
        pfMat4[5] = fOneMinusCos * pfVec[1];
        pfMat4[6] = fOneMinusCos * pfVec[0] ;
        pfMat4[7] = fOneMinusCos * pfVec[1];
        pfMat4[8] = 2.0 * fOneMinusCos * pfVec[2] ;
    };
}

// Returns a matrix composed of rotating 3 angles (in degrees) around 3 vectors

// Returns a matrix composed of rotating 3 angles (in degrees) around 3 vectors

void vConv3Rot3Vec3DMat3D(const double fRot0,
			  const double fRot1,
			  const double fRot2,
			  const double fVec0[3],
			  const double fVec1[3],
			  const double fVec2[3],
			  double *pfMat0,
              double *pfMatDeriv0,   // 1st derivative wrt fRot0
              double *pfMatDeriv1,   // 1st derivative wrt fRot1
              double *pfMatDeriv2    // 1st derivative wrt fRot2              
              )   // The resultant 3x3 matrix
{
  double fTmat0[3][3], fTmat1[3][3],  fTmat2[3][3], fTmat3[3][3];
  double fTDmat0[3][3], fTDmat1[3][3],  fTDmat2[3][3];
  double *pfTD0 = (double*)NULL ,*pfTD1 = (double*)NULL,*pfTD2 = (double*)NULL;
  double fTmatX[3][3],fTmatXX[3][3];

  if (pfMatDeriv0) {
      // Assume that if we are asking for one derivative, we want the other ones as well.
      pfTD0 = &fTDmat0[0][0];
      pfTD1 = &fTDmat1[0][0];
      pfTD2 = &fTDmat2[0][0];
  };

  vConvRotVec3DMat3D(fRot0, fVec0, (double*)fTmat0,pfTD0);
  vConvRotVec3DMat3D(fRot1, fVec1, (double*)fTmat1,pfTD1);
  vConvRotVec3DMat3D(fRot2, fVec2, (double*)fTmat2,pfTD2);

  if (pfMatDeriv0) {
      vMulMat3DMat3D(fTmat2,fTmat1,fTmatX);
      vMulMat3DMat3D(fTmatX,fTDmat0,fTmatXX);
      vCopyMat3D(&fTmatXX[0][0],pfMatDeriv0);
      vMulMat3DMat3D(fTmat2,fTDmat1,fTmatX);
      vMulMat3DMat3D(fTmatX,fTmat0,fTmatXX);
      vCopyMat3D(&fTmatXX[0][0],pfMatDeriv1);
      vMulMat3DMat3D(fTDmat2,fTmat1,fTmatX);
      vMulMat3DMat3D(fTmatX,fTmat0,fTmatXX);
      vCopyMat3D(&fTmatXX[0][0],pfMatDeriv2);
  };

  vMulMat3DMat3D(fTmat2, fTmat1, fTmat3);   // Is this correct order?
  vMulMat3DMat3D(fTmat3, fTmat0, fTmatX);

  double *pfTemp1, *pfTemp2;
  pfTemp1 = pfMat0;
  pfTemp2 = &fTmatX[0][0];

  for (int i = 0; i < 9; i ++)
    *pfTemp1++ = *pfTemp2++;
}



//int nDeconvMat3DVec3DRot3(const double fMatrix[3][3],
int nDeconvMat3DVec3DRot3(double fMatrix[3][3],
			   const double* pfVec0,
			   const double* pfVec1,
			   const double* pfVec2,
			   double *pfRot0, double *pfRot1, double *pfRot2,
			   const double fSignRot1)
{
  // Take a 3x3 rotation matrix and decompose it into rotations around the
  // fVec0, fVec1, fVec2 axes.  This is the inverse of vConv3Rot3...
  // This code is derived from work of David Thomas, while he was at
  // MRC, Cambridge.
  // NEEDS REPAIR: when when vecs are XYZ and Y=+-90, then X and Z are ||
  //               and this algorithm fails.

  int   i, j;
  int   nStat;
  int   nTemp;
  double fSubConst, fDivConst, fZeroConst;
  double fTempMat0[3][3], fTempMat1[3][3], fTempMat2[3][3];
  double fInvTempMat01[3][3];
  double fTempVec0[3];
  double fB[2];
  double fTemp;

  // Setup bogus results in case of early return

  *pfRot0 = -999.0;
  *pfRot1 = -999.0;
  *pfRot2 = -999.0;

  nStat = nCalcGetGonConstants(pfVec0, pfVec1, pfVec2,
			       &fSubConst, &fDivConst, &fZeroConst);
  if (0 == nStat)
    {
      // Constants OK, now try to do the math
      // First try figure out what the second axis angle would be.

      fTemp = 0.0;
      for (j = 0; j < 3; j++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      fTemp = fTemp + pfVec2[i] * fMatrix[j][i] * pfVec0[j];
	    }
	}
      if (0.0 != fDivConst)
	{
	  fTemp = (fTemp - fSubConst) / fDivConst;
	}
      else
	{
	  return (-2);  // Failure
	}
      if (1.0001 < fTemp)
	{
	  return (-2);
	}
      else if (-1.0001 > fTemp)
	return (-2);
      else
	{
	  if (1.0  < fTemp) fTemp =  1.0;
	  if (-1.0 > fTemp) fTemp = -1.0;
	  fTemp = acos(fTemp);
	}

      // fTemp contains the angle difference, now put the requested sign on it
      // and add in the offset

      fTemp = fZeroConst + fTemp * fSignRot1;
      nTemp =  (int) (fTemp / (2.0 * Gs_dPI)); // These two lines: fTemp mod 2pi
      fTemp = fTemp - (double) nTemp * 2.0 * Gs_dPI;

      *pfRot1  = fTemp / Gs_dRADIANS_PER_DEGREE;

      // Get rotation matrix for the second angle

      vConvRotVec3DMat3D(*pfRot1, pfVec1, (double*)fTempMat1);

      // Rotate 3rd axis by this rotation matrix

      vMulMat3DVec3D(fTempMat1, pfVec0, fTempVec0);

      // Evaluate matrix which when operated on (1,cos(Rot2),sin(Rot2)) will
      // be equal to Rot0

      vCross3D(pfVec2, fTempVec0, &fTempMat2[0][0]);       // 1st column
      fTemp = fDot3D(pfVec2, fTempVec0);
      for (i = 0; i < 3; i++)
	{
	  fTempMat2[2][i] = pfVec2[i] * fTemp;               // 3rd column
	  fTempMat2[1][i] = fTempVec0[i] - fTempMat2[2][i];  // 2nd column
	}
      // Rotate pfVec0 through required rotation

      vMulMat3DVec3D(fMatrix, pfVec0, fTempVec0);

      // Now work out the normal vector

      fB[0] = fDot3D(fTempVec0, &fTempMat2[0][0]);
      fB[1] = fDot3D(fTempVec0, &fTempMat2[1][0]);

      if (fB[0]*fB[0] + fB[1]*fB[1] < 1.0E-20)
	return (-2);

      *pfRot2 = atan2((double)fB[0], (double)fB[1]) / Gs_dRADIANS_PER_DEGREE;

      // Get rotation matrix for the third angle

      vConvRotVec3DMat3D(*pfRot2, pfVec2, (double*)fTempMat2);

      // Combine with rotation matrix for second angle

      vMulMat3DMat3D(fTempMat2, fTempMat1, fTempMat0);
      fTemp = fInvMat3D((double*) fTempMat0, (double*)fInvTempMat01);
      if (0.0 == fTemp) return (-2);
      vMulMat3DMat3D(fInvTempMat01, fMatrix, fTempMat2);

      // Work out Rot2 injector matrix (in fTempMat1)
      //                   T
      //  = pfVec2 * pfVec2

      for (i = 0; i < 3; i++)
	{
	  vMulVec3DScalar(pfVec0, pfVec0[i], &fTempMat1[i][0]);
	  // Convert Rot2 injector matrix to minus Rot2 projector matrix
	  fTempMat1[i][i] = fTempMat1[i][i] - 1.0;
	}

      // Dot product of the minus Rot2 projector matrix with the altered matrix
      fB[1] = -fDotND(9, &fTempMat1[0][0], &fTempMat2[0][0]);

      // Dot product of the Rot2-cross matrix with the altered matrix

      fB[0] =  pfVec0[0] * (fTempMat2[1][2] - fTempMat2[2][1])
             + pfVec0[1] * (fTempMat2[2][0] - fTempMat2[0][2])
             + pfVec0[2] * (fTempMat2[0][1] - fTempMat2[1][0]);

      if (fB[0]*fB[0] + fB[1]*fB[1] < 1.0E-20)
	return (-3);

      // Finally solve for the first axis angle
      *pfRot0 = atan2((double)fB[0], (double)fB[1]) / Gs_dRADIANS_PER_DEGREE;
    }

  /*
  // Double check result
  double a3x3Result[3][3];
  vConv3Rot3Vec3DMat3D((double)*pfRot0, (double)*pfRot1, (double)*pfRot2, 
		       pfVec0, pfVec1, pfVec2, (double*)&a3x3Result[0][0]);
  cout << "Result matrix is: ";
  vListMat3D(&a3x3Result[0][0]);
  */
  return (nStat);
}

//void vDeconvMat3D3XYZ(const double fMatrix[3][3],
void vDeconvMat3D3XYZ(double fMatrix[3][3],
		      double *pfRotX, double *pfRotY, double *pfRotZ)
{
  //mrp fMatrix = (RotZ)(RotY)(RotX)
  //     Extract the missetting angles from fMatrix
  //     If X, Y, Z are rotation angles around the x,y,z axes, then fMatrix =
  //
  //    cosYcosZ        sinXsinYcosZ - cosXsinZ     cosXsinYcosZ - sinXsinZ
  //    cosYsinZ        sinXsinYsinZ + cosXcosZ     cosXsinYsinZ - sinXcosZ
  //     -sinY           sinXcosY                    cosXcosY
  //
  //     Watch out for special cases when any angle is near 0, 90, 180, 270!

  double fSinY, fCosY;
  fSinY = -fMatrix[0][2];

  //  vListMat3D(&fMatrix[0][0]);

  if (1.0 < fSinY)
    fSinY = 1.0;
  else if (-1.0 > fSinY)
    fSinY = -1.0;

  fCosY = 1.0  -  fSinY * fSinY;

  if (fCosY < 1.0E-6)
    {
      // Special case, rotation in Y is either -90 or +90 degrees.
      // This means X and Z are coincident (parallel or anti-parallel).
      // Let rotation in X be 0.0, and the missetting all in Z.  For this
      // case:  fMatrix[1][1] = cosZ, and fMatrix[0][1] = -sinZ
      // jwp: actually we should let input rotX stay the same and put
      //      all the rotation in Z!

      *pfRotZ = atan2((double)-fMatrix[1][0], (double)fMatrix[1][1]) / Gs_dRADIANS_PER_DEGREE;
      *pfRotY = 90.0;
      if (0.0 > fSinY) *pfRotY = -90.0;

      if ( (*pfRotX >= -360.0) && (*pfRotX <= +360.0) )
	{
	  // input rotX is believable, so leave it alone and adjust rotZ
	  *pfRotZ = *pfRotZ + (*pfRotX * fSinY);
	}
      else
	{
	  *pfRotX = 0.0;
	}
    }
  else
    {
      //  Let CY be positive so that -90 < Y < +90 degrees.

      double fRotX, fRotY, fRotZ;
      fCosY  = sqrt(fCosY);

      // Try to keep results close to inputs, (i.e. in same octant) if
      // you think inputs are legit.

      fRotX = atan2((double)fMatrix[1][2], (double)fMatrix[2][2]) / Gs_dRADIANS_PER_DEGREE;
      fRotY = atan2((double)fSinY, (double)fCosY)                 / Gs_dRADIANS_PER_DEGREE;
      fRotZ = atan2((double)fMatrix[0][1], (double)fMatrix[0][0]) / Gs_dRADIANS_PER_DEGREE;

      if ( (*pfRotX >= -360.0) && (*pfRotX <= +360.0) )
	{
	  *pfRotX = fRotX;
	}
      else
	{
	  *pfRotX = fRotX;
	}
      if ( (*pfRotY >= -360.0) && (*pfRotY <= +360.0) )
	{
	  *pfRotY = fRotY;
	}
      else
	{
	  *pfRotY = fRotY;
	}
      if ( (*pfRotZ >= -360.0) && (*pfRotZ <= +360.0) )
	{
	  *pfRotZ = fRotZ;
	}
      else
	{
	  *pfRotZ = fRotZ;
	}
    }
}

void vZeroMat3D(double fMatrix[9])
{
  // Set a 3x3 matrix to all zeroes

  int i;
  for (i = 0; i < 9; i++) fMatrix[i] = 0.0;
}

void vZeroMat(const int nDim0, const int nDim1, double fMatrix[])
{
  int i;
  for (i = 0; i < nDim0*nDim1; i++) fMatrix[i] = 0.0;
}


void vIdentMat3D(double fMatrix[3][3])
{
  // Set a 3x3 matrix to the identity matrix

  int i;
  vZeroMat3D(&fMatrix[0][0]);
  for (i = 0; i < 3; i++) fMatrix[i][i] = 1.0;

}

int nCalcGetGonConstants(const double *pfVec1,  // Input rotation vector 1
                         const double *pfVec2,  // Input rotation vector 2
                         const double *pfVec3,  // Input rotation vector 3
                         double *pfSubConst,
                         double *pfDivConst,
                         double *pfZeroConst)
{
  // Compute 3 constants needed by nDeconvMat3DVec3DRot3.

  double fVec1dotVec2;
  double fVec2dotVec3;
  double fVec1dotVec3;
  double fVec2crossVec3[3];
  double fDet;

  // Compute the triple scalar product Vec1 . (Vec2 x Vec3)

  vCross3D(pfVec2, pfVec3, fVec2crossVec3);
  fDet = fDot3D(pfVec1, fVec2crossVec3);

  // Compute some dot products

  fVec1dotVec2 = fDot3D(pfVec1, pfVec2);
  fVec2dotVec3 = fDot3D(pfVec2, pfVec3);
  fVec1dotVec3 = fDot3D(pfVec1, pfVec3);
  *pfSubConst  = fVec1dotVec2 * fVec2dotVec3;
  fVec1dotVec3 = fVec1dotVec3 - *pfSubConst;  // lhs is temporary variable now
  *pfDivConst = sqrt(fDet * fDet + fVec1dotVec3 * fVec1dotVec3);
  if (*pfDivConst > 1e-15) {
    *pfZeroConst = atan2((double)fDet, (double)fVec1dotVec3);
    return (0);
  }
  else
    return (1);
}

double fSign(const double fVar)
{
  // Return the sign of fVar as +1.0 or -1.0

  if (fVar >= 0.0)
    return (1.0);
  else
    return (-1.0);
}

void vBuildUpperRightND(const int nDim, double *pfMat)
{
  // Copy the lower left triangle of a nDim by nDim matrix to the upper right
  // This makes the matrix symmetric and equal to its tranpose

  int   i, j, k, l;
  for (i = 0; i < nDim; i++)
    {
      for (j = i+1; j < nDim; j++)
	{
	  k = i * nDim + j;
	  l = j * nDim + i;
	  pfMat[l] = pfMat[k];
	}
    }
}

void vBuildBasis3D(double a3fVec[3],double a3fVecBasis[3][3])
{
    int nx;
    int nLargest,nSmallest,nOther;

    // Find the largest and smallest values.  These will get swapped.
    nLargest=0;
    nSmallest=1;
    for (nx=0;nx<3;nx++) {
        if (a3fVec[nx]>a3fVec[nLargest])
            nLargest = nx;
        if (a3fVec[nx]<a3fVec[nSmallest])
            nSmallest = nx;
    };
    nOther = (-nLargest -nSmallest + 6) % 3;
    if ((nLargest == nSmallest) ||
        (nLargest == nOther) ||
        (nSmallest == nOther)) {
        printf("Error in vBuildBasis3D()\n");
        return;
    };

    vCopyVec3D(a3fVec,a3fVecBasis[0]);
    a3fVecBasis[1][nSmallest] = a3fVecBasis[0][nLargest];
    a3fVecBasis[1][nLargest] = - a3fVecBasis[0][nSmallest];    
    a3fVecBasis[1][nOther] = 0.0;
    vCross3D(a3fVecBasis[0],a3fVecBasis[1],a3fVecBasis[2]);
    return;
};

void
vAffineTransform(const int n0Start,
                 const int n0End,
                 const int n1Start,
                 const int n1End,
                 const int n2Start,
                 const int n2End,
                 const double *pfMatrix,
                 const double *pfTrans,
                 double *pfResult)
{
  //  Transform one 3D grid to another 3D grid by an affine transformation.
  //  That is:
  //   G2  =   M I1 + T    where M is a 3x3 matrix and T is a 3D vector
  //   --      = --   -          =                     -
  //  The integer indices of the elements of grid 1 are found in I1 and the
  //  result is grid 2 in G2.
  //  Indices are numbered from 0 as in the C language (not from 1 as in
  //  Fortran).
  //
  //  Do it with minimal multiplication for speed.
  // This routine is derived from work of Gerard Bricogne.

  // Does matrix come in row-major or column-major?
  //
  // Set up constants

  double fM00     = pfMatrix[0];
  double fM01     = pfMatrix[1];
  double fM02     = pfMatrix[2];
  double fM10     = pfMatrix[3];
  double fM11     = pfMatrix[4];
  double fM12     = pfMatrix[5];
  double fM20     = pfMatrix[6];
  double fM21     = pfMatrix[7];
  double fM22     = pfMatrix[8];
  double fT0      = pfTrans[0];
  double fT1      = pfTrans[1];
  double fT2      = pfTrans[2];

  double f00Start = fM00 * (double) n0Start;
  double f01Start = fM01 * (double) n0Start;
  double f02Start = fM02 * (double) n0Start;

  double f10Start = fM10 * (double) n1Start;
  double f11Start = fM11 * (double) n1Start;
  double f12Start = fM12 * (double) n1Start;

  double f20Start = fM20 * (double) n2Start;
  double f21Start = fM21 * (double) n2Start;
  double f22Start = fM22 * (double) n2Start;

  double fT00, fT01, fT02, fT10, fT11, fT12, fT20, fT21, fT22;
  double f0Sum, f1Sum, f2Sum;

  fT20           = f20Start + fT0;
  fT21           = f21Start + fT1;
  fT22           = f22Start + fT2;

  double *pfTemp;

  pfTemp = pfResult;

  int i, j, k;   // Loop counters

  for (k = n2Start; k <= n2End; k++)
    {
      fT20 = fT20 + fM20;
      fT21 = fT21 + fM21;
      fT22 = fT22 + fM22;

      fT10 = f10Start;
      fT11 = f11Start;
      fT12 = f12Start;

      for (j = n1Start; j <= n1End; j++)
        {

          fT10 = fT10 + fM10;
          fT11 = fT11 + fM11;
          fT12 = fT12 + fM12;

          f0Sum = fT20 + fT10;
          f1Sum = fT21 + fT11;
          f2Sum = fT22 + fT12;

          fT00 = f00Start;
          fT01 = f01Start;
          fT02 = f02Start;

          for (i = n0Start; i <= n0End; i++)
            {
              fT00      = fT00 + fM00;
              fT01      = fT01 + fM01;
              fT02      = fT02 + fM02;
              *pfTemp++ = f0Sum + fT00;
              *pfTemp++ = f1Sum + fT01;
              *pfTemp++ = f2Sum + fT02;
            } // end i
        }  // end j
    }  // end k
};


int nEigen2D(double a2x2fMat[2][2],double a2x2fEigenVec[2][2],double a2fEigenVals[2]) {
    int nRoots;
    int nVec;
    double f0,f1;

    nRoots = nSolveQuadratic(1.0,-(a2x2fMat[0][0] + a2x2fMat[1][1])/2.0,a2x2fMat[0][0]*a2x2fMat[1][1] - a2x2fMat[0][1]*a2x2fMat[1][0],&a2fEigenVals[0]);

    if (nRoots == 0)
        return 1;
    
    for (nVec = 0; nVec < 2; nVec++) {
        f0 = (a2x2fMat[0][0] - a2fEigenVals[nVec]);
        f0 = f0*f0 + a2x2fMat[1][0]*a2x2fMat[1][0];
        f1 = (a2x2fMat[1][1] - a2fEigenVals[nVec]);
        f1 = f1*f1 + a2x2fMat[0][1]*a2x2fMat[0][1];
        if (f0 > f1) {
            a2x2fEigenVec[nVec][0] = -a2x2fMat[1][0];
            a2x2fEigenVec[nVec][1] = a2x2fMat[0][0] - a2fEigenVals[nVec];
        } else {
            a2x2fEigenVec[nVec][0] = a2x2fMat[1][1] - a2fEigenVals[nVec];
            a2x2fEigenVec[nVec][1] = -a2x2fMat[0][1];
            std::swap(f0,f1);
        };
        
        if (nVec == 0) {
            if (sqrt(f0)/max(max(ABS(a2x2fMat[0][0]),ABS(a2x2fMat[1][1])),max(ABS(a2x2fMat[0][1]),ABS(a2x2fMat[1][0])))< 1e-6) {
                a2x2fEigenVec[0][0] = 1.0;
                a2x2fEigenVec[1][1] = 1.0;
                a2x2fEigenVec[0][1] = 0.0;
                a2x2fEigenVec[1][0] = 0.0;
                if (nRoots == 1)
                    a2fEigenVals[1] = a2fEigenVals[0];
                break;
            };
            
            if (nRoots == 1) {
                return 1;
            };
        };
    };
    
    for (nVec = 0; nVec < 2; nVec++) {
        f0 = sqrt(a2x2fEigenVec[nVec][0]*a2x2fEigenVec[nVec][0] + a2x2fEigenVec[nVec][1]*a2x2fEigenVec[nVec][1]);
        a2x2fEigenVec[nVec][0]/= f0;
        a2x2fEigenVec[nVec][1]/= f0;
    };
    return 0;
};

// Converted from old Fortran.  In this C++ version, please remember
// that the first element of an array has an index of 0.


void
vEigen(const int nDim, const int nMode, double *pfA, double *pfR)

// Routine vEIGEN
//
//        Purpose
//           Compute eigenvalues and eigenvectors of a real symmetric
//           matrix
//
//        Usage
//           int nDim, nMode;
//           double *pfA, *pfR;
//           vEigen(nDim, nMode, pfA, pfR)
//
//        Description of parameters
//         pfA - original matrix (symmetric), destroyed in computation.
//               resultant eigenvalues are developed in diagonal of
//               matrix pfA in descending order.
//  P.S. The equation for the indices of pfA is A[j][i] -> A[i+(j*j+j)/2])
//         pfRr - resultant matrix of eigenvectors (stored columnwise,
//                in same sequence as eigenvalues)
//         nDim - order of matrices pfA and pfR
//         nMode- input code
//                   0   Compute eigenvalues and eigenvectors
//                   1   Compute eigenvalues only (pfR need not be
//                       dimensioned but must still appear in calling
//                       sequence)
//
//        Remarks
//           original matrix pfA must be real symmetric (storage mode=1)
//           matrix pfA cannot be in the same location as matrix pfR
//
//        Subroutines and function subprograms required
//           sqrtf, abs
//
//        Method
//           Diagonalization method originated by Jacobi and adapted
//           by von Neumann for large computers as found in 'Mathematical
//           Methods for Digital Computers', edited by A. Ralston and
//           H.S. Wilf, John Wiley and Sons, New York, 1962, Chapter 7
//          See also Numerical Recipes, Chapter 11, subroutine JACOBI.
//        If a double precision version of this routine is desired, the
//        following variables should be declared double:
//
//     double precision a,r,anorm,anrmx,thr,x,y,sinx,sinx2,cosx,
//    1                 cosx2,sincs,range
//
//        The double precision version of this subroutine must also
//        contain double precision functions.  sqrtf in statements
//        below must be changed to sqrtd.  abs below
//        must be changed to absd. The fRANGE constant should
//        be changed to 1.0d-12.

{

  int   i, j, k, l, m;
  int   nIQ, nJQ, nLL, nMM, nLM, nILR, nIMR, nILQ, nIMQ, nIL, nIM, nLQ, nMQ;
  int   nInd;
  double fX, fY, fSinX, fCosX, fSinX2, fCosX2, fSinCS, fTHR, fANORM, fANRMX;
  double *pfTemp;

  double fRANGE = (double)1.0E-6;

  // Generate identity matrix in pfR

  pfTemp = pfR;
  if (0 == nMode)
    {
      for (j = 0; j < nDim; j++)
	{
	  for (i = 0; i < nDim; i++)
	    {
	      if (i == j)
		*pfTemp++ = 1.0;
	      else
		*pfTemp++ = 0.0;
	    }
	}
    }

  // Compute initial and final norms (anorm and anormx)

  fANORM = 0.0;
  for (i = 0; i < nDim; i++)
    {
      for (j = i; j < nDim; j++)
	{
	  if (i != j)
	    {
	      nIQ = i + (j * j  + j) / 2;
	      fANORM = fANORM + pfA[nIQ] * pfA[nIQ];
	    }
	}
    }
  if (fANORM > 0.0)
    {
      fANORM = 1.41421356 * sqrt(fANORM);
      fANRMX = fANORM * fRANGE / (double) nDim;

      // Initialize indicators and compute threshold, fThr

      fTHR = fANORM;

      do
	{                                   // until (fTHR <= fANRMX)
	  fTHR = fTHR / (double) nDim;
	  do
	    {                                 // until (nInd == 0){
	      nInd = 0;
	      for (l = 0; l < nDim - 1; l++)
		{   // from start to next to last col
		  for (m = l + 1; m < nDim; m++)
		    { // from l+1 to last col

		      // Compute sin and cos

		      nMQ = (m * m  + m) / 2;        //     [m][m]
		      nLQ = (l * l  + l) / 2;        //     [l][l]
		      nLM = l + nMQ;
		      if ( fabs(pfA[nLM]) >= fTHR)
			{
			  nInd = 1;         // Set flag that we did some calcs
			  nLL  = l + nLQ;
			  nMM  = m + nMQ;
			  fX = 0.5 * (pfA[nLL] - pfA[nMM]);
			  fY = -pfA[nLM] / sqrt(pfA[nLM] * pfA[nLM]
						 +  fX * fX);
			  if (fX < 0.0) fY = -fY;
			  fSinX  = fY / sqrt(2.0 * (1.0 +
						     ( sqrt(1.0 - fY * fY))));
			  fSinX2 = fSinX * fSinX;
			  fCosX  = sqrt(1.0 - fSinX2);
			  fCosX2 = fCosX * fCosX;
			  fSinCS = fSinX * fCosX;

			  // Rotate l and m columns
			
			  nILQ = nDim * l;
			  nIMQ = nDim * m;
			  for (i = 0; i < nDim; i++)
			    {
			      nIQ = (i * i  +  i) / 2;
			      if (i != l)
				{
				  if (i != m)
				    {
				      if (i < m)
					nIM = i + nMQ;
				      else
					nIM = m + nIQ;
				      if (i < l)
					nIL = i + nLQ;
				      else
					nIL = l + nIQ;
				      fX       = pfA[nIL] * fCosX
					         -  pfA[nIM] * fSinX;
				      pfA[nIM] = pfA[nIL] * fSinX
					         +  pfA[nIM] * fCosX;
				      pfA[nIL] = fX;
				    }
				}
			      if (0 == nMode)
				{
				  nILR      = nILQ + i;
				  nIMR      = nIMQ + i;
				  fX        = pfR[nILR] * fCosX
				              -  pfR[nIMR] * fSinX;
				  pfR[nIMR] = pfR[nILR] * fSinX
				              +  pfR[nIMR] * fCosX;
				  pfR[nILR] = fX;
				}
			    } // end of i loop
			
			  fX       = 2.0 * pfA[nLM] * fSinCS;
			  fY       =  pfA[nLL] * fCosX2  +  pfA[nMM] * fSinX2
			              -  fX;
			  fX       =  pfA[nLL] * fSinX2  +  pfA[nMM] * fCosX2
			              +  fX;
			  pfA[nLM] = (pfA[nLL] - pfA[nMM]) * fSinCS  +
			              pfA[nLM] * (fCosX2 - fSinX2);
			  pfA[nLL] =  fY;
			  pfA[nMM] =  fX;
			} // end if
		    } // end m loop
		} // end l loop
	    } while (1 == nInd);     // Repeat until inner loops not entered

	  //   Compare threshold with final norm

	} while (fTHR > fANRMX);

      //  Sort eigenvalues and eigenvectors

      nIQ = -nDim;
      for (i = 0; i < nDim; i++)
	{
	  nIQ = nIQ + nDim;
	  nLL  = i + (i * i + i) / 2;          // element [i][i]
	  nJQ  = nDim * (i - 1);
	  for (j = i; j < nDim; j++)
	    {
	      nJQ = nJQ + nDim;
	      nMM = j + (j * j + j) / 2;         // element [j][j]
	      if ( pfA[nLL] < pfA[nMM])
		{
		  // Swap diagonal elements A[nLL} and A[nMM]

		  fX      = pfA[nLL];
		  pfA[nLL] = pfA[nMM];
		  pfA[nMM] = fX;
		  if (0 == nMode)
		    {

		      // Swap off-diagonal elements

		      for (k = 0; k < nDim; k++)
			{
			  nILR      = nIQ + k;
			  nIMR      = nJQ + k;
			  fX        = pfR[nILR];
			  pfR[nILR] = pfR[nIMR];
			  pfR[nIMR] = fX;
			}
		    }
		}
	    } // end j
	} // end i
    }
}

void
vTRED2(const int nDim, double **ppfMat,
		   double *pfDiag, double *pfOffDiag)
{
  //  Householder reduction of a real, symmetric matrix in
  //	**pfMat[1..nDim][1..nDim] (NOTE FORTRAN-style indexing!)
  //  Matrix **pfMat is destroyed in the computation and on return
  //  holds the orthogonal Q matrix effecting the transformation.
  //  pfDiag[1..nDim] holds the diagonal elements
  //  pfOffDiag[1] is set to 0
  //  pfOffDiag[2..n] holds the off-diagonal elements
  //
  // This code is derived from algorithms published in
  // Numerical Recipes in C: The Art of Scientific Programming
  // (c) 1988, Cambridge University Press.
  // See pp. 373ff.

  int    i, j, k, L;          // Loop counters
                              // (L uppercase so not confused with 1
  double  fScale, fHH, fH, fG, fF;

  for (i = nDim; i >=2; i--)
    {
      L = i - 1;
      fH = fScale = 0.0;
      if (L > 1)
        {
          for (k = 1; k <= L; k++)
            {
              fScale += fabs(ppfMat[i][k]);
            }
          if (fScale == 0.0)
            {
              // Skip transformation
              pfOffDiag[i] = ppfMat[i][L];
            }
          else
            {
              for (k = 1; k <= L; k++)
                {
                  ppfMat[i][k] /= fScale;               // Use scale a's f tra
                  fH += ppfMat[i][k] * ppfMat[i][k];  // Form sigma in h
                }
              fF             = ppfMat[i][L];
              fG             = fF > 0.0 ? -sqrt(fH) : sqrt(fH);
              pfOffDiag[i]   = fScale * fG;
              fH            -= fF * fG;
              ppfMat[i][L] = fF - fG;
              fF             = 0.0;
              for (j = 1; j <= L; j++)
                {
                  ppfMat[j][i] = ppfMat[i][j] / fH;  // Store u/H in ith col
                  fG             = 0.0;
                  for (k = 1; k <= j; k++)
                    {
                      fG        += ppfMat[j][k] * ppfMat[i][k];
                    }
                  for (k = j+1; k <= L; k++)
                    {
                      fG        += ppfMat[k][j] * ppfMat[i][k];
                    }
                  pfOffDiag[j]  = fG / fH;// Form element of p in temp unused el
                  fF      += pfOffDiag[j] * ppfMat[i][j];

                } // end j loop
              fHH = fF / (fH + fH);
              for (j = 1; j <= L; j++)
                {
                  fF           = ppfMat[i][j];
                  pfOffDiag[j] = fG = pfOffDiag[j] - fHH * fF;
                  for (k = 1; k <= j; k++)
                    {
                      ppfMat[j][k] -= (fF * pfOffDiag[k] +  fG * ppfMat[i][k]);
                    }
                } // end j loop
            }
        }
      else
        {
          pfOffDiag[i] = ppfMat[i][L];
        }
      pfDiag[i] = fH;
    }
  pfDiag[1]    = 0.0;
  pfOffDiag[1] = 0.0;
  for (i = 1; i <= nDim; i++)
    {
      L = i - 1;
      if (pfDiag[i])
        {               // This skipped when i=0 jwp: I do not like this! :(
          for (j = 1; j <= L; j++)
            {
              fG = 0.0;
              for (k = 1; k <= L; k++)
                {
                  fG += ppfMat[i][k] * ppfMat[k][j];
                }
              for (k = 1; k <= L; k++)
                {
                  ppfMat[k][j] -= fG * ppfMat[k][i];
                }
            }
        }
      pfDiag[i]    = ppfMat[i][i];
      ppfMat[i][i] = 1.0;
      for (j = 1; j <= L; j++)
        {
          ppfMat[j][i] = ppfMat[i][j] = 0.0;
        }
    }
}

void
vTQLI(const int nDim, double *pfDiag, double *pfOffDiag,
      double **ppfEigenVecs)

{
  // QL algorithm with implicit shifts, to determine the eigenvalues and
  // eigenvectors of a real, symmetric matrix previously reduced by vTRED2.
  // On input pfDiag[1..nDim] contains the diagonal elements of the tridiagonal
  // matrix.  On output, it returns the eigenvalues.  The vector
  // pfOffDiag[1..nDim] inputs the subdiagonal elements of the tridiagonal
  // matrix with pfOffDiag[1] arbitrary.  On output pfOffDiag is destroyed.  The
  // matrix ppfEigenVecs is input as the identity matrix OR if the eigenvectors
  // of a matrix that has been reduced by vTRED2 are required, then
  // ppfEigenVecs is the matrix output by vTRED2.  In either case, the kth
  // column of ppfEigenVecs returns the normalized eigenvector corresponding
  // to ppfDiag[k].
  //
  // This code is derived from algorithms published in
  // Numerical Recipes in C: The Art of Scientific Programming
  // (c) 1988, Cambridge University Press.
  // See pp. 380ff.

  int m, L, iter, i, k;   // Loop counters  lowercase l is too easily confused
                          // with digit 1
  double fS, fR, fP, fG, fF, fDD, fC, fB;

  for (i = 2; i <= nDim; i++)
    {
      pfOffDiag[i-1] = pfOffDiag[i]; // Convenient to renumber array offdiag
    }
  pfOffDiag[nDim] = 0.0;

  for (L = 1; L <= nDim; L++)
    {
      iter = 0;
      do
	{
	  for (m = L; m <= nDim-1; m++)
	    {
	      fDD = fabs(pfDiag[m]) + fabs(pfDiag[m+1]);

	      // Look for a single small subdiagonal element to split the matrix

	      if ( (double) (fabs(pfOffDiag[m]) + fDD) == fDD) break;
	    }

	  // If there was no break, then m = n hopefully

	  if (m != L)
	    {
	      if (iter++ == 30)
		{
		  cout << "Too many iterations in vTQLI!\n";
		  return;
		}
	      fG = (pfDiag[L+1] - pfDiag[L]) / (2.0 * pfOffDiag[L]);
	      fR = sqrt( (fG * fG) + 1.0);
	      fG = pfDiag[m] - pfDiag[L] + pfOffDiag[L] / (fG + SIGND(fR, fG));
	      fS = fC = 1.0;
	      fP = 0.0;
	      for (i = m-1; i >= L; i--)
		{
		  fF = fS * pfOffDiag[i];
		  fB = fC * pfOffDiag[i];
		  if (fabs(fF) >= fabs(fG))
		    {
		      fC              = fG / fF;
		      fR              = sqrt( (fC * fC) + 1.0);
		      pfOffDiag[i+1] = fF * fR;
		      fC             *= (fS = 1.0 / fR);
		    }
		  else
		    {
		      fS              = fF / fG;
		      fR              = sqrt( (fS * fS) + 1.0);
		      pfOffDiag[i+1]  = fG * fR;
		      fS             *= (fC = 1.0 / fR);
		    }
		  fG           =  pfDiag[i+1] - fP;
		  fR           = (pfDiag[i]   - fG) * fS  +  2.0 * fC * fB;
		  fP           = fS * fR;
		  pfDiag[i+1]  = fG + fP;
		  fG           = fC * fR  - fB;
		  for (k = 1; k <= nDim; k++)
		    {
		      fF                   = ppfEigenVecs[k][i+1];
		      ppfEigenVecs[k][i+1] = fS * ppfEigenVecs[k][i] + fC * fF;
		      ppfEigenVecs[k][i]   = fC * ppfEigenVecs[k][i] - fS * fF;
		    }
		}
	      pfDiag[L]    = pfDiag[L] - fP;
	      pfOffDiag[L] = fG;
	      pfOffDiag[m] = 0.0;
	    } // endif  m != L
	} while (m != L);
    } // end l loop
}

void
vEigsrt(const int nDim, double *pfEigval, double **ppfEigvec)
{
  // Given the eigenvalues pfEigval[1...nDim] and eigenvectors
  // ppfEigvec[1...nDim][1...nDim] as output from vTQLI, this routine
  // sorts the eigenvalues into ASCENDING order, and rearranges the
  // columns of ppfEigvec accordingly.
  // From Numerical Recipes in C Section 11.1 p. 366.

  int i, j, k;
  double fTemp;

  for (i = 1; i <= nDim; i++)
    {
      fTemp = pfEigval[k=i];
      for (j = i + 1; j <= nDim; j++)
	{
	  if (pfEigval[j] >= fTemp) fTemp = pfEigval[k=j];     // Save largest
	}
      if (k != i)
	{
	  pfEigval[k] = pfEigval[i];                  // Need to swap
	  pfEigval[i] = fTemp;
	  for (j = 1; j <= nDim; j++)
	    {
	      fTemp           = ppfEigvec[j][i];
	      ppfEigvec[j][i] = ppfEigvec[j][k];
	      ppfEigvec[j][k] = fTemp;
	    }
	}
    }
}

int
nMatrixCompress(const int nNumRowsCols, const int *pnDelFlags,
		double *pfMatrix, double *pfVector)
{

  // Delete rows and columns of a square matrix
  //       pfMatrix[nNumRowsCols][nNumRowsCols]
  // according to values of pnDelFlags.  If pnDelFlags[i] is non-zero,
  // then delete row and column i.
  // Return the number of rows and columns
  // (nNum) left in pfMatrix which is now pfMatrix[nNum][nNum];
  // Delete elements pfVector, so that it matches ...

  int i, j;             // Loop counters
  int nNum;
  double *pfTempMat;
  double *pfTempVec;

  nNum = 0;
  pfTempMat = pfMatrix;
  pfTempVec = pfVector;

  for (i = 0; i < nNumRowsCols; i++)
    {
      if (0 == pnDelFlags[i])
	{
	  nNum++;
	  *pfTempVec++ = pfVector[i];   // This element stays in the vector
	}
      for (j = 0; j < nNumRowsCols; j++)
	{
	  if ( (0 == pnDelFlags[i]) && (0 == pnDelFlags[j]) )
	    {
	      // This element stays in matrix

	      *pfTempMat++ = pfMatrix[i * nNumRowsCols + j];
	    }
	}
    }
  return (nNum);
}


int
nRealFFT1D(const int nSize,const int nSign, double *pfData) {
    int i,i1,i2,i3,i4;
    double c1 = 0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;
    int n = nSize;
    theta = 3.141592653589793238/((double) (n>>1));
    if (nSign == 1) {
        c2 = -0.5;
        nFFT1D( nSize/2,1,pfData);
    } else {
        c2 = 0.5;
        theta = -theta;
    };
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    for (i=1;i < (n >> 2); i++) {
        i2 = 1 + (i1 = i + i);
        i4 = 1 + (i3 = n - i1);
        h1r = c1 *(pfData[i1] + pfData[i3]);
        h1i = c1 * (pfData[i2] - pfData[i4]);
        h2r = -c2 *(pfData[i2] + pfData[i4]);
        h2i = c2*(pfData[i1] - pfData[i3]);
        pfData[i1] = h1r+ wr*h2r - wi*h2i;
        pfData[i2] = h1i + wr*h2i + wi*h2r;
        pfData[i3] = h1r - wr*h2r + wi*h2i;
        pfData[i4] = -h1i + wr*h2i + wi*h2r;
        wr = (wtemp = wr)*wpr - wi*wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    if (nSign == 1) {
        pfData[0] = (h1r = pfData[0]) + pfData[1];
        pfData[1] = h1r - pfData[1];
    } else {
        pfData[0] = c1*((h1r = pfData[0]) + pfData[1]);
        pfData[1] = c1*(h1r - pfData[1]);
        nFFT1D(nSize/2,-1,pfData);
    };
    return 0;
};

#define SWAPFFT(a, b) tempr=(a); (a)=(b); (b)=tempr;

int
nFFT1D(const int nSize, const int nSign, double *pfData)
{
    // 1D Fast Fourier Transform routine from "Numerical Recipes in C", pp. 411-
    // nSize MUST be an integer power of 2
    // No checks are made on the input arguments.
    
    
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    int nn;
    int isign;

    // The array convention here is Fortran not C.  Decrement the array pointer to compensate.
    pfData--;
    
    isign = nSign;
    nn = nSize;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i+=2)
    {
        if (j > i)
        {
            SWAPFFT(pfData[j],   pfData[i]);
            SWAPFFT(pfData[j+1], pfData[i+1]);
        }
        m = n >> 1;
        while (m >=2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j+= m;
    }
    mmax = 2;
    while (n > mmax)
    {
        istep = 2 * mmax;
        theta = 6.28318530717959 / (double)(isign * mmax);
        wtemp = sin(0.5 * theta);
        wpr   = -2.0 * wtemp * wtemp;
        wpi   = sin(theta);
        wr    = 1.0;
        wi    = 0.0;
        for (m = 1; m < mmax; m+=2)
        {
            for (i = m; i <=n; i+=istep)
            {
                j          = i + mmax;
                tempr      = (double)(wr * pfData[j]    -  wi * pfData[j+1]);
                tempi      = (double)(wr * pfData[j+1]  +  wi * pfData[j]);
                pfData[j]    = pfData[i]   - tempr;
                pfData[j+1]  = pfData[i+1] - tempi;
                pfData[i]   += tempr;
                pfData[i+1] += tempi;
            }
            wr = (wtemp=wr) * wpr  -  wi * wpi  + wr;
            wi = wi * wpr  + wtemp * wpi  + wi;
        }
        mmax = istep;
    }
    return (0);
}





/*
double* pfCast2double(int nSize,float* pfMat,bool bPush = TRUE);
float* pfCast2float(int nSize,double* pfMat,bool bPush = TRUE);
nCast2doublePop() { pfCast2double(0,NULL,FALSE); return 1; };
nCast2floatPop() { pfCast2float(0,NULL,FALSE); return 1; };
*/

void vCopyMat3D(const float *pfMat1, double *pfMat2) {
    for (int i = 0; i < 9; i++) pfMat2[i] = pfMat1[i];
};

void vCopyMat3D(const double *pfMat1, float *pfMat2) {
    for (int i = 0; i < 9; i++) pfMat2[i] = pfMat1[i];
};

void vCopyVec3D(const float *pfVec1, double *pfVec2) {
    for (int i = 0; i < 3; i++) pfVec2[i] = pfVec1[i];
};

void vCopyVec3D(const double *pfVec1, float *pfVec2) {
    for (int i = 0; i < 3; i++) pfVec2[i] = pfVec1[i];
};

void vCopyVecND(const int nDim, const float *pfVec1, double *pfVec2) {
    for (int i = 0; i < nDim; i++) pfVec2[i] = pfVec1[i];
};

void vCopyVecND(const int nDim, const double *pfVec1, float *pfVec2) {
    for (int i = 0; i < nDim; i++) pfVec2[i] = pfVec1[i];
};

float fDot3D (const double fVec1[3], const float fVec2[3])
{
  return (fVec1[0] * fVec2[0]  +  fVec1[1] * fVec2[1]  + fVec1[2] * fVec2[2]);
}

float fDot3D (const float fVec1[3], const double fVec2[3])
{
  return (fVec1[0] * fVec2[0]  +  fVec1[1] * fVec2[1]  + fVec1[2] * fVec2[2]);
}


float* pfCast3D(const double* pfIn,int nWhich) {
    const int nMaxWhich = 5;
    static float anx3fVecs[nMaxWhich][3];
    if (nWhich>nMaxWhich)
        return (float*)NULL;
    vCopyVec3D(pfIn,anx3fVecs[nWhich]);
    return &anx3fVecs[nWhich][0];
};
float* pfCastND(int nDim,const double* pfIn,int nWhich) {
    const int nMaxWhich = 5;
    const int nMaxDim = 10;
    static float anxnfVecs[nMaxWhich][nMaxDim];
    if (nWhich>nMaxWhich)
        return (float*)NULL;
    vCopyVecND(nDim,pfIn,anxnfVecs[nWhich]);
    return &anxnfVecs[nWhich][0];
};

double* pfCast3D(const float* pfIn,int nWhich) {
    const int nMaxWhich = 5;
    static double anx3fVecs[nMaxWhich][3];
    if (nWhich>nMaxWhich)
        return (double*)NULL;
    vCopyVec3D(pfIn,anx3fVecs[nWhich]);
    return &anx3fVecs[nWhich][0];
};
double* pfCastND(int nDim,const float* pfIn,int nWhich) {
    const int nMaxWhich = 5;
    const int nMaxDim = 10;
    static double anxnfVecs[nMaxWhich][nMaxDim];
    if (nWhich>nMaxWhich)
        return (double*)NULL;
    vCopyVecND(nDim,pfIn,anxnfVecs[nWhich]);
    return &anxnfVecs[nWhich][0];
};


int nCompVecND(const int nDim,double *pfComp1,double *pfComp2) {
    int nx;
    for (nx = 0; nx < nDim; nx++) {
        if (pfComp1[nx] != pfComp2[nx])
            return 1;
    };
    return 0;
};


int nSolveQuadratic(const double fA, const double fB, const double fC,double *pfV)
{
  // This solves the special quadratic equation: AX^2 + 2BX + C = 0
  // used by the nPredictReflns routine.  The (up to 2) solutions are returned
  // in pfV[0] and pfV[1].  The number of solutions is returned by the routine.

  double fD;

  if (fA != 0.0) {
    fD = fB * fB - fA * fC;
    if (fD == 0.0) {
      *pfV = -fB / fA;
      return (1);
    }
    else if (fD > 0.0) {
      fD = sqrt(fD);
      pfV[0] = (-fB + fD)  / fA;
      pfV[1] = (-fB - fD)  / fA;
      return (2);
    }
  }
  else if (fB != 0.0) {

    //  Here if A = 0 and B not = 0

    *pfV = -0.5f * fC / fB;
    return (1);
  }
  return (0);
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
///////////  Sort Routines  /////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

int float_cmp(const void* a,const void* b) {
    float* pfa = (float*) a;
    float* pfb = (float*) b;
    if (*pfa>*pfb) {
        return 1;
    } else if (*pfa==*pfb)
        return 0;
    else
        return -1;
};


int double_cmp(const void* a,const void* b) {
    double* pfa = (double*) a;
    double* pfb = (double*) b;
    if (*pfa>*pfb) {
        return 1;
    } else if (*pfa==*pfb)
        return 0;
    else
        return -1;
};


int int_cmp(const void* a,const void* b) {
    int* pna = (int*) a;
    int* pnb = (int*) b;
    if (*pna>*pnb) {
        return 1;
    } else if (*pna==*pnb)
        return 0;
    else
        return -1;
};


float* g_pfCmpFloats = (float*)NULL;
int float_cmp_rel(const void* a,const void* b) {
    int* pna = (int*) a;
    int* pnb = (int*) b;
    if (g_pfCmpFloats[*pna]>g_pfCmpFloats[*pnb]) {
        return 1;
    } else if (g_pfCmpFloats[*pna]==g_pfCmpFloats[*pnb])
        return 0;
    else
        return -1;
};
double* g_pfCmpDoubles = (double*)NULL;
int double_cmp_rel(const void* a,const void* b) {
    int* pna = (int*) a;
    int* pnb = (int*) b;
    if (g_pfCmpDoubles[*pna]>g_pfCmpDoubles[*pnb]) {
        return 1;
    } else if (g_pfCmpDoubles[*pna]==g_pfCmpDoubles[*pnb])
        return 0;
    else
        return -1;
};


int* g_pnCmpInts = (int*)NULL;

int int_cmp_rel(const void* a, const void* b)
{
    int* pna = (int*) a;
    int* pnb = (int*) b;
    
    if( g_pnCmpInts[*pna] > g_pnCmpInts[*pnb] )
        return 1;
    else if( g_pnCmpInts[*pna] == g_pnCmpInts[*pnb] )
        return 0;
    else
        return -1;
}

//RB this function is using the difference between the values of a and b as a tie-breaker, when g_pnCmpInts is the same for both a and b. 
int int_cmp_rel_2(const void* a, const void* b)
{
    int* pna = (int*) a;
    int* pnb = (int*) b;
    
    if( g_pnCmpInts[*pna] > g_pnCmpInts[*pnb] )
        return 1;
    else if( g_pnCmpInts[*pna] < g_pnCmpInts[*pnb] )
        return -1;
    else if( g_pnCmpInts[*pna] == g_pnCmpInts[*pnb] )
    {
        if( *pna > *pnb )
            return 1;
        else if( *pna < *pnb )
            return -1;
        else
            return 0;
    }
    else
        return 0;
}

int unsigned_short_int_cmp(const void* a,const void* b)
{
    unsigned short int* puia = (unsigned short int*) a;
    unsigned short int* puib = (unsigned short int*) b;
    if (*puia>*puib) {
        return 1;
    } else if (*puia==*puib)
        return 0;
    else
        return -1;
}

Cstring* g_psCmpCstrings = (Cstring*)NULL;

int Cstring_cmp_rel(const void* a,const void* b) {
    int* pna = (int*) a;
    int* pnb = (int*) b;
    if (g_psCmpCstrings[*pna]>g_psCmpCstrings[*pnb]) {
        return 1;
    } else if (g_psCmpCstrings[*pna]==g_psCmpCstrings[*pnb])
        return 0;
    else
        return -1;
};


float g_fCmpFloatsTol;
int float_cmp_tol(const void* a,const void* b) {
    float* pfa = (float*) a;
    float* pfb = (float*) b;
    if (*pfa>*pfb + g_fCmpFloatsTol) {
        return 1;
    } else if (*pfa < *pfb - g_fCmpFloatsTol)
        return -1;
    else
        return 0;
};

double g_fCmpDoublesTol;
int double_cmp_tol(const void* a,const void* b) {
    double* pfa = (double*) a;
    double* pfb = (double*) b;
    if (*pfa>*pfb + g_fCmpDoublesTol) {
        return 1;
    } else if (*pfa < *pfb - g_fCmpDoublesTol)
        return -1;
    else
        return 0;
};



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////////////  Integer math function /////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


int nInvMat3D(int* pnMat1,int* pnMat2) {
    int nDet;

    nDet  =   pnMat1[0] * (pnMat1[4] * pnMat1[8]  -  pnMat1[7] * pnMat1[5])
        - pnMat1[1] * (pnMat1[3] * pnMat1[8]  -  pnMat1[5] * pnMat1[6])
        + pnMat1[2] * (pnMat1[3] * pnMat1[7]  -  pnMat1[4] * pnMat1[6]);
    if ( (0 != nDet) && (NULL != pnMat2) )
    {
        pnMat2[0] =  (pnMat1[4] * pnMat1[8]  -  (pnMat1[5] * pnMat1[7]));
        pnMat2[1] = -(pnMat1[1] * pnMat1[8]  -  (pnMat1[2] * pnMat1[7]));
        pnMat2[2] =  (pnMat1[1] * pnMat1[5]  -  (pnMat1[2] * pnMat1[4]));
        pnMat2[3] = -(pnMat1[3] * pnMat1[8]  -  (pnMat1[5] * pnMat1[6]));
        pnMat2[4] =  (pnMat1[0] * pnMat1[8]  -  (pnMat1[2] * pnMat1[6]));
        pnMat2[5] = -(pnMat1[0] * pnMat1[5]  -  (pnMat1[2] * pnMat1[3]));
        pnMat2[6] =  (pnMat1[3] * pnMat1[7]  -  (pnMat1[4] * pnMat1[6]));
        pnMat2[7] = -(pnMat1[0] * pnMat1[7]  -  (pnMat1[1] * pnMat1[6]));
        pnMat2[8] =  (pnMat1[0] * pnMat1[4]  -  (pnMat1[1] * pnMat1[3]));
    }
    return (nDet);
};

void vMulMat3DVec3D(int* pnMat,int* pnVec1,int* pnVec2) {
    int i,j;
    for (i = 0; i < 3; i++)
    {
        pnVec2[i] = 0;
        for (j = 0; j < 3; j++)
        {
            pnVec2[i] = pnVec2[i] + (pnMat[j*3+i] * pnVec1[j]);
        }
    };
};

bool bDivides(int nDim,int* pnVec,int nFactor) {
    int i;
    if (nFactor<=0)
        return false;
    for (i = 0; i < nDim; i++) {
        if (abs(pnVec[i]) % nFactor)
            return false;
    };
    return true;
};


bool bDivides(int nDim,int* pnVec) {
    int nMax;
    int nx;
    int nFactor;
    nMax = 0;
    for (nx = 0; nx < nDim;nx++)
        nMax = max(abs(pnVec[nx]),nMax);
    nFactor = nMax;
    while (nFactor >1) {
        if (bDivides(nDim,pnVec,nFactor))
            return true;
        nFactor--;
    };
    return false;
};





int nSplineCompute(double* pfX,double* pfY,int nPoints,double fYp0,double fYpN,double* pfYpp) {
    
    int i,k;
    double p,qn,sig,un;
    itr<double> u;
    
    
    if (fYp0 > 0.99e30)         
        pfYpp[0]=u[0]=0.0;
    else { 
        pfYpp[0] = -0.5;
        u[0]=(3.0/(pfX[1]-pfX[0]))*((pfY[1]-pfY[0])/(pfX[1]-pfX[0])-fYp0);
    }
    
    for (i=1;i<nPoints-1;i++) { 
        sig=(pfX[i]-pfX[i-1])/(pfX[i+1]-pfX[i-1]);
        p=sig*pfYpp[i-1]+2.0;
        pfYpp[i]=(sig-1.0)/p;
        u[i]=(pfY[i+1]-pfY[i])/(pfX[i+1]-pfX[i]) - (pfY[i]-pfY[i-1])/(pfX[i]-pfX[i-1]);
        u[i]=(6.0*u[i]/(pfX[i+1]-pfX[i-1])-sig*u[i-1])/p;
    }
    
    if (fYpN > 0.99e30) 
        qn=un=0.0;            
    else {             
        qn=0.5;
        un=(3.0/(pfX[nPoints-1]-pfX[nPoints-2]))*(fYpN-(pfY[nPoints-1]-pfY[nPoints-2])/(pfX[nPoints-1]-pfX[nPoints-2]));
    }    
    pfYpp[nPoints-1]=(un-qn*u[nPoints-2])/(qn*pfYpp[nPoints-2]+1.0);
    for (k=nPoints-2;k>=0;k--) 
        pfYpp[k]=pfYpp[k]*pfYpp[k+1]+u[k];
    return 0;
};

DTREK_EXPORT double fSplineCalc(double fX,double fX1,double fY1,double fYpp1,double fX2,double fY2,double fYpp2) {
    double h;
    double a,b;
    double fReturn;

    if (fX1 > fX2) {
        std::swap(fX1,fX2);
        std::swap(fY1,fY2);
        std::swap(fYpp1,fYpp2);
    };

    h = fX2 - fX1;
        
    a=(fX2-fX)/h;
    b=(fX-fX1)/h; 
        
    fReturn = a*fY1 + b*fY2 + ((a*a*a-a)*fYpp1 + (b*b*b-b)*fYpp2)*(h*h)/6.0;
    return fReturn;

};

/*  Notes on ellipsoids:

    1)  Every ellipsoid is described by the equation xtAx + btx + c
    2)  x is in units of pixels relative to some reference point.  In C3Ddata::nCalcMask() when the ellipse is first calculated, this
        reference point is the shoebox coordiante <0,0>.  Before exiting, however, nCalcMask()
        'normalizes' the reference point for each slice ellipsoid to be the centroid of the shoebox.  Both coordinates are relative to
        the C3Ddata::m_nOrigOffset
    3)  In C3Ddata::nCalcMask() the ellipsoid defines the spot by it's value.  A pixel who's value is <= 1.0 is in the spot... otherwise
        it's outside.
    4)  The btx + c portion of the equation controls the offseting of the ellipsoid.  Thus, we really have a second *implicit* offset.

    5)  To pass to tilt/angle form, we need to figure out the real center of the ellipsoid.  This center should be applied to the
        reference point mentioned in (2)

        *two* offsets that are encoded
        in the xtAx + btx + c


*/

DTREK_EXPORT int nGetEllipsoidTiltAngleCenter(double* a2x2fEllipsoidA,double* a2fEllipsoidb,double fEllipsoidc,double* a2fRelOffset,double* a2fMajorMinor,double* pfRotOffset) {
    // Need to find the true offset of the ellipsoid, and return this to the calling routine.  The calculated ellipse will have a center at that point.

    double a2x2fInvMat[2][2];
    double a2fShift[2];
    double a2x2fShiftEllipsoidA[2][2];
    double fShiftc;
    double a2x2fEigenVecs[2][2];
    double a2fEigenVals[2];

    double a2fTemp[2];
    
    nInvMatND_svd(2,&a2x2fEllipsoidA[0],&a2x2fInvMat[0][0]);
    vMulMatNDVecND(2,&a2x2fInvMat[0][0],&a2fEllipsoidb[0],&a2fShift[0]);
    a2fShift[0]/=2.0;
    a2fShift[1]/=2.0;
    a2fRelOffset[0] = -a2fShift[0];
    a2fRelOffset[1] = -a2fShift[1];

    // Need to reclaculate the 'c' value and scale the matrix accordingly.
    vMulMatNDVecND(2,a2x2fEllipsoidA,&a2fShift[0],&a2fTemp[0]);
    fShiftc = fEllipsoidc -(a2fShift[0]*a2fEllipsoidb[0] + a2fShift[1]*a2fEllipsoidb[1]) + (a2fTemp[0]*a2fShift[0] + a2fTemp[1]*a2fShift[1]);
    if (fShiftc >= 1.0)
        return 1;
    vMulVecNDScalar(4,a2x2fEllipsoidA,1.0/(1.0 - fShiftc),&a2x2fShiftEllipsoidA[0][0]);

    // Find eigenvectors/eigenvalues.
    vEigen(2,&a2x2fShiftEllipsoidA[0][0],&a2fEigenVals[0],&a2x2fEigenVecs[0][0]);
    if ((a2fEigenVals[0] <= 0.0) || (a2fEigenVals[1] <= 0.0))
        return 1;
    a2fMajorMinor[0] = 1.0/sqrt(a2fEigenVals[0]);
    a2fMajorMinor[1] = 1.0/sqrt(a2fEigenVals[1]);

    if (a2fMajorMinor[0] < a2fMajorMinor[1]) {
        std::swap(a2fMajorMinor[0],a2fMajorMinor[1]);
        std::swap(a2x2fEigenVecs[0][0],a2x2fEigenVecs[1][0]);
        std::swap(a2x2fEigenVecs[0][1],a2x2fEigenVecs[1][1]);
    };

    
    *pfRotOffset = acos(a2x2fEigenVecs[0][0]) / Gs_dRADIANS_PER_DEGREE;
    if (a2x2fEigenVecs[0][1]<0.0)
        *pfRotOffset  = 360.0 - *pfRotOffset;


    // Test:  Make a vector to multiply the quadratic by.
    /*
    int nx;
    for (nx = 0; nx < 2; nx++) {
        a2fTemp[0] = a2x2fEigenVecs[nx][0]*a2fMajorMinor[nx];
        a2fTemp[1] = a2x2fEigenVecs[nx][1]*a2fMajorMinor[nx];

        vMulMatNDVecND(2,&a2x2fShiftEllipsoidA[0][0],a2fTemp,a2fTemp2);
        f0 = a2fTemp2[0]*a2fTemp[0] + a2fTemp2[1]*a2fTemp[1];
    };
    */
    


    return 0;
};

DTREK_EXPORT double fEvalEllipse(double fPix0, double fPix1, double* pfA, double* pfb, double* pfc)
{
    double      fValue = fPix0 * fPix0 * (pfA[0]) + 
                         fPix0 * fPix1 * (pfA[1] + pfA[2]) + 
                         fPix1 * fPix1 * (pfA[3]) + 
                         fPix0 * pfb[0] + 
                         fPix1 * pfb[1] + 
                         *pfc;
    
    return fValue;
}

DTREK_EXPORT int nBuildEllipse(double a2fMajorMinor[2], 
                               double a2fRelOffset[2],
                               double fRotOffset,
                               double* a2x2fEllipsoidAOut,
                               double* a2fEllipsoidbOut,
                               double* pfEllipsoidcOut)
{
    const int   c_nPoints = 8;
    
    double      a2fPoint[2];
    double      a2x2lfGMat[2][2];
    double      a2fhVec[2];

    double      pfCoeffA[7][7];
    double      pfCoeffB[7];
    double*     pf0 = &pfCoeffA[0][0];
    double*     pf1 = &pfCoeffB[0];

    // Calculate the new ellipsoid.
    nInitQuadratic(2, FALSE, &pf0, &pf1);
    
    for(int nx = 0; nx < c_nPoints; nx++)
    {
        a2fPoint[0] = a2fRelOffset[0] + a2fMajorMinor[0] * cos(fRotOffset * Gs_dRADIANS_PER_DEGREE) * 
                      cos(nx * 360.0 / c_nPoints * Gs_dRADIANS_PER_DEGREE) - 
                      a2fMajorMinor[1] * sin(fRotOffset*Gs_dRADIANS_PER_DEGREE) *
                      sin(nx*360.0/c_nPoints*Gs_dRADIANS_PER_DEGREE);
        
        a2fPoint[1] = a2fRelOffset[1] + a2fMajorMinor[0]* sin(fRotOffset * Gs_dRADIANS_PER_DEGREE) *
                      cos(nx*360.0 / c_nPoints * Gs_dRADIANS_PER_DEGREE) + 
                      a2fMajorMinor[1] * cos(fRotOffset*Gs_dRADIANS_PER_DEGREE) * 
                      sin(nx * 360.0 / c_nPoints * Gs_dRADIANS_PER_DEGREE);
        
        nAddQuadratic(2,1.0,1.0,&a2fPoint[0],&pfCoeffA[0][0],&pfCoeffB[0]);
    }

    a2fPoint[0] = a2fRelOffset[0];
    a2fPoint[1] = a2fRelOffset[1];
    nAddQuadratic(2,0.0,1.0,&a2fPoint[0],&pfCoeffA[0][0],&pfCoeffB[0]);
    
    double      fConst = 0.0;
    if (nSolveQuadratic(2,&a2x2lfGMat[0][0],&a2fhVec[0],&fConst,&pfCoeffA[0][0],&pfCoeffB[0]))
        return 1;
    
    // Copy the ellipsoid over so that we can use it in the future.
    vCopyVecND(4,&a2x2lfGMat[0][0],a2x2fEllipsoidAOut);
    vCopyVecND(2,&a2fhVec[0],a2fEllipsoidbOut);            
    
    *pfEllipsoidcOut = fConst;
    
    return 0;
}

DTREK_EXPORT int nBuildEllipse(double fRadius,double* a2x2fEllipsoidA, double* a2fEllipsoidb, double* pfEllipsoidc)
{
    a2x2fEllipsoidA[0] = a2x2fEllipsoidA[3] = 1.0 / ( fRadius * fRadius );
    
    a2x2fEllipsoidA[1] = a2x2fEllipsoidA[2] = 0.0;
    
    a2fEllipsoidb[0] = a2fEllipsoidb[1] = 0.0;
    
    *pfEllipsoidc = 0;
    
    return 0;
}

DTREK_EXPORT int nMergeEllipsoids(double* pfOffset1,double* a2x2fEllipsoidA1,double* a2fEllipsoidb1,double fEllipsoidc1,double* pfOffset2,double* a2x2fEllipsoidA2,double* a2fEllipsoidb2,double fEllipsoidc2,double* a2x2fEllipsoidAOut,double* a2fEllipsoidbOut,double* pfEllipsoidcOut) {

    /* Step 1:  Get ellipsoids to be relative to same offset and compute parameteric form */

    double a2fRelOffset1[2];
    double a2fRelOffset2[2];
    double a2fMajorMinor1[2];
    double a2fMajorMinor2[2];
    double fOffset1,fOffset2;
    const double fRotStep = 5.0;

    // We have different centers for each ellipsoid.  Transform ellipsoid 2 to the offset of ellipsoid 1
    vConvertEllipsoidRelative(a2x2fEllipsoidA2,a2fEllipsoidb2,&fEllipsoidc2,pfOffset2[0],pfOffset2[1],pfOffset1[0],pfOffset1[1]);


    // Having converted ellipsoid 2 to be relative to ellipsoid 1 coordinates, we use pfOffset1 as the offset for both
    // ellipsoids.

    // Calculate the parameteric form for each ellipsoid.
    if (nGetEllipsoidTiltAngleCenter(a2x2fEllipsoidA1,a2fEllipsoidb1,fEllipsoidc1,a2fRelOffset1,a2fMajorMinor1,&fOffset1))
        return 1;
    if (nGetEllipsoidTiltAngleCenter(a2x2fEllipsoidA2,a2fEllipsoidb2,fEllipsoidc2,a2fRelOffset2,a2fMajorMinor2,&fOffset2))
        return 1;
    // a2fRelOffset1 and a2fRelOffset2 give the new center (relative to pfOffset1[])

    /* Step 2:  Calculate a set of points which each lie in one but not both of the ellipsoids */

    const int nParamPoints = 20;
    double a2fPoints[nParamPoints*2][2];
    double a2fTestPoint[2];
    double a2fTestPointRot[2];
    int nPoints;
    int nPoint;
    
    nPoints = 0;
    // Step through points on each ellipsoid.  
    // A point is only used if it is outside the other ellipsoid

    for (nPoint = 0; nPoint < nParamPoints; nPoint++) {
        a2fTestPoint[0] = 
            cos(Gs_dRADIANS_PER_DEGREE*(nPoint*360.0/nParamPoints))*a2fMajorMinor1[0];
        a2fTestPoint[1] = 
            sin(Gs_dRADIANS_PER_DEGREE*(nPoint*360.0/nParamPoints))*a2fMajorMinor1[1];
        a2fTestPointRot[0] = a2fRelOffset1[0] + cos(Gs_dRADIANS_PER_DEGREE*fOffset1)*a2fTestPoint[0] - sin(Gs_dRADIANS_PER_DEGREE*fOffset1)*a2fTestPoint[1];
        a2fTestPointRot[1] = a2fRelOffset1[1] + sin(Gs_dRADIANS_PER_DEGREE*fOffset1)*a2fTestPoint[0] + cos(Gs_dRADIANS_PER_DEGREE*fOffset1)*a2fTestPoint[1];
        if (fEvalEllipse(a2fTestPointRot[0],a2fTestPointRot[1],a2x2fEllipsoidA2,a2fEllipsoidb2,&fEllipsoidc2)>1.0) {
            // The point is outside ellipsoid2.  Use it.
            a2fPoints[nPoints][0] = a2fTestPointRot[0];
            a2fPoints[nPoints][1] = a2fTestPointRot[1];
            nPoints++;
        };
    };
    for (nPoint = 0; nPoint < nParamPoints; nPoint++) {
        a2fTestPoint[0] = 
            cos(Gs_dRADIANS_PER_DEGREE*(nPoint*360.0/nParamPoints))*a2fMajorMinor2[0];
        a2fTestPoint[1] = 
            sin(Gs_dRADIANS_PER_DEGREE*(nPoint*360.0/nParamPoints))*a2fMajorMinor2[1];
        a2fTestPointRot[0] = a2fRelOffset2[0] + cos(Gs_dRADIANS_PER_DEGREE*fOffset2)*a2fTestPoint[0] - sin(Gs_dRADIANS_PER_DEGREE*fOffset2)*a2fTestPoint[1];
        a2fTestPointRot[1] = a2fRelOffset2[1] + sin(Gs_dRADIANS_PER_DEGREE*fOffset2)*a2fTestPoint[0] + cos(Gs_dRADIANS_PER_DEGREE*fOffset2)*a2fTestPoint[1];
        if (fEvalEllipse(a2fTestPointRot[0],a2fTestPointRot[1],a2x2fEllipsoidA1,a2fEllipsoidb1,&fEllipsoidc1)>1.0) {
            // The point is outside ellipsoid1.  Use it.
            a2fPoints[nPoints][0] = a2fTestPointRot[0];
            a2fPoints[nPoints][1] = a2fTestPointRot[1];
            nPoints++;
        };
    };


    /* Step 3:  Given this set of points, compute a new ellipsoid */

    double      a2fDir[2];
    double      pfCoeffA[7][7];
    double      pfCoeffB[7];
    double      a2x2lfGMat[2][2];
    double      a2fhVec[2];
    double      fConst;
    double*     pf0,*pf1;

    
    // Calculate the new ellipsoid.
    pf0 = & pfCoeffA[0][0];
    pf1 = & pfCoeffB[0];
    nInitQuadratic(2,FALSE,&pf0,&pf1);

    for (nPoint = 0; nPoint < nPoints; nPoint++)
        nAddQuadratic(2,1.0,1.0,&a2fPoints[nPoint][0],&pfCoeffA[0][0],&pfCoeffB[0]);

    a2fDir[0] = 0.0;
    a2fDir[1] = 0.0;
    nAddQuadratic(2,0.0,1.0,&a2fDir[0],&pfCoeffA[0][0],&pfCoeffB[0]);
    if (nSolveQuadratic(2,&a2x2lfGMat[0][0],&a2fhVec[0],&fConst,&pfCoeffA[0][0],&pfCoeffB[0]))
        return 1;
    // Copy the ellipsoid over so that we can use it in the future.
    vCopyVecND(4,&a2x2lfGMat[0][0],a2x2fEllipsoidAOut);
    vCopyVecND(2,&a2fhVec[0],a2fEllipsoidbOut);            
    *pfEllipsoidcOut = fConst;


    return 0;
};

DTREK_EXPORT void vConvertEllipsoidRelative(float* a2x2fEllipsoidAIn_,float* a2fEllipsoidb_,float* pfEllipsoidc_,float fCentIn0,float fCentIn1,float fCentOut0,float fCentOut1) {
    double a2x2fEllipsoidAIn[2][2];
    double a2fEllipsoidb[2];
    double fEllipsoidc;

    vCopyVecND(4,a2x2fEllipsoidAIn_,&a2x2fEllipsoidAIn[0][0]);
    vCopyVecND(2,a2fEllipsoidb_,&a2fEllipsoidb[0]);
    fEllipsoidc = *pfEllipsoidc_;

    vConvertEllipsoidRelative(&a2x2fEllipsoidAIn[0][0],&a2fEllipsoidb[0],&fEllipsoidc,fCentIn0,fCentIn1,fCentOut0,fCentOut1);
    vCopyVecND(2,&a2fEllipsoidb[0],a2fEllipsoidb_);
    *pfEllipsoidc_ = fEllipsoidc;

};

DTREK_EXPORT void vConvertEllipsoidRelative(double* a2x2fEllipsoidAIn,double* a2fEllipsoidb,double* pfEllipsoidc,double fCentIn0,double fCentIn1,double fCentOut0,double fCentOut1) {
    float a2fEllipsoidbOut[2];
    float a2fShift[2];
    float a2x2fEllipsoidA[2][2];

    a2x2fEllipsoidA[0][0] = a2x2fEllipsoidAIn[0];
    a2x2fEllipsoidA[0][1] = a2x2fEllipsoidAIn[1];
    a2x2fEllipsoidA[1][0] = a2x2fEllipsoidAIn[2];
    a2x2fEllipsoidA[1][1] = a2x2fEllipsoidAIn[3];

    a2fShift[0] = fCentIn0 - fCentOut0;
    a2fShift[1] = fCentIn1 - fCentOut1;
    a2fEllipsoidbOut[0] = -2.0*(a2x2fEllipsoidA[0][0]*a2fShift[0] + a2x2fEllipsoidA[1][0]*a2fShift[1]) + a2fEllipsoidb[0];
    a2fEllipsoidbOut[1] = -2.0*(a2x2fEllipsoidA[1][0]*a2fShift[0] + a2x2fEllipsoidA[1][1]*a2fShift[1]) + a2fEllipsoidb[1];
    *pfEllipsoidc = -a2fEllipsoidb[0]*a2fShift[0] -a2fEllipsoidb[1]*a2fShift[1] + 
        a2x2fEllipsoidA[0][0]*a2fShift[0]*a2fShift[0] +
        a2x2fEllipsoidA[1][0]*a2fShift[1]*a2fShift[0] +
        a2x2fEllipsoidA[0][1]*a2fShift[0]*a2fShift[1] +
        a2x2fEllipsoidA[1][1]*a2fShift[1]*a2fShift[1] + *pfEllipsoidc;
    a2fEllipsoidb[0] = a2fEllipsoidbOut[0];
    a2fEllipsoidb[1] = a2fEllipsoidbOut[1];
};


DTREK_EXPORT double fCalcRotationOffset(double a3fE3[3],double a3fS0[3],double a3fXRC[3],double* a3fXRCR) {
    
    
    // We want to discover the amount of rotation necc. to bring XRC
    // back onto the Ewald sphere.
    // a3fXRO contains the 'ideal' Ewald sphere position (calculated from (x,y) points on plate)
    // a3fXRC contains the 'calculated' Ewald sphere position (probably NOT on the sphere).
    
    double a3fDStarV1[3];
    double a3fDStarV2[3];
    double a3fDStarV3[3];
    double fk1,fk2,fb;
    double a2fSolutions[2];

    double a3fTempVec1[3];
    int nx,ny;
    double f0,f1;
    double fRotErr;
    double a3x3fTempMat1[3][3];
    
    // Find the three vectors v1,v2 and v3.
    // Project the offset d* vector of the spot down into the rotation plane to get v1
    vMulVec3DScalar(a3fE3,fDot3D(a3fE3,a3fXRC),a3fTempVec1);
    vSubVec3DVec3D(a3fXRC,a3fTempVec1,a3fDStarV1);
    // Take cross product of m_a3fRotVector and v1 to get v2.
    // This gives v1 cross v2 == m_a3fRotVector
    vCross3D(a3fE3,a3fDStarV1,a3fDStarV2);
    // Calculate v3
    vSubVec3DVec3D(a3fXRC,a3fDStarV1,a3fDStarV3);
    
    // Calculate constants a (int a3fTemp1), k1,k2 and b
    
    vSubVec3DVec3D(a3fDStarV3,a3fS0,a3fTempVec1);
    fb = fDot3D(a3fDStarV1,a3fDStarV1) + fDot3D(a3fTempVec1,a3fTempVec1) - 1.0;
    fk1 = 2.0* fDot3D(a3fDStarV1,a3fTempVec1);
    fk2 = 2.0* fDot3D(a3fDStarV2,a3fTempVec1);
    
    // Solve quadratic system.  Choose solution with the lowest absolute value.
    nx = nSolveQuadratic(fk1*fk1 + fk2*fk2, fk2*fb,fb*fb - fk1*fk1,a2fSolutions);
    if (nx) {
        f0 = 1000.0;
        for (ny=0;ny<nx;ny++) {
            if ((a2fSolutions[ny]>1.0) || (a2fSolutions[ny]<-1.0)) {
                break;
            };
            f1 = asin(a2fSolutions[ny]);
            if (ABS(f1)<ABS(f0))
                f0 = f1;
        };
    };
    if ((nx==0) || (ny !=nx) || (f0==1000.0)) {
        return -9999.0;
    } else
        fRotErr = f0/fRADIANS_PER_DEGREE;

    if (a3fXRCR) {
        // Calculate the rotated X-calc
        vConvRotVec3DMat3D(fRotErr,&a3fE3[0],&a3x3fTempMat1[0][0]);
        vMulMat3DVec3D(a3x3fTempMat1,a3fXRC,a3fXRCR);

    };
    
    return fRotErr;
    
};


DTREK_EXPORT int nConvexHull(int nNumPoints,double* aa2fPoints,int* anPointsInHull) {
    // This is *not* a very good convex hull algorithm, but it is simple,
    // and works well for well behaved points (i.e. a set of points who's distribution
    // is roughly circular in shape).  "needlelike" distributions will probably not work as well.

    double a2fNorm[2];
    double fMaxNorm;
    double fAngularStep;
    double fStep;
    double f0;
    int    nMaxNorm;
    int    nx;
    int    nPoint;

    // Mark all points as outside the hull
    for (nx=0;nx<nNumPoints;nx++)
        anPointsInHull[nx] = 0;
    fAngularStep = 360.0/(nNumPoints*2);
    for (fStep = 0.0; fStep < 360.0; fStep += fAngularStep) {
        a2fNorm[0] = cos(fStep*Gs_dRADIANS_PER_DEGREE);
        a2fNorm[1] = sin(fStep*Gs_dRADIANS_PER_DEGREE);
        fMaxNorm = -1.0;
        nMaxNorm = -1;
        for (nPoint = 0; nPoint < nNumPoints;nPoint++) {
            f0 = a2fNorm[0]*aa2fPoints[nPoint*2 + 0] + a2fNorm[1]*aa2fPoints[nPoint*2 + 1];
            if (f0 > fMaxNorm) {
                fMaxNorm = f0;
                nMaxNorm = nPoint;
            };
        };
        if (fMaxNorm > 0.0) 
            anPointsInHull[nMaxNorm] = 1;
    };
    return 0;
};

int nPlot(const char* pcName,itr<double>& afY,itr<double>* pafX) {
    FILE* pFOut;
    int nx;
    static Cstring sLastName;
    char* pcNameIn;
    char* pcAttrib;

    if ((pcName == NULL) || (!pcName[0])) {
        pcNameIn = sLastName.string();
        pcAttrib = "at";
    } else {
        pcNameIn = (char*) pcName;
        pcAttrib = "wt";
        sLastName = pcName;
    }

    pFOut = fopen(pcNameIn,pcAttrib);

    if (!pFOut)
        return 1;

    if ((pcName == NULL) || (!pcName[0])) {
        fprintf(pFOut,"\n");
    };
    
    for (nx = 0; nx < afY.size(); nx++) 
        fprintf(pFOut,"%lf %lf\n",(double) ((pafX)?((*pafX)[nx]):nx),afY[nx]);
    fclose(pFOut);
    return 0;

};

int nPlot(const char* pcName,itr<double>& afY,double fStartX,double fStepX) {
    FILE* pFOut;
    int nx;
    static Cstring sLastName;
    char* pcNameIn;
    char* pcAttrib;

    if ((pcName == NULL) || (!pcName[0])) {
        pcNameIn = sLastName.string();
        pcAttrib = "at";
    } else {
        pcNameIn = (char*) pcName;
        pcAttrib = "wt";
        sLastName = pcName;
    }

    pFOut = fopen(pcNameIn,pcAttrib);
    if (!pFOut)
        return 1;
    for (nx = 0; nx < afY.size(); nx++) 
        fprintf(pFOut,"%lf %lf\n",(double)  fStartX + nx*fStepX,afY[nx]);
    fclose(pFOut);
    return 0;

};


int nInvFunc(double fStartX,
             double fStepX,
             itr<double>& afFunctionY,
             double fInvStartY,
             double fInvEndY,
             double fInvStepY,
             itr<double>& afInvFunctionX,
             bool bForceMonotonic)
{
    double fMinY,fMaxY;
    double fXForMinY;
    double fXForMaxY;
    int    nLess,nGreater;
    int    nY;
    int    nPointsInv,nPoints;
    double fThisY,fLastY,fY;
    double* pfInvFunctionX;
    double* pfFunctionY;
    int    nBoundLower,nBoundUpper;

    if (!afFunctionY.size())
        return 1;

    -afInvFunctionX;
    fMinY = fMaxY = afFunctionY[0];
    fXForMinY = fXForMaxY = fStartX;
    if (afFunctionY.last() < fMinY) {
        fMinY = afFunctionY.last();
        fXForMinY = (afFunctionY.size()-1)*fStepX + fStartX;
    };
    if (afFunctionY.last() > fMaxY) {
        fMaxY = afFunctionY.last();
        fXForMaxY = (afFunctionY.size()-1)*fStepX + fStartX;
    };

    
    pfFunctionY = &afFunctionY[0];
    nPoints = afFunctionY.size();
    nLess = 0;
    nGreater = 0;
    fLastY = pfFunctionY[0];
    for (nY = 1; nY < nPoints; nY++) {
        fThisY = pfFunctionY[nY];
        if (fThisY > fLastY)
            nGreater++;
        else if (fThisY < fLastY)
            nLess++;

        fLastY = fThisY;
    };

    // Make sure we are monotonic.
    if ((nGreater > 0) && (nLess > 0)) {
        if (bForceMonotonic) {
            fLastY = pfFunctionY[0];
            for (nY = 1; nY < nPoints; nY++) {
                fThisY = pfFunctionY[nY];
                if (nGreater > nLess) {
                    if (fThisY < fLastY)
                        fThisY = pfFunctionY[nY] = fLastY;
                } else {
                    if (fThisY > fLastY)
                        fThisY = pfFunctionY[nY] = fLastY;
                };
                fLastY = fThisY;
            };            
        } else
            return 1;
    };
    

    // Now verify the output ranges.
    // If fInvStartY > fInvEndY it is a sign that we need to automatically adjust the end range in computation.
    if (fInvStartY > fInvEndY)
        fInvEndY = max(fInvStartY,fMaxY);
    fInvEndY = fInvStartY + fInvStepY*(nPointsInv = (int) floor(0.5 + (fInvEndY - fInvStartY)/fInvStepY));

    if (nGreater) {
        // Monotonic increasing.  
        nBoundLower = 0;
        nBoundUpper = min(nPoints-1,1);

    } else {
        // Monotonic decreasing.
        nBoundLower = nPoints - 1;
        nBoundUpper = max(0,nPoints - 2);
    };

    afInvFunctionX.setsize(nPointsInv);
    pfInvFunctionX = &afInvFunctionX[0];
    for (fY = fInvStartY; ABS(fY - fInvEndY)>fInvStepY/10;fY+= fInvStepY) {
        if (fY <= fMinY) {
            // Just use the minimum value.
            *pfInvFunctionX = fXForMinY;
        } else if (fY >= fMaxY) {
            *pfInvFunctionX = fXForMaxY;
        } else {
            // Find the bounding range.
            // This MUST terminate.
            while ((fY < pfFunctionY[nBoundLower]) || (fY > pfFunctionY[nBoundUpper])) {
                if (nGreater) {
                    if (nBoundUpper == nPoints - 1)
                        return 1;
                    nBoundLower++;
                    nBoundUpper++;
                } else {
                    if (nBoundLower == 0)
                        return 1;
                    nBoundLower--;
                    nBoundUpper--;
                };
            };
            // Compute the averaged value.
            *pfInvFunctionX = fStartX + fStepX*((nBoundLower)*(pfFunctionY[nBoundUpper] - fY) +
                (nBoundUpper)*(fY - pfFunctionY[nBoundLower]))/(pfFunctionY[nBoundUpper] - pfFunctionY[nBoundLower]);
        };
        pfInvFunctionX++;
    };

    return 0;
};



int nCreateMapping(itr<int>& anMap,int nNumEntries) {
    int nx;
    anMap.setsize(nNumEntries + 1);
    anMap[0] = nNumEntries;     // Store the size of the map in the first entry.
    for (nx = 1; nx < nNumEntries + 1; nx++)
        anMap[nx] = -1;
    return 0;
};

int nAddToMapping(itr<int>& anMap,int nFrom,int nTo,bool bPruneDuplicates) {
    int nNext;
    if ((!anMap.size()) || (nFrom<0) || (nFrom >= anMap[0]))
        return 1;
    nNext = anMap[nFrom + 1];
    if (bPruneDuplicates) {
        while (nNext != -1) {
            if (anMap[nNext + 1] == nTo)
                return 0;
            nNext = anMap[nNext];
        };
    };

    anMap[nFrom + 1] = anMap.size();
    anMap + nNext;
    anMap + nTo;
    return 0;
};

int nGetMapping(itr<int>& anMap,int nFrom,itr<int>& anTo,bool bFIFO) {
    int nNext;
    int nx;
    if ((!anMap.size()) || (nFrom<0) || (nFrom >= anMap[0]))
        return 1;
    -anTo;
    nNext = anMap[nFrom + 1];
    while (nNext != -1) {
        anTo + anMap[nNext + 1];
        nNext = anMap[nNext];
    };
    if (bFIFO) {
        for (nx = 0; nx < anTo.size()/2; nx++)
            std::swap(anTo[nx],anTo[anTo.size() - 1 - nx]);
    };
    return 0;
};


int nDeconvMat3DVec3DRot1(double fMatrix[3][3],double* pfVector,double* pfRot) {
    double a3x3fReduced[3][3];
    double a3fVectorBefore[3];
    double a3fVectorAfter[3];
    double a3fVectorCross[3];
    double a3fTempVec1[3];
    double a3fTempVec2[3];
    int nPivot;
    int nx,ny;
    double fMax;
    double fMult;

    // Determine the absolute maximum value.  This is useful so that we don't have a problem pivoting.
    fMax = 0.0;
    for (nx = 0; nx < 3; nx++) {
        for (ny = 0; ny < 3; ny++) {
            fMax = max(ABS(fMatrix[nx][nx]),fMax);
        };
    };
    if (fMax < 1e-10)
        return 1;

    // We know that the matrix has an eignvalue of zero if it is a rotation matrix.
    vCopyMat3D(&fMatrix[0][0],&a3x3fReduced[0][0]);
    for (nx = 0; nx < 3; nx++)
        a3x3fReduced[nx][nx] -= 1.0;
    nPivot = 0;
    for (nx = 1;nx < 3; nx++) {
        if (ABS(a3x3fReduced[0][nx])>ABS(a3x3fReduced[0][nPivot]))
            nPivot = nx;
    };
    if (ABS(a3x3fReduced[0][nPivot])/fMax < 1e-10) {
        // The pivot vector is identically zero.
        pfVector[0] = 1.0;
        pfVector[1] = 0.0;
        pfVector[2] = 0.0;
        a3fVectorBefore[0] = 0.0;
        a3fVectorBefore[1] = 1.0;
        a3fVectorBefore[2] = 0.0;

    } else {
        // We can pivot on nPivot.  Reduce the other two rows.
        fMult = -a3x3fReduced[0][(nPivot + 1) % 3]/a3x3fReduced[0][nPivot];
        for (ny = 0; ny < 3;ny++) {
            a3fTempVec1[ny] = a3x3fReduced[ny][(nPivot + 1) % 3] + a3x3fReduced[ny][nPivot]*fMult;
            a3fTempVec2[ny] = a3x3fReduced[ny][nPivot];
        };
        vCross3D(a3fTempVec1,a3fTempVec2,pfVector);
        fNormVec3D(pfVector);
        vCross3D(a3fTempVec1,pfVector,a3fVectorBefore);
    };
    
    fNormVec3D(a3fVectorBefore);
    vMulMat3DVec3D(fMatrix,a3fVectorBefore,a3fVectorAfter);
    fNormVec3D(a3fVectorAfter);
    vCross3D(a3fVectorBefore,a3fVectorAfter,a3fVectorCross);
    
    *pfRot = acos(max(-1.0,min(1.0,fDot3D(a3fVectorAfter,a3fVectorBefore))))/ Gs_dRADIANS_PER_DEGREE;
    if (fDot3D(a3fVectorCross,pfVector)<0.0)
        *pfRot = 360.0 - *pfRot;

    return 0;
};


/*  This function computes the upper or lower envelope for a function.  It starts at the minimum or maximum value, and builds a montonic
    increasing or decreasing function on that range that bounds the input function.  The range of application can be restricted by the user.
    This function is useful to find envelopes of noise, and therefore to construct an "averaged" function that is monotonic.

*/

int nUpperLowerEnvelope(int nMode,itr<double>& afDataYIn,itr<double>& afDataYOut,int nRangeLower,int nRangeUpper,itr<double>* pafDataXIn) {
    int nPointer,nInterpolate,nStartMax;
    int nLocalMin,nLocalMax,nLastMax;
    int nDir;
    int nModeSign;
    double fInterpolate0,fInterpolate1;
    double* pfDataYIn;
    double* pfDataYOut;
    double* pfDataXIn;

    if (!afDataYIn.size())
        return 0;

    pfDataYIn = &afDataYIn[0];
    afDataYOut.setsize(afDataYIn.size());
    pfDataYOut = &afDataYOut[0];
    pfDataXIn = (pafDataXIn)?(&(*pafDataXIn)[0]):NULL;
  
    // If the upper range was not specified, use the whole array.
    if (nRangeUpper == -1)
        nRangeUpper = afDataYIn.size() - 1;
    nStartMax = nRangeLower;
    // Set the mode sign.  The algorithm below is coded to find the upper envelope starting at the absolute maximum in [nRangeLower,nRangeUpper]
    // So in problems where we want to find the lower envelope starting at the absolute min, we simply apply a sign to each data point using "nModeSign"
    nModeSign = (nMode == g_nFindApplyUpperEnvelope)?1:-1;
    
    // Find the optimum value, which will be the minimum or maximum depending upon which mode we are in.
    for (nPointer = nRangeLower + 1; nPointer <= nRangeUpper;nPointer++) {
        if (pfDataYIn[nPointer]*nModeSign > pfDataYIn[nStartMax]*nModeSign)
            nStartMax = nPointer;
    };
    
    pfDataYOut[nStartMax] = pfDataYIn[nStartMax];
    // We go forward and backwards from the starting point.  nDir sets the direction "forward".
    for (nDir = -1; nDir <=1; nDir +=2 ) {
        // We start one slot to the right (nDir == 1) or to the left (nDir == -1).
        for (nPointer = nStartMax + nDir; (nPointer >=nRangeLower) && (nPointer <= nRangeUpper) ;nPointer += nDir) {
            // Searching for the first value that does not follow the progression.
            // Remember to access pfYDataIn since we have not yet built values for pfYDataOut on this range yet.
            if (pfDataYIn[nPointer]*nModeSign > pfDataYIn[nPointer-nDir]*nModeSign) {
                // We found a local min 
                // Track forward to find the next local max 
                nLocalMin = nPointer - nDir;
                for (nLocalMax = nLocalMin + nDir; (pfDataYIn[nLocalMax]*nModeSign > pfDataYIn[nLocalMax - nDir]*nModeSign) && (nLocalMax>=nRangeLower) && (nLocalMax<= nRangeUpper);nLocalMax+=nDir);
                nLocalMax -= nDir;
                // Copy all intermediate values.
                while (1) {
                    pfDataYOut[nPointer] = pfDataYIn[nPointer];
                    if (nPointer == nLocalMax)
                        break;
                    nPointer+= nDir;
                };
                
                // Now track backward and find the associated point that is greater.
                // We now must use pfDataYOut because past values might have been modified by a previous step.
                // Note that we *should* find a value here becuase by definition, we started at the minimum (or maximum) value, and are at some point lower than that.
                for (nLastMax = nLocalMin - nDir; (pfDataYOut[nLastMax]*nModeSign <  pfDataYOut[nLocalMax]*nModeSign) && (nLastMax>=nRangeLower) && (nLastMax<= nRangeUpper);nLastMax-=nDir);
                
                fInterpolate0 = pfDataYOut[nLastMax];
                fInterpolate1 = pfDataYOut[nLocalMax];
                // Now interplate all of the points between the last maximum and this.
                for (nInterpolate = nLastMax; nInterpolate != nLocalMax + nDir; nInterpolate+= nDir) {
                    if (pfDataXIn) 
                        pfDataYOut[nInterpolate] = (fInterpolate0*ABS(pfDataXIn[nInterpolate] - pfDataXIn[nLocalMax]) + fInterpolate1*ABS(pfDataXIn[nInterpolate] - pfDataXIn[nLastMax]))/(ABS(pfDataXIn[nLocalMax] - pfDataXIn[nLastMax]));
                    else
                        pfDataYOut[nInterpolate] = (fInterpolate0*ABS(nInterpolate - nLocalMax) + fInterpolate1*ABS(nInterpolate - nLastMax))/(ABS(nLocalMax - nLastMax));
                };
                nPointer = nLocalMax;
                
            } else
                pfDataYOut[nPointer] = pfDataYIn[nPointer];
        };
    };

    // Copy the rest of the values that were not in the range of application.
    for (nPointer = 0; nPointer < nRangeLower; nPointer++) 
        pfDataYOut[nPointer] = pfDataYIn[nPointer];
    for (nPointer = afDataYIn.size() - 1; nPointer > nRangeUpper; nPointer--)
        pfDataYOut[nPointer] = pfDataYIn[nPointer];
    
    return 0;
};








#define nCopyNewVector(nBest,fValue) {\
    vSubVecNDVecND(m_nSimplexPoints,&m_afParamPointSum[0],&m_aafParamPointVectors[nBest*m_nSimplexPoints + 0],&m_afParamPointSum[0]);\
    vAddVecNDVecND(m_nSimplexPoints,&m_afParamPointSum[0],&afDelta[0],&m_afParamPointSum[0]);\
    vCopyVecND(m_nSimplexPoints,&afDelta[0],&m_aafParamPointVectors[nBest*m_nSimplexPoints + 0]);\
    m_afParamPointValues[nBest] = fValue; \
    if (nVerbose)\
        printf("Point %2d is better. (%.2lf)\n",nBest,fValue);\
    };\
    bPointWorked  = true;


itr<double> m_aafParamPointVectors; // [nSIMPLEXPARAMS+1][nSIMPLEXPARAMS]
itr<double> m_afParamPointValues;   // [nSIMPLEXPARAMS+1]
itr<double> m_afParamPointSum;      // [nSIMPLEXPARAMS]
itr<double> m_afParamGoodPoint;     // [nSIMPLEXPARAMS]
int m_nSimplexState1;
int m_nSimplexState2;
int m_nSimplexPoints;


int nSimplexRun(double fValue,itr<double>& afDelta,itr<double>* pafStep) {

    int nVerbose = 0;
   
    if (pafStep) {
        m_nSimplexState1 = 0;
        m_nSimplexState2 = -1;
        // Count the number of free variables available.
        m_nSimplexPoints = pafStep->size();
        m_aafParamPointVectors.setsize((m_nSimplexPoints+1)*m_nSimplexPoints);
        m_afParamPointValues.setsize(m_nSimplexPoints+1);
        m_afParamPointSum.setsize(m_nSimplexPoints);
        m_afParamGoodPoint.setsize(m_nSimplexPoints);
    };




    int nx,ny;
    int nBest;
    int nWorst,nSecondWorst;
    double fFac,fFac1,fFac2;
    bool bPointWorked;

    if (m_nSimplexState1==0) {
        if (m_nSimplexState2==-1) {
            // Initialize the points in the simplex.
            for (nx=0;nx<= m_nSimplexPoints; nx++) {
                vZeroMat(m_nSimplexPoints,1,&m_aafParamPointVectors[nx*m_nSimplexPoints + 0]);
                for (ny=0;ny< m_nSimplexPoints; ny++) 
                    m_aafParamPointVectors[nx*m_nSimplexPoints + ny]=((rand() % 2)*2 - 1)*(*pafStep)[ny];
            };
            m_nSimplexState2 = 0;
            if (nVerbose)
                printf("Initializing.\n");
        } 
        if (nVerbose)
            printf("Initializing Vertex %d.\n",m_nSimplexState2);


        m_afParamPointValues[m_nSimplexState2] = fValue;
        m_nSimplexState2++;

        
        if (m_nSimplexState2==m_nSimplexPoints) {
            // Initialize for the computation steps.
            m_nSimplexState1 = 1;
            m_nSimplexState2 = 0;
            // Compute the sum of the points.
            vZeroMat(m_nSimplexPoints,1,&m_afParamPointSum[0]);
            for (nx=0;nx<m_nSimplexPoints;nx++) {
                vAddVecNDVecND(m_nSimplexPoints,&m_aafParamPointVectors[nx*m_nSimplexPoints + 0],&m_afParamPointSum[0],&m_afParamPointSum[0]);
            };
        } else
            vCopyVecND(m_nSimplexPoints,&m_aafParamPointVectors[m_nSimplexState2*m_nSimplexPoints + 0],&afDelta[0]);
    };
    while (m_nSimplexState1==1) {
        // Calculate Worst, SecondWorst and Best
        if (m_afParamPointValues[0]>m_afParamPointValues[1]) {
            nWorst = 0;
            nSecondWorst = 1;
        } else {
            nWorst = 1;
            nSecondWorst = 0;
        };
        nBest = 0;
        for (nx=1;nx<m_nSimplexPoints;nx++) {
            if (m_afParamPointValues[nx]<m_afParamPointValues[nBest])
                nBest = nx;
            if (m_afParamPointValues[nx]>m_afParamPointValues[nWorst]) {
                nSecondWorst = nWorst;
                nWorst = nx;
            };
        };
        bPointWorked = false;

        // Make sure that all the points are not bad.  If so, then we will restore the
        // zero point solution and terminate.
        for (nx=0;nx<m_nSimplexPoints;nx++) {
            if (m_afParamPointValues[nx]<1e10)
                break;
        };
        if (nx==m_nSimplexPoints) {
            vZeroMat(m_nSimplexPoints,1,&afDelta[0]);
            return 2;
        };

        // Check for termination conditions.
        if (ABS(m_afParamPointValues[nWorst] - m_afParamPointValues[nBest])/
            (ABS(m_afParamPointValues[nWorst]) + ABS(m_afParamPointValues[nBest]))<0.001) {
            vCopyVecND(m_nSimplexPoints,&m_aafParamPointVectors[0],&afDelta[0]);
            return 1;
        };
        
        // Start of simplex algorithm.
        
        
        fFac = 0.0;
        switch (m_nSimplexState2) {
        case 0:
            // Try a factor -1.0
            fFac = -1.0;
            m_nSimplexState2 = 1;
            
            break;
        case 1:
            // Just tried -1.0.  
            if (fValue<m_afParamPointValues[nWorst]) {
                nCopyNewVector(nWorst,fValue);
                
            };
            if (fValue<m_afParamPointValues[nBest]) {
                // Did Better than the current best.
                // try -2.0.
                fFac = -2.0;
                m_nSimplexState2 = 2;
            } else if (fValue>m_afParamPointValues[nSecondWorst]) {
                // Did Worse than the second worst point.
                fFac = 0.5;
                m_nSimplexState2 = 3;
            } else {
                // Avoid an infinite loop.
                if (!bPointWorked)
                    return 2;
                m_nSimplexState2 = 0;
                fFac = 0.0;
            };
            
            break;
        case 2:
            // Just tried -2.0.
            if (fValue<m_afParamPointValues[nWorst]) {
                nCopyNewVector(nWorst,fValue);
            };
            m_nSimplexState2 = 0;
            fFac = 0.0;
        case 3:
            // Just tried 0.5
            if (fValue<m_afParamPointValues[nWorst]) {
                nCopyNewVector(nWorst,fValue);
                m_nSimplexState2 = 0;
                fFac = 0.0;
            } else {
                // Can't get rid of the high point.  
                // Do a contraction.
                
                if (nVerbose)
                    printf("Contracting.\n");
                for (nx=0;nx<m_nSimplexPoints;nx++) {
                    if (nx!=nBest) {
                        for (ny=0;ny<m_nSimplexPoints;ny++) {
                            m_aafParamPointVectors[nx*m_nSimplexPoints + ny] = 
                                0.5*(m_aafParamPointVectors[nx*m_nSimplexPoints + ny] + m_aafParamPointVectors[nBest*m_nSimplexPoints + ny]);
                        };
                    };
                };
                // Restart, since we need to recompute values for most vectors.
                m_nSimplexState1 = 0;
                m_nSimplexState2 = 0;
                fFac = 0.0;
                // Copy over first solution.
                vCopyVecND(m_nSimplexPoints,&m_aafParamPointVectors[0],&afDelta[0]);
            };
        };

        if (fFac == 0.0)
            continue;
        else {
            if (nVerbose)
                printf("Trying Factor=%.1lf\n",fFac);

            // Compute the new test point best on fFac.
            fFac1 = (1.0 - fFac)/(m_nSimplexPoints-1);
            fFac2 = fFac1 - fFac;
            for (nx=0;nx<m_nSimplexPoints;nx++) {
                afDelta[nx] = m_afParamPointSum[nx]*fFac1 - m_aafParamPointVectors[nWorst*m_nSimplexPoints + nx]*fFac2;
            };
            return 0;
        };
    };
    return 0;
};

int m_nMultiOptimizeBatchSize;
int m_nMultiOptimizeBestToPickInBatch;
int m_nMultiOptimizeNumTrials;
int m_nMultiOptimizenNumInitialSolutions;
int m_nMultiOptimizeNumParams;
bool m_bMultiOptimizeMinimize;
bool m_bMultiOptimizePrint = false;
bool m_bMultiOptimizeInitialized = false;
int m_nMultiOptimizeCount;
int m_nMultiOptimizeTrialCount;
itr<double> m_afMultiOptimizeLastParamsOut;
itr<double> m_afMultiOptimizeBestValues;
itr<double> m_afMultiOptimizeBestParams;
itr<double> m_afMultiOptimizeLastBestParams;
itr<double> m_afMultiOptimizeDistWeighting;
double m_fMultiOptimizeMaxDist = 0.0;

/*  This function should be called once initialy to load starting point solutions and set parameters.
    It must then be called in a while() loop until it returns 1.  At the last return, the parameter values of 
    the best solution are returned in afParams. and fValue.

    Here is the way the optimization works:
    
    Initialize trial solutions so that we have nNumInitialSolutions in m_afMultiOptimizeLastBestParams
    for (Trial from 1 to nNumTrials) {
        Initialize local "best values" array m_afMultiOptimizeBestParams
        for (each solution we had pending from the last trial) {
            for (Batch from 1 to nBatchSize) {
                Store new solution into m_afMultiOptimizeBestParams if the solution is better than the worst solution.
            }
        }
        Store the best nBestToPickInBatch solutions into m_afMultiOptimizeLastBestParams for the next trial
    }

    Notes:
    m_nMultiOptimizeCount takes on values from 0 to NumPendingTrialSolutions*BatchSize

*/

int nMultiOptimize(double& fValue,itr<double>& afParams,bool bMinimize,int nBatchSize,int nBestToPickInBatch,int nNumTrials,int nNumInitialSolutions) {
    int nx,ny,nz;

    if (nBatchSize >= 1) {
        // Initialize the system.
        m_nMultiOptimizeBatchSize = nBatchSize;
        m_nMultiOptimizeBestToPickInBatch = nBestToPickInBatch;
        m_nMultiOptimizeNumTrials = nNumTrials;
        m_nMultiOptimizenNumInitialSolutions = nNumInitialSolutions;

        // Consistency check.
        if ((m_nMultiOptimizeBestToPickInBatch > nBatchSize) || (m_nMultiOptimizeBatchSize<=0) || (m_nMultiOptimizeBestToPickInBatch<=0))
            return 2;
        m_nMultiOptimizeCount = -1;
        m_nMultiOptimizeTrialCount = 0;

        // The user might have provided an initial weighting vector for distance comparisons.
        m_afMultiOptimizeLastBestParams.clear();
        if (afParams.size()) {
            m_afMultiOptimizeDistWeighting + afParams;
        } else {
            -m_afMultiOptimizeDistWeighting;
        };

        m_afMultiOptimizeBestValues.clear();
        m_afMultiOptimizeBestParams.clear();
        m_bMultiOptimizeMinimize = bMinimize;
        m_bMultiOptimizeInitialized = true;
        return 0;
    };
    if ((!m_bMultiOptimizeInitialized) || (!afParams.size()))
        return 2;

    if (m_afMultiOptimizeLastBestParams.size()/afParams.size() < m_nMultiOptimizenNumInitialSolutions) {
        // We are still loading solutions.
        m_nMultiOptimizeNumParams = afParams.size();
        m_afMultiOptimizeLastBestParams + afParams;
        if (!m_afMultiOptimizeDistWeighting.size()) {
            for (nx = 0; nx < m_nMultiOptimizeNumParams;nx++)
                m_afMultiOptimizeDistWeighting + 1.0;
        };
        m_afMultiOptimizeLastParamsOut.setsize(m_nMultiOptimizeNumParams);
        return 0;
    };

    // On the very very first pass through, we don't have data yet.  This corresponds to the very first
    // while (nMultiOptimize()) call, so we have not gone through the body of the loop yet where values are computed.
    if (m_nMultiOptimizeCount>=0) {
        int nToPrint;
    
        // Make sure afParams is of the correct size.
        if ((m_nMultiOptimizeNumParams<=0) || (afParams.size() % m_nMultiOptimizeNumParams) || (m_afMultiOptimizeDistWeighting.size() % m_nMultiOptimizeNumParams))
            return 2;

        // Save the last solution.  The only time this does not happen is on the
        // very first batch of the very first trial.
        // Search for the worst value.
        // However, we also check if the worst solution is within a tolerable distance to our present solution.
        // If not, then we only over-write the worst solution if no other solutions were within a tolerable distance.
        
        double fThisDist;
        bool   bAcceptDistFound = false;
        int    nWorstPick = 0;   
        int    nWorstPickAndAcceptDist = -1;
        
        for (nx = 0; nx < m_afMultiOptimizeBestValues.size(); nx++) {
            
            
            fThisDist = 0.0;
            for (nz = 0; nz < m_nMultiOptimizeNumParams; nz++)
                fThisDist += 
                (afParams[nz] - m_afMultiOptimizeBestParams[nx*m_nMultiOptimizeNumParams + nz])*
                (afParams[nz] - m_afMultiOptimizeBestParams[nx*m_nMultiOptimizeNumParams + nz])/
                (m_afMultiOptimizeDistWeighting[nz]*m_afMultiOptimizeDistWeighting[nz]);
            fThisDist = sqrt(fThisDist);
            
            if (fThisDist < m_fMultiOptimizeMaxDist)
                bAcceptDistFound = true;
            
            if (nx == 0) {
                nWorstPick = nx;
                if (fThisDist < m_fMultiOptimizeMaxDist)
                    nWorstPickAndAcceptDist = nx;
            } else if (m_bMultiOptimizeMinimize) {
                if (m_afMultiOptimizeBestValues[nx] > m_afMultiOptimizeBestValues[nWorstPick]) {
                    nWorstPick = nx;
                };
                if ((fThisDist < m_fMultiOptimizeMaxDist) && ((nWorstPickAndAcceptDist==-1) || (m_afMultiOptimizeBestValues[nx] > m_afMultiOptimizeBestValues[nWorstPickAndAcceptDist]))) {
                    nWorstPickAndAcceptDist = nx;
                };
            } else {
                if (m_afMultiOptimizeBestValues[nx] < m_afMultiOptimizeBestValues[nWorstPick]) {
                    nWorstPick = nx;
                };
                if ((fThisDist < m_fMultiOptimizeMaxDist) && ((nWorstPickAndAcceptDist==-1) || (m_afMultiOptimizeBestValues[nx] < m_afMultiOptimizeBestValues[nWorstPickAndAcceptDist]))) {
                    nWorstPickAndAcceptDist = nx;
                };
            };
        };
        
        // Act based on results of above loop.
        
        if ((m_afMultiOptimizeBestValues.size() < m_nMultiOptimizeBestToPickInBatch) && (!bAcceptDistFound)) {
            // Case 1)   Our solution is outside the radius of any of the other solutions AND we have not
            //           filled our quota of 'best' solutions.  Use this solution.
            //           This will also take care of adding the very first solution.

            nToPrint = m_afMultiOptimizeBestValues.size();
            m_afMultiOptimizeBestParams + afParams;
            m_afMultiOptimizeBestValues + fValue;
            
        } else if (bAcceptDistFound) {
            // Case 2)  We found at least one solution which was within a refinement distance.  This means
            //          that we only update the worst solution having an acceptible distance.
            
            if (((m_afMultiOptimizeBestValues[nWorstPickAndAcceptDist] > fValue) && (m_bMultiOptimizeMinimize)) ||
                ((m_afMultiOptimizeBestValues[nWorstPickAndAcceptDist] < fValue) && (!m_bMultiOptimizeMinimize))) {
                memcpy(&m_afMultiOptimizeBestParams[nWorstPickAndAcceptDist*m_nMultiOptimizeNumParams + 0],&afParams[0],sizeof(double)*m_nMultiOptimizeNumParams);
                m_afMultiOptimizeBestValues[nWorstPickAndAcceptDist] = fValue;
                nToPrint = nWorstPickAndAcceptDist;
            } else           
                nToPrint = -1;
        } else {
            // Case 3)  No solutions found within the refinement distance.  So, update
            //          the absolute worst solution.
            if (((m_afMultiOptimizeBestValues[nWorstPick] > fValue) && (m_bMultiOptimizeMinimize)) ||
                ((m_afMultiOptimizeBestValues[nWorstPick] < fValue) && (!m_bMultiOptimizeMinimize))){
                memcpy(&m_afMultiOptimizeBestParams[nWorstPick*m_nMultiOptimizeNumParams + 0],&afParams[0],sizeof(double)*m_nMultiOptimizeNumParams);
                m_afMultiOptimizeBestValues[nWorstPick] = fValue;
                nToPrint = nWorstPick;
            } else
                nToPrint = -1;
        } 


        if ((m_bMultiOptimizePrint) && (nToPrint>=0)) {
            printf("Trial [%d] Slot [%d] %8.2lf [",m_nMultiOptimizeTrialCount,nToPrint,m_afMultiOptimizeBestValues[nToPrint]);
            for (nx =0; nx < afParams.size(); nx++) 
                printf("%8.2lf ",afParams[nx]);
            printf("]\n");
        };
    };

    m_nMultiOptimizeCount++;

    if ((m_nMultiOptimizeCount) && (!(m_nMultiOptimizeCount % (m_afMultiOptimizeLastBestParams.size()/m_nMultiOptimizeNumParams*m_nMultiOptimizeBatchSize)))) {

        -m_afMultiOptimizeLastBestParams + m_afMultiOptimizeBestParams;
                
        m_nMultiOptimizeTrialCount++;               
        // Are we at our last iteration?
        if (m_nMultiOptimizeTrialCount == m_nMultiOptimizeNumTrials) {
            // Yes.  We are done, load the best solution into afParams and exit
            for (nx= 1, ny = 0; nx < m_afMultiOptimizeBestValues.size(); nx++) {
            if (((m_afMultiOptimizeBestValues[ny] > m_afMultiOptimizeBestValues[nx]) && (m_bMultiOptimizeMinimize)) ||
                ((m_afMultiOptimizeBestValues[ny] < m_afMultiOptimizeBestValues[nx]) && (!m_bMultiOptimizeMinimize)))
                ny = nx;
            }
            memcpy(&afParams[0],&m_afMultiOptimizeBestParams[ny*m_nMultiOptimizeNumParams + 0],sizeof(double)*m_nMultiOptimizeNumParams);
            fValue = m_afMultiOptimizeBestValues[ny];
            return 1;
        };
        m_nMultiOptimizeCount = 0;
    };
    // Copy over the next solution. solution

    memcpy(&afParams[0],&m_afMultiOptimizeLastBestParams[(m_nMultiOptimizeCount/m_nMultiOptimizeBatchSize)*m_nMultiOptimizeNumParams],sizeof(double)*m_nMultiOptimizeNumParams);
    memcpy(&m_afMultiOptimizeLastParamsOut[0],&afParams[0],sizeof(double)*m_nMultiOptimizeNumParams);

    return 0;
};



/*

      {
        double fValue;
        const int nNumTrials = 100;
        const int nBatchSize = 100;
        const int nBestToPickInBatch = 2;
        itr<double> afParams;
        double fStep;
        int nStat;
        const double fPoint00 = 0.0;
        const double fPoint01 = 0.0;
        const double fPoint10 = 6.0;
        const double fPoint11 = 8.0;
                
        m_bMultiOptimizePrint = true;
        nMultiOptimize(fValue,-afParams,true,nBatchSize,nBestToPickInBatch,nNumTrials,1);
        -afParams + 3.1 + 4.0;
        nMultiOptimize(fValue,afParams);
        
        while (!(nStat = nMultiOptimize(fValue,afParams))) {
            fStep = 1.0/(m_nMultiOptimizeTrialCount + 1);
            m_fMultiOptimizeMaxDist = fStep*2.1;
            
            afParams[0] += ((rand() % 2)*2 - 1.0)*fStep;
            afParams[1] += ((rand() % 2)*2 - 1.0)*fStep;
            fValue = min(
                sqrt((afParams[0] - fPoint00)*(afParams[0] - fPoint00) + (afParams[1] - fPoint01)*(afParams[1] - fPoint01)),
                sqrt((afParams[0] - fPoint10)*(afParams[0] - fPoint10) + (afParams[1] - fPoint11)*(afParams[1] - fPoint11)));
            
        };
        printf("done.\n");
    };

*/


itr<int>        g_anSwapSizes;
itr<void*>      g_apvSwapPointers;

int qsort_swap_arrays(const void  *elem1, const void *elem2,size_t width) {
    unsigned char* pcPoint;
    unsigned char cChar;
    int nSize;
    int nIndex1,nIndex2;
    int nx,ny;

    nIndex1 = ((unsigned char*) elem1 - (unsigned char*) g_apvSwapPointers[0])/g_anSwapSizes[0];
    nIndex2 = ((unsigned char*) elem2 - (unsigned char*) g_apvSwapPointers[0])/g_anSwapSizes[0];
    for (nx = 0; nx < g_anSwapSizes.size(); nx++) {
        pcPoint = (unsigned char*) g_apvSwapPointers[nx];
        nSize = g_anSwapSizes[nx];
        for (ny = 0; ny < nSize; ny++) {
            cChar = pcPoint[nIndex1*nSize + ny];
            pcPoint[nIndex1*nSize + ny] = pcPoint[nIndex2*nSize + ny];
            pcPoint[nIndex2*nSize + ny] = cChar;
        };
    };
    return 0;
};


void qsortswap(void* base,size_t num,size_t width,int (*compare )(const void *elem1, const void *elem2 ),int(*myswap)(const void  *elem1, const void *elem2,size_t width)) 
{
    char *lo, *hi;              /* ends of sub-array currently sorting */
    char *mid;                  /* points to middle of subarray */
    char *loguy, *higuy;        /* traveling pointers for partition step */
    unsigned size;              /* size of the sub-array */
    char *lostk[30], *histk[30];
    int stkptr;                 /* stack for saving sub-array to be processed */

    /* Note: the number of stack entries required is no more than
       1 + log2(size), so 30 is sufficient for any array */

    if (num < 2 || width == 0)
        return;                 /* nothing to do */

    stkptr = 0;                 /* initialize stack */

    lo = (char *)base;
    hi = (char *)base + width * (num-1);        /* initialize limits */

    /* this entry point is for pseudo-recursion calling: setting
       lo and hi and jumping to here is like recursion, but stkptr is
       prserved, locals aren't, so we preserve stuff on the stack */
recurse:

    size = (hi - lo) / width + 1;        /* number of el's to sort */

                                         /* First we pick a partititioning element.  The efficiency of the
                                         algorithm demands that we find one that is approximately the
                                         median of the values, but also that we select one fast.  Using
                                         the first one produces bad performace if the array is already
                                         sorted, so we use the middle one, which would require a very
                                         wierdly arranged array for worst case performance.  Testing shows
                                         that a median-of-three algorithm does not, in general, increase
    performance. */
    
    mid = lo + (size / 2) * width;      /* find middle element */
    myswap(mid, lo, width);               /* swap it to beginning of array */
    
                                          /* We now wish to partition the array into three pieces, one
                                          consisiting of elements <= partition element, one of elements
                                          equal to the parition element, and one of element >= to it.  This
                                          is done below; comments indicate conditions established at every
    step. */
    
    loguy = lo;
    higuy = hi + width;
    
    /* Note that higuy decreases and loguy increases on every iteration,
    so loop must terminate. */
    for (;;) {
    /* lo <= loguy < hi, lo < higuy <= hi + 1,
    A[i] <= A[lo] for lo <= i <= loguy,
        A[i] >= A[lo] for higuy <= i <= hi */
        
        do  {
            loguy += width;
        } while (loguy <= hi && compare(loguy, lo) <= 0);
        
        /* lo < loguy <= hi+1, A[i] <= A[lo] for lo <= i < loguy,
        either loguy > hi or A[loguy] > A[lo] */
        
        do  {
            higuy -= width;
        } while (higuy > lo && compare(higuy, lo) >= 0);
        
        /* lo-1 <= higuy <= hi, A[i] >= A[lo] for higuy < i <= hi,
        either higuy <= lo or A[higuy] < A[lo] */
        
        if (higuy < loguy)
            break;
        
            /* if loguy > hi or higuy <= lo, then we would have exited, so
            A[loguy] > A[lo], A[higuy] < A[lo],
        loguy < hi, highy > lo */
        
        myswap(loguy, higuy, width);
        
        /* A[loguy] < A[lo], A[higuy] > A[lo]; so condition at top
        of loop is re-established */
    }
    
    /*     A[i] >= A[lo] for higuy < i <= hi,
    A[i] <= A[lo] for lo <= i < loguy,
    higuy < loguy, lo <= higuy <= hi
    implying:
    A[i] >= A[lo] for loguy <= i <= hi,
    A[i] <= A[lo] for lo <= i <= higuy,
    A[i] = A[lo] for higuy < i < loguy */
    
    myswap(lo, higuy, width);     /* put partition element in place */
    
                                  /* OK, now we have the following:
                                  A[i] >= A[higuy] for loguy <= i <= hi,
                                  A[i] <= A[higuy] for lo <= i < higuy
    A[i] = A[lo] for higuy <= i < loguy    */
    
    /* We've finished the partition, now we want to sort the subarrays
    [lo, higuy-1] and [loguy, hi].
    We do the smaller one first to minimize stack usage.
    We only sort arrays of length 2 or more.*/
    
    if ( higuy - 1 - lo >= hi - loguy ) {
        if (lo + width < higuy) {
            lostk[stkptr] = lo;
            histk[stkptr] = higuy - width;
            ++stkptr;
        }                           /* save big recursion for later */
        
        if (loguy < hi) {
            lo = loguy;
            goto recurse;           /* do small recursion */
        }
    }
    else {
        if (loguy < hi) {
            lostk[stkptr] = loguy;
            histk[stkptr] = hi;
            ++stkptr;               /* save big recursion for later */
        }
        
        if (lo + width < higuy) {
            hi = higuy - width;
            goto recurse;           /* do small recursion */
        }
    }
    
    
    /* We have sorted the array, except for any pending sorts on the stack.
    Check if there are any, and do them. */
    
    --stkptr;
    if (stkptr >= 0) {
        lo = lostk[stkptr];
        hi = histk[stkptr];
        goto recurse;           /* pop subarray from stack */
    }
    else
        return;                 /* all subarrays done */
}

///////////////////////////////////////////////////////////////////////////
// Round off an input double value keeping an input number of decimal digits
double	dRoundOff(double dIn, int nDecDigits, DTREK_WORD wCtrl)
{
	if( 0.0 == dIn )
        return 0.0;

    bool    bPlus = dIn > 0.0 ? true : false;

    dIn = fabs(dIn);
    
    double		dTenPower = pow(10.0, nDecDigits);
	
	dIn *= dTenPower;
    
    switch(wCtrl)
    {
    case DTREK_VEC_ROUNDOFF_CEILING:
        dIn = ceil(dIn);
        break;
    case DTREK_VEC_ROUNDOFF_FLOOR:
        dIn = floor(dIn);
        break;
    default:
        dIn = ceil(dIn) - dIn <= 0.5 ? ceil(dIn) : floor(dIn);
        break;
    }
	
    dIn /= dTenPower;

    if( !bPlus )
        dIn = -dIn;

	return dIn;
}
////////////////////////////////////////////////////////////////////////////
// Get a nearest multiple integer for an input double value.
// It is assumeed the integer is positive.
double dGetNearestMultipleInteger(double dIn, int nMultInteger, DTREK_WORD wCtrl)
{
    if( nMultInteger <= 0 )
        return dIn; // we require that the integer is positive
    
    dIn /= nMultInteger;

    dIn = dRoundOff(dIn, 0, wCtrl);

    dIn *= nMultInteger;

    return dIn;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Given the input value fIn, find a check integer from vector anCheckValues, closest to the input value.
int nGetNearestIntFromVector(double fIn, std::vector<int>& anCheckValues, DTREK_WORD wCtrl)
{
    if( 0 == (int)anCheckValues.size() )
        return (int)dRoundOff(fIn, 0, wCtrl); // simply return the closest integer
    
    double  fDist_closest = 0.0; 
    double  fDist_ii = 0.0;

    int     ii_Closest = 0;

    for(int ii=0; ii < anCheckValues.size(); ii++)
    {
        fDist_ii = fabs(fIn - anCheckValues[ii]);

        if( 0 == ii || fDist_ii < fDist_closest )
        {
            fDist_closest = fDist_ii;
            ii_Closest = ii;
        }
        else if( fDist_ii == fDist_closest && 0U != wCtrl ) // if there is a preference for floor vs ceiling - make a choice accordingly
        {
            if( wCtrl == DTREK_VEC_ROUNDOFF_CEILING ) // choose the one that is larger
            {
                if( anCheckValues[ii_Closest] < anCheckValues[ii] )
                {
                    fDist_closest = fDist_ii;
                    ii_Closest = ii;
                }
            }
            else if(  wCtrl == DTREK_VEC_ROUNDOFF_FLOOR )
            {
                if( anCheckValues[ii_Closest] > anCheckValues[ii] )
                {
                    fDist_closest = fDist_ii;
                    ii_Closest = ii;
                }
             }
        }
    }
    
    return anCheckValues[ii_Closest];
}
////////////////////////////////////////////////////////////////////////////////////
// A helper function to be used by command line parsers
bool bIsSkipCommandLineArgument(const char* sArg)
{
	if( 0==strcmp(sArg, "\\") || 0==strcmp(sArg, "") )
	{
		printf("\n--------- Skipping argument --------\n");
		return true;
	}
	
	return false;
}
///////////////////////////////////////////////////////////////////////////

