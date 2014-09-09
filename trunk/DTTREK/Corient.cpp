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
// Corient.cc            Initial author: J.W. Pflugrath           26-Jun-1995
//  This file contains the member functions of class Corient which implements
//    crystal orientation encapsulation of d*TREK.
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
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Corient.h"         // Class definition and prototypes

#include <ctype.h>

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
#endif

//+Definitions, constants, and initialization of static member variables

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Corient::Corient()
{
 (void) nInitValues();
}

Corient::Corient(Cimage_header& oHeader)
{
  int     nStat;
  Cstring sTemp;

  (void) nInitValues();

  if (!oHeader.bIsAvailable())
    {
      cout << "Corient ERROR: image header is not valid."
           << "  Cannot construct Corient object!\n";
      nStat = 1;
    }
  else 
    {
      // Try to get required information from the header

      nStat = 0;

      m_poCrystal     = new Ccrystal (oHeader);
      m_bNewCrystal   = TRUE;
      m_poCrysGonio   = new Cgoniometer(oHeader, Ccrystal::ms_sCrystalPrefix);
      m_bNewCrysGonio = TRUE;

      m_poSource      = new Csource(oHeader);
      m_bNewSource    = TRUE;
    }
}

Corient::~Corient() 
{
  if (m_bNewSource && (NULL != m_poSource) ) 
    {
      delete m_poSource;
      m_poSource   = NULL;
      m_bNewSource = FALSE;
    }

  if (m_bNewCrystal && (NULL != m_poCrystal) )
    {
      delete m_poCrystal;
      m_poCrystal   = NULL;
      m_bNewCrystal = FALSE;
    }

  if (m_bNewCrysGonio && (NULL != m_poCrysGonio) ) 
    {
      delete m_poCrysGonio;
      m_poCrysGonio   = NULL;
      m_bNewCrysGonio = FALSE;
    }
}

int Corient::nInitValues(void)
{
  m_poSource         = NULL;
  m_poCrystal        = NULL;
  m_poCrysGonio      = NULL;

  m_bNewSource       = FALSE;
  m_bNewCrystal      = FALSE;
  m_bNewCrysGonio    = FALSE;

  m_nNumSoln         = 0;
  return (0);
}

int Corient::nList() 
{
  int i, j;
  Cstring a3sVecNames[3];
  (void) m_poCrysGonio->nGetNames(3, &a3sVecNames[0]);
  printf("\nOrient listing\n"
         "When the crystal goniometer is at:\n"
         "           %9s %9s %9s\n"
         "             %9.3f %9.3f %9.3f"
	 "\nthe angles between crystal vectors and lab vectors are:\n",
	 a3sVecNames[0].string(), a3sVecNames[1].string(),
	 a3sVecNames[2].string(),
	 m_fGonAnglesUsed[0], m_fGonAnglesUsed[1], m_fGonAnglesUsed[2]);

  //                  123456789 123456789 123456789 123456789 123456789 123456789 
  printf("                                Crystal vectors\n");
  printf("                  a*        b*        c*        a         b         c\n");
  printf("Lab vectors:\n");
  for (i = 0; i < 7; i++)
    {
      printf("%9s    ", m_sVecNames[i].string());
      for (j = 0; j < 6; j++)
	{
	  printf("%9.3f ", m_fVecAngles[i][j]);
	}
      printf("\n");
    }
  return (0);
}

int Corient::nSetup(void)
{
  int   nStat;
  double fDet;

  // Calculate all the matrices and vectors used by Corient

  // Crystal stuff:

  nStat = m_poCrystal->nCalcRotMatrix();
  m_poCrystal->vGetRotMatrix(&m_fCrysRotMatrix[0][0]);
  if (0 != nStat) return (nStat);

  nStat = m_poCrystal->nCalcBMatrix();
  if (0 != nStat) return (nStat);
  m_poCrystal->vGetBMatrix(&m_fCrysBMatrix[0][0]);
  int nx;
  double *pfTemp;
  pfTemp = &m_fCrysBMatrix[0][0];
  for (nx = 0; nx < 9; nx++, pfTemp++)
    if (ABS(*pfTemp) < 0.000001) *pfTemp = 0.0;

  // Invert crystal B matrix

  fDet = fInvMat3D(&m_fCrysBMatrix[0][0], &m_fInvCrysBMatrix[0][0]);
  if (0.0 == fDet) return (4);
  pfTemp = &m_fInvCrysBMatrix[0][0];
  for (nx = 0; nx < 9; nx++, pfTemp++)
    if (ABS(*pfTemp) < 0.000001) *pfTemp = 0.0;

  // Transpose InvCrysBMatrix

  vCopyMat3D((double *)m_fInvCrysBMatrix, (double *)m_fTransInvCrysBMatrix);
  vTranMat3D(m_fTransInvCrysBMatrix);

  // Get crystal goniometer rotation matrix at datum
  
  m_poCrysGonio->vCalcGetRotMatrix(&m_fGonMatrix[0][0],-1);

  // Multiply goniometer rotation matrix by crystal rotation matrix

  vMulMat3DMat3D(m_fGonMatrix, m_fCrysRotMatrix, m_fGonCrysRotMatrix);

  return (0);
}

int Corient::nCalcGetVecAngles(const double *pfVec, double *pfAngles)

{
  // Calculate the angles between a 3D vector pfVec in real space and
  // a*, b*, c*, a, b, c when crystal is at goniometer datum

  int   i, j;
  int   nStat;
  double fVecIn[3];
  double fVecABC[3];
  double fVecABCRotd[3];
  double fMatTemp1[2][3][3];
  double fMatTemp2[3][3];
  double fTemp;
  double fTemp1;

  vCopyVec3D(pfVec, fVecIn);

  // Set resultant angles to bogus values in case of premature return
  for (i = 0; i < 6; i++) 
    {
      pfAngles[i] = -999.0;
    }

  fTemp = fNormVec3D(fVecIn);
  if (0.0 == fTemp) return (4);

  nStat = nSetup();  // This may be redundant but do it to make sure things
                     // are up to date anyways

  vMulMat3DMat3D(m_fGonCrysRotMatrix, m_fCrysBMatrix , fMatTemp2);
  vCopyMat3D((double *)fMatTemp2, &fMatTemp1[0][0][0]); // j=0 below
  vMulMat3DMat3D(m_fGonCrysRotMatrix, m_fTransInvCrysBMatrix , fMatTemp2);
  vCopyMat3D((double *)fMatTemp2, &fMatTemp1[1][0][0]); // j=1 below
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  fVecABC[j] = 0.0;
	  if (i == j)
	    fVecABC[j] = 1.0;
	}

      // fVecABC = 1,0,0 or 0,1,0 or 0,0,1
      //                              -1T
      // Rotate fVecABC by GCB  or GC(B )   (GCB is fMatTemp2)
      //                   ===     == =      ===    =========
      // in order to get direction of a, a* in lab frame

      for (j = 0; j < 2; j++)   // j = 0: real space; j = 1: reciprocal space
	{
	  vCopyMat3D(&fMatTemp1[j][0][0], (double *) fMatTemp2);
	  vMulMat3DVec3D(fMatTemp2, fVecABC, fVecABCRotd);
	  fTemp = fNormVec3D(fVecABCRotd);
	  fTemp = fDot3D(fVecIn, fVecABCRotd);  // Angle = acos (VecIn . VecRotated)
	  if (-1.0 > fTemp) fTemp = -1.0;
	  if ( 1.0 < fTemp) fTemp =  1.0;
	  //	  if ( (fTemp > .99) && (fTemp < 1.0) )
	  //  cout << "fTemp close to, but not 1.0\n";
	  fTemp1 = acos(fTemp) / Gs_dRADIANS_PER_DEGREE;
	  pfAngles[i + 3*j] = fTemp1;
	  //if ( (fTemp1 > 0.01) && (fTemp1 < 90.0) )
	  //  cout << "funny fTemp value: " << fTemp1 << endl;
	}
    }
  return (0);
}

int Corient::nCalcAngles(const double *pfGonAnglesIn)
{
  // Calculate the angles that some vectors make with the principal crystal
  // axes a, b, c, a*, b*, c* when goniometer is at angles *pfGonAngles
  // If pfGonAnglesIn == NULL, then use the goniometer datum values

  int   i, j;
  int   nStat;
  double fMat[3][3];
  double fSave[3];
  double fVecTemp[3];

  if (NULL == pfGonAnglesIn)
    {
      nStat = m_poCrysGonio->nGetDatum(3, (double *) m_fGonAnglesUsed);
    }
  else 
    {
      for (i = 0; i < 3; i++)
	m_fGonAnglesUsed[i] = pfGonAnglesIn[i];
    }

  m_poSource->vCalcGetS0(&m_fVecs[0][0]);
  m_sVecNames[0] = "Source";

  nStat = m_poCrysGonio->nGetVectors(3, &m_fVecs[1][0]);
  if (0 < nStat)
    {
      // This crystal goniometer does not have 3 axes
    }

  nStat = m_poCrysGonio->nGetNames(3, &m_sVecNames[1]);
  nStat = m_poCrysGonio->nGetDatum(3, (double *)fSave);
  nStat = m_poCrysGonio->nSetDatum(3, m_fGonAnglesUsed);
  if (0 < nStat)
    {
      // This crystal goniometer does not have 3 axes
    }

  // Compute laboratory direction of 2nd and 3rd goniometer axes after 
  // rotating goniometer to datum

  m_poCrysGonio->vCalcGetRotMatrix((double *)fMat, -1,
				 m_fGonAnglesUsed[0],
				 m_fGonAnglesUsed[1], 0.0);
  vMulMat3DVec3D(fMat, &m_fVecs[3][0], fVecTemp);
  vCopyVec3D(fVecTemp, &m_fVecs[3][0]);

  m_poCrysGonio->vCalcGetRotMatrix((double *)fMat, -1,
				 m_fGonAnglesUsed[0],
				 0.0, 0.0);
  vMulMat3DVec3D(fMat, &m_fVecs[2][0], fVecTemp);
  vCopyVec3D(fVecTemp, &m_fVecs[2][0]);

  // Create lab vectors along X, Y, Z

  for (i = 0; i < 3; i++) 
    {
      for (j = 0; j < 3; j++)
	{
	  m_fVecs[4+i][j] = 0.0;
	  if (i == j)
	    m_fVecs[4+i][j] = 1.0;
	}    
    }

  m_sVecNames[4] = "Lab X";
  m_sVecNames[5] = "Lab Y";
  m_sVecNames[6] = "Lab Z";

  for (i = 0; i < 7; i++)
    {
      (void) fNormVec3D(&m_fVecs[i][0]);
      nCalcGetVecAngles(&m_fVecs[i][0], &m_fVecAngles[i][0]);
    }  
  nStat = m_poCrysGonio->nSetDatum(3, (double *)fSave);
  return (0);
}

int Corient::nSetOrient(const Cstring &sOrient) 
{
  // Try to set the orientation according to the input string.
  //  Crystal_vector         Lab_vector
  //  [-]a|b|c[*] ||  [-]x|y|z|source|gon1
  //              == 
  //  rot angle vector
  int  i, j;
  int  nStat;

  Cstring sTemp;
  Cstring sInput;

  nStat = nSetup(); // Get all the matrices we'll need

  // Parse the input string

  m_sRequest = sOrient;
  sInput     = sOrient;
  sInput.downcase();

  m_nMode = 0;
  if (sInput.contains("rot"))
    {
      m_nMode = 1;
      nRotate(sInput);
    }
  else
    {
      m_nMode = 0;
      for (j = 0; j < 2; j++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      m_fCrysVec[j][i] = 0.0;
	      m_fLabVec[j][i]  = 0.0;
	    }
	}
    
      // Search for the operator "||"
    
      sTemp = sInput.before("||");
      i = 0;
      while ( ("" != sTemp) && (2 > i) )
	{
	  nStat = nParseVec(sTemp, &m_fCrysVec[i][0]);
      
	  sTemp = sInput.after("||");
      
	  if ("" != sTemp)
	    {
	      if (sTemp.contains(';')) sTemp = sTemp.before(';');
	      nStat = nParseVec(sTemp, &m_fLabVec[i][0]);
	    }

	  // Set sTemp to unparsed part of sInput

	  sInput = sInput.after(';');
	  sTemp  = sInput.before("||");
	  i++;
	}
    
      nStat = 0;
      for (i = 0; i < 2; i++)
	{
	  if (0.0 == fLenVec3D(&m_fCrysVec[i][0]))
	    {
	      nStat++;
	      cout << "Invalid ";
	      cout << "crystal vector " << i+1 << ": "
		   << m_fCrysVec[i][0] << ", "
                   << m_fCrysVec[i][1] << ", "
                   << m_fCrysVec[i][2] << endl;
	    }

	  if (0.0 == fLenVec3D(&m_fLabVec[i][0]))
	    {
	      nStat++;
	      cout << "Invalid ";
	      cout << "lab  vector " << i+1 << ": "
                   << m_fLabVec[i][0] << ", "
                   << m_fLabVec[i][1] << ", "
                   << m_fLabVec[i][2] << endl;
	    }
	}
    }
  return (nStat);
}

int Corient::nCalcOrient(void)
{
  // Calculate new goniometer values for when member variable
  // m_fCrysVec[i] is aligned with m_fLabVec[i], for i=0,1
  // This work follows the work of Phil Evans, MRC found in MADNES
  // The major difference is that the crystal orientation angles in d*TREK
  // are always referenced to a goniostat angles of 0.,0.,0.. This simplifies
  // alot of the code, since if you move the goniometer, the crystal orientation
  // angles do not change.

  int   nStat;
  int   i, j, k;
  double fDet;
  double fMatC[3][3];  // Base vectors in crystal frame for required orientation
  double fInvMatC[3][3];  // Inverse of fMatC
  double fMatL[3][3];  // Base vectors in lab     frame for required orientation
  double fSignPerm[4][3] = { { 1.0,  1.0,  1.0},  // Sign permutations
			    {-1.0, -1.0,  1.0},
			    { 1.0, -1.0, -1.0},
			    {-1.0,  1.0, -1.0} };
  double fMatLSigned[3][3];
  double fNewGonMatrix[3][3];
  double fGonVecs[3][3];
  double fSignTemp;
  double *pfTemp;

  m_nNumSoln = 0;
  nStat = nSetup(); // Calculate most all the matrices we'll need
  fDet = fInvMat3D((double*)m_fCrysRotMatrix, (double*)m_fInvCrysRotMatrix);
  if (0.0 == fDet) return (-1);

  nStat = m_poCrysGonio->nGetVectors(3, &fGonVecs[0][0]);
  if (0 < nStat) 
    {
      // This crystal goniometer does not have 3 axes so we are in trouble
    }

  if (0 == m_nMode)
    {
      // Construct fMatC
      // First vector || crystal vector 1

      vCopyVec3D(&m_fCrysVec[0][0], &fMatC[0][0]);
      fDet = fNormVec3D(&fMatC[0][0]);
      if (0.0 == fDet) return (-1);
    
      // Third vector _|_ to crystal vectors 1 and 2

      vCross3D(&m_fCrysVec[0][0], &m_fCrysVec[1][0], &fMatC[2][0]);
      fDet = fNormVec3D(&fMatC[2][0]);
      if (0.0 == fDet) return (-1);

      // Second vector _|_ to first and third vectors 2 = 3 X 1

      vCross3D(&fMatC[2][0], &fMatC[0][0], &fMatC[1][0]);
      fDet = fNormVec3D(&fMatC[1][0]);
      if (0.0 == fDet) return (-1);

      // Construct fMatL
      // First vector || lab vector 1

      vCopyVec3D(&m_fLabVec[0][0], &fMatL[0][0]);
      fDet = fNormVec3D(&fMatL[0][0]);
      if (0.0 == fDet) return (-1);

      // Third vector _|_ to lab vectors 1 and 2

      vCross3D(&m_fLabVec[0][0], &m_fLabVec[1][0], &fMatL[2][0]);
      fDet = fNormVec3D(&fMatL[2][0]);
      if (0.0 == fDet) return (-1);

      // Second vector _|_ to first and third vectors

      vCross3D(&fMatL[2][0], &fMatL[0][0], &fMatL[1][0]);
      fDet = fNormVec3D(&fMatL[1][0]);
      if (0.0 == fDet) return (-1);

      // Now we have the two frames in fMatC and fMatL.              -1
      // The matrix which superimposes the two frams is fMatL * fMatC  .
      // This is then the orientation matrix GC (m_fGonRotMatrix * m_fCrysRotMatrix).
      //                                     ==
      // The new goniometer rotation matrix is then:
      //                     -1
      //   G' = fMatL * fMatC   * m_fInvCrysRotMatrix,
      //   =
      //   But G' can have several permutations all rotated by 180 degrees.
  
      fDet = fInvMat3D((double*)fMatC, (double*)fInvMatC);
      if (0.0 == fDet) return (-1);
      vMulMat3DMat3D(fInvMatC, m_fInvCrysRotMatrix, fMatC);  // fMatC is temp var here

      // Loop over 4 sign permutations of fMatL

      pfTemp = &m_fAngleSoln[0][0];  // Stupid pointer tricks

      for (k = 0; k < 4; k++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      for (j = 0; j < 3; j++)
		{
		  fMatLSigned[i][j] = fSignPerm[k][i] * fMatL[i][j];
		}
	    }
	  vMulMat3DMat3D(fMatLSigned, fMatC, fNewGonMatrix);
	  fSignTemp = 1.0;
	  for (j = 0; j < 2; j++)
	    {
	      nStat = nDeconvMat3DVec3DRot3(fNewGonMatrix,
					    &fGonVecs[2][0],
					    &fGonVecs[1][0],
					    &fGonVecs[0][0],
					    pfTemp+2, pfTemp+1, 
					    pfTemp, fSignTemp);

	      if (0 == nStat)
		{
		  pfTemp = pfTemp + 3;
		  m_nNumSoln++;
		}
	      fSignTemp = -fSignTemp; // for next time round
	    }
	}
    }
  else if (1 == m_nMode)
    {
      // Rotate mode was used by nSetOrient

      // Construct
      vConvRotVec3DMat3D(m_fLabVec[0][0], &m_fCrysVec[0][0], (double *)fMatC);
      vMulMat3DMat3D(fMatC, m_fGonMatrix, fNewGonMatrix);
      fSignTemp = 1.0;
      pfTemp = &m_fAngleSoln[0][0];  // Stupid pointer tricks
      for (j = 0; j < 2; j++)
	{
	  nStat = nDeconvMat3DVec3DRot3(fNewGonMatrix,
					&fGonVecs[2][0],
					&fGonVecs[1][0],
					&fGonVecs[0][0],
					pfTemp+2, pfTemp+1, pfTemp, fSignTemp);

	  if (0 == nStat)
	    {
	      pfTemp = pfTemp + 3;
	      m_nNumSoln++;
	    }
	  fSignTemp = -fSignTemp; // for next time round
	}
    }
  return (0);
}

int Corient::nListSoln(void)
{
  if (0 >= m_nNumSoln)
    {
      cout << "No solutions available." << endl;
      cout << "Usage:\n"
           << "-EITHER-\n"
	   << "   crystal_vector1 || lab_vector1; crystal_vector2 || lab_vector2\n"
	   << " where\n"
	   << " crystal_vector is one of\n"
	   << "   a, b, c, -a, -b, -c, a*, b*, c*, -a*, -b*, -c*, (h k l), [r s t]\n"
	   << "   (h k l) is a Miller index, while [r s t] is a real direction\n"
	   << "                              along real unit cell axes;\n"
	   << " lab_vector is one of\n"
           << "   x, y, z, -x, -y, -z, source, -source, gon1, -gon1, xl yl zl\n"
	   << "   gon1 is the name of the first goniometer axis (i.e. omega)\n"
           << "   xl yl zl is a sequence of 3 numbers.\n"
           << "-OR-\n"
	   << "   rot fAngle lab_vector\n"
	   << " where\n"
	   << " rot is the word rot\n"
           << " fAngle is a relative amount to rotate in degrees\n"
           << " lab_vector is described above and defaults to source if missing.\n" << endl;
      return (-1);
    }
  else
    {
      int     i;
      double   a3fGonAngles[3];
      Cstring a3sVecNames[3];
      printf("\nSolutions to orient request (>> %s <<):\n", m_sRequest.string());
      (void) m_poCrysGonio->nGetNames(3, &a3sVecNames[0]);
      printf("   Number  %9s %9s %9s\n",
	     a3sVecNames[0].string(), a3sVecNames[1].string(),
	     a3sVecNames[2].string());

      (void) m_poCrysGonio->nGetDatum(3, (double *) a3fGonAngles);
      printf("%9d    %9.3f %9.3f %9.3f <- Current goniometer relative zero\n",
	     (int) 0, a3fGonAngles[0], a3fGonAngles[1], a3fGonAngles[2]);

      for (i = 0; i < m_nNumSoln; i++)
	{
	  printf("%9d    %9.3f %9.3f %9.3f\n",
		 i+1, m_fAngleSoln[i][0], m_fAngleSoln[i][1], m_fAngleSoln[i][2]);
	}
    }
  return (0);
}

int 
Corient::nPickSoln(const int nSolnIn)
{
  // List solutions stored in member variables then prompt
  // for a solution to select.  Update the m_poCrysGonio object
  // with the valid selected angle.
  // Return the solution number - 1.
  // That is a number between 0 and m_nNumSoln-1.  If a negative
  // number is returned, then there was either an error or
  // no solution was chosen.
  // 
  int nStat;
  int nSoln;

  nSoln = nSolnIn;
  nStat = nListSoln();
  if (0 == nStat)
    {
      if ( (0 >= nSoln) || (nSoln > m_nNumSoln) )
	{
	  // If solution is out of range, prompt for the soln number

	  cout << "Enter the number of the selected solution: ";
	  cin >> nSoln;
	}
      if ( (0 < nSoln) && (nSoln <= m_nNumSoln) )
	{
	  nStat = m_poCrysGonio->nSetDatum(3, &m_fAngleSoln[nSoln-1][0]);
	  if (0 == nStat)
	    {
	      nStat = nCalcAngles(&m_fAngleSoln[nSoln-1][0]);
	    }
	}
      else if (0 != nStat)
	nStat = -1;
    }
  return (nStat);
}

int Corient::nUpdateHeader(Cimage_header* poHeader)
{
  // Warning, SCAN_CRYS_RELZERO is not changed!

  return (m_poCrysGonio->nUpdateHeader(poHeader, Ccrystal::ms_sCrystalPrefix));
}

int Corient::nRotate(const Cstring& sRotate)
{
  // Input string is of form "rot fAngle vec"
  // where
  //     rot  is the character string rot
  //     fAngle is a floating point number
  //     is a vector parsable by nParseVec()

  Cstring sInput;
  Cstring sWords[2];
  int    nStat;
  
  sInput = sRotate;
  sInput.downcase();
  if (sInput.contains("rot"))
    {
      sInput = sInput.after("rot");

      // Remove any whitespace at beginning of sInput

      while(isspace((int)sInput.GetAt(0)))
	sInput = sInput.after(0);

      nStat  = split(sInput, sWords, 2, " \n\t");
      if (1 >= nStat)
	{
	  // No vector was specified, default to source

	  sInput = sWords[0] + "source";
	}
      nStat = sscanf(sWords[0].string(), "%f", &m_fLabVec[0][0]);
      if (1 == nStat)
	{
	  sInput = sInput.after(sWords[0]);
	  nStat = nParseVec(sInput, &m_fCrysVec[0][0]);
	}
      else
	nStat = -1;
    }
  return (nStat);
}

int Corient::nParseVec(const Cstring& sVec, double *pfVec)
{
  // Parse a string that contains a vector description and
  // return the real space vector in pfVec.
  // The input Cstring sVec can contain one of the following:
  // Crystal principal vectors:
  //  a   -a   a*  -a*   b  -b   b*   -b*   c   -c   c*   -c  
  // Lab vectors:
  //  x   -x   y   -y    z   -z  source -source  gon1  -gon1
  //   (where gon1 is the name of the 1st crystal goniometer axis)
  // Crystal reciprocal vector (i.e. an hkl)
  // ( 1.0 2.0 3.0 )         
  // Crystal real vector (i.e. along real cell axes)
  // [ 1.0 2.0 3.0 ]
  // Lab vector
  // 1.0 2.0 3.0
  // The string is converted to lower case before parse, so things are
  // case insensitive.
  //

  Cstring sInput;
  int    nType;   // 0 lab basis; 1 real crystal basis; 2 recip crystal basis
  int    nStat;
  double  fVecTemp[3];
  Cstring sGonFirst;
  Cstring sWords[5];
  double  fDet;

  sInput    = sVec;

  // Remove any whitespace at beginning of sInput

  while(isspace((int)sInput.GetAt(0)))
    sInput = sInput.after(0);

  sInput.downcase();
  sGonFirst = "omega";
  nStat     = m_poCrysGonio->nGetNames(1, &sGonFirst);
  sGonFirst.downcase();

  pfVec[0] = 0.0;
  pfVec[1] = 0.0;
  pfVec[2] = 0.0;

  nType = 0;
  if (   ("a" == sInput)
      || ("a*" == sInput)
      || ("-a" == sInput)
      || ("-a*" == sInput) )
    {         // MUST check -n before checking n!
      pfVec[0] = 1.0;
      nType    = 1;
    }
  else if (   ("b" == sInput)
	   || ("b*" == sInput)
	   || ("-b" == sInput)
	   || ("-b*" == sInput) )
    {
      pfVec[1] = 1.0;
      nType    = 1;
    }
  else if (   ("c" == sInput)
	   || ("c*" == sInput)
	   || ("-c" == sInput)
	   || ("-c*" == sInput) )
    {
      pfVec[2] = 1.0;
      nType    = 1;
    }
  else if (   ("x" == sInput)
	   || ("-x" == sInput) )
    pfVec[0] = 1.0;
  else if (   ("y" == sInput)
	   || ("-y" == sInput) )
    pfVec[1] = 1.0;
  else if (   ("z" == sInput)
	   || ("-z" == sInput) )
    pfVec[2] = 1.0;
  else if (sInput.contains("source"))
    m_poSource->vCalcGetS0(pfVec);
  else if (sInput.contains(sGonFirst) || sInput.contains("gon1"))
    nStat = m_poCrysGonio->nGetVectors(1, pfVec);
  else
    {
      // No letter found, so what was found?
      if (sInput.contains('('))
	{
	  nType  = 2;
	  sInput = sInput.after('(');
	  if (sInput.contains(')')) sInput = sInput.before(')');
	}
      else if (sInput.contains('['))
	{
	  nType  = 1;
	  sInput = sInput.after('[');
	  if (sInput.contains(']')) sInput = sInput.before(']');
	}

      // Remove any whitespace at beginning of sInput

      while(isspace((int)sInput.GetAt(0)))
	sInput = sInput.after(0);

      nStat = split(sInput, sWords, 3, " ,\n\t");
      if (3 != nStat) 
	nStat = 4;
      else if (3 == nStat)
	{
	  nStat = 1;
	  int i;
	  for (i = 0; ( (i < 3) && (1 == nStat) ); i++) 
	    nStat = sscanf(sWords[i].string(), "%lf", &pfVec[i]);
	  if (1 == nStat) nStat = 0;  // Reset nStat to success
	}
      sInput = "";  // Do this so checks for * and - yield nothing.
    }
  if (sInput.contains('-'))
    vMulVec3DScalar(pfVec, -1.0, pfVec);
  if (sInput.contains('*') && (1 == nType)) nType = 2;

  if (0 != nType)
    {
      // Orthogonalize crystal vector

      vCopyVec3D(pfVec, fVecTemp);
      if (2 == nType)
	{
	  vMulMat3DVec3D(m_fCrysBMatrix, fVecTemp, pfVec);
	}
      else if (1 == nType)
	{
	  vMulMat3DVec3D(m_fTransInvCrysBMatrix, fVecTemp, pfVec);
	}
    }
  fDet = fNormVec3D(pfVec);  // Normalize it

  if (0.0 == fDet) return (-1);
  return (0);
}

void Corient::vGetAngleSoln(double fAngleSoln[8][3])
{
    for( int i=0; i<8; i++ ) {
        for( int j=0; j<3; j++ ) {
            fAngleSoln[i][j] = m_fAngleSoln[i][j];
        }
    }
}
