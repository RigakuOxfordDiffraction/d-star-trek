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
// Csource.cc            Initial author: J.W. Pflugrath           01-May-1995
//    This file contains the member functions of class Csource.
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

#include "Dtrek.h"
#include "Csource.h"            // Class definition and prototypes
                                //  Csource.h includes others
#include "dtrekdefs.h"
//+Definitions, constants, and initialization of static member variables

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

Cstring Csource::ms_sSourceVectors            = D_K_SourceVectors;
Cstring Csource::ms_sSourceValues             = D_K_SourceValues;
Cstring Csource::ms_sSourceIntensity          = D_K_SourceIntensity;
Cstring Csource::ms_sSourceSize               = D_K_SourceSize;
Cstring Csource::ms_sSourcePolarz             = D_K_SourcePolarz;
Cstring Csource::ms_sSourceCrossfire          = D_K_SourceCrossfire;
Cstring Csource::ms_sSourceKey                = D_K_SourceKey;
Cstring Csource::ms_sSourceSpectralDispersion = D_K_SourceSpectralDispersion;
Cstring Csource::ms_sSourceRefineFlags        = D_K_SourceRefineFlags;
Cstring Csource::ms_sSourcePrefix             = D_K_SourcePrefix;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Csource::Csource()
{
  (void) nInitValues();
}

Csource::Csource(Cimage_header& oHeader)
{
  (void) nInitValues(oHeader);
}

Csource::~Csource()
{
  if (NULL != m_poWavelength)
    {
      delete m_poWavelength;
      m_poWavelength = NULL;
    }
  m_eThe_State = eSource_unknown_state; 
}

int Csource::nInitValues()
{
  m_eThe_State               = eSource_unknown_state; 
  m_sDescription             = "unknown";
  m_sSource_key              = "";

  m_poWavelength             = new Cwavelength();

  m_a2fRotation[0]           = 0.0;   // Goes with m_a3x3fVector[1].
  m_a2fRotation[1]           = 0.0;
  m_a3x3fVector[0][0]        = 0.0;   // m_a3x3fVector[0][] should prob. be m_a3x3fVector[2][]
  m_a3x3fVector[0][1]        = 0.0;   // Then m_a2fRotation[0] would match 
  m_a3x3fVector[0][2]        = 1.0;
  m_a3x3fVector[1][0]        = 0.0;
  m_a3x3fVector[1][1]        = 1.0;
  m_a3x3fVector[1][2]        = 0.0;
  m_a3x3fVector[2][0]        = 1.0;
  m_a3x3fVector[2][1]        = 0.0;
  m_a3x3fVector[2][2]        = 0.0;
  m_a2fSpectralDispersion[0] = 0.0;
  m_a2fSpectralDispersion[1] = 0.0;
  m_a4fCrossfire[0]          = 0.0;
  m_a4fCrossfire[1]          = 0.0;
  m_a4fCrossfire[2]          = 0.0;
  m_a4fCrossfire[3]          = 0.0;
  m_a4fPolarz[0]             = 0.5;
  m_a4fPolarz[1]             = 1.0;
  m_a4fPolarz[2]             = 0.0;
  m_a4fPolarz[3]             = 0.0;
  m_fIntensity               = 0.0;
  m_a4fSize[0]               = 0.0;
  m_a4fSize[1]               = 0.0;
  m_a4fSize[2]               = 0.0;
  m_a4fSize[3]               = 0.0;

  return (0);
}

int Csource::nInitValues(Cimage_header& oHeader)
{
  int    nStat, nTemp;
  Cstring sTemp;

  (void) nInitValues();

  if (!oHeader.bIsAvailable()) 
    {
      cout << "Csource ERROR: image header is not valid."
           << "  Cannot construct source!\n";
      nStat = 1;
    }
  else if (0 == oHeader.nGetValue(ms_sSourceKey, &sTemp)) 
    {
      cout << "Csource ERROR: cannot use database key yet!\n";
      nStat = 2;
  }
  else 
    {
      // Try to get required information from the header

      nStat = 0;

      if (NULL != m_poWavelength)
	delete m_poWavelength;

      m_poWavelength  = new Cwavelength(oHeader, ms_sSourcePrefix);
  
      if (eWavelength_unknown_type == m_poWavelength->m_eThe_Type) 
	{
	  cout << "Csource WARNING: cannot construct wavelength!\n";
	  nStat++;
	}

      nTemp = oHeader.nGetValue(ms_sSourceVectors, 9, &m_a3x3fVector[0][0]);
      nTemp = oHeader.nGetValue(ms_sSourceValues, 2, m_a2fRotation);

      // For the above, when not present, use the defaults, so no error.

      nTemp = oHeader.nGetValue(ms_sSourcePolarz, 4, m_a4fPolarz);
      if (0 != nTemp) 
	{
	  cout << "Csource WARNING: no " << ms_sSourcePolarz 
               << " keyword!\n";
	  nStat++;
	}

      nTemp = oHeader.nGetValue(ms_sSourceSpectralDispersion, 2, 
				m_a2fSpectralDispersion);
      if (0 != nTemp)
	{
	  cout << "Csource WARNING: no " << ms_sSourceSpectralDispersion 
               << " keyword!\n"; 
	  nStat++;
	}

      nTemp = oHeader.nGetValue(ms_sSourceCrossfire, 4, m_a4fCrossfire);
      if (0 != nTemp) 
	{
	  cout << "Csource WARNING: no " << ms_sSourceCrossfire
               << " keyword!\n";
	  nStat++;
	}

      nTemp = oHeader.nGetValue(ms_sSourceSize, 4, m_a4fSize);
      if (0 != nTemp) 
	{
	  cout << "Csource WARNING: no " << ms_sSourceSize 
               << " keyword!\n";
	  nStat++;
	}
    }

  if (0 != nStat)
    {
      // If there are errors, no sense have this around
      m_eThe_State = eSource_unknown_state;
      return (1);
    }
  else
    {
      m_eThe_State = eSource_notrefined_state;
      return (0);
    }
}

int Csource::nList(const int nFlag)
{
  cout << "\nSource listing: \n\n";
  m_poWavelength->nList(nFlag);
  
  if (1 > nFlag)
    {
      cout << "   Direction vector: " << m_a3x3fVector[0][0] << ", "
                                      << m_a3x3fVector[0][1] << ", "
                                      << m_a3x3fVector[0][2] << endl;
      cout << "    Rotation values: " << m_a2fRotation[0] << ", "
                                      << m_a2fRotation[1] << endl;
      cout << "Spectral Dispersion: " << m_a2fSpectralDispersion[0] << ", "
                                      << m_a2fSpectralDispersion[1] << endl;

      cout << "          Crossfire: " << m_a4fCrossfire[0] << ", "
                                      << m_a4fCrossfire[1] << ", "
                                      << m_a4fCrossfire[2] << ", "
                                      << m_a4fCrossfire[3] << endl;
    }

  cout << "       Polarization: " << m_a4fPolarz[0] << ", "
                                  << m_a4fPolarz[1] << ", "
                                  << m_a4fPolarz[2] << ", "
                                  << m_a4fPolarz[3] << endl;

  if (1 > nFlag)
    {
      cout << "               Size: " << m_a4fSize[0] << ", "
                                      << m_a4fSize[1] << ", "
                                      << m_a4fSize[2] << ", "
                                      << m_a4fSize[3] << endl;
    }

  cout << "          Intensity: " << m_fIntensity << endl;

  cout << flush;
  return (0);
}

void Csource::vCalcGetS0(float *pfS0,float *pfDeriv)
{
    double pfS0_[3];
    double a2x3x3fDeriv[2][3][3];
    double* pfDerivPoint;

    if (pfDeriv)
        pfDerivPoint = &a2x3x3fDeriv[0][0][0];
    else
        pfDerivPoint = NULL;

    vCalcGetS0(pfS0_,pfDerivPoint);
    vCopyVec3D(pfS0_,pfS0);
    if (pfDeriv) 
        vCopyVecND(2*3*3,pfDerivPoint,pfDeriv);

};

void Csource::vCalcGetS0(double *pfS0,double *pfDeriv) {

  // Return in *pfS0 a 3D vector which describes the source Vector after
  // applying any rotations specified in m_a3x3fVector[1][*] and [2][*] and
  // normalizing it.
  // If the derivatives are requested, then fill these in as well.

  double fMat[3][3];
  double pfDerivDummy[3][3];    // If we request derivatives, all derivatives must be requested.  
                                // But we don't want the third one, so store it here.
  double m_a3x3fVector_[3][3];

  vCopyMat3D(&m_a3x3fVector[0][0],&m_a3x3fVector_[0][0]);

  vConv3Rot3Vec3DMat3D((double) m_a2fRotation[0], // Usually Rot Y
                       (double) m_a2fRotation[1], // Usually Rot X
                       (double) 0.0,              // Rotation around Z (i.e. source vector) is always 0. 
                       &m_a3x3fVector_[1][0],     // 2nd source vector, usually 0 1 0
                       &m_a3x3fVector_[2][0],     // 3rd source vector, usually 1 0 0
                       &m_a3x3fVector_[0][0],     // 1st source vector, usually 0 0 1  
                       (double*)fMat,
                       (pfDeriv), (pfDeriv)?(pfDeriv+9):NULL,(pfDeriv)?(&pfDerivDummy[0][0]):NULL
               );
  vMulMat3DVec3D(fMat, &m_a3x3fVector_[0][0], pfS0);
  (void) fNormVec3D(pfS0);
};

void Csource::vCalcGetUnrotatedS0(double *pfS0) {
    vCopyVec3D(&m_a3x3fVector[0][0],pfS0);
};

void Csource::vGetPolarz(float *pfPol)
{
  (void) fNormVec3D(&m_a4fPolarz[1]);

  for (int i=0; i < 4; i++) pfPol[i] = m_a4fPolarz[i];
}

int
Csource::nUpdateHeader(Cimage_header* poHeader)
{
  // What about the wavelength?

  if ( (NULL == poHeader) ||  (!poHeader->bIsAvailable()))
    return (-1);
  if (NULL != m_poWavelength)
    {
      (void) m_poWavelength->nUpdateHeader(poHeader, ms_sSourcePrefix);
    }
  (void) poHeader->nReplaceValue(ms_sSourceVectors, 9, &m_a3x3fVector[0][0], 3);
  (void) poHeader->nReplaceValue(ms_sSourceValues, 2, m_a2fRotation, 4);
  (void) poHeader->nReplaceValue(ms_sSourcePolarz, 4, m_a4fPolarz);
  (void) poHeader->nReplaceValue(ms_sSourceSpectralDispersion, 2,
				 m_a2fSpectralDispersion, 5);
  (void) poHeader->nReplaceValue(ms_sSourceCrossfire, 4, m_a4fCrossfire, 4);
  (void) poHeader->nReplaceValue(ms_sSourceSize, 4, m_a4fSize, 4);

  return (0);
}

void
Csource::vSetRotation(const float fRot0, const float fRot1)
{
  m_a2fRotation[0] = fRot0;
  m_a2fRotation[1] = fRot1;
};

int
Csource::nDiff(Csource& oSourceAdd,Csource& oSourceSubtract) {
    int nx;
    for (nx=0;nx<2;nx++)
        m_a2fRotation[nx] += oSourceAdd.m_a2fRotation[nx] - oSourceSubtract.m_a2fRotation[nx];
    return 0;
};
