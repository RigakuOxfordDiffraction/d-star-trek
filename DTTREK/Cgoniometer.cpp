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
// Cgoniometer.cc         Initial author: J.W. Pflugrath           01-May-1995
//    This file contains the member functions of class Cgoniometer.
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
#include "Cgoniometer.h"        // Class definition and prototypes
                                //  Cgoniometer.h includes others
#include "dtrekdefs.h"
//+Definitions, constants, and initialization of static member variables

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

Cstring Cgoniometer::ms_sGonioNumValues   = D_K_GonioNumValues;
Cstring Cgoniometer::ms_sGonioKey         = D_K_GonioKey;
Cstring Cgoniometer::ms_sGonioNames       = D_K_GonioNames;
Cstring Cgoniometer::ms_sGonioUnits       = D_K_GonioUnits;
Cstring Cgoniometer::ms_sGonioVectors     = D_K_GonioVectors;
Cstring Cgoniometer::ms_sGonioValues      = D_K_GonioValues;
Cstring Cgoniometer::ms_sGonioValuesMin   = D_K_GonioValuesMin;
Cstring Cgoniometer::ms_sGonioValuesMax   = D_K_GonioValuesMax;
Cstring Cgoniometer::ms_sGonioCollisionOffsets   = D_K_GonioCollisionOffsets;
Cstring Cgoniometer::ms_sGonioDescription = D_K_GonioDescription;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cgoniometer::Cgoniometer()
{
  (void) nInitValues();
}

Cgoniometer::Cgoniometer(const int nNumValues,
			 const float *pfAxes, const float *pfValues,
			 const Cstring *psNames, const Cstring *psUnits)
{
  (void) nInitValues();
  int i;

  m_nNumValues       = nNumValues;
  m_pfCurrentValue   = new float[m_nNumValues];

  m_pfHardwareMinLimit = new float[m_nNumValues];
  m_pfHardwareMaxLimit = new float[m_nNumValues];

  m_pfDatumCollisionOffset = new float[m_nNumValues];

  m_pfRequestedValue = new float [m_nNumValues];
  m_pfDatumValue     = new float [m_nNumValues];
  m_pfSetValue       = new float [m_nNumValues];
  m_psName           = new Cstring [m_nNumValues];
  m_psUnits          = new Cstring [m_nNumValues];
  m_pfRequestedRate  = new float [m_nNumValues];
  m_pfSetRate        = new float [m_nNumValues];

  for (i = 0; i < m_nNumValues; i++)
    {
      m_pfCurrentValue[i]   = pfValues[i];
      m_pfRequestedValue[i] = pfValues[i];
      m_pfDatumValue[i]     = pfValues[i];
      m_pfSetValue[i]       = pfValues[i];
      m_psName[i]           = psNames[i];
      m_psUnits[i]          = psUnits[i];
      m_pfRequestedRate[i]  = 0.0;
      m_pfSetRate[i]        = 0.0;
    }
  m_pfVector         = new float [3*m_nNumValues];
  for (i = 0; i < 3 * m_nNumValues; i++)
    {
      m_pfVector[i] = pfAxes[i];
    }
  int j = 0;
  for (i = 0; i < m_nNumValues; i++)
    {
      j = 3 * i;
      if (0.0 >= fLenVec3D(&m_pfVector[j]))
	{
	  cout << "ERROR, Cgoniometer: the vector " << m_sRelPrefix 
               << ms_sGonioVectors << " " 
               << i << " length <= 0 [" 
               << m_pfVector[j] << ", " << m_pfVector[j+1] << ", "
               << m_pfVector[j+2] << "]!\n"
               << flush;
	}
    }

  m_eThe_State       = eGoniometer_unknown_state;
  m_eThe_Type        = eGoniometer_other_type;

  vGetNumRotTransValues();
}

Cgoniometer::Cgoniometer(Cimage_header& oHeader, const Cstring& sPre)
{
  (void) nInitValues(oHeader, sPre);
}

Cgoniometer::~Cgoniometer()
{
    if( NULL != m_pfCurrentValue )
    {
        delete [] m_pfCurrentValue;
        m_pfCurrentValue = NULL;
    }
  
    if( NULL != m_pfHardwareMinLimit )
    {
        delete [] m_pfHardwareMinLimit;
        m_pfHardwareMinLimit = NULL;
    }
  
    if( NULL != m_pfHardwareMaxLimit )
    {
        delete [] m_pfHardwareMaxLimit;
        m_pfHardwareMaxLimit = NULL;
    }
  
    if( NULL != m_pfDatumCollisionOffset )
    {
        delete [] m_pfDatumCollisionOffset;
        m_pfDatumCollisionOffset = NULL;
    }
    
    if (NULL != m_pfRequestedValue)
    {
      delete [] m_pfRequestedValue;
      m_pfRequestedValue = NULL;
    }
  if (NULL != m_pfDatumValue)
    {
      delete [] m_pfDatumValue;
      m_pfDatumValue = NULL;
    }
  if (NULL != m_pfSetValue)
    {
      delete [] m_pfSetValue;
      m_pfSetValue = NULL;
    }
  if (NULL != m_psName)
    {
      delete [] m_psName;
      m_psName = NULL;
    }
  if (NULL != m_psUnits)
    {
      delete [] m_psUnits;
      m_psUnits = NULL;
    }
  if (NULL != m_pfRequestedRate)
    {
      delete [] m_pfRequestedRate;
      m_pfRequestedRate = NULL;
    }
  if (NULL != m_pfSetRate)
    {
      delete [] m_pfSetRate;
      m_pfSetRate = NULL;
    }
  if (NULL != m_pfVector)
    {
      delete [] m_pfVector;
      m_pfVector = NULL;
    }
  m_eThe_State = eGoniometer_unknown_state;
  m_eThe_Type  = eGoniometer_unknown_type;
  m_nNumValues = 0;
}

int Cgoniometer::nInitValues()
{
  m_eThe_State       = eGoniometer_unknown_state;
  m_eThe_Type        = eGoniometer_unknown_type;
  m_nNumValues       = 0;
  m_nNumRotValues    = 0;
  m_nNumTransValues  = 0;
  m_pfCurrentValue   = NULL;
  m_pfHardwareMinLimit = NULL;
  m_pfHardwareMaxLimit = NULL;
  m_pfDatumCollisionOffset = NULL;
  m_pfRequestedValue = NULL;
  m_pfDatumValue     = NULL;
  m_pfSetValue       = NULL;
  m_psName           = NULL;
  m_psUnits          = NULL;
  m_pfRequestedRate  = NULL;
  m_pfSetRate        = NULL;
  m_pfVector         = NULL;
  m_sKey             = "";
  m_sDescription     = "UNKNOWN";
  m_sPrefix          = "";
  m_sRelPrefix       = "";
  memset(&m_anGonioOffsetTableEntries[0],0,sizeof(int)*(nMaxGonioRotationOffsets+1));
  memset(&m_anGonioOffsetTableRotAxis[0],0,sizeof(int)*nMaxGonioRotationOffsets);
  memset(&m_afGonioOffsetTableGonio[0],0,sizeof(double)*nMaxGonioRotationOffsets*5);
  memset(&m_afGonioOffsetTableOffset[0],0,sizeof(double)*nMaxGonioRotationOffsets*5);

  return (0);
}

int Cgoniometer::nInitValues(Cimage_header& oHeader, const Cstring& sPre)
{
  int    i;       // Loop counter
  int    nStat, nTemp;
  Cstring sTemp;

  (void) nInitValues();

  m_sPrefix = sPre;
  if (m_sPrefix.contains(D_K_CrystalPrefix))
      m_sRelPrefix = m_sPrefix.from(D_K_CrystalPrefix);
  else
      m_sRelPrefix = sPre;

  if (!oHeader.bIsAvailable())
    {
      cout << "Cgoniometer ERROR: image header is not valid."
           << "  Cannot construct goniometer!\n";
      nStat = 1;
    }
  else if (0 == oHeader.nGetValue(m_sRelPrefix + ms_sGonioKey, &sTemp))
    {
      cout << "Cgoniometer ERROR: cannot use database key yet!\n";
      nStat = 2;
    }
  else
    {
      // Try to get required information from the header

      nStat = 0;

      nTemp = oHeader.nGetValue(m_sRelPrefix + ms_sGonioNumValues, &m_nNumValues);
      if (0 != nTemp)
	{
	  cout << "Cgoniometer ERROR: reading " << m_sRelPrefix
	       <<  ms_sGonioNumValues << " keyword!\n";
	  m_nNumValues = 0;
	  return (1);                // This is fatal!
	}

      // Allocate space for all variables
      m_pfCurrentValue   = new float[m_nNumValues];
      
      m_pfHardwareMinLimit = new float[m_nNumValues];
      m_pfHardwareMaxLimit = new float[m_nNumValues];
      
      m_pfDatumCollisionOffset = new float[m_nNumValues];

      m_pfRequestedValue = new float[m_nNumValues];
      m_pfDatumValue     = new float[m_nNumValues];
      m_pfSetValue       = new float[m_nNumValues];
      m_psName           = new Cstring[m_nNumValues];
      m_psUnits          = new Cstring[m_nNumValues];
      m_pfRequestedRate  = new float[m_nNumValues];
      m_pfSetRate        = new float[m_nNumValues];
      
      for (i = 0; i < m_nNumValues; i++)
	{
	  m_pfCurrentValue[i]   = -9999.0;
      
      m_pfHardwareMinLimit[i] = -9999.0;
      m_pfHardwareMaxLimit[i] = -9999.0;
      
      m_pfDatumCollisionOffset[i] = -9999.0;
	  
      m_pfRequestedValue[i] = 0.0;
	  m_pfDatumValue[i]     = 0.0;
	  m_pfSetValue[i]       = 0.0;
	  m_psName[i]           = "UNKNOWN";
	  m_psUnits[i]          = "UNKNOWN";
	  m_pfRequestedRate[i]  = 0.0;
	  m_pfSetRate[i]        = 0.0;
	}
      m_pfVector         = new float [3*m_nNumValues];

      nTemp = oHeader.nGetValue(m_sRelPrefix + ms_sGonioNames,
				m_nNumValues, m_psName);
      if (0 != nTemp)
	{
	  cout << "Cgoniometer ERROR: reading " << m_sRelPrefix
               << ms_sGonioNames << " keyword!\n";
	  nStat++;
	}

      nTemp = oHeader.nGetValue(m_sRelPrefix + ms_sGonioUnits,
				m_nNumValues, m_psUnits);
      if (0 != nTemp)
	{
	  cout << "Cgoniometer ERROR: reading " << m_sRelPrefix
	       << ms_sGonioUnits << " keyword!\n";
	  nStat++;
	}

      nTemp = oHeader.nGetValue(m_sRelPrefix + ms_sGonioVectors,
				3*m_nNumValues, m_pfVector);
      if (0 != nTemp)
	{
	  cout << "Cgoniometer ERROR: reading " << m_sRelPrefix
	       << ms_sGonioVectors << " keyword!\n";
	  nStat++;
	}

      // It is not required to have the following keyword, so
      // do not flag if not found

      nTemp = oHeader.nGetValue(m_sRelPrefix + ms_sGonioDescription,
				&m_sDescription);

      nTemp = oHeader.nGetValue(m_sPrefix + ms_sGonioValues,
				m_nNumValues, m_pfDatumValue);
      
//+JWP 2011-11-04
      // Special case for when 2theta is wrong
      if (D_K_RxPrefix == m_sPrefix)
	if (-99999.0 == m_pfDatumValue[1]) 
	  m_pfDatumValue[1] = 0.0;
//-JWP 2011-11-04

      if (0 != nTemp)
	    {
	      cout << "ERROR Cgoniometer: reading " << m_sPrefix
	           << ms_sGonioValues << " keyword!\n";
	      nStat++;
	    }
    
      /////////////////////////////////////////////////////////
      // It is not required to have the following keywords, so
      // do not flag if not found

      nTemp = oHeader.nGetValue(m_sPrefix + ms_sGonioValuesMin,
				m_nNumValues, m_pfHardwareMinLimit);
      
      nTemp = oHeader.nGetValue(m_sPrefix + ms_sGonioValuesMax,
				m_nNumValues, m_pfHardwareMaxLimit);
     
      nTemp = oHeader.nGetValue(m_sPrefix + ms_sGonioCollisionOffsets,
				m_nNumValues, m_pfDatumCollisionOffset);
      /////////////////////////////////////////////////////////// 
    }
  int j = 0;
  for (i = 0; i < m_nNumValues; i++)
    {
      if (0.0 >= fLenVec3D(&m_pfVector[j]))
	{
	  cout << "ERROR Cgoniometer: the vector " << m_sRelPrefix 
               << ms_sGonioVectors << " " 
               << i << " length <= 0 [" 
               << m_pfVector[j] << ", " << m_pfVector[j+1] << ", "
               << m_pfVector[j+2] << "]!\n"
               << flush;
	  nStat++;
	}
      j = j + 3;
    }
  vGetNumRotTransValues();

  if (0 == nStat)
      nStat += nLoadOffsets(&oHeader);    

  if (0 != nStat)
    {
      // If there are errors, no sense having this around

      cout << "ERROR Cgoniometer::nInitValues(Cimage_header,...): object not valid!\n" << flush;
      m_eThe_State = eGoniometer_unknown_state;
      return (1);
    }
  else
    {
      m_eThe_State = eGoniometer_not_ready_to_move_state;
      if (1 == m_nNumValues)
	m_eThe_Type = eGoniometer_single_type;
      else if (3 > m_nNumValues)
	m_eThe_Type = eGoniometer_other_type;
      else if (3 == m_nNumValues)
	// To Do: test for kappa or eulerian here:
	m_eThe_Type = eGoniometer_eulerian_type;
      else
	m_eThe_Type = eGoniometer_other_type;
      return (0);
    }
}

int
Cgoniometer::nList(const int nFlag)
{
  printf("\n%s Goniometer listing: \n\n",  m_sPrefix.string());
  printf("        Description: %s\n", m_sDescription.string());
  printf("   Number of values: %d\n\n", m_nNumValues);
  
  printf("         Name    Datum   Current  Units   Vector\n");
  printf("=====================================================================\n");
/*
Name             Datum   Current  Units   Vector
====================================================================
Omega         -360.000  -360.000    deg   (  1.000,   0.000,   0.000)
1234567890123 12345678  12345678 123456 12
*/
  Cstring sTemp;
  for (int i = 0; i < m_nNumValues; i++)
    {
      sTemp = m_psName[i].substr(0, min((int)m_psName[i].length(), 13));
      printf("%13s", sTemp.string());
      printf(" %8.3f  ", m_pfDatumValue[i]);
      if (-9999.0 != m_pfCurrentValue[i])
	printf("%8.3f", m_pfCurrentValue[i]);
      else
	printf("%8s", "Unknown");
      printf(" %6s   (%7.3f, %7.3f, %7.3f)\n", m_psUnits[i].string(),
	     m_pfVector[3*i], m_pfVector[3*i + 1], m_pfVector[3*i + 2]);
    }
  printf("=====================================================================\n");
  fflush(stdout);
  return (0);
}

void
Cgoniometer::vGetNumRotTransValues() {
    int nRots;
    nRots = 0;
    for (nRots = 0; nRots < m_nNumValues; nRots++) {
        if ("deg" != m_psUnits[nRots])
            break;
    };
    m_nNumRotValues = nRots;
    m_nNumTransValues = m_nNumValues - m_nNumRotValues;
};


int
Cgoniometer::nGetDatum(const int nNum, float *pfDatum)
{
  int i, j;
  j = min(nNum, m_nNumValues);
  for (i = 0; i < j; i++)
    {
      pfDatum[i] = m_pfDatumValue[i];
    }
  return (m_nNumValues-nNum);
}

int
Cgoniometer::nGetDatum(const int nNum, double *pfDatum)
{
  int i, j;
  j = min(nNum, m_nNumValues);
  for (i = 0; i < j; i++)
    {
      pfDatum[i] = (double)m_pfDatumValue[i];
    }
  return (m_nNumValues-nNum);
}

int
Cgoniometer::nSetDatum(const int nNum, const float *pfDatum)
{
  int i, j;
  j = min(nNum, m_nNumValues);
  for (i = 0; i < j; i++)
    {
      m_pfDatumValue[i] = pfDatum[i];
    }
  return (m_nNumValues-nNum);
}
int
Cgoniometer::nSetDatum(const int nNum, const double *pfDatum)
{
  int i, j;
  j = min(nNum, m_nNumValues);
  for (i = 0; i < j; i++)
    {
      m_pfDatumValue[i] = (float)pfDatum[i];
    }
  return (m_nNumValues-nNum);
}

int
Cgoniometer::nGetVectors(const int nNum, float *pfVec)
{
  int i, j, k;
  k = min(nNum, m_nNumValues);
  for (i = 0; i < k; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  pfVec[3*i+j] = m_pfVector[3*i+j];
	}
    }
  return (m_nNumValues - nNum);
}
int
Cgoniometer::nGetVectors(const int nNum, double *pdVec)
{
  int i, j, k;
  k = min(nNum, m_nNumValues);
  for (i = 0; i < k; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  pdVec[3*i+j] = (double)m_pfVector[3*i+j];
	}
    }
  return (m_nNumValues - nNum);
}

int
Cgoniometer::nGetNames(const int nNum, Cstring *psNames)
{
  int i, k;
  k = min(nNum, m_nNumValues);
  for (i = 0; i < k; i++)
    {
      psNames[i] = m_psName[i];
    }
  return (m_nNumValues-nNum);
}

void
Cgoniometer::vCalcGetRotMatrix(float *pfMatrix, int nVecNum,
                               const float fRot1, const float fRot2,
                               const float fRot3, const float fRot4)
{
  int nx;
  float  *pfTemp;
  double *pdTemp;
  double a9dTemp[9];
  pdTemp = &a9dTemp[0];
  pfTemp = pfMatrix;
  for (nx = 0; nx < 9; nx++)
    {
      *pfTemp++ = 0.0;
      *pdTemp++ = 0.0;
    }
  

  double a4dTemp[4];
  a4dTemp[0] = fRot1;
  a4dTemp[1] = fRot2;
  a4dTemp[2] = fRot3;
  a4dTemp[3] = fRot4;
  vCalcGetRotMatrix(&a9dTemp[0], nVecNum, a4dTemp[0], a4dTemp[1],
		    a4dTemp[2], a4dTemp[3]);

  pdTemp = &a9dTemp[0];
  pfTemp = pfMatrix;
  for (nx = 0; nx < 9; nx++)
    *pfTemp++ = (float) *pdTemp++;
}


void
Cgoniometer::vCalcGetRotMatrix(double *pfMatrix,int nVecNum,
                               const double fRot1, const double fRot2,
                               const double fRot3, const double fRot4)
{
    // This assumes the first m_nNumRotValues goniometer values are rotations
    // Change later so that it looks at Units (see GetTransVector)
    
    int    i;
    double fRot[4];
    double fVec[4][3];
    double fRotIn[4];
    float  afGonioRotTemp[5];
   
    fRotIn[0] = fRot1;
    fRotIn[1] = fRot2;
    fRotIn[2] = fRot3;
    fRotIn[3] = fRot4;
    
    
    // Get the rotation matrix defined by rotations of fRot1, fRot2 and fRot3
    // around the goniometer axes or if these are -999 by
    // first 3 datum values
    
    int n1 = 0;
    if (4 == m_nNumRotValues)
    {
        n1 = 1;
    }
    int nx = 0;

    for (i = 0; i < m_nNumRotValues;i++)
        afGonioRotTemp[i] = m_pfDatumValue[i];

    if (nVecNum>=0) {
        nApplyAllGonioOffsets(&afGonioRotTemp[0],nVecNum);
    }    

    for (i = n1; i < m_nNumRotValues; i++)
    {
        if (m_nNumValues > i)
        {
            if (-998.0 >= fRotIn[i])
                fRot[nx]  = afGonioRotTemp[i];    // The value at datum
            else
                fRot[nx]  = fRotIn[i];          // The input value
            fVec[nx][0] = m_pfVector[3*i];      // around this axis
            fVec[nx][1] = m_pfVector[3*i+1];
            fVec[nx][2] = m_pfVector[3*i+2];
        }
        else
        {
            fRot[nx]    = 0.0;
            fVec[nx][0] = 1.0;
            fVec[nx][1] = 0.0;
            fVec[nx][2] = 0.0;
        }
        nx++;
    }
    
    // On goniometers, the last axis (i.e. phi)  is applied first.
    
    vConv3Rot3Vec3DMat3D(fRot[2], fRot[1], fRot[0],
        &fVec[2][0], &fVec[1][0], &fVec[0][0],
        pfMatrix);
    if (4 == m_nNumRotValues)
    {
        // Apply the 4th rotation, closest to goniometer baseplate
        // to the pfMatrix
        
        double a3x3fMat0[3][3];
        double a3x3fMat1[3][3];
        double a3x3fMat2[3][3];
        
        fRot[0]    = fRot1;
        fVec[0][0] = m_pfVector[0];      // around this axis
        fVec[0][1] = m_pfVector[1];
        fVec[0][2] = m_pfVector[2];
        vConvRotVec3DMat3D(fRot[0], &fVec[0][0], &a3x3fMat0[0][0]);
        vCopyMat3D(pfMatrix, &a3x3fMat1[0][0]);
        vMulMat3DMat3D(a3x3fMat0, a3x3fMat1, a3x3fMat2);
        vCopyMat3D(&a3x3fMat2[0][0], pfMatrix);
    }
}

int Cgoniometer::nCalcGetTransVector(float *pfVecOut)
{
    int   i;
    float fTemp[3];
    int   nStat;
    
    // Initialize translation vector to 0.0.  If no
    // translation found, then this will be returned
    
    pfVecOut[0] = 0.0;
    pfVecOut[1] = 0.0;
    pfVecOut[2] = 0.0;
    
    nStat = 1;
    for (i = 0; i < m_nNumValues; i++)
    {
        if (m_psUnits[i] == "mm")
        {
            
            // A goniometer value with mm units was found
            
            vMulVec3DScalar(&m_pfVector[3*i], m_pfDatumValue[i], fTemp);
            vAddVec3DVec3D(pfVecOut, fTemp, pfVecOut);
            nStat = 0;      // At least 1 translation was found
        }
    }
    return (nStat);
}

int Cgoniometer::nCalcGetTransVectorDeriv(int nTransComp,double* pfVecOut)
{
  
  int   i;
  int   nStat;
  int   nComp;

  // Initialize translation vector to 0.0.  If no
  // translation found, then this will be returned

  pfVecOut[0] = 0.0;
  pfVecOut[1] = 0.0;
  pfVecOut[2] = 0.0;

  nStat = 1;
  nComp = 0;
  for (i = 0; i < m_nNumValues; i++)  {
      if (m_psUnits[i] == "mm") {
          if (nComp == nTransComp) {
              vCopyVec3D(&m_pfVector[3*i],pfVecOut);
              nStat = 0;      // At least 1 translation was found
          };
          nComp++;
      };
  };
  return (nStat);
};

int Cgoniometer::nSetTransVector(const int nDim, const float *pfVecIn)
{
  int   i, j;

  float a3x3fTransBasis[3][3];
  float a3x3fTransBasisInv[3][3];
  float a3fVecOut[3];
  if (nDim!= 3)
      return 1;
  // Is the translation vector rotated?  If so, we un-rotate it first.

  // Build the matrix of basis vectors.
  for (i = 0, j= 0; i < m_nNumValues; i++) {
      if (m_psUnits[i] == "mm") {
          // A goniometer value with mm units was found
          vCopyVec3D(&m_pfVector[3*i],a3x3fTransBasis[j++]);
      };
  };
  if (j!=3)
      return 1;
  // Can we assume we have ones and zeros?  If so, use vTranMat3D()
  fInvMat3D(&a3x3fTransBasis[0][0],&a3x3fTransBasisInv[0][0]);
  vMulMat3DVec3D(a3x3fTransBasisInv,pfVecIn,a3fVecOut);

  // Stuff the values back into the m_pfDatumValue array.
  for (i = 0, j= 0; i < m_nNumValues; i++) {
      if (m_psUnits[i] == "mm") {
          m_pfDatumValue[i] = a3fVecOut[j++];
      };
  };

  return 0;
}

int
Cgoniometer::nUpdateHeader(Cimage_header* poHeader, const Cstring& sPre)
{
  int     nStat;
  Cstring sPrefix;
  Cstring sRelPrefix;
  Cstring sTemp;
  int     nx;
  sPrefix = sPre;

  if (sPrefix.contains(D_K_CrystalPrefix))
      sRelPrefix = sPrefix.from(D_K_CrystalPrefix);
  else
      sRelPrefix = sPre;

  nStat = poHeader->nReplaceValue(sRelPrefix + ms_sGonioNames,
				  m_nNumValues, m_psName);
  nStat = nStat + poHeader->nReplaceValue(sRelPrefix + ms_sGonioUnits,
					  m_nNumValues, m_psUnits);

  // A little code to write out the gonio vectors in fewer decimals places
  // for values of -1, 0, 1 and more decimal places for other values.

  sTemp = "";
  float *pfTemp;
  pfTemp = m_pfVector;
  for (nx = 0; nx < m_nNumValues*3; nx++, pfTemp++)
    {
      if ( (0.0 == *pfTemp) || (-1.0 == *pfTemp) || (1.0 == *pfTemp) )
	sTemp += Cstring(*pfTemp, 0, 1) + ' ';
      else
	sTemp += Cstring(*pfTemp, 0, 5) + ' ';
    }

  nStat = nStat + poHeader->nReplaceValue(sRelPrefix + ms_sGonioVectors,
					  sTemp);

  nStat = nStat + poHeader->nReplaceValue(sRelPrefix + ms_sGonioDescription,
					  m_sDescription);
  nStat = nStat + poHeader->nReplaceValue(sPrefix + ms_sGonioValues,
					  m_nNumValues, m_pfDatumValue, 4);
  return (nStat);
}

float Cgoniometer::fGetDistance(void)
{
  // Assume the largest translation is the crystal to detector distance

  int   i;
  for (i = 0; i < m_nNumValues; i++)
    {
      if (m_psName[i].contains("ist"))
	{
	  // A goniometer value with "ist" in the name was found

	  return (m_pfDatumValue[i]);
	}
    }

  // One was not found

  return  0.0f;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Get a hardware axis limits by the axis name
bool Cgoniometer::bGetHardwareLimits(const char* ccAxisName, float& fMin, float& fMax)
{
    int   i = 0;
    for(i = 0; i < m_nNumValues; i++)
    {
        if( m_psName[i].contains(ccAxisName) )
        {
            // A goniometer value with "ist" in the name was found
            if( m_pfHardwareMinLimit[i] == -9999.0 || m_pfHardwareMaxLimit[i] == -9999.0 )
                return false; // not supplied
            
            fMin = m_pfHardwareMinLimit[i];     
            fMax = m_pfHardwareMaxLimit[i];     

            return true;
        }
    }

    // Named info was not found
    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////
bool Cgoniometer::bGetHardwareLimits(const char* ccAxisName, CSegment& segLimits)
{
    float   fMin = 0.0f;
    float   fMax = 0.0f;

    if( !bGetHardwareLimits(ccAxisName, fMin, fMax) )
        return false; 
    
    segLimits.vSet((double)fMin, (double)fMax);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////
// Get a hardware axis collision offset by the axis name
bool Cgoniometer::bGetCollisionOffset(const char* ccAxisName, float& fOffset)
{
    int   i = 0;
    for(i = 0; i < m_nNumValues; i++)
    {
        if( m_psName[i].contains(ccAxisName) )
        {
            if( m_pfDatumCollisionOffset[i] == -9999.0 )
                return false; // not supplied
            
            fOffset = m_pfDatumCollisionOffset[i];     

            return true;
        }
    }

    // Named info was not found
    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////


float
Cgoniometer::fGetSwing(void)
{
  int   i;
  for (i = 0; i < m_nNumValues; i++)
    {
      if (m_psName[i].contains("wing"))
	{
	  // A goniometer value with "wing" in the name was found

	  return (m_pfDatumValue[i]);
	}
      else if (m_psName[i].contains("heta"))
	{
	  // A goniometer value with "heta" in the name was found

	  return (m_pfDatumValue[i]);
	}
    }

  // One was not found

  return (0.0);
}

int
Cgoniometer::nSetDistance(const float fDist)
{
  int   i;
  for (i = 0; i < m_nNumValues; i++)
    {
      if (m_psName[i].contains("ist"))
	{
	  // A goniometer value with "ist" in the name was found

	  m_pfDatumValue[i] = fDist;
	  return (0);
	}
    }

  // One was not found

  return (1);
}

int
Cgoniometer::nSetSwing(const float fSwing)
{
  int   i;
  for (i = 0; i < m_nNumValues; i++)
    {
      if (m_psName[i].contains("wing"))
	{
	  // A goniometer value with "wing" in the name was found

	  m_pfDatumValue[i] = fSwing;
	  return (0);
	}
      else if (m_psName[i].contains("heta"))
	{
	  // A goniometer value with "heta" in the name was found

	  m_pfDatumValue[i] = fSwing;
	  return (0);
	}
    }

  // One was not found

  return (1);
}

int
Cgoniometer::nGetRotVector(const int nVecNum, float *pfVec)
{
  // Get the vector specified by the axis number nVecNum when the goniometer
  // is at the datum position.

  if ( (0 > nVecNum) || (m_nNumValues <= nVecNum) || (3 < nVecNum) )
    {
      // Out of range

      return (-1);
    }

  int   i;
  float a3fRot[4], a3x3fRotMat[3][3];

  for (i = 0; i < 4; i++)
    {
      a3fRot[i] = 0.0;                    // Use 0 for all rotations...
      if ( (i < nVecNum) && (i < m_nNumValues) )
	{                                // ... unless the rotation PRECEDES
	  a3fRot[i] = m_pfDatumValue[i]; //    our vector of interest
	
	}
    }

  vCalcGetRotMatrix(&a3x3fRotMat[0][0], nVecNum, a3fRot[0], a3fRot[1], a3fRot[2],
		    a3fRot[3]);
  vMulMat3DVec3D(a3x3fRotMat, &m_pfVector[3*nVecNum], pfVec);
  return (0);
}

int Cgoniometer::nGetRotVector(Cstring& sName,float *pfVec) {
    // Find the goniometer axis.
    int nNum = nGetNum(sName);
    float a3fRotVec[3];

    if (nNum<0)
      {
        printf("ERROR:  Could not find rotation vector '%s'\n",sName.string());
        return -1;
      }
    if (nGetRotVector(nNum,&a3fRotVec[0])) {
        printf("ERROR:  Could not find rotation vector '%s'\n",sName.string());
        return -1;
    };
    vCopyVec3D(&a3fRotVec[0],pfVec);
    return 0;
};

Cstring
Cgoniometer::sGetName(const int nNum)
{
  if ( (0 <= nNum) && (nNum < m_nNumValues) )
    return (m_psName[nNum]);
  else
    return ("Illegal axis");
}

int
Cgoniometer::nGetNum(const Cstring sName)
{
    int i;
    Cstring sName1;
    Cstring sName2;
    sName1 = sName;
    sName1.upcase();
    for (i = 0; i < m_nNumValues; i++)
    {
        sName2 = m_psName[i];
        sName2.upcase();
        if (sName1 == sName2)
            return (i);
    }
    return (-1);
}


int  
Cgoniometer::nDiff(Cgoniometer& oGoniometerAdd,Cgoniometer& oGoniometerSubtract) {
    int nx;

    for (nx=0;nx<m_nNumValues;nx++) {
        m_pfDatumValue[nx] += oGoniometerAdd.m_pfDatumValue[nx] - oGoniometerSubtract.m_pfDatumValue[nx];
    };
    return 0;
};



int
Cgoniometer::nLoadOffsets(Cimage_header* poHeader) {
    int nStat;
    Cstring sTemp;
    int nx,ny;

    nStat = 0;
    
    // Get goniometer offsets.  
    
    sTemp = m_sRelPrefix;
    sTemp += D_K_GonioOffsetTableEntries;
    memset(&m_anGonioOffsetTableEntries[0],0,sizeof(int)*(nMaxGonioRotationOffsets+1));
    poHeader->nGetValue(
        sTemp,nMaxGonioRotationOffsets,
        &m_anGonioOffsetTableEntries[0]);
    if (m_anGonioOffsetTableEntries[0]) {
        for (nx=0,ny=0;m_anGonioOffsetTableEntries[nx];nx++)
            ny += m_anGonioOffsetTableEntries[nx];
        

        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableRotAxis;
        nStat = poHeader->nGetValue(
            sTemp,nx,
            &m_anGonioOffsetTableRotAxis[0]);

        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableOffset;
        nStat += poHeader->nGetValue(
            sTemp,ny,
            &m_afGonioOffsetTableOffset[0]);

        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableGonio;
        nStat += poHeader->nGetValue(
            sTemp,ny,
            &m_afGonioOffsetTableGonio[0]
            );
    } else {
        nStat = 0;
        m_anGonioOffsetTableEntries[0] = 0;
    };
    return nStat;
};


int
Cgoniometer::nSaveOffsets(Cimage_header* poHeader) {
    int nStat;
    Cstring sTemp;
    int nx,ny;


    nStat = 0;

    // Save information on goniometer offsets.
    
    if (m_anGonioOffsetTableEntries[0]) {
        for (nx=0,ny=0;m_anGonioOffsetTableEntries[nx];nx++) 
            ny += m_anGonioOffsetTableEntries[nx];
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableEntries;
        poHeader->nReplaceValue(
            sTemp,nx,
            &m_anGonioOffsetTableEntries[0]);
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableRotAxis;
        poHeader->nReplaceValue(
            sTemp,nx,
            &m_anGonioOffsetTableRotAxis[0]);
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableOffset;
        poHeader->nReplaceValue(
            sTemp,ny,
            &m_afGonioOffsetTableOffset[0]);
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableGonio;
        poHeader->nReplaceValue(
            sTemp,ny,
            &m_afGonioOffsetTableGonio[0]);
    } else {
        // Make sure these fields don't exist in the header!
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableEntries;
        poHeader->nDelete(sTemp);
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableGonio;
        poHeader->nDelete(sTemp);
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableOffset;
        poHeader->nDelete(sTemp);
        sTemp = m_sRelPrefix;
        sTemp += D_K_GonioOffsetTableRotAxis;
        poHeader->nDelete(sTemp);
        
    };

    return nStat;    
};



int 
Cgoniometer::nFindGonioOffsets(float* pfDatumValue,int nRotVector,int& nSet,int& nPoint) {
    int nx;

    if (pfDatumValue == NULL)
        pfDatumValue = m_pfDatumValue;

    nSet = 0;
    nPoint = 0;
    while ((m_anGonioOffsetTableEntries[nSet]) && (nSet < nMaxGonioRotationOffsets)) {
        for (nx = 0; nx < m_anGonioOffsetTableEntries[nSet]; nx++) {
            if (ABS(m_afGonioOffsetTableGonio[nPoint + nx] - pfDatumValue[nx])>0.001)
                break;
        };
        if ((nx == m_anGonioOffsetTableEntries[nSet]) && (nRotVector + 1 == m_anGonioOffsetTableRotAxis[nSet])) {
            return 0;
        };
        nPoint += m_anGonioOffsetTableEntries[nSet];        
        nSet++;
    };
    return 1;
};

int 
Cgoniometer::nApplyAllGonioOffsets(float* pfDatumValue,int nRotVector) {
    int nSet = 0;
    int nPoint = 0;
    int nx;

    if (pfDatumValue == NULL)
        pfDatumValue = m_pfDatumValue;

    // Lookup the goniometer constants.
    if (!nFindGonioOffsets(pfDatumValue,nRotVector,nSet,nPoint))  {
        for (nx = 0; nx < m_anGonioOffsetTableEntries[nSet]; nx++)
            pfDatumValue[nx] += m_afGonioOffsetTableOffset[nPoint + nx];
        return 0;
    };
    return 1;    
};

int 
Cgoniometer::nSaveGonioOffsets(float* pfDatumValue,int nRotVector,int nNumDatum,double* pfRotOffsets) {
    int nx;
    int nSet,nPoint;

    if (pfDatumValue == NULL)
        pfDatumValue = m_pfDatumValue;

    if (!nFindGonioOffsets(pfDatumValue,nRotVector,nSet,nPoint))  {
        for (nx = 0; nx < m_anGonioOffsetTableEntries[nSet]; nx++) {
            if (nx != nRotVector)
                m_afGonioOffsetTableOffset[nPoint + nx] = pfRotOffsets[nx];
        return 0;
        };
    };
    if (nSet < nMaxGonioRotationOffsets) {
        m_anGonioOffsetTableEntries[nSet] = nNumDatum;
        m_anGonioOffsetTableRotAxis[nSet] = nRotVector + 1;
        for (nx = 0; nx < m_anGonioOffsetTableEntries[nSet]; nx++) {
            if  (nx == nRotVector)
                m_afGonioOffsetTableOffset[nPoint + nx] = 0.0;
            else
                m_afGonioOffsetTableOffset[nPoint + nx] = pfRotOffsets[nx];
            m_afGonioOffsetTableGonio[nPoint + nx] = pfDatumValue[nx];
        };
    } else
        return 1;
    return 0;
};

int 
Cgoniometer::nGetGonioRotOffset(float* pfDatumValue,int nRotVector,double& fRotOffset) {
    int nSet,nPoint;

    if (pfDatumValue == NULL)
        pfDatumValue = m_pfDatumValue;

    if (!nFindGonioOffsets(pfDatumValue,nRotVector,nSet,nPoint))  {
        fRotOffset = m_afGonioOffsetTableOffset[nPoint + nRotVector];
        return 0;
    };
    fRotOffset = 0.0;
    return 1;
};

int Cgoniometer::nSaveGonioRotOffset(float* pfDatumValue,int nRotVector,int nNumDatum,double fRotOffset) {
    int nx;
    int nSet,nPoint;

    if (pfDatumValue == NULL)
        pfDatumValue = m_pfDatumValue;

    if (!nFindGonioOffsets(pfDatumValue,nRotVector,nSet,nPoint))  {
        m_afGonioOffsetTableOffset[nPoint + nRotVector] = fRotOffset;
        return 0;
    };
    if (nSet < nMaxGonioRotationOffsets) {
        m_anGonioOffsetTableEntries[nSet] = nNumDatum;
        m_anGonioOffsetTableRotAxis[nSet] = nRotVector + 1;
        for (nx = 0; nx < m_anGonioOffsetTableEntries[nSet]; nx++) {
            if  (nx == nRotVector)
                m_afGonioOffsetTableOffset[nPoint + nx] = fRotOffset;
            else
                m_afGonioOffsetTableOffset[nPoint + nx] = 0.0;
            m_afGonioOffsetTableGonio[nPoint + nx] = pfDatumValue[nx];
        };
    } else
        return 1;
    return 0;
};

int  Cgoniometer::nReconcileHeaderWithImageGoniometer(Cgoniometer* poHeaderGoniometer) {
    if (m_nNumValues == poHeaderGoniometer->m_nNumValues) {
        memcpy(m_pfVector,poHeaderGoniometer->m_pfVector,sizeof(float)*3*m_nNumValues);
        return 0;
    } else
        return 1;
};
