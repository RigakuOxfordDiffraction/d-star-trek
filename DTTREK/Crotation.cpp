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
// Crotation.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Crotation which implements
//    the rotation encapsulation of d*TREK.
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
#include "Crotation.h"         // Class definition and prototypes
#include "Cscan.h"
#include "dtrekdefs.h"
//+Definitions, constants, and initialization of static member variables

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif



Cstring Crotation::ms_sRotation    = D_K_Rotation;
Cstring Crotation::ms_sRotVector   = D_K_RotVector;
Cstring Crotation::ms_sRotAxisName = D_K_RotAxisName;
Cstring Crotation::ms_sRotLimits   = D_K_RotLimits;
int     Crotation::ms_nVerbose     = 1;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Crotation::Crotation()
{
 (void) nInitValues();
}

Crotation::Crotation(Cimage_header& oHeader, const Cstring& sPre)
{
  Cstring sTemp;
  Cstring sPrefix;

  sPrefix = sPre;
  (void) nInitValues();

  if (!oHeader.bIsAvailable())
    {
      if (ms_nVerbose) {
          
          cout << "Crotation ERROR: image header is not valid. "
              << "  Cannot construct rotation object!\n";
      };
    }
  else
    {
      // Try to get required information from the header

      float fInitValues[10];
      int   nTemp;
      float fTemp;
      nTemp = oHeader.nGetValue(sPrefix + ms_sRotation, 10, fInitValues);
      if (0 == nTemp)
        {
          m_fStart               = fInitValues[0];
          m_fEnd                 = fInitValues[1];
          m_fIncrement           = fInitValues[2];
          m_fTime                = fInitValues[3];

	  // Bug for BM19 images
	  nTemp = oHeader.nGetValue("BM19_ROTATION_START", &m_fStart);
	  if (0 == nTemp)
	    {
	      nTemp = oHeader.nGetValue("BM19_ROTATION_END", &m_fEnd);
	      if (0 == nTemp) m_fIncrement = m_fEnd - m_fStart;
	      nTemp = oHeader.nGetValue("BM19_SCAN_EXPOSURE_TIME", &m_fTime);
	    }

	  // Bug for ID19 images
	  nTemp = oHeader.nGetValue("ID19_ROTATION_START", &m_fStart);
	  if (0 == nTemp)
	    {
	      nTemp = oHeader.nGetValue("ID19_ROTATION_END", &m_fEnd);
	      if (0 == nTemp) m_fIncrement = m_fEnd - m_fStart;
	      nTemp = oHeader.nGetValue("ID19_SCAN_EXPOSURE_TIME", &m_fTime);
	    }

          m_nNumOscillations     = (int) fInitValues[4];
          m_a2nDarkImages[0]     = (int) fInitValues[5];
          m_a2nDarkImages[1]     = (int) fInitValues[6];
          m_fDarkChangeLimit     = fInitValues[7];
          m_a2nDCoffsetImages[0] = (int) fInitValues[8];
          m_a2nDCoffsetImages[1] = (int) fInitValues[9];

          // The following are optional rotation keywords:

          nTemp = oHeader.nGetValue(sPrefix + ms_sRotVector, 3, fInitValues);
          if (0 == nTemp)
            {
              if (0.0 >= fLenVec3D(fInitValues))
                {
                  if (ms_nVerbose) {
                      cout << "Crotation ERROR: image header is not valid!"
                          << "\n  BOGUS value for " << sPrefix + ms_sRotVector
                          << "\n     changed to 1 0 0!\n" << endl;
                  };
                  fInitValues[0] = 1.0;
                  fInitValues[1] = 0.0;
                  fInitValues[2] = 0.0;
                }

              vSetVector(fInitValues);
            }
          nTemp = oHeader.nGetValue(sPrefix + ms_sRotAxisName, &sTemp);
          if (0 == nTemp)
            {
              vSetName(sTemp);
            }

          nTemp = oHeader.nGetValue(sPrefix + ms_sRotLimits,3,fInitValues);
          if (0 == nTemp)
              vSetRotMinMaxRange(min(fInitValues[0],fInitValues[1]),max(fInitValues[0],fInitValues[1]),fInitValues[2]);
          m_eThe_State = eRotation_available_state;
        }
      else if (0 == oHeader.nGetValue(sPrefix + "OSC_START", &fTemp))
        {
          // Modify to partially read BNL images.

          m_fStart = fTemp;
          m_fEnd   = fTemp;
          if (0 == oHeader.nGetValue(sPrefix + "OSC_RANGE", &fTemp))
            {
              m_fIncrement = fTemp;
              m_fEnd       = (float)((double)m_fStart + (double)fTemp);  // Large start, small inc?
            }
          m_eThe_State = eRotation_available_state;
        }
      else
        {
          if (ms_nVerbose) {
              cout << "Crotation ERROR: image header is not valid. "
                  << "  Cannot construct rotation object!\n";
          };
          m_eThe_State = eRotation_unknown_state;
        }
    }
}

Crotation::Crotation(const float fInitValues[])
{
  (void) nInitValues();
  m_fStart               = fInitValues[0];
  m_fEnd                 = fInitValues[1];
  m_fIncrement           = fInitValues[2];
  m_fTime                = fInitValues[3];
  m_nNumOscillations     = (int) fInitValues[4];
  m_a2nDarkImages[0]     = (int) fInitValues[5];
  m_a2nDarkImages[1]     = (int) fInitValues[6];
  m_fDarkChangeLimit     = fInitValues[7];
  m_a2nDCoffsetImages[0] = (int) fInitValues[8];
  m_a2nDCoffsetImages[1] = (int) fInitValues[9];
  m_eThe_State           = eRotation_available_state;
}

Crotation::Crotation(const Crotation& oOther)
{
  m_eThe_State           = oOther.m_eThe_State;
  m_fStart               = oOther.m_fStart;
  m_fEnd                 = oOther.m_fEnd;
  m_fIncrement           = oOther.m_fIncrement;
  m_fTime                = oOther.m_fTime;
  m_nNumOscillations     = oOther.m_nNumOscillations;
  m_a2nDarkImages[0]     = oOther.m_a2nDarkImages[0];
  m_a2nDarkImages[1]     = oOther.m_a2nDarkImages[1];
  m_fDarkChangeLimit     = oOther.m_fDarkChangeLimit;
  m_a2nDCoffsetImages[0] = oOther.m_a2nDCoffsetImages[0];
  m_a2nDCoffsetImages[1] = oOther.m_a2nDCoffsetImages[1];
  m_a3fLimits[0]         = oOther.m_a3fLimits[0];
  m_a3fLimits[1]         = oOther.m_a3fLimits[1];
  m_a3fLimits[2]         = oOther.m_a3fLimits[2];

  m_sName                = oOther.m_sName;
  m_a3fVector[0]         = oOther.m_a3fVector[0];
  m_a3fVector[1]         = oOther.m_a3fVector[1];
  m_a3fVector[2]         = oOther.m_a3fVector[2];
}

Crotation::Crotation(Cscan *poScan, const int nSeqNum)
{
  // Create a rotation object from a scan object and a sequence number
  // within the scan

  (void) nInitValues();

  // Copy the rotation object out of the scan,
  // then reset rotation start and end

  *this            = *(poScan->m_poRotation);
  m_fStart         = poScan->fCalcGetRotStart(nSeqNum);
  m_fEnd           = poScan->fCalcGetRotEnd(nSeqNum);
  m_eThe_State     = eRotation_available_state;
}

Crotation& Crotation::operator=(const Crotation& oOther)
{
  if (&oOther == this)    // These 2 statements prevent ...
    return *this;         //   ... assigning an object to itself

  m_eThe_State           = oOther.m_eThe_State;
  m_fStart               = oOther.m_fStart;
  m_fEnd                 = oOther.m_fEnd;
  m_fIncrement           = oOther.m_fIncrement;
  m_fTime                = oOther.m_fTime;
  m_nNumOscillations     = oOther.m_nNumOscillations;
  m_a2nDarkImages[0]     = oOther.m_a2nDarkImages[0];
  m_a2nDarkImages[1]     = oOther.m_a2nDarkImages[1];
  m_fDarkChangeLimit     = oOther.m_fDarkChangeLimit;
  m_a2nDCoffsetImages[0] = oOther.m_a2nDCoffsetImages[0];
  m_a2nDCoffsetImages[1] = oOther.m_a2nDCoffsetImages[1];
  m_a3fLimits[0]         = oOther.m_a3fLimits[0];
  m_a3fLimits[1]         = oOther.m_a3fLimits[1];
  m_a3fLimits[2]         = oOther.m_a3fLimits[2];

  m_sName                = oOther.m_sName;
  m_a3fVector[0]         = oOther.m_a3fVector[0];
  m_a3fVector[1]         = oOther.m_a3fVector[1];
  m_a3fVector[2]         = oOther.m_a3fVector[2];

  return (*this);
}

Crotation::~Crotation()
{
//  cout << "Crotation::destructor called!\n";
}

int Crotation::nInitValues(void)
{
//  cout << "Crotation::nInitValues called\n";

  m_eThe_State           = eRotation_unknown_state;
  m_fStart               = 0.0;
  m_fEnd                 = 0.0;
  m_fIncrement           = 0.0;
  m_fTime                = 0.0;
  m_nNumOscillations     = 0;
  m_a2nDarkImages[0]     = 0;    // No dark images
  m_a2nDarkImages[1]     = 0;
  m_fDarkChangeLimit     = 100.0; // Dark image must change by factor of 100!
  m_a2nDCoffsetImages[0] = 0;
  m_a2nDCoffsetImages[1] = 0;
  m_sName                = "Unknown X axis";
  m_a3fVector[0]         = 1.0;
  m_a3fVector[1]         = 0.0;
  m_a3fVector[2]         = 0.0;
  m_a3fLimits[0]         = 0.0;
  m_a3fLimits[1]         = 360.0;
  m_a3fLimits[2]         = 1000.0;


  return (0);
}

int Crotation::nList(const int nVerbose,Cgoniometer* poGonio)
{
  if (10 > nVerbose)
    printf("\nRotation list:"
           "\n            Start: %8.3f"
           "\n              End: %8.3f"
           "\n        Increment: %8.3f"
           "\n             Time: %8.3f\n",
           m_fStart, m_fEnd, m_fIncrement, m_fTime);
  if (1 > nVerbose)
    {
      cout << "     Oscillations: " << m_nNumOscillations
           << "\n      Axis vector: " << m_a3fVector[0] << ", "
                                      << m_a3fVector[1] << ", "
                                      << m_a3fVector[2];
      cout << "\n        Axis name: " << m_sName << endl;
    }
  if (2 > nVerbose)
    {
      cout << "       Dark start: " << m_a2nDarkImages[0]
           << "\n      Dark update: " << m_a2nDarkImages[1]
           << "\nDark change limit: " << m_fDarkChangeLimit
           << "\n  DCoffset  start: " << m_a2nDCoffsetImages[0]
           << "\n  DCoffset update: " << m_a2nDCoffsetImages[1]
           << "\n    Min Rot limit: " << m_a3fLimits[0]
           << "\n    Max Rot limit: " << m_a3fLimits[1]
           << "\n  Rot Range limit: " << m_a3fLimits[2]
           << endl;
    }
  
    if( 10 <= nVerbose && poGonio )
    {
        int nAxis;
        float a3fTemp[3];

        nAxis = poGonio->nGetNum(sGetName());
        if (nAxis >=0)
        {
            printf("Rotation axis:\n");
            if (poGonio->nGetRotVector(nAxis,&a3fTemp[0])) 
            return 0;

            printf("%13s (%4.0lf - %4.0lf)               (%7.3f, %7.3f, %7.3f)\n",
            m_sName.string(), 
            m_a3fLimits[0],m_a3fLimits[1],
            a3fTemp[0], a3fTemp[1], a3fTemp[2]);

            printf("=====================================================================\n");
        }
            
        if( nVerbose == 11 )
        {
            printf("Rotation range: %4.0lf - %4.0lf\n", m_fStart, m_fEnd);
            printf("=====================================================================\n");
        }
    }

    return 0;
}

void Crotation::vSetVector(const float *pfValues)
{
  m_a3fVector[0] = pfValues[0];
  m_a3fVector[1] = pfValues[1];
  m_a3fVector[2] = pfValues[2];
}

void Crotation::vSetName(const char *pcNameIn)
{
    m_sName = pcNameIn;
}

void Crotation::vSetName(const Cstring& sNameIn)
{
    m_sName = sNameIn;
}

void Crotation::vGetVector(float *pfValues)
{
    fNormVec3D(m_a3fVector);
  
    for(int i = 0; i < 3; i++) 
        pfValues[i] = m_a3fVector[i];
}

int  Crotation::nSetRotRange(const float fStartIn, const float fEndIn)
{
    if( fEndIn < fStartIn ) 
        return 1;

    m_fStart = fStartIn;
    m_fEnd   = fEndIn;

    return 0;
}

void Crotation::vCalcGetRotMatrix(const float fAngle, float *pfMatrix)
{
  // Calculate and return the rotation matrix specified by a rotation
  // by fAngle degrees around the rotation axis.  This is a convenience
  // routine.  Note that the member variable m_a3fVector should
  // specify the rotation axis vector AFTER all crystal goniometer
  // rotations have been applied, that is at DATUM.

  vConvRotVec3DMat3D(fAngle, m_a3fVector, pfMatrix);
}

int
Crotation::nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre)
{
  // Place this rotation object in the header

  float fInitValues[10];

  fInitValues[0] = m_fStart;
  fInitValues[1] = m_fEnd;
  fInitValues[2] = m_fIncrement;
  fInitValues[3] = m_fTime;
  fInitValues[4] = (float) m_nNumOscillations;
  fInitValues[5] = (float) m_a2nDarkImages[0];
  fInitValues[6] = (float) m_a2nDarkImages[1];
  fInitValues[7] = m_fDarkChangeLimit;
  fInitValues[8] = (float) m_a2nDCoffsetImages[0];
  fInitValues[9] = (float) m_a2nDCoffsetImages[1];

  (void) poHeader->nReplaceValue(sPre + ms_sRotation, 10, fInitValues, 3);
  (void) poHeader->nReplaceValue(sPre + ms_sRotVector, 3, m_a3fVector, 3);
  (void) poHeader->nReplaceValue(sPre + ms_sRotAxisName, m_sName);
  (void) poHeader->nReplaceValue(sPre + ms_sRotLimits,3,m_a3fLimits,3);

  return (0);
}

void
Crotation::vDeleteFromHeader(Cimage_header *poHeader, const Cstring &sPre)
{
  // Delete this rotation from the header

  (void) poHeader->nDelete(sPre + ms_sRotation);
  (void) poHeader->nDelete(sPre + ms_sRotVector);
  (void) poHeader->nDelete(sPre + ms_sRotAxisName);
  (void) poHeader->nDelete(sPre + ms_sRotLimits);
}
