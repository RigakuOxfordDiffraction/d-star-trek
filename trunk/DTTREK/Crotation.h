#ifndef DT_CROTATION_H
#define DT_CROTATION_H
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
// Crotation.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Crotation
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

#include "Cstring.h"
#include "Cimage_header.h"
#include "dtrekvec.h"

class Cscan;        // Forward declaration of class Cscan.
class Cgoniometer;  // Forward declaration of class Cgoniometer.

//+Definitions and constants

enum eRotation_states
{
  eRotation_unknown_state,
  eRotation_available_state
};

//+Code begin

class DTREK_EXPORT Crotation {

public:

  static Cstring ms_sRotation;
  static Cstring ms_sRotVector;
  static Cstring ms_sRotAxisName;
  static Cstring ms_sRotLimits;
  DTREK_WIN_DLL_DATA_EXPORT static int     ms_nVerbose;         

protected:
  eRotation_states m_eThe_State;
  float   m_fStart;            // Start angle of rotation in degrees relative to
  float   m_fEnd;              // End angle of rotation in degrees relative to .
  float   m_fIncrement;        // Increment per image in degrees
  float   m_fTime;             // Time per image in seconds
  int     m_nNumOscillations;  // Number of oscillations per image
  int     m_a2nDarkImages[2];  // Number of dark images at start,
                               //    dark image update interval (in images)
  float   m_fDarkChangeLimit;  // Dark image change limit as fraction;
  int     m_a2nDCoffsetImages[2];// Number of DCoffset images at start,
                               //    DCoffset update interval (in images)
  float   m_a3fLimits[3];      // Rotation limits. [0] = min [1] = max [2] = range. (not necc == [1] - [0])
  Cstring m_sName;             // Name of axis rotated
  float   m_a3fVector[3];      // Direction of axis as a vector
                               //   if/when fStart=0.0;
//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

public:
    

Crotation ();                          // Construct an empty rotation object

  // Construct rotation object from header :
Crotation (Cimage_header& oHeader, const Cstring& sPre="");

Crotation (const float fInitValues[]); // Construct rotation object from float

Crotation (const Crotation& oOther);   // Copy constructor

Crotation (Cscan *poScan, const int nSeqNum);

Crotation& operator=(const Crotation& oOther);

~Crotation ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
//

int nInitValues(void);
int nList(const int nVerbose=0,Cgoniometer* poGonio = NULL);
int nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre="");
void vDeleteFromHeader(Cimage_header *poHeader, const Cstring &sPre="");

void vSetVector(const float *pfValues);
void vSetName(const Cstring& sNameIn);
void vSetName(const char *pcNameIn);
inline Cstring& sGetName(void)
  { return (m_sName); }

int  nSetRotRange(const float fStartIn, const float fEndIn);

inline void vSetRotStart(const float fStartIn) {m_fStart = fStartIn;}
inline void vSetRotEnd(const float fEndIn) {m_fEnd = fEndIn;}
inline void vSetIncrement(const float fIncIn) {m_fIncrement = fIncIn;}
inline float fGetMidValue(void) { return ((m_fStart + m_fEnd) * 0.5f); }
inline float fGetIncrement(void) { return (m_fIncrement); }
inline float fGetRotStart(void) { return (m_fStart); }
inline float fGetRotEnd(void) { return (m_fEnd); }
inline float fGetExposureTime(void) { return (m_fTime); }
inline void  vSetExposureTime(const float fTimeIn)
  { m_fTime = fTimeIn; }
inline float fGetRotMax() { return m_a3fLimits[1]; };
inline float fGetRotMin() { return m_a3fLimits[0]; };
inline float fGetRotRange() { return m_a3fLimits[2]; };
inline void  vSetRotMinMaxRange(float fRotMin,float fRotMax,float fRange) { m_a3fLimits[0] = fRotMin; m_a3fLimits[1] = fRotMax; m_a3fLimits[2] = fRange;};
void vGetVector(float *pfValues);
void vCalcGetRotMatrix(const float fAngle, float *pfMatrix);
inline bool bIsAvailable(void)
  { return (m_eThe_State == eRotation_available_state); }

inline void vSetNumOsc(const int nNum)
  { m_nNumOscillations = nNum; }
inline int nGetNumOsc(void)
  { return (m_nNumOscillations); }

inline void vSetDarkIntvl(const int nNum)
  { m_a2nDarkImages[1] = nNum; }
inline int nGetDarkIntvl(void)
  { return (m_a2nDarkImages[1]); }

inline void vSetDarkInit(const int nNum)
  { m_a2nDarkImages[0] = nNum; }
inline int nGetDarkInit(void)
  { return (m_a2nDarkImages[0]); }

 private:

};  // end of class Crotation

#endif   // DT_CROTATION_H

