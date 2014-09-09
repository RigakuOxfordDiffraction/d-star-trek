#ifndef DT_CROTATIONMS_H
#define DT_CROTATIONMS_H
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
// CrotationMS.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Crotation for Microsoft platforms.
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

//+Definitions and constants


enum eRotation_states {
  eRotation_unknown_state,
  eRotation_available_state
};

//+Code begin

class Crotation {

public:

  eRotation_states eThe_State;
  float  fStart;              // Start angle of rotation in degrees relative to
  float  fEnd;                // End angle of rotation in degrees relative to ..
  float  fIncrement;          // Increment per image in degrees
  float  fTime;               // Time per image in seconds
  int    nNumOscillations;    // Number of oscillations per image
  int    nDarkImages[2];      // Number of dark images at start,
                              //    dark image update interval (in images)
  float  fDarkChangeLimit;    // Dark image change limit as fraction;
  int    nDCoffsetImages[2];  // Number of DCoffset images at start,
                              //    DCoffset update interval (in images)
  Cstring sName;               // Name of axis rotated
  float  fVector[3];          // Direction of axis as a vector 
                              //   if/when fStart=0.0;
  static Cstring m_ssRotation;
  static Cstring m_ssRotVector;
  static Cstring m_ssRotAxisName;

//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Crotation ();                          // Construct an empty rotation object

  // Construct rotation object from header :

Crotation (const float fInitValues[]); // Construct rotation object from float

Crotation (const Crotation& oOther);   // Copy constructor

Crotation& operator=(const Crotation& oOther);

~Crotation ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

int nInitValues(void);
int nList(const int nVerbose=0);

void vSetVector(const float fValues[3]);
void vSetName(const Cstring& sNameIn);
inline Cstring& sGetName(void)
  { return (sName); }

int  nSetRotRange(const float fStartIn, const float fEndIn);


inline void vSetRotStart(const float fStartIn) {fStart = fStartIn;}
inline void vSetRotEnd(const float fEndIn) {fEnd = fEndIn;}
inline void vSetIncrement(const float fIncIn) {fIncrement = fIncIn;}
inline float fGetMidValue(void) { return ((fStart + fEnd) * 0.5); }
inline float fGetIncrement(void) { return (fIncrement); }
inline float fGetRotStart(void) { return (fStart); }
inline float fGetRotEnd(void) { return (fEnd); }
inline float fGetExposureTime(void) { return (fTime); }
inline void  vSetExposureTime(const float fTimeIn) 
  { fTime = fTimeIn; }
void vGetVector(float *pfValues);

inline bool bIsAvailable(void) 
  { return (eThe_State == eRotation_available_state); }

inline void vSetNumOsc(const int nNum)
  { nNumOscillations = nNum; }
inline int nGetNumOsc(void)
  { return (nNumOscillations); }

 private:

};  // end of class Crotation

#endif   // DT_CROTATIONMS_H

