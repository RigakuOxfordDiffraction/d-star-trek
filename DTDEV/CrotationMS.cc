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
// CrotationMS.cc          Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Crotation which implements
//    the rotation encapsulation of d*TREK for Microsoft platforms
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

#include "Crotation.h"         // Class definition and prototypes

//+Definitions, constants, and initialization of static member variables

Cstring Crotation::m_ssRotation    = "ROTATION";
Cstring Crotation::m_ssRotVector   = "ROTATION_VECTOR";
Cstring Crotation::m_ssRotAxisName = "ROTATION_AXIS_NAME";

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Crotation::Crotation()
{
 (void) nInitValues();
}

Crotation::Crotation(const float fInitValues[])
{
  (void) nInitValues();
  fStart             = fInitValues[0];
  fEnd               = fInitValues[1];
  fIncrement         = fInitValues[2];
  fTime              = fInitValues[3];
  nNumOscillations   = (int) fInitValues[4];
  nDarkImages[0]     = (int) fInitValues[5];
  nDarkImages[1]     = (int) fInitValues[6];
  fDarkChangeLimit   = fInitValues[7];
  nDCoffsetImages[0] = (int) fInitValues[8];
  nDCoffsetImages[1] = (int) fInitValues[9];
  eThe_State         = eRotation_available_state;
}

Crotation::Crotation(const Crotation& oOther)
{
  eThe_State        = oOther.eThe_State;
  fStart            = oOther.fStart;
  fEnd              = oOther.fEnd;
  fIncrement        = oOther.fIncrement;
  fTime             = oOther.fTime;
  nNumOscillations  = oOther.nNumOscillations;
  nDarkImages[0]    = oOther.nDarkImages[0];
  nDarkImages[1]    = oOther.nDarkImages[1];
  fDarkChangeLimit  = oOther.fDarkChangeLimit;
  nDCoffsetImages[0]= oOther.nDCoffsetImages[0];
  nDCoffsetImages[1]= oOther.nDCoffsetImages[1];

  sName             = oOther.sName;
  fVector[0]        = oOther.fVector[0];
  fVector[1]        = oOther.fVector[1];
  fVector[2]        = oOther.fVector[2];
}

Crotation& Crotation::operator=(const Crotation& oOther)
{
  if (&oOther == this)    // These 2 statements prevent ...
    return *this;         //   ... assigning an object to itself

  eThe_State        = oOther.eThe_State;
  fStart            = oOther.fStart;
  fEnd              = oOther.fEnd;
  fIncrement        = oOther.fIncrement;
  fTime             = oOther.fTime;
  nNumOscillations  = oOther.nNumOscillations;
  nDarkImages[0]    = oOther.nDarkImages[0];
  nDarkImages[1]    = oOther.nDarkImages[1];
  fDarkChangeLimit  = oOther.fDarkChangeLimit;
  nDCoffsetImages[0]= oOther.nDCoffsetImages[0];
  nDCoffsetImages[1]= oOther.nDCoffsetImages[1];

  sName             = oOther.sName;
  fVector[0]        = oOther.fVector[0];
  fVector[1]        = oOther.fVector[1];
  fVector[2]        = oOther.fVector[2];
}

Crotation::~Crotation()
{
//  cout << "Crotation::destructor called!\n";
}

int Crotation::nInitValues(void)
{
//  cout << "Crotation::nInitValues called\n";
  eThe_State = eRotation_unknown_state;
  fStart     = 0.0;
  fEnd       = 0.0;
  fIncrement       = 0.1;
  fTime            = 0.0;
  nNumOscillations = 0;
  nDarkImages[0]   = 0;    // No dark images
  nDarkImages[1]   = 0;
  fDarkChangeLimit = 100.0; // Dark image must change by factor of 100!
  nDCoffsetImages[0]= 0;
  nDCoffsetImages[1]= 0;
  sName             = "Unknown X axis";
  fVector[0]        = 1.0;
  fVector[1]        = 0.0;
  fVector[2]        = 0.0;

  return (0);
}

int Crotation::nList(const int nVerbose)
{
  cout << "\nRotation list:"
       << "\n            Start: " << fStart
       << "\n              End: " << fEnd
       << "\n        Increment: " << fIncrement
       << "\n             Time: " << fTime << endl;
  if (1 > nVerbose)
    {
      cout << "     Oscillations: " << nNumOscillations
           << "\n      Axis vector: " << fVector[0] << ", " 
                                     << fVector[1] << ", " << fVector[2];
      cout << "\n        Axis name: " << sName << endl;
    }
  if (2 > nVerbose)
    {
      cout << "       Dark start: " << nDarkImages[0]
           << "\n      Dark update: " << nDarkImages[1]
           << "\nDark change limit: " << fDarkChangeLimit
           << "\n  DCoffset  start: " << nDCoffsetImages[0]
           << "\n  DCoffset update: " << nDCoffsetImages[1]
           << endl;
    }
  return (0);
}

void Crotation::vSetVector(const float fValues[])
{
  fVector[0] = fValues[0];
  fVector[1] = fValues[1];
  fVector[2] = fValues[2];
}

void Crotation::vSetName(const Cstring& sNameIn)
{
  sName = sNameIn;
}

void Crotation::vGetVector(float *pfValues)
{
//  (void) fNormVec3D(fVector);
  for (int i = 0; i < 3; i++) pfValues[i] = fVector[i];
}

int  Crotation::nSetRotRange(const float fStartIn, const float fEndIn)
{
  if (fEndIn < fStartIn) return (1);
  fStart = fStartIn;
  fEnd   = fEndIn;
  return (0);
};


