#ifndef DT_CSOURCE_H
#define DT_CSOURCE_H
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
// Csource.h        Initial author: J.W. Pflugrath           01-May-1995
//    This file is the header file for class Csource
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
#include "Dtrek.h"
#include "Cwavelength.h"
#include "Cimage_header.h"
#include "dtrekvec.h"

//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eSource_states {
  eSource_unknown_state,
  eSource_notrefined_state,
  eSource_refined_state
};

//+Code begin

// This should probably be DERIVED from the Cwavelength class!
class DTREK_EXPORT Csource {

public:

eSource_states m_eThe_State;  // The state
Cstring        m_sType;

Cwavelength   *m_poWavelength;

float          m_a3x3fVector[3][3];
float          m_a2fSpectralDispersion[2];
float          m_a4fCrossfire[4];
float          m_a4fPolarz[4];
float          m_fIntensity;
float          m_a4fSize[4];

float          m_a2fRotation[2];

Cstring        m_sDescription;
Cstring        m_sSource_key;

// Some static keywords used in image headers

static Cstring  ms_sSourceVectors;
static Cstring  ms_sSourceValues;
static Cstring  ms_sSourceIntensity;
static Cstring  ms_sSourceSize;
static Cstring  ms_sSourcePolarz;
static Cstring  ms_sSourceCrossfire;
static Cstring  ms_sSourceKey;
static Cstring  ms_sSourceSpectralDispersion;
static Cstring  ms_sSourceRefineFlags;
static Cstring  ms_sSourcePrefix;

// float         *pfIntVSTime;
// bool           bScaleNeeded;
//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Csource ();                       // Construct an empty source object
Csource (Cimage_header& oHeader); // Construct a source object from header

~Csource ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:


int    nList(const int nFlag=0);
void   vCalcGetS0(float *pfS0,float *pfDeriv = NULL);
void   vCalcGetS0(double *pfS0,double *pfDeriv = NULL);
void   vCalcGetUnrotatedS0(double *pfS0);
void   vGetPolarz(float *pfPol);
int    nUpdateHeader(Cimage_header* poHeader);
void   vSetRotation(const float fRot0, const float fRot1);
int    nDiff(Csource& oSourceAdd,Csource& oSourceSubtract);

inline bool bIsAvailable(void) {return (m_eThe_State != eSource_unknown_state);}

inline float fGetWavelength(const int nNum=0)
{ return (m_poWavelength->fGetWavelength(nNum)); }

inline float fGetSpectralDispersion()
{ return ((float)sqrt((double)(m_a2fSpectralDispersion[0]*m_a2fSpectralDispersion[0] +
		m_a2fSpectralDispersion[1]*m_a2fSpectralDispersion[1]))); }

inline float fGetRotation(const int nNum=0)
{ return (m_a2fRotation[nNum]); }

//

private:

int nInitValues(void);
int nInitValues(Cimage_header& oHeader);

};  // end of class Csource

#endif   // DT_CSOURCE_H

