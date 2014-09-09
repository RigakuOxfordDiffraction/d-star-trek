#ifndef DT_CNONUNF_H
#define DT_CNONUNF_H
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
// Cnonunf.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Cnonunf
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
#include "Cimage.h"
#include "C3Ddata.h"
#include "Cdetector.h"

//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eNonunf_states {
  eNonunf_unknown_state,
  eNonunf_none_state,
  eNonunf_simple_scale_state,
  eNonunf_simple_mask_state,
  eNonunf_dark_nonunf_state,
  eNonunf_dark_dcoffset_state
};

//+Forward class declarations

class Cdetector;   // Forward declaration of Cdetector class.
class CImagePool;
//+Code begin

class DTREK_EXPORT Cnonunf {

public:

eNonunf_states m_eThe_State;

Cstring        m_sNonunf_key;
int            m_nMethod;
float          m_fScaleFactor;
float          m_fMinFactorLimit;
float          m_fMaxFactorLimit;
float          m_fNumerator;
float          m_fDenominator;
int            m_a4nBadFlag[4];
float          m_fBadFlag;
float          m_fMinRawPixOKValue;

Cimage        *m_poNonunfImage;
Cimage        *m_poDarkImage;
Cimage        *m_poDCoffsetImage;

// Special for speed: a pointer to a member function.  The syntax of this
// is tricky, so watch out. 

float  (Cnonunf::*m_prfNonunfFactor)(const int Px0, const int Px1);

static int ms_nBadPixFlag;
static int ms_nSatPixFlag;

static Cstring ms_sNonunfStateUnknown;
static Cstring ms_sNonunfStateNone;
static Cstring ms_sNonunfStateSimpleScale;
static Cstring ms_sNonunfStateSimpleMask;
static Cstring ms_sNonunfStateDarkNonunf;
static Cstring ms_sNonunfStateDarkDCoffsetNonunf;
static Cstring ms_sNonunfType;
static Cstring ms_sNonunfDenominator;
static Cstring ms_sNonunfNumerator;
static Cstring ms_sNonunfFlag1;
static Cstring ms_sNonunfFlag2;
static Cstring ms_sNonunfFlag3;
static Cstring ms_sNonunfFlag4;
static Cstring ms_sNonunfInfo;
static Cstring ms_sNonunfKey;

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cnonunf ();                    // Construct an empty nonunf object
                               // Construct from disk files:

Cnonunf(const Cstring& sNonunfName); // Construct a simple non-uniformity 
                                    // correction object

Cnonunf(const Cstring& sNonunfName, const Cstring& sDarkName);

Cnonunf(const Cstring& sNonunfName, const Cstring& sDarkName,
        const Cstring& sDCoffsetName);

// Construct from information in the header of an image

Cnonunf(Cimage_header& oHeader, const Cstring& sPrefix = "");

~Cnonunf ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

int nInitValues(void);
int nInitValues(const Cstring& sNonunfName, Cimage_header *poHeader=NULL);
int nInitValues(const Cstring& sNonunfName, const Cstring& sDarkName);
int nInitValues(const Cstring& sNonunfName, const Cstring& sDarkName,
                const Cstring& sDCoffsetName);
int nInitValues(Cimage_header& oHeader, const Cstring& sPrefix = "");

int nList(void);

int nCorrectImage(Cimage *poImage);
int nCorrect3Ddata(C3Ddata *po3Ddata, float *pfPerLayerFactors=NULL);
int nCorrectPixel(const int nPx0, const int nPx1, Cimage *poImage,
                  float *pfValue, const bool bCheckBounds=TRUE);

int nScaleDark(const float fScaleFactor);
int nSubtractDCoffset(Cimage *poImage);
int nSubtractDarkAndDCoffset(Cimage *poImage);

int nSetFactorLimits(const float fMinLimit, const float fMaxLimit);

bool bPixelIsBad(const int nPx0, const int nPx1);
inline bool bIsAvailable(void) { return (m_eThe_State != eNonunf_unknown_state); }
inline float fGetBadFlag(void) { return (m_fBadFlag); }
inline void vSetBadFlag(const float fBad) { m_fBadFlag = fBad; }
inline void vSetMinRawPixOKValue(const float fMinPix) 
             { m_fMinRawPixOKValue = fMinPix; }

private:
CImagePool*       m_pImagePool;


float fSimpleNonunfFactor(const int nPx0, const int nPx1);
float fDarkDenomNonunfFactor(const int nPx0, const int nPx1);
float fDarkNumerNonunfFactor(const int nPx0, const int nPx1);
float fNoneNonunfFactor(const int nPx0, const int nPx1);
float fUnknownNonunfFactor(const int nPx0, const int nPx1);
float fGetMaskValue(const int nPx0, const int nPx1);
};  // end of class Cnonunf

#endif   // DT_CNONUNF_H
