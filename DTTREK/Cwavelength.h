#ifndef DT_CWAVELENGTH_H
#define DT_CWAVELENGTH_H
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
// Cwavelength.h        Initial author: J.W. Pflugrath           01-May-1995
//    This file is the header file for class Cwavelength
//
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

//+Definitions and constants


// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eWavelength_types {
  eWavelength_unknown_type,
  eWavelength_single_type,
  eWavelength_multiple_type,
  eWavelength_spectrum_type
};

//+Code begin

class DTREK_EXPORT Cwavelength {

public:

eWavelength_types  m_eThe_Type;    // The type of wavelength
int                m_nAlloc;       // The number of wavelengths allocated
float             *m_pfWavelength; // The values of the wavelengths in Angstrom
Cstring            m_sInputUnits;  // The input units of the wavelengths

static Cstring  ms_sWavelength;
static Cstring  ms_sWavelengthKey;

private:
  static float     ms_fAngToEv;
  static float     ms_fEvToAng;
  static float     ms_afKnownKa1Ka2Values[];

//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments
 public:

Cwavelength ();                    // Construct an empty wavelength object

Cwavelength (const float fLambda); // Construct simple wavelength object

Cwavelength (Cimage_header& oHeader, const Cstring sPrefix = "", bool bA1A2Check = TRUE); 
                                      // Construct a wavelength object from
                                      // image header

Cwavelength (const int nNum, const float* pfLambda); // Construct multi-
                                                     // wavelength object

~Cwavelength ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

inline float fGetWavelength(const int nNum=0)
{
    return (nNum < m_nAlloc ? (float)*(m_pfWavelength + nNum) : 0.0f); 
}

int   nSetWavelength(const float fLambda);
int   nSetWavelengths(const int nNum, const float *pfLambda);
int   nList(const int nFlag=0);
int   nUpdateHeader(Cimage_header* poHeader, const Cstring& sPre);
int   nVerifyKa1Ka2Information();

private:

int nInitValues(void);

bool m_bA1A2Check;

};  // end of class Cwavelength

#endif   // DT_CWAVELENGTH_H

