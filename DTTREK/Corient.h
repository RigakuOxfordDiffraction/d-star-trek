#ifndef DT_CORIENT_H
#define DT_CORIENT_H
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
// Corient.h        Initial author: J.W. Pflugrath           01-May-1995
//    This file is the header file for class Corient
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

#include "Dtrek.h"
#include "Cimage_header.h"
#include "Csource.h"
#include "Ccrystal.h"
#include "Cgoniometer.h"
#include "Crotation.h"
#include "dtrekvec.h"

//+Definitions and constants
//+Code begin

class DTREK_EXPORT Corient {

public:
private:

  Csource     *m_poSource;     // Pointer to source object
  Ccrystal    *m_poCrystal;    // Pointer to crystal object
  Cgoniometer *m_poCrysGonio;  // Pointer to crystal goniometer;
  Cstring      m_sRequest;     // Copy of requested orient command

  bool m_bNewSource;
  bool m_bNewCrystal;
  bool m_bNewCrysGonio;

  double m_fGonMatrix[3][3];        // Goniometer rotation matrix at datum

  double m_fCrysRotMatrix[3][3];    // Crystal rotation matrix when goniometer
                                   //   is at 0,0,0

  double m_fInvCrysRotMatrix[3][3]; // Inverse crystal rotation matrix 
                                   //   when goniometer is at 0,0,0

  double m_fGonCrysRotMatrix[3][3]; // Previous 2 matrices multiplied together

  double m_fCrysBMatrix[3][3];      // Crystal B matrix when goniometer
                                   //   is at 0,0,0
  double m_fInvCrysBMatrix[3][3];   // Inverse crystal B matrix when 
                                   //   goniometer is at 0,0,0
  double m_fTransInvCrysBMatrix[3][3]; // Transpose of inverse crystal B matrix 
                                   // orientation matrix when goniometer
				   //   is at 0,0,0


  double m_fGonAnglesUsed[3];       // Goniometer used angles for angles 
                                   //   in m_fVecAngles
  Cstring m_sVecNames[7];

  double m_fVecAngles[7][6];
  double m_fVecs[7][3];
  double m_fCrysVec[2][3];
  double m_fLabVec[2][3];
  int   m_nMode; 
  int   m_nNumSoln;
  double m_fAngleSoln[8][3];
//

public:

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Corient ();       // Construct an empty orient object
                  // In the future this should not be allowed to happen.

Corient (Cimage_header& oHeader);

~Corient ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

int nInitValues(void);
int nInitValues(Cimage_header& oHeader);
int nList(void);
int nSetup(void);
int nCalcGetVecAngles(const double *pfVec, double *pfAngles);
int nCalcAngles(const double *pfGonAngles = NULL);
int nSetOrient(const Cstring &sWishList);
int nCalcOrient(void);
int nListSoln(void);
int nPickSoln(const int nSolnIn=-1);
int nUpdateHeader(Cimage_header* poHeader);
int nParseVec(const Cstring &sVec, double *pfVec);
int nRotate(const Cstring& sRotate);
int nGetNumSoln() { return m_nNumSoln; }
void vGetAngleSoln(double fAngleSoln[8][3]);

 private:

};  // end of class Corient

#endif   // DT_CORIENT_H

