#ifndef DT_CTRANSFORM_H
#define DT_CTRANSFORM_H
//
// Copyright (c) 1998 Molecular Structure Corporation
//
// Ctransform.h        Initial author: J.W. Pflugrath        22-December-1998
//    This file is the header file for class Ctransform, a barebones
//    non-uniformity & spatial distortion correction class.  It is a stripped
//    down version of Ccalibrate.h
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
#include "Cimage.h"
#include "Cstring.h"

//+Code begin

class Ctransform {

public:

private:

  Cstring      m_sTransformName;

  Cimage      *m_poImgTransform;
  Cimage      *m_poImgDark;
  Cimage      *m_poImgBadInterp;
  int          m_a2nTOutDim[2];

  int          m_nNumBadInterpOut;
  long m_lNumTime;
  double m_dSumTime;
  double m_dSumTime2;

public:

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Ctransform (void);

~Ctransform ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

inline void vSetTransformName(Cstring &rsName) {m_sTransformName = rsName; }

int nTransform(Cimage *poImageIn, Cimage *poImageOut, Cimage *poImageDark);

int nAddNonunfInfo(Cimage_header &oHeader);
int nAddSpatialInfo(Cimage_header &oHeader);

int nGetTimingInfo(double *pdTime, double *pdSigma);

 private: // Some of the above function will be made private

Cstring sGetDetectorPrefix(Cimage_header &oHeader);

friend int nInterpolationFunction(int nOutDim0, int nOutDim1, int nBadValue,
                           int nNumBadInterpOut, Cimage* poImgBadInterp,
                           Cimage* poImageOut);

};  // end of class Ctransform

#endif   // DT_CTRANSFORM_H
