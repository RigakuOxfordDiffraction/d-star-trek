#ifndef DT_CBFDTREK_H
#define DT_CBFDTREK_H
//
// Copyright (c) 2007 Rigaku Americas Corp.
//
// cbfdtrek.h        Initial author: J.W. Pflugrath           22-Aug-2007
//    This file is the header file for a CBF style image file routines
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

//+Code begin

//  The following header definition comes from ....

/*
*/

/* This number is  written into the byte_order fields in the
   native byte order of the machine writing the file */

#include "cbf.h"

// Forward declaration of classes needed by the function prototypes below

class Cimage_header;

typedef struct _CBFDTREK_header {
  char ac2048text[15360];
  int   nSize1;
  int   nSize2;
  float fExpTime;
  float fPixSize0;
  float fPixSize1;
  float fWavelength;
  float fDet2theta;
  float fBeam0;
  float fBeam1;
  float fDetDist;
  float fOmega;
  float fChi;
  float fKappa;
  float fPhi;
  float fRotStart;
  float fRotInc;
  float fPolarz;
  float fSensorThickness;
} tagCBFDTREK_header;

int nReadCBFHeader(const int* pnFile, const Cstring& rsFilename,
                              char *pcBuffer, Cimage_header *poHeader);

int nReadCBFData(const Cstring &rsFilename, void *pvData,
                 Cimage_header *poHeader);

#endif   // DT_CBFDTREK_H
