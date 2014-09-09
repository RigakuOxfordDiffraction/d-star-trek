#ifndef DT_RAXIS_H
#define DT_RAXIS_H
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
// RAXIS.h        Initial author: J.W. Pflugrath           05-Jul-1995
//    This file is the header file for an RAXIS style image file
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
#include "Cstring.h"

//+Code begin


typedef struct _RAXIS_header {
  char  a10cDevice[10];
  char  a10cVersion[10];
  char  a20cCrystal[20];
  char  a12cCry_system[12];
  float a6fCell[6];
  char  a12cSpace[12];
  float fMosaic;
  char  a80cMemo[80];
  char  a84cReserve1[84];
  char  a12cDate[12];
  char  a20cOperatorname[20];
  char  a4cTarget[4];
  float fWave;                // X-ray wavelength
  char  a20cMonochro[20];     // type of monochromator
  float fMono_2;              // monochromator 2theta (deg)
  char  a20cCollimator[20];   // collimator size and type
  char  a4cFilter[4];         // filter type (Ni, etc.)
  float fCamera;              // Crystal to detector distance
  float fKv;
  float fMa;
  char  a12cFocus[12];        // focus info
  char  a80cXraymemo[80];     // X-ray memo
  INT4  nCylinder;            // something for Weissenberg photo?
  float fWeissenberg;         // something for Weissenberg photo?
  char  a56cReserve2[56];     // Reserved for future use
  char  a4cSpindle[4];        // Crystal mount axis closest to spindle axis
  char  a4cXray_axis[4];      // Crystal mount axis closest to beam axis
  float a3fPhi[3];            // Phi datum, phi start, phi end
  INT4  nOsc;                 // Number oscillations (frame number?)
  float fEx_time;             // Exposure time in minutes?
  float a2fXray[2];           // Direct beam position in pixels
  float a3fCircle[3];         // omega, chi, 2theta
  float fMu;                  // Spindle inclination angle?

#ifndef RAXIS_HEADER_BINARY_STRESS_INFO 
  char  a204cScanTemplate[204]; /* This space is now used for storing the scan
                                   templates information - tlh, 01 Feb 1999 */
#else
  long	posn;
  long	slit;	
  long	side; /* Iso or side Iso:0 Side:1 */	
  char	posz[64];
  char	psis[128];
#endif
  
  INT4  a2nPix_num[2];        // Number of fast, slow pixels
  float a2fPix_size[2];       // Size of fast, slow direction in mm
  INT4  a2nRecord[2];         // Record length in bytes, number of records
  INT4  nRead_start;          // For partial reads, 1st read line
  INT4  nIP_num;              // Which imaging plate 1, 2 ?
  float fRatio;               // Output ratio for high value pixels
  float a2fFading[2];         // Fading time to start of read, end of read
  char  a10cCpu[10];          // Type of computer "IRIS", "VAX", "SUN", etc
  char  a10cIp[10];           // Type of IP
  INT4  a3nDrxz[3];           // IP scanning codes??
  //+3-Apr-2002 modifications start
  float fPixShiftOdd;         // Pixel shift to odd lines
  float fIntRatioOdd;         // Intensity ratio to odd lines
  INT4  nMagicNum;            // Magic number to indicate next values are legit
  INT4  nNumGonAxes;          // Number of goniometer axes
  float a5x3fGonVecs[5][3];   // Goniometer axis vectors
  float a5fGonStart[5];       // Start angles for each of 5 axes
  float a5fGonEnd[5];         // End angles for each of 5 axes
  float a5fGonOffset[5];      // Offset values for each of 5 axes
  INT4  nScanAxisNum;         // Which axis is the scan axis?
  char  a40cAxesNames[40];    // Names of the axes (space or comma separated?)
  //  char  a8cReserve4[8];       // Reserved for future use.
  //  char  a180cReserve4[180];
  //- 3-Apr-2002 mods end
  //  char  a1024cReserve5[1024];
} RAXIS_header;

#define RAXIS_MAGIC_NUM 1

// Forward declaration of classes needed by the function prototypes below

class Cimage;
class Cimage_header;

DTREK_EXPORT extern RAXIS_header g_oRaxisHeader;

DTREK_EXPORT int nReadRAXISheader(const int* pnFile, const Cstring& rsFilename,
		     char *pcBuffer, Cimage_header *poHeader);

DTREK_EXPORT int nReadRAXISdata(const int *pnFile, unsigned short int *puiData,
		   Cimage_header *poHeader);

DTREK_EXPORT int nWriteRAXIS(Cimage &roImage, const Cstring &rsName,RAXIS_header* ptKnownHeader=NULL);

DTREK_EXPORT int nBuildRAXISHeader(RAXIS_header& rtHeader,Cimage& oImage);

#endif   // DT_RAXIS_H
