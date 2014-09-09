#ifndef DT_CSCANMS_H
#define DT_CSCANMS_H
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
// CscanMS.h        Initial author: J.W. Pflugrath           06-Mar-1996
//    This file is the header file for class Cscan for Microsoft systems
//
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

#include "Crotation.h"

//+Definitions and constants


// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eScan_states {
  eScan_unknown_state,
  eScan_ready_to_start_state,
  eScan_not_ready_to_start_state,
  eScan_active_state,
  eScan_finished_state,
  eScan_available_state
};

enum eScan_modes
{
    eScanMode_Unknown,			//An unknown mode
    eScanMode_StillClosed,		//For multiple DARK images
    eScanMode_StillOpen,	        //For Multiple STILLS
    eScanMode_ScanClosed,		//For Testing Purposes only
    eScanMode_ScanOpen			//For X-ray Diffraction Data 
};


//+Code begin

class Cscan {

public:

eScan_states eThe_State;
eScan_modes  eThe_Mode;
Cstring       sScan_key;

Cstring       sFileTemplate;     // A filename template for images
int          nFileStartSeqNum;  // The starting sequence number for the
                                //   first image in the scan
int          nFileStepNum;      // The sequence number step value
int          nFileSeqNum;       // The sequence number for the image to access

int          nNumImages;        // The number of images in the scan

int          nNumImagesAlloc;   // Number of images allocated in poTheImages
int          nNumImagesAvail;   // Number of images available in poTheImages

Crotation    *poRotation;      // Pointer to a Crotation object that describes
                               //  the scan
int          nNumDetectors;    // The number of detectors
int          nNumCounters;     // The number of counters
  
  float        fWavelength;

  int         nNumDetDatum;
  float       *pfDetDatum;


  int         nNumCrysDatum;
  float       *pfCrysDatum;

// Some static keywords used in image headers

  static Cstring  m_ssScanRotPrefix;
  static Cstring  m_ssScanTemplate;
  static Cstring  m_ssScanKey;
  static Cstring  m_ssScanSeqInfo;
  static Cstring  m_ssScanWavelength;
  static Cstring  m_ssScanCrysDatum;
  static Cstring  m_ssScanDetDatum;
  static Cstring  m_ssScanMode;
  static Cstring  m_ssScanModeStillO;
  static Cstring  m_ssScanModeStillC;
  static Cstring  m_ssScanModeScanO;
  static Cstring  m_ssScanModeScanC;
//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cscan ();                       // Construct an empty scan object
                                // Construct from disk files:
Cscan (const Cstring sTemplate, const int nSeqNum, const int nStepNum);

//Cscan (Cscan &oOther);         // Copy constructor
//Cscan& operator=(const Cscan& oOther);

~Cscan ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

int nInitValues(void);
int nList(void);
int nSetRotation(const Crotation& oRotation);
int nSetNonunf(const Cnonunf& oNonunf);
int nSetDatum(const int nCrys, float *pfCrys, const int nDet, float *pfDet);
int nGetDatum(int *pnCrys, float *pfCrys, int *pnDet, float *pfDet);

int nGetImageName(Cstring* sName);
int nGetNumImages(void);

// Inline functions:

inline eGetMode(void)
  {
    return (eThe_Mode);
  }
inline vSetMode(eScan_modes eModeIn = eScanMode_ScanOpen)
  { eThe_Mode = eModeIn; }
inline bool bIsAvailable(void)
  {  return (eThe_State == eScan_available_state); }
inline void vInitSeqNum(void) { nFileSeqNum = nFileStartSeqNum; }
inline void vNextSeqNum(void) { nFileSeqNum = nFileSeqNum + nFileStepNum; }
inline void vPrevSeqNum(void) { nFileSeqNum = nFileSeqNum - nFileStepNum; }
inline void vSetSeqNum(const int nImageNum)
  { nFileSeqNum = nFileStartSeqNum + (nImageNum * nFileStepNum); }
inline int nGetSeqNum(void) { return (nFileSeqNum); }
inline int nGetSeqNum(const int nImageNum) 
  { return (nFileStartSeqNum + (nImageNum * nFileStepNum)); }

inline int nGetSeqInc(void)
  { return (nFileStepNum); }
inline int vSetSeqInc(const int nInc)
  { nFileStepNum = nInc; }
inline int vSetSeqStart(const int nStart)
  { nFileStartSeqNum = nStart; }

inline Cstring& sGetTemplate(void)
  { return (sFileTemplate); }
inline void vSetTemplate(const Cstring &sTemplateIn)
  { sFileTemplate = sTemplateIn; }

inline void vSetWavelength(const float fWave)
  { fWavelength = fWave; }
inline float fGetWavelength(void)
  { return (fWavelength); }

 private:

};  // end of class CscanMS

#endif   // DT_CSCANMS_H

