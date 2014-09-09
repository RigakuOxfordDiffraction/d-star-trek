#ifndef DT_CSCAN_H
#define DT_CSCAN_H
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
// Cscan.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Cscan
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
#include "Cstring.h"
#include "Cimage.h"
#include "Crotation.h"
#include "Cspatial.h"
#include "Cnonunf.h"

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

typedef void (*vProgessCallbackProc) (void *,           // Pointer to object 
		                      const int nNum);  // progress number

//+Code begin

class DTREK_EXPORT Cscan {

public:

//bool         m_bDezingerImage;
//int          m_nDetBinMode;
Cstring      m_sDetectorOptions;

eScan_states m_eThe_State;
eScan_modes  m_eThe_Mode;
Cstring      m_sScan_key;

Cstring      m_sFileTemplate;    // A filename template for images
int          m_nFileStartSeqNum; // The starting sequence number for the
                                 //   first image in the scan
int          m_nFileStepNum;     // The sequence number step value
int          m_nFileSeqNum;      // The sequence number for the image to access
int          m_nNumImages;       // The number of images in the scan

int          m_nNumImagesAlloc;  // Number of images allocated in poTheImages
int          m_nNumImagesAvail;  // Number of images available in poTheImages

Cimage      *m_poTheImages;      // Probably dangerous to hold more than a few
                                 //   images in memory.

int          m_nNumDarkImagesAlloc;  // Number of dark images allocated in
int          m_nNumDarkImagesAvail;  // Number of dark images available
Cimage      *m_poTheDarkImages;      // Pointer to the dark images

int          m_nNumDCoffsetImagesAlloc; // Number of DCoffset images allocated
int          m_nNumDCoffsetImagesAvail; // Number of DCoffset images available
Cimage      *m_poTheDCoffsetImages;     // Pointer to DCoffset images

Crotation   *m_poRotation;       // Pointer to a Crotation object that describes
                                 //  the scan
Cspatial    *m_poSpatial;        // Pointer to a spatial distortion object
Cnonunf     *m_poNonunf;         // Pointer to a nonunformity object
 bool        m_bNewNonunf;

int          m_nNumDetectors;    // The number of detectors
int          m_nNumCounters;     // The number of counters

float        m_fWavelength;
int          m_nWavelengthOption;
int          m_nWavelengthOptimize;
int          m_nNumDetDatum;
float       *m_pfDetDatum;

int          m_nNumCrysDatum;
float       *m_pfCrysDatum;

// Some static keywords used in image headers

  static Cstring  ms_sScanPrefix;
  static Cstring  ms_sScanTemplate;
  static Cstring  ms_sScanKey;
  static Cstring  ms_sScanSeqInfo;
  static Cstring  ms_sScanWavelength;
  static Cstring  ms_sScanWavelengthOpts;
  static Cstring  ms_sScanCrysDatum;
  static Cstring  ms_sScanDetDatum;
  static Cstring  ms_sScanMode;
  static Cstring  ms_sScanModeUnknown;
  static Cstring  ms_sScanModeStillO;
  static Cstring  ms_sScanModeStillC;
  static Cstring  ms_sScanModeScanO;
  static Cstring  ms_sScanModeScanC;
  static Cstring  ms_sScanDezinger;
  static Cstring  ms_sScanDetBinMode;
  static Cstring  ms_sScanDetectorOptions;
//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cscan ();                       // Construct an empty scan object
                                // Construct from disk files:
Cscan (const Cstring sTemplate, const int nSeqNum, const int nStepNum);
Cscan (Cimage_header& oHeader, const Cstring& sPre="");   // Construct scan from an image header

//Cscan (Cscan &oOther);         // Copy constructor
//Cscan& operator=(const Cscan& oOther);

~Cscan ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

int nInitValues(void);
int nList(void);
int nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre="");
void vDeleteFromHeader(Cimage_header *poHeader, const Cstring &sPre="");
int nSetRotation(const Crotation& oRotation);
int nSetNonunf(const Cnonunf& oNonunf);
int nSetDatum(const int nCrys, float *pfCrys, const int nDet, float *pfDet);
int nGetDatum(int *pnCrys, float *pfCrys, int *pnDet, float *pfDet);

int nGetImage(Cimage* poImage);
Cstring sGetImageName(const int nSeqNum = -999999);
int nGetImageName(Cstring* sName);
int nOverlay(Cimage* poImage, const int nNum, 
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg) = NULL,
	     void *pObj=NULL);
int nUnderlay(Cimage* poImage, const int nNum,
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg) = NULL,
	     void *pObj=NULL);
int nAdd(Cimage* poImage, const int nNum,
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg) = NULL,
	     void *pObj=NULL);
int nTile(Cimage* poImage, const int nNum,
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg) = NULL,
	     void *pObj=NULL);
int nAvgSD(Cimage* poImage, const int nNum,
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg) = NULL,
	     void *pObj=NULL);

int nGetNumImages(void);

float fCalcGetRotStart(const int nSeqNum);
float fCalcGetRotEnd(const int nSeqNum);

// Inline functions:

inline eScan_modes eGetMode(void) { return (m_eThe_Mode);  }
inline void vSetMode(eScan_modes eModeIn = eScanMode_ScanOpen)
  { m_eThe_Mode = eModeIn; }
inline bool bIsAvailable(void)
  {  return (eScan_available_state == m_eThe_State); }
inline void vInitSeqNum(void) { m_nFileSeqNum = m_nFileStartSeqNum; }
inline void vNextSeqNum(void) { m_nFileSeqNum = m_nFileSeqNum + m_nFileStepNum; }
inline void vPrevSeqNum(void) { m_nFileSeqNum = m_nFileSeqNum - m_nFileStepNum; }
inline void vSetNumImgs(const int nNum) { m_nNumImages = nNum; }
inline void vSetSeqNum(const int nSeqNum)
  { m_nFileSeqNum = nSeqNum; }
inline void vSetImgNum(const int nImageNum)
  { m_nFileSeqNum = m_nFileStartSeqNum + (nImageNum * m_nFileStepNum); }
inline int nGetSeqNum(void) { return (m_nFileSeqNum); }
inline int nGetSeqNum(const int nImageNum) 
  { return (m_nFileStartSeqNum + (nImageNum * m_nFileStepNum)); }
       int nGetSeqNum(const Cstring& rsImageName);
inline int nGetSeqInc(void)
  { return (m_nFileStepNum); }
inline void vSetSeqInc(const int nInc)
  { m_nFileStepNum = nInc; }
inline void vSetSeqStart(const int nStart)
  { m_nFileStartSeqNum = nStart; m_nFileSeqNum = nStart; }

inline int nGetSeqNumImages(void) { return (m_nNumImages); };
inline int nGetNumWildcards(void) { 
    int nx = -1,nPlaces = 0; 
    while ( (nx=m_sFileTemplate.find('?', nx+1)) != -1) 
        nPlaces++; 
    return nPlaces;
};


inline Cstring& sGetTemplate(void)
  { return (m_sFileTemplate); }
inline void vSetTemplate(const Cstring &sTemplateIn)
  { m_sFileTemplate = sTemplateIn; }

inline void vSetWavelength(const float fWave)
  { m_fWavelength = fWave; }
inline float fGetWavelength(void)
  { return (m_fWavelength); }
inline void vSetWavelengthOpts(const int nOption=0, const int nOptimize=0)
  { m_nWavelengthOption = nOption; m_nWavelengthOptimize = nOptimize; }
inline void vGetWavelengthOpts(int *pnOption, int *pnOptimize)
  { *pnOption = m_nWavelengthOption; *pnOptimize = m_nWavelengthOptimize; }

inline Cstring sGetDetectorOptions(void)
  {  return m_sDetectorOptions; }
inline void vSetDetectorOptions(Cstring s)
  {  m_sDetectorOptions = s; return; }

//inline bool bGetDezinger(void)
//  {  return (m_bDezingerImage); }
//inline void vSetDezinger(const bool bValue)
//  {  m_bDezingerImage = bValue; return; }

//inline int nGetDetBinMode(void)
//  {  return (m_nDetBinMode); }
//inline void vSetDetBinMode(const int nValue)
//  {  m_nDetBinMode = nValue; return; }

 private:

};  // end of class Cscan

#endif   // DT_CSCAN_H

