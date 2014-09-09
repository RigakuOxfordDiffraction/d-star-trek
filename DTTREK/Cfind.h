#ifndef DT_CFIND_H
#define DT_CFIND_H
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
// Cfind.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Cfind
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
#include "Crotation.h"
#include "Cspatial.h"
#include "Cnonunf.h"
#include "Cscan.h"
#include "Creflnlist.h"
#include "C3Ddata.h"
#include "Cdetector.h"
#include "Csource.h"
#include "CImageWait.h"

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))

#include "CXprop.h"

#endif

//+Definitions and constants
// For safety on the enum have the unknown method first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eFind_methods
{
  eFind_unknown_method,
  eFind_2D_method,
  eFind_2D_merge_method
};

//+Code begin

class DTREK_EXPORT Cfind {

public:

  float       m_fSigma;         // Sigma above background threshold
  float       m_fMinimum;       // Minimum threshold
  float       m_fPeakMinimum;   // Minimum threshold (peak max mixel - mean bkg pixel)
  int         m_a4nCircle[4];   // Search Circle limits: center Px0, Px1, min max rad
  int         m_a4nRect[4];     // Search Rect limits (for entire image), min max Px0, min max Px1
  int         m_a2nSpotWindow[2];//C3Ddata spot window in Px0, Px1
  float       m_fNear;          // No neighbors nearer than this
  int         m_a2nPMTInfo[2];  // PMT info parameters: Num highest pixels, min threshold
  int         m_nPeakFilter;    // The number of pixels in 3x3 window > thresh.
  int         m_a2nSeqNum[2];   // Start and ending SeqNum in poScan to search
  int         m_nEveryNthImage; // Argument to specify only some images between seq nums
  int         m_nDumpRefln;     // Write every m_nDump.. refln shoebox to disk

  int         m_nBLoops;        // Num loop to get best local background & sigma
  int         m_nDisplay;       // Flag for sending display messages back to parent process
  bool        m_bIncludeSat;    // Include saturated reflections.
  bool        m_bIncludeBackground; // Include background pixel information in output reflection list.
  eFind_methods m_eMethod;      // 2D or 3D method

  C3Ddata     m_oShoebox;               // Used for nFindCentroid2D()
  Cscan      *m_poScan;                 // Pointer to scan object to search
  Creflnlist *m_poReflnlist;            // Pointer to reflnlist object to put peaks in.
  Cimage_header m_oHeader;              // Default header for various things
  Cgoniometer* m_poGoniometer;          // Stores the CRYSTAL_ goniometer information for the last input image.
  Crotation*   m_poRotation;            // Stores the ROTATION vector for the last input image.

  float       m_a2fDetXYSigmas[2];      // Set to legal values.  Tries to obtain from header somehow.
  int         m_nTotalImages;           // Total images processed.
  int         m_nTotalSpotsProcessed;
  int         m_nTotalSatSpotsProcessed;
  C3DdataDetectorProfile m_oDetectorAreas; // Detector area information.
  C3DdataDetectorProfile m_oDetectorAreasIn; // Detector area information (read in from header).

  bool        m_bPeakInfo; // whether to use the m_oDetectorAreas information

  bool        m_bRank;  // whether to do ranking-related operations

private:
    int       m_nDim0, m_nDim1; // Dimensions of images finding spots in pixels
    int       m_nRadMin2, m_nRadMax2;
    int       m_nFirstReflnAtLast2DFind;  // First reflection added before the last 2D find.
    bool      m_bNewScan;
    bool      m_bNewNonunf;
    bool      m_bNewReflnlist;
    bool      m_bList;
    float     m_fResolutionMin;         // Min resolution limit in A (99999 is min)
    float     m_fResolutionMax;         // Max resolution limit in A (0.1 A is max)
    float     m_fDefaultRockingWidth;   // If we happen to get a rocking width of 0.0 in the input image.
    
    C3Ddata  *m_po3Ddata;
    Cstring m_sOptionString;


public:

    static Cstring ms_sMergeHead;
    static Cstring ms_sMergeNext;
    static Cstring ms_sIndex;
    static Cstring ms_sOnFirstLastImage;
    static Cstring ms_sDtfindOptions;
    static Cstring ms_sDtfindSpotSize;
    static Cstring ms_sDtfindAvgSpotsPerImage;
    static Cstring ms_sDtfindFractionSaturated;
    static Cstring ms_sDtfindPhotonValues;
    static Cstring ms_sDtfindPhotonIoverSig;  

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
    CXprop *m_poXprop;
    Cstring *m_psXpropUpdate;
#endif



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cfind ();                       // Construct an empty find object

Cfind (Cscan *poScan, Creflnlist *poReflnlist);
Cfind (Cimage_header& oHeader, Creflnlist *poReflnlist, Cnonunf *poNonunf=NULL, int nNumberOfImages=-1);

~Cfind ();

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

int nInitValues(void);
int nList(void);
int nGetPeaks(void);

int  nBackgroundAvgSD(Cimage & oImage, const int nPxArea[4],
                      const float fExcludeMin,
                      const float fExcludeMax,
                      float *fAvg, float *fSD);

int nMarkGonioNumbers() { m_nFirstReflnAtLast2DFind = m_poReflnlist->nGetNumReflns(); return 0; };
int nAddGonioNumbers(Cimage& oImage);

int nFind2D(Cimage& oImage, const float fRot, Cdetector *poDetectorIn=NULL,
            Csource *poSourceIn = NULL,bool bAddMergeInformation = false);

void vSetResolution(const float fResoMin, const float fResoMax);

int nMerge2DReflections(Creflnlist* poNewList, bool bRejectOverlappedDifferentHKL=false);
int nMerge2DReflectionsWithTwins(Creflnlist* poNewList, bool bRejectOverlappedDifferentHKL=false);

int nMerge2D_OverlapReject(Creflnlist& oReflnlist,itr<int>& anRefsIn,Crefln& oNewRef,Creflnlist& oNewList);

#define DTREK_DTFIND_FC2D_LOOK_FOR_WEAK_SPOTS		            0x0001
// Reporting weak spots implies looking for them, so DTREK_DTFIND_FC2D_LOOK_FOR_WEAK_SPOTS is a bitwise part of DTREK_DTFIND_FC2D_REPORT_WEAK_SPOTS
#define DTREK_DTFIND_FC2D_REPORT_WEAK_SPOTS		                0x0003
// Calculate peak info even if no spot is found. Useful for background evaluation.

#define DTREK_DTFIND_FC2D_REPORT_BKG_INFO_DISREGARD_ERROR       0x0004                                  

#define DTREK_DTFIND_FC2D_GET_EXTRA_INFO_SAT_COUNT              0x0008
#define DTREK_DTFIND_FC2D_GET_EXTRA_INFO_PEAK_MAX               0x0010
#define DTREK_DTFIND_FC2D_GET_EXTRA_INFO_ELLIPSE_COUNT          0x0020

#define DTREK_DTFIND_FC2D_GET_PEAK_AREA_INFO                    0x0040

#define DTREK_DTFIND_FC2D_GET_EXTRA_INFO   (DTREK_DTFIND_FC2D_GET_EXTRA_INFO_SAT_COUNT |\
                                            DTREK_DTFIND_FC2D_GET_EXTRA_INFO_PEAK_MAX |\
                                            DTREK_DTFIND_FC2D_GET_EXTRA_INFO_ELLIPSE_COUNT)   

int nFindCentroid2D(Cimage& oImage,Crefln* poRefln, Cdetector* poDetector = NULL, DTREK_WORD wCtrl=0U);

int nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre="");
inline void vSetOptionString(Cstring& rsOptionString)
     {m_sOptionString = rsOptionString; }

int  nExpandReflnlist(Creflnlist *poReflnlistIn,bool bExpandExtra,bool bAddEllipsoidInformation);
int  nWrite(Creflnlist& oList,Cstring& sName);


private:


bool bPixelInLimits(const int nPx0, const int nPx1);

int nFind2D(void);

public:
int nDoSearch(float  fResolutionAngstr_LowestTheta,
              float  fResolutionAngstr_HighestTheta,
              eFind_methods eMethod,
	          float fIntensityOverSigmaRatio,
              int   nFirstImage,
              int   nLastImage,
              int   nStep);

int nWriteReflectionList(Cstring& strReflectionFileName);
void vHeaderToDetectorAreas();
};  // end of class Cfind

#endif   // DT_CFIND_H

