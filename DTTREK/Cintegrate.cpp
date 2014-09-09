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
// Cintegrate.cc        Initial author: J.W. Pflugrath           18-Sep-1995
//  This file contains the member functions of class Cintegrate which implements
//    the reflection integrating encapsulation of d*TREK.
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
//+Description
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#if (defined(SSI_PC) || (defined(NO_X_WINDOWS)))
#include <string.h>
#endif


#include "Dtrek.h"
#include "dtreksys.h"
#include "Cintegrate.h"         // Class definition and prototypes
#include "Crefine.h"
#include "Cfind.h"

#ifndef NO_X_WINDOWS
#ifndef SSI_PC

#include "CXprop.h"

#else

#include "CrclHelper.h"
// For Sleep()

#endif
#endif

#ifdef ANL_TIMER
#include "anlTimer.h"
#define ANLSTART(x) anl_start_timer(x);
#define ANLSTOP(x) anl_stop_timer(x);
#define ANLPRINT(x) anl_print_timer(x);
#else
#define ANLSTART(x) 
#define ANLSTOP(x) 
#define ANLPRINT(x)
#endif

using namespace std;

//+Definitions, constants, and initialization of static member variables


Cstring Cintegrate::ms_snDataFlag  = D_K_IntegratenDataFlag;
Cstring Cintegrate::ms_snErrorFlag = D_K_IntegratenErrorFlag;
Cstring Cintegrate::ms_sfDataAddr  = D_K_IntegratefDataAddr;
Cstring Cintegrate::ms_sfWidth0    = D_K_IntegratefWidth0;
Cstring Cintegrate::ms_sfWidth1    = D_K_IntegratefWidth1;


// Maybe the following should be enumerated... // Flag to mark ...
int Cintegrate::ms_nNewFlag        = -1;        // new reflections
int Cintegrate::ms_nActiveFlag     = -2;        // active reflections
int Cintegrate::ms_nDelFlag        = -7;        // reflections to delete


//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cintegrate::Cintegrate(Cimage_header *poHeaderIn, Creflnlist* poReflnlistOut)
{
  (void) nInitValues();
  m_poHeader         = poHeaderIn;                 // Should we copy this?
  m_poScan           = new Cscan(*m_poHeader);
  m_bNewScan         = TRUE;
  //+jwp 12-Jun-2001
  // Use a nonunf from poDetector later on.
  // m_poScan->m_poNonunf = new Cnonunf(*m_poHeader);
  //-jwp 12-Jun-2001
  m_anSeqNum[0]      = m_poScan->nGetSeqNum();
  m_anSeqNum[1]      = m_poScan->nGetSeqNum();
  m_poScan->vInitSeqNum();

  m_poReflnlist      = poReflnlistOut;

  //+2011-06-27 JWP
  // If this is a RAPID image, then make spotchase true
  // RX_DETECTOR_SIZE=460.0 256.0;
  int a2nRAPIDSize[2];
  a2nRAPIDSize[0] = 0;
  a2nRAPIDSize[1] = 0;
  (void) poHeaderIn->nGetValue("RX_DETECTOR_SIZE", 2, a2nRAPIDSize);
  if ( (460.0 == a2nRAPIDSize[0]) && (256.0 == a2nRAPIDSize[1]) )
    m_bSpotChase       = TRUE;
  //-2011-06-27 JWP

  // Look for DTINTEGRATE_OBLIQUE keyword in header and if found, set up
  // oblique incidence correction factors

  int i, nStat;
  m_nNumObliqueCorr = 0;
  nStat = poHeaderIn->nGetValue(D_K_IntegrateOblique, &m_nNumObliqueCorr);
  if ( (0 == nStat) && (0 < m_nNumObliqueCorr) )
    {
      if (NULL != m_pfObliqueCorr)
        delete [] m_pfObliqueCorr;
      if (NULL != m_pfObliqueCorrNorm)
        delete [] m_pfObliqueCorrNorm;
      if (NULL != m_pfObliqueCorrType)
        delete [] m_pfObliqueCorrType;

      m_pfObliqueCorr     = new float [m_nNumObliqueCorr+1];
      m_pfObliqueCorrNorm = new float [m_nNumObliqueCorr+1];
      m_pfObliqueCorrType = new float [m_nNumObliqueCorr+1];

      float *pfTemp;
      pfTemp = new float [2*m_nNumObliqueCorr + 1];
      nStat = poHeaderIn->nGetValue(D_K_IntegrateOblique, 2*m_nNumObliqueCorr+1,
                                    pfTemp);
      int nOff = 1;
      for (i = 0; i < m_nNumObliqueCorr; i++)
        {
          // Shift the values, since the first is really the number of factors

          // 0 means absorp, 1 means Trans

          m_pfObliqueCorrType[i] = pfTemp[nOff++];
          m_pfObliqueCorr[i]     = pfTemp[nOff++];

          if (   (0.0 != m_pfObliqueCorrType[i])
              && (1.0 != m_pfObliqueCorrType[i]))
            {
              cout << "WARNING, oblique incidence correction type is not 0 or 1!" << endl;
            }

          // Calculate 1 - exp(factor), when angle is normal to detector plane.
          // for Transmission types
          // Calculate -exp(factor), when angle is normal to detector plane.
          // for Absorption types

          m_pfObliqueCorrNorm[i] = m_pfObliqueCorrType[i]
                                   - (float)exp(m_pfObliqueCorr[i]);
          if (0.0 == m_pfObliqueCorrNorm[i])
            {
              cout << "ERROR: Oblique incidence correction factor SHOULD NOT be 0!\n";
              if (1 == m_nNumObliqueCorr)
                {
                  m_nNumObliqueCorr = 0;
                  cout << "       NOT USING any OBLIQUE INCIDENCE CORRECTION!\n\n";
                }
              else
                {
                   cout << "       Reset to normal to detector factor to 1.0!\n"
                        << "       ALL RESULTS WILL BE BOGUS!\n\n";
                   m_pfObliqueCorrNorm[i] = 1.0;
                }
            }
        }
      delete [] pfTemp;
    }
}

Cintegrate::Cintegrate(Cscan* poScanIn, Creflnlist* poReflnlistOut)
{
  (void) nInitValues();
  m_poScan       = poScanIn;
  m_bNewScan     = FALSE;
  m_poReflnlist  = poReflnlistOut;
  m_anSeqNum[0]  = m_poScan->nGetSeqNum();   // Get first sequence number
  m_anSeqNum[1]  = m_poScan->nGetSeqNum();   // ?? this is wrong
  m_poScan->vInitSeqNum();

}


Cintegrate::Cintegrate()
{
 (void) nInitValues();
}

Cintegrate::~Cintegrate()
{
  C3Ddata oDummy;

  if (NULL != m_poReflnlistRefine) 
    {
      delete m_poReflnlistRefine;
      m_poReflnlistRefine = NULL;
    }
  if (NULL != m_poReflnlistPartial)
  {
      delete m_poReflnlistPartial;
      m_poReflnlistPartial = NULL;
  }

  if (NULL != m_pfPerImageFactors)
    {
      delete [] m_pfPerImageFactors;
      m_pfPerImageFactors = NULL;
    }

  if (NULL != m_pcBackgroundMask)
  {
      delete [] m_pcBackgroundMask;
      m_pcBackgroundMask = NULL;
  }

  if (NULL != m_pcSelectFields)
    {
      delete [] m_pcSelectFields;
      m_pcSelectFields = NULL;
    }
  if (NULL != m_pcSelectFieldsPartial)
    {
      delete [] m_pcSelectFieldsPartial;
      m_pcSelectFieldsPartial = NULL;
    }


  vDeleteFreeShoeboxes();
  vTermRogue();

  if ( (m_bNewScan) && (NULL != m_poScan) )
    {
      if (NULL != m_poScan->m_poNonunf)
        {
          //+jwp 12-June-2001 probably not needed anymore         
          //delete m_poScan->m_poNonunf;
          m_poScan->m_poNonunf = NULL;
        }
      delete m_poScan;
      m_poScan = NULL;
    }

  if (NULL != m_po3DCentroid)
    {
      delete m_po3DCentroid;
      m_po3DCentroid = NULL;
    }
  if (NULL != m_po3DCentroid2)
    {
      delete m_po3DCentroid2;
      m_po3DCentroid2 = NULL;
    }


  if (NULL != m_pfObliqueCorr)
    {
      delete [] m_pfObliqueCorr;
      m_pfObliqueCorr = NULL;
    }
  if (NULL != m_pfObliqueCorrNorm)
    {
      delete [] m_pfObliqueCorrNorm;
      m_pfObliqueCorrNorm = NULL;
    }
  if (NULL != m_pfObliqueCorrType)
    {
      delete [] m_pfObliqueCorrType;
      m_pfObliqueCorrType = NULL;
    }
  if (NULL != m_pfExcludeReso[0])
    {
      delete [] m_pfExcludeReso[0];
      m_pfExcludeReso[0] = NULL;
    }
  if (NULL != m_pfExcludeReso[1])
    {
      delete [] m_pfExcludeReso[1];
      m_pfExcludeReso[1] = NULL;
    }

  if (NULL != m_poProfit2DPerSlice)
    {
      delete m_poProfit2DPerSlice;
      m_poProfit2DPerSlice = NULL;
      C3Ddata::m_poProfit = NULL;
    }

  if (NULL != m_pnHKLIndex)
    delete[] m_pnHKLIndex;


  if (NULL != m_ptRefineResults)
    {
      delete [] m_ptRefineResults;
      m_ptRefineResults = NULL;
    }

  if (NULL != m_pfNearestCentroidBuf0)
    {
      delete[] m_pfNearestCentroidBuf0;
      m_pfNearestCentroidBuf0 = NULL;
    }
  if (NULL != m_pfNearestCentroidBuf1)
    {
      delete[] m_pfNearestCentroidBuf1;
      m_pfNearestCentroidBuf1 = NULL;
    }
    
    DELETE_THE_IMAGE_WAIT_PTR;
}

int
Cintegrate::nInitValues(void)
{
  // nInitValues initializes the values of the object and allocates arrays
  // It is called by the constructor

  int i;
  int nx,ny;


  m_bScaleShoeboxWithShifts = TRUE;
  m_fFattening       = 1.2f;
  m_nShoeboxLimit    = 50;
  m_eSBdatatype      = e3Ddata_ushort;
  m_fSpotSizeMultiplier = 2.8f;
  m_bWriteScanBitmap = FALSE;
  m_bWriteEllipsoids = FALSE;
  m_bSpotChase       = FALSE;
  m_bKa12            = TRUE;
  m_nMinPeakRad      = 2;
  m_nIntegrateWeakFlag = 1;
  m_bDoFastPrediction = FALSE;
  m_fSigmaAboveBackground = 2.8f;

  m_nProfitNumReflections = 50;
  m_nProfitMaxImageRange = 7;

  m_nVerbose         = 1;                   // Set verbose level
  m_nProfitFlag      = DTI_DO_NOT_PROFIT;   // Do not do any profile fitting.
  m_nPad3D           = 0;                   // No padding
  m_fPad3D           = 1.0;                 // No pad multiplier.
  m_nMinStrongSpotsForIntegrate = 8;
  m_nTwinID          = 1;                // Twin ID for this integration.

  m_nImagesPerScaleBatch = 1;            // A batch for every image.
  m_nImagesPerRefineBatch = 4;           // Refine every four images.

  m_anSeqNum[0]      = 0;
  m_anSeqNum[1]      = 0;
  //m_fWaitTime        = 0.0;
  m_nDisplay         = 0;
#ifdef SSI_PC
  m_nDisplayFreq     = 1;
  m_nDisplayCtr		 = 0;
#endif
  m_sOutputHeader    = "";
  m_sXpropString     = "";
  m_a2fResolution[0] = 0.0;    // 0.0 flags that resolution is not set by user
  m_a2fResolution[1] = 99999.99f;
  m_sBatchPrefix     = "";                  // Default prefix is null (empty, none)
  m_nPrerefineFlag   = 2;                   // Usually do prerefinement with images of 2 refine batches
  
  m_nPrefindFlag         = 0;               // By default do not do prefind
  m_fPrefindSigma        = 5.0;

  m_a2fMosaicityModel[0] = 1.0;
  m_a2fMosaicityModel[1] = 0.0;
  m_bDetermineFixedMosaicity = FALSE;
  m_fFixedMosaicity  = -1.0;                // User must specify a value >0 or -1 for other features.
  m_fMaxIntersectionFraction = 0.03f;// Profile percentage.
  m_a3fPCellVector[0] = 0.0;
  m_a3fPCellVector[1] = 0.0;
  m_a3fPCellVector[1] = 0.0;
  m_fPCellAngle = 0.0;



  m_sReflnFilename           = sTransSymbol(sDtrekGetPrefix() + "dtintegrate.ref");
  m_sReflnPartialsFilename   = sTransSymbol(sDtrekGetPrefix() + "dtintpartials.ref");
  m_sReflnRejectsFilename    = sTransSymbol(sDtrekGetPrefix() + "dtintrejects.ref");
  m_sReflnProfitFilename     = sTransSymbol(sDtrekGetPrefix() + "dtprofit.ref");
  m_fDelHKLThresh    = 0.25;         // Pixels with float delHKLsq more than

  m_poImage          = NULL;
  m_bNewImage        = FALSE;
  m_poScan           = NULL;
  m_bNewScan         = FALSE;
  m_poReflnlist      = NULL;
  m_poReflnlistRefine= NULL;
  m_poReflnlistPartial = NULL;
  m_bNewReflnlist    = FALSE;
  m_poRefine         = NULL;
  m_poHeader         = NULL;

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  m_poXprop          = NULL;
#endif

  for (i = 0; i < 20; i++)
    {
      m_afErrorLimits[i] = 0.0;
    }
  for (i = 0; i < nMaxErrorStatusTypes; i++)
    {
      m_anNumReflns[i] = 0;
    }
  m_pcSelectFields   = NULL;
  m_pcSelectFieldsPartial = NULL;
  m_nReflnFileOutCount = 0;
  m_nReflnPartialOutCount = 0;
  m_pFReflnFileOut   = NULL;
  m_pFReflnPartialOut= NULL;
  m_pFReflnBadOut    = NULL;

  
  m_nNumRefNoErrors     = 0;
  m_nNumRefBadErrors    = 0;
  m_nNumPredictedUnique = 0;
  m_nNumRefSpecial      = 0;
  m_nNumRefWeak         = 0;
  m_nNumRefPartialStart = 0;
  m_nNumRefPartialEnd   = 0;
   
  // Refine results

  m_ptRefineResults    = NULL;
  m_nNumRefineResults  = 0;
  m_nSizeRefineResults = 0;

  // Sorting arrays.
  m_pnHKLIndex = NULL;
  // 2D profile fitting data.
  m_poProfit2DPerSlice = NULL;
  m_sProfitDiagnostic = "";

  // Create a shoebox for temporary use during centroid determinations
  m_po3DCentroid  = new C3Ddata(30, 30, 15,
                                0, 0, 0, e3Ddata_float);
  m_po3DCentroid2 = new C3Ddata(30, 30, 15,
                                0, 0, 0, e3Ddata_float);


  m_nRefineFlag             = 1;       // Default is yes, do refinement
  m_nRefineMinRefln         = 50;      // Default 50 reflections min for refine
  m_nRefineMaxRefln         = 50000;   // Default maximum number of reflections to refine.

  m_nDoIntegrateCount       = 0;       // Only used for printing I think.
  m_nDoRefineCount          = 0;       // Indicates which pass of refinement parameters are currently in use.

  m_fRotIncrement = 0.0;
  m_nDim0 = 0;
  m_nDim1 = 0;
  m_a2nWinMax[0]     = 0;             // Default spot window max size from dtfind
  m_a2nWinMax[1]     = 0;             // Default spot window max size from dtfind
  m_a2nWinMin[0]     = 0;
  m_a2nWinMin[1]     = 0;
  m_fDetectorGain    = 2.5;
  m_nNumBeforePred   = 0;         // Loaded before each new image prediction done
  m_nNumAfterPred    = 0;          // Loaded after each new image prediction done


  m_nCurrentImgNum    = 0;                        // per batch

  m_fOverallRotStart = -999.0;
  m_fOverallRotEnd   = -999.0;
  m_nNonunfMaskBadPix   = 1;    // These are masks for the non-uniformity field of the reflection list.
  m_nNonunfMaskRing     = 8;


  m_pfPerImageFactors       = NULL;
  m_pcBackgroundMask        = NULL;

  // Initialize mosaicity profiles.
  for (nx=0;nx<MAX_SHOEBOX_INTENSITY_PROFILE-2;nx++) {
      for (ny=0;ny<MAX_SHOEBOX_INTENSITY_PROFILE;ny++) {
          m_pfShoeboxIntensityProfilesIntensity[nx][ny] = 0.0;
      };
      m_pnShoeboxIntensityProfilesCounts[nx] = 0;
  };

  // Initialize the mosaicity bins.
  for (ny=0;ny<MAX_MOSAICITY_PROFILES;ny++) {
      for (nx=0;nx<MAX_MOSAICITY_PROFILE_SIZE;nx++) {      
        m_aanMosaicityProfilesCounts[ny][nx] = 0;
      };
      m_afMosaicityValues[ny] = 0.0;
  };

  
  m_pfNearestCentroidBuf0 = NULL;
  m_pfNearestCentroidBuf1 = NULL;
  m_nNearestCentroidSize = 0;



  m_bRogueActive = FALSE;
  m_nRogueMask = 0;
  m_nRogueNegMask = 0;
  m_nRogueImageCount = 0;
  m_nRoguePrintNumber = 1;
  m_nRogueShoeboxImageCount = 0;
  m_nRoguePx0 = 0;
  m_nRoguePx1 = 0;
  m_nRogueDim0 = 0;
  m_nRogueDim1 = 0;
  m_nMaxRoguePx1 = 0;
  m_nRogueBorder = 500;
  m_sRogueTemplate = "test???.osc";
  m_poRogueImage = NULL;
  m_poRogueReflnlist = NULL;
  m_pnRogueReflnlistHKLSort = NULL;
  m_nRoguePlaceInImages = 0;
  m_nRogueImageNumber = 1;  // Image number 1 is the first in a scan specified by -seq 1 ? 
  m_nMaxRogueShoeboxes = 0;
  m_nRogueShoeboxes = 0;
  m_nRogueMinBoxDim = -1;
 
    // Try to keep track of shoeboxes and re-use them
  m_ppoShoeBFree      = NULL;
  m_pnShoeBFreeSize   = NULL;
  m_nNumShoeBAlloc    = 0;
  m_nNumShoeBFree     = 0;

  m_nErrorDoNotPrintMask = 
                          nSetErrorMask(Crefln::ms_nErrorSpecial,     FALSE,NULL) 
                          //nSetErrorMask(Crefln::ms_nErrorOnEdge0,  FALSE,NULL)
                          // nSetErrorMask(Crefln::ms_nErrorOffEdge0,  FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorOnEdge1,  FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorOffEdge1,  FALSE,NULL)
                           | nSetErrorMask(Crefln::ms_nErrorOnEdge2,     FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorOffEdge2,    FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorHKLBound,    FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorRing,        FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorTooDark,     FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorBackground,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorBackgroundSD,FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorNonunfA,     FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorNonunfB,     FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorPartialStart,FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorPartialEnd,  FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorRotTooWide,  FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorIntProblem,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorTooFewPix,   FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorUnaccountedPix, FALSE,NULL)
                          ;
                                                  
  m_nErrorDoNotProcessMask =   
                          //nSetErrorMask(Crefln::ms_nErrorSpecial,     FALSE,NULL)
                          nSetErrorMask(Crefln::ms_nErrorOnEdge0,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorOffEdge0,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorOnEdge1,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorOffEdge1,  FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorOnEdge2,     FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorOffEdge2,    FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorHKLBound,    FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorRing,        FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorTooDark,     FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorBackground,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorBackgroundSD,FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorNonunfA,     FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorNonunfB,     FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorPartialStart,FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorPartialEnd,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorRotTooWide,  FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorIntProblem,  FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorTooFewPix,   FALSE,NULL)
                          //| nSetErrorMask(Crefln::ms_nErrorUnaccountedPix, FALSE,NULL)
                          | nSetErrorMask(Crefln::ms_nErrorOverlap,   FALSE,NULL)
                          ;

  
  // Initialize so there are no oblique incidence factors
  m_nNumObliqueCorr   = 0;
  m_pfObliqueCorr     = NULL;
  m_pfObliqueCorrNorm = NULL;
  m_pfObliqueCorrType = NULL;

  m_nNumExcludeResoRings    = 0;
  m_pfExcludeReso[0]        = NULL;
  m_pfExcludeReso[1]        = NULL;

  m_poSocket = NULL;

  return (0);
}


int
Cintegrate::nIntegrate(const int nNumImages)
{

  // Integrate a fine-sliced 3D scan of single crystal diffraction images.
  // 0.  Read in first image that has info for ...
  //
  // 1.  Get info for source, crystal, crystal goniometer, detector,
  //     detector goniometer.
  //
  // 2.  Loop for each batch:
  //      -C. Get centroids from any completed shoeboxes (but not integrated
  //              shoeboxes from last time around)
  //      -B. Refine on those observed centroids
  //       A. Predict reflections that appear upcoming image and
  //            previous ones that were completed
  //       B. Resolve any overlaps with previous or continuing reflections
  //            (keeping most recently calculated position?)
  //       Now integrate the completed shoeboxes with their new calculated
  //          positions based on refinement of the centroids in the batch
  //
  //       C. Loop over images in the batch
  //          a. Loop over reflections that fall on this images
  //             1. Build 3D shoeboxes with windows from the image
  //                a.  New reflections
  //                b.  Continuing reflections
  //                c.  Done reflections (potentially process here, pass 0)
  //             2. Process Done reflections (pass 1?)
  //                a. Correct for non-uniformity of response and synchrotron
  //                b. Compute local background
  //                c. Generate Profit shoebox, may ellipsoidal shoebox first
  //          Next image
  //       D. End of batch? then Refine or refine with current reflections
  //                        centroids determined from refs in this batch
  //          Keep a header around for each refinement.
  //          Repredict spot centroids?, then integrate the 3D boxes
  //       E. Fit Profit shoeboxes at end of batch, compute Rsym, etc as you
  //          go along if requested
  //   Next batch
  //

  
  if( m_nPrefindFlag != 0 )
  {
        if( 0 != nDoPrefind() )
        {
            cout << endl << "WARNING: Pre-find failed" << endl; 
        }
  }

  int i, j;//, k;                 // Loop counters
  int nStat;                   // Local status
  int nBackupPredict;          // Flag whether to backup prediction or not
  int nPrerefineSave = m_nPrerefineFlag;

  // Get initial info for integration.  All the information is derived from
  // an image header found in member poHeader.

  Cimage       oImage;
  Cdetector   *poDetector;
  Ccrystal    *poCrystal;
  Cgoniometer *poCrysGonio;
  Csource     *poSource;
  Cpredict    *poPredict;
  Crefine     *poRefine;
  float        fRotEnd;
  int          nPack;
  LONG         lImageFileSize = 0;
  

  Cstring      sImageName;
  Cstring      sTemp;
  double       a3fTemp[3];

  float        fSourceCountsStart = 1.0;
  float        fSourceIntensity   = 0.0;

//////////////////////////////////////////////////////////////////////////////
#ifdef ANL_TIMER
  anl_reset_timer(0,"nIntegrate3D setup");
  anl_reset_timer(1,"nUpdateRefine");
  anl_reset_timer(2,"new Cpredict");
  anl_reset_timer(3,"new Crefine");
  anl_reset_timer(4,"nExpandReflnlist");
  anl_reset_timer(6,"Setup for first image");
  anl_reset_timer(7,"nPredict");
  anl_reset_timer(8,"Loop thru reflections to pack");
  anl_reset_timer(9,"vSort");
  anl_reset_timer(10,"nBuildShoeboxes");
  anl_reset_timer(11,"Integrate reflections");
  anl_reset_timer(12,"nDoIntegrate");
  anl_reset_timer(13,"Individual end-of-batch");
  anl_reset_timer(14,"nDoRefine of batch");
  anl_reset_timer(15,"nUpdateRefine");
  anl_reset_timer(16,"end-of-scan");
  anl_reset_timer(17,"Close reflection files");
  anl_reset_timer(18,"vListResults");
  anl_reset_timer(20,"End-of-batch");
  anl_reset_timer(21,"Image Processing Time");
  anl_reset_timer(22,"nGetImage");
#endif

  ANLSTART(0);
  time(&m_nStartTime);

  if (NULL == m_poHeader)
    {
      cout << "ERROR in Cintegrate::nIntegrate, header is NULL!\n";
      return (-1);
    }
  else if (!m_poHeader->bIsAvailable())
    {
      cout << "ERROR in Cintegrate::nIntegrate, header is not available!\n";
      return (-2);
    }

  // Initialize refinement results for later summary table output

  if (NULL != m_ptRefineResults)
    delete [] m_ptRefineResults;
  m_nNumRefineResults  = 0;
  m_nSizeRefineResults = 1;
  if (0 < m_nImagesPerRefineBatch)
    m_nSizeRefineResults = nNumImages / m_nImagesPerRefineBatch + 2;
  m_ptRefineResults = new tagRefineResults [m_nSizeRefineResults];


  
  bool bInvertSourceIntensityInfo = FALSE;
  bInvertSourceIntensityInfo = ("" != sGetEnv("DTREK_INVERT_SOURCEINT"));
  if (bInvertSourceIntensityInfo)
    {
      cout << "INFO, source intensity info will be inverted!\n" << flush;
    }

  // At this point we know we have a legitimate header, so try
  // to create required objects from the header

  nStat = 0;

  // Get the detector info

  poDetector = new Cdetector (*m_poHeader);

  if (!poDetector->bIsAvailable())
    {
      delete poDetector;
      cout << "ERROR in Cintegrate::nIntegrate, detector is not available!\n";
      return (-3);
    }

  // Get the source info

  poSource = new Csource(*m_poHeader);
  if (!poSource->bIsAvailable())
    {
      delete poDetector;
      delete poSource;
      cout << "ERROR in Cintegrate::nIntegrate, source is not available!\n";
      return (-4);
    }

  // Get the crystal info

  poCrystal = new Ccrystal(*m_poHeader);
  if (!poCrystal->bIsAvailable())
    {
      delete poDetector;
      delete poSource;
      delete poCrystal;
      cout << "ERROR in Cintegrate::nIntegrate, crystal is not available!\n";
      return (-4);
    }

  // Get the crystal goniometer info

  poCrysGonio = new Cgoniometer(*m_poHeader, Ccrystal::ms_sCrystalPrefix);
  if (!poCrysGonio->bIsAvailable())
    {
      delete poDetector;
      delete poSource;
      delete poCrystal;
      delete poCrysGonio;
      cout << "ERROR in Cintegrate::nIntegrate, crystal goniometer is not available!\n";
      return (-6);
    }
  ANLSTOP(0);



  // Initialize the profile fitting.  We must do this before pre-refinement, because
  // we will use the first few 1D centroids from pre-refinement to initialize the system.

  //  if (m_nProfitFlag & DTI_DO_2D_PROFIT)
    {
      //      cout << "INFOP: nINitProfit2D \n" << flush;
      if (nInitProfit2D(*poDetector))
        {
          printf("ERROR:  In Cintegrate::nInitProfit2D().\n");
          return (-4);
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////

  // Some local rotation objects for predicting in batches

  Crotation oScanRotation;      // Overall rotation of the scan
  Crotation oImageRotation;     // Rotation of the image
  Crotation oPredRotation;      // Rotation object for prediction

  // Do some preliminaries

  ANLSTART(1);
  nStat = nUpdateRefine(poSource, poDetector, poCrystal, poCrysGonio,m_poScan->m_poRotation);
  ANLSTOP(1);

  if (0 != nStat)
    return (nStat);
  
  Cstring sActiveMaskFile;

  sActiveMaskFile = sGetEnv("CC_ACTIVE_MASK");

#if ((defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  // Doctor up the nonunf image for SSI_PC

  // SSI looks for things in the registry instead of the environment, so
  // don't be suprised if the active mask filename is "" from the environment.

  if (   ("" == sActiveMaskFile) 
      && (eNonunf_simple_mask_state == poDetector->m_poNonunf->m_eThe_State) )
  {
       // Grab location of files from registry
      Cstring sPath = CCrclHelper::GetInstance()->GetInstallDir();

       // Append base filename to above directory

       sPath = sFileBuildName(sPath,"Mask Files");
       sPath = sFileBuildName(sPath,"ActiveMask");

       // Is the simple mask file a NxN binned image?
       // Or can we tell?

       Cimage_header *pHeader = &(poDetector->m_poNonunf->m_poNonunfImage->m_oHeader);
       int nDetectors,nStatus,UnbinnedDims[2],n1,n2,nBinX,nBinY;
       Cstring sKeyword,sDetName,*psDetNames;

       //
       // Get the name of the current detector
       // We are assuming that the "current" detector is the last one in the list.
       // 

       nDetectors = 0;
       sKeyword = D_K_DetectorNumber;
       nStatus = pHeader->nGetValue(sKeyword,&nDetectors);
       if (0 == nStatus)
         {
           if (0 < nDetectors)
             psDetNames = new Cstring [nDetectors];
           else
            nStatus = -1;
      }
      if (0 == nStatus)
        {
          sKeyword = D_K_DetectorNames;
          nStatus = pHeader->nGetValue(sKeyword,nDetectors,psDetNames);
        }
      if (0 == nStatus)
        {
          sDetName = psDetNames[nDetectors-1];
          sKeyword = sDetName;
          sKeyword += D_K_UnbinnedDimensions;
          nStatus = pHeader->nGetValue(sKeyword,2,UnbinnedDims);
        }
      if (0 == nStatus)
        {
          poDetector->m_poNonunf->m_poNonunfImage->nGetDimensions(&n1,&n2);
          nBinX = nBinY = 0;
          if (0 != n1)
            nBinX = UnbinnedDims[0]/n1;
          if (0 != n2)
            nBinY = UnbinnedDims[1]/n2;
        }
       if (0 != nStatus)
         {
           nBinX = nBinY = 0;
         }
       if (0 < nDetectors)
         delete [] psDetNames;
       
      if (0 != nBinX && 0 != nBinY)
        {
          sPath += nBinX;
          sPath += 'x';
          sPath += nBinY;
        }
       sPath += ".img";

       cout << "looking for active mask filename = " << sPath;

      // Does this file exist?  If so, use it as the active mask filename

       if (bFileExists(sPath))
         sActiveMaskFile = sPath;
     }

#endif

  if (   ("" != sActiveMaskFile)
      && (eNonunf_simple_mask_state == poDetector->m_poNonunf->m_eThe_State) )
    {
      // Active mask filename exists!

      Cimage *poImgActiveMask;

      poImgActiveMask = new Cimage(sActiveMaskFile);
      if (poImgActiveMask->bIsAvailable())
        {
          cout << "Combining active mask file with simple mask file.\n" << flush;
          j = poImgActiveMask->nGetDimension();

          // is the active mask the same size as the simple mask?
          
          if (poDetector->m_poNonunf->m_poNonunfImage->nGetDimension() == j)
            {
              unsigned short int uiMask;
              poImgActiveMask->nSetNextPixel(0,0);
              for (i = 0; i < j; i++)
                {
                  uiMask = poImgActiveMask->uiGetNextPixel();
                  if (0 == uiMask)
                    {
                      // BAD programming practice since the member variables are public!
                      // A 0 in the active mask should be propagated to the
                      // simple mask used by the Cnonunf object.

                      poDetector->m_poNonunf->m_poNonunfImage->vSetPixel(i, uiMask);
                    }
                }
            }
          // poDetector->m_poNonunf->m_poNonunfImage->nWrite("simplemask.img");
        }
      delete poImgActiveMask;
    }

//-jwp 28-July-1999

////////////////////////////////////////////////////////////////////////

  // Declare a place to put reflections ...

  Creflnlist *poReflnlistL;
  poReflnlistL = new Creflnlist();

  // Share source, detector, crystal and crystal goniometer objects

  ANLSTART(2);
  poPredict = new Cpredict (poSource, poDetector, poCrystal, poCrysGonio);
  poPredict->vSetNonunfFlag(m_nNonunfMaskBadPix);  // Do not predict refs in bad nonunf regions
  poPredict->vSetImageDriftOn();
  poPredict->vSetMosaicityLinearCoeffs();
  ANLSTOP(2);

  if (0.0 < m_a2fResolution[0])
    {
      // Somewhere the resolution was set, so use it

      poPredict->vSetResolution(m_a2fResolution[0],
                                m_a2fResolution[1]);
    }

  ANLSTART(3);
  //+jwp  12-Jun-2001
  m_poScan->m_poNonunf = poDetector->m_poNonunf;
  //-jwp
  poRefine  = new Crefine  (poSource, poDetector, poCrystal, poCrysGonio,
                            m_poScan->m_poRotation);
  ANLSTOP(3);
//  poRefine  = new Crefine(*m_poHeader, &oReflnlistRefine);

  // Give the reflection list the fields needed for predict and refine objects

  (void) poPredict->nExpandReflnlist(poReflnlistL);
  (void) poRefine->nExpandReflnlist(poReflnlistL);


  // Expand Reflnlist to have fields we need here in Cintegrate

  ANLSTART(4);
  (void) nExpandReflnlist(poReflnlistL);

  ANLSTOP(4);

  //  poReflnlistL->nListFields();

  // A local pointer to a refln for convenience

  Crefln  *poRefln;

////////////////////////////////////////////////////////////////////////

  // If the overwrite environment variable is set, then make sure
  // any existing reflnlist file gets a version number appended, so it is
  // not overwritten

  nStat = nFileAppendVersion(m_sReflnFilename, TRUE);
  if (0 != nStat)
    {
      cout << "ERROR renaming reflnlist output file: " << m_sReflnFilename << "!\n";
#ifdef SSI_PC
      return (2);
#else
      exit (2);
#endif
      /////////
    }

  nStat = nFileAppendVersion(m_sReflnPartialsFilename, TRUE);
  if (0 != nStat)
    {
      cout << "ERROR renaming reflnlist output file: " << m_sReflnPartialsFilename << "!\n";
#ifdef SSI_PC
      return (2);
#else
      exit (2);
#endif
      /////////
    }

  // Select the proper fields to write out (actually mark the
  //   fields to NOT write out).

  m_pcSelectFields = new char [poReflnlistL->m_nTotalFieldsPlus1];
  m_pcSelectFieldsPartial = new char [poReflnlistL->m_nTotalFieldsPlus1];
  for (i = 0; i < poReflnlistL->m_nTotalFieldsPlus1; i++)
    {
      m_pcSelectFields[i] = '\0';  
      m_pcSelectFieldsPartial[i] = char(1);
    }
  m_pcSelectFields[poReflnlistL->m_nTotalFieldsPlus1-1] = '\0';
  m_pcSelectFieldsPartial[poReflnlistL->m_nTotalFieldsPlus1-1] = '\0';

  
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fProfitIntensity, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fProfitSigmaI, char(1),
                                    m_pcSelectFields);
/*
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fProfitErrorModelBackgroundVar, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fProfitErrorModelFittingVar, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fProfitErrorModelCounts, char(1),
*/
  (void) poReflnlistL->nSelectField(eReflnField_int_type,
                                    poReflnlistL->m_nFI_nPartialFlag, char(1),
                                    m_pcSelectFields);

  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fObsXmm, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fObsYmm, char(1),
                                    m_pcSelectFields);
//  (void) poReflnlistL->nSelectField(eReflnField_float_type,
//                                    poReflnlistL->m_nFI_fObsRotStart, char(1),
//                                    m_pcSelectFields);
//  (void) poReflnlistL->nSelectField(eReflnField_float_type,
//                                    poReflnlistL->m_nFI_fObsRotEnd, char(1),
//                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fCalcXmm, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fCalcYmm, char(1),
                                    m_pcSelectFields);
//  (void) poReflnlistL->nSelectField(eReflnField_float_type,
//                                    poReflnlistL->m_nFI_fCalcRotStart, char(1),
//                                    m_pcSelectFields);
//  (void) poReflnlistL->nSelectField(eReflnField_float_type,
//                                    poReflnlistL->m_nFI_fCalcRotEnd, char(1),
//                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fPartiality, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fRecipCoord0, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fRecipCoord1, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fRecipCoord2, char(1),
                                    m_pcSelectFields);
// tjn: add field (m_nFI_nPackedHKL):  when reading the reflection file into scaling or profile fitting,
// it's nice to already have this field in the list to avoid pointer re-allocation.
//  (void) poReflnlistL->nSelectField(eReflnField_int_type,
//                                    poReflnlistL->m_nFI_nPackedHKL, char(1),
//                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_int_type,
                                    m_nFI_nDataFlag, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    m_nFI_fWidth0, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    m_nFI_fWidth1, char(1),
                                    m_pcSelectFields);
 (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                   poReflnlistL->m_nFI_fCalcMosCoeffA, char(1),
                                   m_pcSelectFields);
 (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                   poReflnlistL->m_nFI_fCalcMosCoeffB, char(1),
                                   m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fEllipsoidA00, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fEllipsoidA01, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fEllipsoidA11, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fEllipsoidb0, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fEllipsoidb1, char(1),
                                    m_pcSelectFields);
  (void) poReflnlistL->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fEllipsoidc, char(1),
                                    m_pcSelectFields);
  
  // Create the partials list.
  if (m_poReflnlistPartial)
      delete m_poReflnlistPartial;
  m_poReflnlistPartial = new Creflnlist(*poReflnlistL);

  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fObsRotMid, char(0),
                                    m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                    poReflnlistL->m_nFI_fCalcRotMid, char(0),
                                    m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_int_type,
                                    0, char(0),
                                    m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_int_type,
                                    1, char(0),
                                    m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_int_type,
                                    2, char(0),
                                    m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                    0, char(0),
                                    m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                    1, char(0),
                                    m_pcSelectFieldsPartial);

  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fProfitIntensity, char(0),
                                            m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fProfitSigmaI, char(0),
                                            m_pcSelectFieldsPartial);

  (void) m_poReflnlistPartial->nSelectField(eReflnField_int_type,
                                            poReflnlistL->m_nFI_nPartialFlag, char(0),
                                            m_pcSelectFieldsPartial);

  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fEllipsoidA00, char(0),
                                            m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fEllipsoidA01, char(0),
                                            m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fEllipsoidA11, char(0),
                                            m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fEllipsoidb0, char(0),
                                            m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fEllipsoidb1, char(0),
                                            m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fEllipsoidc, char(0),
                                            m_pcSelectFieldsPartial);
  (void) m_poReflnlistPartial->nSelectField(eReflnField_float_type,
                                            poReflnlistL->m_nFI_fBackgroundSigma, char(0),
                                            m_pcSelectFieldsPartial);
  

  // Make the refinement reflection list.
  if (m_poReflnlistRefine)
      delete m_poReflnlistRefine;
  m_poReflnlistRefine = new Creflnlist(*poReflnlistL);



////////////////////////////////////////////////////////////////////////
  // Initialize the per image factors

  if (NULL != m_pfPerImageFactors)
    delete [] m_pfPerImageFactors;

  m_pfPerImageFactors = new float [nNumImages];
  for (i = 0; i < nNumImages; i++)
    {
      m_pfPerImageFactors[i] = 1.0;
    }

  
  if (poDetector->nCalcPixelShiftArray(*poSource,*m_poScan->m_poRotation,*poCrysGonio))
      return 1;
  //poDetector->nTest(*poSource,*m_poScan->m_poRotation);



 // Copy overall scan rotation to local variable and manipulate it

  //+2009-12-08 JWP  Backdoor method to run pre-refinment multiple times if envvar is > 1
  int nPrerefineLoops = 1;
  if ("" != sGetEnv("DTREK_DTINT_PREREFINELOOPS"))
    {
      nPrerefineLoops = atoi(sGetEnv("DTREK_DTINT_PREREFINELOOPS").string());
      printf("INFO:  Pre-refine loop set by environment: %d\n\n", nPrerefineLoops);
    }
  //-2009-12-08 JWP

pre_refine_rollback:

  // Initialization of counters.
  m_nCurrentImgNum    = 0;                        // per batch
  
  m_nNumRefineResults  = 0;
  m_nDoIntegrateCount       = 0;       // Only used for printing I think.
  m_nDoRefineCount          = 0;       // Indicates which pass of refinement parameters are currently in use.


  m_nNumRefNoErrors     = 0;
  m_nNumRefBadErrors    = 0;
  m_nNumPredictedUnique = 0;
  m_nNumRefSpecial      = 0;
  m_nNumRefWeak         = 0;
  m_nNumRefPartialStart = 0;
  m_nNumRefPartialEnd   = 0;


  // Reflection list creation.
  if (m_pFReflnFileOut) {
      fclose(m_pFReflnFileOut);
      nFileDelete(m_sReflnFilename);
      m_nReflnFileOutCount = 0;
  };
  if (m_pFReflnPartialOut) {
      fclose(m_pFReflnPartialOut);
      nFileDelete(m_sReflnPartialsFilename);
      m_nReflnPartialOutCount = 0;
  };
  if (m_pFReflnBadOut) {
      fclose(m_pFReflnBadOut);
      nFileDelete(m_sReflnRejectsFilename);
  };

  m_pFReflnBadOut       = fopen(m_sReflnRejectsFilename.string(),"wb");
  m_pFReflnPartialOut   = fopen(m_sReflnPartialsFilename.string(),"wb");
  m_pFReflnFileOut      = fopen(m_sReflnFilename.string(),"wb");
  
  if ((!m_pFReflnFileOut) || (!m_pFReflnPartialOut) || (!m_pFReflnBadOut))
  {
      cout << "Error opening " << m_sReflnFilename << "!\n";
#ifdef SSI_PC
      return (2);
#else
      exit (2);
#endif
      /////////
  }

  nDeleteMarkedReflns(poReflnlistL,-999);
  m_poReflnlistPartial->vDeleteAll();

  sTemp = "";
  poReflnlistL->eSetWriteBinary(eReflnoutput_default);
  poReflnlistL->nPutCrystal(*m_poHeader);
  poReflnlistL->nWrite(m_sReflnFilename, NULL, m_pcSelectFields, m_pFReflnFileOut);

  if (NULL != m_pFReflnBadOut)
    poReflnlistL->nWrite(m_sReflnRejectsFilename, NULL, m_pcSelectFields,m_pFReflnBadOut);

  m_poReflnlistPartial->eSetWriteBinary(eReflnoutput_default);
  m_poReflnlistPartial->nPutCrystal(*m_poHeader);
  m_poReflnlistPartial->nWrite(m_sReflnPartialsFilename, NULL, m_pcSelectFieldsPartial, m_pFReflnPartialOut);


  // Open the reflection file and select the fields to write out
  // It must be done here, after the reflnlist has already been expanded
  // to have all the proper fields

  m_poScan->vInitSeqNum();
  m_poScan->vSetSeqNum(m_anSeqNum[0]);
  oScanRotation  = *m_poScan->m_poRotation;
  oScanRotation.vSetRotStart(m_poScan->fCalcGetRotStart(m_anSeqNum[0]));
  oScanRotation.vSetRotEnd(m_poScan->fCalcGetRotEnd(m_anSeqNum[1]));

  if (m_bDoFastPrediction) {
      poPredict->vSetFastPredict(oScanRotation.fGetRotStart(),oScanRotation.fGetRotEnd());
  };


  m_fOverallRotStart = oScanRotation.fGetRotStart();
  m_fOverallRotEnd   = oScanRotation.fGetRotEnd();
  m_fRotIncrement    = oScanRotation.fGetIncrement();

  oImageRotation = oScanRotation;         // Copy scan rotation to image rotation
  
  {
      float a3fTempBuf[3];
      int nStatL;
      nStatL = poCrysGonio->nGetNum(m_poScan->m_poRotation->sGetName());
      if (poCrysGonio->nGetRotVector(nStatL,&a3fTempBuf[0]))
	{
	  cout << "ERROR, problem with rotation axis name: "
               << m_poScan->m_poRotation->sGetName() << "\n";
          return 1;
	}
      vCopyVec3D(a3fTempBuf,m_a3fRotVec);
      // m_poScan->m_poRotation->vGetVector(m_a3fRotVec);// Get rotation axis to memb var
  };

  

  // Initialize ending rotation, so that start gets set correctly later on

  oImageRotation.vSetRotEnd(m_fOverallRotStart);

  if (1 > m_nImagesPerRefineBatch)    {
      cout << "INFO: Images per refinement batch reset to 2\n" << flush;
      m_nImagesPerRefineBatch = 2;
  }

  nStat               = 0;
  nBackupPredict      = 0;

  m_poScan->nGetImageName(&sImageName);

  if( THE_IMAGE_WAIT_PTR->bWaitRequired() && !THE_IMAGE_WAIT_PTR->bDoWaitForFirstImage(sImageName) )
  {
        nStat = 1; // couldn't open the very first image
  }

  // Process images in the scan until status is not 0
  // or until the end of the scan or requested number of images

  //mlw 1/21/98
  //mlw changed while loop end for parallel code
  int nReadStat = 0;
  int nPercentFilesize = 95;
  if ("" != sGetEnv("DTREK_DTINT_PERCENTFILESIZE"))
    {
      nPercentFilesize = atoi(sGetEnv("DTREK_DTINT_PERCENTFILESIZE").string());
      cout << "INFO: Image file size required to be at least "
           << nPercentFilesize << "% of the first image size." << endl << flush;
    }

  while(   (0 == nStat)
         && (oImageRotation.fGetRotEnd() < oScanRotation.fGetRotEnd())
         && (m_nCurrentImgNum < nNumImages)
       )
    {
      ANLSTART(21);
      m_poScan->nGetImageName(&sImageName);
      
      if( THE_IMAGE_WAIT_PTR->bWaitRequired() && !THE_IMAGE_WAIT_PTR->bDoWaitForNextImage(sImageName) )
      {
          break; // We did not get this image, but we still want to return without error.
      }

      ANLSTART(22);
      oImage.m_bReadApplyGonioMask    = TRUE;
      oImage.m_bReadApplyEmbeddedMask = TRUE;
      nStat = m_poScan->nGetImage(&oImage);
      ANLSTOP(22);
      ANLPRINT(22);

      if ( (0 != nStat) || (!oImage.bIsAvailable()) )
        {
          // ERROR: Image is not available!

	  nReadStat = 1;
	  //cout << "ERROR: setting nReadStat to 1\n";
          //cout << flush;
          cout << "ERROR: problem reading image " << sImageName << "!\n";
          if (0 == m_nCurrentImgNum)
            {
              // Die a not so graceful death.

              delete poReflnlistL;
              delete poRefine;
              delete poDetector;
              delete poSource;
              delete poCrystal;
              delete poCrysGonio;
              delete poPredict;

              // Danger! memory leak here?  SSI/CC folks watch out!

              return (-1);
            }
        }
      else
        {
          oImage.nGetDimensions(&m_nDim0, &m_nDim1);

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
          if (NULL != m_poXprop)
            {

              (void) m_poScan->nGetImageName(&sTemp);
              sTemp =  sTemp + " Template: "
                             + m_poScan->sGetTemplate()
                             + " Seq: "
                             + Cstring(m_poScan->nGetSeqNum());
//              poXprop->hSetProperty("DTDISPLAY_IMAGE_UPDATE",
//                                  sTemp + " //dtintegrate//");
            }
#endif

          m_poScan->vNextSeqNum();        // Set next image sequence number

          // Set rotation start to be the end of the previous image.
          // Compute fRotEnd in way to avoid roundoff accumulation

          oImageRotation.vSetRotStart(oImageRotation.fGetRotEnd());

          fRotEnd = m_fOverallRotStart
                   + (m_fRotIncrement
                      * (float)(m_nCurrentImgNum+1));
          oImageRotation.vSetRotEnd(fRotEnd);

          if (0 < m_nVerbose)
            {
              time_t nTime;
              time(&nTime);
              printf("\n==================================================================");
              if (m_nPrerefineFlag)
                {
                  if (m_nCurrentImgNum+1 > (m_nImagesPerRefineBatch * nPrerefineSave))
                    nPrerefineSave++;
                  printf("\nPRE-REFINE IMAGE #: %d (restart expected after image %d of %d)",m_nCurrentImgNum + 1,
                         min(nNumImages,m_nImagesPerRefineBatch * nPrerefineSave), nNumImages);
                }
              else
                {
                  printf("\nIMAGE #: %d (of %d)",m_nCurrentImgNum + 1,
                         nNumImages);
                }
              printf("\n   Name: %s",oImage.sGetName().string());
              printf("\n   Date: %s       Time: %s",sGetDate().string(),
                     sGetTime().string());
              printf("       Elapsed: %.0lf",difftime(nTime,m_nStartTime));
              printf("\n==================================================================\n");
              oImageRotation.nList(3);
            }

          /////////////////////////////////////////////
          //   Deal with variable source intensity   //
          //     Note: to override header use:       //
          //        'setenv SOURCE_INTENSITY 0.0'    //
          /////////////////////////////////////////////

          fSourceIntensity = 0.0;         // Reset

          // Look for special override keyword

          char a10cTemp[10];
//+jwp 1-Aug-2001
          //          sprintf(a10cTemp, "_%03d", m_poScan->nGetSeqNum( )-1);
          sprintf(a10cTemp, "_%04d", m_poScan->nGetSeqNum( )-1);
//-jwp 1-Aug-2001
          (void) oImage.m_oHeader.nGetValue(Csource::ms_sSourceIntensity
                                            + a10cTemp,
                                            &fSourceIntensity);

          if (0.0 >= fSourceIntensity)
            (void) oImage.m_oHeader.nGetValue(Csource::ms_sSourceIntensity,
                                              &fSourceIntensity);

          if (0.0 < fSourceIntensity)
            {
                        if (bInvertSourceIntensityInfo)
                                fSourceIntensity = 1.0 / fSourceIntensity;

              // Found a source intensity keyword in the header, so use it.
              // If not found, since the factors are set to 1.0, nothing happens

                        ANLSTART(6);
              if (0 == m_nCurrentImgNum)
                {

                  // First image, initialize the per image scale factor.
                  // Scale to source intensity given in SOURCE_INTENSITY_INIT
                  // keyword or SOURCE_INTENSITY in first image header.

                  int nStatL;
                  float fSourceL;
                  nStatL = oImage.m_oHeader.nGetValue(
                               Csource::ms_sSourceIntensity + "_INIT",
                                          &fSourceL);
                  if (0 == nStatL)
                    {
                      // Compare to given input SOURCE_INTENSITY_INIT

                      fSourceCountsStart = fSourceL;
                    }
                  else
                    {
                      // Compare to first image

                      fSourceCountsStart = fSourceIntensity;
                    }
                }
              else
                {
                  // Not first image, compute and store per image scale factor

                  m_pfPerImageFactors[m_nCurrentImgNum] = fSourceCountsStart
                                                          / fSourceIntensity;
                }
              if (0 < m_nVerbose)
                {
                  cout << "\nImage source intensity: "
                       << setw(6) << setprecision(5)
                       << fSourceIntensity
                       << "\nImage intensity factor: "
                       << setw(6) << setprecision(5)
                       << m_pfPerImageFactors[m_nCurrentImgNum] << '\n' << flush;
                }
            }

          if (0 == m_nCurrentImgNum)
            {
              // If this is the first image in the scan do some preliminary
              // things that should help with memory management:
              //
              //
	      // 0.  Figure out if images are ushort or (long) int or ???.
	      //
              // 1.  Predict a lot of reflections in order to allocate
              //     a list of reflections next to each other in memory.
              //     Then use reflnlist->vDeleteAll() to put them in the
              //     reflnlist free list for re-use.
              //
              // 2. Create lots of empty shoeboxes, so they are all close
              //    together in memory.
              //
              // 3. Allocate memory for the shoeboxes with an average size.
              //
              // 4. Delete the memory for the shoeboxes so it goes back on the
              //    free list
              //
              // 5. Fill in spot ellipsoid information.
              //
              // 6. Initialize the detector area information.
              //
              // 7. Initialize profile fitting if it's enabled.
              //
              /////////////////////////////////////////////////////////////////

	      // Figure out the internal shoebox data type for this scan integration
	      // based on the pixel value data type in the image

	      m_eSBdatatype      = e3Ddata_ushort;
	      cout << "INFO: eImage_data_type is " << oImage.sGetDataType() << endl;
	      if (eImage_I4 == oImage.m_oHeader.eGetImageDataType())
		{
		  // We are going to make the C3Ddata float to start with
		  // since this will save the expense of an extra conversion

		  m_eSBdatatype      = e3Ddata_float;

		  // Also set the minimum OK value in a raw pixel to 0 unless
		  // overridden by the environment variable.

		  if ("" == sGetEnv("DTREK_NONUNF_OKVALUE"))
		    {
		      poDetector->m_poNonunf->vSetMinRawPixOKValue(0.0);
		      cout << "INFO: Setting min raw pixel OK value to 0.\n";
		    }
		}

              oPredRotation  = oImageRotation; // Copy image rot to pred rot
              oPredRotation.vSetRotStart( m_fOverallRotStart);
              oPredRotation.vSetRotEnd( m_fOverallRotStart + m_fRotIncrement
                                       * (float) (m_nImagesPerScaleBatch+2));
              poPredict->nPredict(&oPredRotation, poReflnlistL);
              m_nNumAfterPred =  poReflnlistL->nGetNumReflns();

              // Get average reflection width

              float fAvgWidth = 0.0;
              for (j = 0; j < m_nNumAfterPred; j++)
                {
                  poRefln    = poReflnlistL->poGetRefln(j);
                  fAvgWidth +=
                    poRefln->fGetField(poReflnlistL->m_nFI_fCalcRotWidth);
                }
              fAvgWidth = fAvgWidth / (float)m_nNumAfterPred;

              m_nNumAfterPred =  poReflnlistL->nGetNumReflns();

              // Delete reflns in list (they go onto free list)

              poReflnlistL->vDeleteAll();

              // Save the image file size
               if( THE_IMAGE_WAIT_PTR->bWaitRequired() && !THE_IMAGE_WAIT_PTR->bSaveImageSize(sImageName) )
                return -1;

              if ((m_a2oDetectorProfiles[0].fGetAverageLevel() == 0.0) &&
                  (m_a2oDetectorProfiles[1].fGetAverageLevel() == 0.0)) {
                  // Initialize peak info if we haven't already.
                  // (remember, we might have pre-refined already)

                  m_nProfileLoading = 1;
                  m_nProfileReading = 0;

                  m_fTargetDetectorProfilesLevel = 2.5;
                  if (m_poHeader) {
                      if (m_a2oDetectorProfiles[0].nInitValues(*m_poHeader))
                          m_nProfileLoading = 0;
                  };
                                    
                  m_a2oDetectorProfiles[m_nProfileLoading].nInitPeakAreas(m_nDim0,m_nDim1);
              } 
                          
            }
          ANLSTOP(6);

          ///////////////////////////////////////////////////////////
          // Predict reflections that fall on this image.          //
          // Append newly predicted reflections to end of the list //
          ///////////////////////////////////////////////////////////

          m_nNumBeforePred = poReflnlistL->nGetNumReflns();

          // Adjust prediction of reflections by the amount of padding
          // BUT do not predict reflections that END  BEFORE the overall rotation start
          // OR                             that BEGIN AFTER the overall rotation end.

          oPredRotation  = oImageRotation; // Copy image rot to pred rot
          oPredRotation.vSetRotStart( max(oPredRotation.fGetRotStart()
                                          - (m_fRotIncrement
                                             * (float) m_nPad3D)
                                          - (m_fRotIncrement
                                             * (float) nBackupPredict),
                                          -9999.0 /*m_fOverallRotStart */)
                                     );
          oPredRotation.vSetRotEnd( min(oPredRotation.fGetRotEnd()
                                        + (m_fRotIncrement
                                           * (float) m_nPad3D),
                                        9999.0 /*m_fOverallRotEnd */)
                                   );

          // Then predict reflections for this range

          if (2 < m_nVerbose) cout << "Predicting reflns...\n" << flush;
          ANLSTART(7);
          poPredict->nPredict(&oPredRotation, poReflnlistL);

          ANLSTOP(7);

          // Set the bounds here (after prediction) for consistency.
          // We only predicted beyond the bounds so that we could get partial reflections that
          // might linger on the scan limits, but are not predicted by the mosaicity.
          oPredRotation.vSetRotStart(max(oPredRotation.fGetRotStart(),m_fOverallRotStart));
          oPredRotation.vSetRotEnd(min(oPredRotation.fGetRotEnd(),m_fOverallRotEnd));

          if (2 < m_nVerbose) cout << "Done predicting reflns...\n" << flush;

          // Mark reflections in excluded rings

          int nRingField
            = poReflnlistL->nGetFieldIndex(poReflnlistL->ms_snNonunfFlag);
          int nRingFlag;
          for (i = 0; i < m_nNumExcludeResoRings; i++)
            {
              for (j = m_nNumBeforePred; j < poReflnlistL->nGetNumReflns(); j++)
                {
                  poRefln = poReflnlistL->poGetRefln(j);
                  if ( (m_pfExcludeReso[0][i]
                        >  poRefln->fGetField(poReflnlistL->m_nFI_fResolution))
                      &&
                      (m_pfExcludeReso[1][i]
                       <  poRefln->fGetField(poReflnlistL->m_nFI_fResolution)) )
                    {
                      // Mark reflection as bad

                      nRingFlag = poRefln->nGetField(nRingField) | m_nNonunfMaskRing;
                      poRefln->vSetField(nRingField, nRingFlag);
                    }
                }
            }

                  if (poReflnlistL->nGetNumReflns() > m_nNearestCentroidSize) {
                          delete[] m_pfNearestCentroidBuf0;
                          m_pfNearestCentroidBuf0 = NULL;
                          delete[] m_pfNearestCentroidBuf1;
                          m_pfNearestCentroidBuf1 = NULL;
                          m_nNearestCentroidSize = poReflnlistL->nGetNumReflns();
                          m_pfNearestCentroidBuf0 = new float[m_nNearestCentroidSize];
                          m_pfNearestCentroidBuf1 = new float[m_nNearestCentroidSize];
                  };
                  poPredict->nCalcNearestCentroids(*poReflnlistL,m_nNumBeforePred,-1,
                          oImageRotation.fGetRotStart(),
                          oImageRotation.fGetRotEnd(),
              30,
                          m_pfNearestCentroidBuf0,
                          m_pfNearestCentroidBuf1);

  
//+2009-12-07 JWP
// Try to have more of these refln files since computers are faster
          static int s_nAB = 0;
          s_nAB++;
	  if (5 < s_nAB) s_nAB = 0;
//-2009-12-07 JWP

          /////////////////////////////////////////////////////////////////////////////////
          // Write out only the newly predicted reflections and send update
          if( 1 == m_nDisplay )
          {
#ifdef SSI_PC  
			  if( m_nDisplay == 0 || (m_nDisplayCtr % m_nDisplayFreq) == 0 ) {
				  char      a255cReffile[255];
	          
				  sprintf(a255cReffile, sDtrekGetPrefix() + "dtintpred%1d.ref",  s_nAB);
	          
				  poReflnlistL->vSetOverwrite(TRUE);
	          
				  poReflnlistL->nWrite((Cstring)a255cReffile, NULL, NULL, NULL, m_nNumBeforePred);
	          
				  poReflnlistL->vSetOverwrite(FALSE);
		          
				  if (m_poSocket)
				  {
					  Cstring reffile = sGetCWD() + a255cReffile;
					  m_poSocket->SendString("dtUpdateDisplay Image = {" + oImage.sGetName() + "} ReflnFile = {" + reffile + "} Module = dtintegrate");
				  }

	              m_sXpropString = a255cReffile;
				  CCrclHelper::GetInstance()->vSendIntegrateUpdateDisplay(oImage.sGetName(), a255cReffile, FALSE, TRUE, FALSE);
			  }
			  m_nDisplayCtr++;
#else
              char      a255cReffile[255];
          
              sprintf(a255cReffile, sDtrekGetPrefix() + "dtintpred%1d.ref",  s_nAB);
          
              poReflnlistL->vSetOverwrite(TRUE);
          
              poReflnlistL->nWrite((Cstring)a255cReffile, NULL, NULL, NULL, m_nNumBeforePred);
          
              poReflnlistL->vSetOverwrite(FALSE);
	          
              if (m_poSocket)
              {
                  Cstring reffile = sGetCWD() + a255cReffile;
                  m_poSocket->SendString("dtUpdateDisplay Image = {" + oImage.sGetName() + "} ReflnFile = {" + reffile + "} Module = dtintegrate");
              }
#endif

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
              if( NULL != m_poXprop )
              {
                  sTemp = sTemp + " Reflnlist: " + sGetCWD() + (Cstring)a255cReffile;
	              m_sXpropString = sGetCWD() + (Cstring)a255cReffile;
                  m_poXprop->hSetProperty("DTDISPLAY_IMAGE_UPDATE", sTemp + " //dtintegrate//");
              }
#endif

/*#ifdef SSI_PC              
              m_sXpropString = a255cReffile;
			  if( m_nDisplay == 0 || (m_nDisplayCtr % m_nDisplayFreq) == 0 ) {
				CCrclHelper::GetInstance()->vSendIntegrateUpdateDisplay(oImage.sGetName(), a255cReffile, FALSE, TRUE, FALSE);
			  }
			  m_nDisplayCtr++;
#endif*/
          }

          // Now loop through all reflections predicted to be on this image
          // Create packed indices for the newly predicted reflections


          ANLSTART(8);
          m_nNumAfterPred =  poReflnlistL->nGetNumReflns();

          if (m_pnHKLIndex)
            delete[] m_pnHKLIndex;
          m_pnHKLIndex = new int [m_nNumAfterPred];

          for (j = m_nNumBeforePred; j < m_nNumAfterPred; j++)
            {
              poRefln = poReflnlistL->poGetRefln(j);
              nPack = poRefln->nPackHKL();
              poRefln->vSetField(poReflnlistL->m_nFI_nPackedHKL, nPack);
              poRefln->vSetField(m_nFI_nDataFlag, (int) ms_nNewFlag);
              poRefln->vSetField(poReflnlistL->m_nFI_nRefineBatch,m_nDoRefineCount);
            }
          ANLSTOP(8);

          // Now that all the new reflections have a packed index,
          // sort the reflections on the packed index.  In this way,
          // reflections with identical indices will be next to each
          // other
          ANLSTART(9);
          poReflnlistL->vSort(eReflnField_int_type,
                              poReflnlistL->m_nFI_nPackedHKL,
                              m_pnHKLIndex);

          ANLSTOP(9);

          // Now go process the list of reflections and build shoeboxes.
          // When shoebox is full, correct for nonunf and find special
          // centroid (all done in nBuildShoeboxes)

          nStat = nBuildShoeboxes(poReflnlistL, &oImage,
                                  &oImageRotation, poDetector
                                  );
          ANLSTOP(10);

          if (0 < nBackupPredict)
            {
              // This is a flag that some refinement was done, and
              // new centroids for done reflections were predicted.
              // This means we can integrate them.

              // Reset the backup of predictions back to 0, so we
              // do not enter this block before it is necessary

              nBackupPredict = 0;
            }


        }  // End of processing for the image

      m_nCurrentImgNum++;           // Increment image number in the scan

#ifdef SSI_PC
      // If it is ever desired that CC be alerted when an image is finished
      // integrating, this is the place to send a message.
#endif

      // Check if end of batch or end of scan

      if ((0 == (m_nCurrentImgNum % m_nImagesPerRefineBatch)) ||
          (oImageRotation.fGetRotEnd() >= oScanRotation.fGetRotEnd())
          || (m_nCurrentImgNum == nNumImages))
        {
          ANLSTART(13);
          ANLSTART(20);
          // Current image is last image in a scan, so do some end of batch processing

          if (2 < m_nVerbose)
            {
              cout << "... end of batch ..." << endl << flush;
            }

          /////////////////////////////////////////////////////
          //                                                 //
          // Refine crystal, detector and source properties  //
          // based on completed shoeboxes with strong spots  //
          //                                                 //
          /////////////////////////////////////////////////////

          if (1 == m_nRefineFlag)
            {
              nBackupPredict = m_nImagesPerRefineBatch;  // Backup the next prediction by 3(?) images
              ANLSTART(14);
              if ( (0 < m_nSizeRefineResults) || (NULL == m_ptRefineResults) )
                m_ptRefineResults[m_nNumRefineResults].nSeq 
                  = m_poScan->nGetSeqNum()-1;

              nStat = nDoRefine(poReflnlistL, poRefine);
              ANLSTOP(14);

              // Need next statement, so that when we calculate hkl for each
              // pixel, we have the correct Source, Detector and Crystal props
              
              if (!nStat) 
              {
                  ANLSTART(15);
                  (void) nUpdateRefine(poSource, poDetector,
                      poCrystal, poCrysGonio,m_poScan->m_poRotation);
                  ANLSTOP(15);
                  
                  // Force use of mosaicity found in poCrystal and not in poPredict
                  
                  poCrystal->vGetMosaicity(a3fTemp);
                  poPredict->vSetCrysMosaicity(a3fTemp);
              }
              else if (m_nPrerefineFlag)
              {
                  // A refinement failed for some reason during pre-refinement
                  // so increase nPrerefineSave
                  
                  nPrerefineSave++;                  
              }
             
              
              if (m_nPrerefineFlag)
              {
                  if (!nStat)
                    m_nPrerefineFlag --;
                  if (m_nCurrentImgNum == nNumImages)
                      m_nPrerefineFlag = 0;
                  // If we are pre-refineing, then we used this first iteration to get data.  Now we are going back.
                  if (m_nPrerefineFlag == 0) 
                  {
		      nPrerefineLoops--;
		      if (0 < nPrerefineLoops)
			{
			  printf("\n==================="
                                 "\nPRE-REFINE AGAIN, restarting ..."
                                 "\n===================\n\n");
			  m_nPrerefineFlag = nPrerefineSave;
			}
                      else
			printf("\n==================="
                          "\nPRE-REFINE FINISHED, restarting ..."
                          "\n===================\n\n");
                      

                      // Take the average (or maximum) mosaicity for the preceeding refinements.
                      if ((m_nNumRefineResults) && (m_fFixedMosaicity <= 0.0) && (m_bDetermineFixedMosaicity)) {
                          int nx;
                          m_fFixedMosaicity = 0.0;
                          for (nx = 0; nx < m_nNumRefineResults; nx++) {
                              m_fFixedMosaicity += m_ptRefineResults[nx].fMos;
                          };
                          m_fFixedMosaicity /= m_nNumRefineResults;
                          vSetFixedMosaicity(m_fFixedMosaicity,poPredict);
                      };
                      
                      goto pre_refine_rollback; 
                      
                  }
                }
                nStat = 0;
            }

            ANLSTOP(20);
            ANLSTOP(13);

//+JWP 2007-09-08
	    // Call dtscaleaverage if needed
	    if ("" != sGetEnv("DTINT_DTSCALE"))
	      {
		Cstring sScale;
		sScale = sGetEnv("DTINT_DTSCALE");
		cout << "INFO: trying DTINT_DTSCALE:\n"
                     << sScale << endl << flush;
		nDoSystemCommand(sScale + '&');
	      }

//-JWP 2007-09-08

        } // end of special end-of-batch processing
         
      ANLSTOP(21);
    }  // while loop

    ANLSTART(17);
    fclose(m_pFReflnFileOut);  // Close reflection file
    m_pFReflnFileOut = NULL;

    ANLSTOP(17);

    
    if (NULL != m_pFReflnBadOut)
    {
        fclose(m_pFReflnBadOut);  // Close reflection file
        m_pFReflnBadOut = NULL;
    }
    
    if (NULL != m_pFReflnPartialOut) 
    {
        fclose(m_pFReflnPartialOut);
        m_pFReflnPartialOut = NULL;
    }


  // Clean-up allocated objects and memory
  
  //Clean_up: ;
    nDeleteMarkedReflns(poReflnlistL,-999);

  // Insert header into dtintegrate reflection list.
  {
    Ccrystal oCrystal(*m_poHeader);
    int nStat;
    if (oCrystal.bIsAvailable()) 
      {
        // tjn:  Do not do this!  It causes the memory to get exhausted for large files.
        // nStat = poReflnlistL->nUpdateCrystal(m_sReflnFilename,&oCrystal);
        // nStat += poReflnlistL->nUpdateCrystal(m_sReflnRejectsFilename,&oCrystal);
        nStat = 0;
      } 
    else
      nStat = 1;    
    
    if (nStat)
      {
        printf("ERROR: Could not insert crystal object into reflection list '%s'\n",
               m_sReflnFilename.string());
      }
  }

  // Do final profile fitting if we are setup to do so.

  if ((m_poProfit2DPerSlice) && (m_bWriteScanBitmap)) {
      sTemp = sTransSymbol(sDtrekGetPrefix() + "dtintegrate.dat");
          m_poProfit2DPerSlice->nWriteScanBitmap(sTemp,m_poHeader);
  };
  
  //+jwp 22-Apr-2003
  // Update the display one last time with the predicted list, 
  // so that users are not confused by the display of the reflns used
  // in the last refinement.
#ifndef NO_X_WINDOWS
#ifndef SSI_PC
  if (NULL != m_poXprop)
    {
      // Write out only the newly predicted reflections and send update
      nWait(1000);  // We have to wait a little bit, say 1 second
      m_poXprop->hSetProperty("DTDISPLAY_REFLN_UPDATE",
			      m_sXpropString);
    }
#endif
#endif
  //-jwp

  
  //+jwp
  // If there was a read image failure, we still want to profile fit what
  // we have already done.
  if ((nStat != 0) && (1 == nReadStat))
    {
      nStat = 0;  // Make sure we can enter the profile-fitting block below
    }
  //-jwp
  
  //  if ((!nStat) && ((DTI_DO_2D_PROFIT & m_nProfitFlag))) {
  // Do the following to create a dtprofit.ref even if profile fitting not
  // called for.

    {
      Cprofit2DSlices oProfit;
      nStat = oProfit.nLoadHeader(*m_poHeader);
      
      int nFileParts = oProfit.nSplitFileProcessSize(*m_poHeader,g_nMaxReflnsToProfit);
      int nFilePart;
      int nFirstPart,nLastPart;
      
      
      nFirstPart = min(nFileParts,0);
      nLastPart = max(0,nFileParts) - 1;
      for (nFilePart = nFirstPart; nFilePart <= nLastPart; nFilePart++) {
          if (!nStat)
              nStat = oProfit.nLoad(m_sReflnFilename,m_sReflnPartialsFilename,(sTemp = ""),nFilePart);
          
          oProfit.vSetIntegrateOutputName(m_sReflnProfitFilename);
          
          oProfit.m_fIntersectMax = m_fMaxIntersectionFraction;
          
          if (!nStat)
              nStat = oProfit.nSyncLists(false);
          if (!nStat)
              nStat = oProfit. nPrintAxial(nFilePart==nLastPart);
          if (!nStat)
              nStat = oProfit.nRefineMosaicityModel();
          if (!nStat)
              nStat = oProfit.nPrintMosaicityModel(nFilePart);
          if (!nStat)
              nStat = oProfit.nFindIntersections();
          if (!nStat)
              nStat = oProfit.nCalcRemoveNoise(nFilePart == nFirstPart,nFilePart == nLastPart);
          if (!nStat)
              nStat = oProfit.nWrite(nFilePart,true,false,m_bWriteEllipsoids);  
      };
      if (!nStat) 
          oProfit.nWriteParts(true,false,m_bWriteEllipsoids);
      if (nStat) {
          printf("ERROR:  During applications of profiles to integrated list.\n");
      };
      for (i=0;i<oProfit.m_nTotalIntersections;i++)
          nSetErrorMask(Crefln::ms_nErrorOverlap, TRUE, NULL);
      
  };

  ANLSTART(18);

  vListIntensityProfiles();
  vListResults();
  ANLSTOP(18);

  nStat = 0;
  if (1 == nReadStat)
    {
      cout << "\nWARNING: failed to read a requested image (see above), \n"
           << "         so perhaps not all the specified images were processed.\n\n" 
           << flush;
    }

  delete poReflnlistL;
  delete poRefine;
  delete poDetector;
  delete poSource;
  delete poCrystal;
  delete poCrysGonio;
  delete poPredict;
#ifdef ANL_TIMER
  anl_print_timer(0);
  anl_print_timer(1);
  anl_print_timer(2);
  anl_print_timer(3);
  anl_print_timer(4);
  anl_print_timer(6);
  anl_print_timer(7);
  anl_print_timer(8);
  anl_print_timer(9);
  anl_print_timer(10);
  anl_print_timer(11);
  anl_print_timer(12);
  anl_print_timer(13);
  anl_print_timer(14);
  anl_print_timer(15);
  anl_print_timer(16);
  anl_print_timer(17);
  anl_print_timer(18);
  anl_print_timer(20);
  anl_print_timer(21);
  anl_print_timer(22);
#endif

  return (nStat);
}

int
Cintegrate::nList(void)
{
  int i;
  cout
       << "\n==============================================="
       << "\nIntegrate object listing:"
       << "\n==============================================="
       << "\nVerbose level:                    " << m_nVerbose
       << "\nScan sequence range:              " << m_anSeqNum[0]
       << ", " << m_anSeqNum[1]
         ;
  if (-999.0 != m_fOverallRotStart)
    {
      // Print out overall rotation start, end only if known.

      cout << "\nOverall scan rotation start, end: " << m_fOverallRotStart
           <<      ", " << m_fOverallRotEnd;
    }
  cout << "\nResolution range:                 " << m_a2fResolution[0]
       <<      ", " << m_a2fResolution[1]
       << "\nImages/ Scale batch:              " << m_nImagesPerScaleBatch
       << "\nImages/ Refine batch:             " << m_nImagesPerRefineBatch
       << "\nBatch prefix:                     " << m_sBatchPrefix
       << "\nImage padding:                    " << m_nPad3D
       << "\nImage padding multiplier:         " << m_fPad3D
       << "\nMinimum peak radius in pixels:    " << m_nMinPeakRad       
       << "\nSpot size multiplier:             " << m_fSpotSizeMultiplier
       << "\nInput max window size:            " << m_a2nWinMax[0] << ", " << m_a2nWinMax[1]
       << "\nInput min window size:            " << m_a2nWinMin[0] << ", " << m_a2nWinMin[1]
       << "\nMaximum refln overlap fraction:   " << m_fMaxIntersectionFraction
       << "\nAlpha1 / Alpha2 splitting option: ";
  if (m_bSpotChase)
    cout << "On";
  else
    cout << "Off";
  if (m_fFixedMosaicity > 0.0)
    cout << "\nMosaicity fixed to value:         " << m_fFixedMosaicity;
  else if (m_a2fMosaicityModel[0] == 0.0)
    cout << "\nMosaicity fixed to value:         " << m_a2fMosaicityModel[1];
  else
    cout << "\nMosaicity model:                  " << m_a2fMosaicityModel[0] << ", " 
                                                   << m_a2fMosaicityModel[1]
         << "\n  mosaicity used = (Refined_mosaicity * " << m_a2fMosaicityModel[0]
         << ") + " << m_a2fMosaicityModel[1];
  cout << "\nNum oblique incidence factors:    " << m_nNumObliqueCorr;
      for (i = 0; i < m_nNumObliqueCorr; i++)
        {
          cout << "\n    Oblique incidence factor " << i+1 << ":   "
               << m_pfObliqueCorr[i] <<  "   Norm: "
               << m_pfObliqueCorrNorm[i]
               << "   Type: "
               << m_pfObliqueCorrType[i];
        }
  cout << "\nNum excluded resolution rings:    " << m_nNumExcludeResoRings;
      for (i = 0; i < m_nNumExcludeResoRings; i++)
        {
          cout << "\n    Ring " << i+1 << " resolution min,max:    "
               << m_pfExcludeReso[0][i] <<  ", "
               << m_pfExcludeReso[1][i];
        }
  cout
       << "\n===============================================\n"
       << endl << flush;

  return (0);
}


int
Cintegrate::nExpandReflnlist(Creflnlist *poReflnlistIn)
{

// Add required fields to the reflection list if they are not already there.
// The fields are named by the static member variables
// (besides H, K, L, Intensity, SigmaI)
//
// Creflnlist::ms_snPackedHKL
// ssnDataFlag

  Creflnlist *poReflnlistL;   // Local reflection list

  if (NULL == poReflnlistIn)
    {
      poReflnlistL = m_poReflnlist;     // Use on member variable
    }
  else
    {
      poReflnlistL = poReflnlistIn;   // Use passed reflection list
    }

  // This field name is a member of Creflnlist class.

  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_snPackedHKL);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPx0);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPx1);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsRotMid);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsRotWidth);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsRotSigma);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfBackground);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfBackgroundSigma);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_ssBatch);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfProfitIntensity);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfProfitSigmaI);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_snPartialFlag);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_snRefineBatch);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfSvec0);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfSvec1);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfSvec2);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfS0vec0);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfS0vec1);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfS0vec2);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidA00);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidA01);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidA11);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidb0);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidb1);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidc);

  // These field names are local to Cintegrate class

  m_nFI_nDataFlag  = poReflnlistL->nExpandGetField(ms_snDataFlag);
  m_nFI_nErrorFlag = poReflnlistL->nExpandGetField(ms_snErrorFlag);
  m_nFI_fWidth0    = poReflnlistL->nExpandGetField(ms_sfWidth0);
  m_nFI_fWidth1    = poReflnlistL->nExpandGetField(ms_sfWidth1);


  if (0 > m_nFI_nDataFlag)
    return (m_nFI_nDataFlag);
  return (0);
}
int
Cintegrate::nBuildShoeboxes(Creflnlist *poReflnlist,
                            Cimage *poImageIn,
                            Crotation *poRotationIn,
                            Cdetector *poDetector
                            )
{
  // During integration, process reflections
  // The Crefln field m_nFI_nDataFlag is used as a flag
  // to designate the status of a reflection and its processing.
  // If this field is equal to:
  //    ms_nNewFlag             This is a brand new reflection with no shoebox
  //    ms_nActiveFlag          This is an active reflection with a shoebox
  //    ms_nDelFlag             This means the reflection should be deleted;
  //                            shoebox data if any existed was already deleted


  int i, j;                      // Loop counters
  int nStatLocal;                // Local status flag
  int nStat;                     // Returned status flag
  int nNumReflns;                // Variable to hold number of reflns in list
  int nPack1, nPack2;            // Packed hkl's for 2 reflns
  int nField1, nField2;          // Field indexes for 2 reflns
  int nToDelete;                 // Counter for number of reflns to delete
  int nFull;                     // Counter for full reflections
  int nNew;                      // Counter for new reflections
  int nActive;                   // Counter for active (new & continuing) reflns
  int nDone;                     // The reflection was integrated.
  int nProblem;                  // The reflection had a problem integrating.

  Crefln  *poRefln1;             // Pointer to a reflection
  Crefln  *poRefln2;             // Pointer to next reflection
  C3Ddata *poShoebox;            // Pointer to a shoebox


  // Initialize counters

  nNumReflns = poReflnlist->nGetNumReflns();
  nToDelete  = 0;
  nFull      = 0;
  nNew       = 0;
  nActive    = 0;
  nStat      = 0;
  nStatLocal = 0;
  nDone      = 0;
  nProblem   = 0;



  // Now loop through rest of reflections

  if (2 < m_nVerbose)
    {
      cout << "Processing reflns ...\n"
           << "  there are " << nNumReflns << " reflns in the list." << endl << flush;
    }

  for (i = 0; (i < nNumReflns) && (0 == nStat); i++)
  {


      poRefln1   = poReflnlist->poGetRefln(m_pnHKLIndex[i]);
      nPack1     = poRefln1->nGetField(poReflnlist->m_nFI_nPackedHKL);
      nField1    = poRefln1->nGetField(m_nFI_nDataFlag);

      
      if (ms_nDelFlag == nField1)
      {
          // Do nothing this reflection was already marked for deletion
      }
      else
      {
          // This reflection not marked for deletion
          // So look to see if the next reflection is a duplicate
          
          if (i+1 < nNumReflns)
          {
              // Peek ahead at the next not deleted reflection...
              
              j = i+1;
              do
              {
                  poRefln2 = poReflnlist->poGetRefln(m_pnHKLIndex[j]);
                  nPack2   = poRefln2->nGetField(poReflnlist->m_nFI_nPackedHKL);
                  nField2  = poRefln2->nGetField(m_nFI_nDataFlag);
                  j++;
              }
              while ( (j < nNumReflns) && (ms_nDelFlag == nField2) );
              j--;
              
              // If it has same index and it has not been marked for deletion...
              
              if ( (nPack1 == nPack2) && (ms_nDelFlag != nField2) )
              {
                  // ... then it is a duplicate, so check sorting order...
                  


                  m_pfNearestCentroidBuf0[min(m_pnHKLIndex[j],m_pnHKLIndex[i])] = m_pfNearestCentroidBuf0[max(m_pnHKLIndex[j],m_pnHKLIndex[i])];
                  m_pfNearestCentroidBuf1[min(m_pnHKLIndex[j],m_pnHKLIndex[i])] = m_pfNearestCentroidBuf1[max(m_pnHKLIndex[j],m_pnHKLIndex[i])];
                  
                  if (m_pnHKLIndex[i] > m_pnHKLIndex[j])
                  {
                      // Keep the newer predicted reflection in poRefln2
                      // and make sure order in poReflnlist is correct
                      
                      Crefln *poTemp = poReflnlist->poGetRefln(m_pnHKLIndex[i]);
                      poReflnlist->m_ppoTheReflns[m_pnHKLIndex[i]]
                          = poReflnlist->poGetRefln(m_pnHKLIndex[j]);
                      poReflnlist->m_ppoTheReflns[m_pnHKLIndex[j]] = poTemp;
                      
                      poRefln1       = poReflnlist->poGetRefln(m_pnHKLIndex[i]);
                      poRefln2       = poReflnlist->poGetRefln(m_pnHKLIndex[j]);
                      nField1        = poRefln1->nGetField(m_nFI_nDataFlag);
                      nField2        = poRefln2->nGetField(m_nFI_nDataFlag);
                  }
                  
                  // Duplicates should always be new and the first of them
                  // should not be new... so it is an error if they are not
                  
                  if ( (ms_nNewFlag == nField1) || (ms_nNewFlag != nField2) )
                  {
                      // Reflection was predicted twice.
                      // Probably is high-lorentz.
                      // To be safe, we will reject both of these.
                      poRefln1->vSetField(m_nFI_nDataFlag, ms_nDelFlag);
                      poRefln2->vSetField(m_nFI_nDataFlag, ms_nDelFlag);
                      nToDelete++;
                      nToDelete++;
                      continue;
                  }
                  
                  // The newer reflection should have the more accurate
                  // calculated stuff, so keep the new one and
                  // transfer the shoebox pointer of the old one to the
                  // new one.  Also transfer any values calculated
                  // in any first pass integration, such as
                  // observed centroid, background avg and sigma,
                  // peak intensity and sigma
                  
                  // TODO: if the newly predicted centroid is too far from
                  //       the previous, then reject this reflection!
                  
                  poRefln2->vSetField(m_nFI_nDataFlag, nField1);
                  poRefln2->vSetField(m_nFI_nErrorFlag,
                      poRefln1->nGetField(m_nFI_nErrorFlag));
                  
                  // It is important to mark a refln as overlapped if it was
                  // overlapped in ANY of the images it was predicted to be
                  // on, so OR together the Nonunf flag field that has overlap
                  // info
                  
                  poRefln2->vSetField(poReflnlist->m_nFI_nNonunfFlag,
                      poRefln1->nGetField(poReflnlist->m_nFI_nNonunfFlag)
                      | poRefln2->nGetField(poReflnlist->m_nFI_nNonunfFlag)  );
                  
                  poRefln2->vSetField(poReflnlist->m_nFI_fObsPx0,
                      poRefln1->fGetField(poReflnlist->m_nFI_fObsPx0));
                  poRefln2->vSetField(poReflnlist->m_nFI_fObsPx1,
                      poRefln1->fGetField(poReflnlist->m_nFI_fObsPx1));
                  poRefln2->vSetField(poReflnlist->m_nFI_fObsRotMid,
                      poRefln1->fGetField(poReflnlist->m_nFI_fObsRotMid));
                  poRefln2->vSetField(poReflnlist->m_nFI_fObsRotWidth,
                      poRefln1->fGetField(poReflnlist->m_nFI_fObsRotWidth));
                  poRefln2->vSetField(poReflnlist->m_nFI_fObsRotSigma,
                      poRefln1->fGetField(poReflnlist->m_nFI_fObsRotSigma));
                  
                  poRefln2->vSetIntensity(poRefln1->fGetIntensity());
                  poRefln2->vSetSigmaI(poRefln1->fGetSigmaI());
                  
                  poRefln2->vSetField(poReflnlist->m_nFI_fBackground,
                      poRefln1->fGetField(poReflnlist->m_nFI_fBackground));
                  poRefln2->vSetField(poReflnlist->m_nFI_fBackgroundSigma,
                      poRefln1->fGetField(poReflnlist->m_nFI_fBackgroundSigma));


                  if ((ABS(poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotStart) - m_fOverallRotStart)< 0.01) ||
                      (ABS(poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotEnd) - m_fOverallRotEnd)< 0.01)) {
                      // This was probably a partial which had it's obs rot start and end adjusted.
                      poRefln2->vSetField(poReflnlist->m_nFI_fCalcRotStart,
                          poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotStart));
                      poRefln2->vSetField(poReflnlist->m_nFI_fCalcRotEnd,
                          poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotEnd));
                      poRefln2->vSetField(poReflnlist->m_nFI_fCalcRotMid,
                          poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotMid));
                      poRefln2->vSetField(poReflnlist->m_nFI_fCalcRotWidth,
                          poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotWidth));
                  };

                  
                  if (NULL != poRefln2->m_pvUserData)
                  {
                      cout << "WARNING! Memory leak in nBuildShoebox!\n";
                  }
                  poRefln2->m_pvUserData = poRefln1->m_pvUserData;
                  poRefln1->m_pvUserData = NULL;

                  // If the old one was a partial, then we copy over the calcrot's as well.
                  
                  // Mark the old one for deletion
                  
                  poRefln1->vSetField(m_nFI_nDataFlag, ms_nDelFlag);
                  nField1 = ms_nDelFlag;
                  nToDelete++;
              }
            }  // end of peek ahead
            
            // Now process this reflection
            
            if (ms_nNewFlag == nField1)
            {
                // Request for new reflection, 
                
                poRefln1->vSetField(m_nFI_nErrorFlag, 0);
                m_nNumPredictedUnique++;
                
                if ((poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotStart)
                    < poRotationIn->fGetRotStart()) && (ABS(poRotationIn->fGetRotStart() - m_fOverallRotStart)>0.001))
                {
                    // Reflection started before this image, so probably
                    // already in the list OR a problem with the prediction
                    // algorithm OR scan partial on first image
                    
                    // Partial only at beginning of image, but not scan partial
                    // This happens when you backup for prediction.
                    
                    m_nNumPredictedUnique--;  // Not UNIQUE, so uncount it
                    nField1 = ms_nDelFlag;
                    poRefln1->vSetField(m_nFI_nDataFlag, ms_nDelFlag);
                    nToDelete++;
                } else {
                    // New reflection, create shoebox for it, and save the
                    // shoebox index.  If this reflection fits on only this
                    // image, it may be new and done at the same time.
                    
                    nStat  = nNewShoebox(poRefln1, poRotationIn);
                    
                    if (0 == nStat)
                    {
                        nField1 = ms_nActiveFlag;
                        poRefln1->vSetField(m_nFI_nDataFlag, ms_nActiveFlag);
                        nNew++;
                    }
                    else
                    {
                        // Some error in the shoebox dimensions.  
                        
                        // Was rejected just because it was a partial?
                        if ((nStat & (1 << Crefln::ms_nErrorPartialStart)) || 
                            (nStat & (1 << Crefln::ms_nErrorPartialEnd))) {
                            
                            
                            if (nStat & (1 << Crefln::ms_nErrorPartialStart))
                                m_nNumRefPartialStart++;
                            else if (nStat & (1 << Crefln::ms_nErrorPartialEnd))
                                m_nNumRefPartialEnd++;
                            
                            
                            // Keep track of the reflection.
                            poRefln1->vSetField(m_nFI_nDataFlag,  ms_nDelFlag);
                            poRefln1->vSetField(m_nFI_nErrorFlag,nStat);
                            nStat = 0;
                            
                            if (8 < m_nVerbose) {
                                
                                cout << "Scan partial:  "
                                    << poRefln1->nGetH() << ", "
                                    << poRefln1->nGetK() << ", "
                                    << poRefln1->nGetL()
                                    << ", Start, End: "
                                    << poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotStart)
                                    << " -> "
                                    << poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotEnd)
                                    << " Width: " << setprecision(4)
                                    << poRefln1->fGetField(poReflnlist->m_nFI_fCalcRotWidth)
                                    << '\n' << flush;
                            }
                        } else {
                            // Some other error in the shoebox dimensions.  
                            poRefln1->vSetField(m_nFI_nDataFlag,  ms_nDelFlag);
                            nToDelete++;
                            nStat = 0;
                        };
                    };
                }
            }
        }  // end not deleted reflection
    } // end of for (i...

  // Second loop to fill shoeboxes.
  // First sort reflnlist on the shoebox data address
  // Copy the image data to the shoebox if needed...

  vInitBackgroundMask();

  for (i = 0; (i < nNumReflns) && (0 == nStat); i++)
  {


      poRefln1   = poReflnlist->poGetRefln(i);
      nField1    = poRefln1->nGetField(m_nFI_nDataFlag);
      
      if (ms_nActiveFlag == nField1)
      {
          // Active reflection with shoebox allocated.
          // Reflection may be new, continuing and done at the same time.
          
          nActive++;
          
          // Debug stuff
          int nH = poRefln1->nGetH();
          int nK = poRefln1->nGetK();
          int nL = poRefln1->nGetL();

          
          // Fill the next layer, unless the shoebox is full
          
          poShoebox = (C3Ddata*)poRefln1->m_pvUserData;
          vAddBackgroundMask(poShoebox);
          if (!poShoebox->bIsFilled())
          {
              // Fill information about competing spots for each slice.              
              poShoebox->vSetCompetingSpot(poShoebox->m_nFilled[2],
                  m_pfNearestCentroidBuf0[i],
                  m_pfNearestCentroidBuf1[i]);

              nStat = poShoebox->nFill2D(poImageIn, -1);
              
              
              if (poShoebox->bIsFilled())
              {
                  nFull++;                  
                  nStatLocal = nPass0Shoebox(poRefln1, poDetector);

                  if (nStatLocal == 0)
                      nDone++;
                  else
                      nProblem++;
                  
                  poRefln1->vSetField(m_nFI_nDataFlag, ms_nDelFlag);
                  vDeleteShoebox(poShoebox);          
                  poRefln1->m_pvUserData = NULL;
                  
              }
          }
      } else if (ms_nDelFlag == nField1) {
          // Do Nothing, this refln was probably deleted because it was
          // partial or duplicate
      }
      else
      {
          // Should never get here, but have error message just in case;
          // unless poRefln1 was deleted above. (Could be leftover full?)
          
          cout << "ERROR: unknown type: "
              << poRefln1->nGetH() << ", "
              << poRefln1->nGetK() << ", "
              << poRefln1->nGetL() << ", "
              << nField1 << '\n';
      }
  }  // end shoebox filling i loop


//jwp
  printf("\nReflection dispositions\n"
         "=================================================\n");
  printf("Status:    New  Active    Full   DNormal Dspecial\n"
         "Number:%7d %7d %7d %7d %7d\n"
         "=================================================\n",
         nNew, nActive, nFull, nDone, nProblem);

  nStat = nDeleteMarkedReflns(poReflnlist, ms_nDelFlag);

  if (m_a2oDetectorProfiles[m_nProfileLoading].fGetAverageLevel()>=m_fTargetDetectorProfilesLevel) {
      m_fTargetDetectorProfilesLevel = 3.5;
      m_a2oDetectorProfiles[m_nProfileLoading].nPrintPeakAreas(false,m_fRotIncrement);
      m_nProfileReading = m_nProfileLoading;
      m_nProfileLoading = ! m_nProfileReading;
      m_a2oDetectorProfiles[m_nProfileLoading].nInitPeakAreas(m_nDim0,m_nDim1);
  };

  if ((nDone>10) && (m_nDoIntegrateCount> m_nProfitMaxImageRange)) {     
      m_nDoIntegrateCount = 0;
  };
  m_nDoIntegrateCount++;

  return (nStat);
}

// This routine calculates the smallest shoebox dimensions that could hold the reflection.

int
Cintegrate::nCalcSmallestShoeboxRange(int nExt2,int nCenterSlice,
                                      int nCent0,int nCent1,
                                      int& nStart0,int& nExt0,int& nStart1,int& nExt1)
{
    
    int a2nCent[2];
    int nMax0,nMax1;
    int nEnd0,nEnd1;
    C3DdataDetectorArea oArea;
    
    // Loop through each slice and increment the observed shifts from the center.
    nEnd0 = nStart0 = nEnd1 = nStart1 = 0;
    
    a2nCent[0] = 0;
    a2nCent[1] = 0;
    nStart0 = nEnd0 = a2nCent[0];
    nStart1 = nEnd1 = a2nCent[1];
    nEnd0 = max(nEnd0,a2nCent[0]+ ((int) 50));
    nStart0 = min(nStart0,a2nCent[0]- ((int) 50));
    nEnd1 = max(nEnd1,a2nCent[1]+ ((int) 50));
    nStart1 = min(nStart1,a2nCent[1]- ((int) 50));
    
    nExt0  = nEnd0 - nStart0 + 1;
    nExt1  = nEnd1 - nStart1 + 1;
    
    nStart0 += nCent0;
    nStart1 += nCent1;
    
    // Maximum limits on spot size calculated in nMax0,nMax1.
    
    m_a2oDetectorProfiles[m_nProfileReading].nGetPeakArea(nStart0,nStart1,oArea);
    nMax0 = (int) (oArea.fSize[0]*m_fSpotSizeMultiplier);
    nMax1 = (int) (oArea.fSize[1]*m_fSpotSizeMultiplier);
    if (m_a2nWinMax[0]>0) {
        nMax0 = min(nMax0,m_a2nWinMax[0]);
        nMax1 = min(nMax1,m_a2nWinMax[1]);
    };
    if (m_a2nWinMin[0]>0) {
        nMax0 = max(nMax0,m_a2nWinMax[0]);
        nMax1 = max(nMax1,m_a2nWinMax[1]);
    };
    
    if (nExt0 > nMax0)
    {
        nStart0 = nCent0 - nMax0/2;
        nExt0   = nMax0;
    }
    if (nExt1 > nMax1)
    {
        nStart1 = nCent1 - nMax1/2;
        nExt1   = nMax1;
    }
    
    if (nExt0 + nStart0 > m_nDim0) 
        nExt0 = m_nDim0 - nStart0;
    if (nExt1 + nStart1 > m_nDim1)
        nExt1 = m_nDim1 - nStart1;
    if (nStart0 < 0) {
        nExt0 += nStart0;
        nStart0 = 0;
    };
    if (nStart1 < 0) {
        nExt1 += nStart1;
        nStart1 = 0;
    };
    
    return 0;
};

int
Cintegrate::nNewShoebox(Crefln *poRefln, Crotation *poRotation)
{
  // Allocate a new shoebox for the input reflection poRefln.
  // poRotation contains the rotation object for the first image of the
  // shoebox.
  // Return  = 0 for success

  int     nStat;      
  int     nx;//,ny;
  double  f0,f1;

  nStat  = 0;

  int   nExt0, nExt1, nExt2;
  int   nStart0, nStart1, nStart2;
  int   nMaxImages;

  float fTemp;
  float fRotMid;
  float fRotWidth;
  bool bPartial;        // A reflection that is a partial at the begin or end of scan.
  bool bPadPartial;     // Even worse, a reflection that is only a partial (possibly) if pad is added.
  int  nPartialStat;
  
  bPartial = false;
  bPadPartial = false;
  nPartialStat = 0;
  nx=0;
  

  if (  poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid) -
      (poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth)*0.5*m_fPad3D)
      < m_fOverallRotStart) {
      bPartial = true;
      if (m_fFixedMosaicity == 0.0) {
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotStart,(float) (m_fOverallRotStart + 0.001));
          // Special case:  Don't reject partials.  This is a user flag.
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotEnd,(float) (m_fOverallRotStart + 0.02));
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid,(float) (m_fOverallRotStart + 0.01));
      } else {
        nPartialStat |=  nSetErrorMask(Crefln::ms_nErrorPartialStart, TRUE, NULL);
        nx++;
      };
  };
  
  if (  poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid) +
      (poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth)*0.5*m_fPad3D)
      > m_fOverallRotEnd) {
      bPartial = true;
      
      if (m_fFixedMosaicity == 0.0) {
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotStart,(float) (m_fOverallRotEnd - 0.001));
          // Special case:  Don't reject partials.  This is a user flag.
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotEnd,(float) (m_fOverallRotEnd - 0.02));
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid,(float) (m_fOverallRotEnd - 0.01));
      } else {
          nPartialStat |=  nSetErrorMask(Crefln::ms_nErrorPartialEnd, TRUE, NULL);
          nx++;
      };
  };
  
  fRotMid = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);
  fRotWidth = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth);

  // Adjust the Mid and Width rotation parameters to agree with the Start and End rotation parameters.
  f0 = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotStart);
  f1 = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotEnd);
     
  poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid,(float)((f1+f0)*0.5));
  poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth,(float)(f1-f0));

  if (!bPadPartial) {
      // In case bPadPartial == true: 
      // These reflections were completely off the image.
      // However, we might still want to use them if the padding intersects with this image.
      // We already set their Start and End values to be in the image range.
      // For our image overlap prediction however, we want to use the original calculated values.
      fRotMid = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);
      fRotWidth = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth);
  };


  f0 =  fRotMid - fRotWidth*0.5*m_fPad3D -  m_fRotIncrement * m_nPad3D;
  f1 =  fRotMid + fRotWidth*0.5*m_fPad3D +  m_fRotIncrement * m_nPad3D;

  nMaxImages = (int) (0.5 + (m_fOverallRotEnd - m_fOverallRotStart)
               / m_fRotIncrement);

  nStart2 = (int) floor((f0 - m_fOverallRotStart)/ m_fRotIncrement);
  nExt2 = (int) floor((f1 -  m_fOverallRotStart)/ m_fRotIncrement);
  nExt2 -=nStart2;
  nExt2++;


  if ((nStart2+nExt2 -1 < m_nCurrentImgNum) || (nStart2 >= nMaxImages)) {
      // The reflection is just plumb off the image!  There is nothing we can do.
      return Crefln::ms_nErrorIntProblem;
  };

  // Otherwise, we should be able to restrict the ranges.

  if (nStart2 < m_nCurrentImgNum) {
      nExt2 += (nStart2 - m_nCurrentImgNum);
      nStart2 = m_nCurrentImgNum;
  };
  if (nStart2 + nExt2 - 1 >= nMaxImages) {
      nExt2 = nMaxImages - nStart2;
  };
           

  // Find the centroid.
  int   nCenterSlice;

  // Calculate the centroid slice.
  fTemp         = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);
  nCenterSlice  = (int) max(0.0,(fTemp - m_fOverallRotStart)/m_fRotIncrement - nStart2); 

  // Get the ellipsoid.
  nStat = nCalcSmallestShoeboxRange(
      nExt2,nCenterSlice,
      (int) poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx0),
      (int) poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx1),
      nStart0,nExt0,nStart1,nExt1);
      
  if (nStat) {
      // Serious error.  We should NEVER have this with a correctly modeled ellipsoid.
      printf("Serious error in nNewShoebox().\n");
      return Crefln::ms_nErrorIntProblem;
  };



  if (m_nShoeboxLimit < nExt2)
    {
      // Reject TOO WIDE in rotation angle reflections silently.
      return Crefln::ms_nErrorRotTooWide;
    }

  

  C3Ddata *poShoebox;


  poShoebox   = poGetShoebox(nExt0, nExt1, nExt2, nStart0, nStart1, nStart2);
  poRefln->m_pvUserData = poShoebox;

  // Save the shoebox's data address as a pseudo address for sorting on later

  unsigned long int ulDataAddr;
  double   dAddr;
  ulDataAddr = (unsigned long int) poShoebox;
  dAddr = (double) ulDataAddr;
  // poRefln->vSetField(m_nFI_fDataAddr, (float)dAddr);

  return (nPartialStat);
}

int
Cintegrate::nPass0Shoebox(Crefln *poRefln, Cdetector *poDetector)
{

  int nStat;
  nStat = 0;

  C3Ddata* poShoebox = (C3Ddata*) poRefln->m_pvUserData;
  
  if (5 < m_nVerbose)
    {
      cout << "Pass 0 reflection: "
           << poRefln->nGetH() << ", "
           << poRefln->nGetK() << ", "
           << poRefln->nGetL() << ", "
           << poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotStart) << ", "
           << poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotEnd) << ", "
           << poShoebox->nGetExt(2) << ", "
           << poShoebox->nGetOffset(2) << '\n' << flush;
    }

  // Convert shoebox to float, but put in new re-usable location, leaving
  // original shoebox untouched (and compressed!)


  int nExt, nFill;
  nExt = poShoebox->nGetExt(2);
  nFill = poShoebox->nGetLayersFilled();
  nStat = poShoebox->nConvertToFloat(NULL, m_po3DCentroid);
  poShoebox = m_po3DCentroid;

  if (0 != nStat)
    {
      poRefln->vSetField(m_nFI_nErrorFlag,
                         poRefln->nGetField(m_nFI_nErrorFlag)
                         | nSetErrorMask(Crefln::ms_nErrorNonunfB, TRUE, NULL));
      cout << "You are really in trouble!\n";
      return (ms_nDelFlag);
    }

  // Check nonunf flags for other types of errors.

  if (0 != (nStat = poRefln->nGetField(poRefln->m_poReflnlist->m_nFI_nNonunfFlag)))
    {
      // Go through each of the types of flag to see what the error is.
      if (nStat & m_nNonunfMaskRing) {
          poRefln->vSetField(m_nFI_nErrorFlag,
              poRefln->nGetField(m_nFI_nErrorFlag)
              | nSetErrorMask(Crefln::ms_nErrorRing,
              TRUE, NULL));
      } else if (nStat & m_nNonunfMaskBadPix) {
          poRefln->vSetField(m_nFI_nErrorFlag,
              poRefln->nGetField(m_nFI_nErrorFlag)
              | nSetErrorMask(Crefln::ms_nErrorNonunfB,
              TRUE, NULL));
      } else {
          cout << "ERROR:  Nonunf flag not recognized: " << nStat << "\n";
      };

      if (NULL != m_pFReflnBadOut)
        {
          // Write out really bad rejected refln to special file

          poRefln->vWrite(m_pFReflnBadOut, m_pcSelectFields);
        }

      // Count as bad since they won't be seen again

      m_nNumRefBadErrors++;
      return (ms_nDelFlag);
    }

  if (NULL != poDetector)
    {
      // Correct the shoebox for non-uniformity
      // if the detector object and its non-uniformity object is available

      // First, find out which image is the first image in the shoebox

      nStat = poDetector->m_poNonunf->nCorrect3Ddata(poShoebox,
                              &m_pfPerImageFactors[poShoebox->nGetOffset(2)]);

      // The above flags for both bad nonunf pixels and saturated pixels
      // But the saturated pixel is looking in the whole shoebox, not
      // just the peak area, so it may include a neighboring saturated pixel
      // We'll treat saturated pixels later, so undo here

      if (0 != (Cnonunf::ms_nSatPixFlag & nStat))     // Maybe a XOR would be better?
        {
          nStat = nStat - Cnonunf::ms_nSatPixFlag;      // Remove saturated pixel error
        }

      if (0 != (Cnonunf::ms_nBadPixFlag & nStat))
        {
          poRefln->vSetField(m_nFI_nErrorFlag,
                             poRefln->nGetField(m_nFI_nErrorFlag)
                             | nSetErrorMask(Crefln::ms_nErrorNonunfA, FALSE, NULL));

          nStat = nStat - Cnonunf::ms_nBadPixFlag;      // Remove bad pixel error
        }
    }
  else
    {
      // No detector object

      nStat = -1;
    }

  if (0 != nStat)
    {
      // There is a severe problem with this reflection so return with severe
      // error

      poRefln->vSetField(m_nFI_nDataFlag, ms_nDelFlag);
      m_nNumRefBadErrors++;
      if (NULL != m_pFReflnBadOut)
        {
          // Write out really bad rejected refln to special file

          poRefln->vWrite(m_pFReflnBadOut, m_pcSelectFields);
        }

      cout << "ERROR BEFORE COG CALCULATION: " << m_nNumRefBadErrors << endl;
      cout << "NONUNF NUM: " << m_anNumReflns[Crefln::ms_nErrorNonunfA] << endl << flush;
      return (ms_nDelFlag);
    }

  // Find center of gravity and intensity of shoebox peak in a special way


  nStat = nCalcMask(poRefln, poDetector, poShoebox);

  if (nStat == 0) {
      m_a2oDetectorProfiles[m_nProfileLoading].nAddPeakArea(
          (int) poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fObsPx0),
          (int) poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fObsPx1),
          poShoebox);
  };

  if (0 == (m_nErrorDoNotProcessMask & nStat)) {
      float a3fVDC[3];
      float a3fS[3];
      float a3fSrot[3];
      float a3fS0rot[3];
      float a3x3fRotMatHKL[3][3];
      float a3x3fTemp[3][3];
      float fRLP;
      int   nObl;
      int   j;
      float fIntensity,fIntensitySigma;
//mrp      char a8cBatch[8];
      char a100cBatch[100]; //mrp
	  memset(a100cBatch, '\0', sizeof(char)*100); //mrp
      Cstring sBatch;

      // Apply LP correction
      
      fRLP = 1.0f / (poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fLorentz)
          * poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fPolarz));
      
      // Apply any oblique incidence correction
      // See for example: Zaleski, Wu & Coppens (1998)
      //                     J. Appl. Cryst. 31, 302-304.
      // Icorr = Iobs * Kobl;
      // Tobl  = (1.0 - exp(Foblnorm)) / (1.0 - exp(Foblnorm/cosA))
      // The refln has 1/cosA in the fCalc_oblique field.
      // m_nNumObliqueCorr may equal 0.
      // Remember (1.0 - exp(Foblnorm)) = m_pfObliqueCorrNorm
      
      // !!!Note: The above is a TRANSMISSION effect.
      //          An ABSORPTION effect would be
      // Aobl = exp(Foblnorm) / exp(Foblnorm/cosA)
      //      = -exp(Foblnorm) / -exp(Foblnorm/cosA)
      
      m_fRefOblique = 1.0;
      for (nObl = 0; nObl < m_nNumObliqueCorr; nObl++)
      {
          m_fRefOblique *= (m_pfObliqueCorrNorm[nObl]
              / (m_pfObliqueCorrType[nObl]
              - (float)exp(m_pfObliqueCorr[nObl]
              * poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fOblique))));
      }
      
      fRLP *= m_fRefOblique;
      
      fIntensity = poShoebox->fGetPeakIntensity();
      fIntensitySigma = poShoebox->fGetPeakIntensitySigma();
      
      if (-999.0 != fIntensitySigma)
          poRefln->vSetSigmaI(fIntensitySigma * fRLP);
      if (-999.0 != fIntensity)
          poRefln->vSetIntensity(fIntensity * fRLP);
      
      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_nPackedHKL,0);
      
      // Figure out S and S0 vectors.
      
      if (poDetector->m_poSpatial->nPixeltoMM(
          poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx0),
          poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx1),
          &a3fVDC[0],&a3fVDC[1],&a3fVDC[2]))
          return 1 << Crefln::ms_nErrorIntProblem;
      
      vConvRotVec3DMat3D(-poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid), m_a3fRotVec, &a3x3fRotMatHKL[0][0]);
      //m_poScan->m_poRotation->vCalcGetRotMatrix(-poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid),
      //    &a3x3fRotMatHKL[0][0]);
      
      vMulMat3DVec3D(m_a3x3fDetDD, a3fVDC, a3fS); // Scattered beam wavevector
      (void) fNormVec3D(a3fS);                    // Normalize sc. beam wavevec.
      vMulVec3DScalar(a3fS, m_fLenS0, a3fS);      // to be same length as S0
      vMulMat3DMat3D(m_a3x3fInvCrysRotGonMat, a3x3fRotMatHKL, a3x3fTemp);
      vMulMat3DVec3D(a3x3fTemp, a3fS,  a3fSrot);
      vMulMat3DVec3D(a3x3fTemp, m_a3fS0, a3fS0rot);
      for (j = 0; j < 3; j++)
      {
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fSvec[j], a3fSrot[j]);
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fS0vec[j], a3fS0rot[j]);
      }
      
      
      
      // Set batch id for this reflection
      // Set batch id based on calculated rotation angle,
      // Unfortunately this bit does not know about minimum
      // number of reflns, average number of reflns, etc 
      // like dtrebatch
      
      float fTemp   = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);
      int nBatchNum;
      fTemp     = ((fTemp - m_fOverallRotStart) / m_fRotIncrement);
      nBatchNum = (int)(fTemp / (float)m_nImagesPerScaleBatch) + 1;
      
      nBatchNum = (nBatchNum * m_nImagesPerScaleBatch) 
          - 1 + m_anSeqNum[0];
      
//mrp      if (10000 > nBatchNum)
      nBatchNum %= 10000; //mrp
//mrp          sprintf(a8cBatch, "%s%04d", m_sBatchPrefix.string(), nBatchNum);
      sprintf(a100cBatch, "%s%04d", m_sBatchPrefix.string(), nBatchNum); //mrp
//mrp      else
//mrp          sprintf(a8cBatch, "%s%d", m_sBatchPrefix.string(), nBatchNum);
//mrp          sprintf(a100cBatch, "%s%d", m_sBatchPrefix.string(), nBatchNum); //mrp
//mrp      a8cBatch[7] = '\0';
//mrp      sBatch      = (Cstring) a8cBatch;
      sBatch      = (Cstring) a100cBatch; //mrp
      
      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_sBatch, sBatch);
               
      
      if ((nStat==0) && (m_poReflnlistRefine->nGetNumReflns()<m_nRefineMaxRefln))
      {
          m_poReflnlistRefine->nInsert(poRefln);
      };


      // Write out refln
      
      poRefln->vWrite(m_pFReflnFileOut, m_pcSelectFields);
      m_nReflnFileOutCount++;
      if (ferror(m_pFReflnFileOut))
          return 1;

          

      // Write out partials.
      int nImageStart,nImageEnd,nImage;      
      double fObsRotStart,fObsRotEnd;
      double fCalcRotMid;
      float  a2x2fEllipsoidA[2][2];
      float  a2fEllipsoidb[2];
      float  fEllipsoidc;
      float  a2fObsCentroidWritten[2];   // The centroid written out might not be the one returned by the shoebox routine.

      a2fObsCentroidWritten[0] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fObsPx0);
      a2fObsCentroidWritten[1] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fObsPx1);

      poShoebox->vGetPeakStartEnd(&nImageStart, &nImageEnd);
      for (nImage = nImageStart; nImage <= nImageEnd;nImage++) {
          fObsRotStart = (float) (nImage + poShoebox->nGetOffset(2))*m_fRotIncrement + m_fOverallRotStart;
          fObsRotEnd = (float) (nImage + 1 + poShoebox->nGetOffset(2))*m_fRotIncrement + m_fOverallRotStart;
          fCalcRotMid = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotMid,(float) (nImage + poShoebox->nGetOffset(2) + 0.5)*m_fRotIncrement + m_fOverallRotStart);
          poRefln->vSetIntensity(fRLP*poShoebox->pfGetPerSliceIntensity()[nImage]);
          poRefln->vSetSigmaI(fRLP*poShoebox->pfGetPerSliceSigma()[nImage]);
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fProfitIntensity,fRLP*poShoebox->pfGetPerSliceProfitIntensity()[nImage]);
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fProfitSigmaI,fRLP*poShoebox->pfGetPerSliceProfitSigma()[nImage]);
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_nPartialFlag,(poShoebox->nGetPerSliceProfitFlag(nImage)));
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fBackgroundSigma,(poShoebox->pfGetPerSliceBackgroundSigma()[nImage]));
          poShoebox->vGetPerSliceEllipsoidRelObsCent(nImage,&a2x2fEllipsoidA[0][0],&a2fEllipsoidb[0],&fEllipsoidc,&a2fObsCentroidWritten[0]);
          poRefln->nPutGetEllipse(true,&a2x2fEllipsoidA[0][0],&a2fEllipsoidb[0],&fEllipsoidc);
          poRefln->vWrite(m_pFReflnPartialOut,m_pcSelectFieldsPartial);
          m_nReflnPartialOutCount++;
      };

      // Update information on "cutpoints"
      if (!((m_nReflnFileOutCount-1) % 10000)) {
          itr<int> anCutPoints;
          if (m_nReflnFileOutCount == 1) {
              anCutPoints[0] = 0;
          } else {
              if (m_poHeader->nGetValue((Cstring) "DTPROFIT_CUTPOINTS",(int) 1,&anCutPoints[0]))
                  return 1;
              anCutPoints.setsize(((int) anCutPoints[0]) + 1);
              if (m_poHeader->nGetValue((Cstring) "DTPROFIT_CUTPOINTS",(int) (1 + anCutPoints[0]),&anCutPoints[0])) {
                  printf("ERROR: Generating cut point information.\n");
                  return 1;
              };
          };
          anCutPoints[(int) (anCutPoints[0]) + 1] = (m_nReflnFileOutCount-1);
          anCutPoints[(int) (anCutPoints[0]) + 2] = (m_nReflnPartialOutCount-(nImageEnd - nImageStart+1));
          anCutPoints[0] += 2;
          m_poHeader->nReplaceValue((Cstring) "DTPROFIT_CUTPOINTS",(int) (1 + anCutPoints[0]),&anCutPoints[0]);
      };

     
  };

  nExt = poShoebox->nGetExt(2);
  nFill = poShoebox->nGetLayersFilled();
  if (   (nFill != poShoebox->nGetLayersFilled())
      || (nExt  != poShoebox->nGetExt(2)) )
    {
      cout << "UNFILLED!: nFill, nExt, nFs, nEs: "
           << nFill << ", " << nExt << ", " << poShoebox->nGetLayersFilled()
           << ", " << poShoebox->nGetExt(2) << "\n";
    }
  return (nStat);
}


#define PRINTSHOEBOX(nMask,poRefln,arg1,arg2) {\
          poShoebox->nLoadPrint(!m_nRoguePlaceInImages);\
          nAddRogue(nMask,poRefln,*poShoebox,poShoebox->pfGetPrintMask());\
          }

#define BPRINTSHOEBOX(_nStat) ((m_bRogueActive) && (((m_nRogueMask==0) && ((_nStat)==0)) || ((_nStat) & m_nRogueMask)) && (!((_nStat) & m_nRogueNegMask)))


int
Cintegrate::nCalcMask(Crefln *poRefln, Cdetector *poDetector,
                         C3Ddata *poShoeboxIn)
{
   

  // This routine uses the static member scratch variables so that
  //   it is not always doing new/delete.

  C3Ddata *poShoebox;

  double f0;
  int   nStat = 0; // Local status
  float a3fWidth[3];
  float a3fCent[3];
  int   a3nCent[3];
  float pfCOG[3];
  float fWidth;
  int   nx,ny;
  int   nStart, nEnd;
  float fTemp;
  C3DdataDetectorArea oArea;

  // Select the shoebox to use, either in the poRefln or in the argument list

  if (NULL == poShoeboxIn)
    {
      poShoebox = (C3Ddata*)poRefln->m_pvUserData;
    }
  else
    {
      poShoebox = poShoeboxIn;
    }

    


  // Get the putative center peak pixel wrt to origin of this shoebox
  // Probably should check input to centroid to see if it makes sense

  a3fCent[0] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx0);
  a3nCent[0] = (int) a3fCent[0];
  a3fCent[1] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx1);
  a3nCent[1] = (int) a3fCent[1];
 
  //  if ((ABS(a3fCent[0] - 526) < 10) && (ABS(a3fCent[1] - 1584) < 10))
  //      nx=nx;    
  //if ((poRefln->nGetH() == 5) && (poRefln->nGetK() == -2) && (poRefln->nGetL() == -10))
  //    nx=nx;
    

  if (m_poProfit2DPerSlice) {
          m_poProfit2DPerSlice->vSetCurrentImage(poShoebox->m_nOrigOffset[2] + m_anSeqNum[0]);
          m_poProfit2DPerSlice->m_nDim0 = m_nDim0;
          m_poProfit2DPerSlice->m_nDim1 = m_nDim1;
  };
  
    
  float fJreso;
  fJreso     = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fResolution);
  fTemp      = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);
  fTemp      = fTemp - m_fOverallRotStart;
  a3fCent[2] = (fTemp / m_fRotIncrement);  // Images from beginning
  a3nCent[2] = (int) a3fCent[2];

  // Calculate width.  The width must be converted from units of degrees to units of "images" by dividing by the m_fRotIncrement (degrees per image)
  fWidth     = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth)/m_fRotIncrement;
  // 

  // Fill in fields that C3Ddata::nCalcMask() needs to know about.
  ((C3DdataInput&) *poShoebox).vInit();
 
  if (m_a2oDetectorProfiles[m_nProfileReading].nGetPeakArea(a3nCent[0],a3nCent[1],oArea,NULL)>=0)  {
  };
  //WAS In Flop(C3Ddata): poShoebox->m_bUseFastBackground = TRUE;
  // Copy values of the ellipsoid constants over.
  poShoebox->m_a2x2fEllipsoidAIn[0][0] = oArea.a2x2fEllipsoidA[0][0];
  poShoebox->m_a2x2fEllipsoidAIn[0][1] = oArea.a2x2fEllipsoidA[0][1];
  poShoebox->m_a2x2fEllipsoidAIn[1][0] = oArea.a2x2fEllipsoidA[1][0];
  poShoebox->m_a2x2fEllipsoidAIn[1][1] = oArea.a2x2fEllipsoidA[1][1];
  poShoebox->m_a2fEllipsoidbIn[0] =   oArea.a2fEllipsoidb[0];
  poShoebox->m_a2fEllipsoidbIn[1] =   oArea.a2fEllipsoidb[1];
  poShoebox->m_fEllipsoidcIn = oArea.fEllipsoidc;
  poShoebox->m_fMinProfitPeakIntensity = oArea.fIntensity;

  
  if (m_nProfitFlag & DTI_DO_2D_PROFIT) {
          poShoebox->m_bProfileFitWeak = true;
  };
  
  if (m_bKa12) {
      double a8fTempBuf[8];
      if (!poDetector->nCalcPixelShift((int) a3fCent[0],(int) a3fCent[1],
          &a8fTempBuf[0],&a8fTempBuf[2])) {
          vCopyVecND(2,&a8fTempBuf[0],&poShoebox->m_a2fShiftKa1[0]);
          vCopyVecND(2,&a8fTempBuf[2],&poShoebox->m_a2fShiftKa2[0]);
      };
  };

  
  // Check for preasure cell setting.
  while (m_fPCellAngle > 0.0) {
      // This 'while' loop is only used so that we can 'break' out when there is an error.
      double        a3fPCellRot[3];
      float         a3x3fRotMatHKL[3][3];
      float         a3x3fTemp[3][3];
      double        a3x3fTemp2[3][3];
      double        fAngle;
      float         a3fVDC[3];
      float         a3fS[3];
  
      
      if (poDetector->m_poSpatial->nPixeltoMM(
          poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx0),
          poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx1),
          &a3fVDC[0],&a3fVDC[1],&a3fVDC[2])) {
          nStat =  1 << Crefln::ms_nErrorIntProblem;     
          break;
      };

      vMulMat3DVec3D(m_a3x3fDetDD, a3fVDC, a3fS); // Scattered beam wavevector
      fNormVec3D(a3fS);
      
      vConvRotVec3DMat3D(poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid), m_a3fRotVec, &a3x3fRotMatHKL[0][0]);
      vMulMat3DMat3D(a3x3fRotMatHKL,m_a3x3fGonMat, a3x3fTemp);
      vCopyMat3D(&a3x3fTemp[0][0],&a3x3fTemp2[0][0]);
      vMulMat3DVec3D(a3x3fTemp2, m_a3fPCellVector,  a3fPCellRot);
      f0 = ABS(fDot3D(a3fPCellRot,a3fS));
      fAngle = acos(min(1.0,f0))/Gs_dRADIANS_PER_DEGREE;
      if (fAngle > m_fPCellAngle) {
          nStat  = 1 << Crefln::ms_nErrorIntProblem;
          break;
      };
      f0 = ABS(fDot3D(a3fPCellRot,m_a3fS0));
      fAngle = acos(min(1.0,f0))/Gs_dRADIANS_PER_DEGREE;
      if (fAngle > m_fPCellAngle) {
          nStat = 1 << Crefln::ms_nErrorIntProblem;
          break;
      };
      nStat = 0;
      break;
  } 

  // nStat is virtually always 0, so !nStat is virtually always true.

  if (!nStat) 
    {
      poShoebox->m_fFattening            = m_fFattening;
      poShoebox->m_fSigmaBelowBackground = poShoebox->m_fSigmaAboveBackground 
                                         = m_fSigmaAboveBackground;
      poShoebox->m_fDetectorGain         = m_fDetectorGain;
      poShoebox->m_bSpotChase            = m_bSpotChase;
      poShoebox->m_nMinPeakRad           = m_nMinPeakRad;
      poShoebox->m_nIntegrateWeakFlag    = m_nIntegrateWeakFlag;
      
      nStat = poShoebox->nCalcMask(a3fCent, fWidth);

      // If we received an error status, 
      // then check if we are printing this type of error.
      // If so, log it.
      
      if (BPRINTSHOEBOX(nStat))
	PRINTSHOEBOX(nStat,poRefln,a3fCent,fWidth);
      
      bool bNoBackgroundError = (0 == (nStat
          & (  nSetErrorMask(Crefln::ms_nErrorBackground, FALSE, NULL)
          | nSetErrorMask(Crefln::ms_nErrorBackgroundSD, FALSE, NULL) ) ) );
      
      if (bNoBackgroundError)
	{
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fBackground,
			     poShoebox->fGetPerPixelBackground());
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fBackgroundSigma,
			     poShoebox->fGetPeakBackgroundSigma());
	}
    }
  
  
  nSetErrorMasks(nStat,TRUE);

  if (nStat & (1 << Crefln::ms_nErrorTooFewPix)) {
      m_nNumRefWeak++;
  } else if (nStat & m_nErrorDoNotProcessMask)
      m_nNumRefBadErrors++;
  else if (nStat == 0)
      m_nNumRefNoErrors++;
  else 
      m_nNumRefSpecial++;

  if (!(nStat & m_nErrorDoNotProcessMask)) {

      // Fill in bitmap information.
      if ((m_bWriteScanBitmap)  && (m_poProfit2DPerSlice))
          m_poProfit2DPerSlice->nAddBitmap(*poShoebox);

      // Set the error mask, since it won't get set any other time.
     
      
      /* Calculation of rotation centroid:  
         This is a tricky operation due to the noise that can be added
         due to the uncertainty in rotation centroid.
         
         We assume that the refined values for the calculated centroid widths are 
         approx. correct.  We also assume that we should use the calculated rotation midpoints
         if they do not disagree with the observed rotation midpoints.
        
         We will distinguish three cases to be processed seperately:
         
         Case 1)  We found peak intensity on 1 image.
         Case 2)  We found peak intensity on 2 images.
         Case 3)  We found peak intensity on 3 or more images.

         Case 1)  Since we have so little data for rotation centroid, it would be incorrect to
                  assume that the centroid is in the middle of the rotation interval. Instead, we
                  will always use the calculated rotation midpoint if it lies completely in the rotation
                  range;  that is, [calc-rot-start,calc-rot-end] is a subset of [rotation-rocking-start,rotation-rocking-end].
                  
                  If not, then we might have to shift it so that it lies (as much as possible) in the rocking range for this image.
                  For this, we choose the minimal shift that causes one of the following to happen:

                    a)  calc-rot-start = rotation-rocking-start
                    b)  calc-rot-end   = rotation-rocking-end
                    c)  calc-rot-mid   = rotation-rocking-mid.

                  Note that (c) will only occur if the calculated rotation width is wider than the rocking width of an
                  image.  

        Case 2)   We have more information, but still very little.  Simply taking a weighted average of the
                  centroids will tend to cause errors.  To see why, imagine that a very skinny spot is located 
                  on the border between the two images.  If this spot is distributed 50/50 between the images, then
                  everything works out okay since the weighted centroid will get correctly calculated to be on the border
                  between the images.  On the other hand, if the spot is distributed 10/90, then the spot will get 
                  disproportionally shifted towards the more intense image;  this happens because we calculate centroids 
                  assuming intensity is located at the midpoint of each image rotation width.  

                  An alternate way to solve this is to assume some intensity distribution of the peak, and to shift the
                  peak in such a way that it stradles the two images.  This is the method that is employed below.

        Case 3)   We have enough information, since the spot is at least as wide as one image width.

      */

      double fObsRot[2];
      double fCalcRot[2];
      double fObsCent;
      double fCalcCent;
      int    nImageCount;
      int    nImageStart;
      int    nImageEnd;


      // Set peak preliminary intensity and sigma calculated in nCalcMask
      
      poRefln->vSetIntensity(poShoebox->fGetPeakIntensity());
      poRefln->vSetSigmaI(poShoebox->fGetPeakIntensitySigma());

      // Get the peak centroid calculated in poShoebox->nCalcMask


      poShoebox->vGetCentroid(pfCOG);

      // Apply the shoebox offsets

      for (nx = 0; nx < 3; nx++)
        {
          pfCOG[nx] += (float) poShoebox->nGetOffset(nx);
        }

      // Set observed centroid to the predicted centroid if we were weak.
      // Otherwise, use the observed centroid from the shoebox.
      
      if (!(nStat & (1 << Crefln::ms_nErrorTooFewPix))) {
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx0, pfCOG[0]);
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx1, pfCOG[1]);
      } else {
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx0,
              poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx0));          
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx1,
              poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx1));
      };

      poShoebox->vGetPeakStartEnd(&nImageStart, &nImageEnd);

      nImageCount = nImageEnd - nImageStart + 1;
      fObsRot[0] = (float) (nImageStart + poShoebox->nGetOffset(2))*m_fRotIncrement + m_fOverallRotStart;
      fObsRot[1] = (float) (nImageEnd + poShoebox->nGetOffset(2)+1)*m_fRotIncrement + m_fOverallRotStart;
      fObsCent   = (float) (pfCOG[2])*m_fRotIncrement + m_fOverallRotStart;
      fCalcRot[0] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid) - 0.5*poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth);
      fCalcRot[1] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid) + 0.5*poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotWidth);
      fCalcCent   = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);

      if (nImageCount == 1) {
          bool bUseSmartEstimation = TRUE;
          // Case 1)
          
          if (!bUseSmartEstimation) {
              // The user disabled smart estimation of centroid.
              fObsCent = fObsCent;
          } else if ((fObsRot[0]<=fCalcRot[0]) && (fObsRot[1]>=fCalcRot[1])) {
              // The whole calculated spot is inside the image limits.
              fObsCent = fCalcCent;
          } else if (fCalcRot[1] - fCalcRot[0] >= fObsRot[1] - fObsRot[0]) {
              // case (c):  The calculated width is greater than the observed width.
              // Just use the center.
              fObsCent = fObsCent;
          } else {
              if (fabs(fObsRot[0] - fCalcRot[0]) < fabs(fObsRot[1] - fCalcRot[1])) {
                  // case (a):  Shift so that starting ranges match.
                  fObsCent = fCalcCent + fObsRot[0] - fCalcRot[0];
              } else {
                  // case (b):  Shift so that ending ranges match.
                  fObsCent = fCalcCent + fObsRot[1] - fCalcRot[1];
              };
          };
      } else if (nImageCount == 2) {
          // Case 2)

          bool bUseTriangular = TRUE;
          bool bUseUniform = TRUE;
          double fPercentData[2];
         
          fPercentData[1] = (fObsCent - (fObsRot[0] + m_fRotIncrement*0.5))/m_fRotIncrement;
          fPercentData[0] = 1.0 - fPercentData[1];

          if (fPercentData[1]<0.0) {
              // If we happen to be using negative intensities, the centroid could have been shifted illegally.
              fObsCent = fObsRot[0] + m_fRotIncrement*0.5;
          } else if (fPercentData[0]<0.0) {
              // If we happen to be using negative intensities, the centroid could have been shifted illegally.
              fObsCent = fObsRot[1] - m_fRotIncrement*0.5;
          } else if (fObsRot[1] - fObsRot[0]<fCalcRot[1] - fCalcRot[0]) {
              // The calculated rotation range is larger than the observed rotation range.
              // Just keep the observed centroid, since it is acurate enough.
              fObsCent = fObsCent;
          } else if (bUseTriangular) {
              // Assume triangular shape.
              // From the centroid, we glean what percentage of intensity was on each image
              if (fPercentData[0]<0.5) 
                  fObsCent = (fObsRot[0] + m_fRotIncrement) + (1.0 - sqrt(fPercentData[0]*2.0))*(fCalcRot[1] - fCalcRot[0])*0.5;
              else
                  fObsCent = (fObsRot[0] + m_fRotIncrement) - (1.0 - sqrt(fPercentData[1]*2.0))*(fCalcRot[1] - fCalcRot[0])*0.5;
          } else if (bUseUniform) {
              // Assume a uniform shape.
              fObsCent = (fObsRot[0] + m_fRotIncrement) + (0.5 - fPercentData[0])*(fCalcRot[1] - fCalcRot[0]);
          } else {
              fObsCent = fObsCent;
          };
      } else {
          // Case 3)
          fObsCent = fObsCent;
      };
     
      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotMid, (float) fObsCent);          
      poShoebox->vGetPeakSize(a3fWidth);
      poRefln->vSetField(m_nFI_fWidth0, a3fWidth[0]);
      poRefln->vSetField(m_nFI_fWidth1, a3fWidth[1]);

      poShoebox->vGetIntensePeakStartEnd(&nStart,&nEnd);

      if ((-1 < nStart) && (poRefln->fGetIntensity()/max(1.0,poRefln->fGetSigmaI())>=5.0)
          && (!(nStat & (1 << Crefln::ms_nErrorTooFewPix))))
         {

          // Add the shoebox intensity profile.  These are used for seeing how well the mosaicity is predicting the spot.

          nx = nImageEnd - nImageStart + 1 - 2;
          if (    (nx>=0) && (nx<MAX_SHOEBOX_INTENSITY_PROFILE)
               && (poRefln->fGetIntensity()>=100.0))
            {
              m_pnShoeboxIntensityProfilesCounts[nx]++;         
              double f0 = 0.0;
              for (ny=nImageStart;ny<=nImageEnd;ny++)
                {
                  f0 += poShoebox->pfGetPerSliceIntensity()[ny];
                }
              if (f0 >= 100.0)
                {
                  for (ny=nImageStart;ny<=nImageEnd;ny++)
                    {
                      m_pfShoeboxIntensityProfilesIntensity[nx][ny-nImageStart] += poShoebox->pfGetPerSliceIntensity()[ny]/f0;
                    }
                }
            }

          // Get the minimum rot width that could account for the observed width
          // assuming a symmetric in rotation angle peak shape.

         fObsRot[0] = (float) (nStart + poShoebox->nGetOffset(2))*m_fRotIncrement + m_fOverallRotStart;
         fObsRot[1] = (float) (nEnd + poShoebox->nGetOffset(2)+1)*m_fRotIncrement + m_fOverallRotStart;

          fWidth = 2.0*min((fObsCent - fObsRot[0]),(fObsRot[1] - fObsCent));

          if (fWidth < 0.0) {
              fObsRot[0] = (float) (nImageStart + poShoebox->nGetOffset(2))*m_fRotIncrement + m_fOverallRotStart;
              fObsRot[1] = (float) (nImageEnd + poShoebox->nGetOffset(2)+1)*m_fRotIncrement + m_fOverallRotStart;
              fWidth = 2.0*min((fObsCent - fObsRot[0]),(fObsRot[1] - fObsCent));
          };

          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotWidth,
              fWidth);
      } else {
          // For whatever reason, this spot does not have width information.  Perhaps it is too weak?
          // Set the width to 0.0 so that it does not affect mosaicity.
          poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotWidth,
              (float) 0.0);          
      }

      // Set the rotation width sigma.
      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotSigma,
          poShoebox->fGetRotationSigma());


      // Set the error mask (still kept in nStat)
      poRefln->vSetField(m_nFI_nErrorFlag,
                         poRefln->nGetField(m_nFI_nErrorFlag)
                         | (nStat));

    
    }
  else
    {
      // Some errors occurred during the preliminary peak processing
      // Unsuccessful determination of observed centroid, leave observed
      // centroid same as calculated one

      poRefln->vSetIntensity(-999.0);
      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx0,
          poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx0));

      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx1,
          poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx1));

      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotMid,
          poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid));

      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotWidth,
                         (float) -999.0);
      poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotSigma, 
                         (float) -999.0);

      poRefln->vSetField(m_nFI_fWidth0, (float)-999.0);
      poRefln->vSetField(m_nFI_fWidth1, (float)-999.0);

      poRefln->vSetField(m_nFI_nErrorFlag,
                         poRefln->nGetField(m_nFI_nErrorFlag)
                         | (nStat));
    }
  return nStat;
}





int
Cintegrate::nUpdateRefine(Csource *poSource, Cdetector *poDetector,
                          Ccrystal *poCrystal, Cgoniometer *poCrysGonio,
                          Crotation *poRotation)
{
  // Update the member variables that
  // depend on the crystal, detector, source objects
  // These may change after every refinement.

  int nStat;
  nStat = poDetector->nCalcGetDDDN(&m_a3x3fDetDD[0][0], m_a3fDetDN);

  float a4fPolNorm[4];

  poSource->vCalcGetS0(m_a3fS0);

  // Get polarization normal and transfer to member variable for use elsewhere
  poSource->vGetPolarz(a4fPolNorm);
  for (int i = 0; i < 3; i++)
    {
      m_a3fPolNorm[i] = a4fPolNorm[i+1];
    }

  m_fLenS0 = fLenVec3D(m_a3fS0);

  //      -1 -1 -1 -1
  //  h = B  C  G  R   x
  //  -   =  =  =  =   -
  //
  //

  float a3x3fTemp[3][3];
  float a3x3fInvCrysOrientMat[3][3];
  float fDet;

  m_fWavelength = poSource->m_poWavelength->fGetWavelength();

  // Calculate and get CB
  //                   ==

  nStat = poCrystal->nCalcOrientMatrix(m_fWavelength);

  if (0 != nStat)
    return (nStat);

  poCrystal->vGetOrientMatrix(&a3x3fTemp[0][0]);

  //                   -1 -1
  // Invert CB  to get B  C
  //        ==         =  =

  fDet = fInvMat3D(&a3x3fTemp[0][0], &a3x3fInvCrysOrientMat[0][0]);
  if (0.0 == fDet)
    return (-1);

  float a3x3fInvGonMatrix[3][3];

  // Calculate and get G,  then invert it
  //                   =

  poCrysGonio->vCalcGetRotMatrix(&a3x3fTemp[0][0],poCrysGonio->nGetNum(m_poScan->m_poRotation->sGetName()));
  vCopyMat3D(&a3x3fTemp[0][0],&m_a3x3fGonMat[0][0]);
  fDet = fInvMat3D(&a3x3fTemp[0][0], &a3x3fInvGonMatrix[0][0]);
  if (0.0 == fDet)
    return (-1);

  // Compute                      -1 -1 -1
  //  m_a3x3InvCrysOrientGonMat = B  C  G
  //                              =  =  =

  vMulMat3DMat3D(a3x3fInvCrysOrientMat, a3x3fInvGonMatrix, m_a3x3fInvCrysOrientGonMat);

  // Calculate and get C and invert it
  //                   =

  poCrystal->nCalcRotMatrix();
  poCrystal->vGetRotMatrix(&a3x3fTemp[0][0]);

  float a3x3fInvCrysRotMatrix[3][3];

  fDet = fInvMat3D(&a3x3fTemp[0][0], &a3x3fInvCrysRotMatrix[0][0]);
  if (0.0 == fDet)
    return (-1);

  //                            -1 -1
  //  m_a3x3InvCrysRotGonMat  = C  G
  //                            =  =

  vMulMat3DMat3D(a3x3fInvCrysRotMatrix, a3x3fInvGonMatrix, m_a3x3fInvCrysRotGonMat);

  //                                        -1
  // One could use m_a3x3InvCrysRotGonMat = G  instead, but ...
  //                                        =
  //vCopyMat3D((float*)a3x3fInvGonMatrix , (float*)m_a3x3fInvCrysRotGonMat);

  return (nStat);
}
int
Cintegrate::nDeleteMarkedReflns(Creflnlist *poReflnlist, const int nFlag)
{
  // Delete reflections whose nDataFlag field is equal to nFlag,
  // since new reflections are
  // at the end of the list, we can loop from the end of the reflnlist
  // back towards the beginning until we have deleted all those flagged.

  int i;
  int nField1;
  int nNumDel = 0;
  int nNumReflns = poReflnlist->nGetNumReflns();
  Crefln  *poRefln;
  C3Ddata *poShoebox;

  for (i = nNumReflns-1; i>=0; i--)
    {
      poRefln = poReflnlist->poGetRefln(i);
      nField1 = poRefln->nGetField(m_nFI_nDataFlag);
      if ((nFlag == nField1) || (nFlag == -999))
        {
          if (NULL != poRefln->m_pvUserData)
            {
              poShoebox = (C3Ddata *)poRefln->m_pvUserData;
              vDeleteShoebox(poShoebox);
              poRefln->m_pvUserData = NULL;
            }

          poReflnlist->nDelete(i);
          nNumDel++;
        }
    }
  if (2 < m_nVerbose)
    {
      cout << "Reflections freed/deleted: " << nNumDel << '\n' << flush;
    }
  return (0);
}

int
Cintegrate::nDoRefine(Creflnlist *poReflnlist, Crefine *poRefine)
{
  //int     i,j;
  int     nStat, nNumReflns;
  //Crefln *poRefln;
  double fPrevMosaicity;

  // Copy reflections flagged as DoneStrong from the input to refine list

  if (5 < m_nVerbose)
    {
      cout << "... about to try to refine\n" << flush;
    }

  nNumReflns =  poReflnlist->nGetNumReflns();

  nStat = 0;
//  cout << "In nDoRefine, nNumReflnsUsable: " << oReflnlistRefine.nGetNumReflns()
//       << endl;
  if (m_nRefineMinRefln <= m_poReflnlistRefine->nGetNumReflns())
    {
      // Refinement code copied straight from dtrefine


      Cstring sRefineOptions;
      nStat = m_poHeader->nGetValue(Crefine::ms_sDtrefineOptions, &sRefineOptions);
      if (0 == nStat)
        {
          cout << "Refinement options found:\n    " << sRefineOptions << endl << flush;

          // TODO: If dtintegrate -display is NOT specified, then remove
          //       any possible '-display' from sRefineOptions
        }
      else
        {
          // Refinement options not in header, so create some default options

          sRefineOptions = "+All -cycles 100 -rej 1 1 2 -verbose 0 -go -verbose 1 -go ";
          cout << "WARNING: no refinement options found, using\n   "
               << sRefineOptions << endl << flush;
        }

      fPrevMosaicity = poRefine->m_ppoCrystals[0]->fGetMosaicity();
      poRefine->nSetNoReflectionMerges(TRUE);
      if (m_nPrerefineFlag)
	poRefine->nSetRotOffsetAdjust(TRUE);
      else
	poRefine->nSetRotOffsetAdjust(FALSE);

      nStat = poRefine->nDoRefinement(sRefineOptions, m_poReflnlistRefine);

      int nx = m_nNumRefineResults;
      if (0 == nStat)
        {
          
          if ( (0 < m_nSizeRefineResults) || (NULL == m_ptRefineResults) )
            {
              poRefine->m_ppoCrystals[0]->vGetCell(
                              &m_ptRefineResults[m_nNumRefineResults].fA,
                              &m_ptRefineResults[m_nNumRefineResults].fB,
                              &m_ptRefineResults[m_nNumRefineResults].fC,
                              &m_ptRefineResults[m_nNumRefineResults].fAlp,
                              &m_ptRefineResults[m_nNumRefineResults].fBet,
                              &m_ptRefineResults[m_nNumRefineResults].fGam);
              poRefine->m_ppoCrystals[0]->vGetOrientAngles(
                              &m_ptRefineResults[m_nNumRefineResults].fRot1,
                              &m_ptRefineResults[m_nNumRefineResults].fRot2,
                              &m_ptRefineResults[m_nNumRefineResults].fRot3);
              double a3fMosaicity[3];
              double a6dDetParams[6];
              poRefine->m_ppoCrystals[0]->vGetMosaicity(a3fMosaicity);
              m_ptRefineResults[m_nNumRefineResults].fMos
                = a3fMosaicity[0];
              
              poRefine->m_ppoDetector[0]->nGetStandardDetectorParams(
                                                 a6dDetParams);
              m_ptRefineResults[m_nNumRefineResults].fT1 = a6dDetParams[3];
              m_ptRefineResults[m_nNumRefineResults].fT2 = a6dDetParams[4];
              m_ptRefineResults[m_nNumRefineResults].fT3 = a6dDetParams[5];
              m_ptRefineResults[m_nNumRefineResults].fDetRot1 = a6dDetParams[0];
              m_ptRefineResults[m_nNumRefineResults].fDetRot2 = a6dDetParams[1];
              m_ptRefineResults[m_nNumRefineResults].fDetRot3 = a6dDetParams[2];
              m_ptRefineResults[m_nNumRefineResults].fSrcRot1 
                = poRefine->m_poSource->fGetRotation(0);
              m_ptRefineResults[m_nNumRefineResults].fSrcRot2 
                = poRefine->m_poSource->fGetRotation(1);

              poRefine->vGetRMS(&m_ptRefineResults[m_nNumRefineResults].fResDeg,
                                &m_ptRefineResults[m_nNumRefineResults].fResMM);

              // Write out 2 lines of refinement results for a plot program to find

              printf("PLOTCRYS:\n"
                  "  Seq    a      b      c     alp    bet    gam    Rot1    Rot2    Rot3    Mos\n");
              printf(" %4d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %7.2f %7.2f %7.2f %6.3f\n",
                  m_ptRefineResults[nx].nSeq,
                  m_ptRefineResults[nx].fA, 
                  m_ptRefineResults[nx].fB, 
                  m_ptRefineResults[nx].fC, 
                  m_ptRefineResults[nx].fAlp,
                  m_ptRefineResults[nx].fBet,
                  m_ptRefineResults[nx].fGam,
                  m_ptRefineResults[nx].fRot1,
                  m_ptRefineResults[nx].fRot2,
                  m_ptRefineResults[nx].fRot3,
                  m_ptRefineResults[nx].fMos);
              
              printf("PLOTDET:\n"
                  "  Seq  TransX TransY  Dist   RotZ   RotX   RotY  SrcRot1 SrcRot2  rmsMM rmsDeg\n");
              printf(" %4d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %7.2f %7.2f %7.3f %6.3f\n",
                  m_ptRefineResults[nx].nSeq,
                  m_ptRefineResults[nx].fT1,
                  m_ptRefineResults[nx].fT2,
                  m_ptRefineResults[nx].fT3,
                  m_ptRefineResults[nx].fDetRot1,
                  m_ptRefineResults[nx].fDetRot2,
                  m_ptRefineResults[nx].fDetRot3,
                  m_ptRefineResults[nx].fSrcRot1,
                  m_ptRefineResults[nx].fSrcRot2,
                  m_ptRefineResults[nx].fResMM,
                  m_ptRefineResults[nx].fResDeg);
              fflush(stdout);
              if (m_nNumRefineResults < m_nSizeRefineResults)
                  m_nNumRefineResults++;
          }
	  m_ptRefineResults[nx].fMosMod = m_ptRefineResults[nx].fMos;
          if (m_fFixedMosaicity > 0.0)
            {
              cout << "\nINFO: Integrate Mosaicity fixed to " << m_fFixedMosaicity
                   << "\n\n";
              poRefine->m_ppoCrystals[0]->vSetMosaicity(m_fFixedMosaicity);
	      m_ptRefineResults[nx].fMosMod = m_fFixedMosaicity;
            }
          else if (0.0 == m_a2fMosaicityModel[0])
            {
              cout << "\nINFO: Integrate Mosaicity fixed to " << m_a2fMosaicityModel[1]
                   << "\n\n";
              poRefine->m_ppoCrystals[0]->vSetMosaicity(m_a2fMosaicityModel[1]);
	      m_ptRefineResults[nx].fMosMod = m_a2fMosaicityModel[1];
            }
          else if ( (0.0 == m_a2fMosaicityModel[1]) && (1.0 == m_a2fMosaicityModel[0]) )
            {
              float fTempMos;
              fTempMos = poRefine->m_ppoCrystals[0]->fGetMosaicity();
              cout << "\nINFO: Integrate Mosaicity accepts the " << fTempMos 
                   << " deg value from refinement.\n\n";
	      m_ptRefineResults[nx].fMosMod = fTempMos;
            }
          else
            {
              float fTempMos;
              fTempMos = poRefine->m_ppoCrystals[0]->fGetMosaicity();
              fTempMos = (m_a2fMosaicityModel[0] * fTempMos) + m_a2fMosaicityModel[1];
              cout << "\nINFO: Integrate Mosaicity reset to " << fTempMos
                   << " by (-mosaicitymodel " << m_a2fMosaicityModel[0] << " "
                   << m_a2fMosaicityModel[1] << ") option.\n\n";
              poRefine->m_ppoCrystals[0]->vSetMosaicity(fTempMos);
	      m_ptRefineResults[nx].fMosMod = fTempMos;
            }

          nStat = poRefine->nUpdateHeader(m_poHeader,m_nTwinID);
          Cimage_header oHeader;
          oHeader = *m_poHeader;
          (void) oHeader.nReplaceValue(Cimage_header::ms_sComment,
                                       "Cintegrate refinement results");

          // Write out results up to now in case something crashes
#ifdef SSI_PC
          Cstring   sOutputHeaderFileName(m_sOutputHeader);
#else
          Cstring   sOutputHeaderFileName(sDtrekGetPrefix() + "dtintegrate.head");
#endif
          (void) oHeader.nWrite(sOutputHeaderFileName);
          m_poReflnlistRefine->vSetOverwrite(TRUE);
          (void) m_poReflnlistRefine->nWrite(sDtrekGetPrefix() + "dtintrefine.ref");

          if( 1 == m_nDisplay )
          {
#ifdef SSI_PC
                CCrclHelper::GetInstance()->vSendIntegrateUpdateDisplay( NULL, "dtintrefine.ref", TRUE, FALSE, TRUE );
#else
          #ifndef NO_X_WINDOWS
                (void) poRefine->nNotifyDisplay(sGetCWD() + sDtrekGetPrefix() + "dtintrefine.ref",
                                                sGetCWD() + sDtrekGetPrefix() + "dtintegrate.head");
          #endif          
#endif          
          }
        }

      // Delete all reflections that have now been refined.
      m_poReflnlistRefine->vDeleteAll();
      m_nDoRefineCount++;
  } 
  else 
  {
      printf("Not enough strong reflections to do a refinement.\n");
      nStat = 1;
  }

  return (nStat);
}

int
Cintegrate::nInitProfit2D(Cdetector& oDetector) {
    Cstring sTemp;

    m_poProfit2DPerSlice = new Cprofit2DPerSlice;
    //WAS In Flop(C3Ddata):m_poProfit2DPerSlice->m_sProfitDiagnostic = m_sProfitDiagnostic;    
    m_poProfit2DPerSlice->vSetInterpConditions(m_nProfitNumReflections,m_nProfitMaxImageRange);
    m_poProfit2DPerSlice->m_nDim0 = m_nDim0;
    m_poProfit2DPerSlice->m_nDim1 = m_nDim1;
    C3Ddata::m_poProfit = m_poProfit2DPerSlice;
    return 0;
};



void
Cintegrate::vSetResolution(const float fResoMin, const float fResoMax)
{
  m_a2fResolution[0] = max(fResoMin, fResoMax);
  m_a2fResolution[1] = min(fResoMin, fResoMax);
}

void
Cintegrate::vListResults(void)
{
  // Output a summary of the results

  int i;

  printf("\n\n*** dtintegrate ***\n");
  vListRefineResults();
  printf("\n"
         "================================================================\n"
         "Summary of results for scan rotation from %9.3f to %9.3f\n"
         "         with image sequence numbers from %5d     to %5d\n"
         "----------------------------------------------------------------\n\n",
         m_fOverallRotStart, m_fOverallRotEnd,
         m_anSeqNum[0], m_anSeqNum[1]);
  printf("Total reflections predicted:          %10d\n", m_nNumPredictedUnique);
  printf("Total reflections with no errors:     %10d\n", m_nNumRefNoErrors);
  printf("Total reflections processed weak:     %10d\n", m_nNumRefWeak);
  printf("Total reflections with bad errors*:   %10d\n", m_nNumRefBadErrors);
  printf("Total reflections partial in scan*:   %10d\n",

            m_nNumRefPartialStart + m_nNumRefPartialEnd);

  printf("\nReflection integration status codes\n"
           "====================================\n");
  printf(" %22s  %10s\n", "Status", "Num reflns");
  printf("------------------------------------\n");
  printf("%22s: %10d\n", "No errors, no warnings", m_nNumRefNoErrors);
  printf("------------------------------------\n");
  for (i = 0; i < nMaxErrorStatusTypes; i++)
  {
      if (!(m_nErrorDoNotPrintMask & (1 << i)))
      {
          if (   ("?" != Crefln::ms_asErrorStrings[i])
              && ("" != Crefln::ms_asErrorStrings[i]) )
          {
              printf("%22s: %10d",  Crefln::ms_asErrorStrings[i].string(),
                  m_anNumReflns[i]);
              if (0 == (m_nErrorDoNotProcessMask & (1 << i) ))
              {
                  printf("  ");
              }
              else
              {
                  printf(" *");
              }
              printf("%s\n",  Crefln::ms_asErrorPhrases[i].string());
          }
      }
  }
  printf("====================================\n"
         "              *Rejected from output.\n");
  fflush(stdout);
}

void
Cintegrate::vListIntensityProfiles(void) {
    int nx,ny,nz;
    int nNormalCounts;
    double fPercent;
    const int nMaxCounts = 60;


    printf("\n\n");
    for (nx=0;nx<MAX_SHOEBOX_INTENSITY_PROFILE-2;nx++) {
        if (m_pnShoeboxIntensityProfilesCounts[nx]>=20) {
            printf("Average rocking curve for strong reflns on %d images (%d contributors).\n",nx+2,m_pnShoeboxIntensityProfilesCounts[nx]);
            printf("==========================================================================\n");
            for (ny=0;ny<nx+2;ny++) {
                fPercent = (m_pfShoeboxIntensityProfilesIntensity[nx][ny]/m_pnShoeboxIntensityProfilesCounts[nx]);
                printf("Image %2d %4.1f%% ",ny,fPercent*100.0);
                nNormalCounts = (int) (fPercent*nMaxCounts);
                for (nz=0;nz<nNormalCounts;nz++)
                    printf("*");
                printf("\n");
            };
            printf("==========================================================================\n\n");
        };
    };
    return;
};

void
Cintegrate::vDeleteFreeShoeboxes(void)
{
  // Delete all shoeboxes in the free list of shoeboxes

  int i;   // Loop counter

  if (8 < m_nVerbose)
    {
      cout << "Cintegrate::vDeleteFreeShoeboxes called with "
           << m_nNumShoeBFree
           << " freed shoeboxes in the list.\n" << flush;
    }

  if (NULL != m_ppoShoeBFree)
    {
      for (i = 0; i < m_nNumShoeBAlloc; i++)
        {
          if (NULL != m_ppoShoeBFree[i])
            {
              delete m_ppoShoeBFree[i];
              m_ppoShoeBFree[i]    = NULL;
              m_pnShoeBFreeSize[i] = -1;
            }
        }
      delete [] m_ppoShoeBFree;
      m_ppoShoeBFree    = NULL;
      delete [] m_pnShoeBFreeSize;
      m_pnShoeBFreeSize = NULL;
      m_nNumShoeBAlloc  = 0;
      m_nNumShoeBFree   = 0;
    }
}

C3Ddata*
Cintegrate::poGetShoebox(const int nExt0, const int nExt1, const int nExt2,
                         const int nOrig0, const int nOrig1, const int nOrig2)
{
    C3Ddata* poShoebox;

    poShoebox = new C3Ddata(nExt0, nExt1, nExt2,
        nOrig0, nOrig1, nOrig2, m_eSBdatatype);

    return poShoebox;
};

void
Cintegrate::vDeleteShoebox(C3Ddata *poShoebox)
{
  if (NULL == poShoebox) return;

  delete poShoebox;
  return;
}


void Cintegrate::vInitBackgroundMask() {
    if (m_pcBackgroundMask==NULL) {
        // Allocate background.
        m_pcBackgroundMask = new unsigned char[m_nDim0*m_nDim1/8+1];
    };
    memset(m_pcBackgroundMask,0,m_nDim0*m_nDim1/8+1);
};

void Cintegrate::vAddBackgroundMask(C3Ddata* poShoebox) {
    int nOffset;
    int nx,ny;
    for (ny=0;ny<poShoebox->m_nExt[1];ny++) {
        nOffset = (poShoebox->m_nOrigOffset[0]) + (poShoebox->m_nOrigOffset[1] + ny)*m_nDim0;
        for (nx=0;nx<poShoebox->m_nExt[0];nx++,nOffset++) {
            m_pcBackgroundMask[nOffset/8] |= (1 << (nOffset % 8));
        };
    };
};

int
Cintegrate::nSetErrorMask(const int nErrorNumber, const bool bCountIt,
                          int *pnMask)
{
  // Computes and sets the error mask for error number nErrorNumber
  // in the input value *pnMask (if pnMask is given).
  // If bCountIt is true, then count the error in the member arrays, too.
  // Returns the mask for the nErrorNumber.  (TODO: or 0 if out of range, etc)

    if (nErrorNumber==Crefln::ms_nErrorIntProblem)
        pnMask = pnMask;
  int nMaskReturn;
  nMaskReturn =  (1 << nErrorNumber);
  if (NULL != pnMask)
    {
      *pnMask = *pnMask | nMaskReturn;
    }
  if (bCountIt)
    {
      m_anNumReflns[nErrorNumber]++;
    }
  return (nMaskReturn);
}

int 
Cintegrate::nSetErrorMasks(const int nErrorNumbersBitwise, const bool bCountIt) {
    int nx,nz;
    int nSeriousErrorNumberBitwise = 0;
    int nNotSeriousErrorNumberBitwise = 0;
    
    for (nx=0;nx<32;nx++) {
        if (nErrorNumbersBitwise & ( ((unsigned int) 1) << nx)) {
            if ((!nSeriousErrorNumberBitwise) && (m_nErrorDoNotProcessMask & ((unsigned int) 1 << nx)))
                nSeriousErrorNumberBitwise = nx+1;
            if ((!nNotSeriousErrorNumberBitwise) && (!(m_nErrorDoNotProcessMask & ((unsigned int) 1 << nx))))
                nNotSeriousErrorNumberBitwise = nx+1;
        };
    }; nz = 0;
    if (nSeriousErrorNumberBitwise)
        nSetErrorMask(nSeriousErrorNumberBitwise-1, bCountIt, &nz); 
    else if (nNotSeriousErrorNumberBitwise)
        nSetErrorMask(nNotSeriousErrorNumberBitwise-1, bCountIt, &nz); 
    return 0;
};


void
Cintegrate::vAddExcludeResoRing(const float fResoMin, const float fResoMax)
{
  // Add a ring of excluded resolution to the array of exclusions

  int i;
  float *pfMin, *pfMax;
  float *pfMinTemp, *pfMaxTemp;

  pfMin     = m_pfExcludeReso[0];
  pfMax     = m_pfExcludeReso[1];
  pfMinTemp = m_pfExcludeReso[0];
  pfMaxTemp = m_pfExcludeReso[1];

  // Allocate new memory to hold all excluded rings

  m_pfExcludeReso[0] = new float [m_nNumExcludeResoRings + 1];
  m_pfExcludeReso[1] = new float [m_nNumExcludeResoRings + 1];

  // Copy previous excluded rings to new location

  for (i = 0; i < m_nNumExcludeResoRings; i++)
    {
      m_pfExcludeReso[0][i] = *pfMinTemp++;
      m_pfExcludeReso[1][i] = *pfMaxTemp++;
    }

  // Add new excluded ring

  m_pfExcludeReso[0][m_nNumExcludeResoRings] = max(fResoMin, fResoMax);
  m_pfExcludeReso[1][m_nNumExcludeResoRings] = min(fResoMin, fResoMax);

  // Increment excluded ring counter

  m_nNumExcludeResoRings++;

  // Free memory of previous

  delete [] pfMin;
  delete [] pfMax;
}


//////////////////////////////////////////////////////////////////
/////////////////// Image drawing routines ///////////////////////
//////////////////////////////////////////////////////////////////

#include <ctype.h>

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


char* g_pcBitText = 
"  ###      #     #####   #####  #       #######  #####  #######  #####   #####   \n"
" #   #    ##    #     # #     # #    #  #       #     # #    #  #     # #     #  \n"
"#     #  # #          #       # #    #  #       #           #   #     # #     #  \n"
"#     #    #     #####   #####  #    #  ######  ######     #     #####   ######  \n"
"#     #    #    #             # #######       # #     #   #     #     #       #  \n"
" #   #     #    #       #     #      #  #     # #     #   #     #     # #     #  \n"
"  ###    #####  #######  #####       #   #####   #####    #      #####   #####   \n"
"                                                                                 \n"
"                                                                                 \n"
"   ##    #####    ####   #####   ######  ######   ####   #    #     #         #  \n"
"  #  #   #    #  #    #  #    #  #       #       #    #  #    #     #         #  \n"
" #    #  #####   #       #    #  #####   #####   #       ######     #         #  \n"
" ######  #    #  #       #    #  #       #       #  ###  #    #     #         #  \n"
" #    #  #    #  #    #  #    #  #       #       #    #  #    #     #    #    #  \n"
" #    #  #####    ####   #####   ######  #        ####   #    #     #     ####   \n"
"                                                                                 \n"
"                                                                                 \n"
" #    #  #       #    #  #    #   ####   #####    ####   #####    ####    #####  \n"
" #   #   #       ##  ##  ##   #  #    #  #    #  #    #  #    #  #          #    \n"
" ####    #       # ## #  # #  #  #    #  #    #  #    #  #    #   ####      #    \n"
" #  #    #       #    #  #  # #  #    #  #####   #  # #  #####        #     #    \n"
" #   #   #       #    #  #   ##  #    #  #       #   #   #   #   #    #     #    \n"
" #    #  ######  #    #  #    #   ####   #        ### #  #    #   ####      #    \n"
"                                                                                 \n"
"                                                                                 \n"
"                                                          # #                 #  \n"
" #    #  #    #  #    #  #    #   #   #  ######           # #    #   #       #   \n"
" #    #  #    #  #    #   #  #     # #       #          #######   # #       #    \n"
" #    #  #    #  #    #    ##       #       #             # #   #######    #     \n"
" #    #  #    #  # ## #    ##       #      #            #######   # #     #      \n"
" #    #   #  #   ##  ##   #  #      #     #               # #    #   #   #       \n"
"  ####     ##    #    #  #    #     #    ######           # #           #        \n"
"\n"
"\n"
;

int textout(Cimage& rc,int p0,int p1,int nValue0,int nValue1,char* str) {
   static bool bInit=FALSE;
   static bool cBits[40][8][8];      // Character x line y is in cBits[x][y][0..7]
   static char* cpChars="0123456789abcdefghijklmnopqrstuvwxyz #*/";

   char* pcPoint;
   char* pcTemp;
   int   nx,ny,nz,nw;

   int   n0,n1;
   int   nPt,nIndex,nVal;
   
   if (bInit==FALSE) {
       pcTemp = new char[strlen(g_pcBitText)+1];
       strcpy(pcTemp,g_pcBitText);
       for (nx=0;nx<40;nx+=10) {
           for (ny=0;ny<8;ny++) {
               if ((nx==0) && (ny==0)) 
                   pcPoint = strtok(pcTemp,"\n");
               else
                   pcPoint = strtok(NULL,"\n");
               for (nz=0;nz<10;nz++) {
                   for (nw=0;nw<8;nw++) {
                       cBits[nx+nz][ny][nw]=((pcPoint[nz*8+nw]=='#')?TRUE:FALSE);
                   };
               };
           };
           bInit=TRUE;
       }
       delete[] pcTemp;
   };
   // Find position of upper left hand corner.
   n0=p0-(strlen(str)/2)*8;
   n1=p1-4;
   if (n0<0) n0=0;
   if (n1<0) n1=0;
   if (n0+strlen(str)*8>=rc.nGetDimension(0)) n0=rc.nGetDimension(0)-strlen(str)*8-2;
   if (n1+8>=rc.nGetDimension(1)) n1=rc.nGetDimension(1)-8-1;
   // Print!
   for (nPt=0;nPt< strlen(str);nPt++) { 
      for (nIndex=0;nIndex<40;nIndex++) 
          if (tolower(str[strlen(str) - 1 - nPt])==cpChars[nIndex]) 
            break;
      if (nIndex==40) 
          return 1;
      for (nx=0;nx<8;nx++) {
          for (ny=0;ny<8;ny++) {
              if (cBits[nIndex][7-ny][7-nx]) 
                  nVal=nValue1; 
              else 
                  nVal=nValue0;
              rc.nSetPixel(n0+8*nPt+nx,n1+ny,(unsigned short int)nVal);
          };
      };

   };
return 0;
};

int Cintegrate::nSetRogueData(int nMask,int nNegMask,Cstring& sTemplate,Cstring& sReflnlist,int nImageNumber,int nMinDimension) { 
    int nx;

    m_nRogueMask = nMask;
    m_nRogueNegMask = nNegMask; 
    m_sRogueTemplate = sTemplate; 
    m_bRogueActive = TRUE; 
    m_nRoguePlaceInImages = (nImageNumber!=-1);
    m_nRogueImageNumber = nImageNumber;
    m_nRogueMinBoxDim = nMinDimension;

    if (sReflnlist.length()) {
        if (m_poRogueReflnlist) 
            delete m_poRogueReflnlist;
        m_poRogueReflnlist = new Creflnlist(sReflnlist);
        if (!m_poRogueReflnlist->bIsAvailable()) {
            printf("ERROR: Could not open rogue input file '%s'\n",sReflnlist.string());
            return 1;
        };
        if (m_pnRogueReflnlistHKLSort)
            delete[] m_pnRogueReflnlistHKLSort;
        m_pnRogueReflnlistHKLSort = new int[m_poRogueReflnlist->nGetNumReflns()];
        m_poRogueReflnlist->nExpandGetField(Creflnlist::ms_snPackedHKL);
        for (nx=0;nx<m_poRogueReflnlist->nGetNumReflns();nx++) {
            (*m_poRogueReflnlist)[nx].vSetField(m_poRogueReflnlist->m_nFI_nPackedHKL,(*m_poRogueReflnlist)[nx].nPackHKL());
        };
        m_poRogueReflnlist->vSort(eReflnField_int_type,m_poRogueReflnlist->m_nFI_nPackedHKL,m_pnRogueReflnlistHKLSort);
    };
    return 0;
};


void Cintegrate::vSetRogueImageCutoff(int nImage) {
    if ((m_poRogueImage) && 
        (m_nRoguePlaceInImages) &&
        (nImage-m_nRogueImageNumber>=2)) {
        m_poRogueImage->nWrite(m_sRogueNextImage);
        m_poRogueImage = NULL;
        printf("Rogue image '%s' written\n",m_sRogueNextImage.string());
    };
};

int Cintegrate::nAddRogue(int nMask,Crefln* poRefln,C3Ddata& oShoebox,float* pfExternData,int* pnMask) {
    int nStat;
    int n0,n1,n2;
    int nx,ny,nz;
    float* pfPointer;
    int* pnMaskPointer;
    double fAverage;
    double fAverageCt;
    double f0,f1;

    
    if (m_poRogueReflnlist) {
        nx = poRefln->nPackHKL();
        ny = m_poRogueReflnlist->nFindFirst(m_poRogueReflnlist->m_nFI_nPackedHKL,nx,m_pnRogueReflnlistHKLSort);        
        if (ny==-1)
            return 1;
        nz = 0;
        while ((ny<m_poRogueReflnlist->nGetNumReflns()) && 
            ((*m_poRogueReflnlist)[m_pnRogueReflnlistHKLSort[ny]].nGetField(m_poRogueReflnlist->m_nFI_nPackedHKL)==nx)) {
            f0 = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcRotMid);
            f1 = (*m_poRogueReflnlist)[m_pnRogueReflnlistHKLSort[ny]].fGetField(m_poRogueReflnlist->m_nFI_fCalcRotMid);
            if (fabs(f0-f1)<5.0) {
                nz = 1;
                break;
            } else
                ny++;
        };
        if (!nz)
            return 1;
    };

    oShoebox.vGetPeakStartEnd(&nx,&ny);
    if ((m_nRogueMinBoxDim!=-1) && 
        ((nx==-1) ||
        (ny - nx + 1<m_nRogueMinBoxDim)))
        return 0;
    // Use the user provided external data if necc.
    if (pfExternData) 
        pfPointer = pfExternData;
    else
        pfPointer = oShoebox.m_pfData;
    pnMaskPointer = pnMask;
    do {
        // If there is not yet an image, create a new one and a new name template.
        if (!m_poRogueImage) {
            // Don't do anything for 'place in image' format if we have already written an image.
            // if (m_nRoguePlaceInImages && m_nRogueImageCount)
            //    return 0;
            m_nRogueImageCount++;
            m_nRogueShoeboxImageCount = 0;
            m_nRoguePx0 = 0;
            m_nRoguePx1 = 0;
            m_nMaxRoguePx1 = 0;
            m_nRogueShoeboxes = 0;
            if ((m_nRoguePlaceInImages) && (m_nRogueDim0==0)) {
                m_nRogueDim0 = m_nDim0;
                m_nRogueDim1 = m_nDim1;
            } else {
                m_nRogueDim0 = 2000;
                m_nRogueDim1 = 2000;
            };


            m_poRogueImage = new Cimage(m_nRogueDim0,m_nRogueDim1,eImage_uI2);
            // Initialize image memory.
            m_poRogueImage->nSetNextPixel(0,0);
            for (n0=0;n0<m_nRogueDim1*m_nRogueDim0;n0++)
                m_poRogueImage->vSetNextPixel((unsigned short int) 0);

            
            Cscan oScan(m_sRogueTemplate,m_nRogueImageCount,1);
            m_sRogueNextImage = oScan.sGetImageName();
        };
        
        if (!m_nRoguePlaceInImages) {
            
            // If an old image existed, is there enough space to fit the next shoebox?
            if (oShoebox.m_nExt[0] + 1 + m_nRoguePx0<m_nRogueDim0) {
                if (((oShoebox.m_nExt[1])*oShoebox.m_nExt[2]) + 2 + m_nRoguePx1 < m_nRogueDim1) {
                    m_nMaxRoguePx1 = max( m_nMaxRoguePx1 , ((oShoebox.m_nExt[1])*oShoebox.m_nExt[2]) + 2 + m_nRoguePx1 + 8*m_nRoguePrintNumber );
                    nStat = 0;
                } else 
                    nStat = 1;
            } else {
                m_nRoguePx1 = m_nMaxRoguePx1;
                m_nRoguePx0 = 0;
                m_nMaxRoguePx1 = 0;
                if ((((oShoebox.m_nExt[1])*oShoebox.m_nExt[2]) + 2 + m_nRoguePx1 < m_nRogueDim1) &&
                    (oShoebox.m_nExt[0] + 1 + m_nRoguePx0 < m_nRogueDim0)) {
                    m_nMaxRoguePx1 = max( m_nMaxRoguePx1 , ((oShoebox.m_nExt[1])*oShoebox.m_nExt[2]) + 2 + m_nRoguePx1 );
                    nStat = 0;
                } else
                    nStat = 1;
            };
            if (nStat) {
                if (!m_nRogueShoeboxImageCount) {
                    printf("WARNING:  Could not paste shoebox of dimensions %dx%dx%d into image.\n",
                        oShoebox.m_nExt[0],
                        oShoebox.m_nExt[1],
                        oShoebox.m_nExt[2]);
                    break;
                } 
                // We write out the previous shoebox.
                m_poRogueImage->nWrite(m_sRogueNextImage);
                delete m_poRogueImage;
                m_poRogueImage = NULL;
            };           
        } else
            nStat = 0;
    } while (nStat);

    // Count the maximum number of shoeboxes per image.
    m_nRogueShoeboxes++;
    m_nMaxRogueShoeboxes = max(m_nRogueShoeboxes,m_nMaxRogueShoeboxes);

    if (!m_nRoguePlaceInImages) {

        fAverage = 0.0;
        fAverageCt = 0.0;        
        if (oShoebox.m_nExt[2]>3)
            oShoebox.m_nExt[2] = oShoebox.m_nExt[2];
        // Now, we paste the shoebox in the image.
        for (n2 = 0; n2< oShoebox.m_nExt[2]; n2++) {
            for (n1 = 0; n1< oShoebox.m_nExt[1]; n1++) {
                for (n0 = 0; n0< oShoebox.m_nExt[0]; n0++) {
                    f0 = max(0.0,(((pnMaskPointer) && (*pnMaskPointer==0))?(0.0):(*pfPointer)));
                    fAverageCt++;
                    fAverage+= f0;
                    (m_poRogueImage->*m_poRogueImage->prnSetPixel)(
                        n0 + m_nRoguePx0,
                        n1 + m_nRoguePx1 + oShoebox.m_nExt[1]*n2,f0);
                    
                    pfPointer++;
                    if (pnMaskPointer)
                        pnMaskPointer++;
                };
            };
        };
        // Print borders.
        for (n0=0;n0< oShoebox.m_nExt[0];n0++) {
            for (n1=0;n1<2;n1++) 
                (m_poRogueImage->*m_poRogueImage->prnSetPixel)(n0 + m_nRoguePx0,n1 + m_nRoguePx1 + oShoebox.m_nExt[1]*oShoebox.m_nExt[2],(float) m_nRogueBorder);
        };
        for (n1=0;n1< oShoebox.m_nExt[1]*oShoebox.m_nExt[2]+2;n1++) {
            (m_poRogueImage->*m_poRogueImage->prnSetPixel)(oShoebox.m_nExt[0] + m_nRoguePx0,n1 + m_nRoguePx1,(float) m_nRogueBorder);
        };
        if (m_nRoguePrintNumber) {
            char pcBuf[20];
            for (nx=0,ny=0;nx<32;nx++) {
                if (nMask & (1 << nx)) {
                    pcBuf[ny++] = 'A' + nx;
                };
            };
            pcBuf[ny] = 0;
            textout(*m_poRogueImage,
                m_nRoguePx0 + oShoebox.m_nExt[0]/2,
                m_nRoguePx1 + oShoebox.m_nExt[1]*oShoebox.m_nExt[2] + 2 +4,
                (int) (fAverage/fAverageCt*0.5),
                (int) (fAverage/fAverageCt), 
                pcBuf);
        };

        // Increment the pixel pointer, and the number of shoeboxes written to this image.        
        m_nRoguePx0 += oShoebox.m_nExt[0] + 1;
        m_nRogueShoeboxImageCount++;


    } else {


        // Place the appropriate slice of the shoebox directly into the image.
        
        if ((m_nRogueImageNumber >= oShoebox.m_nOrigOffset[2] + m_anSeqNum[0]) && (m_nRogueImageNumber <= oShoebox.m_nOrigOffset[2] + m_anSeqNum[0]  + oShoebox.m_nExt[2]- 1)) { 
            if (pnMaskPointer)
                pnMaskPointer += (m_nRogueImageNumber - (oShoebox.m_nOrigOffset[2] + m_anSeqNum[0]))*oShoebox.m_nExt[0]*oShoebox.m_nExt[1];
            pfPointer += (m_nRogueImageNumber - (oShoebox.m_nOrigOffset[2] + m_anSeqNum[0]))*oShoebox.m_nExt[0]*oShoebox.m_nExt[1];
            for (n1 = 0; n1< oShoebox.m_nExt[1]; n1++) {
                for (n0 = 0; n0< oShoebox.m_nExt[0]; n0++) {
                    int nPix0,nPix1;
                    float f0;
                    nPix0 = n0 + oShoebox.m_nOrigOffset[0];
                    nPix1 = n1 + oShoebox.m_nOrigOffset[1];
                    
                    f0 = (m_poRogueImage->*m_poRogueImage->prfGetPixel)(nPix0,nPix1);
                    if (f0!=1.0)
                        f0 = max(0.0,(((pnMaskPointer) && (*pnMaskPointer==0))?(0.0):(*pfPointer)));
                    // Set all 0.0 pixels to 1.0 so that these are not masked out by overlapping shoeboxes.
                    if (f0 == 0.0)
                        f0 = 1.0;
                    
                    (m_poRogueImage->*m_poRogueImage->prnSetPixel)(nPix0, nPix1,f0);
                    pfPointer++;
                    if (pnMaskPointer)
                        pnMaskPointer++;
                };
            };            
        };
    };


    return 0;
};


void Cintegrate::vSetFixedMosaicity(double fFixMosaicity,Cpredict* poPredict) {
    m_fFixedMosaicity = fFixMosaicity;
    Ccrystal oCrystal(*m_poHeader);
    oCrystal.vSetMosaicity(fFixMosaicity);
    oCrystal.nUpdateHeader(m_poHeader);
    if (poPredict) {
        poPredict->vSetCrysMosaicity(fFixMosaicity);
    };

};


void Cintegrate::vTermRogue() {
    if (m_poRogueImage) {
        m_poRogueImage->nWrite(m_sRogueNextImage);
        delete m_poRogueImage;
        m_poRogueImage = NULL;
        printf("NOTICE:  Total of %d Rogue galleries printed with template '%s'\n"
               "Maximum of %d shoeboxes in each gallery.\n"
            ,m_nRogueImageCount,m_sRogueTemplate.string(),m_nMaxRogueShoeboxes);
    };
    if (m_poRogueReflnlist) {
        delete m_poRogueReflnlist;
        m_poRogueReflnlist = NULL;
    };
    if (m_pnRogueReflnlistHKLSort) {
        delete[] m_pnRogueReflnlistHKLSort;
        m_pnRogueReflnlistHKLSort = NULL;
    };
    
    return;
};


void 
Cintegrate::vListRefineResults(void)
{
  int nx,ny;//,nz;
  int nVar1,nVar2;
  const int nVariables = 20;
  const int nTopCorr = 5;
  const double fMinTopCorr = 0.9;
  double afCorr[nTopCorr];
  int    anCorr1[nTopCorr];
  int    anCorr2[nTopCorr];
  tagRefineResults oResults;
  double  aafCorr[nVariables][nVariables];
  double  afAvg[nVariables];
  double  afVar[nVariables];
  float* apfCorr[nVariables] = {
        &oResults.fA, &oResults.fB, &oResults.fC, &oResults.fAlp, &oResults.fBet, &oResults.fGam, 
        &oResults.fRot1, &oResults.fRot2, &oResults.fRot3, &oResults.fMos,
        &oResults.fT1, &oResults.fT2, &oResults.fT3, &oResults.fDetRot1, 
        &oResults.fDetRot2, &oResults.fDetRot3, &oResults.fSrcRot1, &oResults.fSrcRot2,
        &oResults.fResMM, &oResults.fResDeg
  };  
  char* apcCorr[nVariables] = {
      "a","b","c","alp","bet","gam","Rot1","Rot2","Rot3",
          "Mos","TransX","TransY","Dist","RotZ","RotX","RotY",
          "SrcRot1","SrcRot2","rmsMM","rmsDeg"
  };
  
  if ( (0 >= m_nNumRefineResults) || (NULL == m_ptRefineResults) )
  {
      return;
  }

  
  // Initialize the results.
  for (nx=0;nx<nVariables;nx++) {
      afVar[nx] = 0.0;
      afAvg[nx] = 0.0;
      for (ny=0;ny<nVariables;ny++) 
          aafCorr[nx][ny] = 0.0;
  };
  // Take the averages.
  for (nx = 0; nx < m_nNumRefineResults; nx++) {
      oResults = m_ptRefineResults[nx];
      for (nVar1 = 0; nVar1 < nVariables; nVar1++) 
          afAvg[nVar1] += *apfCorr[nVar1];
  };
  for (nVar1 = 0; nVar1 < nVariables; nVar1++) 
      afAvg[nVar1]/=m_nNumRefineResults;

  // Take deviations.
  for (nx = 0; nx < m_nNumRefineResults; nx++) {
      oResults = m_ptRefineResults[nx];
      for (nVar1 = 0; nVar1 < nVariables; nVar1++) 
          afVar[nVar1] += (*apfCorr[nVar1] - afAvg[nVar1])*(*apfCorr[nVar1] - afAvg[nVar1]);
  };
  // Take correlations.
  for (nVar1 = 0; nVar1 < nVariables; nVar1++) {
      for (nVar2 = 0; nVar2 < nVariables; nVar2++) {
          if (nVar1 != nVar2) {
              for (nx = 0; nx < m_nNumRefineResults; nx++) {
                  oResults = m_ptRefineResults[nx];
                  aafCorr[nVar1][nVar2] += (*apfCorr[nVar1] - afAvg[nVar1])*(*apfCorr[nVar2] - afAvg[nVar2]);
              };
              if ((afVar[nVar1]>1e-5) && (afVar[nVar2]>1e-5)) {
                  aafCorr[nVar1][nVar2]/=sqrt(afVar[nVar1]);
                  aafCorr[nVar1][nVar2]/=sqrt(afVar[nVar2]);
              } else {
                  aafCorr[nVar1][nVar2]= 0.0;
                  aafCorr[nVar1][nVar2]= 0.0;
              };
          };
      }
  };

  // Find the top correlations. (maximum ABS(R))
  for (nx=0;nx<nTopCorr;nx++)
      afCorr[nx] = 0.0;
  
  for (nVar1 = 0; nVar1 < nVariables; nVar1++) {
      for (nVar2 = 0; nVar2 < nVariables; nVar2++) {
          if (nVar1 < nVar2) {
              for (nx=0,ny=0;nx<nTopCorr;nx++) {
                  if (ABS(afCorr[nx])<ABS(afCorr[ny]))
                      ny = nx;
              }
              if (ABS(aafCorr[nVar1][nVar2]) > ABS(afCorr[ny])) {
                  afCorr[ny] = aafCorr[nVar1][nVar2];
                  anCorr1[ny] = nVar1;
                  anCorr2[ny] = nVar2;
              };
          };
      };
  };

  
  printf("Summary of crystal refinement results during integration\n"
"------------------------------------------------------------------------------------\n"
"  Seq    a      b      c     alp    bet    gam    Rot1    Rot2    Rot3   Mos  MosMod\n"
"------------------------------------------------------------------------------------\n");
  for (nx = 0; nx < m_nNumRefineResults; nx++)
    {
      printf(" %4d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %7.2f %7.2f %7.2f %5.2f %5.2f\n",
             m_ptRefineResults[nx].nSeq,
             m_ptRefineResults[nx].fA, 
             m_ptRefineResults[nx].fB, 
             m_ptRefineResults[nx].fC, 
             m_ptRefineResults[nx].fAlp,
             m_ptRefineResults[nx].fBet,
             m_ptRefineResults[nx].fGam,
             m_ptRefineResults[nx].fRot1,
             m_ptRefineResults[nx].fRot2,
             m_ptRefineResults[nx].fRot3,
             m_ptRefineResults[nx].fMos,
             m_ptRefineResults[nx].fMosMod);

      //  float fT1, fT2, fT3, fDetRot1, fDetRot2, fDetRot3, fSrcRot1, fSrcRot2;
      //  float fResMM, fResDeg;
      //  int   nSeq;
      if (nx > 0)
        {
          // Compute sums for later averaging

          m_ptRefineResults[0].fA    += m_ptRefineResults[nx].fA;
          m_ptRefineResults[0].fB    += m_ptRefineResults[nx].fB;
          m_ptRefineResults[0].fC    += m_ptRefineResults[nx].fC;
          m_ptRefineResults[0].fAlp  += m_ptRefineResults[nx].fAlp;
          m_ptRefineResults[0].fBet  += m_ptRefineResults[nx].fBet;
          m_ptRefineResults[0].fGam  += m_ptRefineResults[nx].fGam;
          m_ptRefineResults[0].fRot1 += m_ptRefineResults[nx].fRot1;
          m_ptRefineResults[0].fRot2 += m_ptRefineResults[nx].fRot2;
          m_ptRefineResults[0].fRot3 += m_ptRefineResults[nx].fRot3;
          m_ptRefineResults[0].fMos  += m_ptRefineResults[nx].fMos;
          m_ptRefineResults[0].fMosMod  += m_ptRefineResults[nx].fMosMod;
        }
    }
  printf(
"------------------------------------------------------------------------------------\n");
  if (0 < m_nNumRefineResults)
    {
      m_ptRefineResults[0].fA    /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fB    /= (float)m_nNumRefineResults; 
      m_ptRefineResults[0].fC    /= (float)m_nNumRefineResults; 
      m_ptRefineResults[0].fAlp  /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fBet  /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fGam  /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fRot1 /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fRot2 /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fRot3 /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fMos /= (float)m_nNumRefineResults;
      m_ptRefineResults[0].fMosMod /= (float)m_nNumRefineResults;
      printf(" Avg: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %7.2f %7.2f %7.2f %5.2f %5.2f\n\n",
             m_ptRefineResults[0].fA,
             m_ptRefineResults[0].fB,
             m_ptRefineResults[0].fC,
             m_ptRefineResults[0].fAlp,
             m_ptRefineResults[0].fBet,
             m_ptRefineResults[0].fGam,
             m_ptRefineResults[0].fRot1,
             m_ptRefineResults[0].fRot2,
             m_ptRefineResults[0].fRot3,
             m_ptRefineResults[0].fMos,
             m_ptRefineResults[0].fMosMod);

      if (0.0 == m_a2fMosaicityModel[0])
	{
	  printf("INFO: Please note Mosaicity was actually fixed to %.2f"
		 "\n      by action of the -mosaicitymodel %.2f %.2f option.\n",
         m_a2fMosaicityModel[1],
		 m_a2fMosaicityModel[0],
		 m_a2fMosaicityModel[1]);
	}
      else if ( (0.0 != m_a2fMosaicityModel[1]) || (1.0 != m_a2fMosaicityModel[0]) )
	{
	  printf("INFO: Please note Mosaicity was actually adjusted"
		 "\n      by action of the -mosaicitymodel %.2f %.2f option.\n",
		 m_a2fMosaicityModel[0],
		 m_a2fMosaicityModel[1]);
	}

      if (NULL != m_poHeader)
        {
          printf("   Updating header with average crystal results...\n\n");
          float a6fCell[6];
          a6fCell[0] = m_ptRefineResults[0].fA;
          a6fCell[1] = m_ptRefineResults[0].fB;
          a6fCell[2] = m_ptRefineResults[0].fC;
          a6fCell[3] = m_ptRefineResults[0].fAlp;
          a6fCell[4] = m_ptRefineResults[0].fBet;
          a6fCell[5] = m_ptRefineResults[0].fGam;
          (void) m_poHeader->nReplaceValue(Ccrystal::ms_sCrystalXUnitCell,
                                           6, a6fCell);
          a6fCell[0] = m_ptRefineResults[0].fRot1;
          a6fCell[1] = m_ptRefineResults[0].fRot2;
          a6fCell[2] = m_ptRefineResults[0].fRot3;
          (void) m_poHeader->nReplaceValue(Ccrystal::ms_sCrystalXOrientAngles,
                                           3, a6fCell);
	  // Use the fMosMod value!
          (void) m_poHeader->nReplaceValue(Ccrystal::ms_sCrystalXMosaicity,
                                           m_ptRefineResults[0].fMosMod);
          // The header will get written out by dtintegrate.cc
        }
    }
  else
    printf("\n");

  printf("Summary of detector and source refinement results during integration\n"
"------------------------------------------------------------------------------\n"
"  Seq  TransX TransY  Dist   RotZ   RotX   RotY  SrcRot1 SrcRot2  rmsMM rmsDeg\n"
"------------------------------------------------------------------------------\n");
  for (nx = 0; nx < m_nNumRefineResults; nx++)
    {
      printf(" %4d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %7.2f %7.2f %7.3f %6.3f\n",
             m_ptRefineResults[nx].nSeq,
             m_ptRefineResults[nx].fT1,
             m_ptRefineResults[nx].fT2,
             m_ptRefineResults[nx].fT3,
             m_ptRefineResults[nx].fDetRot1,
             m_ptRefineResults[nx].fDetRot2,
             m_ptRefineResults[nx].fDetRot3,
             m_ptRefineResults[nx].fSrcRot1,
             m_ptRefineResults[nx].fSrcRot2,
             m_ptRefineResults[nx].fResMM,
             m_ptRefineResults[nx].fResDeg);
      }
  printf(
"------------------------------------------------------------------------------\n\n");

  // 8 is an arbitrary number in the next line

  if (8 < m_nNumRefineResults)
    {
      for (nx=0;nx<nTopCorr;nx++)
        {
          if (ABS(afCorr[nx])>=fMinTopCorr)
            {
              printf("INFO:  There is a correlation of %.3lf between variables %s and %s\n",
                     afCorr[nx],apcCorr[anCorr1[nx]],apcCorr[anCorr2[nx]]);
            }
        }
    }
  else
    printf("INFO: Not enough refinements in above lists to determine correlations.\n");
  return;
}

// Prefind and prerefine on a specified number of batches
// return 0, if ok, or 1, if error
int  Cintegrate::nDoPrefind()
{
    int         nStat = 0;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Figure out the first and last image for prefind
    if( NULL == m_poScan )
        return 1;

    int         nImageSequenceStep = m_poScan->nGetSeqInc();

    int         nFirstImagePrefind = m_anSeqNum[0];
    int         nLastImagePrefind = m_anSeqNum[0] + nImageSequenceStep * (m_nPrefindFlag * m_nImagesPerRefineBatch - 1);
//+2011-06-16 JWP
// If (m_nPrefindFlag < 0), then use abs(m_nPrefindFlag) images evenly spaced from first image in the scan to the last image
// in the scan
    int         nTotalImagesPrefind = 1;
//-2011-06-16 JWP


    if (0 > m_nPrefindFlag)
      {
        nLastImagePrefind = m_anSeqNum[1];
//+2011-06-16 JWP
	// Adjust nImageSequenceStep
        nTotalImagesPrefind = m_anSeqNum[1] - m_anSeqNum[0] + 1;
        nImageSequenceStep = nTotalImagesPrefind / (ABS(m_nPrefindFlag)-1);
	if (1 > nImageSequenceStep) nImageSequenceStep = 1;
//-2011-06-16 JWP
      }
    else if( nLastImagePrefind > m_anSeqNum[1] )
      {
	nLastImagePrefind = m_anSeqNum[1];
	nTotalImagesPrefind = m_nPrefindFlag * m_nImagesPerRefineBatch;  // total number of images to do peak search on
	nImageSequenceStep = 1;
      }
    else
      {
	nTotalImagesPrefind = m_nPrefindFlag * m_nImagesPerRefineBatch;  // total number of images to do peak search on
	nImageSequenceStep = 1;
      }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    Creflnlist      oPrefindList;	   

    Cfind           oFind(*m_poHeader, 
                          &oPrefindList, 
                          NULL, 
//+2011-06-16 JWP
                          nTotalImagesPrefind);  // total number of images to do peak search on
//-2011-06-16 JWP    

    oFind.nExpandReflnlist(&oPrefindList, false, true);

    
    cout << "\n================================================================\n"
         << "Performing Pre-find + refinement with image(s) " << nFirstImagePrefind 
         << " to " << nLastImagePrefind << " Sequence Step " << nImageSequenceStep << " ..."
         << "\n================================================================\n" << flush;

    if ( 0 != (nStat = oFind.nDoSearch(m_a2fResolution[0],
                                      m_a2fResolution[1],
                                      eFind_2D_method,
                                      m_fPrefindSigma,
                                      nFirstImagePrefind,
                                      nLastImagePrefind,
                                      nImageSequenceStep)) )
        return 1;
    
    
    Cstring     sPrefindReflnlistFilename = sDtrekGetPrefix() + "prefind.ref";
    if( 0 != (nStat = oFind.nWriteReflectionList(sPrefindReflnlistFilename)) )
        return 1;

    Crefine*     poRefine  = new Crefine(*m_poHeader, &oPrefindList);

    poRefine->nExpandReflnlist(&oPrefindList);
    
    Cstring sRefineOptions;
    nStat = m_poHeader->nGetValue(Crefine::ms_sDtrefineOptions, &sRefineOptions);
    if( 0 != nStat )
        return 1;

    cout << "\n==========================\n"
         <<   "Pre-find + refinement ...          \n"
         <<   "==========================\n" << endl << flush;
    nStat = poRefine->nDoRefinement(sRefineOptions, &oPrefindList);
    cout << "============================\n"
         <<   "End of Pre-find + refinement\n"
         <<   "============================\n" << endl << flush;

    delete  poRefine;

    return nStat;
}
