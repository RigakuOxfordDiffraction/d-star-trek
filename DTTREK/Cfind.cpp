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
// Cfind.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Cfind which implements
//    the spot finding encapsulation of d*TREK.
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

#include "Dtrek.h"
#include "dtreksys.h"
#include "Cfind.h"         // Class definition and prototypes
#include "dtrekdefs.h"
#include "Cstat.h"
#include "dtarray.h"


#ifndef NO_X_WINDOWS
#ifndef SSI_PC

#include "CXprop.h"

#else // ifdef SSI_PC

#include "CrclHelper.h"

#endif
#endif
#define nMAXSPOTSIZE  60

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


static bool bListOfReflnIndexesHasDifferentHKLs(itr<int>& anReflnIndexes, Creflnlist& oReflnlist);

//+Definitions, constants, and initialization of static member variables

Cstring Cfind::ms_sDtfindOptions               = D_K_DtfindOptions;
Cstring Cfind::ms_sDtfindSpotSize              = D_K_DtfindSpotSize;
Cstring Cfind::ms_sDtfindAvgSpotsPerImage      = D_K_DtfindAvgSpotsPerImage;
Cstring Cfind::ms_sDtfindFractionSaturated     = D_K_DtfindFractionSaturated;
Cstring Cfind::ms_sDtfindPhotonValues          = D_K_DtfindPhotonValues;
Cstring Cfind::ms_sDtfindPhotonIoverSig        = D_K_DtfindPhotonIoverSig;


//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cfind::Cfind(Cimage_header& oHeader, 
             Creflnlist* poReflnlistOut,
             Cnonunf *poNonunf,
             int nNumberOfImages): 
m_oShoebox(1,1,1,0,0,0)
{

  (void) nInitValues();
  m_oHeader          = oHeader;
  m_poScan           = new Cscan(oHeader);
  m_bNewScan         = TRUE;
  if (NULL == poNonunf)
    {
      m_poScan->m_poNonunf = new Cnonunf(oHeader);
      m_bNewNonunf   = TRUE;
    }
  else
    {
      m_poScan->m_poNonunf = poNonunf;
      m_bNewNonunf   = FALSE;
    }
  m_a2nSeqNum[0]     = m_poScan->nGetSeqNum();
  
  if( nNumberOfImages <= 0 )
      nNumberOfImages = 10;     // Look at the first 10 images by default

  m_poScan->vSetImgNum(nNumberOfImages-1);            
  
  m_a2nSeqNum[1]     = m_poScan->nGetSeqNum();
  
  m_poScan->vInitSeqNum();

  m_poReflnlist      = poReflnlistOut;
  m_bNewReflnlist    = FALSE;


  m_poGoniometer = new Cgoniometer(m_oHeader,(Cstring) D_K_CrystalPrefix);
  m_poRotation = new Crotation(m_oHeader);
  if ((!m_poGoniometer->bIsAvailable()) || (!m_poRotation->bIsAvailable())) {
      delete m_poGoniometer;
      m_poGoniometer = NULL;
      delete m_poRotation;
      m_poRotation = NULL;
  };


  // This gets dimensions of image, if these are 0, then we should
  // defer until first image is available.

  (void) oHeader.nGetDimensions(&m_nDim0, &m_nDim1);

  // If the image file does not contain binary data, then the
  // values of m_a4nRect and m_a4nCircle should be 0!

  if ( (0 < m_nDim0) && (0 < m_nDim1) )
    {
      // Default search rect is whole image less 1% border along each edge.

      m_a4nRect[0]    = m_nDim0 / 100;
      m_a4nRect[1]    = m_nDim1 / 100;
      m_a4nRect[2]    = m_nDim0 - m_a4nRect[0];
      m_a4nRect[3]    = m_nDim1 - m_a4nRect[1];

      // Default search circle is centered at mid-point of image
      // and includes whole image

      m_a4nCircle[0]  = m_nDim0 / 2;
      m_a4nCircle[1]  = m_nDim1 / 2;
      m_a4nCircle[2]  = 0;
      m_a4nCircle[3] = (int) (sqrt( (double)(  (m_a4nCircle[0] * m_a4nCircle[0])
                                             + (m_a4nCircle[1] * m_a4nCircle[1])
                                            ) ) );
      m_nRadMin2         = m_a4nCircle[2] * m_a4nCircle[2];
      m_nRadMax2         = m_a4nCircle[3] * m_a4nCircle[3];

    }

    // RB per JWP: to make sure that the ellipsoids from the header 
    // are not used for the peak search. The rationale: When the search 
    // is run first time, no ellipsoid info is available yet.
    // When the search is run again, and the ellipsoid info is available, 
    // the peak search may  produce a slightly different result, 
    // which may be confusing for the user - the users
    // normally don't change the header files.
    Cimage_header         oTmpHeader;          // an empty header
    m_oDetectorAreasIn.nInitValues(oTmpHeader);

  m_eMethod     = eFind_2D_method;

  Cstring sTemp;
  Cstring sDetName;

  // Try to get the sigmas for x and y pixels.
  if (oHeader.nGetValue(sTemp=D_K_SpatialXYSigma,2,&m_a2fDetXYSigmas[0])) {
      // Didn't work.  Try to get the detector name.
      if (!oHeader.nGetValue(sTemp=D_K_DetectorNames,&sDetName)) {
          // Have the detector name.  Try to get with name prefix.
          if (oHeader.nGetValue(sTemp=sDetName+D_K_SpatialXYSigma,2,&m_a2fDetXYSigmas[0])) 
            {
              // That didn't work either.  So we have failed.

//              cout << "INFO:  Was not able to find "<< D_K_SpatialXYSigma << " or "
//                   << sDetName << D_K_SpatialXYSigma << "\n" 
//                   << "Estimating <x,y> centroid uncertainties as <" 
//                   << m_a2fDetXYSigmas[0] << "," << m_a2fDetXYSigmas[1] <<">\n" << flush;
          }
      } else {
//              cout << "INFO:  Was not able to find "<< D_K_SpatialXYSigma << "\n" 
//                   << "Estimating <x,y> centroid uncertainties as <" 
//                   << m_a2fDetXYSigmas[0] << "," << m_a2fDetXYSigmas[1] <<">\n" << flush;
      };
  };

}

Cfind::Cfind(Cscan* poScanIn, Creflnlist* poReflnlistOut): m_oShoebox(1,1,1,0,0,0)
{
  (void) nInitValues();
  m_poScan       = poScanIn;
  m_bNewScan     = FALSE;
  m_poReflnlist  = poReflnlistOut;
  m_bNewReflnlist= FALSE;
  m_a2nSeqNum[0] = m_poScan->nGetSeqNum();
  m_poScan->vSetImgNum(9);            // Look at first 10 images by default
  m_a2nSeqNum[1] = m_poScan->nGetSeqNum();
  m_poScan->vInitSeqNum();
  m_eMethod      = eFind_2D_merge_method;
}


Cfind::Cfind(): m_oShoebox(1,1,1,0,0,0)
{
 (void) nInitValues();
}

Cfind::~Cfind()
{
    if (NULL != m_po3Ddata)
    {
        delete m_po3Ddata;
        m_po3Ddata = NULL;
    }
    if (m_bNewNonunf && (NULL != m_poScan->m_poNonunf) )
    {
        delete m_poScan->m_poNonunf;
        m_poScan->m_poNonunf = NULL;
    }
    if (m_bNewScan && (NULL != m_poScan))
    {
        delete m_poScan;
        m_poScan = NULL;
    }
    if (m_bNewReflnlist && (NULL != m_poReflnlist))
    {
        delete m_poReflnlist;
        m_poReflnlist = NULL;
    }
    if (m_poGoniometer)
    {
        delete m_poGoniometer;
        m_poGoniometer = NULL;
    }
    if (m_poRotation)
    {
        delete m_poRotation;
        m_poRotation = NULL;
    };
}

int Cfind::nInitValues()
{
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  m_poXprop = NULL;
  m_psXpropUpdate = NULL;
#endif

  m_eMethod          = eFind_unknown_method;
  m_poScan           = NULL;
  m_poReflnlist      = NULL;
  m_sOptionString    = "-h";
  m_poGoniometer     = NULL;
  m_poRotation       = NULL;

  m_bNewScan         = FALSE;
  m_bNewNonunf       = FALSE;

  // This following C3Ddata object is a trick
  // to make the C3Ddata class more efficient

  m_po3Ddata         = new C3Ddata(1, 1, 1, 0, 0, 0);

  m_nDisplay         = 0;
  m_fSigma           = 5.0;
  m_fMinimum         = 0.0;
  m_fPeakMinimum     = 0.0;
  m_a4nCircle[0]     = -1;  // Center Px0
  m_a4nCircle[1]     = -1;  // Center Px1
  m_a4nCircle[2]     = 0;   // Min radius
  m_a4nCircle[3]     = 50;  // Max radius
  m_nRadMin2         = m_a4nCircle[2] * m_a4nCircle[2];
  m_nRadMax2         = m_a4nCircle[3] * m_a4nCircle[3];
  m_a4nRect[0]       = -1;  // Min Px0
  m_a4nRect[1]       = 0;   // Min Px1
  m_a4nRect[2]       = 0;   // Max Px0
  m_a4nRect[3]       = 0;   // Max Px1
  m_a2nSpotWindow[0] = 0;
  m_a2nSpotWindow[1] = 0;
  m_bIncludeSat      = TRUE;
  m_bIncludeBackground = FALSE;
  m_fNear            = 0.0;
  m_nPeakFilter      = 6;
  m_nDumpRefln       = 0;
  m_nBLoops          = 2;
  m_nDim0            = 0;
  m_nDim1            = 0;
  m_a2nSeqNum[0]     = 0;
  m_a2nSeqNum[1]     = 0;
  m_nEveryNthImage   = 1;  // Surprise!  This is not used, DTMainFind.cpp 
  m_fResolutionMin   = 999999.0f;
  m_fResolutionMax   = 0.00001f;
  m_a2fDetXYSigmas[0]= 1; 
  m_a2fDetXYSigmas[1]= 1; 
  m_fDefaultRockingWidth = 0.5;

  m_a2nPMTInfo[0] = 0;
  m_a2nPMTInfo[1] = 0;
  
  m_nTotalImages = 0;
  m_nTotalSpotsProcessed = 0;
  m_nTotalSatSpotsProcessed = 0;


  m_bList     = true;
  
  m_bPeakInfo = false;

  m_bRank = false;
  
  return (0);
}

// Find peaks in the given Image:

int
Cfind::nFind2D(Cimage& oImage, 
               const float fRot,
               Cdetector *poDetectorIn,  
               Csource *poSourceIn,
               bool bAddMergeInformation)
               
{
    int     i, j, k, l, m, n;          // Loop counters
    int     nRad0, nRad1;
    float   fAvg;
    float   fSD;
    float   fMin;
    float   fMax;
    float   fThreshold;
    float   fValue;
    float   fTemp;
    float   fDetReso;
    float   fSatValue;
    float   fNear;
    float   fInc;
    int     nStat;
    int     nArea[4];
    int     a2nBackground[2];
    Crefln *poRefln;
    Creflnlist *poReflnlist;
    Cgoniometer* poGoniometer;
    
    bool    bNewSource   = FALSE;
    bool    bNewDetector = FALSE;
    
    Cstring sReflnlistTmp;
    Cstring sSeqTmp;
    static  int nImgSeq = 0;
    
    
    fInc     = m_poScan->m_poRotation->fGetIncrement();
    if (0.0 == fInc) {
        // Prevent divide by zero later
        fInc = m_fDefaultRockingWidth;
    };
    
    
    if (!oImage.bIsAvailable())
    {
        cout << "ERROR Cfind::nFind2D, Image not available!" << endl;
        return (1);
    }
    
    fSatValue = oImage.fGetSatValue();
    
    nStat = nExpandReflnlist(NULL,bAddMergeInformation,true);
    
    if (0 != nStat)
    {
        cout << "ERROR Cfind::nFind2D, Could not expand reflection list!\n";
        return (nStat);
    }
    
    if (   (0 > m_poReflnlist->m_nFI_fObsPx0)
        || (0 > m_poReflnlist->m_nFI_fObsPx1)
        //|| (0 > m_poReflnlist->m_nFI_fObsPxPeak0)
        //|| (0 > m_poReflnlist->m_nFI_fObsPxPeak1)
        || (0 > m_poReflnlist->m_nFI_fObsSharpness)
        || (0 > m_poReflnlist->m_nFI_fObsRotMid) 
        )
    {
        cout << "ERROR Cfind reflection list missing necessary fields!\n";
        return (-2);
    }
          
    poReflnlist = new Creflnlist(*m_poReflnlist); // Temporary list
    Crefln oRefln (poReflnlist, 0, 0, 0);         // Create refln with hkl = 0 0 0

    poGoniometer = new Cgoniometer(oImage.m_oHeader,(Cstring) D_K_CrystalPrefix);
    if (!poGoniometer->bIsAvailable()) {
        printf("ERROR:  Could not initialize goniometer.\n");
        delete poReflnlist;
        delete poGoniometer;
        return (-3);
    };
    
    // Correct the entire image for non-uniformity of response once and
    // do not correct later any pixels in a 2D search.
    
    nStat = m_poScan->m_poNonunf->nCorrectImage(&oImage);
    if (0 != nStat)
    {
        printf("ERROR Cfind::nFind2D, Could not correct image for non-uniformity!\n");
        delete poReflnlist;
        delete poGoniometer;
        return (nStat);
    }
    
    // Problem if original -scan file had no image in it,
    // then m_a4nRect and m_a4nCircle are probably set wrong, so redo
    
    if ( (0 >= m_nDim0) || (0 >= m_nDim1) )
    {
        (void) oImage.nGetDimensions(&m_nDim0, &m_nDim1);
        m_bList = TRUE;
    }
    
    if (   (0 >= m_a4nRect[0])
        || (0 >= m_a4nRect[1])
        || (0 >= m_a4nRect[2])
        || (0 >= m_a4nRect[3]) )
    {
        // Default search rect is whole image less 1% border along each edge.
        
        m_a4nRect[0]    = m_nDim0 / 100;
        m_a4nRect[1]    = m_nDim1 / 100;
        m_a4nRect[2]    = m_nDim0 - m_a4nRect[0];
        m_a4nRect[3]    = m_nDim1 - m_a4nRect[1];
        m_bList = TRUE;
    }
    
    if (   (0 >= m_a4nCircle[0])
        || (0 >= m_a4nCircle[1]) )
    {
        // Default search circle is centered at mid-point of image
        // and includes whole image
        
        m_a4nCircle[0]  = m_nDim0 / 2;
        m_a4nCircle[1]  = m_nDim1 / 2;
        m_a4nCircle[2]  = 0;
        m_a4nCircle[3]  = (int) (sqrt( (double)( (m_a4nCircle[0] * m_a4nCircle[0])
            + (m_a4nCircle[1] * m_a4nCircle[1])
            ) ) );
        m_bList = TRUE;
    }
    
    m_nRadMin2         = m_a4nCircle[2] * m_a4nCircle[2];
    m_nRadMax2         = m_a4nCircle[3] * m_a4nCircle[3];
    
    //a2nBackground[0]   = max(m_nDim0 / 16, 32); // Default local background
    //a2nBackground[1]   = max(m_nDim1 / 16, 32); // tile is 1/256th of whole image
    a2nBackground[0]   = max(m_nDim0 / 16, 48); // Default local background
    a2nBackground[1]   = max(m_nDim1 / 16, 48); // tile is 1/256th of whole image

    
    // tjn:  Disable listing.  Why should we print the same info over and over again?
    // m_bList = TRUE;
    
   
    // The 'nearness' in this pass takes care of 'really close' spots.
    // A better per-image nearness check is done after we have computed spot shapes.
    fNear = 4.0f;
    
    if (m_bList)
    {
        // Some of the search parameters were updated, so inform user
        
        //      cout << "\nSome search parameters were changed, current values:\n";
        nList();
        m_bList = FALSE;
    }
    
    // Get some detector and source info to see if we can calculate
    // resolution
    
    Cdetector *poDetector;
    Csource   *poSource;
    float     a3fS0[3];
    float     fWavelength;
    
    
    if (NULL == poDetectorIn)
    {
        poDetector = new Cdetector(m_oHeader, "", TRUE, FALSE);
        bNewDetector = TRUE;
    }
    else
    {
        poDetector   = poDetectorIn;
        bNewDetector = FALSE;
    }
    if (NULL == poSourceIn)
    {
        poSource = new Csource(m_oHeader);
        bNewSource = TRUE;
    }
    else
    {
        poSource = poSourceIn;
        bNewSource = FALSE;
    }
    
    if (poSource->bIsAvailable())
    {
        poSource->vCalcGetS0(&a3fS0[0]);  // Required for resolution calc below
        fWavelength = poSource->fGetWavelength();
    }
    
    if (poDetector->bIsAvailable())
    {
        poDetector->vUpdateDDDN(); // Required for resolution calc below
    }

    poDetector->nCalcPixelShiftArray(*poSource,*m_poScan->m_poRotation,*poGoniometer);
    
    nArea[2] = a2nBackground[0];
    nArea[3] = a2nBackground[1];
    for (j = m_a4nRect[1]; j < m_a4nRect[3]; j = j + a2nBackground[1])
    {
        nArea[1] = j;
        for (i = m_a4nRect[0]; i < m_a4nRect[2]; i = i + a2nBackground[0])
        {
            
            fAvg = 0.0;      // Reset fAvg, fSD, fMin and fMax for this area's
            fSD  = 0.0;      //   local background calculation
            fMin = 0.0;
            fMax = fSatValue;
            
            if (0.0 < m_fSigma)
            {
                // Determine local background, if needed
                
                nArea[0] = i;
                
                nStat = nBackgroundAvgSD(oImage, nArea, fMin, fMax,
                    &fAvg, &fSD);
                
                
                for (k = 1; (k < m_nBLoops) && (0 == nStat); k++)
                {
                    // Re compute local background and SD, but exclude pixels
                    // 3 sigma from average (Do this twice more!)
                    // (MAYBE CHECK FOR CONVERGENCE)
                    
                    fMax = fAvg + max(0.01f, 3.0f * fSD);
                    nStat = nBackgroundAvgSD(oImage, nArea, fMin, fMax,
                        &fAvg, &fSD);
                }
                if (0 == nStat)
                {
                    fThreshold = fAvg + m_fSigma * max(0.01f, fSD);
                    fThreshold = max(fThreshold, m_fMinimum);
                }
            }
            else
            {
                nStat = 0;
                fThreshold = m_fMinimum;
            }
            
            if (0 == nStat)
            {
                // At this point we have a threshold we can test pixels against,
                // so ...
                
                for (l = j; l < j + a2nBackground[1]; l++)
                {
                    nRad1 = l - m_a4nCircle[1];
                    nRad1 = nRad1 * nRad1;
                    if (   (l >= m_a4nRect[1])
                        && (l <= m_a4nRect[3])
                        && (nRad1 <= m_nRadMax2) )
                    {
                        // l is within limits ...
                        for (k = i; k < i + a2nBackground[0]; k++)
                        {
                            nRad0 = k - m_a4nCircle[0];
                            nRad0 = nRad0 * nRad0  + nRad1;
                            if (   (k >= m_a4nRect[0])
                                && (k <= m_a4nRect[2])
                                && (nRad0 >= m_nRadMin2)
                                && (nRad0 <= m_nRadMax2) )
                            {
                                fValue = (oImage.*oImage.prfGetPixel)(k, l);
                                
                                // Threshold test, but do not look at saturated pixels
                                
                                if ( (fValue >= fThreshold) && ((fValue < fSatValue) || (m_bIncludeSat)) )
                                {
                                    // Neighbor check (pkcrit check!)
                                    
                                    nStat = 0;           // Used as a counter here
                                    fMax = fValue;
                                    
                                    for (n = l-1; (n < l+2) && (nStat >=0); n++)
                                    {
                                        for (m = k-1; (m < k+2) && (nStat >=0); m++)
                                        {
                                            fTemp = (oImage.*oImage.prfGetPixel)(m, n);
                                            if (fTemp > fValue)
                                            {
                                                // Leave loop
                                                nStat = -100;
                                            }
                                            else if (fTemp >= fThreshold)
                                            {
                                                // Increment the number of pixels above fTh.
                                                nStat++;
                                            }
                                        }
                                    }  // end of m loop
                                }  // end of n loop
                                
                                // In a 3x3 box the peak pixel must be the maximum
                                // value and there must be m_nPeakFilter pixels
                                // above fThreshold:
                                
                                if (nStat >= m_nPeakFilter)
                                {
                                    // Check if this refln is too near a refln
                                    // that has been added to the list in this
                                    // call to nFind2D.
                                    // Use m and nStat as temp vars.
                                    
                                    m      = poReflnlist->nGetNumReflns();
                                    nStat  = 0;
                                    for (n = 0; (n < m) && (0 == nStat); n++)
                                    {
                                        poRefln = poReflnlist->poGetRefln(n);
                                        if (((float)fabs((double)poRefln->fGetField(
                                                                                    poReflnlist->m_nFI_fObsPx0)
                                            - (double)k) < fNear)
                                            &&
                                            ((float)fabs((double)poRefln->fGetField(
                                            poReflnlist->m_nFI_fObsPx1)
                                            - (double)l) < fNear) )
                                        {
                                            // Rather than just flagging as too near,
                                            // save the stronger of the two near ones
                                            
                                            if (poRefln->fGetIntensity() > fValue)
                                            {
                                                // Flag as too near, so not inserted
                                                
                                                nStat = -101;
                                            }
                                            else
                                            {
                                                // Delete the one in the list and
                                                // reset counters since list is shorter
                                                
                                                poReflnlist->nDelete(n);
                                                n--;
                                                m--;
                                            }
                                        }
                                    }
                                    if (0 == nStat)
                                    {
                                        // Insert reflection into Reflnlist
                                        // for later consideration.
                                        
                                        oRefln.vSetField(poReflnlist->m_nFI_fObsPx0,
                                            (float)k + (float)0.5);
                                        oRefln.vSetField(poReflnlist->m_nFI_fObsPx1,
                                            (float)l + (float)0.5);
                                        oRefln.vSetField(poReflnlist->m_nFI_fObsRotWidth,
                                            (float) fInc);
                                        oRefln.vSetIntensity(fValue);
                                        
                                        if (poReflnlist->m_nFI_fObsRotSigma>=0) {
                                            // We don't know the rotation sigma... Just assume standard distribution.
                                            oRefln.vSetField(poReflnlist->m_nFI_fObsRotSigma,
                                                (float) (fInc/sqrt(12.0))
                                                );
                                        };
                                        oRefln.vSetSigmaI(fValue);
                                        
                                        // Remember that fRot is
                                        // EITHER
                                        //  the mid-point of rotation angle
                                        //  of this image.
                                        // OR
                                        //  the image number in a 3D scan
                                        
                                        oRefln.vSetField(poReflnlist->m_nFI_fObsRotMid, fRot);
                                        
                                        fDetReso = 0.0f;
                                        if (   poDetector->bIsAvailable()
                                            && poSource->bIsAvailable() )
                                        {
                                            
                                            fDetReso = fWavelength
                                                *
                                                poDetector->fCalcGetResolution((float)k,
                                                (float)l,
                                                a3fS0);
                                            
                                            oRefln.vSetField(poReflnlist->m_nFI_fDetResolution, fDetReso);
                                        }
                                        
                                        // Resolution test before insertion
                                        
                                        if (   (0.0 != fDetReso)
                                            && (m_fResolutionMin >= fDetReso)
                                            && (m_fResolutionMax <= fDetReso) )
                                            poReflnlist->nInsert(&oRefln);
                                    }
                                }
                                nStat = 0;  // Reset nStat
                            }  // k within limits
                        }  // end of k loop
                    }  // l within limits
        }  // end of l loop
        }
        else
        {
        /*      cout << "ERROR in Cfind::nFind2D, bad local background for area "
        << nArea[0] << ", " << nArea[1] << ", " << nArea[2]
        << ", " << nArea[3] << endl;
            */
            nStat = 0;
        }
    }  // end of i loop
    }  // end of j loop
    
    cout << poReflnlist->nGetNumReflns()
        << " preliminary spots found in 2D search with rotation angle " << fRot << " degs." << endl << flush;
    
#if (!defined(NO_X_WINDOWS))
#ifndef SSI_PC
    
    sSeqTmp = Cstring(nImgSeq++);
    
    if ( (0 == nStat) && (NULL != m_poXprop) )
    {
        // Write temporary list only if requested for updating display
        
        sReflnlistTmp = sGetCWD() + sDtrekGetPrefix() + "dtfind" + sSeqTmp
            + "tmp.ref";
        
        nStat = poReflnlist->nWrite(sReflnlistTmp);
        if (0 == nStat)
        {
            if (NULL == m_psXpropUpdate)
            {
                // No string with image name, so just update reflnlist
                
                m_poXprop->hSetProperty("DTDISPLAY_REFLN_UPDATE",
                    sReflnlistTmp);
            }
            else
            {
                m_poXprop->hSetProperty("DTDISPLAY_IMAGE_UPDATE",
                    *m_psXpropUpdate + " Reflnlist: "
                    + sReflnlistTmp);
            }
        }
    }
    
#else
    if (m_nDisplay) {
        // JME 2/27/01
        sReflnlistTmp = "dtfindtmp.ref";
        nStat = poReflnlist->nWrite(sReflnlistTmp);
        CCrclHelper::GetInstance()->vSendFindUpdateDisplay( oImage.sGetName(), sReflnlistTmp );
    }
    
#endif // SSI_PC
#endif // NO_X_WINDOWS
    
    // Are we using the 2D merge algorithm?
    if (bAddMergeInformation) {
        // Mark the all inserted reflections as NOT 'on edge'.
        int nFI_nMark = m_poReflnlist->nGetFieldIndex(ms_sOnFirstLastImage);
        for (j = 0; j < poReflnlist->nGetNumReflns(); j++) {
            (*poReflnlist)[j].vSetField(nFI_nMark,false);
        };
    };
    
    
    // Insert the list of newly found reflns into our main list.
    nStat = m_poReflnlist->nInsertListFrom(*poReflnlist);
    
    delete poReflnlist;
    delete poGoniometer;
    if (bNewDetector)
        delete poDetector;
    if (bNewSource)
        delete poSource;
           
    return (nStat);
}

int
Cfind::nAddGonioNumbers(Cimage& oImage) {
    int nx,ny;
    int nRots = 0;

    int a2nGonioSetStart[2];        // Starting index in reflection list for goniometer modifications.
    int a2nGonioSetEnd[2];          // Ending index in reflection list for goniometer modifications.
    Cgoniometer* a2poGonio[2];      // Goniometer to apply to the range specified by a2nGonioSetStart and a2nGonioSetEnd
    int a2nGonioRotAxis[2];         // Goniometer rotation axis.
    int nGonioSets = 0;             // Number of goniometer sets to modify.
    int nGonioSet;


    if (!m_poGoniometer) {
        m_poGoniometer = new Cgoniometer(oImage.m_oHeader,(Cstring) D_K_CrystalPrefix);
        m_poRotation = new Crotation(oImage.m_oHeader);
        if ((!m_poGoniometer->bIsAvailable()) || (!m_poRotation->bIsAvailable()))
            return 1;
        return 0;
    };
    Cgoniometer oGonio(oImage.m_oHeader,(Cstring) D_K_CrystalPrefix);
    Crotation oRotation(oImage.m_oHeader);

    if ((!oGonio.bIsAvailable()) || (!oRotation.bIsAvailable()))
        return 1;

    // See if the goniometer constants are in the reflection list
    

        nRots = 0;
        for (nRots = 0; nRots < m_poGoniometer->nGetNumValues(); nRots++)
          {
            if ("deg" != m_poGoniometer->m_psUnits[nRots])
              break;
          };

    if (m_poReflnlist->m_nFI_fGonio1<0)
      {
        for (nx=0;nx< nRots;nx++)
          {
            if (ABS(m_poGoniometer->m_pfDatumValue[nx] - oGonio.m_pfDatumValue[nx])>0.01)
              break;
          }
        // If we made it all the way through, there is nothing to do.
        if (nx == nRots)
            return 0;

        nx = nRots;
        // Add the requested fields to the reflection list.
        if (nx-->0)
            m_poReflnlist->nExpandGetField(Creflnlist::ms_sfGonio1);
        if (nx-->0)
            m_poReflnlist->nExpandGetField(Creflnlist::ms_sfGonio2);
        if (nx-->0)
            m_poReflnlist->nExpandGetField(Creflnlist::ms_sfGonio3);
        if (nx-->0)
            m_poReflnlist->nExpandGetField(Creflnlist::ms_sfGonio4);
        if (nx-->0)
            m_poReflnlist->nExpandGetField(Creflnlist::ms_sfGonio5);

        m_poReflnlist->nExpandGetField(Creflnlist::ms_snGonioRotAxis);

        a2nGonioSetStart[nGonioSets] = 0;
        a2nGonioSetEnd[nGonioSets] = m_nFirstReflnAtLast2DFind-1;
        a2poGonio[nGonioSets] = m_poGoniometer;
        a2nGonioRotAxis[nGonioSets] = m_poGoniometer->nGetNum(m_poRotation->sGetName());
        if (a2nGonioRotAxis[nGonioSets]<0)
            return 1;
        nGonioSets++;
      }

    a2nGonioSetStart[nGonioSets] = m_nFirstReflnAtLast2DFind;
    a2nGonioSetEnd[nGonioSets] = m_poReflnlist->nGetNumReflns()-1;
    a2poGonio[nGonioSets] = &oGonio;
    a2nGonioRotAxis[nGonioSets] = oGonio.nGetNum(oRotation.sGetName());
    if (a2nGonioRotAxis[nGonioSets]<0)
        return 1;
    nGonioSets++;

    
    // If we get here, then we need to add goniometer information to the header.
    for (nGonioSet = 0; nGonioSet < nGonioSets; nGonioSet++) {
        for (nx=a2nGonioSetStart[nGonioSet];nx<=a2nGonioSetEnd[nGonioSet];nx++) {
            ny = 0;
            if (ny++ < nRots)
                (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_fGonio1,(float) a2poGonio[nGonioSet]->m_pfDatumValue[ny-1]);
            if (ny++ < nRots)
                (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_fGonio2,(float) a2poGonio[nGonioSet]->m_pfDatumValue[ny-1]);
            if (ny++ < nRots)
                (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_fGonio3,(float) a2poGonio[nGonioSet]->m_pfDatumValue[ny-1]);
            if (ny++ < nRots)
                (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_fGonio4,(float) a2poGonio[nGonioSet]->m_pfDatumValue[ny-1]);
            if (ny++ < nRots)
                (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_fGonio5,(float) a2poGonio[nGonioSet]->m_pfDatumValue[ny-1]);
            if (m_poReflnlist->m_nFI_nGonioRotAxis>=0)
                (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_nGonioRotAxis,a2nGonioRotAxis[nGonioSet]+1);
        };
    };

    return 0;
};

bool
Cfind::bPixelInLimits(const int nPx0, const int nPx1)
{
  if (nPx0 < m_a4nRect[0]) return (FALSE);
  if (nPx1 < m_a4nRect[1]) return (FALSE);
  if (nPx0 > m_a4nRect[2]) return (FALSE);
  if (nPx1 > m_a4nRect[3]) return (FALSE);
  int nTemp0, nTemp1;
  nTemp0 = nPx0 - m_a4nCircle[0];
  nTemp0 = nTemp0 * nTemp0;
  nTemp1 = nPx1 - m_a4nCircle[1];
  nTemp1 = nTemp1 * nTemp1;
  nTemp0 = nTemp1 + nTemp0;
  if (nTemp0 < m_nRadMin2) return (FALSE);
  if (nTemp0 > m_nRadMax2) return (FALSE);
  return (TRUE);
}

int
Cfind::nBackgroundAvgSD(Cimage& oImage, const int nPxArea[4],
                             const float fExcludeMin,
                             const float fExcludeMax,
                             float *pfAvg, float *pfSD)
{

// Compute background and standard deviation in a portion of an image.
// Exclude from the calculation pixels with values <= fExcludeMin and
// => fExcludeMax.  Look only in the the image area defined by nPxArea.
// Exclude pixels flagged as bad by the non-uniformity object and that are
// outside the find limits (m_a4nRect and m_a4nCircle).
//
  int   i, j;
  int   nStat;
  float fValue;
  float fSum1, fSum2, fNum;
  int   nRad0, nRad1;
  //+2011-Nov-14 JWP
  float fBadValue = 0.0;
  if  ( (eImage_I4 == oImage.nGetDataType()) || (eImage_realIEEE == oImage.nGetDataType()))
    {
      fBadValue = -1.0;
    }
  //-2011-Nov-14 JWP
  fSum1 = 0.0;
  fSum2 = 0.0;
  fNum  = 0.0;

  for (j = nPxArea[1]; j < nPxArea[1]+nPxArea[3]; j++)
    {
      nRad1 = j - m_a4nCircle[1];
      nRad1 = nRad1 * nRad1;
      if (   (j >= m_a4nRect[1])
          && (j <= m_a4nRect[3])
          && (nRad1 <= m_nRadMax2) )
        {
          // j is within limits ...

          for (i = nPxArea[0]; i < nPxArea[0]+nPxArea[2]; i++)
            {
              nRad0 = i - m_a4nCircle[0];
              nRad0 = nRad0 * nRad0  + nRad1;
              if (   (i >= m_a4nRect[0])
                  && (i <= m_a4nRect[2])
                  && (nRad0 >= m_nRadMin2)
                  && (nRad0 <= m_nRadMax2) )
                {
                  fValue = (oImage.*oImage.prfGetPixel)(i, j);
                  //+-2011-Nov-14 JWP next line 
                  if (fBadValue < fValue)
                    {
                      // Pixel is in limits and is not bad

                      if ( (fValue >= fExcludeMin) && (fValue <= fExcludeMax) )
                        {
                          // This pixel passes all tests to include in the
                          //  background and sd calculation.
                          fSum1 = fSum1 + fValue;
                          fSum2 = fSum2 + (fValue * fValue);
                          fNum  = fNum + 1.0f;
                        }
                    }
                }
            }  // end of i loop
        }
    }  // end of j loop

  // For a valid average and sd to be calculated we must have at least 1
  // pixel in the average and more than 1/2 of the pixels in the specified
  // area must contribute!

  if ( (fNum > 1.0) && (fNum > (nPxArea[2] * nPxArea[3] / 2)) )
    {
      *pfAvg = fSum1 / fNum;
      fSum2  = (fSum2 - (fSum1*fSum1 / fNum)) / (fNum - 1.0f);
      if (0.0 <= fSum2)
        {
          *pfSD = (float)sqrt((double)fSum2);
          nStat = 0;
        }
      else
        {
          *pfSD = (float)sqrt((double)-fSum2);
          cout << "Cfind::nBackgroundAvgSD, sqrt of negative!\n";
          nStat = 1;
        }
    }
  else
    {
//    cout << "Cfind::nBackgroundAvgSD, not enough pixels to compute!\n";
      nStat = 2;
    }
  return (nStat);
}

int
Cfind::nExpandReflnlist(Creflnlist *poReflnlistIn,bool bExpandExtra,bool bAddEllipsoidInformation)
{

// Add required fields to the reflection list if they are not already there.
// The fields are named by the static member variables in the Creflnlist class.

  Creflnlist *poReflnlistL;   // Local reflection list

  if (NULL == poReflnlistIn)
    {
      poReflnlistL = m_poReflnlist;     // Use on member variable
    }
  else
    {
      poReflnlistL = poReflnlistIn;   // Use passed reflection list
    }

  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPx0);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPx1);
  //(void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPxPeak0);
  //(void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPxPeak1);
  
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsSharpness);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPeakAreaCount);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsPeakBorderCount);
  
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsRotMid);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsRotWidth);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfObsRotSigma);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfDetResolution);

  if (bAddEllipsoidInformation) {
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidA00);
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidA01);
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidA11);
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidb0);
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidb1);
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfEllipsoidc);
  };
  
  
  if (bExpandExtra) {
      (void) poReflnlistL->nExpandGetField(ms_sMergeHead);           // For the sake of nOverlapCheck()
      (void) poReflnlistL->nExpandGetField(ms_sMergeNext);           // For the sake of nOverlapCheck()
      (void) poReflnlistL->nExpandGetField(ms_sIndex);               // For the sake of nOverlapCheck()
      (void) poReflnlistL->nExpandGetField(ms_sOnFirstLastImage);    // For the sake of nMerge2DReflections();
  };

  if (m_bIncludeBackground) {
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfBackground);
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfBackgroundSigma);
  };

  // Add the Gonio fields if they are already in m_poReflnlist.
  if ((poReflnlistIn) && (m_poReflnlist->m_nFI_fGonio1>=0)) {
      poReflnlistL->nExpandGetField(poReflnlistL->ms_snGonioRotAxis);
      poReflnlistL->nExpandGetField(poReflnlistL->ms_sfGonio1);
  };
  if ((poReflnlistIn) && (m_poReflnlist->m_nFI_fGonio2>=0))
      poReflnlistL->nExpandGetField(poReflnlistL->ms_sfGonio2);
  if ((poReflnlistIn) && (m_poReflnlist->m_nFI_fGonio3>=0))
      poReflnlistL->nExpandGetField(poReflnlistL->ms_sfGonio3);
  if ((poReflnlistIn) && (m_poReflnlist->m_nFI_fGonio4>=0))
      poReflnlistL->nExpandGetField(poReflnlistL->ms_sfGonio4);
  if ((poReflnlistIn) && (m_poReflnlist->m_nFI_fGonio5>=0))
      poReflnlistL->nExpandGetField(poReflnlistL->ms_sfGonio5);


  return (0);
}

Cstring Cfind::ms_sMergeHead                    = "nMergeHead";       
Cstring Cfind::ms_sMergeNext                    = "nMergeNext";       
Cstring Cfind::ms_sIndex                        = "nIndex";               
Cstring Cfind::ms_sOnFirstLastImage             = "nOnFirstLastImage";


int
Cfind::nList(void)
{
  cout << "Find object listing:"
       << "\n     Sigma: " << m_fSigma
       << "\nResolution: " << m_fResolutionMin << " to " << m_fResolutionMax
       << "\n   Minimum: " << m_fMinimum
       << "\nCircle lim: " << m_a4nCircle[0] << ", " << m_a4nCircle[1] << ", "
                           << m_a4nCircle[2] << ", " << m_a4nCircle[3]
       << "\n  Rect lim: " << m_a4nRect[0] << ", " << m_a4nRect[1] << ", "
                           << m_a4nRect[2] << ", " << m_a4nRect[3]
       << "\nSpot wind.: " << m_a2nSpotWindow[0] << ", " << m_a2nSpotWindow[1]
       << "\nPeak filt.: " << m_nPeakFilter
       << "\nBack. tile: " << (m_nDim0/16) << ", " << (m_nDim1/16)
       << "\n Seq. num.: " << m_a2nSeqNum[0] << ", " << m_a2nSeqNum[1]
       << "\n Every Nth: " << m_nEveryNthImage
       << "\nImage dim.: " << m_nDim0 << ", " << m_nDim1
       << "\n   3D dump: " << m_nDumpRefln
       << endl << flush;
  
        if( m_bPeakInfo && m_oDetectorAreasIn.bHasPrintInfo() ) 
          {
            m_oDetectorAreasIn.nPrintPeakAreas(false,1.0);
          }
  
  return (0);
}

int Cfind::nFindCentroid2D(Cimage& roImage, 
                                                   Crefln* poRefln, 
                                                   Cdetector* poDetector, 
                                                   DTREK_WORD wCtrl)
{
    // Find the center of gravity for a spot described by poRefln in
    // image roImage.  This is a 2D method!
    // Saturated and otherwise bad pixels are simply ignored in the calculation.
    // Place centroid, Intensity, SigmaI into *poRefln.

  int   nStat = 0;
  int   nx = 0;
  int   i= 0;
  float fSD = 0.0f;
  float fB = 0.0f;
  float fBSD = 0.0f;
  float fInt = 0.0f;
  
  float a3fCentroid[3] = {0.0f};
  float a3fSize[3] = {0.0f};
  float fInc = 0.0;
  int   nArea[4] = {0};
  int   nDim0 = 0, nDim1 = 0;     // Dimensions of the image
  
  C3DdataDetectorArea oArea;

  fInc     = m_poScan->m_poRotation->fGetIncrement();
  if (0.0 == fInc) {
      // Prevent divide by zero later
      fInc = m_fDefaultRockingWidth;
  };



  if (   (0 > poRefln->m_poReflnlist->m_nFI_fObsPx0)
      || (0 > poRefln->m_poReflnlist->m_nFI_fObsPx1)
      || (0 > poRefln->m_poReflnlist->m_nFI_fObsRotMid) 
      )
    {
      cout << "ERROR Cfind reflection list missing necessary fields!\n";
      return (-2);
    }

  roImage.nGetDimensions(&nDim0, &nDim1);

  if (! m_oDetectorAreas.bIsInitialized()) {
      m_oDetectorAreas.nInitPeakAreas(nDim0,nDim1);
  };

  a3fCentroid[0] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fObsPx0);
  a3fCentroid[1] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fObsPx1);
  a3fCentroid[2] = poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fObsRotMid);

  poRefln->vSetIntensity((float) -999.0); //Signal to delete later

  // See if we can obtain area information.
  

  if ( (0 >= m_a2nSpotWindow[0]) || (0 >= m_a2nSpotWindow[1]) )
  {
      m_oDetectorAreasIn.nGetPeakArea((int) a3fCentroid[0],(int) a3fCentroid[1],oArea);      
      nArea[2] = (int) min(nMAXSPOTSIZE,3.5*oArea.fSize[0]);
      nArea[3] = (int) min(nMAXSPOTSIZE,3.5*oArea.fSize[1]);
  }
  else
  {
      nArea[2] = m_a2nSpotWindow[0];
      nArea[3] = m_a2nSpotWindow[1];
  }
  nArea[0] = (int) a3fCentroid[0] - (nArea[2] / 2);
  nArea[1] = (int) a3fCentroid[1] - (nArea[3] / 2);

  // Limit the image area to within the image

  nArea[0] = max(0, nArea[0]);
  nArea[1] = max(0, nArea[1]);
  nArea[2] = min(nArea[2], nDim0 - nArea[0] + 1);
  nArea[3] = min(nArea[3], nDim1 - nArea[1] + 1);

  a3fCentroid[2] = 0.5;

  if ((nArea[2]<=0) || (nArea[3]<=0) || (nArea[0]>=nDim0) || (nArea[1]>=nDim1))
      return 1;

  C3Ddata o3DdataLocal(nArea[2], nArea[3], 1, nArea[0], nArea[1], 0); 
  
  // Create C3Ddata object

  m_oShoebox.m_nExt[0]        = nArea[2];
  m_oShoebox.m_nExt[1]        = nArea[3];
  m_oShoebox.m_nExt[2]        = 1;
  m_oShoebox.m_nOrigOffset[0] = nArea[0];
  m_oShoebox.m_nOrigOffset[1] = nArea[1];
  m_oShoebox.m_nOrigOffset[2] = 0;
  m_oShoebox.vCheckMemory();

  C3Ddata& o3Ddata = o3DdataLocal; // We once constructed a local shoebox here every time.

  

  // Fill the 3Ddata object

  nStat = o3Ddata.nFill2D(&roImage,0);

  // Don't allow any zero values in the shoebox.  
//+JWP 2009-07-18
//+JWP 2013-05-12  Allow realIEEE to have zero value pixels, too.
// TODO: This causes a problem with Pilates "long int" images 
//       that have a value of 0 as OK and -1, -2 for not OK images.
  float fBadValue = 0.0;
  long nNumZeros;
  double dNumZeroFraction;
  if ( (eImage_I4 == roImage.nGetDataType()) || (eImage_realIEEE == roImage.nGetDataType()) )
    {
      fBadValue = -1.0;
    }
  o3Ddata.nCountZeros(&nNumZeros, fBadValue);
  dNumZeroFraction =  (double)nNumZeros / (double)(nArea[2] * nArea[3]);
  //cout << "Size of window: " << nArea[2] << " x " << nArea[3] << endl;
  //cout << "Num zeros, % : " << nNumZeros << ", " << 100.0 * dNumZeroFraction; << endl;
  if (0.2 < dNumZeroFraction)
    return 1; // Too many zeros 
  if (m_bIncludeSat)
      o3Ddata.m_bIncludeSat = TRUE;

  // Convert to float, was already corrected for nonunf previously
  nStat = o3Ddata.nConvertToFloat();
 
  if (0 == nStat)
    {
      // Conversion was successful, so ...
      // Calculate the centroid of a peak in the shoebox given a starting
      // pixel near the peak center


      a3fCentroid[0] -= (float) o3Ddata.nGetOffset(0);
      a3fCentroid[1] -= (float) o3Ddata.nGetOffset(1);
      a3fCentroid[2] -= (float) o3Ddata.nGetOffset(2);


      
      fSD   = 0.0;
      fInt  = -1.0;

      nStat = o3Ddata.nCalcGetPeakInfo(a3fCentroid, &fInt, &fSD,&fB,&fBSD);

      // RB: reject the spot if the peak maximum is not high enough above the background
      if( m_fPeakMinimum > 0.0 && o3Ddata.fGetPeakMaxValue() - o3Ddata.fGetPerPixelBackground() < m_fPeakMinimum )
          return 1; 

      // cout << "Status from o3Ddata.nCalcGetPeakInfo: " << nStat << '\n';

      // Mask out on edge error in 3rd dimension, since this is a 2D search
      nStat = nStat & ~((1 << Crefln::ms_nErrorOnEdge2) | (1 << Crefln::ms_nErrorUnaccountedPix));

      if( nStat == (1 << Crefln::ms_nErrorTooFewPix) && wCtrl & DTREK_DTFIND_FC2D_LOOK_FOR_WEAK_SPOTS )
      {
          // Try finding the weak spot.  
          // This is important because we might be finding weak spots when used for delta-time strategy.
          // Use the previous peak ellipsoid as a mask.

          // We must provide a standard ellipsoid.

          o3Ddata.m_a2x2fEllipsoidAIn[0][0] = oArea.a2x2fEllipsoidA[0][0];
          o3Ddata.m_a2x2fEllipsoidAIn[0][1] = oArea.a2x2fEllipsoidA[0][1];
          o3Ddata.m_a2x2fEllipsoidAIn[1][0] = oArea.a2x2fEllipsoidA[1][0];
          o3Ddata.m_a2x2fEllipsoidAIn[1][1] = oArea.a2x2fEllipsoidA[1][1];
          o3Ddata.m_a2fEllipsoidbIn[0] =   oArea.a2fEllipsoidb[0];
          o3Ddata.m_a2fEllipsoidbIn[1] =   oArea.a2fEllipsoidb[1];
          o3Ddata.m_fEllipsoidcIn = oArea.fEllipsoidc;
          
          nStat = o3Ddata.nCalcGetPeakInfo(a3fCentroid, &fInt, &fSD,&fB,&fBSD,NULL);
          
          nStat = nStat & (~(1 << Crefln::ms_nErrorTooFewPix));  // mask out the "too few pixels" error
      }

      // If centroid was not really calculated, reject the spot

      if (0.0 == a3fCentroid[0])
          nStat = 1;

      /// RB A slight change: Do not reject a weak spot, if a DTREK_DTFIND_FC2D_REPORT_WEAK_SPOTS flag is set
      
      //// If centroid is below user specified sigma cut off, reject the spot.
      // if ( (fSD==0.0) || ( (fInt/fSD<m_fSigma) && !(wCtrl & DTREK_DTFIND_FC2D_LOOK_FOR_WEAK_SPOTS) ) )
      //    nStat = 1;
     
      if( fSD == 0.0 )
      {
            if ( !(wCtrl & DTREK_DTFIND_FC2D_REPORT_WEAK_SPOTS) )
                nStat = 1;
            else 
                fSD = 1.0; // RB: An arbitrary value to avoid dividing by zero. Later can calculate a real sigma value from bkg info
      }
      else if ( fInt/fSD < m_fSigma && !(wCtrl & DTREK_DTFIND_FC2D_LOOK_FOR_WEAK_SPOTS) )
        nStat = 1;
      
      
      if( 0 == nStat || wCtrl & DTREK_DTFIND_FC2D_REPORT_BKG_INFO_DISREGARD_ERROR )
      {
          if ( m_poReflnlist->m_nFI_fBackground >= 0 )
              poRefln->vSetField(m_poReflnlist->m_nFI_fBackground, (float)fB);

          if ( m_poReflnlist->m_nFI_fBackgroundSigma >= 0 )
              poRefln->vSetField(m_poReflnlist->m_nFI_fBackgroundSigma, (float)fBSD);
      }

      if( 0 == nStat )
      {
          // Transfer results to poRefln
          a3fCentroid[0] += (float) o3Ddata.nGetOffset(0);
          a3fCentroid[1] += (float) o3Ddata.nGetOffset(1);
          a3fCentroid[2] += (float) o3Ddata.nGetOffset(2);

          m_oDetectorAreas.nAddPeakArea((int) a3fCentroid[0],(int) a3fCentroid[1],&o3Ddata);
          
          poRefln->vSetIntensity(fInt);
          poRefln->vSetSigmaI(fSD);

          poRefln->vSetField(m_poReflnlist->m_nFI_fObsPx0,
              a3fCentroid[0]);
          poRefln->vSetField(m_poReflnlist->m_nFI_fObsPx1,
              a3fCentroid[1]);

          poRefln->vSetField(m_poReflnlist->m_nFI_fObsSharpness,
              (float) o3Ddata.fGetPeakSharpness());
          
          if ((m_poReflnlist->m_nFI_fEllipsoidA00>=0) &&
              (m_poReflnlist->m_nFI_fEllipsoidA01>=0) &&
              (m_poReflnlist->m_nFI_fEllipsoidA11>=0) &&
              (m_poReflnlist->m_nFI_fEllipsoidb0>=0) &&
              (m_poReflnlist->m_nFI_fEllipsoidb1>=0) &&
              (m_poReflnlist->m_nFI_fEllipsoidc>=0)) {
              // Transform so that we are using the observed centroid as the reference point.
              float a2x2fEllipsoidA[2][2];
              float a2fEllipsoidb[2];
              float fEllipsoidc;
              o3Ddata.vGetPerSliceEllipsoidRelObsCent(0,&a2x2fEllipsoidA[0][0],&a2fEllipsoidb[0],&fEllipsoidc);
              poRefln->nPutGetEllipse(true,&a2x2fEllipsoidA[0][0],&a2fEllipsoidb[0],&fEllipsoidc);
          };


          
          if (m_poReflnlist->m_nFI_fObsRotSigma>=0) {
              // We don't know the rotation sigma... Just assume standard distribution.
              poRefln->vSetField(m_poReflnlist->m_nFI_fObsRotSigma,
                  (float) (fInc/sqrt(12.0))
                  );
          };        
          o3Ddata.vGetPeakSize(a3fSize);         
          
//+JWP 2011-11-27
// If m_a2nPMTInfo[0] is not 0, then do something special:
//   Sort pixel values in o3Ddata, and return the average of the highest m_a2nPMTInfo[0] pixel values
//   in fIntensity
          if (0 < m_a2nPMTInfo[0])
            {
              // Copy into new C3Ddata before sort?
              C3Ddata o3Dcopy(o3Ddata);
              (void) o3Dcopy.nFill2D(&roImage,0);
              float fTopValues = 0.0;
              float *pfo3DcopyValues = o3Dcopy.m_pfData;
              int   nNumToUse = m_a2nPMTInfo[0];
              int   nNumPoints = o3Dcopy.nGetExt(0) * o3Dcopy.nGetExt(1) * o3Dcopy.nGetExt(2);
              qsort(&pfo3DcopyValues[0], nNumPoints, sizeof(float), float_cmp);
 
              if (nNumPoints < m_a2nPMTInfo[0]) nNumToUse = nNumPoints;
              for (i = 0; i < nNumToUse; i++)
                {
                  // Threshold check here?
                  // Sort is low to high, so go backwards through the array
                  fTopValues += pfo3DcopyValues[nNumPoints-1-i];
                }
              if (0 < nNumToUse)
                fTopValues = fTopValues / (float) nNumToUse;
              // We may wish to put the fTopValues in a different field
              poRefln->vSetField(m_poReflnlist->m_nFI_fObsPeakBorderCount, (float)fTopValues);
            }
//-JWP 2011-11-27

          m_nTotalSpotsProcessed++;
          o3Ddata.vGetSatCount(&nx);
          if (nx)
            m_nTotalSatSpotsProcessed++;

                  if(wCtrl & DTREK_DTFIND_FC2D_GET_EXTRA_INFO)
                  {
                                REFLN_EXTRA_INFO*               pstExtraInfo = new REFLN_EXTRA_INFO;    
                                pstExtraInfo->vInitialize();

                                poRefln->m_pvUserData = (void*)pstExtraInfo;
                          
                if(wCtrl & DTREK_DTFIND_FC2D_GET_EXTRA_INFO_SAT_COUNT)    
                    pstExtraInfo->m_nSaturatedCount = nx;
                
                if(wCtrl & DTREK_DTFIND_FC2D_GET_EXTRA_INFO_PEAK_MAX)    
                                    pstExtraInfo->m_dMaxPixelValue = o3Ddata.fGetPeakMaxValue();
                     
                if(wCtrl & DTREK_DTFIND_FC2D_GET_EXTRA_INFO_ELLIPSE_COUNT)    
                {
                                    // Get number of pixels in the peak ellipse
                                    int         nMaskArrayLength = o3Ddata.nGetExt(0) * o3Ddata.nGetExt(1);
                                    int*        pnMask = o3Ddata.pnGetMask();
                                    for(int ii=0; ii < nMaskArrayLength; ii++)
                                    {
                                            if( pnMask[ii] & C3Ddata::ms_nEllipsoidMask )
                                                    pstExtraInfo->m_nPeakEllipsoidPixelCount++;
                                    }
                }
          }

          if( wCtrl & DTREK_DTFIND_FC2D_GET_PEAK_AREA_INFO )
          {
              int       nPeakAreaPixelCount   = 0;
              int       nPeakEllipsPixelCount = 0;
              int       nPeakBorderPixelCount = 0;
              
              o3Ddata.vGet2DPeakAreaPixelsCount(nPeakAreaPixelCount, nPeakEllipsPixelCount, nPeakBorderPixelCount);

              poRefln->vSetField(m_poReflnlist->m_nFI_fObsPeakAreaCount,   (float)nPeakAreaPixelCount);
              poRefln->vSetField(m_poReflnlist->m_nFI_fObsPeakBorderCount, (float)nPeakBorderPixelCount);
          }

          o3Ddata.vGetPeakSize(a3fSize);
          /*
          // tjn:  These have been removed:  we are now using the ellipsoid constants if they are available.         
          if (m_poReflnlist->m_nFI_fObsPx0Width>=0) {
              poRefln->vSetField(m_poReflnlist->m_nFI_fObsPx0Width,
                  (float) a3fSize[0]);
          };
          if (m_poReflnlist->m_nFI_fObsPx1Width>=0) {
              poRefln->vSetField(m_poReflnlist->m_nFI_fObsPx1Width,
                  (float) a3fSize[1]);
          };
          */
                              
      }
      else
      {
          nStat = 1;
      }
  }
  
  return (nStat);
}

int Cfind::nMerge2D_OverlapReject(Creflnlist& oReflnlist,
                                  itr<int>&   anRefsIn,
                                  Crefln&     oNewRefln,
                                  Creflnlist& oNewList)
{
    if( 0 == anRefsIn.size() )
        return 0; // nothing to do
    
    bool    bDuplicateOverlap;
    int     nRefSort;
    int     nRefSort2;
    int     nRef,nRef2,nRef3;
    int     nMaxContrib = 0;
    double  fRotStep;
    
    double  f0,f1;
    int     nx;
    
    itr<double> afRotMid;
    itr<int>    anRefIndex;
    itr<int>    anCompete;
    itr<int>    anSourceRefs;    // Used as temporary with nGetMapping()
   
    
    fRotStep = m_poScan->m_poRotation->fGetIncrement();
  
    -anRefIndex;
    -anCompete;
        
    for(nRefSort = 0; nRefSort < anRefsIn.size();nRefSort++)
    {
        nRef2 = anRefsIn[nRefSort];
        
        // Check for a 'bad' overlap.  This happens when two 2D centroids having the same observed rotation midpoint,
        // are significantly overlapping (i.e. one ellipsoid covers the other's centroid)
        // Also, check for a 'duplicate' overlap.  This happens when a spot has 'double peaks' that
        // both get picked up in the first pass 2D find.  (i.e. both reflections' ellipsoids cover the others' centroids)
        // This can in fact happen when nFind2D() finds two nearby peaks and nFindCentroid2D() refines them to almost the same position.
        
        bDuplicateOverlap = false;
        for(nRefSort2 = nRefSort+1; nRefSort2 < anRefsIn.size(); nRefSort2++)
        {
            nRef3 = anRefsIn[nRefSort2];
            
            if( oReflnlist[nRef2].fGetField(oReflnlist.m_nFI_fObsRotMid) == oReflnlist[nRef3].fGetField(oReflnlist.m_nFI_fObsRotMid) )
            {
              f0 = oReflnlist[nRef2].fGetIntensity();
              f1 = oReflnlist[nRef3].fGetIntensity();
              bDuplicateOverlap = true;
            }
        }
        
        if( !bDuplicateOverlap )
        {

            anRefIndex + nRef2;
     
            afRotMid + oReflnlist[nRef2].fGetField(oReflnlist.m_nFI_fObsRotMid);

            if( oReflnlist.m_nFI_nSourceRefNum >= 0 )
            {
                if( oReflnlist[nRef2].nGetField(oReflnlist.m_nFI_nSourceRefNum) >= 0 )
                    anCompete + 1;
                else
                    anCompete + 0;
            } 
            else
                anCompete + 0;
        }
    }
    
    if( 0 == anRefIndex.size() )
        return 0;
    
    ///////////////////////
    // DO THE ACTUAL MERGE
    ///////////////////////

    double  a3fCent[3];
    double  fIntensity = 0.0;
    double  fVarianceNumer = 0.0;
    double  fVariance = 0.0;
    double  fMaxRot = 0.0;
    double  fMinRot = 0.0;
    double  fMax = 0.0;
    
    double fSharpness = 0.0;
    double fPeakAreaCount = 0.0;
    double fPeakBorderCount = 0.0;

    double a7fEllipsoid[7],a7fEllipsoidNext[7],a7fEllipsoidTemp[7];
    double a2fCenter[2],a2fCenterNext[2];

    //////////////////////////////////////////////////////////////////////////
    // Sort based on the rotation mid point.
    itr<int>    anRotMidSort;
    -anRotMidSort;
    for(nx=0; nx < afRotMid.size(); nx++)
        anRotMidSort + nx;
    
    g_pfCmpDoubles = &afRotMid[0];
    
    if( anRotMidSort.size() > 0 )
        qsort(&anRotMidSort[0], anRotMidSort.size(), sizeof(int), double_cmp_rel);
    //////////////////////////////////////////////////////////////////////////
    
    int     nContrib = 0;
    bool    bNewRefln = false;

    for(nRefSort = 0; nRefSort <= anRotMidSort.size(); nRefSort++) 
    {
        bNewRefln = false;
        
        if( nRefSort == 0 || nRefSort == anRotMidSort.size() )
            bNewRefln = true;
        else if( ABS(afRotMid[anRotMidSort[nRefSort]] - afRotMid[anRotMidSort[nRefSort - 1]] - fRotStep) > 0.001 )
            bNewRefln = true;
        else if( anCompete[anRotMidSort[nRefSort]] || anCompete[anRotMidSort[nRefSort-1]] )
            bNewRefln = true;

        if( bNewRefln ) 
        {
            // Write out previous reflection.
            if( nContrib > 0 )
            {
                // Convert the ellipsoid to the new centroid coordinates.
                vConvertEllipsoidRelative(&a7fEllipsoid[0],&a7fEllipsoid[4],&a7fEllipsoid[6],a2fCenter[0],a2fCenter[1],a3fCent[0]/fIntensity,a3fCent[1]/fIntensity);
                oNewRefln.nPutGetEllipse(true,&a7fEllipsoid[0],&a7fEllipsoid[4],&a7fEllipsoid[6]);
                
                // Build reflection information
                oNewRefln.vSetIntensity(fIntensity);
                oNewRefln.vSetSigmaI(sqrt(fVariance));
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsPx0,
                    (float) (a3fCent[0]/fIntensity));
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsPx1,
                    (float) (a3fCent[1]/fIntensity));
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsSharpness,
                    (float)fSharpness);
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsPeakAreaCount,
                    (float)fPeakAreaCount);
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsPeakBorderCount,
                    (float)fPeakBorderCount);
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsRotMid,
                    (float) (a3fCent[2]/fIntensity));
                
                oNewRefln.vSetField(oNewList.m_nFI_fDetResolution,
                    (float) oReflnlist[anRefsIn[0]].fGetField(oReflnlist.m_nFI_fDetResolution));
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsRotWidth,
                    (float) (2.0*max(fMaxRot - a3fCent[2]/fIntensity,a3fCent[2]/fIntensity - fMinRot)));
                
                oNewRefln.vSetField(oNewList.m_nFI_fObsRotSigma,
                    (float) (fRotStep * sqrt(fVarianceNumer/(max(1.0,fIntensity) * max(1.0, fIntensity)*12))));
                
                if ((oNewList.m_nFI_nCompeteRefNum >= 0) && (oNewList.m_nFI_nSourceRefNum >= 0) && (oReflnlist.m_nFI_nSourceRefNum >= 0)) 
                {
                    if ((nContrib == 1) && (anCompete[anRotMidSort[nRefSort-1]])) 
                        oNewRefln.vSetField(oNewList.m_nFI_nSourceRefNum,oReflnlist[anRefIndex[anRotMidSort[nRefSort-1]]].nGetField(oReflnlist.m_nFI_nSourceRefNum));
                    else
                        oNewRefln.vSetField(oNewList.m_nFI_nSourceRefNum,-1);
                }

                oNewList.nInsert(&oNewRefln);               
            }

            // Zero variables.
            vZeroMat(3, 1, a3fCent);
            fIntensity = 0.0;
            fVariance = 0.0;
            fVarianceNumer = 0.0;
            nContrib = 0;
            fMaxRot = -9999.0;
            fMinRot = -9999.0;
            fMax = 0.0;
            fSharpness = 0.0;

            fPeakAreaCount = 0.0;
                
            fPeakBorderCount = 0.0;

            a7fEllipsoid[0] = 0.0;
        }

        if (nRefSort < anRotMidSort.size())
            nRef = anRefIndex[anRotMidSort[nRefSort]];
        else
            break;

        // Add in new reflection information.
        f0 = oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsRotMid) + fRotStep / 2.0;

        if( fMaxRot < f0 || fMaxRot <= -9999 )
            fMaxRot = f0;
        
        f0 = oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsRotMid) - fRotStep / 2.0;
        
        if( fMinRot > f0 || fMinRot <= -9999 )
            fMinRot = f0;
        //////////////////////////////////////////////////////////////////////////////////

        f0 = oReflnlist[nRef].fGetSigmaI();
        fVariance += f0*f0;
        
        f0 = oReflnlist[nRef].fGetIntensity();
        fIntensity     += f0;
        fVarianceNumer += f0*f0;
        
        if( fMax < f0 )
        {
            // Assume largest peak is on most intense spot.
            fMax = f0;
            
            fSharpness = oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsSharpness);
        
            fPeakAreaCount   = oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsPeakAreaCount);
            fPeakBorderCount = oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsPeakBorderCount);
        }
        
        a3fCent[0] += oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsPx0)*f0;
        a3fCent[1] += oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsPx1)*f0;
        a3fCent[2] += oReflnlist[nRef].fGetField(oReflnlist.m_nFI_fObsRotMid)*f0;
        
        nContrib++;
        
        nMaxContrib = max(nMaxContrib, nContrib);   


        if( 0 == oReflnlist[nRef].nPutGetEllipse(false,
                                                 &a7fEllipsoidNext[0],
                                                 &a7fEllipsoidNext[4],
                                                 &a7fEllipsoidNext[6],
                                                 &a2fCenterNext[0]) )
        {
            if( a7fEllipsoid[0] == 0.0 )
            {
                vCopyVecND(7, &a7fEllipsoidNext[0], &a7fEllipsoid[0]);
                vCopyVecND(2, &a2fCenterNext[0]   , &a2fCenter[0]);
            } 
            else
            {
                // Merge with the last ellipsoid.
                nMergeEllipsoids(&a2fCenter[0], 
                                &a7fEllipsoid[0],&a7fEllipsoid[4],a7fEllipsoid[6],
                                &a2fCenterNext[0],
                                &a7fEllipsoidNext[0],&a7fEllipsoidNext[4],a7fEllipsoidNext[6],
                                &a7fEllipsoidTemp[0],&a7fEllipsoidTemp[4],&a7fEllipsoidTemp[6]);
                
                vCopyVecND(7, &a7fEllipsoidTemp[0], &a7fEllipsoid[0]);
            }
        }
    }

    return nMaxContrib;
}

int Cfind::nMerge2DReflectionsWithTwins(Creflnlist* poNewList, bool bRejectOverlappedDifferentHKL) 
{
    int         nTwinID;
    int         nRefPred;
    bool        bAnotherTwinLaw;
    bool        bAnotherTwin;
    int         nTwinLaw;
    int         nInserted;
    int         nStat;
    int         nRef;
    int         nMaxOrigRefNum = 0;
    Creflnlist* poOriginalSpotList;
    itr<int>    anMapFromOrig2Calc;
    itr<int>    anTo;
    int         nx;
    -anTo;
    -anMapFromOrig2Calc;

    // Go through the list, and figure out where each reflection in the original file is mapped to.
    if ((m_poReflnlist->m_nFI_nSourceRefNum>=0) && (m_poReflnlist->m_nFI_nCompeteRefNum)) {

        nMaxOrigRefNum = 0;
        for (nRef = 0; nRef < m_poReflnlist->nGetNumReflns(); nRef++) {
            nMaxOrigRefNum = max(nMaxOrigRefNum,(*m_poReflnlist)[nRef].nGetField(m_poReflnlist->m_nFI_nSourceRefNum));
        };
        nCreateMapping(anMapFromOrig2Calc,nMaxOrigRefNum + 1);
        // Add to mapping.
        for (nRef = 0; nRef < m_poReflnlist->nGetNumReflns(); nRef++) {
            nAddToMapping(anMapFromOrig2Calc,(*m_poReflnlist)[nRef].nGetField(m_poReflnlist->m_nFI_nSourceRefNum),nRef);
        };
        for (nRef = 0; nRef < m_poReflnlist->nGetNumReflns();nRef++) {
            (*m_poReflnlist)[nRef].vSetField(m_poReflnlist->m_nFI_nSourceRefNum,-1);
        };

        for (nRef = 0; nRef < nMaxOrigRefNum; nRef++)
          {
            nGetMapping(anMapFromOrig2Calc,nRef,anTo,false);

            if (anTo.size() >= 2)
              {
                int a3nHKL[3],a3nPrevHKL[3];
                
                int         nTwinID = 0;
                int         nPrevTwinID = 0;
                
                // We might have this spot in competition.
                // Check to see if there are two unique HKL's or two unique twin IDs which 'claim' the spot.
                // If so, then the spot is competed for by the various solutions.
                for (nx = 0; nx < anTo.size(); nx++)
                  {
		    if (0 > anTo[nx])
		      {
			cout << "WARNING! anTo[nx] < 0: " << anTo[nx] << '\n';
			break;
		      }
                    a3nHKL[0] = (*m_poReflnlist)[anTo[nx]].nGetH();
                    a3nHKL[1] = (*m_poReflnlist)[anTo[nx]].nGetK();
                    a3nHKL[2] = (*m_poReflnlist)[anTo[nx]].nGetL();
                    nTwinID = (m_poReflnlist->m_nFI_nTwinID>=0)?((*m_poReflnlist)[anTo[nx]].nGetField(m_poReflnlist->m_nFI_nTwinID)):1;
                    if ((nx>0) && ((a3nHKL[0]!=a3nPrevHKL[0]) || (a3nHKL[1]!=a3nPrevHKL[1]) || (a3nHKL[2]!=a3nPrevHKL[2]) || (nTwinID!=nPrevTwinID)))
                      {
                        break;
                      };
                    nPrevTwinID = nTwinID;
                    a3nPrevHKL[0] = a3nHKL[0];
                    a3nPrevHKL[1] = a3nHKL[1];
                    a3nPrevHKL[2] = a3nHKL[2];
                  };
                if (nx < anTo.size())
                  {
                    for (nx = 0; nx < anTo.size(); nx++)
		      if (0 <= anTo[nx])
			(*m_poReflnlist)[anTo[nx]].vSetField(m_poReflnlist->m_nFI_nSourceRefNum,nRef);
                  };
              };
          };
        // Create a new field for the circularly linked 'compete' spot list.
        poNewList->nExpandGetField(poNewList->ms_snCompeteRefNum);
        poNewList->nExpandGetField(poNewList->ms_snSourceRefNum);
        
    } else {
      anMapFromOrig2Calc.clear();
    };

    // Make some mappings so that we are sure of 
    poOriginalSpotList = m_poReflnlist;

    bAnotherTwin = true;
    for (nTwinID = 1; bAnotherTwin; nTwinID++) 
    {
        
        nTwinLaw = 0;
        bAnotherTwin = false;
        do {
            Creflnlist oComponentList(*poOriginalSpotList);
            m_poReflnlist = & oComponentList;
            
            bAnotherTwinLaw = false;
            
            nInserted = 0;
            oComponentList.vDeleteAll();
            for (nRefPred = 0; nRefPred < poOriginalSpotList->nGetNumReflns(); nRefPred++) {
                if ((poOriginalSpotList->m_nFI_nTwinID>=0) && ((*poOriginalSpotList)[nRefPred].nGetField(poOriginalSpotList->m_nFI_nTwinID)==nTwinID + (nTwinLaw+1)*1000))
                    bAnotherTwinLaw = true;
                if ((poOriginalSpotList->m_nFI_nTwinID>=0) && (((*poOriginalSpotList)[nRefPred].nGetField(poOriginalSpotList->m_nFI_nTwinID) % 1000) == nTwinID + 1))
                    bAnotherTwin = true;
                if ((poOriginalSpotList->m_nFI_nTwinID<0) || ((*poOriginalSpotList)[nRefPred].nGetField(poOriginalSpotList->m_nFI_nTwinID)==nTwinID + nTwinLaw*1000)) {
                    oComponentList.nInsert(&(*poOriginalSpotList)[nRefPred]);
                    nInserted++;
                };
            };
            if (nInserted)
                        {
                nStat = nMerge2DReflections(poNewList, bRejectOverlappedDifferentHKL);            
                    if (nStat) return 1;
                        }
            
            nTwinLaw ++;
        } while (bAnotherTwinLaw);
        m_poReflnlist = poOriginalSpotList;
    };
    
    




    
    // Update the competing reflection data.
    if ((nMaxOrigRefNum>0) && (poNewList->m_nFI_nCompeteRefNum >= 0) && (poNewList->m_nFI_nSourceRefNum >= 0)) {

        nCreateMapping(anMapFromOrig2Calc,nMaxOrigRefNum + 1);
        // Add to mapping.
        for (nRef = 0; nRef < poNewList->nGetNumReflns(); nRef++) {
            nx = (*poNewList)[nRef].nGetField(poNewList->m_nFI_nSourceRefNum);
            if (nx >= 0) 
                nAddToMapping(anMapFromOrig2Calc,nx,nRef);
        };

        // Set all compete information to be null.
        for (nRef = 0; nRef < poNewList->nGetNumReflns();nRef++) {
            (*poNewList)[nRef].vSetField(poNewList->m_nFI_nCompeteRefNum,-1);
        };

        for (nRef = 0; nRef < nMaxOrigRefNum; nRef++) {
            nGetMapping(anMapFromOrig2Calc,nRef,anTo,false);
            if (anTo.size() >= 2) {
                // We have this spot in competition.
                // Create a circularly linked list containing all competitors.
                for (nx = 0; nx < anTo.size(); nx++) 
                    (*poNewList)[anTo[nx]].vSetField(poNewList->m_nFI_nCompeteRefNum,anTo[(nx+1) % anTo.size()]);
            };
        };
    };

    return 0;

};

int Cfind::nMerge2DReflections(Creflnlist* poNewList, bool bRejectOverlappedDifferentHKL)
{
    int     nRef  = 0;
    int     nRef2 = 0;

    int     nStat = 0;
    
    int     nContrib = 0;
    
    int     nFI_nHead = -1;
    int     nFI_nNext = -1;
    int     nFI_nOnEdge = -1;
    
    double  fRotStep = 0;

    const int nMaxBinCounts = 10;
    int     nPrevRefs = 0;
    int     nMergedRefs = 0;
    
    int     anCounts[nMaxBinCounts];
    itr<int> anRefs;

    int     nx = 0;

    itr<double> afScanAdjust;
    itr<double> afLastGonio;
    itr<double> afThisGonio;
    
    double      fScanAdjust = 0;
    
    int     nField = -1;
    
    int*        apnFields[] = 
    { 
        & m_poReflnlist->m_nFI_fObsRotMid,
        & m_poReflnlist->m_nFI_fObsRotEnd,
        & m_poReflnlist->m_nFI_fCalcRotStart,
        & m_poReflnlist->m_nFI_fCalcRotEnd,
        & m_poReflnlist->m_nFI_fCalcRotMid,
        NULL
    };

    int*    pnDelFlag = new int[m_poReflnlist->nGetNumReflns()];

    int     nErrorRefs = 0;
    
    int     nTotalReflns = m_poReflnlist->nGetNumReflns();
    printf("\nTotal number of reflns in m_poReflnlist: %d\n", nTotalReflns);

    for(nRef=0; nRef < m_poReflnlist->nGetNumReflns(); nRef++) 
    {
        Crefln*     poRefln = m_poReflnlist->poGetRefln(nRef);
        
        if( poRefln->fGetIntensity() <= 0.0 )
        {
            pnDelFlag[nRef] = -1;
            nErrorRefs++;
        } 
        else
            pnDelFlag[nRef] = 0;
    }
    //printf("\nTotal number of reflns with error in m_poReflnlist: %d\n", nErrorRefs);
    
    m_poReflnlist->nDelete(-1, pnDelFlag);
    
    delete[] pnDelFlag;
    
    // Determine the number of distinct goniometer scans we have.
    // We can only merge goniometer scans that have the same constants in them.
    // Since we already account for the rotation of each reflection in this merge,
    // we can fix this problem by making sure that different scans have vastly different
    // rotation centroids.  This is done by adjusting the rotation centroids accordingly
    // before we do any merging, and the unadjusting the centroids thereafter.  This will
    // garauntee that two reflections from two different scans are not accidentally merged.
    // (In reality, this would probably not happen anyway because two different scans have 
    // different rotation ranges.  Yet there is always the possibility).

    fScanAdjust = 0.0;
    -afScanAdjust;
    -afLastGonio;
    
    for(nRef=0; nRef < m_poReflnlist->nGetNumReflns(); nRef++) 
    {
        -afThisGonio;
        
        if (m_poReflnlist->m_nFI_fGonio1>=0)
            afThisGonio + (double) ((*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio1));
        
        if (m_poReflnlist->m_nFI_fGonio2>=0)
            afThisGonio + (double) ((*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio2));
        
        if (m_poReflnlist->m_nFI_fGonio3>=0)
            afThisGonio + (double) ((*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio3));
        
        if (m_poReflnlist->m_nFI_fGonio4>=0)
            afThisGonio + (double) ((*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio4));
        
        if (m_poReflnlist->m_nFI_fGonio5>=0)
            afThisGonio + (double) ((*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio5));

        if (afLastGonio.size() == afThisGonio.size())
        {
            for (nx=0;nx < afLastGonio.size(); nx++)
            {
                if (ABS(afThisGonio[nx] - afLastGonio[nx]) > 0.01)
                    break;
            }
            
            if (nx != afLastGonio.size())
            {
                // We have a different goniometer setting.
                // Increment the goniometer offset.
                fScanAdjust += 1000.0;
            } 
        }             
        
        // Apply the goniometer offset to this reflection.
        for (nField = 0; apnFields[nField];nField++)
        {
            if (*apnFields[nField]>=0)
                (*m_poReflnlist)[nRef].vSetField(*apnFields[nField],
                                                (float) ( fScanAdjust + (*m_poReflnlist)[nRef].fGetField(*apnFields[nField]) ) 
                                                );
        }
        
        afScanAdjust + fScanAdjust;
        
        afLastGonio.copy(afThisGonio);
    }
    

    // Initialize counting information.
    nMergedRefs = 0;
    for (nx=0;nx<nMaxBinCounts; nx++) 
        anCounts[nx] = 0;


    // This is used for width sigma calculations.
    fRotStep = m_poScan->m_poRotation->fGetIncrement();
    if (fRotStep==0.0)
    {
        fRotStep = m_fDefaultRockingWidth;
        
        printf("WARNING:  Rocking width of image is 0.0 !!?? \n"
               "Setting arbitrarily to width of %.1lf degrees\n",fRotStep);        
    }

    nFI_nHead = m_poReflnlist->nExpandGetField(ms_sMergeHead);
    nFI_nNext = m_poReflnlist->nExpandGetField(ms_sMergeNext);
    nFI_nOnEdge = m_poReflnlist->nExpandGetField(ms_sOnFirstLastImage);
        
    nStat = m_poReflnlist->nOverlapCheck(nFI_nHead, nFI_nNext, fRotStep);       

    if( nStat < 0)
        return 1;

    nExpandReflnlist(poNewList, false, true);
    
    poNewList->nInsertListFrom(*m_poReflnlist,NULL,NULL,0);   // Initializes nInsertReflnFrom().
    
    Crefln      oNewRefln(poNewList);

    float       fTemp = 0.0f;
    // Next we need to merge reflections...
    for(nRef=0; nRef < m_poReflnlist->nGetNumReflns(); nRef++)
    {
        oNewRefln = (*m_poReflnlist)[nRef];
        
        ///////////////////////////////////////
        // RB This code is for debugging only
        //int     nH =  oNewRefln.nGetH();
        //int     nK =  oNewRefln.nGetK();
        //int     nL =  oNewRefln.nGetL();
        //
        //if( nH ==  -9 && nK == 9 && nL == 32 )
        //{
        //    
        //    float   fInt =  oNewRefln.fGetIntensity();
        //    
        //    int     nTrap1 = 0;
        //}
        ///////////////////////////////////////


        // If this reflection is the head reflection then process it.
        if (nRef == (*m_poReflnlist)[nRef].nGetField(nFI_nHead))
        {
            // Build list of reflections.
            -anRefs;
            
            for(nRef2=nRef; (nRef2!=-1); nRef2 = (*m_poReflnlist)[nRef2].nGetField(nFI_nNext)) 
                anRefs + nRef2;

            if( bRejectOverlappedDifferentHKL && anRefs.size() > 1 )
            {
                if( bListOfReflnIndexesHasDifferentHKLs(anRefs, *m_poReflnlist) )
                    continue;  // we reject such reflections
            }
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Adjust the reflection centroid fields to compensate for the scan.
            fScanAdjust = afScanAdjust[anRefs[0]];
            if( 0.0 != fScanAdjust )
            {
                for(nRef2 = 0; nRef2 < anRefs.size(); nRef2++)
                {
                    for (nField = 0; apnFields[nField]; nField++)
                    {
                        fTemp = (*m_poReflnlist)[anRefs[nRef2]].fGetField(*apnFields[nField]);

                        if( *apnFields[nField] >=0 )
                            (*m_poReflnlist)[anRefs[nRef2]].vSetField(*apnFields[nField], (float)(-fScanAdjust + fTemp));
                    }               
                }
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            nPrevRefs = poNewList->nGetNumReflns();

            if( 1 == anRefs.size() )
            {
                poNewList->nInsert(&oNewRefln);  // insert just one reflection itself              
                nContrib = 1;
            }
            else
            {
                // insert one or more reflections after merging them
                nContrib = nMerge2D_OverlapReject(*m_poReflnlist, anRefs, oNewRefln, *poNewList);
            }
            
            nMergedRefs += (poNewList->nGetNumReflns() - nPrevRefs);
            
            anCounts[ min(nMaxBinCounts-1, nContrib-1) ]++;
        }
    }

    printf("\n\n");
    printf("Reflection status after merging\n");
    printf("-------------------------------\n");
    printf("Raw 2D centroids                     %5d\n",m_poReflnlist->nGetNumReflns());
    printf("Merged                               %5d\n",nMergedRefs);
    printf("\n");
    
    for (nx=0;nx<nMaxBinCounts; nx++)
    {
        if (anCounts[nx])
            printf("Reflns. on %3d consecutive image(s)  %5d\n",nx+1,anCounts[nx]);
    }
    
    printf("\n");
    
    return 0;
}

int
Cfind::nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre)
{
    int nx;
  // Place the find results in the header

  poHeader->nReplaceValue(ms_sDtfindOptions, m_sOptionString);
  if ((m_poReflnlist) && (m_nTotalImages))
      poHeader->nReplaceValue(ms_sDtfindAvgSpotsPerImage,m_poReflnlist->nGetNumReflns()/m_nTotalImages);
  else
      poHeader->nReplaceValue(ms_sDtfindAvgSpotsPerImage,0);

  poHeader->nReplaceValue(ms_sDtfindFractionSaturated,m_nTotalSatSpotsProcessed/((double) max(1,m_nTotalSpotsProcessed)));

  m_oDetectorAreas.nUpdateHeader(poHeader);

  double afPhotonIoverSig[10];
  double afPhotonValues[10];
  int    nRef;
  double fIoverSig;
  for (nx=0;nx<10;nx++) {
      afPhotonIoverSig[nx] = exp(nx*0.4);
      afPhotonValues[nx] = 0.0;
  };
  for (nRef = 0; nRef < m_poReflnlist->nGetNumReflns(); nRef++) {
      fIoverSig = (*m_poReflnlist)[nRef].fGetIntensity()/max(1.0,(*m_poReflnlist)[nRef].fGetSigmaI());
      for (nx=0;(nx < 10 - 1) && (fIoverSig > afPhotonIoverSig[nx]); nx++)
        afPhotonValues[nx] += (*m_poReflnlist)[nRef].fGetIntensity();
  };
  poHeader->nReplaceValue(ms_sDtfindPhotonValues,10,&afPhotonValues[0]);
  poHeader->nReplaceValue(ms_sDtfindPhotonIoverSig,10,&afPhotonIoverSig[0]);

  return (0);
}

int  Cfind::nWrite(Creflnlist& oList,Cstring& sName) {
    char* pcSelectFields;
    int nx;
    
//+jwp 13-Mar-2003 
// Put some info here
    if (0 < m_nTotalSpotsProcessed)
      {
        printf("There were %d" 
               " preliminary spots of which %d"
               " were marked as saturated,\n or"
               "%6.2lf%% of them.\n", 
               m_nTotalSpotsProcessed, m_nTotalSatSpotsProcessed,
               100.0 * (double)m_nTotalSatSpotsProcessed / (double)m_nTotalSpotsProcessed);
        fflush(stdout);
      }
//
//-jwp 13-Mar-2003    // Put some info here

    pcSelectFields = new char [oList.m_nTotalFieldsPlus1];
    
    for (nx = 0; nx < oList.m_nTotalFieldsPlus1; nx++)
    {
        pcSelectFields[nx] = '\1';  
    }
    pcSelectFields[oList.m_nTotalFieldsPlus1-1] = '\0';
    
    // Only write out the following:
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fIntensity, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fSigmaI, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fObsPx0, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fObsPx1, char(0),
        pcSelectFields);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    (void) oList.nSelectField(eReflnField_float_type, oList.m_nFI_fObsSharpness,        char(0),    pcSelectFields);
    
    (void) oList.nSelectField(eReflnField_float_type, oList.m_nFI_fObsPeakAreaCount,    char(0),    pcSelectFields);

    (void) oList.nSelectField(eReflnField_float_type,  oList.m_nFI_fObsPeakBorderCount, char(0),    pcSelectFields);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fObsRotMid, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fObsRotWidth, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fObsRotSigma, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fDetResolution, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fEllipsoidA00, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fEllipsoidA01, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fEllipsoidA11, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fEllipsoidb0, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fEllipsoidb1, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fEllipsoidc, char(0),
        pcSelectFields);
    
    (void) oList.nSelectField(eReflnField_int_type,
        oList.m_nFI_nH, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_int_type,
        oList.m_nFI_nK, char(0),
        pcSelectFields);
    (void) oList.nSelectField(eReflnField_int_type,
        oList.m_nFI_nL, char(0),
        pcSelectFields);

        if (oList.m_nFI_fGonio1>=0)
                (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fGonio1, char(0),
        pcSelectFields);
        if (oList.m_nFI_fGonio2>=0)
                (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fGonio2, char(0),
        pcSelectFields);
        if (oList.m_nFI_fGonio3>=0)
                (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fGonio3, char(0),
        pcSelectFields);
        if (oList.m_nFI_fGonio4>=0)
                (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fGonio4, char(0),
        pcSelectFields);
        if (oList.m_nFI_fGonio5>=0)
                (void) oList.nSelectField(eReflnField_float_type,
        oList.m_nFI_fGonio5, char(0),
        pcSelectFields);
        if (oList.m_nFI_nGonioRotAxis>=0)
                (void) oList.nSelectField(eReflnField_int_type,
        oList.m_nFI_nGonioRotAxis, char(0),
        pcSelectFields);
    
    
    oList.nWrite(sName,NULL,pcSelectFields);
    
    delete[] pcSelectFields;
    return 0;
};


void Cfind::vSetResolution(const float fResoMin, const float fResoMax)
{
  // Set the min and max resolution in Angstroms.  Min is closest to
  // beam, max is as far away from beam, thus Min has a higher value
  // in Angstroms than Max

  m_fResolutionMin = max(fResoMin, fResoMax);
  m_fResolutionMax = min(fResoMin, fResoMax);
}


    // Estimate the xy sigmas.      
int nEstimateSigmaXY(int nPx0,int nPx1,float* fSig0,float* fSig1);
    // Estimate the rotation sigmas.
int nEstimateSigmaRot(int nPx0,int nPx1,float* fSig);

////////////////////////////////////////////////////////////////////////////////////////////////////
// RB: adapted from dtfind.cpp
// It is assumed that the Cfind object has the header information set already
int Cfind::nDoSearch(float  fResolutionAngstr_LowestTheta,
                     float  fResolutionAngstr_HighestTheta,
                     eFind_methods eMethod,
                     float fIntensityOverSigmaRatio,
                     int   nFirstImage,
                     int   nLastImage,
                     int   nImageSequenceStep)
{
    m_eMethod = eMethod;
    m_fSigma = fIntensityOverSigmaRatio;
    
    // It is assumed that m_a2nSeqNum is already set
    if( nFirstImage < m_a2nSeqNum[0] )
        return 1;  // inconsistent starting image 
    if( nLastImage > m_a2nSeqNum[1] )
        nLastImage = m_a2nSeqNum[1]; // cannot go beyond available images

    //cout << "Cfind::nDoSearch, SeqFirst, Last: " << m_a2nSeqNum[0] << ", " << m_a2nSeqNum[1] << " \n";
    //cout << "Cfind::nDoSearch, First, Last: " << nFirstImage << ", " << nLastImage << " \n";

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Determine the required resolution
    
    // Check the input.
    
    // Make sure that the input resolution value in angstroms is higher for the lower thetas than for the higher thetas
    if( fResolutionAngstr_LowestTheta < fResolutionAngstr_HighestTheta )
      std::swap(fResolutionAngstr_LowestTheta, fResolutionAngstr_HighestTheta);
    
    // Get the detector resolution

    // Create the source object just once
    Csource*        poSource   = new Csource(m_oHeader);
    if (!poSource->bIsAvailable())
    {
        cout << "ERROR: dtfind - creating source object!\n";
        return 1;
    }
    float       fTempW = poSource->fGetWavelength();
    float       a3fS0[3];
    poSource->vCalcGetS0(&a3fS0[0]);  // Required for resolution calc below

    // Create the detector object just once and do not create
    // its non-uniformity, since we have that in the poScan object.
    Cdetector*      poDetector = new Cdetector(m_oHeader, "", TRUE, FALSE);
    if (!poDetector->bIsAvailable())
    {
        cout << "ERROR: dtfind - creating detector object!\n";
        return 1;
    }
    float       fTemp1=0.0;
    float       fTemp2=0.0;
    poDetector->nGetResolution(a3fS0, &fTemp1, &fTemp2);
    fTemp1 *= fTempW;
    fTemp2 *= fTempW;

    if( fTemp1 < fResolutionAngstr_LowestTheta )
        fResolutionAngstr_LowestTheta = fTemp1;
    
    if (fTemp2 > fResolutionAngstr_HighestTheta )
        fResolutionAngstr_HighestTheta = fTemp2;
    
    vSetResolution(fResolutionAngstr_LowestTheta, fResolutionAngstr_HighestTheta);
    cout << "\nResolution limits of an image are " << fTemp1 << " to " << fTemp2 << '\n';
    cout << "Resolution limits of peak search are " << fResolutionAngstr_LowestTheta << " to " << fResolutionAngstr_HighestTheta << '\n' << endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    Cimage*         poImage = new Cimage();
    Crotation*      poRotation = NULL;
    Crefln*         poRefln = NULL;
    int             nStat = 0;
    
    int             nFirstReflnThisImage = 0;
    bool            bImageAvailable = true;
    
    Cstring         sImagePath("");

    for(int iImage = nFirstImage; iImage <= nLastImage; iImage+=nImageSequenceStep)
    {
      //cout << "Cfind::nDoSearch, iImage: " << iImage << " \n";
      //cout << "Cfind::nDoSearch, First, Last: " << nFirstImage << ", " << nLastImage << " \n";
      poImage->m_bReadApplyGonioMask = true;
      poImage->m_bReadApplyEmbeddedMask = true;
        
      ///////////////////////////////////////////////////////////////
      // Get current image
      m_poScan->vSetSeqNum(iImage);
      
      sImagePath = "";
      m_poScan->nGetImageName(&sImagePath);
        
      //////////////////////////////////////////////////////////////////
      if( THE_IMAGE_WAIT_PTR->bWaitRequired() )
        {    
          if( iImage == nFirstImage )
            bImageAvailable = THE_IMAGE_WAIT_PTR->bDoWaitForFirstImage(sImagePath);
          else
            bImageAvailable = THE_IMAGE_WAIT_PTR->bDoWaitForNextImage(sImagePath);
        }

        if( !bImageAvailable )
            return 1;
        ////////////////////////////////////////////////////////////////////

        printf("...reading image %s\n", sImagePath.string());
        fflush(stdout);

        m_poScan->nGetImage(poImage);
    
        ///////////////////////////////////////////////////////////////
        // Set up rotation object for this image
        if( poRotation )
        {
            delete poRotation;
            poRotation = NULL;
        }
        
        poRotation = new Crotation(poImage->m_oHeader);
        if( !poRotation->bIsAvailable() )
        {
            // Rotation info not in image header, so try to
            // create rotation from the scan info we have
            delete poRotation;
            poRotation = NULL;

            poRotation = new Crotation(m_poScan, m_poScan->nGetSeqNum());
        }
        
        nMarkGonioNumbers();
        
        nStat = nFind2D(*poImage, 
                        poRotation->fGetMidValue(),
                        poDetector, 
                        poSource,
                        eFind_2D_merge_method == m_eMethod);
        
        // RB: not sure why TJN used this...
        // We might need to add goniometer information.
        if( nAddGonioNumbers(*poImage) ) 
        {
            nStat = 1;
            break;
         }
        
        // OK, found spots, now find centroid of each one
        DTREK_WORD        wCtrl = m_bRank ? DTREK_DTFIND_FC2D_GET_PEAK_AREA_INFO : 0U; 

        for(int iRefln = nFirstReflnThisImage; iRefln < m_poReflnlist->nGetNumReflns(); iRefln++)
        {
            poRefln = m_poReflnlist->poGetRefln(iRefln);
            
            nStat = nFindCentroid2D(*poImage, poRefln, poDetector, wCtrl);
            
            if (0 != nStat)
            {
                // Reset local status and mark refln for deletion
                nStat = 0;
                poRefln->vSetIntensity(-999.0);
            }
        }
    
        // Update nFirstReflnThisImage
        nFirstReflnThisImage = m_poReflnlist->nGetNumReflns();
        
        if( THE_IMAGE_WAIT_PTR->bWaitRequired() && !THE_IMAGE_WAIT_PTR->bSaveImageSize(sImagePath) )
            return 1;
//+2011-06-16
        // If the image step is largish and will go over the end point by 1 or 2 images, reset so last image is used
        if ( (3 < nImageSequenceStep) && (nLastImage < (iImage + nImageSequenceStep)))
          {
            if (  2 > ((iImage + nImageSequenceStep) - nLastImage))
              iImage = nLastImage - nImageSequenceStep;
          }
//-2011-06-16
    }
    
    // Delete bad reflections
    for(int jj = m_poReflnlist->nGetNumReflns()-1; jj >= 0; jj--)
    {
        poRefln = m_poReflnlist->poGetRefln(jj);
        
        if( 0.0 > poRefln->fGetIntensity() )
        {
            m_poReflnlist->nDelete(jj);
        }
    }         

    return nStat;
}

int Cfind::nWriteReflectionList(Cstring& strReflectionFileName)
{
    return m_poReflnlist->nWrite(strReflectionFileName);
}

void Cfind::vHeaderToDetectorAreas()
{
    m_oDetectorAreasIn.nInitValues(m_oHeader);
}


/// RB 2/16/04 For possible later use
//void Cfind::vGetPeakAreaSize(int nLocation0_PIX, int nLocation1_PIX, int& nSize0_PIX, int& nSize1_PIX)
//{
//  if ( m_a2nSpotWindow[0] >= 0 || m_a2nSpotWindow[1] >= 0 )
//  {
//      nSize0_PIX = m_a2nSpotWindow[0];
//      nSize1_PIX = m_a2nSpotWindow[1];
//  }
//  else
//  {
//      C3DdataDetectorArea       oArea;
//      m_oDetectorAreasIn.nGetPeakArea(nLocation0_PIX, nLocation1_PIX, oArea);      
//      
//      nSize0_PIX = (int)min(nMAXSPOTSIZE, 3.5 * oArea.fSize[0]);
//      nSize1_PIX = (int)min(nMAXSPOTSIZE, 3.5 * oArea.fSize[1]);
//  }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Go through a list of reflection list indexes and see if they correspond to reflections with different HKLs
bool bListOfReflnIndexesHasDifferentHKLs(itr<int>& anReflnIndexes, Creflnlist& oReflnlist)
{
    if( anReflnIndexes.size() < 2 )
        return false; // there should be at least two reflection list indexes to refer to different HKLs.
    
    int     nFirstHKL = oReflnlist[anReflnIndexes[0]].nPackHKL(); // we trust that anReflnIndexes[0] is a valid index for the reflection list
    int     nAnotherHKL = 0;

    for(int ii = 1; ii < anReflnIndexes.size(); ii++)
    {
        nAnotherHKL = oReflnlist[anReflnIndexes[ii]].nPackHKL(); // we trust that anReflnIndexes[ii] is a valid index for the reflection list
        
        if( nAnotherHKL != nFirstHKL )
            return true;
    }

    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
