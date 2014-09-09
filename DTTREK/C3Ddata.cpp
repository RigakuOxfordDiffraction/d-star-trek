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
// C3Ddata.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class C3Ddata which implements
//    the 3Ddata encapsulation of d*TREK.
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

#include <string.h>          // For memcpy
#include "Dtrek.h"
#include "C3Ddata.h"         // Class definition and prototypes
#include "dtsvd.h"

#include "dtrekdefs.h"
#include "minmax.h"
#include "Cprofit2D.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

//+Definitions, constants, and initialization of static member variables


#define C3DDATA_MAX_FREE 2048

int     C3Ddata::ms_nNumObjects       = 0;
int*    C3Ddata::ms_pnFreeSize        = NULL;
float** C3Ddata::ms_ppfFreeFloat      = NULL;
int     C3Ddata::ms_nFreeCount        = 0;
int     C3Ddata::ms_nSoftCount        = 0;
int     C3Ddata::ms_nHardCount        = 0;

int     C3Ddata::ms_nAllocSmooth      = 0;
int*    C3Ddata::ms_pnMask            = NULL;
float*  C3Ddata::ms_pfMask            = NULL;
float*  C3Ddata::ms_pfPrint           = NULL;

int     C3Ddata::ms_nAllocScratchInt = 0;
int*    C3Ddata::ms_apnScratchInt[2] = { NULL,NULL};

int     C3Ddata::ms_nHighMask =                0x0001;
int     C3Ddata::ms_nEllipsoidMask =           0x0002;
int     C3Ddata::ms_nBorderMask =              0x0004;
int     C3Ddata::ms_nLowMask =                 0x0008;
int     C3Ddata::ms_nSatMask =                 0x0010;
int     C3Ddata::ms_nBadMask =                 0x0020;
int     C3Ddata::ms_nReservedMask =            0x0040;
int     C3Ddata::ms_nCentroidPixelMask =       0x0080;
int     C3Ddata::ms_nContiguousHighMask =      0x0100;
int     C3Ddata::ms_nMinPeakPixels =           0x0200;

int     C3Ddata::ms_nBoundaryPeakPixels =      0x0400;  // a thin (1 pixel thick) pixel layer around the peak area

int		C3Ddata::ms_nSpecialOptions = -1;

int		C3Ddata::ms_nPerSliceIntegrateWeak = 4;
int		C3Ddata::ms_nPerSliceInMosaicity = 1;
int		C3Ddata::ms_nPerSliceCentroid = 2;
int     C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity = 8;
int     C3Ddata::ms_nPerSlicePostRefineNotInRmergeMosaicity = 16;

Cprofit2DPerSlice* C3Ddata::m_poProfit = NULL;

//+Code begin

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
///////////////// C3DdataInput/C3DdataOutput ///////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



void C3DdataInput::vCopy(const C3DdataInput& oOther) { memcpy(this,&oOther,sizeof(C3DdataInput)); };
void C3DdataOutput::vCopy(const C3DdataOutput& oOther) { memcpy(this,&oOther,sizeof(C3DdataOutput)); };
C3DdataInput::C3DdataInput() { vInit(); };
C3DdataOutput::C3DdataOutput() { vInit(); };

void C3DdataInput::vInit() {
  m_a2fShiftKa1[0]		 = 0.0;
  m_a2fShiftKa1[1]		 = 0.0;
  m_a2fShiftKa2[0]		 = 0.0;
  m_a2fShiftKa2[1]		 = 0.0;
  m_fDetectorGain        = 2.5;     // See Cintegrate.cc
  m_a2x2fEllipsoidAIn[0][0] = 0.0;
  m_a2x2fEllipsoidAIn[0][1] = 0.0;
  m_a2x2fEllipsoidAIn[1][0] = 0.0;
  m_a2x2fEllipsoidAIn[1][1] = 0.0;
  m_a2fEllipsoidbIn[0] = 0.0;
  m_a2fEllipsoidbIn[1] = 0.0;
  m_fEllipsoidcIn = 0.0;

  m_fUseBackgroundBorder = 0.5;
  m_fOnEdgeMaxPercent    = 0.00f;   // Set to zero tollerance.
  m_fMinPercentBackground= 0.02f;
  m_fSigmaAboveBackground= 3.0;    // See Cintegrate.cc
  m_fSigmaBelowBackground= 2.50;
  m_nBackgroundBorder    = 2;       // See Cintegrate.cc    
  m_fFattening           = 1.2f;
  m_fFractionIntense     = 0.98f;
  m_fMosaicityLowestSpotIoverSig = 6.0;
  m_nMinPeakRad = 2;
  m_nIntegrateWeakFlag     = 1;
  m_bPrint                 = FALSE;
  m_bFitToPlanar           = TRUE;
  m_bSpotChase             = TRUE;
  m_bUsePerSliceBackground = TRUE;
  m_bProfileFitWeak		   = FALSE;
  m_bIncludeSat            = FALSE;
  m_fMinProfitPeakIntensity = 0.0;
  m_a3nHKL[0]              = -999;
  m_a3nHKL[1]              = -999;
  m_a3nHKL[2]              = -999;
};


void C3DdataOutput::vInit() {

    m_fAvgBackgroundAvg = 0.0;
    m_fAvgBackgroundSigma = 0.0;
    m_a3fAvgBackground[0] = 0.0;
    m_a3fAvgBackground[1] = 0.0;
    m_a3fAvgBackground[2] = 0.0;
    m_a3fCentroid[0]       = 0.0;
    m_a3fCentroid[1]       = 0.0;
    m_a3fCentroid[2]       = 0.0;
    m_fRotSigma            = -999.0;
    m_fPeakIntensity       = -999.0;
    m_fPeakIntensitySigma  = -999.0;
    m_fPeakBackgroundSigma = -999.0;
    m_a3fPeakSize[0]       = 1.0;
    m_a3fPeakSize[1]       = 1.0;
    m_a3fPeakSize[2]       = 1.0;
    m_fPeakMinValue        = -999.0;
    m_fPeakMaxValue        = -999.0;
    
    m_fPeakSharpness       = 0.0;

    m_a3fPeakMaxValue[0]   = 0.0;
    m_a3fPeakMaxValue[1]   = 0.0;
    m_a3fPeakMaxValue[2]   = 0.0;
    m_nPeakStart           = -1;
    m_nPeakEnd             = -1;
    m_nIntensePeakStart    = -1;
    m_nIntensePeakEnd      = -1;
    m_nCentroid            = -1;
    m_nSatCount            = 0;
    m_a2x2nEllipsoidRange[0][0] = 0;
    m_a2x2nEllipsoidRange[0][1] = 0;
    m_a2x2nEllipsoidRange[1][0] = 0;
    m_a2x2nEllipsoidRange[1][1] = 0;
    m_a2x2fEllipsoidA[0][0] = 0.0;
    m_a2x2fEllipsoidA[0][1] = 0.0;
    m_a2x2fEllipsoidA[1][0] = 0.0;
    m_a2x2fEllipsoidA[1][1] = 0.0;
    m_a2fEllipsoidb[0]      = 0.0;
    m_a2fEllipsoidb[1]      = 0.0;
    m_fEllipsoidc           = 0.0;

    m_a3fAbsCentroid[0]       = 0.0;
    m_a3fAbsCentroid[1]       = 0.0;
    m_a3fAbsCentroid[2]       = 0.0;
    m_fRotWidth = 0.0;
};


void
C3DdataOutput::vGetCentroid(float *pfCentroid)
{
  register int i;
  for (i = 0; i < 3; i++)
    {
      pfCentroid[i] = m_a3fCentroid[i];
    }
}


void
C3DdataOutput::vGetPeakSize(float *pfSize)
{
  register int i;
  for (i = 0; i < 3; i++)
    {
      pfSize[i] = m_a3fPeakSize[i];
    }
}

void
C3DdataOutput::vGetPeakStartEnd(int *pnStart, int *pnEnd) {
  *pnStart = m_nPeakStart;
  *pnEnd   = m_nPeakEnd;
}
void
C3DdataOutput::vGetIntensePeakStartEnd(int *pnStart,int *pnEnd) {
    *pnStart = m_nIntensePeakStart;
    *pnEnd = m_nIntensePeakEnd;
};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//////////////////////// C3Ddata ///////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//+Public functions

// Constructors, destructors and assignments

C3Ddata::C3Ddata()
{
  if (0 == ms_nNumObjects) 
      vFreeInit();
  ms_nNumObjects++;
 (void) nInitValues();
}

C3Ddata::C3Ddata(const int nInExt0, const int nInExt1, const int nInExt2,
		 const int nOrig0, const int nOrig1, const int nOrig2,
		 const e3Ddata_types eType)
{
  if (0 == ms_nNumObjects) 
      vFreeInit();

  ms_nNumObjects++;
  (void) nInitValues();            // m_nFilled is set to 0s in here
  m_nExt[0]        = max(nInExt0, 0); // Any error checking to see if less than 0?
  m_nExt[1]        = max(nInExt1, 0); // Any error checking to see if less than 0?
  m_nExt[2]        = max(nInExt2, 0); // Any error checking to see if less than 0?
  m_nOrigOffset[0] = nOrig0;
  m_nOrigOffset[1] = nOrig1;
  m_nOrigOffset[2] = nOrig2;
  m_pfData         = NULL;
  m_puiData        = NULL;
  m_pflData	   = NULL;
  m_pnData         = NULL;
  m_eType          = eType;
  vGetMemory(m_eType);
}

C3Ddata::C3Ddata(const C3Ddata& oOther)
{
  if (0 == ms_nNumObjects) 
      vFreeInit();
  ms_nNumObjects++;
  m_eType          = oOther.m_eType;
  m_nExt[0]        = oOther.m_nExt[0];
  m_nExt[1]        = oOther.m_nExt[1];
  m_nExt[2]        = oOther.m_nExt[2];
  m_nFilled[0]     = 0;
  m_nFilled[1]     = 0;
  m_nFilled[2]     = 0;
  m_nOrigOffset[0] = oOther.m_nOrigOffset[0];
  m_nOrigOffset[1] = oOther.m_nOrigOffset[1];
  m_nOrigOffset[2] = oOther.m_nOrigOffset[2];
  m_puiData = NULL;
  m_pfData  = NULL;
  m_pflData = NULL;
  m_pnData  = NULL;
  vGetMemory(m_eType);
}

C3Ddata::~C3Ddata()
{
  vFreeMemory(&m_pflData);
  m_puiData = NULL;
  m_pfData  = NULL;
  m_pflData = NULL;
  m_pnData  = NULL;
  m_nAlloc  = 0;
  ms_nNumObjects--;
  if (0 >= ms_nNumObjects)
    {
      ms_nNumObjects = 0;
      vFreeDelete();
      vFreeScratch();
    }
}


C3Ddata& C3Ddata::operator+=(C3Ddata& o3DdataIn) {
    int a6nIntersection[6];
    int a3nPos[3];

    if (!bIntersect(o3DdataIn,a6nIntersection)) 
        return (*this);
    if (m_eType != o3DdataIn.m_eType)
        return (*this);

    for (a3nPos[0] = a6nIntersection[0]; a3nPos[0]<=a6nIntersection[0 + 3]; a3nPos[0]++) {
        for (a3nPos[1] = a6nIntersection[1]; a3nPos[1]<=a6nIntersection[1 + 3]; a3nPos[1]++) {
            for (a3nPos[2] = a6nIntersection[2]; a3nPos[2]<=a6nIntersection[2 + 3]; a3nPos[2]++) {
                if (m_eType==e3Ddata_float) {
                    vSetValue(o3DdataIn.fGetValue(a3nPos[0] - o3DdataIn.m_nOrigOffset[0],a3nPos[1] - o3DdataIn.m_nOrigOffset[1],a3nPos[2] - o3DdataIn.m_nOrigOffset[2]),
                        a3nPos[0] - m_nOrigOffset[0], a3nPos[1] - m_nOrigOffset[1], a3nPos[2] - m_nOrigOffset[2]);
                } else if (m_eType==e3Ddata_double) {
                    vSetValue(o3DdataIn.lfGetValue(a3nPos[0] - o3DdataIn.m_nOrigOffset[0],a3nPos[1] - o3DdataIn.m_nOrigOffset[1],a3nPos[2] - o3DdataIn.m_nOrigOffset[2]),
                        a3nPos[0] - m_nOrigOffset[0], a3nPos[1] - m_nOrigOffset[1], a3nPos[2] - m_nOrigOffset[2]);
                    
                } else if (m_eType==e3Ddata_ushort) {
                    vSetValue(o3DdataIn.uiGetValue(a3nPos[0] - o3DdataIn.m_nOrigOffset[0],a3nPos[1] - o3DdataIn.m_nOrigOffset[1],a3nPos[2] - o3DdataIn.m_nOrigOffset[2]),
                        a3nPos[0] - m_nOrigOffset[0], a3nPos[1] - m_nOrigOffset[1], a3nPos[2] - m_nOrigOffset[2]);
                };
            };
        };
    };
    // Automatically set the filled status since it's so easy.
    vSetFilled(a6nIntersection[3]+1 - m_nOrigOffset[0],a6nIntersection[4]+1 - m_nOrigOffset[1],a6nIntersection[5]+1 - m_nOrigOffset[2]);
    return (*this);
};

int 
C3Ddata::nInitValues(void)
{
  for (int i = 0; i < 3; i++)
    {
      m_nExt[i]        = 0;
      m_nFilled[i]     = 0;
      m_nOrigOffset[i] = 0;
    }
  m_nAlloc         = 0;
  m_pfData         = NULL;
  m_pflData        = NULL;
  m_puiData        = NULL;
  m_pnData         = NULL;
  m_eType          = e3Ddata_ushort;
  m_poImage        = NULL;
        
  if (ms_nSpecialOptions == -1) 
    {
      Cstring sTemp = sGetEnv("DTREK_C3DDATA_SPECIAL");
      if (sTemp.length()!=0) 
        {
          if (1 != sscanf(sTemp.string(),"%d",&ms_nSpecialOptions)) 
            {
              printf("Could not read bitmask in DTREK_C3DDATA_SPECIAL");
              ms_nSpecialOptions = 0x000b; // This is 11
            }
          else
            {
              cout << "SPECIAL OPTIONS: " << ms_nSpecialOptions << endl;
            }
        }
      else
        {
          ms_nSpecialOptions = 
            (1 << eSpecial_TrimEdges) |
                                //(1 << eSpecial_DoNotIntegrateGarbage) |
            (1 << eSpecial_PhiHKLLimits);
                                //(1 << eSpecial_UseMosaicityBoundsOnStrong);
        }
    }
        
  return (0);
}

int C3Ddata::nFill2D(Cimage *poImage, const int nWhich)
{
  int nLayer;
  if (0 >= m_nAlloc)
    {
      cout << "ERROR: 3Ddata No space allocated!\n";
      return (1);
    }
  else if ( (nWhich >= m_nExt[2]) || (-1 > nWhich) )
    {
      cout << "ERROR: 3Ddata request to fill beyond allocated space!\n";
      return (2);
    }
  else if (-1 == nWhich)
    {
      nLayer = m_nFilled[2];   // Automatically fill the next layer
      m_nFilled[2]++;          // and keep track of it
    }
  else
    {
      nLayer = nWhich;
    }

  // Always try to point to last image used.  Be careful about using this!

  m_poImage = poImage;

  int i, j;
  // DEBUG
  if (m_puiData !=  (unsigned short *) m_pfData)
    {
      cout << "ERROR in C3Ddata::nFill2D address problems!\n";
      m_puiData = (unsigned short *) m_pfData;
    }
  if (m_pnData !=  (int *) m_pfData)
    {
      cout << "ERROR in C3Ddata::nFill2D address problems!\n";
      m_pnData = (int *) m_pfData;
    }

  int a2nOffset[2];
  static bool s_bFirstTime = TRUE;
  static int  s_a2nSpecial[2] = {0, 0};

  if (s_bFirstTime)
    {
      // Get special offsets from image header for IP2

      s_bFirstTime = FALSE;

      // If the keyword doesn't exist, a2nOffset will remain unchanged.

      if (0 == poImage->m_oHeader.nGetValue("RAXIS_IP2_OFFSETS", 2, &s_a2nSpecial[0]))
	{
	  cout << "INFO: Pixel offsets for IP2 set to " << s_a2nSpecial[0] 
	       << ", " << s_a2nSpecial[1] << endl << flush;
	}
    }

  a2nOffset[0] = m_nOrigOffset[0];
  a2nOffset[1] = m_nOrigOffset[1];

  if (2 == poImage->nGetDetNum())
    {
      // Apply a special additional offset for IP2 images.
      // WARNING: array limits were checked when this C3Ddata object
      //          was instantiated.  A problem could arise if the offsets
      //          cause one to go outside the image boundaries!!!!

      a2nOffset[0] += s_a2nSpecial[0];
      a2nOffset[1] += s_a2nSpecial[1];
    }

  if (e3Ddata_float == m_eType)
    {
      // TODO: make sure image array bounds extents not exceeded

      float *pfTemp;
      i = m_nExt[0] * m_nExt[1] * nLayer;  // the offset into m_pfData
      pfTemp = m_pfData + i;
      for (j = a2nOffset[1]; j < a2nOffset[1]+m_nExt[1]; j++)
	{
	  for (i = a2nOffset[0]; i < a2nOffset[0]+m_nExt[0]; i++)
	    {
	      *pfTemp++ = (poImage->*poImage->prfGetPixel)(i, j);
	    }
	}
    }
  else if (e3Ddata_double == m_eType)
    {
      // TODO: make sure image array bounds extents not exceeded

      double *pflTemp;
      i = m_nExt[0] * m_nExt[1] * nLayer;  // the offset into m_pfData
      pflTemp = m_pflData + i;
      for (j = a2nOffset[1]; j < a2nOffset[1]+m_nExt[1]; j++)
	{
	  for (i = a2nOffset[0]; i < a2nOffset[0]+m_nExt[0]; i++)
	    {
	      *pflTemp++ = (poImage->*poImage->prfGetPixel)(i, j);
	    }
	}
    }
  else if (e3Ddata_ushort == m_eType)
    {
      // Maybe check that image data type supports this!!!!

      unsigned short int *puiTemp;
      i = m_nExt[0] * m_nExt[1] * nLayer;  // the offset into m_puiData
      puiTemp = m_puiData + i;
      for (j = a2nOffset[1]; j < a2nOffset[1]+m_nExt[1]; j++)
	{
	  poImage->nSetNextPixel(a2nOffset[0], j);
	  for (i = a2nOffset[0]; i < a2nOffset[0]+m_nExt[0]; i++)
	    {
	      *puiTemp++ = poImage->uiGetNextPixel();
	    }
	}
    }
  else if (e3Ddata_int == m_eType)
    {
      // TODO: make sure image array bounds extents not exceeded
      cout << "CBF e3Data_int type\n";
      int *pnTemp;
      i = m_nExt[0] * m_nExt[1] * nLayer;  // the offset into m_pfData
      pnTemp = m_pnData + i;
      for (j = a2nOffset[1]; j < a2nOffset[1]+m_nExt[1]; j++)
	{
	  for (i = a2nOffset[0]; i < a2nOffset[0]+m_nExt[0]; i++)
	    {
	      *pnTemp++ = (int) (poImage->*poImage->prfGetPixel)(i, j);
	    }
	}
    }

  return (0);
}


void
C3Ddata::vGetShiftedCentroid(float a3fCent[3],double fIntensity,int nSlice,bool bSet,int a2nAdjCent[2]) {
    const int nMaxSlices = 40;
    static double aa2fObsCent[nMaxSlices][2];
    static double afObsWeight[nMaxSlices];
    static int anObsSlice[nMaxSlices];
    static int nObsSlices;
    int nx=0;
    
    if ((nSlice < 0) && (bSet))
        nObsSlices = 0;
    else if (bSet) {
        
        for (nx = 0; nx < nObsSlices;nx++) {
            if (anObsSlice[nx] == nSlice)
                break;
        };
        if ((nx < 40) && (fIntensity > 0.0)) {
            aa2fObsCent[nx][0] = a3fCent[0];
            aa2fObsCent[nx][1] = a3fCent[1];
            afObsWeight[nx] = max(100,fIntensity);
            anObsSlice[nx] = nSlice;
            if (nx == nObsSlices)
                nObsSlices++;
        };
    } else {
        double  fY[nMaxSlices];
        double  fX[nMaxSlices];
        double  fYSig[nMaxSlices];
        double* apfFuncs[2] = {ms_pfConst,&fX[0] };
        double  a2fArgs[2];
        double  a2fArgsSigma[2];
        int     nCent;
        
        for (nCent = 0; nCent < 2; nCent++) {
            for (nx = 0; nx < nObsSlices;nx++) {
                fY[nx] = aa2fObsCent[nx][nCent];
                fYSig[nx] = 1.0/sqrt(afObsWeight[nx]);
                fX[nx] = anObsSlice[nx];
            };
            if (nObsSlices == 0) {
                a2nAdjCent[nCent] = (int) a3fCent[nCent];
            } else if (nObsSlices == 1) {
                a2nAdjCent[nCent] = (int) aa2fObsCent[0][nCent];
            } else {
                // Do a least squares fit.
                nSolveLeastSquares(2,nObsSlices,&fY[0],&fYSig[0],&apfFuncs[0],&a2fArgs[0],&a2fArgsSigma[0],NULL);
                a2nAdjCent[nCent] = (int) ((nSlice*a2fArgs[1] + a2fArgs[0]));
            };
        };
    };
    return;
};


int
C3Ddata::nBackgroundAvgSD(int nSlice,float* pfBackground,float* pfAverage,float* pfDeviation,float* pfMaxValue,float* pfPeakCentroid,float* pfTotalCentroid) {
    int nx,ny;
    double f0,f1;
    float*  pf0;
    int* pn0;
    int n0,n1;
    
    int nPass;
    double fSumX;
    double  fAverage;
    double  fDeviation;
    double  a3fBackground[3];
    double  fSumX2;
    double  fMaxCentroid;
    int     a2nMaxPixel[2];
    int     nSum;
    int     nSliceSize;
    int     nSliceOffset;
    double  fSatVal,fNominalSatVal;
    double  fBadVal;
    int*    pnMask;
    int*    pnUseInBackground;
    int     nHighCount;
    int     nLowCount;
    int     nBorderCount;
    int     nBackgroundCount;           // Number of pixels actually used in background (<=nSliceCandidatePixels)
    
    const int nStatPasses = 2;  // Number of passes through the satistical loop.


    if (NULL != m_poImage)
        fNominalSatVal = fSatVal = m_poImage->fGetSatValue();
    else
        fNominalSatVal = fSatVal = 999999.0;
    if (m_bIncludeSat)
        fSatVal = 1e20;

    fBadVal = -999.0;


    nSliceSize = m_nExt[0]*m_nExt[1];
    nSliceOffset = nSlice*nSliceSize;
    pnMask = &ms_pnMask[nSliceOffset];



    *pfMaxValue = 0.0;
    a2nMaxPixel[0] = 0;
    a2nMaxPixel[1] = 0;
    
    // Mask out pixels which are in the bounding ellipsoid.
    // We have always traditionally stored the background determination mask in pnUseInBackground,
    // and will continue to do so... but really, it's just the pixels that aren't in the bounding ellipsoid.
    pnUseInBackground = ms_apnScratchInt[0];
    pn0 = pnUseInBackground;
    for (nx = 0; nx < nSliceSize; nx++,pn0++) {
        if (ms_pnMask[nx] & ms_nEllipsoidMask)
            *pn0 = 0;
        else
            *pn0 = 1;
    };

    /*
        Errors that can occur during background estimation:

        1)  Too few pixels to estimate background (after border has been removed).
        2)  The deviations from the planer background do not match the observed detector gain.
        3)  The planar background slope is too high.

    */

  
    
    // Find the statistics.
    // We flag the mask status of each pixel.
    
    for (nPass = 0; nPass< nStatPasses;nPass++) {
        fSumX = 0.0;
        fSumX2 = 0.0;
        nSum = 0;

        pf0 = m_pfData + nSliceOffset;        
        pn0 = pnUseInBackground;
        for (nx = 0; nx<nSliceSize; nx++,pf0++,pn0++) {
            if (*pn0) {
                if (((!nPass) || (!(pnMask[nx] & ms_nHighMask))) && (*pf0<fSatVal) && (*pf0> fBadVal)) {

                    fSumX2+=*pf0 * *pf0;
                    fSumX+= *pf0;
                    nSum++;
                } 
            };
            // Calculate maximum value and peak centroid here.
            // This maximum is over all pixels in shoebox minus a boarder area.
            if (nPass == 0) {
                if (*pf0>*pfMaxValue) {
                    n0 = (nx % m_nExt[0]);
                    n1 = (nx / m_nExt[0]);
                    if ((abs(n0 - m_nExt[0]/2)<m_nExt[0]/4) &&
                        (abs(n1 - m_nExt[1]/2)<m_nExt[1]/4)) {
                        *pfMaxValue = *pf0;
                        a2nMaxPixel[0] = n0;
                        a2nMaxPixel[1] = n1;
                    };
                };
            } 
        };       

        if (nSum>=2)
            fDeviation = max(1.0,sqrt((double)max(0.0,(fSumX2*nSum-fSumX*fSumX)/nSum/(nSum-1))));
        else
            fDeviation = 0.0;
        fAverage = fSumX/max(1,nSum);
        f0 = fAverage + fDeviation*m_fSigmaAboveBackground;
        f1 = fAverage - fDeviation*m_fSigmaBelowBackground;
        
        pf0 = m_pfData + nSliceOffset;
        pn0 = pnUseInBackground;
        nHighCount = 0;
        nLowCount = 0;
        // This code will include saturated pixels.
        for (nx=0;nx<nSliceSize;nx++,pf0++,pn0++) {
            pnMask[nx] &= ~(ms_nHighMask | ms_nLowMask | ms_nSatMask | ms_nBadMask | ms_nBorderMask | ms_nReservedMask);

            // On the first pass (if we are doing two passes) we only need to compute maskes for pixels
            // that are contributing to the background.  On the second pass, we need to compute masks for
            // ALL pixels.

            if ((nPass== nStatPasses-1) || (*pn0)) {
                if (*pf0 >= f0) {
                    if (*pf0>=fSatVal)
                        pnMask[nx] |= ms_nSatMask;
                    if (*pf0>=fNominalSatVal)
                        m_nSatCount ++;
                    pnMask[nx] |= ms_nHighMask;
                    nHighCount++;
                };
                if (*pf0 <= f1) {
                    if (*pf0<=fBadVal)
                        pnMask[nx] |= ms_nBadMask;
                    else {
                        pnMask[nx] |= ms_nLowMask;
                        nLowCount++;
                    };
                }
            };
        };
    };


    // Add the border.  This is done for all pixels.
    nBorderCount = 0;
    for (nPass = 0; nPass< m_nBackgroundBorder;nPass++) {
        int a2nDelta[2];
        for (n1=0,nx=0;n1<m_nExt[1];n1++) {
            for (n0=0;n0<m_nExt[0];n0++,nx++) {
                if (pnMask[nx] & (ms_nHighMask | ms_nBorderMask | ms_nSatMask | ms_nBadMask)) {
                    for (a2nDelta[0] = -1; a2nDelta[0] <= 1; a2nDelta[0]++) {
                        for (a2nDelta[1] = -1; a2nDelta[1] <= 1; a2nDelta[1]++) {
                            if ((a2nDelta[0] + n0 >=0) && (a2nDelta[0] + n0 < m_nExt[0]) &&
                                (a2nDelta[1] + n1 >=0) && (a2nDelta[1] + n1 < m_nExt[1])) {
                                ny = (a2nDelta[0] + n0) + m_nExt[0]*(a2nDelta[1] + n1);
                                if (!(pnMask[ny] & (ms_nHighMask | ms_nBorderMask | ms_nSatMask | ms_nBadMask | ms_nLowMask | ms_nReservedMask))) {
                                    pnMask[ny] |= ms_nReservedMask;
                                    nBorderCount++;
                                };
                            };
                        };
                    };
                };
            };
        };
        for (n1=0,nx=0;n1<m_nExt[1];n1++) {
            for (n0=0;n0<m_nExt[0];n0++,nx++) {
                if (pnMask[nx] & ms_nReservedMask) {
                    pnMask[nx] &= ~(ms_nReservedMask);
                    pnMask[nx] |= ms_nBorderMask;
                };
            };
        };
    };
    /*
    // Remove the 'high' status on all pixels that don't have a minimum # of neighbors.
    for (n1=0,nx=0;n1<m_nExt[1];n1++) {
        for (n0=0;n0<m_nExt[0];n0++,nx++) {
            if (pnMask[nx] & (ms_nHighMask)) {
                int nHighFound = 0;
                int a2nDelta[2];
                for (a2nDelta[0] = -1; a2nDelta[0] <= 1; a2nDelta[0]++) {
                    for (a2nDelta[1] = -1; a2nDelta[1] <= 1; a2nDelta[1]++) {
                        if ((a2nDelta[0] + n0 >=0) && (a2nDelta[0] + n0 < m_nExt[0]) &&
                            (a2nDelta[1] + n1 >=0) && (a2nDelta[1] + n1 < m_nExt[1])) {
                            ny = (a2nDelta[0] + n0) + m_nExt[0]*(a2nDelta[1] + n1);
                            if (pnMask[ny] & ms_nHighMask) 
                                nHighFound++;
                        };
                    };
                };
                if (nHighFound < 3)
                    pnMask[nx] &= (~ms_nHighMask);
            };
        };
    };
    */
    
    ny = ms_nBorderMask | ms_nSatMask | ms_nBadMask | ms_nHighMask;
    
    if (m_bFitToPlanar) {
        // Calculate the planar background. 
        // This only uses pixels that are used in background calculations.
        
        double a3x3fMatA[3][3];
        double a3x3fMatAInv[3][3];
        double a3x3fMatAAdd[3][3];
        double a3fMatb[3];
        vZeroMat(3,3,&a3x3fMatA[0][0]);
        vZeroMat(3,1,&a3fMatb[0]);
        
        nBackgroundCount = 0;
        pn0 = pnUseInBackground;
        
        for (n1=0,nx=0;n1<m_nExt[1];n1++) {
            for (n0=0;n0<m_nExt[0];n0++,nx++,pn0++) {
                if (*pn0) {
                    if (!(pnMask[nx] & ny)) {
                        f0 = ((double) n0)/m_nExt[0];
                        f1 = ((double) n1)/m_nExt[1];
                        nBackgroundCount++;
                        
                        a3x3fMatAAdd[0][0] = f0*f0;
                        a3x3fMatAAdd[1][0] = f1*f0;
                        a3x3fMatAAdd[2][0] = f0;
                        a3x3fMatAAdd[1][1] = f1*f1;
                        a3x3fMatAAdd[2][1] = f1;
                        a3x3fMatAAdd[2][2] = 1;
                        a3x3fMatAAdd[0][1] = a3x3fMatAAdd[1][0];
                        a3x3fMatAAdd[0][2] = a3x3fMatAAdd[2][0];
                        a3x3fMatAAdd[1][2] = a3x3fMatAAdd[2][1];
                        vAddVecNDVecND(9,&a3x3fMatAAdd[0][0],&a3x3fMatA[0][0],&a3x3fMatA[0][0]);
                        a3fMatb[0] += m_pfData[nSliceOffset + nx]*f0;
                        a3fMatb[1] += m_pfData[nSliceOffset + nx]*f1;
                        a3fMatb[2] += m_pfData[nSliceOffset + nx];
                    };
                };
            };
        };
        
        if (0.0 == fInvMat3D(&a3x3fMatA[0][0],&a3x3fMatAInv[0][0])) {
            return 1 << Crefln::ms_nErrorBackground;
        };
        
        vMulMat3DVec3D(a3x3fMatAInv,a3fMatb,a3fBackground);
    } else {
        a3fBackground[0] = 0.0;
        a3fBackground[1] = 0.0;
        a3fBackground[2] = fAverage;
    };


    
    *pfMaxValue = m_pfData[nSliceOffset + a2nMaxPixel[0] + a2nMaxPixel[1]*m_nExt[0]] - (a3fBackground[0]*a2nMaxPixel[0]/m_nExt[0] + a3fBackground[1]*a2nMaxPixel[1]/m_nExt[1] + a3fBackground[2]);


    // Compute the deviations from the planar background. 
    fDeviation = 0.0;
    pn0 = pnUseInBackground;
    fMaxCentroid = 0.0;
    pfPeakCentroid[0] = 0.0;
    pfPeakCentroid[1] = 0.0;

    float a2fCompeteCentroid[2];
    float fTotalCentroidRadius;
    float fMaxTotalCentroid;
    float fTotalCentroid;
    int   nTotalCentroid;

    // Are we spot chasing?  If so, we are interested in the 
    // total spot centroid.
    if ((pfTotalCentroid)) {
        pfTotalCentroid[0] = 0.0;
        pfTotalCentroid[1] = 0.0;
        pfTotalCentroid[2] = 0.0;
        fTotalCentroid = 0.0;
        nTotalCentroid = 0;
        fMaxTotalCentroid = 0.0;
    };

    vGetCompetingSpot(nSlice,&a2fCompeteCentroid[0],&a2fCompeteCentroid[1]);
    fTotalCentroidRadius = 
        ((a2fCompeteCentroid[0] - m_a3fAbsCentroid[0])*(a2fCompeteCentroid[0] - m_a3fAbsCentroid[0]) +
        (a2fCompeteCentroid[1] - m_a3fAbsCentroid[1])*(a2fCompeteCentroid[1] - m_a3fAbsCentroid[1]))/4.0;        

    nBackgroundCount = 0;
    ny = ms_nSatMask | ms_nBadMask | ms_nHighMask;

    for (n1=0,nx=0;n1<m_nExt[1];n1++) {
        for (n0=0;n0<m_nExt[0];n0++,nx++,pn0++) {
            f0 = a3fBackground[0]*n0/m_nExt[0] + a3fBackground[1]*n1/m_nExt[1] + a3fBackground[2];
            f1 = m_pfData[nSliceOffset + nx] - f0;
            if (*pn0) {
                if (!(pnMask[nx] & ny)) {
                    fDeviation += (f1*f1);
                    nBackgroundCount++;
                };
            };
            
            f0 = 
                (n0 + m_nOrigOffset[0] - m_a3fAbsCentroid[0])*(n0 + m_nOrigOffset[0] - m_a3fAbsCentroid[0]) +
                (n1 + m_nOrigOffset[1] - m_a3fAbsCentroid[1])*(n1 + m_nOrigOffset[1] - m_a3fAbsCentroid[1]);

            
            if ((f1 >= 0.8*(*pfMaxValue)) &&
                (abs(n0 - m_nExt[0]/2)<m_nExt[0]/4) &&
                (abs(n1 - m_nExt[1]/2)<m_nExt[1]/4) &&
                (f0 < fTotalCentroidRadius)
                ) {
                
                fMaxCentroid += f1;
                pfPeakCentroid[0] += f1*n0;
                pfPeakCentroid[1] += f1*n1;
            };

            if (ms_pnMask[nx] & ms_nEllipsoidMask) {
                
                if ((pfTotalCentroid) && 
                    (pnMask[nx] & ms_nHighMask) && 
                    (!(pnMask[nx] & ms_nSatMask))) {
                    
                    fTotalCentroid += f1;
                    nTotalCentroid ++;
                    pfTotalCentroid[0] += f1*n0;
                    pfTotalCentroid[1] += f1*n1;
                    fMaxTotalCentroid = max(fMaxTotalCentroid,f1);
                };
            };
        };
    };

    if (nBackgroundCount >= 2)
        fDeviation = sqrt((double)(fDeviation/(nBackgroundCount-1)));
    else
        fDeviation = 0.0;

    // Calculate the peak centroid.  This is used for sharpness calculations.
    if (fMaxCentroid!=0.0)  {
        pfPeakCentroid[0]/= fMaxCentroid;
        pfPeakCentroid[1]/= fMaxCentroid;
    };

    // Calculate the total centroid.  This is used for spot chasing.
    if ((pfTotalCentroid) && (fTotalCentroid>0.0)) {
        // Make sure that our spot centroid is justified.
        // It might not be if we don't have a high enough peak.
        pfTotalCentroid[2] = fTotalCentroid;
    } else if (pfTotalCentroid) {
        pfTotalCentroid[2] = 0.0;
    };

    *pfDeviation = (float) fDeviation;
    *pfAverage = (float) (a3fBackground[0]*0.5 + a3fBackground[1]*0.5 + a3fBackground[2]);
    vCopyVec3D(&a3fBackground[0],pfBackground);

  return 0;    
};


#ifdef DEBUG_C3DDATA 
#define DEBUG_C3DDATA_DIM0 3000
#define DEBUG_C3DDATA_DIM1 3000

#ifndef DEBUG_C3DDATA_DIM0
#error Must define DEBUG_C3DDATA_DIM0
#endif
#ifndef DEBUG_C3DDATA_DIM1
#error Must define DEBUG_C3DDATA_DIM1
#endif

/*
Only use these features when you want to debug Cfind stuff outside of a Cintegrate
context.
*/

#include "Cintegrate.h"

Cintegrate g_oRogueIntegrate;
int g_nRogueCount = 0;
int g_nRogueImage = -1;
float g_a3fRogueCent[3];


int
C3Ddata::nAddRogue(int nRogueImage) {

    // Was the rouge data initialized?

    if ((g_nRogueCount==0) || (g_nRogueImage!=nRogueImage)) {
        g_oRogueIntegrate.vSetRogueData(0,0,"C3Ddata???",nRogueImage,-1);
        g_nRogueImage = nRogueImage;
        g_oRogueIntegrate.vSetRogueImageDim(DEBUG_C3DDATA_DIM0,DEBUG_C3DDATA_DIM1);
    };

    nLoadPrint(TRUE);
    m_nOrigOffset[2] += nRogueImage;
    g_oRogueIntegrate.nAddRogue(2,0,*this,pfGetPrintMask());
    m_nOrigOffset[2] -= nRogueImage;
    g_nRogueCount++;

    return 0;
};

void C3Ddata::vTermRogue() {
    g_oRogueIntegrate.vTermRogue();
    g_nRogueCount = 0;
};
#else

int
C3Ddata::nAddRogue(int nImageNumber) {
    return 0;
};

void C3Ddata::vTermRogue() {
};
#endif



int
C3Ddata::nCalcGetPeakInfo(float *pfCentroid,
			  float *pfInt,  float *pfSigmaI,
			  float *pfAvgOut, float *pfSDOut,
			  float *pfSize)
{
  // Calculate centroid of a peak in the data array
  // Modified input: pfCentroid[0], pfCentroid[1], pfCentroid[2] the
  //                 coordinates of a pixel near the peak centroid.
  //                 The coordinates are relative to the data array origin.
  //                 That is, (0,0,0) is the coordinate of the
  //                 first pixel in the data array.
  //
  // Return new Centroid, Intensity, sigma, background average
  // and background sd.
  //
  // Return 0    if no errors
  //        != 0 some kind of error


  float a3fCent[3];
  int   a3nCent[3];
  int   nStat;
  int i;

  // Get the putative center peak pixel wrt to origin of this shoebox
  // Probably should check input to centroid to see if it makes sense

  for (i = 0; i < 3; i++)
    {
      a3fCent[i] = pfCentroid[i];
      a3nCent[i] = (int) pfCentroid[i];
    }

  // nCalcMask expects an input pixel with origin offsets applied

  for (i = 0; i < 3; i++)
    {
      a3fCent[i] += m_nOrigOffset[i];
      a3nCent[i] += m_nOrigOffset[i];
    }
#ifdef DEBUG_C3DDATA 
  vCopyVec3D(a3fCent,g_a3fRogueCent);
#endif

  // Assuming that the user calling this has not set any of the 
  // variables providing crystalographic data.
  //((C3DdataInput&) *this).vInit();

  ((C3DdataOutput&) *this).vInit();

  m_bSpotChase = FALSE;
  nStat = nCalcMask(a3fCent,m_nExt[2]-0.0001);
  

  // m_a3fCentroid and m_a3fPeakSize are calculated in nCalcMask

  for (i = 0; i < 3; i++)
    {
      pfCentroid[i] = m_a3fCentroid[i];
    }

  if (NULL != pfSize)
    {
      for (i = 0; i < 3; i++)
	{
	  pfSize[i] = m_a3fPeakSize[i];
	}
    }

  if (NULL != pfInt)
    *pfInt = fGetPeakIntensity();
  if (NULL != pfSigmaI)
    *pfSigmaI = fGetPeakIntensitySigma();
  if (NULL != pfAvgOut)
    *pfAvgOut = fGetPerPixelBackground();
  if (NULL != pfSDOut)
    *pfSDOut = fGetPerPixelBackgroundSigma();

  return (nStat);
}

int
C3Ddata::nList(const int nLayer = 0)
{
  int i, j;
  float fTemp;
  for (j = 0; j < m_nExt[1]; j++)
    {
      for (i = 0; i < m_nExt[0]; i++)
	{
	  fTemp = fGetValue(i, j, nLayer);
	  cout << fTemp << " ";
	}
      cout << endl;
    }
  cout << endl << flush;
  return (0);
}

void C3Ddata::vCheckMemory() {
    int nNumPoints= m_nExt[0] * m_nExt[1] * m_nExt[2] + m_nScratchFloatSize*m_nExt[2];

    if (nNumPoints>m_nAlloc) {
	  delete [] m_pflData;
	  m_puiData = NULL;
	  m_pfData  = NULL;
	  m_pflData = NULL;
	  m_pnData  = NULL;
	  m_nAlloc  = 0;
	  vGetMemory(eGetType());
    };
};

int
C3Ddata::nConvertToFloat(Cimage* poImageIn, C3Ddata *po3DOther)
{
  // Convert a shoebox of type e3Ddata_ushort to e3Ddata_float
  // New (2008-02-04): 
  //     Shoebox could be float already and the right thing will be done!
  // We really need to know something about the original image
  // format in order to do this.
  // If (po3DOther != NULL), then place converted results into po3DOther

  if ( (e3Ddata_float == m_eType) && (NULL == po3DOther) )
    {
      // Returns right away if already of float type
      return (0);      
    }
  if (NULL != poImageIn)
    {
      m_poImage = poImageIn;
    }
  if (NULL == m_poImage)
    {
      cout << "C3Ddata ERROR, no image properties available for conversion!\n";
      return (4);
    }

  int i;
  int nAlloc;
  int nNumPoints;
  unsigned short int *puiTemp, *puiData;

  float  *pfTemp, *pfData;
  double *pdTemp;
  int    *pnTemp, *pnData;

  nAlloc = m_nExt[0] * m_nExt[1] * m_nExt[2] + m_nScratchFloatSize*m_nExt[2];
  nNumPoints = m_nExt[0] * m_nExt[1] * m_nExt[2];

  // Save pointer to our original data for looping below and for deleting
  // Do it here, since value of m_puiData may change before we need it.

  puiData = m_puiData;
  puiTemp = puiData;
  pfData  = m_pfData;
  pnData  = m_pnData;

  // Now decide where to put the converted data and proceed accordingly ...

  if (NULL != po3DOther)
    {
      //cout << "C3Ddata: to other float\n";
      
      // We want to put the results into another shoebox and leave this
      // shoebox unchanged, other shoebox must be of type float, so ...
      
      if (   (po3DOther->m_nAlloc < nAlloc)
          || (e3Ddata_float != po3DOther->eGetType()) )
	{
          // Need to allocate more space
          
          delete [] po3DOther->m_pflData;
          po3DOther->m_puiData = NULL;
          po3DOther->m_pfData  = NULL;
          po3DOther->m_pflData = NULL;
          po3DOther->m_nAlloc  = 0;
          po3DOther->m_nExt[0] = m_nExt[0];
          po3DOther->m_nExt[1] = m_nExt[1];
          po3DOther->m_nExt[2] = m_nExt[2];
          po3DOther->vGetMemory(e3Ddata_float);
          po3DOther->m_eType = e3Ddata_float;
	}
      
      // Copy all this's info to the *po3DOther object
      // ::nCopyFrom() for all data types

      (void) po3DOther->nCopyFrom(*this);

      if (NULL != poImageIn)
	po3DOther->m_poImage = poImageIn;
      else
	po3DOther->m_poImage = m_poImage;

      // Set our initial results pointer to the other's address

      pfTemp = po3DOther->m_pfData;
      pdTemp = po3DOther->m_pflData;
      pnTemp = po3DOther->m_pnData;
    }

  // Now convert the unsigned short int to float
  // using the conversion procedure of the image

  if (e3Ddata_float == m_eType)
    {
      if (NULL == po3DOther)
	{
	  // Do nothing because it's copy in place 
	  // and the data_type is already float
	  // But because of logic above, should never get here
	  cout << "WARNING C3Ddata: do nothing!\n";
	}
      else
	{
	  //cout << "C3Ddata: float to new float\n";
	  memcpy(pfTemp, pfData, sizeof(float) * nNumPoints);
	}
    }
  else if (e3Ddata_ushort == m_eType)
    {
      //cout << "C3Ddata: ui to new float\n";
      for (i = 0; i < nNumPoints; i++)
	{
	  *pfTemp++ = (m_poImage->*m_poImage->prfCvtUItoFloat)(*puiTemp++);
	}
      if (NULL == po3DOther)
	{
	  // Now free the unsigned short int memory

	  pfTemp = (float*)puiData;
	  pdTemp = (double*)puiData;
	  vFreeMemory(&pdTemp);
	}
    }
  return (0);
}

void
C3Ddata::vZero(void)
{
    // Zero the data array
    // The binary representation of (float) 0.0 had better be the same as
    // (unsigned short int) 0.
    
    int i;
    int nNumPoints = m_nExt[0]*m_nExt[1]*m_nExt[2];
    if (e3Ddata_ushort == m_eType)
    {
        unsigned short int *puiTemp = m_puiData;
        for (i = 0; i < nNumPoints; i++)
        {
            *puiTemp++ = 0;
        }
    }
    else if (e3Ddata_double == m_eType)
    {
        double *pflTemp = m_pflData;
        for (i = 0; i < nNumPoints; i++)
        {
            *pflTemp++ = 0.0;
        }
    }
    else
    {
        float *pfTemp = m_pfData;
        for (i = 0; i < nNumPoints; i++)
        {
            *pfTemp++ = 0.0;
        }
    }
}

int
C3Ddata::nCalcGetMax(int *pn0Max, int *pn1Max, int *pn2Max,
		     float *pfValue)
{
  int i, j, k;

  *pn0Max = *pn1Max = *pn2Max = 0;
  if (e3Ddata_float == m_eType)
    {
      float *pfTemp = m_pfData;
      float fTest   = *pfTemp;
      for (k = 0; k < m_nExt[2]; k++)
	{
	  for (j = 0; j < m_nExt[1]; j++)
	    {
	      for (i = 0; i < m_nExt[0]; i++)
		{
		  if (fTest < *pfTemp)
		    {
		      *pn0Max = i;
		      *pn1Max = j;
		      *pn2Max = k;
		      fTest   = *pfTemp;
		      *pfValue = fTest;
		    }
		  pfTemp++;
		}
	    }
	}
      return (0);
    }
  else if (e3Ddata_double == m_eType)
    {
      double *pflTemp = m_pflData;
      double fTest   = *pflTemp;
      for (k = 0; k < m_nExt[2]; k++)
	{
	  for (j = 0; j < m_nExt[1]; j++)
	    {
	      for (i = 0; i < m_nExt[0]; i++)
		{
		  if (fTest < *pflTemp)
		    {
		      *pn0Max = i;
		      *pn1Max = j;
		      *pn2Max = k;
		      fTest   = *pflTemp;
		      *pfValue = fTest;
		    }
		  pflTemp++;
		}
	    }
	}
      return (0);
    }
  else if (e3Ddata_ushort == m_eType)
    {
      unsigned short int *puiTemp = m_puiData;
      unsigned short int uiTest = *puiTemp;
      for (k = 0; k < m_nExt[2]; k++)
	{
	  for (j = 0; j < m_nExt[1]; j++)
	    {
	      for (i = 0; i < m_nExt[0]; i++)
		{
		  if (uiTest < *puiTemp)
		    {
		      *pn0Max = i;
		      *pn1Max = j;
		      *pn2Max = k;
		      uiTest  = *puiTemp;
		      *pfValue = (float) uiTest;
		    }
		  puiTemp++;
		}
	    }
	}
      return (0);
    }
  return (-2);
};

void
C3Ddata::vAllocScratch(int nAlgorithm)
{
  // Make certain that member scratch variables are
  // dimensioned large enough for the input shoebox.

    if (nAlgorithm == 1) {
        int i,j;
        j = m_nExt[0]*m_nExt[1];
        if (j > ms_nAllocScratchInt) {
            i = 0;
            if (NULL != ms_apnScratchInt[i])
                delete [] ms_apnScratchInt[i];
            ms_apnScratchInt[i] = new int[j];
            memset(ms_apnScratchInt[i],0,sizeof(int)*j);
            
            i = 1;
            if (NULL != ms_apnScratchInt[i])
                delete [] ms_apnScratchInt[i];
            ms_apnScratchInt[i] = new int[j*16];
            memset(ms_apnScratchInt[i],0,sizeof(int)*j*16);

            ms_nAllocScratchInt = j;

            i = 0;
        };
    };
    
    if (nAlgorithm == 1) {
        if (m_nAlloc > ms_nAllocSmooth)
        {
            // Smooth and Mask are too small, so delete them and make a larger ones.
            // In this way, we don't spend time allocating and freeing
            // these, but keep re-using them, but let
            // them grow to the maximum shoebox size.
            
            if (NULL != ms_pfMask)   delete [] ms_pfMask;
            if (NULL != ms_pfPrint)   delete [] ms_pfPrint;
            if (NULL != ms_pnMask)   delete [] ms_pnMask;
            ms_nAllocSmooth = m_nAlloc;
            ms_pnMask       = new int   [ms_nAllocSmooth];
            ms_pfMask       = new float [ms_nAllocSmooth];
            ms_pfPrint      = new float [ms_nAllocSmooth];

        }
    };

}

void
C3Ddata::vFreeScratch() 
{
    if (NULL != ms_pfMask)
    {
        delete [] ms_pfMask;
        ms_pfMask = NULL;
    }
    if (NULL != ms_pfPrint)
    {
        delete [] ms_pfPrint;
        ms_pfPrint = NULL;
    }
    if (NULL != ms_pnMask)
    {
        delete [] ms_pnMask;
        ms_pnMask = NULL;
    }
    ms_nAllocSmooth = 0;
    

    int i;
    for (i = 0; i <2; i++)
    {
        if (NULL != ms_apnScratchInt[i])
            delete [] ms_apnScratchInt[i];
        ms_apnScratchInt[i] = NULL;
    };
    ms_nAllocScratchInt = 0;

};


void
C3Ddata::vFreeInit(void)
{
  // Make sure all allocated memory is freed before initializing

  vFreeDelete();

  // Now initialize

  if (NULL == ms_pnFreeSize)
    ms_pnFreeSize = new int [C3DDATA_MAX_FREE];
  if (NULL == ms_ppfFreeFloat)
    ms_ppfFreeFloat = new float* [C3DDATA_MAX_FREE];


  register int i;
  for (i = 0; i < C3DDATA_MAX_FREE; i++)
    {
      ms_pnFreeSize[i]    = 0;
      ms_ppfFreeFloat[i]  = NULL;
    }
  ms_nFreeCount = 0;
}

void
C3Ddata::vFreeDelete(void)
{
//  cout << "C3Ddata::vFreeDelete called\n";
//  cout << "C3Ddata re-use stats:  Hard: " << ms_nHardCount
//                     << "  Soft: " << ms_nSoftCount << endl;

  if (NULL != ms_ppfFreeFloat)
    {
      register int i;
      for (i = 0; i < ms_nFreeCount; i++)
	{
	  if (NULL != ms_ppfFreeFloat[i])
	    {
	      delete [] ms_ppfFreeFloat[i];
	      ms_pnFreeSize[i]   = 0;
	      ms_ppfFreeFloat[i] = NULL;
	    }
	}
      delete [] ms_ppfFreeFloat;
      ms_ppfFreeFloat = NULL;
    }
  if (NULL != ms_pnFreeSize)
    {
      delete [] ms_pnFreeSize;
      ms_pnFreeSize = NULL;
    }
  ms_nFreeCount = 0;
}

void
C3Ddata::vGetMemory(const e3Ddata_types eType)
{
  // Look for a piece of memory in the free list
  // When called, nSize is the size of the data required in number of grid
  // points.
    // On return, m_nAlloc is the total number of elements requested.
    
    int nSize;
    int nDataSize;
    int nRequested;
    int nx;

    nRequested = m_nExt[0]*m_nExt[1]*m_nExt[2];
    
    if (0 >= nRequested)
        return;
    
    if (NULL != m_pfData)
        cout << "ERROR is C3Ddata memory leak in vGetMemory!\n";
    
    if ((sizeof(const float)!=4) || (sizeof(const double) != 8) || (sizeof(const unsigned short int)!=2))
        cout << "ERROR:  Sizes are not expected.\n";
    
    if (e3Ddata_ushort == eType)
        nSize = nRequested / 2 + (nRequested % 2) ;
    else if (e3Ddata_double == eType)
        nSize = nRequested * 2;      
    else
        nSize = nRequested;

    nDataSize = nSize + (nSize % 2);
    nSize = nDataSize + m_nScratchFloatSize*m_nExt[2]+ sizeof(double) + sizeof(double);
       
    double *pdTemp;
    pdTemp     = new double [nSize / 2];
    m_pfData   = (float *)pdTemp;
    m_puiData  = (unsigned short int*) pdTemp;
    m_pflData  = (double*) pdTemp;
    m_pnData   = (int *)pdTemp;

    ms_nHardCount++;
    m_nAlloc = m_nExt[0]*m_nExt[1]*m_nExt[2] + m_nExt[2]*m_nScratchFloatSize;

    for (nx =0; nx < m_nScratchFloatSize;nx++) {
        m_apfScratchFloat[nx] = &m_pfData[nDataSize + nx*m_nExt[2]];
    };
    memset(m_apfScratchFloat[0],0,m_nExt[2]*m_nScratchFloatSize*sizeof(float));
    
    // If the free list is full, or if the free list contains
    // more than 1/3 the number of objects
    // delete it all to start fresh
    
    if (   (C3DDATA_MAX_FREE == ms_nFreeCount)
        || (3 * ms_nFreeCount > ms_nNumObjects) )
    {
        // Clear out the free list since nothing was suitable
        
        cout << "just BEFORE vFreeInit new float nsize: " << nSize
            << endl << flush;
        vFreeInit();
        cout << "just AFTER vFreeInit new float nsize: " << nSize
            << endl << flush;
    }
    vAllocScratch(1);    
}

void
C3Ddata::vFreeMemory(double **ppdData)
{
  // Free memory for this object by placing in our free list

  if (NULL == *ppdData) return;

  delete [] *ppdData;
  *ppdData  = NULL;
}

int
C3Ddata::nWrite(const Cstring& rsFilename, const int nDir, Crefln *poRefln)
{
  int nStat;
  Cimage *poImage;

  poImage = new Cimage(this, nDir, poRefln);
  nStat   = poImage->nWrite(rsFilename);

  delete poImage;

  return nStat;
}

int
C3Ddata::nCopyFrom(const C3Ddata& ro3DdataIn)
{
  // Copy all but the m_eType, m_nAlloc and the actual pixel data
  // from another C3Ddata object to this one.
  // Returns the number of points in the input C3Ddata object.

  int   i;                     // Loop counter
  int   nNeeded = 1;
  for (i = 0; i < 3; i++)
    {
      // Copy extents and offsets

      m_nExt[i]         = ro3DdataIn.m_nExt[i];
      m_nFilled[i]      = ro3DdataIn.m_nFilled[i];
      m_nOrigOffset[i]  = ro3DdataIn.m_nOrigOffset[i];
      m_a3fCentroid[i]  = ro3DdataIn.m_a3fCentroid[i];
      m_a3fPeakSize[i]  = ro3DdataIn.m_a3fPeakSize[i];
      nNeeded           = nNeeded * ro3DdataIn.m_nExt[i];  // Get total pixels
    }
  // Copy all parameters over as a block.
  C3DdataInput::vCopy(ro3DdataIn);
  C3DdataOutput::vCopy(ro3DdataIn);
  
  int nRequested,nSize,nDataSize;

  for (i = 0; i < m_nScratchFloatSize; i++) {
      nRequested = (m_nExt[0]*m_nExt[1]*m_nExt[2]);
      if (e3Ddata_ushort == m_eType)
          nSize = nRequested / 2 + (nRequested % 2) ;
      else if (e3Ddata_double == m_eType)
          nSize = nRequested * 2;      
      else
          nSize = nRequested;
      nDataSize = nSize + (nSize % 2);

      m_apfScratchFloat[i] = &m_pfData[nDataSize];

      m_apfScratchFloat[i] += i*m_nExt[2];
      memcpy(m_apfScratchFloat[i],ro3DdataIn.m_apfScratchFloat[i],sizeof(float)*m_nExt[2]);
  };

  return (nNeeded);
}



bool
C3Ddata::bIntersect(const C3Ddata& oDataOther,int a6nIntersection[6]) {
    int nx;
    bool bConstraintFound;

    int a3nMax[3];
    int a3nMin[3];
    
    for (nx=0;nx<3;nx++) {
        bConstraintFound = FALSE;
        a3nMin[nx] = m_nOrigOffset[nx];
        a3nMax[nx] = m_nOrigOffset[nx] + m_nExt[nx] -1;
        if ((oDataOther.m_nOrigOffset[nx]>=a3nMin[nx]) && (oDataOther.m_nOrigOffset[nx]<=a3nMax[nx])) {
            bConstraintFound = TRUE;
            a3nMin[nx] = oDataOther.m_nOrigOffset[nx];
        };
        if ((oDataOther.m_nOrigOffset[nx] + oDataOther.m_nExt[nx] - 1 >=a3nMin[nx]) && 
            (oDataOther.m_nOrigOffset[nx] + oDataOther.m_nExt[nx] - 1 <=a3nMax[nx])) {
            bConstraintFound = TRUE;
            a3nMax[nx] = oDataOther.m_nOrigOffset[nx] + oDataOther.m_nExt[nx] - 1;
        };
        if (!bConstraintFound)
            return FALSE;
    };
    for (nx=0;nx<3;nx++) {
        a6nIntersection[nx]  = a3nMin[nx];
        a6nIntersection[nx + 3] = a3nMax[nx];
    };
    return TRUE;

};


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////// Latest C3Ddata routines added by tjn ////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


// Macros to transform coordinates.

#define C2N(n0,n1,n) n=((n0) + (n1)*m_nExt[0])
#define N2C(n,n0,n1) { (n1) = ((n)/m_nExt[0]); (n0) = ((n) - (n1)*m_nExt[0]); }

#define BELOW_SIG_CUTOFF
#define CENTROID_MASK_PASTE
#define TAIL_REMOVAL
#define INCREASING_TAILS_CHECK

int C3Ddata::nCalcMask(const float a3fAbsCentroid[3],const float fRotWidth) 
{

    int   nx,ny,nz;
    double f0,f1,f2;
    int   n0,n1;
    float a3fTemp1[3];


    // General variables.
    int     nSize;
    int     nSliceSize;
    int     a2nSliceCent[2];   // Center for the given slice.
    int     a3nCent[3];        // Offset center calculated from absolute center.
    float   a3fCent[3];        // Offset center calculated from absolute center.
    int     nPeakStartSearch;  // Set either to nMosPeakStartSearch or 0
    int     nPeakEndSearch;    // Set either to nMosPeakEndSearch or m_nExt[2]-1
    int     nMosPeakStartSearch;//Maximum value that m_nPeakStart could get (based on input rotation width.
    int     nMosPeakEndSearch; // Minimum value that m_nPeakEnd could get (based on input rotation width.
    int     nStatSlice;        // Return status for a particular slice.
    int     nStatSlices;       // Global return status.  Returned if all other errors are 0.
    int     a2x2nSearch[2][2]; // Ranges of search for a given slice.
    int     nWarningMask;      // Error bits that should be treated as warnings.


    // Variables that control looping through slices.

    int nSlices;            // Number of slices added.
    int nPositiveSlices;    // Number of slices whose [2] entry is greater than the center slice.
    int nNegativeSlices;    // Number of slices whose [2] entry is less than the center slice.
    int nSlice;             // The slice under consideration.
    int nLastSlice;         // Used when adding a new slice.
    int nSliceOffset;       // Offset in pixels for the beggining of the slice.
    bool bSliceAdded;       // Was the last slice added?
    bool bContinueLoop;     // Continue looping at finding slices.
    bool bPadOnUpperEnd;    // Should we pad on the upper end of the slice.
    bool bPadOnLowerEnd;    // Should we pad on the lower end of the slice.
    double fPeakMaxValue;
    double fPeakMinValue;
    float* afSliceIntensity;            // Slice intensities for each slice.
    float* afSliceSigma;                        // Slice sigma;
    float* afSliceIntensityProfit;      // Slice intensities for each slice.
    float* afSliceSigmaProfit;          // Slice sigma (profit);
    float* afSliceSigmaBackground;  // Slice sigma from background.
    float* afSliceBackground;
    float* afSliceStatus;
    float* afSliceStat;

    // Variables used in centroid determination/ peak determination.
    double fSumXIntensity;
    double fSumYIntensity;
    double fSumZIntensity;
    double fSumSliceXIntensity;
    double fSumSliceYIntensity;
    double fSumSliceVarianceNumerator;
    double fSumSliceVarianceDenominator;
    double fSumNonnegativeIntensity;
    double fSumSliceIntensity;
    double fSumSliceNonnegativeIntensity;
    double fSumSliceBackgroundVariance;
    double fSumSliceIntensityVariance;
    double fSumBackground;
    double fSum0Dim;
    double fSum1Dim;
    double fSumIntensityVariance;
    double fSumBackgroundVariance;
    double fSumSlicePixelsOnEdge0;
    double fSumSlicePixelsOnEdge1;
    

      // This will contain cumulative errors which do not cause an exit.
    if (e3Ddata_float != m_eType)
        return ( (1<<Crefln::ms_nErrorBackgroundSD) );
    
    if ((ABS(a3fAbsCentroid[0] - 760) < 20) && (ABS(a3fAbsCentroid[1] - 2335) < 20))
    {
    };

    
    // Clear all data in the output fields.  These shall be recomputed.
    ((C3DdataOutput&) *this).vInit();
    // If we are doing bootstraping, get the putative center peak pixel wrt to origin of this shoebox

    for (nx=0; nx<3; nx++)
      {
        a3fCent[nx]          = a3fAbsCentroid[nx] - m_nOrigOffset[nx];
        a3nCent[nx]          = (int) a3fCent[nx];
        m_a3fAbsCentroid[nx] = a3fAbsCentroid[nx];
      }
    
    m_fRotWidth = fRotWidth;        
    

    nStatSlices     = 0;
    nStatSlice      = 0;

    fPeakMaxValue   = 0.0;
    fPeakMinValue   = 99999999;
    nSize           = m_nExt[0]*m_nExt[1]*m_nExt[2];
    nSliceSize      = m_nExt[0]*m_nExt[1];
    nSlices         = 0;
    nPositiveSlices = 0;
    nNegativeSlices = 0;
    nSlice          = a3nCent[2];
    m_nCentroid     = a3nCent[2];
    bContinueLoop   = TRUE;
    bPadOnUpperEnd  = TRUE;
    bPadOnLowerEnd  = TRUE;
    nWarningMask    = ~((  (1 << Crefln::ms_nErrorUnaccountedPix) 
                         | (1 << Crefln::ms_nErrorTooFewPix)));

    afSliceIntensity       = m_apfScratchFloat[m_nScratchFloatIntensity];
    afSliceSigma           = m_apfScratchFloat[m_nScratchFloatSigma];
    afSliceIntensityProfit = m_apfScratchFloat[m_nScratchFloatIntensityProfit];
    afSliceSigmaProfit     = m_apfScratchFloat[m_nScratchFloatSigmaProfit];
    afSliceSigmaBackground = m_apfScratchFloat[m_nScratchFloatSigmaBackground];
    afSliceBackground      = m_apfScratchFloat[m_nScratchFloatBackground];
    afSliceStatus          = m_apfScratchFloat[m_nScratchFloatStatus];
    afSliceStat            = m_apfScratchFloat[m_nScratchFloatStat];

    for (nx=0; nx<m_nExt[2]; nx++)
      {
        afSliceIntensity[nx]       = 0.0;
        afSliceSigma[nx]           = 0.0;
        afSliceIntensityProfit[nx] = 0.0;
        afSliceSigmaProfit[nx]     = 0.0;
        afSliceSigmaBackground[nx] = 0.0;
        afSliceStatus[nx]          = 0;
        afSliceStat[nx]            = 0;
      }

    // Off bounds checking for centroid.

    int a3nBoundsOff[3]    = {3,3,0};
    int a3nBoundsErrors[3] = { Crefln::ms_nErrorOffEdge0, 
                               Crefln::ms_nErrorOffEdge1,
                               Crefln::ms_nErrorOffEdge2 };

    for (nx = 0; nx < 3; nx++)  
      {
        // Check for out-of-bounds
        if ((a3nCent[nx]<0) || (a3nCent[nx]>=m_nExt[nx]))
          {
            // Instant death!
            return (1 << a3nBoundsErrors[nx]) | nStatSlices;
          }

        // Check on-edge       
        if (   (a3nBoundsOff[nx] > a3nCent[nx]) 
            || (m_nExt[nx] - 1 - a3nBoundsOff[nx] < a3nCent[nx]) )
          {
            if (nx < 2) 
              return (1 << a3nBoundsErrors[nx]) | nStatSlices;
          }
      }

    if (   (a3fCent[2] + m_fRotWidth/2.0>=m_nExt[2])
        || (a3fCent[2] - m_fRotWidth/2.0<=0.0))
      {
        return (1 << Crefln::ms_nErrorOffEdge2);
      }

    nMosPeakStartSearch = max(0,(int) (a3fCent[2] - m_fRotWidth/2.0));
    nMosPeakEndSearch   = min(m_nExt[2]-1,(int) (a3fCent[2] + m_fRotWidth/2.0));

    if (m_fRotWidth < 0.1)
      {
        // Allow a little slop here.  We add an extra 10% of an image width to either side.
        // This is not the same as adjusting the mosaicity (has units of degrees), because at this point, units of "images" are used.
        nMosPeakStartSearch = max(0,(int) (a3fCent[2] - m_fRotWidth/2.0 - 0.1));
        nMosPeakEndSearch = min(m_nExt[2]-1,(int) (a3fCent[2] + m_fRotWidth/2.0 + 0.1));
      }

    if (! (ms_nSpecialOptions & (1 << eSpecial_UseMosaicityBoundsOnStrong)))
      {
        // Use this if you want to look over the entire slice.
        // cout << "OFF UseMosaicityBoundsOnStrong\n";
        nPeakStartSearch = 0;
        nPeakEndSearch   = m_nExt[2]-1;
      } 
    else
      {
        // cout << "ON  UseMosaicityBoundsOnStrong\n";
        nPeakStartSearch = nMosPeakStartSearch;
        nPeakEndSearch = nMosPeakEndSearch;
      }
   
    // Initialize mask value.  Since comparisions to neighboring mask pixels
    // are made below, the initial values of the mask pixels are important.
    memset(ms_pnMask,0,sizeof(int)*nSize);

    //Initialize all variables that are computed in this loop.
    fSumXIntensity = 0.0;
    fSumYIntensity = 0.0;
    fSumZIntensity = 0.0;
    fSumSliceVarianceNumerator = 0.0;
    fSumSliceVarianceDenominator = 0.0;
    fSumNonnegativeIntensity = 0.0;
    fSumSliceIntensity = 0.0;
    fSumSliceNonnegativeIntensity = 0.0;
    fSumSliceBackgroundVariance = 0.0;
    fSumSliceIntensityVariance = 0.0;
    fSumBackground = 0.0;
    fSum0Dim = 0.0;
    fSum1Dim = 0.0;
    fSumIntensityVariance = 0.0;
    fSumBackgroundVariance = 0.0;

    // Calculate background for entire shoebox.
    
    // We don't trust all of the background slices.  We will only use
    // the slices with the lowest deviations.
    // Nevertheless, the bit-masks will still 
    //    use per-slice background information,        
    //    since we will not be recalculating these values.
    
    // Variables for background;
    float       a3fBackground[3];
    float       fAverage;
    float       fDeviation;
    float       fMaxValue;
    float       fOverallMaxValue = 0.0;
    float       a2fPeakCentroid[2];
    float       a2fCompeteCentroid[2];
    int         a2nCompeteCentroid[2];
    const int   nBackgroundSlots = 2;
    float       afMinDeviations[nBackgroundSlots];
    float       afMinAverage[nBackgroundSlots];
    float       afMinBackgrounds[nBackgroundSlots][3];
    bool        abUsed[nBackgroundSlots];       
    bool        bSliceWeak;
    bool        bUsedAnotherEllipsoid;
    bool        bStrongFoundAtCentroid;

    // These variables used in the loop below.
    float a2x2fEllipsoidA[2][2];           // Ellipsoid A for each pass.
    float a2fEllipsoidb[2];                // Ellipsoid b for each pass
    float fEllipsoidc;                     // Ellipsoid c for each pass.
    float a2fCenter[2];                    // Center of ellipsoid for each pass.
    float fAverageDim;                     // Average Min/Max bounding dimensions of ellipsoid for each pass.
    float a3fTotalCentroid[3];             // Used for spot chasing in the background determiniation loop.
    float a4fCumulTotalCentroid[4];        // Used for spot chasing in the background determiniation loop.
    int   a2nCentroidSliceCent[2];         // The ellipsoid saved must have a reference point.
    int   a2nOrigSliceCent[2];             // The ellipsoid returned from vGetShiftedCentroid() saved for a check.

    
    // We are not using any of the slots at first.
    for (ny = 0; ny < nBackgroundSlots; ny++)
      abUsed[ny] = FALSE;
    for (ny = 0; ny < 4; ny++) 
      a4fCumulTotalCentroid[ny] = 0;

    // If we have a user supplied ellipsoid, 
    //    use that ellipsoid as an initial 'plunk' down spot.
    // The centroid from this first pass will 
    //    then be used as the starting point if we are spot chasing.

    if ((m_a2x2fEllipsoidAIn[0][0] > 0.0) && (m_bSpotChase))
      {
        nCalcEllipsoidBitmask(0, a3nCent, ms_nEllipsoidMask,
                              m_a2x2fEllipsoidAIn,
                              m_a2fEllipsoidbIn,
                              m_fEllipsoidcIn, a2x2nSearch);
      }
    
    // Initialize the centroid mechanism.

    vGetShiftedCentroid(a3fTemp1,0.0,-1,true,a2nSliceCent);

    // Use the backgrounds with the lowest deviation.
    
    for (nx = 0; nx< m_nExt[2];nx++) 
      {
        nStatSlice = nBackgroundAvgSD(nx, a3fBackground, &fAverage,
                                      &fDeviation, &fMaxValue,
                                      &a2fPeakCentroid[0],
                                      &a3fTotalCentroid[0]);
        vCopyVec3D(a3fBackground,afSliceBackground + nx*3);
        afSliceIntensity[nx] = a3fTotalCentroid[2];
        if (    (m_bSpotChase) 
             && (nx >= nMosPeakStartSearch) 
             && (nx <= nMosPeakEndSearch) 
             && (a3fTotalCentroid[2]!=0.0) )
          {
            a4fCumulTotalCentroid[0] += a3fTotalCentroid[0];
            a4fCumulTotalCentroid[1] += a3fTotalCentroid[1];
            a4fCumulTotalCentroid[2] += nx*a3fTotalCentroid[2];
            a4fCumulTotalCentroid[3] += a3fTotalCentroid[2];

            // Initialize the spot chase mechanism.
            a3fTemp1[0] = a3fTotalCentroid[0] / a3fTotalCentroid[2];
            a3fTemp1[1] = a3fTotalCentroid[1] / a3fTotalCentroid[2];
            vGetShiftedCentroid(a3fTemp1, a3fTotalCentroid[2], nx,
                                TRUE, a2nSliceCent);
          }
        
        if (fMaxValue > fOverallMaxValue)
          {
            m_a3fPeakMaxValue[0]   = a2fPeakCentroid[0];
            m_a3fPeakMaxValue[1]   = a2fPeakCentroid[1];
            m_a3fPeakMaxValue[2]   = (float) nx;
            fOverallMaxValue       = fMaxValue;
          }
        
        if (nStatSlice)
          return nStatSlice;
        // See if the deviation of this slice can be used.
        if (nx < nBackgroundSlots)
            ny = nx;
        else 
          {                
            for (ny=0; ny<nBackgroundSlots; ny++)
              {
                if (afMinDeviations[ny] > fDeviation)
                    break;
              }
          }
        if (ny < nBackgroundSlots)
          {
            vCopyVec3D(&a3fBackground[0],&afMinBackgrounds[ny][0]);
            afMinAverage[ny] = fAverage;
            afMinDeviations[ny] = fDeviation;
            abUsed[ny] = TRUE;
          }
      }
    
    a3nCent[2] = nFindBestCentroidSlice(a3fCent, afSliceIntensity,
                                        nMosPeakStartSearch, nMosPeakEndSearch);
    if (a3nCent[2]<nMosPeakStartSearch)
      nMosPeakStartSearch = a3nCent[2];
    if (a3nCent[2]>nMosPeakEndSearch)
      nMosPeakEndSearch = a3nCent[2];
    
    if ((a4fCumulTotalCentroid[3]) && (m_bSpotChase))
      {
        for (ny = 0,nz = 0; ny < m_nExt[2]; ny++)
          {
            if (afSliceIntensity[ny] > afSliceIntensity[nz])
              nz = ny;
          }
        nx = (int) (a4fCumulTotalCentroid[2]/a4fCumulTotalCentroid[3] + 0.5);
        if ((nx == a3nCent[2]) && (nz == a3nCent[2]))
          {
            // Adjust the spot centroid since we are fairly certain 
            //  it's correct.
            
            a3fCent[0] = a4fCumulTotalCentroid[0]/a4fCumulTotalCentroid[3];
            m_a3fAbsCentroid[0] = a3fCent[0] + m_nOrigOffset[0];
            a3fCent[1] = a4fCumulTotalCentroid[1]/a4fCumulTotalCentroid[3];
            m_a3fAbsCentroid[1] = a3fCent[1] + m_nOrigOffset[1];
          }
      }

    nSlice      = a3nCent[2];
    m_nCentroid = a3nCent[2];

    for (nx=0;nx<m_nExt[2];nx++) 
      {
        afSliceIntensity[nx] = 0.0;
        afSliceSigma[nx]     = 0.0;
      }
    
    // Calculate the average values for the background.

    vZeroMat(3,1,&a3fBackground[0]);
    fDeviation = 0.0;
    fAverage   = 0.0;
    for (ny = 0, nx = 0; ny < nBackgroundSlots; ny++)
      {
        if (abUsed[ny])
          {
            nx++;
            vAddVec3DVec3D(&a3fBackground[0], &afMinBackgrounds[ny][0],
                           &a3fBackground[0]);
            fDeviation += afMinDeviations[ny] * afMinDeviations[ny];
            fAverage   += afMinAverage[ny];
          }
      }
    vMulVec3DScalar(&a3fBackground[0], 1.0/(float)nx, m_a3fAvgBackground);
    m_fAvgBackgroundAvg   = fAverage /(float)nx;
    m_fAvgBackgroundSigma = sqrt((double)(fDeviation/(float)nx));
    if (!m_bUsePerSliceBackground) 
      {
        // Copy the average background to each of the background slots
        for (nx=0;nx< m_nExt[2];nx++)
          vCopyVec3D(m_a3fAvgBackground,afSliceBackground + nx*3);
      }
    bStrongFoundAtCentroid = false;

    // Continue looping over the 2D slices in the shoebox.
    // Depending upon whether this is the bootstrap calculation or not, we will
    // calculate spot masks.

    while (bContinueLoop)
      {
        fSumSliceIntensity = 0.0;
        fSumSliceNonnegativeIntensity = 0.0;
        fSumSliceBackgroundVariance = 0.0;
        fSumSliceIntensityVariance = 0.0;
        fSumSliceXIntensity = 0.0;
        fSumSliceYIntensity = 0.0;
        fSumSlicePixelsOnEdge0 = 0.0;
        fSumSlicePixelsOnEdge1 = 0.0;        

        nStatSlice = 0;
        nSliceOffset = nSlice * nSliceSize;

        // Get the centroid shift.

        vGetShiftedCentroid(a3fCent,0.0,nSlice,false,a2nSliceCent);

        // Save this for checking later.

        a2nOrigSliceCent[0] = a2nSliceCent[0];
        a2nOrigSliceCent[1] = a2nSliceCent[1];

        // Calculate all of the pixels that are in the contiguous pixel
        // region that has the centroid as the seed.
        // we use the stack algorithm to do this.

        if ((a2nSliceCent[0]<0) || (a2nSliceCent[0]>=m_nExt[0]))
          nStatSlice |= (1 << Crefln::ms_nErrorOffEdge0);
        if ((a2nSliceCent[1]<0) || (a2nSliceCent[1]>=m_nExt[1])) 
          nStatSlice |= (1 << Crefln::ms_nErrorOffEdge1);

        // For use with printing, but otherwise ignored.
        if (!nStatSlice)
          ms_pnMask[nSliceOffset + a2nSliceCent[0] 
                    + a2nSliceCent[1]*m_nExt[0]] |= ms_nCentroidPixelMask;

        // Main mask loop:

        if (!nStatSlice)
          {
            /*
             Make two passes.  We are looking for the smallest ellipsoid that
             fits the data.  The first ellipsoid is derived by shifting 
             the proposed ellipsoid over by the per-image shift, 
             and pasting it on the slice.  Since this sometimes does not work
             (or is too large for it's own good), we always compute a secondary
              ellipsoid in competition with the first.

            The secondary ellipsoid is found as follows:
                1)  We take the centroid of all the high pixels in the slice.  
                2)  If this centroid is not high itself, we fail. 
                3)  Otherwise, we find a new contiguous high pixel region 
                       (starting at this centroid).
                4)  From this, we build the secondary ellipse.

            The winning ellipse must have the following properties:

                1)  No on-edge ms_nContiguousHighMask pixels 
                2)  If the ellipsoid is the secondary one, then it cannot 
                        shift over by more than 0.5 
                        of the average diameter of the primary ellipse.  
                        This will keep us from picking up
                        neighboring spots when computing the secondary ellipse.
                2)  If both elipsoids work, then the ellipsoid must have 
                        a larger number of high pixels or a smaller area 
                        with the same number of high pixels.
            */
            
            vZeroMat(2,1,&a2x2fEllipsoidA[0][0]);
            vZeroMat(2,1,&a2fEllipsoidb[0]);
            fEllipsoidc  = 0.0;
            a2fCenter[0] = 0.0;
            a2fCenter[1] = 0.0;
            fAverageDim  = 0;
            a2x2nSearch[0][0] = a2nSliceCent[0];
            a2x2nSearch[0][1] = a2nSliceCent[0];
            a2x2nSearch[1][0] = a2nSliceCent[1];
            a2x2nSearch[1][1] = a2nSliceCent[1];
                     
            bSliceWeak = false;
            bUsedAnotherEllipsoid = false;

            nx = nCalcContiguous(nSlice,ms_nHighMask,
                                 ms_nContiguousHighMask,a2nSliceCent);
            if (nx)
              {
                bUsedAnotherEllipsoid = true;
                
                // Compute the ellispoid to use.
                if (   (nSlice == m_nCentroid)
                    && (m_a2x2fEllipsoidAIn[0][0]!=0.0))
                  {
                    // Use the user supplied ellipsoid for weak spots.
                    a2x2fEllipsoidA[0][0] = m_a2x2fEllipsoidAIn[0][0];
                    a2x2fEllipsoidA[0][1] = m_a2x2fEllipsoidAIn[0][1];
                    a2x2fEllipsoidA[1][0] = m_a2x2fEllipsoidAIn[1][0];
                    a2x2fEllipsoidA[1][1] = m_a2x2fEllipsoidAIn[1][1];
                    a2fEllipsoidb[0]      = m_a2fEllipsoidbIn[0];
                    a2fEllipsoidb[1]      = m_a2fEllipsoidbIn[1];
                    fEllipsoidc           = m_fEllipsoidcIn;
                    afSliceStatus[nSlice] = ((int) afSliceStatus[nSlice])
                                             | ms_nPerSliceIntegrateWeak;
                } 
                else if ((m_a2x2fEllipsoidA[0][0]))
                  {
                    // Use the ellipsoid from the centroid slice if it's available.
                    a2x2fEllipsoidA[0][0] = m_a2x2fEllipsoidA[0][0];
                    a2x2fEllipsoidA[0][1] = m_a2x2fEllipsoidA[0][1];
                    a2x2fEllipsoidA[1][0] = m_a2x2fEllipsoidA[1][0];
                    a2x2fEllipsoidA[1][1] = m_a2x2fEllipsoidA[1][1];
                    a2fEllipsoidb[0]      = m_a2fEllipsoidb[0];
                    a2fEllipsoidb[1]      = m_a2fEllipsoidb[1];
                    a2nSliceCent[0]       = a2nCentroidSliceCent[0];
                    a2nSliceCent[1]       = a2nCentroidSliceCent[1];
                    fEllipsoidc           = m_fEllipsoidc;
                    afSliceStatus[nSlice] = ((int) afSliceStatus[nSlice])
                                              | ms_nPerSliceIntegrateWeak;
                  } 
                else 
                  {
                    // Calculate a generic ellipsoid.from the bits available.
                    nStatSlice |= nCalcEllipsoid(nSlice, a2nSliceCent,
                                                 a2x2fEllipsoidA,
                                                 a2fEllipsoidb,
                                                 fEllipsoidc);
                  }
                
                if (   (nSlice < nMosPeakStartSearch)
                    || (nSlice>nMosPeakEndSearch))
                  {
                    if (m_nIntegrateWeakFlag <= 1)
                      nStatSlice |= 1 << Crefln::ms_nErrorIntProblem;
                  } 
                else if (m_nIntegrateWeakFlag == 0)
                  {
                    if (bStrongFoundAtCentroid)
                      nStatSlice |= 1 << Crefln::ms_nErrorIntProblem;
                  }
                
                if (m_nIntegrateWeakFlag <= 1)
                  {
                    // If we were too weak, then restrict further padding somewhat.
                    if (m_nCentroid == nSlice)
                      {
                        nStatSlice |= 1 << Crefln::ms_nErrorTooFewPix;
                        bPadOnUpperEnd = FALSE;
                        bPadOnLowerEnd = FALSE;
                      } 
                    else if (nSlice < m_nCentroid)
                      bPadOnLowerEnd = FALSE;
                    else
                      bPadOnUpperEnd = FALSE;
                  }
              } 
            else 
              {
                // Calculate a generic ellipsoid.from the bits available.

                nStatSlice |= nCalcEllipsoid(nSlice, a2nSliceCent,
                                             a2x2fEllipsoidA,
                                             a2fEllipsoidb,
                                             fEllipsoidc);
                bStrongFoundAtCentroid = bStrongFoundAtCentroid 
                                         || (nSlice== m_nCentroid);
              }
            // Set the bits in the ellipsoid.
            nCalcEllipsoidBitmask(nSlice, a2nSliceCent, ms_nEllipsoidMask,
                                  a2x2fEllipsoidA,a2fEllipsoidb,fEllipsoidc,
                                  a2x2nSearch,
                                  /* afSliceBackground + nSlice*3 */ NULL);

            if (!(nStatSlice & nWarningMask))
              {
                nStatSlice |= nCheckFuzzySpotShape(nSlice);

                // Save the shift  (This code has been disabled,
                //  but would come in handy if we need to look at the centroid).
                // nInvMatND_svd(2,&a2x2fEllipsoidA[0][0],&a2x2fTempMat[0][0]);
                // vMulMatNDVecND(2,&a2x2fTempMat[0][0],&a2fEllipsoidb[0],&a2fTemp1[0]);
                // a2fCenter[0] =  (-0.5*a2fTemp1[0]) + a2nSliceCent[0];
                // a2fCenter[1] =  (-0.5*a2fTemp1[1]) + a2nSliceCent[1];

                // Fill in the selections based on the "nBestEllipsoid"            
                if (!m_a2x2fEllipsoidA[0][0])
                  {
                    // Fill in the centroid ellipsoid.
                    m_a2x2nEllipsoidRange[0][0] = a2x2nSearch[0][0];
                    m_a2x2nEllipsoidRange[1][0] = a2x2nSearch[1][0];
                    m_a2x2nEllipsoidRange[0][1] = a2x2nSearch[0][1];
                    m_a2x2nEllipsoidRange[1][1] = a2x2nSearch[1][1];
                    vCopyVecND(4,&a2x2fEllipsoidA[0][0],&m_a2x2fEllipsoidA[0][0]);
                    vCopyVecND(2,&a2fEllipsoidb[0],&m_a2fEllipsoidb[0]);
                    m_fEllipsoidc = fEllipsoidc;
                    a2nCentroidSliceCent[0] = a2nSliceCent[0];
                    a2nCentroidSliceCent[1] = a2nSliceCent[1];
                    
                  } 
                else
                  {
                    m_a2x2nEllipsoidRange[0][0] = min(a2x2nSearch[0][0],m_a2x2nEllipsoidRange[0][0]);
                    m_a2x2nEllipsoidRange[1][0] = min(a2x2nSearch[1][0],m_a2x2nEllipsoidRange[1][0]);
                    m_a2x2nEllipsoidRange[0][1] = max(a2x2nSearch[0][1],m_a2x2nEllipsoidRange[0][1]);
                    m_a2x2nEllipsoidRange[1][1] = max(a2x2nSearch[1][1],m_a2x2nEllipsoidRange[1][1]);
                  }

                // Fill in per-slice ellipsoid and convert to a standard
                //   reference point of 0.0 (instead of the a2nSliceCent[] 
                //   reference which is used.

                vCopyVecND(4,&a2x2fEllipsoidA[0][0],&m_apfScratchFloat[m_nScratchFloatEllipsoidA][nSlice*4]);
                vCopyVecND(2,&a2fEllipsoidb[0],&m_apfScratchFloat[m_nScratchFloatEllipsoidb][nSlice*2]);
                m_apfScratchFloat[m_nScratchFloatEllipsoidc][nSlice] = fEllipsoidc;
                vConvertEllipsoidRelative(&m_apfScratchFloat[m_nScratchFloatEllipsoidA][nSlice*4],&m_apfScratchFloat[m_nScratchFloatEllipsoidb][nSlice*2],&m_apfScratchFloat[m_nScratchFloatEllipsoidc][nSlice],a2nSliceCent[0],a2nSliceCent[1],0,0);

                // Save the average dimension.

                fAverageDim = 0.5 * (a2x2nSearch[0][1] + a2x2nSearch[1][1] 
                                   - a2x2nSearch[0][0] - a2x2nSearch[1][0])+1;
              } 
                            
            // Must at least have the centroid in the ellipsoid.

            if ((nSlices == 0) && (!(ms_pnMask[nSliceOffset + a2nOrigSliceCent[0] + a2nOrigSliceCent[1]*m_nExt[0]] & ms_nEllipsoidMask))) 
              nStatSlice |= 1 << Crefln::ms_nErrorIntProblem;

            // Check for competing slice.
            if ((ms_nSpecialOptions & (1 << eSpecial_PhiHKLLimits)))
              {
                vGetCompetingSpot(nSlice,&a2fCompeteCentroid[0],&a2fCompeteCentroid[1]);
                if (   (a2fCompeteCentroid[0] != 0.0)
                    || (a2fCompeteCentroid[1] != 0.0))
                  {
                    a2nCompeteCentroid[0] = (int) (a2fCompeteCentroid[0] - m_nOrigOffset[0]);
                    a2nCompeteCentroid[1] = (int) (a2fCompeteCentroid[1] - m_nOrigOffset[1]);
                    // Do we have a competing spot?
                    
                    if (   (a2nCompeteCentroid[0]>=0) 
                        && (a2nCompeteCentroid[0]<m_nExt[0])
                        && (a2nCompeteCentroid[1]>=0) 
                        && (a2nCompeteCentroid[1]<m_nExt[1])
                        && (ms_pnMask[nSliceOffset + a2nCompeteCentroid[0]
                                       + a2nCompeteCentroid[1]*m_nExt[0]]
                            & ms_nEllipsoidMask))
                      {
                        // Is the shoebox we are integrating supposed to 
                        // have this intensity as well?
                        // If so, this is a full fledged overlap, 
                        // and we might as well return now.

                        if (   (nSlice>=nMosPeakStartSearch) 
                            && (nSlice<=nMosPeakEndSearch))
                          {
                            return 1 << Crefln::ms_nErrorIntProblem;
                          } 
                        else 
                          {
                            nStatSlice |= 1 << Crefln::ms_nErrorIntProblem;
                          }
                      }
                  }
              }

            // Save copies of variables that might be flagged 
            // as errors in the loop below.

            double _fSumXIntensity                = fSumXIntensity;
            double _fSumYIntensity                = fSumYIntensity;
            double _fSumZIntensity                = fSumZIntensity;
            double _fSumNonnegativeIntensity      = fSumNonnegativeIntensity;
            double _fSumSliceIntensity            = fSumSliceIntensity;
            double _fSumSliceNonnegativeIntensity = fSumSliceNonnegativeIntensity;
            double _fSumSliceBackgroundVariance   = fSumSliceBackgroundVariance;
            double _fSumSliceIntensityVariance    = fSumSliceIntensityVariance;
            double _fSumBackground                = fSumBackground;
            double _fSumBackgroundVariance        = fSumBackgroundVariance;
            double _fSumIntensityVariance         = fSumIntensityVariance;
            
            
            // Final processing loop.

            for (n1=0,nx=0;   (n1<m_nExt[1]) 
                           && (!(nStatSlice & nWarningMask)); n1++)
              {
                for (n0=0;(n0<m_nExt[0]) 
                           && (!(nStatSlice & nWarningMask));n0++,nx++)
                  {
                    // if this pixel is in the mask..               
                    
                    if (ms_pnMask[nSliceOffset + nx] & ms_nEllipsoidMask)
                      {
                        // We have a mask pixel.  

                        // Check for saturated or nonunf pixels.
                        if (ms_pnMask[nx + nSliceOffset] & ms_nSatMask)
                          {
                            // If we are bootstrapping, then don't flag as bad 
                            // Otherwise, this reflection should get rejected, no matter what
                            // slice we are on.
                            return 1 << Crefln::ms_nErrorTooDark;
                          } 
                        else if (ms_pnMask[nx + nSliceOffset] & ms_nBadMask)
                          {
                            // Always return an error for this as well.
                            return 1 << Crefln::ms_nErrorNonunfB;
                          } 
                        else
                          {
                            if (m_bFitToPlanar)
                              {
                                f1 = (afSliceBackground + nSlice*3)[0]*(float)n0/(float)m_nExt[0] 
                                    + (afSliceBackground + nSlice*3)[1]*(float)n1/(float)m_nExt[1] 
                                    + (afSliceBackground + nSlice*3)[2];
                                
                              }
                            else
                              f1 = m_a3fAvgBackground[2];
                            f0 = m_pfData[nx + nSliceOffset] - f1;
                            
                            // Also, make sure that the pixel is not on the shoebox edge.
                            if ((n1==0) || (n1==m_nExt[1]-1))
                              {
                                fSumSlicePixelsOnEdge0 += f0;
                              }
                            if ((n0==0) || (n0==m_nExt[0]-1))
                              {
                                fSumSlicePixelsOnEdge1 += f0;
                              }
                            
                            f2 = max(f0,0.0);
                            
                            // Add the pixel in.  Count its statistics.
                            
                            fSumSliceXIntensity           += (n0 - a3nCent[0]) * f2;
                            fSumSliceYIntensity           += (n1 - a3nCent[1]) * f2;
                            fSumXIntensity                += (n0 - a3nCent[0]) * f2;
                            fSumYIntensity                += (n1 - a3nCent[1]) * f2;
                            fSumZIntensity                += (nSlice - m_nCentroid) * f2;
                            fSumNonnegativeIntensity      += f2;
                            fSumSliceIntensity            += f0;
                            fSumSliceNonnegativeIntensity += f2;
                            fSumBackground                += f1;
                            fSumBackgroundVariance        += m_fAvgBackgroundSigma * m_fAvgBackgroundSigma;
                            fSumSliceBackgroundVariance   += m_fAvgBackgroundSigma * m_fAvgBackgroundSigma;
                            fPeakMaxValue                 = max(fPeakMaxValue, f0);
                            fPeakMinValue                 = min(fPeakMinValue, f0);
                          }
                      }
                  }
              }
            fSumSliceIntensityVariance = max(0.0,fSumSliceIntensity)* m_fDetectorGain;
            fSumIntensityVariance += fSumSliceIntensityVariance;

            if (! (nStatSlice & nWarningMask))
              {
                // Set the shifted centroid if it's available.
                a3fTemp1[0] = a3nCent[0] + fSumSliceXIntensity
                                           / (fSumSliceNonnegativeIntensity);
                a3fTemp1[1] = a3nCent[1] + fSumSliceYIntensity 
                                           / (fSumSliceNonnegativeIntensity);
                a3fTemp1[2] = 0.0;
                vGetShiftedCentroid(a3fTemp1,fSumSliceIntensity,nSlice,true,a2nSliceCent);
              }

            if (!(nStatSlice & nWarningMask))
              {
                afSliceIntensity[nSlice] = fSumSliceIntensity;
                afSliceSigma[nSlice] = sqrt((double)(fSumSliceBackgroundVariance + fSumSliceIntensityVariance));
                afSliceSigmaBackground[nSlice] = sqrt((double)fSumSliceBackgroundVariance);
              }
            
            // If this is a non-centroid slice, and the centroid slice had I/sig >= a2fMinSliceReject[*][0] then
            // reject this slice if the (I/sig[centroid])/(I/sig[slice])> a2fMinSliceReject[*][1]
            
            if (  (nSlices)
                && 
                  (ms_nSpecialOptions & (1 << eSpecial_DoNotIntegrateGarbage))
                &&
                  (!(nStatSlice & nWarningMask)))
              {
                
                const float a2fMinSliceReject[][2] = {{3.0,1.0},{5.0,2.0}};
                const int nSliceRejectCount = 2;                    
                
                // cout << "ON  DoNotIntegrateGarbage\n";
                
                for (nx=0;nx<nSliceRejectCount;nx++)
                  {
                    f0 = fSumSliceIntensity/max(1.0,sqrt((double)(fSumSliceBackgroundVariance + fSumSliceIntensityVariance)));
                    f1 = afSliceIntensity[m_nCentroid]/max(1.0,afSliceSigma[m_nCentroid]);                    
                    if ((f1>a2fMinSliceReject[nx][0]) && (f0<a2fMinSliceReject[nx][1]))
                      nStatSlice |= 1 << Crefln::ms_nErrorIntProblem;
                  }                           
              }

            if ((!nSlices) || (!bUsedAnotherEllipsoid))
              {
                if (!nSlices)
                  {
                    m_nIntensePeakStart = nSlice;
                    m_nIntensePeakEnd = nSlice;
                  } 
                else if (!nStatSlice)
                  {
                    m_nIntensePeakStart = min(nSlice,m_nIntensePeakStart);
                    m_nIntensePeakEnd = max(nSlice,m_nIntensePeakEnd);
                  }
              } 


            // Fill in flags as to whether this spot is in mosaicity bounds, or has the correct centroid.
            if ((nSlice>=nMosPeakStartSearch) && (nSlice<=nMosPeakEndSearch))
              afSliceStatus[nSlice] = ((int) afSliceStatus[nSlice]) | ms_nPerSliceInMosaicity;
            if (nSlice == m_nCentroid)
              afSliceStatus[nSlice] = ((int) afSliceStatus[nSlice]) | ms_nPerSliceCentroid;
            
            afSliceStat[nSlice] = nStatSlice;
                       
            if (nStatSlice & nWarningMask)
              {
                // Undo affects of this integration.
                fSumXIntensity = _fSumXIntensity;
                fSumYIntensity = _fSumYIntensity;
                fSumZIntensity =_fSumZIntensity;
                fSumNonnegativeIntensity = _fSumNonnegativeIntensity;
                fSumSliceIntensity = _fSumSliceIntensity;
                fSumSliceNonnegativeIntensity = _fSumSliceNonnegativeIntensity;
                fSumSliceBackgroundVariance = _fSumSliceBackgroundVariance;
                fSumSliceIntensityVariance = _fSumSliceIntensityVariance;
                fSumBackground = _fSumBackground;
                fSumBackgroundVariance = _fSumBackgroundVariance;
                fSumIntensityVariance = _fSumIntensityVariance;
              }           

            // Depending upon how many pixels were 'on edge', we throw errors.
            if (   (fSumSlicePixelsOnEdge0 !=0.0)
                || (fSumSlicePixelsOnEdge1 !=0.0))
              {
                if (   (fSumSliceIntensity<0.0)
                    || (fSumSliceIntensity*m_fOnEdgeMaxPercent < fSumSlicePixelsOnEdge0+fSumSlicePixelsOnEdge1))
                  {
                    if (fSumSlicePixelsOnEdge0!=0.0)
                      nStatSlice |= 1 << Crefln::ms_nErrorOnEdge0;
                    if (fSumSlicePixelsOnEdge1!=0.0)
                      nStatSlice |= 1 << Crefln::ms_nErrorOnEdge1;
                  }
              } 
          }
        
        if (nStatSlice & nWarningMask)
          {
            bSliceAdded = FALSE;
            if (nSlices==0)
              nStatSlices |= nStatSlice;
          } 
        else
          {
            bSliceAdded = TRUE;
            nStatSlices |= nStatSlice;
            
            if (m_nPeakStart == -1)
              {
                m_nPeakStart = nSlice;
                m_nPeakEnd = nSlice;
              } 
            else
              {
                m_nPeakEnd = max(nSlice,m_nPeakEnd);
                m_nPeakStart = min(nSlice,m_nPeakStart);                
              }

            // Copy averaging propreties over.
            fSum0Dim += (a2x2nSearch[0][1] - a2x2nSearch[0][0]+1);
            fSum1Dim += (a2x2nSearch[1][1] - a2x2nSearch[1][0]+1);
            fSumSliceVarianceNumerator   += fSumSliceIntensity*fSumSliceIntensity;
            fSumSliceVarianceDenominator += fSumSliceIntensity;
          }

        // Compute the next slice that needs to get considered.    
        
        if (!bSliceAdded)
          {
            // If we failed to add this slice, then we see what we can try next.
            if (nSlices == 0)
              {
                // Don't attempt looking above or below the centroid slice:  we have failed
                nPositiveSlices = -1;
                nNegativeSlices = -1;
                bContinueLoop = FALSE;
              } 
            else if (nPositiveSlices >= 0) 
              {
                // Don't attempt looking any more above the centroid slice: Look below the centroid slice.
                nPositiveSlices = -1; 
                bContinueLoop = TRUE;
              } 
            else if (nNegativeSlices >= 0) 
              {
                // Don't attempt looking above or below the centroid slice:  we are done.
                nNegativeSlices = -1;
                bContinueLoop = FALSE;
              } 
            else
              // Should not get here.
              bContinueLoop = FALSE;
          } 
        else
          {
            // Just keep on doing what we were doing 
            bContinueLoop = TRUE;
            nSlices++;
          }

        if (   (nPositiveSlices >= 0)
             && 
            (a3nCent[2] + nPositiveSlices + 1<=nPeakEndSearch)
             &&
            ((a3nCent[2] + nPositiveSlices + 1 <= nMosPeakEndSearch) || (bPadOnUpperEnd)))
          {
            // Look at the next slice above the centroid slice. (maybe it is several slices above the centroid)
            nLastSlice = a3nCent[2] + nPositiveSlices;
            nSlice = a3nCent[2] + nPositiveSlices + 1;
            nPositiveSlices++;
          } 
        else if (  (nNegativeSlices >= 0)
                 && 
                   (a3nCent[2] - nNegativeSlices - 1>=nPeakStartSearch)
                 &&
                   ((a3nCent[2] - nNegativeSlices - 1 >= nMosPeakEndSearch)
                     || (bPadOnLowerEnd)))
          {
            // Look at the next slice below the centroid slice. 
            //  (maybe it is several slices below the centroid)
           
            nPositiveSlices = -1;
            nSlice = a3nCent[2] - nNegativeSlices - 1;
            nLastSlice = a3nCent[2] - nNegativeSlices;
            nNegativeSlices++;
          } 
        else
          {
            nPositiveSlices = -1;
            nNegativeSlices = -1;
            bContinueLoop = FALSE;
          }
      }

    volatile int nThad = 0;
    if (nThad==1)
      {
        nPrintMask(true, ms_nBorderMask);
        
        nPrintMask(true, ms_nMinPeakPixels);
        nPrintMask(true, ms_nHighMask);
        nPrintMask(true, ms_nContiguousHighMask);
        nPrintMask(true, ms_nEllipsoidMask);
      }

    // Calculate the average centroid.
    // The centroids sums used the shifted centroid at all places as a reference.
    m_a3fCentroid[0] = a3nCent[0] + 0.5;
    m_a3fCentroid[1] = a3nCent[1] + 0.5;
    m_a3fCentroid[2] = 0.5 + m_nCentroid;
    if (fSumNonnegativeIntensity > 0.0)
      {
        m_a3fCentroid[0] += fSumXIntensity / (fSumNonnegativeIntensity);
        m_a3fCentroid[1] += fSumYIntensity / (fSumNonnegativeIntensity);
        m_a3fCentroid[2] += fSumZIntensity / (fSumNonnegativeIntensity);
      }


    // We might want to use a profile fitted spot.
    if ((m_bProfileFitWeak) && (m_poProfit) && (m_nPeakStart>=0))
      {
        for (nSlice = m_nPeakStart; nSlice <= m_nPeakEnd; nSlice++)
          {
            double a2fCentroid[2];
            a2fCentroid[0] = m_a3fCentroid[0];
            a2fCentroid[1] = m_a3fCentroid[1];
            if (m_poProfit->nFit(nSlice,*this,
                                 a2fCentroid,afSliceIntensityProfit[nSlice],
                                 afSliceSigmaProfit[nSlice]))
              {
                afSliceIntensityProfit[nSlice] = 0.0;
                afSliceSigmaProfit[nSlice] = 0.0;
              } 
          }
      }
  
    // Peak intensity
    m_fPeakIntensity = 0.0;
    m_fPeakIntensitySigma = 0.0;
    m_fPeakBackgroundSigma = 0.0;
    if (m_nPeakStart>=0)
      {
        for (nSlice = m_nPeakStart; nSlice <= m_nPeakEnd; nSlice++)
          {
            if (afSliceIntensity[nSlice]!=0.0)
              {
                if (((int) afSliceStatus[nSlice]) & ms_nPerSliceIntegrateWeak)
                  {
                    m_fPeakIntensity +=afSliceIntensityProfit[nSlice];
                    m_fPeakIntensitySigma += afSliceSigmaProfit[nSlice]*afSliceSigmaProfit[nSlice];
                  } 
                else if ((afSliceIntensityProfit[nSlice]==0.0) || (afSliceSigmaProfit[nSlice]>afSliceSigma[nSlice]))
                  {
                    m_fPeakIntensity +=afSliceIntensity[nSlice];
                    m_fPeakIntensitySigma += afSliceSigma[nSlice]*afSliceSigma[nSlice];
                  } 
                else
                  {
                    m_fPeakIntensity +=afSliceIntensityProfit[nSlice];
                    m_fPeakIntensitySigma += afSliceSigmaProfit[nSlice]*afSliceSigmaProfit[nSlice];
                  }
                m_fPeakBackgroundSigma += afSliceSigmaBackground[nSlice]*afSliceSigmaBackground[nSlice];
              }
          }
      }
    
    m_fPeakIntensitySigma = sqrt(m_fPeakIntensitySigma);
    m_fPeakBackgroundSigma = sqrt(m_fPeakBackgroundSigma);

    if (   (!nStatSlices) 
        && (m_poProfit) && (m_nPeakStart>=0) 
        && (m_fPeakIntensity>=m_fMinProfitPeakIntensity))
      {
        m_poProfit->nAdd(m_nPeakStart,m_nPeakEnd,*this);
      }

    // Check if out-of-bounds of shoebox

    for (nx=0; nx < 3; nx++)
      {
        if (   (0 > m_a3fCentroid[nx]) 
            || (m_a3fCentroid[nx] >= m_nExt[nx]) )
          {
            cout << "OCentroid: " << m_a3fCentroid[0] << ", "
                << m_a3fCentroid[1] << ", " << m_a3fCentroid[2] << endl;
            
            if (nx < 2)
                m_a3fCentroid[nx] = a3nCent[nx] + 0.5;
            else
                m_a3fCentroid[nx] = m_nCentroid  + 0.5;
          }
        if (    (0 > a3nCent[nx]) 
            || (a3nCent[nx] >= m_nExt[nx]) )
          {
            cout << "NOCentroid: " << a3nCent[0] << ", "
                << a3nCent[1] << ", " << a3nCent[2] << endl;
          }
      }
      
    // Copy other parameters over.
    m_fPeakMinValue = fPeakMinValue;
    m_fPeakMaxValue = fPeakMaxValue;

    // Rotation sigma.
    if (fSumSliceVarianceDenominator==0.0)
      m_fRotSigma = 1.0;
    else
      m_fRotSigma = sqrt((double)(fSumSliceVarianceNumerator/(fSumSliceVarianceDenominator*fSumSliceVarianceDenominator*12)));

    // Peak shape
    if (nSlices!=0)
      {
        m_a3fPeakSize[0] = fSum0Dim/nSlices;
        m_a3fPeakSize[1] = fSum1Dim/nSlices;
      } 
    else
      {
        m_a3fPeakSize[0] = 10;
        m_a3fPeakSize[1] = 10;
      }
    m_a3fPeakSize[2] = m_nPeakEnd - m_nPeakStart + 1;

    // Sharpness constant.

    m_fPeakSharpness = 2.0*sqrt(
        (m_a3fPeakMaxValue[0] + 0.5 - m_a3fCentroid[0])*(m_a3fPeakMaxValue[0] + 0.5 - m_a3fCentroid[0]) +
        (m_a3fPeakMaxValue[1] + 0.5 - m_a3fCentroid[1])*(m_a3fPeakMaxValue[1] + 0.5 - m_a3fCentroid[1])
        )/(sqrt(m_a3fPeakSize[0]*m_a3fPeakSize[0] + m_a3fPeakSize[1]*m_a3fPeakSize[1]));


    // Figure out what portion of the found shoebox should be used 
    //     to calculate the observed width.
    // This is very important, since otherwise, the mosaicity tends 
    //     to increase as an integration proceeds.
    // That is, the shoebox algorithm tends to find data where 
    //     perhaps there is none.

    if (    (m_fPeakIntensity>0.0)
         && (m_nIntensePeakStart!=-1)
         && (m_nIntensePeakEnd!=-1)
         && (ms_nSpecialOptions & (1 << eSpecial_TrimEdges)))
    {
      double fCumulativeRemoval;
      bool bLeftCandidate = true;
      bool bRightCandidate = true;

      //        cout << "ON  TrimEdges\n";
      fCumulativeRemoval = 0.0;
      while (   (m_nIntensePeakEnd != m_nIntensePeakStart) 
             && ((bLeftCandidate) || (bRightCandidate)))
        {
          bLeftCandidate = ((afSliceIntensity[m_nIntensePeakStart] + fCumulativeRemoval)/m_fPeakIntensity<1.0 - m_fFractionIntense);
          bRightCandidate = ((afSliceIntensity[m_nIntensePeakEnd] + fCumulativeRemoval)/m_fPeakIntensity<1.0 -  m_fFractionIntense);
          bLeftCandidate = bLeftCandidate || (afSliceIntensity[m_nIntensePeakStart]/afSliceSigma[m_nIntensePeakStart]<m_fMosaicityLowestSpotIoverSig);
          bRightCandidate = bRightCandidate || (afSliceIntensity[m_nIntensePeakEnd]/afSliceSigma[m_nIntensePeakEnd]<m_fMosaicityLowestSpotIoverSig);
          if ((m_nIntensePeakEnd - m_nIntensePeakStart==1) && (bLeftCandidate) && (bRightCandidate))
            {
              // Removing both sides removes the peak. 
              m_nIntensePeakStart = -1;
              m_nIntensePeakEnd = -1;
              bLeftCandidate = false;
              bRightCandidate = false;
            }
          if (bLeftCandidate)
            {
              fCumulativeRemoval+= afSliceIntensity[m_nIntensePeakStart];
              m_nIntensePeakStart++;
            }
          if (bRightCandidate)
            {
              fCumulativeRemoval += afSliceIntensity[m_nIntensePeakEnd];
              m_nIntensePeakEnd--;
            }
        }
    }

    return nStatSlices;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Peak area pixels are all high intensity contiguous area pixels, used to build the ellipsoid. Peak area border pixels are pixels that belong
// to the contiguous area, yet have at least one immediate neighbour that does not belong to that area.   
void C3Ddata::vGet2DPeakAreaPixelsCount(int& nPeakAreaPixels, 
                                        int& nPeakEllipsPixels,
                                        int& nPeakBorderPixels)
{
    nPeakAreaPixels   = 0;
    nPeakEllipsPixels = 0;
    nPeakBorderPixels = 0;
    
    int     nShoeBoxSize0 = m_nExt[0];
    int     nShoeBoxSize1 = m_nExt[1];

    for(int ii=0; ii < nShoeBoxSize1; ii++)
    {
        for(int jj=0; jj < nShoeBoxSize0; jj++)
        {
            if( ms_pnMask[ii * nShoeBoxSize0 + jj] & ms_nEllipsoidMask )
                nPeakEllipsPixels++;

            if( ms_pnMask[ii * nShoeBoxSize0 + jj] & ms_nContiguousHighMask )
            {
                nPeakAreaPixels++;
                
                if( ii == 0 || ii == nShoeBoxSize1 - 1 ||
                    jj == 0 || jj == nShoeBoxSize0 - 1 || // if these contiguous are pixels are shoe-box border pixels, they are automatically peak area border pixels

                   !( ms_pnMask[(ii - 1) * nShoeBoxSize0 + jj] & ms_nContiguousHighMask ) ||
                   !( ms_pnMask[(ii + 1) * nShoeBoxSize0 + jj] & ms_nContiguousHighMask ) ||
                   !( ms_pnMask[ii * nShoeBoxSize0 + jj - 1]   & ms_nContiguousHighMask )   ||
                   !( ms_pnMask[ii * nShoeBoxSize0 + jj + 1]   & ms_nContiguousHighMask ) )
                {    
                    nPeakBorderPixels++;
                    
                    ms_pnMask[ii * nShoeBoxSize0 + jj] |= ms_nBoundaryPeakPixels;
                }
            }

        }
    }
#ifdef RB_DEBUG    
    nPrintMask(true, ms_nBoundaryPeakPixels);
    nPrintMask(true, ms_nContiguousHighMask);
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
int C3Ddata::nPrintMask(bool bPrintFull,int nMaskToPrint) 
{
    C3Ddata*        poShoebox;
    int a3n[3];
    int nStart;
    int nEnd;
    
    //char acChars[] = { '*','x','o','%','@',0 };
    
    int  anCharAssignment[32];
    
    int  nx,ny;

    poShoebox = this;
    
    //////////////////////////////
    Cstring   strMaskName("");  
    if( nMaskToPrint == ms_nHighMask )
        strMaskName = "High Mask";
    else if(nMaskToPrint == ms_nEllipsoidMask )
        strMaskName = "Ellipsoid Mask";
    else if(nMaskToPrint ==ms_nBorderMask )
        strMaskName = "Border Mask";
    else if(nMaskToPrint == ms_nLowMask )
        strMaskName = "Low Mask";
    else if(nMaskToPrint == ms_nSatMask )
        strMaskName = "Saturated Mask";
    else if(nMaskToPrint == ms_nBadMask )
        strMaskName = "Bad Mask";
    else if(nMaskToPrint == ms_nReservedMask )
        strMaskName = "Reserved Mask";
    else if(nMaskToPrint == ms_nCentroidPixelMask )
        strMaskName = "Centroid Pixel Mask";
    else if(nMaskToPrint == ms_nContiguousHighMask )
        strMaskName = "Contiguous Mask";
    else if(nMaskToPrint == ms_nMinPeakPixels )
        strMaskName = "Min Peak Mask";
    else if(nMaskToPrint == ms_nBoundaryPeakPixels )
        strMaskName = "Boundary Peak Mask";
    ///////////////////////////////
    printf("\n%s\n", strMaskName.string());
    printf("\n----\n");

    for (nx=0,ny=0;nx<32;nx++) {
        if (nMaskToPrint & (1 << nx))
            anCharAssignment[nx] = ny++;
    };

    if (bPrintFull) 
    {
        nStart = 0;
        nEnd = m_nExt[2] - 1;
    } 
    else
        poShoebox->vGetPeakStartEnd(&nStart,&nEnd);
    
    for (a3n[2]=nStart; a3n[2] <= nEnd; a3n[2]++) 
    {
        for (a3n[1] = 0; a3n[1]< poShoebox->m_nExt[1];a3n[1]++)
        {
            for(a3n[0]=0; a3n[0] < poShoebox->m_nExt[0]; a3n[0]++)
            {

                if (nMaskToPrint & ms_pnMask[a3n[2]*poShoebox->m_nExt[0]*poShoebox->m_nExt[1] + 
                                             a3n[1]*poShoebox->m_nExt[0] + a3n[0]]) 
                    printf("*");
                else
                    printf(".");

            };
            
            printf("\n");
        }

        printf("\n");
    }
    
    return 0;
}

int C3Ddata::nLoadPrint(bool bPrintHKLBounds) {
    int nx,ny;
    int n0,n1;
    int nSize;
    int nSliceSize;
    int a2nDelta[2];
    int nSlice;

    nSize = m_nExt[0]*m_nExt[1]*m_nExt[2];
    nSliceSize = m_nExt[0]*m_nExt[1];
    
    // Copy all the data over initially.
    for (nx=0;nx<nSize;nx++)
        ms_pfPrint[nx] = m_pfData[nx];

    if (m_nPeakStart == -1)
        return 0;
    for (nSlice = m_nPeakStart; nSlice <= m_nPeakEnd; nSlice++) {
        
        for (nx = 0,n0=0;n0<m_nExt[0];n0++) {
            for (n1=0;n1<m_nExt[1];n1++) {
                nx = n0 + n1*m_nExt[0];
                if (ms_pnMask[nSlice*nSliceSize + nx] &  ms_nCentroidPixelMask)
                    ms_pfPrint[ nSlice*nSliceSize + nx] = 0.0;
                else if (ms_pnMask[nSlice*nSliceSize + nx] &  ms_nEllipsoidMask) {
                    bool bOnEdge = false;
                    for (a2nDelta[0] = -1; a2nDelta[0] <= 1; a2nDelta[0]++) {
                        for (a2nDelta[1] = -1; a2nDelta[1] <= 1; a2nDelta[1]++) {
                            if ((a2nDelta[0] + n0 >=0) && (a2nDelta[0] + n0 < m_nExt[0]) &&
                                (a2nDelta[1] + n1 >=0) && (a2nDelta[1] + n1 < m_nExt[1])) {

                                ny = (a2nDelta[0] + n0) + m_nExt[0]*(a2nDelta[1] + n1);
                                if (!(ms_pnMask[nSlice*nSliceSize + ny] & ms_nEllipsoidMask)) 
                                    bOnEdge = true;
                            }
                        };
                    };
                    if ((a2nDelta[0] + n0 ==0) && (a2nDelta[0] + n0 == m_nExt[0] - 1) &&
                        (a2nDelta[1] + n1 ==0) && (a2nDelta[1] + n1 == m_nExt[1] - 1)) 
                        bOnEdge = true;
                    if (bOnEdge) {
                        ms_pfPrint[nSlice*nSliceSize + nx] = 0.0;
                    };
                };
            };
        };
    };      

    return 0;
};


int  C3Ddata:: nCalcContiguous(int nSlice,int nFlagsIn,int nFlagsOut,int a2nStart[2]) {




    int*    pnStack;           // Used when finding peak for ellipsoidal modeling.
    int     nStackPoint;
    int     nLowCount;          // Number of low pixels.
    int     nPeakCount;        // Number of high pixels starting at centroid.
    int     nPeakCountAdjPeak; // Number of high pixels that were adjacent to another peak.
    int     nPeakCountAdjPeak6;
    int     nPeakCountAdjPeak7;
    int     nMinPeakCount;
    int     nSliceOffset;
    int     n0,n1,nx;
    int     nFlagsInWithPeak;
    bool    bStrongPeak;


    nSliceOffset = nSlice*m_nExt[0]*m_nExt[1];

    // Clear all flags previously set.
    for (nx=0;nx<m_nExt[0]*m_nExt[1];nx++)
        ms_pnMask[nSliceOffset+nx] &= ~nFlagsOut;

    // Set the minimum peak pixels in the mask.
    nFlagsInWithPeak = nFlagsIn | ms_nMinPeakPixels;
	{
		int a2nRange0[2];
		int a2nRange1[2];
		double a2fAxis[2];
		double a2fAxisCent[2];
		double fDistAxisSq;
		double fDistNormAxisSq;
		double fMaxDistAxisSq;
		double fMaxDistNormAxisSq;
		double f0;
		
		if (m_a2fShiftKa1[0] == m_a2fShiftKa2[0]) {
			m_a2fShiftKa1[0] = 0.01f;
			m_a2fShiftKa2[0] = -0.01f;
		};
		if (m_a2fShiftKa1[1] == m_a2fShiftKa2[1]) {
			m_a2fShiftKa1[1] = 0.01f;
			m_a2fShiftKa2[1] = -0.01f;
		};
		
		a2nRange0[0] = (int) min(m_a2fShiftKa1[0],m_a2fShiftKa2[0]) - m_nMinPeakRad;
		a2nRange1[0] = (int) min(m_a2fShiftKa1[1],m_a2fShiftKa2[1]) - m_nMinPeakRad;
		a2nRange0[1] = (int) max(m_a2fShiftKa1[0],m_a2fShiftKa2[0]) + m_nMinPeakRad;
		a2nRange1[1] = (int) max(m_a2fShiftKa1[1],m_a2fShiftKa2[1]) + m_nMinPeakRad;
		a2fAxis[0] = (m_a2fShiftKa1[0] - m_a2fShiftKa2[0])/2.0;
		a2fAxis[1] = (m_a2fShiftKa1[1] - m_a2fShiftKa2[1])/2.0;
		a2fAxisCent[0] = (m_a2fShiftKa1[0] + m_a2fShiftKa2[0])/2.0;
		a2fAxisCent[1] = (m_a2fShiftKa1[1] + m_a2fShiftKa2[1])/2.0;
		f0 = sqrt(a2fAxis[0]*a2fAxis[0] + a2fAxis[1]*a2fAxis[1]);
		fMaxDistAxisSq = (m_nMinPeakRad + f0)*(m_nMinPeakRad + f0);
		fMaxDistNormAxisSq = m_nMinPeakRad*m_nMinPeakRad;
		a2fAxis[0]/=f0;
		a2fAxis[1]/=f0;
        nMinPeakCount = 0;
		
		for (n0 = a2nRange0[0]; n0 <= a2nRange0[1]; n0++) {
			for (n1 = a2nRange1[0]; n1 <= a2nRange1[1]; n1++) {
				if ((n0+a2nStart[0]>=0) && (n0+a2nStart[0]<m_nExt[0]) && (n1+a2nStart[1]>=0) && (n1+a2nStart[1]<m_nExt[1])) {
					fDistAxisSq = a2fAxis[0]*(n0-a2fAxisCent[0]) + a2fAxis[1]*(n1 - a2fAxisCent[1]);
					fDistAxisSq *= fDistAxisSq;
					fDistNormAxisSq = (n0-a2fAxisCent[0])*(n0-a2fAxisCent[0]) + (n1-a2fAxisCent[1])*(n1-a2fAxisCent[1]) - fDistAxisSq;
					if ((fDistAxisSq<fMaxDistAxisSq) && (fDistNormAxisSq < fMaxDistNormAxisSq)) {
						C2N(n0+a2nStart[0],n1+a2nStart[1],nx);
						ms_pnMask[nSliceOffset + nx] |= ms_nMinPeakPixels;
                        nMinPeakCount++;
					};
				};
			};
		};
		// If you want to see the masks for these...
		//nPrintMask(true,ms_nMinPeakPixels);
		//printf("[%d %d]\n",m_nOrigOffset[0],m_nOrigOffset[1]);
	};

    pnStack     = ms_apnScratchInt[1];
    C2N(a2nStart[0],a2nStart[1],pnStack[0]);    
    pnStack[1] = 0;
    nStackPoint = 2; 
    nLowCount = 0;
    nPeakCount = 0;
    nPeakCountAdjPeak = 0;
    nPeakCountAdjPeak6 = 0;
    nPeakCountAdjPeak7 = 0;

    //nPrintMask(true, ms_nHighMask);
    //nPrintMask(true, ms_nContiguousHighMask);

    //$$ shijie yao : April 23, 2008
    //For certain CBF images, the 
    //  nPeakCountAdjThisPeak >= 6 for nPeakCountAdjPeak6++;
    //  nPeakCountAdjThisPeak >= 7 for nPeakCountAdjPeak7++;
    //are too strong restrict to qualify a continuous region. 
    //Relaxing the constraints to (3,4) seems worked well on CDD
    //images, while rejecting less peak regions for "sharp" peaks.
    //Another option will be to check the image type as the following
    //code block, and set the constraints accordingly.

    int nAdj6Count = 3; // = 6; //original count condition for adj6
    int nAdj7Count = 4; // = 7; //original count condition for adj7

    //** OPTION #2 : set constraints according to image type **
    //
    //Cstring tmpName("PILT_DETECTOR_DESCRIPTION"); //looking for "PILATUS conversion"
    //Cstring tmpVal;
    //m_poImage->m_oHeader.nGetValue(tmpName, &tmpVal);
    //if(tmpVal.contains(Cstring("PILATUS")))
    //{
    //    nAdj6Count = 3;
    //    nAdj7Count = 4;
    //}

    while (nStackPoint--) {
        int nPix0;
        int nPix1;
        int nPix;
        int nPrevMask;      // Mask value of the pixel that launched this pixel.
        int nPeakCountAdjThisPeak;
        
        nPrevMask = pnStack[nStackPoint--];
        nPix = pnStack[nStackPoint];

        N2C(nPix,nPix0,nPix1);
        
        // Continue poping points until there are no more points left in the stack to check.
        
        // Perhaps this pixel was already marked.
        if (ms_pnMask[nSliceOffset + nPix] & nFlagsOut)
            continue;
        
        // Don't count the pixel unless it really is in nFlagsIn.  We include
        // a minimum locus of pixels that will be marked with ms_nMinPeakPixels
        if (ms_pnMask[nSliceOffset + nPix] & (nFlagsIn)) {
            nPeakCount++;
            // We are also interested in how many peaks were adjacent to others.
            // A mask that has a large number of disconnected peaks is not desired.
            if (nPrevMask & nFlagsIn)
                nPeakCountAdjPeak++;
        };
        if (ms_pnMask[nSliceOffset + nPix] & ms_nLowMask)
            nLowCount++;

        // Mark the pixel.
        ms_pnMask[nSliceOffset + nPix] |= nFlagsOut;        
        nPeakCountAdjThisPeak = 0;
        
        for (n0 = -1; n0 <=1; n0++) {
            for (n1 = -1; n1 <=1; n1++) {
                C2N(n0+nPix0,n1+nPix1,nx);
                // If the pixel is a high pixel, AND it is in range
                // AND it is not currently marked.
                if ((n0+nPix0>=0) && (n0+nPix0<m_nExt[0]) && (n1+nPix1>=0) && (n1+nPix1<m_nExt[1])) {
                    if ((ms_pnMask[nx + nSliceOffset] & (nFlagsInWithPeak))
                        && ((ms_pnMask[nx + nSliceOffset] & nFlagsOut)==0)) {
                        pnStack[nStackPoint++] = nx;                      
                        pnStack[nStackPoint++] = ms_pnMask[nSliceOffset + nPix];
                    };
                    if ((nx!=nPix) && (ms_pnMask[nx + nSliceOffset] & nFlagsIn) && (ms_pnMask[nPix + nSliceOffset] & nFlagsIn))
                        nPeakCountAdjThisPeak++;
                };
            };
        };
        //if (nPeakCountAdjThisPeak >= 6)
        //    nPeakCountAdjPeak6++;
        //if (nPeakCountAdjThisPeak >= 7)
        //    nPeakCountAdjPeak7++;

        if (nPeakCountAdjThisPeak >= nAdj6Count)
            nPeakCountAdjPeak6++;
        if (nPeakCountAdjThisPeak >= nAdj7Count)
            nPeakCountAdjPeak7++;
    };

    if (nPeakCount>=nMinPeakCount) {
        // Ignore the min peak shape entirely.
        for (n1=0,nx=0;(n1<m_nExt[1]);n1++) {
            for (n0=0;(n0<m_nExt[0]);n0++,nx++) {
                if (ms_pnMask[nSliceOffset + nx] & ms_nMinPeakPixels) {
                    if (!(ms_pnMask[nSliceOffset + nx] & nFlagsIn)) 
                        ms_pnMask[nSliceOffset + nx] &= ~nFlagsOut;
                };
            };
        };
    };

    //nPrintMask(true, ms_nHighMask);
    //nPrintMask(true, ms_nContiguousHighMask);

    bStrongPeak = nPeakCountAdjPeak7 || (nPeakCountAdjPeak6>=2);	   
    return (!bStrongPeak);
};




int C3Ddata::nCalcEllipsoid(int nSlice,int a2nSliceCent[2],float a2x2fEllipsoidA[2][2],float a2fEllipsoidb[2],float& fEllipsoidc)
{
    int nx = 0;

    // Variables used in ellipsoidal modelling.
    const int nSearchPlanes = 16;
    int aa2nSearchPlanes[nSearchPlanes][2] = {
        {1,-1},{1,1},{3,1},{1,3},{1,-3},{3,-1},{1,0},{0,1},
        {-1,1},{-1,-1},{-3,-1},{-1,-3},{-1,3},{-3,1},{-1,0},{0,-1}
    };
    int aa2nSearchPlanesNormal[nSearchPlanes][2] = {
        {1,-1},{1,1},{1,0},{0,1},{0,-1},{1,0},{1,0},{0,1},
        {-1,1},{-1,-1},{-1,0},{0,-1},{0,1},{-1,0},{-1,0},{0,-1}
    };

    int anSearchPoint[16];
    int anSearchMax[16];
    int aa2nSearchMax[16][2];

    double      a2fDir[2];
    double      pfCoeffA[7][7];
    double      pfCoeffB[7];
    double      a2x2lfGMat[2][2];
    double      a2fhVec[2];
    double      fConst;
    double*     pf0,*pf1;
    int         a2nAvgCent[2];

    int         n0,n1;
    int         nSliceSize;
    int         nSliceOffset;
    int         nPlane;
    double      fFattening;

    nSliceSize = m_nExt[0]*m_nExt[1];
    nSliceOffset = nSliceSize*nSlice;
    fFattening = m_fFattening;

    for (nPlane = 0; nPlane < nSearchPlanes; nPlane++) {
        anSearchMax[nPlane] = -999999;
    };
    for (n1=0,nx=0;(n1<m_nExt[1]);n1++) {
        for (nPlane = 0; nPlane < nSearchPlanes; nPlane++) 
            anSearchPoint[nPlane] = n1*aa2nSearchPlanes[nPlane][1];
        for (n0=0;(n0<m_nExt[0]);n0++,nx++) {
            if ((ms_pnMask[nx + nSliceOffset] & (ms_nContiguousHighMask))) {
                for (nPlane = 0; nPlane < nSearchPlanes; nPlane++) {
                    if (anSearchPoint[nPlane]>anSearchMax[nPlane]) {
                        anSearchMax[nPlane] = anSearchPoint[nPlane];
                        aa2nSearchMax[nPlane][0] = n0;
                        aa2nSearchMax[nPlane][1] = n1;
                    };
                };
            };                                
            for (nPlane = 0; nPlane < nSearchPlanes; nPlane++) 
                anSearchPoint[nPlane] += aa2nSearchPlanes[nPlane][0];
        };
    };

    // Find the average centroid.
    a2nAvgCent[0] = 0;
    a2nAvgCent[1] = 0;
    n0 = 0;
    for (nPlane = 0; nPlane < nSearchPlanes;nPlane++) {
        if (anSearchMax[nPlane] != -999999) {
            a2nAvgCent[0] +=aa2nSearchMax[nPlane][0];
            a2nAvgCent[1] +=aa2nSearchMax[nPlane][1];
            n0++;
        };
    };
    if (!n0) {
        return (1 << Crefln::ms_nErrorIntProblem);
    };
    a2nAvgCent[0] /= n0;
    a2nAvgCent[1] /= n0;
    a2nSliceCent[0] = a2nAvgCent[0];
    a2nSliceCent[1] = a2nAvgCent[1];

    // Construct the ellipsoid for the centroid slice.            
    
    pf0 = & pfCoeffA[0][0];
    pf1 = & pfCoeffB[0];
    nInitQuadratic(2,FALSE,&pf0,&pf1);

    for (nPlane = 0; nPlane < nSearchPlanes; nPlane++) {
        if (anSearchMax[nPlane] != -999999) {
            // Note:  Make sure to keep the addition of aa2nSearchPlanes[] since this elliminates the 
            // possibility that there are less that 6 unique search points required for the system to be solved.
            a2fDir[0] = aa2nSearchPlanesNormal[nPlane][0]*fFattening + ( aa2nSearchMax[nPlane][0] - a2nSliceCent[0]);
            a2fDir[1] = aa2nSearchPlanesNormal[nPlane][1]*fFattening + ( aa2nSearchMax[nPlane][1] - a2nSliceCent[1]);
            nAddQuadratic(2,1.0,1.0,&a2fDir[0],&pfCoeffA[0][0],&pfCoeffB[0]);
        };
    };

    // Add one point in the center so that we have a well defined quadratic function.
    a2fDir[0] = 0.0;
    a2fDir[1] = 0.0;
    nAddQuadratic(2,0.0,1.0,&a2fDir[0],&pfCoeffA[0][0],&pfCoeffB[0]);
    if (nSolveQuadratic(2,&a2x2lfGMat[0][0],&a2fhVec[0],&fConst,&pfCoeffA[0][0],&pfCoeffB[0]))
        return (1 << Crefln::ms_nErrorIntProblem);
    // Copy the ellipsoid over so that we can use it in the future.
    vCopyVecND(4,&a2x2lfGMat[0][0],&a2x2fEllipsoidA[0][0]);
    vCopyVecND(2,&a2fhVec[0],&a2fEllipsoidb[0]);            
    fEllipsoidc = fConst;

    return 0;
};

int C3Ddata::nCheckFuzzySpotShape(int nSlice) {
    int nStatSlice = 0;
    int nPeakCount = 0;
    int nAreaCount = 0;
    int n0,n1,nx;
    int nSliceOffset;
    int nSliceSize;

    nSliceSize = m_nExt[0]*m_nExt[1];
    nSliceOffset = nSlice*nSliceSize;

    
    // check for on edge high pixels.  
    // But only do this for strong reflections.
    for (n1=0,nx=0;(n1<m_nExt[1]) ;n1++) {
        for (n0=0;(n0<m_nExt[0]);n0++,nx++) {
            if (ms_pnMask[nx + nSliceOffset] & ms_nContiguousHighMask) {
                nAreaCount++;
                if (ms_pnMask[nx + nSliceOffset] & ms_nEllipsoidMask) 
                    nPeakCount++;
            };
        };
    };
    if (((float) (nAreaCount - nPeakCount))/nAreaCount>0.15)
        nStatSlice |= (1 << Crefln::ms_nErrorUnaccountedPix);
    return nStatSlice;
};

int C3Ddata::nCalcEllipsoidBitmask(int nSlice,int a2nSliceCent[2],int nBitsToSet,float a2x2fEllipsoidA[2][2],float a2fEllipsoidb[2],float fEllipsoidc,int a2x2nSearch[2][2],float* pfEllipsoidalBackground) {
    int n0,n1,nx;
    int nSliceOffset;
    float a2fTemp1[2];
    float a2fTemp2[2];
    double f0,f1;
    double fK2,fK3;
    int nTargetBackgroundPixels;
    int nCurrentPixels;        
    const int nBackgroundPixels = 50;

    nSliceOffset = nSlice*m_nExt[0]*m_nExt[1];
    a2x2nSearch[0][0] = a2nSliceCent[0];
    a2x2nSearch[0][1] = a2nSliceCent[0];
    a2x2nSearch[1][0] = a2nSliceCent[1];
    a2x2nSearch[1][1] = a2nSliceCent[1];

    // Calculate the background ellipse size for background calculation.
    
    if (pfEllipsoidalBackground) {
        double a2fRelOffset[2];
        double a2fMajorMinor[2];
        double fRotOffset;
        double a7fEllipsoid[7];

        vCopyVecND(4,&a2x2fEllipsoidA[0][0],&a7fEllipsoid[0]);
        vCopyVecND(2,&a2fEllipsoidb[0],&a7fEllipsoid[4]);
        a7fEllipsoid[6] = fEllipsoidc;

        nGetEllipsoidTiltAngleCenter(&a7fEllipsoid[0],&a7fEllipsoid[4],a7fEllipsoid[6],&a2fRelOffset[0],&a2fMajorMinor[0],&fRotOffset);
        nCurrentPixels = max((int) (a2fMajorMinor[0]*a2fMajorMinor[1]*Gs_fPI),1);
        nTargetBackgroundPixels = max(nBackgroundPixels,(int) (nCurrentPixels*2.0));
        fK2 = nTargetBackgroundPixels/(a2fMajorMinor[0]*a2fMajorMinor[1]*Gs_fPI) + 1.0;
        nTargetBackgroundPixels = max(nBackgroundPixels,(int) (nCurrentPixels*1.5));
        fK3 = nTargetBackgroundPixels/(a2fMajorMinor[0]*a2fMajorMinor[1]*Gs_fPI) + 1.0;
    } else {
        fK2 = 1e6;
        fK3 = 1e6;
    };

    

    double a3x3fMatA[3][3];
    double a3x3fMatAInv[3][3];
    double a3x3fMatAAdd[3][3];
    double a3fMatb[3];
    vZeroMat(3,3,&a3x3fMatA[0][0]);
    vZeroMat(3,1,&a3fMatb[0]);
        
    for (n1=0,nx=0;n1<m_nExt[1];n1++) {
        for (n0=0;n0<m_nExt[0];n0++,nx++) {
            ms_pnMask[nx + nSliceOffset] &= ~nBitsToSet;
            // if this pixel is in the mask.
            a2fTemp1[0] = n0 - a2nSliceCent[0];
            a2fTemp1[1] = n1 - a2nSliceCent[1];
            vMulMatNDVecND(2,&a2x2fEllipsoidA[0][0],&a2fTemp1[0],&a2fTemp2[0]);
            f1 = a2fTemp1[0]*a2fTemp2[0] + a2fTemp1[1]*a2fTemp2[1] + a2fTemp1[0]*a2fEllipsoidb[0] + a2fTemp1[1]*a2fEllipsoidb[1] + fEllipsoidc;
            if (f1<=1.0) {
                // We have a mask pixel.  
                a2x2nSearch[0][0] = min(n0,a2x2nSearch[0][0]);
                a2x2nSearch[1][0] = min(n1,a2x2nSearch[1][0]);
                a2x2nSearch[0][1] = max(n0,a2x2nSearch[0][1]);
                a2x2nSearch[1][1] = max(n1,a2x2nSearch[1][1]);

                ms_pnMask[nx + nSliceOffset] |= nBitsToSet;
                
            } else if ((f1 < fK2) && (f1 > fK3) && (m_pfData[nSliceOffset + nx] > 0.0) && (!(ms_pnMask[nSliceOffset + nx] & (ms_nHighMask | ms_nBorderMask)))) {
                f0 = ((double) n0)/m_nExt[0];
                f1 = ((double) n1)/m_nExt[1];

                a3x3fMatAAdd[0][0] = f0*f0;
                a3x3fMatAAdd[1][0] = f1*f0;
                a3x3fMatAAdd[2][0] = f0;
                a3x3fMatAAdd[1][1] = f1*f1;
                a3x3fMatAAdd[2][1] = f1;
                a3x3fMatAAdd[2][2] = 1;
                a3x3fMatAAdd[0][1] = a3x3fMatAAdd[1][0];
                a3x3fMatAAdd[0][2] = a3x3fMatAAdd[2][0];
                a3x3fMatAAdd[1][2] = a3x3fMatAAdd[2][1];
                vAddVecNDVecND(9,&a3x3fMatAAdd[0][0],&a3x3fMatA[0][0],&a3x3fMatA[0][0]);
                a3fMatb[0] += m_pfData[nSliceOffset + nx]*f0;
                a3fMatb[1] += m_pfData[nSliceOffset + nx]*f1;
                a3fMatb[2] += m_pfData[nSliceOffset + nx];
            };
        };
    };    
    if ((a3x3fMatA[2][2] > nBackgroundPixels/2) && (pfEllipsoidalBackground)) {
        if (0.0 != fInvMat3D(&a3x3fMatA[0][0],&a3x3fMatAInv[0][0])) {
            double a3fTemp[3];
            vMulMat3DVec3D(a3x3fMatAInv,a3fMatb,a3fTemp);
            vCopyVec3D(a3fTemp,pfEllipsoidalBackground);
        };
    };


    return 0;
};

int C3Ddata::nFindBestCentroidSlice(const float a3fCent[3],const float* pfPerSliceIntensity,int nMosPeakStartSearch,int nMosPeakEndSearch) {
    int nNumSlicesToTest;
    int a3nSlicesToTest[3];
    double f0;
    int n0;
    int nx,ny;
    const double fBorderTol = 0.4;
    
    f0 = a3fCent[2] - floor(a3fCent[2]);
    n0 = (int) floor(a3fCent[2]);
    nNumSlicesToTest = 0;
    if ((n0 >=nMosPeakStartSearch) && (n0 <= nMosPeakEndSearch))
        a3nSlicesToTest[nNumSlicesToTest++] = n0;

    if (n0 - 1>=nMosPeakStartSearch)
        a3nSlicesToTest[nNumSlicesToTest++] = n0 - 1;
    if (n0 + 1 <= nMosPeakEndSearch)
        a3nSlicesToTest[nNumSlicesToTest++] = n0 + 1;
    
    if (nNumSlicesToTest==1)
        return a3nSlicesToTest[0];

    if (nNumSlicesToTest==0)
        return (int) a3fCent[2];
       
    for (nx=1,ny=0;nx < nNumSlicesToTest;nx++) {
        if (pfPerSliceIntensity[a3nSlicesToTest[nx]]>pfPerSliceIntensity[a3nSlicesToTest[ny]])
            ny = nx;
    };
    return a3nSlicesToTest[ny];
};

void C3Ddata::vGetPerSliceEllipsoidRelObsCent(int nSlice,float* a2x2fEllipsoidA,float* a2fEllipsoidb,float* pfEllipsoidc,float* pfAbsObsCentroidAlternate) {
    float a2fRelCentroid[2];
    vCopyVecND(4,&m_apfScratchFloat[m_nScratchFloatEllipsoidA][nSlice*4],a2x2fEllipsoidA);
    vCopyVecND(2,&m_apfScratchFloat[m_nScratchFloatEllipsoidb][nSlice*2],a2fEllipsoidb);
    *pfEllipsoidc = m_apfScratchFloat[m_nScratchFloatEllipsoidc][nSlice];
    if (pfAbsObsCentroidAlternate) {
        a2fRelCentroid[0] = pfAbsObsCentroidAlternate[0] - m_nOrigOffset[0];
        a2fRelCentroid[1] = pfAbsObsCentroidAlternate[1] - m_nOrigOffset[1];
    } else {
        a2fRelCentroid[0] = m_a3fCentroid[0];
        a2fRelCentroid[1] = m_a3fCentroid[1];
    };

    vConvertEllipsoidRelative(a2x2fEllipsoidA,a2fEllipsoidb,pfEllipsoidc,0,0,a2fRelCentroid[0],a2fRelCentroid[1]);

};

int C3Ddata::nCountZeros(long *plNumberOfZeroPixels, const float fBadValue)
{
  //+2009-07-18 JWP
  // In preparation of changing this routine to count bad pixels instead
  // of just pixels with a pixel value of 0.0, I have added the fBadValue
  // argument which is set to 0.0 if missing from the call of the method.
  // Thus it should work as before in legacy code.
  // The problem is that for Pilates images, a value of 0 is not a bad pixel.
  // And pixel values of -1 and -2 are bad.
  //-2009-07-18 JWP

  // Compute the number of zero pixels in a 3Ddata object

  int   i;
  int   nStat;
  long lZero,lNumPix;
  float *pfTemp;

  // Make sure data is converted to floating point before proceeding

  if (e3Ddata_float != m_eType)
    {
      nStat =  nConvertToFloat();
      cout << "WARNING: converting to float in C3Ddata::nCountZeros()!\n";
      if (0 != nStat)
	{
	  return (nStat);
	}
    }

  pfTemp = m_pfData;
  lZero = 0;

  lNumPix = (long)m_nExt[0] * m_nExt[1] * m_nExt[2]; // Tot. num of pixels to loop thru
  for (i = 0; i < lNumPix; i++)
    {
      if(fBadValue >= *pfTemp)
         lZero++;
      pfTemp++;
    }

  *plNumberOfZeroPixels = lZero;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///////////////////// C3DdataDetectorProfile ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void C3DdataDetectorArea::vZero()
{ 
    memset(this, 0, sizeof(C3DdataDetectorArea)); 
};

void C3DdataDetectorArea::vDup() {
    fIntensity = 1000.0;
    fSigma = 100.0;
    fBackground = 10.0;
    fBackgroundSigma = 1.0;
    fSize[0] = 12;
    fSize[1] = 12;
    fSize[2] = 1;
    nBuildEllipse(fSize[0],&a2x2fEllipsoidA[0][0],&a2fEllipsoidb[0],&fEllipsoidc);
    fShoebox[0] = 30;
    fShoebox[1] = 30;
};

int C3DdataDetectorArea::nReadWriteHeader(bool bRead,int nIndex,Cstring& sBuffer,Cstring& sName)
{
    
    char* apcNames[] =
    {
        "EllipsoidA00",
        "EllipsoidA01",
        "EllipsoidA10",
        "EllipsoidA11",
        "Ellipsoidb0",
        "Ellipsoidb1",
        "Ellipsoidc",
        "Intensity",
        "Sigma",
        "Background",
        "BackgroundSigma",
        "Size0",
        "Size1",
        "Size2",
        "Shoebox0",
        "Shoebox1",
        NULL
    };

    double* apfPointers[] = 
    {
            &a2x2fEllipsoidA[0][0],
            &a2x2fEllipsoidA[1][0],
            &a2x2fEllipsoidA[0][1],
            &a2x2fEllipsoidA[1][1],
            &a2fEllipsoidb[0],
            &a2fEllipsoidb[1],
            &fEllipsoidc,
            &fIntensity,
            &fSigma,
            &fBackground,
            &fBackgroundSigma,
            &fSize[0],
            &fSize[1],
            &fSize[2],
            &fShoebox[0],
            &fShoebox[1]
    };

    double      f0 = 0.0;
    char        pcBuf[100];
    
    Cstring     sTemp("");
    
    if( bRead )
    {
        if (!apcNames[nIndex])
            return 1;
        
        if (!sName.length())
        {
            sName= apcNames[nIndex];
            return 0;
        } 
        else if (!sBuffer.length())
        {
            return 1;
        }

        sTemp = sBuffer.before(' ');
        
        if (1 != sscanf(sTemp.string(),"%lf",&f0))
            return 1;
        
        *apfPointers[nIndex] = f0;
        
        sBuffer = sBuffer.after(' ');
    }
    else
    {
        if (!apcNames[nIndex])
            return 1;
        
        if (!sName.length()) 
            sName = apcNames[nIndex];
        
        sprintf(pcBuf,"%.5lf ",*apfPointers[nIndex]);
        
        sBuffer += pcBuf;        
    }
    
    return 0;
}


C3DdataDetectorProfile::C3DdataDetectorProfile(int nDim0,int nDim1) 
{
    m_nDim0 = 0;
    m_nDim1 = 0;
    m_nContrib = 0;
    m_bMeshInitialized = false;
    m_bHasPrintAreaData = false;

    for(int ii=0; ii < C3Ddata_Print_Div_0; ii++)
    {
        for(int jj=0; jj < C3Ddata_Print_Div_1; jj++)
        {
            m_aoDetectorAreasPrint[ii][jj].vZero();
            m_aoDetectorAreasVarPrint[ii][jj].vZero();
        }
    }

    if (m_nDim0>0)
        nInitPeakAreas(m_nDim0,m_nDim1);
}

int  C3DdataDetectorProfile::nAddPeakArea(int nPix0,int nPix1,C3Ddata* poShoebox) {
    C3DdataDetectorArea oData;
    float a3fSize[3];
    //double f0;

    vCopyVecND(4,poShoebox->pfGetEllipsoidA(),&oData.a2x2fEllipsoidA[0][0]);
    if (oData.a2x2fEllipsoidA[0][0] == 0.0)
        return 1;
    oData.a2fEllipsoidb[0] = poShoebox->pfGetEllipsoidb()[0];
    oData.a2fEllipsoidb[1] = poShoebox->pfGetEllipsoidb()[1];
    oData.fEllipsoidc = poShoebox->fGetEllipsoidc();
    oData.fIntensity = poShoebox->fGetPeakIntensity();
    oData.fSigma = poShoebox->fGetPeakIntensitySigma();
    oData.fBackground = poShoebox->fGetPerPixelBackground();
    oData.fBackgroundSigma = poShoebox->fGetPerPixelBackgroundSigma();
    poShoebox->vGetPeakSize(&a3fSize[0]);
    oData.fSize[0] = a3fSize[0];
    oData.fSize[1] = a3fSize[1];
    oData.fSize[2] = a3fSize[2];
    oData.fShoebox[0] = poShoebox->m_nExt[0];
    oData.fShoebox[1] = poShoebox->m_nExt[1];
    m_oDetectorAreas.nAdd(nPix0,nPix1,(double*) &oData);
    return 0;

};

int C3DdataDetectorProfile::nUpdateHeader(Cimage_header *poHeader)
{
    int nArea;
    int nAreas;
    int nVar;
    Cstring sLabel;
    Cstring sData;
    Cstring sTemp;
    int nLoop;
    char* acPrefix[2] = {"DA_","DAV_"};
    if (nPrintPeakAreas(true,1.0)) 
        return 1;

    nAreas = C3Ddata_Print_Div_1 * C3Ddata_Print_Div_0;
    
    C3DdataDetectorArea*        poArea = NULL;
    
    for(nLoop = 0; nLoop < 2; nLoop++)
    {
        for (nVar = 0; 1; nVar++)
        {
            sData = "";
            sLabel = "";
            
            for(nArea = 0; nArea < nAreas; nArea++)
            {
                if( 0 == nLoop )
                    poArea = &m_aoDetectorAreasPrint[0][0]    + nArea;
                else
                    poArea = &m_aoDetectorAreasVarPrint[0][0] + nArea;
                
                if( poArea->nReadWriteHeader(false, nVar, sData, sLabel) )
                    break;
            }
            
            if( nArea == 0 )
                break;
            else if( nArea < nAreas )
                return 1;
            
            sTemp = acPrefix[nLoop];
            
            sTemp += sLabel;
            
            sLabel = sTemp;
            
            if (poHeader->nReplaceValue(sLabel,sData))
                return 1;
        }
    }
    
    poHeader->nReplaceValue((Cstring) "DA_DIM0", m_nDim0);
    poHeader->nReplaceValue((Cstring) "DA_DIM1", m_nDim1);
    
    return 0;
}

int C3DdataDetectorProfile::nInitValues(Cimage_header& oHeader) {
    int nArea;
    int nAreas;
    int nVar;
    int nLoop;
    Cstring sLabel;
    Cstring sData;
    Cstring sTemp;
    char* acPrefix[2] = {"DA_","DAV_"};

    memset(&m_aoDetectorAreasPrint[0][0],0,sizeof(C3DdataDetectorArea)*C3Ddata_Print_Div_0*C3Ddata_Print_Div_1);
    memset(&m_aoDetectorAreasVarPrint[0][0],0,sizeof(C3DdataDetectorArea)*C3Ddata_Print_Div_0*C3Ddata_Print_Div_1);
    
    nAreas = C3Ddata_Print_Div_1*C3Ddata_Print_Div_0;
    m_bHasPrintAreaData = false;
    if (oHeader.nGetValue((Cstring) "DA_DIM0",&m_nDim0))
        return 1;
    if (oHeader.nGetValue((Cstring) "DA_DIM1",&m_nDim1))
        return 1;

    for (nLoop = 0; nLoop < 2; nLoop++) {
        for (nVar = 0; 1; nVar++) {
            sData = "";
            sLabel = "";
            if (m_aoDetectorAreasPrint[0][0].nReadWriteHeader(true,nVar,sData,sLabel))
                break;
            sTemp = (acPrefix[nLoop]);
            sTemp += sLabel;
            sLabel = sTemp;
            if (oHeader.nGetValue(sLabel,&sData)) {
                if ((nVar ==0) && (nLoop ==1))
                    break;
                else
                    return 1;
            }
            
            for (nArea = 0; nArea < nAreas; nArea++) {
                C3DdataDetectorArea* poArea = ((nLoop==0)?((& m_aoDetectorAreasPrint[0][0]) + nArea):((& m_aoDetectorAreasVarPrint[0][0]) + nArea));
                if (poArea->nReadWriteHeader(true,nVar,sData,sLabel))
                    break;
            };
            if (nArea == 0)
                break;
            else if (nArea < nAreas) 
                return 1;
            
        };
    };
    m_bHasPrintAreaData = true;
    return 0;
};

int  C3DdataDetectorProfile::nGetPeakArea(C3DdataDetectorArea& oArea) {
    int nStat;

    nStat = nGetPeakArea(m_nDim0/2,m_nDim1/2,oArea,0);  
    return nStat;
};



int C3DdataDetectorProfile::nGetPeakArea(int nPix0,int nPix1,C3DdataDetectorArea& oArea,C3DdataDetectorArea* poAreaVariance,int nLevelToUse) {
    int nStat;
    
    if (bIsInitialized()) {
        nStat = m_oDetectorAreas.nGet(nPix0,nPix1,(double*) &oArea,(double*) poAreaVariance,nLevelToUse);  
        // If we specified a level to use, then try relaxing the condition.
        if ((nStat==-1) && (nLevelToUse>=0))
            nStat = m_oDetectorAreas.nGet(nPix0,nPix1,(double*) &oArea,(double*) poAreaVariance);  
        if (nStat>=0)
            m_nContrib = m_oDetectorAreas.m_anInterpCounts.last();
    } else if (m_bHasPrintAreaData) {
        int nArea0,nArea1;
        nArea0 = max(0,min(C3Ddata_Print_Div_0-1,(C3Ddata_Print_Div_0*nPix0)/m_nDim0));
        nArea1 = max(0,min(C3Ddata_Print_Div_1-1,(C3Ddata_Print_Div_1*nPix1)/m_nDim1));
        memcpy(&oArea,&m_aoDetectorAreasPrint[nArea0][nArea1],sizeof(oArea));
        if (poAreaVariance)
            memcpy(poAreaVariance,&m_aoDetectorAreasVarPrint[nArea0][nArea1],sizeof(oArea));
        m_nContrib = 1;
        nStat = 0;        
    } else {
        nStat = -1;
        m_nContrib = 0;
    };
    
    if (nStat < 0) {
        oArea.vDup();
        if (poAreaVariance)
            poAreaVariance->vZero();
        return nStat;
    };

    return nStat;
};


/*
int C3DdataDetectorProfile::nGetPeakArea(int nPix0,int nPix1,C3DdataDetectorArea& oArea,C3DdataDetectorArea* poAreaVariance,int nLevelToUse) {
    int nStat;

    m_nContrib = 0;

    nStat = -1;
    if (m_bHasPrintAreaData) {
        int nArea0,nArea1;
        nArea0 = max(0,min(C3Ddata_Print_Div_0-1,(C3Ddata_Print_Div_0*nPix0)/m_nDim0));
        nArea1 = max(0,min(C3Ddata_Print_Div_1-1,(C3Ddata_Print_Div_1*nPix1)/m_nDim1));
        memcpy(&oArea,&m_aoDetectorAreasPrint[nArea0][nArea1],sizeof(oArea));
        if (poAreaVariance)
            memcpy(poAreaVariance,&m_aoDetectorAreasVarPrint[nArea0][nArea1],sizeof(oArea));
        m_nContrib = 1;
        nStat = 0;        
    }     


    if ((bIsInitialized()) && (nStat==-1)) {
        nStat = m_oDetectorAreas.nGet(nPix0,nPix1,(double*) &oArea,(double*) poAreaVariance,nLevelToUse);  
        // If we specified a level to use, then try relaxing the condition.
        if ((nStat==-1) && (nLevelToUse>=0))
            nStat = m_oDetectorAreas.nGet(nPix0,nPix1,(double*) &oArea,(double*) poAreaVariance);  
        if (nStat>=0)
            m_nContrib = m_oDetectorAreas.m_anInterpCounts.last();
    } 
    
    if (nStat < 0) {
        oArea.vDup();
        if (poAreaVariance)
            poAreaVariance->vZero();
        return nStat;
    };

    return nStat;
};
*/

int  C3DdataDetectorProfile::nInitPeakAreas(int nDim0,int nDim1) 
{
    m_nDim0 = nDim0;
    m_nDim1 = nDim1;

    m_oDetectorAreas.vSetInterpConditions(10,0.0 /*8*/);
    m_oDetectorAreas.vInit(10,sizeof(C3DdataDetectorArea)/sizeof(double),0,m_nDim0,0,m_nDim1);    
    
    m_aoDetectorAreasPrint[0][0].fSize[0] = 0.0;
    
    m_bMeshInitialized = true;
    
    return 0;    
};

int  C3DdataDetectorProfile::nPrintPeakAreas(bool bComputeOnly,double fRotIncrement)
{
    int nUsed;
    
    static char acLine[] = "===============================================================================\n";
    int nRegion0,nRegion1;
    int nPix0,nPix1;
    double fDims[2];
    double fShift[2];
    int    nAvgContrib;
    double fOffset;
    
    double a2x2fInvSlopeMat[2][2] = {{-1,2},{0.5,-0.5}};
    
    C3DdataDetectorArea         oArea;
    C3DdataDetectorArea         oAreaVar;
    C3DdataDetectorArea         oAvgArea;
    
    oArea.vZero();
    oAreaVar.vZero();
    oAvgArea.vZero();

    // First, see if we have enough data to print!
    if (nGetPeakArea(oArea)<0) {
        printf("WARNING:  No strong spots found on first pass of integration:\n"
               "          Using estimated defaults!\n");
        if (!bComputeOnly)
            printf("WARNING: Cannot Print strong peak table.\n");
        return 1;
    };
    
    if (!bComputeOnly) {
        printf("\nStrong peak info listing\n");    
        printf(acLine);
        printf("  Pix  Pix  Num  Counts  Sigma Back Back Box Box Major Minor  Rot. Shift Shift\n");
        printf("  [0]  [1]                      Avg  Sig [0] [1]  Axis  Axis  Axis   [0]   [1]\n");
        printf(acLine);
    };

    memset(&oAvgArea,0,sizeof(oAvgArea));
    nAvgContrib = 0;

    if (m_oDetectorAreas.nGetNumContrib()) 
        m_bHasPrintAreaData = false;
    
    
    for (nRegion0=0;nRegion0<C3Ddata_Print_Div_0;nRegion0++) {
        for (nRegion1=0;nRegion1<C3Ddata_Print_Div_1;nRegion1++) {
            nPix0 = (int) ((m_nDim0*(nRegion0 + 0.5))/C3Ddata_Print_Div_0);
            nPix1 = (int) ((m_nDim1*(nRegion1 + 0.5))/C3Ddata_Print_Div_1);
            
            nGetPeakArea(nPix0,nPix1,oArea,&oAreaVar,min(C3Ddata_Print_Div_0,C3Ddata_Print_Div_1) - 1);
            nUsed = m_nContrib;
            // Copy the area over.  This is used in default spot size determination.
            memcpy(&m_aoDetectorAreasPrint[nRegion0][nRegion1],&oArea,sizeof(C3DdataDetectorArea));
            memcpy(&m_aoDetectorAreasVarPrint[nRegion0][nRegion1],&oAreaVar,sizeof(C3DdataDetectorArea));
            
            nGetEllipsoidTiltAngleCenter(&oArea.a2x2fEllipsoidA[0][0],&oArea.a2fEllipsoidb[0],oArea.fEllipsoidc,&fShift[0],&fDims[0],&fOffset);
                        
            if (!bComputeOnly) {
                printf(" %4d %4d %4d %7.0f %6.0f %4.0f %4.0f %3d %3d %5.1f %5.1f %5.1f %5.1f %5.1f\n",                           
                    nPix0,nPix1,
                    nUsed,
                    oArea.fIntensity,
                    oArea.fSigma,
                    oArea.fBackground,
                    oArea.fBackgroundSigma,
                    (int) oArea.fShoebox[0],(int) oArea.fShoebox[1],
                    2.0*fDims[0],2.0*fDims[1],
                    oArea.fSize[2]*fRotIncrement,
                    fShift[0],fShift[1]
                    );
            };
	    vAddVecNDVecND(sizeof(oAvgArea)/sizeof(double),(double*) 
			   &oAvgArea,(double*) &oArea,(double*) &oAvgArea);
            nAvgContrib++;
        };
    };
    if (!bComputeOnly) {
        C3DdataDetectorArea& oArea = oAvgArea;
        printf(acLine);
        vMulVecNDScalar(sizeof(oAvgArea)/sizeof(double),(double*) &oAvgArea,1.0/nAvgContrib,(double*) &oAvgArea);
        nGetEllipsoidTiltAngleCenter(&oArea.a2x2fEllipsoidA[0][0],&oArea.a2fEllipsoidb[0],oArea.fEllipsoidc,&fShift[0],&fDims[0],&fOffset);
        printf(" %4s %4s %4s %7.0f %6.0f %4.0f %4.0f %3d %3d %5.1f %5.1f %5.1f %5.1f %5.1f\n",                           
            " All","area",
            "avg:",
            oArea.fIntensity,
            oArea.fSigma,
            oArea.fBackground,
            oArea.fBackgroundSigma,
            (int) oArea.fShoebox[0],(int) oArea.fShoebox[1],
            2.0*fDims[0],2.0*fDims[1],
            oArea.fSize[2]*fRotIncrement,
            fShift[0],fShift[1]
            );
        printf(acLine);
        fflush(stdout);
    };
    m_bHasPrintAreaData = true;

    
    return 0;
}

int 
C3Ddata::nResetLowValues(const float fMinValue, const float fNewValue)
{
  int nCount = 0;
  
  // Compute the number of zero pixels in a 3Ddata object

  int   i;
  int   nStat;
  long lZero,lNumPix;
  float *pfTemp;

  // Make sure data is converted to floating point before proceeding

  if (e3Ddata_float != m_eType)
    {
      nStat =  nConvertToFloat();
      cout << "WARNING: converting to float in C3Ddata::nResetLowValues()!\n";
      if (0 != nStat)
	{
	  return (nStat);
	}
    }

  pfTemp = m_pfData;
  nCount = 0;

  lNumPix = (long)m_nExt[0] * m_nExt[1] * m_nExt[2]; // Tot. num of pixels to loop thru
  for (i = 0; i < lNumPix; i++)
    {
      if (fMinValue >= *pfTemp)
	{
	  nCount++;
	  *pfTemp = fNewValue;
	}
      pfTemp++;
    }
  return (nCount);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
