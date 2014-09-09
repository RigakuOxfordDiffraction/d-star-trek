#ifndef DT_CBACKGROUND_H
#define DT_CBACKGROUND_H
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
// CBackground.h     Initial author: Thaddeus Niemeyer        15-Nov-2001
//    Background tools
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


#include "Cimage.h"

/*  Recent changes:

	1) Remove nBuildBackgroundMap().  This function used statistical calculations, and that is bad for a spot finder.  
	   We have known that this function is unsatisfactory for some time.
	2) Added nBuildPeakPixels() function.  
	3) Added functions that load contiguous pixel information.
	3) Added nCalcBackground() function.  This uses code that was originally in nBuildBackgroundMap() and can be called at any time to
	   reload the background values for the image.

*/

class Cdetector;
class Csource;

class DTREK_EXPORT Cbackground {

    float*  m_pfBackground;          // Background mask.  Masks out peaks and user determined regions.
    unsigned char*   m_pcBackground;          // 1 = peak 0 = background
    int     m_nBackgroundParams[4];  // [0] = grid spacing in 0, [1] = grid spacing in [1], [2] = size of grid in [0]. [3] = size of grid in [1].

    int m_nDim[2];
    Cimage& m_oCi;
    int     nAllocateBackground();
    static  unsigned char ms_cZero;

public:
    int     m_nMaskValue;
    int     m_nTempValue;

    bool    bInRange(int nx,int ny)              { return ((nx>=0) && (ny>=0) && (nx<m_nDim[0]) && (ny<m_nDim[1]) && (rcGetBackground(nx,ny)!=m_nMaskValue)); };
    float&  rfGetBackground(int nPix0,int nPix1) { return m_pfBackground[(nPix1/m_nBackgroundParams[1])*m_nBackgroundParams[2] + (nPix0/m_nBackgroundParams[0])]; };
    float   fGetSmoothBackground(int nPix0,int nPix1);
    unsigned char&   rcGetBackground(int px0,int px1); 
    unsigned char&   rcGetBackground(int px0,int px1,unsigned char* pcUser);
	int     nCalcBackground(int nNumRegions,float fSigmaReject,Cimage* poBackgroundImage);
    int     nBuildAdjacentToPeaks(int nLevels);
    int     nPrintBackground(Cstring& sName);
    int     nPrintPeaks(Cstring& sName);
    int     nPrintAdjacentToPeaks(Cstring& sName);
	int     nBuildPeakPixels(int nBoxDim,bool bSaveBackground = false);

	int		nFindSpot(int nPix0,int nPix1,itr<int>& anPix0,itr<int>& anPix1,itr<int>& anStack0,itr<int>& anStack1,int& nCount,int& nMaxCount);
	int		nFindSpots(itr<int>& anPix0,itr<int>& anPix1,itr<int>& anDim0,itr<int>& anDim1,itr<int>& anTotal,itr<int>& anCount,itr<int>& anMax,itr<int>* panCumulPix0 = NULL,itr<int>* panCumulPix1 = NULL);

    int     nExcludeZingers(int nMinAdjacent,int nNewValue);  // Sets zinger peaks (masked pixels with less than nMinAdjacent neighbors that are masked) to 
    int     nExcludeMinimum(double fMinimum);
    int     nExcludeAdjacentToExcluded(int nSearchValue);
    int     nFind2D(Creflnlist& oList,const double fRot,double fSignificanceLevel,Cdetector* poDetector,Csource* poSource);
	int		nIntegrate(int nBoxSize,itr<double>& af2Theta,double f2ThetaStep,Cdetector* poDetector,Csource* poSource,int nCompression = 1);
    int     nFatten(int nFattenCount);
    int     nMaskRegion(int n00,int n01,int n10,int n11,bool bInvert);
    int     nClearMasks(int nNewValue = 0);
    int     nClearAll();
    
    Cbackground(Cimage& oImage);
    ~Cbackground();
};
#endif   // DT_CBACKGROUND_H//
