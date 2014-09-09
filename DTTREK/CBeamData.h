#ifndef DT_CBEAMDATA_H
#define DT_CBEAMDATA_H
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CBeamData.h     Initial author: Thaddeus Niemeyer        15-Nov-2001
//    Beam Center and Beam Stop tools
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
#include "Cimage_header.h"
#include "CBackground.h"
#include "Cdetector.h"
#include "Csource.h"

const int nMaxDPoints = 200;

class Cbeamdata {
    int m_nDim[2];                  // Dimensions of image. (I try to avoid mutiple calls to Cimage::nGetDimension())

    Cimage& m_oCi;                  // Reference to image.
    Cimage  m_oCompImage;           // If we compressed the image, this pointer will be non-null, and m_oCi will reference it.
    Cdetector* m_poDetector;        // Detector pointer (out of convenience)
    Csource* m_poSource;            // Source pointer (out of convenience)

    char*  m_pc2thetaProfile;       // Constructed from user supplied peaks.  Contains a '1' at peak locations, and a '0' at non-peak locations.
    
    float  m_f2lambda[2];           // KA1, KA2
    float  m_a3fS0Wave[3];          // S0 vector.
    float  m_fWavelength;

    char*  m_pcInt;                 // Integral of pixels.  Contains a '1' when above m_fIntSigma '2' when invalid and '0' when below m_fIntSigma.
    char*  m_pcCalcInt;             // Used during integration. Same size as m_pcDistInt, but it contains predicted locations.
    char*  m_pcCalcInt2;

    double** m_ppfSpotValue;        // Double array defining the spot.
    double** m_ppfSpotCount;
    int      m_nSpotDim;            // Each dimension in the double array.
    int      m_nCenter;
    int      m_nDataPerPixel;
    double   m_fFWHM;               // Full width at Half max.
    double   m_fFW10M;              // Full width at 10 percent.
    int      m_nCompression;        // The compression factor.


public:
    Cbackground m_oBackground;      // Background for the image.  

    bool    bInRange(int nx,int ny)             { return ((nx>=0) && (ny>=0) && (nx<m_nDim[0]) && (ny<m_nDim[1])); };

    int     nRejectRectangle(int npx00,int npx01,int npx10,int npx11);
    
    int     nBuildAveragePeak(int nDataPerPixel,float fRejectPercent);
    int     nPrintAveragePeak(Cstring& sOutput);
    int     nPrintAveragePeakStats(Cstring& sOutputFileStub);

    int     nFindBeamStopShadow(int nInitCenter0,int nInitCenter1);
    int     nFindResoRings(double d2ThetaMin_DEG,
                           double dMaxExpected2ThetaRingWidth_DEG,
                           int nImageSizeCompressionFactor,
                           itr<double>& adExpectedResoRings_A,
                           Cimage_header* poHeader);
    
    int     nFitParams(double fReso,double f02Theta,double fRad2Theta,double f0Dist,double fRadDist,double fRadBeam);
    int     nInitDetector(double f2Theta = -999.0,double fDistance = 0.0,double fBeam0 = -99999.0,double fBeam1 = -99999.0);
    double  f2dtothetaToReso(int nPix0,int nPix1);
    int     nSimplexRun(double fValue,double* afDelta,double* afStep);
    int     nCompress(Cimage& oOld,int nComp);

    float   fRmerge(int nCenter0,int nCenter1,bool bUseLineRMerge);
    int     nFindLowestRmerge(int& nValue0,int& nValue1);

    Cbeamdata(Cimage& _oCi,int nCompression = 1);
    ~Cbeamdata();

private:
    bool bWriteResoRingsToHeader(itr<double>& afResoRingMid, 
                                 itr<double>& afResoRingStrength, 
                                 double dMaxExpected2ThetaRingWidth_DEG,
                                 Cimage_header* poHeader);

    bool bUpdateInternalHeader(Cimage_header* poHeader);
};
#endif   // DT_CBEAMDATA_H//
