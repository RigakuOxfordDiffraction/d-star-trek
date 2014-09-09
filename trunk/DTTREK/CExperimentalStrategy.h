//
// Copyright (c) 2000 Molecular Structure Corporation
//
// dtprofit.cc     Initial author: T.J. Niemeyer         Jul-2003
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

#ifndef EXPERIMENTAL_STRATEGY_INCLUDE
#define EXPERIMENTAL_STRATEGY_INCLUDE

#include "Cimage_header.h"
#include "Cdetector.h"
#include "Csource.h"

typedef struct _tagSpotSizeInfo
{
    void vInit()
    {
        fMaxSpotSize        = -1.0;    
        fAverageSpotSize    = -1.0;
        fMinSpotSize        = -1.0;    
    }

    _tagSpotSizeInfo& operator=(const _tagSpotSizeInfo& anotherInfo)
    {
        fMaxSpotSize = anotherInfo.fMaxSpotSize;
        fAverageSpotSize = anotherInfo.fAverageSpotSize;
        fMinSpotSize = anotherInfo.fMinSpotSize;

        return *this;
    }
    
    void vToArray(double fSpotSize[3])
    {
        fSpotSize[0] = fMaxSpotSize;
        fSpotSize[1] = fAverageSpotSize;
        fSpotSize[2] = fMinSpotSize;
    }

    void vFromArray(double fSpotSize[3])
    {
        fMaxSpotSize     = fSpotSize[0];  
        fAverageSpotSize = fSpotSize[1];  
        fMinSpotSize     = fSpotSize[2];  
    }

    double  fMaxSpotSize;
    double  fAverageSpotSize;
    double  fMinSpotSize;
}SPOT_SIZE_INFO;

class DTREK_EXPORT CExperimentalStrategy
{

    Cimage_header&  m_oHeader;

    double          m_fDebyeRingCoverage;
    int             m_nMethod;

public:

    DTREK_WIN_DLL_DATA_EXPORT static char* ms_pcHelp;
    
    static double   ms_fDefaultSpotSize[3];

    enum            
    {
        eResoOuter,
        eResoInner,
        eResoChiCoverage
    };
    
    enum
    {
        eComputeDistanceSpotSeparationDefault = -1,
        eComputeDistanceSpotSeparationWorstCase,
        eComputeDistanceSpotSeparationAverageCase,
        eComputeDistanceSpotSeparationBestCase
    };
    
    Cdetector       m_oDetector;
    Csource         m_oSource;
    Cstring         m_sOutputGraph;
    bool            m_bPrint;

    int nCalcReso(double& fResoHigh,double& fResoLow,double f2Theta = -999.0,double fDistance = -999.0);
    
    double fCalcBest2ThetaSwingForTargetHighReso(double f2ThetaMin,double f2ThetaMax,double fTargetHighReso, double& fBestHighReso);
    double fCalcBest2ThetaSwingForTargetLowReso(double f2ThetaMin,double f2ThetaMax,double fTargetLowReso, double& fBestLowReso);

    double fCalcBestDistance(double fDistanceMin,double fDistanceMax,double fTargetReso);
    
    double fCalcBestDistanceSpotSeparation(double fPixelSep,
                                           int nModeAxis=eComputeDistanceSpotSeparationDefault,
                                           int nModeSpot=eComputeDistanceSpotSeparationDefault,
                                           SPOT_SIZE_INFO* pstSpotSizeInfo=NULL);
    
    bool bGetSpotSizeFromHeader(SPOT_SIZE_INFO& pstSpotSizeInfo);  

    double fConvertToDegrees(double fReso);
    int    nCalcBest2ThetaSwings(double fInitSwing,double fDegreesOverlaping,itr<double>& af2Thetas);
    
    int nBuildResoShells(double fResoLow,double fResoHigh,int nNumShells,double* pfShellStart);
    int nParseResoCommandLine(int argc,char** argv,double& fResoHigh,double& fResoLow,double& f2ThetaSwing,double& fDistance);        

    bool bIsAvailable();
    void vSetResoMethod(int nMethod,double fDebyeRingCoverage = 30.0) { m_nMethod = nMethod,m_fDebyeRingCoverage = fDebyeRingCoverage; };
    void vUpdateHeader();

    CExperimentalStrategy(Cimage_header& oHeader);
};

#endif
