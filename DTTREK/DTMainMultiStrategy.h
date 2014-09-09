//
// Copyright (c) 2005 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainMultiStrategy.h    Initial author: RB          21-Mar-2005
// This file contains the definitions of class CDTMainMultiStrategy

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

#ifndef DT_MAIN_MULTI_STRATEGY_H
#define DT_MAIN_MULTI_STRATEGY_H

#include "DTMainStrategy.h"
#include "StrategyScans.h"
#include "Segment.h"
#include "ResoBins.h"
#include "CExperimentalStrategy.h"

class Creflnlist;

class DTREK_EXPORT CDTMainMultiStrategy : public CDTMain
{
public:
    CDTMainMultiStrategy();
    ~CDTMainMultiStrategy();
    virtual int nExecute(unsigned int argc, char* argv[]);

private:
    virtual void vError(const int nErrorNum, const Cstring& sMessage);
    
    bool bSetAuto2Theta(Cstring& sError);
    bool bSetUserInput2Theta(Cstring& sError);
    
    bool bSetAutoDistance(Cstring& sError);
    bool bSetUserInputDistance(Cstring& sError);
    
    bool bGenerate2ThetaPositions(Cstring& sError);
    bool bGenerateValidDistance(Cstring& sError);
    
    bool bSetResolutionLimitsFromCommandLine(Cstring& sResoOption, Cstring& sError);

    bool bSetResolutionLimitsFrom2ThetaPositions(Cstring& sResoOption, Cstring& sError);
    
    void vPrintSummaryScanTable();

    void vUpdateHeader();

    bool bGenerateResoBins(double dResoMin, double dResoMax, int nNumberOfBins);
    void vExpandReflnList(Creflnlist* poList);

    void vSetExecutableNameDTSTRATEGY();

    void vRedistributeDetectorOverlaps();
    bool bLoadPreviousReflections(Cstring& sError);

    bool bParseReflnOverlapCheckParams(Cstring& sError);
    bool bParseMosaicityInput(Cstring& sError);

#define DTREK_MULTISTRATEGY_DET_DISTANCE       0x0001
#define DTREK_MULTISTRATEGY_DET_2THETA         0x0002
    bool bReadSetHardwareLimits(Cstring& sError, DTREK_WORD wCtrl);
    
    void vSetDetectorPositionToHeader(DTREK_WORD wCtrl);
    void vGetDetectorPositionFromHeader(DTREK_WORD wCtrl);
    
    void vRemoveMultiScanInfoFromHeader();
    void vClearMultiStrategyOutputFromHeader();
    
    bool bSetTestScansFromCommandLineToHeader(Cstring& sError);

    bool bInsureTestScansInHeaderAreInNewFormat(Cstring& sError);
    bool bApplyRestrictedRotationLimitsToScansInHeader(Cstring& sError);
    
    void vSaveTestScansFromHeader();

    bool bCheckReflectionOverlap(Cstring& sError);

    bool bGetResoRangeForCurrentDetectorPosition(double& dResoMin, double& dResoMax, Cstring& sError);

    bool bInsureHeaderRotLimitsConsistency(Cstring& sError);

    void vUpdateTargetCompletenessTolerance(double fCompleteness, double fTolerance);
    void vSaveTargetCompletenessTolerance();
    bool bDoExternalResoBinsStatistics(int nNumberOfBins, Creflnlist* pReflnList, Cstring& sError);
    bool bGetReasonableCrystalGonioHardwareLimits(Cstring& sAxisName,
                                                  float&   fAxisHardwareMin,
                                                  float&   fAxisHardwareMax,
                                                  bool     bOptimize,
                                                  Cstring& sError);
    bool bGetImageRotationWidth();

    bool bCalculateAutomaticExposureTime();

private:    //shijie yao 03:26:08 -- helper function
    void vUpdateStrategyEntryWithImageRotWidth(Cstring EntryStr, float imageWidth);

private:
    Cstring         m_sError;
    CDTMainStrategy m_oMainStrategy;
    Creflnlist*     m_poPredictReflnList;
    Creflnlist*     m_poPreviousReflnList;

    CResoStrategyBins*      m_poResoBins;

    double          m_dUserResoMin;
    double          m_dUserResoMax;

    double          m_dDistance;
    
    int             m_nUserSpotSepPix;
    CSegment        m_segUserDistLimits;
    
    
    std::vector<double>          m_ad2ThetaPositions;
    double          m_d2ThetaFromHeader;
    
    int             m_nActive2ThetaPositionIndex; // the one currently used by the algorithm

    std::vector<CstrategyScans*>   m_aoScanSolutionSet; // to keep info on all strategy scans, generated by strategy
    
    CstrategyScans*                m_paoSaveTestScans;  // to save the test scans user inputted

    double          m_dFractionDetectorOverlap;
    
    CSegment        m_segUser2ThetaLimits;
    
    int             m_n2ThetaDirectionSign;
    
    double          m_dUser2ThetaMust;

    double          m_dAllDetectorPositionsAllScansTotalRotWidth;

    bool            m_bScaleSetFromCommandLine;

    bool            m_bResolutionSet;

    STRATEGY_REFLN_OVERLAP_CHECK*    m_pstOverlapCheckParams;
    SPOT_SIZE_INFO                   m_stSpotSizeInfo;

    double          m_fUserInputMosaicity;

    double          m_fUserTargetCompleteness;
    double          m_fUserCompletnessTolerance;

    bool            m_bWriteOverlapReflnList;

    // Since the earliest strategy versions both test and result scans were described by the same keyword prefix SCAN_ (S*_).
    // Beginning with d*TREK 9.6 we are distingushing TEST and RESULT keywords.
    bool            m_bUseOldHeaderScanFormat;
    bool            m_bTestScansGeneratedFromCommandLine;
    bool            m_bTestScansGeneratedFromOldHeaderFormat;

    REFLN_FILE_STRAT_INFO*          m_pstPreviousReflnsInfo;

    double          m_fUserImageRotationWidth;
    double          m_fImageAutoExposureTime;
};
#endif   // !DT_MAIN_MULTI_STRATEGY_H
