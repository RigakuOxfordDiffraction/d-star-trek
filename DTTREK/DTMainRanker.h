//
// Copyright (c) 2004 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainRanker.h     Initial author:                RB 15-Mar-2004

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

#ifndef DT_MAIN_RANKER_H
#define DT_MAIN_RANKER_H

#include "DTMain.h"

class DTREK_EXPORT CDTMainRanker : public CDTMain
{
public:
    CDTMainRanker();
    ~CDTMainRanker();
    virtual int nExecute(unsigned int argc, char* argv[]);

private:
    virtual void vError(const int nErrorNum, const Cstring& sMessage);
    
    Cstring             m_strImageFilePath;
    Cstring             m_strCurrentHeaderFilePath;
    Cstring             m_strInputHeaderFilePath;

    Cstring             m_strResolutionLimits;
    Cstring             m_strResolutionBins;

    int                 m_nTotalNumberOfResoBins; 
    int                 m_nNumberOfHighResoBins;  
    int                 m_nNumberOfLowResoBins;

    Cstring             m_strTemplate;

    Cstring             m_strBeamCenter;
    Cstring             m_strDetectorDistance;
    Cstring             m_strTwoTheta;
    Cstring             m_strWavelength;
    
    Cstring             m_strMinPeakCount;
    Cstring             m_strSigmaFIND;
    Cstring             m_strSigmaREFINE;
    
    Cstring             m_strAnisotropy;
    Cstring             m_strStrategy;
    
    Cstring             m_strFindOptions;
    //Cstring             m_strRankFind;
    Cstring             m_strIndexOptions;
    Cstring             m_strRefineOptions;
    Cstring             m_strRankOptions;

    Cstring             m_strExcludeRankRules;

    Cstring             m_strResultsFileName;
    int                 m_nNumberOfFilePathChars;

    Cstring             m_strSampleNameForReport;
    bool                m_bRulesHelpOnly;

    std::vector< std::pair<int, int> >       m_vecScanImageSequenceRanges;

    bool                m_bUseNonunf;
    Cstring             m_strETE;

    typedef enum
    {
        enProcNone=-1,
        enProcExtractHeader,
        enProcFind,
        enProcIndex,
        enProcRefine,
        enProcRank
    }enumProcessingStep;
    
    bool bCommandArgumentArrayToRankerCommandLineParameters(unsigned int argc, char* argv[]);
    bool bCreateImageHeaderFile();
    
    bool bGetImageTemplate();

    bool bDeriveImageSequence();

    void vGetImageSequenceInformationFromDTREKCommandLine(enumProcessingStep eStep);
    bool bUpdateHeaderFileAccordingToRankerCommandLine();
    enumProcessingStep enDoImageProcessing();
    bool bExecuteDTREKCommand(CDTMainRanker::enumProcessingStep enProcStep, Cstring& strOptions);
    bool bDoFind();
    bool bDoIndex();
    bool bDoRefine();
    bool bDoRank();
    bool bOutputResults(int nNumberOfRankingRules);
    void vSetResoBinsParams();

    int  nGetNumberOfImagesInSeqRanges();
    bool bGenerateSampleNameForReport();
    bool bParseResultsFileOption();
};
#endif//!DT_MAIN_RANKER_H
