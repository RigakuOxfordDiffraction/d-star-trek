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
// DTMainRank.h     Initial author:       RB 15-Mar-2004
// Wrapper-class around former dtrank entry point

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

#ifndef DT_MAIN_RANK_H
#define DT_MAIN_RANK_H

#include "DTMain.h"
#include "Cimage_header.h"

typedef struct _tagRESOLUTION_SHELL_INFO
{
    int         m_nReflnCount;
    double      m_dResoStart_A;
    double      m_dResoEnd_A;
    double      m_dIOverSig;
    double      m_dIOverSigDev;

    void        vInitialize(){m_nReflnCount=0; m_dResoStart_A=0.0; m_dResoEnd_A=0.0; m_dIOverSig=0.0; m_dIOverSigDev=0.0;}

    bool        bSet(Cstring  sInfo)
    {
        std::vector<double>     adTemp;

        if( 0 != sInfo.nListToVector(adTemp, " ") )  // expected string like: "3.51 2.78 201 111.5 5.0"
            return false;

        int     nNumberOfTokens = adTemp.size();
    
        if( nNumberOfTokens < 5 )
            return false;

        m_dResoStart_A = adTemp[0];
        m_dResoEnd_A   = adTemp[1];

        m_nReflnCount  = (int)adTemp[2];

        m_dIOverSig = adTemp[3];

        m_dIOverSigDev = adTemp[4];
        
        return true;
    }
}RESOLUTION_SHELL_INFO;


class DTREK_EXPORT CDTMainRank : public CDTMain
{
public:
    CDTMainRank();
    virtual int nExecute(unsigned int argc, char* argv[]);
    void vGetRulesDescription(Cstring& strDesc);
    void vSetSampleNameForReport(Cstring& strSampleNameForReport){m_strSampleNameForReport=strSampleNameForReport;}

    
private:
    virtual void vError(const int nErrorNum, const Cstring& sMessage);
    
    void vRankResolutionShellInfo(Cstring& strDTREKModuleName,
                                  int    nTotalBins,
                                  int    nHighestBins,
                                  int    nLowestBins,
                                  
                                  int    nMinAcceptableReflnCountPerShell,
                                  int    nPointsAwardForReflnCountPerShell,
                                  int    nMinAcceptableReflnCountPerShellForIOverSigmaCalc,
                                  int    nPointsAwardForIOverSigPerShell,
                                  int    nPointsAwardForIOverSigRotDevPerShell,
                                  double dIOverSigScaleCoefficient,
                                  int    nPointsPenaltyMaxReso);
    
    void vPrintDtRankDescription(FILE* pFile=NULL);
    
    
#define DTREK_RANK_START            0x0001
#define DTREK_RANK_TALLY            0x0002
#define DTREK_RANK_PRINT_SUBSCORES  0x0004
    void vRankToOutput(int nRuleNumber,
                       char* pcDescription=NULL, 
                       int nRulePoints=0,
                       int nRuleMaxPoints=0,
                       DTREK_WORD wCtrl=0U);
    
    void vRankToHeader(DTREK_WORD wCtrl=0U);

    void vRankPeakSharpness();
    void vRankStrategyResults();

    int   m_nCumulativePointCount;
    int   m_nActiveRuleNumber;
    int   m_nActiveRulePoints;

    Cstring m_strActiveRuleSubpoints;

    bool  m_bRefineRankingInfoMustBeInHeader;
    Cstring m_strSampleNameForReport;
};
#endif// !DT_MAIN_RANK_H
