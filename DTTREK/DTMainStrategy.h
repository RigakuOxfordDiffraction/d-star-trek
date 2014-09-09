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
// DTMainStrategy.h    Initial author: RB          21-Mar-2005
// This file contains the definitions of class CDTMainStrategy

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

#ifndef DT_MAIN_STRATEGY_H
#define DT_MAIN_STRATEGY_H

#include "DTMain.h"
#include "Cstrategy.h"
#include "ResoBins.h"
class Cimage_header;
class Creflnlist;
class Cpredict;

class DTREK_EXPORT CDTMainStrategy : public CDTMain
{
public:
    CDTMainStrategy();
    ~CDTMainStrategy();
    
    void vInit();

    virtual int nExecute(unsigned int argc, char* argv[]);
    void vPrintOptionHelp();
    
    double  dGetTotalRotWidth();
    
    void vSetExternalPredictReflnListPtr(Creflnlist* pList);
    
    void vSetExternalResoBinsPtr(CResoStrategyBins* pBins);

    void vSetExternalScanInfoPtr(CstrategyScans* pScans);

    void vDoExternalResoBinsStatistics(CResoStrategyBins* poResoBins, Creflnlist* poPredictReflnList);
    void vDoExternalReflnListCalcSetDStarCubed(Creflnlist* poReflnList, double& dStarCubedMin, double& dStarCubedMax);
    void vDoExternalReflnListSortReduce(Creflnlist* poReflnList);

    double  fGetLengthScaleFactor()const{return m_fLengthScaleFactor;}

    void vSetLengthScaleFactor(const double dScaleFactor){m_fLengthScaleFactor = dScaleFactor;}

    void vSetVerbose(DTREK_WORD wVerbose){m_wVerbose = wVerbose;}

    void vSetUserInputMosaicity(double fMos){m_fUserInputMosaicity=fMos;}

    double fGetBestCompleteness();
    double fGetCompletenessTolerance();

private:
    virtual void vError(const int nErrorNum, const Cstring& sMessage);
    
    void vCreateStrategyObject();
    
    void vDeleteStrategyObject();
    
    void vCreatePredictReflnList();
    void vExpandPredictReflnList();
    
    void vDeletePredictReflnList();

    void vDeleteResoBins();
    
    Cpredict* pSetupPredictObjectForScan(int nScan,
                                         double f2ThetaSwing,
                                         double fDistance,
                                         double* a2fResolution);

    bool m_bExternalPredictReflnListPtr;
    bool m_bExternalResoBinsPtr;

    Cstrategy*      m_poStrategy;
    
    CstrategyScans*  m_poExternScanInfo;

    Creflnlist*     m_poPredictReflnList;
    
    CResoStrategyBins*      m_poResoBins;

    double          m_fLengthScaleFactor;
    
    bool            m_bAnom;

#define     DTREK_DTMAINSTRAT_VERBOSE_PREDICT   0x0001
    DTREK_WORD      m_wVerbose;
    
    double          m_fUserInputMosaicity;
};
#endif   // !DT_MAIN_STRATEGY_H
///////////////////////////////////

