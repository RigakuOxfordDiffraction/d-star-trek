//
// Copyright (c) 2000 Molecular Structure Corporation
//                    9009 New Trails Drive
//                    The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved.
//
// Ccompare.h     Initial author: T.J. Niemeyer  Spring 2000
// Header file for Ccompare.cc and dtcompare.cc
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

#ifndef DTCOMPARE_HEADER
#define DTCOMPARE_HEADER

#include "Dtrek.h"
#include "Creflnlist.h"
#include "Ccrystal.h"
#include <string.h>

enum eComparisonType {APPLES_APPLES,APPLES_ORANGES,APPLES_MAC};

class DTREK_EXPORT Ccompare {
    Creflnlist& m_oList1;
    Creflnlist& m_oList2;
    Ccrystal m_oCrystal;
    int m_nStat;

    Creflnlist* m_poList[2];
    int* m_pnReject[2];
    int* m_pnIndex[2];
    int* m_pnIndexU[2];
    int* m_pnIndexRMerge[2];
    int  m_nNumRefs[2];
    int  m_nFI_nIndex[2];
    int  m_nFI_nIndexU[2];
    int  m_nFI_fIndexRMerge[2];

    double* m_pfIntensityBinsNumer[2];
    double* m_pfIntensityBinsDenom[2];
    int*    m_pnIntensityBinsContrib[2];
    int     m_nIntensityBins;
    double  m_fIntensityMax;

    bool    m_bResortRMergeArray[2];
public:
    bool    m_bSilent;
    double  m_a2fRMerge[2];
    bool bIsAvailable() { return (m_nStat==0); };


    int nEquivocate(eComparisonType eType,bool bEqualRedundancy,int nMinRedundancy,int nRejectFlag = 1,int nUseSigmasIn = -1);
    int nCount(int nDataSet,int nRejectFlag);
    int nReject(double fPercent,int nDataSet);
    int nWriteNonRejected(int nDataSet,Cstring& sName);
    int nComputeRMerge(int nMinRedundancy,bool bDupSigmas);
    int nClearRejects();
    int nInitIntensityBins(double fPercent,int nNumBins);
    int nWriteIntensityBins(Cstring& sFileOut);
    
    int nWriteList(int nList,Cstring& sNewName);
    Ccompare(Creflnlist& oRef1,Creflnlist& oRef2,Cimage_header& oHeader);
    ~Ccompare();
};

#endif
