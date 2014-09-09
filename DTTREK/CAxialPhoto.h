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
// CAxialPhoto.h     Thaddeus J. Niemeyer           17-Apr-2002
//    Computes an axial photograph from images.
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
#ifndef CAXIALPHOTO_HEADER 
#define CAXIALPHOTO_HEADER


#include "Dtrek.h"
#include "dtreksys.h"
#include "Cfind.h"
#include "Crotation.h"
#include "dtarray.h"

class DTREK_EXPORT CAxialPhoto {
    Cimage_header m_oHeader;
    Ccrystal    m_oCrystal;
    Cscan       m_oScan;
    Cdetector*  m_poDetector;
    Csource     m_oSource;
    Cgoniometer m_oGonioIn;
    double      m_a3x3fAxialGonio[3][3];
    double      m_a3fAxialRotVec[3];
    double      m_a360fDegCounts[360];
    double      m_a360fCumulDegCounts[360];

    Cimage*     m_poBackground;     // Contains cumulative background pixels. (ui2)
    Cimage*     m_poPhotons;        // Photon counts for a particular image (ieee)
    Cimage*     m_poAxialPhotons;   // Cumulative photons for completed axial photo. (ieee)
    Cimage*     m_poInput;          // Input image read in here (ui2 or whatever else is used)

    int         m_nStat;
    int         m_nAxisCode;
    itr<int>    m_anSeqNum;
    double      m_fSigma;
    double      m_a2fReso[2];
    int         m_nDim0,m_nDim1;
    bool        m_bPixelAdjust;
    double      m_a2fRotRange[2];
    Cstring     m_sOutput;
    int         m_nHKLRestrict;
    int         m_nHKLValue;
    double      m_fHKLWidth;

public:
    void vSetResolution(double fResoHigh,double fResoLow);
    void vAddSeq(int nSeqStart,int nSeqEnd);
    void vSetSigma(double fSigma);
    void vSetAxisCode(int nAxisCode);
    void vSetRotRange(double fRotStart,double fRotEnd);
    void vSetOutputImageName(Cstring& sName) { m_sOutput = sName; };
    void vSetPixelAdjust(bool bAdjust) { m_bPixelAdjust = bAdjust; };
    void vSetZone(int nHKL,int nValue);
    void vSetZoneWidth(double fHKLWidth);

    void vInit();
    bool bIsAvailable() { return !m_nStat; };

    void vInitDegCounts(double* pfDegCounts);
    void vPrintDegCounts(double* pfDegCounts);

    int nRun();
    int nAddBackground(Cimage& oImage);
    int nFindPeaks(Cimage& oImage);
    int nAddPeaks(double fRot);
    int nSetAxialGonio();
    int nBuildFinalImage();

    CAxialPhoto(Cimage_header& oHeader);
    ~CAxialPhoto();

};

#endif
