//
// Copyright (c) 2001 Molecular Structure Corporation
//
// CinterMesh.h        Initial author: T.J. Niemeyer  Jan-2001
//  This file contains the prototype of class CinterpMesh.h
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

#ifndef CINTERPMESH_HEADER
#define CINTERPMESH_HEADER
#include "Dtrek.h"
#include "Cstring.h"
#include "dtarray.h"

const int g_nMaxInterpLevel = 20;
typedef int (*nInterpCallback)(int nDataPoints,double* pfInData,double* pfOutData);


class DTREK_EXPORT CinterpMesh {

	itr<double> m_afOrigItems;			// Original items if we are saving those. (size depends on # of items added)
	itr<int>	m_anNextItem;			// Next item pointer. (size depends on # of items added)
	itr<int>	m_anOrigItem;			// Index into m_afOrigItems[]

    itr<double> m_afTempBuf;			// Set to two times the width for quick calculations in pfCalcGetItem().
	double*		m_pfTempBuf;			// Points to m_afTempBuf[0];

	double* m_pfItems;					// Averaged items. (NULL if we are saving originals).
    int* m_pnCounts;					// Counts for averaged items.
	int* m_pnFirst;						// Start index into m_anNextItem. (NULL if we are NOT saving originals).
    
    static int ms_nLevelPointers[g_nMaxInterpLevel+2];
    static double ms_fLevelGridStep[g_nMaxInterpLevel+2];
    static int ms_nLevelGridStep[g_nMaxInterpLevel+2];
    static bool ms_bInitialized;
    
    double m_fDim0Range[3];
    double m_fDim1Range[3];
    int m_nNumContrib;
    int m_nMinimumToInterp;
    int m_nWidth;
    int m_nMaxLevel;
	bool m_bSaveOriginals;
	double m_fSigmaReject;				// Only non-zero if m_bSaveOriginals is true.



public:
    // These are set with each call to nGet()
    itr<int> m_anInterpOffset;
    itr<int> m_anInterpCounts;

    void vClear();
    void vInit(int nMaxLevel,int nWidth,double fDim0Min,double fDim0Max,double fDim1Min,double fDim1Max);
    void vSetInterpConditions(int nMinimumToInterp,double fSigmaReject = 0.0);
    int nGetNumContrib() { return m_nNumContrib; };
    int nAdd(double fDim0,double fDim1,double* fItems,nInterpCallback pxCallback = NULL);
    int nGet(double fDim0,double fDim1,double* fItems,double* fItemsSigma,int nLevelToUse = -1);
    int nTest(Cstring& sTestImage);
    int nPrint(Cstring& sTestImage,int nDim0,int nDim1);
    
    double* pfGetItem(int nOffset) { return &m_pfItems[nOffset*m_nWidth]; };	// WARNING:  Only works if we were not saving original data
	double* pfCalcGetItem(int nOffset,double* pfBuffer,double* pfVarBuffer = NULL);
    double fGetAverageLevel();
    int nGetMemorySize() { return ms_nLevelPointers[m_nMaxLevel+1]; };
    double fGetRangeMax(int nVar) { return ((nVar==0)?(m_fDim0Range[1]):(m_fDim1Range[1])); };
    double fGetRangeMin(int nVar) { return ((nVar==0)?(m_fDim0Range[0]):(m_fDim1Range[0])); };

    CinterpMesh();
    ~CinterpMesh();
};

#include "CSphericalHarmonic.h"

class DTREK_EXPORT CSphericalMesh {
    CSphericalHarmonic  m_oHarmonics;

    int                 m_nDim0Azi;
    double              m_fStep0Azi;
    int                 m_nDim1Polar;
    double              m_fStep1Polar;

    int                 m_nNumContrib;
    int                 m_nNumContribAtLastRecompute;
    int                 m_nWidth;
    
    itr<double>         m_afOrigItems;          // User provided points.
    itr<double>         m_afOrigItemsSigma;     // User provided points.
    itr<double>         m_afPolar;              // User provided points.
    itr<double>         m_afAzi;                // User provided points.

    itr<double>         m_afMeshValues;         // Calculated points.

    // Temporary variables used to make things go a little bit faster.
    itr<double>         m_afValues;

    
public:
    int                 nRecompute();

    void vInit(int nHarmonicOrder,int nWidth,double fMaxPolar,double fPolarStep,double fAziStep);
    int nAdd(double fPolar,double fAzi,double fSigma,double* pfItems);
    int nGet(double fPolar,double fAzi,double* pfItems);
    int nGetNumContrib() { return m_afOrigItems.size(); };
    int nGetNumContribSinceLastCompute() { return m_nNumContrib - m_nNumContribAtLastRecompute; };
    bool bReadyToAdd() { return ((nGetNumContrib()==0) ||(nGetNumContribSinceLastCompute()>=1)); };
    void vClear();
    CSphericalMesh();
};

#endif

