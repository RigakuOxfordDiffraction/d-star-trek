///////////////////////////////////////////////////////////////////////////
//
// ImageExpTimeEvaluator.h: interface for the CImageExposureTimeEvaluator class.
//
///////////////////////////////////////////////////////////////////////////

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
 
#ifndef CC_IETE_H
#define CC_IETE_H

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "Creflnlist.h"

class DTREK_EXPORT CImageExposureTimeEvaluator  
{
public:
	CImageExposureTimeEvaluator();
	virtual ~CImageExposureTimeEvaluator();

	bool ProcessImage(const char* const strImageName, bool bCleanPeakList=false);

	bool bGetSaturatedPeaksPercent(double dExposureTime_sec, double& dSaturatedPeaksPercent);
    bool bGetExposureTimeWithSaturationLimit(double dSaturatedPeaksPercent, double& dExposureTime_sec, bool& bZeroPeaksSaturated);
    
    bool bSuggestImageExposureTime(const char* const strImageName, double dMaxSaturatedPeak_Percent, double& dRecommendedExposureTime_sec);
    
    bool bGetExposureTimeWithTargetSignalToNoise(double& dObservedSignalToNoise,
                                                 int& nObservedPeakCount,
                                                 double& dSignalToNoiseOptimizedExposureTime_sec);

	bool bGetExposureTimeWithTargetSignalToNoise(double& dTargetResShell_Begin_A,
                                                 double& dTargetResShell_End_A,
                                                 double& dObservedSignalToNoiseInTargetShell,
                                                 int&    nObservedPeakCountInTargetShell,
                                                 double& dSignalToNoiseOptimizedExposureTime_sec);

	double GetImageSaturatedPeaksPercent()const{return m_dImageSaturatedPeaksPercent;}
	double GetImageExposureTimeSec()const{return m_dImageExposureTime_sec;}
	bool IsCCDImage(){return m_bCCDImage;}
    
    DTREK_WORD  wGetStatus()const{return m_wStatus;}
    bool        bIsError(){return (0U != m_wStatus);}
    
    static Cstring ms_sError;
    static Cstring ms_sMaxPixelValue;
    static Cstring ms_sSignalToNoiseRatio;

    void vStartImageTestSequence(double  dTargetResolution_A,
                                 double  dTargetSignalToNoiseInTargetShell,
                                 bool bResetMaxObservedCountInTargetShell=true);
    
    int nGetTotalNumberOfReflns(){return m_oPeaks.nGetNumReflns();}

    void vSetMaxObservedPeakCountInTargetShell(int nMaxReflnCount){m_nMaxObservedPeakCountInTargetShell=nMaxReflnCount;} 
    int nGetMaxObservedPeakCountInTargetShell(){return m_nMaxObservedPeakCountInTargetShell;}

    void vSetPixelSatValue(double fSatValue){m_dImageSaturatedPixelValue = fSatValue;}
    double fGetPixelSatValue(){return m_dImageSaturatedPixelValue;}

protected:
private:

#define     DTREK_IETE_INIT_ALL     0x0001
    void Initialize(DTREK_WORD wCtrl=0U);

	void CleanPeakList();
    void vSortPeakList(eReflnFieldType eType, int nField);
    
    bool GetTargetResolutionShellLimits(double dTargetResolution_A, 
                                        double& dTargetResShell_Begin_A, 
                                        double& dTargetResShell_End_A);

    bool bSetResolutionShellLimits();
    bool bGetSaveImageResolution(Cimage_header& oHeader);

    int     m_nSortKeyField;
    int*    m_pnSortIndex;

    int     m_nMaxPixelValueFieldIndex;

	double	m_dImageSaturatedPixelValue;  
	int		m_nPedestal;					// "dark current" in CCDs

	double	m_dImageSaturatedPeaksPercent;
	double	m_dImageExposureTime_sec;      
	
	Creflnlist      m_oPeaks;

	bool	m_bCCDImage;

    std::vector<double>          m_daResolutionShellLimits;

    // Target related members
    int         m_nMaxObservedPeakCountInTargetShell;

    float       m_fImageResMin;
    float       m_fImageResMax;

    double      m_dTargetResolution_A;
    double      m_dTargetSignalToNoiseInTargetShell;

    double      dGetTimePowerCoefficient();

#define DTREK_IETE_IMAGE_NOT_PROCESSED              0x0001
#define DTREK_IETE_INVALID_IMAGE_HEADER             0x0002
#define DTREK_IETE_PEAK_SEARCH_FAILURE              0x0004
#define DTREK_IETE_INVALID_RESOLUTION               0x0008
#define DTREK_IETE_RESOLUTION_BINNING_FAILURE       0x0020
#define DTREK_IETE_TOO_FEW_PEAKS_IN_TARGET_SHELL    0x0040

    DTREK_WORD  m_wStatus;
};
//////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CImageExposureTimeCalculator 
{
public:
    CImageExposureTimeCalculator(){m_sError = "";}
public:
#define     DTREK_IETC_FULL_MERGED_REFLNS       0x0001
bool    bCalculate(Cimage_header& oHeader,
                   double   fMinReso,
                   double   fMaxReso,
                   int      nResoBins,
                   double   fImageRotWidth,
                   double   fTargetIoverSigma,
                   double&  fProposedExposureTime,
                   DTREK_WORD wCtrl);
    
    Cstring     sGetLastError();
private:
    bool bValidateResolutionInput(double& fResoMin, double& fResoMax);
private:
    Cstring     m_sError;
};
#endif // !CC_IETE_H
