//
// ImageExpTimeEvaluator.cpp : implementation file
//

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
 
#ifndef VC9
	#ifdef SSI_PC
		#include "stdafx.h"
		#define SSI_PC_ASSERT   ASSERT
	#else
		#define  SSI_PC_ASSERT
	#endif
#else
	#define SSI_PC_ASSERT 
#endif

//#define     IETE_DEBUG_OUTPUT

#include "ImageExpTimeEvaluator.h"
#include "Crefln.h"
#include "Cfind.h"
#include "Crotation.h"

Cstring CImageExposureTimeEvaluator::ms_sError              = "nError";       
Cstring CImageExposureTimeEvaluator::ms_sMaxPixelValue      = "nMaxPixelValue";
Cstring CImageExposureTimeEvaluator::ms_sSignalToNoiseRatio = "fSignalToNoiseRatio";

const int       c_nMinRequiredSpotCountInTargetShell = 10; // so we could do statistics on that shell
const double    c_dMinimumAcceptableSignalToNoisePerPeak = 3.0;
const double    c_dMaximumAcceptablePeakLossCountInTargetShell_percent = 10.0;

// I/Sigma as a function of time changes proportionally to time ^ (1 / Power)
const double    c_dTimePower_IP  = 2.5;      
const double    c_dTimePower_CCD = 5.0;

CImageExposureTimeEvaluator::CImageExposureTimeEvaluator()
{    
    Initialize(DTREK_IETE_INIT_ALL);
}

void CImageExposureTimeEvaluator::Initialize(DTREK_WORD wCtrl)
{
	m_nPedestal = 0;

	m_dImageExposureTime_sec = 0.0;      
	
	m_bCCDImage = FALSE;

    m_nSortKeyField = -999;
    
    m_pnSortIndex = NULL;

    m_nMaxPixelValueFieldIndex = -1;
    
    m_daResolutionShellLimits.clear();
    
    m_wStatus = DTREK_IETE_IMAGE_NOT_PROCESSED;

    if( wCtrl == DTREK_IETE_INIT_ALL )
    {
	    m_dImageSaturatedPeaksPercent = 0.0;
	    
        CleanPeakList();
	    
        m_dImageSaturatedPixelValue = 0.0;  
        
        m_nMaxObservedPeakCountInTargetShell = 0;

        m_dTargetResolution_A = 1.8;
        m_dTargetSignalToNoiseInTargetShell = 10.0;
    }    
    
    m_fImageResMin = -1.0f;
    m_fImageResMax = -1.0f;
}

void CImageExposureTimeEvaluator::vStartImageTestSequence(double  dTargetResolution_A,
                                                          double  dTargetSignalToNoiseInTargetShell,
                                                          bool bResetMaxObservedCountInTargetShell)
{
    m_dTargetResolution_A = dTargetResolution_A;
    m_dTargetSignalToNoiseInTargetShell = dTargetSignalToNoiseInTargetShell;

    if( m_fImageResMax != -1.0f && m_fImageResMin != -1.0f )
        bSetResolutionShellLimits();

    if( bResetMaxObservedCountInTargetShell )
        m_nMaxObservedPeakCountInTargetShell = 0;
}


CImageExposureTimeEvaluator::~CImageExposureTimeEvaluator()
{
	CleanPeakList();
}

void CImageExposureTimeEvaluator::CleanPeakList()
{
	REFLN_EXTRA_INFO*	pstReflnExtraInfo = NULL;
	
	for(int nRef = 0; nRef < m_oPeaks.nGetNumReflns(); nRef++)
	{
		pstReflnExtraInfo = (REFLN_EXTRA_INFO*)(m_oPeaks.poGetRefln(nRef)->m_pvUserData);

		if( pstReflnExtraInfo )
		{
			delete 	pstReflnExtraInfo;
			(m_oPeaks.poGetRefln(nRef))->m_pvUserData = NULL;
		}
	}

	m_oPeaks.vDeleteAll();
}

/////////////////////////////////////////////////////////////////////////////////////////
// Read image header to determine image parameters: exposure time, saturated pixel value, etc.
// Run a peak search to set up a peak list.
// While setting up the peak list, determine the percentage of saturated reflections.
bool CImageExposureTimeEvaluator::ProcessImage(const char* strImageName, bool bCleanPeakList)
{
	Initialize();  // this initialization should not touch the target values !

    if( bCleanPeakList )
    {
        CleanPeakList();
        m_dImageSaturatedPeaksPercent = 0.0;
    }

    m_wStatus |= DTREK_IETE_IMAGE_NOT_PROCESSED;

	Cstring			    sFileName(strImageName);
    if( "" == sFileName )
	{
		m_wStatus |= DTREK_IETE_INVALID_IMAGE_HEADER;
		return false;
	}

	Cimage_header       oHeader(sFileName);
	if( !oHeader.bIsAvailable() ) // Rotation information not available
	{
		m_wStatus |= DTREK_IETE_INVALID_IMAGE_HEADER;
		return false;
	}
	
    if ( !bGetSaveImageResolution(oHeader) )
    {
		m_wStatus |= DTREK_IETE_RESOLUTION_BINNING_FAILURE;
        return false;
    }

	oHeader.nGetValue("DARK_PEDESTAL", &m_nPedestal);

	Cstring			sDetectorNames;			
	oHeader.nGetValue("DETECTOR_NAMES", &sDetectorNames);
	m_bCCDImage = sDetectorNames.contains("CCD");

	Crotation           oRot(oHeader);
       
	if( !oRot.bIsAvailable() ) // Rotation information not available
	{
		m_wStatus |= DTREK_IETE_INVALID_IMAGE_HEADER;
		return false;
	}

	m_dImageExposureTime_sec = oRot.fGetExposureTime();
	
	////////////////////////////////////////////////////////////////////////////////////////////
	Creflnlist      oListStart;	   // initial temporary list of reflections before getting their 2D center

	int			nFI_nError = 0;
	nFI_nError = m_oPeaks.nExpandGetField(ms_sError);    
    nFI_nError = oListStart.nExpandGetField(ms_sError);    
    
	Cfind           oFind(oHeader, &oListStart);
    oFind.m_bIncludeBackground = TRUE;
	
	oFind.nExpandReflnlist(&m_oPeaks,false,true);
    oFind.nExpandReflnlist(&oListStart,false,true);
    
    m_nMaxPixelValueFieldIndex = m_oPeaks.nExpandGetField(ms_sMaxPixelValue, eReflnField_float_type);
    ///////////////////////////////////////////////////////////////////////////////////////

    Cstring		sDetectorName("");    
    oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1, &sDetectorName);
    Cdetector	oDetector(oHeader, sDetectorName);

    Csource		oSource(oHeader); 
	Cimage		oImage(sFileName);
	
	//    if( nSetScanHeader(oHeader,false) )
	//          return FALSE;

	m_dImageSaturatedPixelValue = oImage.fGetSatValue();
	
	oFind.m_eMethod = eFind_2D_method;
	oFind.m_fSigma = 2.0;
	int		nStat = oFind.nFind2D(oImage, 0.0, &oDetector, &oSource);

	//          // Use the sensitive 2D search used in Cbackground
	//          // This tends to generate a lot of spurious spots, but we want this.
	//          Cbackground oBackground(oImage);
	//          oBackground.nBuildPeakPixels(max(50,max(oImage.m_nDim[0],oImage.m_nDim[1])/30));
	//          nStat = oFind.nExpandReflnlist(NULL,false,true);
	//          nStat += oBackground.nFind2D(oListStart,0.0,3.0,m_poDetector,m_poSource);
	//       
	//       
	//  
	if( nStat != 0 )
	{
        m_wStatus |= DTREK_IETE_PEAK_SEARCH_FAILURE;
		return false;
	}                                                                                                  	

	Crefln		oRefln(&m_oPeaks);
	
	REFLN_EXTRA_INFO*	pstReflnExtraInfo = NULL;
	
	for(int nRef = 0; nRef < oListStart.nGetNumReflns(); nRef++)
	{
		oRefln = oListStart[nRef];

		nStat = oFind.nFindCentroid2D(oImage, &oRefln, NULL, DTREK_DTFIND_FC2D_GET_EXTRA_INFO);
		if( 0 != nStat )
			continue;
		
		oRefln.vSetField(nFI_nError,nStat);
		
		pstReflnExtraInfo = (REFLN_EXTRA_INFO*)oRefln.m_pvUserData;
        SSI_PC_ASSERT(pstReflnExtraInfo);
		// Count saturated peaks
		if( pstReflnExtraInfo->m_nSaturatedCount > 0 && pstReflnExtraInfo->m_nPeakEllipsoidPixelCount > 1 ) // the requirement that the ellipsoid has at least one pixel is a sanity requirement
		{
			m_dImageSaturatedPeaksPercent += 1.0;
		
			// Roughly estimate peak max value based on saturated pixel value and number of saturated pixels vs
			// total number of peak pixels
            
            // For weird cases, when the number of pixels in the ellisoid is less or equal 
            // to the number of saturated pixels, re-set the number of saturated pixels to
            // be equal to the number of pixels in the ellipsoid minus one. 

            // 2DO: start using the number of pixels in peak mask as opposed to the number of pixels in ellipsoid

            if( pstReflnExtraInfo->m_nSaturatedCount >= pstReflnExtraInfo->m_nPeakEllipsoidPixelCount )
                pstReflnExtraInfo->m_nSaturatedCount = pstReflnExtraInfo->m_nPeakEllipsoidPixelCount - 1;  
			
            double		dA = sqrt((double)pstReflnExtraInfo->m_nPeakEllipsoidPixelCount);
			double		dB = sqrt((double)pstReflnExtraInfo->m_nSaturatedCount);

			pstReflnExtraInfo->m_dMaxPixelValue = m_dImageSaturatedPixelValue *  dA / ( dA - dB );
        
            float       ff = oRefln.fGetField(m_oPeaks.m_nFI_fObsPx0);
            float       gg = oRefln.fGetField(m_oPeaks.m_nFI_fObsPx1);
        }
		
		//oRefln.vSetField(m_oPeaks.m_nFI_fIntensity, 
        //			    (float)pstReflnExtraInfo->m_dMaxPixelValue);// A trick: save the max pixel intensity value
        //  														// as the total peak intensity, so we could use an
        //  														// existing sort function on the list later.
        //  														// Should be ok, because we are not using the
        //														    // total peak intensity anyway.
        
		oRefln.vSetField(m_nMaxPixelValueFieldIndex, 
        			    (float)pstReflnExtraInfo->m_dMaxPixelValue);
        
        m_oPeaks.nInsert(&oRefln);	// Save reflection in the member list
	}
	
	int			nFoundPeakCount = m_oPeaks.nGetNumReflns();
	
    printf("\nImage %s\n", strImageName);
    printf("\nCurrent total number of reflections examined for saturation on all images: %d\n", nFoundPeakCount);
    printf("\nCurrent total number of saturated reflections found on all images: %d\n", (int)m_dImageSaturatedPeaksPercent);

	if( nFoundPeakCount > 0 )
	{
        m_dImageSaturatedPeaksPercent /= nFoundPeakCount / 100.0; // get the percentage
        printf("\nCurrent total percentage of saturated reflections: %.1f\n", m_dImageSaturatedPeaksPercent);
	}
	else
	{
        SSI_PC_ASSERT( m_dImageSaturatedPeaksPercent==0.0 );
    }

	m_wStatus &= ~DTREK_IETE_IMAGE_NOT_PROCESSED;

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Given exposure time, evaluate the number of saturated peaks 
bool CImageExposureTimeEvaluator::bGetSaturatedPeaksPercent(double dExposureTime_sec, double& dSaturatedPeaksPercent)
{
	if ( 0.0 >= m_dImageExposureTime_sec )
		return false;						// cannot evaluate if no real measurement time is available

	if( m_oPeaks.nGetNumReflns() == 0 )
		return false; // list of peaks does not exist

    
    // Sort found peaks on maximum pixel value
    //vSortPeakList(eReflnField_float_type, m_oPeaks.m_nFI_fIntensity);
    vSortPeakList(eReflnField_float_type, m_nMaxPixelValueFieldIndex);
	
    double		dTimeScale = dExposureTime_sec / m_dImageExposureTime_sec;
	
	// Loop over all found peaks and count the number of peaks that will be saturated
	dSaturatedPeaksPercent = 0.0;
	
	Crefln*			poRefln = NULL;
	REFLN_EXTRA_INFO*	pstReflnExtraInfo = NULL;
	
	int			nSortedIndex = -1;
	double		dEvalMaxPixelValue = 0.0;
	double		dEvalMaxPixelValueError = 0.0;
	
	for(int nRef = m_oPeaks.nGetNumReflns()-1; nRef >= 0; nRef--)
	{
		nSortedIndex = m_pnSortIndex[nRef];

		poRefln = m_oPeaks.poGetRefln(nSortedIndex);

        SSI_PC_ASSERT(poRefln);
		
		pstReflnExtraInfo = (REFLN_EXTRA_INFO*)(poRefln->m_pvUserData);
		SSI_PC_ASSERT(pstReflnExtraInfo);

		dEvalMaxPixelValue = (pstReflnExtraInfo->m_dMaxPixelValue - m_nPedestal) * dTimeScale + m_nPedestal;
		//dEvalMaxPixelValueError = pow(pstReflnExtraInfo->m_dMaxPixelValue, 0.5) * dTimeScale;
		dEvalMaxPixelValueError = 0.00001; // math accuracy?

		if( dEvalMaxPixelValue - m_dImageSaturatedPixelValue > (-1.0) * dEvalMaxPixelValueError )
		{
			dSaturatedPeaksPercent += 1.0;
		}
		else
			break;  // since all peaks are intensity-sorted, there is no point looping further
	}

	if( m_oPeaks.nGetNumReflns() > 0 )
	{

        dSaturatedPeaksPercent /= (double)m_oPeaks.nGetNumReflns() / 100.0;
	}
	else
	{
		SSI_PC_ASSERT(dSaturatedPeaksPercent==0.0);
	}
	
	return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Given percentage of saturated peaks, evaluate the exposure time
bool CImageExposureTimeEvaluator::bGetExposureTimeWithSaturationLimit(double dSaturatedPeaksPercent, 
													                  double& dExposureTime_sec, 
													                  bool& bZeroPeaksSaturated)
{
	if ( 0.0 >= m_dImageExposureTime_sec )
		return false;						// cannot evaluate, if no real measurement time is available
	
	if( dSaturatedPeaksPercent < 0.0 || dSaturatedPeaksPercent > 100.0 )	// DDX should have taken care of that
		return false;	  // meaningless input percentage

	int		nTotalPeakCount = m_oPeaks.nGetNumReflns();
	if( 0 >= nTotalPeakCount )
		return false; // meaningless number of found peaks
	
    // Sort found peaks on maximum pixel value
    //vSortPeakList(eReflnField_float_type, m_oPeaks.m_nFI_fIntensity);
    vSortPeakList(eReflnField_float_type, m_nMaxPixelValueFieldIndex);

	// How many peaks should be saturated to satisfy the input value?
	double		dRequiredSaturatedPeaksCount = nTotalPeakCount * dSaturatedPeaksPercent / 100.0;
	int 		nRequiredSaturatedPeaksCount = ceil(dRequiredSaturatedPeaksCount) - dRequiredSaturatedPeaksCount <= 0.5 ?
												ceil(dRequiredSaturatedPeaksCount) : floor(dRequiredSaturatedPeaksCount);	

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if( nRequiredSaturatedPeaksCount == 0 )
	{
		//		nRequiredSaturatedPeaksCount = nTotalPeakCount - 1;	  // we should have at least one peak to work with!
		nRequiredSaturatedPeaksCount = 1; // we should have at least one peak to work with!
		bZeroPeaksSaturated = true;
	}
	else
	{
		bZeroPeaksSaturated = false;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	int			nThresholdPeakIndex = nTotalPeakCount - nRequiredSaturatedPeaksCount - 1;// the time it takes to saturate this peak is what we need
	int			nThresholdPeakIndex = nTotalPeakCount - nRequiredSaturatedPeaksCount;    // the time it takes to saturate this peak is what we need
	SSI_PC_ASSERT(nThresholdPeakIndex >= 0 && nThresholdPeakIndex < nTotalPeakCount);
	
	int			nSortedIndex = m_pnSortIndex[nThresholdPeakIndex];	
	Crefln*			poRefln	= m_oPeaks.poGetRefln(nSortedIndex); 
	SSI_PC_ASSERT(poRefln);
	REFLN_EXTRA_INFO*	pstReflnExtraInfo = (REFLN_EXTRA_INFO*)(poRefln->m_pvUserData);
	SSI_PC_ASSERT(pstReflnExtraInfo);
	
	if( 0.0 == pstReflnExtraInfo->m_dMaxPixelValue - m_nPedestal )
		return false; // safety check to avoid dividing by zero
	
	
	dExposureTime_sec = m_dImageExposureTime_sec * (m_dImageSaturatedPixelValue - m_nPedestal) / 
												   (pstReflnExtraInfo->m_dMaxPixelValue - m_nPedestal);

	return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Return false if the crystal is so bad that no peaks have been found
bool CImageExposureTimeEvaluator::bSuggestImageExposureTime(const char * strImageName, 
                                                            double dMaxSaturatedPeak_Percent,
                                                            double& dRecommendedExposureTime_sec)
{
    if( '\0' != *strImageName && !ProcessImage(strImageName, true) )  // true means clean peak list
        return false;
    
    printf("\n=======\nRESULTS:\n=======\n");

    printf("Current exposure time: %.1f sec.\n", m_dImageExposureTime_sec);
    
    dRecommendedExposureTime_sec = 0.0;

    double      dSaturationLimitedExposureTime_sec = 0.0;
    bool        bZeroPeaksSaturated = false;
    if( !bGetExposureTimeWithSaturationLimit(dMaxSaturatedPeak_Percent, dSaturationLimitedExposureTime_sec, bZeroPeaksSaturated) )
        return false;

    printf("\nTarget saturated peak percentage: %.1f%%.  Max recommended exposure time: %.1f sec\n", 
                                                    dMaxSaturatedPeak_Percent,        
                                                    dSaturationLimitedExposureTime_sec);

    dRecommendedExposureTime_sec = dRoundOff(dSaturationLimitedExposureTime_sec, 1); // just in case if the rest of the function returns false...

    double      dObservedSignalToNoise = 0.0;
    int         nObservedPeakCount = 0;
    double      dSignalToNoiseOptimizedExposureTime_sec = 0.0;
    
    if ( !bGetExposureTimeWithTargetSignalToNoise(dObservedSignalToNoise, nObservedPeakCount, dSignalToNoiseOptimizedExposureTime_sec) )
        return false;
    
    printf("\nTarget resolution: %.1f. Target I/Sigma: %.1f\n", m_dTargetResolution_A,             
                                                                m_dTargetSignalToNoiseInTargetShell);
    
    printf("\nNumber of peaks found in the target resolution shell: %d\n", nObservedPeakCount);
    printf("\nObserved I/Sigma in the target resolution shell: %.1f\n", dObservedSignalToNoise);

    printf("\nMinimum recommended exposure time: %.1f sec.\n", dSignalToNoiseOptimizedExposureTime_sec);
    
    dRecommendedExposureTime_sec = min(dSaturationLimitedExposureTime_sec, dSignalToNoiseOptimizedExposureTime_sec);
    
    dRecommendedExposureTime_sec = dRoundOff(dRecommendedExposureTime_sec, 1);

    return true;
}

bool CImageExposureTimeEvaluator::bGetExposureTimeWithTargetSignalToNoise(double& dObservedSignalToNoiseInTargetShell,
                                                                          int&    nObservedPeakCountInTargetShell,
                                                                          double& dSignalToNoiseOptimizedExposureTime_sec)
{
    double      d1 = 0.0;
    double      d2 = 0.0;
    
    return bGetExposureTimeWithTargetSignalToNoise(d1,
                                                   d2,
                                                   dObservedSignalToNoiseInTargetShell,
                                                   nObservedPeakCountInTargetShell,
                                                   dSignalToNoiseOptimizedExposureTime_sec);
}

bool CImageExposureTimeEvaluator::bGetExposureTimeWithTargetSignalToNoise(double& dTargetResShell_Begin_A,
                                                                          double& dTargetResShell_End_A,
                                                                          double& dObservedSignalToNoiseInTargetShell,
                                                                          int&    nObservedPeakCountInTargetShell,
                                                                          double& dSignalToNoiseOptimizedExposureTime_sec)
{
    if( m_dTargetResolution_A <= 0.0 )
        return false; // sanity check
    
    //////////////////////////////////
    // Initialize
    m_wStatus &= ~DTREK_IETE_TOO_FEW_PEAKS_IN_TARGET_SHELL;

    dObservedSignalToNoiseInTargetShell = 0.0;
    nObservedPeakCountInTargetShell = 0;
    dSignalToNoiseOptimizedExposureTime_sec = 0.0;

    //double      dObservedBkgToNoisePower = 0.0;

    // Determine which of the resolution shells corresponds to the given target resolution
    //double  dTargetResShell_Begin_A = 0.0;
    //double  dTargetResShell_End_A = 0.0;
    dTargetResShell_Begin_A = 0.0;
    dTargetResShell_End_A   = 0.0;
    
    if ( !GetTargetResolutionShellLimits(m_dTargetResolution_A, dTargetResShell_Begin_A, dTargetResShell_End_A) )
    {
        printf("\n\nCannot determine target resolution shell.\n\n");
        return false; // target resolution cannot be achieved 
    }
    printf("\nTarget resolution shell: %.2f - %.2f A\n", dTargetResShell_Begin_A, dTargetResShell_End_A);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sort the peak list on resolution
    int     nTotalPeaks = m_oPeaks.nGetNumReflns();

    if( 0 == nTotalPeaks )
        return false; // no peaks to work with... bad crystal?

    vSortPeakList(eReflnField_float_type, m_oPeaks.m_nFI_fDetResolution);

    // Loop over all reflections to fish out the ones that are "trapped" in the target shell
    Crefln*     poRefln = NULL;
    float       fRes = 0.0;
    double      dSignalOverNoise = 0.0;
    float       fNoise = 0.0;
    float       fMaxPixelValue = 0.0;
    float       fIntegratedIntensity = 0.0;
    float       fSigmaIntegratedIntensity = 0.0;
    float       fBkg = 0.0;
    int         nSortedIndex = -1;

    Creflnlist      oTargetShellReflnList;
    int             nSignalToNoiseRatioIndex = oTargetShellReflnList.nExpandGetField(ms_sSignalToNoiseRatio, eReflnField_float_type);
    Crefln		    oTargetShellRefln(&oTargetShellReflnList);

#ifdef IETE_DEBUG_OUTPUT
    FILE*       pfOut = fopen("C:\\projects\\ete\\ExposureTimeEval.txt","awt"); // for debugging
    //double      dB2N = 0.0;  // background/noise power
#endif
    int         ii=0;

    double      dSumBkg = 0.0;
    double      dSumSignal = 0.0;
    double      dSumNoise = 0.0;
    double      dIntegratedIntensitySum = 0.0;
    double      dIOverSigmaSum = 0.0;
	
	int		iRef = 0;

    for(iRef=0; iRef < nTotalPeaks; iRef++)
    {
	    nSortedIndex = m_pnSortIndex[iRef];

	    poRefln = m_oPeaks.poGetRefln(nSortedIndex);

        fRes = poRefln->fGetField(m_oPeaks.m_nFI_fDetResolution);

        if( fRes > dTargetResShell_End_A )
            break;  // since peaks are sorted on resolution, no point looping further
        else if( fRes > dTargetResShell_Begin_A )
        {
           fNoise           = poRefln->fGetField(m_oPeaks.m_nFI_fBackgroundSigma);
           //fMaxPixelValue   = poRefln->fGetField(m_oPeaks.m_nFI_fIntensity);
           fMaxPixelValue   = poRefln->fGetField(m_nMaxPixelValueFieldIndex);
           fBkg             = poRefln->fGetField(m_oPeaks.m_nFI_fBackground);

           fIntegratedIntensity = poRefln->fGetField(m_oPeaks.m_nFI_fIntensity);
           fSigmaIntegratedIntensity = poRefln->fGetField(m_oPeaks.m_nFI_fSigmaI);

           //           if( fNoise < 1.0 || fBkg <= 1.0 || fMaxPixelValue <= fBkg || 
           //               fSigmaIntegratedIntensity <= 0.0 || fIntegratedIntensity <= 0.0 )     // safety checks
           //               continue;
           if( fNoise < 1.0 || fBkg <= 1.0 ||
               fSigmaIntegratedIntensity <= 0.0 || fIntegratedIntensity <= 0.0 )     // safety checks
               continue;
           
           //dSignalOverNoise = (double)(fMaxPixelValue - fBkg)/(double)fNoise;

           dSignalOverNoise = (double)(fIntegratedIntensity)/(double)fSigmaIntegratedIntensity;

           if( dSignalOverNoise < c_dMinimumAcceptableSignalToNoisePerPeak )
               continue; // a safety feature: in case dtfind found some spurious peaks here is a way to reject them
           dObservedSignalToNoiseInTargetShell += dSignalOverNoise;
           
           //dB2N = log((double)fBkg) / log((double)fNoise);
           //dObservedBkgToNoisePower += dB2N;

           dSumBkg += fBkg;
           dSumSignal += fMaxPixelValue;
           dSumNoise += fNoise;

           dIntegratedIntensitySum += fIntegratedIntensity;
           dIOverSigmaSum += fIntegratedIntensity/fSigmaIntegratedIntensity;
           
           oTargetShellRefln.vSetField(nSignalToNoiseRatioIndex, (float)dSignalOverNoise);
           oTargetShellReflnList.nInsert(&oTargetShellRefln);           
           
            ++ii;

            //#ifdef IETE_DEBUG_OUTPUT
            //           fprintf(pfOut, "%d\t%f\t%f\t%f\t%f\n", ii, (double)fMaxPixelValue, (double)fBkg, (double)fNoise, dSignalOverNoise);
            //#endif
        }
    }

    dSumBkg = ii > 0 ? dSumBkg/ii : 0.0;
    dSumNoise = ii > 0 ? dSumNoise/ii : 0.0;
    dSumSignal = ii > 0 ? dSumSignal/ii : 0.0;
    double  dSTNTemp = ii > 0 ? dObservedSignalToNoiseInTargetShell/ii : 0.0;
    dIntegratedIntensitySum = ii > 0 ? dIntegratedIntensitySum/ii : 0.0;
    dIOverSigmaSum = ii > 0 ? dIOverSigmaSum/ii : 0.0;

#ifdef IETE_DEBUG_OUTPUT
    fprintf(pfOut, "%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", m_dImageExposureTime_sec, ii, dSumBkg, dSumNoise, dSumSignal, dSTNTemp, dIntegratedIntensitySum, dIOverSigmaSum);

    fclose(pfOut);
#endif

    nObservedPeakCountInTargetShell = oTargetShellReflnList.nGetNumReflns();

    ///////////////////////////////////////////////////////////////////
    // Update m_nMaxObservedPeakCountInTargetShell
    // 
    if( nObservedPeakCountInTargetShell > m_nMaxObservedPeakCountInTargetShell )
         m_nMaxObservedPeakCountInTargetShell = nObservedPeakCountInTargetShell;
    ///////////////////////////////////////////////////////////////////

    
    double      dTimePower = dGetTimePowerCoefficient();

    
    // Now work with the peaks trapped in the target shell
    if( nObservedPeakCountInTargetShell < c_nMinRequiredSpotCountInTargetShell )
    {
        // There is not enough reflections in the shell. Either the crystal is not good enough, or
        // maybe it is just a matter of exposure time.  If peaks are masked by noise, we can guess <I/Sigma> is roughly 1, 
        
        m_wStatus |= DTREK_IETE_TOO_FEW_PEAKS_IN_TARGET_SHELL;
        
        dObservedSignalToNoiseInTargetShell = 1.0;
    
        dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec * pow(m_dTargetSignalToNoiseInTargetShell, dTimePower); // assuming <I/Sigma> grows proportionally to the exp time to the power 1/dTimePower
        
        return true;
    }

    dObservedSignalToNoiseInTargetShell /= nObservedPeakCountInTargetShell;
    //dObservedBkgToNoisePower /= nObservedPeakCountInTargetShell;

    //    if( dObservedBkgToNoisePower > 2.0 ) // it is unrealistic that noise grows slower than the square root of time
    //        dObservedBkgToNoisePower = 2.0;

    // Sort the target shell peaks on their signal-to-noise values
    oTargetShellReflnList.vSort(eReflnField_float_type, nSignalToNoiseRatioIndex, NULL);
    int*        pnTargetShellSortIndices = oTargetShellReflnList.pnGetSortIndex();

    //// The mean signal-to-noise might be positively biased, so we may need to get the median or the mode 
    //// of the distribution of the signal-to-noises and use that metric instead of the mean.
    
    //    // Find the median of the distribution of signal-to-noise values
    //    if( (nObservedPeakCountInTargetShell/2) * 2 != nObservedPeakCountInTargetShell ) // if the peak count is odd 
    //    {
    //        int     nMidSorted = pnTargetShellSortIndices[nObservedPeakCountInTargetShell/2];    
    //        poRefln = oTargetShellReflnList.poGetRefln(nMidSorted);
    //        dObservedSignalToNoiseInTargetShell = (double)poRefln->fGetField(nSignalToNoiseRatioIndex);
    //    }
    //    else
    //    {
    //        int     nMidSorted_1 = pnTargetShellSortIndices[nObservedPeakCountInTargetShell/2];
    //        int     nMidSorted_2 = pnTargetShellSortIndices[nObservedPeakCountInTargetShell/2 - 1];
    //        
    //        poRefln = oTargetShellReflnList.poGetRefln(nMidSorted_1);
    //        dObservedSignalToNoiseInTargetShell = (double)poRefln->fGetField(nSignalToNoiseRatioIndex);
    //        
    //        poRefln = oTargetShellReflnList.poGetRefln(nMidSorted_2);
    //        dObservedSignalToNoiseInTargetShell += (double)poRefln->fGetField(nSignalToNoiseRatioIndex);
    //
    //        dObservedSignalToNoiseInTargetShell /= 2.0;
    //    }
    //    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Get the suggested exposure time change coefficient, based solely on the Signal To Noise target
    double       dSignalToNoiseSuggestedChangeCoeff = m_dTargetSignalToNoiseInTargetShell / dObservedSignalToNoiseInTargetShell;
   
     if( dSignalToNoiseSuggestedChangeCoeff >= 1.0 )
    {
        dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec * pow(dSignalToNoiseSuggestedChangeCoeff, dTimePower); // assuming <I/Sigma> grows proportionally to the exp time to the power 1/dTimePower
        //dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec * pow(dSignalToNoiseSuggestedChangeCoeff, dObservedBkgToNoisePower);
       
       return true; // if the suggested time is higher than this image time - go for it!
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // If the suggested time is lower than this image time, watch out: we don't want to lose too many peaks, as we lower the exposure time...

    // Figure out the projected peak loss (from S/N considerations) to compare it with the ALLOWED loss.
    // If the projected loss is less than allowed - go for what the S/N ratio suggests.
    // Else go by the nCriticalSortedPeakIndexInTargetShell and suggest such an exposure time that would 
    // keep that "critical" peak above the projected noise level. All peaks listed after the critical peak 
    // will automatically survive the lower exposure time, because the peak list is sorted.
    
    int         nTargetShellSortedIndex = -1;
    int         nPojectedPeakLossCountInTargetShell = 0;
    bool        bLosingPeaks = false;
    for(iRef=0; iRef < nObservedPeakCountInTargetShell; iRef++)
    {
        nTargetShellSortedIndex = pnTargetShellSortIndices[iRef];

	    poRefln = oTargetShellReflnList.poGetRefln(nTargetShellSortedIndex);

        dSignalOverNoise = (double)poRefln->fGetField(nSignalToNoiseRatioIndex);
	    
        // Find out if the peak may dissapear as a result of lowering the exposure time
        if( dSignalOverNoise * dSignalToNoiseSuggestedChangeCoeff < c_dMinimumAcceptableSignalToNoisePerPeak ) 
        {
            nPojectedPeakLossCountInTargetShell++;
            bLosingPeaks = true;
        }
        else if( bLosingPeaks )
            break;  // since the peak list is sorted, no reason looking further
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // How many peaks are we allowed to lose anyway?
    int         nAllowedPeakLossCount = m_nMaxObservedPeakCountInTargetShell *  c_dMaximumAcceptablePeakLossCountInTargetShell_percent
                                        / 100.0 
                                        - (m_nMaxObservedPeakCountInTargetShell - nObservedPeakCountInTargetShell);

    if( nAllowedPeakLossCount > nPojectedPeakLossCountInTargetShell )
    {
       dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec * pow(dSignalToNoiseSuggestedChangeCoeff, dTimePower); // assuming <I/Sigma> grows proportionally to the exp time to the power 1/dTimePower
       //dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec * pow(dSignalToNoiseSuggestedChangeCoeff, dObservedBkgToNoisePower);
       return true; // go for the lower exposure time, if the peak loss is within the allowed limits...
    }
    else if( nAllowedPeakLossCount < 1 )
    {
       dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec;
       return true; // keep the current exposure, if cannot lose any peaks
    }
    else
    {
        // Define the critical sorted index
        int         nCriticalSortedPeakIndexInTargetShell = nAllowedPeakLossCount - 1;

        // Find the unsorted index that corresponds to the sorted
        int         nCriticalUnsortedPeakIndexInTargetShell = pnTargetShellSortIndices[nCriticalSortedPeakIndexInTargetShell];

        poRefln = oTargetShellReflnList.poGetRefln(nCriticalUnsortedPeakIndexInTargetShell);
        
        dSignalOverNoise = (double)poRefln->fGetField(nSignalToNoiseRatioIndex);
    
        double      dPeakLossLimitedSignalToNoiseSuggestedChangeCoeff = dSignalOverNoise < c_dMinimumAcceptableSignalToNoisePerPeak ? 
                                                                        1.0 : c_dMinimumAcceptableSignalToNoisePerPeak / dSignalOverNoise; 
        
        dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec
                                                  * pow(dPeakLossLimitedSignalToNoiseSuggestedChangeCoeff, dTimePower);         // assuming <I/Sigma> grows proportionally to the exp time to the power 1/dTimePower
        //dSignalToNoiseOptimizedExposureTime_sec = m_dImageExposureTime_sec
        //                                          * pow(dPeakLossLimitedSignalToNoiseSuggestedChangeCoeff, dObservedBkgToNoisePower);
        
        return true; // go for the lower exposure time with the limited peak loss
    }
   
   return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CImageExposureTimeEvaluator::vSortPeakList(eReflnFieldType eType, int nField)
{
    if( m_nSortKeyField == nField )
        return; // already sorted

    m_oPeaks.vSort(eType, nField, NULL);
    
    m_nSortKeyField = nField;
    m_pnSortIndex = m_oPeaks.pnGetSortIndex();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CImageExposureTimeEvaluator::GetTargetResolutionShellLimits(double dTargetResolution_A, 
                                                                 double& dTargetResShell_Begin_A, 
                                                                 double& dTargetResShell_End_A)
{
    dTargetResShell_Begin_A = -1.0;
    dTargetResShell_End_A   = -1.0;
    
    int     nShellLimitCount = m_daResolutionShellLimits.size();

    // Go through predetermined resolution shells and find where the input resolution falls.
    // Return false if no shell matches the input resolution

    for(int ii=0; ii < nShellLimitCount-1; ii++)
    {
        if( dTargetResolution_A >= m_daResolutionShellLimits[ii] && 
            dTargetResolution_A < m_daResolutionShellLimits[ii+1] )
        {
            dTargetResShell_Begin_A = m_daResolutionShellLimits[ii];
            dTargetResShell_End_A   = m_daResolutionShellLimits[ii+1];
            
            return true;
        }
    }
    
    return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CImageExposureTimeEvaluator::bGetSaveImageResolution(Cimage_header& oHeader)
{
    //////////////////////////////////////////////////////////////////////////
    // Get the max (edge) and min resolution
    if( !oHeader.bIsAvailable() ) // Rotation information not available
		return false;

    Csource              oSource(oHeader);
    if ( !oSource.bIsAvailable() )
        return false;

    float           a3fDirectBeamVector[3];
    oSource.vCalcGetS0(a3fDirectBeamVector);

    Cdetector       oDetector(oHeader);
    if ( !oDetector.bIsAvailable() )
        return false;

    float   fJunk = 0.0;
    if( 0 != oDetector.nGetResolution(a3fDirectBeamVector, &m_fImageResMin, &fJunk, &m_fImageResMax) )
        return false;

    float           fWavelength_A = oSource.fGetWavelength();
    m_fImageResMin *= fWavelength_A;
    m_fImageResMax *= fWavelength_A;

    ///////////////////////////////////////////////////////////////////////////////////////
    // safety check
    if( m_fImageResMin <= 0.0 || m_fImageResMax <= 0.0 || m_fImageResMax > m_fImageResMin )
        return false;
    
    if( !bSetResolutionShellLimits() )
        return false;
    
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CImageExposureTimeEvaluator::bSetResolutionShellLimits()
{
    ////////////////////////////////////////////////////////////////////////
    // Safety check
    if( m_dTargetResolution_A <= 0.0 )
        return false;
    
    const   int     cnresolutionShellCount = 7;  // For LYSO 10 shells would make sense, if we wanted the corner resolution,
                                                 // but we need the edge one, so 7 is sufficient.   
    m_daResolutionShellLimits.clear();
    
    ///////////////////////////////////////////////////////////////////////
    // Define the maximum resolution in such a way that if the volume between the maximum resolution 
    // and the minimum resolution is broken down into 10 parts, then the target resolution will be
    // in the middle of the last (10th) shell
    float       fResMax = pow(m_dTargetResolution_A * (double)m_fImageResMin, 3.);
    fResMax /= cnresolutionShellCount * 2.0 * pow((double)m_fImageResMin, 3.0) - pow(m_dTargetResolution_A, 3.0);
    fResMax = pow((double)(cnresolutionShellCount * 2 - 1) * fResMax, 1./3.);
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Calculate a resolution shell volume (reciprocal space) without the coefficient of 4/3 PI, 
    // because that coefficient will be present on both sides of the equation,
    // when it comes to calculating the shell limits
    double      fResoMaxCubed_reciprocal = 1.0 / pow((double)fResMax, 3.);
    double      fResoMinCubed_reciprocal = 1.0 / pow((double)m_fImageResMin, 3.);
    double      fResoShellVolume_reciprocal = (fResoMaxCubed_reciprocal - fResoMinCubed_reciprocal) / cnresolutionShellCount;
    
    // Now fill the vector of resolution shell limits (direct space)
    double      fResoShellLimit = 0.0;
    m_daResolutionShellLimits.push_back((double)fResMax); 
    for(int nn=1; nn < cnresolutionShellCount; nn++)
    {
        fResoShellLimit = 1.0 / pow(fResoMaxCubed_reciprocal - nn * fResoShellVolume_reciprocal, 1.0/3.0); 
        m_daResolutionShellLimits.push_back(fResoShellLimit);
    }
    m_daResolutionShellLimits.push_back((double)m_fImageResMin);
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CImageExposureTimeEvaluator::dGetTimePowerCoefficient()
{
    return m_bCCDImage ? c_dTimePower_CCD : c_dTimePower_IP;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CImageExposureTimeCalculator::bCalculate(Cimage_header& oHeader,
                                              double   fResoMin,
                                              double   fResoMax,
                                              int      nResoBins,
                                              double   fImageRotWidth,
                                              double   fTargetIoverSigma,
                                              double&  fProposedExposureTime,
                                              DTREK_WORD wCtrl)
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Validate input
    if( !bValidateResolutionInput(fResoMin, fResoMax) )
        return false;

    if( nResoBins <= 0 )
    {
        m_sError = "Input number of resolution bins is zero or negative.";
        return false;
    }

    if( fImageRotWidth <= 0.0 )
    {
        m_sError = "Input image rotation width is zero or negative.";
        return false;
    }

    if( fTargetIoverSigma <= 0.0 )
    {
        m_sError = "Input target I/Sigma is zero or negative.";
        return false;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    fProposedExposureTime = 0.0;
    m_sError = "";

    double  fSaveTargetIOverSigmaInTheLastResoShell = fTargetIoverSigma;
    if( wCtrl == DTREK_IETC_FULL_MERGED_REFLNS )
    {
        fTargetIoverSigma *= 4.0; // Fudge factor. We know that I/Sigma after screening is ca 4 times larger
                                  // than the I/Sigma after dtscaleaverage.  We need to work more on this.
    }

    Cstring              sExposureTimeEvaluationInfo(""); 
    Cstring              sExpTimeHeaderKeyword("DTREFINE_RANK_");
    sExpTimeHeaderKeyword += Cstring(D_K_ExposureTimeEvaluationInfo);
    if ( 0 != oHeader.nGetValue(sExpTimeHeaderKeyword, &sExposureTimeEvaluationInfo) )
    {
        m_sError = "The header file does not have necessary information from dtranker (1)";
        return false; 
    }                    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<double>     afOptions;
    sExposureTimeEvaluationInfo.nListToVector(afOptions, " ");
    if( afOptions.size() < 6 )
    {
        m_sError = "The header file does not have necessary information from dtranker (2)";
        return false;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // We will skip sigmaA and sigmaB for now.
    double      fParam_A = afOptions[0];
    double      fParam_B = afOptions[2];
    double      fParam_ImageRotationWidth  = afOptions[4];
    double      fParam_ImageExposureTime   = afOptions[5];
    double      fParam_IOverSigmaTimePower = afOptions[6];
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // We need to figure out the target resolution from the user resolution range and the number of resolution bins
    // Generate bins for overall resolution range
    CResoStrategyBins*  poResoBins = new CResoStrategyBins(fResoMin, fResoMax, nResoBins);  
    CResoBin*   poBin = poResoBins->pGetBin(nResoBins-1);
    if( !poBin )
    {
        m_sError = "Failed to get the last resolution bin object";
        return false;
    }

    double      fTargetReso = poBin->fGetMiddleReso();
    delete poResoBins;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double      fExtrapolatedIOverSigmaAtTargetResolution = pow(10.0, (fParam_A +
                                                            fParam_B * (0.25/(float)fTargetReso/(float)fTargetReso) ) );

    fProposedExposureTime = pow(fParam_ImageExposureTime, fParam_IOverSigmaTimePower) *
                                fTargetIoverSigma/fExtrapolatedIOverSigmaAtTargetResolution;
    fProposedExposureTime = pow(fProposedExposureTime, 1.0/fParam_IOverSigmaTimePower);
    fProposedExposureTime *= fImageRotWidth/fParam_ImageRotationWidth; // assuming that exposure time 
                                                                       // should grow proportionally 
                                                                       // to the image rotation width
    
    std::vector<int>        anExposureTimeCheckValues;
    anExposureTimeCheckValues.push_back(1);
    anExposureTimeCheckValues.push_back(2);
    anExposureTimeCheckValues.push_back(5);
    anExposureTimeCheckValues.push_back(10);
    anExposureTimeCheckValues.push_back(15);
    anExposureTimeCheckValues.push_back(20);
    anExposureTimeCheckValues.push_back(30);
    anExposureTimeCheckValues.push_back(40);
    anExposureTimeCheckValues.push_back(60);
    anExposureTimeCheckValues.push_back(90);

    fProposedExposureTime = nGetNearestIntFromVector(fProposedExposureTime, 
                                                     anExposureTimeCheckValues, 
                                                     DTREK_VEC_ROUNDOFF_CEILING);

    printf("\nExposure time evaluation.\n");
    printf("Number of resolution shells: %d\n", nResoBins);
    printf("Target resolution (mean resolution of the highest resolution shell): %.3f A\n", fTargetReso);
    printf("Target I/Sigma in the highest resolution shell: %.1f\n", fSaveTargetIOverSigmaInTheLastResoShell);
    printf("Image rotation width: %.2f deg\n", fImageRotWidth);
    printf("Proposed exposure time: %.1f sec\n\n", fProposedExposureTime);

    // Save the exposure time evaluation results in the header
    char    cTemp[256];
    sprintf(cTemp, "%.0f %.2f %.2f %.1f", fProposedExposureTime, 
                                          fImageRotWidth,
                                          fTargetReso,
                                          fSaveTargetIOverSigmaInTheLastResoShell);

    oHeader.nReplaceValueDictionary(Cstring("DTMULTISTRATEGY"), Cstring(D_K_ProposedExposureTime), cTemp);

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CImageExposureTimeCalculator::bValidateResolutionInput(double& fResoMin, double& fResoMax)
{
    const double    cfMaxLowReso = 999.99;
    
    m_sError = "";
    char    cTemp[256];

    double       fReso1 = fResoMin;
    double       fReso2 = fResoMax;

    // If a value is not supplied, treat it as if the low resolution is not supplied.
     if( fReso1 <= 0.0 )
        fReso1 = cfMaxLowReso; // a low resolution default

    if( fReso2 <= 0.0 )
        fReso2 = cfMaxLowReso; // a low resolution default

    // Just a safety check
    if( cfMaxLowReso == fReso1 && cfMaxLowReso == fReso2 )
    {
        sprintf(cTemp, "Wrong resolution input: %f %f", fResoMin, fResoMax);
        m_sError = cTemp;
        return false;
    }

    //OK, now where is the low and where is the high?
    fResoMin = max(fReso1, fReso2);
    fResoMax = min(fReso1, fReso2);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////
Cstring CImageExposureTimeCalculator::sGetLastError()
{
    Cstring     sError("ERROR: Failed to calculate the exposure time. \n");
    sError += m_sError;

    return sError;
}
///////////////////////////////////////////////////////////////////////////////////////////////





















