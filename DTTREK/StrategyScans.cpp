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
// StrategyScans.cpp       Initial author: RB     06-May-2005
// This file contains the member functions of class CStrategyScans

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
#include <algorithm>
#include "dtreksys.h"
#include "StrategyScans.h"
#include "Cimage_header.h"
#include "Crotation.h"
#include "Cgoniometer.h"

static bool bScansCompare(const CstrategyScan& first, const CstrategyScan& second);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyAxis::CstrategyAxis(const CstrategyAxis& anotherAxis)
{
    m_nInfoIndex        = anotherAxis.m_nInfoIndex;   // a reference to the information about this axis in the AxisInfo
    m_dValue            = anotherAxis.m_dValue;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyAxis& CstrategyAxis::operator=(const CstrategyAxis& anotherAxis)
{
    m_nInfoIndex        = anotherAxis.m_nInfoIndex;   // a reference to the information about this axis in the AxisInfo
    m_dValue            = anotherAxis.m_dValue;

    return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyAxisInfo::CstrategyAxisInfo(const CstrategyAxisInfo& anotherAxisInfo)
{
    m_strName           = anotherAxisInfo.m_strName;
    m_bRotAxis          = anotherAxisInfo.m_bRotAxis;
    m_dCollisionOffset  = anotherAxisInfo.m_dCollisionOffset;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyAxisInfo& CstrategyAxisInfo::operator=(const CstrategyAxisInfo& anotherAxisInfo)
{
    m_strName           = anotherAxisInfo.m_strName;
    m_bRotAxis          = anotherAxisInfo.m_bRotAxis;
    m_dCollisionOffset  = anotherAxisInfo.m_dCollisionOffset;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyScan::CstrategyScan()
{
    vInit();
}
//////////////////////////////////////////////////////////////////////////////////
CstrategyScan::CstrategyScan(const CstrategyScan& anotherScan)
{
    m_aoAxes.clear();
    for(int ii=0; ii < (int)anotherScan.m_aoAxes.size(); ii++)
    {
        m_aoAxes.push_back(anotherScan.m_aoAxes[ii]);
    }

    m_dRot_Begin     [enSingleScan]   = anotherScan.m_dRot_Begin[enSingleScan];  
    m_dRot_End       [enSingleScan]   = anotherScan.m_dRot_End[enSingleScan];
    
    m_dCompleteness  [enSingleScan]   = anotherScan.m_dCompleteness[enSingleScan];
    m_dRedundancy    [enSingleScan]   = anotherScan.m_dRedundancy[enSingleScan];
    m_dRedundancyDev [enSingleScan]   = anotherScan.m_dRedundancyDev[enSingleScan];
    
    m_bUsed          [enSingleScan]   = anotherScan.m_bUsed[enSingleScan];      
    
    ///////////////////////////////////////////////////////////////////////////////
    m_dRot_Begin     [enMultipleScan]    = anotherScan.m_dRot_Begin[enMultipleScan];  
    m_dRot_End       [enMultipleScan]   = anotherScan.m_dRot_End[enMultipleScan];
    

    m_dCompleteness  [enMultipleScan]   = anotherScan.m_dCompleteness[enMultipleScan];
    m_dRedundancy    [enMultipleScan]   = anotherScan.m_dRedundancy[enMultipleScan];
    m_dRedundancyDev [enMultipleScan]   = anotherScan.m_dRedundancyDev[enMultipleScan];
    
    m_bUsed          [enMultipleScan]   = anotherScan.m_bUsed[enMultipleScan];      
    
    m_segRotLimits                      = anotherScan.m_segRotLimits;
    m_dMaxRotRange                      = anotherScan.m_dMaxRotRange;

    m_nOriginalScanIndex                = anotherScan.m_nOriginalScanIndex;

    m_dDetRelZero[enStratScanDetRelZero_TwoTheta] = anotherScan.m_dDetRelZero[enStratScanDetRelZero_TwoTheta];
    m_dDetRelZero[enStratScanDetRelZero_Distance] = anotherScan.m_dDetRelZero[enStratScanDetRelZero_Distance];

    m_fImageRotWidth = anotherScan.m_fImageRotWidth;
    m_fImageExpTime  = anotherScan.m_fImageExpTime;


	//shijie yao : 2008.03.12
	//m_dIncrementalWidth = anotherScan.dGetIncrementalWidth();
	m_dCryGonioValMax = anotherScan.dGetCryGonioValMax();
	m_dCryGonioValMin = anotherScan.dGetCryGonioValMin();

	m_enPaddedType[0]	= anotherScan.m_enPaddedType[0];
	m_enPaddedType[1]	= anotherScan.m_enPaddedType[1];
	m_dPaddedVal[0]		= anotherScan.m_dPaddedVal[0];
	m_dPaddedVal[1]		= anotherScan.m_dPaddedVal[1];
	m_enShiftedType[0]	= anotherScan.m_enShiftedType[0];
	m_enShiftedType[1]	= anotherScan.m_enShiftedType[1];
	m_dShiftedVal[0]	= anotherScan.m_dShiftedVal[0];
	m_dShiftedVal[1]	= anotherScan.m_dShiftedVal[1];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyScan::CstrategyScan(Cimage_header* poHeader, int nScanIndex, int nNumAxes, const Cstring& sStrategyPrefix)
{
    vInit();

    // Get gonio values
    
    // Build scan prefix
    Cstring         strScanCrystalPrefix = sBuildStrategyScanPrefix(sStrategyPrefix, D_K_CrystalPrefix, nScanIndex);

    double*         pdGonioValues = new double[nNumAxes];
    int             nStat = poHeader->nGetValue(strScanCrystalPrefix + D_K_GonioValues, nNumAxes, pdGonioValues);
    
    if( 0 != nStat )
    {
        delete [] pdGonioValues;
        return;
    }

    // Add axes objects 
    for(int ii=0; ii < nNumAxes; ii++)
    {
        m_aoAxes.push_back(CstrategyAxis(ii, pdGonioValues[ii]));
    }

    delete [] pdGonioValues;
    ///////////////////////////

    Cstring     sScanPrefix = sBuildStrategyScanPrefix(sStrategyPrefix, D_K_ScanPrefix, nScanIndex);
    Crotation     oRotation(*poHeader, sScanPrefix);
    
    m_dRot_Begin[enSingleScan] = m_dRot_Begin[enMultipleScan] = (double)oRotation.fGetRotStart();
    m_dRot_End  [enSingleScan] = m_dRot_End  [enMultipleScan] = (double)oRotation.fGetRotEnd();
    
    m_dMaxRotRange = (double)oRotation.fGetRotRange();

    m_segRotLimits.vSet((double)oRotation.fGetRotMin(), (double)oRotation.fGetRotMax());
    
    m_nOriginalScanIndex = nScanIndex;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // See if the header uses D_K_ScanDetDatum keyword to specify the detector position. If it does not - no problem,
    // the detector position will be read by the container class CstrategyScans using detector D_K_GonioValues keyword.  
    Cstring     sScanDetRelZero(D_K_ScanDetDatum);

    double      dDetRelZero[3] = {0.0};
    if( 0 == poHeader->nGetValue(sScanPrefix + sScanDetRelZero, 3, dDetRelZero) ) // the first value should be equal to 2, then two theta and distance
    {
        m_dDetRelZero[enStratScanDetRelZero_TwoTheta] = dDetRelZero[1];
        m_dDetRelZero[enStratScanDetRelZero_Distance] = dDetRelZero[2];
    }


	///////////////////////////
	//shijie yao : 2008.03.12
	//Noticed there was a var m_fImageRotWidth defined, which should server the 
	//same as m_dIncrementalWidth. But the original code didn't set its value
	//properly. It's value was set from the STRATEGY_INPUT_S*_SCAN_ROTATION,
	//which only keeps the original copy of the SCAN_ROTATION values, even the
	//-rotimage option is provided.
	//After adding code in CDTMainMultiStrategy::bGetImageRotationWidth() to
	//update the STRATEGY_INPUT_S*_SCAN_ entries together with the SCAN_ROTATION,
	//for the new imagewidth value, m_fImageRotWidth and m_dIncrementalWidth are 
	//sequal now. 
	//
	m_fImageRotWidth = oRotation.fGetIncrement();
	m_dCryGonioValMax = oRotation.fGetRotMax();
	m_dCryGonioValMin = oRotation.fGetRotMin();
	//Crotation  oScanRot(*poHeader, D_K_ScanPrefix);	// use SCAN_ROTATION values for rotImageWidth
	//m_dIncrementalWidth = oScanRot.fGetIncrement();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyScan& CstrategyScan::operator=(const CstrategyScan& anotherScan)
{
    m_aoAxes.clear();
    for(int ii = 0; ii < (int)anotherScan.m_aoAxes.size(); ii++)
    {
        m_aoAxes.push_back(anotherScan.m_aoAxes[ii]);
    }

    m_dRot_Begin     [enSingleScan]    = anotherScan.m_dRot_Begin[enSingleScan];  
    m_dRot_End       [enSingleScan]   = anotherScan.m_dRot_End[enSingleScan];
    
    m_dCompleteness  [enSingleScan]   = anotherScan.m_dCompleteness[enSingleScan];
    m_dRedundancy    [enSingleScan]   = anotherScan.m_dRedundancy[enSingleScan];
    m_dRedundancyDev [enSingleScan]   = anotherScan.m_dRedundancyDev[enSingleScan];
    
    m_bUsed          [enSingleScan]   = anotherScan.m_bUsed[enSingleScan];      
    
    ///////////////////////////////////////////////////////////////////////////////
    m_dRot_Begin     [enMultipleScan]    = anotherScan.m_dRot_Begin[enMultipleScan];  
    m_dRot_End       [enMultipleScan]   = anotherScan.m_dRot_End[enMultipleScan];
    
    m_dCompleteness  [enMultipleScan]   = anotherScan.m_dCompleteness[enMultipleScan];
    m_dRedundancy    [enMultipleScan]   = anotherScan.m_dRedundancy[enMultipleScan];
    m_dRedundancyDev [enMultipleScan]   = anotherScan.m_dRedundancyDev[enMultipleScan];
    
    m_bUsed          [enMultipleScan]   = anotherScan.m_bUsed[enMultipleScan];      
    
    m_segRotLimits                      = anotherScan.m_segRotLimits;
    m_dMaxRotRange                      = anotherScan.m_dMaxRotRange;

    m_nOriginalScanIndex                = anotherScan.m_nOriginalScanIndex;

    m_dDetRelZero[enStratScanDetRelZero_TwoTheta] = anotherScan.m_dDetRelZero[enStratScanDetRelZero_TwoTheta];
    m_dDetRelZero[enStratScanDetRelZero_Distance] = anotherScan.m_dDetRelZero[enStratScanDetRelZero_Distance];

    m_fImageRotWidth = anotherScan.m_fImageRotWidth;
    m_fImageExpTime  = anotherScan.m_fImageExpTime;

	//shijie yao : 2008.03.12
	//m_dIncrementalWidth = anotherScan.dGetIncrementalWidth();
	m_dCryGonioValMax	= anotherScan.dGetCryGonioValMax();
	m_dCryGonioValMin	= anotherScan.dGetCryGonioValMin();
	m_enPaddedType[0]	= anotherScan.m_enPaddedType[0];
	m_enPaddedType[1]	= anotherScan.m_enPaddedType[1];
	m_dPaddedVal[0]		= anotherScan.m_dPaddedVal[0];
	m_dPaddedVal[1]		= anotherScan.m_dPaddedVal[1];
	m_enShiftedType[0]	= anotherScan.m_enShiftedType[0];
	m_enShiftedType[1]	= anotherScan.m_enShiftedType[1];
	m_dShiftedVal[0]	= anotherScan.m_dShiftedVal[0];
	m_dShiftedVal[1]	= anotherScan.m_dShiftedVal[1];
    return *this;
}

//shijie yao : 2008.03.14
// shift the range to make it start with a whole number.
//
//
//
//shijie yao : 2008.03.12
// make the range an integral # of m_incrementalWidth while staying within hardware limits.
// the hardare limits are taking from the first value of 
// CRYSTAL_GONIO_VALUES_MAX
// CRYSTAL_GONIO_VALUES_MIN
// from the header file. 
void CstrategyScan::vPaddingRange(enStrategyScanType type)
{
	//Init
	m_enPaddedType[type]	= enNone;
	m_enShiftedType[type]	= enNone;
	m_dPaddedVal[type]	= 0;
	m_dShiftedVal[type]	= 0;

	double dRotBegin = m_dRot_Begin[type];
	double dRotEnd   = m_dRot_End[type];
	double dWidth    = fabs(dRotBegin - dRotEnd);

	int nScanCnt = (int)(dWidth / m_fImageRotWidth);
	double dDiffAng = nScanCnt * m_fImageRotWidth - dWidth; //is always a + number

	bool bNeedPad = (dDiffAng != 0);
		
	//[1] Pad the ends so that the range will contain an integral # of images of width m_fImageRotWidth,
	//    while staying within the hardware limits:
	if(bNeedPad == true)
	{
		//[1.0] Angle value needs to be pad 
		double dHalfPadWidth = ((nScanCnt+1) * m_fImageRotWidth - dWidth) * 0.5;

		//[1.1] Add half of the dDiffAng to front and end, if possible
		if(	dRotBegin - dHalfPadWidth >= m_dCryGonioValMin &&
			dRotEnd + dHalfPadWidth <= m_dCryGonioValMax)
		{
			m_dRot_Begin[type]  -= dHalfPadWidth;
			m_dRot_End[type]    += dHalfPadWidth;
			m_enPaddedType[type] = enBothEndExtended;
			m_dPaddedVal[type]   = dHalfPadWidth;
		}
		//[1.2] Add all the dDiffAng to front, if possible
		else if(dRotBegin - dDiffAng >= m_dCryGonioValMin)
		{
			m_dRot_Begin[type]   -= dDiffAng;
			m_enPaddedType[type]  = enFrontExtended;
			m_dPaddedVal[type]    = dDiffAng;
		}
		//[1.3] Add all the dDiffAng to the end 
		else if(dRotEnd + dDiffAng <= m_dCryGonioValMax)
		{
			m_dRot_End[type]     += dDiffAng;
			m_enPaddedType[type]  = enEndExtended;
			m_dPaddedVal[type]    = dDiffAng;
		}
		//[1.4] Down size the range on both ends to make it fit
		else
		{
			double dHalfExtraAng = ( dWidth - nScanCnt * m_fImageRotWidth ) * 0.5;
			m_dRot_Begin[type]  += dHalfExtraAng;
			m_dRot_End[type]    -= dHalfExtraAng;	
			m_enPaddedType[type] = enBothEndShrinked;
			m_dPaddedVal[type]   = dHalfExtraAng;
		}
	}

	//[2] Shift the range so that it will start with a whole number,
	//while staying in the hardware limits:
	int intPartOfStart = (int)(m_dRot_Begin[type]);
	double shiftVal = m_dRot_Begin[type] - intPartOfStart;	//positive fractional value of start;
	double shiftBackVal = 1.0 - shiftVal;

	if(shiftVal != 0.0)
	{
		bNeedPad = true;
		//[2.1] Shift to front
		if(intPartOfStart >= m_dCryGonioValMin)
		{
			m_dRot_Begin[type]	-= shiftVal;
			m_dRot_End[type]	-= shiftVal;

			m_enShiftedType[type]	= enShiftedToFront;
			m_dShiftedVal[type]		= shiftVal;
		}
		//[2.2] Shift to back
		else
		{
			m_dRot_Begin[type]	+= shiftBackVal;
			m_dRot_End[type]	+= shiftBackVal;

			m_enShiftedType[type]	= enShiftedToEnd;
			m_dShiftedVal[type]		= shiftBackVal;
		}

		//[2.3] Shift to back and delete one image width from back
		if(m_dRot_End[type] > m_dCryGonioValMax)
		{
			m_enShiftedType[type]	= enShiftedToEndAndChop;
			m_dShiftedVal[type]		= 0;
		}

		while(m_dRot_End[type] > m_dCryGonioValMax)
		{
			m_dRot_End[type]	-= m_fImageRotWidth;
			m_dShiftedVal[type]	+= m_fImageRotWidth;
		}
	}
};

void CstrategyScan::vSetRotation(const double dBegin, const double dEnd, enStrategyScanType eType)
{
	m_dRot_Begin[eType]	= dBegin; 
	m_dRot_End[eType]	= dEnd;
	vPaddingRange(eType);
}
////////////////////////////////////////////////////////////////////////////////////////
void CstrategyScan::vInit()
{
    m_aoAxes.clear();

    m_dRot_Begin     [enSingleScan] = 0.0;  
    m_dRot_End       [enSingleScan] = 0.0;
    
    m_dCompleteness  [enSingleScan] = 0.0;
    
    m_dRedundancy    [enSingleScan] = 0.0;
    m_dRedundancyDev [enSingleScan] = 0.0;
    
    m_bUsed          [enSingleScan] = false;      
    /////////////////////////////////////////

    m_dRot_Begin     [enMultipleScan] = 0.0;  
    m_dRot_End       [enMultipleScan] = 0.0;
    
    m_dCompleteness  [enMultipleScan] = 0.0;
    
    m_dRedundancy    [enMultipleScan] = 0.0;
    m_dRedundancyDev [enMultipleScan] = 0.0;
    
    m_bUsed          [enMultipleScan] = false;
    
    m_segRotLimits.vSet(0.0, 0.0);
    m_dMaxRotRange = 0.0;
    m_nOriginalScanIndex = -1;

    m_dDetRelZero[enStratScanDetRelZero_TwoTheta] = cdStratDetUnrealisticValue;
    m_dDetRelZero[enStratScanDetRelZero_Distance] = cdStratDetUnrealisticValue;

    m_fImageRotWidth = -1.0;
    m_fImageExpTime = -1.0;

	//shijie yao : 2008.03.12
	m_enPaddedType[0]	= m_enPaddedType[1]		= enNone;
	m_enShiftedType[0]	= m_enShiftedType[1]	= enNone;

	m_dPaddedVal[0]		= m_dPaddedVal[1]		= 0.0;
	m_dShiftedVal[0]	= m_dShiftedVal[1]		= 0.0;
}
//////////////////////////////////////////////////////////////////////////////////
double CstrategyScan::dGetAxisValue(int iIndex)
{
    if( iIndex < 0 || iIndex > (int)m_aoAxes.size() - 1 )
        return 0.0;

    return m_aoAxes[iIndex].dGetValue();
}
//////////////////////////////////////////////////////////////////////////////////
bool CstrategyScan::bSetAxisValue(int iIndex, double dValue)
{
    if( iIndex < 0 || iIndex > (int)m_aoAxes.size() - 1 )
        return false;

    m_aoAxes[iIndex].vSetValue(dValue);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////
CstrategyScans::CstrategyScans() :
m_seqResolution(-1.0,-1.0)
{
    m_dTwoTheta = cdStratDetUnrealisticValue;
    m_dDistance = cdStratDetUnrealisticValue;

    m_dTwoThetaCollisionOffset = cdStratDetUnrealisticValue;
    m_dDistanceCollisionOffset = cdStratDetUnrealisticValue;

    m_bMultipleScanSolution = false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyScans::~CstrategyScans()
{
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CstrategyScans::vDeleteScan(int iIndex)
{
    if( iIndex < 0 || iIndex > m_aoScans.size() - 1 )
        return;
    
    std::vector<CstrategyScan>::iterator    it = m_aoScans.begin() + iIndex; 
    m_aoScans.erase(it);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function updates scans in the header object.  Return value - number of scans updated.
int CstrategyScans::nUpdateHeader(Cimage_header* poHeader, 
                                  int nScanNumberOffset, 
                                  const Cstring& sStrategyPrefix,
                                  bool bUpdateDetectorPosition)
{
    if( !poHeader )
        return 0; //safety check

    int     nScansUpdated = 0;
    int     nScanNumberForRecord = -1;
    
    CstrategyScan::enStrategyScanType      eScanType = !m_bMultipleScanSolution ? CstrategyScan::enSingleScan : CstrategyScan::enMultipleScan;
    
    // Let's create a rotation object based on scan 0, then substitute variable fields in a loop. 
    Cstring     sScanPrefix(D_K_ScanPrefix);
    Cstring     sScanDetRelZero(D_K_ScanDetDatum);
    
    Crotation       oRotation(*poHeader, sScanPrefix);
    if( !oRotation.bIsAvailable() )  // The scan must be in the header.
    {
        printf("\nERROR: Cannot create rotation object for writing a header\n");
        return 0;
    }

    CSegment    segRotLimits(0.0, 0.0);
    
    Cstring     sScanCrystalPrefix("");
    int             nNumAxes = nGetNumAxes();
    double*         pdGonioValues = new double[nNumAxes];

    int             nUsed = -1;
    double          dDetRelZero[3] = {2.0};  // 2, because there will be two values: 2Theta and Distance
    
    double          fExposureTime = -1.0;
    for(int ii=0; ii < (int)m_aoScans.size(); ii++)
    {
        if( !bIsScanUsed(ii, eScanType) )
            continue;
        
        nUsed++;
        nScanNumberForRecord = nUsed + nScanNumberOffset; 

        
        oRotation.vSetRotStart((float)dGetScanRotBegin(ii, eScanType));
        oRotation.vSetRotEnd((float)dGetScanRotEnd(ii, eScanType));   

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Check if the the rotation limits have been set. If not, they will be both 0.
        m_aoScans[ii].vGetRotLimits(segRotLimits);

        if( !segRotLimits.bIsEmpty() )
            oRotation.vSetRotMinMaxRange((float)segRotLimits.dGetBeg(), (float)segRotLimits.dGetEnd(), (float)segRotLimits.dGetWidth());
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        fExposureTime = m_aoScans[ii].fGetImageExposureTime();
        if( fExposureTime > 0.0 )
            oRotation.vSetExposureTime(fExposureTime);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        sScanPrefix = sBuildStrategyScanPrefix(sStrategyPrefix, D_K_ScanPrefix, nScanNumberForRecord);

        if( 0 != oRotation.nUpdateHeader(poHeader, sScanPrefix) )
            continue;
        
        // Now update detector position
        if( bUpdateDetectorPosition )
        {
            dDetRelZero[1] = dGetTwoTheta();
            dDetRelZero[2] = dGetDistance();

            poHeader->nReplaceValue(sScanPrefix + sScanDetRelZero, 3, dDetRelZero); 
        }
        /////////////////////////////////////////////////////////////////////////

        // Now update setting angles
        sScanCrystalPrefix = sBuildStrategyScanPrefix(sStrategyPrefix, D_K_CrystalPrefix, nScanNumberForRecord);
       
        for(int jj=0; jj < nNumAxes; jj++)
        {
            pdGonioValues[jj] = dGetScanAxisValue(ii, jj);
        }
        
        poHeader->nReplaceValue(sScanCrystalPrefix + D_K_GonioValues, nNumAxes, pdGonioValues);
    
        nScansUpdated++;
    }
    
    delete [] pdGonioValues;
    
    return nScansUpdated;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CstrategyScans::nPrintResultScanTable(int nScanNumberOffset, 
                                          bool bPrintHeader, 
                                          double& dTotalRotWidth,
                                          DTREK_WORD wCtrl,
                                          REFLN_FILE_STRAT_INFO* pReflnFileStrategyInfo)
{
    dTotalRotWidth = 0.0;
    
    if( 0 == nGetNumScans() )
        return 0; // nothing to do

    if( 1 > nGetNumAxes() )
        return 0; // Should not happen. There should be at least one rotation axis. 
    //////////////////////////////////////////////////////////////////////////////
    
    if( bPrintHeader )
    {
        Cstring     sHorizLine("----------------------------------------------------------------------------\n");
    
        // We need to figure out the axes names for the table header.
        // Assume every scan has the same axes names.
        Cstring     sHeader (" # 2Theta  Dist");
        Cstring     sTemp("");
        Cstring     sTemp2("");

        int                 nAxisNameLength = 0;
        int                 nNumberOfPaddingSpaces = 0;

        const int           c_nMaxAxisColumnWidth = 7;// assuming that an axis name would not be longer than 6 characters + a space for separation
    
        for(int jj=0; jj < nGetNumAxes(); jj++)
        {
            vGetAxisName(jj, sTemp);
        
            nAxisNameLength = sTemp.length();
            nNumberOfPaddingSpaces = c_nMaxAxisColumnWidth - nAxisNameLength;
        
            sTemp2 = "";
            for(int kk=0; kk < nNumberOfPaddingSpaces; kk++)
            {
                sTemp2 += " ";
            }
            sTemp2 += sTemp;
        
            sHeader += sTemp2; 
        }
    
        sHeader += "   Start    End  Width  %%Comp  Redn  +/-\n";

        if( DTREK_STRATSCANS_PRST_ALL != wCtrl )
            printf("\nBest Multiple Scan Solution\n");
        
        printf(sHorizLine.string());
        printf(sHeader.string());
        printf(sHorizLine.string());
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Print fixed scan information
    if( pReflnFileStrategyInfo )
    {
        pReflnFileStrategyInfo->m_sFilePath.vAbbreviateFilePath(58);
        
        int     nNumberOfCharsBeforeStats = 2 + // scan number for record
                                           14 +    // 2Theta and distance
                                           nGetNumAxes() * 7 +  // crystal gonio axes values
                                           22;  // rotation axis values
                                            
        Cstring     sNum(nNumberOfCharsBeforeStats);
        Cstring     sFilePathFormat("%-");
        sFilePathFormat += sNum;
        sFilePathFormat += 's';

        printf(sFilePathFormat.string(), pReflnFileStrategyInfo->m_sFilePath.string());
        printf("%6.1f %5.2f %4.2f\n", pReflnFileStrategyInfo->m_fCompleteness, 
                                      pReflnFileStrategyInfo->m_fRedundancy, 
                                      pReflnFileStrategyInfo->m_fRedundancyDev);
    }

    CstrategyScan::enStrategyScanType      eScanType = !m_bMultipleScanSolution ? CstrategyScan::enSingleScan : CstrategyScan::enMultipleScan;

    bool        bUsed = false;
    
    double      dRotBegin = 0.0;
    double      dRotEnd = 0.0;
    double      dWidth = 0.0;

    Cstring     strYesNo("");
    int         nScanNumberForRecord = -1;
    int         nUsedScanNumber = -1;
    for(int ii=0; ii < nGetNumScans(); ii++)                                 
    {
        bUsed = bIsScanUsed(ii, eScanType);

        // Multiple scan solution should only list used scans
        if( !bUsed && DTREK_STRATSCANS_PRST_ALL != wCtrl )
            continue;
        
        nUsedScanNumber++;
        nScanNumberForRecord = nUsedScanNumber + nScanNumberOffset;

        printf("%2d", nScanNumberForRecord); // scan number
        
        printf(" %6.1f %5.1f ",m_dTwoTheta, m_dDistance); 
        
        for(int jj=0; jj < nGetNumAxes(); jj++)
        {
            printf("%6.1f",  dGetScanAxisValue(ii, jj) );

            ////////////////////////////////////////
            // Specify whether it is a rotation axis
            if( bIsRotationAxis(jj) )
                printf("*");
            else
                printf(" ");
        }
        ////////////////////////////////////////////////////////////////////////////
        
        dRotBegin = dGetScanRotBegin(ii, eScanType);
        dRotEnd   = dGetScanRotEnd  (ii, eScanType);
        dWidth    = fabs(dRotBegin - dRotEnd);
        
        dTotalRotWidth += dWidth;
        
        ///////////////////////////////////////////////////////////////////////////
        // Format Rotation Begin, End, Width, 2Theta and Distance
        if( dWidth > 0.0 )  
        {
            printf(" %6.1f %6.1f %6.1f ", dRotBegin, dRotEnd, dWidth);
        }
        else  // scan not tested?
        {
            printf("%6s %5s %6s %6s %6s",  "   ---", 
                                           "  ---",
                                           "   ---", 
                                           "   ---", 
                                           "   ---"); 
        }
        ///////////////////////////////////////////////////////////////////////////
        
        
        if( dWidth > 0.0 )
        {
            printf("%6.1f %5.2f %4.2f", dGetScanCompleteness (ii, eScanType), 
                                        dGetScanRedundancy   (ii, eScanType), 
                                        dGetScanRedundancyDev(ii, eScanType));
        }
        else // scan not tested?
        {
            printf("%6s %5s %4s", "   ---", 
                                   "  ---",
                                    " ---");
        }
        ////////////////////////////////////////////////////////////////////////////////////

        printf("\n");
    }

    return (nUsedScanNumber+1);
}

int CstrategyScans::nPrintResultScanTable(int nScanNumberOffset, 
                                          bool bPrintHeader, 
                                          double& dTotalRotWidth,
										  std::vector<Cstring>& padInfoVec,
                                          DTREK_WORD wCtrl,
                                          REFLN_FILE_STRAT_INFO* pReflnFileStrategyInfo)
{
	int nScansUsed = nPrintResultScanTable(nScanNumberOffset, 
                                          bPrintHeader, 
                                          dTotalRotWidth,
                                          wCtrl,
                                          pReflnFileStrategyInfo);

	if(nScansUsed > 0)
	{

		CstrategyScan::enStrategyScanType eScanType = 
			!m_bMultipleScanSolution ? CstrategyScan::enSingleScan : CstrategyScan::enMultipleScan;

		bool        bUsed = false;

		//shijie 2008.03.17 : var for padding info
		char msgStr[2048];	//fixed length! ok to make this large because it is local.
		msgStr[0] = 0;
		CstrategyScan::enStrategyScanRangePaddingType enPadType;
		CstrategyScan::enStrategyScanRangePaddingType enShiftType;



		Cstring     strYesNo("");
		int         nScanNumberForRecord = -1;
		int         nUsedScanNumber = -1;
		for(int ii=0; ii < nGetNumScans(); ii++)                                 
		{
			bUsed = bIsScanUsed(ii, eScanType);

			// Multiple scan solution should only list used scans
			if( !bUsed && DTREK_STRATSCANS_PRST_ALL != wCtrl )
				continue;
	        
			nUsedScanNumber++;
			nScanNumberForRecord = nUsedScanNumber + nScanNumberOffset;

	  
			enPadType = enGetPaddingType(ii, eScanType);
			enShiftType = enGetShiftingType(ii, eScanType);

			if(enShiftType != CstrategyScan::enNone)	//shifted
			{
				sprintf(msgStr, "%d : range shifted to start at whole angle value with integral number of images;",
								nScanNumberForRecord);
			}
			else if(enPadType != CstrategyScan::enNone)
			{
				if(enPadType == CstrategyScan::enFrontExtended)
				{
					sprintf(msgStr, "%d : range starting angle extended to contain integral number of images;",
									nScanNumberForRecord);
				}
				else if(enPadType == CstrategyScan::enEndExtended)
				{
					sprintf(msgStr, "%d : range ending angle extended to contain integral number of images;",
									nScanNumberForRecord);
				}
				else if(enPadType == CstrategyScan::enBothEndExtended)
				{
					sprintf(msgStr, "%d : range angles extended to contain integral number of images;",
									nScanNumberForRecord);
				}
				else	//both ends shrinked 
				{
					sprintf(msgStr, "%d : range angles shrinked to contain integral number of images;",
									nScanNumberForRecord);
				}
			}
            else
                msgStr[0] = 0;

			Cstring tmpStr = msgStr;
			if(tmpStr.length() > 0)
				padInfoVec.push_back(tmpStr);
		}
	}

    return nScansUsed;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CstrategyScans::CstrategyScans(Cimage_header* poHeader, const Cstring& sStrategyPrefix) :
m_seqResolution(-1.0,-1.0)
{
    m_dTwoTheta = cdStratDetUnrealisticValue;
    m_dDistance = cdStratDetUnrealisticValue;

    m_dTwoThetaCollisionOffset = cdStratDetUnrealisticValue;
    m_dDistanceCollisionOffset = cdStratDetUnrealisticValue;

    m_bMultipleScanSolution = false;

    bSetAxisInfo(poHeader);
    
    bSetScans(poHeader, sStrategyPrefix);

    bSetDetector(poHeader);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CstrategyScans::bSetAxisInfo(Cimage_header* poHeader)
{
    if( !poHeader )
        return false; // safety

    vClearAxes(); // initialize

    // How many crystal gonio values are there?
    Cstring         strCrystalPrefix(D_K_CrystalPrefix);

    int             nAxes = 0;
    if( 0 != poHeader->nGetValue(strCrystalPrefix + D_K_GonioNumValues, &nAxes) )
        return false;

    if( nAxes <= 0 )
        return false;

    // Get axes names
    Cstring*        psGonioNames = new Cstring[nAxes];
    if( 0 != poHeader->nGetValue(strCrystalPrefix + D_K_GonioNames, nAxes, psGonioNames) )
    {
        delete [] psGonioNames;
        return false;
    }

    // Get rotation axis name
    Cstring         strScanPrefix(D_K_ScanPrefix);

    Cstring         strRotAxis("");
    if( 0 != poHeader->nGetValue(strScanPrefix + D_K_RotAxisName, &strRotAxis) )
        return false;

    int     ii = 0;
    // Get axes collision offsets
    double*         pdGonioCollisionOffsets = new double[nAxes];
    for(ii=0; ii < nAxes; ii++)
        pdGonioCollisionOffsets[ii] = 0.0;
    poHeader->nGetValue(strCrystalPrefix + D_K_GonioCollisionOffsets, nAxes, pdGonioCollisionOffsets);

    // Now add crystal axes objects 
    for(ii=0; ii < nAxes; ii++)
    {
        vAddAxis(psGonioNames[ii], strRotAxis==psGonioNames[ii], pdGonioCollisionOffsets[ii]);
    }

    delete [] psGonioNames;
    delete [] pdGonioCollisionOffsets;
    
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CstrategyScans::bSetScans(Cimage_header* poHeader, const Cstring& sStrategyPrefix)
{
    if( !poHeader )
        return false; // safety

    vClearScans();

    /////////////////////////////////////////////////////////////////////////////////
    // Figure out how many scans are there
    int     iScanMax = -1;
    while(true)
    {
        iScanMax++;
        Crotation     oRotation(*poHeader, sBuildStrategyScanPrefix(sStrategyPrefix, D_K_ScanPrefix, iScanMax));
        
        oRotation.ms_nVerbose = 0;  // so that when it fails at the end of this while loop, it does not say it's an error, because we do want it to fail.

        if( !oRotation.bIsAvailable() )
            break;
    }
    iScanMax--; // last valid scan index
    /////////////////////////////////////////////////////////////////////////////////

    if( iScanMax < 0 )
        return false;
    
    for(int iScan=0; iScan <= iScanMax; iScan++)
        vAddScan(poHeader, iScan, m_aoAxes.size(), sStrategyPrefix);

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////
void CstrategyScan::vGetAxes(std::vector<CstrategyAxis>& aoAxes)
{
    aoAxes.clear();

    for(int ii = 0; ii < (int)m_aoAxes.size(); ii++)
    {
        aoAxes.push_back(CstrategyAxis(m_aoAxes[ii]));
    }
}
//////////////////////////////////////////////////////////////////////////////////////
bool CstrategyScans::bSetDetector(Cimage_header* poHeader)
{
    if( !poHeader )
        return false; // safety

    Cstring              strDetectorPrefix(""); 
    if ( 0 != poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    Cgoniometer         oDetGonio(*poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
        return false;
    
    m_dDistance = (double)oDetGonio.fGetDistance();
    m_dTwoTheta = (double)oDetGonio.fGetSwing();
     
    // Get and set collision offsets
    float   fDistCollOffset = 0.0f;
    float   fTwoThetaCollOffset = 0.0f;

    const char*     c_pcTwoTheta = "2Theta";
    const char*     c_pcDistance = "Dist";

    if( oDetGonio.bGetCollisionOffset(c_pcDistance, fDistCollOffset) )
        m_dDistanceCollisionOffset = (double)fDistCollOffset;
    else
        m_dDistanceCollisionOffset = 0.0;

    if( oDetGonio.bGetCollisionOffset(c_pcTwoTheta, fTwoThetaCollOffset) )
        m_dTwoThetaCollisionOffset = (double)fTwoThetaCollOffset;
    else
        m_dTwoThetaCollisionOffset = 0.0;

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CstrategyScans::bSetScanAxisValue(int iScan, const char* ccAxisName, double dValue)
{
    Cstring     strInputAxisName(ccAxisName);
    strInputAxisName.upcase();

    
    int             iAxis = -1;
    Cstring         strName("");
    
    for(iAxis=0; iAxis < (int)m_aoAxes.size(); iAxis++)
    {
        m_aoAxes[iAxis].vGetName(strName);

        strName.upcase();
        
        if( strName == strInputAxisName )
            break;  // axis name found!
    }
    
    return m_aoScans[iScan].bSetAxisValue(iAxis, dValue);   // if iAxis is ivalid, the return will be false
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CstrategyScans::vSortScans()
{
    std::sort(m_aoScans.begin(), m_aoScans.end(), bScansCompare); 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CstrategyScans::vSetImageExposureTime(double fTime)
{
    for(int iScan=0; iScan < nGetNumScans(); iScan++)
        m_aoScans[iScan].vSetImageExposureTime(fTime);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static bool bScansCompare(const CstrategyScan& first, const CstrategyScan& second)
{
    if( !(first.bIsUsed(CstrategyScan::enMultipleScan) && second.bIsUsed(CstrategyScan::enMultipleScan)) )
    {
        if( first.bIsUsed(CstrategyScan::enMultipleScan) )
            return true; // 1st used, 2nd - not used
        else if( second.bIsUsed(CstrategyScan::enMultipleScan) )
            return false; // 1st not used, 2nd - used
    }

    // If we made it this far, it means that either both or none of the two scans are used
    if( first.dGetCompleteness(CstrategyScan::enMultipleScan) > second.dGetCompleteness(CstrategyScan::enMultipleScan) )
        return true;
    else if( first.dGetCompleteness(CstrategyScan::enMultipleScan) == second.dGetCompleteness(CstrategyScan::enMultipleScan) )
    {
        if( first.dGetRotWidth(CstrategyScan::enMultipleScan) > second.dGetRotWidth(CstrategyScan::enMultipleScan) )
            return true;
        else if( first.dGetRotWidth(CstrategyScan::enMultipleScan) == second.dGetRotWidth(CstrategyScan::enMultipleScan) )
        {
            if( first.dGetRedundancy(CstrategyScan::enMultipleScan) > second.dGetRedundancy(CstrategyScan::enMultipleScan) )
                return true;
        }
    }
    
    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


