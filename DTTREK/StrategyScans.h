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
// Strategyscans.h        Initial author: RB     06-May-2005
// This file contains the definitions of class CStrategyScans

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

#ifndef DT_STRATEGY_SCANS_H
#define DT_STRATEGY_SCANS_H

#include <vector>
#include "dtrekdefs.h"
#include "Cstring.h"
#include "Segment.h"

const double cdStratDetUnrealisticValue = -999.0;

typedef struct _tagREFLN_FILE_STRAT_INFO
{
    Cstring  m_sFilePath;
    double   m_fCompleteness; 
    double   m_fRedundancy; 
    double   m_fRedundancyDev;
    _tagREFLN_FILE_STRAT_INFO(Cstring& sFilePath){m_sFilePath=sFilePath; m_fCompleteness=0.0; m_fRedundancy=0.0; m_fRedundancyDev = 0.0;}
}REFLN_FILE_STRAT_INFO;

class Cimage_header;

class DTREK_EXPORT CstrategyAxisInfo
{
public:
    CstrategyAxisInfo(){m_strName=""; m_bRotAxis=false; m_dCollisionOffset=0.0;}
    CstrategyAxisInfo(Cstring strName, bool bRot, double dCollisionOffset){m_strName=strName;
                                                                           m_bRotAxis=bRot;
                                                                           m_dCollisionOffset=dCollisionOffset;}
    CstrategyAxisInfo& operator=(const CstrategyAxisInfo& anotherAxisInfo);
    CstrategyAxisInfo(const CstrategyAxisInfo& anotherAxisInfo);
public:
    bool bIsRotationAxis()const{return m_bRotAxis;}
    void vGetName(Cstring& strName)const{strName=m_strName;}
    
    double dGetCollisionOffset(){return m_dCollisionOffset;}
private:
    Cstring     m_strName;
    bool        m_bRotAxis;

    double      m_dCollisionOffset;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CstrategyAxis
{
public:
    CstrategyAxis(){m_nInfoIndex=-1; m_dValue=0.0;}
    CstrategyAxis(int nInfoIndex, double dValue){m_nInfoIndex=nInfoIndex; m_dValue=dValue;}
    CstrategyAxis(const CstrategyAxis& anotherAxis);
public:
    double dGetValue()const{return m_dValue;}
    int nGetInfoIndex()const{return m_nInfoIndex;}
    
    void vSetValue(double dValue){m_dValue = dValue;}

    CstrategyAxis& operator=(const CstrategyAxis& anotherAxis);

private:
    int         m_nInfoIndex;   // a reference to the information about this axis in the AxisInfo
    double      m_dValue;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//RB Class CstrategyScan is used to keep information about scans, which is needed for reporting.
class DTREK_EXPORT CstrategyScan
{
public:
    CstrategyScan();

    CstrategyScan(const CstrategyScan& anotherScan);
    
    CstrategyScan(Cimage_header* poHeader, int nScanIndex, int nNumAxes, const Cstring& sStrategyPrefix="");

    ~CstrategyScan(){}

    CstrategyScan& operator=(const CstrategyScan& anotherScan);
    
    enum enStrategyScanType
    {
        enSingleScan=0,
        enMultipleScan
    };

	//shijie yao
	enum enStrategyScanRangePaddingType
	{
		enNone = 0,
		enBothEndExtended,
		enFrontExtended,
		enEndExtended,
		enBothEndShrinked,
		enShiftedToFront,
		enShiftedToEnd,
		enShiftedToEndAndChop
	};
    
    void vSetUsed(const bool bUsed, enStrategyScanType eType){m_bUsed[eType]=bUsed;}
    bool bIsUsed(enStrategyScanType eType)const{return m_bUsed[eType];}

    //void vSetRotation(const double dBegin, const double dEnd, enStrategyScanType eType){m_dRot_Begin[eType]=dBegin; m_dRot_End[eType]=dEnd;}
    //shijie yao : 2008.03.12
	void vSetRotation(const double dBegin, const double dEnd, enStrategyScanType eType);

    double dGetRot_Begin(enStrategyScanType eType)const{return m_dRot_Begin[eType];}
    double dGetRot_End(enStrategyScanType eType)const{return m_dRot_End[eType];}

    double dGetRotWidth(enStrategyScanType eType)const{return fabs(m_dRot_End[eType] - m_dRot_Begin[eType]);}

    int  nGetNumAxes()const{return (int)m_aoAxes.size();}
    double dGetAxisValue(int iIndex);
    bool bSetAxisValue(int iIndex, double dValue);

    void vSetCompleteness (const double dCompleteness, enStrategyScanType eType)  {m_dCompleteness[eType] = dCompleteness;}
    void vSetRedundancy   (const double dRedundancy, enStrategyScanType eType)    {m_dRedundancy[eType] = dRedundancy;}
    void vSetRedundancyDev(const double dRedundancyDev, enStrategyScanType eType) {m_dRedundancyDev[eType] = dRedundancyDev;}

    void vSetRotLimits(double dBegin, double dEnd){m_segRotLimits.vSet(dBegin, dEnd);}
    
    void vSetImageExposureTime(double fTime){m_fImageExpTime=fTime;}

    double dGetCompleteness(enStrategyScanType eType)const{return m_dCompleteness[eType];}
    double dGetRedundancy(enStrategyScanType eType)const{return m_dRedundancy[eType];}
    double dGetRedundancyDev(enStrategyScanType eType)const{return m_dRedundancyDev[eType];}
    
    void   vGetAxes(std::vector<CstrategyAxis>& aoAxes);
    void   vGetRotLimits(CSegment& segLimits){segLimits=m_segRotLimits;}
    double dGetMaxRotRange()const{return m_dMaxRotRange;}

    int    nGetOriginalScanIndex()const{return m_nOriginalScanIndex;}
    
    enum enStrategyScanDetRelZero
    {
        enStratScanDetRelZero_TwoTheta = 0,
        enStratScanDetRelZero_Distance = 1
    };
    double dGetDetRelZero(enStrategyScanDetRelZero eType){return m_dDetRelZero[eType];}
    
    double fGetImageExposureTime()const{return  m_fImageExpTime;}
    double fGetImageRotationWidth()const{return m_fImageRotWidth;}

private:

    std::vector<CstrategyAxis>  m_aoAxes;

    double      m_dRot_Begin[2];       // The first value is used when a scan is a single scan solution, the second is used when a scan is part of a multiple scan solution
    double      m_dRot_End[2];

    CSegment    m_segRotLimits;
    double      m_dMaxRotRange;

    double      m_dCompleteness[2];
    double      m_dRedundancy[2];
    double      m_dRedundancyDev[2];

    bool        m_bUsed[2];
    
    int         m_nOriginalScanIndex;

    double      m_dDetRelZero[2];

    double      m_fImageRotWidth;
    double      m_fImageExpTime;

    void vInit();

	//shijie yao : 2008.03.12
private:
	double	m_dCryGonioValMax;	//CRYSTAL_GONIO_VALUES_MAX
	double	m_dCryGonioValMin;	//CRYSTAL_GONIO_VALUES_MIN
public:
	double dGetCryGonioValMax()		const { return m_dCryGonioValMax;	};
	double dGetCryGonioValMin()		const { return m_dCryGonioValMin;	};
	enStrategyScanRangePaddingType	m_enPaddedType[2];
	double							m_dPaddedVal[2];
	enStrategyScanRangePaddingType	m_enShiftedType[2];
	double							m_dShiftedVal[2];
private:
	void vPaddingRange(enStrategyScanType type);



};
////////////////////////////////////////////////////////////////////////////////////////

class DTREK_EXPORT CstrategyScans
{
public:
    CstrategyScans();
    
    CstrategyScans(Cimage_header* poHeader, const Cstring& sStrategyPrefix="");
    
    ~CstrategyScans();
public:
    int     nGetNumAxes()const{return (int)m_aoAxes.size();}
    void    vClearAxes(){m_aoAxes.clear();}
    void    vAddAxis(Cstring& sName, bool bRot, double dCollisionOffset){m_aoAxes.push_back(CstrategyAxisInfo(sName, bRot, dCollisionOffset));}
    void    vClearScans(){m_aoScans.clear();}
    
    void    vAddScan(Cimage_header* poHeader, 
                     int nScanIndex, 
                     int nNumAxes,
                     const Cstring& sStrategyPrefix="")
                     {m_aoScans.push_back(CstrategyScan(poHeader, nScanIndex, nNumAxes, sStrategyPrefix));}
    
    void    vDeleteScan(int iIndex);
    
    int     nGetNumScans(){return (int)m_aoScans.size();}
    
    void    vGetAxisName(int iAxis, Cstring& strAxisName){m_aoAxes[iAxis].vGetName(strAxisName);}
    bool    bIsRotationAxis(int iAxis){return m_aoAxes[iAxis].bIsRotationAxis();}

    bool    bIsScanUsed(int iScan, CstrategyScan::enStrategyScanType eScanType){return m_aoScans[iScan].bIsUsed(eScanType);}
    bool    bIsScanUsed(int iScan){return m_aoScans[iScan].bIsUsed(m_bMultipleScanSolution ? CstrategyScan::enMultipleScan : CstrategyScan::enSingleScan);}
    
    void    vSetScanUsed(int iScan, CstrategyScan::enStrategyScanType eScanType, bool bUsed){m_aoScans[iScan].vSetUsed(bUsed, eScanType);}
    
    double  dGetScanAxisValue(int iScan, int iAxis){return m_aoScans[iScan].dGetAxisValue(iAxis);}
    bool    bSetScanAxisValue(int iScan, const char* ccAxisName, double dValue);

    double  dGetAxisCollisionOffset(int iAxis){return m_aoAxes[iAxis].dGetCollisionOffset();}

    double  dGetScanCompleteness (int iScan, CstrategyScan::enStrategyScanType eScanType){return m_aoScans[iScan].dGetCompleteness (eScanType);} 
    void    vSetScanCompleteness (int iScan, CstrategyScan::enStrategyScanType eScanType, double dCompl)
                                 {m_aoScans[iScan].vSetCompleteness(dCompl, eScanType);} 
    
    double  dGetScanRedundancy   (int iScan, CstrategyScan::enStrategyScanType eScanType){return m_aoScans[iScan].dGetRedundancy   (eScanType);} 
    void    vSetScanRedundancy   (int iScan, CstrategyScan::enStrategyScanType eScanType, double dRed)
                                 {m_aoScans[iScan].vSetRedundancy(dRed, eScanType);} 
    
    double  dGetScanRedundancyDev(int iScan, CstrategyScan::enStrategyScanType eScanType){return m_aoScans[iScan].dGetRedundancyDev(eScanType);} 
    void    vSetScanRedundancyDev(int iScan, CstrategyScan::enStrategyScanType eScanType, double dDev)
                                 {m_aoScans[iScan].vSetRedundancyDev(dDev, eScanType);} 
    
    double  dGetScanRotBegin(int iScan, CstrategyScan::enStrategyScanType eScanType){return m_aoScans[iScan].dGetRot_Begin(eScanType);}
    double  dGetScanRotEnd  (int iScan, CstrategyScan::enStrategyScanType eScanType){return m_aoScans[iScan].dGetRot_End  (eScanType);}
    
    double  dGetScanRotBegin(int iScan){return m_aoScans[iScan].dGetRot_Begin(m_bMultipleScanSolution ? CstrategyScan::enMultipleScan : CstrategyScan::enSingleScan);}
    double  dGetScanRotEnd  (int iScan){return m_aoScans[iScan].dGetRot_End  (m_bMultipleScanSolution ? CstrategyScan::enMultipleScan : CstrategyScan::enSingleScan);}
    
    void    vSetScanRotation(int iScan, CstrategyScan::enStrategyScanType eScanType, double dStart, double dEnd)
                            {m_aoScans[iScan].vSetRotation(dStart, dEnd, eScanType);}
    
    void    vSetScanRotLimits(int iScan, double dBegin, double dEnd){m_aoScans[iScan].vSetRotLimits(dBegin, dEnd);}
    void    vSetScanRotLimits(int iScan, CSegment segLimits){m_aoScans[iScan].vSetRotLimits(segLimits.dGetBeg(), segLimits.dGetEnd());}
    
    void    vGetScanRotLimits(int iScan, CSegment& segLimits){m_aoScans[iScan].vGetRotLimits(segLimits);}

    double  dGetScanMaxRotRange(int iScan){return m_aoScans[iScan].dGetMaxRotRange();}

    int     nUpdateHeader(Cimage_header* pHeader, 
                          int nScanNumberOffset=0, 
                          const Cstring& sStrategyPrefix="",
                          bool bUpdateDetectorPosition=true);

    void    vSetMultipleScanSolution(bool bIsMSS){m_bMultipleScanSolution=bIsMSS;}
    bool    bIsMultipleScanSolution(){return m_bMultipleScanSolution;}

    double  dGetTwoTheta()const{return m_dTwoTheta;}
    double  dGetDistance()const{return m_dDistance;}
    
    void    vSetTwoTheta(const double dTwoTheta){m_dTwoTheta=dTwoTheta;}
    void    vSetDistance(const double dDistance){m_dDistance=dDistance;}

    double  dGetTwoThetaCollisionOffset()const{return m_dTwoThetaCollisionOffset;}
    double  dGetDistanceCollisionOffset()const{return m_dDistanceCollisionOffset;}
    
    void    vSetTwoThetaCollisionOffset(const double dTwoThetaCollisionOffset){m_dTwoThetaCollisionOffset=dTwoThetaCollisionOffset;}
    void    vSetDistanceCollisionOffset(const double dDistanceCollisionOffset){m_dDistanceCollisionOffset=dDistanceCollisionOffset;}
    
    double  dGetScanDetRelZero(int iScan, CstrategyScan::enStrategyScanDetRelZero eType){return m_aoScans[iScan].dGetDetRelZero(eType);}

#define  DTREK_STRATSCANS_PRST_ALL    0x0001
    int     nPrintResultScanTable(int nScanNumberOffset, 
                                  bool bPrintHeader, 
                                  double& dTotalRotWidth, 
                                  DTREK_WORD wCtrl=0U,
                                  REFLN_FILE_STRAT_INFO* pReflnFileStrategyInfo=NULL);

	//shijie : 2008.03.17 also returns the padding info in the padInfoVec.
    int     nPrintResultScanTable(int nScanNumberOffset, 
                                  bool bPrintHeader, 
                                  double& dTotalRotWidth, 
								  std::vector<Cstring>& padInfoVec,
                                  DTREK_WORD wCtrl=0U,
                                  REFLN_FILE_STRAT_INFO* pReflnFileStrategyInfo=NULL);

    void    vGetAxes(int iScan, std::vector<CstrategyAxis>& aoAxes){m_aoScans[iScan].vGetAxes(aoAxes);}

    void    vSortScans();

    int    nGetOriginalScanIndex(int iScan)const{return m_aoScans[iScan].nGetOriginalScanIndex();}

    void   vSetResoRange(double dResoMin, double dResoMax){m_seqResolution=CSegment(dResoMin,dResoMax);}
    double  dGetResoMax()const{return m_seqResolution.dGetBeg();}
    double  dGetResoMin()const{return m_seqResolution.dGetEnd();}

    void vSetImageExposureTime(double fTime);


	//shijie 2008.03.17
	CstrategyScan::enStrategyScanRangePaddingType 
		enGetPaddingType(int iScan, CstrategyScan::enStrategyScanType eScanType) { return m_aoScans[iScan].m_enPaddedType[eScanType]; };

	CstrategyScan::enStrategyScanRangePaddingType 
		enGetShiftingType(int iScan, CstrategyScan::enStrategyScanType eScanType) { return m_aoScans[iScan].m_enShiftedType[eScanType]; };

private:
    bool bSetAxisInfo(Cimage_header* poHeader);
    bool bSetScans(Cimage_header* poHeader, const Cstring& sStrategyPrefix="");
    bool bSetDetector(Cimage_header* poHeader);

private:
    std::vector<CstrategyScan>      m_aoScans;
    std::vector<CstrategyAxisInfo>  m_aoAxes;

    double    m_dTwoTheta;
    double    m_dDistance;

    double    m_dTwoThetaCollisionOffset;
    double    m_dDistanceCollisionOffset;

    bool m_bMultipleScanSolution;

    CSegment  m_seqResolution;
};
#endif   // !DT_STRATEGY_SCANS_H
