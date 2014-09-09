#ifndef DT_CSTRATEGY_H
#define DT_CSTRATEGY_H
//
// Copyright (c) 1997 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cstrategy.h        Initial author: J.W. Pflugrath           10-May-1997
//    This file is the header file for class Cstrategy
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

//+Include files

#include "Dtrek.h"
#include "Cimage_header.h"
#include "Cdetector.h"
#include "Csource.h"
#include "Ccrystal.h"
#include "Crotation.h"
#include "Creflnlist.h"
#include "dtrekvec.h"
#include "dtarray.h"
#include "ResoBins.h"
#include "StrategyScans.h"

//+Definitions and constants
const int nMAX_TEST_STRATEGY = 5;
const int nMAX_REDUNDANCY = 50;
class Cstrategy;

//+Code begin
class DTREK_EXPORT CstrategyRange
{
public:
  itr<int>      m_anScanNumbers;
  itr<double>   m_afRotStart;
  itr<double>   m_afRotEnd;

  double        m_fRedundancy;
  double        m_fRedundancyDev;
  double        m_fCompleteness;

  bool          m_bIsAvailable;

private:
  Cstrategy*    m_poStrategy;
public:

  int nPrint(Cstrategy& oStrategy);
  int nUpdateHeader(Cimage_header& oHeader, Cstrategy& oStrategy);    // Assumes oHeader contains scans to start with.
  bool operator < (CstrategyRange& oOther);
  int operator -() { -m_anScanNumbers; -m_afRotStart; -m_afRotEnd; return 0; };
  int operator += (CstrategyRange& oOther) { m_anScanNumbers += oOther.m_anScanNumbers; m_afRotStart += oOther.m_afRotStart; m_afRotEnd += oOther.m_afRotEnd; return 0; };
  
  CstrategyRange(const CstrategyRange& oOther)
  {
    m_anScanNumbers += oOther.m_anScanNumbers; 
    m_afRotStart    += oOther.m_afRotStart; 
    m_afRotEnd      += oOther.m_afRotEnd;
      
    m_fRedundancy    = oOther.m_fRedundancy;    
    m_fRedundancyDev = oOther.m_fRedundancyDev; 
    m_fCompleteness  = oOther.m_fCompleteness;  
    m_bIsAvailable   = oOther.m_bIsAvailable;
    m_poStrategy     = oOther.m_poStrategy;
  }
  
  CstrategyRange& operator= (const CstrategyRange& oOther)
  {
    -m_anScanNumbers;
    -m_afRotStart;   
    -m_afRotEnd;     
    
    m_anScanNumbers += oOther.m_anScanNumbers; 
    m_afRotStart    += oOther.m_afRotStart; 
    m_afRotEnd      += oOther.m_afRotEnd;
      
    m_fRedundancy    = oOther.m_fRedundancy;    
    m_fRedundancyDev = oOther.m_fRedundancyDev; 
    m_fCompleteness  = oOther.m_fCompleteness;  
    m_bIsAvailable   = oOther.m_bIsAvailable;
    m_poStrategy     = oOther.m_poStrategy;

    return *this;
  }

  CstrategyRange(Cstrategy* poStrategy) :
  m_poStrategy(poStrategy)
  { 
      m_fRedundancy = 0.0;     
      m_fRedundancyDev = 0.0;  
      m_fCompleteness = 0.0;   
      m_bIsAvailable = FALSE;
  }

  CstrategyRange() :
  m_poStrategy(NULL)
  { 
      m_fRedundancy = 0.0;     
      m_fRedundancyDev = 0.0;  
      m_fCompleteness = 0.0;   
      m_bIsAvailable = FALSE;
  }
  
  void vSetStrategyPtr(Cstrategy* poStrategy){m_poStrategy=poStrategy;}

  double dGetTotalRotWidth();
  
  double dGetScanRotBegin(int iScan){return m_afRotStart[iScan];}
  double dGetScanRotEnd(int iScan){return m_afRotEnd[iScan];}

  void vSetScanRotBeg(int iScan, double dVal){m_afRotStart[iScan]=dVal;}
  void vSetScanRotEnd(int iScan, double dVal){m_afRotEnd[iScan]=dVal;}

  // "Original index" refers to the index a scan had when we just inputted all scans.
  int nGetScanOriginalIndex(int iScan)
  {
      if(iScan < 0 || iScan > m_anScanNumbers.size() - 1)
          return -1;
      else
          return m_anScanNumbers[iScan];
  }
  
  int nGetNumberOfScans(){return m_anScanNumbers.size();}
  
  void vGetScansList(Cstring& strScans, bool bIncludeFixedScans=false);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT Cstrategy 
{
  friend class CstrategyRange;
public:
  static Cstring    ms_sf2STLcu;

private:
  bool      m_bFixRotStart;
  
  int       m_nVerbose;           // Verbosity flag
  
  int       m_nNumResoBins;       // Number of resolution bins for tables
  
  double    m_fCompletenessMin;   // Minimum completeness desired (i.e. target completeness)
  double    m_fCompletenessTol;   // A maximum pad value to add to the completeness to help it reach the target
  double    m_fRotationTol;       // A minimum percentage rotation increase allowed to try and reach the target

  double    m_dMinCompletenessIncreaseOnExpansion; // A minimum completeness increase required 
                                                   // to expand a range trying to reach the target

  double    m_fMinCompPerDegree;  // Minimum increase in percent completeness per degree.  Used only in single scan computations.
  
  double    m_fRedundancyMin;     // Minimum redundancy allowed
  double    m_fCompletenessTolOnFailure;
  double    m_fRedundancyTolOnFailure;
  int       m_a2nFirstScan[2];    // First scan(s) to choose from.
  int       m_a2nSecondScan[2];   // Second scan(s) to choose from.
  itr<int>  m_anFixedScans;       // Always choose reflections from these scans.
  itr<int>  m_anExcludeScans;     // Do not choose reflections from these scans.  This is usefull for high-redundancy strategy.
  bool      m_bAnom;              // Assume I+ != I- during statistical calculations.
  int       m_nMaxTime;           // Used for multiple scan algorithm.
  double    m_fPad;               // Padding.  Used to pad the predictions.

  double    m_fCellLengthFactor;  // Cell length factor for speed
  double    m_fCellVolumeFactor;  // Cell volume divisor for speed

  itr<int>      m_anScanStart;    // In the total predicted reflection list, sorted on scans, we find indices corresponding to the beginning of a new scan. We store these indices in this array.
  itr<int>      m_anScanEnd;      // In the total predicted reflection list, sorted on scans, we find indices corresponding to the end of a new scan. We store these indices in this array.

  itr<double>   m_afScanRotMin;  // Minimum rotation for each scan.
  itr<double>   m_afScanRotMax;  // Maximum rotation for each scan.
  itr<double>   m_afScanRotRange;// Range for each scan.
  itr<double>   m_afScanGonioAngles;    // FOR PRINTING ONLY (see nRememberScanGonioSettings(nScan))
  std::vector<Cstring>  m_asScanGonioAxis;      // FOR PRINTING ONLY (see nRememberScanGonioSettings(nScan))
  
  Ccrystal    *m_poCrystal;      // Pointer to crystal object
  Creflnlist  *m_poReflnlistIn;  // Pointer to input reflnlist object
  Creflnlist  *m_poReflnlist;    // Pointer to working reflnlist object
  int         *m_pnIndex;        // (Packed HKL,Scan)
  int         *m_pnIndexA;       // (Scan,Packed HKL)
  int         *m_pnIndex2;       // (Packed HKL,Scan) (assume triclinic, anom)

 
  int               m_nFI_fRotStart;
  int               m_nFI_fRotEnd;
  
  int               m_nFI_f2STLcu;
  
  int               m_nFI_nPackedHKL;
  int               m_nFI_nPackedHKL2;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  CstrategyScans*                   m_poScanInfo;
  bool                              m_bExternalScanInfoPtr;
  double                            m_dTotalRotWidthAllowed;
  int                               m_nTotalNumberOfScansAllowed;

  double            m_dTotalRotWidth;

  CResoStrategyBins*        m_poResoBins;
  bool              m_bExternalResoBinsPtr;

private:
#define DTREK_DTSTRATEGY_FRTSSA_SET_USED                        0x0001
#define DTREK_DTSTRATEGY_FRTSSA_UPDATE_ROT_RANGE                0x0002
#define DTREK_DTSTRATEGY_FRTSSA_UPDATE_STATS                    0x0004
#define DTREK_DTSTRATEGY_FRTSSA_UPDATE_ALL_RESULTS              (DTREK_DTSTRATEGY_FRTSSA_UPDATE_ROT_RANGE | DTREK_DTSTRATEGY_FRTSSA_UPDATE_STATS)
//#define DTREK_DTSTRATEGY_FRTSSA_RESET_USED                      0x0008

  int  nFromRangeToStrategyScansArray(CstrategyRange& oRange, DTREK_WORD wCtrl=0U);
  bool nFromStrategyScansArrayToScanRotationArrays(bool bCheckRotWidthOnly=true);

public:
  bool              m_bZeroOmegaWithPhiOffset;
  double            m_fZeroOmega;
  CstrategyRange    m_oBestRange;

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments


Cstrategy (Cimage_header *poHeader);

~Cstrategy ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes

public:
    bool bIsValid(Cstring& sError);
#define DTREK_STRATEGY_LC_MSG_USER_ACCEPTED           0x0001
#define DTREK_STRATEGY_LC_MSG_LOW_COMPLETENESS        0x0002
    int   nListCompleteness(CstrategyRange* poRange, DTREK_WORD wCtrl);

    int   nClearRange(CstrategyRange& oRange);
    int   nCalcCompletenessSensitivity(CstrategyRange& oRange,bool bPrint);
    int   nSetup(Cimage_header& oHeader,Creflnlist& oList);
    int   nExpand(Creflnlist& oList);
    int   nUnPackRedundancy();
    int   nCalcScans();
    int   nCalcScanSingle(CstrategyRange& oTestRange,double fMaxSearchStep,int nIndexToChange = 0);
    int   nCalcScanMultiple(CstrategyRange& oTestRange);
    
    int   nAdjustPhiToZero(Cimage_header& oHeader,CstrategyRange& oTestRange);
    void  vShiftReflnRotValueInScan(int iScan, double dRotShift);

    int   nCalcChiCoverage(Cdetector& oDetector,Csource& oSource,double fResoLow,double fResoHigh,double fDebyeRingCoverage, bool bPrint = true);

    bool bScanUsed(int nScan,int nOrder = 0);
    bool bIsFixedScan(int nScan);

    void vRememberScanGonioSettings(Cgoniometer* poCrysGonio,Crotation* poRotation,int nScan);

    int  nUpdateHeader(Cimage_header* poHeader);

    inline void vSetNumResoBins(const int nBins) { m_nNumResoBins = nBins; }

    inline void vSetVerboseLevel(const int nLev) { m_nVerbose = nLev; }
    inline int  nGetVerboseLevel(void) { return (m_nVerbose); }
    inline void vSetCompleteness(const double fComplete) { m_fCompletenessMin = fComplete; }
    //inline void vSetCompletenessSensitivity(const double fMinCompPerDegree) { m_fMinCompPerDegree = fMinCompPerDegree; };
    inline void vSetCompletenessTolerance(const double fTol){m_fCompletenessTol = fTol;};
    inline void vSetRedundancy(const double fRedundancy) { m_fRedundancyMin = fRedundancy; }
    inline double fGetRedundancy()const{return m_fRedundancyMin;}
    inline void vSetCellScale(const double fScaleLength) { m_fCellLengthFactor = fScaleLength; }

    bool bLeaveBestRangeReflectionsOnly(bool bLeaveFixedScans=true);
    int  nWrite(const Cstring& sName = "");

    inline void vSetFirstScan(int nFirstScanStart,int nFirstScanEnd) { m_a2nFirstScan[0] = nFirstScanStart; m_a2nFirstScan[1] = nFirstScanEnd; };
    inline void vSetSecondScan(int nSecondScanStart,int nSecondScanEnd) { m_a2nSecondScan[0] = nSecondScanStart; m_a2nSecondScan[1] = nSecondScanEnd; };

    void vAddFixedScan(int nFixedScanStart, int nFixedScanEnd);

    inline void vSetAnom(bool bAnom) { m_bAnom = bAnom; };
    inline void vSetMaxTime(int nMaxTime) { m_nMaxTime = nMaxTime; };
    inline void vSetPad(double fPad) { m_fPad = fPad; };

    void vSetAxisInfo(Cimage_header* poHeader);
    void vSaveDetectorPosition(Cimage_header* poHeader, double dTwoThetaSwing, double dDistance);

    void vAddStrategyScan(Cimage_header* poHeader, int nScanIndex, int nNumAxes);

    int  nGetNumAxes()const;

    void vGetAxisName(int iIndex, Cstring& strAxis)const;
    bool bIsRotationAxis(int iIndex);

    void vPrintSummary(bool bPrintSingleScanSolution);
    void vPrintSummaryScanTable(bool bSingleScanSummary);
    void vSetFixRotStart(bool bFix){m_bFixRotStart=bFix;}

    double dGetTotalRotWidth(){return m_dTotalRotWidth;}

    void vSetExternalResoBinsPtr(CResoStrategyBins* pBins);

    void vSetExternalWorkingReflnListPtr(Creflnlist* pList);

    void vSetScanInfoObj(CstrategyScans* pScans=NULL);

    void vFillStats(CstrategyRange* poTestRange, bool bMustBeResolutionConsistent=true);

    void vSetMaximumTotalRotation(double dRot){m_dTotalRotWidthAllowed = dRot;}
    void vSetMaximumTotalScans(int nScans){m_nTotalNumberOfScansAllowed = nScans;}

    bool bIsMultipleScanSolution(){return m_poScanInfo->bIsMultipleScanSolution();}
    bool bTryExpandRange(CstrategyRange& oInputRange, double dExtraRotValue);
    void vSortFindBestExpandableRange(std::vector<CstrategyRange>& vRanges);

    double fGetCompletenessTolerance()const{return m_fCompletenessTol;}

#define DTREK_CSTRATEGY_LIST_TITLE	        0x0001
#define DTREK_CSTRATEGY_LIST_CRYSTAL        0x0002
#define DTREK_CSTRATEGY_LIST_REFLN_LIST	    0x0004	 
#define DTREK_CSTRATEGY_LIST_RESO	        0x0008
    int nList(DTREK_WORD wCtrl);
    void vReflnListCalcSetDStarCubed(double& dStarCubedMin, double& dStarCubedMax);
    void vReflnListSortReduce();

private: // Some of the above function will be made private

    int  nInitValues(void);
    int  nInitValues(Cimage_header& oHeader);

    int  nCompleteness(CstrategyRange& oTestRange);
    int  nCompletenessMultiScan(CstrategyRange& oTestRange);
    int  nComputeBestScanToAdd(itr<int>& anScansInUse);
    bool bSetResoBinsObject(double dStarCubedMin, double dStarCubedMax);
};
// end of class Cstrategy
#endif   // DT_CSTRATEGY_H
