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
// ResoBins.h       Initial author: RB               04-Apr-2005

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

#ifndef DT_RESO_BINS_H
#define DT_RESO_BINS_H

#include "Dtrek.h"
#include "dtreksys.h"
#include "Cstat.h"

#include "Segment.h"

class Ccrystal;

class DTREK_EXPORT CResoBin
{

public:
    CResoBin();
    CResoBin(double dLow, double dHigh);
    CResoBin(const CResoBin& anotherBin);
    
    CResoBin& operator=(const CResoBin& anotherBin);

    virtual ~CResoBin(){}
    
    double dGetLowReso() {return m_dLow;}
    double dGetHighReso(){return m_dHigh;}
    double fGetMiddleReso();
    
protected:
    double      m_dLow;
    double      m_dHigh;
};


///////////////////////////////////////////////////////////////////////////////////////////////
#define     DTREK_RESOBIN_LIMITS                        0x0001

#define     DTREK_RESOBIN_SCALE_THEORETICAL_COUNT       0x0002

#define     DTREK_RESOBIN_SCALE_EXPERIMENTAL_COUNT      0x0004

#define     DTREK_RESOBIN_ALL_COUNT (DTREK_RESOBIN_SCALE_THEORETICAL_COUNT | DTREK_RESOBIN_SCALE_EXPERIMENTAL_COUNT)   

#define     DTREK_RESOBIN_ALL (DTREK_RESOBIN_LIMITS | DTREK_RESOBIN_ALL_COUNT)   

struct tagReflectionStats
{
  int   m_nNumRefls;               // Number of reflns (nNumUsed+nNumRejects = nNumRefs)
  int   m_nNumRejs;                // Number of rejected reflections.
  int   m_nNumUsed;                // Number of NON rejected reflections.
  int   m_nNumMults;               // Number of redundancy > 2 groups (minus number of singles)
  int   m_nNumSingles;             // Number of singles (redundancy==1)
  
  int   m_nNumUniqueCalc;          // Number of computed unique reflections.

  void vInit(DTREK_WORD wCtrl=DTREK_RESOBIN_ALL_COUNT);

  void vScaleReflnCount(double dScale, DTREK_WORD wCtrl);

  tagReflectionStats();
  tagReflectionStats(const tagReflectionStats& anotherStats);
  tagReflectionStats& operator += (const tagReflectionStats& anotherStats);
};
///////////////////////////////////////////////////////////////////////////////////////////////

class DTREK_EXPORT CResoStrategyBin : public CResoBin
{
public:
    CResoStrategyBin();
    CResoStrategyBin(double dMinReso, double dMaxReso);
    ~CResoStrategyBin();
    
    CResoStrategyBin(const CResoStrategyBin& anotherBin);
    CResoStrategyBin& operator=(const CResoStrategyBin& anotherBin);

private:

    tagReflectionStats      m_stStats;

public:
    
    void vInit(DTREK_WORD wCtrl=DTREK_RESOBIN_ALL);

    void vAddRefl()         {m_stStats.m_nNumRefls++;}
    void vAddRej()          {m_stStats.m_nNumRejs++;}
    void vAddUsed()         {m_stStats.m_nNumUsed++;}
    void vAddMult()         {m_stStats.m_nNumMults++;}
    void vAddSingle()       {m_stStats.m_nNumSingles++;}
    
    void vAddUniqueCalc() {m_stStats.m_nNumUniqueCalc++;}

    void vSetUniqueCalc(int nCalc){m_stStats.m_nNumUniqueCalc = nCalc;}

    void vAddStats(const tagReflectionStats& anotherStats){m_stStats += anotherStats;}
    tagReflectionStats*    pGetStats(){return &m_stStats;}

    int nGetNumRefls()     {return m_stStats.m_nNumRefls;}
    int nGetNumRejs()      {return m_stStats.m_nNumRejs;}
    int nGetNumUsed()      {return m_stStats.m_nNumUsed;}
    int nGetNumMults()     {return m_stStats.m_nNumMults;}
    int nGetNumSingles()   {return m_stStats.m_nNumSingles;}

    int nGetNumUniqueFound(){return (m_stStats.m_nNumSingles + m_stStats.m_nNumMults);}

    int nGetNumUniqueCalc()  {return m_stStats.m_nNumUniqueCalc;}

    void vScaleReflnCount(double dScale, DTREK_WORD wCtrl){m_stStats.vScaleReflnCount(dScale, wCtrl);}
private:
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CResoRankingBin : public CResoBin
{
public:
    CResoRankingBin()
    {
        m_nMinNumberOfReflnsForStats = 0;
        m_fIOverSigmaMean = -1.0;
        m_fIOverSigmaStandardDeviation = -1.0;
        m_fIOverSigmaStandardDeviationWRTRotation = -1.0;
        m_fIntensityDistributionMean = -1.0;
        m_fIntensityDistributionSigma = -1.0;
        m_fIntensityDistributionMin = -1.0;
        m_fIntensityDistributionMax = -1.0;
        m_nIntensityDistributionReflnCount = -1;

        m_fIOverSigmaDistributionMean = -1.0;
        m_fIOverSigmaDistributionSigma = -1.0;

        m_oIntensityStats.vClearRaw();
        m_oIOverSigmaStats.vClearRaw();
    }
    
    CResoRankingBin(double dMinReso, double dMaxReso) : CResoBin(dMinReso, dMaxReso)
    {
        m_nMinNumberOfReflnsForStats = 0;
        m_fIOverSigmaMean = -1.0;
        m_fIOverSigmaStandardDeviation = -1.0;
        m_fIOverSigmaStandardDeviationWRTRotation = -1.0;
        m_fIntensityDistributionMean = -1.0;
        m_fIntensityDistributionSigma = -1.0;
        m_fIntensityDistributionMin = -1.0;
        m_fIntensityDistributionMax = -1.0;
        m_nIntensityDistributionReflnCount = -1;

        m_fIOverSigmaDistributionMean = -1.0;
        m_fIOverSigmaDistributionSigma = -1.0;

        m_oIntensityStats.vClearRaw();
        m_oIOverSigmaStats.vClearRaw();
    }
    
    CResoRankingBin(const CResoRankingBin& anotherBin);
    
    ~CResoRankingBin(){}
    
    CResoRankingBin& operator=(const CResoRankingBin& anotherBin);

    void vSetMinNumberOfReflnsForStats(int nMinNumberOfReflnsForStats){m_nMinNumberOfReflnsForStats=nMinNumberOfReflnsForStats;}
    
    void vAddIOverSigmaValue(double dRot, double dIOverSigma)
                    {m_vecRotAndIOverSigma.push_back(std::pair<double, double>(dRot, dIOverSigma));}

    void vClearIntensityStats(){m_oIntensityStats.vClearRaw();}
    void vClearIOverSigmaStats(){m_oIOverSigmaStats.vClearRaw();}
    
    void vAddToIntensityStats(double fI){m_oIntensityStats.vAddRaw(fI);}
    void vAddToIOverSigmaStats(double fI){m_oIOverSigmaStats.vAddRaw(fI);}
private:
    std::vector< std::pair<double, double> >        m_vecRotAndIOverSigma; // array of pairs: "Rotation, IOverSigma"
    Cstat                                           m_oIntensityStats;
    Cstat                                           m_oIOverSigmaStats;
    
    int     m_nMinNumberOfReflnsForStats;
    double  m_fIOverSigmaMean;
    double  m_fIOverSigmaStandardDeviation;
    double  m_fIOverSigmaStandardDeviationWRTRotation;
    double  m_fIntensityDistributionMean;
    double  m_fIntensityDistributionSigma;
    double  m_fIntensityDistributionMin;
    double  m_fIntensityDistributionMax;
    int     m_nIntensityDistributionReflnCount;

    double  m_fIOverSigmaDistributionMean;
    double  m_fIOverSigmaDistributionSigma;

public:
    void  vCalcIOverSigma(std::vector<CSegment>& vecRotBins);
    void  vDoIntensityStatistics();

    double fGetIOverSigmaMean()const{return m_fIOverSigmaMean;}
    double fGetIOverSigmaMeanStandardDeviation()const{return m_fIOverSigmaStandardDeviation;}
    double fGetIOverSigmaMeanStandardDeviationWRTRotation()const{return m_fIOverSigmaStandardDeviationWRTRotation;}
    
    double fGetIntensityDistributionMean(){return m_fIntensityDistributionMean; }
    double fGetIntensityDistributionSigma(){return m_fIntensityDistributionSigma; }
    double fGetIOverSigmaDistributionMean(){return m_fIOverSigmaDistributionMean; }
    double fGetIOverSigmaDistributionSigma(){return m_fIOverSigmaDistributionSigma; }
    int fGetIntensityDistributionReflectionCount(){return m_nIntensityDistributionReflnCount; }

    int nGetReflnCount(){return (int)m_vecRotAndIOverSigma.size();}
private:
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CResoRASBin : public CResoBin
{
public:
    CResoRASBin()
    {
        m_oDeltaAStats.vClearRaw();
        m_oDeltaCStats.vClearRaw();

        m_fDeltaAAverage = 0.0;
        m_fDeltaCAverage = 0.0;
    }
    
    CResoRASBin(double dMinReso, double dMaxReso) : CResoBin(dMinReso, dMaxReso)
    {
        m_oDeltaAStats.vClearRaw();
        m_oDeltaCStats.vClearRaw();

        m_fDeltaAAverage = 0.0;
        m_fDeltaCAverage = 0.0;
    }
    
    CResoRASBin(const CResoRASBin& anotherBin)
    {
        m_oDeltaAStats = anotherBin.m_oDeltaAStats;
        m_oDeltaCStats = anotherBin.m_oDeltaCStats;

        m_fDeltaAAverage = anotherBin.m_fDeltaAAverage;
        m_fDeltaCAverage = anotherBin.m_fDeltaCAverage;
    }
    
    ~CResoRASBin(){}
    
    CResoRASBin& operator=(const CResoRASBin& anotherBin)
    {
        m_oDeltaAStats = anotherBin.m_oDeltaAStats;
        m_oDeltaCStats = anotherBin.m_oDeltaCStats;
        
        m_fDeltaAAverage = anotherBin.m_fDeltaAAverage;
        m_fDeltaCAverage = anotherBin.m_fDeltaCAverage;

        return *this;
    }

    void vAddDeltaA(double fdA){m_oDeltaAStats.vAddRaw(fdA);}
    void vAddDeltaC(double fdC){m_oDeltaCStats.vAddRaw(fdC);}

    void vTallyUpStats()
    {
        m_fDeltaAAverage = m_oDeltaAStats.fAverageRaw();
        m_fDeltaCAverage = m_oDeltaCStats.fAverageRaw();
    }

    double  fGetDeltaAAverage(){return m_fDeltaAAverage;}
    double  fGetDeltaCAverage(){return m_fDeltaCAverage;}

    int nGetDeltaACount(){return m_oDeltaAStats.nSize();}
    int nGetDeltaCCount(){return m_oDeltaCStats.nSize();}

private:
    Cstat   m_oDeltaAStats;
    Cstat   m_oDeltaCStats;

    double  m_fDeltaAAverage;
    double  m_fDeltaCAverage;
public:
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T_RESOBIN> class DTREK_EXPORT CResoBins
{
public:
    CResoBins(double dMinReso, double dMaxReso, int nNumberOfBins)
    {
        m_dResoMin = -1.0;
        m_dResoMax = -1.0;

        m_dDStarMinCubed = -1.0;
        m_dDStarMaxCubed = -1.0;           
        m_dDStarCubedConstSlope = -1.0;
    
        if( nNumberOfBins <= 0 )
            nNumberOfBins = 10; // default

        vSet(dMinReso, dMaxReso, nNumberOfBins);
    }
    //////////////////////////////////////////////////////////////
    CResoBins()
    {
        m_dResoMin = -1.0;
        m_dResoMax = -1.0;

        m_dDStarMinCubed = -1.0;
        m_dDStarMaxCubed = -1.0;           
        m_dDStarCubedConstSlope = -1.0;
    }
    ///////////////////////////////////////////////////////////////
    virtual ~CResoBins()
    {
        vDeleteBins();
    }
    ///////////////////////////////////////////////////////////////
protected:
    std::vector<T_RESOBIN*>       m_apoBins;
 
    double      m_dResoMin;
    double      m_dResoMax;

private:
    double      m_dDStarMinCubed;  
    double      m_dDStarMaxCubed;
    double      m_dDStarCubedConstSlope;
    
private:
    void vDeleteBins()
    {
        for(int ii=0; ii < (int)m_apoBins.size(); ii++)
        {
            delete m_apoBins[ii];
            m_apoBins[ii] = NULL;
        }
    
        m_apoBins.clear();
    }
public:
    bool bIsAvailable()
    {
        if( 0 == m_apoBins.size() )
            return false;

        if( 0.0 >= m_dDStarCubedConstSlope )
            return false;

        if( m_dResoMin <= m_dResoMax )
            return false;

        if( m_dDStarMinCubed >= m_dDStarMaxCubed )
            return false;

        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void vInitBins(DTREK_WORD wCtrl=DTREK_RESOBIN_ALL)
    {
        for(int ii=0; ii < m_apoBins.size(); ii++)
        {
            m_apoBins[ii]->vInit(wCtrl);
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void vSet(double dMinReso, double dMaxReso, int nNumberOfBins)
    {
        if( 0 >= nNumberOfBins )
        {
            printf("\n\nError: invalid number of resolution bins.\n\n");
            return;
        }
    
        if( dMinReso == dMaxReso )
        {
            printf("\n\nERROR: cannot set up resolution bins when minimum resolution is equal to maximum.\n\n");
            return;
        }

        vDeleteBins();
    
        if( dMinReso < dMaxReso )
	  std::swap(dMinReso, dMaxReso);

        m_dResoMin = dMinReso;
        m_dResoMax = dMaxReso;

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Calculate a resolution shell volume (reciprocal space) without the coefficient of 4/3 PI, 
        // because that coefficient will be present on both sides of the equation,
        // when it comes to calculating the shell limits
        m_dDStarMinCubed = 1.0 / pow((double)dMinReso, 3.);
        m_dDStarMaxCubed = 1.0 / pow((double)dMaxReso, 3.);
    
        m_dDStarCubedConstSlope = (m_dDStarMaxCubed - m_dDStarMinCubed) / nNumberOfBins;


        double      dVolumeReciprocalShell = (m_dDStarMaxCubed - m_dDStarMinCubed) / nNumberOfBins;
    
        double      dBeginReciprocalVolume = m_dDStarMinCubed;
        double      dEndReciprocalVolume = -1.0;  // just to initialize

        T_RESOBIN*       pBin = NULL;
        for(int ii=0; ii < nNumberOfBins; ii++)
        {
            dEndReciprocalVolume = dBeginReciprocalVolume + dVolumeReciprocalShell;

            pBin = new T_RESOBIN(1.0 / pow(dBeginReciprocalVolume, 1.0/3.0), 1.0 / pow(dEndReciprocalVolume, 1.0/3.0));
        
            m_apoBins.push_back(pBin);
        
            dBeginReciprocalVolume = dEndReciprocalVolume;
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // This function gets a resolution bin by "d* cubed". It uses the fact 
    // that in the d* cubed space all bins have a constant width.
    T_RESOBIN*   pGetBin(double dDStarCubed)
    {
    #ifdef RB_DEBUG
        if( m_apoBins.size() <= 0 || m_dDStarCubedConstSlope == 0.0 )
        {
            printf("Problem in CResoBins::pGetBin().\n");
            return NULL;
        }
    #endif
    
        int     iIndex = -1;

        if( dDStarCubed < m_dDStarMaxCubed && dDStarCubed >= m_dDStarMinCubed )
            iIndex    = (int) ((dDStarCubed - m_dDStarMinCubed) / m_dDStarCubedConstSlope );

        if( -1 == iIndex )
            return NULL; // bin not found
        else 
            return m_apoBins[iIndex];
    }
    /////////////////////////////////////////////////////////////////////
    T_RESOBIN*   pGetBin(int iIndex)
    {
        if( iIndex < 0 || iIndex > m_apoBins.size() - 1 )
            return NULL;

        return m_apoBins[iIndex];
    }
    /////////////////////////////////////////////////////////////////////
    T_RESOBIN*   pGetBinReso(double fReso)
    {
        for(int ii=0; ii < m_apoBins.size(); ii++)
        {
	  //cout << "For bin " << ii << " low: " << m_apoBins[ii]->dGetLowReso()
	  //   << " high:" << m_apoBins[ii]->dGetHighReso() << endl;
            if( fReso <= m_apoBins[ii]->dGetLowReso() && fReso > m_apoBins[ii]->dGetHighReso() )
                return m_apoBins[ii];
        }
        
        return NULL;
    }

    /////////////////////////////////////////////////////////////////////

    int nGetNumberOfBins(){return m_apoBins.size();}

    double  dGetMaxReso()const{return m_dResoMax;}
    double  dGetMinReso()const{return m_dResoMin;}

    double  dGetDStarCubedMin()const{return m_dDStarMinCubed;}
    double  dGetDStarCubedMax()const{return m_dDStarMaxCubed;}

    //bool  bAdd(CResoBins* pBins){return true;}
};  
////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CResoStrategyBins : public CResoBins<CResoStrategyBin>
{
public:
    CResoStrategyBins(double dMinReso, double dMaxReso, int nNumberOfBins);
    ~CResoStrategyBins(){}

private:
    double      m_dCompleteness;
    double      m_dMultiplicity;
    double      m_dMultiplicityDev;

    void vInitTotalStats();
public:
    void  vScaleReflnCount(double dScale, DTREK_WORD wCtrl);
    void vCountTheoreticalUniqueReflectionsInResoShells(Ccrystal* poCrystal, bool bAnom);
    void vListStats();
    void vCalculateTotalStats();
    
    double  dGetCompleteness(){return m_dCompleteness;}
    double  dGetMultiplicity(){return m_dMultiplicity;}
    double  dGetMultiplicityDev(){return m_dMultiplicityDev;}
};
///////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CResoRASBins : public CResoBins<CResoRASBin>
{
public:
    CResoRASBins()
    {
        vInit();
    }

public:
    void vInit()
    {
        m_fDeltaAAverage = 0.0;
        m_fDeltaCAverage = 0.0;

        m_nDeltaACount = 0;
        m_nDeltaCCount = 0;
    }
    
    void vTallyUpStats();
    void vPrintStats();

private:
    double  m_fDeltaAAverage;
    double  m_fDeltaCAverage;

    int     m_nDeltaACount;
    int     m_nDeltaCCount;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif   // !DT_RESO_BINS_H

