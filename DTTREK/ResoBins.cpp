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
// ResoBins.cpp       Initial author: RB               04-Apr-2005

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

#include "ResoBins.h"
#include "Ccrystal.h"
#include "Cstat.h"

///////////////////////////////////////////////////////////////////////
CResoBin::CResoBin()
{
    m_dLow  = -1.0;
    m_dHigh = -1.0;
}
///////////////////////////////////////////////////////////////////////
CResoBin::CResoBin(double dLow, double dHigh)
{
    m_dLow  = dLow; 
    m_dHigh = dHigh;
}
///////////////////////////////////////////////////////////////////////
CResoBin::CResoBin(const CResoBin& anotherBin)
{
    m_dLow   = anotherBin.m_dLow;  
    m_dHigh  = anotherBin.m_dHigh; 
}
///////////////////////////////////////////////////////////////////////
CResoBin& CResoBin::operator=(const CResoBin& anotherBin)
{
    m_dLow   = anotherBin.m_dLow;  
    m_dHigh  = anotherBin.m_dHigh; 

    return *this;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function gets the middle resolution of a resolution bin. It uses the fact 
// that in the d* cubed space all bins have a constant width.
double CResoBin::fGetMiddleReso()
{
    if( 0.0 == m_dLow || 0.0 == m_dHigh ) // a sanity check
        return -1.0; 

    double  fDStarMinCubed = 1.0 / pow((double)m_dLow, 3.);
    double  fDStarMaxCubed = 1.0 / pow((double)m_dHigh, 3.);

    double  fMiddleReso = pow((fDStarMaxCubed + fDStarMinCubed)/2, -1.0/3.0);
    
    return fMiddleReso;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tagReflectionStats::tagReflectionStats()
{
    vInit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tagReflectionStats::tagReflectionStats(const tagReflectionStats& anotherStats)
{
    m_nNumRefls       = anotherStats.m_nNumRefls;
    m_nNumRejs        = anotherStats.m_nNumRejs;
    m_nNumUsed        = anotherStats.m_nNumUsed;
    m_nNumMults       = anotherStats.m_nNumMults;
    m_nNumSingles     = anotherStats.m_nNumSingles;
    
    m_nNumUniqueCalc  = anotherStats.m_nNumUniqueCalc;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize a tagReflectionStats structure to all zeros
void tagReflectionStats::vInit(DTREK_WORD wCtrl)
{
    if( wCtrl & DTREK_RESOBIN_SCALE_EXPERIMENTAL_COUNT )
    {
        m_nNumRefls       = 0;
        m_nNumRejs        = 0;
        m_nNumUsed        = 0;
        m_nNumMults       = 0;
        m_nNumSingles     = 0;
    }

    if( wCtrl & DTREK_RESOBIN_SCALE_THEORETICAL_COUNT )
    {
        m_nNumUniqueCalc  = 0;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tagReflectionStats& tagReflectionStats::operator +=(const tagReflectionStats& anotherStats)
{
    m_nNumRefls       += anotherStats.m_nNumRefls;
    m_nNumRejs        += anotherStats.m_nNumRejs;
    m_nNumUsed        += anotherStats.m_nNumUsed;
    m_nNumMults       += anotherStats.m_nNumMults;
    m_nNumSingles     += anotherStats.m_nNumSingles;
    
    m_nNumUniqueCalc  += anotherStats.m_nNumUniqueCalc;
    
    return (*this);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tagReflectionStats::vScaleReflnCount(double dScale, DTREK_WORD wCtrl)
{
    if( wCtrl == DTREK_RESOBIN_SCALE_EXPERIMENTAL_COUNT )
    {
        m_nNumRefls      = (int)((double)m_nNumRefls * dScale);
    
        m_nNumRejs       = (int)((double)m_nNumRejs * dScale);
        m_nNumUsed       = (int)((double)m_nNumUsed * dScale);
        m_nNumMults      = (int)((double)m_nNumMults * dScale);
        m_nNumSingles    = (int)((double)m_nNumSingles * dScale);
    }

    if( wCtrl == DTREK_RESOBIN_SCALE_THEORETICAL_COUNT )
        m_nNumUniqueCalc    = (int)((double)m_nNumUniqueCalc * dScale);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CResoStrategyBin::CResoStrategyBin()
{
    vInit();
}
////////////////////////////////////////////
CResoStrategyBin::CResoStrategyBin(double dMinReso, double dMaxReso) :
CResoBin(dMinReso, dMaxReso)
{
    // We don't call vInit() here, because that would reset m_dLow and m_dHigh
    m_stStats.vInit(DTREK_RESOBIN_ALL);
}
////////////////////////////////////////////
CResoStrategyBin::~CResoStrategyBin()
{
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CResoStrategyBin::CResoStrategyBin(const CResoStrategyBin& anotherBin) :
CResoBin(anotherBin)
{
    m_stStats     = anotherBin.m_stStats;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CResoStrategyBin& CResoStrategyBin::operator=(const CResoStrategyBin& anotherBin)
{
    CResoBin::operator=(anotherBin);

    m_stStats        = anotherBin.m_stStats;

    return *this;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoStrategyBin::vInit(DTREK_WORD wCtrl)
{
    if( wCtrl & DTREK_RESOBIN_LIMITS )
    {
        m_dLow  = -1.0;  
        m_dHigh = -1.0; 
    }

    m_stStats.vInit(wCtrl);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CResoStrategyBins::CResoStrategyBins(double dMinReso, double dMaxReso, int nNumberOfBins) :
CResoBins<CResoStrategyBin>(dMinReso, dMaxReso, nNumberOfBins)
{
    vInitTotalStats();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoStrategyBins::vInitTotalStats()
{
    m_dCompleteness         = 0.0;
    m_dMultiplicity         = 0.0;
    m_dMultiplicityDev      = 0.0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoStrategyBins::vScaleReflnCount(double dScale, DTREK_WORD wCtrl)
{
    for(int ii=0; ii < m_apoBins.size(); ii++)
    {
        m_apoBins[ii]->vScaleReflnCount(dScale, wCtrl);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoStrategyBins::vCountTheoreticalUniqueReflectionsInResoShells(Ccrystal* poCrystal, bool bAnom)
{
    // All this stuff is done to be able to call CCrystal::nCountUnique()
    int         nBins = nGetNumberOfBins();

    double*     pdLowReso  = new double[nBins];
    double*     pdHighReso = new double[nBins];
    int*        pnNumUnique = new int[nBins];

    CResoStrategyBin*       pBin = NULL;

    int     ii = 0;
    /////////////////////////////////////////////
    for(ii=0; ii < m_apoBins.size(); ii++)
    {
        pBin = pGetBin(ii);
    
        pdLowReso[ii] = pBin->dGetLowReso();
        pdHighReso[ii] = pBin->dGetHighReso();
    
        pnNumUnique[ii] = 0;
    }
    /////////////////////////////////////////////

    // Now, make the CALL!
    poCrystal->nCountUnique(nBins,
                            pdLowReso,
                            pdHighReso,
                            pnNumUnique,
                            bAnom);

    // Save results in bins
    for(ii=0; ii < m_apoBins.size(); ii++)
    {
        pBin = pGetBin(ii);
    
        pBin->vSetUniqueCalc(pnNumUnique[ii]);
    }

    delete  []    pdLowReso;
    delete  []    pdHighReso;
    delete  []    pnNumUnique;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoStrategyBins::vListStats()
{
    //This function will caculate and report bin table statistics
    char*         pcLine = "--------------------------------------------------------------------------------\n";

    printf("\n\nExpected Completeness vs Resolution\n");

    printf(pcLine);

    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
           "Resolution  ", "  Calc", "Num", " Num", " Num", "   Num", "Num", " Avg", "%Comp", "%Comp");

    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n", 
           "   range    ", "unique", "reflns", "rejs", "mults", "single", "unique", "mult", "shell", "cumul");

    printf(pcLine);

    CResoStrategyBin*       pBin = NULL;
    CResoStrategyBin        oCumulBin;

    for (int ii = 0; ii < m_apoBins.size(); ii++)
    {
        pBin = m_apoBins[ii];

        if( !pBin )
          continue;

        oCumulBin.vAddStats(*(pBin->pGetStats()));

        if( 0 != pBin->nGetNumUniqueCalc()   && 
            0 != pBin->nGetNumUsed()         &&
            0 != pBin->nGetNumUniqueFound() )
        {
          printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d %6.2f %6.1f ",
                 min(pBin->dGetLowReso(), 999.99), 
                 pBin->dGetHighReso(),
  
                 pBin->nGetNumUniqueCalc(),
                 pBin->nGetNumRefls(),
                 pBin->nGetNumRejs(),
                 pBin->nGetNumMults(),
                 pBin->nGetNumSingles(),
                 pBin->nGetNumUniqueFound(),
                 (double) (pBin->nGetNumUsed()) / (double) pBin->nGetNumUniqueFound(),
                 min(100., 100. * (double) pBin->nGetNumUniqueFound() / (double) pBin->nGetNumUniqueCalc()));
        }
        else if( 0 < pBin->nGetNumUniqueCalc() )
        {
          printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d %6s %6.1f ",
                  min(pBin->dGetLowReso(), 999.99), 
                  pBin->dGetHighReso(),
                  pBin->nGetNumUniqueCalc(),
                  pBin->nGetNumRefls(), 
                  pBin->nGetNumRejs(),
                  pBin->nGetNumMults(), 
                  pBin->nGetNumSingles(),
                  pBin->nGetNumUniqueFound(),
                  "---",
                  min(100., 100. * (double) pBin->nGetNumUniqueFound() / (double) pBin->nGetNumUniqueCalc()));
        }
        else
        {
          printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d %6s %6s ",
                  min(pBin->dGetLowReso(), 999.99), 
                  pBin->dGetHighReso(),
                  pBin->nGetNumUniqueCalc(),
                  pBin->nGetNumRefls(), 
                  pBin->nGetNumRejs(),
                  pBin->nGetNumMults(),
                  pBin->nGetNumSingles(),
                  pBin->nGetNumUniqueFound(),
                  "---", 
                  "---");
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if( 0 != oCumulBin.nGetNumUniqueCalc() )
          printf("%6.1f\n", 
            min(100., 100. * (double) oCumulBin.nGetNumUniqueFound() / (double)oCumulBin.nGetNumUniqueCalc()));
        else
          printf("%6s\n", "---");
    }

    printf(pcLine);

    printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d ",
           min(m_dResoMin, 999.99), m_dResoMax,
           oCumulBin.nGetNumUniqueCalc(),
	       oCumulBin.nGetNumRefls(), 
           oCumulBin.nGetNumRejs(),
	       oCumulBin.nGetNumMults(), 
           oCumulBin.nGetNumSingles(),
           oCumulBin.nGetNumUniqueFound());

    if( 0.0 != oCumulBin.nGetNumUniqueCalc() )
    {
        double  dMultiplicity = (double) oCumulBin.nGetNumUsed() / (double) max(1, oCumulBin.nGetNumUniqueFound());
        double  dCompleteness = min(100., 100.* (double) oCumulBin.nGetNumUniqueFound() / (double) oCumulBin.nGetNumUniqueCalc());
    
        printf("%6.2f ",  dMultiplicity);
        printf("%6.1f ",  dCompleteness);
        printf("%6.1f\n", dCompleteness);
    }
    else
    {
        printf("%6s %6s %6s\n", "---", "---", "---");
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoStrategyBins::vCalculateTotalStats()
{
    vInitTotalStats();

    Cstat       oStat;
    oStat.vClearRaw();

    CResoStrategyBin*     pBin = NULL;
    CResoStrategyBin      oCumulBin;

    for (int ii = 0; ii < m_apoBins.size(); ii++)
    {
        pBin = m_apoBins[ii];

        oCumulBin.vAddStats(*(pBin->pGetStats()));

        oStat.vAddRaw( (double)pBin->nGetNumUsed() / max(1.0, (double)pBin->nGetNumUniqueFound()) );
    }     

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( oCumulBin.nGetNumUniqueCalc() > 0 )
        m_dCompleteness = min(100.0, 100.0 * ((double)oCumulBin.nGetNumUniqueFound() / (double)oCumulBin.nGetNumUniqueCalc()));
    else
        m_dCompleteness = 0.0;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( oCumulBin.nGetNumUniqueFound() > 0 )
        m_dMultiplicity =   (double) (oCumulBin.nGetNumUsed()) / (double) oCumulBin.nGetNumUniqueFound();
    else
        m_dMultiplicity = 0.0;

    m_dMultiplicityDev = oStat.fStandardDeviationRaw();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CResoRankingBin::CResoRankingBin(const CResoRankingBin& anotherBin)
: CResoBin(anotherBin)
{
    m_vecRotAndIOverSigma.clear();

    for(int ii=0; ii < anotherBin.m_vecRotAndIOverSigma.size(); ii++)
    {
        m_vecRotAndIOverSigma.push_back(anotherBin.m_vecRotAndIOverSigma[ii]);
    }
    
    m_nMinNumberOfReflnsForStats    = anotherBin.m_nMinNumberOfReflnsForStats;
    
    m_fIOverSigmaMean               = anotherBin.m_fIOverSigmaMean;
    m_fIOverSigmaStandardDeviation  = anotherBin.m_fIOverSigmaStandardDeviation;
    m_fIOverSigmaStandardDeviationWRTRotation  = anotherBin.m_fIOverSigmaStandardDeviationWRTRotation;
    
    m_fIntensityDistributionMean = anotherBin.m_fIntensityDistributionMean;
    m_fIntensityDistributionSigma = anotherBin.m_fIntensityDistributionSigma;
    
    m_fIOverSigmaDistributionMean = anotherBin.m_fIOverSigmaDistributionMean;
    m_fIOverSigmaDistributionSigma = anotherBin.m_fIOverSigmaDistributionSigma;
    
    m_fIntensityDistributionMin = anotherBin.m_fIntensityDistributionMin;
    m_fIntensityDistributionMax = anotherBin.m_fIntensityDistributionMax;
    m_nIntensityDistributionReflnCount = anotherBin.m_nIntensityDistributionReflnCount;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CResoRankingBin& CResoRankingBin::operator=(const CResoRankingBin& anotherBin)
{   
    CResoBin::operator=(anotherBin);

    m_vecRotAndIOverSigma.clear();
    
    for(int ii=0; ii < anotherBin.m_vecRotAndIOverSigma.size(); ii++)
    {
        m_vecRotAndIOverSigma.push_back(anotherBin.m_vecRotAndIOverSigma[ii]);
    }

    m_nMinNumberOfReflnsForStats = anotherBin.m_nMinNumberOfReflnsForStats;
    
    m_fIOverSigmaMean               = anotherBin.m_fIOverSigmaMean;
    m_fIOverSigmaStandardDeviation  = anotherBin.m_fIOverSigmaStandardDeviation;
    m_fIOverSigmaStandardDeviationWRTRotation  = anotherBin.m_fIOverSigmaStandardDeviationWRTRotation;

    m_fIntensityDistributionMean = anotherBin.m_fIntensityDistributionMean;
    m_fIntensityDistributionSigma = anotherBin.m_fIntensityDistributionSigma;

    m_fIOverSigmaDistributionMean = anotherBin.m_fIOverSigmaDistributionMean;
    m_fIOverSigmaDistributionSigma = anotherBin.m_fIOverSigmaDistributionSigma;

    m_fIntensityDistributionMin = anotherBin.m_fIntensityDistributionMin;
    m_fIntensityDistributionMax = anotherBin.m_fIntensityDistributionMax;
    m_fIntensityDistributionMax = anotherBin.m_fIntensityDistributionMax;
    m_nIntensityDistributionReflnCount = anotherBin.m_nIntensityDistributionReflnCount;

    return *this;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The purposes of this function are to evaluate IOverSigma, its overall deviation and deviation over rotation range. 
// For the latter purpose, all reflections, contained in this class's m_vecRotAndIOverSigma array will be binned, 
// according to their rotation value, into rotation range bins, specified by the input array vecRotBins.  
// If a bin gets fewer than a m_nMinNumberOfReflnsForStats reflections, the IOverSigma of that bin is considered to be 0.
void CResoRankingBin::vCalcIOverSigma(std::vector<CSegment>& vecRotBins)
{
    // Create a vector such that each element will store the "number of reflns" (int) and their "total IOverSigma" (double) for a rotation range bin.
    std::vector< std::pair<int, double> >    vecRotBinContent; 
    
    int     ii = 0;
    int     jj = 0;

    /////////////////////////////////////////////
    // Allocate bins and initialize bin contents
    for(ii=0; ii < vecRotBins.size(); ii++)
    {
        vecRotBinContent.push_back( std::pair<int, double>(0, 0.0) ); // initialize
    }
    
    ////////////////////////////////////////////
    // Fill rotation bin contents and calculate the OVERALL IOverSigma mean value at the same time
    double      dRot = 0.0;
    m_fIOverSigmaMean = 0.0;
    for(jj=0; jj < m_vecRotAndIOverSigma.size(); jj++)
    {
        dRot = m_vecRotAndIOverSigma[jj].first;
    
        for(ii=0; ii < vecRotBins.size(); ii++)
        {
            if( vecRotBins[ii].bIsInLeft(dRot) )  // if a bin goes from t_beg to t_end, an element with a value equal to t_end belongs to the next bin on the right
            {
                // Bin found!
                (vecRotBinContent[ii].first)++; // up the refln count
                
                (vecRotBinContent[ii].second) += m_vecRotAndIOverSigma[jj].second; // update the IOverSigma sum
                break;
            }
        }
        
        m_fIOverSigmaMean += m_vecRotAndIOverSigma[jj].second;
    }

    m_fIOverSigmaMean /= max(1.0, (double)m_vecRotAndIOverSigma.size());  // safe divide

    /////////////////////////////////////////////////////////
    if( vecRotBinContent.size() == 1 )
    {
        m_fIOverSigmaStandardDeviationWRTRotation = 0.0; // Must be just one image
        return;
    }
    /////////////////////////////////////////////////////////

    ///////////////////////////////////////////
    // Build an object for rotation statistics
    Cstat       oRotStat;
    oRotStat.vClearRaw();

    int     nNumberOfValidRotationBins = 0;

    for(ii=0; ii < vecRotBinContent.size(); ii++)
    {
        if( vecRotBinContent[ii].first > m_nMinNumberOfReflnsForStats )
        {
            oRotStat.vAddRaw(vecRotBinContent[ii].second / vecRotBinContent[ii].first); // a mean IOverSigma value for this rotation bin
            nNumberOfValidRotationBins++;
        }
        else
            oRotStat.vAddRaw(0.0);
    }                                   
    
    if( nNumberOfValidRotationBins == 0 )
    {
        m_fIOverSigmaStandardDeviationWRTRotation = -1.0; // Deviation cannot be estimated
        return;
    }

    ////////////////////////////////////////////
    //Do statistics
    m_fIOverSigmaStandardDeviationWRTRotation = oRotStat.fStandardDeviationRaw();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoRankingBin::vDoIntensityStatistics()
{
    m_fIntensityDistributionMean = m_oIntensityStats.fAverageRaw();
    m_fIntensityDistributionSigma = m_oIntensityStats.fStandardDeviationRaw();
    
    m_nIntensityDistributionReflnCount = m_oIntensityStats.nSize();

    m_fIOverSigmaDistributionMean = m_oIOverSigmaStats.fAverageRaw();
    m_fIOverSigmaDistributionSigma = m_oIOverSigmaStats.fStandardDeviationRaw();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoRASBins::vTallyUpStats()
{
    vInit();

    for(int ii=0; ii < m_apoBins.size(); ii++)
    {
        m_apoBins[ii]->vTallyUpStats();
        
        m_nDeltaACount += m_apoBins[ii]->nGetDeltaACount();  
        m_nDeltaCCount += m_apoBins[ii]->nGetDeltaCCount();

        m_fDeltaAAverage += m_apoBins[ii]->fGetDeltaAAverage() * m_apoBins[ii]->nGetDeltaACount(); 
        m_fDeltaCAverage += m_apoBins[ii]->fGetDeltaCAverage() * m_apoBins[ii]->nGetDeltaCCount(); 
    }
    
    m_fDeltaAAverage = m_nDeltaACount > 0 ? m_fDeltaAAverage/m_nDeltaACount : 0.0;
    m_fDeltaCAverage = m_nDeltaCCount > 0 ? m_fDeltaCAverage/m_nDeltaCCount : 0.0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CResoRASBins::vPrintStats()
{
    printf("\n\nAnomalous scattering signal analysis*\n");
    printf("--------------------------------------------------------\n");
    printf("   Resolution    Num |I+ - I-|    Delta    Delta   R(as)\n");
    printf("        range  acentric  centric  acentric centric shell\n");
    printf("--------------------------------------------------------\n");
    
    double      fRAS = 0.0;
    
    char    cTemp[256];
    
    Cstring     sDeltaCAverage("");
    Cstring     sDeltaAAverage("");
    Cstring     sRAS("");

    for(int ii=0; ii < m_apoBins.size(); ii++)
    {
        ///////////////////////////////////////////////////////////////////
        if( m_apoBins[ii]->nGetDeltaACount() > 0 )
            sprintf(cTemp, "%9.2f", m_apoBins[ii]->fGetDeltaAAverage());
        else
            sprintf(cTemp, "        -");
        sDeltaAAverage = cTemp;
        ////////////////////////////////////////////////////////////////////
        if( m_apoBins[ii]->nGetDeltaCCount() > 0 )
            sprintf(cTemp, "%9.2f", m_apoBins[ii]->fGetDeltaCAverage());
        else
            sprintf(cTemp, "        -");
        sDeltaCAverage = cTemp;
        
        ///////////////////////////////////////////////////////////////////
        if( m_apoBins[ii]->nGetDeltaCCount() > 0 &&
            m_apoBins[ii]->nGetDeltaACount() > 0 && m_apoBins[ii]->fGetDeltaCAverage() != 0.0 )
        {
            fRAS = m_apoBins[ii]->fGetDeltaAAverage() / m_apoBins[ii]->fGetDeltaCAverage();
            sprintf(cTemp, "%7.2f", fRAS);
        }
        else
             sprintf(cTemp, "      -");
        
        sRAS = cTemp;
        /////////////////////////////////////////////////////////////////////

        
        printf("%6.2lf -%5.2lf%9d%9d%9s%9s%7s\n", 
                   m_apoBins[ii]->dGetLowReso(), 
                   m_apoBins[ii]->dGetHighReso(),
                   m_apoBins[ii]->nGetDeltaACount(),
                   m_apoBins[ii]->nGetDeltaCCount(),
                   sDeltaAAverage.string(),
                   sDeltaCAverage.string(),
                   sRAS.string());
    }

    printf("--------------------------------------------------------\n");

    ///////////////////////////////////////////////////////////////////
    if( m_nDeltaACount > 0 )
        sprintf(cTemp, "%9.2f", m_fDeltaAAverage);
    else
        sprintf(cTemp, "        -");
    sDeltaAAverage = cTemp;
    
    ////////////////////////////////////////////////////////////////////
    if( m_nDeltaCCount > 0 )
        sprintf(cTemp, "%9.2f", m_fDeltaCAverage);
    else
        sprintf(cTemp, "        -");
    sDeltaCAverage = cTemp;

    ///////////////////////////////////////////////////////////////////
    if( m_nDeltaCCount > 0 &&
        m_nDeltaACount > 0 && m_fDeltaCAverage != 0.0 )
    {
        fRAS = m_fDeltaAAverage / m_fDeltaCAverage;
        sprintf(cTemp, "%7.2f", fRAS);
    }
    else
         sprintf(cTemp, "      -");
    
    sRAS = cTemp;
    /////////////////////////////////////////////////////////////////////
    printf("%6.2lf -%5.2lf%9d%9d%9s%9s%7s\n", 
                                    dGetMinReso(), 
                                    dGetMaxReso(),
                                    m_nDeltaACount,
                                    m_nDeltaCCount,
                                    sDeltaAAverage.string(),
                                    sDeltaCAverage.string(),
                                    sRAS.string());

    printf("\n*See Fu, Z.-Q.,Rose, J.P., Wang, B.-C. (2004) Acta Cryst. D60, 499-506;"
           "\n   equations 7, 8 and 9. Delta = <|I+ - I-|/sigma(I)>"
           "\n   for acentric (A) and centric (C) reflns.  R(as) = DeltaA / DeltaC\n");
    fflush(stdout);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
