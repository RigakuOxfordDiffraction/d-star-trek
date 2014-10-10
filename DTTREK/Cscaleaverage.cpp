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
// dtscaleaverage.cc     Initial author: T.J. Niemeyer  Spring 2000
//             Based on dtscalemerge 
//    This reads a d*TREK style reflection file, then scales the reflections
//    and averages them.
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

#include "Cscaleaverage.h"
#include "Cstat.h"
#include "Cimage.h"
#include "Csource.h"

#ifdef SSI_PC
#include "CrclHelper.h"
#endif

#include "dtsvd.h"

// This is here only for the sake of dup() and dup2() in one routine.
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

#define USE_SIGMA_WEIGHTS TRUE

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////        tagStatsScaleAverage   ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


void tagStatsScaleAverage::vInitStats() {
  int nx;
  // Initialize a tagStats structure to all zeroes
  nNumRefs    = 0;
  nNumRejects = 0;
  nNumUsed    = 0;
  nNumUsed2   = 0;
  nNumExcluded= 0;
  nNumMults   = 0;
  nNumSingles = 0;

  fNumer        = 0.0;
  fNumer_rim    = 0.0;
  fNumer_pim    = 0.0;
  fNumer_pcv    = 0.0;

  fNumer_rimplus  = 0.0;
  fNumer_rimminus = 0.0;
  fNumer_plus     = 0.0;
  fNumer_minus    = 0.0;

  fDenom        = 0.0;
  fDenomAnom    = 0.0;
  fSumChiSq     = 0.0;
  fSumChiSqA    = 0.0;
  fEAddNumer    = 0.0;
  fEAddDenom    = 0.0;
  nSumChiGroups = 0;
  nSumChiContrib= 0;
  nSumChiAGroups= 0;
  nSumChiAContrib= 0;
  fSumIoverSig  = 0.0;
  fSumIoverSig2 = 0.0;
  fSum          = 0.0;

  nLastRefFlag = -1;
  for (nx = 0; nx<10; nx++)
      anNumInBinMults[nx] = 0;
  return;
};

tagStatsScaleAverage& tagStatsScaleAverage::operator += (tagStatsScaleAverage& oTag) {
  int nx;
  nNumRefs      += oTag.nNumRefs;
  nNumRejects   += oTag.nNumRejects;
  nNumUsed      += oTag.nNumUsed;
  nNumUsed2     += oTag.nNumUsed2;
  nNumExcluded  += oTag.nNumExcluded;
  nNumMults     += oTag.nNumMults;
  nNumSingles   += oTag.nNumSingles;

  fNumer    += oTag.fNumer;

  //shijie yao -- 04:01:08 added members for additonal merging R factors
  fNumer_rim += oTag.fNumer_rim;
  fNumer_pim += oTag.fNumer_pim;
  fNumer_pcv += oTag.fNumer_pcv;

  fNumer_rimplus  += oTag.fNumer_rimplus;
  fNumer_rimminus += oTag.fNumer_rimminus;
  fNumer_plus     += oTag.fNumer_plus;
  fNumer_minus    += oTag.fNumer_minus;

  fDenom        += oTag.fDenom;
  fDenomAnom    += oTag.fDenomAnom;
  fSumChiSq     += oTag.fSumChiSq;
  fSumChiSqA    += oTag.fSumChiSqA;
  fEAddNumer    += oTag.fEAddNumer;
  fEAddDenom    += oTag.fEAddDenom;
  nSumChiGroups += oTag.nSumChiGroups;
  nSumChiContrib+= oTag.nSumChiContrib;
  nSumChiAGroups += oTag.nSumChiAGroups;
  nSumChiAContrib+= oTag.nSumChiAContrib;
  fSumIoverSig  += oTag.fSumIoverSig;
  fSumIoverSig2 += oTag.fSumIoverSig2;
  fSum          += oTag.fSum;

  for (nx=0;nx<10;nx++)
      anNumInBinMults[nx] += oTag.anNumInBinMults[nx];
  return *this;
};

Cscaleaverage::Cscaleaverage(Creflnlist& oListIn,Cimage_header& oHeader) : 
m_oReflnlist(oListIn), 
m_oCrystal(oHeader)
{
    int nx;
    float f0;

    // Get the crystal B matrix into variable
    m_oCrystal.nCalcBMatrix();
    m_oCrystal.vGetBMatrix(&m_a3x3fBMat[0][0]);

    m_nNumRefs = m_oReflnlist.nGetNumReflns();

    m_bExcludeSysAbsences = FALSE;
    m_bComputeOverlap=FALSE;
    m_bUseBackgroundSigma = FALSE;

    // JWP: These next 2 lines are kinda weird 
    //      since one should just look at m_nFI_fIntensity and m_nFI_fSigmaI

    m_nFI_nIntensityI =  oListIn.nExpandGetField(oListIn.ms_sfIntensity);
    m_nFI_nSigmaI     =  oListIn.nExpandGetField(oListIn.ms_sfSigmaI);

    if ((m_nFI_nIntensityI<0) || (m_nFI_nSigmaI<0) || (m_oReflnlist.m_nFI_nPackedHKL<0)) {
        printf("Did not find all reflection list fields:\n");
        cout << oListIn.ms_ssBatch << "\n"
             << oListIn.ms_sfIntensity << "\n"
             << oListIn.ms_sfSigmaI << "\n"
             << oListIn.ms_snPackedHKL << "\n" << flush;
    };

    if (m_oReflnlist.m_nFI_sBatch < 0) {
      Cstring sTemp;
      m_nBatchI = m_oReflnlist.nExpandGetField(m_oReflnlist.ms_ssBatch);
      printf("Reflection list does not have batch field...\n");
      printf("Adding batch field with batch number of '1001'\n");
      sTemp = "1001";
      for (nx=0;nx<m_oReflnlist.nGetNumReflns();nx++) {
          m_oReflnlist[nx].vSetField(m_oReflnlist.m_nFI_sBatch,sTemp);
      };

    } else
        m_nBatchI = m_oReflnlist.m_nFI_sBatch;

    m_fUserEAdd = -1.0;
    m_fUserEMul = -1.0;
    m_fUserEAddMult = 1.0;

    m_fMaxScaleIntensity = 1e20;
    m_fSigma=3.0;
    m_fReject=5.0;
    m_fDevReject=500000.0;
    m_fMaxDevReject = 50.0;
    m_fChiSquared=0.0;
    m_nGroups = 0;
    m_bRejectBatches = FALSE;
    m_sBatchRejectTemplate = "";
    m_pfMergeSigma = new float[m_nNumRefs];
    m_pfChi = new float[m_nNumRefs];
    m_pfDev = new float[m_nNumRefs];
    m_pfDev2 = new float[m_nNumRefs];
    m_pfOrigIntensity = new float[m_nNumRefs];
    m_pnGroup = new int[m_nNumRefs];
    m_pnGroupSize = new int[m_nNumRefs];


    m_nChiBins = 9;
    m_pfChiDistribBins = new float[m_nChiBins];
    m_pnChiDistribGroups = new int[m_nChiBins];
    m_pfChiDistribMultiplicity = new float[m_nChiBins];

    for (nx=0;nx<m_nChiBins;nx++)
      {
        f0 = (float)nx;
        
        f0 = exp((f0/(m_nChiBins-1))*log(10.0));
        
        m_pfChiDistribBins[nx] = (float)( ((int) (f0*10))/10.0 );
      };

    // Hardcode rather than use previous for loop 

    m_pfChiDistribBins[0] = 1.3f;
    m_pfChiDistribBins[1] = 1.6f;
    m_pfChiDistribBins[2] = 2.0f;
    m_pfChiDistribBins[3] = 3.0f;
    m_pfChiDistribBins[4] = 4.0f;
    m_pfChiDistribBins[5] = 5.0f;
    m_pfChiDistribBins[6] = 6.5f;
    m_pfChiDistribBins[7] = 8.0f;
    m_pfChiDistribBins[8] = 10.0f;

    // Table information initialized.

    m_pnReject        = new int[m_nNumRefs];


    for (nx=0;nx<5;nx++) {
        m_pfMin[nx] = NULL;
        m_pfMax[nx] = NULL;
    };

    m_pnIndexBatch    = NULL;

    m_pnIndexIntensity  = NULL;
    m_pfRevIndexReso = NULL;
    m_pnIndexReso = NULL;
    m_pnNumOverlap    = NULL;
    m_pnRevIndex = NULL;

    for (nx=0;nx<m_nNumRefs;nx++)
        m_pnReject[nx]=0;

    m_ptStatsReso     = NULL;
    m_ptStatsInt      = NULL;
    m_ptStatsBatch    = NULL;
    
    m_ptBatchInfo     = NULL;

    m_pfErrorModelEAdd = NULL;
    m_pfErrorModelEMul = NULL;
    m_fErrorModelEAddMinIntensity = 10.0;
    m_fTargetChiSquare = 1.0;

    m_nNumBatches         = 0;
    m_nVerbose            = 2;
    m_nNumBinsReso        = 10;
    m_nNumBinsInt         = 10;
    m_a2fRangeReso[0]     = 0.0;
    m_a2fRangeReso[1]     = 0.0;
    m_a2fRangeInt[0]      = 0.0;
    m_a2fRangeInt[1]      = 20.0;

    m_poScan = new Cscan(oHeader);
    m_pfAbsorbFactors = NULL;
    m_pfAbsorbFactorsSig = NULL;

    m_poPrevRejectList = NULL;
    m_poPartialList = NULL;

    Csource *poSource = NULL;
    poSource  = new Csource(oHeader);
    m_fWavelength = poSource->fGetWavelength();
    delete poSource;
    poSource = NULL;
};

Cscaleaverage::~Cscaleaverage() {
    int nx;

    delete[] m_pfMergeSigma;
    delete[] m_pfChi;
    delete[] m_pfDev;
    delete[] m_pfDev2;
    delete[] m_pfOrigIntensity;
    delete[] m_pnGroup;
    delete[] m_pnGroupSize;

    delete[] m_ptBatchInfo;
    delete[] m_pnNumOverlap;
    delete[] m_pnReject;


    delete[] m_ptStatsReso;
    delete[] m_ptStatsInt;
    delete[] m_ptStatsBatch;

    delete[] m_pfChiDistribBins;
    delete[] m_pnChiDistribGroups;
    delete[] m_pfChiDistribMultiplicity;
        
    delete[] m_pnIndexBatch;
    delete[] m_pnIndexIntensity;
    delete[] m_pfRevIndexReso;
    delete[] m_pnIndexReso;
    delete[] m_pnRevIndex;

    delete m_poScan;
    
    if( m_pfAbsorbFactors )
        delete[] m_pfAbsorbFactors;
    
    if( m_pfAbsorbFactorsSig )
        delete[] m_pfAbsorbFactorsSig;
    
    if( m_pfErrorModelEAdd )
        delete[] m_pfErrorModelEAdd;
    
    if( m_pfErrorModelEMul )
        delete[] m_pfErrorModelEMul;

    for (nx=0;nx<5;nx++)
    {
        delete[] m_pfMin[nx];
        delete[] m_pfMax[nx];
    };

    if (m_poPrevRejectList)
        delete m_poPrevRejectList;
    if (m_poPartialList)
        delete m_poPartialList;
};




double Cscaleaverage::fCalcChiSquared(int nRecomputeRejects)
{
    int nRef,nRefSort;
    int nx,ny,nz;
    int nIdx;               // Batch index.
    int nLastHKL;           // packed HKL of last in reflection in list.
    int nThisHKL;           // packed HKL for current reflection in list.
    int nStartIndex;
    int nEndIndex;
    int nGroup;             // Total number of reflections in the given group.
    int nGroupChi;          // Total number of non chi-rejected reflections in the given group.
    int nGroupChiA;         // Same as previous, but I+ and I- treated separately
    int nGroupsChiA;
    int nGroups;            // Total number of chi-squared groups.
    int nGroupsWithSingles; // Total number of chi-squared groups (plus singles)
    int nNumUsed;           // Number of reflns (participating in some sym-equiv group with multiplicity >=2) Used in Chi-Squared.
    int nBatchRejectCount;  // Number of reflections rejected due to batch rejection.  
    int nBatchRejectPass;   // Number of times we reject batches.

    double fTotalChi;
    double fGroupAvg;        // Average for a symmetrically equiv. group of reflections.
    double fGroupInvVar;     // Used to compute the average.
    double fGroupSigmaMerge; // Computed from the deviations from the fGroupAvg.
    double fGroupSigma;      // Computed from the (unadjusted) sigmas in the group.
    double fGroupSigmaAdj;   // Computed from the (adjusted) sigmas in the group.
    double fGroupChi;        // The chi-squared contribution for the group.
    double fGroupChiA;       // Same as previous, but I+ and I- treated separately
    double fGroupEAdd;       // Eadd used in a group.
    double fIoverSigAvg;     // The average I/sig for unaveraged reflections.
    double fChi;             // The chi-squared for an individual reflection.
    double fChiA;            // The chi-squared for an individual reflection with I+, I- treated separately

    double fGroupAvgPlus;
    double fGroupAvgMinus;
    double fGroupInvVarPlus;
    double fGroupInvVarMinus;
    int    nGroupPlus;
    int    nGroupMinus;
    int    nCentPhase;
    int    nAnomFlag;
//+2010-05-25 JWP
    double fSumGroupInt, fSumGroupSigma, fNumGroupsUsed;
//-2010-05-25 JWP

    double f0,f1;
    double f0InvSquared;
    double f0PlusDelta;
    double f0MinusDelta;

    double fIntensity,fSigma,fRawSigma,fReso;

    double fTemp, fTemp_pcv;

//+2010-07-20
    bool   bUseSigmaWeights = USE_SIGMA_WEIGHTS;
    double fReflnWeight = 1.0;
//-2010-07-20

    int         nBinInt;
    int         nBinReso;
    tagStatsScaleAverage*   ptReso;
    tagStatsScaleAverage*   ptInt;
    tagStatsScaleAverage*   ptBatch;


    fGroupSigma    = 0.0;
    fGroupSigmaAdj = 0.0;

    static Cstat oStat;


    if (nRecomputeRejects) {
        Cstring sTemp;
        Cstring sTemp2;
        char*   pc;

        // Clear all rejection flags.
        for (nx=0;nx<m_nNumRefs;nx++) 
            m_pnReject[nx]=0;
        
        if ((m_fFractionReject) && (m_fChiSquared!=0.0))
            m_fDevReject = fFindDevRejection(m_fFractionReject);
        else
            m_fDevReject= 500000.0;

        sTemp = m_sBatchRejectTemplate;
        pc = strtok(sTemp.string(),",");

        for (nx=0; nx<m_nNumBatches; nx++) 
            m_ptBatchInfo[nx].m_bBatchRejectUser = FALSE;

        while (pc) {
            sTemp2 = pc;
            if (sTemp2.contains('-')) 
                nRejectBatches(sTemp2.before('-'),sTemp2.after('-'));
            else              
                nRejectBatches(sTemp2,(Cstring) "");
            pc = strtok(NULL,",");
        };
        // Copy the user defined rejections over.
        for (nx=0; nx<m_nNumBatches; nx++)
            m_ptBatchInfo[nx].m_bBatchReject = m_ptBatchInfo[nx].m_bBatchRejectUser;        
    };

    nBatchRejectPass = 0;
    //m_oReflnlist.nWrite("jimtest.ref");

    do 
    {    
        // Reject all batches that are flagged with reject flags.
        for (nx=0; nx<m_nNumRefs;nx++)
        {
            int nIdx = m_oReflnlist[m_pnIndex[nx]].nGetField(m_nFI_nBatchIdx);
            if (m_ptBatchInfo[nIdx].m_bBatchReject)
                m_pnReject[nx] = 1;
        }

        //This outer "do" loop is controled by the batch rejection flags.
        //We must loop once through all of the reflections to determine 
        //rejects and statistics.
        //If batch rejection is enabled, we then reject batches based 
        //on the computed statistics.
        //If no batches are rejected we are done.  If one or more batches are 
        //rejected, we need to 
        //make a second pass.  We continue until no new batches are rejected.
        

        
        // Initialize all of the arrays if we are computing statistics.
        // We will assume that these arrays exist.
        
        int* pnTemp;
        // Zero all sums
        if (m_bComputeOverlap)
        {
            pnTemp = m_pnNumOverlap;
            for (nx = 0; nx < m_nNumBatches+1; nx++)
            {
                for (ny = 0; ny < m_nNumBatches; ny++)
                {
                    *pnTemp++ = 0;
                }
            }
        };
        
        for (nx = 0; nx < m_nNumBinsInt+1; nx++)
        {
            m_ptStatsInt[nx].vInitStats();
        }
        
        for (nx = 0; nx < m_nNumBinsReso+1; nx++)
        {
            m_ptStatsReso[nx].vInitStats();
        }
        
        for (nx = 0; nx < m_nNumBatches+1; nx++)
        {
            m_ptStatsBatch[nx].vInitStats();
        }
        
        for (nx=0;nx<10;nx++)
        {
            m_a10nMult[nx] = 0;
            m_a10nMultLE[nx] =  0;
            m_a10nUsedLE[nx] = 0;
            m_a10fSumChi[nx] = 0.0;
        }
        
        
        for (nx=0;nx<m_nChiBins;nx++) {
            m_pnChiDistribGroups[nx] = 0;
            m_pfChiDistribMultiplicity[nx] = 0.0;
        };
        
        
        nGroups = 0;
        nGroupsWithSingles = 0;
        nNumUsed = 0;
        nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        nStartIndex = 0;
        fTotalChi =0.0;
        m_nTotalAllRejects = 0;
        m_nTotalSigmaExcluded = 0;
        m_nTotalAdjustedChiRejects = 0;
        memset(m_pfMergeSigma,0,m_nNumRefs*sizeof(float));
        memset(m_pfChi,0,m_nNumRefs*sizeof(float));
        memset(m_pfDev,0,m_nNumRefs*sizeof(float));
        memset(m_pfDev2,0,m_nNumRefs*sizeof(float));
        
        int nMaxGroupSize = 0;
        
//+2010-05-25 JWP
	// Compute I/Sigma(I)

	fSumGroupInt   = 0.0;
	fSumGroupSigma = 0.0;
	fNumGroupsUsed = 0.0;
//-2010-05-25 JWP

        for (nRefSort=0;nRefSort<m_nNumRefs+1;nRefSort++)
        {
            nRef=m_pnIndex[nRefSort];
            
            if (nRefSort<m_nNumRefs)
                nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
//+-JWP
	    /*
   if ( 1 == (nThisHKL % 2))
     cout << "odd: " << nThisHKL << endl;
   else
     cout << "even: " << nThisHKL << endl;
	    */
//+-JWP
            if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
                continue;

            nEndIndex = nRefSort - 1;
            

            // For each refln in this group,
            // Save the group index for later use 
            // and the number of reflns in this group
            // nGroupsWithSingles here is number of groups 
            //    including groups with singles

            for (nx = nStartIndex; nx <= nEndIndex; nx++)
            {
                m_pnGroup[nx] = nGroupsWithSingles;
                m_pnGroupSize[nx] = (nEndIndex - nStartIndex + 1);
            }
            
            nGroupsWithSingles++; // Increment the group number

            // We have start index and end index, so begin computations.           

            for (nGroup = 0,nx=nStartIndex; nx<=nEndIndex; nx++)
            {
              // for each reflection in the group that is not rejected ...

              if( m_pnReject[nx] == 0 )
              {
                fIntensity = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nIntensityI);
                fRawSigma  = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nSigmaI);
                nIdx       = m_oReflnlist[m_pnIndex[nx]].nGetField(m_nFI_nBatchIdx);

                if (fRawSigma <= 0.0)
                {
                  m_pnReject[nx]=1;           // These must get excluded from the statistics, since they represent bad data.
                } 
                else if ((fIntensity<m_fSigma*fRawSigma) || (fIntensity >= m_fMaxScaleIntensity))
                  {
                    m_pnReject[nx]=2;           // NOT excluded from statistics.
                  }
              }
                
	      if (m_pnReject[nx]!=1)
		nGroup++;                      // Count number of non-excluded reflns in this group
            }
            
            if (nGroup > nMaxGroupSize) 
              nMaxGroupSize = nGroup;            // Save the size of the largest group
            
            // We want to reject reflections that are outside
            // of the rejection limits for symmetry-equivalent reflections.

            f0 = 0.0;
            fGroupAvg = 0.0;
            for (nx=nStartIndex; nx<=nEndIndex; nx++)
            {
                fIntensity = m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                
                fSigma  = fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]);
                
                fRawSigma = max(1.0,m_oReflnlist[m_pnIndex[nx]].fGetSigmaI());
                

                if( m_pnReject[nx] != 1 )
                {
//+2010-07-20
		  if (bUseSigmaWeights) fReflnWeight = max(1.0, (fSigma*fSigma));
		  fGroupAvg += fIntensity / fReflnWeight;
//-2010-07-20
                    f0 += 1.0/max(1.0,(fSigma*fSigma));
                };
            };
            
            if (f0)
              fGroupAvg /= f0;            


            // Reject Chi-squared.
            if (nGroup>1)
            {
                fGroupChi = 0.0;
            
                for (nx=nStartIndex,ny=0; nx<=nEndIndex; nx++)
                {
                    fSigma     = fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]);
                    fIntensity = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nIntensityI);
                    fRawSigma  = m_oReflnlist[m_pnIndex[nx]].fGetSigmaI();                   
                    
                    
                    if( 0.0 != fSigma )
                        fChi = (fIntensity-fGroupAvg)*(fIntensity-fGroupAvg)/max(1.0,(fSigma*fSigma));
                    else
                        fChi = 0.0;      
                    
                    m_pfDev[nx] = (float)fChi;
                    m_pfDev2[nx]= (float)( fIntensity-fGroupAvg );
                    
                    if (m_pnReject[nx]!=1)
                    {
                        fGroupChi += fChi;
                    }
                }
                
                fGroupChi /= (nGroup-1);
                
                for (nx=nStartIndex;nx<=nEndIndex;nx++)
                {
                    m_pfChi[nx] = (float)fGroupChi;
                    
                    if (m_pnReject[nx]!=1)
                    {
                        if (m_pfDev[nx] > m_fDevReject)
                        {
                            m_pnReject[nx] = 1;
                            nGroup--;
                        };
                    };
                };
            } else {
                m_pfDev2[nStartIndex] = 0.0;
                m_pfDev[nStartIndex] = 0.0;
                m_pfChi[nStartIndex] = 0.0;
            };

            
            // Compute average using only non-rejected reflections.
	    fGroupAvg = 0.0;
	    fGroupInvVar = 0.0;
	    fIoverSigAvg = 0.0;

	    fGroupAvgPlus     = 0.0;
	    fGroupAvgMinus    = 0.0;
	    fGroupInvVarPlus  = 0.0;
	    fGroupInvVarMinus = 0.0;
	    nGroupPlus        = 0;
	    nGroupMinus       = 0;
            if (nGroup>0) {
	      for (nx=nStartIndex; nx<=nEndIndex; nx++) {
		if (m_pnReject[nx]!=1) {
		  f1 = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nIntensityI);
		  f0 = max(1.0,fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]));
		  f0InvSquared = 1.0 / (f0 * f0);

		  fGroupAvg += f1 * f0InvSquared;
		  fIoverSigAvg += f1 / f0;
		  fGroupInvVar += f0InvSquared;
			
		  // JWP 2008-06-12
		  // TODO: use temp variables for this
		  nAnomFlag  = m_oReflnlist[m_pnIndex[nx]].nGetField(m_oReflnlist.m_nFI_nFplusminusFlag);
		  nCentPhase = m_oReflnlist[m_pnIndex[nx]].nGetField(m_oReflnlist.m_nFI_nCentPhase);
		  if (nCentPhase > 0) nAnomFlag = 1;  // Force centrics to be counted as I+'s
		  if (1 == nAnomFlag) {
		    // This refln belongs to I+ or centric
		    fGroupAvgPlus    += f1 * f0InvSquared;
		    fGroupInvVarPlus += f0InvSquared;
		    nGroupPlus++;
		    //cout << "+++\n";
		  }
		  else {
		    // This refln belongs to I-
		    fGroupAvgMinus    += f1 * f0InvSquared;
		    fGroupInvVarMinus += f0InvSquared;
		    nGroupMinus++;
		    //cout << "---\n";
		  }
		};
	      };
	      fIoverSigAvg /= nGroup;    
	      
	      // Compute weighted average of the intensity

	      if (fGroupInvVar)
		fGroupAvg/= fGroupInvVar;
	      if (fGroupInvVarPlus)
		fGroupAvgPlus /= fGroupInvVarPlus;
	      if (fGroupInvVarMinus)
		fGroupAvgMinus /= fGroupInvVarMinus;
	      
	      /*** else- block not needed if we init to zero above if-block
		 } else {
                fGroupAvg = 0.0;
                fGroupInvVar = 0.0;
                fIoverSigAvg = 0.0;
	      ***/
            };
           
            if( nGroup > 0 )
            {
                // Compute the group sigma.           
                if( nGroup > 1 )
                {
                    f0 = 0.0;
                    for (nx=nStartIndex; nx<=nEndIndex; nx++)
                    {
                        if (m_pnReject[nx]!=1)
                        {
                            f0 += (m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nIntensityI) -fGroupAvg)*
                                (m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nIntensityI) -fGroupAvg);
                        }
                    }
                    
                    fGroupSigmaMerge = sqrt((double)(f0/(nGroup-1)));
                }
                else
                    fGroupSigmaMerge =0.0;

                for (nx=nStartIndex; nx<=nEndIndex; nx++)
                {
                    if (m_pnReject[nx]==1)
                        m_pfMergeSigma[nx] = 0.0f;
                    else
                        m_pfMergeSigma[nx] = (float)fGroupSigmaMerge;
                }

                // Compute the unadjusted sigma and adjusted sigma for the group.
                fGroupSigma = 0.0;
                fGroupSigmaAdj = 0.0;
                for (nx=nStartIndex; nx<=nEndIndex; nx++) {
                    if (m_pnReject[nx]!=1) {
                        f0 = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nSigmaI);
                        fGroupSigma+= f0*f0;
                        f0 = fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]);
                        fGroupSigmaAdj += f0*f0;
                    };
                };
                // SIG
                fGroupSigmaAdj=sqrt(fGroupSigmaAdj)/nGroup;
                fGroupSigma=sqrt(fGroupSigma)/nGroup;                              
                
            } else {
                nx = nStartIndex;
                f0 = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nSigmaI);
                fGroupSigma= f0;
                f0 = fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]);
                fGroupSigmaAdj = f0;
       
            };

//+2010-05-25 JWP
	    if( nGroup > 1 )
	      {
		fSumGroupInt   += fGroupAvg;
		fSumGroupSigma += fGroupSigma;
		fNumGroupsUsed += 1.0;
	      }
//-2010-05-25 JWP	      

            // Compute the group chi squared.
            // This does its own rejection.
            fGroupChi = 0.0;
            nGroupChi = 0;
	    fGroupChiA = 0.0;
	    int nGroupChiPlusMinus = 0;
	    // There are more unique reflns if I+ and I- are kept separate
	    //
	    nGroupChiA = 0;
	    nGroupsChiA = 0;
            fGroupEAdd = 0.0;
	    bool bFoundPlus  = FALSE;
	    bool bFoundMinus = FALSE;
	    //+JWP 2008-06-17
	    //cout << "nGroup, nGroup+-: " << nGroup << ", " << nGroupPlus << ", " << nGroupMinus << endl;
	    //-JWP 2008-06-17

            for (nx = nStartIndex; nx <= nEndIndex; nx++) {
                if (m_pnReject[nx]!=1) {

		  // This reflection is NOT rejected

		  fIntensity = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nIntensityI);
		  fSigma     = fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]);
		  fRawSigma  = m_oReflnlist[m_pnIndex[nx]].fGetSigmaI();
		  fChi       = 0.0;
		  fChiA      = 0.0;
		  if (0.0 != fSigma) {
		    nAnomFlag  = m_oReflnlist[m_pnIndex[nx]].nGetField(m_oReflnlist.m_nFI_nFplusminusFlag);
		    nCentPhase = m_oReflnlist[m_pnIndex[nx]].nGetField(m_oReflnlist.m_nFI_nCentPhase);
		    if (nCentPhase > 0) nAnomFlag = 1;  // Force centrics to be counted as I+'s

		    // We need to compute differences with different group averages
		    // so we get a reduced ChiSq with I+,I- treated together 
		    // and also treated separately

		    fChi = (fIntensity-fGroupAvg)*(fIntensity-fGroupAvg)/(fSigma*fSigma);
		    //if ( (1 < nGroupMinus) && (-1 == nAnomFlag) )
		    if (-1 == nAnomFlag)
		      {
			// This refln is an I- AND it has more than one refln in the group of I-'s
			fChiA = (fIntensity-fGroupAvgMinus)*(fIntensity-fGroupAvgMinus)/(fSigma*fSigma);
			bFoundMinus = TRUE;
		      }
		    //else if ( (1 < nGroupPlus) && (+1 == nAnomFlag) )
		    else if (+1 == nAnomFlag)
		      {
			// This refln is an I+ AND it has more than one refln in the group of I+'s
			fChiA = (fIntensity-fGroupAvgPlus)*(fIntensity-fGroupAvgPlus)/(fSigma*fSigma);
			bFoundPlus = TRUE;
		      }
		  }
		  if (fChi<m_fDevReject) {
		    fGroupChi += fChi;
		    fGroupEAdd += m_fLastEAddUsed;
		    nGroupChi ++;
		    fGroupChiA += fChiA;
		    nGroupChiA++;
		  };
		  if (fChi>=m_fMaxDevReject) {
		    m_nTotalAdjustedChiRejects++;
		  };
                };
            };
	    if (1 < nGroupPlus) nGroupsChiA++;
	    if (1 < nGroupMinus) nGroupsChiA++;
	    /***
	    cout << "nGroupChi, nGroupChiA, nGroupPlus, nMin: " 
                 << nGroupChi << ", " << nGroupChiA
                 << ", " << nGroupPlus
                 << ", " << nGroupMinus << endl;
	    ***/
            fGroupEAdd/=max(1,nGroupChi);
            
            // Fill in the bins: Batch, Resolution and Intensity.

            fTemp_pcv = 0; //shijie yao : for Rpcv calculation
	    //cout << "Gstart\n";
            for (nx = nStartIndex, ny = 0; nx <= nEndIndex; nx++)
            {
	      /****
	      cout << "HKL: " << m_oReflnlist[m_pnIndex[nx]].nGetH()
                   << ", " << m_oReflnlist[m_pnIndex[nx]].nGetK()
                   << ", " << m_oReflnlist[m_pnIndex[nx]].nGetL()
                   << "  Anom: " 
                   << m_oReflnlist[m_pnIndex[nx]].nGetField(m_oReflnlist.m_nFI_nFplusminusFlag) << endl;
	      ****/
                fIntensity = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_nIntensityI);
                fSigma     = fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]);
                
                
                nIdx     = m_oReflnlist[m_pnIndex[nx]].nGetField(m_nFI_nBatchIdx);
                fReso    = m_oReflnlist[m_pnIndex[nx]].fGetField(m_nFI_f2STLsq);
                fReso    = fReso * sqrt(fReso);
                if ((0.0 != fGroupSigmaAdj) && (m_fSlopeInt != 0.0))
                    nBinInt  = (int) ((fGroupAvg/fGroupSigmaAdj) / m_fSlopeInt);
                else
                    nBinInt = 0;
                nBinInt  = min(nBinInt,m_nNumBinsInt-1);
                nBinInt  = max(nBinInt,0);
                nBinReso = (int) ((fReso - m_a2fRangeReso[0]) / max(1e-10,m_fSlopeReso));
                nBinReso = min(nBinReso,m_nNumBinsReso-1);
                nBinReso = max(nBinReso,0);
                ptReso   = &m_ptStatsReso[nBinReso];
                ptInt    = &m_ptStatsInt[nBinInt];
                ptBatch  = &m_ptStatsBatch[nIdx];
                
                ptReso->nNumRefs++;
                ptInt->nNumRefs++;
                ptBatch->nNumRefs++;
                
                if (m_pnReject[nx] == 1)
                {
                    ptReso->nNumRejects++;
                    ptInt->nNumRejects++;
                    ptBatch->nNumRejects++;       
                    m_nTotalAllRejects++;
                    continue;
                };
                
                if (m_pnReject[nx] == 2) 
                {
                    ptReso->nNumExcluded++;
                    ptInt->nNumExcluded++;
                    ptBatch->nNumExcluded++;
                    m_nTotalSigmaExcluded++;
                };
                
                // Compute I/Sigma(I)
                
                f0                     = fIoverSigAvg;
                ptReso->fSumIoverSig  += f0;
                ptInt->fSumIoverSig   += f0;
                ptBatch->fSumIoverSig += f0;
                if (0.0 != fGroupSigmaAdj)
                    f0                     = fGroupAvg / fGroupSigmaAdj;
                else
                    f0 = 0.0;
                ptReso->fSumIoverSig2 += f0;
                ptInt->fSumIoverSig2  += f0;
                ptBatch->fSumIoverSig2+= f0;
                
                // Want to have fSum for computing the average.
                
                ptReso->fSum          += fIntensity;
                ptInt->fSum           += fIntensity;
                ptBatch->fSum         += fIntensity;
                
                // Add in the multiplicites.
                
                {
                    int nMultBin;
                    if (5 > nGroup)
                    {
                        nMultBin = nGroup;
                    }
                    else if (12 < nGroup)
                    {
                        nMultBin = 7;
                    }
                    else if (8 < nGroup)
                    {
                        nMultBin = 6;
                    }
                    else if (4 < nGroup)
                    {
                        nMultBin = 5;
                    }
                    if (ptReso->nLastRefFlag!=nStartIndex)
                        ptReso->anNumInBinMults[nMultBin]++;
                    if (ptInt->nLastRefFlag!=nStartIndex)
                        ptInt->anNumInBinMults[nMultBin]++;
                    if (ptBatch->nLastRefFlag!=nStartIndex)
                        ptBatch->anNumInBinMults[nMultBin]++;
                };
                
                
                if (nGroup<=1) {
                    ptReso->nNumSingles++;
                    ptInt->nNumSingles++;
                    ptBatch->nNumSingles++;
                    continue;
                } else {
                    ptReso->nNumUsed++;
                    ptInt->nNumUsed++;
                    ptBatch->nNumUsed++;
                };
                
                // Compute number of multiples.
                
		if (nGroup>1) { 
                    // Look back and see when we last referenced the structure
                    if (ptReso->nLastRefFlag!=nStartIndex) {
                        ptReso->nNumMults++;
                        ptReso->nNumUsed2+=nGroup;
                        if (nGroupChi>=2) {
                            ptReso->fEAddNumer += fGroupEAdd*fGroupAvg;
                            ptReso->fEAddDenom += fGroupAvg;
                            ptReso->fSumChiSq+=fGroupChi;
                            ptReso->nSumChiGroups++;
                            ptReso->nSumChiContrib+=nGroupChi;
			    if (nGroupChiA >= 2) {
			      ptReso->fSumChiSqA+=fGroupChiA;
			      ptReso->nSumChiAGroups++;
			      ptReso->nSumChiAGroups += nGroupsChiA;
			      ptReso->nSumChiAContrib+=nGroupChiA;
			    }
                        };
			/*****
                        if (nGroupChiA >= 2) {
			  ptReso->fSumChiSqA+=fGroupChiA;
			  ptReso->nSumChiAGroups++;
			  ptReso->nSumChiAContrib+=nGroupChiA;
			}
			********/
                    };
                    ptReso->nLastRefFlag=nStartIndex;                   
                    
                    if (ptInt->nLastRefFlag!=nStartIndex) {
                        ptInt->nNumMults++;
                        ptInt->nNumUsed2+=nGroup;
                        if (nGroupChi>=2) {
                            ptInt->fEAddNumer += fGroupEAdd*fGroupAvg;
                            ptInt->fEAddDenom += fGroupAvg;
                            ptInt->fSumChiSq+=fGroupChi;
                            ptInt->nSumChiGroups++;
                            ptInt->nSumChiContrib+=nGroupChi;
			    if (nGroupChiA >= 2) {
			      ptInt->fSumChiSqA+=fGroupChiA;
			      ptInt->nSumChiAGroups++;
			      ptInt->nSumChiAContrib+=nGroupChiA;
			    };
                        };
			/**
                        if (nGroupChiA >= 2) {
                            ptInt->fSumChiSqA+=fGroupChiA;
                            ptInt->nSumChiAGroups++;
                            ptInt->nSumChiAContrib+=nGroupChiA;
                        };
			**/
                        
                    };
                    ptInt->nLastRefFlag=nStartIndex;
                    
                    if (ptBatch->nLastRefFlag!=nStartIndex) {
                        ptBatch->nNumMults++;
                        ptBatch->nNumUsed2+=nGroup;
                        if (nGroupChi>=2) {
                            ptBatch->fEAddNumer += fGroupEAdd*fGroupAvg;
                            ptBatch->fEAddDenom += fGroupAvg;
                            ptBatch->fSumChiSq+=fGroupChi;
                            ptBatch->nSumChiGroups++;
                            ptBatch->nSumChiContrib+=nGroupChi;
                        };
			/**
                        if (nGroupChiA>=2) {
                            ptBatch->fSumChiSqA+=fGroupChiA;
                            ptBatch->nSumChiAGroups++;
                            ptBatch->nSumChiAContrib+=nGroupChiA;
                        };
			**/
                    };
                    ptBatch->nLastRefFlag=nStartIndex;
                };
                
                // Compute R-merge statistics.
                //cout << "fDenom+\n";
                ptReso->fDenom+=fabs(fIntensity);
                ptInt->fDenom+=fabs(fIntensity);
                ptBatch->fDenom+=fabs(fIntensity);
                f0 = fabs(fIntensity - fGroupAvg);

                ptReso->fNumer+=f0;
                ptInt->fNumer+=f0;
                ptBatch->fNumer+=f0;

                //shijie yao -- 04:01:08 added members for additonal merging R factors

                //r.i.m algorithm: only difference to Rmerge of above
                //is to scale the group deviation sum (SUM(I - Iagv) by sqrt(N/N-1).
                //We can do this because each item of (I-Iavg) will have the same
                //scaling factor for the sum. SUM {f(N) SUM(I-Iavg) } = SUM { SUM f(N)*(I-Iagv) }.
		//cout << "Numer_rim\n";
                fTemp = f0 * sqrt( (double)nGroup / (double)(nGroup - 1) );
                ptReso->fNumer_rim  += fTemp;
                ptInt->fNumer_rim   += fTemp;
                ptBatch->fNumer_rim += fTemp;

		//+JWP 2008-06-11
		// Use p.i.m. for the I+ != I- statistics

                //p.i.m algorithm only difference to Rmerge() of above
                //is to scale the group deviation sum (SUM(I - Iagv) by sqrt(1/N-1)
                //See the above comments for r.i.m

		fTemp = 0.0;
		nAnomFlag = m_oReflnlist[m_pnIndex[nx]].nGetField(m_oReflnlist.m_nFI_nFplusminusFlag);
		nCentPhase = m_oReflnlist[m_pnIndex[nx]].nGetField(m_oReflnlist.m_nFI_nCentPhase);
		if (nCentPhase > 0) nAnomFlag = 1;  // Force centrics to be counted as I+'s

		// PROBLEM if -scaleanom is NOT active, 
		// then fDenom will have added to it an fIntensity singlet
		// Solution: create and use fDenomAnom

		//cout << "Avg: " << fGroupAvg << "  nGroup+-: " << nGroupPlus << ", " << nGroupMinus << endl;
                if ( (1 < nGroupMinus) && (-1 == nAnomFlag) )
		  {
		    // This refln is an I- AND it has more than one refln in the group of I-'s

		    // fGroupAvgPlus and fGroupAvgMinus may be 0.0 if there
		    // are no reflns in this class

		    //cout << "Avg-: " << fGroupAvgMinus << endl;
		    f0MinusDelta = fabs(fIntensity - fGroupAvgMinus);

		    fTemp = f0MinusDelta * sqrt( (double)nGroupMinus / (double)(nGroupMinus - 1) );
		    //cout << "Numer_pim-\n";
		    ptReso->fNumer_pim  += fTemp;
		    ptInt->fNumer_pim   += fTemp;
		    ptBatch->fNumer_pim += fTemp;
		    ptReso->fDenomAnom  += fabs(fIntensity);
		    ptInt->fDenomAnom   += fabs(fIntensity);
		    ptBatch->fDenomAnom +=fabs(fIntensity);

		  }
                else if ( (1 < nGroupPlus) && (+1 == nAnomFlag) )
		  {
		    // This refln is an I+ AND it has more than one refln in the group of I+'s

		    // fGroupAvgPlus and fGroupAvgMinus may be 0.0 if there
		    // are no reflns in this class
		    //cout << "Avg+: " << fGroupAvgPlus << endl;
		    f0PlusDelta  = fabs(fIntensity - fGroupAvgPlus);
		    fTemp = f0PlusDelta * sqrt( (double)nGroupPlus / (double)(nGroupPlus - 1) );
		    //cout << "Numer_pim+\n";
		    ptReso->fNumer_pim  += fTemp;
		    ptInt->fNumer_pim   += fTemp;
		    ptBatch->fNumer_pim += fTemp;
		    ptReso->fDenomAnom  += fabs(fIntensity);
		    ptInt->fDenomAnom   += fabs(fIntensity);
		    ptBatch->fDenomAnom += fabs(fIntensity);
		  }
		/****
                fTemp = f0 * sqrt( 1.0 / (double)(nGroup - 1) );
                ptReso->fNumer_pim  += fTemp;
                ptInt->fNumer_pim   += fTemp;
                ptBatch->fNumer_pim += fTemp;
		****/

		//-JWP 2008-06-11
		
                //PCV algorithm: the sqrt() need to apply after 
                //the whole group sum. But here is only a single |I-Iavg|.
                /**
                fTemp = sqrt((f0 * f0) / (double)(nGroup - 1) );
                ptReso->fNumer_pcv  += fTemp;
                ptInt->fNumer_pcv   += fTemp;
                ptBatch->fNumer_pcv += fTemp;   
                */
	    };  // end of nx loop
	    //cout << "Gend\n\n";
        
        // Add contribution to this batch, and overlaps in other batches.
        
        if (m_bComputeOverlap) {
            // Since nReduce() sorts secondarily on batch number, we can search through
            // this group to find the start of each "group" of batch numbers.  We will use
            // the oStat to store this information, along with the length of each group.
            // Thus, even elments of oStat contain the batch indexes,
            //       odd elements of oStat contain unique batch lengths.
            
            // In this loop, nz marks the begining of a equi-batch block, and ny marks the
            // end of an equi-batch block.
            oStat.vClear();
            for (nz = nStartIndex,ny = nStartIndex+1; ny<=nEndIndex+1;ny++) {
                if ((ny==nEndIndex+1) ||
                    (m_oReflnlist[m_pnIndex[ny]].nGetField(m_nFI_nBatchIdx)!=
                    m_oReflnlist[m_pnIndex[ny-1]].nGetField(m_nFI_nBatchIdx))) {
                    // Insert the batch index.
                    oStat.vAdd(m_oReflnlist[m_pnIndex[ny-1]].nGetField(m_nFI_nBatchIdx));
                    // Count number of NON-REJECTED reflections in the group in nCount.
                    int nRef;
                    int nCount;
                    for (nCount=0,nRef=nz;nRef<ny;nRef++) {
                        if (m_pnReject[nRef]!=1)
                            nCount++;
                    };
                    // Insert the number of reflections.
                    oStat.vAdd(nCount);
                    nz=ny;
                };
            }
            // We will do a double loop, which will get each pair of batch names.
            // The pair's batches are (oStat[ny],oStat[nz]) and the pair's lengths are (oStat[ny+1],oStat[nz+1])
            for (ny=0;ny<oStat.nSize();ny+=2) {
                for (nz=0;nz<oStat.nSize();nz+=2) {
                    if (nz==ny) {
                        // Only count overlaps here if there really were other batches.
                        if (oStat.nSize()/2>1)
                            m_pnNumOverlap[((int) oStat[ny])*m_nNumBatches+((int) oStat[nz])] += (int) oStat[ny+1];
                    } else
                        m_pnNumOverlap[((int) oStat[ny])*m_nNumBatches+((int) oStat[nz])] += (int) oStat[ny+1];
                };
            };
        };
        
        
        // This accumulation is from earlier code, but it is still taken. 
        // (Recently, it has become more useful)
        
        if (nGroupChi>1) {
            m_a10nMult[min(9,nGroupChi)]++;
            for (nx=1;nx<=min(9,nGroupChi);nx++) {
                m_a10nMultLE[nx]++;
                m_a10nUsedLE[nx] += nGroupChi;
                m_a10fSumChi[nx] += (float)fGroupChi;
            };
            
        };
        
        
        
        if ((nGroupChi!=0) && (nGroupChi!=1)) {
            // Compute which group this chi belongs in.
            for (nx=0;nx<m_nChiBins-1;nx++) {
                if (m_pfChiDistribBins[nx]>fGroupChi/(float)(nGroupChi-1))
                    break;
            };
            
            m_pnChiDistribGroups[nx]++;
            m_pfChiDistribMultiplicity[nx]+=nGroupChi;
            
            fTotalChi+= fGroupChi;
            nNumUsed+=nGroupChi;
            nGroups++;
        } else {
            // ???
            
        };
        
        nStartIndex=nRefSort;
        if (nRefSort<m_nNumRefs)
            nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        
    }; //  end of the BIG: for (nRefSort=0;nRefSort<m_nNumRefs+1;nRefSort++)
    

    for (nx=0;nx<m_nChiBins;nx++) {
        if (m_pnChiDistribGroups[nx])
            m_pfChiDistribMultiplicity[nx]/=m_pnChiDistribGroups[nx];
    };
    
    
    nBatchRejectCount = 0;
    if ((nRecomputeRejects==2) && (m_bRejectBatches) && (nBatchRejectPass==0)) {
        
        bool bReject;
        float fRmergeAverage;
        float fRmergeSigma;
        int   nRmergeRejectsSig = 0;
        int   nRmergeRejectsAbs = 0;
        float fIoverSigAverage;
        float fIoverSigSigma;
        int   nIoverSigRejects = 0;
        float fChiSqAverage;
        float fChiSqSigma;
        int   nChiSqRejectsSig = 0;
        int   nChiSqRejectsAbs = 0;
        int   nPercentBadRejects = 0;
        
        
        // Compute the R-merge sigma and average.
        oStat.vClear();
        for (nx = 0; nx<m_nNumBatches; nx++) {
            if (m_ptStatsBatch[nx].fDenom)
                oStat.vAdd((m_ptStatsBatch[nx].fNumer / m_ptStatsBatch[nx].fDenom));
        };
        oStat.nMarkStat(false,Cstat::S_ABOVEBELOW,5.0);
        fRmergeAverage = (float) oStat.fAverage();
        fRmergeSigma = (float) oStat.fStandardDeviation();
        
        printf("INFO: Distribution of batch Rmerge  average: %6.3f, sigma: %6.3f\n", fRmergeAverage, fRmergeSigma);

        // Compute the I/sig sigma and average.
        oStat.vClear();
        for (nx = 0; nx<m_nNumBatches; nx++) {
            if (m_ptStatsBatch[nx].nNumUsed + m_ptStatsBatch[nx].nNumSingles)
                oStat.vAdd((m_ptStatsBatch[nx].fSumIoverSig2 / (m_ptStatsBatch[nx].nNumUsed + m_ptStatsBatch[nx].nNumSingles)));
        };
        oStat.nMarkStat(false,Cstat::S_ABOVEBELOW,5.0);
        fIoverSigAverage = (float)oStat.fAverage();
        fIoverSigSigma = (float)oStat.fStandardDeviation();
        
        // Compute the Chi-Sq sigma and average.
        oStat.vClear();
        for (nx = 0; nx<m_nNumBatches; nx++) {
            f0 = m_ptStatsBatch[nx].fSumChiSq;
            f1 = (m_ptStatsBatch[nx].nSumChiContrib - m_ptStatsBatch[nx].nSumChiGroups);
            if (f1)
                oStat.vAdd(f0/f1);
        };
        oStat.nMarkStat(false,Cstat::S_ABOVEBELOW,5.0);
        fChiSqAverage = (float)oStat.fAverage();
        fChiSqSigma = (float)oStat.fStandardDeviation();
        
        printf("INFO: Distribution of batch |ChiSq| average: %6.3f, sigma: %6.3f\n", fChiSqAverage, fChiSqSigma);

        // Reject batches based on a large number of rejects.
        // Reject batches who's Rmerge is outside deviation.
        // Reject batches who's <I/sig> value is outside deviation.
        
        if (m_bRejectBatches)
          {
            //printf ("m_fBatchChiSqReject: %.3f\n", m_fBatchChiSqReject);

            if (0.0 > m_fBatchChiSqReject)
              {
                printf("INFO: Reject batches with |ChiSq| > %.3f\n",
                       ABS(m_fBatchChiSqReject));
              }
            else if (1000.0 <= m_fBatchChiSqReject)
              {
                // Print nothing
              }
            else if (2.0 <= m_fBatchChiSqReject)
              {
                printf("INFO: Reject batches with |ChiSq| > %.3f\n      (%0.3f sigma above the average batch |ChiSq|)\n",
                       m_fBatchChiSqReject*fChiSqSigma+fChiSqAverage,
                       m_fBatchChiSqReject);
              }
            if ( (0.0 < m_fBatchRmergeReject) && (1.0 > m_fBatchRmergeReject) )
              {
                printf("INFO: Reject batches with Rmerge > %.3f\n",
                       m_fBatchRmergeReject);
              }
            else if (1000.0 <= m_fBatchRmergeReject)
              {
                // Print nothing
              }
            else if (2.0 <= m_fBatchRmergeReject)
              {
                printf("INFO: Reject batches with Rmerge > %.3f\n      (%0.3f sigma above the average batch Rmerge)\n",
                       m_fBatchRmergeReject*fRmergeSigma+fRmergeAverage,
                       m_fBatchRmergeReject);
              }
            
          }
        
        for (nx = 0; nx<m_nNumBatches; nx++) {
            
            if ((m_ptStatsBatch[nx].nNumRefs) && (((float) m_ptStatsBatch[nx].nNumRejects)/m_ptStatsBatch[nx].nNumRefs>m_fBatchPercentReject)) {
                if (!m_ptBatchInfo[nx].m_bBatchReject)
                    nPercentBadRejects++;
                bReject = TRUE;
            } else if ((m_fBatchRmergeReject>2.0) && (ABS(m_ptStatsBatch[nx].fNumer / max(1.0,m_ptStatsBatch[nx].fDenom) - fRmergeAverage)> fRmergeSigma*m_fBatchRmergeReject)) {
                if (!m_ptBatchInfo[nx].m_bBatchReject)
                    nRmergeRejectsSig++;
                bReject = TRUE;
            } else if ((m_fBatchRmergeReject<1.0) && (m_ptStatsBatch[nx].fNumer / max(1.0,m_ptStatsBatch[nx].fDenom)> m_fBatchRmergeReject)) {
                if (!m_ptBatchInfo[nx].m_bBatchReject)
                    nRmergeRejectsAbs++;
                bReject = TRUE;
            } else if (fabs(m_ptStatsBatch[nx].fSumIoverSig2/ max(1,m_ptStatsBatch[nx].nNumUsed + m_ptStatsBatch[nx].nNumSingles)) > fIoverSigSigma*m_fBatchIoverSigReject) {
                if (!m_ptBatchInfo[nx].m_bBatchReject)
                    nIoverSigRejects++;
                bReject = TRUE;
            } else if ((m_fBatchChiSqReject>2.0) && 
                (0.0 != (f1 = (m_ptStatsBatch[nx].nSumChiContrib - m_ptStatsBatch[nx].nSumChiGroups))) && 
                ((fabs(m_ptStatsBatch[nx].fSumChiSq/f1 - fChiSqAverage)) > fChiSqSigma*m_fBatchChiSqReject)) {
                if (!m_ptBatchInfo[nx].m_bBatchReject)
                    nChiSqRejectsSig++;
                bReject = TRUE;                
            } else if ((m_fBatchChiSqReject<0.0) && 
                (0.0 != (f1 = (m_ptStatsBatch[nx].nSumChiContrib - m_ptStatsBatch[nx].nSumChiGroups))) && 
                ((m_ptStatsBatch[nx].fSumChiSq/f1 > ABS(m_fBatchChiSqReject)))) {
                if (!m_ptBatchInfo[nx].m_bBatchReject)
                    nChiSqRejectsAbs++;
                bReject = TRUE;                
            } else
                bReject = FALSE;
            
            if (bReject) {
                // Reject reflections in this batch. (The actual change to m_pnReject is done at the top of the loop.
                m_ptBatchInfo[nx].m_bBatchReject = TRUE;
                // Make sure to only count reflections that were not previously rejected.
                
                for (nRef=0;nRef<m_nNumRefs;nRef++) {
                    int nIdx = m_oReflnlist[m_pnIndex[nRef]].nGetField(m_nFI_nBatchIdx);
                    if (nIdx == nx) {
                        if (m_pnReject[nRef] != 1) {
                            nBatchRejectCount++;                            
                        };
                    };
                };
            };
        };
        if (m_bRejectBatches)
          if (0 == (nRmergeRejectsAbs + nRmergeRejectsSig + nChiSqRejectsSig + nChiSqRejectsAbs))
            printf("\nINFO: No batches rejected due to Rmerge or |ChiSq| cutoffs.\n");

        if (nRmergeRejectsAbs)
          {
            printf("\nINFO:  %d batches rejected with Rmerge > %.2lf\n",nRmergeRejectsAbs,m_fBatchRmergeReject);
          }
        if (nRmergeRejectsSig)
          {
            printf("\nINFO:  %d batches rejected with Rmerge > %.2lf\n",nRmergeRejectsSig,fRmergeAverage + fRmergeSigma*m_fBatchRmergeReject);
          }
        if (nChiSqRejectsAbs)
          {
            printf("\nINFO:  %d batches rejected with |ChiSq| > %.2lf\n",nChiSqRejectsAbs, ABS(m_fBatchChiSqReject));
          }
        if (nChiSqRejectsSig)
          {
        printf("\nINFO:  %d batches rejected with |ChiSq| > %.2lf\n",nChiSqRejectsSig,fChiSqAverage + fChiSqSigma*m_fBatchChiSqReject);
          }
    }

    nBatchRejectPass++;
    } while ((nRecomputeRejects==2) && (nBatchRejectCount));

    if (nNumUsed - nGroups>0)
        m_fChiSquared = (float)(fTotalChi/(nNumUsed-nGroups));
    else
        m_fChiSquared = 0.0f;
    m_nGroups = nGroups;
    m_nNumUsed = nNumUsed;
    
//+2010-05-25 JWP
    if (0.0 < fNumGroupsUsed)
      {
	fSumGroupInt   = fSumGroupInt   / fNumGroupsUsed;
	fSumGroupSigma = fSumGroupSigma / fNumGroupsUsed;
	
	printf("\nINFO <Intensity> = %.1lf <SigmaInt> = %.1lf,\n  so <Intensity>/<SigmaInt> = %.2lf for %.0lf reflns\n",
	       fSumGroupInt,
	       fSumGroupSigma,
	       fSumGroupInt / fSumGroupSigma,
	       fNumGroupsUsed);
		   
      }
//-2010-05-25 JWP	      

    return m_fChiSquared;
}

bool Cscaleaverage::bCalcAnomalousSignal(Creflnlist *poReflnlist)
{
  // Analyze the anomalous signal in the reflection list output file
  // that contains I+, sigI+, I-, sigI-

  // This routine tries compute an anomalous scattering signal in the data based
  // on the work published here:
  //   Fu, Z.-Q., Rose, J.P. and Wang, B.C. (2004) Acta Crystal. D60, 499-506.
  //
  // There is some question as to what sigmas and intensities to use.

  if (   (0 > poReflnlist->m_nFI_fIntensityPlus)
      || (0 > poReflnlist->m_nFI_fIntensityMinus)
      || (0 > poReflnlist->m_nFI_fSigmaIPlus)
      || (0 > poReflnlist->m_nFI_fSigmaIMinus)
     )
    {
      printf("\nINFO:  no anomalous information in the reflnlist.\n");
      return (FALSE);
    }

  m_oRASBins.vSet
    (
     1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)), 
     1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)),
     10 // bins
     );

  Cspacegroup& oSpace = *m_oCrystal.m_poSpacegroup;
  Crefln *poRefln;

  int nCentric = 0;
  int nRef, nNumRefs;
  double fReso, fResoShell;
  int a3nHKL[3];
  float fIntPlus, fIntMinus;
  float fSigmaPlus, fSigmaMinus, fSigma;
  double fDeltaInt;

  float fISigIThreshold = 0.0001;
  if ("" != sGetEnv("DTREK_RAS_ISIGITHRESHOLD"))
    {
      fISigIThreshold = atof(sGetEnv("DTREK_RAS_ISIGITHRESHOLD").string());
      cout << "INFO:  Refln I/sigI threshold for inclusion in R(as) calculation: " 
           << fISigIThreshold << endl;
    }


  nNumRefs = poReflnlist->nGetNumReflns();
  for (nRef = 0; nRef < nNumRefs; nRef++)
    {
      poRefln     = poReflnlist->poGetRefln(nRef);
      fIntPlus    = poRefln->fGetField(poReflnlist->m_nFI_fIntensityPlus);
      fSigmaPlus  = poRefln->fGetField(poReflnlist->m_nFI_fSigmaIPlus);
      fIntMinus   = poRefln->fGetField(poReflnlist->m_nFI_fIntensityMinus);
      fSigmaMinus = poRefln->fGetField(poReflnlist->m_nFI_fSigmaIMinus);
      fSigma      = poRefln->fGetSigmaI();

      // Only work on this hkl if there is a pair of intensities both + and -
      // and if I/sigI is above a threshhold

      if ( (0.0 < fSigmaPlus) && (0.0 < fSigmaMinus) )
        {
          if (   (fISigIThreshold < (fIntPlus / fSigmaPlus) )
              && (fISigIThreshold < (fIntMinus / fSigmaMinus) )
             )
            {
              // A pair of I+, I- found for analysis

              // Need resolution, need whether centric or acentric

              nCentric = oSpace.nReduceHKL(poRefln, a3nHKL, NULL);

              // if nCentric == 0, then acentric
              // if nCentric >  0, then centric
              // if nCentric <  0, the error (sym absent?)

              fReso         = m_oCrystal.dCalcGetResolution(poRefln);
              //cout << "fReso: " << fReso << endl;
              //fResoShell    = 1.0 / max(1.0e-10, pow(fReso, 0.33333));
              CResoRASBin* poResoBin = m_oRASBins.pGetBinReso(fReso);
              if (NULL != poResoBin)
                {
                  // Inside resolution range 
                  fDeltaInt = fIntPlus - fIntMinus;
                  fDeltaInt = fabs(fDeltaInt);
                  if (0 < nCentric)
                    poResoBin->vAddDeltaC(fDeltaInt / fSigma);
                  else if (0 == nCentric)
                    poResoBin->vAddDeltaA(fDeltaInt / fSigma);
                }
            }
        }
    }

  // Done with all reflns

  m_oRASBins.vTallyUpStats();

  // Use a line similar to the following to print the tables:
  // m_oRASBins.vPrintStats();

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
bool Cscaleaverage::bDoAnomalousStatisticsOnSymEquivGroup(std::vector<Crefln*>& vecSymEquivReflns,
                                                          int nResoField,
                                                          int nPlusMinusFlagField,
                                                          int nCentricityField)
{

    if( vecSymEquivReflns.size() < 2 )
        return false; // need at least two reflections to do the statistics
    
    Cstat       oStatPlusCent;
    oStatPlusCent.vClearRaw();

    Cstat       oStatMinusCent;
    oStatMinusCent.vClearRaw();

    Cstat       oStatPlusAcent;
    oStatPlusAcent.vClearRaw();
    
    Cstat       oStatMinusAcent;
    oStatMinusAcent.vClearRaw();

    Cstat       oStatAll;
    oStatAll.vClearRaw();
    
    double      fGroupStandardDeviation = 0.0;
    
    int     nPlusMinus  = 0;
    int     nCentricity = 0;

    double  fIntensity = 0.0;
    double  fSigma = 0.0;
    double  dWeight = 0.0;
    double  dSigmaSquared = 0.0;

    double  fPlusAcentricSigSq = 0.0;
    double  fPlusCentricSigSq = 0.0;
    double  fMinusAcentricSigSq = 0.0;
    double  fMinusCentricSigSq = 0.0;

    double  fDeltaSigma = 0.0;

    int     nPlusCentSize   = 0;
    int     nPlusAcentSize  = 0;
    int     nMinusCentSize  = 0;
    int     nMinusAcentSize = 0;
    
    Crefln*     poRefln = vecSymEquivReflns[0];
    double  fReso = poRefln->fGetField(nResoField);
    CResoRASBin* poResoBin = m_oRASBins.pGetBinReso(fReso);

    // Average
    for(int ii=0; ii < vecSymEquivReflns.size(); ii++)
    {
        poRefln = vecSymEquivReflns[ii];
        
        nPlusMinus  = poRefln->nGetField(nPlusMinusFlagField);
        nCentricity = poRefln->nGetField(nCentricityField);
        
        fIntensity = poRefln->fGetIntensity();
        fSigma     = poRefln->fGetSigmaI();
        
        if( fIntensity <= 0.0 || fSigma <= 0.0 )
            continue;
        
        //+jwp  Add the weight as 1.0 over the square of sigma as well
        //      so later we can compute a weighted average intensity


        dSigmaSquared =  fSigma * fSigma;
        dWeight    = 1.0 / dSigmaSquared;
        oStatAll.vAddRaw(fIntensity, dWeight);

        if( 1 == nPlusMinus )        // F+
        {
            if( 0 == nCentricity )
            {
                oStatPlusAcent.vAddRaw(fIntensity, dWeight);
                fPlusAcentricSigSq += dSigmaSquared;
            }
            else if( 0 < nCentricity )
            {
                oStatPlusCent.vAddRaw(fIntensity, dWeight);
                fPlusCentricSigSq += dSigmaSquared;
            }
        }
        else if( -1 == nPlusMinus )  // F-
        {
            if( 0 == nCentricity )
            {
                oStatMinusAcent.vAddRaw(fIntensity, dWeight);
                fMinusAcentricSigSq += dSigmaSquared;
            }
            else if( 0 < nCentricity )
            {
                oStatMinusCent.vAddRaw(fIntensity, dWeight);
                fMinusCentricSigSq += dSigmaSquared;
            } 
        }
    }
    
    //+JWP 2007-12-05
    //fGroupStandardDeviation = oStatAll.fStandardDeviationRaw();
    //if( fGroupStandardDeviation == 0.0 )
    //  return false;
    //-JWP 2007-12-05

    nPlusCentSize   = oStatPlusCent.nSize();
    nPlusAcentSize  = oStatPlusAcent.nSize();
    nMinusCentSize  = oStatMinusCent.nSize();
    nMinusAcentSize = oStatMinusAcent.nSize();
    
    if( nPlusCentSize > 0 && nMinusCentSize > 0 )
    {
        //poResoBin->vAddDeltaC(fabs(oStatPlusCent.fAverageRaw() - oStatMinusCent.fAverageRaw())   / fGroupStandardDeviation);
        fDeltaSigma = sqrt(fPlusCentricSigSq/(double)(nPlusCentSize*nPlusCentSize) + fMinusCentricSigSq/(double)(nMinusCentSize*nMinusCentSize));
        if( fDeltaSigma != 0.0 )
            poResoBin->vAddDeltaC(fabs(oStatPlusCent.fAverageRaw() - oStatMinusCent.fAverageRaw())   / fDeltaSigma);
    }
    
    if( nPlusAcentSize > 0 && nMinusAcentSize > 0 )
    {
        //poResoBin->vAddDeltaA(fabs(oStatPlusAcent.fAverageRaw() - oStatMinusAcent.fAverageRaw()) / fGroupStandardDeviation);
        fDeltaSigma = sqrt(fPlusAcentricSigSq/(double)(nPlusAcentSize*nPlusAcentSize) + fMinusAcentricSigSq/(double)(nMinusAcentSize*nMinusAcentSize));
        if( fDeltaSigma != 0.0 )
            poResoBin->vAddDeltaA(fabs(oStatPlusAcent.fAverageRaw() - oStatMinusAcent.fAverageRaw()) / fDeltaSigma);
    }

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Error scaling code.   //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

int Cscaleaverage::nPrintEaddEmul() {
    double fAvgEAdd;
    double fAvgEMul;
    int nx;

    fAvgEAdd = 0.0;
    fAvgEMul = 0.0;
    for (nx=0;nx<m_nNumRefs;nx++) {
        fAvgEAdd += m_pfErrorModelEAdd[nx];
        fAvgEMul += m_pfErrorModelEMul[nx];
    };
    fAvgEAdd /= max(1,m_nNumRefs);
    fAvgEMul /= max(1,m_nNumRefs);
  

    printf("Current error model\n");
    printf("----------------------\n");
    printf(" Emul            %.3f\n", fAvgEMul);
    printf(" Eadd            %.3f\n", fAvgEAdd);
    printf(" |ChiSquared|    %.3f\n", m_fChiSquared);
    printf("----------------------\n");
    return 0;
};


int Cscaleaverage::nChiPlot() {
    int nRefSort;
    int nRef;
    int nx;
    double fReso;
    int nBinReso;
    itr<doublearray> aafPointX;
    itr<doublearray> aafPointY;

    for (nx = 0; nx < m_nNumBinsReso;nx++) {
        -aafPointX[nx];
        -aafPointY[nx];
    };
    
    for (nRefSort=0;nRefSort<m_nNumRefs;nRefSort++) {
        
        nRef=m_pnIndex[nRefSort];
        fReso    = m_oReflnlist[nRef].fGetField(m_nFI_f2STLsq);
        fReso    = fReso * sqrt(fReso);
        nBinReso = (int) ((fReso - m_a2fRangeReso[0]) / max(1e-10,m_fSlopeReso));
        nBinReso = min(nBinReso,m_nNumBinsReso-1);
        nBinReso = max(nBinReso,0);
        if (m_pfDev2[nRefSort]) {
            aafPointX[nBinReso] + (double) ( m_fLastEMulUsed*m_oReflnlist[nRef].fGetSigmaI()); // (double) fCalcSigma(m_oReflnlist[nRef],nRef);
            aafPointY[nBinReso] + (double) ABS(m_pfDev2[nRefSort]);
        };
    };

    for (nx = 0; nx < nBinReso; nx++) {
        nPlot((nx==0)?"scatterplot.txt":"",aafPointY[nx],&aafPointX[nx]);
    };
    return 0;
};



int Cscaleaverage::nReportRejects()
{
    if (!m_nNumRefs)
        return 1;

    printf("\n");
    printf("%8d reflections in data set\n",m_nNumRefs);
    printf("%8d reflections rejected (|ChiSq| > %6.2f)\n",m_nTotalAdjustedChiRejects,m_fMaxDevReject);
    printf("%8d reflections total rejected (%6.2f%%  |Deviation|/sigma > %6.2f)\n",m_nTotalAllRejects,100.0*((float) m_nTotalAllRejects)/max(1,m_nNumRefs),m_fDevReject);
    printf("%8d reflections excluded from scaling/absorption (I/sig <=%5.2f)\n",m_nTotalSigmaExcluded,m_fSigma);
    printf("\n");
    fflush(stdout);
    return 0;
};





double Cscaleaverage::fCalcSigma(Crefln& oRefln, int nRefNumber)
{
    double      fIntensity = oRefln.fGetField(m_nFI_nIntensityI);
    
    if( fIntensity < m_fErrorModelEAddMinIntensity )
        fIntensity = m_fErrorModelEAddMinIntensity;
    
    double      fSigma = oRefln.fGetField(m_nFI_nSigmaI);

    if (0.0 <= m_pfErrorModelEAdd[nRefNumber])
      m_fLastEAddUsed = m_pfErrorModelEAdd[nRefNumber];

    if (0.0 < m_pfErrorModelEMul[nRefNumber])
      m_fLastEMulUsed = m_pfErrorModelEMul[nRefNumber];
    /***
    else
      cout << "MULNEG " << nRefNumber << "E,A: "
           << m_pfErrorModelEAdd[nRefNumber] << ", "
           << m_pfErrorModelEMul[nRefNumber] << endl << flush;
    ***/
    double      fBackgroundSigma = 0.0;
    if( m_bUseBackgroundSigma )
    {
        fBackgroundSigma = oRefln.fGetField(m_oReflnlist.m_nFI_fBackgroundSigma);
    }

    double      dSigma = max(0.0, fSigma * fSigma - fBackgroundSigma * fBackgroundSigma) * m_fLastEMulUsed * m_fLastEMulUsed;
        
    dSigma += fIntensity * fIntensity * m_fLastEAddUsed * m_fLastEAddUsed + fBackgroundSigma * fBackgroundSigma;
    
    dSigma = sqrt(dSigma);

    return dSigma;
}


/*
    This routine attempts to find the value for the chi rejection that rejects a certain percentage of the data.

*/

double Cscaleaverage::fFindDevRejection(double fRejectPercent,double* pfChiSq)
{
    int* pnSort;
    int nx;
    int nRefSort;
    int nRefStart,nRefEnd;
    int nIndex;
    double fDevRejection;
    double fSigma;
    double fIntensity;
    double fGroupAvg;
    double f0;
//+2010-07-20
    bool   bUseSigmaWeights = USE_SIGMA_WEIGHTS;
    double fReflnWeight = 1.0;
//-2010-07-20


    for (nRefSort = 0; nRefSort < m_nNumRefs; nRefSort++)
    {
        nRefStart = nRefEnd = nRefSort;
    
        while ((nRefEnd+1 < m_nNumRefs) && (m_pnGroup[nRefStart]==m_pnGroup[nRefEnd+1]))
            nRefEnd++;


        if( nRefStart < nRefEnd )
        {
            f0 = 0.0;
            fGroupAvg = 0.0;
            
            for (nIndex = nRefStart;nIndex <= nRefEnd; nIndex++)
            {
                fIntensity = m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity();
                
                fSigma  = fCalcSigma(m_oReflnlist[m_pnIndex[nIndex]],m_pnIndex[nIndex]);
                
		if (bUseSigmaWeights) fReflnWeight = max(1.0, (fSigma*fSigma));
                fGroupAvg += fIntensity / fReflnWeight;
                
             
                f0 += 1.0/max(1.0,fSigma*fSigma);
            }
            


            if( f0 )
                fGroupAvg /= f0;            
            
            for(nIndex = nRefStart;nIndex <= nRefEnd; nIndex++)
            {
                fIntensity = m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity();
             
                m_pfDev2[nIndex] = (float)(fIntensity - fGroupAvg);
            }
                
            for(nIndex = nRefStart;nIndex <= nRefEnd; nIndex++)
            {
                fSigma = fCalcSigma(m_oReflnlist[m_pnIndex[nIndex]],m_pnIndex[nIndex]);
            
                if (m_pfDev2[nIndex] == 0.0)
                    m_pfDev[nIndex] = 0.0;
                else
                    m_pfDev[nIndex] = m_pfDev2[nIndex]*m_pfDev2[nIndex]/max(1.0,fSigma*fSigma);
            }
        }

        nRefSort = nRefEnd;
    }


    pnSort = new int[m_nNumRefs];

    for (nx=0;nx<m_nNumRefs;nx++)
        pnSort[nx] = nx;
    
    g_pfCmpFloats = m_pfDev;
    qsort(pnSort,m_nNumRefs,sizeof(int),float_cmp_rel);

    fDevRejection = m_pfDev[pnSort[max(0,(int) (m_nNumRefs - 1 - fRejectPercent*m_nNumRefs))]];
    if (fDevRejection > m_fMaxDevReject) 
        fDevRejection = m_fMaxDevReject;


    delete[] pnSort;

    if( pfChiSq )
    {
        double fChi;
        double fTotalChi;           
        int nChi;
        int nTotalChi;
        
        fTotalChi = 0.0;
        nTotalChi = 0;
        for (nRefSort = 0; nRefSort < m_nNumRefs; nRefSort++) {
            nRefStart = nRefEnd = nRefSort;
            while ((nRefEnd+1 < m_nNumRefs) && (m_pnGroup[nRefStart]==m_pnGroup[nRefEnd+1]))
                nRefEnd++;
            fChi = 0.0;
            nChi = 0;
            for (nIndex = nRefStart;nIndex <= nRefEnd; nIndex++) {
                fSigma = fCalcSigma(m_oReflnlist[m_pnIndex[nIndex]],m_pnIndex[nIndex]);
                f0 = m_pfDev2[nIndex]*m_pfDev2[nIndex]/max(1.0,fSigma*fSigma);
                if ((m_pfDev2[nIndex] != 0.0) && (f0<fDevRejection)) {
                    fChi += f0;
                    nChi ++;
                };
            };
            if (nChi>1) {
                fTotalChi += fChi;
                nTotalChi += nChi - 1;
            };
            nRefSort = nRefEnd;
        };
        if (nTotalChi)
            *pfChiSq = fTotalChi/nTotalChi;
        else
            *pfChiSq = 0;
    };


    return fDevRejection;
};
    



int Cscaleaverage::nRejectBatches(const Cstring& sBatchTemplateStart,const Cstring& sBatchTemplateEnd) {
    int nx;
    for (nx=0; nx<m_nNumBatches; nx++) {
        if (!sBatchTemplateEnd.length()) {
            if (m_ptBatchInfo[nx].m_sBatchName.bRegularMatch(sBatchTemplateStart)) {
                if (!m_ptBatchInfo[nx].m_bBatchRejectUser) {
                    printf("Rejecting batch '%s'\n",m_ptBatchInfo[nx].m_sBatchName.string());
                };
                m_ptBatchInfo[nx].m_bBatchRejectUser = TRUE;
            };
        } else {
            if ((sBatchTemplateStart <= m_ptBatchInfo[nx].m_sBatchName) &&
                (sBatchTemplateEnd >= m_ptBatchInfo[nx].m_sBatchName)) {
                if (!m_ptBatchInfo[nx].m_bBatchRejectUser) {
                    printf("Rejecting batch '%s'\n",m_ptBatchInfo[nx].m_sBatchName.string());
                };
                m_ptBatchInfo[nx].m_bBatchRejectUser = TRUE;
            };
        };
    };

    return 0;
};

int Cscaleaverage::nCompletenessCutoff(double fCutoff,Cstring& sSort1,Cstring& sSort2,bool bReverse) {
    int* pnSort;
    int* pnBeforeRejected;


    Cstring sValue;
    static int nDup = -1;

    int nx,ny;

    pnSort = new int[m_nNumRefs];
    pnBeforeRejected = new int[m_nNumRefs];

    eReflnFieldType eType1;
    eReflnFieldType eType2;
    int nIndex1,nIndex2;
    
    if (sSort1.GetAt(0)=='s')
        eType1 = eReflnField_Cstring_type;
    else if (sSort1.GetAt(0)=='n')
        eType1 = eReflnField_int_type;
    else if (sSort1.GetAt(0) == 'f')
        eType1 = eReflnField_float_type;    
    else
        return 1;
    nIndex1 = m_oReflnlist.nGetFieldIndex(sSort1);
    if (nIndex1<0) {
        printf("Could not find field '%s'\n",sSort1.string());
        return 1;
    };
    
    if (sSort2.length()) {
        if (sSort2.GetAt(0)=='s')
            eType2 = eReflnField_Cstring_type;
        else if (sSort2.GetAt(0)=='n')
            eType2 = eReflnField_int_type;
        else if (sSort2.GetAt(0) == 'f')
            eType2 = eReflnField_float_type;    
        else
            return 1;        
        nIndex2 = m_oReflnlist.nGetFieldIndex(sSort2);
        if (nIndex2<0) {
            printf("Could not find field '%s'\n",sSort2.string());
            return 1;
        };

    } else
        nIndex2 = -1;


    printf("\n\nCalculating cut-off value in data for a completeness of %.2lf percent\n",fCutoff);
    printf("\nSorting on Field(s): %s %s...\n",sSort1.string(),sSort2.string());
    if (nIndex2==-1) 
        m_oReflnlist.vSort(eType1,nIndex1,pnSort);
    else
        m_oReflnlist.vSort2(eType1,nIndex1,eType2,nIndex2,pnSort);
    printf("\ndone.\n\n");

    if (bReverse) {
        for (nx=0;nx<m_nNumRefs/2;nx++) {
            ny = pnSort[nx];
            pnSort[nx] = pnSort[m_nNumRefs - nx - 1];
            pnSort[m_nNumRefs - nx - 1] = ny;
        };
    };

    for (nx=0;nx<m_nNumRefs;nx++)
        pnBeforeRejected[nx] = m_pnReject[nx];

    double fRejectLower;
    double fRejectUpper;
    double fCompleteLower;
    double fCompleteUpper;
    double fReject;
    int    nPass;
    int    nRef;

    nPass = 0;
    fRejectLower = 0.0;
    fRejectUpper = 1.0;

    printf("Finding completeness cut-off point...\n");
    do {
        if (nPass==0)
            fReject = 1.0;
        else if (nPass==1)
            fReject = 0.0;
        else 
            fReject = (fRejectLower + fRejectUpper)/2;
        
        for (nx=0;nx<m_nNumRefs;nx++) {
            nRef = pnSort[nx];
            if (((float) nx)/max(1,m_nNumRefs)>=1.0 - fReject)
                m_pnReject[m_pnRevIndex[nRef]] = 1;
            else
                m_pnReject[m_pnRevIndex[nRef]] = pnBeforeRejected[m_pnRevIndex[nRef]];
        };


        // Duplicate stdout file handle.

        if (nDup == -1) 
        {
#ifdef VC9
            nDup = _dup( 1 );
#else
            nDup = dup( 1 );
#endif
        }
        // Force stdin to point to stdout.  This will disable printing in the nTableCompleteness() function.
#ifdef VC9
        _dup2(0,1);
#else
        dup2(0,1);
#endif

        fCalcChiSquared(false);
        nTableCompleteness();

        if (nPass==0)
            fCompleteLower = m_fOverallPercentComplete*100.0;
        else if (nPass==1)
            fCompleteUpper = m_fOverallPercentComplete*100.0;
        else if (m_fOverallPercentComplete*100.0>fCutoff) {
            fCompleteLower = m_fOverallPercentComplete*100.0;
            fRejectLower = fReject;
        } else {
            fCompleteUpper = m_fOverallPercentComplete*100.0;
            fRejectUpper = fReject;
        };
        nPass++;

    } while ((nPass<2) || (fRejectUpper - fRejectLower>0.001));

    // Reinstate output.
    if (nDup>=0)

#ifdef VC9
        _dup2(nDup,1);
#else
        dup2(nDup,1);
#endif

    printf("Completeness Output.\n");
    nTableCompleteness();
    printf("\n\n");


    for (nx=0;nx<m_nNumRefs;nx++)
        m_pnReject[nx] = pnBeforeRejected[nx];

    printf("Percent Completeness of %.2lf found rejecting %.1lf percent of the data.\n",
        (fCompleteLower+fCompleteUpper)*0.5,100.0*(fRejectLower+fRejectUpper)*0.5);

    int nRefSortAbove;
    int nRefSortBelow;
    int nRefSort;
    bool bRefsEqual;

    nRefSort = nRefSortAbove = nRefSortBelow = min(m_nNumRefs-1,max(0,(int) (m_nNumRefs*(1.0-fReject))));
    sValue = "";
    do {
        nRefSortAbove++;
        if (nRefSortAbove>=m_nNumRefs)
            break;
        bRefsEqual = false;
        if (eType1==eReflnField_Cstring_type) {
            if (!sValue.length())
                sValue = m_oReflnlist[pnSort[nRefSort]].sGetField(nIndex1);
            if (m_oReflnlist[pnSort[nRefSortAbove]].sGetField(nIndex1) ==
                m_oReflnlist[pnSort[nRefSort]].sGetField(nIndex1))
                bRefsEqual = true;
        } else if (eType1==eReflnField_int_type) {
            if (!sValue.length())
                sValue += m_oReflnlist[pnSort[nRefSort]].nGetField(nIndex1);
            if (m_oReflnlist[pnSort[nRefSortAbove]].nGetField(nIndex1) ==
                m_oReflnlist[pnSort[nRefSort]].nGetField(nIndex1))
                bRefsEqual = true;
        } else if (eType1==eReflnField_float_type) {
            if (!sValue.length())
                sValue += m_oReflnlist[pnSort[nRefSort]].fGetField(nIndex1);
            if (m_oReflnlist[pnSort[nRefSortAbove]].fGetField(nIndex1) ==
                m_oReflnlist[pnSort[nRefSort]].fGetField(nIndex1))
                bRefsEqual = true;
        };
    } while (bRefsEqual);
    nRefSortAbove--;

    do {
        nRefSortBelow--;
        if (nRefSortBelow<0)
            break;
        bRefsEqual = false;
        if (eType1==eReflnField_Cstring_type) {
            if (m_oReflnlist[pnSort[nRefSortBelow]].sGetField(nIndex1) ==
                m_oReflnlist[pnSort[nRefSort]].sGetField(nIndex1))
                bRefsEqual = true;
        } else if (eType1==eReflnField_int_type) {
            if (m_oReflnlist[pnSort[nRefSortBelow]].nGetField(nIndex1) ==
                m_oReflnlist[pnSort[nRefSort]].nGetField(nIndex1))
                bRefsEqual = true;
        } else if (eType1==eReflnField_float_type) {
            if (m_oReflnlist[pnSort[nRefSortBelow]].fGetField(nIndex1) ==
                m_oReflnlist[pnSort[nRefSort]].fGetField(nIndex1))
                bRefsEqual = true;
        };
    } while (bRefsEqual);
    nRefSortBelow++;


    printf("\nIn sorted list, Field '%s' has the same value of '%s' \n"
           "in a block of %d reflections around the cut-off point.\n"
           "The cut-off point will reject %.1lf%% of this block\n",
           sSort1.string(),
           sValue.string(),
           nRefSortAbove-nRefSortBelow+1,
           100.0*(nRefSortAbove - nRefSort + 1)/max(1,(nRefSortAbove-nRefSortBelow+1)));

    delete[] pnSort;
    delete[] pnBeforeRejected;
    return 0;
};


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//////////////////////           Tables            //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

Cstring Cscaleaverage::ms_sf2STLsq    = "f2STLsq";
Cstring Cscaleaverage::ms_snBatchIdx  = "nBatchIndex";

int Cscaleaverage::nInitArrays() {
  // Should have reflection list and a few other items, so ready
  // to compute scale factors for batches in the reflnlist

  int      nx,ny,nz;             // Loop variables
  int      nStat;                // Local status
  int      nNumReflns;           // Temp variable holding number of reflns in list
  float    fTemp;                // Temp variable

  Crefln  *poRefln;              // A pointer used for readability of code
  Cstring  sBatchPrev, sBatchCurr; // Some strings to hold batch names

  int      nFI[5];               // Indexes of 5 batch variables

  bool     bNoMatch;
  int      nBatchIndex;
  int      nReflectionsExtinct;

  float    a3fX[3];              // Local var for recip lattice vector Bh = x
  float    a3fH[3];              // hkl of refln
  int      nAlloc;

  nNumReflns   = m_nNumRefs;
  m_nFI_f2STLsq   = m_oReflnlist.nExpandGetField(ms_sf2STLsq);
  m_nFI_nBatchIdx = m_oReflnlist.nExpandGetField(ms_snBatchIdx);

  int *pnDeleteArray;
  pnDeleteArray = new int [nNumReflns];
  
  if (m_pnIndexBatch)
      delete[] m_pnIndexBatch;
  m_pnIndexBatch = new int[nNumReflns];

  printf("Sorting on batch ID ...\n");
  m_oReflnlist.vSort(eReflnField_Cstring_type,m_oReflnlist.m_nFI_sBatch,m_pnIndexBatch);
  printf("...done.\n");

  
  if (m_nNumRefs != 0) {
      m_nNumBatches = 1;
      sBatchPrev = m_oReflnlist[m_pnIndexBatch[0]].sGetField(m_oReflnlist.m_nFI_sBatch);
  } else
      m_nNumBatches = 0;

  // Find number of different batches in the reflection list
  for (nx = 1; nx< m_nNumRefs; nx++) {
      sBatchCurr = m_oReflnlist[m_pnIndexBatch[nx]].sGetField(m_oReflnlist.m_nFI_sBatch);
      if (sBatchCurr != sBatchPrev) {
          m_nNumBatches++;
          sBatchPrev = sBatchCurr;
      };
  };

  if (m_pnNumOverlap)
      delete[] m_pnNumOverlap;
  m_pnNumOverlap    = new int [m_nNumBatches * (m_nNumBatches+1)];


  // Initialize some variables

  nAlloc = m_nNumBatches;
  for (nx= 0; nx< 5; nx++)
    {
      if (m_pfMin[nx])
          delete[] m_pfMin[nx];
      if (m_pfMax[nx])
          delete[] m_pfMax[nx];
      m_pfMin[nx]    = new float [nAlloc];
      m_pfMax[nx]    = new float [nAlloc];
    }

  if (NULL != m_ptBatchInfo)
      delete [] m_ptBatchInfo;

  m_ptBatchInfo = new tagBatch[nAlloc]; 

  // The following fields will be used to divide batches into subbatches
  // Maybe use a member variable and a method to specify these

  nFI[0] = m_oReflnlist.m_nFI_fObsPx0;
  if (0 > nFI[0])
    {
      cout << "No observed pixel 0 field!\n";
      nFI[0] = 0;
//      return (-2);
    }
  nFI[1] = m_oReflnlist.m_nFI_fObsPx1;
  if (0 > nFI[1])
    {
      cout << "No observed pixel 1 field!\n";
      nFI[1] = 0;
//      return (-2);
    }
  nFI[2] = m_oReflnlist.m_nFI_fObsRotMid;
  if (0 > nFI[2])
    {
      cout << "No observed rot mid field!\n";
      nFI[2] = 0;
//      return (-2);
    }

  nFI[3]          = m_nFI_f2STLsq;
  nFI[4]          = m_oReflnlist.m_nFI_fIntensity;

  
  sBatchPrev   = "";
  nBatchIndex  = -1;


  float f2STLsq, f2STLsqMin, f2STLsqMax;
  f2STLsqMin = 0.0;
  f2STLsqMax = 1000.0;
  if (0.0 < m_fResolutionMax)
    f2STLsqMax = 1.0f/(float)(m_fResolutionMax * m_fResolutionMax);
  if (0.0 < m_fResolutionMin)
    f2STLsqMin = 1.0f/(float)(m_fResolutionMin * m_fResolutionMin);


  cout << "Examining batch IDs and reflection resolution limits...\n" << flush;
  m_nNumBatches = 0;
  nReflectionsExtinct = 0;
  for (nx= 0; nx< nNumReflns; nx++)
    {
      poRefln  = m_oReflnlist.poGetRefln(m_pnIndexBatch[nx]);

      // Compute reciprocal lattice vector for potential scaling
      // May be a factor of 4 or 0.25 here?

      a3fH[0] = (float) poRefln->nGetH();
      a3fH[1] = (float) poRefln->nGetK();
      a3fH[2] = (float) poRefln->nGetL();
      vMulMat3DVec3D(m_a3x3fBMat, a3fH, a3fX);

      f2STLsq = fDot3D(a3fX, a3fX);
      poRefln->vSetField(m_nFI_f2STLsq, f2STLsq);
      if ( (f2STLsq > f2STLsqMax) || (f2STLsq < f2STLsqMin) )
        {
          // Mark reflection for deletion?

          pnDeleteArray[m_pnIndexBatch[nx]] = 1;
        }
      else if ( m_bExcludeSysAbsences && 
                m_oCrystal.m_poSpacegroup->bIsExtinct(poRefln->nGetH(),
                                                      poRefln->nGetK(),
                                                      poRefln->nGetL()) ) 
      {
          pnDeleteArray[m_pnIndexBatch[nx]] = 1;
          nReflectionsExtinct++;
      }
      else
        {
          // Passes resolution limits, so work with it
        
          pnDeleteArray[m_pnIndexBatch[nx]] = 0;

          sBatchCurr = poRefln->sGetField(m_oReflnlist.m_nFI_sBatch);
          bNoMatch   = (sBatchPrev != sBatchCurr);
          for (ny = 0; (ny < m_nNumBatches) && bNoMatch; ny++)
            {
              // Check that batchname is in not list of previously found
              // batch names
              // We need to do this faster, maybe with some kind of hashing
              // or other direct searching

              bNoMatch = (m_ptBatchInfo[ny].m_sBatchName != sBatchCurr);
              if (!bNoMatch)
                nBatchIndex = ny;        // Save batch index
            }
          if (bNoMatch)
            {
              // Still no match, so add this batch name to the list of batches

              nBatchIndex = m_nNumBatches;
              m_nNumBatches++;
              if (m_nNumBatches > nAlloc)
              {
                  //This should not happen!
                  printf("Fatal error!\n");
                  exit(1);
              };

              // Initialize this batches info

              for (nz = 0; nz < 5; nz++)
                {
                  m_pfMin[nz][nBatchIndex] = poRefln->fGetField(nFI[nz]);
                  m_pfMax[nz][nBatchIndex] = poRefln->fGetField(nFI[nz]);
                }
              m_ptBatchInfo[nBatchIndex].m_sBatchName = sBatchCurr;
              m_ptBatchInfo[nBatchIndex].m_bBatchReject = FALSE;
              m_ptBatchInfo[nBatchIndex].m_bBatchRejectUser = FALSE;
              m_ptBatchInfo[nBatchIndex].m_fBatchScale = 1.0;
              m_ptBatchInfo[nBatchIndex].m_fBFactorScale = 0.0;
              m_ptBatchInfo[nBatchIndex].m_fBFactorConst = 0.0;
              m_ptBatchInfo[nBatchIndex].m_nNumReflns = 0;
              m_ptBatchInfo[nBatchIndex].m_nBatch = nGetBatchNumber(sBatchCurr,m_ptBatchInfo[nBatchIndex].m_sScan);              

      } // endif bNoMatch

          sBatchPrev = sBatchCurr;

          // Save batch index in refln

          poRefln->vSetField(m_nFI_nBatchIdx, nBatchIndex);

          // Save min and max of 5 refln fields for possible grouping

          for (nz = 0; nz < 5; nz++)
            {
              fTemp = poRefln->fGetField(nFI[nz]);
              m_pfMin[nz][nBatchIndex] = min(m_pfMin[nz][nBatchIndex], fTemp);
              m_pfMax[nz][nBatchIndex] = max(m_pfMax[nz][nBatchIndex], fTemp);
            }
          m_ptBatchInfo[nBatchIndex].m_nNumReflns++;
        }
    }  // i, end of looping through reflns
  cout << "   " << m_nNumBatches << " batch IDs were found in the input reflnlist.\n"
       << "...done.\n" << flush;

  if (nReflectionsExtinct) {
      printf("INFO:  %d reflections found and removed which are extinct in the spacegroup!\n",nReflectionsExtinct);
  };

  int nNumDel =m_oReflnlist.nDelete(1, pnDeleteArray);

  // Reset the number of reflections in the list.

  m_nNumRefs = m_oReflnlist.nGetNumReflns();
  
  if( 0 == m_nNumRefs )
  {
      printf("\n\nWARNING: No reflections left in the reflection list.\n");
      printf("Please check your input resolution limits and other rejection criteria.\n\n");
      return 1;
  }
  
  if (m_pfErrorModelEAdd)
      delete[] m_pfErrorModelEAdd;
  m_pfErrorModelEAdd = NULL;

  if (m_pfErrorModelEMul)
      delete[] m_pfErrorModelEMul;
  m_pfErrorModelEMul = NULL;
  
  vInitErrorModel();

  cout << "\nReflections exceeding resolution limits: " << nNumDel
       << "\n                                   Kept: " << m_nNumRefs
       << " for scaling.\n" << flush;
  

  // Copy the original intensities over into an array.
  // If we are using texsan output, we will want to update with original values at end.
  for (nx=0;nx<m_oReflnlist.nGetNumReflns();nx++) {
      m_pfOrigIntensity[nx] = m_oReflnlist[nx].fGetIntensity();
  };

  printf("Sorting on intensity ...\n");

  m_pnIndexIntensity = new int[m_nNumRefs];
  m_oReflnlist.vSort(eReflnField_float_type,0,m_pnIndexIntensity);

  int       n99PercentOfReflectionCount = (int)ceil(m_oReflnlist.nGetNumReflns() * 0.99); 
  
  int       nSortedIndexOfMaximumIntensityReflection = m_pnIndexIntensity[ max(0, n99PercentOfReflectionCount-1) ];  
  
  m_fMaxScaleIntensity = m_oReflnlist[nSortedIndexOfMaximumIntensityReflection].fGetIntensity();
  
  printf("...done.\n");

  // Check to if nRefineBatch field is present.
  // This is our way of enforcing backwards compatibility with older d*TREK files which used fBackgroundSigma to store
  // Print a warning if this is the case.
  if ((m_oReflnlist.m_nFI_nRefineBatch < 0) || (m_oReflnlist.m_nFI_fBackgroundSigma < 0)) {
      printf("WARNING:  The input file appears to have come from an older version of d*TREK integration.\n");
      printf("WARNING:  Some changes in processing will occur.\n");
      m_bUseBackgroundSigma = FALSE;
  } else 
      m_bUseBackgroundSigma = FALSE; //TRUE;
  

  // Resort the data on batch, because we deleted some reflections.
  if (nNumDel != 0)
    {
      cout << "Resorting on batch ID ...\n" << flush;
      m_oReflnlist.vSort(eReflnField_Cstring_type,m_oReflnlist.m_nFI_sBatch,m_pnIndexBatch);
      cout << "...done.\n" << flush;
    }


  delete [] pnDeleteArray;

  cout << "Sorting and reducing reflnlist to asymmetric unit ..." << endl << flush;
  m_nStat = m_oReflnlist.nReduce(m_oCrystal, m_bScaleAnom);
  cout << "...done.\n" << flush;
  m_pnIndex = m_oReflnlist.pnGetSortIndex();
  if (m_nStat) {
     printf("Could not reduce to asymmetric unit!");
    };
  
  //m_oReflnlist.nWrite("TEST_ReducedHKL.ref", m_pnIndex);

  if (m_pnRevIndex)
      delete[] m_pnRevIndex;
  m_pnRevIndex = new int[m_nNumRefs];
  for (nx=0;nx<m_nNumRefs;nx++)
      m_pnRevIndex[m_pnIndex[nx]] = nx;

  

  if (0 >= m_nNumRefs)
    {
      cout << "ERROR, no reflections pass resolution limits!\n";
#ifdef SSI_PC
      return (-1);
#else
      exit (-1);
#endif
    }

  // Save min, max of reso and intensity for all batches

  m_a2fRangeReso[0] = m_pfMin[3][0];
  m_a2fRangeReso[1] = m_pfMax[3][0];
  m_a2fRangeInt[0]  = m_pfMin[4][0];
  m_a2fRangeInt[1]  = m_pfMax[4][0];
  for (ny = 1; ny < m_nNumBatches; ny++)
    {
      m_a2fRangeReso[0] = min(m_pfMin[3][ny], m_a2fRangeReso[0]);
      m_a2fRangeReso[1] = max(m_pfMax[3][ny], m_a2fRangeReso[1]);
      m_a2fRangeInt[0]  = min(m_pfMin[4][ny], m_a2fRangeReso[0]);
      m_a2fRangeInt[1]  = max(m_pfMax[4][ny], m_a2fRangeReso[1]);
    }

  m_a2fRangeReso[0] = m_a2fRangeReso[0] * (float)sqrt((double)m_a2fRangeReso[0]);
  m_a2fRangeReso[1] = m_a2fRangeReso[1] * (float)sqrt((double)m_a2fRangeReso[1]);

  m_fSlopeReso = (m_a2fRangeReso[1] - m_a2fRangeReso[0])
                  / (float) max(1,m_nNumBinsReso);

  // Make intensity really Intensity / SigmaI
  m_a2fRangeInt[1] = 20.0;
  m_a2fRangeInt[0] = 0.0;
  m_fSlopeInt  = (m_a2fRangeInt[1] - m_a2fRangeInt[0])
                  / (float) max(1,m_nNumBinsInt);


  { 
      if (m_pfRevIndexReso)
          delete[] m_pfRevIndexReso;
      m_pfRevIndexReso = new float[m_nNumRefs];
      float* pfReso;
      pfReso = new float[m_nNumRefs];
      m_pnIndexReso = new int[m_nNumRefs];
            
      for (nx = 0; nx < m_nNumRefs; nx++)
      {
          double fReso;
          fReso    = m_oReflnlist[nx].fGetField(m_nFI_f2STLsq);
          pfReso[nx] = fReso;
          m_pnIndexReso[nx] = nx;
      };

      g_pfCmpFloats = pfReso;
      qsort(m_pnIndexReso,m_nNumRefs,sizeof(int),float_cmp_rel);
      for (nx = 0; nx < m_nNumRefs; nx++)  {
          m_pfRevIndexReso[m_pnIndexReso[nx]] = nx/((float) max(1,m_nNumRefs));
      };

      delete[] pfReso;
  };


  // Clean-up:  Delete new'd variables


  // Initialize other variables.
  if (m_ptStatsReso)
      delete[] m_ptStatsReso;
  m_ptStatsReso     = new tagStatsScaleAverage [m_nNumBinsReso+1];
  if (m_ptStatsInt)
      delete[] m_ptStatsInt;
  m_ptStatsInt      = new tagStatsScaleAverage [m_nNumBinsInt+1];
  if (m_ptStatsBatch)
      delete[] m_ptStatsBatch;
  m_ptStatsBatch    = new tagStatsScaleAverage [m_nNumBatches+1];


  nStat = 0;
  return (nStat);
};

void Cscaleaverage::vInitErrorModel()
{
    if (!m_pfErrorModelEAdd)
        m_pfErrorModelEAdd = new double[m_nNumRefs];
    
    if (!m_pfErrorModelEMul)
        m_pfErrorModelEMul = new double[m_nNumRefs];

    double      dErrAdd = 0.0;
    double      dErrMul = 2.5;
    
    if (0.0 < m_fUserEMul)
      dErrMul = m_fUserEMul;
    
    if (0.0 < m_fUserEAdd)
      dErrAdd = m_fUserEAdd;
 
    for(int nx=0;nx< m_nNumRefs; nx++)
    {
      m_pfErrorModelEMul[nx] = dErrMul;
      m_pfErrorModelEAdd[nx] = dErrAdd;
    }
}


int Cscaleaverage::nTableChiDistribution() {
    int nx,ny,nz;

    printf("\n\nDistribution of reduced ChiSquared (|ChiSq|)");
    printf("\n-----------------------------------------------------------------------------\n");
    printf(" |ChiSq|       ");
    for (nx=0;nx<m_nChiBins;nx++) {
        if (nx==0)
            printf(" <%4.2f ",m_pfChiDistribBins[nx]);
        else if (nx==m_nChiBins-1)
            printf(">%5.2f ",m_pfChiDistribBins[nx]);
        else
            printf("%6.2f ",m_pfChiDistribBins[nx]);
    };
    printf("\n-----------------------------------------------------------------------------\n");
    printf(" Num unique ref");
    for (nx=0;nx<m_nChiBins;nx++)
        printf("%6d ",m_pnChiDistribGroups[nx]);
    printf("\n Num reflns    ");
    for (nx=0;nx<m_nChiBins;nx++)
        printf("%6d ",(int) (m_pnChiDistribGroups[nx]*m_pfChiDistribMultiplicity[nx]));
    printf("\n Avg mult      ");
    for (nx=0;nx<m_nChiBins;nx++)
        printf("%6.1f ",m_pfChiDistribMultiplicity[nx]);

    for (nx=0,nz=0;nx<m_nChiBins;nx++) {
        nz+= m_pnChiDistribGroups[nx];
    };
    printf("\n Shell %%unique ");
    for (nx=0,ny=0;nx<m_nChiBins;nx++) {
        ny = m_pnChiDistribGroups[nx];
        if (nz>0)
            printf("%6.1f ",(100.0*(float)ny)/(float)nz);
        else
            printf("   --- ");
    };
    printf("\n Cumul %%unique ");
    for (nx=0,ny=0;nx<m_nChiBins;nx++) {
        ny+=m_pnChiDistribGroups[nx];
        if (nz>0)
            printf("%6.1f ",(100.0*(float)ny)/(float)nz);
        else
            printf("   --- ");
    };
    printf("\n-----------------------------------------------------------------------------\n\n");
    return 0;
};




int Cscaleaverage::nTableIntensityResolutionPosition()
{
  // Should have reflection list and a few other items, so ready
  // to compute scale factors for batches in the reflnlist

  int      ny,nz;             // Loop variables
  int      nStat;                // Local status

  Cstring  sBatchPrev, sBatchCurr; // Some strings to hold batch names

  int      nFI[5];               // Indexes of 5 batch variables

  // The following fields will be used to divide batches into subbatches
  // Maybe use a member variable and a method to specify these

  nFI[0] = m_oReflnlist.m_nFI_fObsPx0;
  if (0 > nFI[0])
    {
      cout << "No observed pixel 0 field!\n";
      nFI[0] = 0;
//      return (-2);
    }
  nFI[1] = m_oReflnlist.m_nFI_fObsPx1;
  if (0 > nFI[1])
    {
      cout << "No observed pixel 1 field!\n";
      nFI[1] = 0;
//      return (-2);
    }
  nFI[2] = m_oReflnlist.m_nFI_fObsRotMid;
  if (0 > nFI[2])
    {
      cout << "No observed rot mid field!\n";
      nFI[2] = 0;
//      return (-2);
    }

  nFI[3]          = m_nFI_f2STLsq;
  nFI[4]          = m_oReflnlist.m_nFI_fIntensity;

  if (1 < m_nVerbose)
    {
      if (1 < m_nNumBatches)
        {
          printf("\nThere are %d different batches ", m_nNumBatches);
        }
      else
        {
          printf("\nThere is only %d batch ", m_nNumBatches);
        }
      printf("in the input reflection list.\n");
      printf("\nObserved position limits of the Batches\n");
      char *pcShortLine = "--------------------------------------------------------------------------------\n";
      if ( (0 < nFI[0]) || (0 < nFI[1]) || (0 < nFI[2]) )
        {
          // Print these out only if the fields exist in the input reflnlist

          printf(pcShortLine);
          printf("%6s%7s %22s%22s%22s\n", "Batch", "Num",
             (m_oReflnlist.sGetFieldName(nFI[0],eReflnField_float_type)).string(),
             (m_oReflnlist.sGetFieldName(nFI[1],eReflnField_float_type)).string(),
             (m_oReflnlist.sGetFieldName(nFI[2],eReflnField_float_type)).string());
          printf("%6s%7s %11s%11s%11s%11s%11s%11s\n", "name", "refs",
                 "Min", "Max", "Min", "Max", "Min", "Max");
          printf(pcShortLine);
          for (ny = 0; ny < m_nNumBatches; ny++)
            {
              printf("%6s%7d ", m_ptBatchInfo[ny].m_sBatchName.string(), m_ptBatchInfo[ny].m_nNumReflns);
              for (nz = 0; nz < 3; nz++)
                {
                  printf(" %10.1f %10.1f", m_pfMin[nz][ny], m_pfMax[nz][ny]);
                }
              printf("\n");
            }
          printf(pcShortLine);
        }

      printf("\n\nIntensity and Resolution limits of the Batches (before scaling) \n");
      printf(pcShortLine);
      printf("%6s%7s %22s%22s%22s\n", "Batch", "Num", "Intensity ",
             "Resolution ", "[2sinT/lam]^2 ");
      printf("%6s%7s %11s%11s%11s%11s%11s%11s\n", "name", "refs",
             "Min", "Max", "Min", "Max", "Min", "Max");
      printf(pcShortLine);
      for (ny = 0; ny < m_nNumBatches; ny++)
        {
          printf("%6s%7d ", m_ptBatchInfo[ny].m_sBatchName.string(), m_ptBatchInfo[ny].m_nNumReflns);
          nz = 4;
          printf(" %10.1f %10.1f", m_pfMin[nz][ny], m_pfMax[nz][ny]);
          nz = 3;
          printf(" %10.2lf %10.2lf", 1./max(1e-10,sqrt((double)m_pfMin[nz][ny])), 1./max(1e-10,sqrt((double)m_pfMax[nz][ny])));
          printf(" %10.4f %10.4f", m_pfMin[nz][ny], m_pfMax[nz][ny]);
          printf("\n");
        }
      printf(pcShortLine);
      printf("\n");
      fflush(stdout);
    }

  nStat = 0;
  return (nStat);
}

int Cscaleaverage::nTableMultiplicity()
{
  int     nx;        // Loop counters

  char     *pcLine = "--------------------------------------------------------------------------------\n";

  // Print out multiplicity statistics

  printf("\n\nMultiplicity of observed reflections\n");
          printf(pcLine);
          printf("%17s%8s%6d%7d%7d%7d%7d%7d%7d%7s\n", "Mult | ",
                 "1*", (int)2, (int)3, (int)4, (int)5, (int)6,
                 (int)7, (int)8, ">8");
          printf(pcLine);
          printf("%17s", " Refs | ");
          for (nx = 1; nx < 10; nx++)
            {
              printf("%7d", m_a10nMult[nx]);
            }
      printf("\n");
          printf("%17s", "  |ChiSq| | ");
          for (nx = 1; nx < 10; nx++)
            {
          if ((nx==1) || (m_a10nMult[nx]==0))
              printf("%7s"," ----");
          else
                printf("%7.2f", m_a10fSumChi[nx]/max(1,(m_a10nUsedLE[nx]-m_a10nMultLE[nx])));
            }

          printf("\n%s", pcLine);
          printf("*Reflections with a multiplicity of 1 are not used in\n"
                 " scale factor refinement nor in Rmerge calculations.\n\n");

          // Print out number of overlapping reflns in each batch

          if (0 != m_bComputeOverlap)
            nTableOverlap();

          // Print out number of reflns in each category for each batch

          printf("\n\nIout = Iin*Scale_factor*exp(2B_factor*sin^2(Theta)/lambda^2)\n");
          printf("Reflections in input file\n");
          printf(pcLine);
          printf("%13s %7s %7s %9s %9s %9s %9s %9s\n",
                 "Batch", "Scale", "B",     " Num", "  Num", " Num", "  Num", "    Num");
          printf("%13s %7s %7s %9s %9s %9s %9s %9s\n",
                 " name", "factor","factor","refs", "excluded", "rejs", "ovlps", "singles");
          printf(pcLine);

      int nSumReflns   = 0;
      int nSumExcluded = 0;
      int nSumRejected = 0;
      int nSumUsed = 0;
      int nSumSingles = 0;

          for (nx = 0; nx < m_nNumBatches; nx++)
            {
              nSumReflns   += m_ptStatsBatch[nx].nNumRefs;
              nSumExcluded += m_ptStatsBatch[nx].nNumExcluded;
              nSumRejected += m_ptStatsBatch[nx].nNumRejects;
              nSumUsed     += m_ptStatsBatch[nx].nNumUsed;
              nSumSingles  += m_ptStatsBatch[nx].nNumSingles;

              printf("%13s %7.3f %7.3f %9d %9d %9d %9d %9d\n",
                     m_ptBatchInfo[nx].m_sBatchName.string(),
                     m_ptBatchInfo[nx].m_fBatchScale,
                     m_ptBatchInfo[nx].m_fBFactorScale != 0.0f ? -m_ptBatchInfo[nx].m_fBFactorScale : 0.0f,
                     m_ptStatsBatch[nx].nNumRefs,
                     m_ptStatsBatch[nx].nNumExcluded,
                     m_ptStatsBatch[nx].nNumRejects,
                     m_ptStatsBatch[nx].nNumUsed,
                     m_ptStatsBatch[nx].nNumSingles);
            }
          printf(pcLine);
          printf("%13s %7s %7s %9d %9d %9d %9d %9d\n",
                     "All batches","---","---",
             nSumReflns,nSumExcluded,nSumRejected,nSumUsed,nSumSingles);

  return (0);
}


int Cscaleaverage::nTableRMergeVsBatch()
{
    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    int       nTemp;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";


    printf("\n\nIn the tables below Rmerge is defined as:\n");
    printf("    Rmerge = Sum Sum |Ihi - <Ih>| / Sum Sum <Ih>\n");
    printf("              h   i                  h   i\n");
    printf("    where Ihi is the ith used observation for unique hkl h,\n");
    printf("    and  <Ih> is the mean intensity for unique hkl h.\n\n");
    printf("  Num_Obs = Num_Overlaps + Num_Rejects + Num_Singles\n"
        "  Num_Mults counts only sym. related groups with at least\n"
        "  2 reflections.\n\n"
        "  Exclusions are included in statistics.\n\n"
        "  |ChiSq| and Rmerge are taken over all non-rejected\n"
        "  reflections having a redundancy of at least 2.\n\n"
        "  <I/Sig> and Average are taken over all non-rejected\n"
        "  reflections.\n"
        );

    // List batch statistics

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeBatchRow();
#endif

    printf("\n\nRmerge vs Batch\n");
    printf(pcLine);
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
        "Batch",  "Average", "Num", " Num", " Num", "    Num", "I/sig",
        "Rducd", "Rmerge", "Rmerge");
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
        " name", " counts", "obs", "rejs", "ovlps", "single", "unavg",
        "ChiSq", " batch", " cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsBatch;
    ptSCumul = &m_ptStatsBatch[m_nNumBatches];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    for (nx = 0; nx < m_nNumBatches; nx++, ptStats++)
    {
        (*ptSCumul) += (*ptStats);
        nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
        printf("%14s %7.0f %7d %5d %7d %7d ",
            m_ptBatchInfo[nx].m_sBatchName.string(),
            ptStats->fSum / max(1,nTemp),
            ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
            ptStats->nNumSingles);
        if (0.0 != ptStats->fDenom)
        {
            printf("%6.1f %6.2f %6.3f ",
                ptStats->fSumIoverSig / (float) max(1,nTemp),

                ptStats->fSumChiSq
                / max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups)),
                ptStats->fNumer / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeBatchRow( 
                m_ptBatchInfo[nx].m_sBatchName,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
                ptStats->nNumSingles,
                (Cstring)(ptStats->fSumIoverSig / (float) max(1,nTemp)),
                (Cstring)(ptStats->fSumChiSq / max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups))),
                (Cstring)(ptStats->fNumer / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else {
            printf("%6s %6s %6s ", "---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeBatchRow(
                m_ptBatchInfo[nx].m_sBatchName,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
                ptStats->nNumSingles,
                (Cstring)"---", (Cstring)"---", (Cstring)"---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
            printf("%6.3f\n", ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom));
        else
            printf("%6s\n", "---");
    }
#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeBatchRow();
#endif
    printf(pcLine);
    nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("%14s %7.0f %7d %5d %7d %7d ",
        "All batches",
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
        ptStats->nNumSingles);
    if (0.0 != ptStats->fDenom)
    {
        printf("%6.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));

        printf("%6.2f ", ptStats->fSumChiSq/max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups)));
        printf("%6.3f %6.3f\n", ptStats->fNumer / max(1e-10,ptStats->fDenom),
            ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%6s %6s %6s %6s\n ", "---", "---", "---", "---");
    }

    return 0;
}

int Cscaleaverage::nTableRMergeVsIntensity()
{
    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    float     fLow;
    int       nTemp;
    char      cTemp1, cTemp2;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";

    // List intensity statistics

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeISigRow();
#endif

    printf("\n\nRmerge vs Intensity/SigmaI\n");
    printf(pcLine);
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "Int/sigmaI ",  "Average", "Num", " Num", "   Num"," I/sig", " I/sig",
        "Rducd", "Rmerge", "Rmerge");
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "   range   ", " counts", "obs", "rejs", "mults", " unavg","   avg",
        "ChiSq", " shell", " cumul");

    printf(pcLine);
    ptStats  = &m_ptStatsInt[m_nNumBinsInt-1];
    ptSCumul = &m_ptStatsInt[m_nNumBinsInt];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeInt[1];
    pfSlope  = &m_fSlopeInt;
    for (nx = m_nNumBinsInt-1; nx >= 0; nx--, ptStats--)
    {
        (*ptSCumul) += (*ptStats);
        nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
        cTemp1 = ' ';
        cTemp2 = ' ';
        if ((m_nNumBinsInt-1) == nx) cTemp1 = '>';
        if (0 == nx)                 cTemp2 = '<';
        if ( (0.0 != ptStats->fDenom) && (0 != ptStats->nNumUsed) )
        {
            printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d %7.1f %6.1f %6.2f %6.3f ",
                cTemp2, fLow - *pfSlope, cTemp1, fLow,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                ptStats->fSumIoverSig / (float)max(1,nTemp),
                ptStats->fSumIoverSig2/ (float)max(1,nTemp),
                ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups)),
                ptStats->fNumer / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeISigRow( 
                fLow - *pfSlope, fLow,
                (Cstring)(ptStats->fSum / max(1,nTemp)),
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring)(ptStats->fSumIoverSig / (float)max(1,nTemp)),
                (Cstring)(ptStats->fSumIoverSig2/ (float)max(1,nTemp)),
                (Cstring)(ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups))),
                (Cstring)(ptStats->fNumer / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else
        {
            printf("   %c%2.0f -  %c%2.0f%8s%8d%7d %7d %7s %6s %6s %6s ",
                cTemp2, fLow - *pfSlope, cTemp1, fLow,
                "---",
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                "---","---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeISigRow(
                fLow - *pfSlope, fLow,
                (Cstring)"---",
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                (Cstring)"---", (Cstring)"---", (Cstring)"---", (Cstring)"---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
        {
            printf("%6.3f\n", ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom));
        }
        else
        {
            printf("%6s\n", "---");
        }
        fLow = fLow - *pfSlope;
    }

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeISigRow();
#endif

    cTemp1 = '>';
    cTemp2 = '<';
    ptStats = ptSCumul;
    printf(pcLine);
    nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d ",
        cTemp2, m_a2fRangeInt[0], cTemp1, m_a2fRangeInt[1],
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRefs, 
        ptStats->nNumRejects, 
        ptStats->nNumMults);
    if (0.0 != ptStats->fDenom)
    {

        printf("%7.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));
        printf("%6.1f ", ptStats->fSumIoverSig2/ (float)max(1,nTemp));
        printf("%6.2f ", ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups)));
        printf("%6.3f %6.3f\n", ptStats->fNumer / max(1e-10,ptStats->fDenom),
            ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%7s %6s %6s %6s %6s\n", "---", "---", "---", "---", "---");
    }

    return 0;
};


int Cscaleaverage::nTableRMergeVsResolution(bool bSummary) {

    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    double    fLow;
    double    fResoLow, fResoHigh;
    int       nTemp;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeReslnRow();
#endif

    // List resolution statistics

    printf("\nRmerge vs Resolution\n");
    printf(pcLine);
    printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
        "Resolution  ", "Average",  " Num", "   Num", " I/sig", " I/sig",
        "Rducd", "Model","Rmerge", "Rmerge");
    printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
        "   range    ", " counts", "rejs", "mults", " unavg","   avg",
        "ChiSq", "Eadd*"," shell", " cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsReso;
    ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeReso[0];
    pfSlope  = &m_fSlopeReso;
    for (nx = 0; nx < m_nNumBinsReso; nx++, ptStats++)
    {
        fResoLow              = 1.0 / max(1e-10,pow(fLow, 0.333333));
        fResoHigh             = 1.0 / max(1e-10,pow((fLow + *pfSlope), 0.333333));
        (*ptSCumul) += (*ptStats);
        nTemp = ptStats->nNumSingles + ptStats->nNumUsed;
        if ( (0.0 != ptStats->fDenom) && (0 != ptStats->nNumUsed) )
        {
            printf("%6.2lf -%5.2lf%8.0f%7d %7d %7.1f %6.1f %6.2f %6.2f %6.3f ",
                min(fResoLow, 999.99), fResoHigh,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                ptStats->fSumIoverSig / (float)max(1,nTemp),
                ptStats->fSumIoverSig2/ (float)max(1,nTemp),
                ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups)),
                ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom),
                ptStats->fNumer / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeReslnRow( fResoLow, fResoHigh, (Cstring) (ptStats->fSum / max(1,nTemp)),
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring) (ptStats->fSumIoverSig / (float)max(1,nTemp)),
                (Cstring) (ptStats->fSumIoverSig2/ (float)max(1,nTemp)),
                (Cstring) (ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups))),
                (Cstring) (ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom)),
                (Cstring) (ptStats->fNumer / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else
        {
            printf("%6.2f -%5.2f%8s%7d %7d %7s %6s %6s %6s %6s ",
                min(fResoLow,999.99), fResoHigh,
                "---",
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                "---","---","---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeReslnRow( fResoLow, fResoHigh, (Cstring)"---",
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring) "---",(Cstring) "---", (Cstring) "---", (Cstring) "---",(Cstring) "---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
            printf("%6.3f \n", ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom));
        else
            printf("%6s \n", "---");
        fLow = fLow + *pfSlope;
    }
    printf(pcLine);
    nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("%6.2lf -%5.2lf%8.0f%7d %7d ",
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)),
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)),
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRejects, 
        ptStats->nNumMults);
    if (0.0 != ptStats->fDenom)
    {
        printf("%7.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));
        printf("%6.1f ", ptStats->fSumIoverSig2/ (float)max(1,nTemp));
        printf("%6.2f ", ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups)));
        printf("%6.2f ", ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom));
        printf("%6.3f ", ptStats->fNumer / max(1e-10,ptStats->fDenom));
        printf("%6.3f \n", ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%7s %6s %6s %6s %6s %6s \n", "---","---","---", "---", "---", "---");
    }
    fflush(stdout);

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeReslnRow();
#endif

    if(bSummary == true)
    {
        printf("\nI/sig unavg is the mean I/sig for the unaveraged reflections in the input file.\n");
        printf("I/sig avg   is the mean I/sig for the unique reflections in the output file.\n");
        printf(" * When EMul == %.2lf\n",m_fLastEMulUsed);

    ///////// SUMMARY /////////////////
        ptStats--;   // Now points to highest resolution shell

        printf("\n\nSummary of data collection statistics");
        printf("\n-------------------------------------------------------------");
     // printf("\n Wavelength                       %7.5f", m_fWavelength);
        printf("\n Spacegroup                       %s", m_oCrystal.m_poSpacegroup->sGetName().string());
        printf("\n Unit cell dimensions             %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(0), 
               m_oCrystal.fGetCell(1), 
               m_oCrystal.fGetCell(2));
        printf("\n                                  %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(3), 
               m_oCrystal.fGetCell(4), 
               m_oCrystal.fGetCell(5));
        printf("\n Mosaicity                        %.2lf", m_oCrystal.fGetMosaicity());
        printf("\n\n Resolution range                 %.2f - %.2f", 
               1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)), 1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)));
        printf("    (%.2lf - %.2lf)", fResoLow, fResoHigh);
        printf("\n Total number of reflections      %d", 
               ptSCumul->nNumRefs  - ptSCumul->nNumRejects);
        printf("\n Number of unique reflections     %d", 
               ptSCumul->nNumMults + ptSCumul->nNumSingles);

        printf("\n Average redundancy               ");
        if (0 == (ptSCumul->nNumMults  + ptSCumul->nNumSingles))
          {
            printf("---             ");
          }
        else if (10.0 <= ((float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                          / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles))))
          {
            printf("%.2f           ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        else
          {
            printf("%.2f            ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        if (0 == (ptStats->nNumMults  + ptStats->nNumSingles))
          {
            printf("(---)"); 
          }
        else
          {
            printf("(%.2f)",
                   (float)(ptStats->nNumRefs - ptStats->nNumRejects) 
                   / (float)max(1,(ptStats->nNumMults  + ptStats->nNumSingles)));
          }
        printf("\n %% completeness                   %.1f            (%.1f)",
               min(100.0,100.0 * (float)m_nNumUniqueFoundCumul/(float)max(1,m_nNumUniqueCalcCumul)),
               min(100.0, 100.0 * (float)m_nNumUniqueFoundShell/(float)max(1,m_nNumUniqueCalcShell)));
        printf("\n Rmerge                           ");
        if (0.0 != ptSCumul->fDenom)
          printf("%.3f           ", ptSCumul->fNumer / ptSCumul->fDenom);
        else
          printf("---           ");
        if (0.0 != ptStats->fDenom)
          printf("(%.3f)", ptStats->fNumer / ptStats->fDenom);
        else
          printf("(---)");
        printf("\n Reduced ChiSquared               ");
        if (0 != (ptSCumul->nSumChiContrib - ptSCumul->nSumChiGroups))
          printf("%.2f            ",
                 ptSCumul->fSumChiSq / ((float) ptSCumul->nSumChiContrib - (float)ptSCumul->nSumChiGroups));
        else
          printf("---            ");
        if (0 !=  (ptStats->nSumChiContrib - ptStats->nSumChiGroups))
          {
            printf("(%.2f)", 
               ptStats->fSumChiSq / ((float) ptStats->nSumChiContrib - (float)ptStats->nSumChiGroups));
          }
        else
          printf("(---)");
        
        nTemp  = ptSCumul->nNumSingles  + ptSCumul->nNumUsed;
        if (0 != nTemp)
          {
            if (10.0 <= (ptSCumul->fSumIoverSig2 / (float) nTemp) )
              printf("\n Output <I/sigI>                  %.1f ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
            else
              printf("\n Output <I/sigI>                  %.1f  ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
          }
        nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
        if (0 != nTemp)
          {
            printf(                                      "           (%.1f)", 
                   ptStats->fSumIoverSig2 / (float) nTemp);
          }
        printf("\n-------------------------------------------------------------");
        printf("\n  Note: Values in () are for the last resolution shell.\n\n");
    }
///////////////////////////////////////////// SUMMARY /////////////////

    return 0;

};


//shijie yao 04:02:2008
int Cscaleaverage::nTableRrimVsBatch()
{
    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    int       nTemp;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";


    printf("\n\nIn the tables below the redundancy-independent merging\n");
    printf("R factor Rmeas (or Rrim) is defined as:\n");
    printf("    Rmeas = Sum [N/(N-1)]^(1/2) Sum |Ihi - <Ih>| / Sum Sum <Ih>\n");
    printf("             h                   i                  h   i\n");
    printf("    where Ihi is the ith used observation for unique hkl h,\n");
    printf("    and  <Ih> is the mean intensity for unique hkl h.\n\n");
    printf("  Num_Obs = Num_Overlaps + Num_Rejects + Num_Singles\n"
        "  Num_Mults counts only sym. related groups with at least\n"
        "  2 reflections.\n\n"
        "  Exclusions are included in statistics.\n\n"
        "  |ChiSq| and Rmeas are taken over all non-rejected\n"
        "  reflections having a redundancy of at least 2.\n\n"
        "  <I/Sig> and Average are taken over all non-rejected\n"
        "  reflections.\n"
        );

    // List batch statistics

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeBatchRow();
#endif

    printf("\n\nRmeas vs Batch\n");
    printf(pcLine);
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
        "Batch",  "Average", "Num", " Num", " Num", "    Num", "I/sig",
        "Rducd", "Rmeas", "Rmeas");
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
        " name", " counts", "obs", "rejs", "ovlps", "single", "unavg",
        "ChiSq", " batch", " cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsBatch;
    ptSCumul = &m_ptStatsBatch[m_nNumBatches];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    for (nx = 0; nx < m_nNumBatches; nx++, ptStats++)
    {
        (*ptSCumul) += (*ptStats);
        nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
        printf("%14s %7.0f %7d %5d %7d %7d ",
            m_ptBatchInfo[nx].m_sBatchName.string(),
            ptStats->fSum / max(1,nTemp),
            ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
            ptStats->nNumSingles);
        if (0.0 != ptStats->fDenom)
        {
            printf("%6.1f %6.2f %6.3f ",
                ptStats->fSumIoverSig / (float) max(1,nTemp),

                ptStats->fSumChiSq
                / max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups)),
                ptStats->fNumer_rim / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeBatchRow( 
                m_ptBatchInfo[nx].m_sBatchName,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
                ptStats->nNumSingles,
                (Cstring)(ptStats->fSumIoverSig / (float) max(1,nTemp)),
                (Cstring)(ptStats->fSumChiSq / max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups))),
                (Cstring)(ptStats->fNumer_rim / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer_rim/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else {
            printf("%6s %6s %6s ", "---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeBatchRow(
                m_ptBatchInfo[nx].m_sBatchName,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
                ptStats->nNumSingles,
                (Cstring)"---", (Cstring)"---", (Cstring)"---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
            printf("%6.3f\n", ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));
        else
            printf("%6s\n", "---");
    }
#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeBatchRow();
#endif
    printf(pcLine);
    nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("%14s %7.0f %7d %5d %7d %7d ",
        "All batches",
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
        ptStats->nNumSingles);
    if (0.0 != ptStats->fDenom)
    {
        printf("%6.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));

        printf("%6.2f ", ptStats->fSumChiSq/max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups)));
        printf("%6.3f %6.3f\n", ptStats->fNumer_rim / max(1e-10,ptStats->fDenom),
            ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%6s %6s %6s %6s\n ", "---", "---", "---", "---");
    }

    return 0;
}


int Cscaleaverage::nTableRrimVsIntensity()
{
    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    float     fLow;
    int       nTemp;
    char      cTemp1, cTemp2;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";

    // List intensity statistics

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeISigRow();
#endif

    printf("\n\nRmeas vs Intensity/SigmaI\n");
    printf(pcLine);
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "Int/sigmaI ",  "Average", "Num", " Num", "   Num"," I/sig", " I/sig",
        "Rducd", "Rmeas", "Rmeas");
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "   range   ", " counts", "obs", "rejs", "mults", " unavg","   avg",
        "ChiSq", " shell", " cumul");

    printf(pcLine);
    ptStats  = &m_ptStatsInt[m_nNumBinsInt-1];
    ptSCumul = &m_ptStatsInt[m_nNumBinsInt];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeInt[1];
    pfSlope  = &m_fSlopeInt;
    for (nx = m_nNumBinsInt-1; nx >= 0; nx--, ptStats--)
    {
        (*ptSCumul) += (*ptStats);
        nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
        cTemp1 = ' ';
        cTemp2 = ' ';
        if ((m_nNumBinsInt-1) == nx) cTemp1 = '>';
        if (0 == nx)                 cTemp2 = '<';
        if ( (0.0 != ptStats->fDenom) && (0 != ptStats->nNumUsed) )
        {
            printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d %7.1f %6.1f %6.2f %6.3f ",
                cTemp2, fLow - *pfSlope, cTemp1, fLow,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                ptStats->fSumIoverSig / (float)max(1,nTemp),
                ptStats->fSumIoverSig2/ (float)max(1,nTemp),
                ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups)),
                ptStats->fNumer_rim / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeISigRow( 
                fLow - *pfSlope, fLow,
                (Cstring)(ptStats->fSum / max(1,nTemp)),
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring)(ptStats->fSumIoverSig / (float)max(1,nTemp)),
                (Cstring)(ptStats->fSumIoverSig2/ (float)max(1,nTemp)),
                (Cstring)(ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups))),
                (Cstring)(ptStats->fNumer_rim / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer_rim/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else
        {
            printf("   %c%2.0f -  %c%2.0f%8s%8d%7d %7d %7s %6s %6s %6s ",
                cTemp2, fLow - *pfSlope, cTemp1, fLow,
                "---",
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                "---","---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeISigRow(
                fLow - *pfSlope, fLow,
                (Cstring)"---",
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                (Cstring)"---", (Cstring)"---", (Cstring)"---", (Cstring)"---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer_rim/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
        {
            printf("%6.3f\n", ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));
        }
        else
        {
            printf("%6s\n", "---");
        }
        fLow = fLow - *pfSlope;
    }

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeISigRow();
#endif

    cTemp1 = '>';
    cTemp2 = '<';
    ptStats = ptSCumul;
    printf(pcLine);
    nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d ",
        cTemp2, m_a2fRangeInt[0], cTemp1, m_a2fRangeInt[1],
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRefs, 
        ptStats->nNumRejects, 
        ptStats->nNumMults);
    if (0.0 != ptStats->fDenom)
    {

        printf("%7.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));
        printf("%6.1f ", ptStats->fSumIoverSig2/ (float)max(1,nTemp));
        printf("%6.2f ", ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups)));
        printf("%6.3f %6.3f\n", ptStats->fNumer_rim / max(1e-10,ptStats->fDenom),
            ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%7s %6s %6s %6s %6s\n", "---", "---", "---", "---", "---");
    }

    return 0;
};


int Cscaleaverage::nTableRrimVsResolution(bool bSummary) {

    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    double    fLow;
    double    fResoLow, fResoHigh;
    int       nTemp;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeReslnRow();
#endif

    // List resolution statistics

    printf("\nRmeas vs Resolution\n");
    printf(pcLine);
    printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
        "Resolution  ", "Average",  " Num", "   Num", " I/sig", " I/sig",
        "Rducd", "Model","Rmeas", "Rmeas");
    printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
        "   range    ", " counts", "rejs", "mults", " unavg","   avg",
        "ChiSq", "Eadd*"," shell", " cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsReso;
    ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeReso[0];
    pfSlope  = &m_fSlopeReso;
    for (nx = 0; nx < m_nNumBinsReso; nx++, ptStats++)
    {
        fResoLow              = 1.0 / max(1e-10,pow(fLow, 0.333333));
        fResoHigh             = 1.0 / max(1e-10,pow((fLow + *pfSlope), 0.333333));
        (*ptSCumul) += (*ptStats);
        nTemp = ptStats->nNumSingles + ptStats->nNumUsed;
        if ( (0.0 != ptStats->fDenom) && (0 != ptStats->nNumUsed) )
        {
            printf("%6.2lf -%5.2lf%8.0f%7d %7d %7.1f %6.1f %6.2f %6.2f %6.3f ",
                min(fResoLow, 999.99), fResoHigh,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                ptStats->fSumIoverSig / (float)max(1,nTemp),
                ptStats->fSumIoverSig2/ (float)max(1,nTemp),
                ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups)),
                ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom),
                ptStats->fNumer_rim / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeReslnRow( fResoLow, fResoHigh, (Cstring) (ptStats->fSum / max(1,nTemp)),
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring) (ptStats->fSumIoverSig / (float)max(1,nTemp)),
                (Cstring) (ptStats->fSumIoverSig2/ (float)max(1,nTemp)),
                (Cstring) (ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups))),
                (Cstring) (ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom)),
                (Cstring) (ptStats->fNumer / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else
        {
            printf("%6.2f -%5.2f%8s%7d %7d %7s %6s %6s %6s %6s ",
                min(fResoLow,999.99), fResoHigh,
                "---",
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                "---","---","---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeReslnRow( fResoLow, fResoHigh, (Cstring)"---",
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring) "---",(Cstring) "---", (Cstring) "---", (Cstring) "---",(Cstring) "---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
            printf("%6.3f \n", ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));
        else
            printf("%6s \n", "---");
        fLow = fLow + *pfSlope;
    }
    printf(pcLine);
    nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("%6.2lf -%5.2lf%8.0f%7d %7d ",
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)),
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)),
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRejects, 
        ptStats->nNumMults);
    if (0.0 != ptStats->fDenom)
    {
        printf("%7.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));
        printf("%6.1f ", ptStats->fSumIoverSig2/ (float)max(1,nTemp));
        printf("%6.2f ", ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups)));
        printf("%6.2f ", ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom));
        printf("%6.3f ", ptStats->fNumer_rim / max(1e-10,ptStats->fDenom));
        printf("%6.3f \n", ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%7s %6s %6s %6s %6s %6s \n", "---","---","---", "---", "---", "---");
    }
    fflush(stdout);

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeReslnRow();
#endif

    if(bSummary == true)
    {
        printf("\nI/sig unavg is the mean I/sig for the unaveraged reflections in the input file.\n");
        printf("I/sig avg   is the mean I/sig for the unique reflections in the output file.\n");
        printf(" * When EMul == %.2lf\n",m_fLastEMulUsed);

    ///////// SUMMARY /////////////////
        ptStats--;   // Now points to highest resolution shell

        printf("\n\nSummary of data collection statistics");
        printf("\n-------------------------------------------------------------");
     // printf("\n Wavelength                       %7.5f", m_fWavelength);
        printf("\n Spacegroup                       %s", m_oCrystal.m_poSpacegroup->sGetName().string());
        printf("\n Unit cell dimensions             %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(0), 
               m_oCrystal.fGetCell(1), 
               m_oCrystal.fGetCell(2));
        printf("\n                                  %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(3), 
               m_oCrystal.fGetCell(4), 
               m_oCrystal.fGetCell(5));
        printf("\n Mosaicity                        %.2lf", m_oCrystal.fGetMosaicity());
        printf("\n\n Resolution range                 %.2f - %.2f", 
               1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)), 1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)));
        printf("    (%.2lf - %.2lf)", fResoLow, fResoHigh);
        printf("\n Total number of reflections      %d", 
               ptSCumul->nNumRefs  - ptSCumul->nNumRejects);
        printf("\n Number of unique reflections     %d", 
               ptSCumul->nNumMults + ptSCumul->nNumSingles);

        printf("\n Average redundancy               ");
        if (0 == (ptSCumul->nNumMults  + ptSCumul->nNumSingles))
          {
            printf("---             ");
          }
        else if (10.0 <= ((float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                          / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles))))
          {
            printf("%.2f           ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        else
          {
            printf("%.2f            ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        if (0 == (ptStats->nNumMults  + ptStats->nNumSingles))
          {
            printf("(---)"); 
          }
        else
          {
            printf("(%.2f)",
                   (float)(ptStats->nNumRefs - ptStats->nNumRejects) 
                   / (float)max(1,(ptStats->nNumMults  + ptStats->nNumSingles)));
          }
        printf("\n %% completeness                   %.1f            (%.1f)",
               min(100.0,100.0 * (float)m_nNumUniqueFoundCumul/(float)max(1,m_nNumUniqueCalcCumul)),
               min(100.0, 100.0 * (float)m_nNumUniqueFoundShell/(float)max(1,m_nNumUniqueCalcShell)));
        printf("\n Rmeas                            ");
        if (0.0 != ptSCumul->fDenom)
          printf("%.3f           ", ptSCumul->fNumer_rim / ptSCumul->fDenom);
        else
          printf("---           ");
        if (0.0 != ptStats->fDenom)
          printf("(%.3f)", ptStats->fNumer_rim / ptStats->fDenom);
        else
          printf("(---)");
        printf("\n Reduced ChiSquared               ");
        if (0 != (ptSCumul->nSumChiContrib - ptSCumul->nSumChiGroups))
          printf("%.2f            ",
                 ptSCumul->fSumChiSq / ((float) ptSCumul->nSumChiContrib - (float)ptSCumul->nSumChiGroups));
        else
          printf("---            ");
        if (0 !=  (ptStats->nSumChiContrib - ptStats->nSumChiGroups))
          {
            printf("(%.2f)", 
               ptStats->fSumChiSq / ((float) ptStats->nSumChiContrib - (float)ptStats->nSumChiGroups));
          }
        else
          printf("(---)");
        
        nTemp  = ptSCumul->nNumSingles  + ptSCumul->nNumUsed;
        if (0 != nTemp)
          {
            if (10.0 <= (ptSCumul->fSumIoverSig2 / (float) nTemp) )
              printf("\n Output <I/sigI>                  %.1f ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
            else
              printf("\n Output <I/sigI>                  %.1f  ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
          }
        nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
        if (0 != nTemp)
          {
            printf(                                      "           (%.1f)", 
                   ptStats->fSumIoverSig2 / (float) nTemp);
          }
        printf("\n-------------------------------------------------------------");
        printf("\n  Note: Values in () are for the last resolution shell.\n\n");

    ///////////////////////////////////////////// SUMMARY /////////////////
    }
    return 0;

};



int Cscaleaverage::nTableRpimVsBatch()
{
  //+-JWP 2008-06-13  WARNING WARNING: Rpim variables have been usurped. 
  //                  This routine should NOT be used!

    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    int       nTemp;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";


    printf("\n\nIn the tables below the precision-indicating merging\n");
    printf("WARNING WARNING WARNING the follow table should NOT be printed.\n");
    printf("WARNING WARNING WARNING the follow table has invalid numbers.\n");
    printf("R factor Rpim is defined as:\n");
    printf("    Rpim = Sum [1/(N-1)]^(1/2) Sum |Ihi - <Ih>| / Sum Sum <Ih>\n");
    printf("            h                   i                  h   i\n");
    printf("    where Ihi is the ith used observation for unique hkl h,\n");
    printf("    and  <Ih> is the mean intensity for unique hkl h.\n\n");
    printf("  Num_Obs = Num_Overlaps + Num_Rejects + Num_Singles\n"
        "  Num_Mults counts only sym. related groups with at least\n"
        "  2 reflections.\n\n"
        "  Exclusions are included in statistics.\n\n"
        "  |ChiSq| and Rpim are taken over all non-rejected\n"
        "  reflections having a redundancy of at least 2.\n\n"
        "  <I/Sig> and Average are taken over all non-rejected\n"
        "  reflections.\n"
        );

    // List batch statistics

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeBatchRow();
#endif

    printf("\n\nRpim vs Batch\n");
    printf(pcLine);
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
        "Batch",  "Average", "Num", " Num", " Num", "    Num", "I/sig",
        "Rducd", "Rpim", "Rpim");
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
        " name", " counts", "obs", "rejs", "ovlps", "single", "unavg",
        "ChiSq", " batch", " cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsBatch;
    ptSCumul = &m_ptStatsBatch[m_nNumBatches];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    for (nx = 0; nx < m_nNumBatches; nx++, ptStats++)
    {
        (*ptSCumul) += (*ptStats);
        nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
        printf("%14s %7.0f %7d %5d %7d %7d ",
            m_ptBatchInfo[nx].m_sBatchName.string(),
            ptStats->fSum / max(1,nTemp),
            ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
            ptStats->nNumSingles);
        if (0.0 != ptStats->fDenom)
        {
            printf("%6.1f %6.2f %6.3f ",
                ptStats->fSumIoverSig / (float) max(1,nTemp),

                ptStats->fSumChiSq
                / max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups)),
                ptStats->fNumer_pim / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeBatchRow( 
                m_ptBatchInfo[nx].m_sBatchName,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
                ptStats->nNumSingles,
                (Cstring)(ptStats->fSumIoverSig / (float) max(1,nTemp)),
                (Cstring)(ptStats->fSumChiSq / max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups))),
                (Cstring)(ptStats->fNumer_pim / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer_pim/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else {
            printf("%6s %6s %6s ", "---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeBatchRow(
                m_ptBatchInfo[nx].m_sBatchName,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
                ptStats->nNumSingles,
                (Cstring)"---", (Cstring)"---", (Cstring)"---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
            printf("%6.3f\n", ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenom));
        else
            printf("%6s\n", "---");
    }
#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeBatchRow();
#endif
    printf(pcLine);
    nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("%14s %7.0f %7d %5d %7d %7d ",
        "All batches",
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
        ptStats->nNumSingles);
    if (0.0 != ptStats->fDenom)
    {
        printf("%6.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));

        printf("%6.2f ", ptStats->fSumChiSq/max(1,(ptStats->nSumChiContrib - ptStats->nSumChiGroups)));
        printf("%6.3f %6.3f\n", ptStats->fNumer_pim / max(1e-10,ptStats->fDenom),
            ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%6s %6s %6s %6s\n ", "---", "---", "---", "---");
    }

    return 0;
}




int Cscaleaverage::nTableRpimVsIntensity()
{
    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    float     fLow;
    int       nTemp;
    char      cTemp1, cTemp2;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";

    // List intensity statistics

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeISigRow();
#endif

    printf("\n\nRpim vs Intensity/SigmaI\n");
    printf(pcLine);
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "Int/sigmaI ",  "Average", "Num", " Num", "   Num"," I/sig", " I/sig",
        "Rducd", "Rpim", "Rpim");
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "   range   ", " counts", "obs", "rejs", "mults", " unavg","   avg",
        "ChiSq", " shell", " cumul");

    printf(pcLine);
    ptStats  = &m_ptStatsInt[m_nNumBinsInt-1];
    ptSCumul = &m_ptStatsInt[m_nNumBinsInt];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeInt[1];
    pfSlope  = &m_fSlopeInt;
    for (nx = m_nNumBinsInt-1; nx >= 0; nx--, ptStats--)
    {
        (*ptSCumul) += (*ptStats);
        nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
        cTemp1 = ' ';
        cTemp2 = ' ';
        if ((m_nNumBinsInt-1) == nx) cTemp1 = '>';
        if (0 == nx)                 cTemp2 = '<';
        if ( (0.0 != ptStats->fDenom) && (0 != ptStats->nNumUsed) )
        {
            printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d %7.1f %6.1f %6.2f %6.3f ",
                cTemp2, fLow - *pfSlope, cTemp1, fLow,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                ptStats->fSumIoverSig / (float)max(1,nTemp),
                ptStats->fSumIoverSig2/ (float)max(1,nTemp),
                ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups)),
                ptStats->fNumer_pim / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeISigRow( 
                fLow - *pfSlope, fLow,
                (Cstring)(ptStats->fSum / max(1,nTemp)),
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring)(ptStats->fSumIoverSig / (float)max(1,nTemp)),
                (Cstring)(ptStats->fSumIoverSig2/ (float)max(1,nTemp)),
                (Cstring)(ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups))),
                (Cstring)(ptStats->fNumer_pim / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer_pim/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else
        {
            printf("   %c%2.0f -  %c%2.0f%8s%8d%7d %7d %7s %6s %6s %6s ",
                cTemp2, fLow - *pfSlope, cTemp1, fLow,
                "---",
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                "---","---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeISigRow(
                fLow - *pfSlope, fLow,
                (Cstring)"---",
                ptStats->nNumRefs, 
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                (Cstring)"---", (Cstring)"---", (Cstring)"---", (Cstring)"---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer_pim/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
        {
            printf("%6.3f\n", ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenom));
        }
        else
        {
            printf("%6s\n", "---");
        }
        fLow = fLow - *pfSlope;
    }

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeISigRow();
#endif

    cTemp1 = '>';
    cTemp2 = '<';
    ptStats = ptSCumul;
    printf(pcLine);
    nTemp       = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d ",
        cTemp2, m_a2fRangeInt[0], cTemp1, m_a2fRangeInt[1],
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRefs, 
        ptStats->nNumRejects, 
        ptStats->nNumMults);
    if (0.0 != ptStats->fDenom)
    {

        printf("%7.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));
        printf("%6.1f ", ptStats->fSumIoverSig2/ (float)max(1,nTemp));
        printf("%6.2f ", ptStats->fSumChiSq  / max((float) 1.0,((float)ptStats->nSumChiContrib - (float) ptStats->nSumChiGroups)));
        printf("%6.3f %6.3f\n", ptStats->fNumer_pim / max(1e-10,ptStats->fDenom),
            ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%7s %6s %6s %6s %6s\n", "---", "---", "---", "---", "---");
    }

    return 0;
};


int Cscaleaverage::nTableRpimVsResolution(bool bSummary) {

    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    double    fLow;
    double    fResoLow, fResoHigh;
    int       nTemp;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeReslnRow();
#endif

    // List resolution statistics

    printf("\nRpim vs Resolution\n");
    printf(pcLine);
    printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
        "Resolution  ", "Average",  " Num", "   Num", " I/sig", " I/sig",
        "Rducd", "Model","Rpim", "Rpim");
    printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
        "   range    ", " counts", "rejs", "mults", " unavg","   avg",
        "ChiSq", "Eadd*"," shell", " cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsReso;
    ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeReso[0];
    pfSlope  = &m_fSlopeReso;
    for (nx = 0; nx < m_nNumBinsReso; nx++, ptStats++)
    {
        fResoLow              = 1.0 / max(1e-10,pow(fLow, 0.333333));
        fResoHigh             = 1.0 / max(1e-10,pow((fLow + *pfSlope), 0.333333));
        (*ptSCumul) += (*ptStats);
        nTemp = ptStats->nNumSingles + ptStats->nNumUsed;
        if ( (0.0 != ptStats->fDenom) && (0 != ptStats->nNumUsed) )
        {
            printf("%6.2lf -%5.2lf%8.0f%7d %7d %7.1f %6.1f %6.2f %6.2f %6.3f ",
                min(fResoLow, 999.99), fResoHigh,
                ptStats->fSum / max(1,nTemp),
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                ptStats->fSumIoverSig / (float)max(1,nTemp),
                ptStats->fSumIoverSig2/ (float)max(1,nTemp),
                ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups)),
                ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom),
                ptStats->fNumer_pim / max(1e-10,ptStats->fDenom));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeReslnRow( fResoLow, fResoHigh, (Cstring) (ptStats->fSum / max(1,nTemp)),
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring) (ptStats->fSumIoverSig / (float)max(1,nTemp)),
                (Cstring) (ptStats->fSumIoverSig2/ (float)max(1,nTemp)),
                (Cstring) (ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups))),
                (Cstring) (ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom)),
                (Cstring) (ptStats->fNumer / max(1e-10,ptStats->fDenom)),
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        else
        {
            printf("%6.2f -%5.2f%8s%7d %7d %7s %6s %6s %6s %6s ",
                min(fResoLow,999.99), fResoHigh,
                "---",
                ptStats->nNumRejects, 
                ptStats->nNumMults, 
                "---","---","---", "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddRmergeReslnRow( fResoLow, fResoHigh, (Cstring)"---",
                ptStats->nNumRejects, 
                ptStats->nNumMults,
                (Cstring) "---",(Cstring) "---", (Cstring) "---", (Cstring) "---",(Cstring) "---",
                ((ptSCumul->fDenom!=0.0) ? (Cstring)(ptSCumul->fNumer/max(1e-10,ptSCumul->fDenom)) : (Cstring)"---") );
#endif
        }
        if (0.0 != ptSCumul->fDenom)
            printf("%6.3f \n", ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenom));
        else
            printf("%6s \n", "---");
        fLow = fLow + *pfSlope;
    }
    printf(pcLine);
    nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("%6.2lf -%5.2lf%8.0f%7d %7d ",
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)),
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)),
        ptStats->fSum / max(1,nTemp),
        ptStats->nNumRejects, 
        ptStats->nNumMults);
    if (0.0 != ptStats->fDenom)
    {
        printf("%7.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));
        printf("%6.1f ", ptStats->fSumIoverSig2/ (float)max(1,nTemp));
        printf("%6.2f ", ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups)));
        printf("%6.2f ", ptStats->fEAddNumer/max(1e-10,ptStats->fEAddDenom));
        printf("%6.3f ", ptStats->fNumer_pim / max(1e-10,ptStats->fDenom));
        printf("%6.3f \n", ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenom));
    }
    else
    {
        printf("%7s %6s %6s %6s %6s %6s \n", "---","---","---", "---", "---", "---");
    }
    fflush(stdout);

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeReslnRow();
#endif

    if(bSummary == true)
    {
        printf("\nI/sig unavg is the mean I/sig for the unaveraged reflections in the input file.\n");
        printf("I/sig avg   is the mean I/sig for the unique reflections in the output file.\n");
        printf(" * When EMul == %.2lf\n",m_fLastEMulUsed);

    ///////// SUMMARY /////////////////
        ptStats--;   // Now points to highest resolution shell

        printf("\n\nSummary of data collection statistics");
        printf("\n-------------------------------------------------------------");
     // printf("\n Wavelength                       %7.5f", m_fWavelength);
        printf("\n Spacegroup                       %s", m_oCrystal.m_poSpacegroup->sGetName().string());
        printf("\n Unit cell dimensions             %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(0), 
               m_oCrystal.fGetCell(1), 
               m_oCrystal.fGetCell(2));
        printf("\n                                  %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(3), 
               m_oCrystal.fGetCell(4), 
               m_oCrystal.fGetCell(5));
        printf("\n Mosaicity                        %.2lf", m_oCrystal.fGetMosaicity());
        printf("\n\n Resolution range                 %.2f - %.2f", 
               1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)), 1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)));
        printf("    (%.2lf - %.2lf)", fResoLow, fResoHigh);
        printf("\n Total number of reflections      %d", 
               ptSCumul->nNumRefs  - ptSCumul->nNumRejects);
        printf("\n Number of unique reflections     %d", 
               ptSCumul->nNumMults + ptSCumul->nNumSingles);

        printf("\n Average redundancy               ");
        if (0 == (ptSCumul->nNumMults  + ptSCumul->nNumSingles))
          {
            printf("---             ");
          }
        else if (10.0 <= ((float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                          / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles))))
          {
            printf("%.2f           ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        else
          {
            printf("%.2f            ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        if (0 == (ptStats->nNumMults  + ptStats->nNumSingles))
          {
            printf("(---)"); 
          }
        else
          {
            printf("(%.2f)",
                   (float)(ptStats->nNumRefs - ptStats->nNumRejects) 
                   / (float)max(1,(ptStats->nNumMults  + ptStats->nNumSingles)));
          }
        printf("\n %% completeness                   %.1f            (%.1f)",
               min(100.0,100.0 * (float)m_nNumUniqueFoundCumul/(float)max(1,m_nNumUniqueCalcCumul)),
               min(100.0, 100.0 * (float)m_nNumUniqueFoundShell/(float)max(1,m_nNumUniqueCalcShell)));
        printf("\n Rpim                             ");
        if (0.0 != ptSCumul->fDenom)
          printf("%.3f           ", ptSCumul->fNumer_pim / ptSCumul->fDenom);
        else
          printf("---           ");
        if (0.0 != ptStats->fDenom)
          printf("(%.3f)", ptStats->fNumer_pim / ptStats->fDenom);
        else
          printf("(---)");
        printf("\n Reduced ChiSquared               ");
        if (0 != (ptSCumul->nSumChiContrib - ptSCumul->nSumChiGroups))
          printf("%.2f            ",
                 ptSCumul->fSumChiSq / ((float) ptSCumul->nSumChiContrib - (float)ptSCumul->nSumChiGroups));
        else
          printf("---            ");
        if (0 !=  (ptStats->nSumChiContrib - ptStats->nSumChiGroups))
          {
            printf("(%.2f)", 
               ptStats->fSumChiSq / ((float) ptStats->nSumChiContrib - (float)ptStats->nSumChiGroups));
          }
        else
          printf("(---)");
        
        nTemp  = ptSCumul->nNumSingles  + ptSCumul->nNumUsed;
        if (0 != nTemp)
          {
            if (10.0 <= (ptSCumul->fSumIoverSig2 / (float) nTemp) )
              printf("\n Output <I/sigI>                  %.1f ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
            else
              printf("\n Output <I/sigI>                  %.1f  ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
          }
        nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
        if (0 != nTemp)
          {
            printf(                                      "           (%.1f)", 
                   ptStats->fSumIoverSig2 / (float) nTemp);
          }
        printf("\n-------------------------------------------------------------");
        printf("\n  Note: Values in () are for the last resolution shell.\n\n");

    ///////////////////////////////////////////// SUMMARY /////////////////
    }
    return 0;

};


int Cscaleaverage::nTable3RVsResolution(bool bSummary) {

    int nx;

    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    double    fLow;
    double    fResoLow, fResoHigh;
    int       nTemp;
    char     *pcLine =
        "--------------------------------------------------------------------------------\n";

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetRmergeReslnRow();
#endif

    // List resolution statistics

    printf("\nMerging R factors vs Resolution\n");
    printf(pcLine);
    //printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
    //    "Resolution  ", "Average",  " Num", "   Num", " I/sig", " I/sig",
    //    "Rducd", "Model","Rmerge", "Rmerge");
    //printf("%13s %7s %6s %7s %7s %6s %6s %6s %6s %6s \n",
    //    "   range    ", " counts", "rejs", "mults", " unavg","   avg",
    //    "ChiSq", "Eadd*"," shell", " cumul");

    //+JWP 2008-05-23
    // Change order, put Rpim before the others instead of last

    printf("%12s  %5s %5s %6s %6s %6s %5s %6s %6s %6s %6s\n",
        "Resolution ","I/sig", "I/sig", "Rducd", "RmeasA", "RmeasA", "Rducd", "Rmerge", "Rmerge", 
            "Rmeas",  "Rmeas");
    //-JWP 2008-05-23

    printf("%12s  %5s %5s %5s %6s %6s %5s %6s %6s %6s %6s\n",
        "   range  ", "unavg","  avg", "ChiSqA", "shell", "cumul", "ChiSq", 
        "shell", "cumul", "shell", "cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsReso;
    ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeReso[0];
    pfSlope  = &m_fSlopeReso;
    for (nx = 0; nx < m_nNumBinsReso; nx++, ptStats++)
    {
        fResoLow              = 1.0 / max(1e-10,pow(fLow, 0.333333));
        fResoHigh             = 1.0 / max(1e-10,pow((fLow + *pfSlope), 0.333333));
        (*ptSCumul) += (*ptStats);
        nTemp = ptStats->nNumSingles + ptStats->nNumUsed;

        //tmp values

        char acResMin[20], acResMax[20], acIsigUnavg[20], acIsigAvg[20], acRedChiSq[20], acRedChiSqA[20];
        char acRmergeS[20], acRmergeC[20], acRmeasS[20], acRmeasC[20], acRpimS[20], acRpimC[20];

        sprintf(acResMin,    "%6.2lf", min(fResoLow, 999.99));
        sprintf(acResMax,    "%5.2lf", fResoHigh);

        sprintf(acRedChiSq,  "%5.2lf", ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups))); 
	/***
	cout << "chisqa:, sum, contrib, groupsA, groups: " << ptStats->fSumChiSqA << ", " 
             <<  ptStats->nSumChiAContrib  << ", " << ptStats->nSumChiAGroups 
             << ", " << ptStats->nSumChiGroups << endl;
	***/

	sprintf(acRedChiSqA,  "%5.2lf", ptStats->fSumChiSqA / max((float) 1.0,((float) ptStats->nSumChiAContrib  - ((float)ptStats->nSumChiAGroups*0.5)))); 
	
        if ( (0.0 != ptStats->fDenom) && (0 != ptStats->nNumUsed) )
        {
            sprintf(acIsigUnavg, "%6.1f", ptStats->fSumIoverSig / (float)max(1,nTemp)); 
            sprintf(acIsigAvg,   "%6.1f", ptStats->fSumIoverSig2/ (float)max(1,nTemp)); 
	    // Why is this done a second time?

	    //////  sprintf(acRedChiSq,  "%6.2f", ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups))); 

	    if (0.0 != ptStats->fDenomAnom)
	      sprintf(acRpimS,     "%6.3f", ptStats->fNumer_pim / max(1e-10,ptStats->fDenomAnom));  //+-JWP 2008-05-23 order changed
	    else
	      sprintf(acRpimS,   "%6s", "---"); //+-JWP 2008-05-23 order changed

            sprintf(acRmergeS,   "%6.3f", ptStats->fNumer / max(1e-10,ptStats->fDenom));
            sprintf(acRmeasS,    "%6.3f", ptStats->fNumer_rim / max(1e-10,ptStats->fDenom)); 

        }
        else
        {
            sprintf(acIsigUnavg, "%6s", "---");
            sprintf(acIsigAvg,   "%6s", "---"); 
            sprintf(acRpimS,     "%6s", "---"); //+-JWP 2008-05-23 order changed
            sprintf(acRmergeS,   "%6s", "---");
            sprintf(acRmeasS,    "%6s", "---"); 
        }

        if (0.0 != ptSCumul->fDenom)
        {
            sprintf(acRpimC,     "%6.3f",  ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenomAnom)); //+-JWP 2008-05-23 order changed
            sprintf(acRmergeC,   "%6.3f",  ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom)); 
            sprintf(acRmeasC,    "%6.3f",  ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom)); 

        }
        else
        {
            sprintf(acRpimC,     "%6s",  "---");//+-JWP 2008-05-23 order changed
            sprintf(acRmergeC,   "%6s",  "---"); 
            sprintf(acRmeasC,    "%6s",  "---"); 

        }

        printf("%6s -%5s%5s%6s %6s %5s %6s%6s %6s %6s %6s %6s\n",
               acResMin, acResMax, acIsigUnavg, acIsigAvg, 
	       acRedChiSqA, acRpimS, acRpimC, 
	       acRedChiSq, acRmergeS, acRmergeC, acRmeasS, acRmeasC);

        fLow = fLow + *pfSlope;
    }
    printf(pcLine);
    nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
    printf("%6.2lf -%5.2lf",
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)),
        1.0 / (float)max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)));

    if (0.0 != ptStats->fDenom)
    {
        printf("%6.1f ", ptStats->fSumIoverSig / (float)max(1,nTemp));
        printf("%5.1f ", ptStats->fSumIoverSig2/ (float)max(1,nTemp));

	// ChiSqA:
        printf("%6.2f ", ptStats->fSumChiSqA / max((float) 1.0,((float) ptStats->nSumChiAContrib  - ((float)ptStats->nSumChiAGroups*0.5))));
        

        printf("%6.3f ", ptStats->fNumer_pim / max(1e-10,ptStats->fDenomAnom));
        printf("%6.3f ", ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenomAnom));

        printf("%5.2f ", ptStats->fSumChiSq / max((float) 1.0,((float) ptStats->nSumChiContrib  - (float)ptStats->nSumChiGroups)));
        printf("%6.3f ", ptStats->fNumer / max(1e-10,ptStats->fDenom));
        printf("%6.3f ", ptSCumul->fNumer / max(1e-10,ptSCumul->fDenom));

        printf("%6.3f ", ptStats->fNumer_rim / max(1e-10,ptStats->fDenom));
        printf("%6.3f\n", ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));

    }
    else
    {
        printf("%5s%5s%6s %6s %6s %6s %6s %6s %6s %6s\n", "---","---","---", "---", "---", "---", "---", "---", "---", "---");
    }
    fflush(stdout);

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteRmergeReslnRow();
#endif

    if (bSummary)
    {
        printf("\nI/sig unavg is the mean I/sig for the unaveraged reflections in the input file.\n");
        printf("I/sig avg   is the mean I/sig for the unique reflections in the output file.\n");
	printf("Rduced ChiSqA and RmeasA treat I+ and I- reflns separately.\n");
	if (m_bScaleAnom)
	  printf("Note: Since -scaleanom was used for this scaling, RmeasA==Rmeas & ChiSqA==ChiSq.\n");
        printf(" * When EMul == %.2lf\n",m_fLastEMulUsed);


    ///////// SUMMARY /////////////////
        ptStats--;   // Now points to highest resolution shell

        printf("\n\nSummary of data collection statistics");
        printf("\n-------------------------------------------------------------");
     // printf("\n Wavelength                       %7.5f", m_fWavelength);
        printf("\n Spacegroup                       %s", m_oCrystal.m_poSpacegroup->sGetName().string());
        printf("\n Unit cell dimensions             %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(0), 
               m_oCrystal.fGetCell(1), 
               m_oCrystal.fGetCell(2));
        printf("\n                                  %.2lf   %.2lf   %.2lf", 
               m_oCrystal.fGetCell(3), 
               m_oCrystal.fGetCell(4), 
               m_oCrystal.fGetCell(5));
        printf("\n Mosaicity                        %.2lf", m_oCrystal.fGetMosaicity());
        printf("\n\n Resolution range                 %.2f - %.2f", 
               1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)), 1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)));
        printf("    (%.2lf - %.2lf)", fResoLow, fResoHigh);
        printf("\n Total number of reflections      %d", 
               ptSCumul->nNumRefs  - ptSCumul->nNumRejects);
        printf("\n Number of unique reflections     %d", 
               ptSCumul->nNumMults + ptSCumul->nNumSingles);

        printf("\n Average redundancy               ");
        if (0 == (ptSCumul->nNumMults  + ptSCumul->nNumSingles))
          {
            printf("---             ");
          }
        else if (10.0 <= ((float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                          / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles))))
          {
            printf("%.2f           ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        else
          {
            printf("%.2f            ", (float)(ptSCumul->nNumRefs - ptSCumul->nNumRejects) 
                   / (float)max(1,(ptSCumul->nNumMults  + ptSCumul->nNumSingles)));
          }
        if (0 == (ptStats->nNumMults  + ptStats->nNumSingles))
          {
            printf("(---)"); 
          }
        else
          {
            printf("(%.2f)",
                   (float)(ptStats->nNumRefs - ptStats->nNumRejects) 
                   / (float)max(1,(ptStats->nNumMults  + ptStats->nNumSingles)));
          }
        printf("\n %% completeness                   %.1f            (%.1f)",
               min(100.0,100.0 * (float)m_nNumUniqueFoundCumul/(float)max(1,m_nNumUniqueCalcCumul)),
               min(100.0, 100.0 * (float)m_nNumUniqueFoundShell/(float)max(1,m_nNumUniqueCalcShell)));

        printf("\n Rmerge                           ");
        if (0.0 != ptSCumul->fDenom)
          printf("%.3f           ", ptSCumul->fNumer / ptSCumul->fDenom);
        else
          printf("---           ");
        if (0.0 != ptStats->fDenom)
          printf("(%.3f)", ptStats->fNumer / ptStats->fDenom);
        else
          printf("(---)");

        printf("\n Rmeas                            ");
        if (0.0 != ptSCumul->fDenom)
          printf("%.3f           ", ptSCumul->fNumer_rim / max(1e-10,ptSCumul->fDenom));
        else
          printf("---           ");
        if (0.0 != ptStats->fDenom)
          printf("(%.3f)", ptStats->fNumer_rim / max(1e-10, ptStats->fDenom));
        else
          printf("(---)");

        printf("\n RmeasA (I+,I- reflns kept apart) ");
        if (0.0 != ptSCumul->fDenomAnom)
          printf("%.3f           ", ptSCumul->fNumer_pim / max(1e-10,ptSCumul->fDenomAnom));
        else
          printf("---           ");
        if (0.0 != ptStats->fDenomAnom)
          printf("(%.3f)", ptStats->fNumer_pim / max(1e-10, ptStats->fDenomAnom));
        else
          printf("(---)");

        printf("\n Reduced ChiSquared               ");
        if (0 != (ptSCumul->nSumChiContrib - ptSCumul->nSumChiGroups))
          printf("%.2f            ",
                 ptSCumul->fSumChiSq / ((float) ptSCumul->nSumChiContrib - (float)ptSCumul->nSumChiGroups));
        else
          printf("---            ");
        if (0 !=  (ptStats->nSumChiContrib - ptStats->nSumChiGroups))
          {
            printf("(%.2f)", 
               ptStats->fSumChiSq / ((float) ptStats->nSumChiContrib - (float)ptStats->nSumChiGroups));
          }
        else
          printf("(---)");
        
        nTemp  = ptSCumul->nNumSingles  + ptSCumul->nNumUsed;
        if (0 != nTemp)
          {
            if (10.0 <= (ptSCumul->fSumIoverSig2 / (float) nTemp) )
              printf("\n Output <I/sigI>                  %.1f ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
            else
              printf("\n Output <I/sigI>                  %.1f  ", 
                     ptSCumul->fSumIoverSig2 / (float) nTemp);
          }
        nTemp  = ptStats->nNumSingles  + ptStats->nNumUsed;
        if (0 != nTemp)
          {
            printf(                                      "           (%.1f)", 
                   ptStats->fSumIoverSig2 / (float) nTemp);
          }
        printf("\n-------------------------------------------------------------");
        printf("\n  Note: Values in () are for the last resolution shell.\n\n");

    ///////////////////////////////////////////// SUMMARY /////////////////
    }
    return 0;

};





int Cscaleaverage::nTableOverlap(void)
{
  // List matrix of overlaps between scaling batches
  // The matrix may be large, so list out in columns of 10

  int nx,ny,nz;

  int  *pnTemp;
  char *pcDash6 = "------";
  Cstring sTemp;

  {
      // Print the batch overlaps in a file.
      Cimage oBatchImage(m_nNumBatches,m_nNumBatches,eImage_uI2);
      for (nx = 0; nx < m_nNumBatches; nx++) {
          for (ny = 0; ny < m_nNumBatches; ny++) {
             oBatchImage.nSetPixel(nx,ny,(unsigned short int) m_pnNumOverlap[nx * m_nNumBatches + ny]);
          };
      };
      sTemp = "batchoverlaps.img";
      oBatchImage.nWrite(sTemp);
      printf("Writing file '%s' containing batch overlaps.\n",sTemp.string());
  };

  printf("\n\nOverlaps among scaling batches");
  for (nz = 0; nz < m_nNumBatches; nz = nz + 10)
    {
      printf("\n");

      // List batchnames across top

      printf("%6s--", pcDash6);
      for (ny= nz; (ny< m_nNumBatches) && (ny< nz + 10); ny++)
        {
          printf(pcDash6);
        }
      printf("--%6s\n", pcDash6);

      printf("%5s | ", "Batch");
      for (ny= nz; (ny< m_nNumBatches) && (ny< nz + 10); ny++)
        {
          printf("%6s", m_ptBatchInfo[ny].m_sBatchName.string());
        }
      printf ("\n");
      printf("%6s--", pcDash6);
      for (ny= nz; (ny< m_nNumBatches) && (ny< nz + 10); ny++)
        {
          printf(pcDash6);
        }
      printf("--%6s\n", pcDash6);

      for (nx = 0; nx < m_nNumBatches; nx++)
        {
          // List batchnames on left side

//          printf("%5s | ", (const char *) m_psBatchNames[nx]);
          printf("%5s | ", m_ptBatchInfo[nx].m_sBatchName.string());

          pnTemp = &m_pnNumOverlap[nx * m_nNumBatches] + nz;

          for (ny= nz; (ny< m_nNumBatches) && (ny< nz + 10); ny++)
            {
              printf("%6d", *pnTemp++);
            } // end ny(column)

          // List batchnames on right side

          printf(" |%5s\n", m_ptBatchInfo[nx].m_sBatchName.string());
        } // end nx (row)

      // List batchnames across bottom

      printf("%6s--", pcDash6);
      for (ny= nz; (ny< m_nNumBatches) && (ny< nz + 10); ny++)
        {
          printf(pcDash6);
        }
      printf("--%6s\n", pcDash6);

      printf("%5s | ", "Batch");
      for (ny= nz; (ny< m_nNumBatches) && (ny< nz + 10); ny++)
        {
//          printf("%6s", (const char *)m_psBatchNames[ny]);
          printf("%6s", m_ptBatchInfo[ny].m_sBatchName.string());
        }
      printf ("\n");
      printf("%6s--", pcDash6);
      for (ny= nz; (ny< m_nNumBatches) && (ny< nz + 10); ny++)
        {
          printf(pcDash6);
        }
      printf("--%6s\n", pcDash6);
    } // end nz (group of 10)
  printf("\n");
  return 0;
}




int Cscaleaverage::nTableCompleteness(void)
{
    // Calculate and list completeness table

    int i;
    int nNumUniqueCalcCumul;
    int nNumUniqueFoundCumul;
    tagStatsScaleAverage *ptStats;
    tagStatsScaleAverage *ptSCumul;
    float    *pfSlope;
    double   fLow;
    double   fResoLow, fResoHigh;
    char     *pcLine = "--------------------------------------------------------------------------------\n";
    itr<double>     afResoLower;
    itr<double>     afResoUpper;
    itr<int>        anResoCounts;


  // List completeness statistics

//+JWP 2008-05-23
/****
    if (m_bScaleAnom)
    {
        printf("\n\nIn the following completeness and redundancy tables,");
        printf("\nI+ and I- are treated as non-equivalent reflections.\n");
    }
****/
//-JWP 2008-05-23

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vResetCompReslnRow();
#endif

    printf("\n\nCompleteness vs Resolution\n");
    printf(pcLine);
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "Resolution  ", "  Calc", "Num", " Num", " Num", "   Num", "Num",
        " Avg", "%Comp", "%Comp");
    printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
        "   range    ", "unique", "obs", "rejs", "mults", "single",
        "unique", "mult", "shell", "cumul");

    printf(pcLine);
    ptStats  =  m_ptStatsReso;
    ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeReso[0];
    pfSlope  = &m_fSlopeReso;
    int nNumUniqueCalc, nNumUniqueFound;

    nNumUniqueCalcCumul=0;

    afResoLower.clear();
    afResoUpper.clear();
    anResoCounts.clear();
    for (i = 0; i < m_nNumBinsReso; i++)
    {
        fResoLow              = 1.0 / max(1e-10,pow(fLow, 0.333333));
        fResoHigh             = 1.0 / max(1e-10,pow(fLow + *pfSlope, 0.333333));
        afResoLower + fResoLow;
        afResoUpper + fResoHigh;        
        fLow = fLow + *pfSlope;
    };

//+JWP 2008-05-23
    int nScaleAnom = 0;
    nScaleAnom = m_bScaleAnom;
    //if ("" != sGetEnv("DTREK_SCALEANOM"))
    //nScaleAnom = 0;  // Force this to be 0

    m_oCrystal.nCountUnique(m_nNumBinsReso,
			    &afResoLower[0],
			    &afResoUpper[0],
			    &anResoCounts[0],
			    //nScaleAnom);           // Always 0 now
			    m_bScaleAnom);
//-JWP 2008-05-23

    fLow     =  m_a2fRangeReso[0];
    pfSlope  = &m_fSlopeReso;

    for (i = 0; i < m_nNumBinsReso; i++, ptStats++)
    {
        fResoLow              = 1.0 / max(1e-10,pow(fLow, 0.333333));
        fResoHigh             = 1.0 / max(1e-10,pow(fLow + *pfSlope, 0.333333));
        nNumUniqueCalc = anResoCounts[i];

        nNumUniqueFound       = ptStats->nNumMults    + ptStats->nNumSingles;
        (*ptSCumul) += (*ptStats);
        nNumUniqueCalcCumul += nNumUniqueCalc;

        if ( (0 != nNumUniqueCalc) && (0 != ptStats->nNumUsed)
            && (0 != nNumUniqueFound) )
        {
            printf("%6.2lf -%5.2lf%8d%8d%7d%8d%8d%7d%7.2f%7.1f ",
                min(fResoLow, 999.99), fResoHigh,
                nNumUniqueCalc,
                ptStats->nNumRefs, ptStats->nNumRejects,
                ptStats->nNumMults,
                ptStats->nNumSingles,
                nNumUniqueFound,
                (float) (ptStats->nNumUsed + ptStats->nNumSingles)
                / (float) max(1,nNumUniqueFound),
                min(100., 100. * (float) nNumUniqueFound / (float) max(1,nNumUniqueCalc)));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddCompReslnRow( fResoLow, fResoHigh, nNumUniqueCalc,
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumMults, ptStats->nNumSingles,
                nNumUniqueFound, (Cstring) ((float)(ptStats->nNumUsed + ptStats->nNumSingles) / (float) max(1,nNumUniqueFound)),
                (Cstring) min(100., 100. * (float) nNumUniqueFound / (float) max(1,nNumUniqueCalc)),
                ( (nNumUniqueCalcCumul!=0) ? (Cstring)(min(100.0,
                    (100.*(float)(ptSCumul->nNumMults + ptSCumul->nNumSingles)) / max(1,nNumUniqueCalcCumul))) : (Cstring)"---" ) );
#endif
        }
        else if (0 < nNumUniqueCalc)
        {
            printf("%6.2lf -%5.2lf%8d%8d%7d%8d%8d%7d%7s%7.1f ",
                min(fResoLow, 999.99), fResoHigh,
                nNumUniqueCalc,
                ptStats->nNumRefs, ptStats->nNumRejects,
                ptStats->nNumMults, ptStats->nNumSingles,
                nNumUniqueFound,
                "---",
                min(100., 100. * (float) nNumUniqueFound / (float) max(1,nNumUniqueCalc)));
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddCompReslnRow( fResoLow, fResoHigh, nNumUniqueCalc,
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumMults, ptStats->nNumSingles,
                nNumUniqueFound, (Cstring)"---",
                (Cstring) min(100., 100. * (float) nNumUniqueFound / (float) max(1,nNumUniqueCalc)),
                ( (nNumUniqueCalcCumul!=0) ? (Cstring)(min(100.0,
                    (100.*(float)(ptSCumul->nNumMults + ptSCumul->nNumSingles)) / max(1,nNumUniqueCalcCumul))) : (Cstring)"---" ) );
#endif
        }
        else
        {
            printf("%6.2lf -%5.2lf%8d%8d%7d%8d%8d%7d%7s%7s ",
                min(fResoLow, 999.99), fResoHigh,
                nNumUniqueCalc,
                ptStats->nNumRefs, ptStats->nNumRejects,
                ptStats->nNumMults, ptStats->nNumSingles,
                nNumUniqueFound,
                "---", "---");
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vAddCompReslnRow( fResoLow, fResoHigh, nNumUniqueCalc,
                ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumMults, ptStats->nNumSingles,
                nNumUniqueFound, (Cstring)"---", (Cstring)"---",
                ( (nNumUniqueCalcCumul!=0) ? (Cstring)(min(100.0,
                    (100.*(float)(ptSCumul->nNumMults + ptSCumul->nNumSingles)) / max(1,nNumUniqueCalcCumul))) : (Cstring)"---" ) );
#endif
        }
        if (nNumUniqueCalcCumul != 0)
            printf("%6.1f\n", min(100.0,
     
                (100.* (float)(  ptSCumul->nNumMults
                    + ptSCumul->nNumSingles))
                / nNumUniqueCalcCumul));
        else
            printf("%6s\n", "---");
        fLow = fLow + *pfSlope;
    }
    nNumUniqueFound = ptSCumul->nNumMults + ptSCumul->nNumSingles;
    printf(pcLine);
    printf("%6.2lf -%5.2lf%8d%8d%7d%8d%8d%7d ",
        1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)),
        1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)),
        nNumUniqueCalcCumul,
        ptStats->nNumRefs, ptStats->nNumRejects,
        ptStats->nNumMults, ptStats->nNumSingles, ptStats->nNumSingles+ptStats->nNumMults);
    if (0 != nNumUniqueCalcCumul)
    {
        printf("%6.2f ", (float) (ptStats->nNumUsed + ptStats->nNumSingles)
            / (float) max(1,nNumUniqueFound));
        printf("%6.1f ",  min(100., 100.* (float) nNumUniqueFound / max(1,nNumUniqueCalcCumul)));
        printf("%6.1f\n", min(100., 100.* (float) nNumUniqueFound / max(1,nNumUniqueCalcCumul)));
    }
    else
    {
        printf("%7s%7s%7s\n", "---", "---", "---");
    }

    // Save some calculations of the number of unique reflns for
    // use in ::nTableRMergeVsResolution()

    ptStats--;
    m_nNumUniqueFoundCumul = nNumUniqueFound;
    m_nNumUniqueCalcCumul  = nNumUniqueCalcCumul;
    m_nNumUniqueFoundShell = ptStats->nNumMults    + ptStats->nNumSingles;
    m_nNumUniqueCalcShell  = nNumUniqueCalc;

#ifdef SSI_PC
    CCrclHelper::GetInstance()->vCompleteCompReslnRow();
#endif

    //+jwp
    printf("\n\nRedundancy vs Resolution\n");
    printf(pcLine);
    printf("%13s%8s%47s %5s %5s\n",
        "Resolution  ", "  Calc",
        "Percent of reflections measured N times, N = ", "%Comp", "%Comp");
    printf("%13s%8s%5s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n",
        "   range    ", "unique", "0", "1", "2", "3",
        "4", "5-8", "9-12", ">12", "shell", "cumul" );

    printf(pcLine);
    ptStats  =  m_ptStatsReso;
    ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
    ptSCumul->vInitStats();                          // Initialize cumulative stats
    fLow     =  m_a2fRangeReso[0];
    pfSlope  = &m_fSlopeReso;
    int j;
    int nNumMultBinsM1 = 7;

    nNumUniqueCalcCumul=0;
    nNumUniqueFoundCumul=0;
    for (i = 0; i < m_nNumBinsReso; i++, ptStats++)
    {
        fResoLow              = 1.0 / max(1e-10,pow(fLow, 0.333333));
        fResoHigh             = 1.0 / max(1e-10,pow(fLow + *pfSlope, 0.333333));
        nNumUniqueCalc = anResoCounts[i];
        nNumUniqueFound        = ptStats->nNumMults    + ptStats->nNumSingles;
        ptStats->anNumInBinMults[0]   = nNumUniqueCalc-nNumUniqueFound;
        nNumUniqueCalcCumul         += nNumUniqueCalc;
        nNumUniqueFoundCumul        += nNumUniqueFound;
        printf("%6.2lf -%5.2lf%8d", fResoLow, fResoHigh, nNumUniqueCalc);

        for (j = 0; j < nNumMultBinsM1 + 1; j++)
        {
            ptSCumul->anNumInBinMults[j] += ptStats->anNumInBinMults[j];
            if (0 == j)
                printf("%5.1f",
                    max(0.0, min(100.0,
                        (float) ptStats->anNumInBinMults[j]
                        / (float) max(1,nNumUniqueCalc) * 100.0))
                    );
            else
                printf("%6.1f",
                    max(0.0, min(100.0,
                        (float) ptStats->anNumInBinMults[j]
                        / (float) max(1,nNumUniqueCalc) * 100.0))
                    );
        }
        printf("%6.1f%6.1f",
            max(0.0, min(100., 100. * (float) nNumUniqueFound / (float) max(1,nNumUniqueCalc))),
            max(0.0, min(100., 100.* ((float)nNumUniqueFoundCumul) / max(1,nNumUniqueCalcCumul)))
            );


        printf("\n");
        fLow = fLow + *pfSlope;
    }
    ptStats = ptSCumul;
    printf(pcLine);
    printf("%6.2lf -%5.2lf%8d%5.1f",
        1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)),
        1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)),
        (int) nNumUniqueCalcCumul,
        (float)ptStats->anNumInBinMults[0] / (float) max(1,nNumUniqueCalcCumul) * 100.0);
    for (j = 1; j < nNumMultBinsM1+1; j++)
    {
        printf("%6.1f",
            (float)ptStats->anNumInBinMults[j] / (float) max(1,nNumUniqueCalcCumul) * 100.0);
    }
    m_fOverallPercentComplete = 
        max(0.0, min(1.0, ((float) nNumUniqueFoundCumul) / max(1,nNumUniqueCalcCumul)));
    printf("%6.1f%6.1f",100.0*m_fOverallPercentComplete,100.0*m_fOverallPercentComplete);   

    printf("\n");

    
    if (0 < nNumUniqueCalcCumul)
    {
        // Form some cumulative statistics

        printf("\n");
        printf(pcLine);
        printf("%13s          %s", "Resolution  ",
            "Percent of reflections measured AT LEAST N times, N =\n");
        printf("%13s%8s%5s %5s %5s %5s %5s %5s %5s %5s\n",
            "   range    ", " ", " ",  "13", "9", "5",
            "4", "3", "2", "1" );
        printf(pcLine);
        printf("%6.2lf -%5.2lf%13s",
            1.0 / max(1e-10,pow((double)m_a2fRangeReso[0], 0.333333)),
            1.0 / max(1e-10,pow((double)m_a2fRangeReso[1], 0.333333)),
            " ");
        for (j = nNumMultBinsM1; j > 0; j--)
        {
            ptStats->anNumInBinMults[j] += ptStats->anNumInBinMults[j+1];            
            printf("%6.1f", 100.0*((float)ptStats->anNumInBinMults[j]) / max(1,nNumUniqueCalcCumul));

        }
        printf("\n");
        printf("\n");
    }
    fflush(stdout);    
    return (0);
}

int dtscaleaverage_float_cmp(const void* a,const void* b) {
    float* pfa = (float*) a;
    float* pfb = (float*) b;
    if (*pfa>*pfb) {
        return 1;
    } else if (*pfa==*pfb)
        return 0;
    else
        return -1;
};


int Cscaleaverage::nTableAnomSignal() {
    float* pfPlus;                    // For each group of sym-equiv reflections, this will contain the average Plus/Minus value for the group.
    float* pfMinus;                   
    int    nPlusMinusCount;           
    int    nPlusCount;
    int    nMinusCount;

    int    nRefSort,nStartIndex,nEndIndex;
    int    nLastHKL,nThisHKL;
    int    nRef;
    int    a3nHKL[3];
    int    nx,ny;
    double  f0,f1;
    float  fPlusFound,fMinusFound;
    int nPlusFound,nMinusFound;
    Cspacegroup& oSpace = *m_oCrystal.m_poSpacegroup;
    
    pfPlus = new float[m_nNumRefs];
    pfMinus = new float[m_nNumRefs];
    
    nPlusMinusCount = 0;
    nPlusCount = 0;
    nMinusCount = 0;
    
    nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    nStartIndex = 0;
    
    for (nRefSort=0;nRefSort<m_nNumRefs+1;nRefSort++)
    {
        nRef=m_pnIndex[nRefSort];
        
        if (nRefSort<m_nNumRefs)
            nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        if ((nRefSort != m_nNumRefs) && ((nThisHKL >> 1)== (nLastHKL >> 1)))
            continue;
        nEndIndex = nRefSort - 1;
        
        fPlusFound = 0.0;
        fMinusFound = 0.0;
        nMinusFound = 0;
        nPlusFound = 0;
        for (nx = nStartIndex; nx <= nEndIndex; nx++) {
            if (m_pnReject[nx]!=1) {
                oSpace.nReduceHKL(&m_oReflnlist[m_pnIndex[nx]],a3nHKL,&ny);
                if (ny==-1) {
                    nMinusFound++;
                    fMinusFound += m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                };
                if (ny==1) {
                    nPlusFound++;
                    fPlusFound += m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                };
            };
        };
        if ((nPlusFound) && (nMinusFound)) {
            if (nPlusFound) 
                pfPlus[nPlusCount++] = fPlusFound/nPlusFound;
            if (nMinusFound)
                pfMinus[nMinusCount++] = fMinusFound/nMinusFound;
            
            nPlusMinusCount++;
        };

        nStartIndex=nRefSort;
        if (nRefSort<m_nNumRefs)
            nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    }   

    if (nPlusMinusCount<20) {
        printf("Too few Friedel pairs available for anomalous signal test.\n");
        goto exit_place;
    };

    

    // Sort both data arrays.
    qsort(pfPlus,nPlusCount,sizeof(float),dtscaleaverage_float_cmp);
    qsort(pfMinus,nMinusCount,sizeof(float),dtscaleaverage_float_cmp);

    double fMax,fMin;
    fMax = min(pfPlus[nPlusCount-1],pfMinus[nMinusCount-1]);
    fMin = max(pfPlus[0],pfMinus[0]);

    
    // Restore the data so that correct cdf values are computed.
    for (nx=0,ny=0;nx<nPlusCount;nx++) {
        if ((nx+1>=nPlusCount) || (pfPlus[nx+1]>=fMin))
            pfPlus[ny++] = pfPlus[nx];
    };
    nPlusCount = ny;

    for (nx=0,ny=0;nx<nMinusCount;nx++) {
        if ((nx+1>=nMinusCount) || (pfMinus[nx+1]>=fMin))
            pfMinus[ny++] = pfMinus[nx];
    };
    nMinusCount = ny;
    

    // Find the maximum difference between the two distributions.
    // For every entry in the pfMinus array between fMax and fMin,
    // we lookup the nearest value in the pfMinus array.  The difference in indices then 
    // gives the cdf function value differences.    

    double fMaxDiff;
    int nMinusPoint;
    int nPlusPoint;
    fMaxDiff = 0.0;
    nMinusPoint = 0;
    for (nPlusPoint=0;nPlusPoint<nPlusCount;nPlusPoint++) {
        if ((pfPlus[nPlusPoint]<fMax) && (pfPlus[nPlusPoint]>fMin)) {
            // We can be sure to find two points in pfMinus that bound the point pfPlus[nPlusPoint].
            while ((nMinusPoint+1<nMinusCount) && (pfPlus[nPlusPoint]>pfMinus[nMinusPoint+1]))
                nMinusPoint++;
            while ((nMinusPoint>0) && (pfPlus[nPlusPoint]<pfMinus[nMinusPoint]))
                nMinusPoint--;
            // DEBUG:  Should NOT have:
            if ((pfPlus[nPlusPoint]>pfMinus[nMinusPoint+1]) || (pfPlus[nPlusPoint]<pfMinus[nMinusPoint])) {
                printf("Error in anomalous test.\n");
                goto exit_place;
            };
            f0 = fabs(((float) nPlusPoint)/max(1,nPlusCount) - ((float) nMinusPoint)/max(1,nMinusCount));
            if (f0>fMaxDiff) {
                fMaxDiff = f0;
            };
        };
    };
    
   
    double fQSStat;
    double fStat;

    f0 = sqrt( (double) (nPlusCount*nMinusCount)/max(1.0,(double)(nPlusCount + nMinusCount)));
    fQSStat = (f0 + 0.12 + 0.11/max(1e-10,f0))*fMaxDiff;
    f0 = 0.0;
    f1 = 1.0;
    for (nx=0; nx<100;nx++) {
        f0 += 2.0*f1*exp(-2.0*(nx+1)*(nx+1)*fQSStat*fQSStat);
        f1*=-1.0;
    };
    fStat = f0;
        
    printf("Calculated Kolmogorov-Smirnov comparison of I+ and I- intensity distributions\n");
    if (fStat<0.5) {
           printf("indicates that I+ != I- at a significance level of %6.2f%%\n"
           "The test used %d Friedel pairs for comparison.\n",fStat*100.0,nPlusMinusCount);
    } else {
           printf("indicates that I+ == I-.\n"
           );
    };
    printf("The Kolmogorov-Smirnov statistic and it's significance level were calculated\n");


exit_place:
    delete[] pfPlus;
    delete[] pfMinus;
    return 0;
};

#ifdef SSI_PC
int Cscaleaverage::nWriteReflns(Cstring& sName,
                                Cstring& sUnavgName,
                                Cstring& sRejectsName,
                                Cstring& sScaleName,
                                bool bWriteHeader,
                                bool bWriteUnavgHeader, 
                                bool bTexsan,
                                bool bTexsan2,
                                bool bShelx,
                                bool bAnom,
                                                              bool bTexsan2andShelx) 
#else
int Cscaleaverage::nWriteReflns(Cstring& sName,
                                Cstring& sUnavgName,
                                Cstring& sRejectsName,
                                Cstring& sScaleName,
                                bool bWriteHeader,
                                bool bWriteUnavgHeader, 
                                bool bTexsan,
                                bool bTexsan2,
                                bool bShelx,
                                bool bAnom) 
#endif
{
    
    int nStat;
    int nRefSort,nRef;
    int nx,ny;
    int nStartIndex,nEndIndex;
    int nLastHKL;
    int nThisHKL;
    float fSigma,fIntensity;
    
    int a3nHKL[3];
    int nPlusFound;
    int nMinusFound;
    int nGroup;
    int nGroupPlus,nGroupMinus;
    int nGroupSigmasUsed,nMergeSigmasUsed;  // Computed when averaging.
    
    float fGroupAverage;
    float fAveragePlus;
    float fAverageMinus;
    float fGroupSigma;
    float fSigmaPlus;
    float fSigmaMinus;
    float fInvVarPlus;
    float fInvVarMinus;
    float fGroupInvVar;
    float fGroupSigmaMerge;
    
//+2010-07-20
    bool   bUseSigmaWeights = USE_SIGMA_WEIGHTS;
    double fReflnWeight = 1.0;
//-2010-07-20

    Cspacegroup& oSpace = *m_oCrystal.m_poSpacegroup;
    Creflnlist oOutList;
    
    // Consistency checks.

#ifdef SSI_PC 
      Cstring sUnavgName1;
      Cstring sUnavgName2;
      sUnavgName1 = "f2plus.dat";
      sUnavgName2 = "shelxl.hkl";
      if( bTexsan2andShelx ) {
              bTexsan2 = TRUE;
              bShelx = FALSE;
      }
      else {
#endif
    // Not both Texsan flags.
    if ((bTexsan) && (bTexsan2)) 
      {
        printf("WARNING:  Only one of -texsan, -texsan2, and -shelx can be specifed!\n          Choosing -texsan2 format.\n");
        bTexsan=FALSE;
      }
    if ((bTexsan) && (bShelx)) 
      {
        printf("WARNING:  Only one of -texsan, -texsan2, and -shelx can be specifed!\n          Choosing -shelx format.\n");
        bTexsan=FALSE;
      }
    if ((bTexsan2) && (bShelx)) 
      {
        printf("WARNING:  Only one of -texsan, -texsan2, and -shelx can be specifed!\n          Choosing -shelx format.\n");
        bTexsan2=FALSE;
      }
#ifdef SSI_PC 
      }
#endif
    
    // scaleanom ==> anom.
    
    if (m_bScaleAnom)
      {
        // Set the packed HKL entry to be I- for all reflections in the list.
        // Since the nPackedHKL is even for I- and odd (has 1 added) for I+, turn the odds into evens

        for (nRefSort=0;nRefSort<m_nNumRefs;nRefSort++)
            m_oReflnlist[nRefSort].vSetField(m_oReflnlist.m_nFI_nPackedHKL,
            (int) ((((unsigned int) m_oReflnlist[nRefSort].nGetField(m_oReflnlist.m_nFI_nPackedHKL))/2)*2));
        bAnom = TRUE;
      }

#ifdef SSI_PC        
      if( !bTexsan2andShelx ) {
#else
    if (((bShelx) || (bTexsan) || (bTexsan2)) && (sUnavgName.length()==0)) 
      {
        if (bShelx)
            sUnavgName = "shelxl.hkl";
        else
            sUnavgName = "f2plus.dat";
        printf("\nWARNING!  -texsan -texsan2 and -shelx require -ref command.\n");
        printf("INFO:  Setting file name to %s\n",sUnavgName.string());
      }
#endif
#ifdef SSI_PC 
      }
#endif
    
    // Do we need to write out a reflection list?
    
    if (*sName.string() != 0)
      {
        if (bAnom)
          {
            m_nFI_fIntensityPlus  = oOutList.nExpandGetField(Creflnlist::ms_sfIntensityPlus);
            m_nFI_fSigmaIPlus     = oOutList.nExpandGetField(Creflnlist::ms_sfSigmaIPlus);
            m_nFI_fIntensityMinus = oOutList.nExpandGetField(Creflnlist::ms_sfIntensityMinus);
            m_nFI_fSigmaIMinus    = oOutList.nExpandGetField(Creflnlist::ms_sfSigmaIMinus);
          }
        oOutList.nAddExtra(m_oReflnlist);
        {
          // This is a kludge
          Cimage_header oTHead;
          float a2fTemp[2];
          a2fTemp[0] = 1.0;
          a2fTemp[1] = m_fWavelength;
          oTHead.nReplaceValue(D_K_SourceWavelength, 2, a2fTemp);
          oOutList.nPutHeaderInfo(oTHead, D_K_SourceWavelength, false);
        }
        
        Crefln oRefln(&oOutList);
        
        nGroupSigmasUsed = 0;
        nMergeSigmasUsed = 0;
        nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        nStartIndex = 0;

        for (nRefSort=0;nRefSort<m_nNumRefs+1;nRefSort++)
          {
            nRef=m_pnIndex[nRefSort];
            
            if (nRefSort<m_nNumRefs)
                nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
            if ((nRefSort != m_nNumRefs) && (nThisHKL == nLastHKL))
                continue;
            nEndIndex = nRefSort - 1;
            
            fGroupAverage=0.0;
            fAveragePlus=0.0;
            fAverageMinus=0.0;
            fGroupSigma=0.0;
            fSigmaPlus=0.0;
            fSigmaMinus=0.0;
            fInvVarPlus=0.0;
            fInvVarMinus=0.0;
            fGroupInvVar=0.0;
            fGroupSigmaMerge=0.0;
            
            nPlusFound=0;
            nMinusFound=0;
            
            nGroup=0;
            nGroupPlus=0;
            nGroupMinus=0;
            
            // Count the starting and ending locations for I+ and I-.
            // We have built this so that it would be portable, even if we decide
            // to change whether I+ or I- reflections come first in the list.
            
            int nCentric = 0;

            for (nx = nStartIndex; nx <= nEndIndex; nx++)
                if (m_pnReject[nx]!=1)
                  {                

		    //+-JWP The following is probably not needed since
		    //      the reflnlist probably already has the fields
		    //      for the nAnomFlag and nCentPhase
		    // TODO: See if next line is redundant!

                    nCentric = oSpace.nReduceHKL(&m_oReflnlist[m_pnIndex[nx]],a3nHKL,&ny);
		    //+JWP 2007-06-10
		    if (0 > nCentric)
		      {
			// These are usually systematically absent reflections
			// that are in the input file.
			// They should be excluded the next time through ...
			// or perhaps rejected

			m_pnReject[nx] = 2;

			/***
		      cout << "ERROR after nReduceHKL, hkl: "
                           << a3nHKL[0] << ", "
                           << a3nHKL[1] << ", "
                           << a3nHKL[2] << endl;
			***/
		      }
		    else
		      {
			nCentric = 0;  // Comment out this line if 
		                     // you do NOT want the centric I+ and I- treated separately
		      }
                    nGroup++;

		    //+JWP 2007-06-10
		    // This bit will merge so-called I+ and I- of centric reflns
                    if ( (ny==1) || (nCentric > 0) )
                        nPlusFound=1;
                    else if (ny==-1)
                        nMinusFound=1;
		    /* This original code will keep I+ and I- separate for centric reflns
                    if (ny==1)
                        nPlusFound=1;
                    else if (ny==-1)
                        nMinusFound=1;
		    */
		    //-JWP 2007-06-10

                    fGroupSigmaMerge = m_pfMergeSigma[nx];
                    fIntensity = m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                    fSigma     = fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]);
                    
                    
                    if (0.0 >= fSigma) 
                      {
                        // Should not reach this code.
                        cout << "WARNING: sigma <= 0: " << fSigma << " for hkl: "
                            << a3nHKL[0] << " " << a3nHKL[1] << " " << a3nHKL[2]
                            << "  ... reset to 9.0e10." << endl;
                        fSigma = 9.0E10f;
                      }
                    
		    if (bUseSigmaWeights) fReflnWeight = max(1.0, (fSigma*fSigma));
                    fGroupAverage += fIntensity / fReflnWeight;
                    fGroupSigma   +=  1.0/ fReflnWeight;
                    fGroupInvVar  += 1.0 / fReflnWeight;
		    //+JWP 2007-06-10
		    // This next if will treat centric measurements as an I+ if centric > 0
                    if ( (ny==1) || (nCentric > 0) )
		      //-JWP if (1 == ny)
                      {
                        fAveragePlus += fIntensity / max(1.0,(fSigma*fSigma));
                        fSigmaPlus   += 1.0/ max(1.0,(fSigma * fSigma));
                        fInvVarPlus  += 1.0 / max(1.0,(fSigma*fSigma));
                        nGroupPlus++;
                      }
                    if (-1 == ny)
                      {
                        fAverageMinus += fIntensity/max(1.0,(fSigma*fSigma));
                        fSigmaMinus   +=  1.0/ max(1.0,(fSigma * fSigma));
                        fInvVarMinus  += 1.0/max(1.0,(fSigma*fSigma));
                        nGroupMinus++;
                      }
                  }
                if (0 != nGroup)
                  {
                    fGroupAverage /= max(1e-10,fGroupInvVar);
                    fGroupSigma    = 1.0/sqrt((double)max(1e-30,fGroupSigma));
                    
                    if (nPlusFound==1)
                      {
                        fAveragePlus /= max(1e-10,fInvVarPlus);
                        fSigmaPlus    = 1.0/sqrt((double)max(1e-30,fSigmaPlus));
                      }
                    if (nMinusFound==1)
                      {
                        fAverageMinus /= max(1e-10,fInvVarMinus);
                        fSigmaMinus    = 1.0/sqrt((double)max(1e-30,fSigmaMinus));
                      }
                    
                    // Add the reflection(s).
                    
                    // Build the HKL entries.
                    
                    oRefln.vSetIntensity(fGroupAverage);
                    oRefln.vSetSigmaI(fGroupSigma);
                    if (fGroupSigma>fGroupSigmaMerge)
                      nGroupSigmasUsed++;
                    else
                      nMergeSigmasUsed++;
                    
                    if (bAnom)
                      {
                        oRefln.vSetH(a3nHKL[0]);
                        oRefln.vSetK(a3nHKL[1]);
                        oRefln.vSetL(a3nHKL[2]);
                        if (nPlusFound==1)
                          {
                            oRefln.vSetField(m_nFI_fIntensityPlus,fAveragePlus);
                            oRefln.vSetField(m_nFI_fSigmaIPlus,fSigmaPlus);
                          }
                        else
                          {
                            oRefln.vSetField(m_nFI_fIntensityPlus,(float) -1.0);
                            oRefln.vSetField(m_nFI_fSigmaIPlus,(float) -1.0);
                          }
                        
                        if (nMinusFound==1)
                          {
                            oRefln.vSetField(m_nFI_fIntensityMinus,fAverageMinus);
                            oRefln.vSetField(m_nFI_fSigmaIMinus,fSigmaMinus);
                          }
                        else
                          {
                            oRefln.vSetField(m_nFI_fIntensityMinus,(float) -1.0);
                            oRefln.vSetField(m_nFI_fSigmaIMinus,(float) -1.0);
                          }
                      }
                    else
                      {
                        // This first case is no longer reached,
                        // since anom is always turned on when scaleanom is set.
                        // I have left it here so that we can move back to the old way if necc.
                        
                        if ((m_bScaleAnom) && (nMinusFound==1))
                          {
                            oRefln.vSetH(-a3nHKL[0]);
                            oRefln.vSetK(-a3nHKL[1]);
                            oRefln.vSetL(-a3nHKL[2]);
                          }
                        else
                          {
                            oRefln.vSetH(a3nHKL[0]);
                            oRefln.vSetK(a3nHKL[1]);
                            oRefln.vSetL(a3nHKL[2]);
                          }
                      }
                    oOutList.nInsert(&oRefln);
                  }
                
                nStartIndex=nRefSort;
                if (nRefSort<m_nNumRefs)
                  nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
          }
        
        if (   (nMergeSigmasUsed * 10 < (nGroupSigmasUsed + nMergeSigmasUsed))
            || (nGroupSigmasUsed * 10 < (nGroupSigmasUsed + nMergeSigmasUsed)))
          {
            //            printf ("\nWARNING!  Check your error model!\n");
            printf ("\n");
            printf("%10d of %d sigmas computed from merging statistics.\n",
                nMergeSigmasUsed, nGroupSigmasUsed + nMergeSigmasUsed);
            printf("%10d of %d sigmas computed from error model.\n",
                nGroupSigmasUsed, nGroupSigmasUsed + nMergeSigmasUsed);
          }
        
        if (!bWriteHeader) 
          oOutList.vSetNoWriteHeader();
        
        // Always write averaged list in TEXT mode.
        oOutList.eSetWriteBinary(eReflnoutput_text);
        
        cout << "\nWriting " << sName << " ..." << endl << flush;
        nStat = oOutList.nWrite(sName);
        if (0 == nStat)
          {
            cout << "... done writing.\n" << flush;
            bCalcAnomalousSignal((Creflnlist*)&oOutList);
          }
        else
          {
            cout << "ERROR writing reflnlist " << sName << "!\n"
                << endl;
          }
      }
    
    if (sRejectsName.length() !=0) 
      {
        Creflnlist oRejects(m_oReflnlist);
        for (nx=0;nx<m_oReflnlist.nGetNumReflns();nx++) 
          {
            if (m_pnReject[nx] == 1)
              {
                oRejects.nInsert(&m_oReflnlist[m_pnIndex[nx]]);
              }
          }
        oRejects.nWrite(sRejectsName);
      }

    if (*sScaleName.string() != 0) 
      {
        int nFI_fScale  = m_oReflnlist.nExpandGetField("fScale");
        char* pcSelect;
        int nFields;
        int* pnSortIndex;

        nFields  = m_oReflnlist.nGetNumFields();
        pcSelect = new char [nFields + 1];
        for (nx = 0; nx < nFields; nx++)
          pcSelect[nx] = m_oReflnlist.m_pcSelect[nx] + 1;
        pcSelect[nFields] = m_oReflnlist.m_pcSelect[nFields];

        (void) m_oReflnlist.nSelectField(eReflnField_int_type,
            m_oReflnlist.m_nFI_nH, char(0), pcSelect);
        (void) m_oReflnlist.nSelectField(eReflnField_int_type,
            m_oReflnlist.m_nFI_nK, char(0),        pcSelect);
        (void) m_oReflnlist.nSelectField(eReflnField_int_type,
            m_oReflnlist.m_nFI_nL, char(0),        pcSelect);
        (void) m_oReflnlist.nSelectField( eReflnField_float_type,
            m_oReflnlist.m_nFI_fIntensity, char(0),pcSelect);
        (void) m_oReflnlist.nSelectField( eReflnField_float_type,
            m_oReflnlist.m_nFI_fSigmaI, char(0), pcSelect);
        (void) m_oReflnlist.nSelectField( eReflnField_float_type,
            nFI_fScale, char(0),pcSelect);

        for (nx = 0; nx < m_oReflnlist.nGetNumReflns(); nx++) 
          {
            m_oReflnlist[nx].vSetField(nFI_fScale,(float) (m_oReflnlist[nx].fGetIntensity() / m_pfOrigIntensity[nx]*((m_pnReject[m_pnRevIndex[nx]]==1)?(-1.0):(1.0))));
          }        
        pnSortIndex = m_oReflnlist.pnGetSortIndex();
        m_oReflnlist.m_pnSortIndex = NULL;
        nStat = m_oReflnlist.nWrite(sScaleName, NULL, pcSelect);
        m_oReflnlist.m_pnSortIndex = pnSortIndex;

        delete[] pcSelect;
      };
#ifdef SSI_PC
      if (*sUnavgName1.string() != 0)
      {
#else
    if (*sUnavgName.string() != 0)
      {
        char *pcSelect = NULL;
        int  *pnSelectFortran = NULL;
        int  *pnSelectFortran2 = NULL;
        bool bSuccess = FALSE;
#endif
        
        // Set all sigma values to the value they have in the error model.
        
        for (nx = 0; nx < m_oReflnlist.nGetNumReflns(); nx++) 
          {
            if (m_pnReject[nx] == 1)
                m_oReflnlist[m_pnIndex[nx]].vSetSigmaI(-fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]));
            else
              m_oReflnlist[m_pnIndex[nx]].vSetSigmaI(fCalcSigma(m_oReflnlist[m_pnIndex[nx]],m_pnIndex[nx]));
          }

#ifdef SSI_PC
Recycle:
        char *pcSelect = NULL;
        int  *pnSelectFortran = NULL;
        int  *pnSelectFortran2 = NULL;
        bool bSuccess = FALSE;
#endif
        
        if ((!bTexsan) && (!bTexsan2) && (!bShelx)) 
          {
            if (!bWriteUnavgHeader)
                m_oReflnlist.vSetNoWriteHeader();
            bSuccess = TRUE;
          } 
        else
          {
            int nFields;
            
            //why are all fields printing?
            //    modifiy sigma and intensity data here.
           
            m_oReflnlist.vSetNoWriteHeader();
            
            nFields  = m_oReflnlist.nGetNumFields();
            pcSelect = new char [nFields + 1];
            pnSelectFortran = new int[nFields + 1];
            pnSelectFortran2 = new int[nFields + 1];
            
            // "Deselect" all fields
            
            for (nx = 0; nx < nFields; nx++)
                pcSelect[nx] = m_oReflnlist.m_pcSelect[nx] + 1;
            
            // Keep the global selection the same
            
            pcSelect[nFields] = m_oReflnlist.m_pcSelect[nFields];
            
            // Select only the 5 fields we want to keep
            
            (void) m_oReflnlist.nSelectField(eReflnField_int_type,
                m_oReflnlist.m_nFI_nH, char(0), pcSelect);
            (void) m_oReflnlist.nSelectField(eReflnField_int_type,
                m_oReflnlist.m_nFI_nK, char(0),        pcSelect);
            (void) m_oReflnlist.nSelectField(eReflnField_int_type,
                m_oReflnlist.m_nFI_nL, char(0),        pcSelect);
            (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                m_oReflnlist.m_nFI_fIntensity, char(0),pcSelect);
            (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                m_oReflnlist.m_nFI_fSigmaI, char(0), pcSelect);

            if (bShelx)
              {
                m_oReflnlist.nExpandGetField(m_oReflnlist.ms_snNonunfFlag);
                
                // Swap some numbers.  This is very important, since the ordering
                // of absorption information for shelx is different.
                int a6nSInfo[6];
                double a6fSInfo[6];
                a6nSInfo[0] = m_oReflnlist.m_nFI_fS0vec[0];
                a6nSInfo[1] = m_oReflnlist.m_nFI_fS0vec[1];
                a6nSInfo[2] = m_oReflnlist.m_nFI_fS0vec[2];
                a6nSInfo[3] = m_oReflnlist.m_nFI_fSvec[0];
                a6nSInfo[4] = m_oReflnlist.m_nFI_fSvec[1];
                a6nSInfo[5] = m_oReflnlist.m_nFI_fSvec[2];
                qsort(&a6nSInfo[0],6,sizeof(int),int_cmp);
                for (nx = 0; nx < m_oReflnlist.nGetNumReflns(); nx++)
                  {
                    m_oReflnlist[nx].vSetField(m_oReflnlist.m_nFI_fCalcRotMid,(float) 1.0);
                    a6fSInfo[0] = m_oReflnlist[nx].fGetField(m_oReflnlist.m_nFI_fS0vec[0]);
                    a6fSInfo[1] = m_oReflnlist[nx].fGetField(m_oReflnlist.m_nFI_fSvec[0]);
                    a6fSInfo[2] = m_oReflnlist[nx].fGetField(m_oReflnlist.m_nFI_fS0vec[1]);
                    a6fSInfo[3] = m_oReflnlist[nx].fGetField(m_oReflnlist.m_nFI_fSvec[1]);
                    a6fSInfo[4] = m_oReflnlist[nx].fGetField(m_oReflnlist.m_nFI_fS0vec[2]);
                    a6fSInfo[5] = m_oReflnlist[nx].fGetField(m_oReflnlist.m_nFI_fSvec[2]);
                    for (ny = 0; ny < 6; ny++) 
                      m_oReflnlist[nx].vSetField(a6nSInfo[ny],(float) a6fSInfo[ny]);
                  };

                pnSelectFortran[m_oReflnlist.m_nFI_nH] = 4;
                pnSelectFortran[m_oReflnlist.m_nFI_nK] = 4;
                pnSelectFortran[m_oReflnlist.m_nFI_nL] = 4;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fIntensity] = 8;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSigmaI] = 8;                
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fCalcRotMid] = 4;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fS0vec[0]] = 8;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSvec[0]] = 8;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fS0vec[1]] = 8;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSvec[1]] = 8;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fS0vec[2]] = 8;
                pnSelectFortran[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSvec[2]] = 8;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fIntensity] = 2;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSigmaI] = 2;                
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fCalcRotMid] = 0;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fS0vec[0]] = 5;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSvec[0]] = 5;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fS0vec[1]] = 5;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSvec[1]] = 5;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fS0vec[2]] = 5;
                pnSelectFortran2[m_oReflnlist.m_nIntReflnFields + m_oReflnlist.m_nFI_fSvec[2]] = 5;

                (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                    m_oReflnlist.m_nFI_fCalcRotMid, char(0), pcSelect);
                (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                    m_oReflnlist.m_nFI_fSvec[0], char(0), pcSelect);
                (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                    m_oReflnlist.m_nFI_fSvec[1], char(0), pcSelect);
                (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                    m_oReflnlist.m_nFI_fSvec[2], char(0), pcSelect);
                (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                    m_oReflnlist.m_nFI_fS0vec[0], char(0), pcSelect);
                (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                    m_oReflnlist.m_nFI_fS0vec[1], char(0), pcSelect);
                (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                    m_oReflnlist.m_nFI_fS0vec[2], char(0), pcSelect);
              };
            
            if (bTexsan2) 
              {
                float fInt, fSigma, fScale;
                float fLorentz, fPolarz, fLP;
                
                int nFI_fLorentz = m_oReflnlist.m_nFI_fLorentz;
                int nFI_fPolarz  = m_oReflnlist.m_nFI_fPolarz;
                int nFI_fTgeos  = m_oReflnlist.nGetFieldIndex("fTgeos");
                int nFI_fScale  = m_oReflnlist.nExpandGetField("fScale");
                
                // If LP information is not avail., 
                // then we can still run assuming that LP=1.0
                
                if ((0 > nFI_fLorentz) || (0 > nFI_fPolarz))
                  {
                    printf("WARNING: Lorentz or Polarz factor not found.\n");
                  }
                
                // NO LONGER TRUE: We cannot run if we don't have the fTgeos data.
                // STILL TRUE:     But we need the fTgeos column!
                
                if (0 > nFI_fTgeos) 
                  {
                    printf("WARNING: transmisson or scale factor not found.\n");
                  }  
                else
                  {
                    bSuccess = TRUE;
                    
                    // Use the 2nd field, whatever it is, to hold the pre-scaled
                    // Intensity, also with no LP correction
                    
                    (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                        2, char(0), pcSelect);
                    
                    (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                        nFI_fTgeos, char(0), pcSelect);
                    (void) m_oReflnlist.nSelectField( eReflnField_float_type,
                        nFI_fScale, char(0), pcSelect);
                    
                    Crefln *poRefln;
                    
                    for (nx = 0; nx < m_oReflnlist.nGetNumReflns(); nx++)
                      {
                        poRefln = m_oReflnlist.poGetRefln(nx);
                        //+-jwp 28-Feb-2001     fTgeos  = poRefln->fGetField(nFI_fTgeos);
                        
                        if ((0 < nFI_fLorentz) && ( 0 < nFI_fPolarz)) 
                          {
                            fLorentz= poRefln->fGetField(nFI_fLorentz);
                            fPolarz = poRefln->fGetField(nFI_fPolarz);
                            fLP     = fLorentz * fPolarz;
                          } 
                        else
                          {
                            fLP = 1.0f;
                          }

                        
                        fInt   = poRefln->fGetIntensity();
                        fSigma = poRefln->fGetSigmaI();

                        
                        // Compute the scale factor that was applied.
                        
                        if (m_pfOrigIntensity[nx] != 0.0) 
                          {
                            fScale  = fInt / m_pfOrigIntensity[nx];
                          }
                        
                        if (   (0.0 < fSigma)
                            && (0.0 < fLP)
                            && (0.0 < fScale)
                            //+-jwp 28-Feb-2001 && (0.0 < fTgeos) 
                            ) 
                          {
                            //old                         poRefln->vSetField(nFI_fScale,fScale);
                            poRefln->vSetField(nFI_fScale,(float)1.0);
                            poRefln->vSetField(nFI_fTgeos, (float)(1.0/max(1e-10,fScale)));
                            poRefln->vSetIntensity(m_pfOrigIntensity[nx]);
                            // +-tjn 22-Mar-2001 poRefln->vSetSigmaI(m_pfOrigSigmas[nx]);
                            poRefln->vSetSigmaI(fSigma/max(1e-10,fScale));
                            fInt   = m_pfOrigIntensity[nx] * fLP;
                            
                            poRefln->vSetField(2, fInt);
                          }  
                        else if (0.0 < fSigma) 
                          {
                            cout << "WARNING -texsan2 problem: ";
                            poRefln->nList(2);
                          }
                      } // end of for (nx)
                  } 
              }
            else
              {
                bSuccess = TRUE;
              }

            if ((bTexsan2) || (bShelx))
              {
                float fMaxFloatFound   = 0.0f;
                float fMaxFloatAllowed = 999999.f;
                float fInt,fSigma;
                Crefln* poRefln;
                for (nx = 0; nx < m_oReflnlist.nGetNumReflns(); nx++)
                  {
                    poRefln = m_oReflnlist.poGetRefln(nx);
                    fInt = poRefln->fGetIntensity();
                    fSigma = poRefln->fGetSigmaI();
                    fMaxFloatFound = max(fMaxFloatFound, fInt);
                    fMaxFloatFound = max(fMaxFloatFound, fSigma);
                  };
            
                if (fMaxFloatFound > fMaxFloatAllowed)
                  {
                    float fRescaleFactor = 1.0f;
                    float fAddValue = 1.0;
                    do 
                      {
                        fRescaleFactor = fMaxFloatAllowed / max(1.0,(fMaxFloatFound + fAddValue));
                        fAddValue     *= 10.0f;
                      } while (fRescaleFactor * fMaxFloatFound > fMaxFloatAllowed);
                
                printf("INFO:  Re-scaling intensities and sigmas by %f\n",
                    fRescaleFactor);
                for (nx = 0; nx < m_oReflnlist.nGetNumReflns(); nx++)
                  {
                    poRefln = m_oReflnlist.poGetRefln(nx);
                    fInt    = poRefln->fGetIntensity();
                    fSigma  = poRefln->fGetSigmaI();
                    poRefln->vSetIntensity(fInt * fRescaleFactor);
                    poRefln->vSetSigmaI(fSigma * fRescaleFactor);
                    fInt    = poRefln->fGetField(2);
                    poRefln->vSetField(2, fInt * fRescaleFactor);
                  }
                  }
              };
          }
        if (bSuccess)
          {
            m_oReflnlist.nSelect("-fSigmaI<=0.0");
#ifdef SSI_PC
                      if( bTexsan2andShelx ) {
                              if(bTexsan2) {
                                      cout << "Writing " << sUnavgName1 << " ..." << endl << flush;
                              }
                              else if(bShelx) { 
                                      cout << "Writing " << sUnavgName2 << " ..." << endl << flush;
                              }
                      }
                      else {
                              cout << "Writing " << sUnavgName << " ..." << endl << flush;
                      }
#else
            cout << "Writing " << sUnavgName << " ..." << endl << flush;
#endif        
            if ((bTexsan) || (bTexsan2) || (bShelx) || (!bWriteUnavgHeader))
              m_oReflnlist.eSetWriteBinary(eReflnoutput_text);
            if (bShelx)
              {
                m_oReflnlist.vSetFortranFormat(pnSelectFortran,pnSelectFortran2);
              };
#ifdef SSI_PC
                      if( bTexsan2andShelx ) {
                              if(bTexsan2) {
                                      nStat = m_oReflnlist.nWrite(sUnavgName1, NULL, pcSelect);
                              }
                              else if(bShelx) { 
                                      nStat = m_oReflnlist.nWrite(sUnavgName2, NULL, pcSelect);
                              }
                      }
                      else {
                              nStat = m_oReflnlist.nWrite(sUnavgName, NULL, pcSelect);
                      }
#else                         
            nStat = m_oReflnlist.nWrite(sUnavgName, NULL, pcSelect);
#endif 
            if (0 == nStat)
              {
                cout << "... done writing.\n" << flush;
              }
            else
              {
#ifdef SSI_PC
                              if( bTexsan2andShelx ) {
                                      if(bTexsan2) {
                                              cout << "ERROR writing reflnlist " << sUnavgName1 << "!\n"
                                                       << endl;
                                      }
                                      else if(bShelx) { 
                                              cout << "ERROR writing reflnlist " << sUnavgName2 << "!\n"
                                                       << endl;
                                      }
                              }
                              else {
                                      cout << "ERROR writing reflnlist " << sUnavgName << "!\n"
                                               << endl;
                              }
#else
                cout << "ERROR writing reflnlist " << sUnavgName << "!\n"
                     << endl;
#endif
              }
          }
        delete [] pcSelect;
        delete [] pnSelectFortran;
        delete [] pnSelectFortran2;
      }
#ifdef SSI_PC
        if( bTexsan2andShelx && bTexsan2 ) {
                bTexsan2 = FALSE;
                bShelx = TRUE;
                goto Recycle;
        }
#endif
    return (0);
}

void Cscaleaverage::vExpandReflnlist(Creflnlist *poReflnlist)
{
  // Add required fields to the working reflection list if they are not
  // already there..

  poReflnlist->nExpandGetField(Cscaleaverage::ms_sf2STLsq);
  poReflnlist->nExpandGetField(Cscaleaverage::ms_snBatchIdx);
  poReflnlist->nExpandGetField(Creflnlist::ms_snPackedHKL);
  poReflnlist->nExpandGetField(Creflnlist::ms_ssBatch);
  poReflnlist->nExpandGetField("fTgeos");
  poReflnlist->nExpandGetField("fScale");
}


int Cscaleaverage::nFitChiSquared(bool bPrint) 
{
    itr<int> anFirstSubGroupInResoGroup;        // [#ResoBins] For each resolution bin, the index of the first intensity sub group
    itr<int> anSubGroupReso;                    // [#TotalSubgroups] Resolution bin for each sub group
    itr<int> anSubGroupCount;                   // [#TotalSubgroups] Number of contributors to the sub group.
    
    itr<double> afSubGroupMaxIntensity;         // [#TotalSubgroups] Maximum intensity.
    itr<double> afSubGroupAvgIntensity;         // [#TotalSubGroups] Average intensity (important).
    itr<double> afSubGroupMinIntensity;         // [#TotalSubgroups] Minimum intensity
    
    itr<double> afSubGroupEAdd;                 // [#TotalSubgroups] Eadd
    itr<double> afSubGroupEMul;                 // [#TotalSubgroups] Emul

    int nTotalAvgs;

    itr<double> afResoGroupMin;                 // [#ResoBins] maximum m_pfRevIndexReso index for the resolution group.
    itr<double> afResoGroupMax;
    
    int      nResoGroups;
    int      nResoGroup;
    
    double fAvgEMulNumer;
    double fAvgEMulDenom;
    double fAvgEMul = -1.0;

    const int nMinAvgsPerSubGroup = 50;         // Important number. Indicates the minimum number of reflections needed to determine Eadd and/or Emul
    const int nMinAvgsPerGroup = 300;           // Important number. Indicates the minimum number of reflections needed in a resolution group (which is broken into sub-groups).


    int      nLastHKL,nThisHKL;
    int      nStartIndex,nEndIndex,nIndex;
    int      nRefSort;
    int      nRef;
    int      nStat;
    int      nx,ny;
    double   f0,f1;

    vInitErrorModel();
    
    m_fDevReject = fFindDevRejection(m_fFractionReject);

    nResoGroups = g_nErrorModelBins;
    // Determine if we need less resolution bins.  
    // We get the total # of averages.
    nTotalAvgs = 0;
    

    for(nRefSort=0;nRefSort<m_nNumRefs;nRefSort++)
    {
        nRef=m_pnIndex[nRefSort];
        

        if( m_pfDev[nRefSort] <= m_fDevReject && m_pfMergeSigma[nRefSort] > 0.0 )
        {
            nTotalAvgs++;
        }
    }
    

    // If we break the data up into the given resolution bins, will we still have enough reflections per group?
    if (nTotalAvgs/max(1,nResoGroups) < nMinAvgsPerGroup)
    {
        // No?  Decrease the # of groups.
        nResoGroups = nTotalAvgs/max(1,nMinAvgsPerGroup);
    }
    
    if (nResoGroups<1)
        return 1;
    
    // Determine cutoffs for each resolution bin.
    -afResoGroupMin;
    -afResoGroupMax;
    nx = 0;
    for (nRefSort=0;nRefSort<m_nNumRefs;nRefSort++)
    {
        nRef=m_pnIndexReso[nRefSort];
        ny = m_pnRevIndex[nRef];
        if ((nx == 0) && (afResoGroupMin.size() == afResoGroupMax.size()))
            afResoGroupMin + m_pfRevIndexReso[nRef];
        
        if ((m_pfDev[nRefSort] <= m_fDevReject) && (m_pfMergeSigma[ny] > 0.0))
        {
            if (nx == 0) 
                afResoGroupMax + m_pfRevIndexReso[nRef];

            afResoGroupMax.last() = m_pfRevIndexReso[nRef];
            nx++;
            if ((nx >= nTotalAvgs/max(1,nResoGroups)) && (afResoGroupMin.size()!=nResoGroups))
                nx = 0;
        }
    }
    
    int nRetryCount = -1;
retry:
    nRetryCount++;
    if (1 < nRetryCount) cout << "retry: " << nRetryCount << endl << flush;

    anFirstSubGroupInResoGroup.clear();
    anSubGroupReso.clear();
    anSubGroupCount.clear();
    afSubGroupMaxIntensity.clear();
    afSubGroupAvgIntensity.clear();
    afSubGroupMinIntensity.clear();
    afSubGroupEAdd.clear();
    afSubGroupEMul.clear();
    
    int    nSizeOfDeviationsArray = 0;
   
    // Fit models in each resolution group.
    for (nResoGroup = 0; nResoGroup < nResoGroups; nResoGroup++)
    {
        itr<double> afIntensity;
        itr<int> anStartIndex;
        itr<int> anEndIndex;
        itr<int> anIntensitySort;
        int nAvgsPerSubGroup;
        double fResoRank;
        double fGroupAvg;
        double fGroupCount;

        itr<double> afSigmaSquared;
        itr<double> afIntensityFunc;
        itr<double> afDeviations;
        itr<double> afDeviationSigmas;
               
        anFirstSubGroupInResoGroup + anSubGroupReso.size();
        -afIntensity;
        -anStartIndex;
        -anEndIndex;
        
        nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        nStartIndex = 0;
        for (nRefSort=0;nRefSort<=m_nNumRefs;nRefSort++) 
        {
            nRef = m_pnIndex[nRefSort];

            
            if (nRefSort<m_nNumRefs)
                nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
            
            if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
                continue;
            
            nEndIndex = nRefSort-1;
            

            fResoRank = m_pfRevIndexReso[m_pnIndex[nStartIndex]];
            

            fGroupAvg = 0.0;
            fGroupCount = 0.0;
            for (nIndex = nStartIndex; nIndex <= nEndIndex; nIndex++)
            {
                if (m_pfDev[nIndex] <= m_fDevReject)
                {
                    fGroupAvg += m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity();
                    fGroupCount++;
                }
            }
            

            if (fGroupCount>=2.0)
            {
                fGroupAvg /= fGroupCount;
                
                if ((fResoRank >= afResoGroupMin[nResoGroup]) && (fResoRank < afResoGroupMax[nResoGroup]))
                {
                    anStartIndex + nStartIndex;
                    anEndIndex + nEndIndex;
                    afIntensity + fGroupAvg;
                
                
                
                }
            }
            
            nStartIndex=nRefSort;
            
            if (nRefSort<m_nNumRefs)
            {
                nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
                
            }
        }

        
        
        // Decide number of sub-groups and associated intensity cutoffs defining those sub groups.
        g_pfCmpDoubles = &afIntensity[0];
        -anIntensitySort;
        
        for(nx=0; nx < afIntensity.size(); nx++)
            anIntensitySort + nx;
        
        if( anIntensitySort.size() > 0 )    // safety check
            qsort(&anIntensitySort[0], anIntensitySort.size(), sizeof(int), double_cmp_rel);
        
        nAvgsPerSubGroup = max(1,afIntensity.size()/max(1,nMinAvgsPerSubGroup));
        nAvgsPerSubGroup = afIntensity.size()/max(1,nAvgsPerSubGroup);

        for (ny = 0;ny < afIntensity.size();ny++)
        {
            anSubGroupReso + nResoGroup;
            anSubGroupCount + min(nAvgsPerSubGroup,afIntensity.size() - ny);
            
            if (afIntensity.size() - (anSubGroupCount.last() + ny)< nAvgsPerSubGroup/2)
                anSubGroupCount.last() = afIntensity.size() - ny;
        
            afSubGroupMinIntensity + afIntensity[anIntensitySort[ny]];
            afSubGroupMaxIntensity + afIntensity[anIntensitySort[ny + anSubGroupCount.last() - 1]];
            
            ny += anSubGroupCount.last();
        }


        // Now we fit a model to each sub-group
        while( afSubGroupEAdd.size() < anSubGroupReso.size() )
        {
            int nSubGroup;
            int nAvgGroup;
            double fIntensity;
            double fSigma;
            double fDeviation;
            double fDeviationContrib;
            double fGroupAvg;
            double fGroupAvgAvg;
            double fGroupAvgAvgContrib;
            
            -afSigmaSquared;
            -afIntensityFunc;
            -afDeviations;
            -afDeviationSigmas;
            
            nSubGroup = afSubGroupEAdd.size();
            fGroupAvgAvg = 0.0;
            fGroupAvgAvgContrib = 0.0;
            for (nAvgGroup = 0; nAvgGroup < anStartIndex.size(); nAvgGroup++)
            {
                fGroupAvg = afIntensity[nAvgGroup];
                if ((fGroupAvg >= afSubGroupMinIntensity[nSubGroup]) &&
                    (fGroupAvg <= afSubGroupMaxIntensity[nSubGroup])) {

                    fDeviation = 0.0;
                    fDeviationContrib = 0.0;
                    for (nIndex = anStartIndex[nAvgGroup]; nIndex <= anEndIndex[nAvgGroup]; nIndex++)
                    {
                        if (m_pfDev[nIndex] <= m_fDevReject)
                        {
                            fIntensity = m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity();
                        
                            fSigma     = m_oReflnlist[m_pnIndex[nIndex]].fGetSigmaI();

                            fDeviation += (fIntensity-fGroupAvg)*(fIntensity-fGroupAvg);

                            fDeviationContrib += 1.0;
                        }
                    }
            
                    fDeviation /= max(1.0, fDeviationContrib-1);

                    
                    for (nIndex = anStartIndex[nAvgGroup]; nIndex <= anEndIndex[nAvgGroup]; nIndex++)
                    {
                        fIntensity = m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity();
                        fSigma     = m_oReflnlist[m_pnIndex[nIndex]].fGetSigmaI();
                        
                        afSigmaSquared + fSigma*fSigma;
                        afIntensityFunc + fIntensity*fIntensity;
                        afDeviations + fDeviation; 
                        afDeviationSigmas + (1.0/max(1.0,fDeviationContrib)); /* (fIntensity*fIntensity) */
                        fGroupAvgAvg += fGroupAvg;
                        fGroupAvgAvgContrib += 1.0;                        
                    }                    
                }
            }
            
            afSubGroupAvgIntensity[nSubGroup] = fGroupAvgAvg/max(1.0,fGroupAvgAvgContrib);
            
            double    afArgs[2]       = {0.0};
            
            if( (nSizeOfDeviationsArray = afDeviations.size()) > 0 )
            {
                int       anArgsUsed[2]   = {0};
                double    afArgsSigma[2]  = {0.0};
                double*   ppfFuncs[2] = {NULL};
                double    fTotalVariance = 0.0;
                
                ppfFuncs[0] = &afSigmaSquared[0];
                ppfFuncs[1] = &afIntensityFunc[0];
                
                // Enable all variables.
                anArgsUsed[0] = 1;
                anArgsUsed[1] = 1;
                
                // Disable Emul/Eadd if the user requested this.
                if (m_fUserEMul>0.0)
                {
                    anArgsUsed[0] = 2;
                    afArgs[0] = (m_fUserEMul*m_fUserEMul);
                }
                if (fAvgEMul > 0.0)
                {
                    anArgsUsed[0] = 2;
                    afArgs[0] = fAvgEMul*fAvgEMul;
                }
                if (m_fUserEAdd>0.0)
                {
                    anArgsUsed[1] = 2;
                    afArgs[1] = (m_fUserEAdd*m_fUserEAdd);
                }

                // try solving with both Eadd and Emul active.
                nStat = nSolveLeastSquares(2,
                                        nSizeOfDeviationsArray,
                                        &afDeviations[0],
                                        &afDeviationSigmas[0],
                                        ppfFuncs,
                                        afArgs,
                                        afArgsSigma,
                                        anArgsUsed,
                                        &fTotalVariance);
                
                if( afArgs[0] <= 0.0 || nStat )
                {
                    anArgsUsed[0] = 2;
                    afArgs[0] = 0.0;
                    nSolveLeastSquares(2,
                                     nSizeOfDeviationsArray,
                                     &afDeviations[0],
                                     &afDeviationSigmas[0],
                                     ppfFuncs,
                                     afArgs,
                                     afArgsSigma,
                                     anArgsUsed,
                                     &fTotalVariance);
                } 
                else if (afArgs[1] <= 0.0)
                {
                    anArgsUsed[1] = 2;
                    afArgs[1] = 0.0;
                    nSolveLeastSquares(2,
                                     nSizeOfDeviationsArray,
                                     &afDeviations[0],
                                     &afDeviationSigmas[0],
                                     ppfFuncs,
                                     afArgs,
                                     afArgsSigma,
                                     anArgsUsed,
                                     &fTotalVariance);
                }
            }

            f0 = sqrt(max(0.0,afArgs[0]));
            f1 = sqrt(max(0.0,afArgs[1]));            

            afSubGroupEMul + f0;
            afSubGroupEAdd + f1;
        }
    }

    fAvgEMulNumer = 0.0;
    fAvgEMulDenom = 0.0;

 
    // Next, we apply all of the sigma values.
 
    for (nRef = 0; nRef < m_nNumRefs;nRef++)
    {
        int nSubGroup;
        double fIntensity;

        // Determine resolution group.
        // There really is no easier way to do it, since the resolution ranges are of an indeterminant size.
        for (nResoGroup = 0; nResoGroup < nResoGroups - 1; nResoGroup++)
        {
            if ((m_pfRevIndexReso[nRef] >= afResoGroupMin[nResoGroup]) && (m_pfRevIndexReso[nRef] <= afResoGroupMax[nResoGroup]))
                break;
        }

        fIntensity = m_oReflnlist[nRef].fGetIntensity();
        for (nSubGroup = anFirstSubGroupInResoGroup[nResoGroup];
            (nSubGroup < anSubGroupReso.size() - 1) && (anSubGroupReso[nSubGroup] == nResoGroup); 
            nSubGroup++)
        {
               if ((fIntensity >= afSubGroupMinIntensity[nSubGroup]) &&
                    (fIntensity <= afSubGroupMaxIntensity[nSubGroup])) 
                    break;
               else if ((nSubGroup == anFirstSubGroupInResoGroup[nResoGroup]) && (fIntensity < afSubGroupMinIntensity[nSubGroup])) 
                   break;
        }
        
        // Use the Eadd/Emul for this subgroup.
        m_pfErrorModelEAdd[nRef] = afSubGroupEAdd[nSubGroup];
        m_pfErrorModelEMul[nRef] = afSubGroupEMul[nSubGroup];

        fAvgEMulNumer += afSubGroupEMul[nSubGroup]*fIntensity;
        fAvgEMulDenom += fIntensity;
        
    }

    if ((fAvgEMul <= 0.0) && (m_fUserEMul<=0))
    {
        
        fAvgEMul = fAvgEMulNumer/max(1.0,fAvgEMulDenom);
//+JWP 2010-07-08
// It is dangerous to retry too many times, so do not retry after 4 times
	if  (4 > nRetryCount) goto retry;
//-JWP 2010-07-08	
    }

    // Try to determine the intensity, below which the variance needs be estimated using an average I^2 rather than a per-refleciton I^2
    double fIntensityRank;
    const double fMinIntensityRank = 0.03; 
    const double fMaxIntensityRank = 0.4;
    const double fIntensityRankStep = 0.03;
    double fVariance;
    double fBestVariance;
    double fBestIntensityCutoff;
    double fBestChiSquared;
    double fIntensityCutoff;
    double fIntensityCutoffLower;
    double fChiSquared;
    double fChiSquaredContrib;
    double fObsChiSq = 0.0;

    fBestVariance = 0.0;
    for (fIntensityRank = fMinIntensityRank; fIntensityRank <= fMaxIntensityRank; fIntensityRank += fIntensityRankStep)
    {
        fIntensityCutoff = m_oReflnlist[m_pnIndexIntensity[min(m_nNumRefs-1,(int) (m_nNumRefs*fIntensityRank))]].fGetIntensity();
        fIntensityCutoffLower = m_oReflnlist[m_pnIndexIntensity[min(m_nNumRefs-1,(int) (m_nNumRefs*max(0.0,fIntensityRank-fIntensityRankStep)))]].fGetIntensity();
        m_fErrorModelEAddMinIntensity = fIntensityCutoff;
        fVariance = 0.0;
        fChiSquared = 0.0;
        fChiSquaredContrib = 0.0;
        for (nRefSort = 0; nRefSort < m_nNumRefs;nRefSort++)
        {
            double fIntensity;
            double fSigma;
            nRef = m_pnIndex[nRefSort];
            fIntensity = m_oReflnlist[nRef].fGetIntensity();
            
            
            if ((m_pnGroupSize[nRefSort]>=2) && (m_pfDev2[nRefSort]>0.0))
            {
                fSigma   = fCalcSigma(m_oReflnlist[nRef],nRef);
                f0 = (fSigma - ABS(m_pfDev2[nRefSort]));
                fVariance += f0*f0;                
                f0 = (m_pfDev2[nRefSort]*m_pfDev2[nRefSort]/max(1.0,(fSigma*fSigma)));
                
                if (f0 < m_fDevReject)
                {
                    fChiSquared += f0;
                    fChiSquaredContrib += 1.0;
                }
            }
        }            
        
        if((fVariance < fBestVariance) || (fBestVariance == 0.0) || (fIntensityCutoffLower <= 0.0))
        {
            fBestVariance = fVariance;
            fBestIntensityCutoff = fIntensityCutoff;
            fBestChiSquared = fChiSquared/max(1.0,fChiSquaredContrib);
        }
    }
    
    m_fErrorModelEAddMinIntensity = fBestIntensityCutoff;

    
    // Adjust chi-squared so that it's roughly 1.0 over the entire data set.
    // But DON'T adjust Eadd if it was fixed by user, or Emul if it was fixed by user.
    // Also, don't adjust if the user explicitly requested that we not adjust.
    if ((m_fTargetChiSquare>0.0) && ((m_fUserEAdd<0.0) || (m_fUserEMul<0.0)))
    {
        int nState = 0;
        double fMultiplier;        
        while (!(nStat = nBisection(nState,0.2,5.0,fMultiplier,fObsChiSq,m_fTargetChiSquare,0.03)))
        {
            for (nx=0;nx<m_nNumRefs;nx++)
            {            
                if (m_fUserEAdd<=0.0)
                    m_pfErrorModelEAdd[nx]*= fMultiplier;
                if (m_fUserEMul<=0.0)
                    m_pfErrorModelEMul[nx]*= fMultiplier;
            }
            
            fFindDevRejection(m_fFractionReject,&fObsChiSq);
            
            for (nx=0;nx<m_nNumRefs;nx++)
            {
                if (m_fUserEAdd<=0.0)
                    m_pfErrorModelEAdd[nx]/= fMultiplier; // Should be okay to divide by this.
                if (m_fUserEMul<=0.0)
                    m_pfErrorModelEMul[nx]/= fMultiplier;
            }
        }   
        
        if (nStat == 1)
        {
            // We succeeded in finding a chi^2 near 1.0
            for (nx=0;nx<m_nNumRefs;nx++)
            {
                m_pfErrorModelEAdd[nx]*= fMultiplier;
                m_pfErrorModelEMul[nx]*= fMultiplier;
            }
        }
    }
    
    m_fDevReject = fFindDevRejection(m_fFractionReject,&fObsChiSq);
    
   
    if ((bPrint) || ("" != sGetEnv("DTREK_VERBOSE")) )
    {
        itr<double> m_afResoBins;
        itr<double> m_afIntensityBins;

        itr<double> m_afResoEadd;
        itr<double> m_afResoEaddContrib;
        itr<double> m_afResoEmul;
        itr<double> m_afResoEmulContrib;
        itr<double> m_afResoChiSq;
        itr<int>    m_afResoContrib;
        itr<double> m_afIntensityEadd;
        itr<double> m_afIntensityEaddContrib;
        itr<double> m_afIntensityEmul;
        itr<double> m_afIntensityEmulContrib;
        itr<double> m_afIntensityChiSq;
        itr<int>    m_afIntensityContrib;
        itr<double> m_afIntensityCounts;

        
        float  fResoLow     =  m_a2fRangeReso[0];
        float* pfResoSlope  = &m_fSlopeReso;
        float  fIntensityLow     =  m_a2fRangeInt[1];
        float* pfIntensitySlope  = &m_fSlopeInt;
        int    nBin;
        
        for (nBin = 0; nBin <= m_nNumBinsReso + 1; nBin++)
        {
            m_afResoBins + fResoLow;
            fResoLow = fResoLow + *pfResoSlope;        
            m_afResoEadd[nBin] = 0;
            m_afResoEaddContrib[nBin] = 0;
            m_afResoEmul[nBin] = 0;
            m_afResoEmulContrib[nBin] = 0;
            m_afResoChiSq[nBin] = 0;
            m_afResoContrib[nBin] = 0;
        }
        
        for (nBin = 0; nBin <= m_nNumBinsInt + 1; nBin++)
        {
            m_afIntensityBins + fIntensityLow;
            fIntensityLow = fIntensityLow - *pfIntensitySlope;
            m_afIntensityEadd[nBin] = 0;
            m_afIntensityEaddContrib[nBin] = 0;
            m_afIntensityEmul[nBin] = 0;
            m_afIntensityEmulContrib[nBin] = 0;
            m_afIntensityChiSq[nBin] = 0;
            m_afIntensityContrib[nBin] = 0;
            m_afIntensityCounts[nBin] = 0;
        }
        
        for (nRefSort = 0; nRefSort < m_nNumRefs;nRefSort++)
        {
            double fIntensity;
            double fReso;
            double fSigma;
            double fIoverSig;
            nRef = m_pnIndex[nRefSort];
            fIntensity = m_oReflnlist[nRef].fGetIntensity();
            fReso    = m_oReflnlist[nRef].fGetField(m_nFI_f2STLsq);
            fReso    = fReso * sqrt(fReso);
            fSigma   = fCalcSigma(m_oReflnlist[nRef],nRef);
            fIoverSig = fIntensity/max(1.0,fSigma);
            f0 = m_pfDev2[nRefSort]*m_pfDev2[nRefSort]/max(1.0,fSigma*fSigma);

            if ((m_pnGroupSize[nRefSort]>=2) && (m_pfDev2[nRefSort]>0.0) && (f0 < m_fDevReject))
            {
                for (nBin = 0; nBin < m_nNumBinsInt - 1;nBin++)
                {
                    if (fIoverSig > m_afIntensityBins[nBin])
                        break;
                }
                
                m_afIntensityEadd[nBin] += m_fLastEAddUsed;
                m_afIntensityEaddContrib[nBin] += 1.0;
                m_afIntensityEmul[nBin] += m_fLastEMulUsed;
                m_afIntensityEmulContrib[nBin] += 1.0;
                m_afIntensityChiSq[nBin] += m_pfDev2[nRefSort]*m_pfDev2[nRefSort]/max(1.0,fSigma*fSigma)*m_pnGroupSize[nRefSort]/max(1,(m_pnGroupSize[nRefSort]-1));

                m_afIntensityContrib[nBin] ++;
                m_afIntensityCounts[nBin] += fIntensity;
                
                for (nBin = 0; nBin < m_nNumBinsReso;nBin++)
                {
                    if (fReso < m_afResoBins[nBin+1])
                        break;
                }
                
                if (nBin < m_nNumBinsReso)
                {
                    m_afResoEadd[nBin] += m_fLastEAddUsed;
                    m_afResoEaddContrib[nBin] += 1.0;
                    m_afResoEmul[nBin] += m_fLastEMulUsed;
                    m_afResoEmulContrib[nBin] += 1.0;
                    m_afResoChiSq[nBin] += m_pfDev2[nRefSort]*m_pfDev2[nRefSort]/max(1.0,fSigma*fSigma)*m_pnGroupSize[nRefSort]/max(1,(m_pnGroupSize[nRefSort]-1));
                    m_afResoContrib[nBin] ++;
                }
            }
        }
        
        printf("ChiSq fitting parameters vs Resolution\n");
        printf("-----------------------------------------------------------\n");
        printf(" Resolution Range Contrib  |ChiSq|    Emul   Eadd Eadd/Emul\n");
        printf("-----------------------------------------------------------\n");
        for (nBin = 0; nBin < m_nNumBinsReso;nBin++) {
            f0 = m_afResoEmul[nBin]/max(1.0,m_afResoEmulContrib[nBin]);
            f1 = m_afResoEadd[nBin]/max(1.0,m_afResoEaddContrib[nBin]);
            printf("   %6.2lf -%6.2lf %7d  %7.2lf %7.2lf %6.3lf %9.3lf\n",
                1.0 / max(1e-10,pow(m_afResoBins[nBin], 0.333333)),
                1.0 / max(1e-10,pow(m_afResoBins[nBin + 1], 0.333333)),
                (int) m_afResoContrib[nBin],
                m_afResoChiSq[nBin]/max(1.0,m_afResoContrib[nBin]),
                f0,f1,
                (f0 == 0.0)?1.0:f1/f0
                );
        }
        printf("-----------------------------------------------------------\n\n");
        
        printf("ChiSq fitting parameters vs I/sigI\n");
        printf("------------------------------------------------------------------\n");
        printf(" Int/sigmaI Range Counts Contrib  |ChiSq|    Emul   Eadd Eadd/Emul\n");
        printf("------------------------------------------------------------------\n");
        for (nBin = 0; nBin < m_nNumBinsInt;nBin++) {
            f0 = m_afIntensityEmul[nBin]/max(1.0,m_afIntensityEmulContrib[nBin]); 
            f1 = m_afIntensityEadd[nBin]/max(1.0,m_afIntensityEaddContrib[nBin]);
            printf("   %c%4.0f -  %c%4.0f %6.0lf %7d %8.2lf %7.2lf %6.3lf %9.3lf\n",
                (nBin==0)?'>':' ',                
                m_afIntensityBins[nBin],
                (nBin==m_nNumBinsInt-1)?'<':' ',
                m_afIntensityBins[nBin+1],
                m_afIntensityCounts[nBin]/max(1.0,m_afIntensityContrib[nBin]),
                (int) m_afIntensityContrib[nBin],
                m_afIntensityChiSq[nBin]/max(1.0,m_afIntensityContrib[nBin]),
                f0,f1,
                (f0 == 0.0)?1.0:f1/f0
                );

        }
        printf("------------------------------------------------------------------\n\n");
        
    }

    // Print information about Eadd and Emul.

    
    
    return 0;


}
////////////////////////////////////////////////////////////////////////////////////////
// Mode 0 - remember, 
// Mode 1 - get scale and unapply  
// Mode 2 - apply scale only 
// Mode 3 - unapply scale only
int Cscaleaverage::nRememberIntensities(int nScaleType, int nMode) 
{
    float**     ppfFactorIntensity = NULL;
    float**     ppfFactorSigma = NULL;
    
    int         nRef = 0;
    double      f0 = 0.0;
    double      f1 = 0.0;

    if( nScaleType == 0 )
    {        
        ppfFactorIntensity = &m_pfAbsorbFactors;
        ppfFactorSigma     = &m_pfAbsorbFactorsSig;
    } 
    else
    {
	return 1;
    }
        
    if( nMode == 0 )
    {
        if( *ppfFactorIntensity )
            delete[] *ppfFactorIntensity;
        
        if( *ppfFactorSigma )
            delete[] *ppfFactorSigma;
  
        *ppfFactorIntensity = new float[m_oReflnlist.nGetNumReflns()];
        *ppfFactorSigma     = new float[m_oReflnlist.nGetNumReflns()];
        
        for(nRef = 0; nRef < m_oReflnlist.nGetNumReflns();nRef++)
        {
            (*ppfFactorIntensity)[nRef] = m_oReflnlist[nRef].fGetIntensity();
            (*ppfFactorSigma)[nRef]     = m_oReflnlist[nRef].fGetSigmaI();
        }
    } 
    else if( nMode == 1 )
    {
        for(nRef = 0; nRef < m_oReflnlist.nGetNumReflns(); nRef++)
        {
            f0 = (*ppfFactorIntensity)[nRef];
            f1 = (*ppfFactorSigma)[nRef];
        
            if( f1 > 0.0 )
            {
                (*ppfFactorIntensity)[nRef] = m_oReflnlist[nRef].fGetSigmaI()/f1;
            } 
            else if( f0 > 0.0 )
            {
                (*ppfFactorIntensity)[nRef] = m_oReflnlist[nRef].fGetIntensity()/f0;
            } 
            else
                (*ppfFactorIntensity)[nRef] = 1.0;

            m_oReflnlist[nRef].vSetIntensity((float) f0);
            m_oReflnlist[nRef].vSetSigmaI((float) f1);
        }
    }
    else if (nMode == 2)  // apply scale factors 
    {
        for (nRef = 0; nRef < m_oReflnlist.nGetNumReflns();nRef++)
        {
            m_oReflnlist[nRef].vSetIntensity((float) (m_oReflnlist[nRef].fGetIntensity() * (*ppfFactorIntensity)[nRef]));
            m_oReflnlist[nRef].vSetSigmaI((float) (m_oReflnlist[nRef].fGetSigmaI() * (*ppfFactorIntensity)[nRef]));
        }
    }
    else if (nMode == 3)
    {
        for (nRef = 0; nRef < m_oReflnlist.nGetNumReflns();nRef++)
        {
            m_oReflnlist[nRef].vSetIntensity((float) (m_oReflnlist[nRef].fGetIntensity() / max(1.0,(*ppfFactorIntensity)[nRef])));
            m_oReflnlist[nRef].vSetSigmaI((float) (m_oReflnlist[nRef].fGetSigmaI() / max(1.0,(*ppfFactorIntensity)[nRef])));
        }
    } 
    
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cscaleaverage::nExcludePriorReflections()
{
    if ((!m_poPrevRejectList) || (!m_poPrevRejectList->bIsAvailable()))
        return 1;

    int nRefRej;
    int nPackedHKL;
    int nStart;
    double fRotMid,fRejRotMid;
    int nNumNonExcluded,nNumExcluded;
    
    if (!m_poPrevRejectList->pnGetSortIndex()) {
        if (m_poPrevRejectList->nReduce(m_oCrystal, m_bScaleAnom))
            return 1;
    };

    nNumNonExcluded = 0;
    nNumExcluded = 0;

    for (nRefRej = 0; nRefRej < m_poPrevRejectList->nGetNumReflns(); nRefRej++) {
        nPackedHKL = (*m_poPrevRejectList)[nRefRej].nGetField(m_poPrevRejectList->m_nFI_nPackedHKL);
        fRejRotMid = (*m_poPrevRejectList)[nRefRej].fGetField(m_poPrevRejectList->m_nFI_fCalcRotMid);
        nStart = m_oReflnlist.nFindFirst(m_oReflnlist.m_nFI_nPackedHKL,
                        nPackedHKL,m_pnIndex);

        nNumExcluded = 0;
        while ((nStart < m_oReflnlist.nGetNumReflns()) && (m_oReflnlist[m_pnIndex[nStart]].nGetField(m_oReflnlist.m_nFI_nPackedHKL)==nPackedHKL)) {
            fRotMid = m_oReflnlist[m_pnIndex[nStart]].fGetField(m_oReflnlist.m_nFI_fCalcRotMid);
            if (ABS(fRotMid - fRejRotMid) < 5.0) {
                // Reject.
                m_pnReject[nStart] = 1;
                nNumExcluded++;

            };
            nStart++;
        };
        if (!nNumExcluded)
            nNumNonExcluded++;

    };
    printf("INFO:  %d reflections excluded from list of exclusions.\n",nNumExcluded);
    if (nNumNonExcluded) 
        printf("WARNING:  %d reflections in rejection list not found in data set!\n",nNumNonExcluded);

    return 0;
};

