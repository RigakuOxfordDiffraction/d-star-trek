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
//
// Cscale.cc        Initial author: T.J. Niemeyer           Spring 2000
//  This file contains more member functions of class Cscaleaverage
//    which implements the scaling and merging encapsulation of d*TREK.
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

#include <map>

#include "Cscaleaverage.h"
#include "Cabsorb.h"
#include "CSphericalHarmonic.h"

#ifdef SSI_PC
    #include "CrclHelper.h"
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


double Cscaleaverage::fRmerge(bool bIgnore1,bool bIgnore2) {
    double fNumer;
    double fDenom;
    int nLastHKL,nThisHKL;
    int nStartIndex,nEndIndex;
    int nRef,nRefSort;
    int nx;

    double fGroupAverage;
    int   nGroupSize;
    double f0;
    bool a3bAccept[3];
    int nNumUsed;

    fNumer =0.0;
    fDenom =0.0;

    a3bAccept[0]=TRUE;
    a3bAccept[1]=!bIgnore1;
    a3bAccept[2]=!bIgnore2;



    nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    nStartIndex = 0;
    nNumUsed =0;

    for (nRefSort=0;nRefSort<=m_nNumRefs;nRefSort++)  {
            nRef=m_pnIndex[nRefSort];
                
            if (nRefSort<m_nNumRefs)
                nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
            if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
                continue;
            nEndIndex=nRefSort-1;

            // Calculate the average.
            for (nx=nStartIndex,nGroupSize=0,fGroupAverage=0.0;nx<=nEndIndex;nx++) {
                if (a3bAccept[m_pnReject[nx]]) {
                    nNumUsed++;
                    nGroupSize++;
                    fGroupAverage += m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                };
            };
            if (nGroupSize>1) {
                
                fGroupAverage/=nGroupSize;

                // Add in Rmerge statistics.
                for (nx=nStartIndex;nx<=nEndIndex;nx++) {
                    if (a3bAccept[m_pnReject[nx]]) {
                        f0 = m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                        //printf("%7.2lf comp to avg %7.2lf with %7.2lf\n",f0,fGroupAverage,fabs(f0-fGroupAverage));
                        fNumer+=fabs(f0 - fGroupAverage);
                        fDenom+=fGroupAverage;
                    };
                };
            };
                           
      
            nStartIndex=nRefSort;
            if (nRefSort<m_nNumRefs)
                nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    };
    return fNumer/max(1e-10,fDenom);

};

//shijie yao : 03:26:08 the other merge factors
// Redundance-independent merging R factor [JAC(2001) 34, 130-135].
// Based on the fRmerge() algorithm of above.
double Cscaleaverage::fRmerge_RIM(bool bIgnore1, bool bIgnore2)
{
    double fNumer;
    double fNumberTemp;
    double fDenom;
    double fGroupAverage;
    double f0;

    int nLastHKL, nThisHKL;
    int nStartIndex, nEndIndex;
    int nRef, nRefSort;
    int nx;
    int nGroupSize;
    int nNumUsed;

    bool a3bAccept[3];

    fNumer = 0.0;
    fDenom = 0.0;

    a3bAccept[0]=TRUE;
    a3bAccept[1]=!bIgnore1;
    a3bAccept[2]=!bIgnore2;

    nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    nStartIndex = 0;
    nNumUsed =0;

    for (nRefSort=0;nRefSort<=m_nNumRefs;nRefSort++)  {
        nRef=m_pnIndex[nRefSort];
            
        if (nRefSort<m_nNumRefs)
            nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
            continue;
        nEndIndex=nRefSort-1;

        // Calculate the average.
        for (nx=nStartIndex,nGroupSize=0,fGroupAverage=0.0;nx<=nEndIndex;nx++) {
            if (a3bAccept[m_pnReject[nx]]) {
                nNumUsed++;
                nGroupSize++;
                fGroupAverage += m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
            };
        };
        if (nGroupSize>1) {
            
            fGroupAverage/=nGroupSize;

            // Add in Rmerge statistics.
            fNumberTemp = 0.0;
            for (nx=nStartIndex;nx<=nEndIndex;nx++) {
                if (a3bAccept[m_pnReject[nx]]) {
                    f0 = m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                    fNumberTemp += fabs(f0 - fGroupAverage);
                    fDenom += fGroupAverage;
                };
            };

            //r.i.m algorithm: only difference to Rmerge() of above
            //is to scale the group deviation sum (SUM(I - Iagv) by sqrt(N/N-1)
            fNumberTemp *= sqrt((double)nGroupSize / (double)(nGroupSize - 1));
            fNumer += fNumberTemp;
        };

        nStartIndex=nRefSort;
        if (nRefSort<m_nNumRefs)
            nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    };
    return fNumer/max(1e-10,fDenom);
}; 

//shijie yao : 03:26:08 the other merge factors
// Precision-indicating merging R factor [JAC(2001) 34, 130-135].
// Based on the fRmerge() algorithm of above.
double Cscaleaverage::fRmerge_PIM(bool bIgnore1, bool bIgnore2)
{
    double fNumer;
    double fNumberTemp;
    double fDenom;
    double fGroupAverage;
    double f0;

    int nLastHKL, nThisHKL;
    int nStartIndex, nEndIndex;
    int nRef, nRefSort;
    int nx;
    int nGroupSize;
    int nNumUsed;

    bool a3bAccept[3];

    fNumer = 0.0;
    fDenom = 0.0;

    a3bAccept[0]=TRUE;
    a3bAccept[1]=!bIgnore1;
    a3bAccept[2]=!bIgnore2;

    nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    nStartIndex = 0;
    nNumUsed =0;

    for (nRefSort=0;nRefSort<=m_nNumRefs;nRefSort++)  {
        nRef=m_pnIndex[nRefSort];
            
        if (nRefSort<m_nNumRefs)
            nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
            continue;
        nEndIndex=nRefSort-1;

        // Calculate the average.
        for (nx=nStartIndex,nGroupSize=0,fGroupAverage=0.0;nx<=nEndIndex;nx++) {
            if (a3bAccept[m_pnReject[nx]]) {
                nNumUsed++;
                nGroupSize++;
                fGroupAverage += m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
            };
        };
        if (nGroupSize>1) {
            
            fGroupAverage/=nGroupSize;

            // Add in Rmerge statistics.
            fNumberTemp = 0.0;
            for (nx=nStartIndex;nx<=nEndIndex;nx++) {
                if (a3bAccept[m_pnReject[nx]]) {
                    f0 = m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                    fNumberTemp += fabs(f0 - fGroupAverage);
                    fDenom += fGroupAverage;
                };
            };

            //p.i.m algorithm only difference to Rmerge() of above
            //is to scale the group deviation sum (SUM(I - Iagv) by sqrt(1/N-1)
            fNumberTemp *= sqrt(1.0 / (double)(nGroupSize - 1));
            fNumer += fNumberTemp;
        };

        nStartIndex=nRefSort;
        if (nRefSort<m_nNumRefs)
            nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    };
    return fNumer/max(1e-10,fDenom);
}; 


//shijie yao : 03:26:08 the other merge factors
// Precision-indicating merging R factor [JAC(2001) 34, 130-135].
// Based on the fRmerge() algorithm of above.
double Cscaleaverage::fRmerge_PCV(bool bIgnore1, bool bIgnore2)
{
    double fNumer;
    double fNumberTemp;
    double fDenom;
    double fGroupAverage;
    double f0;

    int nLastHKL, nThisHKL;
    int nStartIndex, nEndIndex;
    int nRef, nRefSort;
    int nx;
    int nGroupSize;
    int nNumUsed;

    bool a3bAccept[3];

    fNumer = 0.0;
    fDenom = 0.0;

    a3bAccept[0]=TRUE;
    a3bAccept[1]=!bIgnore1;
    a3bAccept[2]=!bIgnore2;

    nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    nStartIndex = 0;
    nNumUsed =0;

    for (nRefSort=0;nRefSort<=m_nNumRefs;nRefSort++)  {
        nRef=m_pnIndex[nRefSort];
            
        if (nRefSort<m_nNumRefs)
            nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
            continue;
        nEndIndex=nRefSort-1;

        // Calculate the average.
        for (nx=nStartIndex,nGroupSize=0,fGroupAverage=0.0;nx<=nEndIndex;nx++) {
            if (a3bAccept[m_pnReject[nx]]) {
                nNumUsed++;
                nGroupSize++;
                fGroupAverage += m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
            };
        };
        if (nGroupSize>1) {
            
            fGroupAverage/=nGroupSize;

            // Add in Rmerge statistics.
            fNumberTemp = 0.0;
            for (nx=nStartIndex;nx<=nEndIndex;nx++) {
                if (a3bAccept[m_pnReject[nx]]) {
                    f0 = m_oReflnlist[m_pnIndex[nx]].fGetIntensity();
                    //PCV algorithm
                    fNumberTemp += fabs(f0 - fGroupAverage) * fabs(f0 - fGroupAverage);
                    fDenom += fGroupAverage;
                };
            };

            //PCV algorithm
            fNumberTemp *= (1.0 / (double)(nGroupSize - 1));
            fNumer += sqrt(fNumberTemp);
        };

        nStartIndex=nRefSort;
        if (nRefSort<m_nNumRefs)
            nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    };
    return fNumer/max(1e-10,fDenom);
}; 

int Cscaleaverage::nBatchScale(Cstring& sFixedBatch,bool bPrintResults) {
    
    int nLastHKL,nThisHKL;
    int nRef,nRefSort;
    int nStat;
    int nStartIndex,nEndIndex;
    int nIndex;
    int nx,ny;
    int nIter;
    double f0;
    bool bMatrixInversion = true;

    const int nNumIter = 1;

    itr<double> afSumNumer;
    itr<double> afSumDenom;
    itr<double> afSumVarianceNumer;
    itr<double> afSumVarianceDenom;
    itr<double> afCorrMat;
    itr<double> afCorrMatVariance;
    itr<int>    anBatches;
    itr<int>    anBatchStart;
    itr<int>    anBatchEnd;
    itr<double> afScale;


    //shijie yao : 03:26:08 testing the other R-factor calculation code
    double fInitRMerge;
    double fFinalRMerge;
    double fFinalRMergeNoErrorModel;

    double fInitRMerge_RIM;
    double fFinalRMerge_RIM;
    double fFinalRMergeNoErrorModel_RIM;

    double fInitRMerge_PIM;
    double fFinalRMerge_PIM;
    double fFinalRMergeNoErrorModel_PIM;

    double fInitRMerge_PCV;
    double fFinalRMerge_PCV;
    double fFinalRMergeNoErrorModel_PCV;




    int nBatch;
    int nBatchi1,nBatchi2;
    int nBatch1,nBatch2;


    fInitRMerge = fRmerge();   

    //shijie yao : 03:26:08 testing the other R-factor calculation code
    fInitRMerge_RIM = fRmerge_RIM();    
    fInitRMerge_PIM = fRmerge_PIM();    
    fInitRMerge_PCV = fRmerge_PCV();    

    nStat = 0;

    afSumNumer.setsize(m_nNumBatches*m_nNumBatches);
    afSumDenom.setsize(m_nNumBatches*m_nNumBatches);
    afSumVarianceNumer.setsize(m_nNumBatches*m_nNumBatches);
    afSumVarianceDenom.setsize(m_nNumBatches*m_nNumBatches);
    afCorrMat.setsize(m_nNumBatches*m_nNumBatches);
    afCorrMatVariance.setsize(m_nNumBatches*m_nNumBatches);
    afScale.setsize(m_nNumBatches);
    memset(&afSumNumer[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
    memset(&afSumDenom[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
    memset(&afSumVarianceNumer[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
    memset(&afSumVarianceDenom[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
    memset(&afCorrMat[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
    memset(&afCorrMatVariance[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
    memset(&afScale[0],0,sizeof(double)*m_nNumBatches);

    for (nBatch = 0; nBatch < m_nNumBatches; nBatch++)
        afScale[nBatch] = 1.0;
    
    for(nIter = 0; nIter < nNumIter && !nStat; nIter++)
    {
        nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        
        nStartIndex = 0;
        int     nLinearArrayIndex = 0;
        for(nRefSort=0; nRefSort <= m_nNumRefs; nRefSort++)
        {
            nRef = m_pnIndex[nRefSort];
            
            if( nRefSort < m_nNumRefs)
                nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
            
            if( nRefSort != m_nNumRefs && nThisHKL == nLastHKL )
                continue;
            
            nEndIndex = nRefSort - 1;
            
            -anBatches;
            -anBatchStart;
            -anBatchEnd;
            for (nIndex = nStartIndex; nIndex <= nEndIndex; nIndex++)
            {
                nx = m_oReflnlist[m_pnIndex[nIndex]].nGetField(m_nFI_nBatchIdx);
                if( anBatches.size() && nx < anBatches.last() )
                {
                    printf("ERROR:  In batch scaling.\n");
                    return 1;
                }
                
                if( !anBatches.size() || nx > anBatches.last() )
                {
                    anBatches + nx;
                    anBatchStart + nIndex;
                    anBatchEnd + nIndex;
                } 
                else
                    anBatchEnd.last() = nIndex;
            }



            for(nBatchi1 = 0; nBatchi1 < anBatches.size(); nBatchi1++)
            {                
                for(nBatchi2 = 0; nBatchi2 < anBatches.size(); nBatchi2++)
                {
                    if( nBatchi1 == nBatchi2 )
                        continue;

                    double  fBatch1Intensity = 0.0;
                    double  fBatch2Intensity = 0.0;
                    double  fBatch1Variance  = 0.0;
                    double  fBatch2Variance  = 0.0;
                    int     nCount1          = 0;
                    int     nCount2          = 0;
                    double  fIntensity       = 0.0;
                    double  fSigma           = 0.0;
                    
                    nBatch1 = anBatches[nBatchi1];
                    nBatch2 = anBatches[nBatchi2];
                    
                    // Average and get sigma of all elements in each batch.
                    for (nIndex = anBatchStart[nBatchi1]; nIndex <= anBatchEnd[nBatchi1];nIndex++)
                    {
                        if (m_pnReject[nIndex] == 0)
                        {
                            fIntensity = m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity() * afScale[nBatch1];
                            
                            fSigma = fCalcSigma(m_oReflnlist[m_pnIndex[nIndex]],m_pnIndex[nIndex]) * afScale[nBatch1];
                            
                            fBatch1Intensity += fIntensity/max(1.0,(fSigma*fSigma));
                            
                            fBatch1Variance += 1.0/max(1.0,(fSigma*fSigma));
                            
                            nCount1 ++;
                        }
                    }
                    
                    if( !nCount1 )
                        continue;

                    //fBatch1Intensity /= fBatch1Variance;
                    fBatch1Intensity /= max(1e-20,fBatch1Variance);

                    fBatch1Variance = 1.0/max(1e-20,fBatch1Variance);
                    
                    
                    for (nIndex = anBatchStart[nBatchi2]; nIndex <= anBatchEnd[nBatchi2];nIndex++)
                    {
                        if (m_pnReject[nIndex] == 0 )
                        {
                            fIntensity = m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity()*afScale[nBatch2];
                            fSigma = fCalcSigma(m_oReflnlist[m_pnIndex[nIndex]],m_pnIndex[nIndex])*afScale[nBatch2];
                            fBatch2Intensity += fIntensity/max(1.0,(fSigma*fSigma));
                            fBatch2Variance += 1.0/max(1.0,(fSigma*fSigma));
                            nCount2 ++;
                        };
                    };
                    if (!nCount2)
                        continue;
                    fBatch2Intensity /= max(1e-20,fBatch2Variance);
                    fBatch2Variance = 1.0/max(1e-20,fBatch2Variance);
                    
                    
                    f0 = 1.0 / ( fBatch1Variance + fBatch2Variance );
                    
                    nLinearArrayIndex = nBatch1 * m_nNumBatches + nBatch2;
                    
                    afSumNumer[nLinearArrayIndex]          += fBatch1Intensity * fBatch2Intensity * f0;
                    afSumDenom[nLinearArrayIndex]          += fBatch1Intensity * fBatch1Intensity * f0;
                    
                    afSumVarianceNumer[nLinearArrayIndex]  += (fBatch1Variance * fBatch2Intensity * fBatch2Intensity + 
                                                            fBatch2Variance * fBatch1Intensity * fBatch1Intensity) * 
                                                            f0 * f0;

                    afSumVarianceDenom[nLinearArrayIndex]  += 2.0 * fBatch1Variance * fBatch1Intensity * fBatch1Intensity * f0 * f0;
                }
            }
            
            nStartIndex=nRefSort;
            if (nRefSort<m_nNumRefs)
                nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
        // RB restrain groups of batches
        if( 0 != m_vecScaleRestrInfo.size() )
        {
            if( 1 == m_vecScaleRestrInfo.size() && "-1" == m_vecScaleRestrInfo[0].m_sBatchBegin
                                                && "-1" == m_vecScaleRestrInfo[0].m_sBatchEnd )  // see if restr info implies some defaults to be set
                vGenerateDefaultBatchRestrainInfo(m_vecScaleRestrInfo[0].m_dWeight);                
            
            vConvertBatchRestrainInfoToBatchIndices();

            int         nBatchBegin = 0;
            int         nBatchEnd = 0;
            double      dRestrWeight = 0.0;
            double      dRestrain_coefficient = 0.0;

            for(int ii=0; ii < m_vecScaleRestrInfo.size(); ii++)
		    {
			    nBatchBegin             = atoi(m_vecScaleRestrInfo[ii].m_sBatchBegin.string());
			    nBatchEnd               = atoi(m_vecScaleRestrInfo[ii].m_sBatchEnd.string());
                dRestrWeight            = m_vecScaleRestrInfo[ii].m_dWeight;
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // OUTPUT AND ERROR CHECKING
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if( nBatchEnd   <  0             || nBatchBegin <  0             ||
                    nBatchEnd   >= m_nNumBatches || nBatchBegin >= m_nNumBatches || 
                    nBatchBegin >= nBatchEnd )
                {
                    cout << "\nCannot restrain batch scalefactors. Invalid batch index range: " << nBatchBegin << " " << nBatchEnd << endl << flush;
                    cout << "This may be due to wrong input or some batches have been already rejected.\n";
                    continue;
                }
                
                cout << "\nRestraining batch scalefactors. Batch Names(Indices): " << m_ptBatchInfo[nBatchBegin].m_sBatchName.string() 
                                                                                   << "(" << nBatchBegin << ")" 
                                                                                   << " - " 
                                                                                   << m_ptBatchInfo[nBatchEnd].m_sBatchName.string() 
                                                                                   << "(" << nBatchEnd << ")"                                                                                                           
                                                                                   << " Weight: " << dRestrWeight << endl << flush;

                if( dRestrWeight >= 1.0 || dRestrWeight <= 0.0 )
                {
                    cout << "\nCannot restrain batch scalefactors. Invalid restraining weight: " << dRestrWeight << endl << flush;
                    continue;
                }
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                dRestrain_coefficient   = 0.5 * (1.0 / dRestrWeight) * (1.0 / dRestrWeight);

			    for(int jj=nBatchBegin; jj <= nBatchEnd-1; jj++)		// go over all batches that have a right-hand neighbor
			    {
				    for(nBatch1 = jj; nBatch1 <= jj+1; nBatch1++)	  // defines the two batches that have a common "dummy reflection"
				    {                
					    for(nBatch2 = jj; nBatch2 <= jj+1; nBatch2++)
                        {
						    if( nBatch1 == nBatch2 )
							    continue;

                            nLinearArrayIndex = nBatch1 * m_nNumBatches + nBatch2;
                    
						    afSumNumer[nLinearArrayIndex]          += dRestrain_coefficient;
						    afSumDenom[nLinearArrayIndex]          += dRestrain_coefficient;
                    
						    afSumVarianceNumer[nLinearArrayIndex]  += dRestrain_coefficient;
						    afSumVarianceDenom[nLinearArrayIndex]  += dRestrain_coefficient;

                        }  // summation batch 2
				    } // summation batch 1
			    } // all batches that have a left and right neighbor to be restrained
		    } // all "runs" to be restrained
        }////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Now calculate correlation constants K[i,j] where i < j and K[i,j]*Batch(i) = Batch(j) 
        // For Matrix X we assume X[i,j] = X[i*m_nNumBatches + j]
        
        double      fNumer = 0.0;
        double      fDenom = 0.0;
        double      fVarNumer = 0.0;
        double      fVarDenom = 0.0;

        for (nBatch1 = 0; nBatch1 < m_nNumBatches; nBatch1++)
        {
            for (nBatch2 = 0; nBatch2 < m_nNumBatches; nBatch2++)
            {
                nLinearArrayIndex = nBatch1*m_nNumBatches + nBatch2;
                
                fDenom = afSumDenom[nLinearArrayIndex];
                
                if( fDenom )
                {
                    fNumer = afSumNumer[nLinearArrayIndex];
                    
                    afCorrMat[nLinearArrayIndex] = fNumer /max(1e-20, fDenom);
                    
                    fVarNumer = afSumVarianceNumer[nLinearArrayIndex];
                    fVarDenom = afSumVarianceDenom[nLinearArrayIndex];
                    
                    afCorrMatVariance[nLinearArrayIndex] = fVarNumer / max(1e-20, fDenom * fDenom) + 
                                                          fVarDenom * fNumer * fNumer / max(1e-20, fDenom * fDenom * fDenom * fDenom);
                }
            }
        }
        
        // Calculate the final matrix describing equations using scale factors.
        
        itr<double>     afAMat;
        itr<double>     afAMatInv;
        itr<double>     afXVec;
        itr<double>     afEigenVecs;
        itr<double>     afEigenValues;
        double          fAvgAbsSum;
        
        afAMat.setsize(m_nNumBatches*m_nNumBatches);
        afAMatInv.setsize(m_nNumBatches*m_nNumBatches);
        afEigenVecs.setsize(m_nNumBatches*m_nNumBatches);
        afXVec.setsize(m_nNumBatches);
        afEigenValues.setsize(m_nNumBatches);
        memset(&afAMat[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
        memset(&afAMatInv[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
        memset(&afEigenVecs[0],0,sizeof(double)*m_nNumBatches*m_nNumBatches);
        memset(&afXVec[0],0,sizeof(double)*m_nNumBatches);
        memset(&afEigenValues[0],0,sizeof(double)*m_nNumBatches);
        
        if( !bMatrixInversion )
        {
            // Iteratively construct the scale factor values (non-matrix inversion method).
            for( nBatch1 = 0; nBatch1 < m_nNumBatches; nBatch1++ )
            {
                double fScaleNumer = 0.0;
                double fScaleDenom = 0.0;
                double fScale = 0.0;
                for(nBatch2 = 0; nBatch2 < m_nNumBatches; nBatch2++)
                {
                    nLinearArrayIndex = nBatch2*m_nNumBatches + nBatch1;
                    
                    if( nBatch2 != nBatch1 && afCorrMat[nLinearArrayIndex] )
                    {
                        fScaleNumer += afCorrMat[nLinearArrayIndex] / max(1e-20, afCorrMatVariance[nLinearArrayIndex]);
                        fScaleDenom += 1.0/max(1e-20, afCorrMatVariance[nLinearArrayIndex] );
                    } 
                }
                
                fScale = fScaleNumer / max(1e-20, fScaleDenom);
                
                afScale[nBatch1] *= 1.0 / fScale;
            }
        } 
        else
        { 
            fAvgAbsSum = 0.0;
            for (nBatch1 = 0; nBatch1 < m_nNumBatches; nBatch1++)
            {  
                double      fSumInvVariance = 0.0;
                // We will fill the equation describing Batch1's scale factor.
                for (nBatch2 = 0; nBatch2 < m_nNumBatches; nBatch2++)
                {
                    nLinearArrayIndex = nBatch2*m_nNumBatches + nBatch1;

                    if (nBatch1 != nBatch2)
                    {
                        // afCorrMat[nBatch2,nBatch1] contains the conversion from a batch2 scale factor to a batch1 scale factor.
                        if (afCorrMatVariance[nLinearArrayIndex] > 0.0)
                        {
                            fSumInvVariance += 1.0 / afCorrMatVariance[nLinearArrayIndex];
                            afAMat[nLinearArrayIndex] = afCorrMat[nLinearArrayIndex] / afCorrMatVariance[nLinearArrayIndex];
                            fAvgAbsSum += ABS(afAMat[nLinearArrayIndex]);
                        }
                    }
                }
                
                // All of these average together to get the batch1 scale factor.  This is the diagonal element of the matrix.
                // However, if this is the fixed batch, then we fill in the diagonal element with '0' and but the constant in the B-vector
                afAMat[nBatch1 * m_nNumBatches + nBatch1] = -fSumInvVariance;
                fAvgAbsSum += ABS(afAMat[nBatch1 * m_nNumBatches + nBatch1]);
            }
            
            // Invert the matrix.
            if( 0 != nInvMatND_svd(m_nNumBatches, &afAMat[0], &afAMatInv[0], &afEigenVecs[0], &afEigenValues[0]) )
                nStat = 1;
            
            if( !nStat )
            {
                int         nMinEigenValue = 0;
                double      fMinEigenValue = 0.0;
                
                // Look through the eigen values and find one that is nearest to zero.
                for(nx = 0; nx < m_nNumBatches; nx++)
                {
                    if( nx == 0  || fMinEigenValue > afEigenValues[nx] )
                    {
                        fMinEigenValue = afEigenValues[nx];
                        nMinEigenValue = nx;
                    }
                }
                
                
                double      fAvgScale = 0.0;
                for(nx = 0; nx < m_nNumBatches; nx++)
                {
                    afXVec[nx] = afEigenVecs[nMinEigenValue * m_nNumBatches + nx];
                    fAvgScale += afXVec[nx];
                }
                
                fAvgScale /= max(1, m_nNumBatches);
                
                for (nx = 0; nx < m_nNumBatches; nx++)
                {
                    afScale[nx] *= (afXVec[nx]==0.0) ? 0.0 : ( fAvgScale / afXVec[nx] );
                    
                    if( afScale[nx] <= 0.0 || afXVec[nx]==0.0 )
                    {
                        printf("ERROR:  Batch scale factors look very unreasonable.  They will not be used.\n");
                        nStat = 1;
                    }
                }
            }
        }
    }
    

    if( sFixedBatch.length() )
    {
        int     nFixedBatch = 0;;
        // Adjust the fixed batch so that it is equal to 1.0
        for(nx=0; nx < m_nNumBatches; nx++)
        {
            if( m_ptBatchInfo[nx].m_sBatchName == sFixedBatch )
            {
                nFixedBatch = nx;
                break;
            }
        }
        
        if( nx == m_nNumBatches )
        {
            printf("WARNING:  Could not find requested fixed batch '%s'\n", sFixedBatch.string());
        }
        
        f0 = afScale[nFixedBatch];
        for(ny=0; ny < m_nNumBatches; ny++)
        {
            if( ny == nFixedBatch )
                afScale[ny] = 1.0;  // Make it exactly equal to 1.0
            else
                afScale[ny] /= max(1e-10, f0);
        }
    }
    
    if( bPrintResults )
    {
        int nFirstBatch;
        int nBatchGroup = 10;
        char *pcShortLine = "--------------------------------------------------------------------------------\n";
        printf("\nScale factors for %d batches  (variables in brackets were fixed.)\n",
            m_nNumBatches);
        printf(pcShortLine);
        printf(" Batch id    ...     ...     ...     ...     ...     ...     ...     ...     ...\n");
        printf(" Scale factor...     ...     ...     ...     ...     ...     ...     ...     ...\n");
        printf(pcShortLine);
        for (nFirstBatch=0;nFirstBatch<m_nNumBatches;nFirstBatch+=nBatchGroup)
        {
            if (0 < nFirstBatch)
                printf("\n");
            //            printf(" Batch number ");
            for (nx=nFirstBatch;nx<min(m_nNumBatches,nFirstBatch+nBatchGroup);nx++) 
                printf("%8s",m_ptBatchInfo[nx].m_sBatchName.string());
            printf("\n");
            //            printf(" Batch scale  ");
            for (nx=nFirstBatch;nx<min(m_nNumBatches,nFirstBatch+nBatchGroup);nx++) {
                bool bFixed;
                bFixed =  (afScale[nx] == 1.0);
                f0 = afScale[nx];
                if (bFixed)
                    printf(" [%5.3lf]",m_ptBatchInfo[nx].m_fBatchScale * f0);
                else
                    printf(" %7.3lf",m_ptBatchInfo[nx].m_fBatchScale * f0);
            };
            printf("\n");
        };    
        printf(pcShortLine); 
    };

   // Multiply the reflections by their new values.
    
    if (!nStat)
    {
        for (nRef=0;nRef<m_nNumRefs;nRef++)
        {
            nx = m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx);
            f0 = afScale[nx];
            m_oReflnlist[nRef].vSetIntensity(m_oReflnlist[nRef].fGetIntensity()*f0);
            m_oReflnlist[nRef].vSetSigmaI(m_oReflnlist[nRef].fGetSigmaI()*f0);
        }
        
        fFinalRMergeNoErrorModel = fRmerge();

        //shijie yao : 03:26:08 testing the other R-factor calculation code
        fFinalRMergeNoErrorModel_RIM = fRmerge_RIM();    
        fFinalRMergeNoErrorModel_PIM = fRmerge_PIM();    
        fFinalRMergeNoErrorModel_PCV = fRmerge_PCV();    

        nFitChiSquared(false);
        fCalcChiSquared(1);
        
        fFinalRMerge = fRmerge();

        //shijie yao : 03:26:08 testing the other R-factor calculation code
        fFinalRMerge_RIM = fRmerge_RIM();    
        fFinalRMerge_PIM = fRmerge_PIM();    
        fFinalRMerge_PCV = fRmerge_PCV();    
    } 
    else
    {
        fFinalRMergeNoErrorModel = 1.0;
        fFinalRMerge = 1.0;
        fFinalRMerge_RIM = 1.0; //mrp    
        fFinalRMergeNoErrorModel_RIM = 1.0; //mrp    
        fFinalRMerge_PCV = 1.0; //mrp    
        fFinalRMergeNoErrorModel_PCV = 1.0; //mrp    
    }
    

    printf("Batch scale results\n");
    printf("--------------------------------------------------------------------------\n");
    printf("\nRmerge before applying scale factors:                         %.4lf\n", fInitRMerge);
    printf(  "Rmerge after  applying scale factors with    adjusted sigmas: %.4lf\n", fFinalRMerge);
    printf(  "Rmerge after  applying scale factors without adjusted sigmas: %.4lf\n\n", fFinalRMergeNoErrorModel);


    //shijie yao : 03:26:08 testing the other R-factor calculation code
    printf("\nRmeas before applying scale factors:                          %.4lf\n", fInitRMerge_RIM);
    printf(  "Rmeas after  applying scale factors with    adjusted sigmas:  %.4lf\n", fFinalRMerge_RIM);
    printf(  "Rmeas after  applying scale factors without adjusted sigmas:  %.4lf\n\n", fFinalRMergeNoErrorModel_RIM);
    /****+JWP 2008-06-13
    printf("\nRpim before applying scale factors:                           %.4lf\n", fInitRMerge_PIM);
    printf(  "Rpim after  applying scale factors with    adjusted sigmas:   %.4lf\n", fFinalRMerge_PIM);
    printf(  "Rpim after  applying scale factors without adjusted sigmas:   %.4lf\n\n", fFinalRMergeNoErrorModel_PIM);
    -JWP 2008-06-13****/

    printf("\nRpcv before applying scale factors:                           %.4lf\n", fInitRMerge_PCV);
    printf(  "Rpcv after  applying scale factors with    adjusted sigmas:   %.4lf\n", fFinalRMerge_PCV);
    printf(  "Rpcv after  applying scale factors without adjusted sigmas:   %.4lf\n\n", fFinalRMergeNoErrorModel_PCV);


    printf("(All merging R factors are computed for reflections having I/sig(I) >= %.1lf)\n",m_fSigma);
    printf("--------------------------------------------------------------------------\n\n");

	if( fFinalRMerge > fInitRMerge || fFinalRMerge <= 0.0 || 0 != nStat )
    {
		printf("Batch scale failed to decrease R-merge.  Scale factors will not be applied.\n");

        // Un-Multiply the reflections.
        if( !nStat )
        {
            for (nRef=0;nRef<m_nNumRefs;nRef++)
            {
                nx = m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx);
                f0 = max(1e-10,afScale[nx]);
                m_oReflnlist[nRef].vSetIntensity(m_oReflnlist[nRef].fGetIntensity()/f0);
                m_oReflnlist[nRef].vSetSigmaI(m_oReflnlist[nRef].fGetSigmaI()/f0);
            }
            
            nFitChiSquared(false);
            
            fCalcChiSquared(1);
        }        		
    } 
    else
    {
		// Copy the batch scale factors over.
		for (nx=0;nx<m_nNumBatches;nx++)
        {
			m_ptBatchInfo[nx].m_fBatchScale *= afScale[nx];
		}
    }

    fflush(stdout);


    return 0;
}


int Cscaleaverage::nBFactorScale(Cstring& sFixedBatch,bool bPrintResults)
{
    int         nStat = 0;

    int         nBatchIndex=0;
    int         nRefSort=0;
    int         nRefSort2=0;
    int         nRef=0;

    itr<double>         afValues;
    itr<double>         afValuesSigma;
    itr<double>         afBFactorFunc;

    if( !m_nNumRefs )
        return 1;
    

    double    fInitRMerge = fRmerge();    

    int         nFirstRefSortInBatch = 0;
    int         nLastRefSortInBatch=0;

    for(nRefSort = 0; nRefSort <= m_nNumRefs; nRefSort++)
    {
        nRef = nRefSort < m_nNumRefs ? m_pnIndexBatch[nRefSort] : 0;

        if( nRefSort != 0 && ( nRefSort == m_nNumRefs || m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx) != nBatchIndex ) ) 
        {
            nLastRefSortInBatch = nRefSort - 1;

            double      fIntensity = 0.0;
            double      fBFactorFunc = 0.0;
            double*     apfFuncs[2]={NULL};
            double      afArgs[2]={0.0};
            double      afArgsSigma[2]={0.0};

            -afValues;
            -afValuesSigma;
            -afBFactorFunc;
            
            for(nRefSort2 = nFirstRefSortInBatch; nRefSort2 <= nLastRefSortInBatch; nRefSort2++) 
            {
                nRef = m_pnIndexBatch[nRefSort2];
                
                fIntensity = m_oReflnlist[nRef].fGetIntensity();

                if( fIntensity <= 0.0 )
                    continue; // to prevent erroneous measurements from entering refinement

                fBFactorFunc = m_oReflnlist[nRef].fGetField(m_nFI_f2STLsq);
                fBFactorFunc = fBFactorFunc / 2.0;

                if( fIntensity > 0.0 )
                {
                    afValues + log(fIntensity);
                    afValuesSigma + 1.0;
                    afBFactorFunc + fBFactorFunc;
                }
            }

            apfFuncs[0] = &afBFactorFunc[0];
            apfFuncs[1] = ms_pfConst;

            if( 0 == afValues.size() )
                return 1; // safety check

            nStat = nSolveLeastSquares(2, 
                                       afValues.size(),
                                       &afValues[0],
                                       &afValuesSigma[0],
                                       &apfFuncs[0],
                                       &afArgs[0],
                                       &afArgsSigma[0],
                                       NULL);

            if( 0 != nStat )
            {
                return nStat;
            }

            printf("%5d %6s %.2lf %.2lf\n", nBatchIndex+1, m_ptBatchInfo[nBatchIndex].m_sBatchName.string(), afArgs[0], afArgs[1] );
            fflush(stdout);

            ///////////////////////////////////////////////////////
            // RB save B-factor in the m_ptBatchInfo array
            m_ptBatchInfo[nBatchIndex].m_fBFactorScale = afArgs[0]; 
            m_ptBatchInfo[nBatchIndex].m_fBFactorConst = afArgs[1]; 
            ///////////////////////////////////////////////////////

#if 0
            FILE*           pFOut = fopen("wilson.txt","wt");
            
            for(int nPass = 0; nPass < 2; nPass++)
            {
                for(nRefSort2 = nFirstRefSortInBatch; nRefSort2 <= nLastRefSortInBatch; nRefSort2++)
                {
                    nRef = m_pnIndexBatch[nRefSort2];
                    fIntensity = m_oReflnlist[nRef].fGetIntensity();

                    fBFactorFunc = m_oReflnlist[nRef].fGetField(m_nFI_f2STLsq);
                    
                    if (nPass == 0)
                        fprintf(pFOut,"%.4lf %.2lf\n", fBFactorFunc / 2.0, exp(afArgs[1] + fBFactorFunc*afArgs[0]));
                    else
                        fprintf(pFOut,"%.4lf %.2lf\n", fBFactorFunc / 2.0, fIntensity);
                }
                
                fprintf(pFOut,"\n");
                fflush(pFOut);
            }
            
            fclose(pFOut);
#endif    

            nFirstRefSortInBatch = nRefSort;
        }  
        
        if( nRefSort < m_nNumRefs ) 
            nBatchIndex = m_oReflnlist[m_pnIndexBatch[nRefSort]].nGetField(m_nFI_nBatchIdx);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Adjust the fixed batch so that it's B is equal to 0.0
    int     nx = 0;
    if( sFixedBatch.length() )
    {
        for(nx=0; nx < m_nNumBatches; nx++)
        {
            if( m_ptBatchInfo[nx].m_sBatchName == sFixedBatch )
                break;
        }
        
        if( nx == m_nNumBatches )
        {
            printf("WARNING:  Could not find requested fixed batch '%s'\n", sFixedBatch.string());
            nx = 0;
        }
        
        double      f0 = m_ptBatchInfo[nx].m_fBFactorScale;
        double      f1 = m_ptBatchInfo[nx].m_fBFactorConst;
        
        for(int ny=0; ny < m_nNumBatches; ny++)
        {
            m_ptBatchInfo[ny].m_fBFactorScale -= f0;
            
            m_ptBatchInfo[ny].m_fBFactorConst -= f1;
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////
    // Map real batch indices to the indices in the m_ptBatchInfo array
    std::map<Cstring, int>          mapBatchIndices;
    
    for(nx=0; nx < m_nNumBatches; nx++)
    {
        mapBatchIndices.insert(std::pair<Cstring,int>(m_ptBatchInfo[nx].m_sBatchName, nx));
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    
    std::map<Cstring, int>::iterator        oIt;
    
    double    dBFactor = 0.0;

    double    dSTL=0.0;

    Cstring     sBatchName("");

    for(nRef=0; nRef < m_nNumRefs; nRef++)
    {
        sBatchName = m_oReflnlist[nRef].sGetField(m_oReflnlist.m_nFI_sBatch); // get this reflection's batch name

        oIt = mapBatchIndices.find(sBatchName);                  // find this batch's index in m_ptBatchInfo array
        
        if( oIt != mapBatchIndices.end() )
        {
            nBatchIndex = (*oIt).second;
        }
        else
        {
            cout << "\nCould not find batch with ID " << sBatchName << "in the batch list\n";
            return 1;
        }

        dSTL = m_oReflnlist[nRef].fGetField(m_nFI_f2STLsq);
        dSTL /= 2.0;

        dBFactor = exp( (-1.0) * (m_ptBatchInfo[nBatchIndex].m_fBFactorConst + m_ptBatchInfo[nBatchIndex].m_fBFactorScale * dSTL) );
        
        m_oReflnlist[nRef].vSetIntensity(m_oReflnlist[nRef].fGetIntensity()*dBFactor);
        m_oReflnlist[nRef].vSetSigmaI(m_oReflnlist[nRef].fGetSigmaI()*dBFactor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Now check to see if we are doing better on the R-factors...
    double      fFinalRMergeNoErrorModel = fRmerge();
    
    nFitChiSquared(false);
    
    fCalcChiSquared(1);
        
    double      fFinalRMerge = fRmerge();
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("B-factor scale results\n");
    printf("--------------------------------------------------------------------------\n");
    printf("\nRmerge before applying B factors:                             %.4lf\n", fInitRMerge);
    printf(  "Rmerge after  applying B factors with    adjusted sigmas:     %.4lf\n", fFinalRMerge);
    printf(  "Rmerge after  applying B factors without adjusted sigmas:   %.4lf\n\n", fFinalRMergeNoErrorModel);
    printf("(All R-merge values are computed for reflections having I/sig(I) >= %.1lf)\n",m_fSigma);
    printf("--------------------------------------------------------------------------\n\n");

	if( fFinalRMerge > fInitRMerge || fFinalRMerge <= 0.0 || nStat )
    {
		printf("B-factor scale failed to decrease R-merge. B factors will not be applied.\n\n");

        // Un-Multiply the reflections.
        for(nRef=0; nRef < m_nNumRefs; nRef++)
        {
            sBatchName = m_oReflnlist[nRef].sGetField(m_oReflnlist.m_nFI_sBatch); // get this reflection's batch name

            oIt = mapBatchIndices.find(sBatchName);                  // find this batch's index in m_ptBatchInfo array
        
            if( oIt != mapBatchIndices.end() )
            {
                nBatchIndex = (*oIt).second;
            }
            else
            {
                cout << "\nCould not find batch with ID " << nx << "in the batch list\n";
                return 1;
            }

            dSTL = m_oReflnlist[nRef].fGetField(m_nFI_f2STLsq);
            dSTL /= 2.0;

            dBFactor = exp(m_ptBatchInfo[nBatchIndex].m_fBFactorConst + m_ptBatchInfo[nBatchIndex].m_fBFactorScale * dSTL);
        
            m_oReflnlist[nRef].vSetIntensity(m_oReflnlist[nRef].fGetIntensity() * dBFactor);
            m_oReflnlist[nRef].vSetSigmaI(m_oReflnlist[nRef].fGetSigmaI() * dBFactor);
        }
        
        /////////////////////////////////////////////        
        // Erase information from the Batch info
        for(int ny=0; ny < m_nNumBatches; ny++)
        {
            m_ptBatchInfo[ny].m_fBFactorScale = 0.0;
            m_ptBatchInfo[ny].m_fBFactorConst = 0.0;
        }
        /////////////////////////////////////////////
        
        nFitChiSquared(false);
        
        fCalcChiSquared(1);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    fflush(stdout);

    return 0;
}



/*  This routine will attempt to uniformly rebatch the data so that specified overlap
    requirements are achieved.  A uniform batching consists of a partition of the
    rotation range so that each batch gets the same rotation width.

    The routine will assume a format of incoming batch strings:
    
     a)  Either the batch string is 4 digits, OR
     b)  The batch string has >4 characters.
      
    For each unique scan, the routine will sort all reflections in the scan based on increasing
    rotation angle.  It will then attempt to uniformly batch the reflections so that the input
    requirements are met.  It's first try will utilize the full number of input batches.  It 
    proceeds to decrease the scan interval until either the requirements are met, or there is 
    a single batch left.  

*/

int Cscaleaverage::nReBatch(int nMinBatchOverlaps,
                            bool bRebatchImageBoundry,
                            bool bRebatchBookends,
                            bool bRebatch)
{
    int nx;
    double f0,f1;
    int nStat;
    const int nMaxBatchTest = 50;  // Maximum rebatching to be used.
    Cstring sTemp, sTemp2;
    
    Crefln* poRefln;
    int nRef;
    
    Cstring sScanID;
    int nFirstSortRef;
    int nLastSortRef;
    int nLastBatch;
    int nBatches;
    int nScans;
    float fImageWidth;
    /*  Allowable input flags:


    bRebatchImageBoundry=false bRebatchBookends=false bRebatch=true/false
    bRebatchImageBoundry=true bRebatchBookends=false bRebatch=true/false
    bRebatchImageBoundry=false bRebatchBookends=true bRebatch=true/false

    */
    
    int*    pnOverlaps;         // Used in the inside loop below.
    float*  pfBatchStart;       // Starting rotation for batch.
    float*  pfBatchEnd;         // Ending rotation for batch.
    int*    pnBatchI;           // Batch index (w.r.t. overlap table).
    bool*   pbInScan;           // When processing a scan, this indicates which reflections are in the scan.
    float*  pfSuggestedRebatch; // Suggested Rebatch interval length.
    float*  pfPrevBatch;        // Previous average interval length.
    float*  pfMinRot;           // Minimum rotation in the scan
    int*    pnRebatchCount;     // Number of new batches in the scan. (important to avoid 'widow' batches).
    Cstring*    psScanID;       // ScanID for each scan.
    
    
    nStat = 2;  // Return an error if we exit before we are done.
    
    // Per Reflection.
    pbInScan  = new bool[m_nNumRefs];
    
    // Per Batch.
    pnOverlaps = new int[m_nNumBatches];
    pfBatchStart = new float[m_nNumBatches];
    pfBatchEnd = new float[m_nNumBatches];
    pnBatchI    = new int[m_nNumBatches];
    
    // Per Scan.
    pfSuggestedRebatch = new float[m_nNumBatches];    
    pfPrevBatch        = new float[m_nNumBatches];
    pfMinRot           = new float[m_nNumBatches];
    psScanID           = new Cstring[m_nNumBatches];
    pnRebatchCount     = new int[m_nNumBatches];
    if (m_oReflnlist.m_nFI_fObsRotMid<0) {
        printf("Could not find rotation midpoint datum in reflection file.\n");
        goto exit_place;
    };
    
    // First scan starts at first sorted reflection.
    sTemp = m_oReflnlist[m_pnIndexBatch[0]].sGetField(m_oReflnlist.m_nFI_sBatch);
    nLastBatch = nGetBatchNumber(sTemp, sScanID);
    if (nLastBatch<0) {
        goto exit_place;
    };
    nFirstSortRef = 0;
    nBatches = 1;
    nScans = 0;
    
    
    if (bRebatchImageBoundry) {
        if (m_poScan->bIsAvailable()) {
            fImageWidth = m_poScan->m_poRotation->fGetIncrement();
        } else {
            printf("WARNING:  Unable to rebatch on image boundries.  No scan info. available.\n");
            fImageWidth = 0.0;
        };
    };
    
    
    printf("Scan/Batch Information\n");
    printf("-----------------------------------------------------------------------\n");
    printf("   Scan  First    Last   Num   Rot   Rot  Rot/  Avg   Min   Max >=%4d  \n",nMinBatchOverlaps);
    printf("         batch   batch batch start   end batch ovlp  ovlp  ovlp  ovlps \n");
    printf("-----------------------------------------------------------------------\n");
    
    nStat = 0;
    
    for (nRef = 0;nRef<= m_nNumRefs; nRef++) {
        if (nRef<m_nNumRefs) {
            if (nRef==nFirstSortRef) {
                pfBatchEnd[0] = pfBatchStart[0] = m_oReflnlist[m_pnIndexBatch[nFirstSortRef]].fGetField(m_oReflnlist.m_nFI_fObsRotMid);
                pnBatchI[0] = m_oReflnlist[m_pnIndexBatch[nFirstSortRef]].nGetField(m_nFI_nBatchIdx);
            };
            poRefln = m_oReflnlist.poGetRefln(m_pnIndexBatch[nRef]);
            sTemp2 = poRefln->sGetField(m_oReflnlist.m_nFI_sBatch);
            nx = nGetBatchNumber(sTemp2,sTemp);
            if (nx<0) {
                goto exit_place;
            };
            if (sTemp == sScanID) {
                if (nx != nLastBatch) {
                    pfBatchEnd[nBatches] = pfBatchStart[nBatches] = m_oReflnlist[m_pnIndexBatch[nRef]].fGetField(m_oReflnlist.m_nFI_fObsRotMid);
                    pnBatchI[nBatches] = m_oReflnlist[m_pnIndexBatch[nRef]].nGetField(m_nFI_nBatchIdx);

                    
                    nBatches++;
                    nLastBatch = nx;
                    continue;
                } else {
                    pfBatchEnd[nBatches-1] = max(pfBatchEnd[nBatches-1],m_oReflnlist[m_pnIndexBatch[nRef]].fGetField(m_oReflnlist.m_nFI_fObsRotMid) +
                        m_oReflnlist[m_pnIndexBatch[nRef]].fGetField(m_oReflnlist.m_nFI_fObsRotWidth)/2
                        );
                    pfBatchStart[nBatches-1] =  min(pfBatchStart[nBatches-1],m_oReflnlist[m_pnIndexBatch[nRef]].fGetField(m_oReflnlist.m_nFI_fObsRotMid) - 
                        m_oReflnlist[m_pnIndexBatch[nRef]].fGetField(m_oReflnlist.m_nFI_fObsRotWidth)/2
                        );
                    continue;
                };
            };
            nLastSortRef = nRef - 1;
        } else
            nLastSortRef = nRef - 1;
        // Now we are ready to process all reflections in this scan.
        
        for (nx=0;nx<nBatches;nx++) {
            if ((nx==0) || (pfBatchStart[nx]<pfMinRot[nScans]))
                pfMinRot[nScans] = pfBatchStart[nx];
        };
        
        // Mark all reflections that are in this scan with a true.
        for (nx=0;nx<m_nNumRefs;nx++)
            pbInScan[nx] = FALSE;
        for (nx=nFirstSortRef;nx<=nLastSortRef;nx++)
            pbInScan[m_pnIndexBatch[nx]] = TRUE;
        
        // Compute the overlaps for each batch.
        // We save the overlaps for the batches in the "pnOverlaps" array
        
        // These variables are printed.
        int     nOverlaps;
        int     nAngle1,nAngle2;
        int     nNumAboveSpec;
        double  fAverageOverlaps;
        int     nMinOverlaps;
        int     nMaxOverlaps;
        double  fMinAngle;
        double  fMaxAngle;
        double  fRangeAngle;
        double  fAverageAngleStep1; // Calculated from the width of each batch.
        double  fAverageAngleStep2; // Caluclated from the overall scan width.
        double  fPercentPassed;
        
        fAverageOverlaps = 0;
        nMinOverlaps = 0;
        nMaxOverlaps = 0;
        fMinAngle = 0;
        fMaxAngle = 0;
        fAverageAngleStep1 = 0.0;
        fAverageAngleStep2 = 0.0;
        nNumAboveSpec = 0;
        
        for (nAngle1 = 0; nAngle1<nBatches;nAngle1++) {
            
            nOverlaps = 0;
            nOverlaps = m_pnNumOverlap[m_nNumBatches*pnBatchI[nAngle1] + pnBatchI[nAngle1]];

            if ((!nAngle1) || (nOverlaps>nMaxOverlaps)) 
                nMaxOverlaps = nOverlaps;
            if ((!nAngle1) || (nOverlaps<nMinOverlaps)) 
                nMinOverlaps = nOverlaps;
            fAverageOverlaps += nOverlaps;
            fAverageAngleStep1 += pfBatchEnd[nAngle1] - pfBatchStart[nAngle1];
            if ((!nAngle1) || (pfBatchEnd[nAngle1]>fMaxAngle)) 
                fMaxAngle = pfBatchEnd[nAngle1];
            if ((!nAngle1) || (pfBatchStart[nAngle1]<fMinAngle)) 
                fMinAngle = pfBatchStart[nAngle1];
            if (nOverlaps>=nMinBatchOverlaps)
                nNumAboveSpec++;
            pnOverlaps[nAngle1] = nOverlaps;
        };      
        
        fAverageOverlaps/=max(1,nBatches);
        fRangeAngle = fMaxAngle - fMinAngle;
        fPercentPassed = ((float) nNumAboveSpec)/max(1,nBatches);
        fAverageAngleStep1/=max(1,nBatches);
        fAverageAngleStep2 = fRangeAngle/max(1,nBatches);
        
        
        // Print the information for this scan.
        
        
        printf(" %3s???%7s %7s %5d",
            sScanID.string(),
            m_oReflnlist[m_pnIndexBatch[nFirstSortRef]].sGetField(m_oReflnlist.m_nFI_sBatch).string(),
            m_oReflnlist[m_pnIndexBatch[nLastSortRef]].sGetField(m_oReflnlist.m_nFI_sBatch).string(),
            nBatches
            );
        
        printf(" %5.1f %5.1f %5.1f",fMinAngle,fMaxAngle,fAverageAngleStep2);
        printf(" %4d %5d %5d",(int) fAverageOverlaps,nMinOverlaps,nMaxOverlaps);            
        printf(" %6.1f\n",fPercentPassed*100.0);
        
        if ((bRebatch) && (fPercentPassed<0.95) && (!bRebatchBookends)) {
            // Compute what the batching ratio should be.
            Cstat oStat;
            double fBatchWidth;
            
            // Find the average ratio of batch width to overlaps.
            oStat.vClear();
            for (nAngle1 = 0; nAngle1<nBatches;nAngle1++) {
                f0 = (pfBatchEnd[nAngle1] - pfBatchStart[nAngle1]);
                if (f0 >0.001) {
                    oStat.vAdd(max(1,pnOverlaps[nAngle1])/max(1e-10,f0));
                };
            };
            fBatchWidth = nMinBatchOverlaps/max(1e-10,(oStat.fAverage() + 3.0*oStat.fStandardDeviation()));
            // If the suggested rebatch is *smaller* keep the old batch width.
            if (fBatchWidth<fAverageAngleStep2) 
                fBatchWidth = fAverageAngleStep2;
            
            // If we require that batch boundries occur on image rotation width boundries, adjust the width.
            if ((bRebatchImageBoundry) && (fImageWidth!=0.0)) {
                fBatchWidth = fImageWidth*ceil(fBatchWidth/max(1e-10,fImageWidth));
                pnRebatchCount[nScans] = (int) ceil(fRangeAngle/max(1e-10,fBatchWidth));
                // Possibly, the last batch needs to get 'rounded' into the previous batch.
                // This will not occur unless we are doing image boundry batching.
                if ((pnRebatchCount[nScans]>1) && (pnRebatchCount[nScans]*fBatchWidth - fRangeAngle<fImageWidth*0.3))
                    pnRebatchCount[nScans]--;
                
            } else {
                // Adjust for an equally distributed rebatch (no widow batches).
                pnRebatchCount[nScans] = (int)ceil(fRangeAngle/max(1e-10,fBatchWidth));
                fBatchWidth = fRangeAngle/max(1,(pnRebatchCount[nScans]));
            };
            
            
            pfSuggestedRebatch[nScans] = fBatchWidth;
            pfPrevBatch[nScans] = fAverageAngleStep2;
        } else {
            // Tells the rebatch below that we are not interested in rebatching.
            pfSuggestedRebatch[nScans] = 0.0;
        };

        if ((bRebatchBookends) && (bRebatch)) {
            // Rebatch the bookends NOW.

            Cstat oStat;
            double fAvgOverlaps;
            Cstring sBatchToFind;
            Cstring sBatchToReplace;
            int nRef2;


            // Find the average ratio of batch width to overlaps.
            oStat.vClear();
            for (nAngle1 = 0; nAngle1<nBatches;nAngle1++) {
                oStat.vAdd(pnOverlaps[nAngle1]);
            };
            fAvgOverlaps = oStat.fAverage() - 3.0*oStat.fStandardDeviation();
            for (nAngle1 = 0; nAngle1<nBatches;nAngle1++) {
                if ((pnOverlaps[nAngle1] < fAvgOverlaps) && (pnOverlaps[nAngle1] < oStat.fAverage()*0.5)) {
                    // Remove this batch.
                    sBatchToFind = m_ptBatchInfo[pnBatchI[nAngle1]].m_sBatchName;
                    sBatchToReplace = "";
                    f1 = 1000.0;
                    // Find the nearest batch that has enough overlaps already.
                    for (nAngle2= 0; nAngle2 < nBatches; nAngle2++) {
                        if ((nAngle2 != nAngle1) && 
                            ((f0 = ABS((pfBatchStart[nAngle1] + pfBatchEnd[nAngle1]) - 
                            (pfBatchStart[nAngle2] + pfBatchEnd[nAngle2]))) < f1)) {
                            f1 = f0;

                            sBatchToReplace = m_ptBatchInfo[pnBatchI[nAngle2]].m_sBatchName;
                        };
                    };
                    if (f1 != 1000.0) {
                        for (nRef2 = 0;nRef2< m_nNumRefs; nRef2++) {
                            poRefln = m_oReflnlist.poGetRefln(m_pnIndexBatch[nRef2]);
                            if (sBatchToFind == poRefln->sGetField(m_oReflnlist.m_nFI_sBatch)) {
                                poRefln->vSetField(m_oReflnlist.m_nFI_sBatch,sBatchToReplace);
                            };
                        };
                    };
                    nStat = 1;
                };
            };
        };

        psScanID[nScans] = sScanID;
        nScans++;
        
        // End of scan processing.  Set the new parameters for the next scan.
        if (nRef != m_nNumRefs) {
            poRefln = m_oReflnlist.poGetRefln(m_pnIndexBatch[nRef]);
            sTemp = poRefln->sGetField(m_oReflnlist.m_nFI_sBatch);
            nLastBatch = nGetBatchNumber(sTemp,sScanID);
            nFirstSortRef = nRef+1;
            nBatches = 1;
        };
    };
    printf("-----------------------------------------------------------------------\n");
    
    
    
    if ((bRebatch) && (!bRebatchBookends)) {
        int nScan;
        int nNewBatch;
        char pcNewBatch[20];
        Cstring sNewBatch;
        Cstring sTemp;
        
        // See if we need to rebatch.
        for (nScan=0;nScan<nScans;nScan++) {
            
            if (pfSuggestedRebatch[nScan]!=0.0) {
                nStat = 1;
                printf("Rebatching scan %3s???.\n",psScanID[nScan].string());
                printf("Old width = %7.2f New width = %7.2f\n",pfPrevBatch[nScan],pfSuggestedRebatch[nScan]);
                
                nNewBatch = 0;
                for (nRef = 0;nRef< m_nNumRefs; nRef++) {
                    poRefln = m_oReflnlist.poGetRefln(m_pnIndexBatch[nRef]);
                    sTemp = poRefln->sGetField(m_oReflnlist.m_nFI_sBatch);
                    nx = nGetBatchNumber(sTemp,sScanID);
                    if (nx<0) {
                        goto exit_place;
                    };
                    if (sScanID == psScanID[nScan]) {
                        nNewBatch = min(pnRebatchCount[nScan]-1, 
                            (int)
                            ((poRefln->fGetField(m_oReflnlist.m_nFI_fObsRotMid) - pfMinRot[nScan])/max(1e-10,pfSuggestedRebatch[nScan])));
                        sprintf(pcNewBatch,"%s%4.4d",psScanID[nScan].string(),nNewBatch);
                        sNewBatch = pcNewBatch;
                        poRefln->vSetField(m_oReflnlist.m_nFI_sBatch,sNewBatch);
                    };
                };
            } 
        };
    };


    delete[] pbInScan;
    delete[] pnOverlaps;
    delete[] pfBatchStart;
    delete[] pfBatchEnd;
    delete[] pnBatchI;
    delete[] pfSuggestedRebatch;
    delete[] pfPrevBatch;
    delete[] pfMinRot;
    delete[] psScanID;
    delete[] pnRebatchCount;
exit_place:
    return nStat;
};



int Cscaleaverage::nSphericalHarmonics(tScaleProcessOptions& oOption)
{
    int nx;
    double f0,f1;
    double a3fVec1[3];
    double a3fVec2[3];
    double a3x3fMat[3][3];

    int nNumWithProblems;
    int nNumWithSeriousProblems;
    Cstring sNextBatchOrient;
    Cstring sLastBatchOrient;
    Cstring sBatchOrient;
    int nNumOrients;
    bool bOrientFound;
    int nRef,nRefSort;
    int nBatch;
    int nRefsInOrient;
    int nFirstBatch,nLastBatch;

    double a3fStartRot[3];      // S0 at start of rotation scan.
    double a3fEndRot[3];        // S0 at end of rotation scan.
    double a3fOppositeRot[3];   // S0 far away from a3fStartRot.
    double a3fAvgRot[3];        // Average S0 vector.
    double a3fRot[3];           // Temporary S0 vector
    double a3x3fRotMat[3][3];   // Temporary.
    double fRot;                // Temporary.
    double a3x3fOffsetRot[3][3];// Rotation matrix that takes S0 at fStartRot to <0,0,1>
    double fStartRot;           // Rot_Mid Start.
    double fEndRot;             // Rot_Mid End.
    double fOppositeRot;        // Rot_Mid Opposite.
    double fLengthRot;          // Range of Rot_Mid between 0 and 180.
    double fLengthRot360;       // Range of Rot_Mid between 0 and 360.
    double fMaxLength;          // Maximum length of a chord.
    double fCircleRadius;       // The tips of S0 vectors lie in a plane.  This is the radius of the intersection of plane with unit sphere.
    double a3fPerpS0[3];        // Perpendicular vector for plane of S0 vectors.
    double fPerpVecSign = 1.0;  // Might need to multiply the perpendicular vector by a sign.
    double fPolar;              // Polar coordinate.
    double fAzi;                // Azimuthal coordinate.
    double fMinPolar = 360.0;   // Minimum polar coordinate found.
    double fMaxPolar = -360.0;  // Maximum polar coordinate found.
    
    std::vector<Cstring>            asBatchPrefixes;

    // Arrays for actual spherical harmonic code.
    itr<int> anGroup;
    itr<double> afScanRot;
    itr<double> afValues;
    itr<double> afValuesSigma;
    itr<double> afPolar;
    itr<double> afAzi;
    itr<double> afScaleFactorsOut;

    
    nNumWithProblems = 0;
    nNumWithSeriousProblems = 0;
    sNextBatchOrient = "";
    sLastBatchOrient = "";
    asBatchPrefixes.clear();
    afAzi.setsize(m_oReflnlist.nGetNumReflns());
    afPolar.setsize(m_oReflnlist.nGetNumReflns());
    memset(&afAzi[0],0,sizeof(double)*afAzi.size());
    memset(&afPolar[0],0,sizeof(double)*afPolar.size());
           
    printf("\n\n");
    printf("Spherical Harmonic correction scans.\n");
    printf("-------------------------------------------------------------------\n");
    printf("   Scan First  Last    Rot    Rot    Rot     S0    S0    S0    S0\n");
    printf("        batch batch  start    end  range radius   [x]   [y]   [z]\n");
    printf("-------------------------------------------------------------------\n");
    nNumOrients = 0;
    do {
        sLastBatchOrient = sNextBatchOrient;
        bOrientFound = FALSE;
        for (nRef=0;nRef<m_oReflnlist.nGetNumReflns();nRef++) {
            nBatch = m_ptBatchInfo[m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx)].m_nBatch;
            sBatchOrient = m_ptBatchInfo[m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx)].m_sScan;
            if (((!bOrientFound) || (sBatchOrient<sNextBatchOrient)) &&
                ((nNumOrients==0) || (sBatchOrient>sLastBatchOrient))) {
                sNextBatchOrient = sBatchOrient;
                asBatchPrefixes.push_back(sBatchOrient);
                bOrientFound = TRUE;
            };
        };
        
        if ((nNumOrients==0) || (sNextBatchOrient>sLastBatchOrient)) {
            // A new orientation was found.  
            
            // Find a maximum chord on the circle described by
            // the tips of the S0 vectors.  (The tips will lie in a plane,
            // but not necc. a plane passing through the center of the unit sphere).
            nRefsInOrient = 0;
            nFirstBatch = 1000;
            nLastBatch = -1;
            vZeroMat(3,1,a3fAvgRot);
            for (nRef=0;nRef<m_oReflnlist.nGetNumReflns();nRef++) {
                nBatch = m_ptBatchInfo[m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx)].m_nBatch;
                sBatchOrient = m_ptBatchInfo[m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx)].m_sScan;
                if (sBatchOrient==sNextBatchOrient) {
                    fRot = m_oReflnlist[nRef].fGetField(m_oReflnlist.m_nFI_fObsRotMid);
                    for (nx=0;nx<3;nx++)
                        a3fRot[nx] = m_oReflnlist[nRef].fGetField(m_oReflnlist.m_nFI_fS0vec[nx]);
                    
                    if ((nLastBatch==-1) || (fRot<fStartRot)) {
                        vCopyVec3D(a3fRot,a3fStartRot);
                        fStartRot = fRot;                          
                    };
                    if ((nLastBatch==-1) || (fRot>fEndRot)) {
                        vCopyVec3D(a3fRot,a3fEndRot);
                        fEndRot = fRot;                          
                    };
                    if ((nLastBatch==-1) || ((fRot>fOppositeRot) && (fRot-fStartRot<180.0))) {
                        vCopyVec3D(a3fRot,a3fOppositeRot);
                        fOppositeRot = fRot;
                    };
                    nFirstBatch = min(nFirstBatch,nBatch);
                    nLastBatch = max(nLastBatch,nBatch);
                    
                    
                    vAddVec3DVec3D(a3fRot,a3fAvgRot,a3fAvgRot);
                    nRefsInOrient++;
                };
            };
            // Compute average S0 vector
            vMulVec3DScalar(a3fAvgRot,1.0/max(1,nRefsInOrient),a3fAvgRot);
            
            vSubVec3DVec3D(a3fOppositeRot,a3fStartRot,a3fRot);
            fMaxLength = fLenVec3D(a3fRot);
            fLengthRot = fOppositeRot - fStartRot;
            fLengthRot360 = min(360.0,fEndRot-fStartRot);
            
            fCircleRadius = sqrt(fMaxLength*fMaxLength/2.0/max(1e-10,(1.0 - cos(Gs_dRADIANS_PER_DEGREE*fLengthRot))));
            
           
            // Compute perpendicular S0 vector.
            vSubVec3DVec3D(a3fOppositeRot,a3fAvgRot,a3fVec1);
            vSubVec3DVec3D(a3fStartRot,a3fAvgRot,a3fVec2);                
            if ((fLenVec3D(a3fVec1)>0.001) && (fLenVec3D(a3fVec2)>0.001)) 
                vCross3D(a3fVec1,a3fVec2,a3fPerpS0);
            else
                vCopyVec3D(a3fAvgRot,a3fPerpS0);
            fNormVec3D(a3fPerpS0);

            // Compute the rotation matrix that would take an S0 vector at fStartRot to <0,0,1>
            vBuildBasis3D(a3fStartRot,a3x3fRotMat);
            vNormMat3D(a3x3fRotMat);
            vTranMat3D(a3x3fRotMat);
            vZeroMat3D(&a3x3fMat[0][0]);
            a3x3fMat[1][0] = 1.0;
            a3x3fMat[2][1] = 1.0;
            a3x3fMat[0][2] = 1.0;
            vMulMat3DMat3D(a3x3fMat,a3x3fRotMat,a3x3fOffsetRot);
                       
            // Print out information on this scan.
            printf(" %3s???   %3.3d   %3.3d %6.2f %6.2f %6.2f %6.2f %5.2f %5.2f %5.2f\n",
                sNextBatchOrient.string(),
                nFirstBatch,
                nLastBatch,
                fStartRot,fEndRot,fLengthRot360,
                fCircleRadius,a3fPerpS0[0],a3fPerpS0[1],a3fPerpS0[2]
                );

            // Determine Azimuthal and Polar coordinates for each reflection.
            for (nRefSort=0;nRefSort<m_oReflnlist.nGetNumReflns();nRefSort++) {
                nRef = m_pnIndex[nRefSort];
                nBatch = m_ptBatchInfo[m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx)].m_nBatch;
                sBatchOrient = m_ptBatchInfo[m_oReflnlist[nRef].nGetField(m_nFI_nBatchIdx)].m_sScan;
                if (sBatchOrient==sNextBatchOrient) {
                    
                    // Find the rotation matrix that
                    // transforms the S0 vector so that it points down <0,0,1>
                    // This is composed of a matrix taking a3fRot to a3fStartRot, followed
                    // by a3x3fOffsetRot which takes it to <0,0,1>
                    f0 = 5.0;
                    for (nx=0;nx<3;nx++)
                        a3fRot[nx] = m_oReflnlist[nRef].fGetField(m_oReflnlist.m_nFI_fS0vec[nx]);

                    for (nx=0;nx<3;nx++) {
                        vConvRotVec3DMat3D(fPerpVecSign*(fStartRot-m_oReflnlist[nRef].fGetField(m_oReflnlist.m_nFI_fObsRotMid)),&a3fPerpS0[0],&a3x3fMat[0][0]);
                        vMulMat3DMat3D(a3x3fOffsetRot,a3x3fMat,a3x3fRotMat);
                        
                        // Make sure that the S0 vector is getting transformed correctly.
                        vMulMat3DVec3D(a3x3fRotMat,a3fRot,a3fVec1);
                        f1 = fabs(a3fVec1[0]) + fabs(a3fVec1[1]) + fabs(a3fVec1[2]-1.0);
                        
                        // This complicated logic results in the smallest residual getting chossen.
                        if (f1<0.05)
                            break;                            
                        else if (f1<f0) {
                            f0 = f1;
                            fPerpVecSign *= -1.0;
                            if (nx==1)
                                break;
                        } else if (f1==f0) 
                            break;
                        else
                            fPerpVecSign *= -1.0;
                    };                
                    if (nx>=2) {
                        nNumWithProblems++;
                    };                     
                    for (nx=0;nx<3;nx++)
                        a3fRot[nx] = -m_oReflnlist[nRef].fGetField(m_oReflnlist.m_nFI_fSvec[nx]);
                    
                    // Transform the S vector accordingly.  It should now lie in the correct hemisphere. (i.e. [2]>=0.0)
                    
                    vMulMat3DVec3D(a3x3fRotMat,a3fRot,a3fVec1);                
                    // In case we want to bypass the rotation correction: vCopyVec3D(a3fRot,a3fVec1);
                    CSphericalHarmonic::vVec3DToSpherical(&a3fVec1[0],fPolar,fAzi);                
                    fMinPolar = min(fPolar,fMinPolar);
                    fMaxPolar = max(fPolar,fMaxPolar);
                    afPolar[nRefSort] = fPolar;
                    afAzi[nRefSort] = fAzi;    
                };
            };
            
        } else
            break;

        nNumOrients++;
    } while (1);
    printf("---------------------------------------------------------------------\n");
    if (nNumWithSeriousProblems) {
        printf("ERROR: %d bad reflections detected!  Algorithm terminating.\n",nNumWithSeriousProblems);
        return 1;
    };
    if (nNumWithProblems) {
        printf("WARNING: %d reflections had suspicious S0/S vectors!\n",nNumWithProblems);
    };


    // Now build the arrays to send to the absorbtion correction.

    int nLastHKL,nThisHKL;
    int nStartIndex,nEndIndex,nIndex;
    int nGroup;
    int nCode;
    int nScan;
    double fInitRMerge,fFinalRMerge,fFinalRMergeNoErrorModel;

    CSphericalHarmonic oAbsorb;

    oAbsorb.vSetHarmonics(oOption.nHarmonicOrderDifrac,oOption.nHarmonicOrderIncident,2);
    oAbsorb.vSetScaleInfo(m_fSigma,5.0,oOption.nMaxIter,0.0);
    nLastHKL = m_oReflnlist[m_pnIndex[0]].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    nStartIndex = 0;
    nGroup= 1;
    
    for (nRefSort=0;nRefSort<=m_nNumRefs;nRefSort++)  {
        nRef=m_pnIndex[nRefSort];
        
        if (nRefSort<m_nNumRefs)
            nThisHKL = m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
        if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
            continue;
        nEndIndex=nRefSort-1;

        for (nIndex = nStartIndex; nIndex <= nEndIndex; nIndex++) {
            nBatch = m_ptBatchInfo[m_oReflnlist[m_pnIndex[nIndex]].nGetField(m_nFI_nBatchIdx)].m_nBatch;
            sBatchOrient = m_ptBatchInfo[m_oReflnlist[m_pnIndex[nIndex]].nGetField(m_nFI_nBatchIdx)].m_sScan;
            for (nScan = 0; nScan < asBatchPrefixes.size(); nScan++) {
                if (asBatchPrefixes[nScan] == sBatchOrient)
                    break;
            };
            
            CSphericalHarmonic::vSetGroupInfo(nCode,nGroup,nScan,(m_pnReject[nIndex]==1)?g_nSphericalFlagUserReject:0);

            anGroup + nCode;
            afScanRot + m_oReflnlist[m_pnIndex[nIndex]].fGetField(m_oReflnlist.m_nFI_fObsRotMid);
            afValues + m_oReflnlist[m_pnIndex[nIndex]].fGetIntensity();
            afValuesSigma + fCalcSigma(m_oReflnlist[m_pnIndex[nIndex]],m_pnIndex[nIndex]);

            // We already have afPolar and afAzi filled in.
        };
              
        nGroup ++;
        nStartIndex=nRefSort;
        if (nRefSort<m_nNumRefs)
            nLastHKL =  m_oReflnlist[nRef].nGetField(m_oReflnlist.m_nFI_nPackedHKL);
    };

    fInitRMerge = fRmerge();    

    
    if (oAbsorb.nScale(afValues,afValuesSigma,afPolar,afAzi,anGroup,afScanRot,afScaleFactorsOut))
        return 1;
	nNumWithSeriousProblems = 0;
    for (nRefSort=0;nRefSort<m_nNumRefs;nRefSort++)  {
        nRef=m_pnIndex[nRefSort];
		if (afScaleFactorsOut[nRefSort]<=0.0) {
			nNumWithSeriousProblems++;
		};
        m_oReflnlist[nRef].vSetIntensity(m_oReflnlist[nRef].fGetIntensity()/max(1e-10,afScaleFactorsOut[nRefSort]));
        m_oReflnlist[nRef].vSetSigmaI(m_oReflnlist[nRef].fGetSigmaI()/max(1e-10,afScaleFactorsOut[nRefSort]));
    };
	if (nNumWithSeriousProblems) {
		printf("WARNING:  %d reflections had suspicious scale factors.  Factors will not be applied.\n",nNumWithSeriousProblems);
		return 1;
	};

	// Print table of transmission factors in each resolution shell.
	{
		
		double fReso;
		int nBinReso;
		double fAvgTrans;
		double fMinTrans;
		double fMaxTrans;
		double fTrans;
		float  fResoLow,fResoHigh,fLow,fLowStart,*pfSlope;
		int    nContrib;
		int    nUsed,nOutliers;
		int    nRefGroup,nRefScan,nRefReject,nRefBinReso;
		int    nTotalScans;
		int    nScanStart,nScanEnd;
		
		nTotalScans = asBatchPrefixes.size();

		printf("\nSpherical correction factors vs Resolution\n");
		printf("========================================================================\n");
		printf("    Resolution   Scan    Avg    Max    Min   Max/  Total  Total Outliers\n");
		printf("         range         trans  trans  trans    Min  avail   used  removed\n");
		printf("========================================================================\n");
		fLowStart = fLow     =  m_a2fRangeReso[0];
		pfSlope  = &m_fSlopeReso;

		for (nBinReso = 0; nBinReso < m_nNumBinsReso + 1; nBinReso++) {
			
			fResoLow              = 1.0f / max(1e-10,pow((double)fLow, 0.333333));
			fResoHigh             = 1.0f / max(1e-10,pow((double)(fLow + *pfSlope), 0.333333));

			// For the last reso bin, we are using all reflections.
			if (nBinReso == m_nNumBinsReso) {
				nScanStart = nTotalScans;
				nScanEnd = nTotalScans;
				fResoHigh = fResoLow;
				fResoLow = 1.0f / max(1e-10,pow((double)fLowStart, 0.333333));
			} else  if (nTotalScans == 1) {
				nScanStart = 0;
				nScanEnd = 0;
			} else {
				nScanStart = 0;
				nScanEnd = nTotalScans;
			};

			for (nScan = nScanStart; nScan <= nScanEnd; nScan++) {
				nContrib = 0;
				nOutliers = 0;
				nUsed = 0;
				fMinTrans = 0.0;
				fMaxTrans = 0.0;
				fAvgTrans = 0.0;
				for (nRefSort = 0; nRefSort < m_nNumRefs;nRefSort++) {
					CSphericalHarmonic::vGetGroupInfo(anGroup[nRefSort],nRefGroup,nRefScan,nRefReject);
					if ((nRefScan == nScan) || (nScan == nTotalScans)) {
						nRef = m_pnIndex[nRefSort];
						fReso    = m_oReflnlist[nRef].fGetField(m_nFI_f2STLsq);
                                                fReso    = fReso * sqrt(fReso);
						nRefBinReso = (int) ((fReso - m_a2fRangeReso[0]) / max(1e-10,m_fSlopeReso));
						if ((nRefBinReso == nBinReso) || (nBinReso == m_nNumBinsReso)) {
							fTrans = 1.0/afScaleFactorsOut[nRefSort];
							fMaxTrans = (nContrib==0)?fTrans:max(fTrans,fMaxTrans);
							fMinTrans = (nContrib==0)?fTrans:min(fTrans,fMinTrans);
							fAvgTrans += fTrans;
							if (nRefReject & g_nSphericalFlagOutlierReject)
								nOutliers++;
							if (!(nRefReject & (g_nSphericalFlagOutlierReject | g_nSphericalFlagTooLowRedundancy | g_nSphericalFlagExclude)))
								nUsed++;
							nContrib++;
						};
					};
				};
				if (nBinReso == m_nNumBinsReso) {
					printf("========================================================================\n");
				};
				if ((nScan == 0) || (nBinReso == m_nNumBinsReso))
					printf(" %6.2f -%5.2f",min((float) 999.0,fResoLow),fResoHigh);
				else
					printf(" %13s","");
				if (nScan == nTotalScans)
					printf(" %6s","All");
				else
					printf(" %3s???",asBatchPrefixes[nScan].string());

				printf(" %6.2lf %6.2lf %6.2lf %6.2lf %6d %6d %8d\n",fAvgTrans/max(1,nContrib),fMaxTrans,fMinTrans,fMaxTrans/max(0.00001,fMinTrans),nContrib,nUsed,nOutliers);
				
			};
			printf("\n");
			fLow = fLow + *pfSlope;
		};
		

	};
	

    fFinalRMergeNoErrorModel = fRmerge();
    nFitChiSquared(false);
    fCalcChiSquared(1);

    fFinalRMerge = fRmerge();

    printf("Absorption correction results\n");
    printf("--------------------------------------------------------------------------\n");
    printf("\nRmerge before applying scale factors:                         %.4lf\n", fInitRMerge);
    printf(  "Rmerge after  applying scale factors with    adjusted sigmas: %.4lf\n", fFinalRMerge);
    printf(  "Rmerge after  applying scale factors without adjusted sigmas: %.4lf\n\n", fFinalRMergeNoErrorModel);
    printf("(All R-merge values are computed for reflections having I/sig(I) >= %.1lf)\n",m_fSigma);
    printf("--------------------------------------------------------------------------\n\n");


    
    return 0;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
void Cscaleaverage::vSetScaleRestrainInfo(Cstring sInfo)
{
    const double        cdDefaultWeight = 0.1;
    Cstring             sMinusOne("-1");

    m_vecScaleRestrInfo.clear();

    if( "" == sInfo )
    {
        m_vecScaleRestrInfo.push_back(SCALE_RESTRAIN_INFO(sMinusOne, 
                                                          sMinusOne,
                                                          cdDefaultWeight)); // restrain all batches in every scan;
        return;
    }
    else if( !sInfo.contains('-') ) // see if some explicit batch ranges are given
    {
        double      dWeight = atof(sInfo.string());
        m_vecScaleRestrInfo.push_back(SCALE_RESTRAIN_INFO(sMinusOne,
                                                          sMinusOne,
                                                          dWeight)); // restrain all batches in every scan with weight dWeight
        return;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Cstring     sTemp(sInfo);
    std::vector<Cstring>        saRanges;
    sInfo.nListToVector(saRanges, ", ");


    std::vector<Cstring>    saRestrInfo;
    double                  dWeight = 0.0;
    for(int ii=0; ii < saRanges.size(); ii++)
    {
        saRanges[ii].nListToVector(saRestrInfo, "- ");

        if( saRestrInfo.size() != 2 && saRestrInfo.size() != 3 ) //expecting input like 10001-10360 or 10001-10360-0.05
            continue;
        
        dWeight = saRestrInfo.size() == 3 ? atof(saRestrInfo[2]) : cdDefaultWeight;
        
        m_vecScaleRestrInfo.push_back(SCALE_RESTRAIN_INFO(saRestrInfo[0], saRestrInfo[1], dWeight));
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Make sure ALL batches within EACH scan are restrained WITHIN THE SCAN with the given WEIGHT.
void Cscaleaverage::vGenerateDefaultBatchRestrainInfo(double dWeight)
{
    m_vecScaleRestrInfo.clear();
    
    if( m_nNumBatches < 2 )   // should have at least two batches to scale!
        return;

    int         nIndexCurrentScanBegin = 0;

    Cstring     sCurrentScan = m_ptBatchInfo[0].m_sScan;

    Cstring     sCurrentScanSave(sCurrentScan); 

    int         iBatch = 1; // loop variable

    while( iBatch < m_nNumBatches + 1 )
    {
        if( iBatch < m_nNumBatches )
        {
            sCurrentScan = m_ptBatchInfo[iBatch].m_sScan;
        }
        // else do not update the scan number, because there are no more batches

        if( sCurrentScanSave != sCurrentScan || iBatch == m_nNumBatches ) // next scan started or no more batches
        {
            m_vecScaleRestrInfo.push_back(SCALE_RESTRAIN_INFO(m_ptBatchInfo[nIndexCurrentScanBegin].m_sBatchName,
                                                              m_ptBatchInfo[iBatch-1].m_sBatchName,
                                                              dWeight));
            sCurrentScanSave = sCurrentScan;

            nIndexCurrentScanBegin = iBatch;
        }
        
        iBatch++;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The purpose of this function is to convert batch names in the batch restrain info to the indices in the list of the available batches.
// For example, the batch restrain info might say that batches 10001-10360 and 20010-20360 need to be restrained. At the same time, in the batch 
// list, batches 10001-10360 are actually indexed as 0-359, while batches 20010-20360 are indexed as 369-719. 
// Another thing to keep in mind: some batches might have already been rejected, so restrain info must be corrected for that as well.
void Cscaleaverage::vConvertBatchRestrainInfoToBatchIndices()
{
    // Map batch names to batch list indices
    std::map<Cstring, int>                  mapBatchIndices;
    std::map<Cstring, int>::iterator        oIt;

    Cstring                                 sBatchName("");
    Cstring                                 sScanID("");
    
    for(int ii=0; ii < m_nNumBatches; ii++)
    {
        sBatchName = m_ptBatchInfo[ii].m_sBatchName;

        mapBatchIndices.insert(std::pair<Cstring,int>(sBatchName, ii));
    }

    // Now go through the array of restrain info and convert the Begins and Ends to corresponding indices

    // Since some of the batches could have been rejected, there could be cases, when a Begin or an End or both cannot be found.
    
    // If Begin is not found, set it at the next index, unless the next index corresponds to a batch number that is larger than End.
    // If End is not found, set it at the previous index, unless the previous index corresponds to a batch number that is smaller than Begin.

    // Invalidate restrain info that cannot be meaningfully transferred into indices.

    std::vector<int>        anRestrInfoInvalid;

    int     nBatchBeginIndex = 0;
    int     nBatchEndIndex   = 0;

    int     nBegin = 0;
    int     nEnd = 0;

    bool    bFound = false;

    for(int jj=0; jj < m_vecScaleRestrInfo.size(); jj++)
    {
        ///////////////////////////////////////////////////////////////////////////////////////
        // Figure out the begin index
        
        bFound = false;
        
        oIt = mapBatchIndices.find(m_vecScaleRestrInfo[jj].m_sBatchBegin);

        if( oIt != mapBatchIndices.end() )
        {
            nBatchBeginIndex = (*oIt).second;
            m_vecScaleRestrInfo[jj].m_sBatchBegin = nBatchBeginIndex;
            bFound = true;
        }
        else  // The batch is not found. Try to find the closest available batch.
        {
            nBegin = nGetBatchNumber(m_vecScaleRestrInfo[jj].m_sBatchBegin, sScanID); 
            nEnd   = nGetBatchNumber(m_vecScaleRestrInfo[jj].m_sBatchEnd, sScanID);

            for(int iii=nBegin + 1; iii < nEnd; iii++)   // + 1, because the begin batch itself is not found
            {
                vMakeBatchName(iii, sScanID.string(), sBatchName);

                oIt = mapBatchIndices.find(sBatchName);

                if( oIt != mapBatchIndices.end() )
                {
                    nBatchBeginIndex = (*oIt).second;
                    m_vecScaleRestrInfo[jj].m_sBatchBegin = nBatchBeginIndex;
                    bFound = true;
                    break;
                }
            }
        }

        if( !bFound )
            anRestrInfoInvalid.push_back(jj);
    
        ////////////////////////////////////////////////////////////////////////////////////////
        
        ///////////////////////////////////////////////////////////////////////////////////////
        // Figure out the end index
        bFound = false;
        
        oIt = mapBatchIndices.find(m_vecScaleRestrInfo[jj].m_sBatchEnd);

        if( oIt != mapBatchIndices.end() )
        {
            nBatchEndIndex = (*oIt).second;
            m_vecScaleRestrInfo[jj].m_sBatchEnd = nBatchEndIndex;
            bFound = true;
        }
        else  // The batch is not found. Try to find the closest available batch.
        {
            nBegin = nGetBatchNumber(m_vecScaleRestrInfo[jj].m_sBatchBegin, sScanID); 
            nEnd   = nGetBatchNumber(m_vecScaleRestrInfo[jj].m_sBatchEnd,   sScanID);

            for(int iii=nEnd - 1; iii > nBegin; iii--)  // -1, because the end batch itself is not found
            {
                vMakeBatchName(iii, sScanID.string(), sBatchName);

                oIt = mapBatchIndices.find(sBatchName);

                if( oIt != mapBatchIndices.end() )
                {
                    nBatchEndIndex = (*oIt).second;
                    m_vecScaleRestrInfo[jj].m_sBatchEnd = nBatchEndIndex;
                    bFound = true;
                    break;
                }
            }
        }

        if( !bFound )
            anRestrInfoInvalid.push_back(jj);
   
        ////////////////////////////////////////////////////////////////////////////////////////
        // Final check
        nBegin = nGetBatchNumber(m_vecScaleRestrInfo[jj].m_sBatchBegin, sScanID); 
        nEnd   = nGetBatchNumber(m_vecScaleRestrInfo[jj].m_sBatchEnd, sScanID);

        if( nBegin >= nEnd ||
            nBegin < 0  || nBegin > m_nNumBatches - 1 ||
            nEnd   < 0  || nEnd   > m_nNumBatches - 1)
        {
            anRestrInfoInvalid.push_back(jj);
        }
    }

    // Invalidate restrain info that could not be meaningfully transferred to the batch indices
    for(int kk=0; kk < anRestrInfoInvalid.size(); kk++)
    {
        m_vecScaleRestrInfo[anRestrInfoInvalid[kk]].m_sBatchBegin = -2;
        m_vecScaleRestrInfo[anRestrInfoInvalid[kk]].m_sBatchEnd   = -2;
    }
}

