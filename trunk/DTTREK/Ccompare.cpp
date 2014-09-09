#include "Ccompare.h"
#include "Cstat.h"

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
#endif

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
// Ccompare.cc     Initial author: T.J. Niemeyer  Spring 2000
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



Ccompare::Ccompare(Creflnlist& oRef1,Creflnlist& oRef2,Cimage_header& oHeader):
m_oList1(oRef1),m_oList2(oRef2),m_oCrystal(oHeader) {
    int nx,ny;
    Cstring sTemp;

    m_poList[0] = NULL;
    m_poList[1] = NULL;
    m_pnReject[0] = NULL;
    m_pnReject[1] = NULL;
    m_pnIndex[0] = NULL;
    m_pnIndex[1] = NULL;
    m_pnIndexU[0] = NULL;
    m_pnIndexU[1] = NULL;
    m_pnIndexRMerge[0] = NULL;
    m_pnIndexRMerge[1] = NULL;
    m_pfIntensityBinsNumer[0] = NULL;
    m_pfIntensityBinsNumer[1] = NULL;
    m_pfIntensityBinsDenom[0] = NULL;
    m_pfIntensityBinsDenom[1] = NULL;
    m_pnIntensityBinsContrib[0] = NULL;
    m_pnIntensityBinsContrib[1] = NULL;

    m_nStat = 0;

    for (nx=0;nx<2;nx++) {
        m_pfIntensityBinsNumer[nx] = NULL;
        m_pfIntensityBinsDenom[nx] = NULL;
        m_pnIntensityBinsContrib[nx] = NULL;
    };


    m_bResortRMergeArray[0]= TRUE;
    m_bResortRMergeArray[1]= TRUE;
    m_bSilent = FALSE;
    if (!m_oCrystal.bIsAvailable()) {
        printf("ERROR:  Could not load crystal object!\n");
        m_nStat=1;
    } else {
        m_oList1.nPutCrystal(m_oCrystal);
        m_oList2.nPutCrystal(m_oCrystal);
    };

    
    m_poList[1] = &oRef2;
    m_poList[0] = &oRef1;
    // Build the reflection sort arrays.
    
    for (nx=0;(nx<2) && (!m_nStat);nx++) {
        m_nNumRefs[nx] = m_poList[nx]->nGetNumReflns();
        m_pnReject[nx]= new int[m_nNumRefs[nx]];
        m_pnIndexU[nx] = new int[m_nNumRefs[nx]];
        m_pnIndexRMerge[nx] = new int[m_nNumRefs[nx]];
            
            
        // Build unique HKL sort arrays.  
        m_nFI_nIndexU[nx] = m_poList[nx]->nExpandGetField(sTemp = "nTemp");
        for (ny=0;ny<m_nNumRefs[nx];ny++)
            (*m_poList[nx])[ny].vSetField(m_nFI_nIndexU[nx],(*m_poList[nx])[ny].nPackHKL());
        m_poList[nx]->vSort(eReflnField_int_type,m_nFI_nIndexU[nx],m_pnIndexU[nx]);

        // Build R-Merge sort arrays.
        m_nFI_fIndexRMerge[nx] = m_poList[nx]->nExpandGetField(sTemp = "fRMerge");
        
        cout << "Sorting and reflnlist to asymmetric unit ..." << endl << flush;
        if (!m_nStat) {
            if ((nx==0) || (m_poList[0]!=m_poList[1])) {
                m_nStat += m_poList[nx]->nReduce(m_oCrystal, FALSE);
                cout << "...done.\n" << flush;
            };
            m_pnIndex[nx] = m_poList[nx]->pnGetSortIndex();    
            m_nFI_nIndex[nx] = m_poList[nx]->nExpandGetField(Creflnlist::ms_snPackedHKL);
        };
        
    };
};

Ccompare::~Ccompare() {
    int nx;
    for (nx=0;nx<2;nx++) {
        if (m_pnReject[nx])
            delete[] m_pnReject[nx];
        if (m_pnIndexU[nx])
            delete[] m_pnIndexU[nx];
        if (m_pnIndexRMerge[nx])
            delete[] m_pnIndexRMerge[nx];
        if (m_pfIntensityBinsNumer[nx])
            delete[] m_pfIntensityBinsNumer[nx];
        if (m_pfIntensityBinsDenom[nx])
            delete[] m_pfIntensityBinsDenom[nx];
        if (m_pnIntensityBinsContrib[nx])
            delete[] m_pnIntensityBinsContrib[nx];
    };
};

int Ccompare::nClearRejects() {
    int nPass;
    int nx;

    if (!m_bSilent)
        printf("Clearing rejected reflections.\n");
    // Mark all reflections as non-rejected.
    for (nPass=0;nPass<2;nPass++) {
        for (nx=0;nx<m_nNumRefs[nPass];nx++)
            m_pnReject[nPass][nx] = 0;
    };
    return 0;
};

int Ccompare::nEquivocate(eComparisonType eType,
                          bool bEqualRedundancy,
                          int nMinRedundancy,
                          int nUseSigmasIn,
                          int nRejectFlag) 
{
    int nLastHKL;
    int nThisHKL;
    int nStartSource,nEndSource,nSourceFound;
    int nStartComp,nEndComp,nCompFound;
    int nRefSort;
    int nNumRefs;
    int nRef;
    int* pnHKLComp;

    int nPass;
    int nLoop;
    int nx,ny;

    bool bDoCCP4Comparision = (
        (m_poList[0]->m_nFI_fIntensityPlus>=0) &&
        (m_poList[0]->m_nFI_fIntensityMinus>=0) && 
        (m_poList[0]->m_nFI_fSigmaIPlus>=0) && 
        (m_poList[0]->m_nFI_fSigmaIMinus>=0) && 
        (m_poList[1]->m_nFI_fIntensityPlus>=0) &&
        (m_poList[1]->m_nFI_fIntensityMinus>=0) && 
        (m_poList[1]->m_nFI_fSigmaIPlus>=0) && 
        (m_poList[1]->m_nFI_fSigmaIMinus>=0));
    


    if (eType != APPLES_ORANGES) {
        if (!m_bSilent)
            printf("Trimming lists (%s) assuming (%s)\n"
            "with all terms in R-merge containing at least %d reflections.\n",
            (eType==APPLES_APPLES)?"apples to apples":"apples to macintoshes",
            bEqualRedundancy?"Equal Redundancy":"Any Redundancy",
            nMinRedundancy               
            );
        // If we are doing APPLES_APPLES, we need to loop twice, since we don't have any way
        // of detecting the redundancy of the data specified with nMinRedundancy.
        for (nLoop=0;nLoop<2;nLoop++) {

            // No need of a second loop if we are comparing APPLES_MAC
            if ((nLoop==1) && (eType == APPLES_APPLES))
                eType = APPLES_MAC;
            else if (nLoop==1)
                break;

            for (nPass=0;nPass<2;nPass++) {
                Creflnlist& oSource = *m_poList[nPass];
                Creflnlist& oComp = *m_poList[!nPass];
                int* pnIndexSource = (eType==APPLES_APPLES)?(m_pnIndexU[nPass]):(m_pnIndex[nPass]);
                int* pnIndexComp = (eType==APPLES_APPLES)?(m_pnIndexU[!nPass]):(m_pnIndex[!nPass]);
                int  nFieldSource = (eType==APPLES_APPLES)?(m_nFI_nIndexU[nPass]):(m_nFI_nIndex[nPass]);
                int  nFieldComp = (eType==APPLES_APPLES)?(m_nFI_nIndexU[!nPass]):(m_nFI_nIndex[!nPass]);
                
                
                nLastHKL = oSource[pnIndexSource[0]].nGetField(nFieldSource);
                nStartSource = 0;   
                nNumRefs = oSource.nGetNumReflns();
                for (nRefSort=0;nRefSort<nNumRefs+1;nRefSort++) {
                    nRef=pnIndexSource[nRefSort];
                    
                    if (nRefSort<nNumRefs)
                        nThisHKL = oSource[nRef].nGetField(nFieldSource);
                    if ((nRefSort!=nNumRefs) && (nThisHKL==nLastHKL))
                        continue;
                    nEndSource=nRefSort-1;
                                        
                    // Discover the number of non-rejected source reflections in this group.
                    nSourceFound = 0;
                    for (nx=nStartSource;nx<=nEndSource;nx++) {
                        if (m_pnReject[nPass][pnIndexSource[nx]]==0)
                            nSourceFound++;
                    };
                    
                    nStartComp = oComp.nFindFirst(nFieldComp,oSource[pnIndexSource[nStartSource]].nGetField(nFieldSource),pnIndexComp);
                    if (nStartComp<0)
                        nCompFound = 0;
                    else {
                        nCompFound = 0;
                        for (nEndComp=nStartComp;nEndComp<oComp.nGetNumReflns();nEndComp++) {
                            if (oComp[pnIndexComp[nEndComp]].nGetField(nFieldComp)!=oComp[pnIndexComp[nStartComp]].nGetField(nFieldComp))
                                break;
                            if (m_pnReject[!nPass][pnIndexComp[nEndComp]]==0)
                                nCompFound++;                    
                        };
                        nEndComp--;
                        if (oComp[pnIndexComp[nStartComp]].nGetField(nFieldComp) != oSource[pnIndexSource[nStartSource]].nGetField(nFieldSource)) {
                            printf("Computational error made during comparison.\n");
                            return 1;
                        };
                        pnHKLComp = oComp[pnIndexComp[nStartComp]].pnGetHKL();
                    };
                    
                    if ((!nCompFound) || (!nSourceFound) || ((eType == APPLES_MAC) && ((nSourceFound<nMinRedundancy) || (nCompFound<nMinRedundancy)))) {
                        // Mark all of the reflections in the source as rejected, 
                        if (nSourceFound) {
                            for (nx=nStartSource;nx<=nEndSource;nx++) {
                                if (m_pnReject[nPass][pnIndexSource[nx]]==0) 
                                    m_pnReject[nPass][pnIndexSource[nx]] = nRejectFlag;
                            };
                        };
                        // Mark all reflections in the comparison as rejected.
                        if (nCompFound) {
                            for (nx=nStartComp;nx<=nEndComp;nx++) {
                                if (m_pnReject[!nPass][pnIndexComp[nx]]==0) 
                                    m_pnReject[!nPass][pnIndexComp[nx]] = nRejectFlag;
                            };
                        };
                    } else if (bEqualRedundancy) {
                        if (nCompFound>nSourceFound) {
                            for (nx=0,ny=0;ny<nCompFound-nSourceFound;nx++)
                                if (m_pnReject[!nPass][pnIndexComp[nStartComp+nx]]==0) {
                                    m_pnReject[!nPass][pnIndexComp[nStartComp+nx]] = nRejectFlag;
                                    ny++;
                                };
                        } else if (nCompFound<nSourceFound) {
                            for (nx=0,ny=0;ny<nSourceFound - nCompFound;nx++)
                                if (m_pnReject[nPass][pnIndexSource[nStartSource+nx]]==0) {
                                    m_pnReject[nPass][pnIndexSource[nStartSource+nx]] = nRejectFlag;
                                    ny++;
                                };
                        };
                    };

                    // Do further rejections if we have CCP4 output, AND we want I+/I- data to agree.
                    if ((bDoCCP4Comparision) && (nCompFound) && (nSourceFound)) {
                        int nCompContrib,nSourceContrib;

                        nCompContrib = 0;
                        nSourceContrib = 0;
                        for (nx=nStartComp;nx<=nEndComp;nx++) {
                            if (m_pnReject[!nPass][pnIndexComp[nx]]==0) {
                                if (oComp[pnIndexComp[nx]].fGetField(oComp.m_nFI_fSigmaIPlus)>0.0)
                                    nCompContrib++;
                                if (oComp[pnIndexComp[nx]].fGetField(oComp.m_nFI_fSigmaIMinus)>0.0)
                                    nCompContrib++;
                            };
                        };
                        for (nx=nStartSource;nx<=nEndSource;nx++) {
                            if (m_pnReject[nPass][pnIndexSource[nx]]==0) {
                                if (oSource[pnIndexSource[nx]].fGetField(oSource.m_nFI_fSigmaIPlus)>0.0)
                                    nSourceContrib++;
                                if (oSource[pnIndexSource[nx]].fGetField(oSource.m_nFI_fSigmaIMinus)>0.0)
                                    nSourceContrib++;
                            };
                        };
                        if (nSourceContrib != nCompContrib) {
                            // Delete all contributors since they did not have identical I+ and I-.
                            
                            for (nx=nStartComp;nx<=nEndComp;nx++) 
                                m_pnReject[!nPass][pnIndexComp[nx]]=nRejectFlag;
                            for (nx=nStartSource;nx<=nEndSource;nx++) 
                                m_pnReject[nPass][pnIndexSource[nx]]=nRejectFlag;

                        };
                    };


                    if ((nPass == nUseSigmasIn-1) && (nCompFound) && (nSourceFound)) {
                        double a3fAverageSigmaSource[3];
                        int    a3nAverageSigmaSource[3];
                        double fAverageIntensitySource;
                        double fAverageIntensityComp;
                        int    nAverageIntensityComp;
                        double fScale;
                        double f0;
                        // Use the sigmas in oSource.  Copy them to oComp.
                        for (nx=0;nx<3;nx++) {
                            a3fAverageSigmaSource[nx] = 0.0;
                            a3nAverageSigmaSource[nx] = 0;
                        };
                        fAverageIntensitySource = 0.0;
                        fAverageIntensityComp = 0.0;
                        nAverageIntensityComp = 0;
                        
                        for (nx=nStartSource;nx<=nEndSource;nx++) {
                            if (m_pnReject[nPass][pnIndexSource[nx]]==0) {
                                a3fAverageSigmaSource[0] += oSource[pnIndexSource[nx]].fGetSigmaI();
                                a3nAverageSigmaSource[0]++;
                                if (bDoCCP4Comparision) {
                                    f0 = oSource[pnIndexSource[nx]].fGetField(oSource.m_nFI_fSigmaIPlus);
                                    if (f0>0.0) {
                                        a3fAverageSigmaSource[1] += f0;
                                        a3nAverageSigmaSource[1]++;
                                        
                                    };
                                    
                                    f0 = oSource[pnIndexSource[nx]].fGetField(oSource.m_nFI_fSigmaIMinus);
                                    if (f0>0.0) {
                                        a3fAverageSigmaSource[2] += f0;
                                        a3nAverageSigmaSource[2]++;
                                    };
                                };
                                fAverageIntensitySource += oSource[pnIndexSource[nx]].fGetIntensity();
                            };
                        };

                        if (a3nAverageSigmaSource[0]) {
                            for (nx=nStartComp;nx<=nEndComp;nx++) {                                
                                fAverageIntensityComp += oComp[pnIndexComp[nx]].fGetIntensity();
                                nAverageIntensityComp++;
                            };
                            fScale = (fAverageIntensityComp/nAverageIntensityComp)/(fAverageIntensitySource/a3nAverageSigmaSource[0]);
                            for (nx=nStartComp;nx<=nEndComp;nx++) {
                                oComp[pnIndexComp[nx]].vSetSigmaI(a3fAverageSigmaSource[0]/a3nAverageSigmaSource[0]*fScale);
                                if (bDoCCP4Comparision) {
                                    if ((a3fAverageSigmaSource[1]>0.0) && (oComp[pnIndexComp[nx]].fGetField(oComp.m_nFI_fSigmaIPlus)>0.0)) 
                                        oComp[pnIndexComp[nx]].vSetField(oComp.m_nFI_fSigmaIPlus,(float) (a3fAverageSigmaSource[1]/a3nAverageSigmaSource[1]*fScale));
                                    if ((a3fAverageSigmaSource[2]>0.0) && (oComp[pnIndexComp[nx]].fGetField(oComp.m_nFI_fSigmaIMinus)>0.0)) 
                                        oComp[pnIndexComp[nx]].vSetField(oComp.m_nFI_fSigmaIMinus,(float) (a3fAverageSigmaSource[2]/a3nAverageSigmaSource[2]*fScale));

                                }
                            };
                        };

                    };
                    
                    nStartSource=nRefSort;
                    if (nRefSort<nNumRefs)
                        nLastHKL =  oSource[nRef].nGetField(nFieldSource);                
                };
            };
        };
    } else {
        if (!m_bSilent)
            printf("Will not trim lists in (apples to oranges) comparison!\n");
    };
    if (nRejectFlag == 1) {
        m_bResortRMergeArray[0] = TRUE;
        m_bResortRMergeArray[1] = TRUE;
    };

    if (!m_bSilent) {
        printf("%d Reflections used in first  data set.\n",m_nNumRefs[0] - nCount(0,1));
        printf("%d Reflections used in second data set.\n",m_nNumRefs[1] - nCount(1,1));
    };
    return 0;
};

int Ccompare::nComputeRMerge(int nMinRedundancy,bool bDupSigmas)
{
    int nLastHKL;
    int nThisHKL;
    int nStartSource,nEndSource,nSourceFound0,nSourceFound02;
    int nRefSort;
    int nNumRefs;
    int nRef;

    int nPass;
    int nx,ny;
    double f0;
    static Cstat oStat0;
    static Cstat oStat02;

    bool   bUse02;
    double fGroupAvg;
    double fGroupSigma;
    double fNumer0,fNumer02;
    double fDenom0,fDenom02;
    int nNumContributors;

    if (nMinRedundancy<2)
        nMinRedundancy = 2;

    for (nPass=0;nPass<2;nPass++) {
        Creflnlist& oSource = *m_poList[nPass];
        int* pnIndexSource = m_pnIndex[nPass];
        int  nFieldSource = m_nFI_nIndex[nPass];       

        // If we are recording intensity information, initialize.
        if (m_pfIntensityBinsNumer[nPass]) {
            for (nx=0;nx<m_nIntensityBins;nx++) {
                m_pfIntensityBinsNumer[nPass][nx] = 0.0;
                m_pfIntensityBinsDenom[nPass][nx] = 0.0;
                m_pnIntensityBinsContrib[nPass][nx] = 0;
            };
        };
        

        fNumer0 =0.0;
        fNumer02 = 0.0;
        fDenom0 = 0.0;
        fDenom02 = 0.0;
        nNumContributors=0;
        
        nLastHKL = oSource[pnIndexSource[0]].nGetField(nFieldSource);
        nStartSource = 0;   
        nNumRefs = oSource.nGetNumReflns();
        for (nRefSort=0;nRefSort<nNumRefs+1;nRefSort++) {
            nRef=pnIndexSource[nRefSort];
            
            if (nRefSort<nNumRefs)
                nThisHKL = oSource[nRef].nGetField(nFieldSource);
            if ((nRefSort!=nNumRefs) && (nThisHKL==nLastHKL))
                continue;
            nEndSource=nRefSort-1;

            oStat0.vClear();
            oStat02.vClear();
            // Discover the number of non-rejected source reflections in this group.
            for (nx=nStartSource;nx<=nEndSource;nx++) {
                f0 = max(0.0,(*m_poList[nPass])[pnIndexSource[nx]].fGetIntensity());
                if (m_pnReject[nPass][pnIndexSource[nx]]==0) {
                    oStat0.vAdd(f0);                    
                    oStat02.vAdd(f0);
                } else if (m_pnReject[nPass][pnIndexSource[nx]]==2) {
                    oStat02.vAdd(f0);
                };
            };
            nSourceFound0 = oStat0.nSize();
            nSourceFound02 = oStat02.nSize();

            if (nSourceFound02>=nMinRedundancy) {
                bUse02 = (nSourceFound0<nMinRedundancy);

                fGroupAvg = oStat02.fAverage();
                fGroupSigma = oStat02.fStandardDeviation();
                for (nx=nStartSource;nx<=nEndSource;nx++) {
                    if ((m_pnReject[nPass][pnIndexSource[nx]]==0) || ((m_pnReject[nPass][pnIndexSource[nx]]==2) && (bUse02))) {
                        f0 = max(0.0,(*m_poList[nPass])[pnIndexSource[nx]].fGetIntensity());
                        (*m_poList[nPass])[pnIndexSource[nx]].vSetField(m_nFI_fIndexRMerge[nPass],(float) (fabs(fGroupAvg - f0)/max(1.0,f0)));
                        fNumer02 += fabs(fGroupAvg - f0);
                        fDenom02 += f0;
                    };
                };
            } else {
                for (nx=nStartSource;nx<=nEndSource;nx++) 
                    (*m_poList[nPass])[pnIndexSource[nx]].vSetField(m_nFI_fIndexRMerge[nPass],(float) -1.0);
            };

            // Find the R-Merge for all non-rejected reflections in this group.
            if (nSourceFound0>=nMinRedundancy) {
                fGroupAvg = oStat0.fAverage();
                fGroupSigma = oStat0.fStandardDeviation();
                for (nx=nStartSource;nx<=nEndSource;nx++) {
                    if (m_pnReject[nPass][pnIndexSource[nx]]==0) {
                        f0 = max(0.0,(*m_poList[nPass])[pnIndexSource[nx]].fGetIntensity());
                        fNumer0 += fabs(fGroupAvg - f0);
                        fDenom0 += f0;
                        nNumContributors++;
                        if (m_pfIntensityBinsNumer[nPass]) {
                            // Add log intensity bin information.
                            ny = min(m_nIntensityBins-1,(int) (max(1.0,f0)/m_fIntensityMax*m_nIntensityBins));
                            m_pfIntensityBinsNumer[nPass][ny] += fabs(fGroupAvg - f0);
                            m_pfIntensityBinsDenom[nPass][ny] += f0;
                            m_pnIntensityBinsContrib[nPass][ny] ++;
                        };
                    };
                };
                if (bDupSigmas) {
                    for (nx=nStartSource;nx<=nEndSource;nx++) {
                        (*m_poList[nPass])[pnIndexSource[nx]].vSetSigmaI(fGroupSigma);
                    };
                };
            } 



            
            nStartSource=nRefSort;
            if (nRefSort<nNumRefs)
                nLastHKL =  oSource[nRef].nGetField(nFieldSource);

        };
        m_a2fRMerge[nPass] = (fDenom0==0.0)?0.0:(fNumer0/fDenom0);

        if (!m_bSilent) {
            printf("R-merge for %s data set (%d contributors): %7.4f\n",
                (nPass)?"second":"first",
                nNumContributors,
                m_a2fRMerge[nPass]);            
        };        
    };
    return 0;    
};

int Ccompare::nCount(int nDataSet,int nRejectFlag) {
    int nx;
    int nCount;
    nCount = 0;
    for (nx=0;nx<m_nNumRefs[nDataSet];nx++) {
        if (m_pnReject[nDataSet][nx]==nRejectFlag)
            nCount++;
    };
    return nCount;
};

int Ccompare::nReject(double fPercent,int nDataSet) {
    // First, we must count the number of non-rejected reflections.
    // We reserved field a value of '2' for reflections that were rejected with a call to nReject(),
    // and a value of '1' for reflections that were rejected in nEquivocate()
    // We assume that nComputeRMerge() has already been called.

    int nNumEquivocated;
    int nNumRejected;
    int nMaxRejected;

    // Sort the list.
    if (m_bResortRMergeArray[nDataSet]) {
        m_poList[nDataSet]->vSort(eReflnField_float_type,m_nFI_fIndexRMerge[nDataSet],m_pnIndexRMerge[nDataSet]);
        // m_bResortRMergeArray[nDataSet] = FALSE;
    };
    nNumEquivocated = nCount(nDataSet,1);
    nMaxRejected = (int) (fPercent*(m_nNumRefs[nDataSet] - nNumEquivocated));
    nNumRejected = 0;

    for(int nx=m_nNumRefs[nDataSet]-1;nx>=0;nx--) 
    {
        int nRejectFlag;
        float fRMerge;
        nRejectFlag = m_pnReject[nDataSet][m_pnIndexRMerge[nDataSet][nx]];
        if (nRejectFlag!=1) {
            fRMerge = (*m_poList[nDataSet])[m_pnIndexRMerge[nDataSet][nx]].fGetField(m_nFI_fIndexRMerge[nDataSet]);
            
            if (nNumRejected<nMaxRejected) 
                m_pnReject[nDataSet][m_pnIndexRMerge[nDataSet][nx]] = 2;
            else
                m_pnReject[nDataSet][m_pnIndexRMerge[nDataSet][nx]] = 0;
            nNumRejected++;
        };        
    };
    if (!m_bSilent)
        printf("Rejected %d reflections from %s data set.\n",min(nNumRejected,nMaxRejected),nDataSet?"second":"first");
    return 0;
};


int Ccompare::nWriteNonRejected(int nDataSet,Cstring& sName) {
    int* pnSort;
    int nNumRefsOut;
    int nx;

    pnSort = new int[m_nNumRefs[nDataSet]];

    nNumRefsOut = 0;
    for (nx=0;nx<m_nNumRefs[nDataSet];nx++) {
        if (m_pnReject[nDataSet][nx]==0) 
            pnSort[nNumRefsOut++] = nx;
    };
    if (nNumRefsOut!=0) {
        m_poList[nDataSet]->nWrite(sName,pnSort,NULL,NULL,0,nNumRefsOut-1);
        printf("%d reflections written out to file '%s'\n",nNumRefsOut,sName.string());
    } else 
        printf("0 reflections!  File '%s' not generated\n",sName.string());
    delete[] pnSort;
    return 0;
};

/*  Initializes arrays for intensity bin experiment.

*/
int Ccompare::nInitIntensityBins(double fPercent,int nNumBins)
{
    int nPass;
    double fIntensity;
    double f0;
    int* pnTemp;
    // Sort on intensity.
    fIntensity = 0.0;
    for (nPass=0;nPass<2;nPass++) {
        pnTemp = new int[m_nNumRefs[nPass]];
        m_poList[nPass]->vSort(eReflnField_float_type,0,pnTemp);
        f0 = (*m_poList[nPass])[pnTemp[(int) (fPercent*(m_nNumRefs[nPass] - 1))]].fGetIntensity();
        delete[] pnTemp;
        fIntensity = max(f0,fIntensity);
    };
    m_fIntensityMax = 1.2*fIntensity;
    m_nIntensityBins = nNumBins;
    for (nPass=0;nPass<2;nPass++) {
        if (m_pfIntensityBinsNumer[nPass])
            delete[] m_pfIntensityBinsNumer[nPass];
        m_pfIntensityBinsNumer[nPass] = new double[nNumBins];
        if (m_pfIntensityBinsDenom[nPass])
            delete[] m_pfIntensityBinsDenom[nPass];
        m_pfIntensityBinsDenom[nPass] = new double[nNumBins];
        if (m_pnIntensityBinsContrib[nPass])
            delete[] m_pnIntensityBinsContrib[nPass];
        m_pnIntensityBinsContrib[nPass] = new int[nNumBins];
    };
    return 0;
};

int Ccompare::nWriteIntensityBins(Cstring& sFileOut)
{
    FILE* pFOut;
    int nx;

    pFOut = fopen(sFileOut.string(),"w+t");
    if (!pFOut) {
        printf("Could not open file '%s'\n",sFileOut.string());
        return 1;
    };
    for (nx=0;nx<m_nIntensityBins;nx++) {
        if ((m_pfIntensityBinsDenom[0][nx]!=0.0) && (m_pnIntensityBinsContrib[0][nx]>10) &&
            (m_pfIntensityBinsNumer[1][nx]!=0.0) && (m_pnIntensityBinsContrib[1][nx]>10)) {
            fprintf(pFOut,"%f %f %f\n",
                (nx+0.5)*m_fIntensityMax/m_nIntensityBins,
                m_pfIntensityBinsNumer[0][nx]/m_pfIntensityBinsDenom[0][nx],
                m_pfIntensityBinsNumer[1][nx]/m_pfIntensityBinsDenom[1][nx]
                );
        };
    };
    fclose(pFOut);
    return 0;
};
