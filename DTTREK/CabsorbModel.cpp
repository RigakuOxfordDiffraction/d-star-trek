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
// dtscaleaverage.cc     Initial author: T.J. Niemeyer  Spring 2001
//             Based on dtscalemerge 
//    This is a new royalty-less absorption algorithm 
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



#include "Cabsorb.h"


Cabsorbmodel::Cabsorbmodel(Creflnlist& oListIn,int* pnHKLSort,int* pnReject):m_oListIn(oListIn) {   
    int nRef;
    
    m_nMinRedundancy = 2;
    m_nMaxFoundRedundancy = 0;
    m_nDegreesPerS0 = 10;
    m_fMinCoeffRange = 0.6;
    m_fMaxCoeffRange = 1.0/m_fMinCoeffRange;    
        
    m_nStat = 1;
    m_pfSS0Coeffs = NULL;
    m_pfScale = NULL;
    
    m_pnIndex = pnHKLSort;
    m_pnReject = pnReject;
    m_pfScale  = new double[m_oListIn.nGetNumReflns()];     
    // Set all scale factors to 1.0;
    for (nRef=0;nRef<m_oListIn.nGetNumReflns();nRef++) {
        m_pfScale[nRef] = 1.0;
    };
    return;
};

int Cabsorbmodel::nInit(Cscaleaverage& oScale) {
    int nx;
    int nMode;

    
    m_nStat = 1;
    if (nInitReflectionList(m_oListIn))
        return 1;
      
    fComputeRMerge();
    for (nx=180-1;(nx>=0) && (an2ThetaPoints[nx]==0);nx--);
    if (nx)
        m_oSurface.m_nMax2Theta =  nx+2;   

    m_oSurface.m_nDegreesPerS0 = m_nDegreesPerS0;
    m_nStat  = m_oSurface.nLoad(m_oListIn,m_a3nFI_fSVec,m_nFI_nBatchIdx,m_a3nFI_fS0Vec,oScale.poGetBatchInfo());
    if (m_nStat)
        return 1;

    printf("Building S vector mesh.\n");    
    nMode = TOTAL_MODES;

    do {
        nMode--;
        if (nMode==-1)
            break;
        
        m_oSurface.nBuildPoints("S-Vector surface",nMode); 
        if (m_oSurface.nBuildTriangles()) 
            m_nStat++;
        else
            m_nStat += m_oSurface.nFindTriangleCoeffs(m_oListIn);
        if (m_nStat)
            break;        
    } while (m_oSurface.nCheckMeshDimensions(m_oListIn,m_pnIndex,m_pnReject,m_nMinRedundancy));

//    Cstring sTemp;
//    m_oSurface.nPrint(0,sTemp = "k:\\tjn\\sphered.dat","test",TRUE);
    
    if (m_nStat) {
        printf("ERROR:  Unable to determine an appropriate mesh!\n");
        return 1;
    };

    if (m_pfSS0Coeffs)
        delete[] m_pfSS0Coeffs;
    m_pfSS0Coeffs = new double[m_oSurface.m_nNumSPoints*m_oSurface.m_nNumS0Points];

    for (nx=0;nx<m_oSurface.m_nNumSPoints*m_oSurface.m_nNumS0Points;nx++) {
        m_pfSS0Coeffs[nx] = 0.0;
    };
       
    return m_nStat;
};


   
    
 
Cabsorbmodel::~Cabsorbmodel() {
    if (m_pfSS0Coeffs)
        delete[] m_pfSS0Coeffs;
    if (m_pfScale)
        delete[] m_pfScale;
};


int Cabsorbmodel::nInitReflectionList(Creflnlist& oList) {
    Cstring sTemp;
    int nx;
    int nContributors;
    double fContributorsIoverSig;

    // Find the indices for the S0 vector and S vector.
    
    m_a3nFI_fS0Vec[0] = oList.nGetFieldIndex(sTemp = D_K_IntegratefS0vec0);
    m_a3nFI_fS0Vec[1] = oList.nGetFieldIndex(sTemp = D_K_IntegratefS0vec1);
    m_a3nFI_fS0Vec[2] = oList.nGetFieldIndex(sTemp = D_K_IntegratefS0vec2);
    m_a3nFI_fSVec[0] = oList.nGetFieldIndex(sTemp = D_K_IntegratefSvec0);
    m_a3nFI_fSVec[1] = oList.nGetFieldIndex(sTemp = D_K_IntegratefSvec1);
    m_a3nFI_fSVec[2] = oList.nGetFieldIndex(sTemp = D_K_IntegratefSvec2);
    m_nFI_fScale = oList.nExpandGetField(sTemp = "fScale");
    m_nFI_nPackedHKL = oList.nGetFieldIndex(oList.ms_snPackedHKL);
    m_nFI_sBatch = oList.nGetFieldIndex(oList.ms_ssBatch);
    m_nFI_nBatchIdx = oList.nGetFieldIndex(Cscaleaverage::ms_snBatchIdx);
    m_nFI_fSigmaI = oList.nGetFieldIndex(sTemp = "fCalcSigmaI");
    if (m_nFI_fSigmaI<0)
        m_nFI_fSigmaI = 1;

    if (m_nFI_sBatch==-1) {
        printf("WARNING: Could not find batch number.\n");
    };
    
    for (nx=0;nx<3;nx++) {
        if (m_a3nFI_fSVec[nx]<0) 
            break;
        if (m_a3nFI_fS0Vec[nx]<0)
            break;
    };
    if (m_nFI_nPackedHKL<0) {
        printf("Could not find Packed HKL entry.\n");
        return 1;
    };
    if (nx!=3) {
        printf("Could not find all S and S0 field entries.\n");
        return 1;
    };

    return 0;
};



int Cabsorbmodel::nBuildReflectionScaleFactors() {
    int nRef;
    int nCoeff;
    int nx,ny,nz;
    double f0,f1;
    double fAbsorbance;
    double fWeight;
    double fFactor;
    double fMaxScale = -1e10;
    double fMinScale = 1e10;

    
    // Next, apply these coefficients to all non-excluded reflections.
    for (nRef=0;nRef<m_oListIn.nGetNumReflns();nRef++) {
                
        // Remember that the fast direction encodes S points.
        // We are taking a weighted average of 6 coefficients for each reflection.
        fAbsorbance = 0.0;
        fWeight = 0.0;
        
        for (nx=0;nx<3;nx++) {
            nCoeff = m_oSurface.m_a3pnSPoints[nx][nRef];
            
            for (ny=0;ny<2;ny++) {                
                f0 = m_oSurface.m_a2pfS0Coeffs[ny][nRef]*m_oSurface.m_a3pfSCoeffs[nx][nRef];
                f1 = m_pfSS0Coeffs[m_oSurface.m_a2pnS0Points[ny][nRef]*m_oSurface.m_nNumSPoints +
                    m_oSurface.m_a3pnSPoints[nx][nRef]];
                if (f1<0.0)
                    nx=nx;
                fAbsorbance += f0*f1;
                fWeight += f0;
            };
        };
        if (nx==3) {
            fFactor = (fAbsorbance/=fWeight);
            fMaxScale = max(fMaxScale,fFactor);
            fMinScale = min(fMinScale,fFactor);
            
            if (fAbsorbance>0.0) {
                m_pfScale[nRef] = fFactor;
            } 
        };
        
            
    };
    return 0;
};

double Cabsorbmodel::fComputeRMerge() {
    int nx,ny;
    double f0,f1;

    int nLastHKL;
    int nThisHKL;
    int nNumRefs;
    int nRef;
    int nRefSort;
    int nStart;
    int nEnd;
    int nFound;
    int nPass;
    static Cstat oStat;


    double a3fVec[3];
    double fNumer0;
    double fDenom0;
    double fGroupNumer0;
    double fGroupDenom0;
    double fGroupAvg;
    double fGroupSigma;
    double f2Theta;
    int nNumGroupContributors;

    for (nx=0;nx<180;nx++) 
        an2ThetaPoints[nx] = 0;

    nNumRefs = m_oListIn.nGetNumReflns();
    
    fNumer0 = 0.0;
    fDenom0 = 0.0;
    m_nContributors = 0;
    nLastHKL = m_oListIn[m_pnIndex[0]].nGetField(m_nFI_nPackedHKL);    
    m_nMaxFoundRedundancy = 0;
    nStart=0;
    
    for (nRefSort=0;nRefSort<nNumRefs+1;nRefSort++) {
        nRef=m_pnIndex[nRefSort];        
        if (nRefSort<nNumRefs)
            nThisHKL = m_oListIn[nRef].nGetField(m_nFI_nPackedHKL);
        if ((nRefSort!=nNumRefs) && (nThisHKL==nLastHKL))
            continue;
        nEnd=nRefSort-1;

           
        oStat.vClear();
        
        // Discover the number of non-rejected source reflections in this group.
        for (nx=nStart;nx<=nEnd;nx++) {
            f0 = max(0.0,m_oListIn[m_pnIndex[nx]].fGetIntensity());
            if (m_pnReject[nx]==0) {
                oStat.vAdd(f0*fGetScale(m_pnIndex[nx]));
            };
        };


        nx = nStart;
        for (ny=0;ny<3;ny++) {
            a3fVec[ny] = m_oListIn[m_pnIndex[nx]].fGetField(m_a3nFI_fS0Vec[ny]) +
                m_oListIn[m_pnIndex[nx]].fGetField(m_a3nFI_fSVec[ny]);
        };
        f2Theta = acos(max(-1.0,min(1.0,1.0 - 0.5*fDot3D(a3fVec,a3fVec))))/Gs_dRADIANS_PER_DEGREE;
        an2ThetaPoints[max(0,min(180-1,(int) f2Theta))]+=oStat.nSize();

        
        if (oStat.nSize()>=m_nMinRedundancy) {
            
            
           
            nFound = 0;
            for (nx=0;nx<oStat.nSize();nx++) {
                if (oStat.bIsMarked(nx))
                    nFound++;
            };
            
            // Find the R-Merge for all non-rejected reflections in this group.
            
            fGroupAvg = oStat.fAverage();
            fGroupSigma = oStat.fStandardDeviation();
            nNumGroupContributors = 0;
            fGroupNumer0 = 0.0;
            fGroupDenom0 = 0.0;
            for (nx=nStart;nx<=nEnd;nx++) {
                if (oStat.bIsMarked(nx-nStart)) {
                    
                    f0 = max(0.0,m_oListIn[m_pnIndex[nx]].fGetIntensity())*fGetScale(m_pnIndex[nx]);
                    fGroupNumer0 += fabs(fGroupAvg - f0);
                    fGroupDenom0 += f0;
                    nNumGroupContributors++;
                };
            };
        } else
            nFound = 0;
        
        if (nFound>=m_nMinRedundancy) {
            fNumer0 += fGroupNumer0;
            fDenom0 += fGroupDenom0;
            m_nContributors+=nNumGroupContributors;
        } 
        m_nMaxFoundRedundancy = max(nFound,m_nMaxFoundRedundancy);
        
        
        nStart=nRefSort;
        if (nRefSort<nNumRefs)
            nLastHKL =  m_oListIn[nRef].nGetField(m_nFI_nPackedHKL);        
    };

    

    // Clear all '2' flags from the last pass.
    for (nx=0;nx<nNumRefs;nx++) {
        if (m_pnReject[nx]==2)
            m_pnReject[nx] = 0;
    };


    if (m_nContributors==0)
        return 0.0;
    else
        return fNumer0/fDenom0;
};


int Cabsorbmodel::nApply() {
    int nRef;
    int nx;
    double f0;
    for (nRef=0;nRef<m_oListIn.nGetNumReflns();nRef++) {
        f0 = fGetScale(nRef);
        m_oListIn[nRef].vSetIntensity(m_oListIn[nRef].fGetIntensity()*f0);
        m_oListIn[nRef].vSetSigmaI(m_oListIn[nRef].fGetSigmaI()*f0);
    };
    return 0;
};

double Cabsorbmodel::fGetScale(int nRef) {
    return m_pfScale[nRef];
};


int Cabsorbmodel::nRelativeRingScale() {
    int nS0Surface;
    int nPoint;
    int nStart,nEnd;
    int nRing;
    int nx;
    int nPass;
    double f0,f1;
    Cstat oStat;  
    
    for (nRing=0;1;nRing++) {
        m_oSurface.nGetRingStartEnd(nRing,nStart,nEnd);
        if ((nStart==-1) || (nRing==-1))
            break;
        for (nPass=0;nPass<2;nPass++) {
            oStat.vClear();
            nx  = 0;
            for (nS0Surface = 0; nS0Surface< m_oSurface.m_nNumS0Points;nS0Surface++) {            
                for (nPoint=nStart;(nPoint<=nEnd);nPoint++) {
                    f0=m_pfSS0Coeffs[nS0Surface*m_oSurface.m_nNumSPoints + nPoint];
                    if (((nPass==0) && (f0>0.0)) ||
                        ((nPass==1) && (f0<=m_fMaxCoeffRange) && (f0>=m_fMinCoeffRange))) {
                        oStat.vAdd(f0);
                        nx++;
                    }                         
                };
            };                        
            nx -= oStat.nMarkStat(FALSE,Cstat::S_ABOVEBELOW,5.0);
            if (nx) {
                f0 = oStat.fAverage();
                for (nS0Surface = 0; nS0Surface< m_oSurface.m_nNumS0Points;nS0Surface++) {            
                    for (nPoint=nStart;nPoint<=nEnd;nPoint++) {
                        f1 = m_pfSS0Coeffs[nS0Surface*m_oSurface.m_nNumSPoints + nPoint];
                        if (((nPass==0) && (f1>0.0)) ||
                            ((nPass==1) && (f1<=m_fMaxCoeffRange) && (f1>=m_fMinCoeffRange))) 
                            m_pfSS0Coeffs[nS0Surface*m_oSurface.m_nNumSPoints + nPoint]/=f0;
                        else
                            m_pfSS0Coeffs[nS0Surface*m_oSurface.m_nNumSPoints + nPoint] = -1.0;
                    };
                };
            };
        };
    };
    return 0;
};

// All S point coeffs that do not have connections, are set equal to the average of neighboring S point coeffs.
int Cabsorbmodel::nPropogateRingScale() {
    int nCoeff1,nCoeff2;
    int nS0Point1,nS0Point2;
    int nCoeff1i,nCoeff2i;
    int nAdjCt;
    int nStart,nEnd;
    int nPass;
    int nIter;
    int nMaxConnects;
    double fAverage;
    int    nAverage;
    
    nMaxConnects = 1;
    for (nIter=0;nMaxConnects;nIter++) {
        nMaxConnects = 0;
        for (nPass=0;nPass<2;nPass++) {
            for (nCoeff1 = 0; nCoeff1<m_oSurface.m_nNumSPoints;nCoeff1++) {
                for (nS0Point1=0;nS0Point1< m_oSurface.m_nNumS0Points;nS0Point1++) {
                    if (m_pfSS0Coeffs[nCoeff1i = (nS0Point1*m_oSurface.m_nNumSPoints+ nCoeff1)]<0.0) {
                        fAverage = 0.0;
                        nAverage = 0;
                        for (nAdjCt=0;nAdjCt<m_oSurface.m_poSPoints[nCoeff1].nAdjCt;nAdjCt++) {
                            nCoeff2 = m_oSurface.m_poSPoints[nCoeff1].anAdj[nAdjCt];
                            for (nS0Point2 = max(0,nS0Point1-1); nS0Point2<min(m_oSurface.m_nNumS0Points,nS0Point1+1);nS0Point2++) {
                                if (m_oSurface.m_poS0Points[nS0Point2].nOrient==m_oSurface.m_poS0Points[nS0Point1].nOrient) {
                                    if ((nCoeff2!=nCoeff1) || (nS0Point2!=nS0Point1)) {
                                        if (m_pfSS0Coeffs[nCoeff2i = (nS0Point2*m_oSurface.m_nNumSPoints + nCoeff2)]>0.0) {
                                            nAverage++;
                                            fAverage += m_pfSS0Coeffs[nCoeff2i];
                                        };
                                    };
                                };
                            };
                        };
                        // In the first pass, we are only interested in finding the points that have the most neighbors.
                        if (nPass==0)
                            nMaxConnects = max(nMaxConnects,nAverage);
                        else if ((nAverage) && (nAverage==nMaxConnects))
                            m_pfSS0Coeffs[nCoeff1i] = fAverage/nAverage;
                    };
                };
            };    
        };
    };
    return 0;
};

int Cabsorbmodel::nPrintRingScale(bool bPrintCounts) {
    int nx,ny;
    int nRing;
    int nStart,nEnd;
    double f0,f1;
    int nS0Point;
    int nSPoint;
    int nS0PointStart;
    const int nWidth = 6;

    printf("\n");
    for (nRing=1;1;nRing++) {
        m_oSurface.nGetRingStartEnd(nRing,nStart,nEnd);
        if ((nStart==-1) || (nRing==-1))
            break;
        
        printf("%s for ring %d (2theta = %5.1lf)\n",(bPrintCounts)?"Contributor Counts":"Coefficients",nRing,m_oSurface.m_poSPoints[nStart].f2Theta);
        printf("---------------------------------------------------------------------\n");
        nS0PointStart = 0;
        do {
            if (nS0PointStart!=0)
                printf("\n\n");
            printf("S0 Rotation/scan ");
            for (nS0Point=nS0PointStart;nS0Point< min(nWidth+nS0PointStart,m_oSurface.m_nNumS0Points);nS0Point++) {
                printf("%4.0lf.0/%d ",m_oSurface.m_poS0Points[nS0Point].fRot,m_oSurface.m_poS0Points[nS0Point].nOrient+1);
            };
            for (nSPoint=nStart;nSPoint<=nEnd;nSPoint++) {
                printf("\n");
                printf("S  Chi     %3.0lf.0   ",m_oSurface.m_poSPoints[nSPoint].fChi);
                for (nS0Point=nS0PointStart;nS0Point< min(nWidth+nS0PointStart,m_oSurface.m_nNumS0Points);nS0Point++) {
                    if (bPrintCounts) 
                        printf("%6d   ",m_oSurface.m_pnC[nS0Point*m_oSurface.m_nNumSPoints + nSPoint]);
                    else                        
                        printf("%6.3lf   ",m_pfSS0Coeffs[nS0Point*m_oSurface.m_nNumSPoints + nSPoint]);
                };
            };
            nS0PointStart+=nWidth;
        } while (nS0PointStart<m_oSurface.m_nNumS0Points);
        printf("\n");
        printf("----------------------------------------------------------------------\n");
        printf("\n");
    };

    return 0;
};
    

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int Cabsorbmodel::nLinearMethod(int nResoRing) {
    int nx,ny,nz,nw;
    double f0;
    double f1;

    int nSS0Coeffs;
    int nStat;

    double fDeviation;
    double fOrigR;
    double fFinalR; 
    double* pfA;
    double* pfb;
    double* pfX;
    int*    pnC;
    int*    pnUseHarmonic;

    static Cstat oStat;
    int*    pnCoeffAverageChi;
    double* pfCoeffAverageChi;
    double* pfCoeffBuffer;
    double* pfCoeffs;
    double* pfCoeffs2;
    int*    pnCoeffBuffer;
    int*    pnStatus;
    int     nCoeffBuffer;
    int     nHarmonicTerms;

    int     nCoeffs;
    int     nCoeff1,nCoeff2;
    int     nRef2,nRef1;
    double  fCumulScale;
    int     nMappedCoeff;
    int     nUnMappedCoeff;
    int     nRingStart,nRingEnd,nRingCoeffs;


    int nRefSort;
    int nLastHKL;
    int nThisHKL;
    int nNumRefs;
    int nRef;
    int nStart;
    int nEnd;
    int nFound;
    int nPass;
    bool bPrint = FALSE;


    nNumRefs = m_oListIn.nGetNumReflns();

   
    if (nResoRing==-1) {
        // If the user called this function, then all resolution rings must be fit.
        fOrigR = fComputeRMerge();

        // Set all coefficients to -1.0 for now.
        for (nx=0;nx<m_oSurface.m_nNumS0Points*m_oSurface.m_nNumSPoints;nx++)
            m_pfSS0Coeffs[nx] = -1.0;

        // Run the scaling on different rings.
        for (nx=1;(0==nLinearMethod(nx));nx++);


        nRelativeRingScale();        
        nPropogateRingScale();
        nPrintRingScale(FALSE);
        nPrintRingScale(TRUE);

        nBuildReflectionScaleFactors();
        fFinalR = fComputeRMerge();

        printf("Linear Absorption Info\n");
        printf("--------------------------------------\n");
        printf("Total Reflections       %6d\n",nNumRefs);
        printf("Total Contributors      %6d\n",m_nContributors);
        printf("Original R-merge        %6.3f\n",fOrigR);
        printf("Final   R-merge         %6.3f\n",fFinalR);
        printf("Total S points in Mesh  %6d\n",m_oSurface.m_nNumSPoints);

        for (nx=1,nw=1;m_oSurface.nGetRingStartEnd(nx,ny,nz),ny!=-1;nx++) {
            if (nz-ny+1>10)
                nw = nw*100+(nz-ny+1);
            else
                nw = nw*10 + (nz-ny+1);
        };
        printf("S surface topology      %6d\n",nw);
        printf("Total S0 segments       %6d\n",m_oSurface.m_nNumS0Points);
        printf("Degrees per S0          %6d\n",m_nDegreesPerS0);
        printf("--------------------------------------\n");
        return 0;
    } else {
        m_oSurface.nGetRingStartEnd(nResoRing,nRingStart,nRingEnd);
        if ((nRingStart==-1) || (nRingEnd==-1))
            return 1;
        nRingCoeffs = nRingEnd - nRingStart + 1;
        if (m_nContributors==0)
            return 1;
    };

    nSS0Coeffs = nRingCoeffs * m_oSurface.m_nNumS0Points;
    pnStatus      = new int[nSS0Coeffs];    

    // Initialize matrices.
    pfA = new double[nSS0Coeffs*nSS0Coeffs];
    pfb = new double[nSS0Coeffs];
    pfX = new double[nSS0Coeffs];
    pnC = new int[nSS0Coeffs];
    pfCoeffs = new double[nSS0Coeffs];
    pfCoeffs2= new double[m_oSurface.m_nNumSPoints];
    pfCoeffBuffer = new double[m_nMaxFoundRedundancy*16];
    pnCoeffBuffer = new int[m_nMaxFoundRedundancy*16];
    pfCoeffAverageChi = new double[nSS0Coeffs];
    pnCoeffAverageChi = new int[nSS0Coeffs];
    pnUseHarmonic = new int[m_oSurface.m_nNumS0Points];

    for (nPass=-1;nPass<=-1;nPass++) {
        vZeroMat(nSS0Coeffs,nSS0Coeffs,pfA);
        vZeroMat(nSS0Coeffs,1,pfb);   
        vZeroMat(nSS0Coeffs,1,pfCoeffAverageChi);
        for (nx=0;nx<nSS0Coeffs;nx++) {
            pnC[nx] = 0;
            pnCoeffAverageChi[nx] = 0;
        };
        
        for (nx=0;nx<m_oSurface.m_nNumS0Points;nx++) {
            pnUseHarmonic[nx] = 0;

            for (ny=nRingStart;ny<=nRingEnd;ny++) {
                if (m_oSurface.m_pnC[nx*m_oSurface.m_nNumSPoints + ny]<m_oSurface.m_fMinContrib)
                    pnUseHarmonic[nx] = 3;
            };
        };

        
        nLastHKL = m_oListIn[m_pnIndex[0]].nGetField(m_nFI_nPackedHKL);   
        nStart = 0;
        for (nRefSort=0;nRefSort<nNumRefs+1;nRefSort++) {
            
            nRef=m_pnIndex[nRefSort];
            
            if (nRefSort<nNumRefs)
                nThisHKL = m_oListIn[nRef].nGetField(m_nFI_nPackedHKL);
            if ((nRefSort!=nNumRefs) && (nThisHKL==nLastHKL))
                continue;
            nEnd=nRefSort-1;
            
            nFound = 0;
            oStat.vClear();
            for (nx=nStart;nx<=nEnd;nx++) {
                int nSPoint0 = m_oSurface.m_a2pnClosestS[0][m_pnIndex[nx]];
                int nSPoint1 = m_oSurface.m_a2pnClosestS[1][m_pnIndex[nx]];
                if (((nSPoint0>=nRingStart) && (nSPoint0<=nRingEnd)) ||
                    ((nSPoint1>=nRingStart) && (nSPoint1<=nRingEnd))) {
                    f0 = max(0.0,m_oListIn[m_pnIndex[nx]].fGetIntensity());
                    oStat.vAdd(f0*m_pfScale[m_pnIndex[nx]]);
                    if (m_pnReject[nx]==0)
                        nFound++;
                    else
                        oStat.bMark(FALSE,nx-nStart);
                } else {
                    oStat.vAdd(0.0);
                    oStat.bMark(FALSE,nx-nStart);
                };
            };
            // Find out how many remain.

            nFound -= oStat.nMarkStat(FALSE,Cstat::S_ABOVEBELOW,5.0);
            
            if (nFound>=m_nMinRedundancy) {
                // Add coefficients to normal matrix and vector:
                
                for (nRef1=nStart;nRef1<=nEnd;nRef1++) {
                    if (oStat.bIsMarked(nRef1-nStart)) {
                        fDeviation = 0.0;
                        // No coefficients to begin with.
                        nCoeffBuffer = 0;                                    
                        for (nRef2=nStart;nRef2<=nEnd;nRef2++) {
                            if (oStat.bIsMarked(nRef2-nStart)) {
                                double fCoeff;
                                int    nCoeffS;
                                int    nCoeffS0;
                                int    nCoeffSS0;

                                if (nRef1==nRef2)
                                    f0 = 1.0 - 1.0/nFound;
                                else
                                    f0 = -1.0/nFound;
                                fCoeff = f0*max(0.0,m_oListIn[m_pnIndex[nRef2]].fGetIntensity());

                                nCoeffS0 = m_oSurface.m_pnClosestS0[m_pnIndex[nRef2]];
                                if ((m_oSurface.m_a2pnClosestS[0][m_pnIndex[nRef2]]<=nRingEnd) &&
                                    (m_oSurface.m_a2pnClosestS[0][m_pnIndex[nRef2]]>=nRingStart))
                                    nCoeffS = m_oSurface.m_a2pnClosestS[0][m_pnIndex[nRef2]];
                                else if ((m_oSurface.m_a2pnClosestS[1][m_pnIndex[nRef2]]<=nRingEnd) &&
                                    (m_oSurface.m_a2pnClosestS[1][m_pnIndex[nRef2]]>=nRingStart))
                                    nCoeffS = m_oSurface.m_a2pnClosestS[1][m_pnIndex[nRef2]];
                                else {
                                    printf("WARNING:  Error in linear absorption!\n");
                                    nCoeffS = m_oSurface.m_a2pnClosestS[0][m_pnIndex[nRef2]];
                                };        
                                
                                f0 = m_oSurface.m_pfChi[m_pnIndex[nRef2]];
                                nCoeffSS0 = nCoeffS0* nRingCoeffs + nCoeffS - nRingStart;
                                if (pnCoeffAverageChi[nCoeffSS0]) {
                                    f1 = pfCoeffAverageChi[nCoeffSS0]/pnCoeffAverageChi[nCoeffSS0];
                                    if (fabs(f0-f1)>fabs(f0 + Gs_dPI*2.0 -f1))
                                        f0 += Gs_fPI*2.0;
                                    if (fabs(f0-f1)>fabs(f0 - Gs_dPI*2.0 - f1))
                                        f0 -= Gs_dPI*2.0;
                                };
                                pnCoeffAverageChi[nCoeffSS0]++;
                                pfCoeffAverageChi[nCoeffSS0] += f0;

                                
                                if (pnUseHarmonic[nCoeffS0]) {
                                    nHarmonicTerms = min(pnUseHarmonic[nCoeffS0],nRingCoeffs);
                                    // We use as many harmonics as are available as S points.
                                    for (nz=0;nz<nHarmonicTerms;nz++) {
                                        if (nz<=nHarmonicTerms/2)
                                            f0 = cos(nz*m_oSurface.m_pfChi[m_pnIndex[nRef2]]);
                                        else
                                            f0 = sin((nz-nHarmonicTerms/2)*m_oSurface.m_pfChi[m_pnIndex[nRef2]]);
                                        nCoeffSS0 = nCoeffS0* nRingCoeffs + nz;

                                        for (nx=0;nx<nCoeffBuffer;nx++) {
                                            if (pnCoeffBuffer[nx]==nCoeffSS0)
                                                break;
                                        };
                                        if (nx==nCoeffBuffer) {
                                            nCoeffBuffer++;
                                            pnCoeffBuffer[nx] = nCoeffSS0;
                                            pfCoeffBuffer[nx] = 0.0;
                                        };
                                        pfCoeffBuffer[nx] += fCoeff*f0;
                                    };
                                } else {
                                    // Search for the coefficient.
                                    nCoeffSS0 = nCoeffS0* nRingCoeffs + nCoeffS - nRingStart;
                                    for (nx=0;nx<nCoeffBuffer;nx++) {
                                        if (pnCoeffBuffer[nx]==nCoeffSS0)
                                            break;
                                    };
                                    if (nx==nCoeffBuffer) {
                                        nCoeffBuffer++;
                                        pnCoeffBuffer[nx] = nCoeffSS0;
                                        pfCoeffBuffer[nx] = 0.0;
                                    };
                                    pfCoeffBuffer[nx] += fCoeff;
                                };

                                f0 = m_oListIn[m_pnIndex[nRef2]].fGetField(m_nFI_fSigmaI);
                                fDeviation += f0*f0;
                            };
                        };                                                            
                        
                        fDeviation = 1.0;
                        if ((fDeviation>0.0) && (nCoeffBuffer !=0)) {
                            fDeviation = sqrt(fDeviation/nCoeffBuffer);
                            // We have an equation describing the difference between nRef1 and the 
                            // Average of the other reflections. Now add these coeffs. to the normal equations.
                            for (nx=0;nx<nCoeffBuffer;nx++) {
                                nCoeff1 = pnCoeffBuffer[nx];
                                pnC[nCoeff1]++;
                                for (ny=0;ny<nCoeffBuffer;ny++) {                                    
                                    f0 = pfCoeffBuffer[nx]*pfCoeffBuffer[ny]/(fDeviation*fDeviation);
                                    nCoeff2 = pnCoeffBuffer[ny];                                    
                                    pfA[nCoeff1*nSS0Coeffs + nCoeff2] += f0;
                                };                            
                            };
                        };
                    };
                };
            };
            
            nStart=nRefSort;
            if (nRefSort<nNumRefs)
                nLastHKL =  m_oListIn[nRef].nGetField(m_nFI_nPackedHKL);        
        };
        
        nSmartMatrixCompress(nSS0Coeffs,pfA,pfb,pnStatus,TRUE,NULL,0.000);               
        nStat = nSolveAxB_svd(nSS0Coeffs,pfA,pfX,pfb);
        
        for (nx=0,ny=0;nx<nRingCoeffs * m_oSurface.m_nNumS0Points;nx++) {
            if (pnStatus[nx]==1)
                pfCoeffs[nx] = pfX[ny++];
            else if (pnC[nx]==0)
                pfCoeffs[nx] = -1.0;
            else
                pfCoeffs[nx] = 1.0;
        };
        if (bPrint) {
            printf("\n");
            for (ny=0;ny<m_oSurface.m_nNumS0Points;ny++) {
                printf("\n");
                for (nx=0;nx<nRingCoeffs;nx++) {
                    nz = ny*nRingCoeffs + nx;
                    printf("%2d %6d %10f\n",pnStatus[nz],pnC[nz],pfCoeffs[nz]);
                };
            };
        };
        


        for (ny=0;ny<m_oSurface.m_nNumS0Points;ny++) {
            if (pnUseHarmonic[ny]) {
                nHarmonicTerms = min(pnUseHarmonic[ny],nRingCoeffs);
                // If we are using the harmonics, these must be translated.
                for (nCoeff1=0;nCoeff1<nRingCoeffs;nCoeff1++) {
                    if (pnCoeffAverageChi[ny*nRingCoeffs + nCoeff1]) {
                        // Use the average chi angle for reflections that had this coefficient.
                        f0 = 0.0;
                        f1 = pfCoeffAverageChi[ny*nRingCoeffs + nCoeff1]/pnCoeffAverageChi[ny*nRingCoeffs + nCoeff1];
                        
                        for (nCoeff2=0;nCoeff2<nHarmonicTerms;nCoeff2++) {
                            if (nCoeff2<=nHarmonicTerms/2)
                                f0 += pfCoeffs[ny*nRingCoeffs + nCoeff2]*cos(nCoeff2*f1);
                            else
                                f0 += pfCoeffs[ny*nRingCoeffs + nCoeff2]*sin((nCoeff2-nHarmonicTerms/2)*f1);
                        };
                        pfCoeffs2[nCoeff1] = f0;
                    } else
                        pfCoeffs2[nCoeff1] = -1.0;
                };
                for (nCoeff1=0;nCoeff1<nRingCoeffs;nCoeff1++) {
                    pfCoeffs[ny*nRingCoeffs + nCoeff1] = pfCoeffs2[nCoeff1];
                };
            } 
        };
        if (bPrint) {
            printf("\n");
            for (ny=0;ny<m_oSurface.m_nNumS0Points;ny++) {
                printf("\n");
                for (nx=0;nx<nRingCoeffs;nx++) {
                    nz = ny*nRingCoeffs + nx;
                    printf("%2d %6d %10f\n",pnStatus[nz],pnC[nz],pfCoeffs[nz]);
                };
            };
        };        
        
    };

    // Put the completed scale factors in the correct places.
    for (ny=0;ny<m_oSurface.m_nNumS0Points;ny++) {
        for (nx=0;nx<nRingCoeffs;nx++) {                    
            f0 = pfCoeffs[ny*nRingCoeffs + nx];
            m_pfSS0Coeffs[ny*m_oSurface.m_nNumSPoints+ nx + nRingStart] = f0;
        };
    };
        

    delete[] pfA;
    delete[] pfb;
    delete[] pfX;
    delete[] pnC;
    delete[] pfCoeffs;
    delete[] pfCoeffs2;
    delete[] pfCoeffAverageChi;
    delete[] pnCoeffAverageChi;
    delete[] pfCoeffBuffer;
    delete[] pnCoeffBuffer;
    delete[] pnStatus;
    delete[] pnUseHarmonic;
    return 0;
};


