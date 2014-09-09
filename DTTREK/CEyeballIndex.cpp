//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CEyeballIndex.h            Initial author: Thaddeus Niemeyer 24-August-01
//  This file contains the member functions of class CEyeballIndex
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


#include "CEyeballIndex.h"

#include "Cdetector.h"
#include "Csource.h"
#include "Crotation.h"

#include "memory.h"
#ifdef DOSEM_STATUSBAR
#include "Crapid.h"
#include "dosem.h"
#endif
#include "string.h"

#include "dtsvd.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::cerr;
#endif

CEyeballIndex::CEyeballIndex(Cimage_header& oHeader,Creflnlist& oList):m_oList(oList),m_oHeader(oHeader) {
    double f0;
    m_nNumRefs = m_oList.nGetNumReflns();
    m_bIsAvailable = FALSE;
    m_nPrint = 0;

    m_pn3DGridSums = NULL;
    m_pnPrintNeighbors = NULL;
    m_pfPrintResidBuffer = NULL;
    m_pfRefineErrorPercent = NULL;
    m_pfRefineErrorPercentDiffs = NULL;
    m_pfFFTArray = NULL;
    m_pnSort = NULL;

    m_poDiffList = NULL;

    (void) m_oList.nExpandGetField(m_oList.ms_sfRecipCoordD0);
    (void) m_oList.nExpandGetField(m_oList.ms_sfRecipCoordD1);
    (void) m_oList.nExpandGetField(m_oList.ms_sfRecipCoordD2);

    // This will not fill in twinID fields  and the X vectors.
    if (m_oList.nCalcRecipLatticePoints(m_oHeader,f0 = -1.0)) 
        return;

    vInit();

    m_bIsAvailable = TRUE;
};

void CEyeballIndex::vInit() {
    int nRef;
    int nx;
    double a3fX[3];
    double fDet;

    // Must be called in this order:
    // 1) Constructor calls nCalcRecipLatticePoints() to load X vectors.
    // 2) Calculate wavelength.
    // 3) Expand the nIndexSelect field.
    // 4) Expand other nOverlapCheck related fields.
    // 5) Reject reflections for I/sig or reso using nIndexSelect field.

    vCalcGetMaxCellAndWavelength();
    m_oList.nExpandGetField(m_oList.ms_snTwinID);
    m_oList.nExpandGetField(m_oList.ms_sfDetResolution);
    m_nFI_nIndexSelect = m_oList.nExpandGetField("nIndexSelect");
    
    m_nFI_nNeighborFirst = m_oList.nExpandGetField("nNeighborFirst");
    m_nFI_nNeighborNext = m_oList.nExpandGetField("nNeighborNext");
    m_oList.nExpandGetField("nIndex");

    m_nTwinID = -1;
    m_nVerbose = 1;
    m_nPrint = 0;

    m_f3DAngleStep = 5.0;
    m_n3DGridSpacing = 200;
    m_fMax3DAngleBetweenSolns = 10;
    m_f3DAngleStepRefineMax = 1.0;
    m_f3DAngleStepRefineMin = 0.01;
    m_n3DGridSpacingRefine = 400;

    m_f2DAngleStep = 2.0;           
    m_fMax2DAngleBetweenSolns = 10;
    m_n2DGridSpacing = 1000;        
    m_n2DGridSpacingRefine = 10000; 
    
    m_fMaxPruneReject = 0.50;     
    m_nMinRefsAfterPrune = 30;

    m_n3DGridSize = 0;
    m_nFFTArraySize = 0;

    m_nNumBest3DVectors = 0;
    m_nMaxBest3DVectors = m_nBest3DVectors;
    m_nNumBest2DVectors = 0;
    m_nMaxBest2DVectors = m_nBest2DVectors;

    for (nx=0;nx<3;nx++) {
        m_a3fMinX[nx] = 0.0;
        m_a3fMaxX[nx] = 0.0;
    };
    m_fMaxX = 0.0;

    for (nRef=0;nRef<m_oList.nGetNumReflns();nRef++) {
        // Check all reflections regardless of twin status.
        a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
        a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
        a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
        fDet = fLenVec3D(a3fX);
        if (0.0 < fDet) {
            fDet = m_fWavelength / fDet;
            m_oList[nRef].vSetField(m_oList.m_nFI_fDetResolution,(float) fDet);
        } else
            m_oList[nRef].vSetField(m_oList.m_nFI_fDetResolution,(float) -1.0);

        m_oList[nRef].vSetField(m_nFI_nIndexSelect,1);

        for (nx=0;nx<3;nx++) {
            m_a3fMinX[nx] = min(a3fX[nx],m_a3fMinX[nx]);
            m_fMaxX = max(fLenVec3D(a3fX),m_fMaxX);
            m_a3fMaxX[nx] = max(a3fX[nx],m_a3fMaxX[nx]);
        };

    };    

    nRejectIoverSigReso(-1.0,99999999.0,0.0);
    
    // RB 5/19/04  Need to reset the random number generator, so that we could have reproducible results when 
    // dtindex is run repeatedly.
    srand(1);
};

CEyeballIndex::~CEyeballIndex() {
    if (m_pn3DGridSums)
        delete[] m_pn3DGridSums;   
    m_pn3DGridSums = NULL;

    if (m_pnSort)
        delete[] m_pnSort;
    m_pnSort = NULL;

    if (m_pfPrintResidBuffer)
        delete[] m_pfPrintResidBuffer;    
    m_pfPrintResidBuffer = NULL;

    if (m_pfRefineErrorPercent)
        delete[] m_pfRefineErrorPercent;
    m_pfRefineErrorPercent = NULL;

    if (m_pfRefineErrorPercentDiffs)
        delete[] m_pfRefineErrorPercentDiffs;
    m_pfRefineErrorPercentDiffs = NULL;

    if (m_pfFFTArray)
        delete[] m_pfFFTArray;
    m_pfFFTArray = NULL;
    
    if (m_pnPrintNeighbors)
        delete[] m_pnPrintNeighbors;
    m_pnPrintNeighbors = NULL;

    if (m_poDiffList)
        delete m_poDiffList;
    m_poDiffList = NULL;
};

void CEyeballIndex::vPercentComplete(int nStatus,char* pcMessage) {
    static int nLastStatus;
#ifdef DOSEM_STATUSBAR
    if ((nStatus<0) && (!pcMessage))
        df.percent_bar(PB_HIDE);
    else if (pcMessage) {
        df.percent_bar(PB_SHOW,nStatus,100,pcMessage);
    } else 
        df.percent_bar(PB_UPDATE,nStatus);
#else
    if ((nStatus<0) && (!pcMessage))
        printf("\n");
    else if (pcMessage) {
        printf("%*s: %s",strlen(pcMessage),"","00---10---20---30---40---50---60---70---80---90---100\n");
        printf("%s: ",pcMessage);
        nLastStatus = 0;
    } else {
        while (nLastStatus<nStatus) {
            if (nLastStatus % 2)
                printf(".");
            nLastStatus++;
        };
    };
#endif
};



void
CEyeballIndex::vCalcGetMaxCellAndWavelength() {
    // Try to calculate a reasonable default value for the maximum cell
    // length that can/should be indexed based on the minimum spot
    // separation and the crystal to detector distance.
    
    Cstring sFirstDetectorName;
    int nStat;
    double fNum;
    const double fMMSeparationMin = 0.65;

    m_fWavelength = 1.0;
    m_fCellLengthMax = 1000.0;

    Csource oSource(m_oHeader);
    if (!oSource.bIsAvailable()) 
        return;
    m_fWavelength = oSource.fGetWavelength();

    nStat = m_oHeader.nGetValue(Cdetector::ms_sDetectorNames,1, &sFirstDetectorName);
    if (0 == nStat) {
        Cdetector oFirstDetector(m_oHeader,sFirstDetectorName,TRUE,FALSE);
        if (oFirstDetector.bIsAvailable()) {
            fNum = oFirstDetector.fGetDistance();
            if (0.0 != fNum)
            {
                m_fCellLengthMax = (fabs(fNum) * m_fWavelength / fMMSeparationMin);
                m_fCellLengthMax = 1000.0;
                return;
            }
        };
    };        
    
}



int CEyeballIndex::nRejectIoverSigReso(double fMinIoverSig,double fMinReso,double fMaxReso) {

    int nRef;
    double f0;

    if (fMinReso<fMaxReso) {
        f0 = fMinReso;
        fMinReso = fMaxReso;
        fMaxReso = f0;
    };
        
    m_fMinIoverSig = fMinIoverSig;
    m_fMinReso = fMinReso;
    m_fMaxReso = fMaxReso;

    for (nRef=0;nRef<m_nNumRefs;nRef++) {
        double fIntensity = m_oList[nRef].fGetIntensity();
        double fSigmaI = m_oList[nRef].fGetSigmaI();
        double fIoverSig = fIntensity/max(1.0,fSigmaI);

        m_oList[nRef].vSetField(m_nFI_nIndexSelect,0);
        if (fIoverSig>=fMinIoverSig) {
            if ((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
                f0 = m_oList[nRef].fGetField(m_oList.m_nFI_fDetResolution);
                if ((f0>=fMaxReso) && (f0<=fMinReso))
                    m_oList[nRef].vSetField(m_nFI_nIndexSelect,1);
            };
        };
    };
    return 0;
};

int CEyeballIndex::nGetBestPlane2D(int nSolution,double* pfPlaneMat,double* pfDSpacing,double* pfRealSpacing) {

    if (nSolution>=m_nNumBest2DVectors)
        return 1;

    if (pfPlaneMat)
        vCopyMat3D(&m_aa3fBest2DVecs[m_aanBest2DVectorSortOrder[nSolution]][0][0],pfPlaneMat);
    if (pfDSpacing)
        *pfDSpacing = m_aafBest2DVecsDSpacing[m_aanBest2DVectorSortOrder[nSolution]];
    if (pfRealSpacing)
        *pfRealSpacing = m_aafBest2DVecsRealSpacing[m_aanBest2DVectorSortOrder[nSolution]];
    return 0;

};

int CEyeballIndex::nGetBestViewDirections3D(int nBufSize,double* pfVecs) {
    int anDirectionsStored[m_nBest3DVectors];
    int nVectorsStored;
    int nx,ny;

    memset(anDirectionsStored,0,sizeof(int)*m_nBest3DVectors);
    nVectorsStored = 0;
    do {
        for (nx=0,ny=-1;nx<m_nNumBest3DVectors;nx++) {
            if (anDirectionsStored[nx]==0) {
                if ((ny==-1) || (m_aanBest3DVectorCounts[nx]<m_aanBest3DVectorCounts[ny]))
                    ny = nx;
            };
        };
        if (ny!=-1)  {
            vCopyVec3D(&m_aa3fBest3DVecs[ny][0],pfVecs+nVectorsStored*3);
            nVectorsStored++;
            anDirectionsStored[ny] = 1;
        };
    } while ((ny!=-1) && (nVectorsStored<nBufSize));
    return nVectorsStored;
};

int CEyeballIndex::nCalcHeuristic2D(double a3fProjVector[3],int nGridDim,bool bPrune) {
    double a3x3fProjBasis[3][3];
    int nx,ny;
    double f0;
    double a3fTemp[3];
    
    double a3fX[3];    
    double a3fProj[3];
    double a2fProj[2];
    int    a2nBin[2];
    int    nBin;
    
    int a20nHist[20];
    int nTotalContrib;
    int nIgnoreBin;
    int nUsedBins;
    int nRef;
    int nUsedRefs;
    int nTotalGrid;

    nTotalGrid = nGridDim*nGridDim;
    
    if (nTotalGrid>m_n3DGridSize) {
        m_n3DGridSize = nTotalGrid;
        if (m_pn3DGridSums)
            delete[] m_pn3DGridSums;
        m_pn3DGridSums = new int[nTotalGrid];
    };

    if (bPrune) {
        nRejectIoverSigReso(m_fMinIoverSig,m_fMinReso,m_fMaxReso);
    };
    

    // Establish two other vectors perpendicular to this one.
    vBuildBasis3D(a3fProjVector,a3x3fProjBasis);
    vNormMat3D(a3x3fProjBasis);
    
    memset(m_pn3DGridSums,0,sizeof(int)*nTotalGrid);
    nUsedRefs = 0;
    for (nRef=0;nRef<m_nNumRefs;nRef++) {
        // Only use reflections that have a twin id of 0 (indicating that they have not
        // yet been indexed.
        if ((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
            if (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==1) {
                a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
                a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
                a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
                f0 = fDot3D(a3fProjVector,a3fX);
                vMulVec3DScalar(a3fProjVector,f0,a3fTemp);
                vSubVec3DVec3D(a3fX,a3fTemp,a3fProj);
                for (nx=0;nx<2;nx++) {
                    a2fProj[nx] = fDot3D(a3x3fProjBasis[nx+1],a3fProj);
                    a2nBin[nx] = (int) (nGridDim/2 + a2fProj[nx]*nGridDim/2/m_fMaxX);
                    a2nBin[nx] = min(nGridDim-1,max(0,a2nBin[nx]));
                };
                nBin = a2nBin[1]*nGridDim + a2nBin[0];
                m_pn3DGridSums[nBin]++;
                nUsedRefs++;
            };
        };
    };
    
    memset(a20nHist,0,sizeof(int)*20);
    nTotalContrib = 0;
    for (nx=0;nx<nTotalGrid;nx++) {
        a20nHist[min(m_pn3DGridSums[nx],20-1)]+=m_pn3DGridSums[nx];
        nTotalContrib += m_pn3DGridSums[nx];
    };
    for (nIgnoreBin=1,ny=0;nIgnoreBin<20;nIgnoreBin++) {
        ny += a20nHist[nIgnoreBin];
        if (ny>=nTotalContrib*0.2)
            break;
    };
    
    nUsedBins = 0;
    for (nx=0;nx<nTotalGrid;nx++) {
        if (m_pn3DGridSums[nx]>=nIgnoreBin)
            nUsedBins++;

    };
    
    if (bPrune) {
        int nPrunned = 0;
        int nPruneSize = 1;
        do {
            for (nRef=0;nRef<m_nNumRefs;nRef++) {
                if ((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
                    if (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==1) {
                        a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
                        a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
                        a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
                        f0 = fDot3D(a3fProjVector,a3fX);
                        vMulVec3DScalar(a3fProjVector,f0,a3fTemp);
                        vSubVec3DVec3D(a3fX,a3fTemp,a3fProj);
                        for (nx=0;nx<2;nx++) {
                            a2fProj[nx] = fDot3D(a3x3fProjBasis[nx+1],a3fProj);
                            a2nBin[nx] = (int) (nGridDim/2 + a2fProj[nx]*nGridDim/2/m_fMaxX);
                            a2nBin[nx] = min(nGridDim-1,max(0,a2nBin[nx]));
                        };
                        nBin = a2nBin[1]*nGridDim + a2nBin[0];
                        if (m_pn3DGridSums[nBin]==nPruneSize) {
                            if ((nPrunned<nUsedRefs*m_fMaxPruneReject) && (nUsedRefs-nPrunned>=m_nMinRefsAfterPrune)) {
                                m_oList[nRef].vSetField(m_nFI_nIndexSelect,0);
                                nPrunned++;
                            };
                        };
                    };
                };
            };
            nPruneSize++;
        } while (nPruneSize<20);
    };
    return nUsedBins;
};

int CEyeballIndex::nCalcHeuristic1D(double a3fCompressVector[3],int nGridDim,bool bPrune) {
    int nx,ny;
    double f0;
    
    double a3fX[3];    
    int    nBin;
    
    int a20nHist[20];
    int nTotalContrib;
    int nIgnoreBin;
    int nUsedBins;
    int nRef;
    int nTotalGrid;
    int nUsedRefs;
    
    nTotalGrid = nGridDim;



    // Use the 3D grid for our calculations.
    if (nTotalGrid>m_n3DGridSize) {
        m_n3DGridSize = nTotalGrid;
        if (m_pn3DGridSums)
            delete[] m_pn3DGridSums;
        m_pn3DGridSums = new int[nTotalGrid];
    };
    fNormVec3D(a3fCompressVector);
       
    memset(m_pn3DGridSums,0,sizeof(int)*nTotalGrid);
    nUsedRefs = 0;
    for (nRef=0;nRef<m_nNumRefs;nRef++) {
        // Only use reflections that have a twin id of 0 (indicating that they have not
        // yet been indexed.
        if ((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
            if (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==1) {
                a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
                a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
                a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
                // Project the vector onton the compression direction.
                f0 = fDot3D(a3fCompressVector,a3fX);
                nBin = (int) (nGridDim/2 + f0*nGridDim/2/m_fMaxX);
                m_pn3DGridSums[nBin]++;
                nUsedRefs++;
            };
        };
    };
    
    memset(a20nHist,0,sizeof(int)*20);
    nTotalContrib = 0;
    for (nx=0;nx<nTotalGrid;nx++) {
        a20nHist[min(m_pn3DGridSums[nx],20-1)]+=m_pn3DGridSums[nx];
        nTotalContrib += m_pn3DGridSums[nx];
    };
    for (nIgnoreBin=1,ny=0;nIgnoreBin<20;nIgnoreBin++) {
        ny += a20nHist[nIgnoreBin];
        if (ny>=nTotalContrib*0.2)
            break;
    };
    
    if (m_nPrint) {
        FILE* pFOut;
        static int nFile = 0;
        Cstring sFile;
        sFile = "1d";
        sFile += nFile++;
        sFile += ".txt";
        pFOut = fopen(sFile.string(),"wt");
        if (pFOut) {
            for (nx=0;nx<nTotalGrid;nx++) {
                fprintf(pFOut,"%d %d\n",nx,m_pn3DGridSums[nx]);
            };
            fclose(pFOut);
        };
    };
    
    if (bPrune) {
        int nPrunned;
        int nPruneSize = 1;

        nPrunned = 0;          
        do {
            for (nRef=0;nRef<m_nNumRefs;nRef++) {
                if ((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
                    if (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==1) {
                        a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
                        a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
                        a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
                        // Project the vector onton the compression direction.
                        f0 = fDot3D(a3fCompressVector,a3fX);
                        nBin = (int) (nGridDim/2 + f0*nGridDim/2/m_fMaxX);

                        if (m_pn3DGridSums[nBin]==nPruneSize) {
                            if ((nPrunned<nUsedRefs*m_fMaxPruneReject) && (nUsedRefs-nPrunned>=m_nMinRefsAfterPrune)) {
                                m_oList[nRef].vSetField(m_nFI_nIndexSelect,0);
                                nPrunned++;
                            };
                        };
                    };
                };
            };
            nPruneSize++;
        } while (nPruneSize<40);
    };
    
    nUsedBins = 0;
    for (nx=0;nx<nTotalGrid;nx++) {
        if (m_pn3DGridSums[nx]>=nIgnoreBin)
            nUsedBins++;
    };
    return nUsedBins;
};


int CEyeballIndex::nCalcRefineBestViewDirections3D() {
    double fPhi,fTheta;
    double fThetaSteps;
    double a3fProjVector[3];

    int nTotalTests;
    int nTotalGrid;
    int nGridDim;
    int nUsedBins;
    int nWorstSoln;
    int nx;
    int nVector;
    
    nGridDim = m_n3DGridSpacing;
    nTotalGrid = nGridDim*nGridDim;
    m_nNumBest3DVectors = 0;
    nTotalTests = 0;

    vPercentComplete(0,"Calculating 3D vectors");

    for (fPhi=0.0;fPhi<=90;fPhi+=m_f3DAngleStep) {
        fThetaSteps = cos(fPhi*Gs_dRADIANS_PER_DEGREE)*2.0*Gs_dPI/(m_f3DAngleStep*Gs_dRADIANS_PER_DEGREE);
        vPercentComplete((int) (100.0*(fPhi)/90.0));
        for (fTheta=0.0;fTheta<360.0;fTheta+=360.0/fThetaSteps) {
            // Get the projection vector.
            a3fProjVector[0] = cos(fTheta*Gs_dRADIANS_PER_DEGREE)*cos(fPhi*Gs_dRADIANS_PER_DEGREE);
            a3fProjVector[1] = sin(fTheta*Gs_dRADIANS_PER_DEGREE)*cos(fPhi*Gs_dRADIANS_PER_DEGREE);
            a3fProjVector[2] = sin(fPhi*Gs_dRADIANS_PER_DEGREE);

            nUsedBins = nCalcHeuristic2D(a3fProjVector,nGridDim);

            // We now know how many bins were used.

            nWorstSoln = 0;
            for (nx=0;nx<m_nNumBest3DVectors;nx++) {
                if ((fabs(m_aa2fBest3DVecsAngles[nx][0]-fPhi)<m_fMax3DAngleBetweenSolns) &&
                    (fabs(m_aa2fBest3DVecsAngles[nx][1]-fTheta)<m_fMax3DAngleBetweenSolns)) {
                    if (m_aanBest3DVectorCounts[nx]>nUsedBins) {

                        // The new solution is better.
                        m_aanBest3DVectorCounts[nx] = nUsedBins;
                        vCopyVec3D(a3fProjVector,m_aa3fBest3DVecs[nx]);
                        m_aa2fBest3DVecsAngles[nx][0] = fPhi;
                        m_aa2fBest3DVecsAngles[nx][1] = fTheta;

                    };
                    // In any case, break since this bin 'covers' the solution.
                    break;
                };
                if (m_aanBest3DVectorCounts[nx]>m_aanBest3DVectorCounts[nWorstSoln]) 
                    nWorstSoln = nx;
            };
            // If we made it all the way through, and either
            // 1) We have not filled all the bins yet, or 
            // 2) The worst solution is worse than this solution.
            if (nx>=m_nNumBest3DVectors) {

                if (m_nNumBest3DVectors<m_nMaxBest3DVectors) {
                    nWorstSoln = m_nNumBest3DVectors;
                    m_aanBest3DVectorCounts[nWorstSoln] = nUsedBins+1;
                    m_nNumBest3DVectors++;
                };                
                if (m_aanBest3DVectorCounts[nWorstSoln]>nUsedBins) {
                    
                    m_aanBest3DVectorCounts[nWorstSoln] = nUsedBins;
                    vCopyVec3D(a3fProjVector,m_aa3fBest3DVecs[nWorstSoln]);
                    m_aa2fBest3DVecsAngles[nWorstSoln][0] = fPhi;
                    m_aa2fBest3DVecsAngles[nWorstSoln][1] = fTheta;               
                };
            };
            nTotalTests++;
        };
    };
    vPercentComplete();

    vPercentComplete(0,"Refining 3D vectors");
    for (nVector=0;nVector<m_nNumBest3DVectors;nVector++) {
        m_aanBest3DVectorCounts[nVector] = nRefineBestViewDirections(&m_aa3fBest3DVecs[nVector][0],m_nRefine3);
        vPercentComplete(((nVector+1)*100)/m_nNumBest3DVectors);
    };

    int nVector2;
    int nNumBest;
    double f0;

    // Check for redundant vectors.
    nNumBest = 1;
    for (nVector=1;nVector<m_nNumBest3DVectors;nVector++) {
        for (nVector2=0;nVector2<nNumBest;nVector2++) {
            f0 = fMinVectorDifference(m_aa3fBest3DVecs[nVector2],m_aa3fBest3DVecs[nVector]);
            if (f0<0.01) 
                break;
        };
        if (nVector2==nNumBest) {
            vCopyVec3D(m_aa3fBest3DVecs[nVector],m_aa3fBest3DVecs[nNumBest++]);
        };

    };
    m_nNumBest3DVectors = nNumBest;
    vPercentComplete();

    return 0;
};



int CEyeballIndex::nCalcRefineBestViewDirections2D(double* pfPerpVector,bool bKeepPrune) {

    double a3fCompressVector[3];
    double a3fNormalVector[3];
    double a3x3fBasis[3][3];
    double fTheta;
    int nPass;

    double a3fTemp[3];
    int nx;

    double          afBest2DVecsAngles[m_nBest2DVectors];


    
    int             nUsedBins;
    int             nWorstSoln;
    int             nBestSoln;
    int             nVector;
    
    
    vBuildBasis3D(pfPerpVector,a3x3fBasis);
    vNormMat3D(a3x3fBasis);
    
    nCalcHeuristic2D(pfPerpVector,m_n3DGridSpacingRefine,TRUE); 
    
    vPercentComplete(0,"Calculating 2D planes");

    // On the first pass, we will find the best planes.
    // On the second pass, using the best found plane, we prune spots, and re-find the planes.
    for (nPass=0;nPass<2;nPass++) {
        
        m_nNumBest2DVectors = 0;
        for (fTheta=0.0;fTheta<180.0;fTheta+=m_f2DAngleStep) {
            vPercentComplete((int) (100*(fTheta/180.0)/2.0 + nPass*50));
            
            // Build the compression vector.
            vMulVec3DScalar(a3x3fBasis[1],cos(fTheta*Gs_dRADIANS_PER_DEGREE),a3fCompressVector);
            vMulVec3DScalar(a3x3fBasis[2],sin(fTheta*Gs_dRADIANS_PER_DEGREE),a3fTemp);
            vAddVec3DVec3D(a3fTemp,a3fCompressVector,a3fCompressVector);
            vMulVec3DScalar(a3x3fBasis[1],-sin(fTheta*Gs_dRADIANS_PER_DEGREE),a3fNormalVector);
            vMulVec3DScalar(a3x3fBasis[2],cos(fTheta*Gs_dRADIANS_PER_DEGREE),a3fTemp);
            vAddVec3DVec3D(a3fTemp,a3fNormalVector,a3fNormalVector);
            
            // The 'compression' vector is passed to the heuristic. However, we are saving in terms.
            // of a 'normal' vector.
            nUsedBins = nCalcHeuristic1D(a3fCompressVector,m_n2DGridSpacing); 
            
            // We now know how many bins were used.
            
            nWorstSoln = 0;
            nBestSoln = 0;
            for (nx=0;nx<m_nNumBest2DVectors;nx++) {
                if ((fabs(afBest2DVecsAngles[nx]-fTheta)<m_fMax2DAngleBetweenSolns) ||
                    (fabs(afBest2DVecsAngles[nx]-fTheta)>180.0 - m_fMax2DAngleBetweenSolns)) {
                    
                    if (m_aanBest2DVectorCounts[nx]>nUsedBins) {                    
                        // The new solution is better.
                        m_aanBest2DVectorCounts[nx] = nUsedBins;
                        vCopyVec3D(a3fCompressVector,m_aa3fBest2DVecs[nx][0]);
                        afBest2DVecsAngles[nx] = fTheta;                   
                    };
                    // In any case, break since this bin 'covers' the solution.
                    break;
                };
                if (m_aanBest2DVectorCounts[nx]>m_aanBest2DVectorCounts[nWorstSoln]) 
                    nWorstSoln = nx;
                if (m_aanBest2DVectorCounts[nx]<m_aanBest2DVectorCounts[nBestSoln])
                    nBestSoln = nx;
            };
            // If we made it all the way through, and either
            // 1) We have not filled all the bins yet, or 
            // 2) The worst solution is worse than this solution.
            if (nx>=m_nNumBest2DVectors) {
                
                if (m_nNumBest2DVectors<m_nMaxBest2DVectors) {
                    nWorstSoln = m_nNumBest2DVectors;
                    m_aanBest2DVectorCounts[nWorstSoln] = nUsedBins+1;
                    m_nNumBest2DVectors++;
                };                
                if (m_aanBest2DVectorCounts[nWorstSoln]>nUsedBins) {
                    
                    m_aanBest2DVectorCounts[nWorstSoln] = nUsedBins;
                    vCopyVec3D(a3fCompressVector,m_aa3fBest2DVecs[nWorstSoln][0]);
                    afBest2DVecsAngles[nWorstSoln]= fTheta;
                };
            };
        };

        // Prune spots that don't fit well to the best solution.
        if (nPass == 0) 
            nCalcHeuristic1D(m_aa3fBest2DVecs[nBestSoln][0],m_n2DGridSpacing,TRUE); 
    };
    vPercentComplete();
    
    vPercentComplete(0,"Refining 2D vectors");
   
    for (nVector=0;nVector<m_nNumBest2DVectors;nVector++) {

        m_aanBest2DVectorCounts[nVector] = nRefineBestViewDirections(&m_aa3fBest2DVecs[nVector][0][0],m_nRefine2,pfPerpVector);
        // The other two plane vectors are used for display, but should be based on the "pfPerpVector" passed to this routine.
        
        vCopyVec3D(pfPerpVector,&m_aa3fBest2DVecs[nVector][1][0]);
        vCross3D(&m_aa3fBest2DVecs[nVector][0][0],&m_aa3fBest2DVecs[nVector][1][0],&m_aa3fBest2DVecs[nVector][2][0]);
        // Get the D spacing from FFT.  
        // [0] contains the 'Normal' vector, (The vector down which the FFT is taken)
        // [1] contains the 'Perpendicular' vector
        // [2] contains the Contains the cross product of the first two.

        m_aafBest2DVecsRealSpacing[nVector] = fDPSFFT(&m_aa3fBest2DVecs[nVector][0][0]);
        m_aafBest2DVecsDSpacing[nVector] = m_fWavelength/m_aafBest2DVecsRealSpacing[nVector];
        vPercentComplete(((nVector+1)*100)/m_nNumBest2DVectors);
    };
    g_pnCmpInts = &m_aanBest2DVectorCounts[0];
    for (nx=0;nx<m_nNumBest2DVectors;nx++)
        m_aanBest2DVectorSortOrder[nx] = nx;
    qsort(&m_aanBest2DVectorSortOrder[0],m_nNumBest2DVectors,sizeof(int),int_cmp_rel);
    vPercentComplete();    

    // Reset the rejection flags (unless the user wants us to keep the prunning settings.
    if (!bKeepPrune)
        nRejectIoverSigReso(m_fMinIoverSig,m_fMinReso,m_fMaxReso);

    return 0;
};



int CEyeballIndex::nRefineBestViewDirections(double* pfVector,int nMethod,double* pfRestriction)
{
    int nTotalGrid;
    int nGridDim;
    double a3fTemp[3];
    double f0;

    double fStep;
    double a3fVec[3];
    int    nValue;
    double aa3fDeltaVecs[8][3];
    double afDeltaVecValues[8];
    double a3x3fRotVecs[3][3];
    double a3x3fRotMat[3][3];

    int nRot;
    int nMinRot;
    double a8fRot1[8] = {-1.0, 0.0, 1.0,-1.0,1.0,-1.0,0.0,1.0};
    double a8fRot2[8] = {-1.0,-1.0,-1.0, 0.0,0.0, 1.0,1.0,1.0};


    if (nMethod == m_nRefine3)
        nGridDim = m_n3DGridSpacingRefine;
    else if (nMethod == m_nRefine2)
        nGridDim = m_n2DGridSpacingRefine;
    else 
        return -1;

    nTotalGrid = nGridDim*nGridDim;

   
    fStep = m_f3DAngleStepRefineMax;
    vCopyVec3D(pfVector,a3fVec);
    if (nMethod == m_nRefine3)
        nValue = nCalcHeuristic2D(a3fVec,nGridDim);
    else 
        nValue = nCalcHeuristic1D(a3fVec,nGridDim);
    
    do {
        // For the current vector, get 2 rotation vectors.
        vBuildBasis3D(a3fVec,a3x3fRotVecs);
        vNormMat3D(a3x3fRotVecs);
        for (nRot=0;nRot<8;nRot++) {
            
            vConv3Rot3Vec3DMat3D(
                0.0,fStep*a8fRot1[nRot],fStep*a8fRot2[nRot],
                &a3x3fRotVecs[0][0],
                &a3x3fRotVecs[1][0],
                &a3x3fRotVecs[2][0],
                &a3x3fRotMat[0][0]);
            vMulMat3DVec3D(a3x3fRotMat,a3fVec,&aa3fDeltaVecs[nRot][0]);

            if (pfRestriction) {
                // Keep the solution normal to the restriction.
                f0 = fDot3D(&aa3fDeltaVecs[nRot][0],pfRestriction);
                vMulVec3DScalar(pfRestriction,f0,a3fTemp);
                vSubVec3DVec3D(&aa3fDeltaVecs[nRot][0],a3fTemp,&aa3fDeltaVecs[nRot][0]);
            };
            if (nMethod == m_nRefine3)
                afDeltaVecValues[nRot] = nCalcHeuristic2D(&aa3fDeltaVecs[nRot][0],nGridDim);
            else 
                afDeltaVecValues[nRot] = nCalcHeuristic1D(&aa3fDeltaVecs[nRot][0],nGridDim);

            
        };
        // See if the value has decreased for any of our tests.
        for (nRot=0,nMinRot=0;nRot<8;nRot++) {
            if (afDeltaVecValues[nRot]<afDeltaVecValues[nMinRot])
                nMinRot = nRot;
        };
        if (afDeltaVecValues[nMinRot]<nValue) {
            vCopyVec3D(&aa3fDeltaVecs[nMinRot][0],a3fVec);
            nValue = (int) (afDeltaVecValues[nMinRot]);
        } else {
            fStep/=2.0;
        };
        
        
    } while (fStep>m_f3DAngleStepRefineMin);
    
    vCopyVec3D(&a3fVec[0],pfVector);
    return nValue;
};


bool CEyeballIndex::bFindLC(double* pfErrPercent,double* pfBasisVecs,int nNumVecs,double* pfVec,double* pfVecD,int* pnCoeffs,int nNormalize) {
    double a3fVec[3];        // Place the LC vector here.
    double a3fDiff[3];       // a3fVec - pfVec
    double a3fHKLDiff[3];    // fabs(Int(HKL) - float(HKL))
    double fError;
    double fDet;

    double a3x3fMat[3][3];
    double a3x3fInvMat[3][3];
    double a3fLC[3];
    double a3fPrevLC[3];
    double f0;
    int nx;

    if (nNumVecs==1) {
        f0 = fLenVec3D(pfVec)/fLenVec3D(pfBasisVecs);
        if (f0-floor(f0)<0.5)
            pnCoeffs[0] = (int) floor(f0);
        else
            pnCoeffs[0] = (int) (1+floor(f0));
        a3fHKLDiff[0] = a3fHKLDiff[1] = a3fHKLDiff[2] = fabs((double) pnCoeffs[0] - f0);

        if (fDot3D(pfVec,pfBasisVecs)<0.0)
            pnCoeffs[0]*= -1;
        vMulVec3DScalar(pfBasisVecs,(double) pnCoeffs[0],a3fVec);      
        fDet = 10000.0;

        vSubVec3DVec3D(a3fVec,pfVec,a3fDiff);
    } else if (nNumVecs==2) {
        for (nx=0;nx<2;nx++)
            vCopyVec3D(pfBasisVecs+3*nx,a3x3fMat[nx]);
        vCross3D(a3x3fMat[0],a3x3fMat[1],a3x3fMat[2]);
        if (0.0 == (fDet = fInvMat3D(& a3x3fMat[0][0],& a3x3fInvMat[0][0])))
            return FALSE;
        vMulMat3DVec3D(a3x3fInvMat,pfVec,a3fLC);
        vCopyVec3D(a3fLC,a3fPrevLC);
        for (nx=0;nx<3;nx++) {
            if (a3fLC[nx]-floor(a3fLC[nx])<0.5)
                a3fLC[nx] = pnCoeffs[nx] = (int) floor(a3fLC[nx]);
            else
                a3fLC[nx] = pnCoeffs[nx] = (int) (1+floor(a3fLC[nx]));
            a3fHKLDiff[nx] = fabs(pnCoeffs[nx] - a3fPrevLC[nx]);
        };
        if (pnCoeffs[2]!=0)
            return FALSE;
        vMulMat3DVec3D(a3x3fMat,a3fLC,a3fVec);
        vSubVec3DVec3D(a3fVec,pfVec,a3fDiff);
    } else if (nNumVecs==3) {
        for (nx=0;nx<3;nx++)
            vCopyVec3D(pfBasisVecs+3*nx,a3x3fMat[nx]);

        if (0.0 == (fDet= fInvMat3D(& a3x3fMat[0][0],& a3x3fInvMat[0][0])))
            return FALSE;
        vMulMat3DVec3D(a3x3fInvMat,pfVec,a3fLC);
        vCopyVec3D(a3fLC,a3fPrevLC);
        for (nx=0;nx<3;nx++) {
            if (a3fLC[nx]-floor(a3fLC[nx])<0.5)
                a3fLC[nx] = pnCoeffs[nx] = (int) floor(a3fLC[nx]);
            else
                a3fLC[nx] = pnCoeffs[nx] = (int) (1+floor(a3fLC[nx]));
            a3fHKLDiff[nx] = fabs(pnCoeffs[nx] - a3fPrevLC[nx]);
           
        };
        vMulMat3DVec3D(a3x3fMat,a3fLC,a3fVec);
        vSubVec3DVec3D(a3fVec,pfVec,a3fDiff);


        // If the user provided a derivative vector (assumed normalized to 1.0)
        if (pfVecD) {
            
            double a3fTemp[3];
            double a3fTemp2[3];

            // Find the projection of the difference amount onto this derivative vector
            // and subtract the projection.
            vMulVec3DScalar(pfVecD,fDot3D(pfVecD,a3fDiff),a3fTemp);
            vSubVec3DVec3D(a3fDiff,a3fTemp,a3fDiff);

            vMulMat3DVec3D(a3x3fInvMat,a3fTemp,a3fTemp2);
            for (nx=0;nx<3;nx++) {
                if (fabs(a3fTemp2[nx])>=0.5) {
                    *pfErrPercent = 1.0;
                    return TRUE;
                };
            };
        };
    };
    



    if (nNormalize)
        fError = 0.5*fabs(fLenVec3D(a3fDiff))/(fLenVec3D(pfVec) + fLenVec3D(a3fVec));
    else
        fError = fLenVec3D(a3fHKLDiff);    

    if (fError<0.0)
        fError = 1.0;
    *pfErrPercent = fError;
    return TRUE;
};


int CEyeballIndex::nCalcPrintableEdges(double* pfIndexEdges,
                                       int nNumIndexEdges,
                                       double fMaxResidual,
                                       int nMaxEdgesToGenerate,
                                       double* pfEdgesOut0,
                                       double* pfEdgesOut1,
                                       int nTwinID) {
    int nx;
    double f0;
    
    int nFI_nEdge[3];
    int nFI_fResid;
    int nFI_nPackedHKL;
    int a3nHKL[3];
    int a3nHKL2[3];
    int a3nHKLD[3];
    int a3x2nDelta[3][2];
    int a3nDelta[3];
    
    int nRef,nRef2;
    int nPackedHKL;
    int nNeighborsFound;
    
    double a3fX[3];
    double a3fX2[3];
    double a3fXDiff[3];
    double fErrorPercent;
    bool bIndexed;
    
    nFI_nEdge[0] = m_oList.nExpandGetField("nIndexH");
    nFI_nEdge[1] = m_oList.nExpandGetField("nIndexK");
    nFI_nEdge[2] = m_oList.nExpandGetField("nIndexL");
    nFI_fResid = m_oList.nExpandGetField("fIndexResidual");
    nFI_nPackedHKL = m_oList.nExpandGetField(Creflnlist::ms_snPackedHKL);
    
    // Try to index all reciprocal lattice points.
    // place the 'HKL' values in the reflection list, but not in the 'HKL' entry place.
    
    for (nRef=0;nRef<m_nNumRefs;nRef++) {
        if ((nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==nTwinID)) {
            a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
            a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
            a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
            bIndexed = bFindLC(&fErrorPercent,pfIndexEdges,nNumIndexEdges,a3fX,NULL,a3nHKL);
            if (!bIndexed)
                fErrorPercent = 100.0;
            m_oList[nRef].vSetField(nFI_fResid,(float) fErrorPercent);
            for (nx=nNumIndexEdges;nx<3;nx++)
                a3nHKL[nx] = 0;
            for (nx=0;nx<3;nx++)
                m_oList[nRef].vSetField(nFI_nEdge[nx],a3nHKL[nx]);
            // Set the packed HKL value.
            m_oList[nRef].vSetField(nFI_nPackedHKL,m_oList[nRef].nPackHKL(a3nHKL));
        };
    };
    
    // Sort the reflection list on packed HKL.
    if (!m_pnSort)
        m_pnSort = new int[m_nNumRefs];
    if (!m_pfPrintResidBuffer) {
        m_pfPrintResidBuffer = new double[m_nNumRefs*8];
    };
    if (!m_pnPrintNeighbors) {
        m_pnPrintNeighbors = new int[m_nNumRefs*8];
    };
    
    
    m_oList.vSort(eReflnField_int_type,nFI_nPackedHKL,m_pnSort);
    
    for (nx=0;nx<3;nx++) {
        if (nx<nNumIndexEdges) {
            a3x2nDelta[nx][0] = -1;
            a3x2nDelta[nx][1] = 1;
        } else {
            a3x2nDelta[nx][0] = 0;
            a3x2nDelta[nx][1] = 0;
        };
    };
    
    for (nRef=0;nRef<m_nNumRefs;nRef++) {
        if ((nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==nTwinID)) {
            if (m_oList[nRef].fGetField(nFI_fResid)<fMaxResidual) {
                a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
                a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
                a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
                
                nPackedHKL = m_oList[nRef].nGetField(nFI_nPackedHKL);
                for (nx=0;nx<3;nx++)
                    a3nHKL[nx] = m_oList[nRef].nGetField(nFI_nEdge[nx]);
                for (nx=0;nx<8;nx++)
                    m_pfPrintResidBuffer[nRef*8+nx] = -1.0;
                nNeighborsFound = 0;
                for (a3nDelta[0] = a3x2nDelta[0][0]; a3nDelta[0]<=a3x2nDelta[0][1]; a3nDelta[0]++) {
                    for (a3nDelta[1] = a3x2nDelta[1][0]; a3nDelta[1]<=a3x2nDelta[1][1]; a3nDelta[1]++) {
                        for (a3nDelta[2] = a3x2nDelta[2][0]; a3nDelta[2]<=a3x2nDelta[2][1]; a3nDelta[2]++) {
                            if ((a3nDelta[0]!=0) || (a3nDelta[1]!=0) || (a3nDelta[2]!=0)) {
                                for (nx=0;nx<3;nx++)
                                    a3nHKL2[nx] = a3nDelta[nx] + a3nHKL[nx];
                                // We have a tweaked HKL.
                                // See if we can find it's packed value.
                                nPackedHKL = m_oList[nRef].nPackHKL(a3nHKL2);
                                nRef2 = m_oList.nFindFirst(nFI_nPackedHKL,nPackedHKL,m_pnSort);
                                if (nRef2!=-1) {
                                    nRef2 = m_pnSort[nRef2];
                                    // We found an adjacent spot.  See what the residual for the delta vector is.
                                    if (nRef2>nRef) {
                                        // The spot is further along in the list (don't want to grab same edges twice!)
                                        a3fX2[0] = m_oList[nRef2].fGetField(m_oList.m_nFI_fRecipCoord0);
                                        a3fX2[1] = m_oList[nRef2].fGetField(m_oList.m_nFI_fRecipCoord1);
                                        a3fX2[2] = m_oList[nRef2].fGetField(m_oList.m_nFI_fRecipCoord2);
                                        vSubVec3DVec3D(a3fX,a3fX2,a3fXDiff);
                                        // Try to index the difference vector.  Place the residual in the buffer.
                                        bIndexed = bFindLC(&fErrorPercent,pfIndexEdges,nNumIndexEdges,a3fXDiff,NULL,a3nHKLD);
                                        if ((bIndexed) && (fErrorPercent<=fMaxResidual) && (abs(a3nHKLD[0]) + abs(a3nHKLD[1]) + abs(a3nHKLD[2]) == 1)) {
                                            m_pnPrintNeighbors[nRef*8+nNeighborsFound] = nRef2*4 + abs(a3nHKLD[0]) + 2*abs(a3nHKLD[1]) + 3*abs(a3nHKLD[2]);
                                            m_pfPrintResidBuffer[nRef*8+nNeighborsFound] = fErrorPercent;
                                            nNeighborsFound++;                                            
                                        };
                                    };                    
                                };
                            };
                        };
                    };
                };
            };
        };
    };

    // Now we must select a cut-off to get the required # of edges.

	double a2fCutOff[2] = {0, 0};
    double fCutOff;
    int nPass;
    int nPassingRefs;

    for (nPass=0;(nPass<2) || (a2fCutOff[1]-a2fCutOff[0]>0.001);nPass++) {
        if (nPass==0)
            fCutOff = 0.0;
        else if (nPass==1)
            fCutOff = 1.0;
        else 
            fCutOff = (a2fCutOff[0] + a2fCutOff[1])/2;

        nPassingRefs = 0;
        for (nRef=0;nRef<m_nNumRefs;nRef++) {
            for (nx=0;nx<8;nx++) {
                f0 = m_pfPrintResidBuffer[nRef*8+nx];
                if ((f0>=0.0) && (f0<fCutOff))
                    nPassingRefs++;
            };
        };
        if (nPass==0)
            a2fCutOff[0] = fCutOff;
        else if (nPass==1)
            a2fCutOff[1] = fCutOff;
        else if (nPassingRefs<nMaxEdgesToGenerate)
            a2fCutOff[0] = fCutOff;
        else 
            a2fCutOff[1] = fCutOff;
    };

    // Finally, we generate the edges themselves.
    nPassingRefs = 0;
    for (nRef=0;nRef<m_nNumRefs;nRef++) {
        if ((nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==nTwinID)) {
            for (nx=0;nx<8;nx++) {
                f0 = m_pfPrintResidBuffer[nRef*8+nx];
                if ((f0>=0.0) && (f0<fCutOff) && (nPassingRefs<nMaxEdgesToGenerate)) {
                    
                    a3fX[0] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord0);
                    a3fX[1] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord1);
                    a3fX[2] = m_oList[nRef].fGetField(m_oList.m_nFI_fRecipCoord2);
                    a3fX[0] += m_oList[m_pnPrintNeighbors[nRef*8+nx]/4].fGetField(m_oList.m_nFI_fRecipCoord0);
                    a3fX[1] += m_oList[m_pnPrintNeighbors[nRef*8+nx]/4].fGetField(m_oList.m_nFI_fRecipCoord1);
                    a3fX[2] += m_oList[m_pnPrintNeighbors[nRef*8+nx]/4].fGetField(m_oList.m_nFI_fRecipCoord2);
                    vMulVec3DScalar(a3fX,0.5,a3fX);
                    
                    vMulVec3DScalar(pfIndexEdges + ((m_pnPrintNeighbors[nRef*8+nx] % 4)-1)*3,0.5,a3fX2);
                    vSubVec3DVec3D(a3fX,a3fX2,pfEdgesOut0 + nPassingRefs*3);
                    vAddVec3DVec3D(a3fX,a3fX2,pfEdgesOut1 + nPassingRefs*3);
                    
                    nPassingRefs++;
                };
            };
        };
    };


    return nPassingRefs;
};

double CEyeballIndex::fMinVectorDifference(double a3fVec1[3],double a3fVec2[3]) {
    double a3fVec3[3];
    double f0;
    vAddVec3DVec3D(a3fVec1,a3fVec2,a3fVec3);
    f0 = fLenVec3D(a3fVec3);
    vSubVec3DVec3D(a3fVec1,a3fVec2,a3fVec3);
    f0 = min(fLenVec3D(a3fVec3),f0);
    return f0;
};

double CEyeballIndex::fDPSFFT(double a3fVecT[3])
{
  int           ni = 0;
  int           nj = 0;

  double        fTemp = 0.0;
  double        *pfTemp = NULL;
  double        *pfData = NULL;

  int           nNumBins;
  int           nMinBins      =   64; // Must be a power of 2
  int           nMaxBins      = 2048; // Must be a power of 2

  double                fP = 0.0;
  double                fPmin = 0.0;
  double                fPmax = 0.0;

  double                a3fX[3];
  
  double                fDeltaBin = 0.0;
  
  int                   nBin = 0;
  
  unsigned short int    uiFreq = 0U;
  
  int                   nNeeded = 0;
  
  int                   nFI_fFFTValue = 0;


  nFI_fFFTValue = m_oList.nExpandGetField("fFFTValue");

  // Loop through the input reflnlist list and form the dot product
  // p = d* . t ,  then put p in a histogram or frequency array
  //     --   -
  
  // Initialize fPmin and fPmax with first dot product
  for (ni = 0, nj = 0; ni < m_nNumRefs; ni++) {
      if ((m_nTwinID==-1) || (m_oList[ni].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
          if (m_oList[ni].nGetField(m_nFI_nIndexSelect)==1) {
              
              a3fX[0] = m_oList[ni].fGetField(m_oList.m_nFI_fRecipCoord0);
              a3fX[1] = m_oList[ni].fGetField(m_oList.m_nFI_fRecipCoord1);
              a3fX[2] = m_oList[ni].fGetField(m_oList.m_nFI_fRecipCoord2);
              fP         = fDot3D(a3fVecT, a3fX);
              
              // Store P in the fFloatH field for now     
              
              m_oList[ni].vSetField(nFI_fFFTValue, (float) fP);
              if ((nj==0) || (fPmin > fP)) 
                  fPmin = fP;
              if ((nj==0) || (fPmax < fP))
                  fPmax = fP;
              nj++;
          };
      };
  }
  
  // Now build frequency distribution (histogram)

  nNumBins  = nint((fPmax - fPmin) * 5.0 * m_fCellLengthMax);

  // Make nNumBins a power of 2        
  nNumBins = nint(log((double)nNumBins) / log(2.0));
  nNumBins = (int)pow(2.0, (double)nNumBins);

  // cout << "nNumBins " << nNumBins << endl;

  if (nMinBins > nNumBins)
    {
      nNumBins = nMinBins;
      if (6 <= m_nVerbose)
        cout << "WARNING, minimum frequency bins is " << nMinBins << endl;
    }
  if (nMaxBins < nNumBins)
    {
      nNumBins = nMaxBins;
      if (6 <= m_nVerbose)
        cout << "WARNING, maximum frequency bins is " << nMaxBins << endl;
    }

  fDeltaBin = (fPmax - fPmin) / (double)(nNumBins-1);

  if (0.0 >= fDeltaBin)
    {
      cout << "WARNING, minimum delta bin is 0.1" << endl;
      fDeltaBin = 0.1;
    }
//  cout << "deltabin, nNumBins: " << fDeltaBin << ", " << nNumBins << endl;

  int nBinTestMin, nBinTestMax;
  
  //  Minimum test bin of 1/10 the cell max?  I don't like it!
  //  nBinTestMin = nint(0.1f * m_fCellLengthMax * (float)(nNumBins-1) * fDeltaBin);

  nBinTestMin = 3;
  nBinTestMin = max(1, nBinTestMin);
  nBinTestMax = nint(m_fCellLengthMax * (double)(nNumBins-1) * fDeltaBin);
  nBinTestMax = min(nNumBins, nBinTestMax);

  // Make sure enough space exists for the histogram and FFT arrays
        
  nNeeded = nNumBins * 2 + 100;

  if (nNeeded > m_nFFTArraySize)
    {
      if (NULL != m_pfFFTArray)
        delete [] m_pfFFTArray;
      m_nFFTArraySize = nNeeded;
      m_pfFFTArray = new double [nNeeded];
    }      

  pfTemp = m_pfFFTArray;
  for (ni = 0; ni < nNeeded; ni++, pfTemp++)   {
      *pfTemp = 0.0;
  }

  // Build histogram
  
  pfData = m_pfFFTArray;
  for (ni = 0; ni < m_nNumRefs; ni++) {
      if ((m_nTwinID==-1) || (m_oList[ni].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
          if (m_oList[ni].nGetField(m_nFI_nIndexSelect)==1) {
              
              nj       = ni;
              fP      = m_oList[nj].fGetField(nFI_fFFTValue);
              nBin    = nint((fP - fPmin) / fDeltaBin);
              if ((0 <= nBin) && (nNumBins > nBin))
              {
                  // Note we use 2*nBin because of FFT below expects complex numbers
                  pfData[2*nBin] += 1.0;
              }
          }
      };
  };

  // Frequency histogram done, so FFT it
  // The FFT routine expects complex number pairs, so set the imaginary component to 0.
  // This FFT routine is from Numerical Recipes and expects pfData[1] to be the first
  // real value, so pass pfData-1 to achieve this.

  nFFT1D(nNumBins, 1, pfData);

  // Everytime we do an FFT, write the vector out for viewing with O


  // Convert real and imaginary parts to a single square of magnitude,
  // save largest values. pfData[1] should be 0.0

  pfTemp    = pfData;
  for (ni = 0; ni < nNumBins/2; ni++)
    {
      fTemp       = *pfTemp * *pfTemp++;
      fTemp      += *pfTemp * *pfTemp++;
      pfData[ni]   = fTemp;
    }

  // Find a minimum away from origin

  float fMax;
  
  if (m_nPrint) {
      FILE* pFOut;
      static int nPass;
      int nx;
      Cstring sFileName;
      sFileName = "dps";
      sFileName += nPass++;
      sFileName += ".txt";
      pFOut = fopen(sFileName.string(),"wt");
      if (pFOut) {
          for (nx=0;nx<nNumBins/2;nx++) {
              fprintf(pFOut,"%d %lf\n",nx,pfData[nx]);
          };
          fclose(pFOut);
      };
  };

  fMax = pfData[0];
  ni    = nBinTestMin;
  while ( (pfData[ni] <= fMax) && (ni < nNumBins/2) )
    {
      fMax = pfData[ni];
      ni++;
    }

  // Now find maximum in rest of Fourier data

  int nIMax = ni;
  while ( (ni < nBinTestMax) && (ni < nNumBins/2) ) {
      if (pfData[ni] > fMax)
      {
          fMax  = pfData[ni];
          nIMax = ni;
      }
      ni++;
  }

  // This gives approximate percentage of reflns indexed

  fMax = (float)sqrt((double)fMax) /  (float)m_nNumRefs * 100.0f;

  return ((double)nIMax / (double)(nNumBins-1) /  fDeltaBin * m_fWavelength);
}


int CEyeballIndex::nConvertPlanes2OrientMatrix(
        double* pfNorm1,double fSpace1,
        double* pfNorm2,double fSpace2,
        double* pfNorm3,double fSpace3,
        double a3x3fOrientMat[3][3]) {

    double a3x3fNormMat[3][3];
    double a3x3fScaleMat[3][3];
    double a3x3fNormMatInv[3][3];
    vCopyVec3D(pfNorm1,&a3x3fNormMat[0][0]);
    fNormVec3D(&a3x3fNormMat[0][0]);
    vCopyVec3D(pfNorm2,&a3x3fNormMat[1][0]);
    fNormVec3D(&a3x3fNormMat[1][0]);
    vCopyVec3D(pfNorm3,&a3x3fNormMat[2][0]);
    fNormVec3D(&a3x3fNormMat[2][0]);

    if (0.0 == fInvMat3D(&a3x3fNormMat[0][0],&a3x3fNormMatInv[0][0]))
        return 1;
    vZeroMat3D(&a3x3fScaleMat[0][0]);
    a3x3fScaleMat[0][0] = fSpace1;
    a3x3fScaleMat[1][1] = fSpace2;
    a3x3fScaleMat[2][2] = fSpace3;
    vMulMat3DMat3D(a3x3fScaleMat,a3x3fNormMatInv,a3x3fOrientMat);
    vTranMat3D(a3x3fOrientMat);

    if (fDetMat3D(&a3x3fOrientMat[0][0])<0.0) {
        // Matrix has negative determinant:  negate first vector.
        //vMulVec3DScalar(&a3x3fOrientMat[0][0],-1.0,&a3x3fOrientMat[0][0]);
    };
    return 0;
};


int CEyeballIndex::nAutoIndex(double a3x3fNonReduced[3][3]) {
    double a3fFirstView[3];
    double a3fFirstViews[2][3];
    double a3fSecondView[3];
    const int nNumSecondViews = 5;
    int nNumFoundSecondViews;
    double aa3fSecondViews[nNumSecondViews][3];
    double a3x3fPlaneNormals[3][3][3];
    double a3fDSpace[3];
    int nMinRefsAfterPrune;
    int nx;
    int nStat;
    double f0;

    // We only want the 'best' view.  Set the maximum number of 3D views to 1.
    m_nMaxBest3DVectors = 1;
    m_nMaxBest2DVectors = 2;

    // Modify minimum number of reflections after a prune?
    nMinRefsAfterPrune = m_nMinRefsAfterPrune;

    // Find the best view directions.
    nCalcRefineBestViewDirections3D();
    // Get the best view direction.
    nStat = nGetBestViewDirections3D(2,&a3fFirstViews[0][0]);
    vCopyVec3D(&a3fFirstViews[1][0],&a3fFirstView[0]);
    if (nStat != 1)
        return 1;

    // Calculate the best planes in this view.
    nCalcRefineBestViewDirections2D(a3fFirstView,TRUE);
    // Get the two best planes.
    for (nx=0;nx<2;nx++) {
        if (nGetBestPlane2D(nx,&a3x3fPlaneNormals[nx][0][0],&a3fDSpace[nx],NULL))
            return 1;
    };

    // Find the best view directions again.  
    // Presently, we are prunning spots that did not fit on the first view and best two planes.
    m_nMaxBest3DVectors = nNumSecondViews;
    nCalcRefineBestViewDirections3D();
    // Get the best view direction.
    nNumFoundSecondViews = nGetBestViewDirections3D(nNumSecondViews,&aa3fSecondViews[0][0]);
    // Get the first views that is significantly different from a3fFirstView.
    for (nx=0; nx<nNumFoundSecondViews; nx++) {
        f0 = fMinVectorDifference(&aa3fSecondViews[nx][0],&a3fFirstView[0]);
        if (f0>0.05)
            break;
    };
    if (nx==nNumSecondViews)
        return 1;

    // Calculate the best plane for this view.
    vCopyVec3D(&aa3fSecondViews[nx][0],a3fSecondView);
    nCalcRefineBestViewDirections2D(a3fSecondView);

    // Get the best plane.
    if (nGetBestPlane2D(nx,&a3x3fPlaneNormals[2][0][0],&a3fDSpace[2],NULL))
        return 1;

    // Copy the solution to the output array.
    nStat = nConvertPlanes2OrientMatrix(
        &a3x3fPlaneNormals[0][0][0],a3fDSpace[0],
        &a3x3fPlaneNormals[1][0][0],a3fDSpace[1],
        &a3x3fPlaneNormals[2][0][0],a3fDSpace[2],a3x3fNonReduced);
    if (nStat)
        return 1;
        
    // Reset these variables.
    m_nMaxBest3DVectors = m_nBest3DVectors;
    m_nMaxBest2DVectors = m_nBest2DVectors;
    m_nMinRefsAfterPrune = nMinRefsAfterPrune;

    return 0;
};


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

bool CEyeballIndex::bSameOrientMatrix(Ccrystal& oCrystal1,double fDetError) {
    int nTwin;

    oCrystal1.nCalcOrientMatrix(m_fWavelength);
    for (nTwin=1;1;nTwin++) {
        Ccrystal oCrystal(m_oHeader,nTwin);
        if (!oCrystal.bIsAvailable())
            break;
        oCrystal.nCalcOrientMatrix(m_fWavelength);
        if (oCrystal.fCalcTransForOrientMatrix(oCrystal1)<fDetError)
            return true;
    };
    return false;

};



double CEyeballIndex::fRefineSolution(double fPercentToReject,int nNumEdges,double a3x3fEdgesIn[3][3],double a3x3fEdgesOut[3][3]) {
    int nx,ny;

    double   fTestErrPercent;          // For a given test solution, the error percent.
    double   a9fRefineParams[9];       // Used to refine the edge.
    double   a81fRefineParams[81];     // Normal matrix (could be from 9 to 81 elements large).
    int     a3nCoeffs[3];             // Solution coefficients returned by bFindLC
    int     nRef;
    double   fWeight;
	double   fUnRefinedError;			// This function can be called so that it only returns the un-refined error residual.
    int      nNumValid;
    int      nNormalize;

    int     nNumNodesIndexed;
    int     nNumEdgesIndexed;
    int     nNumEdgesFound;

    Creflnlist* poList;
    double*     pfRefineErrorPercent;
    
    nNumEdgesIndexed = 0;
    nNumNodesIndexed = 0;
    nNumEdgesFound = 0;
    fWeight = 0.0;

    nNormalize = 1;



    {
        poList = &m_oList;
        if (!m_pfRefineErrorPercent)
            m_pfRefineErrorPercent = new double[poList->nGetNumReflns()];
        pfRefineErrorPercent = m_pfRefineErrorPercent;
    };

    // Figure out an appropriate error percent.

    nNumValid = 0;
	fUnRefinedError = 0.0;
    for (nRef=0;nRef<poList->nGetNumReflns();nRef++) {
        if (((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) &&
            (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==1)) {
    
            
            double a3fVec[3];
            double a3fVecD[3];
            
            a3fVec[0] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoord0);
            a3fVec[1] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoord1);
            a3fVec[2] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoord2);
            a3fVecD[0] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoordD0);
            a3fVecD[1] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoordD1);
            a3fVecD[2] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoordD2);
            fNormVec3D(a3fVecD);
            
            if (bFindLC(&fTestErrPercent,& a3x3fEdgesIn[0][0],nNumEdges,a3fVec,a3fVecD,& a3nCoeffs[0],nNormalize) && (fTestErrPercent<0.8)) {
                pfRefineErrorPercent[nRef] = fTestErrPercent;
				fUnRefinedError += 1.0/max(0.0001,fTestErrPercent);
    
            } else 
                pfRefineErrorPercent[nRef] = 0.5;
            nNumValid++;                
            
        } else 
            pfRefineErrorPercent[nRef] = 0.0;
    };
    if ((!nNumValid) || (fPercentToReject>=1.0))
        return fUnRefinedError;
    qsort(&pfRefineErrorPercent[0],poList->nGetNumReflns(),sizeof(double),double_cmp);
    nx = min(poList->nGetNumReflns() - 1,max(4,poList->nGetNumReflns() - 1 - (int) (fPercentToReject*nNumValid)));
    do {
        m_fTopErrPercent = pfRefineErrorPercent[nx];
        nx--;
    } while ((m_fTopErrPercent==0.5) && (nx>=0));
    if (nx==-1)
        return 0.0;
    ny = 0;
    do {
        nx--;
        ny++;
    } while ((nx>=0) && (ny < 5));
    if ((nx==-1) || (pfRefineErrorPercent[nx]==0.0))
        return 0.0;
    
    
    // Initialize the refinement.
    vZeroMat(9,1,&a9fRefineParams[0]);
    vZeroMat(9,9,&a81fRefineParams[0]);
    
    for (nx=0;nx<3;nx++) {
        m_a3nOdd[nx] = 0;
        m_a3nEven[nx] = 0;
        m_a3fOdd[nx] = 0.0;
        m_a3fEven[nx] = 0.0;
    };
    
    m_a2nRefined[0] = 0;
    m_a2nRefined[1] = 0;
    
    for (nRef=0;nRef<poList->nGetNumReflns();nRef++) {
        if (((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) &&
            (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==1)) {
            
            double a3fVec[3];
            double a3fVecD[3];
            
            a3fVec[0] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoord0);
            a3fVec[1] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoord1);
            a3fVec[2] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoord2);
            a3fVecD[0] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoordD0);
            a3fVecD[1] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoordD1);
            a3fVecD[2] = (*poList)[nRef].fGetField(poList->m_nFI_fRecipCoordD2);
            fNormVec3D(a3fVecD);

            m_a2nRefined[0]++;
            
            if ((bFindLC(&fTestErrPercent,& a3x3fEdgesIn[0][0],nNumEdges,a3fVec,a3fVecD,& a3nCoeffs[0],nNormalize)) && 
                (fTestErrPercent<m_fTopErrPercent)) {
                
                m_a2nRefined[1]++;
                for (nx=0;nx<3;nx++) {
                    if (a3nCoeffs[nx] % 2) {
                        m_a3nOdd[nx]++;
                        m_a3fOdd[nx] += (m_fTopErrPercent-fTestErrPercent)/m_fTopErrPercent;
                    } else {
                        m_a3nEven[nx]++;
                        m_a3fEven[nx] += (m_fTopErrPercent-fTestErrPercent)/m_fTopErrPercent;
                    };
                };
                
                fWeight += (m_fTopErrPercent-fTestErrPercent)/m_fTopErrPercent;

                /*  
                This refinement will try to get the best edges using a metric that effectively ignores the phi centroid of each spot.
                This is done by only looking at the portion of the residual vector that is orthogonal to the derivative vector.
                (Remember, for each reflection the derivative vector (D vector) denotes the change in X vector as the reflection passesthrough the difracting condition.
                 Thus, the component of the X vector along the D vector is the most uncertain, and should be ignored in this refinement).

                We are trying to deterinine a new set of basis vectors (VV=[v1 v2 v3]) given the Xobs vector (a3fVec), 
                the D vector (a3fDVec) and it's HKL index (a3nCoeffs).  Using the HKL coefficients, we can determine Xcalc
                as VV*HKL. Thus, the residual vector R is Xcalc - Xobs = VV*HKL - Xobs.
                We are interested in minimizing only the portion of R that is orthogonal to D.  Thus, computing two
                arbitrary orthgonal vectors D1 and D2, we must minimize
                
                  dot(D1,VV*HKL - Xobs) ^ 2 + dot(D2,VV*HKL - Xobs) ^ 2

                Fortunatly, this is a quadratic equation in the 9 unknown VV terms.
                This system is solved below.

                */                
                double a3x3fDBasis[3][3];
                const double fOrthoFac = 0.1;

                // Build two orthogonal vectors.
                vBuildBasis3D(a3fVecD,a3x3fDBasis);
                vNormMat3D(a3x3fDBasis);

                // Loop over each derivative.
                for (nx=0;nx<nNumEdges*3;nx++) {
                    for (ny=0;ny<nNumEdges*3;ny++) {
                        a81fRefineParams[nx + nNumEdges*3*ny] +=
                            a3x3fDBasis[1][ny % 3]*a3nCoeffs[ny/3]*
                            a3x3fDBasis[1][nx % 3]*a3nCoeffs[nx/3]+
                            a3x3fDBasis[2][ny % 3]*a3nCoeffs[ny/3]*
                            a3x3fDBasis[2][nx % 3]*a3nCoeffs[nx/3]+
                            fOrthoFac*
                            a3x3fDBasis[0][ny % 3]*a3nCoeffs[ny/3]*
                            a3x3fDBasis[0][nx % 3]*a3nCoeffs[nx/3];
                    };
                    a9fRefineParams[nx] += 
                        fDot3D(a3fVec,a3x3fDBasis[1])*a3x3fDBasis[1][nx % 3]*a3nCoeffs[nx/3] +
                        fDot3D(a3fVec,a3x3fDBasis[2])*a3x3fDBasis[2][nx % 3]*a3nCoeffs[nx/3] +
                        fOrthoFac*fDot3D(a3fVec,a3x3fDBasis[0])*a3x3fDBasis[0][nx % 3]*a3nCoeffs[nx/3]
                        ;
                };
            }; 
        };
    };    

    nx = nSolveAxB_svd(nNumEdges*3,& a81fRefineParams[0],& a3x3fEdgesOut[0][0],& a9fRefineParams[0]);
    if (fDetMat3D(&a3x3fEdgesOut[0][0])<0.0) {
        // Invert one of the edges.
        vMulVec3DScalar(&a3x3fEdgesOut[0][0],-1.0,&a3x3fEdgesOut[0][0]);
    };
    if (nx)
        return 0.0;
    else
        return fWeight/m_fTopErrPercent;
};

double CEyeballIndex::fRefineSolution(double fErrPercent,double a3x3fOrientMatrixIn[3][3],double a3x3fOrientMatrixOut[3][3],int* pnPass) {
    double fLastValue;
    double fThisValue;
    double fTopErrPercent;
    bool bConverged;
    const int nMaxPasses = 10;
    const double fConvergePercent = 0.01;
    int nPass;
    double a3x3fOrientMatrixLast[3][3];
    
    
    fLastValue = -1000;
    nPass = 0;
    
    do {
        fThisValue = fRefineSolution(fErrPercent,3,a3x3fOrientMatrixIn,a3x3fOrientMatrixOut);
        if (fThisValue == 0.0) 
            return 0.0;
        if (fThisValue<fLastValue) {
            // Possibly, the refined value is *worse*
            bConverged = true;
            vCopyMat3D(&a3x3fOrientMatrixLast[0][0],&a3x3fOrientMatrixOut[0][0]);
            break;
        };

        vCopyMat3D(&a3x3fOrientMatrixIn[0][0],&a3x3fOrientMatrixLast[0][0]);
        vCopyMat3D(&a3x3fOrientMatrixOut[0][0],&a3x3fOrientMatrixIn[0][0]);
        fTopErrPercent = m_fTopErrPercent;
        bConverged = (2.0*fabs(fLastValue-fThisValue)/(fabs(fLastValue) + fabs(fThisValue)))<fConvergePercent;
        nPass++;
        fLastValue = fThisValue;
    } while ((nPass<nMaxPasses) && (!bConverged));

    m_fTopErrPercent = fTopErrPercent;
    if (pnPass)
        *pnPass = nPass;

    if (!bConverged)
        return (m_fLastRefined = 0.0);
    else
        return (m_fLastRefined = fLastValue);
};
    
//  This routine is called if the user has a known cell that needs to be refined.  


double CEyeballIndex::fRefineSolution(double fErrPercent,
                                      Ccrystal& oCrystalIn,
                                      Ccrystal& oCrystalOut,
                                      int* pnNumberIndexedReflns) 
{
    double a3x3fOrientMatrixIn[3][3];
    double a3x3fOrientMatrixOut[3][3];
    double fValue;
    double fTopErrPercent;
    int nPass;
    int nDoubledDetected;
    int nTwinDetected;
    int nRestarts;
    

    nDoubledDetected = 0;
    nTwinDetected = 0;

    nRestarts = 0;

    while (nRestarts<=3) {
        
        
        oCrystalIn.nCalcOrientMatrix(m_fWavelength);
        oCrystalIn.vGetOrientMatrix(&a3x3fOrientMatrixIn[0][0]);
        oCrystalOut.nSetOrientMatrix(&a3x3fOrientMatrixIn[0][0],m_fWavelength);
        
        if (0.0 == (fValue = fRefineSolution(fErrPercent,a3x3fOrientMatrixIn,a3x3fOrientMatrixOut,&nPass)))
            return 0.0;
        fTopErrPercent = m_fTopErrPercent;
        if (fTopErrPercent>=0.5)
            return 0.0;
        
        
        oCrystalOut.nSetOrientMatrix(&a3x3fOrientMatrixOut[0][0],m_fWavelength);
        
        // Check for doubled axis automatically.
        
        
        if (nDoubledAxes() && (nRestarts<=3)) {
            int nPrint;
            nPrint = m_nPrint;
            m_nPrint = 0;
            nComputeDoubledAxes(fErrPercent,oCrystalOut,oCrystalIn);
            m_nPrint = nPrint;
            nDoubledDetected=1;
            nRestarts++;
            continue;
        } else
            break;
    };

    if (bSameOrientMatrix(oCrystalOut,0.05)) {
        fValue = 0.0;
        nTwinDetected = 1;
    };

    
    

    if (m_nPrint>=2) {
        printf("\n%6s %6s %6s %6s %6s %6s %8s %8s %8s %9s\n",
            "a","b","c","alpha","beta","gamma","Volume","Remarks","#Indexed","%Residual");
        printf("==============================================================================\n");
        m_nPrint = 1;
    };

    
    if (m_nPrint>=1)
    {
        printf("%6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %8.0lf %8s %8d %9.3lf\n",
               oCrystalOut.fGetCell(0),
               oCrystalOut.fGetCell(1),
               oCrystalOut.fGetCell(2),
               oCrystalOut.fGetCell(3),
               oCrystalOut.fGetCell(4),
               oCrystalOut.fGetCell(5),
               oCrystalOut.fCalcVolume(),
               nTwinDetected?"Twin":(nDoubledDetected?"Doubled?":"Okay"),
               (int)(fValue * fTopErrPercent), 
               fTopErrPercent*100.0);
    }
    
    if ( pnNumberIndexedReflns )
        *pnNumberIndexedReflns = (int)(fValue * fTopErrPercent);

    return fValue;
}


   


int CEyeballIndex::nComputeDoubledAxes(double fErrPercent,
                                       Ccrystal& oCrystalIn,
                                       Ccrystal& oCrystalOut)
{
    double a3x3fOrientMatrixIn[3][3];
    double a3x3fOrientMatrixTest[3][3];
    double a3x3fOrientMatrixBuf[3][3];
    double fCurrentResid;
    double fTestResid;
    double fBestResid;
    int nCurrentDoubledAxes;
    int nPass;
    bool bChanged;
    bool bCellChanged;

    
    oCrystalIn.nCalcOrientMatrix(m_fWavelength);
    oCrystalIn.vGetOrientMatrix(&a3x3fOrientMatrixIn[0][0]);

    bCellChanged = FALSE;
    do {
        fCurrentResid= fRefineSolution(fErrPercent,3,a3x3fOrientMatrixIn,a3x3fOrientMatrixBuf);
        vCopyMat3D(&a3x3fOrientMatrixBuf[0][0],&a3x3fOrientMatrixIn[0][0]);
        nCurrentDoubledAxes = nDoubledAxes();
        fBestResid = fCurrentResid;
        bChanged = FALSE;

        if (nCurrentDoubledAxes) {
            // Need to attempt to get halve one of the axes.

            for (nPass=0;nPass<3;nPass++) {
                vCopyMat3D(&a3x3fOrientMatrixIn[0][0],&a3x3fOrientMatrixTest[0][0]);
                // Try doubling one of the input axes. (that is, halving the orient matrix length.
                vMulVec3DScalar(&a3x3fOrientMatrixTest[nPass][0],2.0,&a3x3fOrientMatrixTest[nPass][0]);
                
                fTestResid= fRefineSolution(fErrPercent,a3x3fOrientMatrixTest,a3x3fOrientMatrixBuf);
                if ((fTestResid!=0.0) && (fTestResid>0.9*fBestResid) && (nDoubledAxes()<nCurrentDoubledAxes)) {
                    vCopyMat3D(&a3x3fOrientMatrixBuf[0][0],&a3x3fOrientMatrixIn[0][0]);
                    bChanged = TRUE;
                    if (m_nPrint)
                        printf("INFO:  Found doubled axis. Halving it.\n");
                    break;
                };
            };           
            
        } else {
            // See if a doubled axes would give a better solution.

            /*
            for (nPass=0;nPass<3;nPass++) {
                vCopyMat3D(&a3x3fOrientMatrixIn[0][0],&a3x3fOrientMatrixTest[0][0]);
                // Try doubling one of the input axes. (that is, halving the orient matrix length.
                vMulVec3DScalar(&a3x3fOrientMatrixTest[nPass][0],0.5,&a3x3fOrientMatrixTest[nPass][0]);
                
                fTestResid= fRefineSolution(fErrPercent,a3x3fOrientMatrixTest,a3x3fOrientMatrixBuf);
                if ((fTestResid!=0.0) && (fTestResid>fBestResid) && (nDoubledAxes()==0)) {
                    bChanged = TRUE;
                    vCopyMat3D(&a3x3fOrientMatrixBuf[0][0],&a3x3fOrientMatrixBest[0][0]);
                    break;
                };
            };
            if (bChanged) {
                // We found a better solution.
                vCopyMat3D(&a3x3fOrientMatrixBest[0][0],&a3x3fOrientMatrixIn[0][0]);
                bChanged = TRUE;
                if (m_nPrint)
                    printf("INFO:  Found halved axis. Doubling the axis\n");
            };
            */
        };

        bCellChanged = bCellChanged || bChanged;

    } while (bChanged);

    oCrystalOut.nSetOrientMatrix(&a3x3fOrientMatrixIn[0][0],m_fWavelength);
    
    return bCellChanged;
};



int CEyeballIndex::nDoubledAxes() {
    int nx;
    int nDoubled;

    nDoubled = 0;
    for (nx=0;nx<3;nx++) {
        if ((m_a3nOdd[nx] + m_a3nEven[nx]>=5) && (m_a3nEven[nx])) {
            if ((m_a3fOdd[nx]/((double) m_a3fOdd[nx] + m_a3fEven[nx])<=0.1) &&
                (m_a3fOdd[nx]/max(1,m_a3nOdd[nx])<m_a3fEven[nx]/m_a3nEven[nx]))
                nDoubled++;
        };
    };
    return nDoubled;
};



bool CEyeballIndex::bRefineBeam(double fRejectFraction,
                                int nBeamSearchRad,
                                int nMaxBeamMovement,
                                bool bSetNewCenter,
                                Ccrystal& oCrystalIn,
                                Ccrystal& oCrystalOut,
                                double *pdBeamPixChange,
                                double& fBestValue,     // RB added two new variables
                                Cstring& sError)
{
    fBestValue = 0.0;
    sError = ""; // initialize

    
    Crotation       oRotation(m_oHeader);
    Csource         oSource(m_oHeader);
    Cgoniometer     oGoniometer(m_oHeader,Ccrystal::ms_sCrystalPrefix);
    double          a3x3fOrientMat[3][3];
    double          a3x3fOrientMatBest[3][3];

    float           fDD[3][3];
    float           fDN[3];
    float           a3fS0[3];
    float           fGonMatrix[3][3];
    float           fInvGonMatrix[3][3];
    float           fX[3];
    int             nStat;
    int             nRef,nx;
    Crefln*         poRefln;
    Cstring         sDetectorName;
    float           f1PX,f2PX,fRot;
    float           f1MM,f2MM,f3MM;
    float           fRotVec[3];
    float           fInvRot[3][3];
    float           fVDCO[3];
    float           fTS[3];
    float           fXR[3];

    *pdBeamPixChange = 0.0;

    if( m_oList.m_nFI_fGonio1 >= 0 )
    {
        //mrp sError = "Beam check will not execute with multi-scan.\n Please use only one scan.";        
        //mrp return false;
        printf("WARNING: Beam check will not work with reflections from multiple scans.\nYou may need to rerun dtfind with an input header having scan angles that match those from the input image(s).\n"); //mrp
    }

    nStat = m_oHeader.nGetValue(Cdetector::ms_sDetectorNames, &sDetectorName);
    if( 0 != nStat )
    {
        sError = "Could not obtain detector information for beam computation.";
        return false;
    }
    
    Cdetector       oDetector(m_oHeader,sDetectorName);

    if ((!oRotation.bIsAvailable()) ||  (!oDetector.bIsAvailable()) || (!oSource.bIsAvailable()) || (!oGoniometer.bIsAvailable()))
    {
        sError = "All objects are not available for beam center computation.";
        
        return false;
    }

    if ((m_oList.m_nFI_fObsPx0<0) || (m_oList.m_nFI_fObsPx1<0) || (m_oList.m_nFI_fObsRotMid<0))
    {
        sError = "All fields are not available for beam center computation.";

        return false;
    }

    oDetector.nCalcGetDDDN(&fDD[0][0], &fDN[0]);

    oSource.vCalcGetS0(a3fS0);
    
    oGoniometer.vCalcGetRotMatrix(&fGonMatrix[0][0],oGoniometer.nGetNum(oRotation.sGetName()));
    
    if (oGoniometer.nGetRotVector(oRotation.sGetName(),&fRotVec[0]))
    {
        nStat =1;
        
        sError = "Beam check will not execute because of a problem with crystal goniometer.";

        return false;
    }   
    
    fInvMat3D(&fGonMatrix[0][0],&fInvGonMatrix[0][0]);


    // Beam shift variables.
    double      fValue;
    float       a2fBeam[2];
    
    float       a2fBestBeam[2];
    
    float       a2fOrigBeam[2];
    
    int         nPass;
    double      a2fMaxSearch[2];
    double      a2fMinSearch[2];
    double      a2fSearch[2];
    int         a2nSearchCount[2];
    double      fStep;
    bool        bOutputImage = true;
    bool        bMoved = FALSE;
    const int   nMaxRefs = 100;
    int         nAvailableRefs;
    Cimage      oImage(nBeamSearchRad*2+1,nBeamSearchRad*2+1,eImage_realIEEE);
    int         nCumulCount,nCumulCountRange;

   
    oDetector.m_poSpatial->nGetBeamPosition(&a2fOrigBeam[0], &a2fOrigBeam[1]);

    printf("Original (input header) beam center: [%5.1lf %5.1lf]\n", a2fOrigBeam[0],a2fOrigBeam[1]);

    // Do a count of the number of vectors used.  We will reduce this number
    // if it is too large.
    nx = 0;
    nAvailableRefs = 0;
    
    for (nRef = 0; nRef<m_nNumRefs;nRef++)
    {
        m_oList[nRef].vSetField(m_nFI_nIndexSelect,0);
        
        if( m_nTwinID ==-1 || m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID) == m_nTwinID ) 
            nAvailableRefs++;
    }
    
    do 
    {
        int     nRandomNumber = rand();
        
        nRef = nRandomNumber % m_nNumRefs;
        
        if ((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) {
            if (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==0) {
                m_oList[nRef].vSetField(m_nFI_nIndexSelect,1);
                nx++;
            };
        }
    } while (nx<min(nAvailableRefs,nMaxRefs));

    for (nPass = 0; nPass < 2; nPass++) {
        
        printf("\n");

        if (nPass == 0) {
            fStep = 1.0;
            a2fMinSearch[0] = max(0,((int) a2fOrigBeam[0]) - nBeamSearchRad);
            a2fMinSearch[1] = max(0,((int) a2fOrigBeam[1]) - nBeamSearchRad);
            a2fMaxSearch[0] = ((int) a2fOrigBeam[0]) + nBeamSearchRad;
            a2fMaxSearch[1] = ((int) a2fOrigBeam[1]) + nBeamSearchRad;
        } else {
            fStep = 0.1;
            a2fMinSearch[0] = a2fBestBeam[0] - 1.0;
            a2fMinSearch[1] = a2fBestBeam[1] - 1.0;
            a2fMaxSearch[0] = a2fBestBeam[0] + 1.0;
            a2fMaxSearch[1] = a2fBestBeam[1] + 1.0;
        };
            printf("%s Pass.  Search dim0 in [%d %d] dim1 in [%d %d]\n",
                (nPass==0)?"First":"Second",
                (int) a2fMinSearch[0],(int) a2fMaxSearch[0],(int) a2fMinSearch[1],(int) a2fMaxSearch[1]);

        
        if (bOutputImage) {
            int nx,ny;            
            for (nx=0;nx<nBeamSearchRad*2+1;nx++) {
                for (ny=0;ny<nBeamSearchRad*2+1;ny++) {
                    oImage.nSetPixel(nx,ny,(float) 0.0);
                };
            };
        };
        
        fBestValue = -1.0;
        nCumulCount=0;
        nCumulCountRange = (int) ((a2fMaxSearch[0] - a2fMinSearch[0])*(a2fMaxSearch[1] - a2fMinSearch[1]));
        for (a2fSearch[0] = a2fMinSearch[0],a2nSearchCount[0] = 0; a2fSearch[0]<=a2fMaxSearch[0]; a2fSearch[0] += fStep,a2nSearchCount[0]++) {
            for (a2fSearch[1] = a2fMinSearch[1],a2nSearchCount[1] = 0; a2fSearch[1]<=a2fMaxSearch[1]; a2fSearch[1] += fStep, a2nSearchCount[1]++,nCumulCount++) { 
                if (!(nCumulCount % max(1,(nCumulCountRange/40))) && (!nPass))
                    printf(".");

                a2fBeam[0] = a2fSearch[0];
                a2fBeam[1] = a2fSearch[1];
                oDetector.m_poSpatial->nSetBeamPosition(
                    a2fBeam[0],
                    a2fBeam[1]
                    );           
                
                
                // Recompute normalized X vectors.
                
                for (nRef = 0; nRef < m_nNumRefs; nRef++)
                {
                    if (((m_nTwinID==-1) || (m_oList[nRef].nGetField(m_oList.m_nFI_nTwinID)==m_nTwinID)) &&
                        (m_oList[nRef].nGetField(m_nFI_nIndexSelect)==1)) {
                        
                        poRefln = m_oList.poGetRefln(nRef);       
                        f1PX = poRefln->fGetField(m_oList.m_nFI_fObsPx0);
                        f2PX = poRefln->fGetField(m_oList.m_nFI_fObsPx1);
                        fRot = poRefln->fGetField(m_oList.m_nFI_fObsRotMid);       
                        nStat = oDetector.m_poSpatial->nPixeltoMM(f1PX, f2PX, &f1MM, &f2MM, &f3MM);       
                        vConvRotVec3DMat3D(-fRot, fRotVec, &fInvRot[0][0]);       
                        fVDCO[0] = f1MM;
                        fVDCO[1] = f2MM;
                        fVDCO[2] = f3MM;
                        vMulMat3DVec3D(fDD, fVDCO, fTS);
                        (void) fNormVec3D(fTS);
                        vAddVec3DVec3D(a3fS0, fTS, fXR);
                        
                        vMulMat3DVec3D(fInvRot, fXR, fTS);      // fTS is a temp var here
                        vMulMat3DVec3D(fInvGonMatrix, fTS, fX); // fTS is a temp var here                
                        
                        poRefln->vSetField(m_oList.m_nFI_fRecipCoord0, fX[0]);
                        poRefln->vSetField(m_oList.m_nFI_fRecipCoord1, fX[1]);
                        poRefln->vSetField(m_oList.m_nFI_fRecipCoord2, fX[2]);
                        
                        float a3fTemp1[3];
                        float a3fTemp2[3];
                        
                        // Get the portion of the X vector perpendicular to the rotation vector.
                        vMulVec3DScalar(fRotVec,fDot3D(fX,fRotVec),a3fTemp1);
                        vSubVec3DVec3D(fX,a3fTemp1,a3fTemp1);
                        vCross3D(fRotVec,a3fTemp1,a3fTemp2);
                        poRefln->vSetField(m_oList.m_nFI_fRecipCoordD0, a3fTemp2[0]);
                        poRefln->vSetField(m_oList.m_nFI_fRecipCoordD1, a3fTemp2[1]);
                        poRefln->vSetField(m_oList.m_nFI_fRecipCoordD2, a3fTemp2[2]);
                    };
                };
                
                
                // Save the orientation matrix, since fRefineSolution might modifiy oCrystalIn.
                oCrystalIn.nCalcOrientMatrix(m_fWavelength);
                oCrystalIn.vGetOrientMatrix(&a3x3fOrientMat[0][0]);
                
                // Calculate the heuristic value.           
                fValue = fRefineSolution(fRejectFraction,oCrystalIn,oCrystalOut);
                
                
                if ((fBestValue==-1.0) || (fValue>fBestValue)) {
                    fBestValue = fValue;
                    a2fBestBeam[0] = a2fBeam[0];
                    a2fBestBeam[1] = a2fBeam[1];
                    oCrystalOut.nCalcOrientMatrix(m_fWavelength);
                    oCrystalOut.vGetOrientMatrix(&a3x3fOrientMatBest[0][0]);
                };
                
                oCrystalIn.nSetOrientMatrix(&a3x3fOrientMat[0][0],m_fWavelength);
                if (bOutputImage)
                    oImage.nSetPixel(a2nSearchCount[0],a2nSearchCount[1],(float) fValue);
            };
        };
        oCrystalOut.nSetOrientMatrix(&a3x3fOrientMatBest[0][0],m_fWavelength);
        
        if (nPass == 0) 
        {
            
            if ((a2fOrigBeam[0] - a2fBestBeam[0])*(a2fOrigBeam[0] - a2fBestBeam[0]) +
                ((a2fOrigBeam[1] - a2fBestBeam[1])*(a2fOrigBeam[1] - a2fBestBeam[1])) < 
                (nMaxBeamMovement*nMaxBeamMovement))
            {
                continue;
            } 
            else
            {
                printf("\nBeam center is not within %d pixels of input beam center [%d %d]\n"
                    "The index solution will not work.\n",nMaxBeamMovement,(int) a2fOrigBeam[0],(int) a2fOrigBeam[1]);
                fBestValue = 0.0;
                break;
            }
        } 
        else
        {
            printf("Calculated pre-reduced cell solution is in agreement with detector beam center!\n");
            printf("Original (input header) beam center: [%5.1lf %5.1lf]\n",a2fOrigBeam[0],a2fOrigBeam[1]);
            printf("New      (calculated)   beam center: [%5.1lf %5.1lf]\n",a2fBestBeam[0],a2fBestBeam[1]);            

            double dX, dY;
            dX = a2fOrigBeam[0]-a2fBestBeam[0];
            dY = a2fOrigBeam[1]-a2fBestBeam[1];
            *pdBeamPixChange =  sqrt(fabs(dX * dX  + dY * dY));
            if (bSetNewCenter)
            {
                printf("Header updated to reflect beam center change.\n");
                oDetector.m_poSpatial->nSetBeamPosition(a2fBestBeam[0],a2fBestBeam[1]);
                oDetector.m_poSpatial->nUpdateHeader(&m_oHeader);
            }
        }
    }

    // Select all reflections again.
    for (nRef = 0; nRef<m_nNumRefs;nRef++)
    {
        m_oList[nRef].vSetField(m_nFI_nIndexSelect,1);
    }

    if( fBestValue == 0.0 )
    {
        printf("WARNING: beam center refinement failed.\n");
        printf("Restoring the original beam center position in the header.\n");
        
        oDetector.m_poSpatial->nSetBeamPosition(a2fOrigBeam[0], a2fOrigBeam[1]);
        oDetector.m_poSpatial->nUpdateHeader(&m_oHeader);
        
        if( bOutputImage )
        {
            printf("Beam center map written to 'beamcenter.img'\n");
          
        }
    
        return false;
    }

    return true;
}


int CEyeballIndex::nAddTweaks() {
	// +jwp 19-Dec-2002
    // Examine the reflns and if they appear to come from a single image,
    // let us 'randomize' the fObs_rot_mid values
	int nFI_fObsRotMid;
	int nx;
	const bool bPrintTweaked = false;

	nFI_fObsRotMid = m_oList.m_nFI_fObsRotMid;
	if (0 < nFI_fObsRotMid)
	{
		float fObsRot, fWidth, fTemp;
		Crefln *poRefln;
		int   ny, nStart, nEnd;
		nx = 0;
		nStart = 0;
		while (nx < m_oList.nGetNumReflns())
		{
			//printf("while nx: %d\n", nx);
			poRefln = m_oList.poGetRefln(nx);
			fObsRot = poRefln->fGetField(nFI_fObsRotMid);
			ny     = 1;
			nStart = nx;
			nx++;
			for ( ; nx < m_oList.nGetNumReflns(); nx++)
			{
				poRefln = m_oList.poGetRefln(nx);
				fTemp   = poRefln->fGetField(nFI_fObsRotMid);
				if (fTemp != fObsRot)
					break;
				ny++;
			}
			if (100 < ny)
			{
				// At least 100 reflns have same fObsRotMid value
				
				poRefln = m_oList.poGetRefln(nStart);
				fWidth  = poRefln->fGetField(m_oList.m_nFI_fObsRotWidth);
				fWidth  = min(1.0, fWidth * 1.2);
				nEnd    = min(nStart+ny, m_oList.nGetNumReflns());
				printf("INFO: reflns %d to %d appear to come from a single image, tweaking\n"
					"      estimated observed_rot_mid values (-notweakobs disables this).\n", nStart, nEnd-1);
				for (nx = nStart; nx < nEnd; nx++)
				{
					//			printf("tweak nx: %d\n", nx);
					poRefln = m_oList.poGetRefln(nx);
					fTemp   = poRefln->fGetField(nFI_fObsRotMid);
					fTemp   = fTemp - (fWidth * 0.5) 
						+ (float) (nx % 5) * 0.25 * fWidth;
					poRefln->vSetField(nFI_fObsRotMid, fTemp);
				}
			}
		} // endwhile
	} // if fobsrotmid exists
	if (bPrintTweaked)
	{
		m_oList.nWrite(sDtrekGetPrefix() + "tweaked.ref");
	}
	return 0;
};



int CEyeballIndex::nSearchForTwinLaw(Ccrystal& oCrystalIn,
                                     Ccrystal& oCrystalOut,
                                     int nMaxIndex,
                                     double fPercentReject,
                                     double fMinTwinFraction)
{
	double a3fTempVec[3];
	int nx;
	
	double a3x3fOrientIn[3][3];
	double a3x3fRealIn[3][3];
    double fMinOrientRealVec;
	double a3fRealRot[3];

    double f0;
	const double fMinRot = 2.0;
	const double fRotStep = 0.5;

    
	double a3x3fRotOrientIn[3][3];
	double a3x3fRotOrientOut[3][3];
	
	double fResid,fBestResid,fBestResidThisHKL;
	double a3x3fRotMat[3][3],a3x3fBestRotMat[3][3],a3x3fBestRotMatThisHKL[3][3]; 
    int    a3nRotHKL[3],a3nBestRotHKL[3],a3nBestRotHKLThisHKL[3];
	double fRot,fBestRot,fBestRotThisHKL;
    int nTwin;
    int nTwinLaw;

    Cimage_header oTempHeader(m_oHeader);
    Ccrystal oTempCrystal;

    nTwin = oCrystalIn.m_nTwin;
    nTwinLaw = oCrystalIn.m_nTwinLaws + 1;
    if (nTwinLaw > g_nMaxTwinLaws) {
        printf("INFO: Crystal component #%d already has %d twin laws (the maximum allowed)\n",nTwin,nTwinLaw);
        return 1;
    };

    
	// Grab the orientation matrix for the known crystal.
	oCrystalIn.nCalcOrientMatrix(1.0);
    oCrystalIn.vGetOrientMatrix(&a3x3fOrientIn[0][0]);
    fInvMat3D(&a3x3fOrientIn[0][0],
        &a3x3fRealIn[0][0]);
    vTranMat3D(a3x3fRealIn);
   
    

	
	// Calculate the minimum reciprocal vector.
	fMinOrientRealVec = 0.0;
	for (nx =0; nx < 3; nx++) {
		if ((nx==0) || (fLenVec3D(a3x3fOrientIn[nx])<fMinOrientRealVec))
			fMinOrientRealVec = fLenVec3D(a3x3fRealIn[nx]);
	};
    
    // The first count (a3nRotHKL[0]) only starts at 0 since we can rotate about a vector or
    // it's negative.
    fBestResid = 0.0;
    for (a3nRotHKL[0] = 0; a3nRotHKL[0] <= nMaxIndex; a3nRotHKL[0]++) {
        for (a3nRotHKL[1] = -nMaxIndex; a3nRotHKL[1] <= nMaxIndex; a3nRotHKL[1]++) {
            for (a3nRotHKL[2] = -nMaxIndex; a3nRotHKL[2] <= nMaxIndex; a3nRotHKL[2]++) {
                fBestResidThisHKL = 0.0;
                if (!bDivides(3,&a3nRotHKL[0])) {
                    // Build the rotation vector.
                    vZeroMat(3,1,&a3fRealRot[0]);
                    for (nx = 0; nx < 3; nx++) {
                        vMulVec3DScalar(a3x3fRealIn[nx],a3nRotHKL[nx],a3fTempVec);
                        vAddVec3DVec3D(a3fRealRot,a3fTempVec,a3fRealRot);
                    };
                    // We must avoid a situation where the vector we created is identicaly zero.
                    if (fLenVec3D(a3fRealRot)>0.98*fMinOrientRealVec) {
                        // Normalize this vector 
                        fNormVec3D(a3fRealRot);
                        
                        // Run through each rotation, searching for the best solution.
                        for (fRot = fMinRot; fRot <= 360.0 - fMinRot; fRot += fRotStep) {
                            // Refine this new twin law, which is a rotation about a vector
                            vConvRotVec3DMat3D(fRot,&a3fRealRot[0],&a3x3fRotMat[0][0]);
                            vMulMat3DMat3D(a3x3fRotMat,a3x3fOrientIn,a3x3fRotOrientIn);
                            fResid = fRefineSolution(1.0,3,a3x3fRotOrientIn,a3x3fRotOrientOut);
                            if (fResid > fBestResidThisHKL) {
                                fBestResidThisHKL = fResid;
                                fBestRotThisHKL = fRot;
                                vCopyMat3D(&a3x3fRotMat[0][0],&a3x3fBestRotMatThisHKL[0][0]);
                                a3nBestRotHKLThisHKL[0] = a3nRotHKL[0];
                                a3nBestRotHKLThisHKL[1] = a3nRotHKL[1];
                                a3nBestRotHKLThisHKL[2] = a3nRotHKL[2];
                            };
                            if (fResid > fBestResid) {
                                fBestResid = fResid;
                                fBestRot = fRot;
                                vCopyMat3D(&a3x3fRotMat[0][0],&a3x3fBestRotMat[0][0]);
                                a3nBestRotHKL[0] = a3nRotHKL[0];
                                a3nBestRotHKL[1] = a3nRotHKL[1];
                                a3nBestRotHKL[2] = a3nRotHKL[2];
                            };
                        };
                    };
                };
                if (fBestResidThisHKL > 0.0) {
                    printf("For HKL = [ %d %d %d ]\n",a3nRotHKL[0],a3nRotHKL[1],a3nRotHKL[2]);
                    printf("Resid = %8.3lf Rot = %8.2lf\n",fBestResidThisHKL,fBestRot);
                };
            };
        };
    };
    
    
    double a3x3fOrientInInv[3][3];
    int nRefineCycles = 5;
    int nCycle;
    // Not sure whether this really does any good.  It might do some harm actually.

    
    for (nCycle = 0; nCycle < nRefineCycles; nCycle++) {
        
        // Refine this solution.
        vMulMat3DMat3D(a3x3fBestRotMat,a3x3fOrientIn,a3x3fRotOrientIn);
        fResid = fRefineSolution(fPercentReject,3,a3x3fRotOrientIn,a3x3fRotOrientOut);
        fInvMat3D(&a3x3fOrientIn[0][0],&a3x3fOrientInInv[0][0]);
        vMulMat3DMat3D(a3x3fRotOrientOut,a3x3fOrientInInv,a3x3fBestRotMat);
        vNormMat3D(a3x3fBestRotMat);
    };
    

    oCrystalIn.nUpdateHeader(&oTempHeader);
    oCrystalOut.nInitValues(oTempHeader,oCrystalIn.m_nTwin);
    oCrystalOut.vSetTwinLaw(oCrystalOut.m_nTwinLaws + 1,&a3x3fBestRotMat[0][0],true);
    oCrystalOut.nUpdateHeader(&oTempHeader);

    m_oList.nCalcRecipLatticePoints(oTempHeader,f0 = 10.0);

      
    int nWorked;


    // Discover how many reflections this new solution actually indexed with the new twin law.
    nWorked = 0;
    for (nx = 0; nx < m_oList.nGetNumReflns(); nx++) {
        if (m_oList[nx].nGetField(m_oList.m_nFI_nTwinID) == nTwin + 1000*nTwinLaw)
            nWorked++;
    };

    if (nWorked > m_oList.nGetNumReflns()*fMinTwinFraction) {
        printf("INFO:  %d reflection out of %d (%.1lf%%) worked with new indexing solution\n",nWorked,m_oList.nGetNumReflns(),
            100.0*nWorked/m_oList.nGetNumReflns());
        m_oList.nList(oTempHeader);
        return 0;
    } else
	    return 1;
};




// The following function is specifically for use in the Lattice viewer.  
// Here is a list of arguments and their purpose.
// nPoints -- number of points in dPoints
// aa3fPointsSelectedInPlate -- reciprocal coordinates from the selected points.  This turns out to be easier than some kind of selection list since I already create this array.
// Output:
// a3fPlaneNormal[] -- <normvector> 
// a3fPlaneOffset[] -- <offsetvector> for the best fit plane (<offsetvector> should be the average of the points that are within some tolerance of the plane).
// fInterPlateDistance -- distance to the next parallel plane.



int CEyeballIndex::nCalcFitPlane(int nPoints, double aa3fPointsSelectedInPlate[][3], double a3fPlaneNormal[3],double a3fPlaneCentroid[3],double* pfInterPlaneDistance) {
    itr<double> afConst;
    itr<double> aafFuncs[4];
    double*     apfFuncs[4];
    itr<int>    anOutlierReject;
    double      a3fArgs[3];
    double      a3fArgsSigma[3];
    double      f0;
    double      fContribLength;
    int         nContrib;
    int         nPoint;
    int         nx;

    
    // Determine the 'constant' value (this in theory can be any #, but we want the
    // least squares to be stable.
    nContrib = 0;
    fContribLength = 0.0;
    for (nPoint = 0; nPoint < nPoints; nPoint++) {
        double a3fTemp[3];
        for (nx = 0; nx < 3; nx++)
            a3fTemp[nx] = aa3fPointsSelectedInPlate[nPoint][nx];
        fContribLength += fLenVec3D(a3fTemp);
        nContrib++;
    };
    fContribLength/=nContrib;


    for (nPoint = 0; nPoint < nPoints; nPoint++) {
        afConst + fContribLength;
        for (nx = 0; nx < 3; nx++)
            aafFuncs[nx] + aa3fPointsSelectedInPlate[nPoint][nx];        
    };
    for (nx = 0; nx < 3; nx++)
        apfFuncs[nx] = &aafFuncs[nx][0];
    
    f0 = ms_fLeastSquaresOutlierRejectIoverSig;
    ms_fLeastSquaresOutlierRejectIoverSig = 3.0;
    if (nSolveLeastSquares(3,nPoints,&afConst[0],NULL,&apfFuncs[0],a3fArgs,a3fArgsSigma,NULL,NULL,&anOutlierReject))
        return 1;
    ms_fLeastSquaresOutlierRejectIoverSig = f0;

    vZeroMat(3,1,&a3fPlaneCentroid[0]);
    nContrib = 0;
    for (nPoint = 0; nPoint < nPoints; nPoint++) {
        if (!anOutlierReject[nPoint]) {
            for (nx = 0; nx < 3; nx++)
                a3fPlaneCentroid[nx] += aa3fPointsSelectedInPlate[nPoint][nx];
            nContrib++;
        };
    };
    if (nContrib) {
        for (nx = 0; nx < 3; nx++)
            a3fPlaneCentroid[nx]/=nContrib;
    };
    
    vCopyVec3D(a3fArgs,a3fPlaneNormal);
    fNormVec3D(a3fPlaneNormal);

    if (pfInterPlaneDistance)
        *pfInterPlaneDistance = fDPSFFT(a3fPlaneNormal);

    return 0;
};

