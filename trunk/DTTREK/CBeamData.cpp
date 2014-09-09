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
// CBeamData.cc     Initial author: Thaddeus Niemeyer        15-Nov-2001
//    Beam Center and Beam Stop tools
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

#include "CBeamData.h"

#include "Cstat.h"
#include "raxis.h"
#include "dtrekvec.h"
#include "dtsvd.h"
#include "dtarray.h"
#include "CExperimentalStrategy.h"

#define WRITETEST(pcTitle,nMax,fX,fY) { FILE* pFTest; int nx; pFTest = fopen(pcTitle,"w+t"); for (nx=0;nx<nMax;nx++) fprintf(pFTest,"%f %f\n",(float) (fX),(float) (fY)); fclose(pFTest); }

Cbeamdata::Cbeamdata(Cimage& _oCi,int nCompression): 
m_oCompImage((nCompression==1)?1:(_oCi.nGetDimension(0)/nCompression),
             (nCompression==1)?1:(_oCi.nGetDimension(1)/nCompression),
             eImage_uI2),
m_oCi((nCompression==1)?_oCi:m_oCompImage),
m_oBackground(m_oCi)
{ 
    m_pcInt = NULL;
    m_pcCalcInt = NULL;
    m_pcCalcInt2 = NULL;

    m_nDim[0]=m_oCi.nGetDimension(0);
    m_nDim[1]=m_oCi.nGetDimension(1);

    m_ppfSpotValue = NULL;
    m_ppfSpotCount = NULL;

    m_poSource = NULL;
    m_poDetector = NULL;
    m_nCompression = nCompression;

    nCompress(_oCi,nCompression);
}

int Cbeamdata::nCompress(Cimage& oOld,int nComp) {
    int nn0,nn1;
    int no0,no1;
    int nDim[2];
    double fAverage;
    int nAverage;

    nDim[0] = oOld.nGetDimension(0);
    nDim[1] = oOld.nGetDimension(1);

    if (nComp>1) {
        // Average data from original image into the new image.
        for (nn0 = 0; nn0 < m_nDim[0]; nn0++) {
            for (nn1 = 0; nn1 < m_nDim[1]; nn1++) {
                fAverage = 0.0;
                nAverage = 0;
                for (no0 = 0; no0 < nComp; no0++) {
                    for (no1 = 0; no1 < nComp; no1++) {
                        if ((no0 + nn0*nComp < nDim[1]) && (no1 + nn1*nComp < nDim[1])) {
                            fAverage +=(oOld.*oOld.prfGetPixel)(no0 + nn0*nComp, no1 + nn1*nComp);
                            nAverage++;
                        }
                    }
                }
                (m_oCi.*m_oCi.prnSetPixel)(nn0,nn1,fAverage/max(1,nAverage));
            }
        }
        m_oCi.m_oHeader = oOld.m_oHeader;
    }
    return 0;
}


Cbeamdata::~Cbeamdata() { 
    int nx;

    if (m_pcInt)
        delete[] m_pcInt;
    if (m_pcCalcInt)
        delete[] m_pcCalcInt;
    if (m_pcCalcInt2)
        delete[] m_pcCalcInt2;
    if (m_ppfSpotValue) {
        for (nx=0;nx<m_nSpotDim;nx++) {
            delete[] m_ppfSpotValue[nx];
        }
        delete[] m_ppfSpotValue;
    }
    if (m_ppfSpotCount) {
        for (nx=0;nx<m_nSpotDim;nx++) {
            delete[] m_ppfSpotCount[nx];
        }
        delete[] m_ppfSpotCount;
    }

    if (m_poSource)
        delete m_poSource;
    if (m_poDetector)
        delete m_poDetector;

}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

float Cbeamdata::fRmerge(int nCenter0,int nCenter1,bool bUseLineRMerge) {
    int npx0;
    int npx1;
    int ncpx0;
    int ncpx1;

    float fDiffSum;
    float fTotalSum;
    float fVal,fValc;

    // We will do all pixels ABOVE the estimated center. 
    fDiffSum=0;
    fTotalSum=0;

    for (npx1=0;npx1<nCenter1;npx1++) {
        m_oCi.nSetNextPixel(0,npx1);
        for (npx0=0;npx0<m_nDim[0];npx0++) {
            fVal = (m_oCi.*m_oCi.prfGetPixel)(npx0,npx1);
            if (bUseLineRMerge) {
                ncpx0 = npx0;
                ncpx1 = 2*nCenter1 - npx1;
            } else {
                ncpx0 = 2*nCenter0 - npx0;
                ncpx1 = 2*nCenter1 - npx1;
            }
            // Make sure the complementary pixel is within range.
            if ((ncpx0>=0) && (ncpx0<m_nDim[0]) && (ncpx1>=0) && 
                (ncpx1<m_nDim[1]) && (fVal) && (m_oBackground.rcGetBackground(npx0,npx1)==0) &&
                (abs(npx0-nCenter0)>=100) && (abs(npx1-nCenter1)>=100) 
                ) {
                fValc = (m_oCi.*m_oCi.prfGetPixel)(ncpx0,ncpx1);
                // See if either pixel is masked out or zero.
                if ((fValc) && (m_oBackground.rcGetBackground(ncpx0,ncpx1)==0)) {
                    // Record the difference.
                    fTotalSum+=fVal+fValc;
                    fDiffSum+=ABS(fVal - fValc);

                }
            }
        }
    }
return ((float) fDiffSum)/((float) fTotalSum);
}

int Cbeamdata::nRejectRectangle(int npx00,int npx01,int npx10,int npx11) {
    int nx,ny;

    for (nx=min(npx00,npx10);nx<=max(npx00,npx10);nx++)
        for (ny=min(npx01,npx11);ny<=max(npx01,npx11);ny++) 
            m_oBackground.rcGetBackground(nx,ny)=1;
return 0;
}


int Cbeamdata::nFindLowestRmerge(int& nValue0,int& nValue1) {

    Cstring sTemp;
    float f0;
    int nx,ny;
    int nTemp0,nTemp1;
    int nBoxCenter0,nBoxCenter1;
    int nVerbose = 1;

    
    

    float fSearchVal[20][20];

    float fSearchMem[21*21];
    int   nSearchMem0[21*21];
    int   nSearchMem1[21*21];
    int   nMem;
    int   nSearchDelta0 = 2;
    int   nSearchDelta1 = 2;
    int   nSearchDim0;
    int   nSearchDim1;
    int   nSearchMem;

    if ((!m_oCi.m_oHeader.nGetValue(D_K_DetectorType, 1, &sTemp)) && (sTemp.contains("R-AXIS-CS"))) {
        // It's a rapid image.  Restrict the 0 pixel search.
        nSearchDelta0 = 0;
    }

    nSearchDim0 = nSearchDelta0*2 + 1;
    nSearchDim1 = nSearchDelta1*2 + 1;
    nSearchMem = nSearchDim0*nSearchDim1;


    printf("Finding R-merge...\n");

    if (nValue0!=-1)
        nBoxCenter0=nValue0;
    else
        nBoxCenter0=m_nDim[0]/2;
    if (nValue1!=-1)
        nBoxCenter1=nValue1;
    else
        nBoxCenter1=m_nDim[1]/2;
    for (nx=0;nx<nSearchMem;nx++) {
        nSearchMem0[nx] = -1;
        nSearchMem1[nx] = -1;
    }

    while (1) {
        
        for (nx = -nSearchDelta0;nx <= nSearchDelta0;nx++) {
            for (ny = -nSearchDelta1;ny <= nSearchDelta1;ny++) {
                
                for (nMem=0;nMem<nSearchMem;nMem++) {
                    if ((nSearchMem0[nMem]==nBoxCenter0+nx) && (nSearchMem1[nMem]==nBoxCenter1+ny)) {
                        fSearchVal[nx+nSearchDelta0][ny+nSearchDelta1] = fSearchMem[nMem];
                        break;
                    }
                }
                if (nMem==nSearchMem) {
                    f0 = fRmerge(nBoxCenter0+nx,nBoxCenter1+ny,nSearchDelta0==0);
                    fSearchVal[nx+nSearchDelta0][ny+nSearchDelta1]=f0;
                    if (nVerbose>=2)
                        printf("[%d %d] = %f\n",nBoxCenter0+nx,nBoxCenter1+ny,f0);
                }
            }
            if (nVerbose>=2)
                printf("\n");
        }
        
        // Print the R-merge value table.
        if (nVerbose>=1) {
            printf("R-merge Values\n       ");
            for (nx=-nSearchDelta0;nx<=nSearchDelta0;nx++)
                printf("%6d ",nBoxCenter0+nx);
            printf("\n");
            for (ny=-nSearchDelta1;ny<=nSearchDelta1;ny++) {
                printf("%6d ",nBoxCenter1+ny);
                for (nx=-nSearchDelta0;nx<=nSearchDelta0;nx++)
                    printf("%6.4f ",fSearchVal[nx+nSearchDelta0][ny+nSearchDelta1]);
                printf("\n");
            }
        }

        // Look for the lowest R-merge value.
        f0=1e20f;
        for (nx=0;nx<nSearchDim0;nx++) {
            for (ny=0;ny<nSearchDim1;ny++) {
                if (fSearchVal[nx][ny]<f0) {
                    nTemp0=nx;
                    nTemp1=ny;
                    f0=fSearchVal[nx][ny];
                }
            }
        }
        
        
        // If the pixel in the center of the box is the lowest one,
        // then we are finished.
        
        if ((nTemp0==nSearchDelta0) && (nTemp1==nSearchDelta1))
            break;
        // We are not done


        
        double* pfA = NULL;
        double* pfb = NULL;
        double pfVec[2];
        double pfQuadb[2];
        double pfQuadA[2][2];
        double pfQuadAInv[2][2];
        double pfShift[2];
        double fConst;
        nInitQuadratic(2,true,&pfA,&pfb);
        for (nx=0;nx<nSearchDim0;nx++) {
            for (ny=0;ny<nSearchDim1;ny++) {
                pfVec[0] = nx- nSearchDelta0;
                pfVec[1] = ny- nSearchDelta1;
                nAddQuadratic(2,fSearchVal[nx][ny],1.0,pfVec,pfA,pfb);
            }
        }
        nSolveQuadratic(2,&pfQuadA[0][0],&pfQuadb[0],&fConst,pfA,pfb);
        pfQuadAInv[0][0] = pfQuadA[1][1];
        pfQuadAInv[1][1] = pfQuadA[0][0];
        pfQuadAInv[0][1] = -pfQuadA[0][1];
        pfQuadAInv[1][0] = -pfQuadA[1][0];
        pfShift[0] = -0.5*(pfQuadb[0]*pfQuadAInv[0][0] + pfQuadb[1]*pfQuadAInv[1][0])/(pfQuadA[0][0]*pfQuadA[1][1] - pfQuadA[0][1]*pfQuadA[1][0]);
        pfShift[1] = -0.5*(pfQuadb[0]*pfQuadAInv[0][1] + pfQuadb[1]*pfQuadAInv[1][1])/(pfQuadA[0][0]*pfQuadA[1][1] - pfQuadA[0][1]*pfQuadA[1][0]);
        delete[] pfA;
        delete[] pfb;



  

        // Recompute new box center in nTemp0 and nTemp1
        nTemp0+=nBoxCenter0-nSearchDelta0;
        nTemp1+=nBoxCenter1-nSearchDelta1;

        // Copy over the values into the menirt save array.
        nMem = 0;
        for (nx = -nSearchDelta0;nx <= nSearchDelta0;nx++) {
            for (ny = -nSearchDelta1;ny <= nSearchDelta1;ny++) {
                nSearchMem0[nMem] = nBoxCenter0 + nx;
                nSearchMem1[nMem] = nBoxCenter1 + ny;
                fSearchMem[nMem] = fSearchVal[nx+nSearchDelta0][ny+nSearchDelta1];
                nMem++;
            }
        }
               
        // New values.
        nBoxCenter0 = nTemp0;
        nBoxCenter1 = nTemp1;
            
            
    }
    printf("Lowest R-merge computed at [%d %d]\n",nBoxCenter0,nBoxCenter1);
    nValue0=nBoxCenter0;
    nValue1=nBoxCenter1;

    return 0;

}






////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


int Cbeamdata::nBuildAveragePeak(int nDataPerPixel,float fRejectPercent) {

    int nPass;
    int npx[2];

    const double fSpotSizeSigmaReject = 2.0;            // Sigma reject for spot sizes.  (i.e., don't use spots that deviate from average spot size by this amount).
    const double fNormSpotIntensity =  1000.0;  // All integrated spot intensites normalized to this value.
    Cstat oStat;
    float fSpotSizeAverage;
    float fSpotSizeSigma;

    const int nMaxInSpot = 1000;
    static int nSpotStack0[nMaxInSpot+10];
    static int nSpotStack1[nMaxInSpot+10];
    static int nSpot0[nMaxInSpot];
    static int nSpot1[nMaxInSpot];
    int nStackPoint;
    int nPixelsInSpot;

    int nSpotsUsed0;
    int nSpotsUsed1;        // Spots used in first integration pass (had spot size within deviation, and had some non-zero intensity)
    int nSpotsUsed2;
    int nSpotsUsed3;        // Spots used in second integration pass (were used in first integration pass, and fit average profile within deviation).
    int nSpotCount;
    bool bSpotRejected;     // Used in last pass to unset pixels for rejected spots.

    double** ppfTestCount;  // These are used as a temporary buffer.
    double** ppfTestValue;
    int*     pnReject;      // Reject flag used in fourth pass.
    float*   pfDev;         // Deviations used to compute the reject flags.


    int nx,ny,nz;
    double f0,f1;

    // Note on mapping:  The centroid is always mapped to m_ppfSpotXXXX[m_nCenter][m_nCenter].  All other points are mapped to the "nearest" bin.

    // On the first pass, (nPass==0) we are just trying to determine the number of pixels (on average) in a spot.
    // Only spots that are within deviation are considered on the second pass.
    // On the second pass, we actually integrate these spots. The centroid is first found, and the spot is
    // normalized to fNormSpotIntensity.
    // On the third pass, we reject the worst fitting peaks.  Some peaks might be near the boarder of the image, or near
    // bad pixels (with bad background estimation).  These should be rejected.
    // The final pass integrates only non-rejected peaks.
    // In each pass, we consider each spot as we find it, and then mark all pixels in it with a value of "2".
    // This way, we consider each spot only once in a pass, and then reset values for the next pass.

    m_nDataPerPixel = nDataPerPixel;
    nSpotsUsed0 = 0;
    nSpotsUsed1 = 0;
    nSpotsUsed2 = 0;
    nSpotsUsed3 = 0;

    // oStat is used to compute statistics on the spot size.
    // After nPass==1 however, it can be used for other things.
    oStat.vClear();
    
    for (nPass=0;nPass<4;nPass++) {

        // Print label to user to tell what we are up to.
        if (nPass==0) 
            printf("Searching for peaks.\n");
        else if (nPass==1)
            printf("Integrating peaks (first pass).\n");
        else if (nPass==2)
            printf("Rejecting %5.1f%% of peaks that have worst fit.\n",100.0*fRejectPercent);
        else if (nPass==3)
            printf("Integrating peaks (second pass).\n");

        if (nPass==1) {
            // Discover average spot size.
            oStat.nMarkStat(false,Cstat::S_ABOVEBELOW,fSpotSizeSigmaReject);
            oStat.nMarkStat(false,Cstat::S_ABOVEBELOW,fSpotSizeSigmaReject);
            fSpotSizeSigma = oStat.fStandardDeviation();
            fSpotSizeAverage = oStat.fAverage();
            m_nSpotDim = (int) (sqrt((double)fSpotSizeAverage) + fSpotSizeSigma*fSpotSizeSigmaReject);
            m_nSpotDim*=nDataPerPixel;
            m_nCenter = m_nSpotDim/2;
            // Allocate spot memory.
            m_ppfSpotCount = new double*[m_nSpotDim];
            m_ppfSpotValue = new double*[m_nSpotDim];
            ppfTestCount = new double*[m_nSpotDim];
            ppfTestValue = new double*[m_nSpotDim];

            for (nx=0;nx<m_nSpotDim;nx++) {
                m_ppfSpotCount[nx] = new double[m_nSpotDim];
                m_ppfSpotValue[nx] = new double[m_nSpotDim];
                ppfTestCount[nx] = new double[m_nSpotDim];
                ppfTestValue[nx] = new double[m_nSpotDim];
            }
        }
        
        // We are integrating on the second and fourth passes.
        if ((nPass==1) || (nPass==3)) {
            for (nx=0;nx<m_nSpotDim;nx++) {
                for (ny=0;ny<m_nSpotDim;ny++) {
                    m_ppfSpotCount[nx][ny] = 0.0;
                    m_ppfSpotValue[nx][ny] = 0.0;
                    ppfTestCount[nx][ny] = 0.0;
                    ppfTestValue[nx][ny] = 0.0;
                }
            }
        }

        // On the third pass, we allocate memory for the deviation and reject arrays.
        if (nPass==2) {
            pnReject = new int[nSpotsUsed1];
            pfDev = new float[nSpotsUsed1];
        }

        // On the fourth pass, we use the deviations calculated on the 3rd pass to reject spots.
        if (nPass==3) {
            // Sort spots so that we can do rejections.
            g_pfCmpFloats = pfDev;
            int* pnSortIndex;
            pnSortIndex = new int[nSpotsUsed1];
            for (nx=0;nx<nSpotsUsed1;nx++) 
                pnSortIndex[nx] = nx;
            qsort(pnSortIndex,nSpotsUsed1,sizeof(int),float_cmp_rel);
            for (nx=0;nx<nSpotsUsed1;nx++)
                pnReject[nx] = 0;
            for (nx=0;nx<min(nSpotsUsed1*fRejectPercent,nSpotsUsed1);nx++) 
                pnReject[pnSortIndex[nSpotsUsed1 - 1 - nx]] = 1;
            delete[] pnSortIndex;
        }

        nSpotCount = 0;
        for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
            for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
                if (m_oBackground.rcGetBackground(npx[0],npx[1])==1) {
                    nStackPoint = 1;
                    nSpotStack0[0] = npx[0];
                    nSpotStack1[0] = npx[1];
                    m_oBackground.rcGetBackground(npx[0],npx[1]) = 2;
                    nPixelsInSpot = 0;
                    while (nStackPoint) {
                        int nPoint0,nPoint1;
                        nStackPoint--;
                        nPoint0 = nSpotStack0[nStackPoint];
                        nPoint1 = nSpotStack1[nStackPoint];
                        if (nPixelsInSpot+1>=nMaxInSpot) {
                            nPixelsInSpot = 0;
                            nStackPoint = 0;
                            break;
                        } else {
                            nSpot0[nPixelsInSpot] = nPoint0;
                            nSpot1[nPixelsInSpot] = nPoint1;
                            nPixelsInSpot++;
                        }
                        if (nStackPoint>=nMaxInSpot) {
                            nPixelsInSpot = 0;
                            nStackPoint = 0;
                            break;
                        }
                        for (nx=-1;nx<=1;nx++) {
                            for (ny=-1;ny<=1;ny++) {
                                if (m_oBackground.bInRange(nPoint0+nx,nPoint1+ny) && 
                                    (m_oBackground.rcGetBackground(nPoint0+nx,nPoint1+ny)==1)) {
                                    nSpotStack0[nStackPoint] = nPoint0+nx;
                                    nSpotStack1[nStackPoint] = nPoint1+ny;
                                    nStackPoint++;
                                    m_oBackground.rcGetBackground(nPoint0+nx,nPoint1+ny) = 2;
                                }
                            }
                        }
                    }
                    if (nPixelsInSpot) {

                        bSpotRejected = TRUE;

                        if (nPass==0) {
                            nSpotsUsed0++;
                            oStat.vAdd(nPixelsInSpot);
                        } else if ((nPass==1) || (nPass==2) || ((nPass==3)  && (!pnReject[nSpotCount]))) {
                            // We only consider spot who's dimensions are within the deviations required.
                            if ((nPixelsInSpot>=fSpotSizeAverage-fSpotSizeSigma*fSpotSizeSigmaReject) &&
                                (nPixelsInSpot<=fSpotSizeAverage+fSpotSizeSigma*fSpotSizeSigmaReject)) {
                                // Find the x and y centroids.
                                double fSumX;
                                double fSumY;
                                double fSumXV;
                                double fSumYV;
                                double fSum;
                                double fCentX;
                                double fCentY;
                                double fNorm;
                                double fAverageBackground;
                                int a2nPix[2];

                                bSpotRejected = FALSE;
                                fSumX = 0.0;
                                fSumY = 0.0;
                                fSumXV = 0.0;
                                fSumYV = 0.0;
                                fSum = 0.0;

                                fAverageBackground = 0.0;
                                for (nx=0;nx<nPixelsInSpot;nx++)
                                    fAverageBackground += m_oBackground.rfGetBackground(nSpot0[nx],nSpot1[nx]);
                                fAverageBackground/=nPixelsInSpot;

                                for (nx=0;nx<nPixelsInSpot;nx++) {
                                    f0 = (m_oCi.*m_oCi.prfGetPixel)(nSpot0[nx],nSpot1[nx]);
                                    f0 = max(0.0,f0 - fAverageBackground);
                                    fSumX += nSpot0[nx];
                                    fSumY += nSpot1[nx];
                                    fSumXV += nSpot0[nx]*f0;
                                    fSumYV += nSpot1[nx]*f0;
                                    fSum += f0;
                                }
                                if (fSum!=0.0) {

                                    // Clear the test buffer.
                                    for (nx=0;nx<m_nSpotDim;nx++) {
                                        for (ny=0;ny<m_nSpotDim;ny++) {
                                            ppfTestValue[nx][ny] = 0.0;
                                            ppfTestCount[nx][ny] = 0.0;
                                        }
                                    }

                                    fCentX = fSumXV/fSum;
                                    fCentY = fSumYV/fSum;
                                    fNorm = fNormSpotIntensity/fSum; 
                                    // For each pixel in the spot, we loop over a nDataPerPixel x nDataPerPixel grid.
                                    // These points are mapped onto the m_ppfSpotValue and m_ppfSpotCount arrays.
                                    for (nz=0;nz<nPixelsInSpot;nz++) {
                                        f0 = (m_oCi.*m_oCi.prfGetPixel)(nSpot0[nz],nSpot1[nz]);
                                        f0 = max(0.0,f0 - fAverageBackground);
                                        f0*=fNorm;
                                        for (a2nPix[0]=0;a2nPix[0]<nDataPerPixel;a2nPix[0]++) {
                                            for (a2nPix[1]=0;a2nPix[1]<nDataPerPixel;a2nPix[1]++) {
                                                nx = (int) ( m_nCenter + (a2nPix[0] - nDataPerPixel/2) + (nSpot0[nz] - fCentX)*nDataPerPixel);
                                                ny = (int) ( m_nCenter + (a2nPix[1] - nDataPerPixel/2) + (nSpot1[nz] - fCentY)*nDataPerPixel);
                                                if ((nx>=0) && (ny>=0) && (nx<m_nSpotDim) && (ny<m_nSpotDim)) {
                                                    ppfTestValue[nx][ny] += f0;
                                                    ppfTestCount[nx][ny] ++;
                                                }
                                            }
                                        }
                                    }

                                    // We could have nPass =  1,2 or 3 .  (These are the 2nd 3rd and 4th passes)
                                    if ((nPass==1) || (nPass==3)) {
                                        if (nPass==1) 
                                            nSpotsUsed1++;
                                        else
                                            nSpotsUsed3++;
                                        // Simply copy the imformation over.
                                        for (nx=0;nx<m_nSpotDim;nx++) {
                                            for (ny=0;ny<m_nSpotDim;ny++) {
                                                m_ppfSpotValue[nx][ny] += ppfTestValue[nx][ny];
                                                m_ppfSpotCount[nx][ny] += ppfTestCount[nx][ny];
                                            }
                                        }
                                    } else if (nPass==2) {
                                        // Compute deviations of spot from average.
                                        pfDev[nSpotsUsed2] = 0.0;
                                        for (nx=0;nx<m_nSpotDim;nx++) {
                                            for (ny=0;ny<m_nSpotDim;ny++) {
                                                if (m_ppfSpotCount[nx][ny]==0.0)
                                                    f0 = 0.0;
                                                else
                                                    f0 = m_ppfSpotValue[nx][ny]/m_ppfSpotCount[nx][ny];
                                                if (ppfTestCount[nx][ny]==0.0)
                                                    f1 = 0.0;
                                                else
                                                    f1 = ppfTestValue[nx][ny]/ppfTestCount[nx][ny];
                                                pfDev[nSpotsUsed2] += fabs(f0-f1);
                                            }
                                        }
                                        nSpotsUsed2++;
                                    }
                                }
                            }
                        }
                        nSpotCount++;

                        if ((nPass==3) && (bSpotRejected)) {
                            // Set all of the pixels in the spot to zero.  We can only do this on the last pass.  It is used to
                            // show the user which spots were utilized.
                            for (nx=0;nx<nPixelsInSpot;nx++) {
                                m_oBackground.rcGetBackground(nSpot0[nx],nSpot1[nx]) = 0;
                            }
                        }
                    }
                }
             }
        }
        // Change all pixels from 2 back to 1
        for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
            for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
                if (m_oBackground.rcGetBackground(npx[0],npx[1]) == 2)
                    m_oBackground.rcGetBackground(npx[0],npx[1]) = 1;
            }
        }
    }

    printf("\n");
    printf("%5d spots found in image.\n",nSpotsUsed0);
    printf("%5d spots used with sizes within %6.2f sigma of average spot size \n",nSpotsUsed1,fSpotSizeSigmaReject);
    printf("%5d spots used in final spot profile.\n",nSpotsUsed3);


    // Delete the temporary memory buffers.
    for (nx=0;nx<m_nSpotDim;nx++) {
        delete[] ppfTestCount[nx];
        delete[] ppfTestValue[nx];
    }
    delete[] ppfTestCount;
    delete[] ppfTestValue;
    delete[] pnReject;
    delete[] pfDev;


    return 0;
}

int Cbeamdata::nPrintAveragePeak(Cstring& sOutput) {
    FILE* pFOut;
    int nx,ny;
    double f0;
    bool bPointsWritten;

    pFOut = fopen(sOutput.string(),"w+t");
    if (!pFOut)
        return 1;

    printf("Writing average peak information to '%s'\n",sOutput.string());
    
    for (nx=0;nx<m_nSpotDim;nx++) {
        bPointsWritten = FALSE;
        for (ny=0;ny<m_nSpotDim;ny++) {
            if (m_ppfSpotCount[nx][ny]==0.0)
                f0 = 0.0;
            else
                f0 = m_ppfSpotValue[nx][ny]/m_ppfSpotCount[nx][ny];
            if (((fabs((double) nx-m_nCenter)/m_nDataPerPixel<2.0*m_fFW10M) && (fabs((double) ny-m_nCenter)/m_nDataPerPixel<2.0*m_fFW10M)) || (!m_fFW10M)) {
                fprintf(pFOut,"%lf %lf %lf\n",((double) nx-m_nCenter)/m_nDataPerPixel,((double) ny-m_nCenter)/m_nDataPerPixel,f0);
                bPointsWritten = TRUE;
            }
        }
        if (bPointsWritten)
            fprintf(pFOut,"\n");
    }
    fclose(pFOut);
    return 0;
}

int Cbeamdata::nPrintAveragePeakStats(Cstring& sOutputFileStub) {

    int nx,ny;
    double f0,f1;
    
    double fSumXV;
    double fSumYV;
    double fSum;
    double fCentX;
    double fCentY;

    const double fPixelStep = 0.5;
    double* pfProfileValues;    
    double* pfProfileCounts;
    double* pfIntProfile;
    int nProfileLength;

    Cstring sTemp;
    FILE* pFOut;

    // Find the centroid, 
    fSumXV = 0.0;
    fSumYV = 0.0;
    fSum = 0.0;
    for (nx=0;nx<m_nSpotDim;nx++) {
        for (ny=0;ny<m_nSpotDim;ny++) {
            if (m_ppfSpotCount[nx][ny]==0.0)
                f0 = 0.0;
            else
                f0 = m_ppfSpotValue[nx][ny]/m_ppfSpotCount[nx][ny];
            fSumXV += nx*f0;
            fSumYV += ny*f0;
            fSum += f0;
        }
    }
    fCentX = fSumXV/fSum;
    fCentY = fSumYV/fSum;

    // Allocate arrays.
    nProfileLength = (int) (m_nSpotDim/fPixelStep+1);
    pfProfileValues = new double[nProfileLength];
    pfProfileCounts = new double[nProfileLength];
    pfIntProfile = new double[nProfileLength];
    for (nx=0;nx<nProfileLength;nx++) {
        pfProfileValues[nx] = 0;
        pfProfileCounts[nx] = 0;
    }

    // Integrate the profile. (Maps peak profile into 1D plot w.r.t. distance to centroid.
    for (nx=0;nx<m_nSpotDim;nx++) {
        for (ny=0;ny<m_nSpotDim;ny++) {
            // Calculate distance to the centroid.
            if (m_ppfSpotCount[nx][ny]==0.0)
                f0 = 0.0;
            else
                f0 = m_ppfSpotValue[nx][ny]/m_ppfSpotCount[nx][ny];
            f1 = sqrt((double)((nx - fCentX)*(nx - fCentX) + (ny - fCentY)*(ny - fCentY)));
            if (((int) (f1/fPixelStep))<nProfileLength) {
                pfProfileValues[(int) (f1/fPixelStep)]+=f0;
                pfProfileCounts[(int) (f1/fPixelStep)]++;
            }
        }
    }

    // Integrate the profile.  This allows the user to get the peak width.
    f0 = 0.0;
    for (nx=0;nx<nProfileLength;nx++) {
        pfIntProfile[nx] = 0.0;
        if (pfProfileCounts[nx])
            f0 += pfProfileValues[nx]/pfProfileCounts[nx];
        pfIntProfile[nx] += f0;
    }
    for (nx=0;nx<nProfileLength;nx++) {
        pfIntProfile[nx]*=1000.0/f0;
    }

    // Figure out m_fFWHM and m_fFW10M
    f0 = 0.0;
    m_fFWHM = 0.0;
    m_fFW10M = 0.0;
    f1 = -1;
    for (nx=0;nx<nProfileLength;nx++) {
        if (pfProfileCounts[nx]!=0.0) {
            if (f1<0.0)
                f1 = pfProfileValues[nx]/pfProfileCounts[nx];
            if ((m_fFWHM==0.0) && (pfProfileValues[nx]/pfProfileCounts[nx]<0.5*f1))
                m_fFWHM = 2.0*nx*fPixelStep/m_nDataPerPixel;
            if ((m_fFW10M==0.0) && (pfProfileValues[nx]/pfProfileCounts[nx]<0.1*f1))
                m_fFW10M = 2.0*nx*fPixelStep/m_nDataPerPixel;
        }
    }

    // Output files.
    sTemp = sOutputFileStub+"_peak.txt";
    pFOut = fopen(sTemp.string(),"w+t");
    printf("Writing peak intensity as a function of pixel radius to '%s'\n",sTemp.string());
    for (nx=0;nx<nProfileLength;nx++) {
        if (pfProfileCounts[nx] != 0.0)
            f0 = pfProfileValues[nx]/pfProfileCounts[nx];
        else
            f0 = 0.0;
        fprintf(pFOut,"%lf %lf\n",nx*fPixelStep/m_nDataPerPixel,f0);
    }
    fclose(pFOut);

    sTemp = sOutputFileStub+"_int.txt";
    pFOut = fopen(sTemp.string(),"w+t");
    printf("Writing integrated peak intensity as a function of pixel radius to '%s'\n",sTemp.string());
    for (nx=0;nx<nProfileLength;nx++) {
        fprintf(pFOut,"%lf %lf\n",nx*fPixelStep/m_nDataPerPixel,pfIntProfile[nx]);
    }
    fclose(pFOut);

    printf("Width of peak at 50%% intensity = %6.2f pixels\n",m_fFWHM);
    printf("Width of peak at 10%% intensity = %6.2f pixels\n",m_fFW10M);

    delete[] pfProfileValues;
    delete[] pfProfileCounts;
    delete[] pfIntProfile;

    return 0;
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


int Cbeamdata::nFindBeamStopShadow(int nInitCenter0,int nInitCenter1) {
    int  a2x2nSearch[2][2]; // Ranges of consideration for beam stop finder.  
    int  a2nPos[2];
    int  nNewPoints;
    int  a2nFoundCenter[2];
    int nBeamStopRadius;    // Radius of beam stop.

    float fVal;
    int nx,ny;
    double f0,f1;
    Cstat oStat;

    const float fImagePercent = 0.5f;            // Use this percentage of the center part of the image.
    const float fBackgroundFlagSigma = 2.0f;     // How much below background must a pixel be before it is flagged as a beamstop pixel.
    const float fPercentRejectsInBeastop = 0.10f;// How many rejects in the "circle" drawn in the center.
    const float fPercentProfileRad = 0.1f;       // What percent of the image length should be discarced when determining background.

    int nImageRegion = min((int) (m_nDim[0]*fImagePercent/2),(int) (m_nDim[1]*fImagePercent/2));
    Cstring sTemp;

    float* pfBackgroundAverage;
    float* pfBackgroundSigma;


    /* Step 1)  Find profile of background as a function of distance from the nInitCenter.  This will be used in the
                next step were we find the masked region.  Pixels are placed in the region based on whether they are above or below
                this background.  Note that values below a certain cutoff are not included.
    */

    const int nInitProfileRad = (int) (m_nDim[0]*fPercentProfileRad);

    pfBackgroundAverage = new float[m_nDim[0]];
    pfBackgroundSigma = new float[m_nDim[0]];
    

    for (nx=nInitProfileRad;nx<m_nDim[0]/2-10;nx++) {
        double fRadStep = 0.5/(nx);
        double fRad;
        int nPix[2];

        oStat.vClear();
        // Step around a circle
        for (fRad=0.0;fRad<Gs_dPI*2.0;fRad+=fRadStep) {
            nPix[0] = (int) (nInitCenter0 + nx*cos(fRad));
            nPix[1] = (int) (nInitCenter1 + nx*sin(fRad));
            if (bInRange(nPix[0],nPix[1])) {
                fVal = (m_oCi.*m_oCi.prfGetPixel)(nPix[0],nPix[1]);
                oStat.vAdd(fVal);
            }
        }
        oStat.nMarkStat(false,Cstat::S_ABOVEBELOW,2.0);
        pfBackgroundSigma[nx] = oStat.fStandardDeviation();
        pfBackgroundAverage[nx] = oStat.fAverage();
    }
    
    for (nx=0;nx<nInitProfileRad;nx++) {
        pfBackgroundSigma[nx] = pfBackgroundSigma[nInitProfileRad];
        pfBackgroundAverage[nx] = pfBackgroundAverage[nInitProfileRad];

    }


    /* Step 2)  Find masked region.  We start with the beam-center as the branching out point for a region that
                will be "grown" iteratively. Pixels in the region are marked in the background map.  
    */

    m_oBackground.rcGetBackground(nInitCenter0,nInitCenter1)=2;

    // Set default ranges.
    a2x2nSearch[0][0] = nInitCenter0 - nImageRegion;
    a2x2nSearch[0][1] = nInitCenter0 + nImageRegion;
    a2x2nSearch[1][0] = nInitCenter1 - nImageRegion;
    a2x2nSearch[1][1] = nInitCenter1 + nImageRegion;

    do {

        // Add Points.
        nNewPoints = 0;
        for (a2nPos[0]=a2x2nSearch[0][0]+1;a2nPos[0]<=a2x2nSearch[0][1]-1;a2nPos[0]++) {
            for (a2nPos[1]=a2x2nSearch[1][0]+1;a2nPos[1]<=a2x2nSearch[1][1]-1;a2nPos[1]++) {
                if ((m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1])==2) ||
                    (m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1])==m_oBackground.m_nMaskValue)) {
                    // Compute distance to center.  This gives us our fAverage and fSigma values.
                    double fAverage;
                    double fSigma;
                    double fDist;
                    fDist = sqrt((double)( (a2nPos[0]-nInitCenter0)*(a2nPos[0]-nInitCenter0) + (a2nPos[1]-nInitCenter1)*(a2nPos[1]-nInitCenter1) ) );
                    fAverage = pfBackgroundAverage[min(m_nDim[0]-1,(int) fDist)];
                    fSigma = pfBackgroundSigma[min(m_nDim[0]-1,(int) fDist)];

                    for (nx=-1;nx<=1;nx++) {
                        for (ny=-1;ny<=1;ny++) {
                            if (m_oBackground.rcGetBackground(a2nPos[0]+nx,a2nPos[1]+ny)==0) {
                                    fVal = (m_oCi.*m_oCi.prfGetPixel)(a2nPos[0] + nx,a2nPos[1] + ny);
                                    if (fVal<fAverage-fSigma*fBackgroundFlagSigma) {
                                        m_oBackground.rcGetBackground(a2nPos[0]+nx,a2nPos[1]+ny)=3;
                                        nNewPoints++;
                                    }
                            }
                        }
                    }
                }
            }
        }
        // Change all 3 pixels to 2,
        // Change all 2 pixels to 1
        // Compute stats of pixels in oStat.

       
        oStat.vClear();
        for (a2nPos[0]=a2x2nSearch[0][0];a2nPos[0]<=a2x2nSearch[0][1];a2nPos[0]++) {
            for (a2nPos[1]=a2x2nSearch[1][0];a2nPos[1]<=a2x2nSearch[1][1];a2nPos[1]++) {
                if ((nx=m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1]))!=0) {
                    if (nx<=3) {
                        fVal = (m_oCi.*m_oCi.prfGetPixel)(a2nPos[0],a2nPos[1]);
                        oStat.vAdd(fVal);
                    }
                    if (nx==3)
                        m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1])=2;
                    else if (nx==2)
                        m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1])=1;

                }
            }
        }

    } while (nNewPoints);

    // m_oBackground.nPrintPeaks(sTemp = "test.img");

    /* Step 3)  We find the "center" of the probably circular beam stop and place the result in a2nFoundCenter[2].  The beam
                stop is probably circular.  To find this center, we compute the point that has the largest minimum distance
                to a 0.90 percent marked spot.  This is found by sending out "feelers" in 16 directions.
    */

    //m_oBackground.nFatten(1);

    
    a2nFoundCenter[0] = nInitCenter0;
    a2nFoundCenter[1] = nInitCenter1;
    

    float fMaxMinDist = 0.0;
    for (a2nPos[0]=a2x2nSearch[0][0]+1;a2nPos[0]<=a2x2nSearch[0][1]-1;a2nPos[0]++) {
        for (a2nPos[1]=a2x2nSearch[1][0]+1;a2nPos[1]<=a2x2nSearch[1][1]-1;a2nPos[1]++) {
            if (m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1])==1)
            {
                int a16x2nFeelers[16][2] = { {1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,1},{1,-1},{-1,-1},
                                             {2,1},{-2,1},{2,-1},{-2,-1},{1,2},{1,-2},{-1,2},{-1,-2} };
                
                float fMinDist = m_nDim[0];
                for (nx=0;nx<16;nx++) {
                    for (ny=0;ny<m_nDim[0]/4;ny++) {
                        int nTestPos0,nTestPos1;
                        nTestPos0 = a2nPos[0]+a16x2nFeelers[nx][0]*ny;
                        nTestPos1 = a2nPos[1]+a16x2nFeelers[nx][1]*ny;
                        if (!m_oBackground.rcGetBackground(nTestPos0,nTestPos1)) {
                            fMinDist = min(fMinDist,
					   (float)sqrt((double)(ny*ny*(a16x2nFeelers[nx][0]*a16x2nFeelers[nx][0] + a16x2nFeelers[nx][1]*a16x2nFeelers[nx][1]))));
                            break;
                        }
                    }
                }
                if (fMaxMinDist<fMinDist) {
                    fMaxMinDist = fMinDist;
                    a2nFoundCenter[0] = a2nPos[0];
                    a2nFoundCenter[1] = a2nPos[1];
                }
            }
        }
    }

    /* Step 4)  We find the "radius" of the center of the beam stop.  This will require expanding the circle until we get a dropoff
                below a specified number of pixels.    
    */

    
    for (ny=2;ny<m_nDim[0]/2;ny++) {
        int nFoundInBackground;
        int nFound;

        nFound = 0;
        nFoundInBackground = 0;
        // Compute all points in the beam stop center.
        for (a2nPos[0]=a2x2nSearch[0][0]+1;a2nPos[0]<=a2x2nSearch[0][1]-1;a2nPos[0]++) {
            for (a2nPos[1]=a2x2nSearch[1][0]+1;a2nPos[1]<=a2x2nSearch[1][1]-1;a2nPos[1]++) {
                f0 = sqrt((double)((a2nPos[0]-a2nFoundCenter[0])*(a2nPos[0]-a2nFoundCenter[0]) + (a2nPos[1]-a2nFoundCenter[1])*(a2nPos[1]-a2nFoundCenter[1])));
                if (f0<=ny) {
                    nFound++;
                    if (m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1])==1) {
                        nFoundInBackground++;
                    }
                }
            }
        }
        if ((((float) nFoundInBackground)/nFound)<1.0-fPercentRejectsInBeastop)
            break;
    }
    nBeamStopRadius = ny;
       

    /* Step 5)  Compute points along "rings".  We already have a bounding box specified in a2x2nSearch.
       These points will correspond to locations (from the right and from the left) where we intersect with
       the beam stop rod.  

    */

    int* a2x2pnPoints[2][2];
    bool* pbReject;

    int nRing;

    a2x2pnPoints[0][0] = new int[nImageRegion];
    a2x2pnPoints[0][1] = new int[nImageRegion];
    a2x2pnPoints[1][0] = new int[nImageRegion];
    a2x2pnPoints[1][1] = new int[nImageRegion];
    pbReject         = new bool[nImageRegion];
    

    for (nRing=nBeamStopRadius;nRing<nImageRegion;nRing++) {
        double fMaxRad;
        double fMinRad;
        double fRad;
        double fRadStep;
        int    nPix[2];

        fRadStep = 0.5/(nRing+1);
        fMaxRad = 0.0;
        fMinRad = Gs_dPI*2.0;

        // Step around a circle
        for (fRad=0.0;fRad<Gs_dPI*2.0;fRad+=fRadStep) {
            nPix[0] = (int) (a2nFoundCenter[0] + nRing*cos(fRad));
            nPix[1] = (int) (a2nFoundCenter[1] + nRing*sin(fRad));
            if (m_oBackground.rcGetBackground(nPix[0],nPix[1])==1) {
                fVal = (m_oCi.*m_oCi.prfGetPixel)(nPix[0],nPix[1]);
                if (fVal!=0.0) {
                    if (fMaxRad<fRad) 
                        fMaxRad = fRad;
                    if (fMinRad>fRad) 
                        fMinRad = fRad;
                }
            }
        }
        // We might (inconveniently) have the beam mask located on the theta=0 axis that we have defined.
        // If so, just swap the two so that the "min" and "max" have the same meaning.
        if (fMaxRad-fMinRad>Gs_dPI) {
            f0 = fMaxRad-Gs_dPI*2.0;
            fMaxRad = fMinRad;
            fMinRad = f0;
        }
        a2x2pnPoints[0][0][nRing] = (int) (a2nFoundCenter[0] + nRing*cos(fMinRad));
        a2x2pnPoints[0][1][nRing] = (int) (a2nFoundCenter[1] + nRing*sin(fMinRad));
        a2x2pnPoints[1][0][nRing] = (int) (a2nFoundCenter[0] + nRing*cos(fMaxRad));
        a2x2pnPoints[1][1][nRing] = (int) (a2nFoundCenter[1] + nRing*sin(fMaxRad));
    }


    /*  Step 6)  Fit the points found in these lines.  We are finding two lines.  The lines might
        be described as (0 as a function of 1) or (1 as a function of 0).
        Also, we only want to use the points furtherest away from the beam stop.

    */

    int nIndepVar;       // 0 if [0] is a function of [1]. 1 if [1] is a function of [0].
    int nLine;
    int nPass;
    double a2x2fLineParams[2][2];    // Parameters for lines.
    double a2fAverageDir[2];          // For points found on the beam stop holder, what is the average vector to the found center.
    int nRingLinearStart;
    int nRingLinearEnd;

    nRingLinearStart = nBeamStopRadius+10;
    nRingLinearEnd = (int) (nImageRegion*0.9);

    // The variable with the larger standard deviation is used.
    oStat.vClear();
    for (nRing=nRingLinearStart;nRing<nRingLinearEnd;nRing++)
        oStat.vAdd(a2x2pnPoints[0][0][nRing]);
    f0 = oStat.fStandardDeviation();
    oStat.vClear();
    for (nRing=nRingLinearStart;nRing<nRingLinearEnd;nRing++)
        oStat.vAdd(a2x2pnPoints[0][1][nRing]);
    f1 = oStat.fStandardDeviation();
    if (f0>f1)
        nIndepVar = 0;
    else
        nIndepVar = 1;


    a2fAverageDir[0] = 0;
    a2fAverageDir[1] = 0;

    double fSumXY = 0.0;
    double fSumX2 = 0.0;
    double fSumX = 0.0;
    double fSumY = 0.0;
    double fSumN = 0.0;


    for (nLine=0;nLine<2;nLine++) {
        for (nPass=0;nPass<2;nPass++) {
            
            
            
            if (nPass==0) {
                // On the first pass, we use all points in our linear least squares.
                for (nRing=nRingLinearStart;nRing<nRingLinearEnd;nRing++)
                    pbReject[nRing] = FALSE;
            } else {
                // On the second pass, we reject all points that did not fit very well.
                oStat.vClear();
                for (nRing=nRingLinearStart;nRing<nRingLinearEnd;nRing++) {
                    f0 = fabs(a2x2pnPoints[nLine][!nIndepVar][nRing] - ( a2x2fLineParams[nLine][1]  + a2x2fLineParams[nLine][0]*a2x2pnPoints[nLine][nIndepVar][nRing]));
                    oStat.vAdd(f0);
                }
                oStat.nMarkStat(FALSE,oStat.S_ABOVEBELOW,3.0);
                for (nRing=nRingLinearStart;nRing<nRingLinearEnd;nRing++) {
                    if (!oStat.bIsMarked(nRing-nRingLinearStart)) {
                        pbReject[nRing] = 1;
                    }
                }
            }
            fSumXY = 0.0;
            fSumX2 = 0.0;
            fSumX = 0.0;
            fSumY = 0.0;
            fSumN = 0.0;
            
            for (nRing=nRingLinearStart;nRing<nRingLinearEnd;nRing++) {
                if (!pbReject[nRing]) {
                    fSumXY += a2x2pnPoints[nLine][0][nRing]*a2x2pnPoints[nLine][1][nRing];
                    fSumX2 += a2x2pnPoints[nLine][nIndepVar][nRing]*a2x2pnPoints[nLine][nIndepVar][nRing];
                    fSumX += a2x2pnPoints[nLine][nIndepVar][nRing];
                    fSumY += a2x2pnPoints[nLine][!nIndepVar][nRing];
                    fSumN++;
                    a2fAverageDir[0] += (a2x2pnPoints[nLine][0][nRing] - a2nFoundCenter[0]);
                    a2fAverageDir[1] += (a2x2pnPoints[nLine][1][nRing] - a2nFoundCenter[1]);
                }
            }
            a2x2fLineParams[nLine][0] = (fSumXY*fSumN - fSumX*fSumY)/(fSumN*fSumX2 - fSumX*fSumX);
            a2x2fLineParams[nLine][1] = (-fSumXY*fSumX + fSumX2*fSumY)/(fSumN*fSumX2 - fSumX*fSumX);

        }
    }



    a2fAverageDir[0]/=fSumN;
    a2fAverageDir[1]/=fSumN;

    m_oBackground.nClearAll();


    /* Step 7)  Fill in the lines.  We want to fill in all pixels in between the lines down to the circle center.
                Also, fill in the beam stop center.
    */

    for (a2nPos[0]=0;a2nPos[0]<m_nDim[0];a2nPos[0]++) {
        for (a2nPos[1]=0;a2nPos[1]<m_nDim[1];a2nPos[1]++) {
            // Is the point on the correct side of the circle center.
            // Use nLine = 0 for this test
            if ((((a2nPos[nIndepVar] - a2nFoundCenter[nIndepVar])*a2fAverageDir[nIndepVar])) +
            (((a2x2fLineParams[0][1]  + a2x2fLineParams[0][0]*a2nPos[nIndepVar]) - a2nFoundCenter[!nIndepVar])*a2fAverageDir[!nIndepVar])>=0.0) {
                // Is the point between the two lines?
                // We should have one positive and one negative deviation.
                f0 = (a2x2fLineParams[0][1]  + a2x2fLineParams[0][0]*a2nPos[nIndepVar]) - a2nPos[!nIndepVar];
                f1 = (a2x2fLineParams[1][1]  + a2x2fLineParams[1][0]*a2nPos[nIndepVar]) - a2nPos[!nIndepVar];
                if (((f0>=0.0) && (f1<=0.0)) || ((f0<=0.0) && (f1>=0.0))) {
                    m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1]) = 1;
                }
            }
            f0 = sqrt((double)((a2nPos[0]-a2nFoundCenter[0])*(a2nPos[0]-a2nFoundCenter[0]) + (a2nPos[1]-a2nFoundCenter[1])*(a2nPos[1]-a2nFoundCenter[1])));
            if (f0<nBeamStopRadius) {
                m_oBackground.rcGetBackground(a2nPos[0],a2nPos[1]) = 1;
            }
        }
    }      

    delete[] a2x2pnPoints[0][0];
    delete[] a2x2pnPoints[0][1];
    delete[] a2x2pnPoints[1][0];
    delete[] a2x2pnPoints[1][1];
    delete[] pbReject;
    delete[] pfBackgroundAverage;
    delete[] pfBackgroundSigma;



    return 0;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

int Cbeamdata::nInitDetector(double f2Theta,double fDistance,double fBeam0,double fBeam1)
{
    Cstring sTemp;
    Cstring m_sDetectorName;
    float m_a3fS0[3];

    m_sDetectorName = "";    
    
    m_oCi.m_oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1, &m_sDetectorName);
    
    if( m_poDetector ) 
        delete m_poDetector;
    m_poDetector = new Cdetector(m_oCi.m_oHeader, m_sDetectorName);
    if (!m_poDetector->bIsAvailable())
        return 1;

    if (fDistance>0.0)
        m_poDetector->m_poGoniometer->nSetDistance(fDistance);
    if (f2Theta>-999.0)
        m_poDetector->m_poGoniometer->nSetSwing(f2Theta);
    if (fBeam0>-9999.0)
        m_poDetector->m_poSpatial->nSetBeamPosition(fBeam0,fBeam1);

    if (m_poSource)
        delete m_poSource;
    m_poSource = new Csource(m_oCi.m_oHeader);
    if (!m_poSource->bIsAvailable())
        return 1;

    m_poSource->vCalcGetS0(m_a3fS0);
    m_fWavelength = m_poSource->fGetWavelength();
    vMulVec3DScalar(m_a3fS0,1.0/m_fWavelength,m_a3fS0Wave);
    m_poDetector->vUpdateDDDN();

    return 0;

}

double Cbeamdata::f2dtothetaToReso(int nPix0,int nPix1) {

    double fReso;
    double fchi;

    if (!m_poDetector)
        return 1;
    fReso = m_poDetector->fCalcGetResolution(nPix0,nPix1,m_a3fS0Wave, NULL,&fchi);

    return fReso;

    
    // Convert resolution into 2theta degrees.
    //f2theta = (asin(min(1.0,m_fWavelength/2.0/fReso)))*2.0;
    //fchi = m_fChiSign*(fchi + m_fChiOffset);
    //if (fchi<0.0)
    //    fchi += 2.0*Gs_dPI;
}


int Cbeamdata::nFitParams(double fReso,double f02Theta,double fRad2Theta,double f0Dist,double fRadDist,double fRadBeam) 
{
    double      fResoMin = 0.0;
    double      fResoMax = 0.0;
    double      f0Beam0  = 0.0;
    double      f0Beam1  = 0.0;

    double      f0 = 0.0;
    double      f1 = 0.0;
    
    float       ff0 = 0.0;
    float       ff1 = 0.0;
    
    double      fSum = 0.0;
    
    int         nx = 0;
    int         ny = 0;


    int             nSimplexPass = 0;
    const int       nMaxSimplexPass = 20;
    double          a4fSteps[4];
    double          a4fDelta[4];
    double          a4fBest[4];
    double          fBestSum = 0.0;
    
    itr<double>     afValues;
    double*         pfValues;
    
    const int       nStep = 6;
    
    int             nDim0 = 0;
    int             nDim1 = 0;
    
    int             n0 = 0;
    int             n1 = 0;

    if( nInitDetector() )
        return 1;

    if (f0Dist <= 0.0) {
        f0Dist = m_poDetector->fGetDistance();
        fRadDist = ABS(f0Dist*fRadDist)/2.0;
    }    
    m_poDetector->m_poSpatial->nGetBeamPosition(&ff0,&ff1);
    f0Beam0 = ff0;
    f0Beam1 = ff1;

    fResoMin = fReso*1.01;
    fResoMax = fReso*0.99;

    printf("Searching by simplex method for detector parameters, resolution = %.3lf\n",fReso);
    printf("Steping 2theta   from %6.2lf to %6.2lf\n",f02Theta - fRad2Theta,f02Theta + fRad2Theta);
    printf("Steping distance from %6.2lf to %6.2lf\n",f0Dist - fRadDist,f0Dist + fRadDist);
    printf("Percent Complete 00---10---20---30---40---50---60---70---80---90---100\n");
    printf("                 +");


    a4fSteps[0] = fRadBeam;
    a4fSteps[1] = fRadBeam;
    a4fSteps[2] = fRad2Theta;
    a4fSteps[3] = fRadDist;

    // Bin the values in the image.
    afValues.setsize((m_nDim[0]/nStep+1)*(m_nDim[1]/nStep+1));
    pfValues = &afValues[0];
    nDim0 = 0;
    nDim1 = 0;
    for (n0=0,nDim0=0;n0<m_nDim[0];n0+=nStep,nDim0++) {
        for (n1=0,nDim1=0;n1<m_nDim[1];n1+=nStep,nDim1++) {
            fSum = 0.0;
            for (nx = 0; nx < nStep; nx++) {
                for (ny = 0; ny < nStep; ny++) {
                    if ((nx+n0 < m_nDim[0]) && (ny+n1 < m_nDim[1])) {
                        fSum += m_oCi.fGetPixel(nx + n0,ny + n1);
                    }                    
                }                
            }
            *pfValues = fSum;
            pfValues++;
        }
    }

    pfValues = &afValues[0];   


    for (nSimplexPass = 0; nSimplexPass < nMaxSimplexPass; nSimplexPass++) {
        nSimplexRun(0.0,a4fDelta,a4fSteps);


        do {

            nInitDetector(a4fDelta[2] + f02Theta,a4fDelta[3] + f0Dist,a4fDelta[0] + f0Beam0,a4fDelta[1] + f0Beam1);

            fSum = 0.0;

            for (n0=0;n0<nDim0;n0++) {
                for (n1=0;n1<nDim1;n1++) {
                    f0 = f2dtothetaToReso(n0*nStep + nStep/2,n1*nStep + nStep/2);
                    if ((f0 >= fResoMax) && (f0 <= fResoMin)) 
                        fSum += pfValues[nDim0*n1 + n0];
                }
            }
            fSum *= -1.0;
                
        } while (!nSimplexRun(fSum,a4fDelta,NULL));

        fSum *= -1.0;
        
        if( fSum > fBestSum )
        {
            for (nx=0; nx< 4; nx++)
                a4fBest[nx] = a4fDelta[nx];
            fBestSum = fSum;
        }
    }
       

    return 0;
}


#define nSIMPLEXPARAMS 4

#define nCopyNewVector(nBest,fValue) {\
    vSubVecNDVecND(nSIMPLEXPARAMS,&m_afParamPointSum[0],&m_aafParamPointVectors[nBest][0],&m_afParamPointSum[0]);\
    vAddVecNDVecND(nSIMPLEXPARAMS,&m_afParamPointSum[0],&afDelta[0],&m_afParamPointSum[0]);\
    vCopyVecND(nSIMPLEXPARAMS,&afDelta[0],&m_aafParamPointVectors[nBest][0]);\
    m_afParamPointValues[nBest] = fValue; \
    if (nVerbose)\
        printf("Point %2d is better. (%.2lf)\n",nBest,fValue);\
    };\
    bPointWorked  = true;


int Cbeamdata::nSimplexRun(double fValue,double* afDelta,double* afStep) {

    int nVerbose = 0;
    
    static double m_aafParamPointVectors[nSIMPLEXPARAMS+1][nSIMPLEXPARAMS];
    static double m_afParamPointValues[nSIMPLEXPARAMS+1];
    static double m_afParamPointSum[nSIMPLEXPARAMS];
    static double m_afParamGoodPoint[nSIMPLEXPARAMS];
    static int m_nSimplexState1;
    static int m_nSimplexState2;
    static int m_nSimplexPoints;

    if (afStep) {
        m_nSimplexState1 = 0;
        m_nSimplexState2 = -1;
        // Count the number of free variables available.
        m_nSimplexPoints = 4;
    }




    int nx,ny;
    int nBest;
    int nWorst,nSecondWorst;
    double fFac,fFac1,fFac2;
    bool bPointWorked;

    if (m_nSimplexState1==0) {
        if (m_nSimplexState2==-1) {
            // Initialize the points in the simplex.
            vZeroMat(nSIMPLEXPARAMS,1,&m_aafParamPointVectors[0][0]);
            for (nx=0;nx<= nSIMPLEXPARAMS; nx++) {
                vZeroMat(nSIMPLEXPARAMS,1,&m_aafParamPointVectors[nx][0]);
                for (ny=0;ny< nSIMPLEXPARAMS; ny++) 
                    m_aafParamPointVectors[nx][ny]=((rand() % 2)*2 - 1)*afStep[ny];
            }
            m_nSimplexState2 = 0;
            if (nVerbose)
                printf("Initializing.\n");
        } 
        if (nVerbose)
            printf("Initializing Vertex %d.\n",m_nSimplexState2);


        m_afParamPointValues[m_nSimplexState2] = fValue;
        m_nSimplexState2++;

        
        if (m_nSimplexState2==m_nSimplexPoints) {
            // Initialize for the computation steps.
            m_nSimplexState1 = 1;
            m_nSimplexState2 = 0;
            // Compute the sum of the points.
            vZeroMat(nSIMPLEXPARAMS,1,&m_afParamPointSum[0]);
            for (nx=0;nx<m_nSimplexPoints;nx++) {
                vAddVecNDVecND(nSIMPLEXPARAMS,&m_aafParamPointVectors[nx][0],&m_afParamPointSum[0],&m_afParamPointSum[0]);
            }
        } else
            vCopyVecND(nSIMPLEXPARAMS,&m_aafParamPointVectors[m_nSimplexState2][0],&afDelta[0]);
    }
    while (m_nSimplexState1==1) {
        // Calculate Worst, SecondWorst and Best
        if (m_afParamPointValues[0]>m_afParamPointValues[1]) {
            nWorst = 0;
            nSecondWorst = 1;
        } else {
            nWorst = 1;
            nSecondWorst = 0;
        }
        nBest = 0;
        for (nx=1;nx<m_nSimplexPoints;nx++) {
            if (m_afParamPointValues[nx]<m_afParamPointValues[nBest])
                nBest = nx;
            if (m_afParamPointValues[nx]>m_afParamPointValues[nWorst]) {
                nSecondWorst = nWorst;
                nWorst = nx;
            }
        }
        bPointWorked = false;

        // Make sure that all the points are not bad.  If so, then we will restore the
        // zero point solution and terminate.
        for (nx=0;nx<m_nSimplexPoints;nx++) {
            if (m_afParamPointValues[nx]<1e10)
                break;
        }
        if (nx==m_nSimplexPoints) {
            vZeroMat(nSIMPLEXPARAMS,1,&afDelta[0]);
            return 2;
        }

        // Check for termination conditions.
        if (ABS(m_afParamPointValues[nWorst] - m_afParamPointValues[nBest])/
            (ABS(m_afParamPointValues[nWorst]) + ABS(m_afParamPointValues[nBest]))<0.001) {
            vCopyVecND(nSIMPLEXPARAMS,&m_aafParamPointVectors[0][0],&afDelta[0]);
            return 1;
        }
        
        // Start of simplex algorithm.
        
        
        fFac = 0.0;
        switch (m_nSimplexState2) {
        case 0:
            // Try a factor -1.0
            fFac = -1.0;
            m_nSimplexState2 = 1;
            
            break;
        case 1:
            // Just tried -1.0.  
            if (fValue<m_afParamPointValues[nWorst]) {
                nCopyNewVector(nWorst,fValue);
                
            }
            if (fValue<m_afParamPointValues[nBest]) {
                // Did Better than the current best.
                // try -2.0.
                fFac = -2.0;
                m_nSimplexState2 = 2;
            } else if (fValue>m_afParamPointValues[nSecondWorst]) {
                // Did Worse than the second worst point.
                fFac = 0.5;
                m_nSimplexState2 = 3;
            } else {
                // Avoid an infinite loop.
                if (!bPointWorked)
                    return 2;
                m_nSimplexState2 = 0;
                fFac = 0.0;
            }
            
            break;
        case 2:
            // Just tried -2.0.
            if (fValue<m_afParamPointValues[nWorst]) {
                nCopyNewVector(nWorst,fValue);
            }
            m_nSimplexState2 = 0;
            fFac = 0.0;
        case 3:
            // Just tried 0.5
            if (fValue<m_afParamPointValues[nWorst]) {
                nCopyNewVector(nWorst,fValue);
                m_nSimplexState2 = 0;
                fFac = 0.0;
            } else {
                // Can't get rid of the high point.  
                // Do a contraction.
                
                if (nVerbose)
                    printf("Contracting.\n");
                for (nx=0;nx<m_nSimplexPoints;nx++) {
                    if (nx!=nBest) {
                        for (ny=0;ny<nSIMPLEXPARAMS;ny++) {
                            m_aafParamPointVectors[nx][ny] = 
                                0.5*(m_aafParamPointVectors[nx][ny] + m_aafParamPointVectors[nBest][ny]);
                        }
                    }
                }
                // Restart, since we need to recompute values for most vectors.
                m_nSimplexState1 = 0;
                m_nSimplexState2 = 0;
                fFac = 0.0;
                // Copy over first solution.
                vCopyVecND(nSIMPLEXPARAMS,&m_aafParamPointVectors[0][0],&afDelta[0]);
            }
        }

        if (fFac == 0.0)
            continue;
        else {
            if (nVerbose)
                printf("Trying Factor=%.1lf\n",fFac);

            // Compute the new test point best on fFac.
            fFac1 = (1.0 - fFac)/(m_nSimplexPoints-1);
            fFac2 = fFac1 - fFac;
            for (nx=0;nx<nSIMPLEXPARAMS;nx++) {
                afDelta[nx] = m_afParamPointSum[nx]*fFac1 - m_aafParamPointVectors[nWorst][nx]*fFac2;
            }
            return 0;
        }
    }
    return 0;
}


int Cbeamdata::nFindResoRings(double d2ThetaMin_DEG,
                              double dMaxExpected2ThetaRingWidth_DEG,
                              int nImageSizeCompressionFactor,
                              itr<double>& adExpectedResoRings_A,
                              Cimage_header* poHeader)
{
    // RB 7/13/05  Here we used to clean the ring information in the header first, but now
    // we shouldn't do that, because this function may be called repeatedly on several images,
    // so we want to keep in the header all (unique) rings from all images.
    
    itr<double>         af2Theta;               // The array for integrated 2theta values.
    itr<double>         afRejectRangesOut;      // The array for rejected ranges.
    itr<double>         afRejectResoMidOut;
    itr<double>         afRejectRangesIntensity;

    /////////////////////////////////////////////////////////
    // Determine the resolution range to search over.
	double      f0 = 0.0;
    double      f1 = 0.0;
    
    bUpdateInternalHeader(poHeader);

    CExperimentalStrategy       oStrat(m_oCi.m_oHeader);
	
    if( !oStrat.bIsAvailable() )
        return 1;
    
    oStrat.vSetResoMethod(oStrat.eResoChiCoverage,50.0);
	oStrat.nCalcReso(f0,f1);
    /////////////////////////////////////////////////////////

    nInitDetector();
    
    double              fStep = 0.1;
    if( 0 != m_oBackground.nIntegrate(5, af2Theta, fStep, m_poDetector, m_poSource, m_nCompression) )
        return 1;
    
    //Cstring sTestTemp = "test.img";
    //m_oBackground.nPrintPeaks((Cstring&)sTestTemp);
    
    int         nBins = min(af2Theta.size(),((int) (oStrat.fConvertToDegrees(f0)/fStep)));
    
    double              fMin2ThetaWidth = 1.0;
    int         nBinsToCheck = (int) (fMin2ThetaWidth/(fStep*2.0));
    
    double      fWavelength = m_poSource->fGetWavelength();  
    
    
    int         nTitlePrinted = 0;
    bool        bTableHeaderPrinted = false;
    //printf("Resolution rings found.\n");
    //printf("--------------------------------------------------\n");
    //printf(" Peak 2Theta 2Theta 2Theta    Reso    Reso    Peak \n");
    //printf("   #     Mid  Start    End   Start     End Percent \n");
    //printf("--------------------------------------------------\n");
    

    // Next, we detect the peaks in 2Theta,
	int         nx = 0;
    int         ny = 0;
    int         nIntValue = 0;
    double      fStart = 0.0;
    double      fEnd = 0.0;
    double      fStartAng = 0.0;
    double      fEndAng = 0.0;
    double      fMidAng = 0.0;
    double      fValidAng = 0.0;
    double      fStartDeg = 0.0;
    double      fEndDeg = 0.0;
    double      fMidDeg = 0.0;
    double      fValueStart = 0.0;
    double      fValueEnd = 0.0;
    double      fIntValue = 0.0;
    double      fBaseValue = 0.0;
    double      fPercentPeak = 0.0;
    
    bool        bFail = true;
    int         nBin = 0;

    for(nBin = max(nBinsToCheck+1,(int) (d2ThetaMin_DEG/fStep)); nBin < nBins - nBinsToCheck - 1;nBin++) 
    {
        if( (af2Theta[nBin+1] > 0.0) && (af2Theta[nBin-1] > 0.0) && 
            (af2Theta[nBin-1] <= af2Theta[nBin]) && 
            (af2Theta[nBin+1] <= af2Theta[nBin]) )     // rb: a test whether af2Theta[nBin] is a local (3-bin) maximum 
        {

			// Try to fit this peice to a quadratic curve.
			itr<double>         afFitPoints;
			double*             apfFuncs[3] = {ms_pfConst, ms_pfX, ms_pfX2};
			double              afQuad[3];
            
            fIntValue = 0.0;
            nIntValue = 0;     //rb: nIntValue is basically the number of bins, used for a peak fit: 2*nBinsToCheck+1
			
            for (int nBinShift = -nBinsToCheck; nBinShift <= nBinsToCheck; nBinShift++)
            {
				afFitPoints + af2Theta[nBinShift + nBin];
                fIntValue += af2Theta[nBinShift + nBin];
                nIntValue++;
			}
            
            fValueStart = af2Theta[nBin - nBinsToCheck];
            fValueEnd = af2Theta[nBin + nBinsToCheck];
			
            if( afFitPoints.size() > 0 &&
                !nSolveLeastSquares(3, afFitPoints.size(), &afFitPoints[0], NULL, apfFuncs, afQuad, NULL, NULL) ) 
            {
				// We must have a concave down fit, and a maximum somewhere inbetween
				bFail = (afQuad[2] >= 0.0);
				if (!bFail)
                {
					f0 = - afQuad[1] / ( 2.0*afQuad[2]);
					bFail =  (f0 <= 0.3) || (f0 >= 0.7);
				}
			} 
            else
				bFail = true;

		
            if (!bFail)
            {
                // Also need to check the "percent peak" condition.
                fBaseValue = 0.5 * (fValueStart + fValueEnd) * nIntValue;
                fPercentPeak = 100.0 * (fIntValue - fBaseValue) / max(1.0, fBaseValue);
                
                // We use this point.
                fMidDeg = fStep*(0.5 + nBin);
                fStartDeg = fStart = fStep*(0.5 + nBin) - dMaxExpected2ThetaRingWidth_DEG/2.0;
                fEndDeg = fEnd = fStep*(0.5 + nBin) + dMaxExpected2ThetaRingWidth_DEG/2.0;
				
                fMidAng = fWavelength/(max(0.00001,sin(Gs_dRADIANS_PER_DEGREE*fMidDeg/2.0)*2));
                
                fStartAng = fWavelength/(max(0.00001,sin(Gs_dRADIANS_PER_DEGREE*fStart/2.0)*2)); 
                
                fEndAng = fWavelength/(max(0.00001,sin(Gs_dRADIANS_PER_DEGREE*fEnd/2.0)*2)); 
			}

            // Make sure we have a strong enough peak.
            if (!bFail)
            {
                if (fPercentPeak < 0.2)
                    bFail = true;
            }

			if ((!bFail) && (adExpectedResoRings_A.size()))
            {
				// Check to see that the d-spacing is within acceptable ranges.
				fValidAng = -1.0;
				for (nx = 0; nx < adExpectedResoRings_A.size(); nx++)
                {
					if ((adExpectedResoRings_A[nx]>=fEndAng) && (adExpectedResoRings_A[nx]<=fStartAng))
                    {
						if ((fValidAng == -1.0) || (ABS(adExpectedResoRings_A[nx] - fMidAng) < ABS(fValidAng - fMidAng)))
                        {
							fValidAng = adExpectedResoRings_A[nx];
						}
					}
				}
				
                if (fValidAng == -1.0)
					bFail = true;
			}

			if (!bFail)
            {
				// Add exclude regions.
				for (nx = 0; nx + 1 < afRejectRangesOut.size(); nx+=2)
                {

					// This rejects if there is *any* overlap.  This is probably not what we want.
					//if (ABS(fStartAng - fEndAng) + ABS(afRejectRangesOut[nx] - afRejectRangesOut[nx+1]) >
					//	ABS(min(fEndAng,afRejectRangesOut[nx+1]) - max(fStartAng,afRejectRangesOut[nx])))
					//	bFail = true;
					if ((fMidAng >= afRejectRangesOut[nx+1]) && (fMidAng <= afRejectRangesOut[nx]))
						bFail = true;
				}
				
                if (!bFail)
                {
					afRejectRangesOut + fStartAng;
					afRejectRangesOut + fEndAng;
                    afRejectRangesIntensity + fPercentPeak;
                    afRejectResoMidOut + fMidAng;
				}
			}
				
			if( !bFail ) 
            {
                fStart = fStartAng;
                fEnd = fEndAng;

                if( !bTableHeaderPrinted )
                {
                    printf("Resolution rings found.\n");
                    printf("--------------------------------------------------\n");
                    printf(" Peak 2Theta 2Theta 2Theta    Reso    Reso    Peak \n");
                    printf("   #     Mid  Start    End   Start     End Percent \n");
                    printf("--------------------------------------------------\n");
                    
                    bTableHeaderPrinted = true;
                }
                
                printf(" %4d %6.2lf %6.2lf %6.2lf %7.4lf %7.4lf %7.2lf\n",
                       nTitlePrinted, fMidDeg, fStartDeg, fEndDeg, fStartAng, fEndAng, fPercentPeak);
                nTitlePrinted++;
            }
        }
    }

    if (!nTitlePrinted) 
    {
        printf("INFO:  No powder rings detected.\n");
    }
    printf("--------------------------------------------------\n");

    // DEBUG ONLY
    if( false ) 
    {
        FILE*           pFOut = fopen("ResoRingPlot.txt","wt");
        
        for (nBin = 0; nBin < nBins; nBin++)
        {
            fprintf(pFOut, "%.2lf %.2lf\n", fStep * (0.5 + nBin), af2Theta[nBin]);
        }
        
        fclose(pFOut);
    }
    
    // Update header information.
    if( poHeader )
    {
        bWriteResoRingsToHeader(afRejectResoMidOut, 
                                afRejectRangesIntensity, 
                                dMaxExpected2ThetaRingWidth_DEG,
                                poHeader);
    }

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cbeamdata::bWriteResoRingsToHeader(itr<double>& afResoRingMid, 
                                        itr<double>& afResoRingStrength, 
                                        double dMaxExpected2ThetaRingWidth_DEG,
                                        Cimage_header* poHeader)
{
    ////////////////////////////////////////////////////////
    // Establish the keyword name
    Cstring     strCallingModuleName = sDtrekGetModuleName();
    strCallingModuleName.upcase();
    strCallingModuleName += "_RANK";

    if( afResoRingMid.size() != afResoRingStrength.size() )
        return false; // safety

    if( 0 == afResoRingMid.size() )
        return true; // nothing to do

    /////////////////////////////////////////////
    // Read rings already in the header
    std::vector<double>     adExistingRings;

    int     nx = -1;
    double  f0 = 0.0;

    while( 0 == poHeader->nGetValueDictionary(strCallingModuleName, "RESO_RING_MID_RESO_*", f0, ++nx) )
    {
        adExistingRings.push_back(f0);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Loop through newly found rings and see if they match the rings already in the header 
    int     iRingIndexInHeader = adExistingRings.size() - 1;
    bool    bRingMatch = false;
    for(int ii=0; ii < afResoRingMid.size(); ii++)
    {
        bRingMatch = false;
        for(int jj=0; jj < adExistingRings.size(); jj++)
        {
            if( fabs(afResoRingMid[ii] - adExistingRings[jj]) < dMaxExpected2ThetaRingWidth_DEG / 2. )
            {
                bRingMatch = true;
                break;
            }
        }
        
        if( !bRingMatch )
        {
            iRingIndexInHeader++;
            
            poHeader->nReplaceValueDictionary(strCallingModuleName, "RESO_RING_MID_RESO_*", afResoRingMid[ii], iRingIndexInHeader);
            poHeader->nReplaceValueDictionary(strCallingModuleName, "RESO_RING_PERCENT_*", afResoRingStrength[ii], iRingIndexInHeader);
        }
    }
    
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RB I am adding this function as a patch to the following problem, m_oCi.m_oHeader originally comes from the header of the image, not
// from the header which is an input file. The input header file may have the correct beam position, wavelength etc., which
// must be passed to m_oCi.m_oHeader
bool Cbeamdata::bUpdateInternalHeader(Cimage_header* poHeader)
{
    if( !poHeader )
        return false; // just a safety check;

    if( !poHeader->bIsAvailable() )
        return false; // ditto

    Cstring              strDetectorPrefix(""); 
    if( 0 != poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    //////////////////////////////////////////////////////////////////////////////////////
    Cspatial            oSpatial(*poHeader, strDetectorPrefix);
    if( !oSpatial.bIsAvailable() )
        return false;

    // RB For now we should not adjust the spatial to the compressed image, because that breaks Thad's code
    // for finding rings. We need to investigate further why it does.

    // One thing I have noiticed is that in nCompress() function he has the following statement:
    // m_oCi.m_oHeader = oOld.m_oHeader;

    // So even though the image is compressed, he wants to use the original image header with the old pixels size etc.
    // So it looks like we should not call bChangeBinning() after all.
    
    //oSpatial.bChangeBinning(m_nCompression, m_nCompression);
    
    oSpatial.nUpdateHeader(&(m_oCi.m_oHeader), strDetectorPrefix);
    //////////////////////////////////////////////////////////////
    
    Cgoniometer         oDetGonio(*poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
        return false;
    
    oDetGonio.nUpdateHeader(&(m_oCi.m_oHeader), strDetectorPrefix);
    ///////////////////////////////////////////////////////////////
    
    Cstring        strSourcePrefix("SOURCE_");
    Cwavelength    oWavelength(*poHeader, strSourcePrefix, false);   // 'false' means no alpha1-alpha2 check

    oWavelength.nUpdateHeader(&(m_oCi.m_oHeader), strSourcePrefix);

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   

