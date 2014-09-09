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
// CAxialPhoto.cc     Thaddeus J. Niemeyer           17-Apr-2002
//    Computes an axial photograph from images.
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

#include "CAxialPhoto.h"
#include "CBackground.h"

void CAxialPhoto::vInit() {
    m_nStat = 0;
    m_anSeqNum.clear();
    m_fSigma = 2.0;
    m_nAxisCode = 100;
    m_a2fReso[0] = 999999.99;
    m_a2fReso[1] = 0.1;
    m_nDim0 = 0;
    m_nDim1 = 0;
    m_sOutput = "dtaxial.img";
    m_bPixelAdjust = false;
    m_a2fRotRange[0] = 0.0;
    m_a2fRotRange[1] = 360.0;
    m_nHKLRestrict = -1;
    m_nHKLValue = 0;
    m_fHKLWidth = 1.0;
    

    m_poBackground = NULL;
    m_poPhotons = NULL;
    m_poAxialPhotons = NULL;   
    m_poInput = NULL;

};

void CAxialPhoto::vSetResolution(double fResoHigh,double fResoLow) {
    if (fResoHigh>fResoLow)
        std::swap(fResoHigh,fResoLow);
    m_a2fReso[0] = fResoLow;
    m_a2fReso[1] = fResoHigh;
};
void CAxialPhoto::vAddSeq(int nSeqStart,int nSeqEnd) {
    if (nSeqStart<= nSeqEnd) {
        m_anSeqNum + nSeqStart;
        m_anSeqNum + nSeqEnd;
    };
};
void CAxialPhoto::vSetSigma(double fSigma) {
    m_fSigma = fSigma;
};
void CAxialPhoto::vSetAxisCode(int nAxisCode) {
    m_nAxisCode = nAxisCode;
};
void CAxialPhoto::vSetRotRange(double fRotStart,double fRotEnd) {
    m_a2fRotRange[0] = fRotStart;
    m_a2fRotRange[1] = fRotEnd;
};

void CAxialPhoto::vSetZone(int nHKL,int nValue) {
    m_nHKLRestrict = nHKL;
    m_nHKLValue = nValue;
};
void CAxialPhoto::vSetZoneWidth(double fHKLWidth) {
    m_fHKLWidth = fHKLWidth;
};


CAxialPhoto::CAxialPhoto(Cimage_header& oHeader):
m_oHeader(oHeader),
m_oCrystal(oHeader),
m_oSource(oHeader),
m_oGonioIn(oHeader, Ccrystal::ms_sCrystalPrefix),
m_oScan(oHeader)
 { 
    Cstring sTemp;

    vInit();
    sTemp = "";
    (void) oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1, &sTemp);
    m_poDetector = new Cdetector(oHeader,sTemp);

    if ((!m_oScan.bIsAvailable()) || 
        (!m_oHeader.bIsAvailable()) || 
        (!m_oCrystal.bIsAvailable()) ||
        (!m_oSource.bIsAvailable()) ||
        (!m_oGonioIn.bIsAvailable()) ||
        (!m_poDetector->bIsAvailable()))
        m_nStat = 1;
    else
        m_nStat = 0;
};

CAxialPhoto::~CAxialPhoto() {
    if (m_poDetector)
        delete m_poDetector;
    m_poDetector = NULL;
    if (m_poBackground)
        delete m_poBackground;
    m_poBackground = NULL;
    if (m_poPhotons)
        delete m_poPhotons;
    m_poPhotons = NULL;
    if (m_poAxialPhotons)
        delete m_poAxialPhotons;
    m_poAxialPhotons = NULL;
    if (m_poInput)
        delete m_poInput;
    m_poInput = NULL;  
};


int CAxialPhoto::nSetAxialGonio()
{
    double a3fHKL[3];
    double a3fRecip[3];
    double a3fScanVec[3];
    double a3x3fUBMat[3][3];
    double a3x3fT[3][3];
    double a3x3fTP[3][3];
    double a3x3fRealOrient[3][3];

    int nx,ny;
    float a3fTemp[3];
    double a3x3fTemp[3][3];
    


    // Build the HKL value that we will project down.
    for (nx=0,ny=1;nx<3;nx++,ny*=10)
        a3fHKL[2-nx] = (m_nAxisCode/ny) % 10;
    
    // Determine the reciprocal space vector for that point.    
    m_oCrystal.nCalcOrientMatrix();
    m_oCrystal.vGetOrientMatrix(&a3x3fUBMat[0][0]);
    
    // Build 'real' cell axes as collumn vectors.
    // We use the fact that trans[UB] * inv[trans[UB]] = I
    // So, each collumn in inv[trans[UB]] is perpendicular to two collumns in UB.
    vCopyMat3D(&a3x3fUBMat[0][0],&a3x3fTemp[0][0]);
    vTranMat3D(a3x3fTemp);
    fInvMat3D(&a3x3fTemp[0][0],&a3x3fRealOrient[0][0]);
    vMulMat3DVec3D(a3x3fRealOrient,a3fHKL,a3fRecip);

    // Normalize the reciprocal vector,and build two other vectors perpendicular to it.
    // Call this matrix T
    fNormVec3D(a3fRecip);
    vBuildBasis3D(a3fRecip,a3x3fT);    
    vNormMat3D(a3x3fT);

    // Get the rotation vector associated with the scan.  Build two other vectores perpendicular to it.
    // Call this matrix T'
    m_oScan.m_poRotation->vGetVector(&a3fTemp[0]); vCopyVec3D(&a3fTemp[0],&a3fScanVec[0]);
    vBuildBasis3D(a3fScanVec,a3x3fTP);
    vNormMat3D(a3x3fTP);
    
    // We want to find a goniometer matrix G to map all three orthonormal vectors T to T'
    // Thus, we want GT=T'
    vTranMat3D(a3x3fT);
    vMulMat3DMat3D(a3x3fTP,a3x3fT,m_a3x3fAxialGonio);
    vCopyVec3D(a3x3fTP[0],m_a3fAxialRotVec);
    return 0;
};

int CAxialPhoto::nAddBackground(Cimage& oImage) {
    int n0,n1;
    double f0,f1;
    bool bInit;

    
    if (!m_poBackground) {
        m_poBackground = new Cimage(oImage);
        bInit = true;
    } else
        bInit = false;
    for (n0 = 0; n0 < m_nDim0; n0++) {
        for (n1 = 0; n1 < m_nDim1; n1++) {
            f0 = (oImage.*oImage.prfGetPixel)(n0,n1);

            if (bInit)
                f1 = 1e20;
            else
                f1 = (m_poBackground->*m_poBackground->prfGetPixel)(n0,n1);
            (m_poBackground->*m_poBackground->prnSetPixel)(n0,n1,min(f0,f1));

        };
    };    


    return 0;
};

int CAxialPhoto::nFindPeaks(Cimage& oImage) {
    Cbackground oBackground(oImage);
    int n0,n1;
    double f0;
    
    oBackground.nBuildPeakPixels(30);
    oBackground.nFatten(2);


    // If we don't have the photon image yet, construct it.

    if (!m_poPhotons) {
        // Build a photon image.  This should use floating point numbers.
        m_poPhotons = new Cimage(m_nDim0,m_nDim1,eImage_realIEEE);
    };

    for (n0 = 0; n0 < m_nDim0; n0++) {
        for (n1 = 0; n1 < m_nDim1; n1++) {
            (m_poPhotons->*m_poPhotons->prnSetPixel)(n0,n1,0.0);
        };
    };


    // Extract all pixels that seem to have significant contribution.
    // Remove background before adding them to m_poPhotons.
    for (n0 = 0; n0 < m_nDim0; n0++) {
        for (n1 = 0; n1 < m_nDim1; n1++) {  
            if (oBackground.rcGetBackground(n0,n1)) { 
                f0 = (oImage.*oImage.prfGetPixel)(n0,n1) - (m_poBackground->*m_poBackground->prfGetPixel)(n0,n1);
                if (f0 >0.0) {
                    f0 += (m_poPhotons->*m_poPhotons->prfGetPixel)(n0,n1);
                    (m_poPhotons->*m_poPhotons->prnSetPixel)(n0,n1,f0);
                };
            };         
        };
    };


    return 0;
};

int CAxialPhoto::nBuildFinalImage() {
    int n0,n1;
    double fValue;
    
    for (n0 = 0; n0 < m_nDim0; n0++) {
        for (n1 = 0; n1 < m_nDim1; n1++) {
            // Use m_poInput to do this.
            fValue = (m_poAxialPhotons->*m_poAxialPhotons->prfGetPixel)(n0,n1);
            fValue += (m_poBackground->*m_poBackground->prfGetPixel)(n0,n1);            
            (m_poInput->*m_poInput->prnSetPixel)(n0,n1,fValue);
        };
    };
    m_poInput->nWrite(m_sOutput);


    return 0;
};

int CAxialPhoto::nAddPeaks(double fRot)
{ 
    int nx=0;
    double f0=0.0;    
    double a3fTempVec[3];
    double a3x3fTempMat[3][3];
    float fPx0,fPx1;
    float f1MM,f2MM,f3MM;
    float a3x3fFloatBuf[3][3];
    float a3fFloatBuf[3];
    int nStat;

    if (!m_poAxialPhotons) {
        int n0,n1;
        // Build a photon image.  This should use floating point numbers.
        m_poAxialPhotons = new Cimage(m_nDim0,m_nDim1,eImage_realIEEE);
        for (n0 = 0; n0 < m_nDim0; n0++) {
            for (n1 = 0; n1 < m_nDim1; n1++) {
                (m_poAxialPhotons->*m_poAxialPhotons->prnSetPixel)(n0,n1,0.0);
            };
        };
    };
    
    int n0,n1;
    int n0Axial,n1Axial;
    int n0Frac,n1Frac;
    int nPix;
    int nRotAngle;
    double f0Axial,f1Axial;
    double f1Px,f2Px;    
    double a3fS0[3];
    double a3x3fDD[3][3];
    double a3x3fDDInv[3][3];
    double a3fDN[3];
    double fValue,fPrevValue;
    double fFraction;


    double a3x3fGonioMatIn[3][3];   
    double a3x3fInvGonioMatIn[3][3];
    double a3x3fInvUBMat[3][3];
    double a3fVDCOIn[3];
    double a3fTSIn[3];    
    double a3fXRIn[3];
    double a3fXIn[3];
    double a3fHKLIn[3];
    double a3fRotVecIn[3];
    double a3x3fInvRotIn[3][3];
    double a3fE3CrossS0In[3];
    double fLorentzIn;

    double a3fXAxial[3];
    double a3fXRotAxial[3];
    double a3fVDCOAxial[3];
    double a3fVDCOScaleAxial[3];
    double a3fE3CrossS0Axial[3];
    double fLorentzAxial;


    m_oGonioIn.vCalcGetRotMatrix(&a3x3fFloatBuf[0][0],m_oGonioIn.nGetNum(m_oScan.m_poRotation->sGetName()));vCopyMat3D(&a3x3fFloatBuf[0][0],&a3x3fGonioMatIn[0][0]);
    fInvMat3D(&a3x3fGonioMatIn[0][0],&a3x3fInvGonioMatIn[0][0]);
    m_oScan.m_poRotation->vGetVector(&a3fFloatBuf[0]);  vCopyVec3D(&a3fFloatBuf[0],&a3fRotVecIn[0]);    
    m_poDetector->nCalcGetDDDN(&a3x3fDD[0][0], &a3fDN[0]);
    m_poDetector->vUpdateDDDN();
    fInvMat3D(&a3x3fDD[0][0], &a3x3fDDInv[0][0]);
    m_oSource.vCalcGetS0(a3fS0);
    m_oCrystal.nCalcOrientMatrix();
    m_oCrystal.vGetOrientMatrix(&a3x3fTempMat[0][0]);
    fInvMat3D(&a3x3fTempMat[0][0],&a3x3fInvUBMat[0][0]);
    
    vConvRotVec3DMat3D(-fRot, a3fRotVecIn, &a3x3fInvRotIn[0][0]);
    vCross3D(a3fRotVecIn, a3fS0, a3fE3CrossS0In);

    vCross3D(m_a3fAxialRotVec,a3fS0,a3fE3CrossS0Axial);

    double fMassCalc0;
    double fMassCalc1;
    double fMassDenom;
    itr<int> anPixArray0;
    itr<int> anPixArray1;
    itr<int> anInvestigatePixArray0;
    itr<int> anInvestigatePixArray1;
    itr<double> afPixValue;
    
    bool bScatteredVecIntersectsDetector = false;


    for (n0 = 0; n0 < m_nDim0; n0++) {
        for (n1 = 0; n1 < m_nDim1; n1++) {
            fValue = (m_poPhotons->*m_poPhotons->prfGetPixel)(n0,n1);
            if (fValue>0.0) {
                -anPixArray0 + n0;
                -anPixArray1 + n1;
                -afPixValue + fValue;
                fMassCalc0 = n0*fValue;
                fMassCalc1 = n1*fValue;
                fMassDenom = fValue;
                (m_poPhotons->*m_poPhotons->prnSetPixel)(n0,n1,0.0);

                if (m_bPixelAdjust) {
                    anInvestigatePixArray0.push(n0);
                    anInvestigatePixArray1.push(n1);
                    // Find all pixels that are adjacent having this value.
                    while (!anInvestigatePixArray0.isempty()) {
                        int nPixLook0,nPixLook1;
                        int nPixDisp0,nPixDisp1;
                        nPixLook0 = anInvestigatePixArray0.pop();
                        nPixLook1 = anInvestigatePixArray1.pop();
                        for (nPixDisp0=-1;nPixDisp0<=1;nPixDisp0++) {
                            for (nPixDisp1=-1;nPixDisp1<=1;nPixDisp1++) {
                                if ((nPixLook0+nPixDisp0>=0) && (nPixLook0+nPixDisp0<m_nDim0) && 
                                    (nPixLook1+nPixDisp1>=0) && (nPixLook1+nPixDisp1<m_nDim1) &&
                                    ((nPixDisp0==0) || (nPixDisp1==0))) {
                                    fValue = (m_poPhotons->*m_poPhotons->prfGetPixel)(nPixLook0 + nPixDisp0,nPixLook1 + nPixDisp1);
                                    if (fValue>0) {
                                        anInvestigatePixArray0.push(nPixLook0 + nPixDisp0);
                                        anInvestigatePixArray1.push(nPixLook1 + nPixDisp1);
                                        fMassCalc0 += (nPixLook0 + nPixDisp0)*fValue;
                                        fMassCalc1 += (nPixLook1 + nPixDisp1)*fValue;
                                        fMassDenom += fValue;
                                        anPixArray0 + (nPixLook0 + nPixDisp0);
                                        anPixArray1 + (nPixLook1 + nPixDisp1);
                                        afPixValue + fValue;
                                        (m_poPhotons->*m_poPhotons->prnSetPixel)(nPixLook0 + nPixDisp0,nPixLook1 + nPixDisp1,0.0);
                                    };
                                };
                            };
                        };
                    };                        
                };
                // Calculate the center of mass for the pixels we obtained.
                // This will give us the single pixel coordinate to translate.


                fMassCalc0/=fMassDenom;
                fMassCalc1/=fMassDenom;
                                

                // Compute the unrotated X reciprocal vector.
                f1Px = fMassCalc0;
                f2Px = fMassCalc1;
                nStat = m_poDetector->m_poSpatial->nPixeltoMM(f1Px, f2Px,
                    &f1MM, &f2MM, &f3MM);
                a3fVDCOIn[0] = f1MM;
                a3fVDCOIn[1] = f2MM;
                a3fVDCOIn[2] = f3MM;
                vMulMat3DVec3D(a3x3fDD, a3fVDCOIn, a3fTSIn);
                (void) fNormVec3D(a3fTSIn);
                vAddVec3DVec3D(a3fS0, a3fTSIn, a3fXRIn);
                fLorentzIn = 1.0f / max(ABS(fDot3D(a3fE3CrossS0In, a3fXRIn)), 1e-10);
                vMulMat3DVec3D(a3x3fInvRotIn, a3fXRIn, a3fTSIn);      // fTS is a temp var here
                vMulMat3DVec3D(a3x3fInvGonioMatIn, a3fTSIn, a3fXIn); // fTS is a temp var here

                // Did the user specify a zone restriction?
                if (m_nHKLRestrict != -1) {
                    // Multiply by inverse goniometer matrix.
                    vMulMat3DVec3D(a3x3fInvUBMat,a3fXIn,a3fHKLIn);
                    if ((a3fHKLIn[m_nHKLRestrict]>m_nHKLValue + m_fHKLWidth) ||
                        (a3fHKLIn[m_nHKLRestrict]<m_nHKLValue - m_fHKLWidth))
                        continue;
                };
               

                // We have a3fXIn.  This completes the calculations which involve the input pixel coordinates.
                
                // Next, we apply our calculated goniometer matrix to get the new goniometer rotated X vector
                vMulMat3DVec3D(m_a3x3fAxialGonio,a3fXIn,a3fXAxial);
                
                
                // We want to discover the amount of rotation necc. to bring XRC
                // back onto the Ewald sphere.
                // a3fXRO contains the 'ideal' Ewald sphere position (calculated from (x,y) points on plate)
                // a3fXRC contains the 'calculated' Ewald sphere position (probably NOT on the sphere).
                
                double a3fDStarV1[3];
                double a3fDStarV2[3];
                double a3fDStarV3[3];
                double fk1,fk2,fb;
                double a2fSolutions[2];
                
                // Find the three vectors v1,v2 and v3.
                // Project the offset d* vector of the spot down into the rotation plane to get v1
                vMulVec3DScalar(m_a3fAxialRotVec,fDot3D(m_a3fAxialRotVec,a3fXAxial),a3fTempVec);
                vSubVec3DVec3D(a3fXAxial,a3fTempVec,a3fDStarV1);
                // Take cross product of m_a3fRotVector and v1 to get v2.
                // This gives v1 cross v2 == m_a3fRotVector
                vCross3D(m_a3fAxialRotVec,a3fDStarV1,a3fDStarV2);
                // Calculate v3
                vSubVec3DVec3D(a3fXAxial,a3fDStarV1,a3fDStarV3);
                
                // Calculate constants a (int a3fTemp1), k1,k2 and b
                
                vSubVec3DVec3D(a3fDStarV3,a3fS0,a3fTempVec);
                fb = fDot3D(a3fDStarV1,a3fDStarV1) + fDot3D(a3fTempVec,a3fTempVec) - 1.0;
                fk1 = 2.0* fDot3D(a3fDStarV1,a3fTempVec);
                fk2 = 2.0* fDot3D(a3fDStarV2,a3fTempVec);
                
                // Solve quadratic system.  Choose solution with the lowest absolute value.
                nx = nSolveQuadratic(fk1*fk1 + fk2*fk2, fk2*fb,fb*fb - fk1*fk1,a2fSolutions);
                if (nx) {
                    int nAngleSolution;
                    int nSolution;
                    double fRotAngle;

                    for (nSolution=0;nSolution<nx;nSolution++)
                    {
                        if ((a2fSolutions[nSolution] >= -1.0) && (a2fSolutions[nSolution]<=1.0))
                        {
                            // Since a2fSolutions[] contains the sin of the rot angle, we need to look at both solutions.
                            for (nAngleSolution = 0; nAngleSolution < 2;nAngleSolution++) 
                            {
                                fRotAngle = asin(a2fSolutions[nSolution])/fRADIANS_PER_DEGREE;
                                if (nAngleSolution == 1)
                                    fRotAngle = 180.0 - fRotAngle;
                                
                                // Build a rotation matrix to take a3fXAxial to it's new position.
                                vConvRotVec3DMat3D(fRotAngle, m_a3fAxialRotVec, &a3x3fTempMat[0][0]);
                                vMulMat3DVec3D(a3x3fTempMat,a3fXAxial,a3fXRotAxial);
                                
                                // Now we can continue with the processing to project this pixel onto the image.
                                vSubVec3DVec3D(&a3fXRotAxial[0],&a3fS0[0],&a3fVDCOAxial[0]);
                                f0 = fLenVec3D(a3fVDCOAxial);
                                
                                // Verify that
                                if (ABS(f0 - 1.0) < 0.0001) 
                                {
                                    if ((fRotAngle>=0.0) && (fRotAngle<360))
                                    {
                                        nRotAngle = (int) fRotAngle;
                                    }
                                    else if ((fRotAngle<0.0) && (fRotAngle+360.0>=0.0))
                                    {
                                        nRotAngle = (int) (fRotAngle + 360.0);
                                    }
                                    else if ((fRotAngle>360) && (fRotAngle-360.0<360.0))
                                    {
                                        nRotAngle = (int) (fRotAngle - 360.0);
                                    }
                                    
                                    m_a360fDegCounts[nRotAngle] += anPixArray0.size();
                                    m_a360fCumulDegCounts[nRotAngle] += anPixArray0.size();


                                    if ((nRotAngle>=m_a2fRotRange[0]) && (nRotAngle<=m_a2fRotRange[1]))
                                    {
                                    } 
                                    else if ((nRotAngle - 360.0>=m_a2fRotRange[0]) && (nRotAngle - 360.0 <=m_a2fRotRange[1]))
                                    {
                                    } 
                                    else if ((nRotAngle + 360.0>=m_a2fRotRange[0]) && (nRotAngle + 360.0 <=m_a2fRotRange[1]))
                                    {
                                    } 
                                    else
                                        continue;


                                    fLorentzAxial = 1.0f / max(ABS(fDot3D(a3fE3CrossS0Axial, a3fXRotAxial)), 1e-10);

                                    // We have a solution.  
                                    // Transform it to detector pixel coordinates.
                                    bScatteredVecIntersectsDetector = false;
                                    if( 0 == m_poDetector->nScaleSOnPlate(a3fVDCOAxial, a3fVDCOScaleAxial) )                                    
                                    {
                                        vMulMat3DVec3D(a3x3fDDInv, a3fVDCOScaleAxial, a3fVDCOAxial);
                                        
                                        if( 0 == m_poDetector->m_poSpatial->nMMtoPixel(a3fVDCOAxial[0],
                                                                                       a3fVDCOAxial[1],
                                                                                       a3fVDCOAxial[2],
                                                                                       &fPx0, 
                                                                                       &fPx1) )
                                            bScatteredVecIntersectsDetector = true;
                                    }
                                    
                                    if( bScatteredVecIntersectsDetector )
                                    {
                                        for (nPix=0;nPix<anPixArray0.size();nPix++)
                                        {
                                            fValue = (afPixValue[nPix] - (m_poBackground->*m_poBackground->prfGetPixel)(anPixArray0[nPix],anPixArray1[nPix]))/fLorentzIn;
                                            if (fValue>0.0)
                                            {
                                                f0Axial = (anPixArray0[nPix] +  fPx0 - fMassCalc0);
                                                f1Axial = (anPixArray1[nPix] +  fPx1 - fMassCalc1);
                                                n0Axial = (int) f0Axial;
                                                n1Axial = (int) f1Axial;
                                                
                                                
                                                if ((n0Axial>=0) && (n1Axial>=0) && (n0Axial + 1<m_nDim0) && (n1Axial + 1<m_nDim1)) 
                                                {
                                                    for (n0Frac = 0; n0Frac < 2; n0Frac++)
                                                    {
                                                        for (n1Frac =0; n1Frac < 2; n1Frac++)
                                                        {
                                                            fFraction = ABS(n0Axial + !n0Frac - f0Axial)*ABS(n1Axial + !n1Frac - f1Axial);
                                                            fPrevValue = (m_poAxialPhotons->*m_poAxialPhotons->prfGetPixel)((int) n0Axial + n0Frac,(int) n1Axial + n1Frac);
                                                            (m_poAxialPhotons->*m_poAxialPhotons->prnSetPixel)((int) n0Axial + n0Frac,(int) n1Axial + n1Frac,fPrevValue + fValue*fFraction);
                                                            
                                                        };
                                                    };
                                                };
                                            };
                                        };
                                    };
                                    // Don't look at the other solution, even if it's available.
                                    break;
                                };
                            };
                        };
                    };
                };             
            };
        };
    };
    return 0;
};

void CAxialPhoto::vInitDegCounts(double* pfDegCounts) {
    int nx;
    for (nx=0;nx<360;nx++)
        pfDegCounts[nx] = 0.0;
};
void CAxialPhoto::vPrintDegCounts(double* pfDegCounts)
{
    int nx,ny;
    int nWidth;
    double fCounts;
    double fMaxStart;
    double fMaxCounts;
    const double afWidths[] = { 5,10,20,40,0.0};

    printf("------------------------------------------------\n");
        
    for (nWidth = 0; afWidths[nWidth] != 0.0; nWidth++) {
        fMaxStart = 0.0;
        fMaxCounts = 0.0;
        for (nx=0;nx<360;nx++) {
            fCounts = 0.0;
            for (ny = 0; ny< (int) afWidths[nWidth];ny++) 
                fCounts += afWidths[(nx + ny) % 360];
            if (fCounts>fMaxCounts) {
                fMaxCounts = fCounts;
                fMaxStart = nx;
            };
        };
        printf("For coverage with width of %.2lf restrict rotation in interval [%.2lf %.2lf]\n",
            afWidths[nWidth],fMaxStart,fMaxStart + afWidths[nWidth]);
    };    
    printf("------------------------------------------------\n");
};

int CAxialPhoto::nRun()
{
    int nImageSeqPair;
    int nImageSeqNum;
    int nImageCount;
    double fRotStart,fRotEnd;
    bool   bPrintPerImageWidthInfo = false;

    if (!m_anSeqNum.size())
        return 1;

    nSetAxialGonio();

    nImageCount = 0;

    for (nImageSeqPair = 0; nImageSeqPair < m_anSeqNum.size(); nImageSeqPair+=2) {
        for (nImageSeqNum = m_anSeqNum[nImageSeqPair]; nImageSeqNum <= m_anSeqNum[nImageSeqPair+1]; nImageSeqNum++) {
            if (!m_poInput) 
                m_poInput = new Cimage();

            m_oScan.vSetSeqNum(nImageSeqNum);
            if (m_oScan.nGetImage(m_poInput)) 
                return 1;
            m_poInput->nGetDimensions(&m_nDim0, &m_nDim1);            

            nAddBackground(*m_poInput);

            
        };
    };
    m_poBackground->nWrite("dtaxial_background.img");
        

    vInitDegCounts(&m_a360fCumulDegCounts[0]);
    for (nImageSeqPair = 0; nImageSeqPair < m_anSeqNum.size(); nImageSeqPair+=2)
    {
        for (nImageSeqNum = m_anSeqNum[nImageSeqPair]; nImageSeqNum <= m_anSeqNum[nImageSeqPair+1]; nImageSeqNum++) 
        {
            if (!m_poInput) 
                m_poInput = new Cimage();

            m_oScan.vSetSeqNum(nImageSeqNum);
            m_oScan.nGetImage(m_poInput);
            
            nFindPeaks(*m_poInput);
            fRotStart = m_oScan.fCalcGetRotStart(nImageSeqNum);
            fRotEnd   = m_oScan.fCalcGetRotEnd(nImageSeqNum);

            vInitDegCounts(&m_a360fDegCounts[0]);
            nAddPeaks(0.5*(fRotStart + fRotEnd));
            
            if( bPrintPerImageWidthInfo )
            {
                printf("Width Info for sequence #%d\n",nImageSeqNum);
                vPrintDegCounts(&m_a360fDegCounts[0]);
            }

            nImageCount++;
        }
    }
    printf("Cumulative Width Info\n");
    vPrintDegCounts(&m_a360fCumulDegCounts[0]);

    nBuildFinalImage();


    return 0;
}
