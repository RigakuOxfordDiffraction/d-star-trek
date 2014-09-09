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


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//////////////////////  CSpoint /////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

CSpoint::CSpoint() {
    int nx;
    nAdjCt=0;
    nRow=0;
    nCol=0;
    nContrib = 0;
};

CSpoint::~CSpoint() { 
};

int CSpoint::nAddAdj(int nTo) {
    int nx;
    for (nx=0;nx<nAdjCt;nx++) {
        if (anAdj[nx]==nTo)
            return 0;
    };
    if (nAdjCt==NEIGHBORVECS) {
        cout << "CSpoint array dimension exceeded.\n" << flush;
        EXIT(0);
    };
    anAdj[nAdjCt++] = nTo;
    return 1;
};

int CSpoint::nRemoveAdj(int nTo) {
    int nx,ny;
    for (nx=0,ny=0;nx<nAdjCt;) {
        if (anAdj[nx]==nTo) 
            nx++;
        else
            anAdj[ny++]=anAdj[nx++];
    };
    nAdjCt=ny;
    return 0;
};
bool CSpoint::bIsAdj(int nTo) {
    int nx;
    for (nx=0;nx<nAdjCt;nx++) {
        if (anAdj[nx]==nTo)
            return TRUE;
    };
    return FALSE;   
};

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//////////////////////  CS0point ///////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

int CS0point::nPrint() {
    printf(" [%d,%d] Rot(%5f) Contrib(%d)\n",nOrient,nPoint,fRot,nContrib);
    return 0;
};



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////  Cabsorbsurface ///////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


Cabsorbsurface::Cabsorbsurface() {
    int nx;

    m_poSPoints = NULL;
    m_poS0Points = NULL;
    m_pnC = NULL;
    for (nx=0;nx<3;nx++) {
        m_a3pfSCoeffs[nx] =  NULL;   
        m_a3pnSPoints[nx] =  NULL;
        m_a3pnTriangles[nx] = NULL;
    };
    for (nx=0;nx<2;nx++) {
        m_a2pfS0Coeffs[nx] = NULL;
        m_a2pnS0Points[nx] = NULL;
    };
    m_a2pnClosestS[0] = NULL;
    m_a2pnClosestS[1] = NULL;
    m_pnClosestS0 = NULL;
    m_pnOrientation = NULL;
    m_pfOrientAngle = NULL;
    m_pfSVectors = NULL;
    m_pfChi = NULL;

    m_nNumRefs = 0;
    m_pnRingStart = NULL;
    m_pnRingEnd = NULL;
    m_nRings = 0;

    m_nDegreesPerS0 = 15;
    m_nMax2Theta = 90;
    m_nNumSPoints = 0;
    m_fMinRingContribAverage = 60;
    m_fMinContrib = 50;    
    m_fMinNonHarmonicRings = 0.5;
};

Cabsorbsurface::~Cabsorbsurface() {
    int nx;
    if (m_poSPoints)
        delete[] m_poSPoints;
    if (m_poS0Points)
        delete[] m_poS0Points;
    if (m_pnC)
        delete[] m_pnC;
    for (nx=0;nx<3;nx++) {
        delete[] m_a3pfSCoeffs[nx];
        delete[] m_a3pnSPoints[nx];
        delete[] m_a3pnTriangles[nx];
    };
    for (nx=0;nx<2;nx++) {
        delete[] m_a2pfS0Coeffs[nx];
        delete[] m_a2pnS0Points[nx];
    };
    delete[] m_a2pnClosestS[0];
    delete[] m_a2pnClosestS[1];
    delete[] m_pnClosestS0;
    delete[] m_pnOrientation;
    delete[] m_pfOrientAngle;
    delete[] m_pfSVectors;
    delete[] m_pfChi;

    delete[] m_pnRingStart;
    delete[] m_pnRingEnd;
};



int Cabsorbsurface::nPrint(int nMode,Cstring& sFilename,char* cpLabel,bool bPrintPoints) {
    
    int nx,ny;
    int nAtomCt;
    double a3fTemp[3];



    for (nx=0;nx<m_nNumSPoints;nx++) {
        printf("#%d Level %d  Neighbors = [ ",nx,m_poSPoints[nx].nRow);
        for (ny=0;ny<m_poSPoints[nx].nAdjCt;ny++) 
            printf("#%d(%d) ",m_poSPoints[nx].nGetAdj(ny),m_poSPoints[m_poSPoints[nx].nGetAdj(ny)].nRow);
        printf(" ] \n");
    };




    if (nMode == 0) {
        FILE* pFOut;
        pFOut = fopen(sFilename.string(),"w+t");
        if (!pFOut)
            return 0;
        fprintf(pFOut,"begin %s \ncolour green \n",cpLabel);
        nAtomCt=0;
        
        for (nx=0;nx<m_nNumSPoints;nx++) {
            
            // Travel out to each adjacency, and back again.
            
            fprintf(pFOut,"m %8.3f%8.3f%8.3f \n",
                m_poSPoints[nx].a3fCoords[0]*5.0,
                m_poSPoints[nx].a3fCoords[1]*5.0,
                m_poSPoints[nx].a3fCoords[2]*5.0);
            
            
            for (ny=0;ny<m_poSPoints[nx].nAdjCt;ny++) {
                fprintf(pFOut,"l  %8.3f%8.3f%8.3f \n",
                    m_poSPoints[m_poSPoints[nx].nGetAdj(ny)].a3fCoords[0]*5.0,
                    m_poSPoints[m_poSPoints[nx].nGetAdj(ny)].a3fCoords[1]*5.0,
                    m_poSPoints[m_poSPoints[nx].nGetAdj(ny)].a3fCoords[2]*5.0);
                fprintf(pFOut,"l %8.3f%8.3f%8.3f \n",
                    m_poSPoints[nx].a3fCoords[0]*5.0,
                    m_poSPoints[nx].a3fCoords[1]*5.0,
                    m_poSPoints[nx].a3fCoords[2]*5.0);
                vSubVec3DVec3D(
                    m_poSPoints[nx].a3fCoords,
                    m_poSPoints[m_poSPoints[nx].nGetAdj(ny)].a3fCoords,
                    a3fTemp);
                if (fLenVec3D(a3fTemp)>1.0) {
                    printf("Length too long (%f) between (%d,%d)\n",fLenVec3D(a3fTemp),nx,m_poSPoints[nx].nGetAdj(ny));
                };
                
            };
            
        };
        if (bPrintPoints) {
            fprintf(pFOut,"colour white \n");
            for (nx=0;nx<m_nNumRefs;nx++) 
                fprintf(pFOut,"dot %f %f %f \n",5.0*m_pfSVectors[nx*3+0],5.0*m_pfSVectors[nx*3+1],5.0*m_pfSVectors[nx*3+2]);
        };
        fprintf(pFOut,"end_object \n");
        fclose(pFOut);
    };

    if (nMode == 1) {
        FILE* pFX;
        FILE* pFY;
        FILE* pFZ;
        double a3fVec[3];
        int a2nGrid[2];
        int a2nPoint[2];
        a2nGrid[0] = 100;
        a2nGrid[1] = 100;
        pFX = fopen("c:\\temp\\xdata.txt","w+t");
        pFY = fopen("c:\\temp\\ydata.txt","w+t");
        pFZ = fopen("c:\\temp\\zdata.txt","w+t");
        if ((!pFX) || (!pFY) || (!pFZ))
            return 1;
        for (a2nPoint[0]=0;a2nPoint[0]<a2nGrid[0];a2nPoint[0]++) {
            double fTheta = (Gs_dPI*2.0*a2nPoint[0])/a2nGrid[0];
            for (a2nPoint[1]=0;a2nPoint[1]<a2nGrid[1];a2nPoint[1]++) {
                double fPsi = (Gs_dPI*a2nPoint[1])/a2nGrid[1]-Gs_dPI*0.5;
                double a3x3fPoints[3][3];
                double a3fCoeffs[3];
                int nTriangle;
                a3fVec[0] = cos(fTheta)*cos(fPsi);
                a3fVec[1] = sin(fTheta)*cos(fPsi);
                a3fVec[2] = sin(fPsi);
                for (nTriangle=0;nTriangle<m_nNumTriangles;nTriangle++) {
                    
                    // Load the triangle points.
                    for (nx=0;nx<3;nx++)
                        vCopyVec3D(m_poSPoints[m_a3pnTriangles[nx][nTriangle]].a3fCoords,a3x3fPoints[nx]);
                    if (bGetPointCoeffs(a3x3fPoints,a3fVec,a3fCoeffs)) 
                        break;
                };
                if (nTriangle==m_nNumTriangles) {
                    printf("WARNING: Could not find printing vector %d %d\n",a2nPoint[0],a2nPoint[1]);
                } else {
                    fprintf(pFX,"%f\t",a3fVec[0]);
                    fprintf(pFY,"%f\t",a3fVec[1]);
                    fprintf(pFZ,"%f\t",a3fVec[2]);                   
                };
            };
            fprintf(pFX,"\n");
            fprintf(pFY,"\n");
            fprintf(pFZ,"\n");
        };
        fclose(pFX);
        fclose(pFY);
        fclose(pFZ);
    };

    return 0;
};

int Cabsorbsurface::nGetRingStartEnd(int nResoRing,int& nRingStart,int& nRingEnd) {
    int nx,ny;
    nRingStart = -1;
    nRingEnd = -1;
    for (nx=0;nx<m_nNumSPoints;nx++) {
        ny = m_poSPoints[nx].nRow;
        if ((ny==nResoRing) && (nRingStart==-1))
            nRingStart = nx;
        if (ny==nResoRing)
            nRingEnd = nx;
    };    
    return 0;
};


double Cabsorbsurface::fGetVectorAngle(double* pfVector) {
    double f0;

    f0 = sqrt(pfVector[0]*pfVector[0] + pfVector[1]*pfVector[1]);
    if (pfVector[1]>=0.0)
        return acos(pfVector[0]/f0);
    else
        return (Gs_dPI*2.0 - acos(pfVector[0]/f0));
};



int Cabsorbsurface::nBuildPoints(char* cpTitle,int nMode) {
    int nx;

    int pnMode0[] = { 4,8,-1};
    int pnMode1[] = { 4,8,8,-1};
    int pnMode2[] = { 4,8,8,8,-1};
    int pnMode3[] = { 4,8,16,16,-1};
    double pfOffset0[] = { 0.0,0.0};
    double pfOffset1[] = { 0.0,0.0,12.5};
    double pfOffset2[] = {0.0,0.0,24.5,0.0};
    double pfOffset3[] = {0.0,0.0,0.0,12.25};
    int* pnModes[] = { pnMode1,pnMode2,pnMode3};
    double* pfOffsets[] = {pfOffset1,pfOffset2,pfOffset3};
    Cstring sTemp;

    nBuildPoints(pnModes[nMode],pfOffsets[nMode]);
    
    return 0;    

return 0;
};



int Cabsorbsurface::nBuildPoints(int* pnRingCounts,double* pfChiOffsets) {
    int nx,ny,nz;
    double f0,f1,f2;
    
    double f2Theta;
    double fChi;
    int    nRing;
    int    nPoint1,nPoint2;
    int    nAdj;

    int  nTotalPoints;
    int  nTotalRings;
    int* anPointsRow;
    int* aanConnections[2];
    double* aa3fPoints[3];
    double* afChi;
    double* af2Theta;

    
    int nConnects;
    int nPoints;
    int nNewPoints;
    int nLastPoints;
    int nLastNewPoints;

    // Count total number of points.
    nTotalRings = 0;
    for (nTotalPoints=0,nx=0;pnRingCounts[nx]>0;nx++) {
        nTotalPoints += pnRingCounts[nx];
        nTotalRings++;
    };

    nTotalPoints += 2;
    for (nx=0;nx<2;nx++)
        aanConnections[nx] = new int[nTotalPoints*8];
    for (nx=0;nx<3;nx++) 
        aa3fPoints[nx] = new double[nTotalPoints];
    afChi = new double[nTotalPoints];    
    af2Theta = new double[nTotalPoints];
    anPointsRow = new int[nTotalPoints];
    m_nNumSPoints = nTotalPoints;

    // Add first point.
    nPoints = 1;
    aa3fPoints[0][0] = 0.0;
    aa3fPoints[1][0] = 0.0;
    aa3fPoints[2][0] = 1.0;
    afChi[0] = 0.0;
    af2Theta[0] = 0.0;
    anPointsRow[0] = 0;

    nConnects = 0;
    nLastNewPoints = 1;
    nLastPoints = 0;
    nPoints = 1;
    
    for (nRing=0;pnRingCounts[nRing]>0;nRing++) {
        nNewPoints = pnRingCounts[nRing];
        // Make connections in the ring.
        // Build points.
        for (nx=0;nx<nNewPoints;nx++) {
            aanConnections[0][nConnects] = nPoints + nx;
            aanConnections[1][nConnects] = nPoints + ((nx + 1) % nNewPoints);
            nConnects++;
            fChi = pfChiOffsets[nRing] + nx*360/nNewPoints;
            if (fChi>=360.0)
                fChi-=360.0;
            afChi[nPoints+nx] = fChi;            
            f2Theta = ((nRing+1)*90.0/nTotalRings);
            f2Theta *= m_nMax2Theta/90.0;
            af2Theta[nPoints+nx] = f2Theta;
            
            aa3fPoints[0][nPoints + nx] = cos(Gs_dRADIANS_PER_DEGREE*fChi)*sin(f2Theta*Gs_dRADIANS_PER_DEGREE);
            aa3fPoints[1][nPoints + nx] = sin(Gs_dRADIANS_PER_DEGREE*fChi)*sin(f2Theta*Gs_dRADIANS_PER_DEGREE);
            aa3fPoints[2][nPoints + nx] = cos(f2Theta*Gs_dRADIANS_PER_DEGREE);
            anPointsRow[nPoints+nx] = nRing+1;
        };
        // Make connetions to previous ring.
        for (nPoint1=0;nPoint1<nNewPoints;nPoint1++) {
            int nMaxConnects = 2;
            int nClosestAbove = 0;
            int nClosestBelow = 1;
            int nClosest = 0;
            // See if there is some point in the previous ring with 'almost' the exact same angle.
            // If so, we will be connecting exactly to it.
            for (nPoint2=0;nPoint2<nLastNewPoints;nPoint2++) {
                if (fabs(afChi[nLastPoints+nPoint2]-afChi[nPoints+nPoint1])<0.1) {
                    nMaxConnects = 1;
                    nClosest = nPoint2;
                };
                f0 = afChi[nLastPoints+nPoint2]-afChi[nPoints+nPoint1];
                if ((f0>0.0) && (f0<0.8*360.0/nLastNewPoints))
                    nClosestAbove = nPoint2;
                if ((f0<0.0) && (-f0<0.8*360.0/nLastNewPoints))
                    nClosestBelow = nPoint2;
            };
            if ((nMaxConnects==1) || (nLastNewPoints==1)) {
                aanConnections[0][nConnects] = nLastPoints + nClosest;
                aanConnections[1][nConnects] = nPoints + nPoint1;
                nConnects++;
            } else {
                aanConnections[0][nConnects] = nLastPoints + nClosestAbove;
                aanConnections[1][nConnects] = nPoints + nPoint1;
                nConnects++;
                aanConnections[0][nConnects] = nLastPoints + nClosestBelow;
                aanConnections[1][nConnects] = nPoints + nPoint1;
                nConnects++;
            };
        };
        nLastNewPoints = nNewPoints;
        nLastPoints = nPoints;
        nPoints += nNewPoints;
    };
    // Build the antartic pole.
    aa3fPoints[0][nPoints] = 0.0;
    aa3fPoints[1][nPoints] = 0.0;
    aa3fPoints[2][nPoints] = -1.0;
    anPointsRow[nPoints] = 200;
    nPoints++;
    // Connect to antartic point.
    for (nPoint1=0;nPoint1<nLastNewPoints;nPoint1++) {
        aanConnections[0][nConnects] = nLastPoints + nPoint1;
        aanConnections[1][nConnects] = nPoints-1;
        nConnects++;
    };
    if (m_pnC)
        delete[] m_pnC;
    m_pnC = new int[m_nNumSPoints*m_nNumS0Points];

    if (m_poSPoints)
        delete[] m_poSPoints;
    m_poSPoints = new CSpoint[m_nNumSPoints];
    

    // Fill in equilateral points.
    for (nPoint1=0;nPoint1<m_nNumSPoints;nPoint1++) {
        for (nx=0;nx<3;nx++)
            m_poSPoints[nPoint1].a3fCoords[nx] = aa3fPoints[nx][nPoint1];
        nx = 0;
        for (nAdj = 0;nAdj<nConnects;nAdj++) {
            if (aanConnections[0][nAdj]==nPoint1) 
                m_poSPoints[nPoint1].anAdj[nx++] = aanConnections[1][nAdj];
            else if (aanConnections[1][nAdj]==nPoint1) 
                m_poSPoints[nPoint1].anAdj[nx++] = aanConnections[0][nAdj];
        };
        m_poSPoints[nPoint1].nAdjCt = nx;
        m_poSPoints[nPoint1].nRow = anPointsRow[nPoint1];
        m_poSPoints[nPoint1].nCol = nPoint1;
        m_poSPoints[nPoint1].fChi = afChi[nPoint1];
        m_poSPoints[nPoint1].f2Theta = af2Theta[nPoint1];
    };

    for (nx=0;nx<2;nx++)
        delete[] aanConnections[nx];
    for (nx=0;nx<3;nx++)
        delete[] aa3fPoints[nx];
    delete[] afChi;
    delete[] af2Theta;
    delete[] anPointsRow;
    return 0;

};



int Cabsorbsurface::nBuildTriangles() {

    int nVert,nPoint;
    int a3nPoints[3];
    int a3nPointsBuf[3];
    bool bFound;
    bool bDivided;
    int nNumVerts;
    int nPoint1i,nPoint2i,nPoint3i,nPoint3;
    int nTriangleCt;
    double a3fVec[3];
    double a3fCoeffs[3];
    double a3x3fPoints[3][3];
    int nx;

    // Allocate the triangles array.  It has an exact number of entries.
    // We must however, count the number of degrees of vertices in the system.
    for (nPoint=0,m_nNumTriangles=0;nPoint<m_nNumSPoints;nPoint++)
        m_nNumTriangles += m_poSPoints[nPoint].nAdjCt;
    if (m_nNumTriangles % 3) {
        cout << "ERROR: triangularization";
        return 1;
    };
    m_nNumTriangles/=3;                  
    for (nx=0;nx<3;nx++) {
        if (m_a3pnTriangles[nx])
            delete[] m_a3pnTriangles[nx];
        m_a3pnTriangles[nx] = new int[m_nNumTriangles];
    };
    nTriangleCt = 0;

    for (nPoint=0;nPoint<m_nNumSPoints;nPoint++) {
        a3nPoints[0] = nPoint;
        nNumVerts = m_poSPoints[nPoint].nAdjCt;
        for (nPoint1i=0;(nPoint1i<nNumVerts);nPoint1i++) {
            a3nPoints[1] = m_poSPoints[nPoint].nGetAdj(nPoint1i);
            for (nPoint2i=0;(nPoint2i<nNumVerts);nPoint2i++) {
                a3nPoints[2] = m_poSPoints[nPoint].nGetAdj(nPoint2i);
                if ((a3nPoints[1] != a3nPoints[2] ) && (m_poSPoints[a3nPoints[1]].nRow==m_poSPoints[a3nPoints[2]].nRow) && (m_poSPoints[a3nPoints[0]].nRow!=m_poSPoints[a3nPoints[1]].nRow)) {
                    // Build the three points in a3x3fPoints.
                    for (nx=0;nx<3;nx++)
                        vCopyVec3D(m_poSPoints[a3nPoints[nx]].a3fCoords,a3x3fPoints[nx]);

                    // Make sure that the two vertices on the same level are adjacent.
                    bDivided = ! m_poSPoints[a3nPoints[1]].bIsAdj(a3nPoints[2]);

                    if (!bDivided) {
                        for (nx=0;nx<3;nx++)
                            a3nPointsBuf[nx]= a3nPoints[nx];
                        // Order the proposed triangle vertices.
                        if (a3nPointsBuf[0]>a3nPointsBuf[1])
                            std::swap(a3nPointsBuf[0],a3nPointsBuf[1]);
                        if (a3nPointsBuf[1]>a3nPointsBuf[2])
                            std::swap(a3nPointsBuf[1],a3nPointsBuf[2]);
                        if (a3nPointsBuf[0]>a3nPointsBuf[1])
                            std::swap(a3nPointsBuf[0],a3nPointsBuf[1]);
                        // Search for the proposed triangle amoung the triangles we have already stored.
                        for (nx=0;nx<nTriangleCt;nx++) {
                            if ((a3nPointsBuf[0]==m_a3pnTriangles[0][nx]) &&
                                (a3nPointsBuf[1]==m_a3pnTriangles[1][nx]) &&
                                (a3nPointsBuf[2]==m_a3pnTriangles[2][nx])) 
                                break;
                        };
                        if (nx==nTriangleCt) {
                            if (nTriangleCt==m_nNumTriangles) {
                                cout << "ERROR: triangularization (too many triangles).\n" << flush;
                                return 1;
                            };
                            for (nx=0;nx<3;nx++)
                                m_a3pnTriangles[nx][nTriangleCt] = a3nPointsBuf[nx];
                            nTriangleCt++;
                        };                       
                    };
                };
            };
        };
    };

    if (m_nNumTriangles != nTriangleCt) {
        cout << "ERROR: triangularization (could not find enough triangles).\n" << flush;
        return 1;
    };
    return 0;

};


// Function returns proportional coefficients for each of four bounding vertices.
// It also tests to make sure the given point is inside the triangle determined by the three points.
bool Cabsorbsurface::bGetPointCoeffs(double a3x3fPoints[3][3],double a3fPoint[3],double a3fCoeffs[3]) {
    double a3fLens[3];
    int a3x2nVecs[3][2] = {{0,1},{1,2},{2,0}};
    int a3nIndex[3] = { 1,2,0};     // Gives which a3fLens[] element used to compute the coefficient.
    double a2x3fVecs[2][3];
    double a3x3fCross[3][3];
    double fTotalLengths;


    int nx,ny;
    double a3fTemp[3];

    
    // Check to see that the vector has the same sign when dotted with each of the 3 cone surface normals.
    int nNegSigns=0;
    int nPosSigns=0;
    int nTested=0;
    for (nx=0;nx<3;nx++) {
        vCross3D(a3x3fPoints[a3x2nVecs[nx][0]],a3x3fPoints[a3x2nVecs[nx][1]],a3fTemp);
        if (fDot3D(a3fTemp,a3fPoint)<0.0)
            nNegSigns++;
        else
            nPosSigns++;
    };
    if ((3 != nPosSigns) && (3 != nNegSigns))
        return FALSE;

    // Get the volume of each of the three triangles
    // Also compute the cross product vectors.
    fTotalLengths =0.0;
    for (nx=0;nx<3;nx++) {
        for (ny=0;ny<2;ny++) 
            vSubVec3DVec3D(a3x3fPoints[a3x2nVecs[nx][ny]],a3fPoint,a2x3fVecs[ny]);
        vCross3D(a2x3fVecs[0],a2x3fVecs[1],a3x3fCross[nx]);
        fTotalLengths += (a3fLens[nx] = fLenVec3D(a3x3fCross[nx]));
    };


    for (nx=0;nx<3;nx++) {
        a3fCoeffs[nx] = a3fLens[a3nIndex[nx]]/fTotalLengths;
    };
    return TRUE;
};

int Cabsorbsurface::nFindTriangleCoeffs(Creflnlist& oReflnlist) {
    int nTriangle;
    int nx,ny;
    int nStat;
    int nRef;
    double a3x3fPoints[3][3];
    double a3fS[3];
    double a3fCoeffs[3];

    nStat = 0;

    for (nRef=0;(nRef<m_nNumRefs) && (!nStat);nRef++) {
        vCopyVec3D(m_pfSVectors+nRef*3,a3fS);
        
        // Get the bounding triangles.
        for (nTriangle=0;nTriangle<m_nNumTriangles;nTriangle++) {
            
            // Load the triangle points.
            for (nx=0;nx<3;nx++)
                vCopyVec3D(m_poSPoints[m_a3pnTriangles[nx][nTriangle]].a3fCoords,a3x3fPoints[nx]);
            if (bGetPointCoeffs(a3x3fPoints,a3fS,a3fCoeffs)) 
                break;
        };
        // For each reflection in the list, we must find a bounding triangle.
        // Otherwise, the point was invalid.
        if (nTriangle==m_nNumTriangles)  {
            if (!nStat)
                printf("ERROR: Could not find reflection vector #%d\n",nRef);
            nStat++;
        } else {
            int a2nRow[2];
            int a2nPoint[2];
            double a2fDist[2];
            int nRow;
            
            // Place the coefficients and points in the array.
            
            a2fDist[0] = -1000;
            a2fDist[1] = -1000;
            a2nRow[0] = 100;
            a2nRow[1] = 100;
            
            for (nx=0;nx<3;nx++) {
                m_a3pnSPoints[nx][nRef] = m_a3pnTriangles[nx][nTriangle];
                m_a3pfSCoeffs[nx][nRef] = a3fCoeffs[nx];
                m_poSPoints[m_a3pnSPoints[nx][nRef]].nContrib++;
                nRow = m_poSPoints[m_a3pnSPoints[nx][nRef]].nRow;
                for (ny=0;ny<2;ny++) {
                    if ((nRow==a2nRow[ny]) || (a2nRow[ny]==100)) {
                        if (a2fDist[ny]<a3fCoeffs[nx]) {
                            a2fDist[ny] = m_a3pfSCoeffs[nx][nRef];
                            a2nPoint[ny] = m_a3pnSPoints[nx][nRef];
                            a2nRow[ny] = nRow;
                            break;
                        };
                    };
                };                                
            };
            if ((a2nRow[0]==100) || (a2nRow[1]==100)) {
                if (!nStat)
                    printf("ERROR:  Bad input reflection ranges!\n");
                nStat++;
                
                a2nPoint[1] = a2nPoint[0];
            };
            
            m_a2pnClosestS[0][nRef] = a2nPoint[0];
            m_a2pnClosestS[1][nRef] = a2nPoint[1];
        };
    };
    return nStat;
};

int Cabsorbsurface::nCheckMeshDimensions(Creflnlist& oList,int* pnIndex,int* pnReject,int nMinRedundancy) {
    int nPoints;

    int nLastHKL;
    int nThisHKL;
    int nRef;
    int nRefSort;

    int nRing;
    int nStart,nEnd;
    Cstat oStat;
    int nRingSuccess;
    int nRingFailure;

    int nx,ny,nz;
    double f0;

    nPoints = m_nNumSPoints*m_nNumS0Points;

    for (nx=0;nx<nPoints;nx++)
        m_pnC[nx] = 0;

    // Go through all of the points, and determine how many contributors each
    // triangle coeff. has.   
    nLastHKL = oList[pnIndex[0]].nGetField(oList.m_nFI_nPackedHKL);   
    nStart = 0;
    for (nRefSort=0;nRefSort<m_nNumRefs+1;nRefSort++) {
        
        nRef=pnIndex[nRefSort];        
        if (nRefSort<m_nNumRefs)
            nThisHKL = oList[nRef].nGetField(oList.m_nFI_nPackedHKL);
        if ((nRefSort!=m_nNumRefs) && (nThisHKL==nLastHKL))
            continue;
        nEnd=nRefSort-1;

        ny = 0;
        for (nx=nStart;nx<=nEnd;nx++) {
            if (pnReject[nx]==0)
                ny++;
        };
        
        if (ny>=nMinRedundancy) {
            for (nx=nStart;nx<=nEnd;nx++) {
                if (pnReject[nx]==0) {
                    for (nz=0;nz<2;nz++) 
                        m_pnC[m_pnClosestS0[pnIndex[nx]]*m_nNumSPoints + m_a2pnClosestS[nz][pnIndex[nx]]]+=ny;
                };
            };
        };
        
        nStart=nRefSort;
        if (nRefSort<m_nNumRefs)
            nLastHKL =  oList[nRef].nGetField(oList.m_nFI_nPackedHKL);        
    };

    // Make sure that *most* triangles have the minimum number of points.  


    nRingSuccess = 0;
    nRingFailure = 0;
    for (nRing=1;1;nRing++) {
        nGetRingStartEnd(nRing,nStart,nEnd);
        if ((nStart==-1) || (nRing==-1))
            break;   
        for (nx=0;nx<m_nNumS0Points;nx++) {
            oStat.vClear();
            for (ny=nStart;ny<=nEnd;ny++)
                oStat.vAdd((double) m_pnC[nx*m_nNumSPoints + ny]);
            oStat.nMarkStat(false,Cstat::S_ABOVEBELOW,3.0);
            f0 = oStat.fAverage();
            if (f0<m_fMinRingContribAverage) {
                nRingFailure++;
            } else {  
                for (ny=0;ny<nEnd-nStart+1;ny++) {
                    f0 = oStat[ny];
                    if ((f0<m_fMinContrib) && (nEnd-nStart+1>4))
                        break;
                };
                if (ny!=nEnd-nStart+1)
                    nRingFailure++;
                else
                    nRingSuccess++;
            };
        };
    };
    if (nRingSuccess/((double) nRingSuccess+nRingFailure)>=m_fMinNonHarmonicRings)
        return 0;   
    else
        return 1;
};



int Cabsorbsurface::nLoad(Creflnlist& oReflnlist,int a3nS[3],int nBatchIdx,int a3nS0[3],tagBatch* poBatches) {
    double f0,f1;
    int nx,ny;
    int nRef;
    int nStat;

    double a3fS[3];
    double a3fVec1[3];
    double a3fVec2[3];
    double a3x3fMat[3][3];

    int nFI_fObsRotMid;
    int nFI_nBatchIdx;
    
    // Initialize the reflection index lists.
    for (nx=0;nx<3;nx++) {
        delete[] m_a3pfSCoeffs[nx];
        delete[] m_a3pnSPoints[nx];        
        m_a3pfSCoeffs[nx] = new double[oReflnlist.nGetNumReflns()];
        m_a3pnSPoints[nx] = new int[oReflnlist.nGetNumReflns()];
    };
    for (nx=0;nx<2;nx++) {
        delete[] m_a2pfS0Coeffs[nx];
        delete[] m_a2pnS0Points[nx];
        m_a2pfS0Coeffs[nx] = NULL;
        m_a2pnS0Points[nx] = NULL;
        if (m_nDegreesPerS0>0) {
            // Only allocate these if we are not using harmonics.
            m_a2pfS0Coeffs[nx] = new double[oReflnlist.nGetNumReflns()];
            m_a2pnS0Points[nx] = new int[oReflnlist.nGetNumReflns()];
        };
    };
    delete[] m_a2pnClosestS[0];
    m_a2pnClosestS[0] = new int[oReflnlist.nGetNumReflns()];
    delete[] m_a2pnClosestS[1];
    m_a2pnClosestS[1] = new int[oReflnlist.nGetNumReflns()];
    delete[] m_pnClosestS0;
    m_pnClosestS0 = new int[oReflnlist.nGetNumReflns()];

    delete[] m_pnOrientation;
    m_pnOrientation = new int[oReflnlist.nGetNumReflns()];
    delete[] m_pfOrientAngle;
    m_pfOrientAngle = new double[oReflnlist.nGetNumReflns()];

    delete[] m_pfSVectors;
    m_pfSVectors = new float[oReflnlist.nGetNumReflns()*3];
    delete[] m_pfChi;
    m_pfChi = new float[oReflnlist.nGetNumReflns()];

    m_nNumRefs=oReflnlist.nGetNumReflns();

    nFI_fObsRotMid = oReflnlist.nGetFieldIndex(Creflnlist::ms_sfObsRotMid);
    
    nStat = 1;
    if ((a3nS0[0]>=0) && (a3nS0[1]>=0) && (a3nS0[2]>=0) && (nFI_fObsRotMid>=0)) {

        Cstring sNextBatchOrient = "";
        Cstring sLastBatchOrient = "";
        Cstring sBatchOrient;

        double a3fStartRot[3];      // S0 at start of rotation scan.
        double a3fEndRot[3];        // S0 at end of rotation scan.
        double a3fOppositeRot[3];   // S0 far away from a3fStartRot.
        double a3fAvgRot[3];        // Average S0 vector.
        double a3fRot[3];           // Temporary S0 vector
        double a3fScatt[3];         // Temporary S vector.
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
        double fPerpVecSign;        // Might need to multiply the perpendicular vector by a sign.

        int nRefsInOrient;          // Number of reflections in this particular orientation.
        int nFirstBatch;            // First batch in orientation.
        int nLastBatch;             // Last batch in orientation.
        int nBatch;                 // Temporary used to get batch.
        int nNumCoeffs;             // Coeffs. induced. in orientation.
        bool bOrientFound;          // Flag used in new orientation finding loop.
        bool bCompleteRotation;     // Did the S0 points go all the way around 360 degrees?

        double fMaxDot;             // Used to error check the 'perpendicular' vector (a3fPerpS0)
        double fMinDot;             // Used to error check the 'perpendicular' vector (a3fPerpS0)
        CS0point* poS0;             // S0 temporary pointer.

        int nNumWithProblems;       // Reflections that had Suspecious S/S0 vectors.
        int nNumWithSeriousProblems;// Reflections that had BAD S/S0 vectors.

        nNumWithProblems = 0;
        nNumWithSeriousProblems = 0;


        printf("\n\n");
        printf("Triangular Mesh Scans/Coefficients\n");
        printf("--------------------------------------------------------------------------\n");
        printf("   Scan First  Last    Rot    Rot    Rot     S0     S0    S0    S0    S0\n");
        printf("        batch batch  start    end  range coeffs radius   [x]   [y]   [z]\n");
        printf("--------------------------------------------------------------------------\n");

        m_nNumS0Points = 0;       
        m_nNumOrients = 0;
        do {
            sLastBatchOrient = sNextBatchOrient;
            bOrientFound = FALSE;
            for (nRef=0;nRef<oReflnlist.nGetNumReflns();nRef++) {
                nBatch = poBatches[oReflnlist[nRef].nGetField(nBatchIdx)].m_nBatch;
                sBatchOrient = poBatches[oReflnlist[nRef].nGetField(nBatchIdx)].m_sScan;
                if (((!bOrientFound) || (sBatchOrient<sNextBatchOrient)) &&
                    ((m_nNumOrients==0) || (sBatchOrient>sLastBatchOrient))) {
                    sNextBatchOrient = sBatchOrient;
                    bOrientFound = TRUE;
                };
            };

            if ((m_nNumOrients==0) || (sNextBatchOrient>sLastBatchOrient)) {
                // A new orientation was found.  
                
                // Find a maximum chord on the circle described by
                // the tips of the S0 vectors.  (The tips will lie in a plane,
                // but not necc. a plane passing through the center of the unit sphere).
                nRefsInOrient = 0;
                nFirstBatch = 1000;
                nLastBatch = -1;
                vZeroMat(3,1,a3fAvgRot);
                for (nRef=0;nRef<oReflnlist.nGetNumReflns();nRef++) {
                    nBatch = poBatches[oReflnlist[nRef].nGetField(nBatchIdx)].m_nBatch;
                    sBatchOrient = poBatches[oReflnlist[nRef].nGetField(nBatchIdx)].m_sScan;
                    if (sBatchOrient==sNextBatchOrient) {
                        fRot = oReflnlist[nRef].fGetField(nFI_fObsRotMid);
                        for (nx=0;nx<3;nx++)
                            a3fRot[nx] = -oReflnlist[nRef].fGetField(a3nS0[nx]);

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
                vMulVec3DScalar(a3fAvgRot,1.0/nRefsInOrient,a3fAvgRot);

                vSubVec3DVec3D(a3fOppositeRot,a3fStartRot,a3fRot);
                fMaxLength = fLenVec3D(a3fRot);
                fLengthRot = fOppositeRot - fStartRot;
                fLengthRot360 = min(360.0,fEndRot-fStartRot);
                

                // Compute the arc length from this chord.
                if (m_nDegreesPerS0<=0) {
                    nNumCoeffs = -2*m_nDegreesPerS0 + 1;
                    bCompleteRotation = FALSE;
                } else if (fLengthRot>0.0001) {
                    fCircleRadius = sqrt(fMaxLength*fMaxLength/2.0/(1.0 - cos(Gs_dRADIANS_PER_DEGREE*fLengthRot)));
                    nNumCoeffs = 1 + max(0,(int) (fLengthRot360*Gs_dRADIANS_PER_DEGREE*fCircleRadius/(m_nDegreesPerS0*Gs_dRADIANS_PER_DEGREE)));
                    // If we went around in a complete circle, then reduce the number of coeffs by 1.
                    if (fabs(fLengthRot360-360.0)<1.0) {
                        nNumCoeffs = max(1,nNumCoeffs-1);
                        bCompleteRotation = TRUE;
                    };
                } else
                    nNumCoeffs = 1;

                // Compute perpendicular S0 vector.
                vSubVec3DVec3D(a3fOppositeRot,a3fAvgRot,a3fVec1);
                vSubVec3DVec3D(a3fStartRot,a3fAvgRot,a3fVec2);                
                if ((fLenVec3D(a3fVec1)>0.001) && (fLenVec3D(a3fVec2)>0.001)) 
                    vCross3D(a3fVec1,a3fVec2,a3fPerpS0);
                else
                    vCopyVec3D(a3fAvgRot,a3fPerpS0);
                fNormVec3D(a3fPerpS0);


                // Print out information on this scan.
                printf(" %3s???   %3.3d   %3.3d %6.2f %6.2f %6.2f %6d %6.2f %5.2f %5.2f %5.2f\n",
                    sNextBatchOrient.string(),
                    nFirstBatch,
                    nLastBatch,
                    fStartRot,fEndRot,fLengthRot360,nNumCoeffs,
                    fCircleRadius,a3fPerpS0[0],a3fPerpS0[1],a3fPerpS0[2]
                    );
                
                if (m_nDegreesPerS0>0) {
                    // Expand S0 points array to include the new points.
                    poS0 = new CS0point [ nNumCoeffs + m_nNumS0Points ] ;
                    if (m_poS0Points) {
                        for (nx=0;nx< m_nNumS0Points;nx++)
                            poS0[nx] = m_poS0Points[nx];
                        delete[] m_poS0Points;
                    };
                    m_poS0Points = poS0;
                    
                    // Compute S0 fields and vectors
                    for (nx=0;nx<nNumCoeffs;nx++) {
                        poS0 = &m_poS0Points[nx + m_nNumS0Points];
                        poS0->nOrient = m_nNumOrients;
                        poS0->nPoint = nx;
                        poS0->nCumulPoint = nx + m_nNumS0Points;
                        poS0->nContrib = 0;
                        if (bCompleteRotation)
                            poS0->fDispRot = nx*fLengthRot360/nNumCoeffs;
                        else if (nNumCoeffs!=1)
                            poS0->fDispRot = nx*fLengthRot360/(nNumCoeffs-1);
                        else
                            poS0->fDispRot = fLengthRot360/2.0;
                        poS0->fRot = poS0->fDispRot + fStartRot;
                        vConvRotVec3DMat3D(poS0->fDispRot,&a3fPerpS0[0],&a3x3fRotMat[0][0]);
                        vMulMat3DVec3D(a3x3fRotMat,a3fStartRot,poS0->a3fCoords);
                    };
                };


                // Compute the rotation matrix that would take an S0 vector at fStartRot to <0,0,1>
                vBuildBasis3D(a3fStartRot,a3x3fRotMat);
                vNormMat3D(a3x3fRotMat);
                vTranMat3D(a3x3fRotMat);
                vZeroMat3D(&a3x3fMat[0][0]);
                a3x3fMat[1][0] = 1.0;
                a3x3fMat[2][1] = 1.0;
                a3x3fMat[0][2] = 1.0;
                vMulMat3DMat3D(a3x3fMat,a3x3fRotMat,a3x3fOffsetRot);

                // Assign coeffs contributions to all reflections.
                fMaxDot = -1.0;
                fMinDot = 1.0;
                fPerpVecSign = 1.0;
                for (nRef=0;nRef<oReflnlist.nGetNumReflns();nRef++) {
                    nBatch = poBatches[oReflnlist[nRef].nGetField(nBatchIdx)].m_nBatch;
                    sBatchOrient = poBatches[oReflnlist[nRef].nGetField(nBatchIdx)].m_sScan;
                    if (sBatchOrient==sNextBatchOrient) {
                        m_pnOrientation[nRef] = m_nNumOrients;
                        m_pfOrientAngle[nRef] = oReflnlist[nRef].fGetField(nFI_fObsRotMid);
                        for (nx=0;nx<3;nx++)
                            a3fRot[nx] = -oReflnlist[nRef].fGetField(a3nS0[nx]);
                        for (nx=0;nx<3;nx++)
                            a3fScatt[nx] = oReflnlist[nRef].fGetField(a3nS[nx]);        
                        f0 = fDot3D(a3fPerpS0,a3fRot);
                        fMaxDot = max(f0,fMaxDot);
                        fMinDot = min(f0,fMinDot);

                        if (m_nDegreesPerS0>0) {
                            double fLargest1;
                            double fLargest2;
                            int    nLargest1;
                            int    nLargest2;
                            // Scan through all of the points, and find the two that have
                            // the largest dot products with the original.
                            nLargest1 = -1;
                            nLargest2 = -1;
                            for (nx=0;nx<nNumCoeffs;nx++) {                                                        
                                f0 = fDot3D(a3fRot,m_poS0Points[m_nNumS0Points + nx].a3fCoords);
                                if (nLargest1 == -1) {
                                    nLargest1 = nx;
                                    fLargest1 = f0;
                                } else if (nLargest2 == -1) {
                                    nLargest2 = nx;
                                    fLargest2 = f0;
                                } else if ((fLargest1<fLargest2) && (f0>fLargest1)) {
                                    nLargest1 = nx;
                                    fLargest1 = f0;
                                } else if ((fLargest2<fLargest1) && (f0>fLargest2)) {
                                    nLargest2 = nx;
                                    fLargest2 = f0;
                                };
                            };
                            if (nNumCoeffs == 1) {
                                // There is only one point.
                                m_a2pfS0Coeffs[0][nRef] = 1.0;
                                m_a2pfS0Coeffs[1][nRef] = 0.0;
                                m_a2pnS0Points[0][nRef] = m_nNumS0Points + nLargest1;
                                m_a2pnS0Points[1][nRef] = m_nNumS0Points + nLargest1;
                                m_poS0Points[m_nNumS0Points + nLargest1].nContrib++;
                                m_pnClosestS0[nRef] = m_a2pnS0Points[0][nRef];
                            } else {
                                // Find the distance to each of the two bounding points.
                                vSubVec3DVec3D(a3fRot,m_poS0Points[m_nNumS0Points + nLargest1].a3fCoords,a3fVec1);
                                vSubVec3DVec3D(a3fRot,m_poS0Points[m_nNumS0Points + nLargest2].a3fCoords,a3fVec2);
                                f0 = fLenVec3D(a3fVec1);
                                f1 = fLenVec3D(a3fVec2);
                                m_a2pfS0Coeffs[0][nRef] = f1/(f0+f1);
                                m_a2pfS0Coeffs[1][nRef] = f0/(f0+f1);
                                m_a2pnS0Points[0][nRef] = m_nNumS0Points + nLargest1;
                                m_a2pnS0Points[1][nRef] = m_nNumS0Points + nLargest2;
                                m_poS0Points[m_nNumS0Points + nLargest1].nContrib++;
                                m_poS0Points[m_nNumS0Points + nLargest2].nContrib++;
                                if (m_a2pfS0Coeffs[0][nRef]>m_a2pfS0Coeffs[1][nRef])
                                    m_pnClosestS0[nRef] = m_a2pnS0Points[0][nRef];
                                else
                                    m_pnClosestS0[nRef] = m_a2pnS0Points[1][nRef];
                            };
                        };

                        // Find the rotation matrix that
                        // transforms the S0 vector so that it points down <0,0,1>
                        // This is composed of a matrix taking a3fRot to a3fStartRot, followed
                        // by a3x3fOffsetRot which takes it to <0,0,1>
                        f0 = 5.0;
                        for (nx=0;nx<3;nx++) {
                            vConvRotVec3DMat3D(fPerpVecSign*(fStartRot-m_pfOrientAngle[nRef]),&a3fPerpS0[0],&a3x3fMat[0][0]);
                            vMulMat3DMat3D(a3x3fOffsetRot,a3x3fMat,a3x3fRotMat);
                            
                            // Make sure that the S0 vector is getting transformed correctly.
                            vMulMat3DVec3D(a3x3fRotMat,a3fRot,a3fVec1);
                            f1 = fabs(a3fVec1[0]) + fabs(a3fVec1[1]) + fabs(a3fVec1[2]-1.0);

                            // This complicated logic results in the smallest vector getting chossen.
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


                        // Transform the S vector accordingly.  It should now lie in the correct hemisphere. (i.e. [2]>=0.0)
                        vMulMat3DVec3D(a3x3fRotMat,a3fScatt,a3fVec1);

                        // Save the transformed S vector.
                        vCopyVec3D(a3fVec1,m_pfSVectors + nRef*3);
                        // Save the transformed S vector angle (important if we use harmonics)
                        m_pfChi[nRef] = fGetVectorAngle(a3fVec1);
                    };
                };                

                if (0) {
                    // Debug:  Print out points.
                    for (nx=0;nx<nNumCoeffs;nx++) 
                        m_poS0Points[m_nNumS0Points + nx].nPrint();
                };
                if (fabs(fMaxDot - fMinDot)>0.01) {
                    printf("WARNING:  Orientation %s??? has inconsistent S0 vectors!\n",sNextBatchOrient.string());
                    nNumWithSeriousProblems++;
                    break;
                };
                m_nNumS0Points += nNumCoeffs;


            } else
            break;
            
            m_nNumOrients++;
        } while (1);
        printf("--------------------------------------------------------------------------\n");
        if (nNumWithSeriousProblems) {
            printf("ERROR: %d bad reflections detected!  Algorithm terminating.\n",nNumWithSeriousProblems);
            goto exit_place;
        };
        if (nNumWithProblems) {
            printf("WARNING: %d reflections had suspicious S0/S vectors!\n",nNumWithProblems);
        };

    };

    nStat = 0;
exit_place:
    return nStat;
};

