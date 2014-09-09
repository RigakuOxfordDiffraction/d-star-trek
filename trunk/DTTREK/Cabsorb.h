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

#ifndef CABSORB_INCLUDE
#define CABSORB_INCLUDE
#include "Dtrek.h"
#include "dtrekvec.h"
#include "Cstring.h"
#include "Creflnlist.h"
#include "dtsvd.h"
#include "Cstat.h"
#include "Cscaleaverage.h"

enum { NEIGHBORVECS = 16 };

struct CSpoint {
    double a3fCoords[3];
    int anAdj[NEIGHBORVECS];
    int nAdjCt;
    int nRow;           // The points go around the sphere in rows.  Row 0 (south pole) has 1 collumn in it.
    int nCol;           // The points in each row are enumerated.  The equatorial row has the most columns.
    int nContrib;       // Number of contributors in the reflection list.

    double fChi;        // This variable changes from 0 to 2Pi as you traverse a row.
    double f2Theta;     // This variable changes as the angle increases.
    int nIndex;         // Index w.r.t. the whole array of points.  All entries in anAdj[] use these numbers.
    
    int nGetAdj(int iIndex) { return anAdj[iIndex]; };
    int nAddAdj(int nTo);
    bool bIsAdj(int nVert);
    int nRemoveAdj(int nTo);

    CSpoint();
    ~CSpoint();
};

struct CS0point {
    double a3fCoords[3];
    double fRot;
    double fDispRot;
    int nOrient;
    int nPoint;
    int nCumulPoint;
    int nContrib;

    int nPrint();
};

#ifdef CONIO
#include <conio.h>
#endif

#ifdef CONIO
#define EXIT(x) { getch(); exit(x); }
#else
#define EXIT(x) exit(x)
#endif

enum {MODE_488,MODE_4888,MODE_481616,TOTAL_MODES};

struct Cabsorbsurface {
    
    CSpoint* m_poSPoints;
    CS0point* m_poS0Points;

    int* m_pnRingStart;       // The points are ordered on rings.  A ring is a block of reflections with the same nRow value.
    int* m_pnRingEnd;
    int  m_nRings;
    int* m_a3pnTriangles[3];  // Triangles.  Filled in by nBuildTriangles()
    int m_nNumTriangles;


    int*    m_a3pnSPoints[3];  // For each reflection this gives the points defining an enclosing triangle of the S vector.
    double* m_a3pfSCoeffs[3];  // For each reflection this gives the coefficient for the above triangular entries.
    int*    m_a2pnS0Points[2]; // For each reflection this gives the 2 points defining the S0 vector
    double* m_a2pfS0Coeffs[2]; // For each reflection this gives the coefficient for the above entries.
    int*    m_pnOrientation;   // Orientation for each reflection.
    int*    m_a2pnClosestS[2]; // Closest of m_a3pnSPoints for each reflection.  One for upper, and one for lower row.
    int*    m_pnClosestS0;     // Closest of m_a3pnS0Points for each reflection.
    double* m_pfOrientAngle;   // Angle for each reflection (we might merge orrientations).
    int m_nNumRefs;            // Number of reflections  (computed in nLoad())
    int m_nNumOrients;         // Number of orientations (computed in nLoad())
    float* m_pfSVectors;       // S vectors pre refln.   (computed in nLoad())
    float* m_pfChi;            // Chi angle.             (computed in nLoad())
    int* m_pnC;                // Counts per segment     (computed in nCheckMeshDimensions())

       
    int nPrint(int nMode,Cstring& sFilename,char* cpLabel,bool bPrintPoints); 


    int m_nNumSPoints;
    int m_nMax2Theta;
    int m_nDegreesPerS0;
    int m_nNumS0Points;
    double m_fMinNonHarmonicRings;
    double m_fMinRingContribAverage;
    double m_fMinContrib;    

    
    int nBuildPoints(int* pnPoints,double* pfChiOffsets);
    int nBuildPoints(char* cpTitle,int nMode); // Builds the absorption point coordinates (Set of normalized vectors).  Orients based on quadratic form vectors.
    int nBuildTriangles();
    int nLoad(Creflnlist& oReflnlist,int a3nS[3],int nBatchIdx,int a3nS0[3],tagBatch* poBatches);                     // Loads coefficients for each reflection.

    int  nCheckMeshDimensions(Creflnlist& oList,int* pnIndex,int* pnReject,int nMinRedundancy);
    int  nFindTriangleCoeffs(Creflnlist& oReflnlist);
    bool bGetPointCoeffs(double a3x3fPoints[3][3],double a3fPoint[3],double a3fCoeffs[3]);
    int  nGetRingStartEnd(int nRing,int& nRingStart,int& nRingEnd);
    double fGetVectorAngle(double* pfVector);


    Cabsorbsurface();
   ~Cabsorbsurface();
};


class Cabsorbmodel {
    double m_fDegreesPerS0;
    Creflnlist& m_oListIn;
    Cabsorbsurface m_oSurface;
    int an2ThetaPoints[180];

    int m_a3nFI_fS0Vec[3];
    int m_a3nFI_fSVec[3];
    int m_nFI_nPackedHKL;
    int m_nFI_sBatch;
    int m_nFI_fScale;
    int m_nFI_nBatchIdx;
    int m_nFI_fSigmaI;
    int m_nStat;

    int* m_pnIndex;             // User provided sort of HKL entries.
    int* m_pnReject;            // Reflections that get rejected.
    double* m_pfScale;          // Scale factor for the reflections.
    double* m_pfSS0Coeffs;      // Coefficients (faster dimension encodes S points).

    int nInitReflectionList(Creflnlist& oList);
    double fGetScale(int nRef);
    int nRelativeRingScale();
    int nPropogateRingScale();
    int nPrintRingScale(bool bPrintCounts);

    // These are loaded after calling fComputeRMerge();
    int m_nContributors;

public:
    bool bIsAvailable() { return (m_nStat == 0); };


    int    m_nMinRedundancy;                    // Minimum redundancy for the R-merge.
    int    m_nMaxFoundRedundancy;               // Loaded after calling fComputeRMerge();
    int    m_nDegreesPerS0;
    double m_fMaxCoeffRange;
    double m_fMinCoeffRange;

    double fComputeRMerge();
    int nBuildReflectionScaleFactors();   
    int nLinearMethod(int nResoRing = -1);
    int nApply();
    int nInit(Cscaleaverage& oScale);

    Cabsorbmodel(Creflnlist& oListIn,int* pnHKLSort,int* pnReject);
    ~Cabsorbmodel();
};
#endif
