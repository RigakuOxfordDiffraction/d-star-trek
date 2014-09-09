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
// CEyeballIndex.cc            Initial author: Thaddeus Niemeyer 24-August-01
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

#include "Dtrek.h"
#include "Creflnlist.h"
#include "Cimage_header.h"

class DTREK_EXPORT CEyeballIndex {
public:
    enum            {m_nBest3DVectors = 10};
    enum            {m_nBest2DVectors = 10};
    enum            {m_nRefine2 = 1,m_nRefine3 = 2};
    int             m_nMaxBest3DVectors;      // Maximum number of 3D vectors (by default equals m_nBest3DVectors)
    int             m_nMaxBest2DVectors;      // Maximum number of 2D plane directions (by default equals m_nBest2DVectors)

private:
    double          m_f3DAngleStep;           // Step used when finding 3D vectors
    int             m_n3DGridSpacing;         // Number of steps used for 3D grid method (actually it's a 2D grid)
    int             m_n3DGridSpacingRefine;   // Number of steps used for 3D grid method (refinement step).
    double          m_f3DAngleStepRefineMax;  // Initial step for refining 3D vectors.
    double          m_f3DAngleStepRefineMin;  // Minimum step for vector refinement.
    double          m_fMax3DAngleBetweenSolns;// Maximum angle between candidate 3D solutions.


    double          m_f2DAngleStep;           // Step used when finding 2D vectors.
    double          m_fMax2DAngleBetweenSolns;// Maximum angle between candidate 2D solutions.
    int             m_n2DGridSpacing;         // Number of bins used in 2D array (it's 1D)
    int             m_n2DGridSpacingRefine;   // Number of bins used in 2D array (it's 1D)
    double          m_f2DAngleStepRefineMax;  // Initial step for refining 3D vectors.
    double          m_f2DAngleStepRefineMin;  // Minimum step for vector refinement.

    
    double          m_fMaxPruneReject;        // When doing a second step (already having some vector), we like to reject the worst data.
    int             m_nMinRefsAfterPrune;     // Mininum number of reflections that should be left after pruning.

    Creflnlist&     m_oList;                  // Input reflection list.
    Cimage_header&  m_oHeader;                // Input header file.
    int*            m_pn3DGridSums;           // Array used in 3D grid sums
    int             m_n3DGridSize;            // Size of 3D Grid.
    double          m_a3fMinX[3];
    double          m_a3fMaxX[3];
    double          m_fMaxX;
    bool            m_bIsAvailable;
    int             m_nNumRefs;

    int*            m_pnSort;                       // Used to sort packed HKL for printing.
    double*         m_pfRefineErrorPercent;         // Used when refining edges.
    double*         m_pfRefineErrorPercentDiffs;    // Used when refining edges (difference vectors).
    double*         m_pfPrintResidBuffer;           // Used when printing.
    int*            m_pnPrintNeighbors;             // Used when printing.

    double*         m_pfFFTArray;             // Used when compute the FFT array.
    int             m_nFFTArraySize;          // Size of array (last stored).

    int             m_nNumBest3DVectors;
    int             m_aanBest3DVectorCounts[m_nBest3DVectors];
    double          m_aa3fBest3DVecs[m_nBest3DVectors][3];
    double          m_aa2fBest3DVecsAngles[m_nBest3DVectors][2];    

    int             m_nNumBest2DVectors;
    double          m_aa3fBest2DVecs[m_nBest2DVectors][3][3];
    double          m_aafBest2DVecsDSpacing[m_nBest2DVectors];
    double          m_aafBest2DVecsRealSpacing[m_nBest2DVectors];
    int             m_aanBest2DVectorCounts[m_nBest2DVectors];
    int             m_aanBest2DVectorSortOrder[m_nBest2DVectors];
    
    int             m_nFI_nIndexSelect;       // Used for I/sig and resolution rejections.
    double          m_fMinIoverSig;           // Set after calling nRejectIoverSigReso()
    double          m_fMinReso;               // Set after calling nRejectIoverSigReso()
    double          m_fMaxReso;               // Set after calling nRejectIoverSigReso()

    int             m_a3nOdd[3];              // Number of odd refined (fRefineSolution())
    int             m_a3nEven[3];             // Number of even refined (fRefineSolution())
    double          m_a3fOdd[3];              // Number of odd refined (I/sig) (fRefineSolution())
    double          m_a3fEven[3];             // Number of even refined (I/sig) (fRefineSolution())
    int             m_a2nRefined[2];          // [0] == Total considered. [1] == Total passed. (fRefineSolution())
    double          m_fLastRefined;           // Last refined value (fRefineSolution())
    double          m_fTopErrPercent;         // Maximum error percent for refinement (fRefineSolution())

    int             m_nFI_nNeighborFirst;     // When nOverlapCheck() is called to get difference vectors, this is used.
    int             m_nFI_nNeighborNext;      // When nOverlapCheck() is called to get difference vectors, this is used.
    Creflnlist*     m_poDiffList;             // Difference vector list (loaded if nComputeDiffVecs() is called).

    int nCalcHeuristic2D(double a3fProjVector[3],int nGridDim,bool bPrune = false);
    int nCalcHeuristic1D(double a3fCompressVector[3],int nGridDim,bool bPrune = false);
    void vPercentComplete(int nStatus=-1,char* pcMessage = NULL);
    bool bFindLC(double* pfErrPercent,double* pfBasisVecs,int nNumVecs,double* pfVec,double* pfVecD,int* pnCoeffs,int nNormalize = 1);
    double fDPSFFT(double a3fVecT[3]);
    double fMinVectorDifference(double a3fVec1[3],double a3fVec2[3]);
    void vInit();
    void vCalcGetMaxCellAndWavelength();


public:

    bool bIsAvailable() { return m_bIsAvailable; };

    int             m_nPrint;                 // Used to do diagnostics.
    int             m_nVerbose;               // Verbosity.
    int             m_nTwinID;                // Which twin componant to use (default: -1==ignore).
    double          m_fWavelength;            // Wavelength;
    double          m_fCellLengthMax;         // Maximum possible cell length (calculated from detector params if available).
    
    int nRejectIoverSigReso(double fMinIoverSig,double fResoMin,double fResoMax);
    int nGetBestViewDirections3D(int nBufSize,double* pfVecs);      // Loads the single perpendicular 'view' vector.
    int nGetBestPlane2D(int nSolution,double* pfPlaneMat,double* pfDSpacing,double* pfRealSpacing);
    int nCalcRefineBestViewDirections3D();
    int nRefineBestViewDirections(double* pfVector,int nMethod,double* pfRestriction = NULL);
    int nCalcRefineBestViewDirections2D(double* pfPerpVector,bool bKeepPrune = false);
    bool bSameOrientMatrix(Ccrystal& oCrystal1,double fErrorPercent);

    int nComputeDiffVecs();
    int nAutoIndex(double a3x3fNonReduced[3][3]);   // Auto index. Be sure to set m_nTwinID.
    int nComputeDoubledAxes(double fErrPercent,Ccrystal& oCrystalIn,Ccrystal& oCrystalOut);
	int nAddTweaks();
	int nSearchForTwinLaw(Ccrystal& oCrystalIn,Ccrystal& oCrystalOut,int nMaxIndex,double fPercentReject,double fMinTwinFraction);

    int    nDoubledAxes();
    double fRefineSolution(double fErrPercent,double a3x3fOrientMatrixIn[3][3],double a3x3fOrientMatrixOut[3][3],int* pnPass = NULL);
    double fRefineSolution(double fRejectPercent,int nNumEdges,double a3x3fEdgesIn[3][3],double a3x3fEdgesOut[3][3]);
    double fRefineSolution(double fRejectPercent,
                           Ccrystal& oCrystalIn,
                           Ccrystal& oCrystalOut,
                           int* pnNumberIndexedReflns=NULL);
        
    bool bRefineBeam(double fRejectFraction,
                     int nBeamSearchRad,
                     int nMaxBeamMovement,
                     bool bSetNewCenter,
                     Ccrystal& oCrystalIn,
                     Ccrystal& oCrystalOut, 
                     double *pdBeamPixChange,
                     double& fBestValue,
                     Cstring& sError);

  

    // Use this routine to move from plane description to edge description.

    int nConvertPlanes2OrientMatrix(
        double* pfNorm1,double fSpace1,
        double* pfNorm2,double fSpace2,
        double* pfNorm3,double fSpace3,
        double a3x3fOrientMat[3][3]);


    // Use this routine only for printing.

    int nCalcPrintableEdges(double* pfIndexEdges,
                        int nNumIndexEdges,
                        double fMaxResidual,
                        int nMaxEdgesToGenerate,
                        double* pfEdgesOut0,
                        double* pfEdgesOut1,
                        int nTwinID = -1);
    int nCalcFitPlane(int nPoints, double aa3fPointsSelectedInPlate[][3], double a3fPlaneNormal[3],double a3fPlaneCentroid[3],double* pfInterPlaneDistance);

    CEyeballIndex(Cimage_header& oHeader,Creflnlist& oList);
    ~CEyeballIndex();
};

