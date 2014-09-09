#ifndef DT_CINDEX_H
#define DT_CINDEX_H
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
// Cindex.h        Initial author: J.W. Pflugrath           26-May-1995
//    This file is the header file for class Cindex
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

//+Include files

#include "Dtrek.h"
#include "Cimage_header.h"
#include "Cdetector.h"
#include "Csource.h"
#include "Ccrystal.h"
#include "Crotation.h"
#include "Creflnlist.h"
#include "C3Ddata.h"
#include "Cimage.h"
#include "dtrekvec.h"
#include "CCellReduce.h"

//+Definitions and constants

//+Code begin

const int DT_NUMSOL_MIN_INIT = 10;
const int DT_NUMSOL_MAX_INIT = 50;


class TestRecord;
class DTREK_EXPORT Cindex {
    friend class TestRecord;
    friend DTREK_EXPORT int dtcell_main(unsigned int argc,char*argv[]);

public:

  float    m_a4fRejectLimits[4];      // 2 mm limits and 2 degree limits
  float    m_fWavelength;
  float    m_a3fS0[3];
  float    m_a3x3fBasisVecs[3][3];
  float    m_a6fInputCell[6];
  float    m_a6fInputErrors[6];
  float    m_a3x3fUBMat[3][3];
  float    m_a3x3fInvUBMat[3][3];

  float    m_a3x3fUBBest[3][3];
  float    m_a3x3fInvUBBest[3][3];

  float    m_a26x6x6fBMat[26][6][6];

  bool     m_bKnownCell;
  float    m_fSubPeak;
  float    m_fCellLengthMax;
  Ccrystal m_oCrystal;   // A place for the solution

  int         m_nVerbose;
  int         m_nNumDetectors;
  int         m_nWhichDetector;
  Cstring    *m_psDetectorNames;

  Cdetector  **m_ppoDetector;  // Pointer to pointer of detector object(s)
  Csource     *m_poSource;     // Pointer to source object
  Cgoniometer *m_poCrysGonio;  // Pointer to crystal goniometer;
  Crotation   *m_poRotation;   // Pointer to rotation object;

  Creflnlist  *m_poReflnlist;  // Pointer to working reflnlist object

  Creflnlist  *m_poDiffVecs;   // Pointer to list of sum and difference vectors
  Creflnlist  *m_poSolutions;  // Pointer to list of answers

  float        m_a44x6fLattCells[45][6];    // These hold the best lattice cells for all 44 lattice 
                                            // characters after a call to the cell reduction.
  float        m_a44fLattResid[45];         // These hold the lattice residuals.
  static int   ms_a44nOrder[];              // Ordering of 44 lattice chars. w.r.t. Int. Tables.

  int          m_nIndex0;
  int          m_nIndex1;
  int          m_nIndex2;
  int          m_a15x3nBIndices[15][3];

  int         *m_pnIndex;
  int         *m_pnSoln;

  int          m_nSolnPick1;
  int          m_nSolnPick2;
  int          m_nLattNum;
  float        m_fMMSeparationMin;
  float        m_fDiffVecAvg;
  float        m_fDiffVecMax;
  float        m_fMinFFTResidual;
  float        m_fSumFFTResidual;

  int   m_nNumDiffVecs;
  int   m_nNumRowVecs;

private:


  C3Ddata  *m_po3DTrick;
// These are values that get indexed:

  float m_a6fCrysUnitCell[6];  // 6 real unit cell values
  float m_a6fCrysRecipCell[6]; // 6 reciprocal unit cell values
  float m_a3fCrysRot[3];       // 3 crystal rotation (missetting/orient) angles

  float m_a3x3fInvSMat[3][3];

  float m_fAngleEps;
  float m_fLengthEps;

  float m_fFourierMax;
  float m_afCos2PI[1025];

  float m_fMinAllowedResidual;
  float m_fMinAllowedResidualStart; // = 0.95;
  float m_fMinAllowedResidualEnd;   // = 0.50;
  float m_fMinAllowedResidualDelta; // = 0.01;
  float        m_fResolutionMin;    // Min resolution limit in A (99999 is min)
  float        m_fResolutionMax;    // Max resolution limit in A (0.01 A is max)


// Some other things

  bool m_bNewDetNames;
  bool m_bNewDetector;
  bool m_bNewSource;
  bool m_bNewCrysGonio;
  bool m_bNewReflnlist;
  bool m_bNewRotation;
  bool m_bReIndex;

  static Cstring ms_sfA;
  static Cstring ms_sfB;
  static Cstring ms_sfC;
  static Cstring ms_sfAlpha;
  static Cstring ms_sfBeta;
  static Cstring ms_sfGamma;
  static Cstring ms_sfVolume;
  static Cstring ms_sfRot1;
  static Cstring ms_sfRot2;
  static Cstring ms_sfRot3;
  static Cstring ms_sfDiffLength;
  static Cstring ms_ssLattice;
  static Cstring ms_ssLattSymm;
  static Cstring ms_snLattNum;
  static Cstring ms_snSortOrder;
  static Cstring ms_snIndexSelect;
  static Cstring ms_snTwinID;
  static Cstring ms_sfSortItem;
  static Cstring ms_snNumIndexed;

  tagLattChar m_a44tLattChar[45];

//

public:

// And the field indexes for each of the above:

  int    m_nFI_fA;
  int    m_nFI_fB;
  int    m_nFI_fC;
  int    m_nFI_fAlpha;
  int    m_nFI_fBeta;
  int    m_nFI_fGamma;
  int    m_nFI_fVolume;
  int    m_nFI_fRot1;
  int    m_nFI_fRot2;
  int    m_nFI_fRot3;
  int    m_nFI_sLattice;
  int    m_nFI_sLattSymm;
  int    m_nFI_nLattNum;
  int    m_nFI_nSortOrder;
  int    m_nFI_fDiffLength;
  int    m_nFI_fSortItem;
  int    m_nFI_nIndexSelect;
  int    m_nFI_nTwinID;

  int    m_nFI_nBins;
  int    m_nFI_fPhi;
  int    m_nFI_fPsi;
  int    m_nFI_fTvec0;
  int    m_nFI_fTvec1;
  int    m_nFI_fTvec2;
  int    m_nFI_fDelta;

  int    m_nFI_nNumIndexed;

  static Cstring ms_sDtindexOptions;

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cindex ();       // Construct an empty index object
                 // In the future this should not be allowed to happen.

Cindex (Cimage_header& oHeader);

~Cindex ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:


void vSetReflectionList(Creflnlist* poList);
void vSetResolution(const float fResoMin, const float fResoMax);
int nInitValues(void);
int nInitValues(Cimage_header& oHeader);
int nList(const int nFlag);
int nCalcGetMaxCell(float *pfMaxCell);

// Below should be private after testing

int nCopyReflns(void);
int nCalcSumDiffs(bool bDoDiffVecs,const int nMaxNumRef=300, const float fCellLengthMax=0.0,int nTwinID = 0);
int nGroupVecs(const float fTest = 0.02);
int nFindIndices(const float fCellMin, const float fCellMax,
                 const int nIndexAbsMin, const int nIndexAbsMax,
                 const float fResCutoff = 0.90);
float fIntegerResidual(const float fVec[3], const int nTopRefs = 0);
int nLoadSolution(const int nSolnNum, Ccrystal *poCrystalIn=NULL);
int nPutSolution(const int nSolnNum, Ccrystal *poCrystalInOut,double fResidual);
int nPruneSolutions(const int nNumToKeep);

int nSetup(void);
int nIndex(void);


int nReduceCell(float *pfCell, float *pfRot=NULL);

int nIndexCos(const float fMinRadius=3.0, const float fMaxRadius=200.0,
              const int nGrid=100, const int nDiffs=250, float *pfAvg=NULL);
int nIndexCos(const float *pfMinAngstroms, const float *pfMaxAngstroms,
              const int nGrid=32, const int nDiffs=250, const int nFlag=0,
              float *pfCOG=NULL);

void vInitSolutions(void);
int  nInsertSolution(const float *pfCell, const float *pfOrient = NULL,
                     const int nSpacegroup = 0);

int nLoadSolutionArray(void);
int nGetPeaks(C3Ddata *po3Ddata, const float fStep, float *pfPeakMin=NULL);
int nListResults(const int nSoln=-1, int *pnBest=NULL);
int nListResultsB(const int nSoln=-1);

int nSelectMatch(const int nSpacegroupNum);
void vSetKnownCell(const float *pfCell);
void vSetKnownErrors(const float *pfErrors);
float fCalcGetGrid(void);
double fGetKnownCellResidual(double a6fKnownCell[6],Ccrystal& oCrystal,double fPromptingResidual = -100.0);

inline void vSetVerboseLevel(const int nLev) { m_nVerbose = nLev; }
inline int  nGetVerboseLevel(void) { return (m_nVerbose); }
inline void vSetCellLengthMax(const float fValue) { m_fCellLengthMax = fValue; }
inline float fGetCellLengthMax(void) { return (m_fCellLengthMax); }
inline float fGetDiffVecAvg(void) { return (m_fDiffVecAvg); }
inline float fGetDiffVecMax(void) { return (m_fDiffVecMax); }
inline void vSetMMSeparationMin(const float fMin) { m_fMMSeparationMin = fMin; }


inline void vSetSubPeak(const float fSubPeak)
            { m_fSubPeak = fSubPeak; }
int nIndexDPS(const int nNumRefsUsed,bool bCoarse);

int nDPSFFT(const float a3fVecT[3], 
            const int nRow, 
            const int nNumRefsUsed, 
            int* pnDataSize, 
            float** ppfData,
            Crefln* poReflnAnswer, 
            Cimage* poImgHist, 
            Cimage* poImgFFT);

int nDPSVecRefine(const float fPhiInc, const float fPsiInc,
                  const int nRow, 
                  const int nNumRefsUsed, int *pnDataSize, float **ppfData,
                  Crefln *poReflnAnswer, Cimage *poImgHist, Cimage *poImgFFT);

int nCalcReflnResolution(Creflnlist *poReflnlist);
 int nDeleteReflns(const double fSigmaCutoff, const double fSharpnessCutoff);
int nReduceCell(Ccrystal& oInputCell,Ccrystal* poOutputCrystals,double fResid = 5.0);

#ifdef SSI_PC
    //void BuildIndexChooseSolutionTclCommand(Cstring &csCommand);
void vBuildTableCmd(Cstring &csCommand);
#endif

//private: // Some of the above functions will be made private

int nExpandReflnlist(void);

public:
    void vSetNumIndexedBySolution(int nSolution, int nIndexed);
    void vWriteRankingInfoToHeader(Cimage_header* poHdr, int nNumUsed);

///////////////////////////////////////////////////////////////////////////////////////////////
// RB: 5/19/04  The member variables below used to be static variables in the member functions.
// They would only get initialized when an object was created for the first time. So we need to
// turn them into initializable class members.
private:
  int       m_nAtom_DPSFFT;
  bool      m_bFirstTime_IndexCos;
  bool      m_bFirstTime_IndexDPS;
  int       m_nJ_IndexCos;
  int       m_nDepth_DPSVecRefine;
///////////////////////////////////////////////////////////////////////////////////////////////
};  // end of class Cindex

#endif   // DT_CINDEX_H
