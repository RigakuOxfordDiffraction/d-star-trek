#ifndef DT_CREFINE_H
#define DT_CREFINE_H
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
// Crefine.h        Initial author: J.W. Pflugrath           12-May-1995
//    This file is the header file for class Crefine
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
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
#include "CXprop.h"
#endif
#include "dtrekvec.h"
#include "dtsvd.h"

//+Definitions and constants

const int nMAXPARAMS = 43;

//+Code begin

struct CrefineTwinData;

class DTREK_EXPORT Crefine {

  enum {nMaxTwins = 10};
  enum eTransTypes {
      DET_TRANS1,    //   0    Det translation 1
      DET_TRANS2,    //   1    Det translation 2
      DET_TRANS3,    //   2    Det translation 3 (aka crystal to detector distance)
      DET_ROT1,      //   3    Det rotation 1
      DET_ROT2,      //   4    Det rotation 2
      DET_ROT3,      //   5    Det rotation 3
      GONIO_MISS1,   //   6    Goniometer miss-setting 1
      GONIO_MISS2,   //   7    Goniometer miss-setting 2
      SOURCE_ROT1,   //   8    Source rotation 1
      SOURCE_ROT2,   //   9    Source rotation 2
      SOURCE_WAVE,   //   10   Source wavelength
      CRYS_ROT1,     //   11   Crystal rotation 1
      CRYS_ROT2,     //   12   Crystal rotation 1
      CRYS_ROT3,     //   13   Crystal rotation 1
      CELL_A,        //   14   a*
      CELL_B,        //   15   b*
      CELL_C,        //   16   c*
      CELL_ALP,      //   17   alpha*
      CELL_BET,      //   18   beta*
      CELL_GAM,      //   19   gamma*
	  CELL_MOSAICITY,//   20
      CELL_CONS,     //   21   Refining with constraints?
	  CELL_RECIP00,	 //   22   Reciprocal shift vector for 0th twin law component.  (identically the UB matrix itself without twin ops).
	  CELL_RECIP01,	 //   23
	  CELL_RECIP02,	 //   24
	  CELL_RECIP10,	 //   25   Recirpocal shift vector for 1st twin law component.
	  CELL_RECIP11,	 //   26
	  CELL_RECIP12,	 //   27
      CELL_RECIP20,  //   28   Recirpocal shift vector for 2nd twin law component.
      CELL_RECIP21,  //   29
      CELL_RECIP22,  //   30
      CELL_RECIP30,  //   31   Recirpocal shift vector for 3rd twin law component.
      CELL_RECIP31,  //   32
      CELL_RECIP32,  //   33
	  CELL_TWINLAW10,//   34   Twin law for 1st twin law component  
	  CELL_TWINLAW11,//   35
	  CELL_TWINLAW12,//   36
      CELL_TWINLAW20,//   37   Twin law for 2nd twin law component  
      CELL_TWINLAW21,//   38   
      CELL_TWINLAW22,//   39
      CELL_TWINLAW30,//   40   Twin law for 3rd twin law component  
      CELL_TWINLAW31,//   41
      CELL_TWINLAW32,//   42
      CELL_RECIP_MAX        =CELL_RECIP32,
      CELL_TWINLAW_MAX      =CELL_TWINLAW32,
      ARRAYFLAG = 1000
  };
  enum eWeightTypes {
      eWeightRotSigma,
      eWeightIntensity,
      eWeightIoverSig,
      eWeightLorentz,
      eWeightRotWidth,
      eWeightNone,
      eWeightLAST
  };
public:

  enum eRefinementReflectionNumbers {
      eReflectionNumberTotal,
      eReflectionNumberAccept,
      eReflectionNumberRejects,
      eReflectionNumberIgnored,
      eReflectionNumberDegRejects,
      eReflectionNumberMMRejects,
      eReflectionNumberOtherComp};



  bool   m_bIsAvailable;
  bool   m_bPrintTables;                        // Do we print data tables at the very end of refinement.
  
  bool  m_bUseOverlapped;                       // Do not reject overlapped different HKLs, when getting reflections from images

  int    m_nDisplay;                            // Flag to update display of calling process
  float  m_a4fRejectLimits[4];                  // 2 mm limits and 2 degree limits
  int    m_a7nReflnNumbers[7];                  // Total, Accept, Rejects, Ignored, DegRejects MMRejects Other-comp
  float  m_a7fReflnNumbers[7];                  // Same, but containing total photon counts.
  float  m_fRMSDEG;
  float  m_fRMSMM;
  float  m_fRMSANG;                             // r.m.s difference between obs & calc 
  float  m_a4fAvgDistanceToNearest[4];          // Average distance to nearest spot. [0],[1] = values [2],[3] = counts
  float  m_fFindPad;                            // Padding to use during sequence refinements.
  float  m_a2fResoUsed[2];                      // Min/Max resolution limits used.
  float  m_a3x10fIndexHKLErrorRefs[3][10];      // Tables describing the indexing residual in each 
  float  m_a3x10fIndexHKLErrorIntensity[3][10];
  Cstring m_sRejectFile;						// Rejection file name (output after nGetReflnsFromFind()).
  Cstring m_sMapFile;                           // Map file for sequence finds.
  Cstring m_sLogRefine;                         // Refinement Log.  
  Cstring m_sLogRefineLast;                     // Last Refinement Log line.
  Cstring m_sLogTwinReject;                     // Important information for twin rejection.
  CrefineTwinData* m_poComp;                    // Important information for each component.
                                        // spot angles and millimeters
  float  m_fSigma;                      // Sigma cutoff.
  float  m_fSharpness;                  // Sharpness cutoff.
  float  m_fLambda;                     // Used in nOpLevenbergMarquardt()
  float  m_fTweakSensitivity;           // Tweak sensitivity.  10=large.  0.1 = small
  float  m_fMaxCorrCoeff;               // Maximum correlation coeff. between two variables for rejection.
  float  m_fAvergeDiagVal;              // Used to print out the "newtonness" of the iteration
  float  m_fResidMin;                   // Stop iterating if we are below this residual.
  float  m_fLorentzMax;                 // Maximum Lorentz factor.
  float* m_pfResidProfile;              // Loaded with the residual profile for all tests.
  float* m_pfLambdaProfile;             // Loaded with the lambda profile for all test.
  int    m_nProfileCt;                  // Size of the above.
  int    m_anRefineFlags[nMAXPARAMS]; // an extra one for crystal mosaicity
  int    m_anSigmaFlags[nMAXPARAMS];  // Which parameters did we refine sigma on?
  int    m_nNumTwins;                   // Number of twins available.
  int    m_nMaxTwinLaws;				// Maximum # of twin laws for any given component.
  int    m_nPad;                        // Padding to use for finding spots.
  int    m_nTwinID;
  itr<int> m_anReflnOverlapMat;             // Overlap data for reflections used in find list.
  itr<double> m_afReflnOverlapIntensities;  // Overlap intensities for each reflection in the find list
  itr<double> m_afReflnMinResid;			// Arrays containing the minimal residual for that reflection.	
  itr<int>    m_anReflnBestTwin;			// Corresponding twin associated with best residual.

  eWeightTypes m_eWeightType;           // Weighting type.
  int          m_nAutoCorrCheck;        // Should we run the auto-correlation check when computing sigma values?
  int          m_nVerbose;
  int          m_nNumDetectors;
  int          m_nCycles;
  int          m_nWhichDetector;
  int          m_nAdjustRotation;       // Should we adjust the crystal orrientation matrix during mosaicity refinement?
  Cstring     *m_psDetectorNames;
  float        m_fResolutionMin;        // Min resolution limit in A (99999 is min)
  float        m_fResolutionMax;        // Max resolution limit in A (0.5 A is max)

  float        m_fTestMosaicityRangeMax;
  float        m_fTestMosaicityRangeMin;
  float        m_fLowestSig;        // Place a cutoff on the lowest sigma weight used for refinement (calculated by program).

  Cdetector  **m_ppoDetector;       // Pointer to pointer of detector object(s)
  Csource     *m_poSource;          // Pointer to source object
  Ccrystal   **m_ppoCrystals;       // Pointer to crystal object
  Cgoniometer *m_poCrysGonio;       // Pointer to crystal goniometer;
  Crotation   *m_poRotation;        // Pointer to rotation object;
  Cimage_header *m_poHeader;        // Pointer to header.  This is required if we are doing twin indexing.

  Creflnlist  *m_poReflnlistIn;     // Pointer to input reflnlist object
  Creflnlist  *m_poReflnlist;       // Pointer to working reflnlist object
  Cstring      m_sOutHeader;        // Output header filename;

  Cstring      m_a6sUserDetNames[6];
  Cstring      m_asWeightNames[eWeightLAST];
  float        m_a6x6fCrystalVecs[6][6];     // Basis for degrees of freedom on the crystal.

  double m_aafParamPointVectors[nMAXPARAMS][nMAXPARAMS];
  double m_afParamPointValues[nMAXPARAMS];
  double m_afParamPointSum[nMAXPARAMS];
  double m_afParamGoodPoint[nMAXPARAMS];


  static double  ms_pfRefineConvergence[nMAXPARAMS];
  static Cstring ms_sDtrefineFiles;
  static Cstring ms_sDtrefineOptions;
  static Cstring ms_sDtrefineRmsMm;
  static Cstring ms_sDtrefineRmsDeg;
  static Cstring ms_sDtrefineReflectionNumbers;
  static Cstring ms_sDtrefineReflectionNumbersPhotons;
  static Cstring ms_sDtrefineTwinsRefined;
  static Cstring ms_sDtrefineTwinRejected;
  static Cstring ms_sDtrefineLambda;
  static Cstring ms_sDtrefineLambdaProfile;
  static Cstring ms_sDtrefineResidProfile;
  static Cstring ms_sDtrefineProfileCount;
  static Cstring ms_sDtrefineResoUsed;

private:

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  CXprop      *m_poXprop;
#endif

// These are values that get refined:

  float m_a6fCrysUnitCell[6];  // 6 real unit cell values
  float m_a6fCrysTurds[6];
  float m_a6fCrysUnitCellSigs[6];
  float m_a3fCrysRot[3];       // 3 crystal rotation (missetting/orient) angles
  float m_a3fCrysRotSigs[3];
  float m_a3fCrysRotShifts[3];
  float m_a3fCrysMosaicity[3];      // 3 crystal mosaicity.
  float m_a3fCrysMosaicitySig[3];
  float m_a3fCrysMosaicityShift[3];
  float m_aa3fCrysRecipShift[g_nMaxTwinLaws+1][3];
  float m_aa3fCrysRecipShiftSig[g_nMaxTwinLaws+1][3];
  float m_aa3fCrysRecipShiftShift[g_nMaxTwinLaws+1][3];
  float m_aa3fCrysTwinLaw[g_nMaxTwinLaws + 1][3];
  float m_aa3fCrysTwinLawSig[g_nMaxTwinLaws + 1][3];
  float m_aa3fCrysTwinLawShift[g_nMaxTwinLaws + 1][3];
  float m_afCrysTwinFraction[g_nMaxTwinLaws + 1];
  float m_afCrysTwinFractionSig[g_nMaxTwinLaws + 1];
  float m_afCrysTwinFractionShift[g_nMaxTwinLaws + 1];
  float m_a10x3fDetTrans[10][3];  // 3 translation vectors for up to 10 detectors
  float m_a10x3fDetTransSigs[10][3];
  float m_a10x3fDetTransShifts[10][3];
  float m_a10x3fDetRot[10][3];    // 3 rotation angles for up to 10 detectors
  float m_a10x3fDetRotSigs[10][3];
  float m_a10x3fDetRotShifts[10][3];
  float m_a2fSourceRot[2];     // 2 rotation angles for the source
  float m_a2fSourceRotSigs[2];
  float m_a2fSourceRotShifts[2];
  float m_fSourceWave;       // 1 source wavelength
  float m_fSourceWaveSig;
  float m_fSourceWaveShift;

// Some other things

  int   m_nFI0, m_nFI1, m_nFI2, m_nFI3, m_nFI4;    // Field indices for reflection centroids and rotation width
  int   m_nFIMM0,m_nFIMM1,m_nFIMM2;

  bool m_bCheckConvergence;
  bool m_bUseWideSlice;
  bool m_bNotUseWideSlice;
  bool m_bNewDetNames;
  bool m_bNewDetector;
  bool m_bNewSource;
  bool m_bNewCrystal;
  bool m_bNewCrysGonio;
  bool m_bNewRotation;
  bool m_bReIndex;
  bool m_bIsRhomboAsHex;
  bool m_bPrintDetectorMosaicity;     // Print detector position dependent mosaicity.
  bool m_bNoReflectionMerges;         // Do we allow reflection merges to take place?

  Cstring       *m_psVarNames;
  static char*  m_pcVarNames[];            // Variables names.
  static eTransTypes m_aeTransWeighting[]; // Used when rejecting correlated variables.
  static float m_afDerivTweakRef[];        // Used in least squared optimization.
  
  static int ms_nDetTrans1;
  static int ms_nDetTrans2;
  static int ms_nDetTrans3;
  static int ms_nDetRot1;
  static int ms_nDetRot2;
  static int ms_nDetRot3;
  static int ms_nGonioMiss1;
  static int ms_nGonioMiss2;
  static int ms_nSourceRot1;
  static int ms_nSourceRot2;
  static int ms_nSourceWave;
  static int ms_nCrysRot1;
  static int ms_nCrysRot2;
  static int ms_nCrysRot3;
  static int ms_nCrysAstar;
  static int ms_nCrysBstar;
  static int ms_nCrysCstar;
  static int ms_nCrysAlps;
  static int ms_nCrysBets;
  static int ms_nCrysGams;
  static int ms_nCrysMosaicity;
  static int ms_nCrysConstraints;
  static int ms_nCrysRecip00;
  static int ms_nCrysRecip01;
  static int ms_nCrysRecip02;
  static int ms_nCrysRecip10;
  static int ms_nCrysRecip11;
  static int ms_nCrysRecip12;
  static int ms_nCrysRecip20;
  static int ms_nCrysRecip21;
  static int ms_nCrysRecip22;
  static int ms_nCrysRecip30;
  static int ms_nCrysRecip31;
  static int ms_nCrysRecip32;
  static int ms_nCrysTwinLaw10;
  static int ms_nCrysTwinLaw11;
  static int ms_nCrysTwinLaw12;
  static int ms_nCrysTwinLaw20;
  static int ms_nCrysTwinLaw21;
  static int ms_nCrysTwinLaw22;
  static int ms_nCrysTwinLaw30;
  static int ms_nCrysTwinLaw31;
  static int ms_nCrysTwinLaw32;

  static int ms_nFirstCrysParam;
  static int ms_nLastCrysParam;

  static Cstring ms_sDetTrans1;
  static Cstring ms_sDetTrans2;
  static Cstring ms_sDetTrans3;
  static Cstring ms_sDetTrans;
  static Cstring ms_sDetRot1;
  static Cstring ms_sDetRot2;
  static Cstring ms_sDetRot3;
  static Cstring ms_sDetRot;
  static Cstring ms_sDetAll;
  static Cstring ms_sGonioMiss1;
  static Cstring ms_sGonioMiss2;
  static Cstring ms_sGonioMiss;
  static Cstring ms_sSourceRot1;
  static Cstring ms_sSourceRot2;
  static Cstring ms_sSourceRot;
  static Cstring ms_sSourceWave;
  static Cstring ms_sCrysRot1;
  static Cstring ms_sCrysRot2;
  static Cstring ms_sCrysRot3;
  static Cstring ms_sCrysRot;
  static Cstring ms_sCrysAstar;
  static Cstring ms_sCrysBstar;
  static Cstring ms_sCrysCstar;
  static Cstring ms_sCrysLengths;
  static Cstring ms_sCrysAlps;
  static Cstring ms_sCrysBets;
  static Cstring ms_sCrysGams;
  static Cstring ms_sCrysAngles;
  static Cstring ms_sCrysCell;
  static Cstring ms_sCrysMosaicity;
  static Cstring ms_sCrysRecipShift;
  static Cstring ms_sCrysRecipShift0;
  static Cstring ms_sCrysRecipShift1;
  static Cstring ms_sCrysRecipShift2;
  static Cstring ms_sCrysRecipShift3;
  static Cstring ms_sCrysTwinLaw;
  static Cstring ms_sCrysTwinLaw1;
  static Cstring ms_sCrysTwinLaw2;
  static Cstring ms_sCrysTwinLaw3;
  static Cstring ms_sCrysAll;
  static Cstring ms_sCrysConstraints;
  static Cstring ms_sAll;
  static Cstring* ms_psFlags[nMAXPARAMS];

  int       m_nNumTotalParams;                             // NPT Total poss. num. of refined params
  int       m_nNumRefinedCellParams;                       // NP
  int       m_nUserRefineFlags[nMAXPARAMS];                // extra one for crystal mosaicity plus nMAXOPTIONS extra for refinement technique.
  int       m_nUserRefineCrysConstraints;                  // Are we refining the crystal system with constraints imposed?
  Cstring   m_sUserRefineNonCrysFlags;                     // Flags for non-crystal related parameters.
  Cstring   m_sUserRefineCrysFlags;                        // Flags for crystal related parameters.
  
  int m_nNumberOfResoBins;

//

public:

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

    Crefine ();       // Construct an empty refine object
                       // In the future this should not be allowed to happen.

    Crefine (Csource *poSourceIn, Cdetector *poDetectorIn,
	     Ccrystal *poCrystalIn, Cgoniometer *poCrysGonioIn,
	     Crotation *poRotationIn, Creflnlist *poReflnlistIn=NULL);

    Crefine (Cimage_header& oHeader, Creflnlist *poReflnlistIn);

    ~Crefine ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
public:

#if ((defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
    void ReportResults(long Total, long Rejected, long Accepted, long Excluded);
#endif


    int
    nComputeAbMats(double* ppfDeriv[nMAXPARAMS],
                   double* pfVal,
                   float* pfSigs,
                   int   anRefineFlags[nMAXPARAMS],
                   int* pnRejectFlag,
                   int nNumParams,
                   int nNumRefs,
                   double* afA,
                   double* afB,
                   int* pnNumParamsInc);

    //void vUsed2RealVarNum(int x, int nNumParams, int *py);
    void vUsed2RealVarNum(int x, int *py);
    int
    nOpLevenbergMarquardt(bool bComputeSigmas,
                                   int nNumParams,
                                   int nNumRefs,
                                   double *ppfDeriv[nMAXPARAMS],
                                   double *pfVal,
                                   float *pfSigs,
                                   int   *pnRejectFlag,
                                   float afDelta[nMAXPARAMS],
                                   float afSigma[nMAXPARAMS]
                                   );
    int nSimplexInit(int nNumParams);
    int nSimplexRun(double fValue,double afDelta[nMAXPARAMS]);


    int nForceCell(double a6fCell[6]);   // Forces cell to obey spacegroup specified in the m_poCrystal object.
    int nUpdateDisplay();               // Display output stuff.  Removed from refiment loops.
    int nInitValues(void);
    int nInitValues(Cimage_header& oHeader);
    int nList(const int nFlag = 1);
    int nUpdateHeader(Cimage_header *poHeader, int nTwinID,const Cstring &sPre="");
    inline void vSetVerboseLevel(const int nLev) { m_nVerbose = nLev; }
    int nExpandReflnlist(Creflnlist *poReflnlistIn = NULL);
    int nCopyReflns(Creflnlist *poReflnlistOther = NULL);
    int nPrintTables(int* pnReject);
    int nSetup(void);
    int nRefineLoop(const int nLoop = 0);
    double fRefineMosaicityLoop(double fOrigMosaicity,float* pfMosCoeffA,float* pfMosCoeffB,float* pfCalcRotMid,float* pfObsRotMid,float* pfObsRotStart,float* pfObsRotWidth,int* pnReject);
    int nSetRefineFlags(const Cstring& sFlags);
    void vBuildRefineFlags();
    int nSetNoReflectionMerges(const bool bSet) { m_bNoReflectionMerges = bSet;return 0; };
    int nSetRejectReflectionFile(Cstring sFile) { m_sRejectFile = sFile; return 0; };
    int nSetRotOffsetAdjust(bool bSet) { m_nAdjustRotation = bSet; return 0;};
    inline void vSetSigma(const float fSigma) { m_fSigma = fSigma; }
    inline void vSetSharpness(const float fSharpness) { m_fSharpness = fSharpness; };
    void vSetResolution(const float fResoMin, const float fResoMax);
    void vSetRejectLimits(const float fRejLimX, const float fRejLimY,
		          const float fRejLimR);
    int nDoRefinement(const Cstring& sRefineOptions,
		      Creflnlist *poReflnlistOther=NULL);
    int nTwinReject(double fMinFractionIntensityToIndex,double fMaxRMSMM);
    int nDoubledAxes(int nTwinID,double fMinFractionInOdd);

    void vLogResiduals(int nNumIter);
    void vPrintResidualLog(int nPrintCode);

    bool bIsAvailable() { return m_bIsAvailable; };
    int nListResults(void);
    inline Cstring sGetOutHeader(void) { return (m_sOutHeader); }

#define DTREK_REFINE_GRFI_PREDICT_MOSAICITY     0x0001
#define DTREK_REFINE_GRFI_FOR_RANKING           0x0002
    int nGetReflnsFromImages(const int nNumImages, const int *pnSeqNum,
			 Cimage_header *poHeader,
			 const float fResoMin = 0.0f,
			 const float fResoMax = 0.0f,
             DTREK_WORD  wCtrl = 0U);


    int nGetReflnsFromFind(Creflnlist& poReflnlistIn,
			 const float fResoMin, const float fResoMax);
    int nGetMosaicity(const int nNumImages, const int *pnSeqNum,
			     Cimage_header *poHeader,
			     const float fResoMin = 0.0f,
			     const float fResoMax = 0.0f,
                 const int   nGonioSet = 0,
                 const float* pfCalcRotMid = NULL,
                 float* pfRotOffset = NULL,
                 int* pnRejectFlag = NULL,
                 int* pnGonioSet = NULL
                 );

    float  fGetMosaicity()const{return m_a3fCrysMosaicity[0];}

    void vSetMosaicity(double fMosaicity);

    int nRoundIndices(double fMatrix[3][3],const double a3fR0[3],const double a3fR1[3],double& fLowest,double a3fLowestResid[3],int*pnHKL);

    int nNotifyDisplay(const Cstring sReflnlist = "", const Cstring sHeader = "");
    void vGetRMS(float *pfRMSDeg, float *pfRMSMM, float *pfRMSInvAngstrom=NULL); 
 
private: // Some of the above function will be made private
     
    int nDeleteHeaderEntries();

    bool m_bRanking;

public:
     void vSetRanking(const bool bRanking){m_bRanking=bRanking;}
     bool bIsRanking()const{return m_bRanking;}

     void vSetNumberOfResoBins(int nBins){m_nNumberOfResoBins=nBins;}

     bool bIsDisplay();
};  // end of class Crefine


struct CrefineTwinData {
    int    m_a7nReflnNumbers[7];      // Above (for each twin)
    float  m_a7fReflnNumbers[7];
    float  m_fRMSDEG;
    float  m_fRMSMM;
    float  m_fRMSANG;

    CrefineTwinData() {  m_fRMSDEG = 0.0; m_fRMSMM = 0.0; m_fRMSANG = 0.0; memset(m_a7nReflnNumbers,0,sizeof(int)*7);memset(m_a7fReflnNumbers,0,sizeof(float)*7); };
};


#endif   // DT_CREFINE_H

