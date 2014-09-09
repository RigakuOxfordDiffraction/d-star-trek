#ifndef DT_CINTEGRATE_H
#define DT_CINTEGRATE_H
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
// Cintegrate.h        Initial author: J.W. Pflugrath           18-Sep-1995
//    This file is the header file for class Cintegrate
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

#include <time.h>

#include "Dtrek.h"
#include "Cimage.h"
#include "Cscan.h"
#include "Creflnlist.h"
#include "Cdetector.h"
#include "Ccrystal.h"
#include "Cgoniometer.h"
#include "Csource.h"
#include "Cpredict.h"
#include "Crefine.h"
#include "C3Ddata.h"
#include "Cspacegroup.h"
#include "Cstring.h"
#include "Cprofit2D.h"
#include "DTCoreSocket.h"

#include "CImageWait.h"

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
#include "CXprop.h"
#endif

//+Definitions and constants

#define FLOAT_HKL_0 10      // Float HKL Regions.  Used for determining the HKL ellipsoids on different regions of the detector.
#define FLOAT_HKL_1 10      // Float HKL Regions.  Used for determining the HKL ellipsoids on different regions of the detector.
#define PRINT_REGION_0 5    // Print Regions:  Used when printing out spot sizes.
#define PRINT_REGION_1 5
#define MAX_SHOEBOX_INTENSITY_PROFILE 50
#define MAX_MOSAICITY_PROFILE_SIZE 100
#define MAX_MOSAICITY_PROFILES 4
#define MOSAICITY_PROFILE_STEP 0.05



typedef struct _tagRefineResults
{
  float fA, fB, fC, fAlp, fBet, fGam, fRot1, fRot2, fRot3, fMos, fMosMod;
  float fT1, fT2, fT3, fDetRot1, fDetRot2, fDetRot3, fSrcRot1, fSrcRot2;
  float fResMM, fResDeg;
  int   nSeq;
} tagRefineResults;



// For safety on the enum have the unknown method first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum {DTI_DO_2D_PROFIT = 1,
      DTI_DO_NOT_PROFIT = 0 };


//+Code begin

class DTREK_EXPORT Cintegrate {

public:

  bool        m_bScaleShoeboxWithShifts;// Should the size of the showbox be adjusted taking centroid shifts into account?
  float       m_fFattening;             // Peak fattening. (passed down to C3Ddata::nCalcMask())
  int         m_a2nWinMax[2];           // The maximum spot size to be used.  (Restricts large sizes that use HKL ellipsoid limits as a bootstrap.
  int         m_a2nWinMin[2];           // The minimum spot size to be used.  (Forces a shoebox size of a certain dimension).
  int         m_nShoeboxLimit;          // Maximum size of a shoebox in the [2] direction.
  e3Ddata_types m_eSBdatatype;          // Default shoebox datatype, normally ushort, but could be int or float.
  float       m_fSpotSizeMultiplier;    // dtfind/dtintegrate spot size multiplier.
  bool        m_bWriteScanBitmap;       // Do we write out a scan bitmap when we are done with profile fitting?
  bool        m_bWriteEllipsoids;       // Do we write out ellipsoid information?
  bool        m_bSpotChase;             // Do we chase the starting spot centroids?
  bool        m_bKa12;                  // Do we look for Ka12 split?
  float       m_fSigmaAboveBackground;  // Sigma above background for a pixel to be considered 'non-weak'
  int         m_nMinPeakRad;            // Minimum peak radius that will be used.  This is also 'smeared' by the ka1/ka2 characteristics.  
  int         m_nIntegrateWeakFlag;     // See C3Ddata for descriptions.  
  bool        m_bDoFastPrediction;      // 'Fast' prediction algorithm enabled.

  int         m_nProfitNumReflections;  // Number of reflections for profile fitting. (passed to Cprofit2D via nProfileFit()).
  int         m_nProfitMaxImageRange;   // Maximum image range for profile fitting (passed to Cprofit2D via nProfileFit()).

  int         m_nVerbose;               // Verbosity level for debugging and such
  int         m_nProfitFlag;            // Flag to write profiles or not
  int         m_nPad3D;                 // Image padding needed for 3D method 
  double      m_fPad3D;                 // Image padding needed for 3D method.
  int         m_nMinStrongSpotsForIntegrate; // Minimum number of strong reflections we must have before we call nDoIntegrate().
  int         m_n3rdPassPad;            // 3rd pass padding.  (set to 0 to disable third pass integration).

  int         m_nTwinID;                // Twin ID for this integration.
  int         m_nImagesPerScaleBatch;   // Images per scale batch
  int         m_nImagesPerRefineBatch;  // Images per refinement.
  int         m_anSeqNum[2];            // Start and ending SeqNum in poScan to integrate  
  time_t      m_nStartTime;             // Time at start of Cintegrate::nIntegrate().
  int         m_nDisplay;               // Flag to update colling processes display (set with -display flag)
#ifdef SSI_PC
  int		  m_nDisplayFreq;			// Update frequency for calling processes display (set with -display flag argument)
  int		  m_nDisplayCtr;
#endif
  Cstring     m_sOutputHeader;          // Name of output header resulting from integration.
  float       m_a2fResolution[2];       // Resolution bounds.
  Cstring     m_sBatchPrefix;           // Batch prefix
  int         m_nPrerefineFlag;         // Flag for pre-refinement type
  
  int         m_nPrefindFlag;            // Flag for prefind option, number of batches to use prefind on.
  double      m_fPrefindSigma;           // A value to be passed to the Cfind object when doing prefind
  
  bool        m_bDetermineFixedMosaicity;// Do we fix mosaicity after pre-refinement?
  float       m_fFixedMosaicity;        // Allows user to fix the mosaicity in integration.
  float       m_a2fMosaicityModel[2];   // Allows user to modify refined mosaicity
  float       m_fMaxIntersectionFraction;       // Passed to profile fitting (see description in Cprofit2D)
  double      m_a3fPCellVector[3];
  double      m_fPCellAngle;

  Cstring     m_sReflnFilename;
  Cstring     m_sReflnPartialsFilename;
  Cstring     m_sReflnRejectsFilename;
  Cstring     m_sReflnProfitFilename;

  int         m_nOverlapFlag;           // Flag to check overlaps or not
  float       m_fDelHKLThresh;          // Threshold for delta HKL cutoff

inline void vSetPrerefineFlag(const int nFlag) { m_nPrerefineFlag = nFlag; }
inline void vSetPrefindFlag(const int nFlag) { m_nPrefindFlag = nFlag; }
void vSetPrefindSigma(const double fSigma){m_fPrefindSigma=fSigma;}
inline void vSetReflnFilename(const Cstring& sNameIn) { m_sReflnFilename = sNameIn; }

inline void vSetReflnProfitFilename(const Cstring& sNameIn){m_sReflnProfitFilename = sNameIn;}
inline Cstring sGetReflnProfitFilename(void){return (m_sReflnProfitFilename);}

inline void vSetVerboseLevel(const int nLev) { m_nVerbose = nLev; }
inline void vSetBatch(const int nImagesPerScaleBatch,const int nImagesPerRefineBatch) { m_nImagesPerScaleBatch = nImagesPerScaleBatch; m_nImagesPerRefineBatch = nImagesPerRefineBatch; }
inline void vSetBatchPrefix(const Cstring& rsTemp) { m_sBatchPrefix = rsTemp; }

   void vSetWaitTime(const double fWaitTime){THE_IMAGE_WAIT_PTR->vSetWaitTime(fWaitTime);}

  // Basic structures.

  Cimage     *m_poImage;      // Pointer to image object		
  bool        m_bNewImage;
  Cscan      *m_poScan;       // Pointer to scan object to integrate
  bool        m_bNewScan;
  Creflnlist *m_poReflnlist;        // Pointer to reflnlist object to put results in.
  Creflnlist *m_poReflnlistRefine;  // Stores reflections used in the refinement.
  Creflnlist *m_poReflnlistPartial; // Reflection list containing partials.
  bool        m_bNewReflnlist;
  Crefine    *m_poRefine;     // A Crefine  object for refining spots

  Cimage_header *m_poHeader;       // Pointer to header object

  CDTCoreSocket* m_poSocket;

  // Error related arrays and structures. 

  float       m_afErrorLimits[nMaxErrorStatusTypes];        // nMaxErrorStatusTypes possible error limits
  int         m_anNumReflns[nMaxErrorStatusTypes];			// Error's during integration. (see ms_asErrorStrings)


  int         m_nNumRefPartialStart;
  int         m_nNumRefPartialEnd;
  int         m_nNumPredictedUnique;
  int         m_nNumRefBadErrors;
  int         m_nNumRefNoErrors;
  int         m_nNumRefSpecial;
  int         m_nNumRefWeak;

  static int ms_nNewFlag;
  static int ms_nActiveFlag;
  static int ms_nDelFlag;

  static Cstring ms_snDataFlag;
  static Cstring ms_snErrorFlag;
  static Cstring ms_sfWidth0;
  static Cstring ms_sfWidth1;
  static Cstring ms_sfDataAddr;
 
  Cstring m_sProfitDiagnostic;

private:

  // Sorting arrays.

  int        *m_pnHKLIndex;

  // 2D profile information.
#ifdef NEWPROFIT
  Cprofit2DSliceModel *m_poProfit2DPerSlice;
#else
  Cprofit2DPerSlice *m_poProfit2DPerSlice;
#endif

  // The following variables are used in centroid determination.
  C3Ddata    *m_po3DCentroid;
  C3Ddata    *m_po3DCentroid2;

  // Misc Constants. (nDoRefine() nDoIntegrate())

  int         m_nRefineFlag;       // Refinement flag 0= do not ref; 1= do refine
  int         m_nRefineMinRefln;   // Minimum reflections needed for refinement
  int         m_nRefineMaxRefln;   // Maximum reflections allowed in refinement
  int         m_nDoRefineCount;    // Indicates which pass of refinement parameters are currently in use.
  int         m_nDoIntegrateCount; // Only used for printing I think.


  // Calculated Constants.  Calculated once, (or once per image), and then used.

  float       m_fRotIncrement;          // Rotation increment for the scan.
  int         m_nDim0, m_nDim1;         // Dimensions of current image
  float       m_a3fRotVec[3];           // Rotation axis vector
  float       m_fDetectorGain;          // Gain of detector.
  int         m_nNumBeforePred;         // Loaded before each new image prediction done
  int         m_nNumAfterPred;          // Loaded after each new image prediction done


  // Counter variables. 
  int         m_nCurrentImgNum;     // Image number in a scan.  Incremented as we traverse through images.

  // The following variables are updated during nUpdateRefine()

  float       m_a3x3fDetDD[3][3];               // Detector projector matrix
  float       m_a3fDetDN[3];                    // Detector projector vector
  float       m_a3fS0[3];                       // Source vector
  float       m_a3fPolNorm[3];                  // Normal to polarization vector
  float       m_fLenS0;                         // Source vector length
  float       m_fWavelength;                    // Source wavelength
  float       m_a3x3fInvCrysOrientGonMat[3][3]; // Inverse matrix used to convert recip lattice point  to hkl.
  float       m_a3x3fInvCrysRotGonMat[3][3];
  float       m_a3x3fGonMat[3][3];

  tagRefineResults *m_ptRefineResults;
  int         m_nNumRefineResults;
  int         m_nSizeRefineResults;

  // Misc. calculated variables.

  float       m_fOverallRotStart;
  float       m_fOverallRotEnd;
  int         m_nNonunfMaskBadPix;  // These are masks for the non-uniformity field of the reflection list. 
  int         m_nNonunfMaskRing;    // These are masks for the non-uniformity field of the reflection list.

  C3DdataDetectorProfile m_a2oDetectorProfiles[2];      // We are reading from one profile, and outputing to the other.
  int                    m_nProfileLoading;             // Which profile is currently getting loaded.
  int                    m_nProfileReading;             // Which profile is being read (not necc. different m_nProfileLoading)
  double                 m_fTargetDetectorProfilesLevel;// CinterpMesh level to which we must fill.
  
  float      *m_pfPerImageFactors; // Factors per image in a shoebox
  float      m_pfShoeboxIntensityProfilesIntensity[MAX_SHOEBOX_INTENSITY_PROFILE-2][MAX_SHOEBOX_INTENSITY_PROFILE];
  int        m_pnShoeboxIntensityProfilesCounts[MAX_SHOEBOX_INTENSITY_PROFILE];
  unsigned char* m_pcBackgroundMask;  
  int        m_aanMosaicityProfilesCounts[MAX_MOSAICITY_PROFILES][MAX_MOSAICITY_PROFILE_SIZE];
  double     m_afMosaicityValues[MAX_MOSAICITY_PROFILES];
  float*     m_pfNearestCentroidBuf0;
  float*     m_pfNearestCentroidBuf1;
  int		 m_nNearestCentroidSize;

  Cstring m_sXpropString;

  // Printing of shoebox data. (Rogue's gallery)
  bool       m_bRogueActive;
  int        m_nRogueMask;
  int        m_nRogueNegMask;
  int        m_nRoguePrintNumber;
  int        m_nRogueImageCount;
  int        m_nRogueShoeboxImageCount;
  int        m_nRogueDim0;
  int        m_nRogueDim1;
  int        m_nRoguePx0;
  int        m_nRoguePx1;
  int        m_nMaxRoguePx1;
  int        m_nRogueBorder;
  Cstring    m_sRogueTemplate;
  Cstring    m_sRogueNextImage;
  Creflnlist*m_poRogueReflnlist;
  Cimage*    m_poRogueImage;
  int*       m_pnRogueReflnlistHKLSort;
  int        m_nRoguePlaceInImages;
  int        m_nRogueImageNumber;
  int        m_nMaxRogueShoeboxes;
  int        m_nRogueShoeboxes;
  int        m_nRoguePass;
  int        m_nRogueMinBoxDim;

  // Variables related to output of profiles and reflections.

  char       *m_pcSelectFields;             // Selected fields in output reflection list.
  char       *m_pcSelectFieldsPartial;      // Selected fields for the partial list.
  int         m_nReflnFileOutCount;         // Used to update important information for size-constrained profile fitting.
  int         m_nReflnPartialOutCount;      // Used to update important information for size-constrained profile fitting.
  FILE*       m_pFReflnFileOut;             // Output file stream for reflection file
  FILE*       m_pFReflnPartialOut;          // Output file stream for partials.
  FILE*       m_pFReflnBadOut;              // Output file stream for rejected reflns file

  // Shoebox related variables.


  C3Ddata   **m_ppoShoeBFree;               // Array of freed shoeboxes
  int        *m_pnShoeBFreeSize;            // Allocated size of each freed shoebox
                                            //  use this to check shoebx size so paging minimized
  int         m_nNumShoeBAlloc;             // Number of slots allocated in free list
  int         m_nNumShoeBFree;              // Number of shoeboxes in free list;

  int         m_nErrorDoNotProcessMask;     // Mask for what not to process
  int         m_nErrorDoNotPrintMask;       // Mask for what should not be printed in summary table.


  // Some things for oblique incidence correction

  int    m_nNumObliqueCorr;
  float *m_pfObliqueCorr;
  float *m_pfObliqueCorrNorm;
  float *m_pfObliqueCorrType;
  float  m_fRefOblique;

  // Rings excluded with -ring command.

  int    m_nNumExcludeResoRings;
  float *m_pfExcludeReso[2];


  // Reflection indices.

  int    m_nFI_nDataFlag;
  int    m_nFI_nErrorFlag;
  int    m_nFI_fWidth0;
  int    m_nFI_fWidth1;

public:

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cintegrate ();                       // Construct an empty integrate object

Cintegrate (Cscan *poScan, Creflnlist *poReflnlist);  // Integrate a scan
Cintegrate (Cimage_header *poHeaderIn,
	    Creflnlist *poReflnlist); // Integrate an image

~Cintegrate ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

#if ( (!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)) )
    CXprop* m_poXprop;
#endif

// 2D profile fitting.

int nInitProfit2D(Cdetector& oDetector);

int nInitValues(void);
int nList(void);
int nIntegrate(const int nNumImages = 99999);

int nExpandReflnlist(Creflnlist* poReflnlistIn = NULL);
int nBuildShoeboxes(Creflnlist *poReflnlistActive,
		    Cimage *poImageIn, Crotation *poRotation,
		    Cdetector *poDetector);

int nNewShoebox(Crefln *poRefln, Crotation *poRotation);


// nDoIntegrate() related functions.

int  nCalcFloatHKL(Crefln *poRefln,
		   Cdetector *poDetector, const int nModeIn,
		   float *pfMax);
int  nCalcSmallestShoeboxRange(int nExt2,int nCenterSlice,
                               int nCent0,int nCent1,
                               int& nStart0,int& nExt0,int& nStart1,int& nExt1);

int nPass0Shoebox(Crefln *poRefln, Cdetector *poDetector);
int nDeleteMarkedReflns(Creflnlist *poReflnlist, const int nFlag);
int nDoRefine(Creflnlist *poReflnlist, Crefine *poRefine);

void vSetFixedMosaicity(double fFixMosaicity,Cpredict* poPredict = NULL);
void vSetGain(double fGain) { m_fDetectorGain = fGain; };
void vAddExcludeResoRing(const float fResoMin, const float fResoMax);
void vSetResolution(const float fResoMin, const float fResoMax);
inline void vSetMosaicityModel(const float fMosMul, const float fMosAdd)
  { m_a2fMosaicityModel[0] = fMosMul; m_a2fMosaicityModel[1] = fMosAdd; }
void vListResults(void);
void vListRefineResults(void);
void vListIntensityProfiles(void);
void vSetRogueImageDim(int nDim0,int nDim1) { m_nRogueDim0 = nDim0; m_nRogueDim1 = nDim1; };
void vSetRogueImageCutoff(int nImageNumber);
int nSetRogueData(int nMask,int nNegMask,Cstring& sTemplate,Cstring& sReflnlist,int nImageNumber=-1,int nMinDimension=-1);

int nAddRogue(int nMask,Crefln* pnHKL,C3Ddata& oShoebox,float* pfExternData = NULL,int* pnMask = NULL);
void vTermRogue();


 private:

int  nSetErrorMask(const int nErrorNumber, const bool bCountIt, int *pnMask);
int  nSetErrorMasks(const int nErrorNumbersBitwise, const bool bCountIt);
bool bPixelInLimits(const int nPx0, const int nPx1);
int  nExpandShoeboxes(const int nAdditional = 0);

int  nCalcMask(Crefln *poRefln, Cdetector *poDetector, C3Ddata *poShoeboxIn);
		  
int  nUpdateRefine(Csource *poSource, Cdetector *poDetector,
		   Ccrystal *poCrystal, Cgoniometer *poCrysGonio,
           Crotation *poRotation);

C3Ddata* poGetShoebox(const int nExt0, const int nExt1, const int nExt2,
		      const int nOrig0, const int nOrig1, const int nOrig2);

void vDeleteFreeShoeboxes(void);
void vDeleteShoebox(C3Ddata *poShoebox);

// Functions dealing with background mask.
void vInitBackgroundMask();
void vAddBackgroundMask(C3Ddata* poShoebox);

int  nDoPrefind();

};  // end of class Cintegrate

#endif   // DT_CINTEGRATE_H
