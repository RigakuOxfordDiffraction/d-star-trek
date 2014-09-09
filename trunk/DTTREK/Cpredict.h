#ifndef DT_CPREDICT_H
#define DT_CPREDICT_H
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
// Cpredict.h        Initial author: J.W. Pflugrath           01-May-1995
//    This file is the header file for class Cpredict
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
#include "Cnonunf.h"
#include "dtrekvec.h"

//+Definitions and constants

enum ePredict_methods {
  ePredict_unknown_method,
  ePredict_2D_method,
  ePredict_3D_method
};

//+Code begin

class DTREK_EXPORT Cpredict {

public:

  int          m_nDisplay;     // Flag to update display of calling process
  int          m_nNumDetectors;
  int          m_nNumCrystals;
  int		   m_nMaxTwinLaws;
  int          m_nWhichDetector;
  int          m_nWhichCrystal;
  int		   m_nWhichCrystalTwinLaw;
  Cstring     *m_psDetectorNames;
  bool        *m_pbActive;     // Which crystal componants are active?
  bool         m_bUseEdgeReso; // Use the resolution at edge of detector, (otherwise at corners (i.e. max bounds))

  Cdetector  **m_ppoDetector;  // Pointer to pointer of detector object(s)
  Csource     *m_poSource;     // Pointer to source object
  Ccrystal   **m_ppoCrystal;   // Pointer to crystal object(s).  
  Cgoniometer *m_poCrysGonio;  // Pointer to crystal goniometer;

  Crotation   *m_poRotation;        // Pointer to rotation object to predict reflns for
  Creflnlist  *m_poReflnlist;       // Pointer to reflnlist object


  double m_fReflnWidthMax;
  double m_fReflnWidthMaxRad;
  double m_fResolutionMin;      // Min resolution limit in A (99999 is min)
  double m_fResolutionMax;      // Max resolution limit in A (0.5 A is max)
  int    m_nStat;               // Initialization status.
  int    m_nScanNumber;         // Which scan number


ePredict_methods m_eMethod;      //  2D or 3D method

private:
  Cstring       m_sStrategyPrefix;
  bool m_bNewDetNames;
  bool m_bNewDetector;
  bool m_bNewSource;
  bool m_bNewCrystal;
  bool m_bNewCrysGonio;
  bool m_bNewRotation;
  bool m_bNewReflnlist;

  double m_afDD[3][3];             // Detector DD matrix
  double m_afDN[3];                // Detector DN vector
  double m_afInvDD[3][3];          // Inverse of detector DD matrix
  double m_fDet;                  // General matrix determinant
  double m_fDetDistance;          // Distance from crystal to detector
  double m_fMax2Theta;
  float  m_fDetResolMin, m_fDetResolMax;  // Min and Max resolution on detector

  double m_afGonMatrix[3][3];      // Goniometer rotation matrix at datum

  double m_afCrysMatrix[3][3];     // Crystal orientation matrix when goniometer
                               //   is at 0,0,0
  double m_afGonCrysMatrix[3][3];  // Previous 2 matrices multiplied together
  double m_afRecipShift[3];			// Reciprocal shift.
  double m_afRecipShiftRotated[3];  // Reciprocal shift vector rotated by goniometer matrix and phi rotation (calculated per-reflection)
  int   m_anSortedAxes[3];        // Indices for which cell direction to loop over

  double m_afS0[3];                // Source direction vector (normalized?)
  double m_fWavelength;           // Source wavelength in Angstroms
  double m_fDstarSqMin;           // d* squared minimum
  double m_fDstarSqMax;           // d* squared maximum
  double m_afPol[4];               // Polarization from poSource
  double m_afSN[3];                // S0 X Polarization
  double m_fSourceDispersion;     // Source spectral dispersion
  double m_fS0dotDN;              // Result of fDot3D(fS0, fDN);
  double m_fDNdotDN;              // Result of fDot3D(fDN, fDN);
  double m_fLenDN;                // |DN|;
  double m_afXRminusS0[3];         // XR - S0,
  double m_afE3crossS0[3];         // Result of vCross3D(fRotVec, fS0,...)
  double m_afRotVec[3];            // Rotation vector extracted from poRotation
  double m_afRotStartMatrix[3][3]; // Composite rotation matrix at start of Rotation
  double m_afRotEndMatrix[3][3];   // Composite rotation matrix at end of Rotation
  double m_afRotPastEndMatrix[3][3]; // Composite rotation matrix past end of Rotn.
  double m_fRotOverallStart, m_fRotOverallEnd;  // Over-all start and end of rotation
  double m_fMaxLorentz;   // Max allow Lorentz factor before flagging as bad

  double m_a3fLocalDN[3];


// Some things calculated for a reflection:

  double m_fMosCoeffA;    // Mosaicity can be refined using these two per-reflection parameters.
  double m_fMosCoeffB;
  double m_fReflnDstarSq; // d* * d* for the reflection
  double m_fReflnStart;   // Starting rotation value of refln in degrees
  double m_fReflnEnd;     // Ending   rotation value of refln in degrees
  double m_fReflnWidth;   // Rotation width of refln in degrees
  double m_fReflnRot0;    // Predicted rot value when ref on Ewald sphere
  double m_fReflnResolution;  // Resolution in Angstrom
  double m_fLorentz;      // Lorentz factor: Itrue = Iobs/Lorentz)
  double m_fPx0;          // 0 Pixel coordinate.
  double m_fPx1;          // 1 Pixel coordinate.
  double m_fPxNorm0;      // 0 Pixel coordinate (normalized).
  double m_fPxNorm1;      // 1 Pixel coordinate (normalized).  
  double m_fXmm;
  double m_fYmm;
  double m_fZRelativeCoordinate;
  double m_fPx0Drift;      // Drift parameter (dpixel/dRot) (w.r.t. a positive change in rotation angle)
  double m_fPx1Drift;      // Drift parameter (dpixel/dRot) (w.r.t. a positive change in rotation angle)
  double m_afXR[3];        // Reciprocal lattice coordinate at diffracting condition
  int   m_nNonunfFlag;     // A flag to tell what to do with reflns with
                           // predicted pixel in bad nonunf region

  double m_fScaleLength;                // Scale factor for crystal cell lengths
  bool m_bCalcImageDrift;               // (loads CalcPx0Drift CAlcPx1Drift) Pixel drift per degree.
  bool m_bCalcMosaicityLinearCoeffs;    // (loads MosCoeffA,MosCoeffB) Loads two mosaicity coeffs for each refln.
  bool m_bCalcMinMosaicityForRefln;     // (loads MosCoeffA) Loads minimum mosaicity required to predict this refln.

  bool m_bDetectorIsCylindrical;        // Bool to flag whether detector is cylindrical or flat

  itr<float>  m_afPrePredictRotMid;     // Pre-predict Rot-mid sorted position.
  itr<int>    m_anPrePredictPackedHKL;  // Pre-predict Packed HKL array.
  double      m_fPrePredictRotWidth;    // Pre-predict 95% rotation width cutoff.
  bool        m_bDoFastPrediction;      // Are we using the 'fast' predition method.  Default should be 'no'.
  double      m_fWholeRotStart;         // What is the start of the entire rotation wedge?
  double      m_fWholeRotEnd;           // What is the end of the entire rotation wedge?

#define     DTREK_PREDICT_VERBOSE_MOSAICITY     0x0001
  DTREK_WORD m_wVerbose;

public:
  void vSetVerbose(DTREK_WORD wVerbose){m_wVerbose=wVerbose;}


////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cpredict ();       // Construct an empty predict object
                   // In the future this should not be allowed to happen.

Cpredict (Csource *poSourceIn, Cdetector *poDetectorIn,
	  Ccrystal *poCrystalIn, Cgoniometer *poCrysGonioIn,
	  Crotation *poRotationIn=NULL, Creflnlist *poReflnlistIn=NULL);

Cpredict (Cimage_header& oHeader, Creflnlist *poReflnlistIn,int nScanNumber = 0, const Cstring& sStrategyPrefix="");

~Cpredict ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

int nInitValues(void);
int nInitValues(Cimage_header& oHeader);
int nList(void);
int nPredict(Crotation *poRotationIn=NULL, Creflnlist *poReflnlistIn=NULL);
int nPredict(Creflnlist& oListIn,double fPad);
int nPredictWholeDataSet(double fRotStart,double fRotEnd);
int nPredictInsertRefln(int nHKL[3],double fStartRad,double fEndRad,double fX12[3],Crefln& oRefln);
int nPredictReflnsFast(const double fRotStart, const double fRotEnd);
int nPredictReflns(const double fRotStart, const double fRotEnd);

int nSetupDetector(const int nWhich = 0);
int nSetupSourceCrystal(const int nWhichTwin = 0,const int nWhichTwinLaw = 0);
int nSetupRotation(const double fRotStart, const double fRotEnd);
int nGetIndexLargest3D(const double fVec[3]);
int nSortAxes(void);

double fDetermineRockingWidth(double fFractionUsed);
void vSetResolution(const double fResoMin, const double fResoMax);
 double dGetResolution(const int nMinMax);
int  nExpandReflnlist(Creflnlist *poReflnlistIn = NULL);
Cnonunf* poGetNonunf(void);
inline void vSetNonunfFlag(const int nNonunfFlag = 0)
{ m_nNonunfFlag = nNonunfFlag; }
inline void vSetScaleLength(const double fScaleLength)
{ m_fScaleLength = fScaleLength; }
inline void vSetMaxLorentz(const double fMaxLorentz) { m_fMaxLorentz = fMaxLorentz; }
inline void vSetImageDriftOn() { m_bCalcImageDrift = TRUE; };
inline void vSetMosaicityLinearCoeffs() { m_bCalcMosaicityLinearCoeffs = TRUE; };
inline void vSetMosaicityMinMosaicity() { m_bCalcMinMosaicityForRefln  = TRUE; };
inline void vUseEdgeReso() {m_bUseEdgeReso = TRUE; };
inline void vSetFastPredict(double fRotStart,double fRotEnd) { m_bDoFastPrediction = TRUE; m_fWholeRotStart = fRotStart; m_fWholeRotEnd = fRotEnd;};
  
inline void vSetCrysMosaicity(const float fMosaicity) 
{ int nx;
    for (nx=0;nx<m_nNumCrystals;nx++) 
        m_ppoCrystal[nx]->vSetMosaicity(fMosaicity);
};
inline void vSetCrysMosaicity(const double* pfMosaicity) 
{ int nx;
    for (nx=0;nx<m_nNumCrystals;nx++) 
        m_ppoCrystal[nx]->vSetMosaicity(pfMosaicity[0],pfMosaicity[1],pfMosaicity[2]);
};

inline bool bIsAvailable() { return (m_nStat == 0); };
int	  nCalcNearestCentroids(Creflnlist& oList,int nReflnStart,int nReflnEnd,double fObsRotStart,double fObsRotEnd,int nMinSearchRadius,float* pfCentroid0,float* pfCentroid1);

 private: // Some of the above function will be made private

void vSwap2Values(int* pnInt1, int* pnInt2);
void vSwap2Values(float* pfFloat11, float* pfFloat2);
void vSwap2Values(double* pfFloat11, double* pfFloat2);

int  nSolveQuadratic(const double fA, const double fB, const double fC,
		     double *pfV);

int nOutsideLimits(const int nMin, const int nMax,
		  int* pnBegin,  int* pnEnd);

void vAdjustLoopLimits(const int nLoop,
		       const double fQ1, const double fQ2,
		       const double fR1, const double fR2,
		       int* pnBegin, int* pnEnd);

double fCalcGetPolarz(void);

int   nCalcDetCoords(void);
int   nCalcDriftDetCoords(void);

int   nCalcRecipCoords(const double fX[3], const double fRotStartRad,
		       const double fRotEndRad);


bool bPixelInLimits(const int nPx0, const int nPx1);
inline int nIntMinusOne(const double fFloat) { return ( (int)(fFloat-1.0));}
inline int nIntPlusOne(const double fFloat)  { return ( (int)(fFloat+1.0));}

};  // end of class Cpredict

#endif   // DT_CPREDICT_H

