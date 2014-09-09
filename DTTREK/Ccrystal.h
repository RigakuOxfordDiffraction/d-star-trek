#ifndef DT_CCRYSTAL_H
#define DT_CCRYSTAL_H
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
// Ccrystal.h        Initial author: J.W. Pflugrath           01-May-1995
//    This file is the header file for class Ccrystal
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
#include "Cstring.h"
#include "Cspacegroup.h"
#include "Cimage_header.h"
#include "dtrekvec.h"

//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eCrystal_states {
  eCrystal_unknown_state,
  eCrystal_notrefined_state,
  eCrystal_refined_state
};
const int g_nMaxTwinLaws = 3;
const int g_nResidInfoRMSMm = 0;
const int g_nResidInfoRMSHKL = 1;
const int g_nResidInfoRMSRot = 2;
const int g_nResidInfoNumIndexed = 3;
const int g_nResidInfoFraction = 4;
const int g_nResidInfoIntensity = 5;
const int g_nResidInfoIoverSig = 6;
const int g_nResidInfoRMSEntries = 7;

//+Code begin

class DTREK_EXPORT Ccrystal {

public:

eCrystal_states m_eThe_State;       // The state of the
double           m_fCell[6];
double           m_fOrientAngles[3];

double           m_fOrientMatrix[3][3];      // Rotation matrix times B matrix (UB matrix)
double           m_fOrientVectors[3][3];
double           m_fRotMatrix[3][3];         // The matrix defined by fOrientAngles and

int				 m_nTwinLaws;								// The # of twin laws available.                     
int				 m_nMaxTwinLaws;							// Number of twin laws available.  This is only meaningful if the crystal was obtained from a header file.
double			 m_afTwinFraction[g_nMaxTwinLaws];			// Twin fraction for twin law (if specified).                       
double			 m_aa3fTwinLaw[g_nMaxTwinLaws][3];	        // Twin law (if specified)
double			 m_aa3fRecipShift[g_nMaxTwinLaws + 1][3];	// Reciprocal shift vectors (refined when twins are active).
													        // The first 3vec is relative to the main component.  
													        // The second is relative to the component induced by first twin law  etc ...
// Note on arrays twin components:  The "twin-law index" should be 0 for the (default) identity twin, and 1 for the first twin law component.
// However, in m_afTwinFraction and m_aa3fTwinLaw we store component 1's stuff in [0] since component 0 does not have data.
// In m_aa3fRecipShift we store 0 in [0], 1 in [1] etc...

double           m_fVolume;
double           m_fExposureTime;
Cstring         m_sDescription;
Cstring         m_sCrystal_key;

double           m_fRecipCell[6];
double           m_fRecipVol;
double           m_fBMatrix[3][3];   // Orthogonalized reciprocal cell matrix

double           m_fMosaicity;
double           m_fEffictiveSpectralDispersion;
double           m_fMosaicityOffset;

int              m_nTwins;           // Number of twin componants.  This is only meaningful if the crystal was obtained from a header file.
int              m_nTwin;            // Which twin is this?  (1 .. m_nTwins) Different twins can be loaded from a cyrstal.


Cspacegroup     *m_poSpacegroup;

DTREK_WIN_DLL_DATA_EXPORT static Cstring   ms_sCrystalPrefix;
static Cstring   ms_sCrystalKey;
static Cstring   ms_sCrystalDescription;
static Cstring   ms_sCrystalXUnitCell;
static Cstring   ms_sCrystalXMosaicity;
static Cstring   ms_sCrystalXOrientAngles;
static Cstring   ms_sCrystalXOrientVectors;
static Cstring   ms_sCrystalXSpacegroup;
static Cstring   ms_sCrystalXNumTwinLaws;
static Cstring   ms_sCrystalXTwinLaws;
static Cstring   ms_sCrystalXTwinFraction;
static Cstring   ms_sCrystalXRecipShift;


////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Ccrystal ();                       // Construct an empty crystal object

Ccrystal (Cimage_header& oHeader,int nTwin = 1); // Construct crystal object from image header

Ccrystal(Ccrystal& oOther);

~Ccrystal ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

int nList(const int nFlag = 0);


double fCalcVolume(void);
int nCalcRecipCell(void);
int nCalcBMatrix(const double fWavelength = 1.0);
int nCalcGetRealMatrix(float *pfMatrix);
int nCalcGetRealMatrix(double* pfMatrix);

int nCalcGetBDeriv(double fMat[7][3][3], int *pnNumberRefinedCellParameters,
		   const double fWavelength = 1.0);
int nCalcGetBDeriv(float fMat[7][3][3], int *pnNumberRefinedCellParameters, 
           const float fWavelength = 1.0);

int nGetDerivVecs(int& nNumberRefinedCellParameters,float fNormDirVecs[6][6]);
int nCalcGetRotDeriv(double fMat[3][3][3]);

int nCalcOrientMatrix(const double fWavelength = 1.0);
int nCalcRotMatrix(void);

double fCalcTransForOrientMatrix(double a3x3fOrientMatrixOther[3][3],double* pfTransMat = NULL);
double fCalcTransForOrientMatrix(Ccrystal& oCrystalOther,double* pfTransMat = NULL);

void vGetRecipCell(float *pfRecip);
void vGetRecipCell(double *pfRecip);
void vGetBMatrix(float *pfMat);
void vGetBMatrix(double *pfMat);
void vGetOrientMatrix(float *pfMat);
void vGetOrientMatrix(double *pfMat);

int  nSetOrientMatrix(double *pfOrientMat, const double fWavelength = 1.0f);
int  nSetOrientMatrix(float *pfOrientMat, const double fWavelength = 1.0f);
void vGetRotMatrix(double *pfMat);
void vGetRotMatrix(float *pfMat);

void vGetOrientAngles(float *pfOrientAngles);
void vGetOrientAngles(double *pfOrientAngles);
void vGetOrientAngles(float *pfOr1, float *pfOr2, float *pfOr3);
void vGetOrientAngles(double *pfOr1, double *pfOr2, double *pfOr3);
inline double fGetOrientAngle(const int nDim) { return(m_fOrientAngles[nDim]); }
void vSetOrientAngles(const float *pfOrientAngles);
void vSetOrientAngles(const double *pfOrientAngles);
void vSetOrientAngles(const double fOr1, const double fOr2, const double fOr3);

void vGetCell(float *pfCell);
void vGetCell(double *pfCell);
void vGetCell(float *pfA, float *pfB, float *pfC,
	      float *pfAlpha, float *pfBeta, float *pfGamma);
void vSetCell(const float *pfCell);
void vSetRecipCell(const float* pfRecipCell);
void vSetRecipCell(const double* pfRecipCell);

void vSetCell(const double *pfCell);
void vSetCell(const double fA, const double fB, const double fC,
	      const double fAlpha, const double fBeta, const double fGamma);

double fGetMosaicity();
void  vSetMosaicity(const double fMosaicity);
void  vSetMosaicity(const double fMosaicity,double fEffSpectralDispersion,double fMosaicityOffset);
void  vGetMosaicity(double* pfMosConst);
void  vGetMosaicity(float* pfMosConst);

int nCalcGetRealCell(const float *pfRecipCell, 
		     float *pfRealCell,
		     const double fWavelength = 1.0,
		     const float *pfRecipCellSig = NULL,
		     float *pfRealCellSig = NULL);
int nCountUnique(
        int nNumReso,
        double* pfResoLower,
        double* pfResoUpper,
        int*    pnResoUnique,        
        const int nAnomFlag=0);

double dCalcGetResolution(Crefln *poRefln); 
bool bLinearScaleToVolumeScale(double dLinearScale, double& dVolumeScale);

inline double fGetCell(const int nDim) { return(m_fCell[nDim]); }
inline bool bIsAvailable(void)
{  return (m_eThe_State != eCrystal_unknown_state); }

int nUpdateHeader(Cimage_header* poHeader);
int nUpdateHeaderSigmaValues(Cimage_header* poHeader,
							 float* pfCellSigmas,
							 float* pfRotSigmas,
							 float* pfMosaicitySigmas,
							 float* pfRecipShiftSigmas,
							 float* pfTwinLaw1Sigmas,
							 float* pfTwinFraction1Sigmas);

// These set twin law related stuff.  Please note that 
int nCalcTwinLawDeriv(int nTwinLawIndex,double fMat[3][3][3]);
void vSetTwinLaw(int nTwinLawIndex,double* pfTwinLaw,bool bIsMatrix);      
int nGetTwinLaw(int nTwinLawIndex,double* pfTwinLaw,bool bIsMatrix);
int nClearTwinLaw(int nTwinLawIndex);
int nMaxTwinLaws() { return m_nMaxTwinLaws; };
int nGetRecipShift(int nTwinLawIndex,double* pfRecipShift);
void vSetRecipShift(int nTwinLawIndex,double* pfRecipShift);
int nGetTwinFraction(int nTwinLawIndex,double* pf0) { if (nTwinLawIndex>0) *pf0 = m_afTwinFraction[nTwinLawIndex-1]; else *pf0 = 1.0; return 0; }; 
void vSetTwinFraction(int nTwinLawIndex,double f0) {  if (nTwinLawIndex>0) m_afTwinFraction[nTwinLawIndex-1] = f0;  };
int nGetTwinLawRotAndVector(int nTwinLawIndex,int nMaxHKLToSearch,double* pfRot,double* pfRotSig,double* pfVector,int* pnHKL);


int nFindNumTwins(Cimage_header& oHeader);
int nWhichTwin() { return m_nTwin; };
void vSetTwin(int nTwin) { m_nTwin = nTwin; return;};
int nNumTwins()  { return m_nTwins; };
int nLoadKeywords(int nTwin);       // Loads the ms_sCrystalX???? keywords to reflect the twin componant requested.

int nInitValues(Cimage_header& oHeader,int nTwin);
static int nGetSetResidValues(Cimage_header& oHeader,int nTwinID,int nTwinLaw,double* pfValues,int nAddDelete);
private:

int nInitValues(void);



};  // end of class Ccrystal

#endif   // DT_CCRYSTAL_H

