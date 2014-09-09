#ifndef DT_CDETECTOR_H
#define DT_CDETECTOR_H
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
// Cdetector.h        Initial author: J.W. Pflugrath           01-May-1995
//    This file is the header file for class Cdetector
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
#include "Cgoniometer.h"
#include "Cspatial.h"
#include "Cnonunf.h"
#include "Cimage_header.h"
#include "dtrekvec.h"

//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eDetector_states {
  eDetector_unknown_state,
  eDetector_notrefined_state,
  eDetector_refined_state
};

enum eDetector_data_types {
  eDetector_other_type,
  eDetector_byte,
  eDetector_ubyte,
  eDetector_I2,
  eDetector_I4,
  eDetector_uI2,
  eDetector_uI4,
  eDetector_realIEEE,
  eDetector_compressed

};

//+Forward class declarations

class Cspatial;   // Forward declaration of Cspatial class.
class Cnonunf;    // Forward declaration of Cnonunf  class.
class Csource;	  // Forward declaration of Csource  class.	
class Crotation;  // Forward declaration of Crotation class.

//+Code begin

class DTREK_EXPORT Cdetector {

public:

static Cstring       ms_sDetectorNames;
static Cstring       ms_sDetectorNumber;
static Cstring       ms_sDetectorDescription;
static Cstring       ms_sDetectorKey;
static Cstring       ms_sDetectorVectors;
static Cstring       ms_sDetectorDimensions;
static Cstring       ms_sDetectorSize;
static itr<double>   ms_afChiCoverage;
static double        ms_fResoMin;
static double        ms_fResoStep;

eDetector_states     m_eThe_State;           // The state
eDetector_data_types m_eThe_Type;            // The type of detector
Cstring              m_sType;

int                  m_a2nDim[2];
float                m_a2fSize[2];
float                m_a2x3fVector[2][3];
float                m_a3x3fDD[3][3];        // Detector DD matrix (has m_a3x3fDDRot applied)
float                m_a3x3fDDRot[3][3];     // Detector DD rotation matrix.
float                m_a3fDN[3];             // Detector DN vector
float                m_a3fLocalDN[3];	     // Detector DN vector, local to incident location calculated in nScaleSOnPlate().
float                m_a3fS[3];              // Scattered S vector (computed in some places and placed here whenever).

bool                 m_bNewGonio;
Cgoniometer         *m_poGoniometer;         // Pointer to detector goniometer
                                             // which has detector position

Cspatial            *m_poSpatial;
Cnonunf             *m_poNonunf;

Cstring              m_sDescription;
Cstring              m_sDetector_key;
Cstring              m_sPrefix;

double*			     m_pfPixelShiftArray;	   // For wavelength dependent shifts.
int		             m_a2nPixelShiftArrayDiv[2];
double*              m_pfProfileArray;
int                  m_a2nProfileArrayDiv[2];

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cdetector ();                       // Construct an empty detector object

// Construct a detector object from header

Cdetector (Cimage_header& oHeader, const Cstring& sPre = "",
	   const bool bSpatial = TRUE, const bool bNonunf = TRUE);

~Cdetector ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

// These routines should be used to get detector parameters instead of a call go m_poGoniometer->nGetDatum().
int nSetStandardDetectorParams(double a6fParams[6],int nStart = 0,int nCount = 6);
int nSetModifiedDetectorParams(double a6fParams[6],int nStart = 0,int nCount = 6);
int nGetStandardDetectorParams(double a6fParams[6],int nStart = 0,int nCount = 6);
int nGetModifiedDetectorParams(double a6fParams[6],int nStart = 0,int nCount = 6);


// Modified and standard derivative matrices.
int nCalcGetDDDNDerivStandard(double a6x3x3fDeriv[6][3][3]);
int nCalcGetDDDNDerivModified(double a6x3x3fDeriv[6][3][3]);


int nList(const int nFlag=0);
int nCalcDDDNRotMat(float* pfRotMat,float* pfRotMatDeriv0=NULL,float* pfRotMatDeriv1=NULL,float* pfRotMatDeriv2=NULL);
int nCalcGetDDDN(float *pfDD, float *pfDN);
int nScaleSOnPlate(const float fS[3],float fTS[3]);
int nScaleSOnPlateItr(float fS[3],float fTS[3],float fShift[3]);

int nCalcDDDNRotMat(double* pfRotMat,double* pfRotMatDeriv0=NULL,double* pfRotMatDeriv1=NULL,double* pfRotMatDeriv2=NULL);
int nCalcGetDDDN(double *pfDD, double *pfDN);
int nScaleSOnPlate(const double fS[3],double fTS[3]);
int nScaleSOnPlateItr(double fS[3],double fTS[3],double fShift[3]);


void vUpdateDDDN(void);
int nGetResolution(const float fS0[3], float *pfResoMin, float *pfResoMax,
                   float *pfResoMaxEdge = NULL);
int nGetResolution(const float fS0[3], int nBorder,int nStep,float fMinFraction,
                   float *pfResoMin,float *pfResoMax,
                   int nChiTestArraySize = 0,double* pfResolutions = NULL,double* pfChiCoverageArray = NULL
                   );

double fCalcGetResolution(const double fPx0, const double fPx1, 
			 const float fS0[3], double *pfX = NULL,double *pfChi = NULL);
int nCalcGetPixelFromThetaChi(double f2Theta,double fChi,double& fPix0,double& fPix1,
                                       const float fS0[3],bool bTruncateOutOfBounds);

int nCalcDetCoords(const float fXR[3],
		   const float fS0[3],
		   float *pfXmm, float *pfYmm,
		   float *pfPx0, float *pfPx1);

Cstring sGetParamName(int nParam) {  return m_poGoniometer->sGetName(nParam); };

inline float fGetDistance(void) { return (m_poGoniometer->fGetDistance()); }
inline float fGetSwing(void) { return (m_poGoniometer->fGetSwing()); }
inline bool bIsAvailable(void)
{  return (eDetector_unknown_state != m_eThe_State); }

inline Cstring sGetPrefix(void) { return (m_sPrefix); };

int nDiff(Cdetector& oDetectorAdd,Cdetector& oDetectorSubtract);
int nUpdateHeader(Cimage_header* poHeader, const Cstring& sPre);
int nUpdateFromHeader(Cimage_header& roHeader, const Cstring& rsPre,
		      const bool bSpatial=TRUE, const bool bNonunf=TRUE);

float fGetPixelDistance(float fPxA0=0.0, float fPxA1=0.0, 
			float fPxB0=0.0, float fPxB1=0.0);
int nCalcPixelShift(const float _fS0[3],double a3fRotVec[3],double fLambda1,double fLambda2,int nPix0,int nPix1,float* pfPixWave);
int nCalcPixelShiftBin(double fKa1[2],double fKa2[2],double a3x2x2fShiftBin[3][2][2],int a3x2x2nShiftBin[3][2][2],int a3x2nShiftBin[3][2]);
int nCalcPixelShiftArray(Csource& oSource,Crotation& oRotation,Cgoniometer& oCrysGonio);
int nCalcPixelShift(double fPix0,double fPix1,double* pfWaveShift0,double* pfWaveShift1);
int nCheckReflnCharacteristics(double fCentPix0,double fCentPix1,double* a2x2fA,double* a2fb,double fc);


int nTest(Csource& oSource,Crotation& oRotation);

void vGetLocalNormal(double* pfLocalNormal){vCopyVec3D(m_a3fLocalDN, pfLocalNormal);}
//

private:

int nInitValues(void);
int nInitValues(Cimage_header& oHeader, const Cstring& sPre,
		const bool bSpatial=TRUE, const bool bNonunf=TRUE);


};  // end of class Cdetector

#endif   // DT_CDETECTOR_H

