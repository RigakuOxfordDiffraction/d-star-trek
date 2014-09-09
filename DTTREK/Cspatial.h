#ifndef DT_CSPATIAL_H
#define DT_CSPATIAL_H
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
// Cspatial.h        Initial author: J.W. Pflugrath           03-Mar-1995
//    This file is the header file for class Cspatial
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

#include "Cstring.h"
#include "Cimage.h"
#include "Cdetector.h"

//+Definitions and constants

// For safety on the enum have the unknown state first.  This way if things
//  are initialized to 0, then the program will sort of know.

enum eSpatial_states {
  eSpatial_unknown_state,
  eSpatial_available_state
};

enum eSpatial_types {
  eSpatial_unknown_type,
  eSpatial_simple_type,
  eSpatial_complex_type,
  eSpatial_interp_table
};


enum eSpatial_geo_types {
    eSpatial_flat_geometry,
    eSpatial_cylindrical_geometry
};


//+Forward class declarations

class Cdetector;   // Forward declaration of Cdetector class.
class CImagePool;

//+Code begin

class DTREK_EXPORT Cspatial {

public:

static Cstring ms_sSpatialKey;
static Cstring ms_sSpatialDistortionType;
static Cstring ms_sSpatialDistortionInfo;
static Cstring ms_sSpatialDistortionVectors;
static Cstring ms_sSpatialBeamPosn;
static Cstring ms_sSpatialTypeUnknown;
static Cstring ms_sSpatialTypeSimple;
static Cstring ms_sSpatialTypeComplex;
static Cstring ms_sSpatialTypeNone;
static Cstring ms_sSpatialTypeInterp;

static Cstring ms_sPixelSize;
static Cstring ms_sPscale;
static Cstring ms_sXintStart;
static Cstring ms_sYintStart;
static Cstring ms_sXinvStart;
static Cstring ms_sYinvStart;
static Cstring ms_sXintStep;
static Cstring ms_sYintStep;
static Cstring ms_sXinvStep;
static Cstring ms_sYinvStep;
static Cstring ms_sXSize;
static Cstring ms_sYSize;
static Cstring ms_sXBeam;
static Cstring ms_sYBeam;
static Cstring ms_sBadFlag;


private:
eSpatial_states     m_eThe_State;
eSpatial_types      m_eThe_Type;
eSpatial_geo_types  m_eThe_Geometry;

float        m_a2fBeamPx[2];      // Pixel of primary hitting the detector
                                  //  when detector is normal to beam. The mm
                                  //  coordinates of this position are 0,0.
float        m_a2fNominalPxSize[2];//Nominal pix size in mm in the fast and slow
                                  // directions.  This will be the actual pixel 
                                  // size when eThe_Type is eSpatial_simple_type
float        m_a2x2fDirVec[2][2]; // Direction vector of the 1st (fast) and 
                                  //   2nd (slow) pixel directions.
float        m_a2x2fInvVec[2][2]; // Direction vector of the 1st (fast) and 

float        m_a2fMaxPixel[2];    // Maximum pixel number in the fast and slow
                                  //   directions.  The minimum pixel number is
                                  //   assumed to be 0.

Cstring      m_sDescription;

// Some fileset information when using Marty Stanton's calibration routines:

Cstring         m_sFileset_name;

static Cstring  ms_sCalpar;
static Cstring  ms_sX_int;
static Cstring  ms_sY_int;
static Cstring  ms_sInv_x_int;
static Cstring  ms_sInv_y_int;

Cimage_header *m_poCalpar;
Cimage        *m_poX_int;
Cimage        *m_poY_int;
Cimage        *m_poInv_x_int;
Cimage        *m_poInv_y_int;

float          m_fPxScale, m_fPxSize;
int            m_n0Start,  m_n1Start, m_n0InvStart, m_n1InvStart;
int            m_n0Step,   m_n1Step,  m_n0InvStep,  m_n1InvStep;
float          m_a2fMMOriginPx[2];
int            m_n0ImageSize, m_n1ImageSize;
float          m_fXCenter, m_fYCenter;

float          m_fRadius;  // cylindrical detector radius
float          m_fCylindricalDetectorShift; // this value is equal to (D-Rc), where D is the distance between
                                            // the center of the lab. coord. system and the center of the detector.
                                            // Rc is the radius of the cylindrical detector.
                                            // In the ideal case D=RC, i.e. m_fCylindricalDetectorShift=0.0  

int            m_nBadFlag;

int            m_nNumErrorOut;
int            m_nVerbose;
//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

public:

Cspatial ();                 // Construct an empty spatial distortion object
Cspatial (const Cstring& sFileset);  // Construct a spatial distortion object by 
                                   //    reading from fileset sFileset
Cspatial(const float fBeam0, const float fBeam1,
         const float fPxSize0, const float fPxSize1,
         const int nDim0=511, const int nDim1= 511,
	 const float fDir00=1.0, const float fDir01=0.0,
	 const float fDir10=0.0, const float fDir11=1.0);

 // Construct from information in the header 
 //    of an image
Cspatial(Cimage_header& oHeader, const Cstring& sPrefix = "");

~Cspatial ();

float fGetCylindricalDetectorRadius()const{return m_fRadius;}
float fGetCylindricalDetectorShift()const{return m_fCylindricalDetectorShift;}

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
//

// int nRead(void);
int nRead(const Cstring& sName);
int nRead_BS(const Cstring& sName);

// int nWrite(void);
// int nWrite(const Cstring& sName);

int nList(void);

int nPixeltoMM(const float fPx0, const float fPx1, float *pfXmm, float *pfYmm, float *pfZRelativeCoordinate);
int nMMtoPixel(float fXmm, float fYmm, float fZRelativeCoordinate, float *pfPx0, float *pfPx1,bool bTrucateOutOfBounds = TRUE);

int nSetBeamPosition(const float fBeamPx0, const float fBeamPx1);
int nGetBeamPosition(float *pfBeamPx0, float *pfBeamPx1);
int nSetScaleFactors(const float fScalePx0, const float fScalePx1);
int nSetMaxPixel(const int nMax0, const int nMax1);
int nGetMaxPixel(int nDim) { return (int) m_a2fMaxPixel[nDim]; };
int nDiff(Cspatial& oSpatialAdd,Cspatial& oSpatialSubtract);

int nDistortImage(Cimage *poImage);
int nUndistorImage(Cimage *poImage);
inline bool bIsAvailable(void)
                 { return (m_eThe_State == eSpatial_available_state); }
inline eSpatial_geo_types eSpatialType(void)
{ return (m_eThe_Geometry); };

inline void vSetRadius(const float fRadius) { m_fRadius = fRadius; }; // Not normally used.
inline void vSetVerbose(const int nVerbose) { m_nVerbose = nVerbose; }
inline float fGetNominalPixelSize(const int nWhich) 
  { return (m_a2fNominalPxSize[nWhich]); }
int nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre = "");

bool bChangeBinning(int nFastRatio, int nSlowRatio);

private:

int nInitValues(void);
int nInitValues(const Cstring& sFileset);
int nInitValues(const float fBeam0, const float fBeam1,
                const float fPxsize0, const float fPxsize1,
                const int nDim0=511, const int nDim1= 511,
		const float fDir00=1.0, const float fDir01=0.0,
		const float fDir10=0.0, const float fDir11=1.0);
int nStantonPixeltoMM(const float fPx0, const float fPx1,
                      float *pfXmm, float *pfYmm);
int nStantonMMtoPixel(const float fXmm, const float fYmm,
                      float *pfPx0, float *pfPx1);

void vInvert(void);

CImagePool*       m_pImagePool;
};  // end of class Cspatial

#endif   // DT_CSPATIAL_H

