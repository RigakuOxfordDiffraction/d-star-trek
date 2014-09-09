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
// Cdetector.cc            Initial author: J.W. Pflugrath           01-May-1995
//    This file contains the member functions of class Cdetector.
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
//+Description

//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Cdetector.h"          // Class definition and prototypes
                                // Cdetector.h includes others
#include "Csource.h"
#include "Crotation.h"
#include "dtrekdefs.h"
#include "dtsvd.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

//+Definitions, constants, and initialization of static member variables

Cstring Cdetector::ms_sDetectorNames       = D_K_DetectorNames;
Cstring Cdetector::ms_sDetectorNumber      = D_K_DetectorNumber;
Cstring Cdetector::ms_sDetectorKey         = D_K_DetectorKey;
Cstring Cdetector::ms_sDetectorDescription = D_K_DetectorDescription;
Cstring Cdetector::ms_sDetectorVectors     = D_K_DetectorVectors;
Cstring Cdetector::ms_sDetectorDimensions  = D_K_DetectorDimensions;
Cstring Cdetector::ms_sDetectorSize        = D_K_DetectorSize;


//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cdetector::Cdetector()
{
  (void) nInitValues();
}

Cdetector::Cdetector(Cimage_header& oHeader, const Cstring& sPre,
		     const bool bSpatial, const bool bNonunf)
{
  (void) nInitValues(oHeader, sPre, bSpatial, bNonunf);
}

Cdetector::~Cdetector()
{
    m_eThe_State = eDetector_unknown_state;

    if( m_bNewGonio && (NULL != m_poGoniometer) )
    {
        delete m_poGoniometer;
        m_poGoniometer = NULL;
        m_bNewGonio    = FALSE;
    }
  
    if( m_poNonunf )
    {
	    delete m_poNonunf;
	    m_poNonunf = NULL;
    }
    
    if( m_poSpatial )
    {
	    delete m_poSpatial;
	    m_poSpatial = NULL;
    }

    if (m_pfPixelShiftArray) 
    {
      delete[] m_pfPixelShiftArray;
      m_pfPixelShiftArray = NULL;
    }
    if (m_pfProfileArray)
    {
      delete[] m_pfProfileArray;
      m_pfProfileArray = NULL;
    }

  // Should we delete poGoniometer?, poSpatial?, poNonunf? YES, if constructed
  //  by this class.
}

int
Cdetector::nInitValues()
{
  int i;       // Loop counter

  m_eThe_State      = eDetector_unknown_state;
  m_eThe_Type       = eDetector_other_type;  
  m_sDescription    = "UNKNOWN";
  m_sDetector_key   = "";
  m_sPrefix         = "";
  m_a2nDim[0]         = 0;
  m_a2nDim[1]         = 0;
  for (i = 0; i < 3; i++)
    {
      m_a2x3fVector[0][i] = 0.0;
      m_a2x3fVector[1][i] = 0.0;
      m_a3fDN[i]          = 0.0;
      m_a3fLocalDN[i] = 0.0;
  }
  vIdentMat3D(m_a3x3fDD);
  vIdentMat3D(m_a3x3fDDRot);

  m_a2x3fVector[0][0]   = 1.0;
  m_a2x3fVector[1][1]   = 1.0;

  m_a2fSize[0]        = 0.0;
  m_a2fSize[1]        = 0.0;

  m_poGoniometer    = NULL;
  m_poSpatial       = NULL;
  m_poNonunf        = NULL;
  m_bNewGonio       = FALSE;
  m_pfPixelShiftArray = NULL;
  m_a2nPixelShiftArrayDiv[0] = 0;
  m_a2nPixelShiftArrayDiv[1] = 0;
  m_pfProfileArray = NULL;
  m_a2nProfileArrayDiv[0] = 0;
  m_a2nProfileArrayDiv[1] = 0;

  return (0);
}

int
Cdetector::nInitValues(Cimage_header& oHeader, const Cstring& sPre,
	   const bool bSpatial, const bool bNonunf)
{
  (void) nInitValues();

  return (nUpdateFromHeader(oHeader, sPre, bSpatial, bNonunf));
}

int
Cdetector::nList(const int nFlag)
{
  cout << "\n" << m_sPrefix << "Detector listing: \n\n";

  cout << "   Pixel dimensions:     " << m_a2nDim[0] << ", "
                                      << m_a2nDim[1] << endl;
  cout << " Nominal size in mm:     " << m_a2fSize[0] << ", "
                                      << m_a2fSize[1] << endl;
  cout << "        Description:     " << m_sDescription << endl;

  if (1 > nFlag)
    {

      cout << "  Fast pixel vector: " << m_a2x3fVector[0][0] << ", "
                                      << m_a2x3fVector[0][1] << ", "
                                      << m_a2x3fVector[0][2] << endl;

      cout << "  Slow pixel vector: " << m_a2x3fVector[1][0] << ", "
                                      << m_a2x3fVector[1][1] << ", "
                                     << m_a2x3fVector[1][2] << endl;
    }
  if (NULL != m_poSpatial)
    {
      m_poSpatial->nList();
/*
      float a2fBeam[2];
      m_poSpatial->nGetBeamPosition(&a2fBeam[0], &a2fBeam[1]);
      cout << " Direct beam posn in px: " << a2fBeam[0] 
                                  << ", " << a2fBeam[1] << endl;
*/
    }
  cout << flush;
  m_poGoniometer->nList(nFlag);
  return (0);
}


int 
Cdetector::nSetStandardDetectorParams(double a6fParams[6],int nStart,int nCount) {
    float a6fOrigParams[6];
    int nx;
    m_poGoniometer->nGetDatum(6,a6fOrigParams);
    for (nx=nStart;nx<nStart+nCount;nx++) {
        a6fOrigParams[nx] = a6fParams[nx];
    };
    return m_poGoniometer->nSetDatum(6,a6fOrigParams);
};

int Cdetector::nSetModifiedDetectorParams(double a6fParams[6],int nStart,int nCount)
{
    double a3x3fRotMat[3][3];
    double a3fVec0[3];
    double a3fVec1[3];
    double a6fOrigParams[6];
    float  a3fFloatBuf[3];
    int nx;


    nGetModifiedDetectorParams(a6fOrigParams);   
    for (nx=nStart;nx<nStart+nCount;nx++) {
        a6fOrigParams[nx] = a6fParams[nx];
    };
   
    // Set the rotations.  These will remain the same.
    vCopyVec3D(a6fOrigParams,a3fFloatBuf);
    m_poGoniometer->nSetDatum(3,a3fFloatBuf);

    // Get the inverse rotation
    nCalcDDDNRotMat(&a3x3fRotMat[0][0]);
    vTranMat3D(a3x3fRotMat);

    // Back-translate the modified translation vector by multiplying by the inverse rotation matrix.
    vCopyVec3D(&a6fOrigParams[3],a3fVec0);
    vMulMat3DVec3D(a3x3fRotMat,a3fVec0,a3fVec1);    
    vCopyVec3D(a3fVec1,a3fFloatBuf); 

    // Represent the translated vector in terms of the goniometer basis vectors using nSetTransVector().
    m_poGoniometer->nSetTransVector(3,a3fFloatBuf);
    
    return 0;
};

int 
Cdetector::nGetStandardDetectorParams(double a6fParams[6],int nStart,int nCount) {
    float a6fOrigParams[6];
    int nx;
    m_poGoniometer->nGetDatum(6,a6fOrigParams);
    for (nx=nStart;nx<nStart+nCount;nx++) {
        a6fParams[nx] = a6fOrigParams[nx];
    };
    return 0;
};

int 
Cdetector::nGetModifiedDetectorParams(double a6fParams[6],int nStart,int nCount) {
    float a6fOrigParams[6];
    double a3x3fDD[3][3];
    double a3fDN[3];
    int nx;

    // The first three modified parameters are exactly the same.
    m_poGoniometer->nGetDatum(6,a6fOrigParams);
    // The three trans values are obtained from the last
    // row of the DDDN matrix.
    nCalcGetDDDN(&a3x3fDD[0][0],&a3fDN[0]);
    for (nx=0;nx<3;nx++)
        a6fOrigParams[nx+3] = a3x3fDD[2][nx];

    for (nx=nStart;nx<nStart+nCount;nx++) {
        a6fParams[nx] = a6fOrigParams[nx];
    };
    return 0;
};



void
Cdetector::vUpdateDDDN(void)
{
  if (0 != nCalcGetDDDN(&m_a3x3fDD[0][0], m_a3fDN))
    {
      vIdentMat3D(m_a3x3fDD);
      m_a3fDN[0] = 0.0;
      m_a3fDN[1] = 0.0;
      m_a3fDN[2] = 0.0;
    }
}
/////////////////////////////////////////////////////////////////////
int Cdetector::nCalcDDDNRotMat(double* pfRotMat,
                               double* pfRotMatDeriv0,
                               double* pfRotMatDeriv1,
                               double* pfRotMatDeriv2) 
{
    double  fTemp1[3];
    double  fTemp2[3*3];
    
    float   fTemp1_[3];
    float   fTemp2_[3*3];
    
    int     nx = 0;

    if( m_poGoniometer->bIsAvailable() )
    {
        m_poGoniometer->nGetDatum(3, &fTemp1_[0]);
        
        m_poGoniometer->nGetVectors(3, &fTemp2_[0]);
        
        for(nx = 0; nx < 3; nx++)
            fTemp1[nx] = fTemp1_[nx];
        
        for(nx = 0; nx < 9; nx++)
            fTemp2[nx] = fTemp2_[nx];
        
        vConv3Rot3Vec3DMat3D(fTemp1[0], 
                             fTemp1[1], 
                             fTemp1[2],
                            &fTemp2[0], 
                            &fTemp2[3], 
                            &fTemp2[6],
                             pfRotMat,
                             pfRotMatDeriv0,
                             pfRotMatDeriv1,
                             pfRotMatDeriv2);
        
        vCopyMat3D(pfRotMat, &m_a3x3fDDRot[0][0]);
    
        return 0;
    }
    else
    {
        return 1;
    }
}
//////////////////////////////////////////////////////////////////////
int Cdetector::nCalcGetDDDNDerivStandard(double a6x3x3fDeriv[6][3][3])
{
    double a3x3fRotMat[3][3];
    double a3x3fTemp0[3][3];
    double a3x3x3fRotDeriv[3][3][3];
    float  fDgds_[3];
    int nx;
    

   
    (void) m_poGoniometer->nCalcGetTransVector(fDgds_); // in MADNES DGDS is
    
    // For the three rotation parameters, we compute the derivative for the rotation matrix.
    nCalcDDDNRotMat(&a3x3fRotMat[0][0],&a3x3x3fRotDeriv[0][0][0],&a3x3x3fRotDeriv[1][0][0],&a3x3x3fRotDeriv[2][0][0]);
    vCopyVec3D(&m_a2x3fVector[0][0],&a3x3fTemp0[0][0]);
    vCopyVec3D(&m_a2x3fVector[1][0],&a3x3fTemp0[1][0]);
    vCopyVec3D(&fDgds_[0],&a3x3fTemp0[2][0]);
    for (nx=0;nx<3;nx++) {
        vMulMat3DMat3D(a3x3x3fRotDeriv[nx],a3x3fTemp0,a6x3x3fDeriv[nx]);        
    };
    // For the three translation values, we simply multiply the appropriate 
    // derivative vector by the rotation matrix.

    vZeroMat3D(&a3x3fTemp0[0][0]);
    for (nx=3;nx<6;nx++) {
        m_poGoniometer->nCalcGetTransVectorDeriv(nx-3,&a3x3fTemp0[2][0]);
        vMulMat3DMat3D(a3x3fRotMat,a3x3fTemp0,a6x3x3fDeriv[nx]);
    };
    return 0;
};

int
Cdetector::nCalcGetDDDNDerivModified(double a6x3x3fDeriv[6][3][3]) {
    int nx,ny;

    nCalcGetDDDNDerivStandard(a6x3x3fDeriv);

    // The rotation derivatives are the same, except for the last col. of each
    // which is set to zero.
    for (nx=0;nx<3;nx++) {
        for (ny=0;ny<3;ny++)
            a6x3x3fDeriv[nx][2][ny] = 0.0;
    };
    // The translation derivatives are easy:  They are set to the correct matrix elements.
    for (nx=3;nx<6;nx++) {
        vZeroMat3D(&a6x3x3fDeriv[nx][0][0]);
        a6x3x3fDeriv[nx][2][nx-3] = 1.0;
    };
    return 0;
};


int
Cdetector::nCalcGetDDDN(double *pfDD, double *pfDN)
{

// Calculates and returns the 3x3 matrix *pfDD and the 3D vector *pfDN
// which describes the detector position 0.0, 0.0 mm in the
// laboratory coordinate system AFTER all rotations and translations
// of the detector goniostat have been applied.
// See "Computational Aspects of Protein Crystal Data Analysis", DL/SCI/R25
// (the Daresbury Study Weekend 23-24 January 1987) pp. 126.... and
// Proceedings of the EEC Cooperative Workshop on Position-Sensitive Detector
// Software (Phases I& II) p.104... and (Phase III) pp. 57-64.
//
//  pfDD = | d       d       d   |
//         |  xx      yx      ox |
//         |                     |
//         | d       d       d   |
//         |  xy      yy      oy |
//         |                     |
//         | d       d       d   |
//         |  xz      yz      oz |
//
//  dx = gon_rot * fVector[0][*]
//  dy = gon_rot * fVector[1][*]
//  do = gon_rot * dettranslation (that is parallel to source)
//
  double  fDx[3];
  double  fDy[3];
  double  fDo[3];

  double  fPo[3];
  double  fDgds[3];
  float   fDgds_[3];
  double  fMatrix[3][3];
  double  m_a2x3fVector_[2][3];
  int nx,ny;

  (void) m_poGoniometer->nCalcGetTransVector(fDgds_); // in MADNES DGDS is

  for (nx=0;nx<3;nx++)
      fDgds[nx] = fDgds_[nx];
  for (nx=0;nx<2;nx++) {
      for (ny=0;ny<3;ny++)
        m_a2x3fVector_[nx][ny] = m_a2x3fVector[nx][ny];
  };

// Shouldn't fPo be defined as the vector that is parallel to source?
// or perpendicular to fVector[0][*] and fVector[1][*]?
// It IS the crystal to detector distance, but doesn't have to be along
// laboratory Z (if the rest of this is correct!)

  fPo[0] = 0.0;
  fPo[1] = 0.0;
  fPo[2] = fDgds[2];  // This is different than MADNES!, just had -dist
//  fPo[2] = -fLenVec3D(fDgds);  // Shouldn't it be this?

  // Explicitly form the rotation matrix since, the order we want is not
  // what the goniometer class gives us.

//  m_poGoniometer->vCalcGetRotMatrix((double*)fMatrix);

  if (0 == nCalcDDDNRotMat(&fMatrix[0][0]))
    {
      vMulMat3DVec3D(fMatrix, &m_a2x3fVector_[0][0], fDx);
      vMulMat3DVec3D(fMatrix, &m_a2x3fVector_[1][0], fDy);
      vMulMat3DVec3D(fMatrix,  fDgds, fDo);
      for (int i = 0; i < 3; i++)
	{
	  pfDD[i]   = fDx[i];
	  pfDD[3+i] = fDy[i];
	  pfDD[6+i] = fDo[i];
	}
      vMulMat3DVec3D(fMatrix,  fPo,    pfDN);
      return (0);
    }
  else
    {
      return (1);
    }
}



int
Cdetector::nCalcDDDNRotMat(float* pfRotMat,float* pfRotMatDeriv0,float* pfRotMatDeriv1,float* pfRotMatDeriv2) 
{
  float fTemp1[3];
  float fTemp2[3*3];

  if (m_poGoniometer->bIsAvailable())
    {
      (void) m_poGoniometer->nGetDatum(3, &fTemp1[0]);
      (void) m_poGoniometer->nGetVectors(3, &fTemp2[0]);
      vConv3Rot3Vec3DMat3D(fTemp1[0], fTemp1[1], fTemp1[2],
			   &fTemp2[0], &fTemp2[3], &fTemp2[6],
			   (float *)pfRotMat,
               pfRotMatDeriv0,pfRotMatDeriv1,pfRotMatDeriv2);
      return (0);
    }
  else
    {
      return (1);
    }

}


int
Cdetector::nCalcGetDDDN(float *pfDD, float *pfDN)
{

// Calculates and returns the 3x3 matrix *pfDD and the 3D vector *pfDN
// which describes the detector position 0.0, 0.0 mm in the
// laboratory coordinate system AFTER all rotations and translations
// of the detector goniostat have been applied.
// See "Computational Aspects of Protein Crystal Data Analysis", DL/SCI/R25
// (the Daresbury Study Weekend 23-24 January 1987) pp. 126.... and
// Proceedings of the EEC Cooperative Workshop on Position-Sensitive Detector
// Software (Phases I& II) p.104... and (Phase III) pp. 57-64.
//
//  pfDD = | d       d       d   |
//         |  xx      yx      ox |
//         |                     |
//         | d       d       d   |
//         |  xy      yy      oy |
//         |                     |
//         | d       d       d   |
//         |  xz      yz      oz |
//
//  dx = gon_rot * fVector[0][*]
//  dy = gon_rot * fVector[1][*]
//  do = gon_rot * dettranslation (that is parallel to source)
//
  float  fDx[3];
  float  fDy[3];
  float  fDo[3];

  float  fPo[3];
  float  fDgds[3];
  float  fMatrix[3][3];

  (void) m_poGoniometer->nCalcGetTransVector(fDgds); // in MADNES DGDS is

// Shouldn't fPo be defined as the vector that is parallel to source?
// or perpendicular to fVector[0][*] and fVector[1][*]?
// It IS the crystal to detector distance, but doesn't have to be along
// laboratory Z (if the rest of this is correct!)

  fPo[0] = 0.0;
  fPo[1] = 0.0;
  fPo[2] = fDgds[2];  // This is different than MADNES!, just had -dist
//  fPo[2] = -fLenVec3D(fDgds);  // Shouldn't it be this?

  // Explicitly form the rotation matrix since, the order we want is not
  // what the goniometer class gives us.

//  m_poGoniometer->vCalcGetRotMatrix((float*)fMatrix);

  if (0 == nCalcDDDNRotMat(&fMatrix[0][0]))
    {
      vMulMat3DVec3D(fMatrix, &m_a2x3fVector[0][0], fDx);
      vMulMat3DVec3D(fMatrix, &m_a2x3fVector[1][0], fDy);
      vMulMat3DVec3D(fMatrix,  fDgds, fDo);
      for (int i = 0; i < 3; i++)
	{
	  pfDD[i]   = fDx[i];
	  pfDD[3+i] = fDy[i];
	  pfDD[6+i] = fDo[i];
	}
      vMulMat3DVec3D(fMatrix,  fPo,    pfDN);
      return (0);
    }
  else
    {
      return (1);
    }
}


double
Cdetector::fCalcGetResolution(const double fPx0, const double fPx1,
			      const float _fS0[3], double *pfX,double *pfChi)
{
  // Calculate the resolution in Angstroms of the given detector pixel.
  // For this to work, the routine Cdetector::nCalcGetDDDN must have been
  // called and worked so that the member variables m_a3x3fDD, m_a3fDN are valid

  // Get mm coordinates of this pixel, transform it by
  // the detector goniometer system to get the scattered
  // beam wavevector for a reflection at that pixel,
  // take the dot product with the source vector to get
  // the cosine of 2theta, then compute resolution from
  // Bragg's Law.

  int   nStat;
  float  fVDC[3];     // Virtual Detector Coordinate
  double fLenS0;      // Length of source vector (should be wavelength or 1.0)
  double fReso;       // Resolution in Angstroms
  float  fS_[3];      // Scattered beam wavevector
  double fS[3];       
  double fS0[3];

  vCopyVec3D(_fS0,fS0);

  fLenS0  = fLenVec3D(fS0);                 // This may also be the wavelength.
  nStat   = m_poSpatial->nPixeltoMM(fPx0, fPx1, &fVDC[0], &fVDC[1], &fVDC[2]);
  if (0 == nStat)
    {
      // If on the detector then make the calculations.
      //  x = s + s0   Reso = 1 / |x|
      //  -   -   --               -

      vMulMat3DVec3D(m_a3x3fDD, fVDC, fS_);  // Calc scattered beam wavevector
      vCopyVec3D(fS_,fS);
      (void) fNormVec3D(fS);                 // Normalize the scattered beam
      vCopyVec3D(fS,m_a3fS);

      if (pfChi) {
          double a3fBeamVec[3];
          double a3fPerpBeamVec[3];
          double a3fPerpPerpBeamVec[3];
          double a3fTempVec1[3];
          double f0,f1,fLen; 

          vCopyVec3D(fS0,a3fBeamVec);
          fNormVec3D(a3fBeamVec);

          // Take x vector (points down) and project onto the S0 vector.
          a3fPerpBeamVec[0] = 1.0;
          a3fPerpBeamVec[1] = 0.0;
          a3fPerpBeamVec[2] = 0.0;
          // Subtract projection of a3fPerpBeamVec on a3fBeamVec (i.e. orthonormalize)
          f1 = fDot3D(a3fPerpBeamVec,a3fBeamVec);
          vMulVec3DScalar(a3fBeamVec,f1,a3fTempVec1);
          vSubVec3DVec3D(a3fPerpBeamVec,a3fTempVec1,a3fPerpBeamVec);
          fNormVec3D(a3fPerpBeamVec);
          // Create other perpendicular vector. Now we should have three orthonormal vectors.
          vCross3D(a3fBeamVec,a3fPerpBeamVec,a3fPerpPerpBeamVec);

          // Project vector down to plane containing Perp vector and PerpPerp vector.
          f1 = fDot3D(fS,a3fBeamVec);
          vMulVec3DScalar(a3fBeamVec,f1,a3fTempVec1);
          vSubVec3DVec3D(fS,a3fTempVec1,a3fTempVec1);
          fLen = fLenVec3D(a3fTempVec1);

          f1 = fDot3D(a3fTempVec1,a3fPerpBeamVec)/max(1e-20,fLen);
          if (f1>=1.0) 
              f0=0.0; 
          else if (f1<=-1.0) 
              f0=Gs_dPI; 
          else 
              f0=acos(f1);
          
          if (fDot3D(a3fTempVec1,a3fPerpPerpBeamVec)<0.0) 
              *pfChi=2.0*Gs_dPI-f0; 
          else 
              *pfChi=f0;

      };

      vMulVec3DScalar(fS, fLenS0, fS);       // wavevec. to be same length as S0
      vAddVec3DVec3D(fS, fS0, fS);           // S + S0 = X
      if (NULL != pfX)
	{
	  // Return the reciprocal lattice vector if asked for
	  pfX[0] = fS[0];
	  pfX[1] = fS[1];
	  pfX[2] = fS[2];
	}
      fReso =  fNormVec3D(fS);               // 1 / |X| is resolution
      if (0.0 == fReso)                      // but prevent divide by zero
	{
	  fReso = DTREK_DET_RESO_MAX;                   // This is close enough to infinity
	}
      else
	{
	  fReso = 1.0f / fReso;               // The answer
	}
      return (fReso);
    }
  else
    {
      return (-1.0);                         // Error
    }
}


int    
Cdetector::nCalcGetPixelFromThetaChi(double f2Theta,double fChi,double& fPix0,double& fPix1,
                                       const float fS0[3],bool bTruncateOutOfBounds) {

    float  _fPix0,_fPix1;
    double f0,f1; 
    double a3fTempVec1[3];
    double a3fTempVec2[3];


    double a3fBeamVec[3];
    double a3fPerpBeamVec[3];
    double a3fPerpPerpBeamVec[3];
    double fS[3];
    
    bool bGreaterThan90;
    
    // Should not need this.  Make powder software call this function.
    // (void) nCalcGetDDDN(&m_a3x3fDD[0][0], m_a3fDN);    
    
    // Take x vector (points down) and project onto the S0 vector.
    a3fPerpBeamVec[0] = 1.0;
    a3fPerpBeamVec[1] = 0.0;
    a3fPerpBeamVec[2] = 0.0;
    vCopyVec3D(fS0,a3fBeamVec);
    fNormVec3D(a3fBeamVec);
    f1 = fDot3D(a3fPerpBeamVec,a3fBeamVec);
    
    vMulVec3DScalar(a3fBeamVec,f1,a3fTempVec1);
    // Subtract projection.
    vSubVec3DVec3D(a3fPerpBeamVec,a3fTempVec1,a3fPerpBeamVec);
    // Create other perpendicular vector.
    vCross3D(a3fBeamVec,a3fPerpBeamVec,a3fPerpPerpBeamVec);

   // Compute the length of the beam vector using the 2theta value.
   if (f2Theta>=Gs_dPI*0.5) {
      f0=ABS(tan(Gs_dPI-f2Theta));
      bGreaterThan90=TRUE;
   } else {
      f0=ABS(tan(f2Theta));
      bGreaterThan90=FALSE;
   };
   // Correction for 2theta near 90 degrees.
   if (f0>1e20) f0=1e20;   

   // Use vector parameterization to figure out the direction.  We are using <cos(chi).sin(chi)> as a 
   // parameterization tuple.

   vMulVec3DScalar(a3fPerpBeamVec,f0*cos(fChi),a3fTempVec1);
   vMulVec3DScalar(a3fPerpPerpBeamVec,f0*sin(fChi),a3fTempVec2);
   vAddVec3DVec3D(a3fTempVec1,a3fTempVec2,fS);
   if (bGreaterThan90) 
      vAddVec3DVec3D(fS,a3fBeamVec,fS);
   else
      vSubVec3DVec3D(fS,a3fBeamVec,fS);

   double a3fVTS[3];
   double a3fVDCC[3];
   double a3x3fDDInv[3][3];
   double a3x3fDD[3][3];

   if( 0 != nScaleSOnPlate(fS,a3fVTS) )
       return 1;
   
   // Calculate the virtual detector coordinates by using the inverse DD matrix.
   vCopyMat3D(&m_a3x3fDD[0][0],&a3x3fDD[0][0]);
   fInvMat3D(&a3x3fDD[0][0],&a3x3fDDInv[0][0]);
   
   vMulMat3DVec3D(a3x3fDDInv, a3fVTS, a3fVDCC);
   if( 0 != m_poSpatial->nMMtoPixel(a3fVDCC[0], a3fVDCC[1],a3fVDCC[2],&_fPix0, &_fPix1,bTruncateOutOfBounds) )
       return 1;
   
   fPix0 = _fPix0;
   fPix1 = _fPix1;

   return 0;
};


int Cdetector::nGetResolution(const float fS0[3], 
                              int nBorder,
                              int nStep,
                              float fMinFraction,
                              float *pfResoMin,
                              float *pfResoMax,
                              int nChiTestArraySize,
                              double* pfResolutions,
                              double* pfChiCoverageArray)
{
    int i,j;
    int nx,ny;
    int nLocalChiMin,nLocalChiMax;
    double fReso,fLocalResoMax,fLocalResoMin;
    double fChi,fLocalChiMax,fLocalChiMin;
    double a2x2fReso[2][2];
    double a2x2fChi[2][2];
    const double fMinBinReso = 10.0;
    const double fStepBinReso = 0.01;
    const int    nBinReso = 1000;
    itr<char> aacResoBins[360];
    itr<double> afResoBins;


    (void) nCalcGetDDDN(&m_a3x3fDD[0][0], m_a3fDN);    

    *pfResoMax = DTREK_DET_RESO_MAX;
    *pfResoMin = 0.0;    
    
    for (ny = 0; ny < 360; ny++) {
        aacResoBins[ny].setsize(nBinReso);
        for (nx=0;nx<nBinReso;nx++)
            aacResoBins[ny][nx] = (char) 0;
    };
    afResoBins.setsize(nBinReso);
    

    for (i = nBorder; (i < m_a2nDim[0] - nBorder) && (i + nStep < m_a2nDim[0]); i += nStep)  {
        for (j = nBorder; (j < m_a2nDim[1] - nBorder) && (j + nStep < m_a2nDim[1]); j += nStep)  {
            
            a2x2fReso[0][0] = fCalcGetResolution((float)i, (float) j, fS0,NULL,&a2x2fChi[0][0]);
            a2x2fReso[0][1] = fCalcGetResolution((float)i, (float) j + nStep, fS0,NULL,&a2x2fChi[0][1]);
            a2x2fReso[1][0] = fCalcGetResolution((float)i + nStep, (float) j, fS0,NULL,&a2x2fChi[1][0]);
            a2x2fReso[1][1] = fCalcGetResolution((float)i + nStep, (float) j + nStep, fS0,NULL,&a2x2fChi[1][1]);            
            if ((m_poNonunf) && 
                (m_poNonunf->bPixelIsBad(i,j) || m_poNonunf->bPixelIsBad(i + nStep,j) || m_poNonunf->bPixelIsBad(i,j + nStep) || m_poNonunf->bPixelIsBad(i + nStep,j + nStep)))
                continue;
            
            fLocalResoMax = DTREK_DET_RESO_MAX;
            fLocalResoMin = 0.0;    
            fLocalChiMax = 0.0;
            fLocalChiMin = Gs_dPI*2.0;

            for (nx=0;nx<2;nx++) {
                for (ny=0;ny<2;ny++) {                    
                    fReso = a2x2fReso[nx][ny];
                    fChi = a2x2fChi[nx][ny];

                    if (fReso > fLocalResoMin) fLocalResoMin = fReso;
                    if (fReso < fLocalResoMax) fLocalResoMax = fReso;
                    if (fChi > fLocalChiMax) fLocalChiMax = fChi;
                    if (fChi < fLocalChiMin) fLocalChiMin = fChi;
                };
            };

            if (fLocalChiMax - fLocalChiMin > Gs_dRADIANS_PER_DEGREE*180.0) {
                fLocalChiMin += Gs_dPI*2.0;
                std::swap(fLocalChiMin,fLocalChiMax);
            };
            nLocalChiMin = max(0,(int) floor(fLocalChiMin / Gs_dRADIANS_PER_DEGREE));
            nLocalChiMax = min(360-1,(int) ceil(fLocalChiMax / Gs_dRADIANS_PER_DEGREE));
            
            for (nx = (int) (fLocalResoMax/fStepBinReso);nx<= min(nBinReso - 1,(int) (fLocalResoMin/fStepBinReso));nx++) {
                for (ny = nLocalChiMin; ny <= nLocalChiMax;ny++) 
                    aacResoBins[ny][nx] = (char) 1;
            };
        };
    };

    for (nx = 0; nx < nBinReso; nx++) {
        afResoBins[nx] = 0.0;
        for (ny = 0; ny < 360; ny++) {
            if (aacResoBins[ny][nx])
                afResoBins[nx] += 1.0;
        };
    };


    // Look for top contributor.
    int nMaxBin = 0;
    double fAvgReso = 0.0;
    int    nAvgReso = 0;
    for (nx=0;nx<nBinReso;nx++) {
        if (afResoBins[nx]>afResoBins[nMaxBin])
            nMaxBin = nx;
        if (afResoBins[nx]>0.0) {
            fAvgReso += afResoBins[nx];
            nAvgReso++;
        };
    };
    fAvgReso/=max(1,nAvgReso);

    // Choose resolution cutoffs where we have a minimum fraction of coverage.

    // Find high resolution by starting at beginning of array.
    for (nx=0;(nx<nBinReso) && (afResoBins[nx]<fAvgReso*fMinFraction);nx++);
    if (nx==nBinReso)
        return 1;
    *pfResoMax = nx*fStepBinReso;

    // Find low resolution by starting at end of array.      

    for (nx=nBinReso-1;(nx>=0) && (afResoBins[nx]<fAvgReso*fMinFraction);nx--);
    if (nx==0) 
        return 1;
    *pfResoMin = nx*fStepBinReso;

    
	
    if (0)
	{
		FILE* pFData;
		static int nPrintNumber = 0;
		Cstring sPrint;
		sPrint = "reso";
		sPrint += nPrintNumber;
		sPrint += ".txt";
		pFData = fopen(sPrint.string(),"wt");
		for (nx=0;nx<nBinReso;nx++) 
			fprintf(pFData,"%f\n",afResoBins[nx]);
		fclose(pFData);
		nPrintNumber++;
	};
	

    // Write out user specified chi coverages.
    for (nx=0; nx< nChiTestArraySize;nx++) {
        int n1,n2;
        double f0;

        f0 = pfResolutions[nx]/fStepBinReso;
        n1 = max(0,min(nBinReso-1,(int) floor(f0)));
        n2 = max(0,min(nBinReso-1,(int) ceil(f0)));
        pfChiCoverageArray[nx] = ABS(n1 - f0)*afResoBins[n2] + ABS(n2 - f0)*afResoBins[n1];
    };
    
   


    float fDum1, fDum2, fDum3, fDum4, fTemp[3];
    int   nStat;
    
    for (i = 0; i < 3; i++)
    {
        fTemp[i] = 0.0;
    }
    
    // Check to see that the beam hits the image plate.
    nStat = nCalcDetCoords(fTemp, fS0, &fDum1, &fDum2, &fDum3, &fDum4);
    
    if (0 == nStat)
    {
        *pfResoMin = DTREK_DET_RESO_MAX;
    }

    return 0;
};



int
Cdetector::nGetResolution(const float fS0[3], 
                          float *pfResoMin,
			              float *pfResoMax,
                          float *pfResoMaxEdge)
{
  // fS0[3] is the vector describing the direct beam from Csource
  // If fS0 is normalized, then the returned resolutions need to be
  // multiplied by fWavelength.
  //
  // Maximum resolution has smallest d-plane spacing,
  // i.e. 1 Angstrom is more max than 3.0 Angstrom.

  // Compute the minimum and maximum resolution on the detector in Angstroms
  //
  // +----+----+   Look at only a few pixels to determine
  // |.   .   .|   their mm positions and then their recip lattice vectors.
  // | .  .  . |   Resolution is from 1/length of reciprocal lattice vector x:
  // |  . . .  |
  // |   ...   |
  // +....+....+   x =  s0 + s          Resolution = fWavelength / |x|
  // |   ...   |   -    --   -
  // |  . . .  |
  // | .  .  . |   (Remember in this program s0 goes from xtal towards source!)
  // |.   .   .|
  // +----+----+
  // Look at pixels around the 4 edges of the detector.

  int   i, j;               // Loop counters
  float fReso;
  float a4fMinResoEdge[4];
  float a4fMaxResoEdge[4];

  *pfResoMax = DTREK_DET_RESO_MAX;
  *pfResoMin = 0.0;
  for (i = 0; i < 4; i++) {
      a4fMinResoEdge[i] = 0.0;
      a4fMaxResoEdge[i] = DTREK_DET_RESO_MAX;
  };

  if (!m_poGoniometer->bIsAvailable())
    {
      // No detector goniometer available, so we cannot calculate resolution.

      return (1);
    }

  (void) nCalcGetDDDN(&m_a3x3fDD[0][0], m_a3fDN);

  // Look at the 4 edges of the detector.

  for (j = 0; (j < m_a2nDim[1]); j++)
  {
      // Look along the j = 0 edge, near i = 0
      
      fReso = -1.0;
      for (i = 0; (i < m_a2nDim[0]) && (fReso < 0.0); i++)
      {
          fReso = fCalcGetResolution((float)i, (float)j, fS0);
          if (fReso > 0.0)
          {
              if (fReso > a4fMinResoEdge[0])
                  a4fMinResoEdge[0] = fReso;
              if (fReso < a4fMaxResoEdge[0])
                  a4fMaxResoEdge[0] = fReso;
          }
          
      }
      
      // Look along the j = 0 edge, near i = m_a2nDim[0]
      
      fReso = -1.0;
      for (i = m_a2nDim[0]-1; (i >= 0) && (fReso < 0.0); i--)
      {
          fReso = fCalcGetResolution((float)i, (float)j, fS0);
          if (fReso > 0.0)
          {
              if (fReso > a4fMinResoEdge[1])
                  a4fMinResoEdge[1] = fReso;
              if (fReso < a4fMaxResoEdge[1])
                  a4fMaxResoEdge[1] = fReso;
          }
      }
  }

  // Next two edges joined at m_a2nDim[0], m_a2nDim[1], see previous loops
  // for comments

  for (i = 0; i < m_a2nDim[0]; i++)
  {
      // Look along the i = 0 edge, near j = 0
      
      fReso = -1.0;
      for (j = 0; (j < m_a2nDim[1]) && (fReso < 0.0); j++)
      {
          fReso = fCalcGetResolution((float)i, (float)j, fS0);
          if (fReso > 0.0)
          {
              if (fReso > a4fMinResoEdge[2])
                  a4fMinResoEdge[2] = fReso;
              if (fReso < a4fMaxResoEdge[2])
                  a4fMaxResoEdge[2] = fReso;
          }
          
      }
      
      // Look along the i = 0 edge, near j = m_a2nDim[1]
      
      fReso = -1.0;
      for (j = m_a2nDim[1]-1; (j >= 0) && (fReso < 0.0); j--)
      {
          fReso = fCalcGetResolution((float)i, (float)j, fS0);
          if (0.0 < fReso)
          {
              if (fReso > a4fMinResoEdge[3])
                  a4fMinResoEdge[3] = fReso;
              if (fReso < a4fMaxResoEdge[3])
                  a4fMaxResoEdge[3] = fReso;
          }
      }
  }

  // Check if direct beam is on detector, if so, then ResoMin is DTREK_DET_RESOMAX
  // fDum1-4 are dummy variables

  float fDum1, fDum2, fDum3, fDum4, fTemp[3];
  int   nStat;

  for (i = 0; i < 3; i++)
    {
      fTemp[i] = -fS0[i];
    }
  nStat = nCalcDetCoords(fTemp, fS0, &fDum1, &fDum2, &fDum3, &fDum4);

  *pfResoMin = max(max(a4fMinResoEdge[0],a4fMinResoEdge[1]),max(a4fMinResoEdge[2],a4fMinResoEdge[3]));
  *pfResoMax = min(min(a4fMaxResoEdge[0],a4fMaxResoEdge[1]),min(a4fMaxResoEdge[2],a4fMaxResoEdge[3]));

  if (0 == nStat)
    {
      *pfResoMin = DTREK_DET_RESO_MAX;
    }

  if (NULL != pfResoMaxEdge) {
      *pfResoMaxEdge = min(min(a4fMinResoEdge[0],a4fMinResoEdge[1]),min(a4fMinResoEdge[2],a4fMinResoEdge[3]));
  };
  if (*pfResoMin < *pfResoMax) return (1);  // Didn't change them, so error

  return (0);
}

int Cdetector::nScaleSOnPlateItr(double fS[3],double fTS[3],double fShift[3])
{
    float _fS[3];
    float _fTS[3];
    float _fShift[3];

    int nStat;
    
    vCopyVec3D(fS,_fS);
 
    vCopyVec3D(fShift,_fShift);
    
    nStat = nScaleSOnPlateItr(_fS,_fTS,_fShift);
    
    vCopyVec3D(_fTS,fTS);
    
    return nStat;
}

int Cdetector::nScaleSOnPlateItr(float fS[3], float fTS[3], float fShift[3])
{
	float a3fTempVec[3];
	float f0;
    int nStat;

    nStat = nScaleSOnPlate(fS,fTS);
    
    if( nStat )
        return nStat;
	
    f0 = fDot3D(m_a3fLocalDN, fShift) / fDot3D(m_a3fLocalDN,fTS);
	
    vMulVec3DScalar(fTS, f0, a3fTempVec);
	
    vSubVec3DVec3D(a3fTempVec, fShift, a3fTempVec);
	
    vAddVec3DVec3D(a3fTempVec, fTS, fS);
	
    fNormVec3D(fS);
    
    nStat = nScaleSOnPlate(fS, fTS);

    if( nStat )
        return nStat;
   
    return 0;
}


int Cdetector::nScaleSOnPlate(const double fS[3],double fTS[3])
{
    float _fS[3];
    float _fTS[3];

    vCopyVec3D(fS,_fS);
    
    int     nStat = nScaleSOnPlate(_fS,_fTS);
    
    vCopyVec3D(_fTS,fTS);
    
    return nStat;
}

////////////////////////////////////////////////////////////////////////////////////////
// Scale the scattered vector S on the detector surface using equation
//
// T S = (XR - S0) * DIST**2 / (XR.DN - S0.DN)   
//   -    --   --               -- --   -- --
// where T is a scalar, XR is a reciprocal-lattice vector, S0 is the primary beam vector,
// DIST is the crystal-to-detector distance, DN is the normal to detector vector.
//
// See Proceedings Phase III, pp. 57-64

// For the cylindrical detector the equation for determining T is
//
//   T**2 * (Sy**2 + Sz**2) = Rc**2, 
//   where Rc is the cylinder radius.
// 
// If the cylidrical detector is shifted by a "shift" along Z-axis, the equation becomes:
//
//   T**2 * (Sy**2 + Sz**2) + (2 * shift * T * Sz) + (shift**2 - Rc**2) = 0
//
// Notice: for the cylindrical detector, the length of the local normal vector, m_a3fLocalDN, 
// in this routine is Rc/T, regardless whether there is a shift or not.
int Cdetector::nScaleSOnPlate(const float fS[3], float fTS[3])
{
    double      fDetDistance = 0.0;
    double      fRMFAC = 0.0;
    double      fDNdotS = 0.0;
    double      f0 = 0.0;
    float       a3x3fRotMatInv[3][3];
    
    float       fSRot[3];        // The fS vector inverse rotated to account for plate rotation.
                                 // This is only done on the RAPID, where m_a3fDN is not relevant.

    fDetDistance = (double)fGetDistance();
    
    if( m_poSpatial->eSpatialType() == eSpatial_flat_geometry )
    {
        fDNdotS = fDot3D(fS,m_a3fDN);
	
        // Get the local detector normal (it's just equal to m_a3fDN)
		vCopyVec3D(m_a3fDN, m_a3fLocalDN);

        if (fDNdotS == 0.0)
            return 1;
        
        fRMFAC = (fDetDistance * fDetDistance) / (fDNdotS);
    } 
    else if( m_poSpatial->eSpatialType() == eSpatial_cylindrical_geometry )
    {
        //RB: We will treat the shifted cylinder case separately from Thad's code
        // for the non-shifted cylinder
        float       fShift = m_poSpatial->fGetCylindricalDetectorShift();
        
        double      dShift  = (double)fShift;
        double      dRadius = (double)m_poSpatial->fGetCylindricalDetectorRadius();
        
        if( fShift != 0.0f )
        {
            // RB: The commented code below may need to be re-visited. For now I am not
            // applying any "rotation correction" to the vector fS, because I simply don't 
            // understand what that correction is for.
   
            //vCopyMat3D(&m_a3x3fDDRot[0][0], &a3x3fRotMatInv[0][0]);
            //vTranMat3D(a3x3fRotMatInv);
            //vMulMat3DVec3D(a3x3fRotMatInv, fS, fSRot);
            //// Use the line below to disable the rotation correction.  (just for testing).
            //// vCopyVec3D(fS,fSRot);
            //f0 = sqrt( fSRot[1] * fSRot[1] + fSRot[2] * fSRot[2] );
            
            double      dSy_squared_plus_Sz_squared = fS[1] * fS[1] + fS[2] * fS[2];
            if( 0.0 == dSy_squared_plus_Sz_squared )
                return 1; // This means that a particular scattered vector is parallel to the cylinder axis 
                          // and therefore cannot intersect the cylindrical surface 
                          // no matter how tall the cylinder is and no matter how long the scattered vector is.
                          // Therefore, we simply cannot find the scale coefficient for vector fS.  
            
            double      dShift_times_Sz = dShift * fS[2];
            double      dShift_times_Sz_squared = dShift_times_Sz * dShift_times_Sz;
            double      dShift_squared = dShift * dShift;
            double      dRc_squared = dRadius * dRadius;

            double      dSquareRootArg = dShift_times_Sz_squared + (dRc_squared - dShift_squared) * dSy_squared_plus_Sz_squared;

            if( dSquareRootArg < 0.0 )
                return 1; // This means that a particular scattered vector cannot intersect the cylindrical surface 
                          // no matter how tall the cylinder is and no matter how long the scattered vector is.
                          // Therefore, we simply cannot find the scale coefficient for vector fS.  

            double      dScaleCoeff = ( sqrt(dSquareRootArg) - dShift_times_Sz ) / dSy_squared_plus_Sz_squared;

            fRMFAC = (float)dScaleCoeff;
            
            // Get the local detector normal
            m_a3fLocalDN[0] = 0.0;
            m_a3fLocalDN[1] = fS[1];
            m_a3fLocalDN[2] = fS[2] + (float) dShift / fRMFAC;
        }
        else
        {
            /// RB This is the old Thad's code for handling the case, when the shift is zero.
            ////////////////////////////////////////////////////////////////
            vCopyMat3D(&m_a3x3fDDRot[0][0], &a3x3fRotMatInv[0][0]);
            vTranMat3D(a3x3fRotMatInv);
            /////////////////////////////////////////////////////////////////

            vMulMat3DVec3D(a3x3fRotMatInv, fS, fSRot);
        
            // Use the line below to disable the rotation correction.  (just for testing).
            // vCopyVec3D(fS,fSRot);
        
            f0 = sqrt( fSRot[1] * fSRot[1] + fSRot[2] * fSRot[2] );
        
            if( f0 == 0.0 )
                return 1;

		    // Get the local detector normal (it's always in the [1]x[2] plane.
		    m_a3fLocalDN[0] = 0.0;
		    m_a3fLocalDN[1] = fSRot[1];
		    m_a3fLocalDN[2] = fSRot[2];

            fRMFAC = fDetDistance / f0;
        }
    } 
    else
        return 1;

    vMulVec3DScalar(fS, fRMFAC, fTS);
    
    return 0; 
}


int Cdetector::nCalcDetCoords(const float fXR[3], 
                              const float fS0[3],
			                  float *pfXmm, 
                              float *pfYmm,
			                  float *pfPx0, 
                              float *pfPx1)
{
  // Calculate the detector coordinates of a reflection, given
  // dimensionless reciprocal lattice coordinates at diffracting
  // condition, detector info and fS0;

  //  This routine is not meant to be run in a loop as it is
  //  also calculates the prerequisites (including an inverted matrix).
  //
  //  Calculates spot coordinates on the detector and tests
  //  whether spot is within detector limits.
  //  Return 0 on success, non-zero on failure

  //  This subroutine created using code written by Albrecht
  //  Messerschmidt.    November 1986 see Proceedings Phase III pp. 57...

  //
  // fXR[3]           Dimensionless reciprocal lattice coordinate of this
  //                      reflection in the laboratory frame.
  // fInvDD[3][3]     Matrix used in calculation of detector coords
  //
  // *pfPx0, *pfPx1  Reflection detector pixel coordinates in FAST, SLOW directions
  // *pfXmm, *pfYmm  Reflection mm coordinates in lab X, Y directions
  //

  //  Local variables

  float     fXRdotDN = 0.0f;
  
  float     fS0dotDN = 0.0f;
  
  float     fXRminusS0[3] = {0.0f};
  
  float     fTS[3] = {0.0f};     // RB: this is the scattered vector S multiplied by the scalar T
  
  float     fVDCC[3] = {0.0f};
  
  int       nStat = 0;
  
  float     fInvDD[3][3];
  
  float     fDet = 0.0f;
  
  float     fZRelativeCoordinate = 0.0f;

  // Get detector DD matrix and DNvec

  vUpdateDDDN();
  
  fDet = fInvMat3D(&m_a3x3fDD[0][0], &fInvDD[0][0]);    // Invert DD matrix
  
  if( 0.0 == fDet )
      return (2);

  // Convert the coordinates of the vector XR in reflection
  // position in diffractometer coordinate system into
  // mm-detector-coordinates.

  // fVDCC - calculated detector coordinates on inclined
  //         detector.
  //
  //  Now we need  T S = (XR - S0) * DIST**2 / (XR.DN - S0.DN)   
  //                 -    --   --               -- --   -- --
  // See Proceedings Phase III, pp. 57-64

  // RB 2/3/2006 The equation above is actually needed to calculate the scale factor T,
  // to be applied to the scattered vector S. This is done in function nScaleSOnPlate()

  // Notice that for the cylindrical detector the equation for determining T is
  //
  //   T**2 * (Sy**2 + Sz**2) = Rc**2, 
  //   where Rc is the cylinder radius.
  // 
  // If the cylidrical detector is shifted by a "shift" along Z-axis, The equation becomes:
  //
  //   T**2 * (Sy**2 + Sz**2) + (2 * shift * T * Sz) + (shift**2 - Rc**2) = 0
  //

  fXRdotDN     = fDot3D(fXR, m_a3fDN);
  fS0dotDN     = fDot3D(fS0, m_a3fDN);

  // There cannot be an intersection if these have a positive dot product.
  // Otherwise, we have the beam hitting the image.
  if( m_poSpatial->eSpatialType() == eSpatial_flat_geometry && fXRdotDN < fS0dotDN )
      return (1);

  vSubVec3DVec3D(fXR, fS0, fXRminusS0);   // Also used by oblique incidence calc
  
  nStat = nScaleSOnPlate(fXRminusS0, fTS);

  if( 0 == nStat )
  {
      // The calculated detector coordinates VDCC are derived
      // by    VDCC = fInvDD * TS
      //       ----   =====     -

      vMulMat3DVec3D(fInvDD, fTS, fVDCC);

     *pfXmm = fVDCC[0];
     *pfYmm = fVDCC[1];
     fZRelativeCoordinate   = fVDCC[2];

     // Convert mm to pixels to see if reflection falls on the detector:

     nStat = m_poSpatial->nMMtoPixel(*pfXmm, *pfYmm, fZRelativeCoordinate, pfPx0, pfPx1);
  }

  return nStat;
}

int
Cdetector::nUpdateHeader(Cimage_header* poHeader, const Cstring& sPre)
{
  // Preliminary, does not update all possible detector items

  Cstring sPrefix;
  sPrefix = sPre;
  if ("" == sPrefix)
    {
      // If the input prefix is blank, then use the member variable as the
      // prefix.

      sPrefix = m_sPrefix;
    }
  return (m_poGoniometer->nUpdateHeader(poHeader, sPrefix));
}

int
Cdetector::nUpdateFromHeader(Cimage_header &roHeader, const Cstring& rsPre,
			     const bool bSpatial, const bool bNonunf)
{
  int     nStat, nTemp;
  Cstring sTemp;

  m_sPrefix = rsPre;

  if ("" == m_sPrefix)
    {
      // Null prefix string, see if the ms_sDetectorNames keyword exists
      // in the header.  If so, change the prefix to the first name
      // in the list of names.  This call is OK, since if there is no
      // ms_sDetectorNames keyword, then sPrefix is set to "" anyways.

      (void) roHeader.nGetValue(ms_sDetectorNames, 1, &m_sPrefix);
    }

  if (!roHeader.bIsAvailable())
    {
      cout << "Cdetector ERROR: image header is not valid."
	   << "  Cannot update detector!\n";;
      nStat = 1;
    }
  else if (0 == roHeader.nGetValue(m_sPrefix + ms_sDetectorKey, &sTemp))
    {
      cout << "Cdetector ERROR: cannot use database key yet!\n";
      nStat = 2;
    }
  else
    {
      // Try to get required information from the header

      nStat = 0;

      // Get detector goniometer from header
      // Detector position always come from the detector goniometer object

      if (NULL != m_poGoniometer)
	{
	  delete m_poGoniometer;
	}
      m_poGoniometer = new Cgoniometer(roHeader, m_sPrefix);
      m_bNewGonio    = TRUE;
      if (!m_poGoniometer->bIsAvailable()) nStat++;

    ///////////////////////////////////////////////////////////////////////////////////////
    if( bSpatial && NULL == m_poSpatial ) 
    {
        // Get spatial distortion info only if requested
        m_poSpatial    = new Cspatial(roHeader, m_sPrefix);

        if ( !m_poSpatial->bIsAvailable() )
            nStat++;
	}
    else if( NULL == m_poSpatial )
	{
        // Otherwise construct a minimal spatial distortion object
        // if one does not already exist
        m_poSpatial   = new Cspatial();
    }
    else
	{
	    // Try to get the direct beam position and update it.
	    // This does not update the pixel size.
        float         a2fBeamPx[2] = {0.0, 0.0};
        nTemp = roHeader.nGetValue(m_sPrefix + Cspatial::ms_sSpatialBeamPosn, 2, a2fBeamPx);
        m_poSpatial->nSetBeamPosition(a2fBeamPx[0], a2fBeamPx[1]);
	}
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    if( bNonunf && NULL == m_poNonunf )
	{
        // Get nonunf info only if requested
        m_poNonunf = new Cnonunf(roHeader, m_sPrefix);

        if ( !m_poNonunf->bIsAvailable() ) 
            nStat++;
	}
    else if( NULL == m_poNonunf )
	{
        // Otherwise construct a minimal nonunf object
        // if one does not already exist
        m_poNonunf = new Cnonunf();
	}

    // Get direction of fast pixel and slow pixel on the detector face

    //23-Oct-1996 Remove this, since the default of 1 0 0 0 1 0 is what we
    //            always want nowadays
    //      nTemp = roHeader.nGetValue(m_sPrefix + ms_sDetectorVectors, 6,
    //				&m_a2x3fVector[0][0]);

    nTemp = roHeader.nGetValue(m_sPrefix + ms_sDetectorDimensions, 2, m_a2nDim);
    
    if (0 != nTemp)
      {
	cout << "Cdetector WARNING no: " << m_sPrefix
             << ms_sDetectorDimensions << " keyword!\n";
	nStat++;
      }
    else if (m_poSpatial->bIsAvailable())
      {
	// Adjust the maximum allowed pixel in the spatial distortion object
	(void) m_poSpatial->nSetMaxPixel(m_a2nDim[0]-1, m_a2nDim[1]-1);
      }

      nTemp = roHeader.nGetValue(m_sPrefix + ms_sDetectorSize, 2, m_a2fSize);
      if (0 != nTemp)
	{
	  cout << "Cdetector WARNING: no " << m_sPrefix
               << ms_sDetectorSize << " keyword!\n";
	  nStat++;
	}

      // It is not required to have the following keyword, so
      // do not flag if not found.

      nTemp = roHeader.nGetValue(m_sPrefix + ms_sDetectorDescription,
				&m_sDescription);
    }

  if (0 != nStat)
    {
      // If there are errors, no sense having this around
      m_eThe_State = eDetector_unknown_state;
      return (1);
    }
  else
    {
      m_eThe_State = eDetector_notrefined_state;
      return (0);
    }
}

float
Cdetector::fGetPixelDistance(float fPxA0, float fPxA1,
			     float fPxB0, float fPxB1)
{
  // Return the distance between 2 pixels in mm.  Return -1.0 if error.
  // If fPxA0 < 0.0, then the first pixels is set to be near the
  // detector center, and the second pixel is a relative offset from that
  // position.

  float a2fAmm[2], a2fBmm[2];
  float fDist;
  float fZRelativeCoordinate = 0.0f;
  int   nStatA, nStatB;

  fDist = -1.0;
  if ( (NULL == m_poSpatial) || (!m_poSpatial->bIsAvailable()) )
    {
      return (fDist);
    }

  // These computations will assume that the Z coordinate can be ignored, since we are
  // only using differential computations.

  if (0.0 <= fPxA0)
    {
      nStatA   = m_poSpatial->nPixeltoMM(fPxA0, fPxA1, &a2fAmm[0], &a2fAmm[1],&fZRelativeCoordinate);
      nStatB   = m_poSpatial->nPixeltoMM(fPxB0, fPxB1, &a2fBmm[0], &a2fBmm[1],&fZRelativeCoordinate);
    }
  else
    {
      float fPxC0, fPxC1, fPxD0, fPxD1;
      nStatA = 1;
      fPxC0 = 0.4f * (float) m_a2nDim[0] + 1.0f;
      fPxC1 = 0.4f * (float) m_a2nDim[1];
      while ( (0 != nStatA) && (0.0 < fPxC0) )
	{
	  fPxC0  -= 1.0;
	  nStatA = m_poSpatial->nPixeltoMM(fPxC0, fPxC1, &a2fAmm[0], &a2fAmm[1],&fZRelativeCoordinate);
	}
      fPxD0 = fPxC0 + fPxB0;
      fPxD1 = fPxC1 + fPxB1;
      nStatB   = m_poSpatial->nPixeltoMM(fPxD0, fPxD1, &a2fBmm[0], &a2fBmm[1],&fZRelativeCoordinate);
    }
  if ( (0 == nStatA) && (0 == nStatB) )
    {
      float fTemp;
      fTemp = (a2fAmm[0] - a2fBmm[0]);
      fDist = fTemp * fTemp;
      fTemp = (a2fAmm[1] - a2fBmm[1]);
      fDist = sqrt((double)(fDist +  fTemp * fTemp));
    }
  return (fDist);
}


int Cdetector::nDiff(Cdetector& oDetectorAdd,Cdetector& oDetectorSubtract) {
    int nStat;

    if ((!m_poGoniometer) || (!oDetectorAdd.m_poGoniometer) || (!oDetectorSubtract.m_poGoniometer) ||
        (!m_poSpatial) || (!oDetectorAdd.m_poSpatial) || (!oDetectorSubtract.m_poSpatial))
        return 1;

    nStat = m_poSpatial->nDiff(*oDetectorAdd.m_poSpatial,*oDetectorSubtract.m_poSpatial);
    return nStat + m_poGoniometer->nDiff(*oDetectorAdd.m_poGoniometer,*oDetectorSubtract.m_poGoniometer);
};




int Cdetector::nCalcPixelShift(const float _fS0[3],
                               double a3fRotVec[3],
                               double fLambda1,
                               double fLambda2,
                               int nPix0,
                               int nPix1,
                               float* pfPixWave)
{

	int nStat;
    int nTotalStat;
	float  fPix0,fPix1;
	float  fVDC[3];     // Virtual Detector Coordinate
	float  fS_[3];      // Scattered beam wavevector
	double fS[3];       
	double fS0[3];
	double fX[3];
    double fXLambdaAdj[3];
    double fXRot[3];
	double a3fVDCC[3];
	double a3x3fDD[3][3];
	double a3x3fDDInv[3][3];
	double a3fVTS[3];

	vCopyVec3D(_fS0,fS0);
	fNormVec3D(fS0);
	
	// Determine the S vector at the current pixel position.
	nStat   = m_poSpatial->nPixeltoMM(nPix0, nPix1, &fVDC[0], &fVDC[1], &fVDC[2]);
	if (nStat) 
		return nStat;
	
	vMulMat3DVec3D(m_a3x3fDD, fVDC, fS_);  // Calc scattered beam wavevector
	vCopyVec3D(fS_,fS);
	fNormVec3D(fS);
	vAddVec3DVec3D(fS, fS0, fX);           // S + S0 = X

    nStat = 0;
    nTotalStat = 0;
    
    // This calculates the lambda modulated X vector
    vMulVec3DScalar(fX, fLambda2/fLambda1, fXLambdaAdj);
    
    if (pfPixWave) {
        if (fCalcRotationOffset(a3fRotVec,fS0,fXLambdaAdj,&fXRot[0]) <= -9999.0)
            nStat = 2;
        if (!nStat) {
            vSubVec3DVec3D(fXRot,fS0,fS);
            
            if (nScaleSOnPlate(fS,a3fVTS))
                nStat = 2;
            if (!nStat) {
                
                vCopyMat3D(&m_a3x3fDD[0][0],&a3x3fDD[0][0]);
                fInvMat3D(&a3x3fDD[0][0],&a3x3fDDInv[0][0]);
                
                vMulMat3DVec3D(a3x3fDDInv, a3fVTS, a3fVDCC);
                if (m_poSpatial->nMMtoPixel(a3fVDCC[0], a3fVDCC[1],a3fVDCC[2],&fPix0, &fPix1,false))
                    nStat = 2;
                
                pfPixWave[0] = fPix0;
                pfPixWave[1] = fPix1;
            }
        }
    };
    
    nTotalStat += nStat;
    
	return nTotalStat;		
};

const int nShiftArrayMaxPixelsPerEntry = 20;
const int nShiftArrayEntries = 4;
const int nProfileArrayMaxPixelsPerEntry = 20;
const int nProfileArrayDim = 60;

int Cdetector::nCalcPixelShiftArray(Csource& oSource,Crotation& oRotation,Cgoniometer& oCrysGonio) {

	float	a3fS0[3];
	float	afTemp[nShiftArrayEntries];
    float   a4fTemp[4];
    double  alfTemp[nShiftArrayEntries];
    double  a3fRotVec[3];
	double  fWavelength0;
	double	fWavelength1,fWavelength2;
    bool bDebug = true;

    

	int nPix0,nPix1,nPix0Shift,nPix1Shift,nShift0,nShift1;
	int nPos0,nPos1;
	int nx,ny;
    int nCount;
    int nAxis;
    int nError;


    
    nAxis = oCrysGonio.nGetNum(oRotation.sGetName());
    if (0 > nAxis)
        nAxis = 0;
    if (oCrysGonio.nGetRotVector(nAxis,&a4fTemp[0])) 
        return 1;
    vCopyVec3D(a4fTemp,a3fRotVec);



    oSource.vCalcGetS0(a4fTemp);            
    vCopyVec3D(a4fTemp,a3fS0);    
	(void) nCalcGetDDDN(&m_a3x3fDD[0][0], m_a3fDN);    
	fWavelength0 = oSource.m_poWavelength->fGetWavelength(0);
	fWavelength1 = oSource.m_poWavelength->fGetWavelength(1);
	fWavelength2 = oSource.m_poWavelength->fGetWavelength(2);

	// We might not have Ka1/Ka2 information.
	if ((fWavelength1==0.0) || (fWavelength2==0.0)) {
		fWavelength1 = fWavelength0;
		fWavelength2 = fWavelength0;
	};
	if (fWavelength0 == 0.0)
		return 1;

    nx = m_a2nDim[0]/nShiftArrayMaxPixelsPerEntry;
    ny = m_a2nDim[1]/nShiftArrayMaxPixelsPerEntry;

    if ((m_pfPixelShiftArray) && ((m_a2nPixelShiftArrayDiv[0] != nx) || (m_a2nPixelShiftArrayDiv[1] != ny))) {
        delete[] m_pfPixelShiftArray;
        m_pfPixelShiftArray = NULL;
    };
    m_a2nPixelShiftArrayDiv[0] = nx;
    m_a2nPixelShiftArrayDiv[1] = ny;

    if (!m_pfPixelShiftArray) {                
        m_pfPixelShiftArray = new double[nShiftArrayEntries*m_a2nPixelShiftArrayDiv[0]*m_a2nPixelShiftArrayDiv[1]];
        memset(m_pfPixelShiftArray,0,sizeof(double)*nShiftArrayEntries*m_a2nPixelShiftArrayDiv[0]*m_a2nPixelShiftArrayDiv[1]);
    };

	nCount = 0;
	for (nPos0 = 0; nPos0 < m_a2nPixelShiftArrayDiv[0];nPos0++) {
		nPix0 = (int) ((0.5*m_a2nDim[0])/m_a2nPixelShiftArrayDiv[0] + (1.0*nPos0*m_a2nDim[0])/m_a2nPixelShiftArrayDiv[0]);
		for (nPos1 = 0; nPos1 < m_a2nPixelShiftArrayDiv[1];nPos1++) {
            nPix1 = (int) ((0.5*m_a2nDim[1])/m_a2nPixelShiftArrayDiv[1] + (1.0*nPos1*m_a2nDim[1])/m_a2nPixelShiftArrayDiv[1]);

            // This loop was added to perhaps help get over the 'bad pixel' problems.

            nx = nPos0 + nPos1*m_a2nPixelShiftArrayDiv[0];

            nError = 1;
            for (nShift0 = -1; (nShift0 <=1) && (nError); nShift0++) {
                for (nShift1 = -1; (nShift1 <=1) && (nError); nShift1++) {
                    nPix0Shift = nPix0 + nShift0;
                    nPix1Shift = nPix1 + nShift1;
                                        
                    ny = nCalcPixelShift(a3fS0,a3fRotVec,fWavelength0,fWavelength1,nPix0,nPix1,&afTemp[0]);
                    ny |= nCalcPixelShift(a3fS0,a3fRotVec,fWavelength0,fWavelength2,nPix0,nPix1,&afTemp[2]);

                    alfTemp[0] = afTemp[0] - nPix0;
                    alfTemp[1] = afTemp[1] - nPix1;
                    alfTemp[2] = afTemp[2] - nPix0;
                    alfTemp[3] = afTemp[3] - nPix1;

                    if ((!ny) && (!nCompVecND(nShiftArrayEntries,&alfTemp[0],&m_pfPixelShiftArray[nShiftArrayEntries*nx]))) {
                        // We are almost certainly calculating the identical values.  Break out.
                        goto bypass_calculation;
                    };
                    
                    if (!ny) {
                        nError = 0;
                        m_pfPixelShiftArray[nShiftArrayEntries*nx + 0] = alfTemp[0];
                        m_pfPixelShiftArray[nShiftArrayEntries*nx + 1] = alfTemp[1];
                        m_pfPixelShiftArray[nShiftArrayEntries*nx + 2] = alfTemp[2];
                        m_pfPixelShiftArray[nShiftArrayEntries*nx + 3] = alfTemp[3];
                        nCount++;
                    };
                };
            };

    	};
	};
	if (nCount == 0) {
		return 1;
	};

bypass_calculation:
    
	return 0;
};

int Cdetector::nCalcPixelShift(double fPix0,double fPix1,double* pfWaveShift0,double* pfWaveShift1) {
	int a2nPos0[2],a2nPos1[2];
	int nx,ny,nz;
    int nPix0,nPix1;
    double afTemp[nShiftArrayEntries];

    if (!m_pfPixelShiftArray)
        return 1;

    nPix0 = (int) fPix0;
    nPix1 = (int) fPix1;

	a2nPos0[0] = max(0,min(m_a2nPixelShiftArrayDiv[0]-1,(nPix0*m_a2nPixelShiftArrayDiv[0])/m_a2nDim[0]));
	a2nPos1[0] = max(0,min(m_a2nPixelShiftArrayDiv[1]-1,(nPix1*m_a2nPixelShiftArrayDiv[1])/m_a2nDim[1]));
    if (a2nPos0[0] + 1 < m_a2nPixelShiftArrayDiv[0])
        a2nPos0[1] = a2nPos0[0] + 1;
    else
        a2nPos0[1] = a2nPos0[0] - 1;
    if (a2nPos1[0] + 1 < m_a2nPixelShiftArrayDiv[1])
        a2nPos1[1] = a2nPos1[0] + 1;
    else
        a2nPos1[1] = a2nPos1[0] - 1;


    nGridInterpClear(nShiftArrayEntries);
    for (nx = 0; nx < 2; nx++) {
        for (ny = 0; ny < 2; ny++) {
            nPix0 = (int) ((0.5*m_a2nDim[0])/m_a2nPixelShiftArrayDiv[0] + (1.0*a2nPos0[nx]*m_a2nDim[0])/m_a2nPixelShiftArrayDiv[0]);
            nPix1 = (int) ((0.5*m_a2nDim[1])/m_a2nPixelShiftArrayDiv[1] + (1.0*a2nPos1[ny]*m_a2nDim[1])/m_a2nPixelShiftArrayDiv[1]);
            nz = a2nPos0[nx] + a2nPos1[ny]*m_a2nPixelShiftArrayDiv[0];
            if (m_pfPixelShiftArray[nShiftArrayEntries*nz]!=0.0)
                nGridInterpAdd(nPix0,nPix1,nShiftArrayEntries,&m_pfPixelShiftArray[nShiftArrayEntries*nz]);
        };
    };
    if (nGridInterpSolve(fPix0,fPix1,nShiftArrayEntries,&afTemp[0]))
        return 1;
	pfWaveShift0[0] = afTemp[0];
	pfWaveShift0[1] = afTemp[1];
	pfWaveShift1[0] = afTemp[2];
	pfWaveShift1[1] = afTemp[3];

    /* But we are not done.  A small error in the values of our shifts is probably noticible.
       We will try scaling both quantities according to the formula:

        (Length(Ka1) + delta)/(Length(Ka2) + delta) = 0.5.

        So, delta = Length(Ka2) - Length(Ka1)/0.5;

    
    double fLengthKa1,fLengthKa2,fDelta,fMult;

    fLengthKa1 = sqrt(pfWaveShift0[0]*pfWaveShift0[0] + pfWaveShift0[1]*pfWaveShift0[1]);
    fLengthKa2 = sqrt(pfWaveShift1[0]*pfWaveShift1[0] + pfWaveShift1[1]*pfWaveShift1[1]);
    fDelta =  fLengthKa2 - fLengthKa1*2.0;
    fMult = 1.0 + fDelta/fLengthKa1;
    pfWaveShift0[0]*=fMult;
    pfWaveShift0[1]*=fMult;
    fMult = 1.0 + fDelta/fLengthKa2;
    pfWaveShift1[0]*=fMult;
    pfWaveShift1[1]*=fMult;
    */
	return 0;
};


/*
#include "CSphericalHarmonic.h"


int Cdetector::nTest(Csource& oSource,Crotation& oRotation) {
	int nPix0,nPix1;
	int nPos0,nPos1;
	double a3fS0[3];
	double fS[3];
	double a3x3fBasis[3][3];
	int nStat;
	double a2fWaveShift0[2];
	double a2fWaveShift1[2];

	float a4fTemp[4];
	float fVDC[3];
	float _fS[3];
	double fPolar,fAzi;
	double f0;
	int nPass;
	int nx;

	CSphericalHarmonic oHarmonics;
	itr<double> afValues;
	itr<double> afValuesSigma;
	itr<double> afPolar;
	itr<double> afAzi;
	itr<int> anGroup;
	itr<double> afScanRot;
	itr<double> afScaleFactorsOut;


	oSource.vCalcGetS0(a4fTemp);            
    vCopyVec3D(a4fTemp,a3fS0);    
	(void) nCalcGetDDDN(&m_a3x3fDD[0][0], m_a3fDN);    
	vBuildBasis3D(a3fS0,a3x3fBasis);
	for (nx=0;nx<3;nx++)
		std::swap(a3x3fBasis[0][nx],a3x3fBasis[2][nx]);
	for (nx=0;nx<3;nx++)
		a3x3fBasis[0][nx] *= -1.0;

	double fMaxDiff;
	double fAvgDiff;
	int    nCountDiff;
	double fAvgRef;
    double fValue,fNewValue;
    int nMaxPos0,nMaxPos1;
    Cimage* poImage = NULL;

	-afValues;
	-afValuesSigma;
	-afPolar;
	-afAzi;
	-anGroup;
	-afScanRot;
	fMaxDiff = 0.0;
	fAvgDiff = 0.0;
	nCountDiff = 0;
	for (nPass = 0; nPass < 2; nPass++) {
		for (nPos0 = 0; nPos0 < m_a2nPixelShiftArrayDiv[0];nPos0++) {
			nPix0 = (int) ((0.5*m_a2nDim[0])/m_a2nPixelShiftArrayDiv[0] + (1.0*nPos0*m_a2nDim[0])/m_a2nPixelShiftArrayDiv[0]);
			for (nPos1 = 0; nPos1 < m_a2nPixelShiftArrayDiv[1];nPos1++) {
				nPix1 = (int) ((0.5*m_a2nDim[1])/m_a2nPixelShiftArrayDiv[1] + (1.0*nPos1*m_a2nDim[1])/m_a2nPixelShiftArrayDiv[1]);
				
				// Determine the S vector at the current pixel position.
				nStat   = m_poSpatial->nPixeltoMM(nPix0, nPix1, &fVDC[0], &fVDC[1], &fVDC[2]);
				if (nStat) 
					return nStat;
				
				vMulMat3DVec3D(m_a3x3fDD, fVDC, _fS);  // Calc scattered beam wavevector
				vCopyVec3D(_fS,fS);
				fNormVec3D(fS);
                fValue = a2fWaveShift0[0] + 20;
				nCalcPixelShift(nPix0,nPix1,&a2fWaveShift0[0],&a2fWaveShift1[0]);
				if (nPass == 0) {
					afValues + fValue;
					afValuesSigma + 1.0;
					CSphericalHarmonic::vVec3DToSpherical(&fS[0],fPolar,fAzi,&a3x3fBasis[0][0]);
					afPolar + fPolar;
					afAzi + fAzi;
					anGroup + 0;
					afScanRot + 0.0;
				} else {
                    fNewValue = fAvgRef*afScaleFactorsOut[nCountDiff];
					f0 = (fNewValue - fValue);
					fAvgDiff += ABS(f0);
					fMaxDiff = max(fMaxDiff,ABS(f0));
					nCountDiff++;
                    poImage->nSetPixel(nPos0,nPos1,(unsigned short int) (100.0*(fValue)));
                    poImage->nSetPixel(nPos0,nPos1 + nMaxPos1,(unsigned short int) (100.0*(fNewValue)));

				};
				
			};
            nMaxPos1 = nPos1;
		};
        nMaxPos0 = nPos0;

        if (!poImage)
            poImage = new Cimage(nMaxPos0,nMaxPos1*2,eImage_uI2);


		oHarmonics.vSetPrint(1);
		oHarmonics.vSetScaleInfo(-10000000.0,-1000000.0,1,0.1);
		oHarmonics.vSetHarmonics(4,0,2);
		oHarmonics.vSetCentroSymmetric(0);

		oHarmonics.nScale(afValues,afValuesSigma,afPolar,afAzi,anGroup,afScanRot,afScaleFactorsOut);

        fAvgRef = 0.0;
        for (nx= 0; nx< afValues.size();nx++) {
			fAvgRef += afValues[nx];
        };
		fAvgRef/=nx;

	};
    poImage->nWrite("comp.img");
	
	return 0;
};
*/



int Cdetector::nCalcPixelShiftBin(double a2fKa1[2],double a2fKa2[2],double a3x2x2fShiftBin[3][2][2],int a3x2x2nShiftBin[3][2][2],int a3x2nShiftBin[3][2]) {
    int nKa12;
    int nx;


    for (nKa12 = 0; nKa12 < 3; nKa12++) {
        double pfShift[2];
        if (nKa12==0) {
            pfShift[0] = a2fKa1[0];
            pfShift[1] = a2fKa1[1];
        } else if (nKa12==1) {
            pfShift[0] = a2fKa2[0];
            pfShift[1] = a2fKa2[1];
        } else {
            pfShift[0] = a2fKa2[0] - a2fKa1[0];
            pfShift[1] = a2fKa2[1] - a2fKa1[1];
        };
        for (nx=0; nx < 2; nx++) {
            a3x2x2nShiftBin[nKa12][nx][0] = (int) floor(pfShift[nx]);
            a3x2x2nShiftBin[nKa12][nx][1] = a3x2x2nShiftBin[nKa12][nx][0] + 1;
            a3x2x2fShiftBin[nKa12][nx][0] = ABS(pfShift[nx] - a3x2x2nShiftBin[nKa12][nx][1]);
            a3x2x2fShiftBin[nKa12][nx][1] = ABS(pfShift[nx] - a3x2x2nShiftBin[nKa12][nx][0]);
            a3x2nShiftBin[nKa12][nx] = ((a3x2x2fShiftBin[nKa12][nx][0]>a3x2x2fShiftBin[nKa12][nx][1])?(a3x2x2nShiftBin[nKa12][nx][0]):(a3x2x2nShiftBin[nKa12][nx][1]));
        };
    };
    return 0;
};
 

int Cdetector::nCheckReflnCharacteristics(double fCentPix0,
                                          double fCentPix1,
                                          double* a2x2fA,
                                          double* a2fb,
                                          double fc)
{
    double afWaveShift[2][2];
    double fOffsetPix0,fOffsetPix1;
    int nx;
    int nSlicesWorked;


    nSlicesWorked = 0;

    // Get the shifts for this spot from the detector object.
    if (nCalcPixelShift(fCentPix0,fCentPix1,&afWaveShift[0][0],&afWaveShift[1][0]))
        return 1;

    if (fEvalEllipse(fCentPix0,fCentPix1,a2x2fA,a2fb,&fc) > 1.0) {
        printf("ERROR:  Problem in C3Ddata::nCheckReflnCharacteristics().\n");
        return 0;
    };
    // We must have the two wave shift vectors inside the ellipsoid to pass this test.
    for (nx = 0; nx < 2 ;nx++) {
        fOffsetPix0 = fCentPix0 + afWaveShift[nx][0] ;
        fOffsetPix1 = fCentPix1 + afWaveShift[nx][1] ;
        if (fEvalEllipse(fOffsetPix0,fOffsetPix1,a2x2fA,a2fb,&fc))
            break;
    };
    if (nx ==2)
        return 0;
    return 1;

};





