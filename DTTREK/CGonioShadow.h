//
// Copyright (c) 2004 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of Argonne National Laboratory.
//
// CGonioShadow.h           Initial author: RB Aug 2004
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


#ifndef GONIO_SHADOW_INCLUDE
#define GONIO_SHADOW_INCLUDE

#include "Dtrek.h"
#include "dtrekdefs.h"
#include "Cstring.h"

// Defines
typedef struct _tagXY
{
    double   dX;
    double   dY;
}XY;

class Cimage;
////////////////////////////////////////////////////////////////////////////////////////
// This class is a helper class for CGonioShadow.
class DTREK_EXPORT CGonioShadowGeometry
{
public:
    CGonioShadowGeometry(std::vector<double>* pdaShadowInfo, double daPixelSize[2]);
    virtual ~CGonioShadowGeometry(){}
    virtual  bool bRotate(double dRotationDeg, double dImageWidth)=0;
    virtual  bool bIsAvailable()=0;
    virtual  bool bErasePixels(Cimage* poImage)=0;
    
protected:
    double  m_dOmegaOffset;
    double  m_dPixelsPerDegree;
    double  m_dPixelScaleFactor_0;
    double  m_dPixelScaleFactor_1;
};

#define DT_PGSG_MIN_NUMBER_OF_VERTICES      3
class DTREK_EXPORT CPolygonGonioShadowGeometry : public CGonioShadowGeometry
{
public:
    CPolygonGonioShadowGeometry(std::vector<double>* pdaShadowInfo, double daPixelSize[2]);
    ~CPolygonGonioShadowGeometry();
    virtual  bool bRotate(double dRotationDeg, double dImageWidth);
    virtual  bool bIsAvailable();
    virtual  bool bErasePixels(Cimage* poImage);
    
private:
    // coordinates as defined by the header file
    double*     m_pdVertice_X_Coords;
    double*     m_pdVertice_Y_Coords;
    
    // coordinates after "rotating" the coordinates, defined in the header file, according to the rotation offset
    double*     m_pdRotated_Vertice_X_Coords;
    
    int         m_nNumberOfVertices;
};

class DTREK_EXPORT CCircularGonioShadowGeometry : public CGonioShadowGeometry
{
public:
    CCircularGonioShadowGeometry(std::vector<double>* pdaShadowInfo, double daPixelSize[2]);
    virtual  bool bRotate(double dRotationDeg, double dImageWidth);
    virtual  bool bIsAvailable();
    virtual  bool bErasePixels(Cimage* poImage);

private:
    // coordinates as defined by the header file
    double      m_dCenter_X_Coord;   
    double      m_dCenter_Y_Coord;
    
    // X coordinate after "rotating" the coordinates, defined in the header file, according to the rotation offset
    double      m_dRotated_Center_X_Coord;
    
    double      m_dInnerRadius;
    double      m_dOuterRadius;
};
////////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CGonioShadow
{
friend class CGonioShadowGeometry;
friend class CPolygonGonioShadowGeometry;
friend class CCircularGonioShadowGeometry;
public:
    CGonioShadow(Cstring& strShadowInfo, double daPixelSize[2]);
    virtual ~CGonioShadow();
    
    bool bIsAvailable();

#define     DT_GH_HARDWARE_LIMITS_OMEGA                 0x0001
#define     DT_GH_HARDWARE_LIMITS_CHI_KAPPA             0x0002
#define     DT_GH_HARDWARE_LIMITS_PHI                   0x0004
#define     DT_GH_HARDWARE_LIMITS_DET_DISTANCE          0x0008
#define     DT_GH_HARDWARE_LIMITS_DET_SWING             0x0010
    bool bIsHardwareLimitsApply(double dParameterValue, DTREK_WORD wParameterType);

    bool bRotate(double dRotationDeg, double dImageWidth){return m_poShadowGeometry->bRotate(dRotationDeg, dImageWidth);}
    bool bErasePixels(Cimage* poImage){return m_poShadowGeometry->bErasePixels(poImage);}

private:
    bool bReadData(Cstring strShadowInfo, double daPixelSize[2]);
    
    std::vector<double>    m_daShadowInfo;  // array with shadow information, read from an input text string

    enum enShadowInfoPosition
    {
        enShadow_Type=0,
        
        enShadow_omega_min,             // Applies only if image Omega is between dOmegaMin and dOmegaMax
        enShadow_omega_max,

        enShadow_chi_kappa_min,         // Applies only if image ChiKappa is between dChiKappaMin and dChiKappaMax
        enShadow_chi_kappa_max,
        
        enShadow_phi_min,               // Applies only if image Phi is between dPhiMin and dPhiMax
        enShadow_phi_max,
        
        enShadow_omega_offset,          // The shadow was defined at this Omega value
        enShadow_pixels_per_degree,     // Shift this object by this pixels per Omega degrees.  Use 0 for an omega-independent (i.e. non-moving) shadow.

        enShadow_det_dist_min,          // Applies only if image detector distance is between dDetDistMin and dDetDistMax
        enShadow_det_dist_max,

        enShadow_det_two_theta_min,     // Applies only if image detector two-theta swing angle is between dDetTwoThetaMin and dTwoThetaMax
        enShadow_det_two_theta_max,

        enShadow_det_pix_size_0,        // The shadow was defined at these pixel sizes
        enShadow_det_pix_size_1,

        enShadow_coordinates_start      // Start of the shadow geometrical description (polygon vertice coordinates or circular area center and two radii)
    };

    enum enShadowType
    {
         enShadowTypePolygon=0,
         enShadowTypeCircle 
    };
    
    CGonioShadowGeometry*        m_poShadowGeometry;
};
#endif  // GONIO_SHADOW_INCLUDE





