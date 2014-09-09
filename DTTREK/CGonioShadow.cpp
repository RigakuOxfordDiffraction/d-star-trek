//
// Copyright (c) 2004 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of Argonne National Laboratory.
//
// CGonioShadow.cpp           Initial author: RB Aug 2004
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

#include "CGonioShadow.h"
#include "Cimage.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


static int nPointInPolygon(int nPoints, double *xp, double *yp, double x, double y);
static void vCheckTransformAngleValue(double& dParameterValue);

CGonioShadow::CGonioShadow(Cstring& strShadowInfo, double daPixelSize[2])
{
    m_poShadowGeometry = NULL;
    
    bReadData(strShadowInfo, daPixelSize);
}
///////////////////////////////////////////////////////////////////////////////////
CGonioShadow::~CGonioShadow()
{
    if( m_poShadowGeometry )
    {
        delete m_poShadowGeometry;
        m_poShadowGeometry = NULL;
    }
}
///////////////////////////////////////////////////////////////////////////////////
bool CGonioShadow::bIsAvailable()
{
    if( !m_poShadowGeometry )
        return false;
    
    return m_poShadowGeometry->bIsAvailable();
}
///////////////////////////////////////////////////////////////////////////////////
bool CGonioShadow::bReadData(Cstring strShadowInfo, double daPixelSize[2])
{
    strShadowInfo.nListToVector(m_daShadowInfo, ' ');

    if( 0 == m_daShadowInfo.size() )
        return false;

    if( enShadowTypePolygon == (int)m_daShadowInfo[enShadow_Type] )
        m_poShadowGeometry = new CPolygonGonioShadowGeometry(&m_daShadowInfo, daPixelSize);
    else if( enShadowTypeCircle == (int)m_daShadowInfo[enShadow_Type] )
        m_poShadowGeometry = new CCircularGonioShadowGeometry(&m_daShadowInfo, daPixelSize);

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
bool CGonioShadow::bIsHardwareLimitsApply(double dParameterValue, DTREK_WORD wParameterType)
{
    double      dTestMin = 0.0;
    double      dTestMax = 0.0;

    switch( wParameterType )
    {
    case DT_GH_HARDWARE_LIMITS_OMEGA:
        dTestMin = m_daShadowInfo[enShadow_omega_min];
        dTestMax = m_daShadowInfo[enShadow_omega_max];
        vCheckTransformAngleValue(dParameterValue);
        break;
    case  DT_GH_HARDWARE_LIMITS_CHI_KAPPA:
        dTestMin = m_daShadowInfo[enShadow_chi_kappa_min];
        dTestMax = m_daShadowInfo[enShadow_chi_kappa_max];
        vCheckTransformAngleValue(dParameterValue);
        break;
    case  DT_GH_HARDWARE_LIMITS_PHI:
        dTestMin = m_daShadowInfo[enShadow_phi_min];
        dTestMax = m_daShadowInfo[enShadow_phi_max];
        vCheckTransformAngleValue(dParameterValue);
        break;
    case  DT_GH_HARDWARE_LIMITS_DET_DISTANCE:
        dTestMin = m_daShadowInfo[enShadow_det_dist_min];
        dTestMax = m_daShadowInfo[enShadow_det_dist_max];
        break;
    case  DT_GH_HARDWARE_LIMITS_DET_SWING:
        dTestMin = m_daShadowInfo[enShadow_det_two_theta_min];
        dTestMax = m_daShadowInfo[enShadow_det_two_theta_max];
        vCheckTransformAngleValue(dParameterValue);
        break;
    default:
        return false; // unknown parameter type
    }
    
    if( dParameterValue < dTestMin || dParameterValue > dTestMax )
        return false;

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////
CGonioShadowGeometry::CGonioShadowGeometry(std::vector<double>* pdaShadowInfo, double daImagePixelSize[2])
{
    m_dOmegaOffset = (*pdaShadowInfo)[CGonioShadow::enShadow_omega_offset];

    m_dPixelsPerDegree = (*pdaShadowInfo)[CGonioShadow::enShadow_pixels_per_degree];

    m_dPixelScaleFactor_0 = 1.;
    m_dPixelScaleFactor_1 = 1.;
        
    if( 0.0 !=  daImagePixelSize[0] )
        m_dPixelScaleFactor_0 = (*pdaShadowInfo)[CGonioShadow::enShadow_det_pix_size_0] / daImagePixelSize[0];
    if( 0.0 !=  daImagePixelSize[1] )
        m_dPixelScaleFactor_1 = (*pdaShadowInfo)[CGonioShadow::enShadow_det_pix_size_1] / daImagePixelSize[1];
}
/////////////////////////////////////////////////////////////////////////////////////////
CPolygonGonioShadowGeometry::CPolygonGonioShadowGeometry(std::vector<double>* pdaShadowInfo, double daImagePixelSize[2]) :
CGonioShadowGeometry(pdaShadowInfo, daImagePixelSize),
m_pdVertice_Y_Coords(NULL),
m_pdVertice_X_Coords(NULL),
m_pdRotated_Vertice_X_Coords(NULL),
m_nNumberOfVertices(0)
{
    if( pdaShadowInfo )
    {
        int     nCoordinatesVectorLength = pdaShadowInfo->size() - CGonioShadow::enShadow_coordinates_start;
    
        if( nCoordinatesVectorLength/2*2==nCoordinatesVectorLength &&
            nCoordinatesVectorLength >= DT_PGSG_MIN_NUMBER_OF_VERTICES * 2 ) // the polygon should be at least a triangle, i.e. 3 pairs of coords!    
        {
            m_nNumberOfVertices = nCoordinatesVectorLength/2;
            m_pdVertice_Y_Coords = new double[m_nNumberOfVertices];
            m_pdVertice_X_Coords = new double[m_nNumberOfVertices];
            
            m_pdRotated_Vertice_X_Coords = new double[m_nNumberOfVertices];
        
            for(int ii=0; ii < m_nNumberOfVertices; ii++)
            {
                m_pdVertice_X_Coords[ii] = (*pdaShadowInfo)[CGonioShadow::enShadow_coordinates_start + 2 * ii] * m_dPixelScaleFactor_0;
                m_pdVertice_Y_Coords[ii] = (*pdaShadowInfo)[CGonioShadow::enShadow_coordinates_start + 2 * ii + 1] * m_dPixelScaleFactor_1;
            }
        }
    }
}

CPolygonGonioShadowGeometry::~CPolygonGonioShadowGeometry()
{
    if( m_pdVertice_Y_Coords )
    {
        delete [] m_pdVertice_Y_Coords;
        m_pdVertice_Y_Coords = NULL;
    }
    
    if( m_pdVertice_X_Coords )
    {
        delete [] m_pdVertice_X_Coords;
        m_pdVertice_X_Coords = NULL;
    }
    
    if( m_pdRotated_Vertice_X_Coords )
    {
        delete [] m_pdRotated_Vertice_X_Coords;
        m_pdRotated_Vertice_X_Coords = NULL;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
bool CPolygonGonioShadowGeometry::bIsAvailable()
{
    if( !m_pdVertice_X_Coords || !m_pdVertice_Y_Coords ||  !m_pdRotated_Vertice_X_Coords )
        return false;
    
    if( m_nNumberOfVertices < DT_PGSG_MIN_NUMBER_OF_VERTICES )
        return false;
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
bool CPolygonGonioShadowGeometry::bRotate(double dRotationOffset, double dImageWidth)
{
    double      dPx0Offset = m_dOmegaOffset - dRotationOffset;

    dPx0Offset *= m_dPixelsPerDegree * m_dPixelScaleFactor_0;
    
    double      dPx0Min =  999999.0;
    double      dPx0Max = -999999.0;

    double      dPx0 = 0.0;

    for(int nPoint = 0; nPoint < m_nNumberOfVertices; nPoint++)
    {
	      dPx0 = m_pdVertice_X_Coords[nPoint] + dPx0Offset;
	      
          if( dPx0Min > dPx0 ) 
              dPx0Min = dPx0;
	      
          if( dPx0Max < dPx0 ) 
              dPx0Max = dPx0;

	      m_pdRotated_Vertice_X_Coords[nPoint] = dPx0;
    }
    //	  //cout << "\n\n";
    if( 20.0 < m_dPixelsPerDegree && 23.0 > m_dPixelsPerDegree )
    {
        // It is probably a shadow on a RAPID, so make sure it loops
        // around to other edge of image if required

        // RAPID pixels usually go from 0 to 4599 along Px0. 
        // A full 360 degrees will be 7980 pixels.
        // That is 22.16667 pixels per degree
        // So an object must fit between 0 and 4599 pixels

        //  There are 4 possibilities:
        //  a  Both entirely outside the image
        //  b. Polygon min/max are both entirely within the image
        //  c. Polygon min is in the image, but quad max is out of the image (on left)
        //  d. Polygon max is in the image, but quad min is out of the image (on right)

        // There is a periodicity of 7980 pixels (per 360 degrees) for
        // RAPID images

	      double        dMaxDim0    = 4600.0 * m_dPixelScaleFactor_0;
	      double        dAdd = 0.0;
	      double        dSub = 0.0;

	      // Shadows are defined for the full RAPID horizontal 
	      // readout of 4600 pixels.  If this is a narrow readout of 
	      // 2800 pixels then shift the shadows horizontal by 1800 pixels as the first 
	      // thing.

	      if( 2800 * m_dPixelScaleFactor_0 == dImageWidth )
          {
		       dAdd = -1800.0 * m_dPixelScaleFactor_0;
		  }

	    
          while( 0.0 > dPx0Min && 0.0 > dPx0Max )
          {
                // If both outside the image on the low side, try to see if
                //    inside the image due to periodicity
                dAdd    += 7980.0 * m_dPixelScaleFactor_0;
                dPx0Min += 7980.0 * m_dPixelScaleFactor_0;
                dPx0Max += 7980.0 * m_dPixelScaleFactor_0;
          }
		  
          while( dMaxDim0 <= dPx0Min && dMaxDim0 <= dPx0Max )
          {
		  // If both outside the image on the high side, try to see if
		  //    inside the image due to periodicity
            dSub    -= 7980.0 * m_dPixelScaleFactor_0;
            dPx0Min -= 7980.0 * m_dPixelScaleFactor_0;
            dPx0Max -= 7980.0 * m_dPixelScaleFactor_0;
          }
	    
          dAdd = dAdd + dSub;
        //cout << "dAdd is " << dAdd << endl;

        // We do not have to worry about a shadow or points being off 
          // the image.
        // The ::bErasePixels*() routines will make sure we do not 
          // go out-of-bounds

        for (int nPoint = 0; nPoint < m_nNumberOfVertices; nPoint++)
		{
		  // Shift first coordinate of each point by dAdd
		  m_pdRotated_Vertice_X_Coords[nPoint] += dAdd;
		}
    }
    
    
    return true;
}	  
///////////////////////////////////////////////////////////////////////////////////////////////
bool CPolygonGonioShadowGeometry::bErasePixels(Cimage* poImage)
{
    if( !poImage )
        return false;

    // Find the minimum and maximum pixel coords, so we could define a rectangle that holds the polygon
    double      dMinX = m_pdRotated_Vertice_X_Coords[0];
    double      dMaxX = m_pdRotated_Vertice_X_Coords[0];

    double      dMinY = m_pdVertice_Y_Coords[0];
    double      dMaxY = m_pdVertice_Y_Coords[0];

    for(int ii = 1; ii < m_nNumberOfVertices; ii++)
    {
        if (dMinX > m_pdRotated_Vertice_X_Coords[ii])
            dMinX = m_pdRotated_Vertice_X_Coords[ii];
        
        if (dMaxX < m_pdRotated_Vertice_X_Coords[ii])
            dMaxX = m_pdRotated_Vertice_X_Coords[ii];
        
        if (dMinY > m_pdVertice_Y_Coords[ii])
            dMinY = m_pdVertice_Y_Coords[ii];	  
        
        if (dMaxY < m_pdVertice_Y_Coords[ii])
            dMaxY = m_pdVertice_Y_Coords[ii];	  
    }

    int       nMin0 = (int) dMinX;
    int       nMax0 = (int) dMaxX;
    int       nMin1 = (int) dMinY;
    int       nMax1 = (int) dMaxY;

    if (0 > nMin0)
        nMin0 = 0;
    if (0 > nMax0)
        nMax0 = 0;
    if (0 > nMin1)
        nMin1 = 0;
    if (0 > nMax1)
        nMax1 = 0;

    if (nMin0 >= poImage->nGetDimension(0))
        nMin0 = poImage->nGetDimension(0) - 1;
    if (nMax0 >= poImage->nGetDimension(0))
        nMax0 = poImage->nGetDimension(0) - 1;
    if (nMin1 >= poImage->nGetDimension(1))
        nMin1 = poImage->nGetDimension(1) - 1;
    if (nMax1 >= poImage->nGetDimension(1))
        nMax1 = poImage->nGetDimension(1) - 1;
    /////////////////////////////////////////////////////////////////////////
  
    double      dX = 0.0;
    double      dY = 0.0;

    float       fErasedValue = 0.0;

    if( poImage->nBytesPerPixel() != sizeof(unsigned short int) )
    {
        cout << "WARNING! Erasing a polygon in non-two-byte pixels!\n";
    }

    // Loop over pixels in the image in the rectangle that holds the polygon
    for(int ny = nMin1; ny <= nMax1; ny++)
    {
        dY = (double)ny;
        for (int nx = nMin0; nx <= nMax0; nx++)
        {
            dX = (double)nx;
            if(0 != nPointInPolygon(m_nNumberOfVertices, m_pdRotated_Vertice_X_Coords, m_pdVertice_Y_Coords, dX, dY) )
            {
                // Point is inside, so "erase" the value
                (void) (poImage->*poImage->prnSetPixel)(nx, ny, fErasedValue);
            }
        }
    }

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
CCircularGonioShadowGeometry::CCircularGonioShadowGeometry(std::vector<double>* pdaShadowInfo, double daPixelSize[2]) :
CGonioShadowGeometry(pdaShadowInfo, daPixelSize),
m_dCenter_X_Coord(0.0),
m_dCenter_Y_Coord(0.0),
m_dRotated_Center_X_Coord(0.0),
m_dInnerRadius(0.0),
m_dOuterRadius(0.0)
{
    if( pdaShadowInfo )
    {
        int     nCircularAreaInfoVectorLength = pdaShadowInfo->size() - CGonioShadow::enShadow_coordinates_start;
    
        if( nCircularAreaInfoVectorLength >= 4 ) // 2 center coordinates and 2 radii  
        {
            m_dCenter_X_Coord = (*pdaShadowInfo)[CGonioShadow::enShadow_coordinates_start] * m_dPixelScaleFactor_0;
            m_dCenter_Y_Coord = (*pdaShadowInfo)[CGonioShadow::enShadow_coordinates_start + 1] * m_dPixelScaleFactor_1;
            
            if( m_dPixelScaleFactor_0 != m_dPixelScaleFactor_1 ) // assuming that pixels must be squares
            {
                cout << "\n Error applying circular shadow. Problem with scaling" << endl;
            }
            
            double  dFirstInputRadius = (*pdaShadowInfo)[CGonioShadow::enShadow_coordinates_start + 2] * m_dPixelScaleFactor_0;
            double  dSecondInputRadius = (*pdaShadowInfo)[CGonioShadow::enShadow_coordinates_start + 3] * m_dPixelScaleFactor_0; // assuming m_dPixelScaleFactor_0 == m_dPixelScaleFactor_1
            
            m_dInnerRadius = min(dFirstInputRadius, dSecondInputRadius);
            m_dOuterRadius = max(dFirstInputRadius, dSecondInputRadius);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CCircularGonioShadowGeometry::bIsAvailable()
{
    if( m_dInnerRadius >= m_dOuterRadius )
        return false;

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CCircularGonioShadowGeometry::bErasePixels(Cimage* poImage)
{
    if( !poImage )
        return false;
	      
    // Zero out the pixels within the circular band between the inner and outer radii.
    // The center does not need to be within the image

    int     nInnerRadius = (int)m_dInnerRadius;
    int     nOuterRadius = (int)m_dOuterRadius;
    
    int     nPxCen0 = (int)m_dRotated_Center_X_Coord;
    int     nPxCen1 = (int)m_dCenter_Y_Coord;

    double      dOuterRadSq = (double)(nOuterRadius * nOuterRadius);
    double      dInnerRadSq = (double)(nInnerRadius * nInnerRadius);

    // Calculate some pixel limits to prevent looking at pixels not on the image
    /////////////////////////////////////////////////
    int         nJmin = nPxCen1 - nOuterRadius - 1;
    if(0 > nJmin) 
        nJmin = 0;
    if(poImage->nGetDimension(1) <= nJmin)
        nJmin = poImage->nGetDimension(1)-1;
    /////////////////////////////////////////////////

    int         nJmax = nPxCen1 + nOuterRadius + 1;
    if(0 > nJmax) 
        nJmax = 0;
    if(poImage->nGetDimension(1) <= nJmax)
        nJmax = poImage->nGetDimension(1)-1;
    ////////////////////////////////////////////////

    int         nImin = nPxCen0 - nOuterRadius - 1;
    if(0 > nImin) 
        nImin = 0;
    if(poImage->nGetDimension(0) <= nImin) 
        nImin = poImage->nGetDimension(0)-1;
    /////////////////////////////////////////////////

    int         nImax = nPxCen0 + nOuterRadius + 1;
    if(0 > nImax)
        nImax = 0;
    if(poImage->nGetDimension(0) <= nImax)
        nImax = poImage->nGetDimension(0)-1;
    /////////////////////////////////////////////////

    double      dJSq = 0.0;
    double      dISq = 0.0;
    double      dTemp = 0.0;
    
    float       fErasedValue = 0.0;

    for(int jj = nJmin; jj <= nJmax; jj++)
    {
        dJSq = (double)jj + 0.5 - m_dCenter_Y_Coord;  // Use pix j center: (j+0.5)
        dJSq = dJSq * dJSq;
        
        for (int ii = nImin; ii <= nImax; ii++)
        {
	        dISq = (double)ii + 0.5 - m_dRotated_Center_X_Coord;  // Use pix i center: (i+0.5)
	        dTemp = dJSq + dISq * dISq;
	        
            if( dOuterRadSq >= dTemp && dInnerRadSq <= dTemp )
                (void) (poImage->*poImage->prnSetPixel)(ii, jj, fErasedValue);
	    }
    }

    // And last, but not least: always include the center ...
    (void)(poImage->*poImage->prnSetPixel)(nPxCen0, nPxCen1, fErasedValue);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CCircularGonioShadowGeometry::bRotate(double dRotationOffset, double dImageWidth)
{
    double      dPx0Offset = m_dOmegaOffset - dRotationOffset;
    
    dPx0Offset *= m_dPixelsPerDegree * m_dPixelScaleFactor_0;

    double      dPx0Min =  999999.0;
    double      dPx0Max = -999999.0;

    /////////////////////////////////////////////////////////////
    double      dPx0 = m_dCenter_X_Coord + dPx0Offset;

    if( dPx0Min > dPx0 ) 
      dPx0Min = dPx0;

    if( dPx0Max < dPx0 ) 
      dPx0Max = dPx0;

    m_dRotated_Center_X_Coord = dPx0;
    /////////////////////////////////////////////////////////////
    
    //	  //cout << "\n\n";
    if( 20.0 < m_dPixelsPerDegree && 23.0 > m_dPixelsPerDegree )
    {
        // It is probably a shadow on a RAPID, so make sure it loops
        // around to other edge of image if required

        // RAPID pixels usually go from 0 to 4599 along Px0. 
        // A full 360 degrees will be 7980 pixels.
        // That is 22.16667 pixels per degree
        // So an object must fit between 0 and 4599 pixels

        //  There are 4 possibilities:
        //  a  Both entirely outside the image
        //  b. Polygon min/max are both entirely within the image
        //  c. Polygon min is in the image, but quad max is out of the image (on left)
        //  d. Polygon max is in the image, but quad min is out of the image (on right)

        // There is a periodicity of 7980 pixels (per 360 degrees) for
        // RAPID images

	      double        dMaxDim0 = 4600.0 * m_dPixelScaleFactor_0;
	      double        dAdd = 0.0;
	      double        dSub = 0.0;

	      // Shadows are defined for the full RAPID horizontal 
	      // readout of 4600 pixels.  If this is a narrow readout of 
	      // 2800 pixels
	      // then shift the shadows horizontal by 1800 pixels as the first 
	      // thing.

	      if( 2800 * m_dPixelScaleFactor_0 == dImageWidth )
          {
                dAdd = -1800.0 * m_dPixelScaleFactor_0;
          }

        //////////////////////////////////////////////////////////////
        // Get the min/max bounds of the circle (not just the center)
        
        // Ask JIM why!
          
        dPx0Min = dPx0Min - m_dInnerRadius;
        dPx0Max = dPx0Min + m_dInnerRadius;
	    ///////////////////////////////////////////////////////////////

        while( 0.0 > dPx0Min && 0.0 > dPx0Max )
		{
            // If both outside the image on the low side, try to see if
            //    inside the image due to periodicity
            dAdd    += 7980.0 * m_dPixelScaleFactor_0;
            dPx0Min += 7980.0 * m_dPixelScaleFactor_0;
            dPx0Max += 7980.0 * m_dPixelScaleFactor_0;
		}
		  
        while( dMaxDim0 <= dPx0Min && dMaxDim0 <= dPx0Max )
		{
		  // If both outside the image on the high side, try to see if
		  //    inside the image due to periodicity
            dSub    -= 7980.0 * m_dPixelScaleFactor_0;
            dPx0Min -= 7980.0 * m_dPixelScaleFactor_0;
            dPx0Max -= 7980.0 * m_dPixelScaleFactor_0;
		}
	    
        dAdd = dAdd + dSub;
        //cout << "dAdd is " << dAdd << endl;

        // We do not have to worry about a shadow or points being off 
          // the image.
        // The ::vErasePixels*() routines will make sure we do not 
          // go out-of-bounds

		  // Shift first coordinate of each point by dAdd
		  m_dRotated_Center_X_Coord += dAdd;
    }

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
static int nPointInPolygon(int nPoints, double *xp, double *yp, double x, double y)
/*
//    int pnpoly(int npol, float *xp, float *yp, float x, float y)

License to Use the pnpoly routine from Wm. Randolph Franklin
(This applies to the CGonioMask::nPointInPolygon() routine only - jwp)

Copyright (c) 1970-2003, Wm. Randolph Franklin 

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

Redistributions of source code must retain the above
copyright notice, this list of conditions and the following
disclaimers.

Redistributions in binary form must reproduce the above
copyright notice in the documentation and/or other materials
provided with the distribution.

The name of Wm. Randolph Franklin may not be used to
endorse or promote products derived from this Software without
specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
{
    int i, j, c = 0;
    
    for(i = 0, j = nPoints-1; i < nPoints; j = i++)
    {
        if (  (   ((yp[i]<=y) && (y<yp[j]))
	     || ((yp[j]<=y) && (y<yp[i])))
	    &&
	    (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i])
          )
	    c = !c;
    }
    
    return c;
}

////////////////////////////////////////////////////////////////////////
// Make sure the angle is between 0 and 360
void  vCheckTransformAngleValue(double& dAngleValue)
{
    if( dAngleValue >= 0.0 && dAngleValue <= 360.0 )
        return;

    if( dAngleValue < 0.0 && dAngleValue >= -360.0 )
    {
        
        dAngleValue += 360.0;

        return;
    }
    
    cout << "\nInvalid shadow angle value" << dAngleValue << endl;
}
//////////////////////////////////////////////////////////////////////////
