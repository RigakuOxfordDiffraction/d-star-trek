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
// Cspatial.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Cspatial which implements
//    the spatial distortion encapsulation of d*TREK.
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
//   Cspatial.cc implements a simple pixel scaling spatial distortion and
//   a more complex spatial distortion developed by Dr. Marty Stanton of
//   Brandeis University.
//
//+ToDo
//
//   Error messages need implementing
//   m_a2x2fDirVec needs proper treatment.
//
//+Include files

#include "Cspatial.h"         // Class definition and prototypes
                              //  Cspatial.h includes Cimage.h, Cimage_header.h
#include "dtrekdefs.h"
#include "ImagePool.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


//+Definitions, constants, and initialization of static member variables

Cstring Cspatial::ms_sSpatialKey               = D_K_SpatialKey;
Cstring Cspatial::ms_sSpatialDistortionType    = D_K_SpatialDistortionType ;
Cstring Cspatial::ms_sSpatialDistortionInfo    = D_K_SpatialDistortionInfo;
Cstring Cspatial::ms_sSpatialDistortionVectors = D_K_SpatialDistortionVectors;
Cstring Cspatial::ms_sSpatialBeamPosn          = D_K_SpatialBeamPosn ;
Cstring Cspatial::ms_sSpatialTypeUnknown       = D_K_SpatialTypeUnknown;
Cstring Cspatial::ms_sSpatialTypeSimple        = D_K_SpatialTypeSimple ;
Cstring Cspatial::ms_sSpatialTypeComplex       = D_K_SpatialTypeComplex;
Cstring Cspatial::ms_sSpatialTypeNone          = D_K_SpatialTypeNone;
Cstring Cspatial::ms_sSpatialTypeInterp        = D_K_SpatialTypeInterp;
Cstring Cspatial::ms_sCalpar                   = ".calpar";
Cstring Cspatial::ms_sX_int                    = ".x_int";
Cstring Cspatial::ms_sY_int                    = ".y_int";
Cstring Cspatial::ms_sInv_x_int                = ".inv_x_int";
Cstring Cspatial::ms_sInv_y_int                = ".inv_y_int";

Cstring Cspatial::ms_sPixelSize                = D_K_CalibPixelSize;
Cstring Cspatial::ms_sPscale                   = D_K_CalibPscale;
Cstring Cspatial::ms_sXintStart                = D_K_CalibXintStart;
Cstring Cspatial::ms_sYintStart                = D_K_CalibYintStart;
Cstring Cspatial::ms_sXinvStart                = D_K_CalibXinvStart;
Cstring Cspatial::ms_sYinvStart                = D_K_CalibYinvStart;
Cstring Cspatial::ms_sXintStep                 = D_K_CalibXintStep;
Cstring Cspatial::ms_sYintStep                 = D_K_CalibYintStep;
Cstring Cspatial::ms_sXinvStep                 = D_K_CalibXinvStep;
Cstring Cspatial::ms_sYinvStep                 = D_K_CalibYinvStep;
Cstring Cspatial::ms_sXSize                    = D_K_CalibXSize;
Cstring Cspatial::ms_sYSize                    = D_K_CalibYSize;
Cstring Cspatial::ms_sXBeam                    = D_K_CalibXBeam;
Cstring Cspatial::ms_sYBeam                    = D_K_CalibYBeam;
Cstring Cspatial::ms_sBadFlag                  = D_K_CalibBadFlag;

const double cdCylindricalDetectorRadius_MM = 127.4;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cspatial::Cspatial()
{
 (void) nInitValues();
}

Cspatial::Cspatial(const Cstring& sFileset)
{
  (void) nInitValues();   // Direction vectors
  (void) nRead(sFileset);
}

Cspatial::Cspatial(const float fBeam0, const float fBeam1,
		   const float fPxSize0, const float fPxSize1,
		   const int nDim0, const int nDim1,
		   const float fDir00, const float fDir01,
		   const float fDir10, const float fDir11)
{
  (void) nInitValues(fBeam0, fBeam1, fPxSize0, fPxSize1, nDim0, nDim1,
		     fDir00, fDir01, fDir10, fDir11);
}

Cspatial::Cspatial(Cimage_header& oHeader, const Cstring& sPre)
{
  Cstring sTemp;
  Cstring sPrefix;
  int nStat;
  int nTemp;
  float fValues[4];

  (void) nInitValues();
  nStat   = 0;
  sPrefix = sPre;
  if ("" == sPrefix)
    {
      // Null prefix string, see if the Cdetector::ms_sDetectorNames keyword
      // exists in the header.  If so, change the prefix to the first name
      // in the list of names.  This call is OK, since if there is no
      // keyword, then sPrefix is set to "" anyways.

      (void) oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1, &sPrefix);
    }

  if (!oHeader.bIsAvailable())
    {
      cout << "Cspatial ERROR: image header is not valid. "
	   << " Cannot construct spatial distortion object!\n";
      nStat = 1;
    }
  else if (0 == oHeader.nGetValue(sPrefix + ms_sSpatialKey, &sTemp))
    {
      cout << "Cspatial ERROR: cannot use database key yet!\n";
      nStat = 2;
    }
  	
  // Try to get required information from the header
  int nDim0, nDim1;  // Dimensions of the image in the header

  nDim0 = 1;
  nDim1 = 1;
  nStat = nStat + oHeader.nGetValue(Cimage_header::ms_sSize1, &nDim0);
  nDim0--;
  nStat = nStat + oHeader.nGetValue(Cimage_header::ms_sSize2, &nDim1);
  nDim1--;

  // Look for direction vectors, but if they do not exist, just ignore
  // This MUST be done before a call to Cspatial::nRead.

  (void) oHeader.nGetValue(sPrefix + ms_sSpatialDistortionVectors, 4,
			   &m_a2x2fDirVec[0][0]);

/*
  cout << "Spatial distortion vectors after header read: " << nStat
       << '\n'
       << m_a2x2fDirVec[0][0] << ", "
       << m_a2x2fDirVec[0][1] << ", "
       << m_a2x2fDirVec[1][0] << ", "
       << m_a2x2fDirVec[1][1] << endl;
*/
  vInvert();

  nTemp = oHeader.nGetValue(sPrefix + ms_sSpatialDistortionType, &sTemp);
  //cout << "In Cspatial, spatial type is : " << sTemp << endl;
  sTemp = sTransSymbol(sTemp);
  if (0 != nTemp)
    {
      cout << "Cspatial ERROR: reading " << sPrefix
           << ms_sSpatialDistortionType << " keyword!\n";
      nStat++;
    }
  else
    {
      if (ms_sSpatialTypeUnknown == sTemp)
	{
	  nStat++;
	}
      else if (ms_sSpatialTypeSimple == sTemp)
	{
	  nTemp = oHeader.nGetValue(sPrefix + ms_sSpatialDistortionInfo,
				    4, fValues);
/*
	  cout << "In Cspatial, spatial values are : "
               << fValues[0] << ", " << fValues[1] << ", "
               << fValues[2] << ", " << fValues[3] << endl;
*/
	  if (0 == nTemp)
	    {
	      (void) nInitValues(fValues[0], fValues[1],
				 fValues[2], fValues[3],
				 nDim0, nDim1,
				 m_a2x2fDirVec[0][0],
				 m_a2x2fDirVec[0][1],
				 m_a2x2fDirVec[1][0],
				 m_a2x2fDirVec[1][1]
				 );
	    }
	  else
	    {
	      cout << "Cspatial: error reading simple spatial distortion from header!\n";
	      nStat++;
	    }	
	}
      else if (ms_sSpatialTypeComplex == sTemp)
	{
	  cout << "Cspatial: do not know how to create complex spatial "
	       << " object!\n";
	  nStat++;
	}	
      else if (ms_sSpatialTypeInterp == sTemp)
	{
	  nTemp = oHeader.nGetValue(sPrefix + ms_sSpatialDistortionInfo,
				    &sTemp);
	  if (0 == nTemp)
	    {
	      nTemp = nRead(sTemp);
	      if (0 != nTemp)
		{
		  cout << "Cspatial: error reading spatial distortion file(s)!\n";
		  nStat++;
		}
	    }
	  else
	    {
	      cout << "Cspatial: error reading spatial distortion fileset "
                   << "name from header!\n";
	      nStat++;
	    }	
	}
      else
	{
	  cout << "Cspatial: invalid " << sPrefix
               << ms_sSpatialDistortionInfo << " keyword in header!\n";
	  nStat++;
	}
    }

  // Now look for direct beam information
  
  nTemp = oHeader.nGetValue(sPrefix + ms_sSpatialBeamPosn, 2, m_a2fBeamPx);
  if (0 == nTemp)
  {
      nTemp = nSetBeamPosition(m_a2fBeamPx[0], m_a2fBeamPx[1]);
  }
  
  if (0 != nStat)
  {
      // This is the spatial distortion fallback.
      nInitValues((float)nDim0 * 0.5f, 
                  (float)nDim1 * 0.5f, 
                  0.1f, 
                  0.1f,
                  nDim0,
                  nDim1);
  }
  else
    {
      m_eThe_State          = eSpatial_available_state;
    }
  // Look for the geometry type.  We are looking at the detector name.
  // or detector type
  // If the detector name is D19_, then this is the cylindrical
  // detector at ILL D19.

  nTemp = oHeader.nGetValue(D_K_DetectorType,&sTemp);
  
  if( 0 == nTemp && sTemp.contains("R-AXIS-CS") )
  {
      // It's a rapid.  We need to get the distance.
      
      Cstring       strPrefix("");
      if( 0 != oHeader.nGetValue(D_K_DetectorNames, 1, &strPrefix) )
            strPrefix = "";

      Cgoniometer   oGon(oHeader, strPrefix);
      if( oGon.bIsAvailable() )
      {
          //m_fRadius = oGon.fGetDistance();
          m_fCylindricalDetectorShift = oGon.fGetDistance() - m_fRadius;
          
          // sanity check
          if( m_fRadius <= 0.0f )
              nStat++;
          
          if( m_fCylindricalDetectorShift < 0.0f && fabs((double)m_fCylindricalDetectorShift) > m_fRadius )
              nStat++; // The cylinder is shifted towards the crystal way too much.
      }
      else
          nStat++;
      
      m_eThe_Geometry = eSpatial_cylindrical_geometry;
//      cout << "GEOMETRY is CYLIND\n" << flush;
  }

}

Cspatial::~Cspatial()
{
   if (NULL != m_poCalpar)
     {
       delete m_poCalpar;
       m_poCalpar = NULL;
     }
   
   if (NULL != m_poX_int)
     {
       m_pImagePool->bFreeImage(m_poX_int->sGetName());
       m_poX_int = NULL;
     }
   
   if (NULL != m_poY_int)
     {
       m_pImagePool->bFreeImage(m_poY_int->sGetName());
       m_poY_int = NULL;
     }
   
   if (NULL != m_poInv_x_int)
     {
       m_pImagePool->bFreeImage(m_poInv_x_int->sGetName());
       m_poInv_x_int = NULL;
     }
   
   if (NULL != m_poInv_y_int)
     {
       m_pImagePool->bFreeImage(m_poInv_y_int->sGetName());
       m_poInv_y_int = NULL;
     }
}

int Cspatial::nInitValues(const Cstring& sFileset)
{
//  (void) nInitValues();  // Danger! Changes direction vectors
  cout << "ERROR in Cspatial::nInitValues(const Cstring& sFileset)! DO NOT USE!\n";
  (void) nRead(sFileset);
  return 0;
}

int Cspatial::nInitValues(const float fBeam0, const float fBeam1,
			  const float fPxsize0, const float fPxsize1,
			  const int nDim0, const int nDim1,
			  const float fDir00, const float fDir01,
			  const float fDir10, const float fDir11)
{
  (void) nInitValues();
  m_a2x2fDirVec[0][0]   = fDir00;
  m_a2x2fDirVec[0][1]   = fDir01;
  m_a2x2fDirVec[1][0]   = fDir10;
  m_a2x2fDirVec[1][1]   = fDir11;
  vInvert();

  // Can beam position be off of detector?  
  // Can it be anywhere if Simple_spatial?

//  if (   (0.0 <= fBeam0)
//      && (0.0 <= fBeam1)
//      && (0.0 < fPxsize0)
  if (   (0.0 < fPxsize0)
      && (0.0 < fPxsize1) )
    {
      m_a2fBeamPx[0]        = fBeam0;
      m_a2fBeamPx[1]        = fBeam1;
      m_a2fNominalPxSize[0] = fPxsize0;
      m_a2fNominalPxSize[1] = fPxsize1;
      (void) nSetMaxPixel(nDim0, nDim1);
      m_eThe_Type           = eSpatial_simple_type;
      m_eThe_State          = eSpatial_available_state;
      m_sDescription        = "Simple pixel scaling spatial distortion";
    }
  else
    {
      m_eThe_Type  = eSpatial_unknown_type;
      m_eThe_State = eSpatial_unknown_state;
    }
  return (0);
}

int Cspatial::nInitValues(void)
{
  m_eThe_State          = eSpatial_unknown_state;
  m_eThe_Type           = eSpatial_unknown_type;
  m_eThe_Geometry       = eSpatial_flat_geometry;
//
  m_a2fBeamPx[0]        = 0.0;
  m_a2fBeamPx[1]        = 0.0;
  m_a2fNominalPxSize[0] = 0.1f;
  m_a2fNominalPxSize[1] = 0.1f;
  m_a2fMaxPixel[0]      = 3071.;
  m_a2fMaxPixel[1]      = 3071.;
  
  //m_fRadius             = 0.0;
  double        dTemp = 0.0;
  if( bGetEnv(D_K_CylindricalDetectorRadius, dTemp) )
    m_fRadius = (float)dTemp;
  else 
    m_fRadius = (float)cdCylindricalDetectorRadius_MM;

  m_fCylindricalDetectorShift = 0.0f;
  //
  
  m_a2x2fDirVec[0][0]   = 1.0;
  m_a2x2fDirVec[0][1]   = 0.0;
  m_a2x2fDirVec[1][0]   = 0.0;
  m_a2x2fDirVec[1][1]   = 1.0;
  m_a2x2fInvVec[0][0]   = 1.0;
  m_a2x2fInvVec[0][1]   = 0.0;
  m_a2x2fInvVec[1][0]   = 0.0;
  m_a2x2fInvVec[1][1]   = 1.0;
  m_sFileset_name       = "distor";
  m_sDescription        = "No spatial distortion information";
  m_poCalpar            = NULL;
  m_poX_int             = NULL;
  m_poY_int             = NULL;
  m_poInv_x_int         = NULL;
  m_poInv_y_int         = NULL;


  m_nNumErrorOut       = 0;
  m_nVerbose           = 0;
  
  m_pImagePool =  CImagePool::poGetInstance();
  
  return (0);

}

int Cspatial::nRead_BS(const Cstring& sFileset)
{
    // Read the P4P file
    Cstring sTemp;
    int nx;
    float f0;
    char* pcTemp;
    
    FILE* pFIn;
    int nStat;
    const int nBufLen = 200;
    char pcBuf[nBufLen];
    
    m_poCalpar      = NULL;
    m_poX_int       = NULL;
    m_poY_int       = NULL;
    m_poInv_x_int   = NULL;
    m_poInv_y_int   = NULL;
    
    
    pFIn = fopen(sFileset.string(),"rt");
    if (!pFIn)
        return 1;
    nStat = 1;
    while (fgets(&pcBuf[0],nBufLen-1,pFIn)) {
        sTemp = pcBuf;
        if ((sTemp.find("DATA")==0) && (sTemp.find("SPATIAL") != -1)) {
            for (pcTemp = strtok(sTemp.string(),"\t\n\r "),nx=0;(pcTemp) && (nx<10-1); pcTemp = strtok(NULL,"\t\n\r "),nx++);            
            if (!pcTemp)
                pcTemp = "";
            sTemp = pcTemp;
            if (sTemp.bParseFloat(f0)) {
                m_n0ImageSize = (int) f0;
                m_n1ImageSize = (int) f0;
                nStat = 0;
            } else
                nStat = 1;                            
            break;
        };
    };
    if (nStat) {
        fclose(pFIn);
        return 1;
    };
    
    int         nInterpSize = 56;
    int         nFile;
    int         nInterp0,nInterp1;
    int         nLinePoint = 0;
    int         nNum,nSign;
    int         nDigits;
    Cimage*     poImage = NULL;
    
    
    pcBuf[0] = 0;
    m_poCalpar      = new Cimage_header;
    
    nStat = 0;
    for (nFile = 0; (nFile < 4) && (!nStat); nFile++) {
        switch (nFile) {
            
        case 2:
            poImage = m_poX_int       = m_pImagePool->poGetImage(nInterpSize, nInterpSize, eImage_uI4, ms_sX_int);
            break;
        case 3:
            poImage = m_poY_int       = m_pImagePool->poGetImage(nInterpSize, nInterpSize, eImage_uI4, ms_sY_int);
            break;
        case 0:
            poImage = m_poInv_x_int   = m_pImagePool->poGetImage(nInterpSize, nInterpSize, eImage_uI4, ms_sInv_x_int);
            break;
        case 1:
            poImage = m_poInv_y_int   = m_pImagePool->poGetImage(nInterpSize, nInterpSize, eImage_uI4, ms_sInv_y_int);
            break;
        };        
        
        for (nInterp1 = 0;(!nStat) && (nInterp1 < nInterpSize); nInterp1++) {
            for (nInterp0 = 0;nInterp0 < nInterpSize; nInterp0++) {
                nNum = 0;
                nSign = 1;
                nDigits = 0;
                while (1) {
                    if (!pcBuf[nLinePoint]) {
                        if (!fgets(&pcBuf[0],nBufLen-1,pFIn)) {
                            nStat = 1;
                            break;
                        };
                        nLinePoint = 0;
                    };
                    
                    if (pcBuf[nLinePoint]=='+') {
                    } else if (pcBuf[nLinePoint]=='-') {
                        nSign = -1;
                    } else if ((pcBuf[nLinePoint]>='0') && (pcBuf[nLinePoint]<='9')) {
                        nNum = nNum*10 + (pcBuf[nLinePoint]-'0');
                        nDigits++;
                    } else if  (
                        (pcBuf[nLinePoint]=='\n') || 
                        (pcBuf[nLinePoint]=='\t') || 
                        (pcBuf[nLinePoint]=='\r') ||
                        (pcBuf[nLinePoint]==' ')) {
                        if (nDigits)
                            break;
                    } else {
                        nStat = 1;
                        break;
                    };
                    nLinePoint++;
                } 
                poImage->nSetPixel(nInterp0,nInterp1,(LONG) (nNum*nSign));
                
            };
        };     
    };

    fclose(pFIn);
    return 0;   
};

int
Cspatial::nRead(const Cstring& sFileset)
{
    int nStat;
    Cstring sTemp;
    m_eThe_State    = eSpatial_unknown_state;   // Default state is unknown
    m_eThe_Type     = eSpatial_interp_table;
    m_sDescription  = "Read from fileset: " + sFileset;
    
    cout << "...reading spatial distortion tables..." << endl;
    
    m_sFileset_name = sFileset;
    sTemp = sFileset;
    sTemp.upcase();
    if (sTemp.contains("P4P") || (sTemp.contains("SPIN"))) {
        nStat = nRead_BS(sFileset);
        
        switch (m_n0ImageSize) {
        case 512:
            m_n0Step = 10;
            m_n1Step = 10;
            break;
        case 1024:
            m_n0Step = 20;
            m_n1Step = 20;
            break;
        case 2048:
            m_n0Step = 40;
            m_n1Step = 40;
            break;
        default:
            nStat = 1;
        };
        
        m_n0InvStep = m_n0Step;
        m_n1InvStep = m_n1Step;
        m_n0Start = -m_n0Step*2;
        m_n1Start = -m_n1Step*2;
        m_n0InvStart = -m_n0InvStep*2;
        m_n1InvStart = -m_n1InvStep*2;
        m_fPxScale = 0.001f;
        // m_fPxSize = 0.1094;
        m_fPxSize = 0.122f;
	if ("" != sGetEnv("BS_PIXEL_SIZE"))
	  m_fPxSize = atof(sGetEnv("BS_PIXEL_SIZE").string());
	cout << "BS pixel size is: " << m_fPxSize << endl << flush;
        m_a2fMMOriginPx[0] = 100;
        m_a2fMMOriginPx[1] = 400;

    }
    else
    { 
        m_sFileset_name = sFileset;
        
        m_poCalpar      = new Cimage_header(sFileset+ms_sCalpar);
        
        m_poX_int       = m_pImagePool->poGetImage(sFileset+ms_sX_int);
        m_poY_int       = m_pImagePool->poGetImage(sFileset+ms_sY_int);
        m_poInv_x_int   = m_pImagePool->poGetImage(sFileset+ms_sInv_x_int);
        m_poInv_y_int   = m_pImagePool->poGetImage(sFileset+ms_sInv_y_int);
        
        if (   (!m_poCalpar->bIsAvailable())
            || (!m_poX_int->bIsAvailable())
            || (!m_poY_int->bIsAvailable())
            || (!m_poInv_x_int->bIsAvailable())
            || (!m_poInv_y_int->bIsAvailable()) )
            nStat = -5;
        else
            nStat = 0;
        
        cout << endl;
        if (!nStat) {
            // Get from the .calpar file some necessary items
            
            int nStatL;
            Cstring sError = "ERROR reading keyword: ";
            nStat = 0;
            nStatL = m_poCalpar->nGetValue(ms_sPixelSize, &m_fPxSize);
            if (0 != nStatL)
                cout << sError << ms_sPixelSize << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sPscale,    &m_fPxScale);
            if (0 != nStatL)
                cout << sError << ms_sPscale << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sXintStart, &m_n0Start);
            if (0 != nStatL)
                cout << sError << ms_sXintStart << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sYintStart, &m_n1Start);
            if (0 != nStatL)
                cout << sError << ms_sYintStart << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sXinvStart, &m_n0InvStart);
            if (0 != nStatL)
                cout << sError << ms_sXinvStart << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sYinvStart, &m_n1InvStart);
            if (0 != nStatL)
                cout << sError << ms_sYinvStart << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sXintStep,  &m_n0Step);
            if (0 != nStatL)
                cout << sError << ms_sXintStep << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sYintStep,  &m_n1Step);
            if (0 != nStatL)
                cout << sError << ms_sYintStep << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sXinvStep,  &m_n0InvStep);
            if (0 != nStatL)
                cout << sError << ms_sXinvStep << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sYinvStep,  &m_n1InvStep);
            if (0 != nStatL)
                cout << sError << ms_sYinvStep << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sXSize,     &m_n0ImageSize);
            if (0 != nStatL)
                cout << sError << ms_sXSize << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sYSize,     &m_n1ImageSize);
            if (0 != nStatL)
                cout << sError << ms_sYSize << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sXBeam,
                &m_a2fMMOriginPx[0]);
            if (0 != nStatL)
                cout << sError << ms_sXBeam << endl;
            nStat += nStatL;
            m_a2fNominalPxSize[0] = m_fPxSize;
            m_a2fNominalPxSize[1] = m_fPxSize;
            
            
            nStatL = m_poCalpar->nGetValue(ms_sYBeam,
                &m_a2fMMOriginPx[1]);
            if (0 != nStatL)
                cout << sError << ms_sYBeam << endl;
            nStat += nStatL;
            nStatL = m_poCalpar->nGetValue(ms_sBadFlag,   &m_nBadFlag);
            if (0 != nStatL)
                cout << sError << ms_sBadFlag << endl;
            nStat += nStatL;
        };
    }
    
    //      cout << "Spatial nStat is: " << nStat << endl;
    // Some of the above can be transferred to fMaxPixel, fPxBeam, m_a2fNominal...
    
    // Check that certain important items do not have invalid values
    if (   (0.0  >= m_fPxSize)
        || (-100 >= m_n0Step)
        || (-100 >= m_n1Step)
        || (-100 >= m_n0InvStep)
        || (-100 >= m_n1InvStep)
        )
    {
        cout << "ERROR Cspatial::nRead() with pixel_size, or step(s)!" << endl;
        nStat++;
    }
    
    if (0 == nStat)
    {
        m_a2fBeamPx[0] = m_a2fMMOriginPx[0];
        m_a2fBeamPx[1] = m_a2fMMOriginPx[1];

        nStat = nSetBeamPosition(m_a2fBeamPx[0], m_a2fBeamPx[1]);
        if (0 != nStat)
        {
            cout << "ERROR Cspatial::nRead() setting direct beam position." << endl;
        }
    }
    
    // Do some clean-up if there were any errors.
    
    if (0 != nStat)
    {
        m_eThe_Type  = eSpatial_unknown_type;

        if( m_poCalpar )
            delete  m_poCalpar;
        
        if( m_poX_int )
            m_pImagePool->bFreeImage(m_poX_int->sGetName());
        
        if( m_poY_int )
            m_pImagePool->bFreeImage(m_poY_int->sGetName());

        if( m_poInv_x_int )
            m_pImagePool->bFreeImage(m_poInv_x_int->sGetName());
      
        if( m_poInv_y_int )
            m_pImagePool->bFreeImage(m_poInv_y_int->sGetName());
        
        m_poCalpar    = NULL;
        m_poX_int     = NULL;
        m_poY_int     = NULL;
        m_poInv_x_int = NULL;
        m_poInv_y_int = NULL;
    }
    else
    {
        m_eThe_State = eSpatial_available_state;
    }
    return (nStat);
}

int
Cspatial::nSetBeamPosition(const float fBeamPx0, const float fBeamPx1)
{
    int nStat = 0;
    if (  (eSpatial_simple_type == m_eThe_Type)
        || (    (0.0 <= fBeamPx0)
        && (0.0 <= fBeamPx1) ) )
    {
        m_a2fBeamPx[0] = fBeamPx0;
        m_a2fBeamPx[1] = fBeamPx1;
        if (eSpatial_interp_table == m_eThe_Type)
        {
            // Compute new m_fXCenter and m_fYCenter ...
            
            m_fXCenter   = 0.0;
            m_fYCenter   = 0.0;
            float f0tmp, f1tmp;
            
            // Input spd fileset has mm (i.e. beam) origin off by 1 pixel,
            // so subtract 1 from each
            
            m_a2fMMOriginPx[0] = m_a2fBeamPx[0];
            m_a2fMMOriginPx[1] = m_a2fBeamPx[1];
            
            // Not sure the -1.0's are needed in next function call
            
            nStat              = nStantonPixeltoMM(
                m_a2fMMOriginPx[0]-1.0f,
                m_a2fMMOriginPx[1]-1.0f,
                &f0tmp, &f1tmp);

              
            if (0 == nStat)
            {
                // OK, successful, so set new m_fXCenter and m_fYCenter
                
                m_fXCenter = f0tmp;
                m_fYCenter = f1tmp;
                
                //	      cout << "m_fX,YCenter: " << m_fXCenter << ", " << m_fYCenter << endl;
                // Test the reverse transformation
                
                nStat = nStantonMMtoPixel(0.0, 0.0, &f0tmp, &f1tmp);

                float fEPS = 0.05f;
                if (   (0 == nStat)
                    && ((float)fabs((double)(f0tmp - (m_a2fMMOriginPx[0]-1.0))) <= fEPS)

                    && ((float)fabs((double)(f1tmp - (m_a2fMMOriginPx[1]-1.0))) <= fEPS) )
                    m_eThe_State = eSpatial_available_state;
                
                // Set some nominal values based on ...
                
                else
                {

                    m_eThe_State = eSpatial_unknown_state;
                    cout << "ERROR: Bad spatial distortion tables: "
                        << "PX -> MM -> PX operations are inconsistent!\n";
                    if ("" != sGetEnv("PX210_SPD_BADTABLE_OVERRIDE"))
                    {
                        cout << "WARNING: You requested to override the ERROR!\n";
                        m_eThe_State = eSpatial_available_state;
                    }
                }
            }
        }
    }
    else
    {
        m_eThe_State = eSpatial_unknown_state;
        nStat = -1;
    }
    
    return (nStat);
}

int Cspatial::nSetScaleFactors(const float fScalePx0, const float fScalePx1)
{
  if (   (0.0 < fScalePx0)  // Do not allow this if interp_table
      && (0.0 < fScalePx1) )
    {
      m_a2fNominalPxSize[0] = fScalePx0;
      m_a2fNominalPxSize[1] = fScalePx1;
      m_fPxSize             = fScalePx0;
      return (0);
    }
  else
    {
      m_eThe_State = eSpatial_unknown_state;
      return (1);
    }
}

int Cspatial::nPixeltoMM(const float fPx0, 
                         const float fPx1, 
			             float *pfXmm, 
                         float *pfYmm, 
                         float *pfZRelativeCoordinate)
{
  //RB: The ZRelativeCoordinate variable represents a scale factor in front of the vector Do, 
  // see  Messerschmidt.    November 1986, Proceedings Phase III page 59
  //
  // For the flat detector geometry ZRelativeCoordinate is simply unity.
  //
  // For the cylindrical geometry, assuming that Do is a vector connecting the crystal with the center of the detector,
  // 
  // ZRelativeCoordinate = (1 - (Rc - sqrt(Rc**2-Ymm**2)) / (Rc+shift)  ),
  // where Rc is the cylinder radius and "shift" is a shift of the detector along the Z-axis.
  // In the ideal case, when "shift" is zero, ZRelativeCoordinate is also a cosine of the angle between Do and the projection
  // of the scattered vector on plane YZ.
  //
    
  int   nStat;
  float f0tmp, f1tmp;

  *pfXmm = 0.1f; // The default Xmm position returned for all errors
  *pfYmm = 0.1f; // The default Ymm position returned for all errors
  *pfZRelativeCoordinate = 0.1f; // The default Ymm position returned for all errors

  if (eSpatial_available_state != m_eThe_State)
    {
      if (100 > m_nNumErrorOut)
	{
	  cout << "     Error in Cspatial::nPixeltoMM, spatial distortion info"
               << " not available!\n";
	  m_nNumErrorOut++;
	}
      return (1);

    }
  else if (   (0.0 > fPx0) || (0.0 > fPx1)
	   || (fPx0 > m_a2fMaxPixel[0]) || (fPx1 > m_a2fMaxPixel[1]) )
    {
      // Input pixels are out-of-bounds
      return (1);
    }
  else if (eSpatial_simple_type == m_eThe_Type)
    {
      //  Flat plate geometry computations.
      
      f0tmp = (fPx0 - m_a2fBeamPx[0]) * m_a2fNominalPxSize[0];
      f1tmp = (fPx1 - m_a2fBeamPx[1]) * m_a2fNominalPxSize[1];
      
      // Correct for m_a2x2fDirVec next (we are already in mm's)
      
      *pfXmm = m_a2x2fDirVec[0][0] * f0tmp  +  m_a2x2fDirVec[1][0] * f1tmp;
      *pfYmm = m_a2x2fDirVec[0][1] * f0tmp  +  m_a2x2fDirVec[1][1] * f1tmp;
      *pfZRelativeCoordinate = 1.0;

      // *pfXmm = *pfXmm - xcenter;
      // *pfYmm = *pfYmm - ycenter;

      if( eSpatial_cylindrical_geometry == m_eThe_Geometry )
      {
          // Adjust above for cylindrical (rapid) geometry.
          f0tmp = *pfYmm / m_fRadius;
          
          // X coordinate stays the same. (it is perpendicular to gonio base plate)
          
          // Y coordinate changes
          *pfYmm = m_fRadius * sin((double) f0tmp);          
      
          if( 0.0f != m_fCylindricalDetectorShift )
          {
              float   fSign = fabs(f0tmp) > Gs_fPI / 2.0f ? -1.0f : 1.0f;
              
              float   fSquareRootArgument = m_fRadius * m_fRadius - (*pfYmm) * (*pfYmm);
              if( fSquareRootArgument < 0.0f )
                  fSquareRootArgument = 0.0f;  // a safety feature against a crazy Y-value
              
              *pfZRelativeCoordinate = m_fCylindricalDetectorShift + fSign * sqrt((double)fSquareRootArgument);
              
              *pfZRelativeCoordinate /= (m_fRadius + m_fCylindricalDetectorShift);
          }
          else
          {
              // Z coordinate represents fraction of radius.
              *pfZRelativeCoordinate = cos((double) f0tmp);
          }
      }
      
      return (0);
    }
  else if (eSpatial_interp_table == m_eThe_Type)
    {
      nStat = nStantonPixeltoMM(fPx0, fPx1, &f0tmp, &f1tmp);

      // Correct for m_a2x2fDirVec next (we are already in mm's)
      // (Probably should do this only if nStat = 0)

      *pfXmm = f0tmp;
      *pfYmm = f1tmp;
      *pfZRelativeCoordinate = 1.0;

//      *pfXmm = m_a2x2fDirVec[0][0] * f0tmp  +  m_a2x2fDirVec[1][0] * f1tmp;
//      *pfYmm = m_a2x2fDirVec[0][1] * f0tmp  +  m_a2x2fDirVec[1][1] * f1tmp;

      return (nStat);
    }
  else
    {
      cout << "     Error in Cspatial::nPixeltoMM, unsupported spd type\n";
      return (2);
    }
}

int Cspatial::nMMtoPixel(float fXmm, 
                         float fYmm, 
                         float fZRelativeCoordinate,
			             float *pfPx0,
                         float *pfPx1,
                         bool bTruncateOutOfBounds)
{
      float f0tmp, f1tmp;

      *pfPx0 = 0.0; // The default 1st pixel coordinate returned for all errors
      *pfPx1 = 0.0; // The default 2nd pixel coordinate returned for all errors

      if (eSpatial_available_state != m_eThe_State)
        {
          if (100 > m_nNumErrorOut)
	    {
	      cout << "     Error in Cspatial::nMMtoPixel, spatial distortion info"
                   << " not available!\n";
	      m_nNumErrorOut++;
	    }
          return (1);
        }
      else if ( eSpatial_simple_type == m_eThe_Type )
        {
          if( eSpatial_cylindrical_geometry == m_eThe_Geometry )
          {
              // Adjust above for cylindrical (rapid) geometry.

              // X coordinate stays the same.
              // Y coordinate is back transformed into flat plate dimensions.
              if( 0.0f == m_fCylindricalDetectorShift )
              {
                  fYmm = ((fYmm>0.0)?(1.0):(-1.0)) * acos( max(-1.0,min(1.0,fZRelativeCoordinate))) * m_fRadius; 
              }
              else
              {
                  double    fSinus = fabs((double)fYmm) / (double)m_fRadius; 
                  double    fArcSinus = asin(min(1.0, fSinus));

                  if( fZRelativeCoordinate < m_fCylindricalDetectorShift / (m_fRadius + m_fCylindricalDetectorShift) )
                    fArcSinus += (Gs_dPI / 2.0);

                  float     fSign = fYmm > 0.0f ? 1.0f : -1.0f;

                  fYmm = fSign * (float)fArcSinus * m_fRadius; 
              }

              // Z coordinate is ignored.
          }

          // Convert mm to mm along pixel axes (is there a translation involved?)

          f0tmp = m_a2x2fInvVec[0][0] * fXmm  +  m_a2x2fInvVec[1][0] * fYmm;
          f1tmp = m_a2x2fInvVec[0][1] * fXmm  +  m_a2x2fInvVec[1][1] * fYmm;

          f0tmp =  (f0tmp / m_a2fNominalPxSize[0]) + m_a2fBeamPx[0];
          f1tmp =  (f1tmp / m_a2fNominalPxSize[1]) + m_a2fBeamPx[1];

          // Check for out-of-bounds conditions

          int nOut = 0;

          if (0.0 > f0tmp)
          {
              if (bTruncateOutOfBounds) {
                  f0tmp = 0.0;
                  nOut++;
              }
          }
          if (f0tmp > m_a2fMaxPixel[0])
          {
              if (bTruncateOutOfBounds) {
                  f0tmp = m_a2fMaxPixel[0];
                  nOut++;
              }
          }
          if (0.0 > f1tmp)
          {
              if (bTruncateOutOfBounds) {
                  f1tmp = 0.0;
                  nOut++;
              }
          }
          if (f1tmp > m_a2fMaxPixel[1])
          {
              if (bTruncateOutOfBounds) {
                  f1tmp = m_a2fMaxPixel[1];
                  nOut++;
              }
          }
          *pfPx0 = f0tmp;
          *pfPx1 = f1tmp;
          if (0 == nOut)
	    return (0);
          else
	    return (1);
    }
  else if (eSpatial_interp_table == m_eThe_Type)
    {
      return (nStantonMMtoPixel(fXmm, fYmm, pfPx0, pfPx1));
    }
  else
    {
      cout << "     Error in Cspatial::nMMtoPixel, unsupported spd type\n";
      return (2);
    }
}

int Cspatial::nList(void)
{
  cout << "    Spatial distortion descriptive text:\n";
  cout << "    " << m_sDescription << endl;

  if (eSpatial_available_state != m_eThe_State)
    {
      cout << "     Error in Cspatial::nList, spatial distortion info"
           << " not available!\n";
      return (1);
    }
  else
    {
      if (eSpatial_simple_type == m_eThe_Type)
	{
	  cout << "    Simple spatial distortion:\n";
	  cout << "    Center of primary beam: " << m_a2fBeamPx[0]
               << ", " << m_a2fBeamPx[1]
               << endl;

	  cout << "           Pixel size (mm): " << m_a2fNominalPxSize[0]
	                                 << ", " << m_a2fNominalPxSize[1]
               << endl;
	  cout << "         Direction vectors: "
	       << m_a2x2fDirVec[0][0] << ", " << m_a2x2fDirVec[0][1] << ", "
	       << m_a2x2fDirVec[1][0] << ", " << m_a2x2fDirVec[1][1] << endl;
	  return (0);
	}
      else if (eSpatial_interp_table == m_eThe_Type)
	{
	  cout << "    Simple spatial distortion:\n";
	  cout << "    Center of primary beam: " << m_a2fBeamPx[0]
               << ", " << m_a2fBeamPx[1]
               << endl;

	  cout << "           Pixel size (mm): " << m_a2fNominalPxSize[0]
	                                 << ", " << m_a2fNominalPxSize[1]
               << endl;
	  cout << " Interp Table Px size (mm): " << m_fPxSize
               <<"\nInterp Table Px scale:      " << m_fPxScale << endl;
	  cout << "         Direction vectors: "
	       << m_a2x2fDirVec[0][0] << ", " << m_a2x2fDirVec[0][1] << ", "
	       << m_a2x2fDirVec[1][0] << ", " << m_a2x2fDirVec[1][1] << endl;
	  return (0);  // Not yet implemented
	}
    else
      {
	cout << "     Error in Cspatial::nList, unsupported spd type\n";
	return (2);
      }
    }
}

int Cspatial::nSetMaxPixel(const int nMax0, const int nMax1)
{
  if ( (0 <= nMax0) && (0 <= nMax1) )
    {
      m_a2fMaxPixel[0] = (float) nMax0;
      m_a2fMaxPixel[1] = (float) nMax1;
      return (0);
    }
  else
    return (1);
}

//
// Implement spatial distortion interpolation tables designed by
//   Marty Stanton, Brandeis University
//

int Cspatial::nStantonPixeltoMM(const float fPx0, 
                                const float fPx1,
                                float* pfXmm, 
                                float* pfYmm)
{
  float f0tmp, f1tmp, f0look, f1look, f0delta, f1delta;
  int   n0look, n1look;
  int   n0Interp[4], n1Interp[4];
  float f0Interp[4], f1Interp[4];
  int   nStat;

  // Set results to something in case we return with error before getting
  // valid results.

  *pfXmm = 0.0;
  *pfYmm = 0.0;

  // *In the Cimage class, elements are numbered
  //   from 0, not 1!





  f0look = ((fPx0+1.0f) - (float) m_n0Start) / (float) m_n0Step; //* + 1.0;
  f1look = ((fPx1+1.0f) - (float) m_n1Start) / (float) m_n1Step; //*  + 1.0;

  n0look = (int) f0look;
  n1look = (int) f1look;

  f0delta = f0look - (float) n0look;
  f1delta = f1look - (float) n1look;

  // Check that the interpolation points are all in bounds and not flagged bad

  nStat = m_poX_int->nGetPixel( n0look, n1look, &n0Interp[0]);
  if ( (0 != nStat) || (n0Interp[0] == m_nBadFlag) ) return (1);

  nStat = m_poX_int->nGetPixel( n0look+1, n1look, &n0Interp[1]);
  if ( (0 != nStat) || (n0Interp[1] == m_nBadFlag) ) return (1);

  nStat = m_poX_int->nGetPixel( n0look, n1look+1, &n0Interp[2]);
  if ( (0 != nStat) || (n0Interp[2] == m_nBadFlag) ) return (1);

  nStat = m_poX_int->nGetPixel( n0look+1, n1look+1, &n0Interp[3]);
  if ( (0 != nStat) || (n0Interp[3] == m_nBadFlag) ) return (1);

  nStat = m_poY_int->nGetPixel( n0look, n1look, &n1Interp[0]);
  if ( (0 != nStat) || (n1Interp[0] == m_nBadFlag) ) return (1);

  nStat = m_poY_int->nGetPixel( n0look+1, n1look, &n1Interp[1]);
  if ( (0 != nStat) || (n1Interp[1] == m_nBadFlag) ) return (1);

  nStat = m_poY_int->nGetPixel( n0look, n1look+1, &n1Interp[2]);
  if ( (0 != nStat) || (n1Interp[2] == m_nBadFlag) ) return (1);

  nStat = m_poY_int->nGetPixel( n0look+1, n1look+1, &n1Interp[3]);
  if ( (0 != nStat) || (n1Interp[3] == m_nBadFlag) ) return (1);

  // Now we have the coordinates of the 4 corners to interpolate from,
  //   so convert to float and do the interpolation

  for (int i = 0; i < 4; i++)
    {
      f0Interp[i] = (float) n0Interp[i];
      f1Interp[i] = (float) n1Interp[i];
    }

  f0tmp = f0Interp[0] +
            f0delta * (f0Interp[1] - f0Interp[0])  +
	    f1delta * (f0Interp[2] - f0Interp[0])  +
	    f0delta * f1delta * (f0Interp[0] + f0Interp[3] -
                                 f0Interp[1] - f0Interp[2]);

  f1tmp = f1Interp[0] +
            f0delta * (f1Interp[1] - f1Interp[0])  +
	    f1delta * (f1Interp[2] - f1Interp[0])  +
	    f0delta * f1delta * (f1Interp[0] + f1Interp[3] -
                                 f1Interp[1] - f1Interp[2]);


  // Apply fDirVec here

  // Correct for m_a2x2fDirVec next

  *pfXmm = m_a2x2fDirVec[0][0] * f0tmp  +  m_a2x2fDirVec[1][0] * f1tmp;
  *pfYmm = m_a2x2fDirVec[0][1] * f0tmp  +  m_a2x2fDirVec[1][1] * f1tmp;

  //  *pfXmm = f0tmp;
  //  *pfYmm = f1tmp;

  // Then scale to proper mm with appropriate beam center offset

  *pfXmm = *pfXmm * m_fPxScale * m_fPxSize - m_fXCenter;
  *pfYmm = *pfYmm * m_fPxScale * m_fPxSize - m_fYCenter;

  if (2 < m_nVerbose)
    {
      cout <<   "Px in:  " << fPx0 << ", " << fPx1
           << "\nInt vl: " << f0tmp << ", " << f1tmp
           << "\nMM out: " << *pfXmm << ", " << *pfYmm << endl;
    }
  return (0);
}

int Cspatial::nStantonMMtoPixel(const float fXmm, const float fYmm,
                      float *pfPx0, float *pfPx1)
{
  float f0tmp, f1tmp, f0look, f1look, f0delta, f1delta;
  int   n0look, n1look;
  int   n0Interp[4], n1Interp[4];
  float f0Interp[4], f1Interp[4];
  int   nStat;

  // Set results to something in case we return with error before getting
  // valid results.

  *pfPx0 = 0.0;
  *pfPx1 = 0.0;

  f0tmp = fXmm + m_fXCenter;
  f1tmp = fYmm + m_fYCenter;

  // Convert mm to mm along pixel axes (is there a translation involved?)

  f0look = m_a2x2fInvVec[0][0] * f0tmp  +  m_a2x2fInvVec[1][0] * f1tmp;
  f1look = m_a2x2fInvVec[0][1] * f0tmp  +  m_a2x2fInvVec[1][1] * f1tmp;

  // Apply fDirVec (not implemented yet)

//  f0look = f0tmp;
//  f1look = f1tmp;

  //

  f0look = f0look / m_fPxSize;
  f1look = f1look / m_fPxSize;

  // The interpolation table is only step resolution

  // *In the Cimage class, elements are numbered
  //   from 0, not 1!
  f0look = ((f0look - (float) m_n0InvStart) / (float) m_n0InvStep); //* + 1.0;
  f1look = ((f1look - (float) m_n1InvStart) / (float) m_n1InvStep); //* + 1.0;

  n0look = (int) f0look;
  n1look = (int) f1look;

  f0delta = f0look - (float) n0look;
  f1delta = f1look - (float) n1look;

  // Check that the interpolation points are all in bounds and not flagged bad

  nStat = m_poInv_x_int->nGetPixel( n0look, n1look, &n0Interp[0]);
  if ( (0 != nStat) || (n0Interp[0] == m_nBadFlag) ) return (1);

  nStat = m_poInv_x_int->nGetPixel( n0look+1, n1look, &n0Interp[1]);
  if ( (0 != nStat) || (n0Interp[1] == m_nBadFlag) ) return (1);

  nStat = m_poInv_x_int->nGetPixel( n0look, n1look+1, &n0Interp[2]);
  if ( (0 != nStat) || (n0Interp[2] == m_nBadFlag) ) return (1);

  nStat = m_poInv_x_int->nGetPixel( n0look+1, n1look+1, &n0Interp[3]);
  if ( (0 != nStat) || (n0Interp[3] == m_nBadFlag) ) return (1);

  nStat = m_poInv_y_int->nGetPixel( n0look, n1look, &n1Interp[0]);
  if ( (0 != nStat) || (n1Interp[0] == m_nBadFlag) ) return (1);

  nStat = m_poInv_y_int->nGetPixel( n0look+1, n1look, &n1Interp[1]);
  if ( (0 != nStat) || (n1Interp[1] == m_nBadFlag) ) return (1);

  nStat = m_poInv_y_int->nGetPixel( n0look, n1look+1, &n1Interp[2]);
  if ( (0 != nStat) || (n1Interp[2] == m_nBadFlag) ) return (1);

  nStat = m_poInv_y_int->nGetPixel( n0look+1, n1look+1, &n1Interp[3]);
  if ( (0 != nStat) || (n1Interp[3] == m_nBadFlag) ) return (1);

  // Now we have the coordinates of the 4 corners to interpolate from,
  //   so convert to float and do the interpolation

  for (int i = 0; i < 4; i++)
    {
      f0Interp[i] = (float) n0Interp[i];
      f1Interp[i] = (float) n1Interp[i];
    }

  f0tmp = f0Interp[0] +
            f0delta * (f0Interp[1] - f0Interp[0])  +
	    f1delta * (f0Interp[2] - f0Interp[0])  +
	    f0delta * f1delta * (f0Interp[0] + f0Interp[3] -
                                 f0Interp[1] - f0Interp[2]);

  f1tmp = f1Interp[0] +
            f0delta * (f1Interp[1] - f1Interp[0])  +
	    f1delta * (f1Interp[2] - f1Interp[0])  +
	    f0delta * f1delta * (f1Interp[0] + f1Interp[3] -
                                 f1Interp[1] - f1Interp[2]);


  *pfPx0 = (f0tmp * m_fPxScale) - 1.0f;  // Pixels are numbered from 0 on ...
  *pfPx1 = (f1tmp * m_fPxScale) - 1.0f;

  // Check that fPx0 and fPx1 are within with image boundaries
  if (   (*pfPx0 < 0.0)
      || (*pfPx1 < 0.0)
      || (*pfPx0 >= (float) m_n0ImageSize)
      || (*pfPx1 >= (float) m_n1ImageSize) )
    return (2);

  return (0);
}


const int       c_nSpatialBeamPosHeaderDecPrecision = 2;
const int       c_nSpatialSimpleInfoHeaderDecPrecision = 5;

int Cspatial::nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre)
{
  // Update the spatial distortion info in an image header
  //
  // Cimage_header *poHeader      The image header to update
  // Returns   0   Success
  //           1   Failure
  //

  Cstring sPrefix;

  sPrefix = sPre;
  if ("" == sPrefix)
    {
      // Null prefix string, see if the  Cdetector::ms_sDetectorNames keyword
      // exists in the header.  If so, change the prefix to the first name
      // in the list of names.  This call is OK, since if there is no
      // keyword, then sPrefix is set to "" anyways.

      (void) poHeader->nGetValue(Cdetector::ms_sDetectorNames, 1, &sPrefix);
    }
  
  poHeader->nReplaceValue(sPrefix + ms_sSpatialBeamPosn, 2, m_a2fBeamPx, c_nSpatialBeamPosHeaderDecPrecision);
  
  if (eSpatial_simple_type == m_eThe_Type)
    {
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionType,
			      ms_sSpatialTypeSimple);
      float fTemp[4];
      fTemp[0] = m_a2fBeamPx[0];
      fTemp[1] = m_a2fBeamPx[1];
      fTemp[2] = m_a2fNominalPxSize[0];
      fTemp[3] = m_a2fNominalPxSize[1];
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionInfo, 4, fTemp, c_nSpatialSimpleInfoHeaderDecPrecision);
    }
  else if (eSpatial_interp_table == m_eThe_Type)
    {
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionType,
			      ms_sSpatialTypeInterp);
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionInfo,
			      m_sFileset_name);
    }
  else if (eSpatial_complex_type == m_eThe_Type)
    {
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionType,
			      ms_sSpatialTypeComplex);
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionInfo,
			      "Undefined");
    }
  else
    {
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionType,
			      ms_sSpatialTypeUnknown);
      poHeader->nReplaceValue(sPrefix + ms_sSpatialDistortionInfo, "None");
    }

  return (0);
}

void
Cspatial::vInvert(void)
{
  // Invert the 2x2 matrix for mm direction transformation
  float fTemp;
  fTemp =  m_a2x2fDirVec[0][0] * m_a2x2fDirVec[1][1]
         - m_a2x2fDirVec[0][1] * m_a2x2fDirVec[1][0];
  if ( (1.0 != fTemp) && (-1.0 != fTemp) )
    {
      cout << "Cspatial WARNING: problem with spatial distortion vectors!\n";
    }
  if (0.0 != fTemp)
    {
      m_a2x2fInvVec[0][0] =  m_a2x2fDirVec[1][1] / fTemp;
      m_a2x2fInvVec[0][1] = -m_a2x2fDirVec[0][1] / fTemp;
      m_a2x2fInvVec[1][0] = -m_a2x2fDirVec[1][0] / fTemp;
      m_a2x2fInvVec[1][1] =  m_a2x2fDirVec[0][0] / fTemp;
    }
  else
    {
      cout << "Cspatial ERROR: problem with spatial distortion vectors!\n";
    }
}

int
Cspatial::nGetBeamPosition(float *pfBeamPx0, float *pfBeamPx1)
{
  *pfBeamPx0 = m_a2fBeamPx[0];
  *pfBeamPx1 = m_a2fBeamPx[1];
  return (0);
}


int Cspatial::nDiff(Cspatial& oSpatialAdd, Cspatial& oSpatialSubtract)
{
    float fBeamPx0Add,fBeamPx1Add;
    float fBeamPx0Sub,fBeamPx1Sub;
    float fBeamPx0,fBeamPx1;

    oSpatialAdd.nGetBeamPosition(&fBeamPx0Add,&fBeamPx1Add);
    oSpatialSubtract.nGetBeamPosition(&fBeamPx0Sub,&fBeamPx1Sub);
    nGetBeamPosition(&fBeamPx0,&fBeamPx1);
    fBeamPx0 += fBeamPx0Add - fBeamPx0Sub;
    fBeamPx1 += fBeamPx1Add - fBeamPx1Sub;
    nSetBeamPosition(fBeamPx0,fBeamPx1);
    return 0;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Change the state of the spatial object, if the pixel size is increased by a ratio.
// RB 8/17/05 This function still needs to be tested. 
bool Cspatial::bChangeBinning(int nFastRatio, int nSlowRatio)
{
    if( nFastRatio <= 0 || nSlowRatio <= 0 )
        return false;
    
    if( eSpatial_simple_type == m_eThe_Type )
    {
        // Change the beam position in pixels
        m_a2fBeamPx[0] /= nFastRatio; 
        m_a2fBeamPx[1] /= nSlowRatio; 
    
        // Change the pixel size in mm
        m_a2fNominalPxSize[0] *= nFastRatio;
        m_a2fNominalPxSize[1] *= nSlowRatio;
    }
    else if(  eSpatial_interp_table == m_eThe_Type )
    {
        if( nFastRatio != nSlowRatio )
            return false;   // since there is just one scale factor in the current implementation
                            // of interp tables, i.e. m_fPxScale, we cannot support two different ratios
            
        m_a2fBeamPx[0] /= nFastRatio; 
        m_a2fBeamPx[1] /= nSlowRatio; 
        
        m_fPxScale /= nFastRatio;
        m_fPxSize *= nFastRatio;
    }
    else
        return false; // unsupported 

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


