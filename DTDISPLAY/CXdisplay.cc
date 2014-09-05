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
// CXdisplay.cc            Initial author: J.W. Pflugrath           18-Jul-1995
//    This file contains the member functions of class CXdisplay.
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

#include "CXdisplay.h"              // Class definition and prototypes

using namespace std;

//+Code begin


//+Definitions, constants, and initialization of static member variables

//+Public functions

// Constructors, destructors and assignments

CXdisplay::CXdisplay(char *pcCimageName, Widget hParentIn, Cimage *poImageIn,
		     CXcolormap *poColormapIn,
		     tagImagePropsDisplay *ptPropsIn,
		     tagReflnPropsDisplay *ptReflnPropsIn) :
UIComponent(pcCimageName)
{
  // Construct a CXdisplay object

  int nStat;
  m_pucData       = NULL;
  m_nDataSize     = 0;
  m_pucScale      = NULL;
  m_nScaleSize    = 0;
  m_hPixmap       = 0;
  m_phXImage      = NULL;
  m_phXImageScale = NULL;
  m_hGC           = NULL;
  m_hGCtext       = NULL;
  m_hGCrefln      = NULL;
  m_hParent       = hParentIn;
  m_poImage       = poImageIn;
  m_pXFont        = NULL;
  m_poDetector    = NULL;

  m_pObj          = NULL;
  m_prvRubberbandCallback = NULL;

  String pcClass, pcProgname;
  XtGetApplicationNameAndClass(XtDisplay(m_hParent), &pcProgname, &pcClass);

  m_bHasFloatValues = FALSE;

  m_bLineInRefln   = TRUE;
  if ("" != sGetEnv("DTDISPLAY_NO_REFLN_LINE"))
    m_bLineInRefln = FALSE;

  m_bPixelSwapBytes = FALSE;
  if (MSBFirst == XImageByteOrder(XtDisplay(m_hParent)))
    m_bPixelSwapBytes =  (eCPU_little_endian == nGetByteOrderCPU());
  else
    m_bPixelSwapBytes =  (eCPU_big_endian == nGetByteOrderCPU());

  if (NULL == ptReflnPropsIn)
    {
      String pcColor;

      m_tReflnProps.a6fReflnSize[0] = 25.0;
      m_tReflnProps.a6fReflnSize[1] = 25.0;
      m_tReflnProps.a6fReflnSize[2] = 0.0;
      m_tReflnProps.a6fReflnSize[3] = 0.0;
      m_tReflnProps.a6fReflnSize[4] = 0.0;
      m_tReflnProps.a6fReflnSize[5] = 0.0;
      m_tReflnProps.sReflnColor     = "red";
      m_tReflnProps.nReflnPredObs   = 0;
      m_tReflnProps.nReflnSymbol    = 4;
      m_tReflnProps.fImageRotStart  = 0.0;
      m_tReflnProps.fImageRotEnd    = 0.0;
      pcColor = XGetDefault (XtDisplay(m_hParent), pcProgname, "reflnPlotColor");
      if (NULL != pcColor)
	{
	  // user has chosen a default color, use that instead
	  m_tReflnProps.sReflnColor = pcColor;
	}
    }
  else
    {
      (void) nSetReflnProps(ptReflnPropsIn, FALSE);
    }

  if (NULL == ptPropsIn)
    {
      m_tImageProps.fZoom             = 0.9;    // Allow room for profiles
      m_tImageProps.nOrient           = -1;      // -1 means not set
      m_tImageProps.fSatPixValue      = 65535.0;
      m_tImageProps.fTooLowPixValue   = -2.0;
      m_tImageProps.sSatPixColor      = "RED";   // UPPERCASE is important here
      m_tImageProps.sTooLowPixColor   = "YELLOW"; // UPPERCASE is important here
      m_tImageProps.nResoCircles      = -1;
      m_tImageProps.nPlotRockingCurve = 0;
      m_tImageProps.fScaleMinSD       = 1.5;
      m_tImageProps.fScaleMaxSD       = 5.0;
      m_tImageProps.fSD               = -1.0;
      m_tImageProps.fImageMin         = 0.0;
      m_tImageProps.fImageMax         = 65535.0;
      m_tImageProps.fAvg              = 100.0;
      m_tImageProps.fScaleMin         = 0.0;
      m_tImageProps.fScaleMax         = 65535.0;
    }
  else
    {
      m_tImageProps.fZoom             = ptPropsIn->fZoom;
      m_tImageProps.nOrient           = ptPropsIn->nOrient;
      m_tImageProps.nResoCircles      = ptPropsIn->nResoCircles;
      m_tImageProps.nPlotRockingCurve = ptPropsIn->nPlotRockingCurve;
      m_tImageProps.fSatPixValue      = ptPropsIn->fSatPixValue;
      m_tImageProps.fTooLowPixValue   = ptPropsIn->fTooLowPixValue;
      m_tImageProps.sSatPixColor      = ptPropsIn->sSatPixColor;
      m_tImageProps.sTooLowPixColor   = ptPropsIn->sTooLowPixColor;
      m_tImageProps.nIntegrateBox[0]  = ptPropsIn->nIntegrateBox[0];
      m_tImageProps.nIntegrateBox[1]  = ptPropsIn->nIntegrateBox[1];
      m_tImageProps.fSD               = ptPropsIn->fSD;
      m_tImageProps.fAvg              = ptPropsIn->fAvg;
      m_tImageProps.fImageMin         = ptPropsIn->fImageMin;
      m_tImageProps.fImageMax         = ptPropsIn->fImageMax;
      m_tImageProps.fScaleMin         = ptPropsIn->fScaleMin;
      m_tImageProps.fScaleMax         = ptPropsIn->fScaleMax;
    }
  
  // The following should be initialized from some X resource, not like this

  m_a10sReflnColors[0] = "blue";
  m_a10sReflnColors[1] = "red";
  m_a10sReflnColors[2] = "green";
  m_a10sReflnColors[3] = "magenta";
  m_a10sReflnColors[4] = "red";
  m_a10sReflnColors[5] = "yellow";
  m_a10sReflnColors[6] = "yellow";
  m_a10sReflnColors[7] = "yellow";
  m_a10sReflnColors[8] = "yellow";
  m_a10sReflnColors[9] = "yellow";

  m_tImageProps.fScaleMinSD       = 1.5;
  m_tImageProps.fScaleMaxSD       = 5.0;
  m_tImageProps.fAspectRatio      = 1.0;
  m_tImageProps.nAspectMode       = 1;       // Default is true x/y perspective

  m_tImageProps.nOrig[0]          = 0;
  m_tImageProps.nOrig[1]          = 0;
  m_tImageProps.nDisplayReflns    = 1;
  m_tImageProps.nDisplayPixValues = 1;
  m_tImageProps.nAutoRescale      = 0;

  m_fScaleMin     = m_tImageProps.fScaleMin;
  m_fScaleMax     = m_tImageProps.fScaleMax;
  m_nScaleMode    = 0;
  m_nXPrev        = 0;
  m_nYPrev        = 0;
  m_fPx0Prev      = 0;
  m_fPx1Prev      = 0;
  m_nDir0         = 1;
  m_nDir1         = 1;
  m_poRubberband  = NULL;
  m_bSpacebarHeldDown = FALSE;
  m_bControlKeyHeldDown = FALSE;

  m_poReflnlist   = NULL;
  m_poColormap    = poColormapIn;
  m_poCursor      = new CXcursor(m_hParent);
  m_nCursorMode   = 0;
  m_nBeamstopShadowMode = 0;
  m_tSpotInfo.nStat = -1;  // Indicates spot info is not available

  int a4nPxArea[4];

  // Compute initial scaling factors for converting Cimage data to X11 Pixmap.
  // Base scaling factors on center 80% of Cimage pixels.

  if (NULL == m_poImage) return;           // Failure
  if (!m_poImage->bIsAvailable()) return;  // Failure

  m_bHasFloatValues = (eImage_realIEEE == m_poImage->m_oHeader.eGetImageDataType());

  m_poImage->nGetDimensions(&m_tImageProps.nDim[0], &m_tImageProps.nDim[1]);
  a4nPxArea[0] = 1;
  a4nPxArea[1] = 1;
  m_poImage->m_oHeader.nGetValue(D_K_DtdisplayTile, 2, a4nPxArea);
  m_tImageProps.nPlotRockingCurve = a4nPxArea[0] * a4nPxArea[1]; 

  m_tImageProps.nExt[0]    = m_tImageProps.nDim[0];
  m_tImageProps.nExt[1]    = m_tImageProps.nDim[1];
  a4nPxArea[0] = m_tImageProps.nExt[0] / 10;
  a4nPxArea[1] = m_tImageProps.nExt[1] / 10;
  a4nPxArea[2] = m_tImageProps.nExt[0] - a4nPxArea[0];
  a4nPxArea[3] = m_tImageProps.nExt[1] - a4nPxArea[1];

  if (0.0 > m_tImageProps.fSD)
    {
      // Compute Average and Standard deviation, only when input SD is
      // less than 0

//      cout << "...computing average, sd, min, max ..." << endl;
      float fMin = 0.0;

      if (m_poImage->nGetDataType() == (int) eImage_realIEEE)
	{
	  fMin = -65536.0;
	}
      else if (m_poImage->nGetDataType() == (int) eImage_I2)
	fMin = -32768.0;
      else if (m_poImage->nGetDataType() == (int) eImage_I4)
	fMin = -65536.0;

      nStat = nAvgSD(*m_poImage, a4nPxArea, 
		     fMin, min(65534.0, m_poImage->fGetSatValue() - 1.0),
		     &m_tImageProps.fImageMin, &m_tImageProps.fImageMax,
		     &m_tImageProps.fAvg, &m_tImageProps.fSD);
    }

  // Add the callbacks to the parent widget

    _clientDataStructs = new UICallbackStruct;
    _clientDataStructs->object = this;
    _clientDataStructs->client_data = (XtPointer)0;
    XtAddCallback(m_hParent,
        XmNexposeCallback,
        CXdisplay::vEventCallbackCallback,
        (XtPointer)_clientDataStructs);
    XtAddCallback(m_hParent,
        XmNinputCallback,
        CXdisplay::vEventCallbackCallback,
        (XtPointer)_clientDataStructs);
    XtAddCallback(m_hParent,
        XmNresizeCallback,
        CXdisplay::vEventCallbackCallback,
        (XtPointer)_clientDataStructs);

  // Get depth of screen for use in depth of images

  m_ulDepthImage = DefaultDepthOfScreen(XtScreen(m_hParent));

  // Get a graphic context.  Replace this later with a Cgc class.

  m_hGC     = XCreateGC(XtDisplay(m_hParent),
			XtWindow(m_hParent), 0, NULL);

  // Get default font in the GC

  m_pXFont = XQueryFont(XtDisplay(m_hParent), XGContextFromGC(m_hGC));

  // Get a graphics context for the drawing of text.
  // Get the text color from the resource database or set to blue if blue
  // exists in the colormap or set to foreground pixel.

  m_ulColorBlue  = m_poColormap->ulGetColorIndex("blue");
  m_ulColorWhite = m_poColormap->ulGetColorIndex("white");

  String pcColor;
  pcColor = XGetDefault (XtDisplay(m_hParent), pcProgname, "pixelValueColor");
  if (pcColor != NULL)
    {
      m_hGCValues.foreground = m_poColormap->ulGetColorIndex(pcColor);	  
    }
  else
    {
      m_hGCValues.foreground = m_ulColorBlue;
    }
  if (m_hGCValues.foreground == 65535)
    {
      m_hGCValues.foreground = BlackPixel(XtDisplay(m_hParent),
					  DefaultScreen(XtDisplay(m_hParent)));
    }

  m_hGCValues.background = m_ulColorWhite;

//best        m_hGCValues.function = GXorReverse;
//bad         m_hGCValues.function = GXorInverted;
//bad         m_hGCValues.function = GXandReverse;
//badwblack   m_hGCValues.function = GXnand;
//secondbest  m_hGCValues.function = GXequiv;
//nope        m_hGCValues.function = GXnor;
//worst       m_hGCValues.function = GXandInverted;
//ok          m_hGCValues.function = GXinvert;
//allblue      m_hGCValues.function = GXcopy;
//mostlygrey      m_hGCValues.function = GXxor;

  m_hGCValues.function = GXcopy;
  m_hGCtext = XCreateGC(XtDisplay(m_hParent), XtWindow(m_hParent),
			GCFunction | GCForeground, &m_hGCValues);

  // Get reflection list color, default is red

  /*pcColor = XGetDefault (XtDisplay(m_hParent), pcProgname, "reflnPlotColor");

  if (NULL != pcColor)
    {
	  m_hGCValues.foreground = m_poColormap->ulGetColorIndex(pcColor);	
	  //m_hGCValues.foreground = m_poColormap->ulGetColorIndex("red");
    }
  else
    {
      m_hGCValues.foreground = m_poColormap->ulGetColorIndex("red");
    }*/
  if (m_hGCValues.foreground == 65535)
    {
      m_hGCValues.foreground = BlackPixel(XtDisplay(m_hParent),
					DefaultScreen(XtDisplay(m_hParent)));
    }

  m_hGCrefln = XCreateGC(XtDisplay(m_hParent), XtWindow(m_hParent),
		      GCFunction | GCForeground, &m_hGCValues);

  //

  m_a5fQuad[0][0] = -1.0f;
  m_a5fQuad[0][1] = -1.0f;
  m_a5fQuad[1][0] = -1.0f;
  m_a5fQuad[1][1] = -1.0f;
  m_a5fQuad[2][0] = -1.0f;
  m_a5fQuad[2][1] = -1.0f;
  m_a5fQuad[3][0] = -1.0f;
  m_a5fQuad[3][1] = -1.0f;
  m_a5fQuad[4][0] = -1.0f;
  m_a5fQuad[4][1] = -1.0f;

  m_fLineWidth = 0.0;
  //

  m_nMin =   0;
  m_nMax = 127;
  m_nMax = min(127, m_poColormap->nNumColors()-1);
  //cout << "m_nMax is " << m_nMax << endl;
  int nColors = m_nMax - m_nMin + 1;

  //  cout << "INFO color for Saturated pixel value color: " 
  //           << m_tImageProps.sSatPixColor << "\n" << flush;
  m_ulSatColorIndex = m_poColormap->ulGetColorIndex(
				m_tImageProps.sSatPixColor.string());
  if (m_ulSatColorIndex == CXcolormap::ms_ulCOLORERROR) m_ulSatColorIndex = 1;
  m_ulTooLowColorIndex = m_poColormap->ulGetColorIndex(
				m_tImageProps.sTooLowPixColor.string());
  if (m_ulTooLowColorIndex == CXcolormap::ms_ulCOLORERROR)
    {
      //      cout << "ERROR bad color for TooLow pixel value color: " 
      //           << m_tImageProps.sTooLowPixColor << "\n" << flush;
      m_ulTooLowColorIndex = 1;
    }
  if (m_poColormap->nNumColors() >= nColors) 
    {
      // Convert Cimage to a byte array, scaled between ...

      m_tImageProps.fScaleMin  = max(m_tImageProps.fImageMin,
                                    (m_tImageProps.fAvg 
				     - m_tImageProps.fScaleMinSD 
				       * m_tImageProps.fSD));
      m_tImageProps.fScaleMax  = min(m_tImageProps.fImageMax,
                                 (m_tImageProps.fAvg 
				     + m_tImageProps.fScaleMaxSD 
				       * m_tImageProps.fSD));
      vRefresh();
    }
}

CXdisplay::~CXdisplay() 
{
  // Clear the window so users no that there is no image display anymore

  XClearWindow(XtDisplay(m_hParent), XtWindow(m_hParent));

  // Remove callbacks added to the parent widget

  XtRemoveCallback(m_hParent,
		   XmNexposeCallback,
		   CXdisplay::vEventCallbackCallback,
		   (XtPointer)_clientDataStructs);
  XtRemoveCallback(m_hParent,
		   XmNinputCallback,
		   CXdisplay::vEventCallbackCallback,
		   (XtPointer)_clientDataStructs);
  XtRemoveCallback(m_hParent,
		   XmNresizeCallback,
		   CXdisplay::vEventCallbackCallback,
		   (XtPointer)_clientDataStructs);

  m_pObj                  = NULL;
  m_prvRubberbandCallback = NULL;
  vDeleteXStuff();

  m_poReflnlist = NULL;

  if (NULL != m_poCursor)
    {
      delete m_poCursor;
      m_poCursor = NULL;
    }
  if (NULL != m_hGC)
    XFreeGC(XtDisplay(m_hParent), m_hGC);

  if (NULL != m_hGCtext)
    XFreeGC(XtDisplay(m_hParent), m_hGCtext);

  if (NULL != m_hGCrefln)
    XFreeGC(XtDisplay(m_hParent), m_hGCrefln);

  // Must free font AFTER freeing of GCs!?

//  if (NULL != m_pXFont)
//    XFreeFont(XtDisplay(m_hParent), m_pXFont);
}

//+Public functions

void
CXdisplay::vDeleteXStuff(void) 
{
  if (NULL != m_phXImage)
    {
      m_phXImage->data = NULL;        // Could only point to m_pucData
      XDestroyImage(m_phXImage);
      m_phXImage = NULL;
    }

  if (NULL != m_phXImageScale)
    {
      m_phXImageScale->data = NULL;   // Could only point to m_pucScale 
      XDestroyImage(m_phXImageScale);
      m_phXImageScale = NULL;
    }

  if (0 != m_hPixmap)
    XFreePixmap(XtDisplay(m_hParent), m_hPixmap);

/*
Try to reuse the allocated memory rather than delete/new
  if (NULL != m_pucData) 
    {
      delete [] m_pucData;
      m_pucData = NULL;
      m_nDataSize = 0;
    }

  if (NULL != m_pucScale) 
    {
      delete [] m_pucScale;
      m_pucScale = NULL;
      m_nScaleSize = 0;
    }
*/

}

int
CXdisplay::nAvgSD(Cimage& oImage, const int nPxArea[4], 
		  const float fExcludeMin,
		  const float fExcludeMax,
		  float *pfMin, float *pfMax,
		  float *pfAvg, float *pfSD)
{
  
// Compute background and standard deviation in a portion of an image.
// Exclude from the calculation pixels with values <= fExcludeMin and 
// => fExcludeMax.  Look only in the the image area defined by nPxArea.
// 
  int   i, j;
  int   nStat;
  float fValue;
  double fSum1, fSum2, fNum;

  fSum1 = 0.0;
  fSum2 = 0.0;
  fNum  = 0.0;

  *pfMin = (oImage.*oImage.prfGetPixel)(nPxArea[0], nPxArea[1]);
  *pfMax = (oImage.*oImage.prfGetPixel)(nPxArea[0], nPxArea[1]);

  for (j = nPxArea[1]; j < nPxArea[1]+nPxArea[3]; j++) {
    for (i = nPxArea[0]; i < nPxArea[0]+nPxArea[2]; i++) {
//jwp      if ( bPixelInLimits(i, j) && !poScan->poNonunf->bPixelIsBad(i, j) ) {
//jwp
//jwp	// Pixel is in the limits and is not bad, so get it and apply nonunf:

	fValue = (oImage.*oImage.prfGetPixel)(i, j);
//jwp               * (poScan->poNonunf->*prfNonunfFactor)(i, j);

	if ( (fValue >= fExcludeMin) && (fValue <= fExcludeMax) ) {

	  // This pixel passes all tests to include in the background and sd
	  //    calculation.
	  fSum1 = fSum1 + fValue;
	  fSum2 = fSum2 + (fValue * fValue);
	  fNum  = fNum + 1.0;
	  if (fValue < *pfMin) *pfMin = fValue;
	  if (fValue > *pfMax) *pfMax = fValue;
	}
//jwp      }
    }  // end of i loop
  }  // end of j loop

  // For a valid average and sd to be calculated we must have at least 1
  // pixel in the average and more than 1/2 of the pixels in the specified
  // area must contribute!

  if ( (fNum > 1.0) && (fNum > (nPxArea[2] * nPxArea[3] / 2)) ) {
    *pfAvg = fSum1 / fNum;
    fSum2  = (fSum2 - (fSum1*fSum1 / fNum)) / (fNum - 1.0);
    if (fSum2 >= 0.0) {
      *pfSD = (float)sqrt(fSum2);
      nStat = 0;
    }
    else {
      *pfSD = (float)sqrt(-fSum2);
      cerr << "nAvgSD, sqrt of negative!\n";
      nStat = 1;
    }
  }
  else {
//    cerr << "nAvgSD, not enough pixels to compute!\n";
    nStat = 2;
  }

  return (nStat);
}

void
CXdisplay::vImageToCLUT(Cimage& oImage, const float fMin, const float fMax,
                        const int nIndexMin, const int nIndexMax)
{
  // Convert Cimage data to indices in a color lookup table

  int   i, j;       // Loop counters
  int   nValue;
  float fValue;
  float fScale;
  int   nPx0, nPx1;
  float fPixRatio;

  unsigned char *pucTemp;
  unsigned short int *puiTemp;

//+2011-10-12 JWP
// Later below want to apply some kind of resolution modification of the contrast/color
// so figure out pixel which is closest to the 
// of current beam position (if on detector)

  float fBeamPx0 = -999.0;
  float fBeamPx1 = -999.0; 
  int   nBeamPx0, nBeamPx1;
  float fResoContrastFactor = 1.0;
  float fContrastFactorApplied = 1.0;
  float fMaxDistSquared = 0.0;
  if (NULL != m_poDetector)
    {
      float a3fX[3], a3fXR[3];
      float fMm0, fMm1;
      float a3x3fMat[3][3];
      float fDistSq, fTemp0, fTemp1;
      int nStat1;

      a3fX[0] = 0.0;
      a3fX[1] = 0.0;
      a3fX[2] = 0.0;

      nStat1 = m_poDetector->nCalcDetCoords(a3fX, m_a3fS0,
					    &fMm0, &fMm1, &fBeamPx0, &fBeamPx1);
      if (0 == nStat1)
	{
	  // So if both fBeamPx0 and fBeamPx1 are greater than 0, we have a beam position
	  if (   (0.0 < fBeamPx0) 
              && (0.0 < fBeamPx1)
	      && ("" != sGetEnv("DTDISPLAY_RESO_CONTRAST")) 
		 ) 
	    fResoContrastFactor = atof(sGetEnv("DTDISPLAY_RESO_CONTRAST").string()); 
	  nBeamPx0 = nint(fBeamPx0);
	  nBeamPx1 = nint(fBeamPx1);
	  fTemp0  = fBeamPx0 - ((float)oImage.nGetDimension(0) - 1.0);
	  fTemp1  = fBeamPx1 - ((float)oImage.nGetDimension(1) - 1.0);
	  // Look at "distance squared" to the 4 corners of the image
	  fDistSq = fBeamPx0*fBeamPx0 + fBeamPx1*fBeamPx1;  // (0, 0)
	  fMaxDistSquared = max(fMaxDistSquared, fDistSq);

	  fDistSq = fBeamPx0*fBeamPx0 + fTemp1 * fTemp1; // (0, nDim1)
	  fMaxDistSquared = max(fMaxDistSquared, fDistSq);

	  fDistSq = fTemp0 * fTemp0 +  fBeamPx1 * fBeamPx1; //(nDim0, 0)
	  fMaxDistSquared = max(fMaxDistSquared, fDistSq);

	  fDistSq = fTemp0 * fTemp0 +  fTemp1 * fTemp1; // (nDim0, nDim1)
	  fMaxDistSquared = max(fMaxDistSquared, fDistSq);
	  if (0.0 >= fMaxDistSquared) fMaxDistSquared = 1.0;
	}
    }
//-2011-10-12 JWP

  // We only allocate up to 255 colors for the image regardless of the
  // number of colors supported by the display.  We do this because then we
  // only need 1 byte per pixel.

  if (   (nIndexMin < 0) || (nIndexMin > 255) 
      || (nIndexMax > m_poColormap->nNumColors())
      || (nIndexMin > m_poColormap->nNumColors())
      || (nIndexMax < 0) || (nIndexMax > 255) )
    return;

  m_nBytesPerPixel = 1;
  if (8 < m_ulDepthImage)
    m_nBytesPerPixel++;
  if (16 < m_ulDepthImage)
    m_nBytesPerPixel++;
  if (24 < m_ulDepthImage)
    m_nBytesPerPixel++;

  //cout << "Depth is " << m_ulDepthImage << " bytes per pixel is: " 
  // << m_nBytesPerPixel << endl;
  if ((XImage*)NULL == m_phXImage)
    {
      // You are in big trouble!

      cout << "ERROR NULL image in vImageToCLUT!\n" << flush;
    }

  int nSizeNeeded = m_phXImage->bytes_per_line * m_phXImage->height;
  if ( (NULL != m_pucData) && (m_nDataSize < nSizeNeeded) )
    {
      delete [] m_pucData;
      m_pucData = NULL;
    }
  if (NULL == m_pucData)
    {
      // Allocate the memory

      m_pucData = new unsigned char [nSizeNeeded];
      m_nDataSize = nSizeNeeded;
    }

  m_phXImage->data = (char*) m_pucData;

//  cout << "XImage bytes_per_line, height: " << m_phXImage->bytes_per_line 
//       << "  " << m_phXImage->height << endl;
//  cout << "XImage wu, hu: " << m_nWidthUsed << "  " << m_nHeightUsed << endl;

  unsigned long int ulColor;
  pucTemp   = m_pucData;                      // Set up pointers for looping
  puiTemp   = (unsigned short int*)m_pucData;
  
  // Change looping because *m_puCdata is not full size, but size of drawingarea

  if ((fMax - fMin) != 0.0)
    fScale = (float) (nIndexMax - nIndexMin) / (fMax - fMin);
  else
    {
      fScale = (float) (nIndexMax - nIndexMin);
    }

  // Figure out looping start, end, increment in order to 

  ///////////////////////////////////////////////////////////////////////

//  cout << "Enter Orient, Zoom, Mode and PixelSizeRatio: " << endl;
//  cin >> nOrient >> fZoom >> nMode >> fSizeRatio;

  // Compute Cimage pixels per Pixmap pixel
  // 1. Get ratio of pixmap edges

  fPixRatio = (float) m_nWidthUsed / (float) m_nHeightUsed;

  // m_fStep* is the number of Cimage pixels per Pixmap pixels

  if (4 > m_tImageProps.nOrient) 
    {
      m_fStep0 = (float) m_tImageProps.nExt[0] / (float) m_nWidthUsed;
      m_fStep1 = (float) m_tImageProps.nExt[1] / (float) m_nHeightUsed;
    }
  else
    {
      m_fStep0 = (float) m_tImageProps.nExt[1] / (float) m_nWidthUsed;
      m_fStep1 = (float) m_tImageProps.nExt[0] / (float) m_nHeightUsed;
    }

  if (1 == m_tImageProps.nAspectMode)
    {
      // Find which direction is longer in real (i.e. mm) terms

      if (4 > m_tImageProps.nOrient) 
	m_fStep1 = m_fStep1 * m_tImageProps.fAspectRatio;
      else
	m_fStep0 = m_fStep0 * m_tImageProps.fAspectRatio;

      float fLen0 = m_fStep0 * (float) m_nWidthUsed;  // fLen is length in 
      float fLen1 = m_fStep1 * (float) m_nHeightUsed; // arbitary units

      if (fLen0 > fLen1)
	m_fStep1 = m_fStep0;
      else
	m_fStep0 = m_fStep1;

      if (4 > m_tImageProps.nOrient) 
	m_fStep1 = m_fStep1 * m_tImageProps.fAspectRatio; // Be sure to re-apply this
      else
	m_fStep0 = m_fStep0 * m_tImageProps.fAspectRatio;
    }
  
  // Recalculate and set nExt

  if (4 > m_tImageProps.nOrient) 
    {
      m_tImageProps.nExt[0]  = (int) (m_fStep0 * (float) m_nWidthUsed   + 0.5);
      m_tImageProps.nExt[1]  = (int) (m_fStep1 * (float) m_nHeightUsed  + 0.5);
    }
  else
    {
      m_tImageProps.nExt[1]  = (int) (m_fStep0 * (float) m_nWidthUsed   + 0.5);
      m_tImageProps.nExt[0]  = (int) (m_fStep1 * (float) m_nHeightUsed  + 0.5);
    }

  if (m_tImageProps.nExt[0] > m_tImageProps.nDim[0])
    m_tImageProps.nExt[0] = m_tImageProps.nDim[0];
  if (m_tImageProps.nExt[1] > m_tImageProps.nDim[1])
    m_tImageProps.nExt[1] = m_tImageProps.nDim[1];

  ///////////////////////////////////////////////////////////////////////

  if (0 >= m_tImageProps.nOrient)
    {
      // 0 -> 0
      // 1 -> 1
      m_nStart0 = m_tImageProps.nOrig[0];
      m_nStart1 = m_tImageProps.nOrig[1];
    }
  else if (1 == m_tImageProps.nOrient)
    {
      // 0 -> 0
      // 1 -> -1
      m_nStart0 = m_tImageProps.nOrig[0];
      m_nStart1 = m_tImageProps.nOrig[1] + m_tImageProps.nExt[1] - 1;
      m_fStep1  = -m_fStep1;
    }
  else if (2 == m_tImageProps.nOrient)
    {
      // 0 -> -0
      // 1 -> 1
      m_nStart0 = m_tImageProps.nOrig[0] + m_tImageProps.nExt[0] - 1;
      m_nStart1 = m_tImageProps.nOrig[1];
      m_fStep0  = -m_fStep0;
    }
  else if (3 == m_tImageProps.nOrient)
    {
      // 0 -> -0
      // 1 -> -1
      m_nStart0 = m_tImageProps.nOrig[0] + m_tImageProps.nExt[0] - 1;
      m_nStart1 = m_tImageProps.nOrig[1] + m_tImageProps.nExt[1] - 1;
      m_fStep0  = -m_fStep0;
      m_fStep1  = -m_fStep1;
    }
  else if (4 == m_tImageProps.nOrient)
    {
      // 0 -> 1
      // 1 -> 0
      m_nStart0 = m_tImageProps.nOrig[1];
      m_nStart1 = m_tImageProps.nOrig[0];
    }
  else if (5 == m_tImageProps.nOrient)
    {
      // 0 -> -1
      // 1 -> 0
      m_nStart0 = m_tImageProps.nOrig[1];
      m_nStart1 = m_tImageProps.nOrig[0] + m_tImageProps.nExt[0] - 1;
      m_fStep1  = -m_fStep1;
    }
  else if (6 == m_tImageProps.nOrient)
    {
      // 0 -> 1
      // 1 -> -0
      m_nStart0 = m_tImageProps.nOrig[1] + m_tImageProps.nExt[1] - 1;
      m_nStart1 = m_tImageProps.nOrig[0];
      m_fStep0  = -m_fStep0;
    }
  else if (7 == m_tImageProps.nOrient)
    {
      // 0 -> -1
      // 1 -> -0
      m_nStart0 = m_tImageProps.nOrig[1] + m_tImageProps.nExt[1] - 1;
      m_nStart1 = m_tImageProps.nOrig[0] + m_tImageProps.nExt[0] - 1;
      m_fStep0  = -m_fStep0;
      m_fStep1  = -m_fStep1;
    }

  m_nDir0   = 1;
  m_nDir1   = 1;
  if (0.0 > m_fStep0) m_nDir0 = -1;
  if (0.0 > m_fStep1) m_nDir1 = -1;

//+2010-10-25 jwp, some tests of rendering the "in-between" pixels
// which may be more useful with something like a Pilatus.

  float fMaxPix, fTestVal;
  bool bGetMaxNeighbor = FALSE;
  bool bTestNeighbor0, bTestNeighbor1;

  // If the step size causes one to skip over more than 3 pixels,
  // then try to display the max image pixel that is within 2 pixels
  // of the desired image pixel.  We hope this will prevent small (in area) spots from
  // not being shown (that is, they will be shown) and not skipped by the mapping from
  // display pixels back to image pixels.
  //
  // Test limits once outside of the inner loop
  // ?? 3.5, 3.5, 2, 3  // For big step sizes
  // ?? 2.5, 2.5, 1, 2  // For med step sizes
  // Do not worry about small step sizes
  // ?? 1.5  1.5  0  1
 
  float fStepLimit    = 1.5;
  int   nLoLimit      = 1;
  int   nHiLimitDelta = 0;

  if (  (2.5 < ABS(m_fStep0)) || (2.5 < ABS(m_fStep1)) )
    {
      nLoLimit      = 2;
      nHiLimitDelta = 1;
      fStepLimit    = 2.5;
    }
  if (  (3.5 < ABS(m_fStep0)) || (3.5 < ABS(m_fStep1)) )
    {
      nLoLimit      = 3;
      nHiLimitDelta = 2;
      fStepLimit    = 3.5;
    }
//+2012-09-26
  if ("" != sGetEnv("DTREK_DTDISPLAY_TOOLOW"))
    {
      m_tImageProps.fTooLowPixValue = atof(sGetEnv("DTREK_DTDISPLAY_TOOLOW").string());
    }
//-2012-09-26
      
  bGetMaxNeighbor = (fStepLimit < ABS(m_fStep0)) || (fStepLimit < ABS(m_fStep1));
  int nn0, nn1;

  for (j = 0; j < m_nHeightUsed; j++) 
    {
      nPx1 = m_nStart1 + (int) ((float) j * m_fStep1);
///////
      bTestNeighbor1 = TRUE;
      // Preliminary boundary checks
      if (nPx1 <= nLoLimit) bTestNeighbor1 = FALSE;
      if ((nPx1 + nHiLimitDelta) > (m_tImageProps.nDim[1]-1)) bTestNeighbor1 = FALSE;
      bTestNeighbor1 = bTestNeighbor1 && bGetMaxNeighbor; // This should be fast
/////////

      for (i = 0; i < m_nWidthUsed; i++) 
	{
	  nPx0 = m_nStart0 + (int) ((float) i * m_fStep0);
/////////
	  bTestNeighbor0 = TRUE;
	  if (bTestNeighbor1)
	    {
	      if (nPx0 <= nLoLimit) bTestNeighbor0 = FALSE;
	      if ((nPx0 + nHiLimitDelta) > (m_tImageProps.nDim[0]-1)) bTestNeighbor0 = FALSE;
	    }
	  bTestNeighbor0 = bTestNeighbor0 && bTestNeighbor1; // This should be fast
/////////
  	  if (4 > m_tImageProps.nOrient)
	    {
	      fValue = (oImage.*oImage.prfGetPixel)(nPx0, nPx1);
/////////
	      if (bTestNeighbor0)
		{
		  // Find maximum nearest neighbor to display

		  // cout << "TN: " << m_fStep0 << ", " << m_fStep1 << ": " << nPx0 << ", " << nPx1 << endl;
		  for (nn0 = -nHiLimitDelta; nn0 < nLoLimit; nn0++)
		    {
		      for (nn1 = -nHiLimitDelta; nn1 < nLoLimit; nn1++)
			{
			  fTestVal = (oImage.*oImage.prfGetPixel)(nPx0+nn0, nPx1+nn1);
			  if (fTestVal > fValue) fValue = fTestVal;
			}
		    }
		}
/////////
	    }
	  else
	    {
	      fValue = (oImage.*oImage.prfGetPixel)(nPx1, nPx0);
/////////
	      if (bTestNeighbor0)
		{
		  // Find maximum nearest neighbor to display 

		  for (nn0 = -nHiLimitDelta; nn0 < nLoLimit; nn0++)
		    {
		      for (nn1 = -nHiLimitDelta; nn1 < nLoLimit; nn1++)
			{
			  fTestVal = (oImage.*oImage.prfGetPixel)(nPx1+nn0, nPx0+nn1);
			  if (fTestVal > fValue) fValue = fTestVal;
			}
		    }
		}
/////////
	    }

	  if (fValue >= m_tImageProps.fSatPixValue)
	    {
	      // The pixel is saturated, so flag it with the saturated color
	      // if one is available.

	      ulColor = m_ulSatColorIndex;
	    }
//+2012-09-26 JWP
	  else if (fValue <= m_tImageProps.fTooLowPixValue)
	    {
	      // The pixel is below a certain value, 
	      // so flag it with the 'too low' color
	      // if one is available.

	      ulColor = m_ulTooLowColorIndex;
	    }
//-2012-09-26 JWP
	  else
	    {

//+2011-10-12 JWP
// We want to apply some kind of resolution modification of fScale here
// If fResoContrastFactor is not 1? 
// nPx0, nPx1 is the image pixel coords (not the display coords)
// nBeamPx0, nBeamPx1 
// Compute "distance"
	      fContrastFactorApplied = 1.0;
	      if (1.0 != fResoContrastFactor)
		{
		  int nDistSq;
		  nDistSq = 1 + (nPx0 - nBeamPx0) * (nPx0 - nBeamPx0) 
                             +  (nPx1 - nBeamPx1) * (nPx1 - nBeamPx1); 
		  fContrastFactorApplied = (float) nDistSq / fMaxDistSquared; // Should go 0 to 1
		  fContrastFactorApplied = 1.0 + (fContrastFactorApplied * fResoContrastFactor);
		}


//-2011-10-12 JWP
	      nValue = (int) ((fValue - fMin) * fScale * fContrastFactorApplied + 0.5);  // Round-off
	      if (nValue < nIndexMin) 
		nValue = nIndexMin;
	      else if (nValue > nIndexMax)
		nValue = nIndexMax;
	      ulColor = m_poColormap->ulColorIndex(nValue);
	    }
	  // DANGER if byte_order of Xserver and Xclient are different!
	  switch (m_nBytesPerPixel)
	    {
	    case 1:
	      *pucTemp++ = (unsigned char)ulColor;
	      break;
	    case 2:
	      *puiTemp++ = (unsigned short int)ulColor;
	      break;
	    default:
	      // Must use relatively slow
	      XPutPixel(m_phXImage, i, j, ulColor);
	      break;
	    }
	}

      // Skip over pad bytes if necessary
      
      switch(m_nBytesPerPixel)
	{
	case 1:
	  for (i = 0; i < (m_phXImage->bytes_per_line - m_nWidthUsed); i++)
	    pucTemp++;
	  break;
	case 2:
	  for (i = 0; i < (m_phXImage->bytes_per_line/sizeof(unsigned short int)
			   - m_nWidthUsed); i++)
	    puiTemp++;
	  break;
	default:
	  break;
	}
    }
}

void CXdisplay::vEventCallback(Widget w,
			       XtPointer client_data, XtPointer call_data)
{
  XmDrawingAreaCallbackStruct *cbs =
    (XmDrawingAreaCallbackStruct *) call_data;

  XEvent *pXevent = cbs->event;

  if (cbs->reason == XmCR_EXPOSE)
    {
      Display *dpy = pXevent->xany.display;
      
      XExposeEvent *pXexpose = (XExposeEvent *)cbs->event;

      XCopyArea (dpy, m_hPixmap, pXevent->xany.window, m_hGC, 
		 pXexpose->x, pXexpose->y, 
		 pXexpose->width, pXexpose->height, 
		 pXexpose->x, pXexpose->y);
#ifdef DEBUG
  XSync(XtDisplay(m_hParent), False);
  XmUpdateDisplay(m_hParent);
#endif

      vPixelInfo();                           // Refresh pixel info

      vDrawColorScale(1);

//      cout << "DA expose:" << pXexpose->x << ", " << pXexpose->y << ", " 
//	   << pXexpose->width << ", " << pXexpose->height << endl;
    }
  else if (cbs->reason == XmCR_INPUT)
    {
      int nX = pXevent->xbutton.x;
      int nY = pXevent->xbutton.y;
      if (pXevent->xany.type == ButtonPress) 
	{
	  //	  cout << "DA button pressed, (x, y):"
	  //	       << nX << ", " << nY << endl;

	  if ( (Button1 == pXevent->xbutton.button) && (RB_NONE_MODE == m_nCursorMode) )
	    {
	      vPixelInfo(nX, nY, 1, 1);
	    }
	  else if ( (Button2 == pXevent->xbutton.button) && (NULL != m_poRubberband) )
	    {
	      // Will do panning on motion_notify
	    }
	  else if ( (  (Button2 == pXevent->xbutton.button)
		   && (0 == (pXevent->xbutton.state
			     & (Button1Mask | Button3Mask)) ) )
		    || (RB_NONE_MODE != m_nCursorMode) )
	    {
	      // Create a new rubberband object. Rectangle, Line or Circle
	      // depending on whether none, ctrl or shift key pressed

	      if ((nX < m_nWidthUsed) && (nY < m_nHeightUsed) )
		{
		  if (   (0 != (pXevent->xbutton.state & Mod1Mask))
		      || (RB_CIR2_MODE == m_nCursorMode)
		      ||  (   (ShiftMask | ControlMask)
			  == (pXevent->xbutton.state & (ShiftMask | ControlMask)) ) )
		    {
		      // Alt (or Shift+Ctrl) pressed, a slightly different circle, cursor stays at center

		      m_poCursor->vSetFont(XC_target, TRUE);
		      m_poRubberband = new CXrubberband(m_hParent, &m_hPixmap,
							nX, nY, m_nWidthUsed,
							m_nHeightUsed, 
							RB_CIR2_MODE);
		    }
		  else if (   (0 != (pXevent->xbutton.state & ShiftMask) )
                           || (RB_TLIN_MODE == m_nCursorMode) )
		    {
		      // Only shift pressed
		  
		      m_poCursor->vSetFont(XC_hand2, TRUE);
		      m_poRubberband = new CXrubberband(m_hParent, &m_hPixmap,
							nX, nY, m_nWidthUsed,
							m_nHeightUsed, 
							RB_TLIN_MODE);
		      vSetLineWidth(m_fLineWidth);
		    }
		  else if (   (0 != (pXevent->xbutton.state & ShiftMask) )
                           || (RB_LINE_MODE == m_nCursorMode) )
		    {
		      // Only shift pressed
		  
		      m_poCursor->vSetFont(XC_hand2, TRUE);
		      m_poRubberband = new CXrubberband(m_hParent, &m_hPixmap,
							nX, nY, m_nWidthUsed,
							m_nHeightUsed, 
							RB_LINE_MODE);
		      vSetLineWidth(0.0);
		    }
		  else if (   (0 != (pXevent->xbutton.state & ControlMask) )
                           || (RB_CIRC_MODE == m_nCursorMode) )
		    {
		      // Ctrl pressed, a circle, origin stays at center

		      m_poCursor->vSetFont(XC_circle, TRUE);
		      m_poRubberband = new CXrubberband(m_hParent, &m_hPixmap,
							nX, nY, m_nWidthUsed,
							m_nHeightUsed, 
							RB_CIRC_MODE);
		    }
		  else
		    {
		      // No modifier buttons pressed or unknown m_nCursorMode,
		      // initialize zoom box drag 
		  
		      m_poCursor->vSetFont(XC_sizing, TRUE);
		      m_poRubberband = new CXrubberband(m_hParent, &m_hPixmap,
							nX, nY, m_nWidthUsed,
							m_nHeightUsed, 
							RB_RECT_MODE);
		    }
		}
	    }
	}
      else if (pXevent->xany.type == ButtonRelease) 
	{
	  vPixelInfo(nX, nY, 2, pXevent->xbutton.button);
	  if (NULL != m_poRubberband)
	    {
	      if (   ( (Button2 == pXevent->xbutton.button) && (RB_NONE_MODE == m_nCursorMode) )
                  || ( (Button1 == pXevent->xbutton.button) && (RB_NONE_MODE != m_nCursorMode) ) )
		{
		  // End of rubberbanding mode...

		  if (   (RB_RECT_MODE == m_nCursorMode) 
                      || (RB_CIR2_MODE == m_nCursorMode) )
		    m_nCursorMode = RB_NONE_MODE;

		  int nX1, nY1, nX2, nY2;
		  int nMode = m_poRubberband->nGetMode();
		  m_poRubberband->vGetCoords(&nX1, &nY1, &nX2, &nY2);
		  delete m_poRubberband;
		  m_poRubberband = NULL;

		  if (RB_LINE_MODE == m_nCursorMode)
		    {
/*
		      //Temp kludge for drawing quadrilateral

		      float fQuad0, fQuad1, fQuad3, fQuad4;
		      if (0 == nDpyPixToImgPix( nX1, nY1, &fQuad0, &fQuad1))
			{
			  if (0 == nDpyPixToImgPix( nX2, nY2, &fQuad3, &fQuad4))
			    {
			      // Add to quadrilateral
			      int iq;
			      for (iq = 1; iq < 4; iq++)
				{
				  if ( (0.0 > m_a5fQuad[iq][0]) || (0.0 > m_a5fQuad[iq][1]) )
				    {
				      m_a5fQuad[iq][0] = fQuad3;
				      m_a5fQuad[iq][1] = fQuad4;
				      break;
				    }
				}
			      if (1 == iq)
				{
				  m_a5fQuad[0][0] = fQuad0;
				  m_a5fQuad[0][1] = fQuad1;
				}
			      vDrawResol();
			    }
			}
*/
		    }
		  else if (RB_CIRC_MODE == m_nCursorMode)
		    {
		      // Reset quadrilateral

		      m_a5fQuad[0][0] = -1.0f;
		      m_a5fQuad[0][1] = -1.0f;
		      m_a5fQuad[1][0] = -1.0f;
		      m_a5fQuad[1][1] = -1.0f;
		      m_a5fQuad[2][0] = -1.0f;
		      m_a5fQuad[2][1] = -1.0f;
		      m_a5fQuad[3][0] = -1.0f;
		      m_a5fQuad[3][1] = -1.0f;
		      m_a5fQuad[4][0] = -1.0f;
		      m_a5fQuad[4][1] = -1.0f;
		    }

		  //	      m_poCursor->vRestore();
		  if (RB_NONE_MODE == m_nCursorMode)
		    m_poCursor->vReset();

		  if (0 == (pXevent->xbutton.state & Button3Mask))
		    {
		      // If button 3 is held down when button 2 released, then
		      // that is abort rubberbanding, so do not enter this section
		  
		      int   nNewOrig[2], nNewExt[2];
		      float fNewOrig[2], fNewExt[2];

		      // TO DO:  What if nDpyPixToImgPix fails in the following!

		      nDpyPixToImgPix( nX1, nY1, &fNewOrig[0], &fNewOrig[1]);
		      nDpyPixToImgPix( nX2, nY2, &fNewExt[0], &fNewExt[1]);
		      if ( (NULL != m_pObj) && (NULL != m_prvRubberbandCallback) )
			{
			  m_prvRubberbandCallback(m_pObj, nMode,
						  fNewOrig[0], fNewOrig[1],
						  fNewExt[0], fNewExt[1]);
			}
		      if (RB_RECT_MODE == nMode)
			{
			  // Rectangle mode 
			  
			  nNewOrig[0] = (int) fNewOrig[0];
			  nNewOrig[1] = (int) fNewOrig[1];
			  nNewExt[0]  = (int) fNewExt[0];
			  nNewExt[1]  = (int) fNewExt[1];

			  // Because of the 8 different orientations, the origin and
			  // extents must be adjusted, so that the origin comes
			  // before the extents.
			  // THIS SHOULD BE IN THE SetRegion ROUTINE!

			  if (nNewExt[0] > nNewOrig[0])
			    {
			      nNewExt[0] = nNewExt[0] - nNewOrig[0];
			    }
			  else
			    {
			      nX1         = nNewOrig[0];
			      nNewOrig[0] = nNewExt[0];
			      nNewExt[0]  = nX1 - nNewOrig[0];
			    }
			  if (nNewExt[1] > nNewOrig[1])
			    {
			      nNewExt[1] = nNewExt[1] - nNewOrig[1];
			    }
			  else
			    {
			      nX1         = nNewOrig[1];
			      nNewOrig[1] = nNewExt[1];
			      nNewExt[1]  = nX1 - nNewOrig[1];
			    }

			  (void) nSetRegion(nNewOrig[0], nNewOrig[1],
					    nNewExt[0], nNewExt[1]);
			  //	      m_tImageProps.fZoom = 1.0;
			  vRefresh();
			}
		    }
		}
	    }
          //+2012-12-08 JWP
          else if (Button4 == pXevent->xbutton.button)
	    {
	      // If the mouse has a scroll wheel, use that to zoom/unzoom
	      //cout << "Button4 release\n";
	      // Try to unzoom the view area
	      int nCent[2];
	      int nNewExt[2], nNewOrig[2];
	      float fFactor = 0.2;
		  
	      if (    (m_tImageProps.nExt[0] >= m_tImageProps.nDim[0]) 
                   && (m_tImageProps.nExt[1] >= m_tImageProps.nDim[1]) )
		{
		  // If already all the way unzoomed, don't do anything
		  //cout << "Button4 max unzoomed\n";		  
		}
	      else
		{
		  nNewExt[0] = m_tImageProps.nExt[0]  + nint((float)m_tImageProps.nExt[0] * fFactor);
		  nNewExt[1] = m_tImageProps.nExt[1]  + nint((float)m_tImageProps.nExt[1] * fFactor);
		  if (nNewExt[0] == m_tImageProps.nExt[0]) nNewExt[0]++;
		  if (nNewExt[1] == m_tImageProps.nExt[1]) nNewExt[1]++;

		  nCent[0]   = m_tImageProps.nOrig[0] + nint((float)m_tImageProps.nExt[0] * 0.5);  // Find the center
		  nCent[1]   = m_tImageProps.nOrig[1] + nint((float)m_tImageProps.nExt[1] * 0.5);
		  nNewOrig[0] = nCent[0] - nNewExt[0] / 2;
		  nNewOrig[1] = nCent[1] - nNewExt[1] / 2;
		  (void) nSetRegion(nNewOrig[0], nNewOrig[1], nNewExt[0], nNewExt[1]);
		  vRefresh();
		}
	    }
	  else if (Button5 == pXevent->xbutton.button)
	    {
	      //cout << "Button5 release\n";
	      // Try to zoom the view area

	      int nCent[2];
	      int nNewExt[2], nNewOrig[2];
	      float fFactor = 0.2;
	      nCent[0] = -1;
	      nCent[1] = -1;
	      if (    (m_tImageProps.nExt[0] <= 1)
                   && (m_tImageProps.nExt[1] <= 1) )
		{
		  // If already all the way zoomed, don't do anything
		  //cout << "Button5 max zoomed\n";
		}
	      else
		{
		  // If starting out almost completely unzoomed, make the first zoom centered at the cursor position

		  //cout << "Image Ext,Dim 0: " << m_tImageProps.nExt[0] << ", " << m_tImageProps.nDim[0] 
		  //     << "       1: " << m_tImageProps.nExt[1] << ", " <<  m_tImageProps.nDim[1] << endl;

		  if (    (10 > ABS(m_tImageProps.nExt[0] - m_tImageProps.nDim[0])) 
                       && (10 > ABS(m_tImageProps.nExt[1] - m_tImageProps.nDim[1])) )
		    {
		      // ... All the way unzoomed, so change center to this image pixel
		      //cout << "Mx unzoomed: x, y: " << nX << ", " << nY << endl;

		      float fPx0, fPx1;
		      int nStat;
		      nStat = nDpyPixToImgPix( nX, nY, &fPx0, &fPx1);
		      if (   (nX >= 0) && (nY >=0) 
                          && (nX < m_nWidthUsed) && (nY < m_nHeightUsed) 
			  && (nStat == 0) )
			{
			  // Legit fPx0, fPx1
			  nCent[0] = nint(fPx0);
			  nCent[1] = nint(fPx1);
			  fFactor = 0.5;
			}
		    }
		  
		  nNewExt[0] = m_tImageProps.nExt[0]  - nint((float)m_tImageProps.nExt[0] * fFactor);
		  nNewExt[1] = m_tImageProps.nExt[1]  - nint((float)m_tImageProps.nExt[1] * fFactor);
		  if (nNewExt[0] == m_tImageProps.nExt[0]) nNewExt[0]--;
		  if (nNewExt[1] == m_tImageProps.nExt[1]) nNewExt[1]--;
		  if (0 >= nNewExt[0]) nNewExt[0] = 1;
		  if (0 >= nNewExt[1]) nNewExt[1] = 1;

		  if (-1 == nCent[0]) // Center was not changed previously, so change it now
		    nCent[0]   = m_tImageProps.nOrig[0] + nint((float)m_tImageProps.nExt[0] * 0.5);  // Find the center
		  if (-1 == nCent[1]) // Center was not changed previously, so change it now
		    nCent[1]   = m_tImageProps.nOrig[1] + nint((float)m_tImageProps.nExt[1] * 0.5);
		  if (1 < nNewExt[0])
		    nNewOrig[0] = nCent[0] - nNewExt[0] / 2;
		  else
		    nNewOrig[0] = m_tImageProps.nOrig[0];
		  if (1 < nNewExt[1])
		    nNewOrig[1] = nCent[1] - nNewExt[1] / 2;
		  else
		    nNewOrig[1] = m_tImageProps.nOrig[1];

		  (void) nSetRegion(nNewOrig[0], nNewOrig[1], nNewExt[0], nNewExt[1]);
		  vRefresh();
		}
	    }
          //-2012-12-08 JWP
        }
      else if (pXevent->xany.type == MotionNotify) 
	{
	  XEvent hXLookAheadEvent;
	  
	  // Try to get only the last MotionNotify event in the queue in
	  // order to speed things up.

	  while (XCheckTypedWindowEvent (XtDisplay(m_hParent), XtWindow(m_hParent),
                MotionNotify, &hXLookAheadEvent))
	    {
	      nX = hXLookAheadEvent.xmotion.x;
	      nY = hXLookAheadEvent.xmotion.y;
	    }

	  if (NULL != m_poRubberband)
	    {
	      // Rubber band object exists

	      if ( ((Button1Mask | Button2Mask)
		    == (pXevent->xmotion.state & (Button1Mask | Button2Mask)))
		   || (0 != (pXevent->xmotion.state & (ShiftMask | ControlMask | Mod1Mask))) )
		{
		  if (   (pXevent->xmotion.state & Button2Mask)
                     && !(pXevent->xmotion.state & Button1Mask)
			 && !m_bSpacebarHeldDown && !m_bControlKeyHeldDown)
		    {
		      m_poRubberband->vDraw(nX, nY);			  
		    }
		  else
		    {
		      // Pan if two buttons or space or control key held down
		      m_poRubberband->vPan(nX, nY);
		    }
		}
	      else
		{
		  m_poRubberband->vDraw(nX, nY);
		}
	      vPixelInfo(nX, nY, 0, 2);
	    }
	  else if ((Button1Mask | Button2Mask)
		   == (pXevent->xmotion.state & (Button1Mask | Button2Mask)))
	    {
	      // This means pan image (button1 and button2 both pressed)
	      m_poCursor->vSetFont(XC_fleur, TRUE);
	      vPan(nX, nY);
	      //	      m_poCursor->vRestore();
	      m_poCursor->vReset();
	      vPixelInfo(nX, nY, 0, 12);
	    }
	  else if (m_bSpacebarHeldDown || m_bControlKeyHeldDown)
	    {
	      m_poCursor->vSetFont(XC_fleur, TRUE);
	      vPan(nX, nY);
	      m_poCursor->vReset();
	      vPixelInfo(nX, nY, 0, 1);
	    }
	  else
	    {
	      vPixelInfo(nX, nY, 0, 1);
	    }
	}
      else if (pXevent->xany.type == KeyPress) 
	{
	  // A a keyboard command
	  XKeyEvent *pXKeyevent = (XKeyEvent *) cbs->event;
	  Modifiers hM2;
	  KeySym  hKeysym;
	  XtTranslateKeycode(XtDisplay(m_hParent), pXKeyevent->keycode,
			     pXKeyevent->state, &hM2, &hKeysym);
	  
	  // Mask out high order bits (why?  It was empirically discovered!)

	  hKeysym = hKeysym & 0xFFFF;
//	  cout << "keycode, keysym: " <<  pXKeyevent->keycode << ", " 
//	       << hKeysym << endl;


	  // All the following shifts need modification by the ORIENTATION.
	  // They make sense ONLY when nOrient = 0!

	  if (   (hKeysym == (KeySym)XK_Delete)
	      || (hKeysym == (KeySym)XK_KP_Delete)
	      || (hKeysym == (KeySym)XK_KP_Subtract)
	      || (hKeysym == (KeySym)XK_minus)
	      || (hKeysym == (KeySym)XK_BackSpace) )
	    {
	      if (   (NULL != m_poReflnlist) 
                  && (0 < m_tImageProps.nDisplayReflns) )
		{
		  float fPx0, fPx1;
		  int nX = pXKeyevent->x;
		  int nY = pXKeyevent->y;
		  if (0 == nDpyPixToImgPix(nX, nY, &fPx0, &fPx1))
		    {
		      // nX is temp
		      nX = m_poReflnlist->nDeleteSpot(fPx0, fPx1);
		      if (0 <= nX)
			{
			  int i = m_poReflnlist->m_nFI_nNonunfFlag;
			  if (0 > i)
			    i = m_poReflnlist->m_nFI_nH;
			  m_poReflnlist->vSort(eReflnField_int_type, i, NULL);
			  vRefresh();
			  vPlotReflnlist();
			}
		    }
		}
	    }
	  else if (hKeysym == (KeySym)XK_space)
	    {
	      // Wheneverthe spacebar is pressed down, reset this private variable to true.

	      //cout << "Spacebar down = true.\n" << flush;
	      m_bSpacebarHeldDown = TRUE;
	    }
	  else if (   (hKeysym == (KeySym)XK_Control_L)
		   || (hKeysym == (KeySym)XK_Control_R))
	    {
	      // Whenever the left or right control key is pressed down, reset this is set to true.

	      //cout << "Control key down = true.\n" << flush;
	      m_bControlKeyHeldDown = TRUE;
	    }

	  else if (   (hKeysym == (KeySym)XK_plus)
		   || (hKeysym == (KeySym)XK_Insert)
		   || (hKeysym == (KeySym)XK_KP_Insert)
		   || (hKeysym == (KeySym)XK_KP_Add) )
	    {
	      if (   (NULL != m_poReflnlist) 
                  && (0 < m_tImageProps.nDisplayReflns) )
		{
		  if (0 == m_tSpotInfo.nStat)
		    {
		      m_poReflnlist->vAddSpot(m_tSpotInfo.fCentPx0, 
					      m_tSpotInfo.fCentPx1,
					      m_tSpotInfo.fIntensity,
					      m_tSpotInfo.fSigmaI,
					      m_tSpotInfo.fResolution,
					      m_tSpotInfo.fRotMid,
					      m_tSpotInfo.fRotWidth);
		      m_tSpotInfo.nStat = -1;
		      int i = m_poReflnlist->m_nFI_nNonunfFlag;
		      if (0 > i)
			i = m_poReflnlist->m_nFI_nH;
		      m_poReflnlist->vSort(eReflnField_int_type, i, NULL);
		      vRefresh();
		      vPlotReflnlist();
		    }
		}
	    }
          else if (   (hKeysym == (KeySym)XK_L) || (hKeysym == (KeySym)XK_l)
		   || (hKeysym == (KeySym)XK_less) 
                   || (hKeysym == (KeySym)XK_comma) )
	    {
	      // Try to enlarge the view area

	      (void) nSetRegion(m_tImageProps.nOrig[0] - m_tImageProps.nExt[0]/2, 
				m_tImageProps.nOrig[1] - m_tImageProps.nExt[1]/2,
				m_tImageProps.nExt[0] * 2,
				m_tImageProps.nExt[1] * 2);
	      vRefresh();
	    }
	  else if ( (hKeysym == (KeySym)XK_G) || (hKeysym == (KeySym)XK_g)
		   || (hKeysym == (KeySym)XK_greater)
		   || (hKeysym == (KeySym)XK_period) )
	    {
	      // Try to enlarge the view area

	      (void) nSetRegion(m_tImageProps.nOrig[0] + m_tImageProps.nExt[0]/4, 
				m_tImageProps.nOrig[1] + m_tImageProps.nExt[1]/4,
				m_tImageProps.nExt[0] / 2,
				m_tImageProps.nExt[1] / 2);
	      vRefresh();
	    }
	  else if ( (hKeysym == (KeySym)XK_J) || (hKeysym == (KeySym)XK_j) 
		   || (hKeysym == (KeySym)XK_KP_Left)
		   || (hKeysym == (KeySym)XK_Left) )
	    {
	      // Shift view 1/10th of way to left

	      (void) nSetRegion(m_tImageProps.nOrig[0] - max(1, m_tImageProps.nExt[0]/10), 
				m_tImageProps.nOrig[1],
				m_tImageProps.nExt[0],
				m_tImageProps.nExt[1]);
	      vRefresh();
	    }
	  else if ( (hKeysym == (KeySym)XK_K) || (hKeysym == (KeySym)XK_k)
		   || (hKeysym == (KeySym)XK_KP_Right)
		   || (hKeysym == (KeySym)XK_Right) )
	    {
	      // Shift view 1/10th of way to right

	      (void) nSetRegion(m_tImageProps.nOrig[0] + max(1, m_tImageProps.nExt[0]/10), 
				m_tImageProps.nOrig[1],
				m_tImageProps.nExt[0],
				m_tImageProps.nExt[1]);
	      vRefresh();
	    }
	  else if ( (hKeysym == (KeySym)XK_I) || (hKeysym == (KeySym)XK_i) 
		   || (hKeysym == (KeySym)XK_KP_Up)
		   || (hKeysym == (KeySym)XK_Up) )

	    {
	      // Shift view 1/10th of way up

	      (void) nSetRegion(m_tImageProps.nOrig[0],
				m_tImageProps.nOrig[1] - max(1, m_tImageProps.nExt[1]/10),
				m_tImageProps.nExt[0],
				m_tImageProps.nExt[1]);
	      vRefresh();
	    }
	  else if ( (hKeysym == (KeySym)XK_M) || (hKeysym == (KeySym)XK_m)
		   || (hKeysym == (KeySym)XK_KP_Down)
		   || (hKeysym == (KeySym)XK_Down) )
	    {
	      // Shift view 1/10th of way down

	      (void) nSetRegion(m_tImageProps.nOrig[0],
				m_tImageProps.nOrig[1] + max(1, m_tImageProps.nExt[1]/10),
				m_tImageProps.nExt[0],
				m_tImageProps.nExt[1]);
	      vRefresh();
	    }
	}
      else if (pXevent->xany.type == KeyRelease) 
	{
	  // Do nothing
	  //+JWP 3009-03-14
	  // We need to know when the space key is released
	  // A a keyboard command
	  XKeyEvent *pXKeyevent = (XKeyEvent *) cbs->event;
	  Modifiers hM2;
	  KeySym  hKeysym;
	  XtTranslateKeycode(XtDisplay(m_hParent), pXKeyevent->keycode,
			     pXKeyevent->state, &hM2, &hKeysym);
	  
	  // Mask out high order bits (why?  It was empirically discovered!)

	  hKeysym = hKeysym & 0xFFFF;
	  if (hKeysym == (KeySym)XK_space)
	    {
	      // Wheneverthe spacebar is released, reset this private variable to false.

	      m_bSpacebarHeldDown = FALSE;
	    }
	  else if (   (hKeysym == (KeySym)XK_Control_L)
		   || (hKeysym == (KeySym)XK_Control_R))
	    {
	      m_bControlKeyHeldDown = FALSE;
	    }
	  //-JWP 3009-03-14
	}
      else 
	{
	  cout << "Unknown input event (x, y):"
	       << nX << ", " << nY << endl;
	}
    }
  else if (cbs->reason == XmCR_RESIZE)
    {
      vRefresh();
    }
}

int
CXdisplay::nCreatePixmap(void)
{
  // This converts a Cimage into a X11 Pixmap.

  // Get width and height of parent widget

  XtVaGetValues(m_hParent, XmNwidth, &m_nWidth, XmNheight, &m_nHeight, NULL);

  // Set used area, zoom factor could be less than 1.0

  if (NULL == m_poImage) return (2);

  m_nWidthUsed  = (int) (m_tImageProps.fZoom * (float)m_nWidth);
  if (m_nWidthUsed > m_nWidth) m_nWidthUsed = m_nWidth;
  m_nHeightUsed = (int) (m_tImageProps.fZoom * (float)m_nHeight);
  if (m_nHeightUsed > m_nHeight) m_nHeightUsed = m_nHeight;

  // Delete any previously used objects (m_phXImage, m_hPixmap)

  vDeleteXStuff();

  // Create an XImage

  m_phXImage = XCreateImage(XtDisplay(m_hParent),
			    DefaultVisual(XtDisplay(m_hParent),
					  DefaultScreen(XtDisplay(m_hParent))),
			    m_ulDepthImage, ZPixmap, 0, (char *)NULL,
			    m_nWidthUsed, m_nHeightUsed, XBitmapPad(XtDisplay(m_hParent)),
			    0);
  if ((XImage*)NULL == m_phXImage)
    {
      // You are in big trouble!

      cout << "ERROR creating image with XCreateImage!\n" << flush;
    }

  // Convert Cimage to a byte array of color indices, scaled by ...
  // vImageToCLUT also allocates and set XImage->data pointer.

  vImageToCLUT(*m_poImage,
	       m_tImageProps.fScaleMin, m_tImageProps.fScaleMax,
	       m_nMin, m_nMax);

  // Create pixmap for the image for quick refresh on expose

  m_hPixmap = XCreatePixmap(XtDisplay(m_hParent),
			    XtWindow(m_hParent),
			    m_nWidthUsed, m_nHeightUsed, 
			    DefaultDepthOfScreen(XtScreen(m_hParent)));

/*
  cout << "Before putting XImage to pixmap\n";
  int jj;
  unsigned short int *puiJtest;
  puiJtest = (unsigned short int*) m_phXImage->data;
  for (jj=0; jj < 10; jj++) 
    cout << "jj, uiJ: " << jj << ": " << *puiJtest++ << endl;
  cout << endl << flush;
*/

  if (   m_bPixelSwapBytes && (2 == m_nBytesPerPixel))
    {
      // Swap bytes on the 2-byte pixels
      vSwapBytes(m_nDataSize, m_pucData);
    }

  // Put the image into the pixmap

  XPutImage(XtDisplay(m_hParent),
	    m_hPixmap,
	    m_hGC,
	    m_phXImage,
	    0, 0, 0, 0, m_nWidthUsed, m_nHeightUsed);
/*
  //+test
  cout << "about to getimage\n" << flush;
  XImage      *m_phXNewImage;
  m_phXNewImage = XGetImage(XtDisplay(m_hParent), m_hPixmap, 0, 0, 
			    m_nWidthUsed, m_nHeightUsed, 0xFFFFFF, ZPixmap);
  cout << "props of gotten image:"
       << "\ndepth:          " <<  m_phXNewImage->depth
       << "\nbytes_per_line: " << m_phXNewImage->bytes_per_line 
       << "\nbitmap_pad:     " <<  m_phXNewImage->bitmap_pad
       << "\nbits_per_pixel: " <<  m_phXNewImage->bits_per_pixel
       << endl << flush;

  puiJtest = (unsigned short int*) m_phXNewImage->data;
  for (jj=0; jj < 10; jj++) 
    cout << "jj, uiJ: " << jj << ": " << *puiJtest++ << endl;
  cout << endl << flush;
  XDestroyImage(m_phXNewImage);
  //-test
*/
  // Plot reflection list

  if (m_tImageProps.nDisplayReflns > 0)    vPlotReflnlist();

  // Plot pixel value list

  if (m_tImageProps.nDisplayPixValues > 0) vPlotPixelValues();

  vDrawResol();  // Draw resolution circles onto the pixmap

  return (0);
}

int
CXdisplay::nSetRegion(const int nOrigIn0, const int nOrigIn1,
		      const int nExtIn0,  const int  nExtIn1)
{
  // Set the origin and the extents of the input Cimage to display
  // Coordinate system is the system of the input Cimage.
  // This function does not resize or changed the display image.

  // Do some error checking...

  if ( (nExtIn0 <= 0) || (nExtIn1 <= 0) )
    {
      // Reset to origin and extents to full size
      m_tImageProps.nOrig[0] = 0;
      m_tImageProps.nOrig[1] = 0;
      m_tImageProps.nExt[0]  = m_tImageProps.nDim[0];
      m_tImageProps.nExt[1]  = m_tImageProps.nDim[1];
      return (0);
    }
/*
  else if (   (nOrigIn0 < 0)
      || (nOrigIn1 < 0)
      || (nExtIn0 > m_tImageProps.nDim[0])
      || (nExtIn1 > m_tImageProps.nDim[1]) )
    {
      return (-1);
    }
*/
  else
    {
      m_tImageProps.nOrig[0] = max(nOrigIn0, 0);
      m_tImageProps.nOrig[1] = max(nOrigIn1, 0);
      m_tImageProps.nExt[0]  = min(nExtIn0, (m_tImageProps.nDim[0] - m_tImageProps.nOrig[0]));
      m_tImageProps.nExt[1]  = min(nExtIn1, (m_tImageProps.nDim[1] - m_tImageProps.nOrig[1]));

      return (0);
    }
};

int
CXdisplay::nDpyPixToImgPix(const int nX, const int nY, float *pf1, float *pf2)
{
  // Convert the display pixel coords to the image coordinates
  // Return error if the coords are not in the image.
  // Do image pixels as float.

  float fPx0, fPx1;

//  fPx0 = (float) m_nStart0 + (((float) nX + 0.5) * m_fStep0);
//  fPx1 = (float) m_nStart1 + (((float) nY + 0.5) * m_fStep1);
  fPx0 = (float) (m_nStart0)  +  (float) nX * m_fStep0; 
  fPx1 = (float) (m_nStart1)  +  (float) nY * m_fStep1;

  if (0 > m_nDir0)
    fPx0 = fPx0 + 1.0;
  if (0 > m_nDir1)
    fPx1 = fPx1 + 1.0;
  if (4 > m_tImageProps.nOrient)
    {
      *pf1 = fPx0;
      *pf2 = fPx1;
    }
  else
    {
      *pf1 = fPx1;
      *pf2 = fPx0;
    }

  if (    (0.0 <= *pf1) && (0.0 <= *pf2)
      && (*pf1 < (float) m_tImageProps.nDim[0])
      && (*pf2 < (float) m_tImageProps.nDim[1]))
    return (0);
  else
    return (-1);
}

int
CXdisplay::nImgPixToDpyPix(const float fPx0In, const float fPx1In,
			   int *pnX, int *pnY)
{
  // Convert the image pixel coords to the display coordinates
  // Return error if the coords are not in the display.
  // Do image pixels as float.

  float fPx0, fPx1;

  if (4 > m_tImageProps.nOrient)
    {
      fPx0 = fPx0In;
      fPx1 = fPx1In;
    }
  else
    {
      fPx0 = fPx1In;
      fPx1 = fPx0In;
    }

  if (0 > m_nDir0)
    fPx0 = fPx0 - 1.0;
  if (0 > m_nDir1)
    fPx1 = fPx1 - 1.0;

  *pnX = (int) ((fPx0 - (float) m_nStart0) / m_fStep0 + 0.5);
  *pnY = (int) ((fPx1 - (float) m_nStart1) / m_fStep1 + 0.5);

  if (   (0 <= *pnX) && (0 <= *pnY) 
      && (*pnX < m_nWidthUsed) && (*pnY < m_nHeightUsed))
    return (0);
  else
    return (-1);
}

void
CXdisplay::vRefresh(void)
{
  // Refresh the pixmap from Cimage and the window from the pixmap.

  // Re-create the pixmap.

  if (0 != nCreatePixmap()) return;

  m_poCursor->vSetWait();

  // Install colormap just in case something happened
  // Should no longer be needed since XtSetWMColormapWindows
  // have been called.

  //cout << "CXcolormap::vInstall() about to be called." << endl;

  if ("" != sGetEnv("DTDISPLAY_REFRESH_CMAP"))
    m_poColormap->vInstall();

  // Copy the entire Pixmap to the window
  // (This is done despite expose elsewhere which also does it!)

  if ( (0 < m_nWidthUsed) && (0 < m_nHeightUsed) )
    {
      XCopyArea (XtDisplay(m_hParent), m_hPixmap, 
		 XtWindow(m_hParent), m_hGC, 
		 0, 0, m_nWidthUsed, m_nHeightUsed, 0, 0);
#ifdef DEBUG
  XSync(XtDisplay(m_hParent), False);
  XmUpdateDisplay(m_hParent);
#endif

    }

  if (1.0 > m_tImageProps.fZoom)
    {
      vDrawProfile();
    }
  vDrawColorScale();

  if (0 <= m_tImageProps.nResoCircles)
    vDrawResol();
  //  m_poCursor->vRestore();
  if (0 != m_nBeamstopShadowMode)
    m_poCursor->vSetFont(XC_gumby, FALSE);
  else
    m_poCursor->vReset();
}

void
CXdisplay::vPixelInfo(const int nXIn, const int nYIn, const int nInfoMode,
		      const int nButton)
{
  // List out x,y and intensity at that point
  // Update cursor display
  // 
  // nInfoMode = 0  motion_notify, default
  //           = 1  button press
  //           = 2  button release

  int   nStat;
  float fValue;
  float fPx0, fPx1;
  int   nXT = nXIn;
  int   nYT = nYIn;

  // Are static for this object or all objects???  Ans: For all objects, so bug!

  static int nTextHeight = 1;
  static int nTextWidth  = 1;
  static int nTextX = 1;
  static int nTextY = 1;

  int nX, nY;
  int nTextLen;

  if (-1 == nXT) nXT = m_nXPrev;
  if (-1 == nYT) nYT = m_nYPrev;

  //  Clear previous string at old position
      
  XCopyArea (XtDisplay(m_hParent), m_hPixmap, 
	     XtWindow(m_hParent), m_hGC, 
	     nTextX, nTextY-nTextHeight, nTextWidth, nTextHeight+5,
	     nTextX, nTextY-nTextHeight);
#ifdef DEBUG
  XSync(XtDisplay(m_hParent), False);
  XmUpdateDisplay(m_hParent);
#endif

  if (NULL == m_poImage) return;

  // Clear pixmap to where previous cursor might have been,
  // except if rubberbanding (but leave displayed if button release!)

  if ( (2 != nInfoMode) && (NULL == m_poRubberband) )
    {
      XCopyArea (XtDisplay(m_hParent), m_hPixmap, 
		 XtWindow(m_hParent), m_hGC, 
		 0, m_nYPrev, m_nWidthUsed+1, 1,
		 0, m_nYPrev);
#ifdef DEBUG
  XSync(XtDisplay(m_hParent), False);
  XmUpdateDisplay(m_hParent);
#endif

      XCopyArea (XtDisplay(m_hParent), m_hPixmap, 
		 XtWindow(m_hParent), m_hGC, 
		 m_nXPrev, 0, 1, m_nHeightUsed+1,
		 m_nXPrev, 0);
#ifdef DEBUG
  XSync(XtDisplay(m_hParent), False);
  XmUpdateDisplay(m_hParent);
#endif

      // If you wanted to leave an extra mark at the previous cursor
      // position, this is where you would put it.  You would have to
      // save the previous,previous position, so you could clear it though.

    }

  nStat = nDpyPixToImgPix( nXT, nYT, &fPx0, &fPx1);

  if (   (nXT >= 0) && (nYT >=0) 
      && (nXT < m_nWidthUsed) && (nYT < m_nHeightUsed) 
      && (nStat == 0) )
    {
      nX = 0;
      nY = 0;

      m_fPx0Prev = fPx0;   // Save image pixel coordinates for use in other
      m_fPx1Prev = fPx1;   // methods such a vPan

      fValue = (m_poImage->*m_poImage->prfGetPixel)((int)fPx0, (int)fPx1);
      char cTemp[20];

      if (!m_bHasFloatValues)
	sprintf(cTemp, "%.1f, %.1f, %.0f", fPx0, fPx1, fValue);
      else if (fValue < 10.0f)
	sprintf(cTemp, "%.1f, %.1f, %.4f", fPx0, fPx1, fValue);
      else
	sprintf(cTemp, "%.1f, %.1f, %.0f", fPx0, fPx1, fValue);

      nTextLen = strlen(cTemp);
      
      if (nInfoMode != 2)
	{
	  // Draw new string at new position
	  
	  nTextHeight = m_pXFont->ascent + m_pXFont->descent;
	  nTextWidth  = XTextWidth(m_pXFont, cTemp, nTextLen);

	  //  Then next statement is for when the text does NOT track
	  //  with the mouse, but is fixed at 0,0

	  if ( (nXT < nTextWidth) && (nYT < nTextHeight) )
	    nX = m_nWidthUsed - nTextWidth;

	  if ( (nX + nTextWidth) > m_nWidthUsed)
	    nX = m_nWidthUsed - nTextWidth;
	  if ( (nY - nTextHeight) < 0)
	    nY = nTextHeight;
	  if (nY > m_nHeightUsed)
	    nY = m_nHeightUsed - nTextHeight;
	  XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		      nX, nY, cTemp, nTextLen);
	}
      // Save position of this string so it can be cleared next time
      
      nTextX = nX;
      nTextY = nY;
      
      if ( (2 > nInfoMode) && (NULL == m_poRubberband) )
	{
	  // Draw cursor lines on the WINDOW, not on the PIXMAP!!
	  //  (but only if rubberbanding is not active)
	  
	  XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		    0, nYT, m_nWidthUsed, nYT);
	  
	  XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		    nXT, 0, nXT, m_nHeightUsed);
	}

      // Save position of cursor so it can be cleared next time

      m_nXPrev = nXT;
      m_nYPrev = nYT;
    }

  if (   (nXT >= 0) && (nYT >=0) 
      && (nXT < m_nWidthUsed) && (nYT < m_nHeight) )
    {
      if ( (1.0 > m_tImageProps.fZoom) && (nYT < m_nHeightUsed) )
	{
	  vDrawProfile();
	  m_nScaleMode = 0;
	  if (   (1 == nButton)
              && (2 == nInfoMode)
              && (1 == m_nBeamstopShadowMode)
	      )
		{
		  // Button release, in image area and beamstop shadow mode
		  nCalcDrawBeamstopShadow(fPx0, fPx1);
		  //Do NOT turn off, let the calling program turn it off:
                  //vSetBeamstopShadowMode(0);
		}

	}
      else if (1 == nButton)
	{
	  if (2 == nInfoMode)
	    {
	      // Button release with color scale
	      vDrawColorScale(3, nXT, nYT);
	    }
	  else
	    {
	      // Cursor 
	      vDrawColorScale(2, nXT, nYT);
	    }
	}
    }
}

int
CXdisplay::nSetImageProps(const tagImagePropsDisplay* ptInputProps)
{
  m_tImageProps.fZoom             = ptInputProps->fZoom;
  m_tImageProps.fScaleMin         = ptInputProps->fScaleMin;
  m_tImageProps.fScaleMax         = ptInputProps->fScaleMax;
  m_tImageProps.fScaleMinSD       = ptInputProps->fScaleMinSD;
  m_tImageProps.fScaleMaxSD       = ptInputProps->fScaleMaxSD;

  m_tImageProps.fImageMin         = ptInputProps->fImageMin;
  m_tImageProps.fImageMax         = ptInputProps->fImageMax;
  m_tImageProps.fAvg              = ptInputProps->fAvg;
  m_tImageProps.fSD               = ptInputProps->fSD;
  m_tImageProps.fAspectRatio      = ptInputProps->fAspectRatio;

//
// These should not be set with this routine.  They are read-only
//  m_tImageProps.nDim[0]           = ptInputProps->nDim[0];
//  m_tImageProps.nDim[1]           = ptInputProps->nDim[1];
//
// These should be set with a call to nSetRegion()
//  m_tImageProps.nOrig[0]          = ptInputProps->nOrig[0];
//  m_tImageProps.nOrig[1]          = ptInputProps->nOrig[1];
//  m_tImageProps.nExt[0]           = ptInputProps->nExt[0];
//  m_tImageProps.nExt[1]           = ptInputProps->nExt[1];

  m_tImageProps.nIntegrateBox[0]  = ptInputProps->nIntegrateBox[0];
  m_tImageProps.nIntegrateBox[1]  = ptInputProps->nIntegrateBox[1];
  m_tImageProps.nAutoRescale      = ptInputProps->nAutoRescale;
  m_tImageProps.nOrient           = ptInputProps->nOrient;
  m_tImageProps.nAspectMode       = ptInputProps->nAspectMode;
  m_tImageProps.nDisplayReflns    = ptInputProps->nDisplayReflns;
  m_tImageProps.nDisplayPixValues = ptInputProps->nDisplayPixValues;
  m_tImageProps.nResoCircles      = ptInputProps->nResoCircles;
  m_tImageProps.nPlotRockingCurve = ptInputProps->nPlotRockingCurve;
  m_tImageProps.fSatPixValue      = ptInputProps->fSatPixValue;
  m_tImageProps.sSatPixColor      = ptInputProps->sSatPixColor;

  vSetSatPixColor(m_tImageProps.sSatPixColor.string());

  // We should perform error checking on the above values

  vRefresh();
  return (0);
}

int
CXdisplay::nGetImageProps(tagImagePropsDisplay *ptOutputProps)
{
  ptOutputProps->fZoom             = m_tImageProps.fZoom;
  ptOutputProps->fScaleMin         = m_tImageProps.fScaleMin;
  ptOutputProps->fScaleMax         = m_tImageProps.fScaleMax;
  ptOutputProps->fScaleMinSD       = m_tImageProps.fScaleMinSD;
  ptOutputProps->fScaleMaxSD       = m_tImageProps.fScaleMaxSD;
  ptOutputProps->fImageMin         = m_tImageProps.fImageMin;
  ptOutputProps->fImageMax         = m_tImageProps.fImageMax;

  ptOutputProps->fAvg              = m_tImageProps.fAvg;
  ptOutputProps->fSD               = m_tImageProps.fSD;
  ptOutputProps->fAspectRatio      = m_tImageProps.fAspectRatio;
  ptOutputProps->nDim[0]           = m_tImageProps.nDim[0];
  ptOutputProps->nDim[1]           = m_tImageProps.nDim[1];
  ptOutputProps->nOrig[0]          = m_tImageProps.nOrig[0];
  ptOutputProps->nOrig[1]          = m_tImageProps.nOrig[1];
  ptOutputProps->nExt[0]           = m_tImageProps.nExt[0];
  ptOutputProps->nExt[1]           = m_tImageProps.nExt[1];
  ptOutputProps->nIntegrateBox[0]  = m_tImageProps.nIntegrateBox[0];
  ptOutputProps->nIntegrateBox[1]  = m_tImageProps.nIntegrateBox[1];
  ptOutputProps->nAutoRescale      = m_tImageProps.nAutoRescale;
  ptOutputProps->nOrient           = m_tImageProps.nOrient;
  ptOutputProps->nAspectMode       = m_tImageProps.nAspectMode;
  ptOutputProps->nDisplayReflns    = m_tImageProps.nDisplayReflns;
  ptOutputProps->nDisplayPixValues = m_tImageProps.nDisplayPixValues;
  ptOutputProps->nPlotRockingCurve = m_tImageProps.nPlotRockingCurve;
  ptOutputProps->nResoCircles      = m_tImageProps.nResoCircles;
  ptOutputProps->fSatPixValue      = m_tImageProps.fSatPixValue;
  ptOutputProps->fTooLowPixValue   = m_tImageProps.fTooLowPixValue;
  ptOutputProps->sSatPixColor      = m_tImageProps.sSatPixColor;
  ptOutputProps->sTooLowPixColor   = m_tImageProps.sTooLowPixColor;
  return (0);
}

int
CXdisplay::nSetReflnProps(const tagReflnPropsDisplay* ptInputProps, bool bRefresh)
{
  m_tReflnProps.a6fReflnSize[0] = ptInputProps->a6fReflnSize[0];
  m_tReflnProps.a6fReflnSize[1] = ptInputProps->a6fReflnSize[1];
  m_tReflnProps.a6fReflnSize[2] = ptInputProps->a6fReflnSize[2];
  m_tReflnProps.a6fReflnSize[3] = ptInputProps->a6fReflnSize[3];
  m_tReflnProps.a6fReflnSize[4] = ptInputProps->a6fReflnSize[4];
  m_tReflnProps.a6fReflnSize[5] = ptInputProps->a6fReflnSize[5];
  m_tReflnProps.sReflnColor     = ptInputProps->sReflnColor;
  m_tReflnProps.nReflnPredObs   = ptInputProps->nReflnPredObs;
  m_tReflnProps.nReflnSymbol    = ptInputProps->nReflnSymbol;
  m_tReflnProps.fImageRotStart  = ptInputProps->fImageRotStart; 
  m_tReflnProps.fImageRotEnd    = ptInputProps->fImageRotEnd; 
  
  if (bRefresh)
    {
      vRefresh();
      vPlotReflnlist();
    }
  return (0);
}

int
CXdisplay::nGetReflnProps(tagReflnPropsDisplay *ptOutputProps)
{
  ptOutputProps->a6fReflnSize[0] = m_tReflnProps.a6fReflnSize[0];
  ptOutputProps->a6fReflnSize[1] = m_tReflnProps.a6fReflnSize[1];
  ptOutputProps->a6fReflnSize[2] = m_tReflnProps.a6fReflnSize[2];
  ptOutputProps->a6fReflnSize[3] = m_tReflnProps.a6fReflnSize[3];
  ptOutputProps->a6fReflnSize[4] = m_tReflnProps.a6fReflnSize[4];
  ptOutputProps->a6fReflnSize[5] = m_tReflnProps.a6fReflnSize[5];
  ptOutputProps->sReflnColor     = m_tReflnProps.sReflnColor;
  ptOutputProps->nReflnPredObs   = m_tReflnProps.nReflnPredObs;
  ptOutputProps->nReflnSymbol    = m_tReflnProps.nReflnSymbol;
  ptOutputProps->fImageRotStart  = m_tReflnProps.fImageRotStart;
  ptOutputProps->fImageRotEnd    = m_tReflnProps.fImageRotEnd;
  return(0);
}

void
CXdisplay::vPlotReflnlist(Creflnlist *poReflnlistIn)
{
  int         i, nStat, nNumReflns;
  int         nDpyPx0, nDpyPx1;
  int         nDpyCalcPx0, nDpyCalcPx1;
  int         nSize0, nSize1;
  int         nSize2, nSize3;
  int         nSize4, nSize5;
  int         nFI_fPx0, nFI_fPx1;
  int         nFI_fCalcPx0, nFI_fCalcPx1;
  int         nFI_fCalcRotStart, nFI_fCalcRotEnd;
  int         nFI_nNonunfFlag;
  bool        bDoColors;
  float       fObsPx0, fObsPx1, fCalcPx0, fCalcPx1;
  float       fMinRawPixOKValue = 1.0;
  Crefln     *poRefln;
  static int  nErrorMsg = 0;

  int nFI_fEllipseAxisMajor;
  int nFI_fEllipseAxisMinor;
  int nFI_fEllipseAxisMajorOffset;
  int nFI_fObsRotMid = -1;
  int nFI_fObsRotWidth = -1;
  float fAxisMajor, fAxisMinor, fAxisMajorOffset;
  
  int         nSymbol;
  float fImageRotStart, fImageRotEnd;
  float fRefRotStart, fRefRotEnd;

  if (NULL != poReflnlistIn)
    {
      // This is a pointer to a new reflection list which we should
      // keep around.

      m_poReflnlist   = poReflnlistIn;
    }
  else if (NULL == m_poReflnlist)
    // No input reflection list, see if there is an old one available...
    {
      return;  // The stored reflection list has gone away;
    }
  
  if (!m_poReflnlist->bIsAvailable())
    return;

  nFI_fPx0          = m_poReflnlist->m_nFI_fObsPx0;
  nFI_fPx1          = m_poReflnlist->m_nFI_fObsPx1;
  nFI_fCalcPx0      = m_poReflnlist->m_nFI_fCalcPx0;
  nFI_fCalcPx1      = m_poReflnlist->m_nFI_fCalcPx1;
  nFI_fCalcRotStart = m_poReflnlist->m_nFI_fCalcRotStart;
  nFI_fCalcRotEnd   = m_poReflnlist->m_nFI_fCalcRotEnd;
  nFI_nNonunfFlag   = m_poReflnlist->m_nFI_nNonunfFlag;

  //+2009-12-07 JWP Sometimes the reflnlist is still bogus
  if (   (100 < nFI_fPx0) 
      && (100 < nFI_fCalcPx0) )
    return; // Suspect something is wrong 
  // The above does prevent core dumps for now, so keep it.
  //-2009-12-07 JWP Sometimes the reflnlist is still bogus
 
  nFI_fEllipseAxisMajor        = m_poReflnlist->nGetFieldIndex("fEllipseAxisMajor");
  nFI_fEllipseAxisMinor        = m_poReflnlist->nGetFieldIndex("fEllipseAxisMinor");
  nFI_fEllipseAxisMajorOffset  = m_poReflnlist->nGetFieldIndex("fEllipseAxisMajorOffset");

  bDoColors         = (0 <= nFI_nNonunfFlag);
  if ("" != sGetEnv("DTREK_NONUNF_OKVALUE"))  // Should come from image header?
    {
      fMinRawPixOKValue  =  atof(sGetEnv("DTREK_NONUNF_OKVALUE").string());
    }

  // Try to do something about the fact the 360 == 0 when dealing with angles

  float fMaxAngle = 360.0;

  fImageRotStart = m_tReflnProps.fImageRotStart + 720.0;
  while (fMaxAngle <= fImageRotStart) fImageRotStart -= 360.0;

  if ((fImageRotStart < 20.0) || (fImageRotStart > 340.0))
    {
      // Make a shift away from 0 degrees

      fMaxAngle = 180.0;
      fImageRotStart = m_tReflnProps.fImageRotStart + 720.0;
      while (fMaxAngle <= fImageRotStart) fImageRotStart -= 360.0;
    }

  fImageRotEnd = m_tReflnProps.fImageRotEnd + 720.0;
  while (fMaxAngle <= fImageRotEnd) fImageRotEnd -= 360.0;

  if (fImageRotEnd < fImageRotStart) 
    fImageRotEnd += 360.0;

/*
  cout << "A. Image rot start: " << m_tReflnProps.fImageRotStart 
       << "\nImage rot end: " << m_tReflnProps.fImageRotEnd << endl;

  cout << "B. Image rot start: " << fImageRotStart 
       << "\nImage rot end: " << fImageRotEnd << endl;
*/

  if ( (0 > nFI_fPx0) || (0 > nFI_fPx1) )
    {
      // Error, the requested fields do not exist in the Reflnlist, so look
      // for calculated fields

      if (0 == nErrorMsg)
	{
	  cerr << "Observed pixel fields not found, looking for "
               << "calculated fields ..."
               << endl;
	}

      nFI_fPx0 = m_poReflnlist->m_nFI_fCalcPx0;
      nFI_fPx1 = m_poReflnlist->m_nFI_fCalcPx1;

      if ( (0 > nFI_fPx0) || (0 > nFI_fPx1) )
	{
	  if (0 == nErrorMsg)
	    {
	      cerr << "  ... not found, no reflections plotted." << endl;
	    }
	  return;
	}
      else if (0 == nErrorMsg)
	{
	  cerr << "  ... found and plotted." << endl;
	}
      nErrorMsg = 1;
    }
  nNumReflns = m_poReflnlist->nGetNumReflns();
  if (0 >= nNumReflns)
    return;  // Nothing to plot, so ...

  m_poCursor->vSetWait();



  nSize0     = nint(ABS(m_tReflnProps.a6fReflnSize[0] * 0.5 / m_fStep0));
  nSize1     = nint(ABS(m_tReflnProps.a6fReflnSize[1] * 0.5 / m_fStep1));

  nSize2 = nSize3 = nSize4 = nSize5 = 0;

  bool bDrawSpot = FALSE;
  bool bDrawBkg  = FALSE;
  bool bDrawEllipse = FALSE;

  if (   (0.0 < m_tReflnProps.a6fReflnSize[2]) 
      && (0.0 < m_tReflnProps.a6fReflnSize[3]) )
    {
      nSize2    = nint(ABS(m_tReflnProps.a6fReflnSize[2] * 0.5));// / m_fStep0));
      nSize3    = nint(ABS(m_tReflnProps.a6fReflnSize[3] * 0.5));// / m_fStep1));
      nSize4    = nint(ABS(m_tReflnProps.a6fReflnSize[4]));
	bDrawSpot = TRUE;
    }

  if (   (0.0 < m_tReflnProps.a6fReflnSize[4]) 
      && (0.0 < m_tReflnProps.a6fReflnSize[5]) )
    {
      nSize4   = nint(ABS(m_tReflnProps.a6fReflnSize[4] * 0.5 / m_fStep0));
      nSize5   = nint(ABS(m_tReflnProps.a6fReflnSize[5] * 0.5 / m_fStep1));
      bDrawBkg = TRUE;
    }

  if (   (0 < nFI_fEllipseAxisMajor) && (0 < nFI_fEllipseAxisMinor)
      && (0 < nFI_fEllipseAxisMajorOffset) )
    {
      bDrawEllipse = TRUE;
    }

  // Do something to the color

  XGCValues hGCtemp;
	
  hGCtemp.foreground = m_poColormap->ulGetColorIndex(
		                 m_tReflnProps.sReflnColor.string());
  XChangeGC(XtDisplay(m_hParent), m_hGCrefln, GCForeground,
	    &hGCtemp);


  fCalcPx0 = -999.0;
  fCalcPx1 = -999.0;
  int *pnIndex;
  
  // Sort index is guaranteed to exist and be valid, but check anyways and
  // return if it does not exist.

  pnIndex = m_poReflnlist->pnGetSortIndex();
  if (NULL == pnIndex)
    return;

  poRefln  = m_poReflnlist->poGetRefln(pnIndex[0]);
  
  // Get nNonunfFlag of first reflection

  int nNewColor = -1;
  int nReflnSymbol;

  nReflnSymbol = m_tReflnProps.nReflnSymbol;

  if (  bDrawEllipse
	&& (0 > nFI_fCalcRotStart) || (0 > nFI_fCalcRotEnd) )
    {
      nFI_fObsRotMid   = m_poReflnlist->m_nFI_fObsRotMid;
      nFI_fObsRotWidth = m_poReflnlist->m_nFI_fObsRotWidth;
    }
  if (    ((0 > nFI_fCalcRotStart) || (0 > nFI_fCalcRotEnd))
	  &&  ( (4 == nReflnSymbol) || (5 == nReflnSymbol) )
      &&  !bDrawEllipse )
    {
      // There is no information in the reflnlist about the rotation
      // start and/or end, so using circles/squares to plot intersection 
      // with the image is invalid, so use circles
      // Exception: if drawing ellipses!

      nReflnSymbol = 2;

    }

  for (i = 0; i < nNumReflns; i++)
    {
      poRefln  = m_poReflnlist->poGetRefln(pnIndex[i]);
      fObsPx0  = poRefln->fGetField(nFI_fPx0);
      fObsPx1  = poRefln->fGetField(nFI_fPx1);
      if ( (0 <= nFI_fCalcPx0) && (0 <= nFI_fCalcPx1) )
	{
	  fCalcPx0 = poRefln->fGetField(nFI_fCalcPx0);
	  fCalcPx1 = poRefln->fGetField(nFI_fCalcPx1);
	}

      if ( (0 > fObsPx0) || (0 > fObsPx1) )
	{
	  // If the observed pixel list has bogus pixel values,
	  // try to use the calculated pixel values

	  fObsPx0 = fCalcPx0;
	  fObsPx1 = fCalcPx1;
	}

      nStat = nImgPixToDpyPix( fObsPx0, fObsPx1, &nDpyPx0, &nDpyPx1);
      if (0 == nStat)
	{
	  // Draw lines on the WINDOW, not on the PIXMAP!!

	  int nLeft, nTop, nRight, nBottom;

	  int nPoints;
	  XPoint a37hPoints[37];
	  nPoints = 25;

/*	  nLeft   = max(0,nDpyPx0 - nSize0);
	  nTop    = max(0, nDpyPx1 - nSize1);
	  nRight  = min(m_nWidthUsed-1, nDpyPx0 + nSize0);
	  nBottom = min(m_nHeightUsed-1, nDpyPx1 + nSize1);
*/

	  nLeft   = nDpyPx0 - nSize0;
	  nRight  = nDpyPx0 + nSize0;
	  nTop    = nDpyPx1 - nSize1;
	  nBottom = nDpyPx1 + nSize1;

	  if (bDoColors)
	    {
	      // If nNonunfFlag field exists, plot non-zero ones as
	      // a diamond in addition to any other symbol.
	      int nTestColor;
	      nTestColor = poRefln->nGetField(nFI_nNonunfFlag);
	      if ( (NULL != m_poImage) && (1.0 < fMinRawPixOKValue) )
		{
		  // If the underlying pixel value in the image is less
		  // than fMinRawPixOKValue, then this refln is flagged
		  // as potentially ignored.

		  float fValue;
		  fValue = (m_poImage->*m_poImage->prfGetPixel)
		    ((int)fObsPx0, (int)fObsPx1);
		  if (fValue < fMinRawPixOKValue)
		    // Flag as ignored
		    nTestColor = 2;
		}
	      if (nNewColor != nTestColor)
		{
		  // Not same, so need new color.

		  nNewColor = nTestColor;
		  if ( (0 <= nNewColor) && (10 > nNewColor) )
		    {
		      hGCtemp.foreground = m_poColormap->ulGetColorIndex(
					 m_a10sReflnColors[nNewColor]);

		      XChangeGC(XtDisplay(m_hParent), m_hGCrefln, GCForeground,
				&hGCtemp);
		    }
		}
	    }

	  // Select refln symbol to display from properties or
	  // from rotation range

	  nSymbol = nReflnSymbol;
	  if ( (0 < nFI_fCalcRotStart) && (0 < nFI_fCalcRotEnd) )
	    {
	      fRefRotStart = poRefln->fGetField(nFI_fCalcRotStart) + 720.0;
	      while (fRefRotStart >= fMaxAngle) fRefRotStart -= 360.0;
	      fRefRotEnd = poRefln->fGetField(nFI_fCalcRotEnd) + 720.0;
	      while (fRefRotEnd >= fMaxAngle) fRefRotEnd -= 360.0;
	      if (fRefRotEnd < fRefRotStart) 
		{
		  fRefRotEnd += 360.0;
		}
	    }
	  else if (   bDrawEllipse 
		      && (0 < nFI_fObsRotMid) && (0 < nFI_fObsRotWidth) )
	    {
	      fRefRotStart = poRefln->fGetField(nFI_fObsRotMid) 
                           - 0.5 * poRefln->fGetField(nFI_fObsRotWidth) + 720.0;
	      while (fRefRotStart >= fMaxAngle) fRefRotStart -= 360.0;
	      fRefRotEnd = poRefln->fGetField(nFI_fObsRotMid) 
                           + 0.5 * poRefln->fGetField(nFI_fObsRotWidth) + 720.0;
	      while (fRefRotEnd >= fMaxAngle) fRefRotEnd -= 360.0;
	      if (fRefRotEnd < fRefRotStart) 
		{
		  fRefRotEnd += 360.0;
		}
	    }
	  if (4 == nSymbol)
	    {
	      // Refln rot start/end has been transformed to be between
	      // 0 and 360, but the m_tReflnProps.fImageRot* has not been;
	      // but fImageRotStart and fImageRotEnd has been, so use them.

	      //	      if (fRefRotStart > m_tReflnProps.fImageRotEnd)
	      if (fRefRotStart > fImageRotEnd)
		{
		  nSymbol = 3;   // Use box
		}
	      //else if (fRefRotEnd < m_tReflnProps.fImageRotStart)
	      else if (fRefRotEnd < fImageRotStart)
		{
		  nSymbol = 3;   // Use box
		}
	      else
		{
		  // Use circle

		  nSymbol = 2;
		}
	    }
	  else if (5 == nSymbol)
	    {
	      if (fRefRotStart > fImageRotEnd)
		{
		  nSymbol = -1;   // Do not plot
		}
	      else if (fRefRotEnd < fImageRotStart)
		{
		  nSymbol = -1;   // Do not plot
		}
	      else
		{
		  // Use circle

		  nSymbol = 2;
		}
	    }

	  if ( (0 < nNewColor) && (-1 != nSymbol) )
	    {
	      // Draw a triangle
/*
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nLeft, nBottom, nRight, nBottom);
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nRight, nBottom, (nLeft + nRight)/2, nTop);
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			(nLeft + nRight)/2, nTop, nLeft, nBottom);
*/
	      // Draw a diamond

	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			(nLeft + nRight)/2, nBottom, 
			nRight, (nTop+nBottom)/2);
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nRight, (nTop+nBottom)/2, 
			(nLeft + nRight)/2, nTop);
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			(nLeft + nRight)/2, nTop, 
			nLeft,  (nTop+nBottom)/2);
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nLeft,  (nTop+nBottom)/2,
			(nLeft + nRight)/2, nBottom);
	    }
	  int nJim;
	  nJim = nSymbol;
	  switch(nSymbol)
	    {
	    case 0: // plus symbol
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nDpyPx0, nTop, nDpyPx0, nBottom);
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nLeft, nDpyPx1, nRight, nDpyPx1);
	      break;
	    case 1: // cross symbol
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nLeft, nTop, nRight, nBottom);
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			nRight, nTop, nLeft, nBottom);
	      break;
	    case 2: // circle symbol
	      XDrawArc(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
		       nLeft, nTop, (nRight-nLeft), (nBottom-nTop), 0, 23040);
	      break;
	    case 3: // box symbol
	      XDrawRectangle(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			     nLeft, nTop, (nRight-nLeft), (nBottom-nTop));
	      break;
	    case 6: // Draw the tilted ellipsoid
	      // Note that the tilt angle is negated here
	      nStat = nCalcEllipse(fObsPx0, fObsPx1, float(nSize0), 
				   float(nSize1), -float(nSize2),
			   nPoints, a37hPoints);
	      if (0 == nStat)
		XDrawLines(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			   a37hPoints, nPoints, CoordModeOrigin);
	      break;
	    case -1:
	      break;  // Do not draw a symbol (could be ignored refln)
	    default:
	      cerr << "Error: Unknown symbol\n";
	      break;
	    }
	  nSymbol = nJim;

	  // Draw line from observed pixel position to calculate pixel position
	  //  but only if they differ

	  if (    (0 <= nSymbol) && (m_bLineInRefln)

	      && ( (fObsPx0 != fCalcPx0) || (fObsPx1 != fCalcPx1) ) )
	    {
	      nStat = nImgPixToDpyPix( fCalcPx0, fCalcPx1, 
				      &nDpyCalcPx0, &nDpyCalcPx1);	  
	      if ( (0 == nStat) && (0 < nDpyCalcPx0) )
		{
		  XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			    nDpyPx0, nDpyPx1, nDpyCalcPx0, nDpyCalcPx1);
		}
	    }
	  if (0 <= nSymbol)
	    {
	      if (bDrawSpot)
		{
		  // Note that the tilt angle is negated here

		  nStat = nCalcEllipse(fObsPx0, fObsPx1, float(nSize2), 
				       float(nSize3), -float(nSize4),
				       nPoints, a37hPoints);
		  //if (0 == nStat)
		    XDrawLines(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			       a37hPoints, nPoints, CoordModeOrigin);
	    /************
		  nLeft   = nDpyPx0 - nSize2;
		  nRight  = nDpyPx0 + nSize2;
		  nTop    = nDpyPx1 - nSize3;
		  nBottom = nDpyPx1 + nSize3;
		  XDrawArc(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			   nLeft, nTop, (nRight-nLeft), (nBottom-nTop),
			   0, 23040);
   **************/
		}
	      if (bDrawBkg)
		{
		  nLeft   = nDpyPx0 - nSize4;
		  nRight  = nDpyPx0 + nSize4;
		  nTop    = nDpyPx1 - nSize5;
		  nBottom = nDpyPx1 + nSize5;
		  XDrawArc(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			   nLeft, nTop, (nRight-nLeft), (nBottom-nTop),
			   0, 23040);
		}
	      if (bDrawEllipse)
		{
		  fAxisMajor       = poRefln->fGetField(nFI_fEllipseAxisMajor);
		  fAxisMinor       = poRefln->fGetField(nFI_fEllipseAxisMinor);
		  fAxisMajorOffset = poRefln->fGetField(nFI_fEllipseAxisMajorOffset);

		  nStat = nCalcEllipse(fObsPx0, fObsPx1, fAxisMajor, fAxisMinor, 
				       fAxisMajorOffset, nPoints, a37hPoints);
		  XDrawLines(XtDisplay(m_hParent), m_hPixmap, m_hGCrefln,
			     a37hPoints, nPoints, CoordModeOrigin);
		}
	    }
	}
    }
  //  m_poCursor->vRestore();
  m_poCursor->vReset();
}

void
CXdisplay::vEventCallbackCallback(Widget w, 
		XtPointer clientData, XtPointer callData)
{
    UICallbackStruct *data = (UICallbackStruct *) clientData;
    CXdisplay *obj = (CXdisplay *)data->object;
    
    obj->vEventCallback(w, data->client_data, callData);
}

void
CXdisplay::vPlotPixelValues(void)
{

  // Plot pixel values in the image, if the pixels are zoomed big enough

  int   i, j;
  int   nPx0, nPx1;
  int   nDpyPx0, nDpyPx1;
  float fValue;
  int   nStat, nLen;
  char  cText[20];
  int   nVertOffset;
  int   nTextHeight;
  int   nTextWidth;

  if (NULL == m_poImage) return;

  if ( (0.02 >= ABS(m_fStep0)) && (0.02 >= ABS(m_fStep1)) )
    {
//      nVertOffset = m_pXFont->ascent + m_pXFont->descent + (int) (0.5/m_fStep1);
      nVertOffset = (int) ABS(0.5 / m_fStep1);
      nTextHeight = m_pXFont->ascent + m_pXFont->descent;
      nPx1 = m_nStart1;
      for (j = 0; j <= (m_tImageProps.nExt[1]+1); j++) 
	{
	  nPx0 = m_nStart0;
	  for (i = 0; i <= (m_tImageProps.nExt[0]+1); i++) 
	    {
	      if (4 > m_tImageProps.nOrient)
		{
		  fValue = (m_poImage->*m_poImage->prfGetPixel)(nPx0, nPx1);
		  nStat = nImgPixToDpyPix((float)nPx0+0.5, (float)nPx1+0.5,
					  &nDpyPx0, &nDpyPx1);
		}
	      else
		{
		  fValue = (m_poImage->*m_poImage->prfGetPixel)(nPx1, nPx0);
		  nStat = nImgPixToDpyPix((float)nPx1+0.5, (float)nPx0+0.5,
					  &nDpyPx0, &nDpyPx1);
		}
//+JWP 2009-01-08
	      if (!m_bHasFloatValues)
		nLen = sprintf(cText, "%.0f", fValue);
	      else if (fValue >= 10.0f)
		nLen = sprintf(cText, "%.0f", fValue);
	      else
		nLen = sprintf(cText, "%.4f", fValue);
//-JWP 2009-01-08
	      if (0 == nStat)
		{
		  // Try to draw centered text

		  nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
		  XDrawString(XtDisplay(m_hParent), m_hPixmap, m_hGCtext,
			      nDpyPx0-nTextWidth/2, nDpyPx1+nTextHeight/2,
			      cText, nLen);
		}
	      nPx0 += m_nDir0;
	    }
	  nPx1 += m_nDir1;
	}
    }
}

void
CXdisplay::vSetSatPixColor(const char *pcColor)
{
  // Change the saturated pixel color, but only if the color is already in the
  // colormap

  unsigned long int ulStat;

  ulStat = m_poColormap->ulGetColorIndex(pcColor);
  if (ulStat != CXcolormap::ms_ulCOLORERROR)
    {
      m_ulSatColorIndex = ulStat;
      m_tImageProps.sSatPixColor = pcColor;
    }
}

void
CXdisplay::vSetTooLowPixColor(const char *pcColor)
{
  // Change the TooLow pixel color, but only if the color is already in the
  // colormap

  unsigned long int ulStat;

  ulStat = m_poColormap->ulGetColorIndex(pcColor);
  if (ulStat != CXcolormap::ms_ulCOLORERROR)
    {
      m_ulTooLowColorIndex = ulStat;
      m_tImageProps.sTooLowPixColor = pcColor;
    }
}

int
CXdisplay::nCreatePS(const Cstring& sPSFile, const Cstring& sComment,
		     const Boolean& bPrintName)
{
  int i, j;
  int nStat = 1;
  Cstring sTemp;
  XImage *phImage;
  unsigned long ulPlaneMask = 0xFFFFFF;  // Get ALL bit planes

  float fAspect;
  int   nWidth   = 450;
  int   nHeight  = 450;

  m_poCursor->vSetWait();

  // Get the contents of the pixmap which has all things drawn on it

  phImage = XGetImage(XtDisplay(m_hParent),
		      m_hPixmap,
		      0, 0, m_nWidthUsed, m_nHeightUsed,
		      ulPlaneMask, ZPixmap);

  if (NULL == phImage)
    {
      //      m_poCursor->vRestore();
      m_poCursor->vReset();
      return (2);
    }

  fAspect = (float) m_nHeightUsed / (float) m_nWidthUsed;

  if (1.0 < fAspect)
    {
      // Height greater than width, reduce width

      nWidth = (int) ( (float) nWidth / fAspect);
    }
  else
    {
      // Width greater than height, reduce height
      nHeight = (int) ((float) nHeight * fAspect);
    }
  

  ofstream oOut(sPSFile.string());
  if (oOut.rdbuf()->is_open()) 
    {
      oOut << "%!\n"
	   << "newpath\n";

	  // Write filename if specified.
      if (bPrintName)
	{
	  sTemp = sPSFile;
	}
      else
	{
	  sTemp = m_poImage->sGetName();
	}
      oOut << "/Helvetica findfont 20 scalefont setfont\n"
	   << "80 210 moveto\n"
	   << "(" << sTemp << ") show ";

      oOut << "80 180 moveto\n"
	   << "/Helvetica findfont 12 scalefont setfont\n";

    // Write comment

      sTemp = sComment;
      i = 0;
      while (0 <= sTemp.find('\n'))
	{
	  i++;                           // Count newlines in sTemp
	  oOut << "(" << sTemp.before('\n') << ") show"
	       << " 80 180 -15 " << i << " mul add moveto\n";
	  sTemp = sTemp.after('\n');
	}
      oOut << "(" << sTemp << ") show\n";
      oOut << "80 270 moveto "
    // Draw box
           << nWidth << " 0 rlineto 0 "
           << nHeight << " rlineto " << -nWidth << " 0 rlineto 0 "
           << -nHeight << " rlineto stroke\n"
           << "0 0 moveto save 80 270 translate gsave " << nWidth << " "
           << nHeight << " scale\n"
           << "/picstr " << m_nWidthUsed << " string def\n"
           << m_nWidthUsed << " " << m_nHeightUsed << " 8\n"
           << "[ " << m_nWidthUsed << " 0 0 " << -m_nHeightUsed 
	   << " 0 " << m_nHeightUsed << " ]\n"
           << "{ currentfile picstr readhexstring pop }\nimage\n";

      // This needs better encapsulation of the colormap.

      if (8 == phImage->depth)
	{
	  // This is faster than using XGetPixel() when it's an 8-bit image

	  // The member variable m_pucData contains the indices into the
	  // the colormap for pixels in the displayed pixmap.  The indices are NOT
	  // the colors!  So we need a way to get the colors from the indices.
      
	  unsigned char *pucTemp;
	  int           nTemp;
	  char          a3cTemp[3];
	  pucTemp = (unsigned char *) phImage->data; 
	  int a256nIndex[256];
	  for (i = 0; i < 256; i++)
	    {
	      a256nIndex[i] = 255;  // All white to start with
	    }

	  // Force a XQueryColors (trick!):
	  //   you have to know the source code to CXcolormap to know this
	  //   call XQueryColors and fills in m_poColormap->m_phXcolor.

	  (void) m_poColormap->ulGetColorIndex(0, 0, 0, 0);

	  for (i = 0; i < m_poColormap->m_nAllocated; i++)
	    {
	      // We need to reverse engineer the color lookup table indexing
	      // Base B&W on the color red
	  
	      nTemp = m_poColormap->m_phXcolor[i].red / 256;// 65356 / 256 = 256
	      if (0   > nTemp) 
		{
//	          cerr << "Warning! 0 > nTemp in CreatePS!\n";
		  nTemp = 0;
		}
	      if (255 < nTemp)
		{
//	      cerr << "Warning! 255 < nTemp in CreatePS!\n";
		  nTemp = 255;
		}
	      a256nIndex[m_poColormap->m_phXcolor[i].pixel] = nTemp;
	    }

	  for (j = 0; j < m_nHeightUsed; j++)
	    {
	      for (i = 0; i < m_nWidthUsed; i++)
		{
		  // 0 is black; 255 is white

		  nTemp = a256nIndex[*pucTemp++];

		  sprintf(a3cTemp, "%02X", nTemp);
		  oOut << a3cTemp;
		  if (39 == i % 40)
		    oOut << "\n";
		}

	      // Skip extra unused bytes in *phImage at the end of each line

	      for (i = m_nWidthUsed; i < phImage->bytes_per_line; i++)
		pucTemp++;
	      oOut << "\n";
	    }
	}
      else
	{
	  // For images with a depth of more than 8 bits, we need
	  // to use a slower method to get the colors:  XGetPixel & XQueryColor

	  unsigned long hPixel;
	  XColor        hXcolor;

	  // Need to convert from Xserver representation of pixels to a single
	  // character or byte 0 to 255 for the grayscale

	  int           nTemp;
	  char          a3cTemp[3];

	  hXcolor.flags = DoRed | DoBlue | DoGreen;

	  for (j = 0; j < m_nHeightUsed; j++)
	    {
	      for (i = 0; i < m_nWidthUsed; i++)
		{
		  // In PostScript 0 is black; 255 is white

		  hPixel = XGetPixel(phImage, i, j);
	      
		  // Get the RGB values for the hPixel

		  hXcolor.pixel = hPixel;
		  //		  cout << hPixel << endl;
		  XQueryColor(XtDisplay(m_hParent), m_poColormap->m_hColormap,
			      &hXcolor);

		  // Use the red value to map to PS gray

		  nTemp = (int)hXcolor.red / 256;
		  if (0 > nTemp) 
		    {
//		      cout << "Warning! 0 > nTemp in CreatePS!\n";
		      nTemp = 0;
		    }
		  else if (255 < nTemp) 
		    {
//		      cout << "Warning! 255 < nTemp in CreatePS!\n";
		      nTemp = 255;
		    }
		  sprintf(a3cTemp, "%02X", nTemp);
		  oOut << a3cTemp;
		  if (39 == i % 40)
		    oOut << "\n";
		}
	      oOut << "\n";
	    }
	}
      oOut << "grestore 1 setlinewidth restore showpage\n";
      oOut.close();
      nStat = 0;
    }
  if (NULL != phImage) XDestroyImage(phImage);
  //  m_poCursor->vRestore();
  m_poCursor->vReset();
  return (nStat);
}

void
CXdisplay::vDrawProfile(const int nXIn, const int nYIn, const int nInfoMode)
{
  // Draw a profile along right-hand and bottom-edge of window.
  // nInfoMode = 0; lines are parallel to window edges
  // nInfoMode = 1; line goes from previous pixel to current pixel (not impl)

  int   i;
  int   nStat;
  int   nOffset, nExt0, nExt1;
  int   nX, nY;
  float fPx0, fPx1;
  int   nPx0, nPx1, nPx0Prev, nPx1Prev;
  float fValue, fMin, fMax, fScale;
  int   nLen, nTextWidth, nTextHeight;
  char  cText[20];
  
  if (0 == nInfoMode)
    {
      // Clear areas for plots.  Do this first, so any errors leave a display
      // cleared of the profile

      // Try to get profiles non-gray, fill rectangle with white since
      // it seems that XClearArea does not alway work correctly.
      // Use m_hGCrefln for foreground since it's foreground is always
      // reset in vPlotReflnlist().

      XGCValues hGCtemp;
	
      hGCtemp.foreground = m_ulColorWhite;
      XChangeGC(XtDisplay(m_hParent), m_hGCrefln, GCForeground,
		&hGCtemp);

      XFillRectangle(XtDisplay(m_hParent), XtWindow(m_hParent),
		     m_hGCrefln, m_nWidthUsed+1, 0, 
		     m_nWidth-m_nWidthUsed, m_nHeightUsed+1);

      XFillRectangle(XtDisplay(m_hParent), XtWindow(m_hParent),
		     m_hGCrefln, 0, m_nHeightUsed+1,
		     m_nWidthUsed+1, m_nHeight-m_nHeightUsed);

      // This is corner which may/should be cleared by another routine:

      XFillRectangle(XtDisplay(m_hParent), XtWindow(m_hParent),
		     m_hGCrefln, m_nWidthUsed+1, m_nHeightUsed+1,
//		     m_nHeight-m_nHeightUsed, m_nWidth-m_nWidthUsed);
		     m_nWidth-m_nWidthUsed, m_nHeight-m_nHeightUsed);
      // Enclose image area in box

      hGCtemp.foreground = m_ulColorBlue;
      XChangeGC(XtDisplay(m_hParent), m_hGCrefln, GCForeground,
		&hGCtemp);

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		0, m_nHeightUsed, m_nWidth, m_nHeightUsed);
  
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		m_nWidthUsed, 0, m_nWidthUsed, m_nHeight);
      
      if (NULL == m_poImage) return;

      if ( (nXIn >= m_nWidthUsed) || (nYIn >= m_nHeightUsed) )
	return;  // Out-of-bounds

      nX = nXIn;
      nY = nYIn;

      if (-1 == nX) nX = m_nXPrev;
      if (-1 == nY) nY = m_nYPrev;

      nStat  = nDpyPixToImgPix( nX, nY, &fPx0, &fPx1);

      if (0 != nStat)
	return; // Out-of-bounds

      nTextHeight = m_pXFont->ascent + m_pXFont->descent;

      fValue = (m_poImage->*m_poImage->prfGetPixel)((int)fPx0, (int)fPx1);
      fMin   = fValue;
      fMax   = fValue;

      nExt0  = m_nWidth  - m_nWidthUsed;
      nExt1  = m_nHeight - m_nHeightUsed;
  
      // Draw extension of cross-hairs through profile plot

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		m_nWidthUsed, nY, m_nWidth,  nY);

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		nX, m_nHeightUsed, nX, m_nHeight);

      if (20 < nExt0)
	{
	  // There is room for a plot

	  // Find min, max along line

	  for (i = 0; i < m_nHeightUsed; i++)           // Loop top to bottom
	    {
	      nStat = nDpyPixToImgPix(nX, i, &fPx0, &fPx1);
	      fValue = (m_poImage->*m_poImage->prfGetPixel)((int)fPx0, (int)fPx1);
	      if (fValue < fMin) fMin = fValue;
	      if (fValue > fMax) fMax = fValue;
	    }

	  if (fMax == fMin)
	    {
	      fMax = fMin + 1.0;
	    }

	  // Scale fValue to fit between m_nWidth and m_nWidthUsed and place
	  // in nPx0
      
	  fScale   =  0.9 * (float) nExt0 / (fMax - fMin);
	  nOffset  = (int) (0.05 * (float) nExt0);
	  nPx0Prev = m_nWidth - (int) ((fValue - fMin) * fScale + nOffset);
	  nPx1Prev = 0;
	  for (i = 0; i < m_nHeightUsed; i++)           // Loop top to bottom
	    {
	      nStat = nDpyPixToImgPix( nX, i, &fPx0, &fPx1);
	      fValue = (m_poImage->*m_poImage->prfGetPixel)((int)fPx0, (int)fPx1);
	      
	      nPx0 = m_nWidth - (int) ((fValue - fMin) * fScale + nOffset);
	      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			nPx0Prev, nPx1Prev, nPx0, i);
	      nPx0Prev = nPx0;
	      nPx1Prev = i;
	    }

	  // Draw min and max values on the plot if there is room

	  nLen        = sprintf(cText, "%.0f", fMax);
	  nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
	  nPx0        = m_nWidthUsed+1;
	  nPx1        = nTextHeight;
	  XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		      nPx0, nPx1, cText, nLen);
	  nLen        = sprintf(cText, "%.0f", fMin);
	  nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
	  nPx0        = m_nWidth - nTextWidth;
	  nPx1        = nTextHeight + nTextHeight;
	  XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		      nPx0, nPx1, cText, nLen);
	}

      if (20 < nExt1)
	{
	  // There is room for a plot

	  // Find min, max along line

	  for (i = 0; i < m_nWidthUsed; i++)           // Loop left to right
	    {
	      nStat = nDpyPixToImgPix(i, nY, &fPx0, &fPx1);
	      fValue = (m_poImage->*m_poImage->prfGetPixel)((int)fPx0, (int)fPx1);
	      if (fValue < fMin) fMin = fValue;
	      if (fValue > fMax) fMax = fValue;
	    }

	  if (fMax == fMin)
	    {
	      fMax = fMin + 1.0;
	    }
	  
	  // Scale fValue to fit between m_nWidth and m_nWidthUsed and place
	  // in nPx0
      
	  fScale   =  0.9 * (float) nExt1 / (fMax - fMin);
	  nOffset  = (int) (0.05 * (float) nExt1);
	  nPx1Prev = m_nHeight - (int) ((fValue - fMin) * fScale + nOffset);
	  nPx0Prev = 0;
	  for (i = 0; i < m_nWidthUsed; i++)           // Loop left to right
	    {
	      nStat = nDpyPixToImgPix( i, nY, &fPx0, &fPx1);
	      fValue = (m_poImage->*m_poImage->prfGetPixel)((int)fPx0, (int)fPx1);
	      
	      nPx1 = m_nHeight - (int) ((fValue - fMin) * fScale + nOffset );
	      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			nPx0Prev, nPx1Prev, i, nPx1);
	      nPx0Prev = i;
	      nPx1Prev = nPx1;
	    }

	  if ( (6 * nTextHeight) < nExt1)
	    {
	      // Draw min and max values on the plot if there is room

	      nLen        = sprintf(cText, "%.0f", fMin);
	      nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
	      nPx0        = 0;
	      nPx1        = m_nHeight;
	      XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			  nPx0, nPx1, cText, nLen);
	      nLen        = sprintf(cText, "%.0f", fMax);
	      nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
	      nPx0        = 0;
	      nPx1        = m_nHeightUsed + nTextHeight;
	      XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			  nPx0, nPx1, cText, nLen);
	    }
	}

      // Plot rocking curve if desired in lower right corner

      if (9 == m_tImageProps.nPlotRockingCurve)
	{
	  // For now, support only a 3x3 tiled image

	  int   a2nOffsets[2][9];
	  int   nDim0, nDim1;
	  int   nOff0, nOff1;
	  float a9fValue[9];

	  (void) m_poImage->nGetDimensions(&nDim0, &nDim1);
	  nOff0 = nDim0 / 3;
	  nOff1 = nDim1 / 3;

	  a2nOffsets[0][0] = -nOff0;
	  a2nOffsets[1][0] = -nOff1;
	  a2nOffsets[0][1] = 0;
	  a2nOffsets[1][1] = -nOff1;
	  a2nOffsets[0][2] = +nOff0;
	  a2nOffsets[1][2] = -nOff1;
	  a2nOffsets[0][3] = -nOff0;
	  a2nOffsets[1][3] = 0;
	  a2nOffsets[0][4] = 0;
	  a2nOffsets[1][4] = 0;
	  a2nOffsets[0][5] = +nOff0;
	  a2nOffsets[1][5] = 0;
	  a2nOffsets[0][6] = -nOff0;
	  a2nOffsets[1][6] = +nOff1;
	  a2nOffsets[0][7] = 0;
	  a2nOffsets[1][7] = +nOff1;
	  a2nOffsets[0][8] = +nOff0;
	  a2nOffsets[1][8] = +nOff1;

	  nX = nXIn;
	  nY = nYIn;

	  if (-1 == nX) nX = m_nXPrev;
	  if (-1 == nY) nY = m_nYPrev;

	  nStat  = nDpyPixToImgPix( nX, nY, &fPx0, &fPx1);	  
	  nX = (int)fPx0;
	  nY = (int)fPx1;

	  // Move to pixel in middle panel of image

	  if (nX <   nOff0) nX += nOff0;
	  if (nX > 2*nOff0) nX -= nOff0;
	  if (nY <   nOff1) nY += nOff1;
	  if (nY > 2*nOff1) nY -= nOff1;

	  fValue = (m_poImage->*m_poImage->prfGetPixel)(nX, nY);
	  fMin   = fValue;
	  fMax   = fValue;
	  for (i = 0; i < 9; i++)
	    {
	      nOff0 = nX + a2nOffsets[0][i];
	      nOff1 = nY + a2nOffsets[1][i];
	      fValue = (m_poImage->*m_poImage->prfGetPixel)(nOff0, nOff1);

	      a9fValue[i] = fValue;
	      if (fValue < fMin) fMin = fValue;
	      if (fValue > fMax) fMax = fValue;
	    }

	  if (fMax == fMin)
	    {
	      fMax = fMin + 1.0;
	    }
	  nExt0  = m_nWidth  - m_nWidthUsed;
	  nExt1  = m_nHeight - m_nHeightUsed;
	  fScale   =  0.9 * (float) nExt1 / (fMax - fMin);

	  nOffset  = (int) (0.05 * (float) nExt1);
	  fValue   = a9fValue[0];
	  nPx1Prev = m_nHeight - (int) ((fValue - fMin) * fScale + nOffset);
	  nPx0Prev = m_nWidthUsed;
	  for (i = m_nWidthUsed; i < m_nWidth; i++)
	    {
	      nPx0 = (i - m_nWidthUsed) / (nExt0 / 9);
	      if (0 > nPx0)  nPx0 = 0;
	      if (9 <= nPx0) nPx0 = 8;
	      fValue = a9fValue[nPx0];	      

	      nPx1 = m_nHeight - (int) ((fValue - fMin) * fScale + nOffset );

	      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			nPx0Prev, nPx1Prev, i, nPx1);
	      nPx0Prev = i;
	      nPx1Prev = nPx1;
	    }

	  if ( (6 * nTextHeight) < nExt1)
	    {
	      // Draw min and max values on the plot if there is room

	      nLen        = sprintf(cText, "%.0f", fMin);
	      nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
	      nPx0        = m_nWidth - nTextWidth;
	      nPx1        = m_nHeight;
	      XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			  nPx0, nPx1, cText, nLen);

	      nLen        = sprintf(cText, "%.0f", fMax);
	      nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
	      nPx0        = m_nWidth - nTextWidth;
	      nPx1        = m_nHeightUsed + nTextHeight;
	      XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			  nPx0, nPx1, cText, nLen);

	      nLen        = sprintf(cText, "Rock.", fMin);
	      nTextWidth  = XTextWidth(m_pXFont, cText, nLen);
	      nPx0        = (m_nWidthUsed + m_nWidth - nTextWidth)/2 ;
	      nPx1        = m_nHeight;
	      XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
			  nPx0, nPx1, cText, nLen);
	    }
	}
    }

// What are these doing here? //////////////////////////////////////////////////
  m_fScaleMin     = m_tImageProps.fScaleMin;
  m_fScaleMax     = m_tImageProps.fScaleMax;
////////////////////////////////////////////////////////////////////////////////
}

void
CXdisplay::vPan(const int nX, const int nY)
{
  // Pan the image by some small amount

  int nDelX, nDelY;
  float fPx0, fPx1;

  if (   (0 <= nX)
      && (0 <= nY)
      && (nX < m_nWidthUsed)
      && (nY < m_nHeightUsed))
    {
      // Still within the image area

      if (0 == nDpyPixToImgPix(nX, nY, &fPx0, &fPx1))
	{
	  nDelX = (int) (m_fPx0Prev - fPx0) + m_tImageProps.nOrig[0];
	  nDelY = (int) (m_fPx1Prev - fPx1) + m_tImageProps.nOrig[1];

	  int nStatSR = 0;
	  nStatSR == nSetRegion(nDelX, nDelY, m_tImageProps.nExt[0],
				m_tImageProps.nExt[1]);
	    if (0 == nStatSR)
	    {
	      vRefresh();
	    }
	    //cout << "::vPan nStat, nDelX, nDelY: " << nStatSR << ", " << nDelX << ", " << nDelY << endl;
	}
    }
}

void
CXdisplay::vDrawColorScale(const int nMode, const int nX, const int nY)
{
  // Draw a pixmap in a part of the drawing area with an image of the color
  // scale.  Allow control of the color scale with the mouse.
  // Mode = 0  -- Initial draw
  // Mode = 1  -- Refresh, Expose
  // Mode = 2  -- from vPixelInfo, ButtonPress or ButtonMotion probably unless
  // Mode = 3  -- from vPixelInfo, ButtonRelease within colorscale graph

  // The color scale is implemented as two linear scales with separate slopes.
  // The left-hand side goes from fMin to (fAvg + 5 SD)
  // The right-hand side goes from (fAvg +5SD) to fMax.
  // (fAvg + 5SD)  is mapped to 4/5ths of the total scalewidth
  // It could have been exponential, but this was thought to be faster??!

  bool bMinChanged = FALSE;
  bool bMaxChanged = FALSE;
  if (1.0 <= m_tImageProps.fZoom) return;
  if (3 == nMode)
    {
      // Release of button1, so new scaling should be applied to main image area

      m_tImageProps.fScaleMin = m_fScaleMin;
      m_tImageProps.fScaleMax = m_fScaleMax;
      m_nScaleMode = 0;
      vRefresh();
      return; 
    }

  // Create image for color scale graphic

  m_nScaleWidth  = m_nWidthUsed;
  m_nScaleHeight = min(15, (m_nHeight - m_nHeightUsed) / 2);

  if ( (0 >= m_nScaleWidth) || (0 >= m_nScaleHeight) )
      return;                    //No need to redraw


  if (-1 < nMode)
    {
      // Create an XImage if necessary

      if (NULL == m_phXImageScale)
	{
	  m_phXImageScale = XCreateImage(XtDisplay(m_hParent),
					 DefaultVisual(XtDisplay(m_hParent),
                                           DefaultScreen(XtDisplay(m_hParent))),
					 m_ulDepthImage, ZPixmap, 0, 
					 (char *)NULL,
					 m_nScaleWidth, m_nScaleHeight, 8, 0);

	}
      else if (NULL == m_phXImageScale->data)
	{
	  // Image existed, but data was NULL, so probably wrong size

	  m_phXImageScale = XCreateImage(XtDisplay(m_hParent),
					 DefaultVisual(XtDisplay(m_hParent),
                                           DefaultScreen(XtDisplay(m_hParent))),
					 m_ulDepthImage, ZPixmap, 0, 
					 (char *)NULL,
					 m_nScaleWidth, m_nScaleHeight, 
					 XBitmapPad(XtDisplay(m_hParent)), 0);
	}

      int nSizeNeeded = m_phXImageScale->bytes_per_line * m_phXImageScale->height;

      // Check if enough memory was allocated already for the items needed

      if ( (NULL != m_pucScale) && (m_nScaleSize < nSizeNeeded) )
	{
	  delete [] m_pucScale;
	  m_pucScale = NULL;
	}
      if (NULL == m_pucScale)
	{
	  // Allocate the memory

	  m_pucScale = new unsigned char [nSizeNeeded];
	  m_nScaleSize = nSizeNeeded;
	}
      m_phXImageScale->data = (char *)m_pucScale;
    }

  int   i, j;
  float fMin, fMax, fAvg, fSD, fAvg5SD;
  fMin = m_tImageProps.fImageMin;
  fMax = m_tImageProps.fImageMax;
  fAvg = m_tImageProps.fAvg;
  fSD  = m_tImageProps.fSD;
  fAvg5SD = min(fMax-1, (fAvg +  m_tImageProps.fScaleMaxSD * fSD));

  if (fMax == fMin) fMax += 1.0;     // Prevent potential divide by zero later
  
  int nIndexMin = 0;
  int nIndexMax = m_poColormap->nNumColors();
  
  int nSlopeChange = 4 * m_nScaleWidth / 5;
  if (fAvg5SD == (fMax-1))
    nSlopeChange = (m_nScaleWidth-1);

  float fScaleL  = (fAvg5SD - fMin) / (float) nSlopeChange;
  float fScaleR  = (fMax - fAvg5SD) / (float) (m_nScaleWidth-nSlopeChange);

  float fTemp;
  int   nTemp;

  // Figure out which pixels the values of m_fScaleMin and m_fScaleMax are
  // set to.

  int   nScaleMinPix, nScaleMaxPix;

  if (m_fScaleMin < fAvg5SD)
    {
      nScaleMinPix = (int) ((m_fScaleMin - fMin) / fScaleL);
    }
  else
    {
      nScaleMinPix = (int) ((m_fScaleMin - fAvg5SD) / fScaleR) + nSlopeChange;
    }
  nScaleMinPix = max(0, nScaleMinPix);
  nScaleMinPix = min(m_nScaleWidth, nScaleMinPix);
    
  if (m_fScaleMax < fAvg5SD)
    {
      nScaleMaxPix = (int) ((m_fScaleMax - fMin) / fScaleL);
    }
  else
    {
      nScaleMaxPix = (int) ((m_fScaleMax - fAvg5SD) / fScaleR + nSlopeChange);
    }
  nScaleMaxPix = min(m_nScaleWidth, nScaleMaxPix);
  nScaleMaxPix = max(0, nScaleMaxPix);  

  nTemp = m_nWidthUsed / 20;
  if (-1 < nMode)
    {
      unsigned char *pucTemp;
      unsigned char *pucTemp2;
      unsigned short int *puiTemp;
      unsigned short int *puiTemp2;
      float          fScale1;
      int nDel[2];
      nDel[0] = ABS(nX-nScaleMinPix);    // Is cursor close to Scale min or max?
      nDel[1] = ABS(nX-nScaleMaxPix);
      i = 0;
      if (nDel[0] > nDel[1]) i = 1;      // Cursor is closer to scalemax
//      if (nTemp > nDel[i])               // Is it close enough?
      if ( (2 == nMode) || (nTemp > nDel[i]) )    // Is it close enough?
	{
	  if (0 == i)
	    {
	      // Cursor is closer to scaleMin, so try to drag it

	      bMinChanged = TRUE;
	      if (nX < nSlopeChange)
		m_fScaleMin = min(m_fScaleMax, (float)nX * fScaleL + fMin);
	      else
		m_fScaleMin = min(m_fScaleMax, (float)(nX-nSlopeChange)
				                * fScaleR + fAvg5SD);
	    }
	  else
	    {
	      // Cursor is closer to scaleMax, so try to drag it

	      bMaxChanged = TRUE;
	      if (nX < nSlopeChange)
		m_fScaleMax = max(m_fScaleMin, (float)nX * fScaleL + fMin);
	      else
		{
		  m_fScaleMax = max(m_fScaleMin, (float)(nX-nSlopeChange)
				    * fScaleR + fAvg5SD);
		}
	    }
	}

      // Fill m_pucScale with a color

      unsigned long int ulColor;
      pucTemp = m_pucScale;
      if (2 == m_nBytesPerPixel)
	{
	  // WARNING! This might have alignment problems on OSF1

	  puiTemp = (unsigned short int*) m_pucScale;
	}

      if (0.001f > (m_fScaleMax - m_fScaleMin))
	{
	  m_fScaleMax = m_fScaleMin + 1.0f;
	}

      fScale1 = (float) (nIndexMax - nIndexMin) / (m_fScaleMax - m_fScaleMin);
      for (i = 0; i < m_nScaleWidth; i++)
	{
	  if (i < nSlopeChange)
	    {
	      // On left hand side

	      fTemp = ((float)i * fScaleL) + fMin;  // Scale k to image value
	    }
	  else
	    {
	      // On right hand side
      	      fTemp = ((float)(i - nSlopeChange) * fScaleR) + fAvg5SD;
	    }
	  nTemp = nIndexMin;
	  if ( (fTemp >= m_tImageProps.fSatPixValue)
	      && (CXcolormap::ms_ulCOLORERROR != m_ulSatColorIndex) )
	    {
	      ulColor = m_ulSatColorIndex;
	    }
	  else if ( (fTemp <= m_tImageProps.fTooLowPixValue)
	      && (CXcolormap::ms_ulCOLORERROR != m_ulSatColorIndex) )
	    {
	      ulColor = m_ulTooLowColorIndex;
	    }
	  else
	    {
	      if (fTemp >= m_fScaleMax)
		nTemp = nIndexMax-1;
	      else if (fTemp > m_fScaleMin)
		{
		  // Scale image value to color index
		  nTemp = nTemp + (int) ((fTemp-m_fScaleMin) * fScale1);
		}
	      ulColor = m_poColormap->ulColorIndex(nTemp);
	    }

	  // In the following, no need to worry about XBitmapPad because
	  // the height is in the inner loop 

	  switch(m_nBytesPerPixel)
	    {
	    case 1:
	      pucTemp2 = pucTemp;
	      for (j = 0; j < m_nScaleHeight; j++)
		{
		  *pucTemp2 = (unsigned char)ulColor;
		  pucTemp2 += m_nScaleWidth;
		}
	      pucTemp++;
	      break;
	    case 2:
	      puiTemp2 = puiTemp;
	      for (j = 0; j < m_nScaleHeight; j++)
		{
		  *puiTemp2 = (unsigned short int)ulColor;
		  puiTemp2 += m_nScaleWidth;
		}
	      puiTemp++;
	      break;
	    default:
	      // Must use the relatively slow XPutPixel
	      for (j = 0; j < m_nScaleHeight; j++)
		{
		  XPutPixel(m_phXImageScale, i, j, ulColor);
		}
	    }
	}
    }

  // Put the image directly on the window and do not use a pixmap.
  // Put only the part of the image that we think has changed.
  
  int nStart, nEnd;
  nStart = 0;
  nEnd   = m_nScaleWidth;
  if ( (2 > nMode) || (0 == m_nScaleMode) )
    {
      // On refresh, expose or first time, redraw entire colorscale

      bMinChanged = TRUE;
      bMaxChanged = TRUE;
    }
  if (!bMinChanged)
    {
      nStart = nScaleMinPix;
    }
  else
    {
      if (m_fScaleMin < fAvg5SD)
	{
	  nScaleMinPix = (int) ((m_fScaleMin - fMin) / fScaleL);
	}
      else
	{
	  nScaleMinPix = (int) ((m_fScaleMin - fAvg5SD) / fScaleR  + nSlopeChange);
	}
      nScaleMinPix = max(0, nScaleMinPix);
      nScaleMinPix = min(m_nScaleWidth, nScaleMinPix);
    }
    
  if (!bMaxChanged)
    {
      nEnd   = nScaleMaxPix;
    }
  else
    {
      if (m_fScaleMax < fAvg5SD)
	{
	  nScaleMaxPix = (int) ((m_fScaleMax - fMin) / fScaleL);
	}
      else
	{
	  nScaleMaxPix = (int) ((m_fScaleMax - fAvg5SD) / fScaleR + nSlopeChange);
	}
      nScaleMaxPix = max(0, nScaleMaxPix);
      nScaleMaxPix = min(m_nScaleWidth, nScaleMaxPix);
    }

  if (bMinChanged || bMaxChanged)
    {
      if (   m_bPixelSwapBytes && (2 == m_nBytesPerPixel))
	{
	  // Swap bytes on the 2-byte pixels
	  vSwapBytes(m_nScaleSize, m_pucScale);
	}

      XPutImage(XtDisplay(m_hParent),
		XtWindow(m_hParent),
		m_hGC,
		m_phXImageScale,
		nStart, 0, nStart, m_nHeightUsed, nEnd, m_nScaleHeight);
      // Make sure image is flushed before drawing
      // lines below, so they are not overwritten

//      XFlush(XtDisplay(m_hParent));
//      XSync(XtDisplay(m_hParent), False);
      XmUpdateDisplay(m_hParent);


      // Draw box around colorscale

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		0, m_nHeightUsed, m_nScaleWidth, m_nHeightUsed);

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		m_nScaleWidth, m_nHeightUsed,
		m_nScaleWidth, m_nHeightUsed+m_nScaleHeight);

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		m_nScaleWidth, m_nHeightUsed+m_nScaleHeight,
		0, m_nHeightUsed+m_nScaleHeight);

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		0, m_nHeightUsed+m_nScaleHeight,
		0, m_nHeightUsed);

      // Draw lines and triangles at nScaleMinPix & nScaleMaxPix
      // Draw filled triangles or wedges, rather than outlined

      char  cText[20];
      int   nLen, nTextWidth, nTextHeight, nForeground;
      XPoint a3hPoints[3];

      // Save foreground color for later

      nForeground = m_hGCValues.foreground;
      nTextHeight = m_pXFont->ascent + m_pXFont->descent;
      a3hPoints[0].x = 0;
      a3hPoints[0].y = m_nHeightUsed;
      a3hPoints[1].x = max(5, nScaleMinPix);  // Do not let triangle disappear
      a3hPoints[1].y = m_nHeightUsed  +  m_nScaleHeight / 2;
      a3hPoints[2].x = 0;
      a3hPoints[2].y = m_nHeightUsed  +  m_nScaleHeight;
      if (bMinChanged)
	{
	  XFillPolygon(XtDisplay(m_hParent),  XtWindow(m_hParent), m_hGCtext,
		       a3hPoints, 3, Convex, CoordModeOrigin);

	}
      nLen = sprintf(cText, "%.0f", m_fScaleMin);
      nTextWidth  = XTextWidth(m_pXFont, cText, nLen);      
      m_hGCValues.foreground = m_poColormap->ulGetColorIndex("black");
      XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCForeground,
		&m_hGCValues);
      XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		  max(5, nScaleMinPix - nTextWidth - 5),
		  a3hPoints[1].y + nTextHeight/2,
		  cText, nLen);
      m_hGCValues.foreground = nForeground;
      XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCForeground, &m_hGCValues);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		nScaleMinPix, m_nHeightUsed+m_nScaleHeight,
		nScaleMinPix, m_nHeightUsed);

      a3hPoints[0].x = m_nScaleWidth;
      a3hPoints[0].y = m_nHeightUsed;
      a3hPoints[1].x = min(m_nScaleWidth-5, nScaleMaxPix); // Do not let tri dis
      a3hPoints[1].y = m_nHeightUsed  +  m_nScaleHeight / 2;
      a3hPoints[2].x = m_nScaleWidth;
      a3hPoints[2].y = m_nHeightUsed  +  m_nScaleHeight;
      if (bMaxChanged)
	{
	  XFillPolygon(XtDisplay(m_hParent),  XtWindow(m_hParent), m_hGCtext,
		       a3hPoints, 3, Convex, CoordModeOrigin);
	}
      nLen        = sprintf(cText, "%.0f", m_fScaleMax);
      nTextWidth  = XTextWidth(m_pXFont, cText, nLen);      
      m_hGCValues.foreground = m_ulColorWhite;
      XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCFunction | GCForeground,
		&m_hGCValues);
      XDrawString(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		  min(m_nScaleWidth-5-nTextWidth, nScaleMaxPix + 5),
		  a3hPoints[1].y + nTextHeight/2,
		  cText, nLen);
      m_hGCValues.foreground = nForeground;
      XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCForeground, &m_hGCValues);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		nScaleMaxPix, m_nHeightUsed+m_nScaleHeight,
		nScaleMaxPix, m_nHeightUsed);

      // The average (which is always less than fAvg5SD!)

      nTemp = (int) ((fAvg - fMin) / fScaleL);
      nTemp = max(0, nTemp);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		nTemp, m_nHeightUsed,
		nTemp, m_nHeightUsed+m_nScaleHeight);

      // The suggested fScaleMax

      nTemp = (int) ((fAvg5SD - fMin) / fScaleL);
      nTemp = min(m_nScaleWidth, nTemp);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		nTemp, m_nHeightUsed,
		nTemp, m_nHeightUsed+m_nScaleHeight/2);

      // The suggested fScaleMin which is always less than fAvg
      
      nTemp = (int) ((fAvg - m_tImageProps.fScaleMinSD * m_tImageProps.fSD 
		      - fMin) / fScaleL);
      nTemp = max(0, nTemp);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGCtext,
		nTemp, m_nHeightUsed+m_nScaleHeight,
		nTemp, m_nHeightUsed+m_nScaleHeight/2);
    }
  m_nScaleMode = 1;
}

void
CXdisplay::vSetResol(Cdetector *poDetector, const float *pfS0,
		      const float fWavelength, const int nRings)
{
  // This routine needs to be called at least once before calling vDrawResol.

  m_poDetector = poDetector;
  for (int i = 0; i < 3; i++)
    {
      m_a3fS0[i] = pfS0[i];
    }
  m_fWavelength = fWavelength;
  m_tImageProps.nResoCircles = nRings;
}
void
CXdisplay::vSetBeamstopShadowMode(const int nMode)
{
  m_nBeamstopShadowMode = nMode;
  if (0 == nMode)
    {
      m_poCursor->vReset();
    }
  else
    {
      // TO DO: Change contrast here
      m_poCursor->vSetFont(XC_gumby, FALSE);
      //m_poCursor->vReset();
    }
  return;
}

int 
CXdisplay::nCalcDrawBeamstopShadow(float fPixIn0, float fPixIn1,
				   float fCutoffValue)
{
  // Calculate and Draw (i.e. mark) the pixels 
  // that are in the beam stop shadow

  int nStat1;

  float fBeamPx0, fBeamPx1;

  //m_nBeamstopShadowMode = 0;
  fBeamPx0 = fPixIn0;
  fBeamPx1 = fPixIn1;

  int nx, ny, nStat;
  int nWidth = 6;

  // Algorithm
  // Find low pixel values that are connected through other low pixel values
  // to the input pixel position.

  // 1. Start at beam center pixel.
  // 2. Average NxN pixel area and get average and SD.
  // 3. Determine a cutoff value; this cutoff value will change as 
  //    resolution changes.
  // 4. Set pixel values less than or equal to this cutoff to blue.

  if (-1.0 < fBeamPx0)
    {
      // Beam center on the image
      int nEdgeSize = 5;
      C3Ddata *po3Ddata;
      po3Ddata = new C3Ddata(nEdgeSize, nEdgeSize, 1, nint(fBeamPx0)-3, nint(fBeamPx1)-3, 0);
      nStat1 = po3Ddata->nFill2D(m_poImage, -1);

      float a3fBkg[3], fAvg, fSD, fMaxValue, a3fPeakCent[3], a3fTotCent[3];
      nStat1 = po3Ddata->nBackgroundAvgSD(0, &a3fBkg[0], &fAvg, &fSD, 
					  &fMaxValue,
					  &a3fPeakCent[0], &a3fTotCent[0]);

      delete po3Ddata;
      //cout << "Average is: " << fAvg << endl;

      float fThreshold = max(fAvg + 3.0 * fSD, 3.0 * fAvg);
      //cout << "Enter new threshold: ";
      //cin >> fThreshold;
      //fThreshold = m_tImageProps.fScaleMax;
      if (0.0f >= fCutoffValue)
	{
	  cout << "Threshold is " << fThreshold << endl << flush;
	}

      cout << "INFO: Automatically trying to flag beamstop shadow pixels...\n"
           << flush;
      // Now 
      int nCent0, nCent1;
      nCent0 = nint(fBeamPx0);
      nCent1 = nint(fBeamPx1);
      
      m_nRecurseCount   = 0;
      m_nShadowPixCount = 0;
      float fSlope = 50.0;
      Cimage oImage(*m_poImage, m_poImage->pcTheData());
      nEdgeSize = 5;
      nStat1 = nRecurseShadow(nCent0, nCent1, nEdgeSize, fAvg, fSD, fThreshold,
			      fSlope, &oImage);

      cout << "INFO: Number of pixels added to the beamstop shadow: "
           <<  m_nShadowPixCount << endl << flush;
    }
  vRefresh();
  return (0);
}

int  
CXdisplay::nCDBeamstopShadow(void)
{
  // This should probably go in the Cimage class.

  // Algorithm (note: this is not actually implemented below!)
  // 1. Predict pixel position of beam center.
  // 2. Computer a radius of 1/6 the edge of the image.
  // 3. At that radius, for angle = 0 to 360, step 10, check an area of 
  //    5 x 5 pixels and get an average background value.
  // 4. Find the angle of the lowest average pixel value.
  // 5. Use a simplex refinement to go to lowest average angle.
  // 6. Repeat steps 3-5 for a radius of 1/4 the edge of the image.
  // 7. Connect the two points for a line.
  // 8. Move out from beam center along the line and do a perpendicular
  //    line at each pixel, determine
  // 9. Set pixel values along the line below a threshold to 0.

  int i, j;
  int nStat1, nStat2;

  if (NULL == m_poDetector)  // No detector info, so do not draw anything
    return (1);


  float fPx0, fPx1, fMm0, fMm1;
  int   nDpyPx0, nDpyPx1, nDpyPx2, nDpyPx3;

  // Get current beam position (if on detector)

  float a3fX[3], a3fXR[3];
  float a3x3fMat[3][3];
  a3fX[0] = 0.0;
  a3fX[1] = 0.0;
  a3fX[2] = 0.0;
  float fBeamPx0, fBeamPx1;

  nStat1 = m_poDetector->nCalcDetCoords(a3fX, m_a3fS0,
				       &fMm0, &fMm1, &fBeamPx0, &fBeamPx1);
  if (0 != nStat1)
    return (1);

  int nDim0, nDim1;
  int a5nRadius[5];
  nDim0 = m_poImage->nGetDimension(0);
  nDim1 = m_poImage->nGetDimension(1);
  for (i = 0; i < 5; i++)
    {
      a5nRadius[i] = (nDim0 + nDim1) / 2;
    }

  return (0);
}
float
CXdisplay::fLineShadow(float fPx0, float fPx1, float fSlope, float fIntercept, 
		      float fWidth, int nNumPoints, int nSmooth, int nDir)
{
  // TODO: If Abs(fSlope) > 1.0, then loop over pixels a different way
  //       to find the diffuse edge

  // JWP:  Maybe this routine can be used to get the threshold value
  //       in a local area for the other routine?
  //
  // Given a line: Y = fSlope * X + fIntercept, 
  // a point on/near the line: fPx0, fPx1
  // and a direction to go in: nDir
  // Do the following:
  // 1. Copy pixel values in an image to a temporary variable
  // 2. Smooth the values of the points
  // 3. Compute 1st and 2nd derivatives of the values
  // 4. Determine a threshold for the values, where point at fPx0, fPx1
  //    is considered low, and something at the end is considered high.
  // 5. Set pixel values in the original image to 0 from fPx0, fPx1 to
  //    to cutoff point.

  // Cimage is m_poImage

  int nx, ny;
  int np, nxx;
  int nDim0, nDim1;
  float *pfValues;
  float *pfSmoothValues;
  float *pfSmSortedValues;
  float *pfFirstDeriv;
  float *pfSecondDeriv;

  pfValues       = new float [nNumPoints+1];
  pfSmoothValues = new float [nNumPoints+1];
  pfSmSortedValues = new float [nNumPoints+1];
  pfFirstDeriv   = new float [nNumPoints+1];
  pfSecondDeriv  = new float [nNumPoints+1];

  // 1. Copy pixel values in an image to pfValues

  nDim0 = m_poImage->nGetDimension(0);
  nDim1 = m_poImage->nGetDimension(1);

  if (0 == nDir) nDir  = 1;

  // 1. Copy nNumPoints image values into pfValues;

  // Find limits within the image for 

  nx    = nint(fPx0);
  ny    = nint(fPx1);
  int nw;
  float fSum;
  int nCount;
  int nWidth = nint(fWidth);
  if (0.0 == fSlope)
    {
      for (np = 0; np < nNumPoints; np++)
	{
	  fSum = 0.0;
	  nCount = 0;
	  pfValues[np] = 0.0;
	  // ny does not change
	  ny = nint(fPx1);
	  if ( (0 <= nx) && (nDim0 > nx) )
	    {
	      for (nw = 0; nw < nWidth; nw++)
		{
		  if ( (0 <= ny) && (nDim1 > ny) )
		    {
		      fSum += (m_poImage->*m_poImage->prfGetPixel)(nx, ny);
		      nCount++;
		    }
		  ny++; // Iterate over width which changes the Y intercept
		}
	    }
	  if (0 < nCount) pfValues[np] = fSum / (float) nCount;
	  nx += nDir;
	}
    }
  else // Look for an edge along a vertical line, so look along ny and not nx
    {
      for (np = 0; np < nNumPoints; np++)
	{
	  fSum = 0.0;
	  nCount = 0;
	  pfValues[np] = 0.0;
	  nx = nint(fPx0);
	  if ( (0 <= ny) && (nDim1 > ny) )
	    {
	      for (nw = 0; nw < nWidth; nw++)
		{
		  if ( (0 <= nx) && (nDim0 > nx) )
		    {
		      fSum += (m_poImage->*m_poImage->prfGetPixel)(nx, ny);
		      nCount++;
		    }
		nx++; // Iterate over width which changes the X intercept
		}
	    }
	  if (0 < nCount) pfValues[np] = fSum / (float) nCount;
	  ny += nDir;
	}
    }

  // 2. Smooth the pfValues array into pfSmoothValues;
  //    Use 1/20th of the total points, but 5 points minimum for smoothing
  
  int nSmPoints = 0;
  nCount = 0;
  fSum = 0.0f;
  
  nSmPoints = max(5, nNumPoints / 20);
  for (np = 0; np < nNumPoints; np++)
    {
      fSum   = 0.0f;
      nCount = 0;
      pfSmoothValues[np] = 0.0;
      for (nx = max(0, (np - nSmPoints /2)); 
	   nx < min(nNumPoints, np + (nSmPoints/2)); nx++)
	{
	  if (0.0f < pfValues[nx])
	    {
	      nCount++;
	      fSum += pfValues[nx];
	    }
	}
      if (nCount > 0) pfSmoothValues[np] = fSum / (float) nCount;
      pfSmSortedValues[np] = pfSmoothValues[np];
    }

  // Find the median smoothed pixel value 

  qsort(&pfSmSortedValues[0], nNumPoints, sizeof(float), float_cmp);
  float fMedianValue;
  fMedianValue = pfSmSortedValues[nNumPoints/2];
  
  // Remove the top 15% of pixel values

  float fTopValueLimit;

  fTopValueLimit = pfSmSortedValues[nNumPoints*85/100];
  
  // Set all pixel values above the fTopValueLimit to fTopValueLimit
  for (np = 0; np < nNumPoints; np++)  
    {
      if (fTopValueLimit < pfSmoothValues[np])
	pfSmoothValues[np] = fTopValueLimit;
    }

  // Find the location of the first pixel above the median value

  int nMedi = 0;
  while (pfSmoothValues[nMedi] < fMedianValue)
    {
      nMedi++;
    }
  //cout << "nMedi = " << nMedi << endl;

  // 3. Compute 1st and 2nd derivatives

  pfFirstDeriv[0]  = 0.0f;
  pfSecondDeriv[0] = 0.0f;
  for (np = 1; np < nNumPoints; np++)
    {
      if (0.0 < pfValues[np])
	{
	  pfFirstDeriv[np]  = pfSmoothValues[np] - pfSmoothValues[np - 1];
	  pfSecondDeriv[np] = pfFirstDeriv[np]   - pfFirstDeriv[np - 1];
	}
      else
	{
	  pfFirstDeriv[np]  = -9999.0f;
	  pfSecondDeriv[np] = 0.0f;
	}
    }

  //  4. Determine a threshold for the values.
  //    a. Find point maximum 1st derivative starting at the expected high value edge
  //       but only look up to the median pixel value
  //    b. Find the average of the 1st derivative at the nSmPoints high end (should be low)
  //    c. Find the point between (a) and then that is at least 10% heigher than value (b).

  int nMax;
  int nMin;
  nMax  = 0;
  nMin  = 0;
  for (np = 1; np < nMedi; np++)
    {
      if (pfFirstDeriv[nMax] < pfFirstDeriv[np])
	{
	  nMax = np;
	}
      if (pfFirstDeriv[nMin] > pfFirstDeriv[np])
	{
	  nMin = np;
	}
    }

  // If the minimum first derivative is hugely negative, this suggests a peak
  // in the array that we do not want to include in the mask.

  //if (nMin > nMax) cout << " WARNING peak for fpixin1: " << fPx1 << endl;

  np = nMedi;
  float a2fAvg[2];
  for (nx = 0; nx < 2; nx++)
    {
      fSum   = 0.0f;
      nCount = 0;
      while ( (0 <= np) && (nCount < nSmPoints) )
	{
	  if (-9999.0f != pfFirstDeriv[np])
	    {
	      fSum += pfFirstDeriv[np];
	      nCount++;
	    }
	  np--;
	}
      a2fAvg[nx] = fSum / (float) nCount;
    }

  // Of these two averaged derivatives, select the slope closest to 0;

  if (ABS(a2fAvg[0]) > ABS(a2fAvg[1]))
    {
      a2fAvg[0] = a2fAvg[1];
    }
  fSum = a2fAvg[0];

  // fSum contains the average 1st derivative closest to "flat"

  int nFlat;
  np = nMax;
  float fThreshold = -999.0;
  while (   (pfFirstDeriv[np] != -9999.0)
         && (pfFirstDeriv[np] > fSum)
	 && (np < nNumPoints) )
    {
      fThreshold = pfSmoothValues[np];
      np++;
    }
  nFlat = np;
  nFlat = min(nMax + 2, nNumPoints);

  // 4.5
  // Report the max value to use as a threshold or cutoff
  // elsewhere:

  if (pfSmoothValues[nFlat] > fThreshold)
    fThreshold = pfSmoothValues[nFlat];

  delete [] pfFirstDeriv;
  delete [] pfSecondDeriv;
  delete [] pfValues;
  delete [] pfSmoothValues;
  delete [] pfSmSortedValues;
  return (fThreshold);
}

int
CXdisplay::nRecurseShadow(int nCent0, int nCent1, 
			  int nEdgeSize, float fAvg, float fSD,
			  float fThreshold, float fSlopePixVal, 
			  Cimage *poImageOut)
{
  // fSlopePixVal is not used just yet.
  
  // Look for pixels in m_poImage which have a value below fThreshold
  // in a region centered on pixel coordinate (nCent0, nCent1) 
  // of area nEdgeSize by nEdgeSize.  Any pixel values below fThreshold
  // will get set to a value of 0 in *poImageOut.
  // If there are no new pixels below the threshold, then exit, otherwise
  // if there are new pixels below the threshold, call this routine again
  // with an adjacent value of (nCent0, nCent1).  For the adjacent area 
  // go out in 4 different directions.

  int nLDebug = 1;
  int nStat;
  int  nCountNewPixelsInShadow;
  int  nCountNewPixelsInShadowDir;
  int nDir;
  int a4x2nEdge[4][2];
  int a2nCent[2];
  int nHalfEdge;
  a2nCent[0] = nCent0;
  a2nCent[1] = nCent1;
  int a5x2nDir[5][2] = { {0, 0}, {1, 0}, {-2, 0}, {0, 1}, {0, -2} };
      
  nHalfEdge = nEdgeSize / 2;

  m_nRecurseCount++;
  //cout << "nRecurseShadow " << nCent0 << ", " << nCent1 
  // << "   count: " << m_nRecurseCount << endl << flush;
  int nTooMany;
  nTooMany = m_poImage->nGetDimension() / 50;  // About 2% of pixels

  // Do not let recursion go forever
  if (   (nTooMany < m_nRecurseCount) 
	 || (m_nShadowPixCount > (nTooMany * 5) ))
    {
      cout << "WARNING: Too many calls to CXdisplay::nRecurseShadow(...)\n"
           << "         Perhap the threshold is too high.\n"
           << endl;
      return (1);
    }
  // Try to calculate new average threshold here?
  C3Ddata *po3Ddata;
  po3Ddata = new C3Ddata(nEdgeSize, nEdgeSize, 1, 
			 a2nCent[0] - nHalfEdge,
			 a2nCent[1] - nHalfEdge,
			 0);
  nStat = po3Ddata->nFill2D(m_poImage, -1);
  float a3fBkg[3], fAvg2, fAvg3, fSD2, fMaxValue, a3fPeakCent[3], a3fTotCent[3];
  nStat = po3Ddata->nBackgroundAvgSD(0, &a3fBkg[0], &fAvg2, &fSD2, 
				     &fMaxValue,
				     &a3fPeakCent[0], &a3fTotCent[0]);
  delete po3Ddata;
  fAvg3 = fAvg2;
  po3Ddata = new C3Ddata(nEdgeSize, nEdgeSize, 1, 
			 a2nCent[0] - 20 * nHalfEdge,
			 a2nCent[1] - 20 * nHalfEdge,
			 0);
  nStat = po3Ddata->nFill2D(m_poImage, -1);
  if (0 == nStat)
    nStat = po3Ddata->nBackgroundAvgSD(0, &a3fBkg[0], &fAvg3, &fSD2, 
				       &fMaxValue,
				       &a3fPeakCent[0], &a3fTotCent[0]);
  delete po3Ddata;

  // Proceed only if the AVERAGE is less than the threshold

  if (fAvg2 > fThreshold)
    {
      //      cout << "Avg > Thresh: " << fAvg2 << " > " << fThreshold << endl;
      return (0);
    }

  // Try to find a suitable threshold value in this region


  float fThresholdUsed;  // Maybe threshold should be lowered as one goes to 
                         // to higher resolution?

 //                    pix0      pix1   slope  intercept
  float fPx0, fPx1;
  fPx0 = nint(nCent0);
  fPx1 = nint(nCent1);
  int nWidth = 12;
  float a4fThresh[4] = {0.0, 0.0, 0.0, 0.0};
  // y = slope * x   + y-intercept
  // 2 horizontal lines
  //                           x    y   slope intercept
  a4fThresh[0] = fLineShadow(fPx0, fPx1, 0.0, fPx1, 
		      (float)nWidth, 151, 7, 1);
  a4fThresh[1] = fLineShadow(fPx0, fPx1, 0.0, fPx1, 
		      (float)nWidth, 151, 7, -1);

  // 2 more vertical lines

  a4fThresh[2] = fLineShadow(fPx0, fPx1, 1000.0, fPx0,
		      (float)nWidth, 151, 7, 1);
  a4fThresh[3] = fLineShadow(fPx0, fPx1, 1000.0, fPx0, 
		      (float)nWidth, 151, 7, -1);
  int np;
  int nCount = 0;
  float fSum = 0.0;

  // Average the top two thresholds

  qsort(&a4fThresh[0], 4, sizeof(float), float_cmp);

  for (np = 2; np < 4; np++)
    {
      if (0.0 < a4fThresh[np])
	{
	  fSum += a4fThresh[np];
	  nCount++;
	}
    }
  if (0 < nCount)
    fSum = fSum / (float)nCount;

  fThresholdUsed = 0.6 * fSum;

  if (0.0 >= fThresholdUsed)
    {
      // fLineShadow did not do a good job;
      fThresholdUsed = fThreshold;
    }
  //cout << "nCent0,1: " << nCent0 << ", " <<  nCent1 << endl;
  //cout << "tholdused, hold: " << fThresholdUsed << "   " << fThreshold << endl;
  nCountNewPixelsInShadow = 0;

  for (nDir = 0; nDir < 5; nDir++)
    {
      nCountNewPixelsInShadowDir = 0;

      // Go look in 5 directions (really 4, but include the center)

      // Shift to a new center to look at

      a2nCent[0] += (18 * nEdgeSize * a5x2nDir[nDir][0] / 20);
      a2nCent[1] += (18 * nEdgeSize * a5x2nDir[nDir][1] / 20);

      int nx, ny;
      int nxx, nyy;
      bool bSetToZero;
      float fPix1;
      float fPix2;
      float fPixPrev = 0.0f;
      nyy = a2nCent[1] - nHalfEdge;
      for (ny = 0; ny < nEdgeSize; ny++)
	{
	  nyy++;
	  nxx = a2nCent[0] - nHalfEdge;
	  for (nx = 0; nx < nEdgeSize; nx++)
	    {
	      nxx++;
	      fPix1 = (m_poImage->*m_poImage->prfGetPixel)(nxx, nyy);
	      fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx, nyy);	      
	      bSetToZero = FALSE;
	      if ( (fPix1 < fThresholdUsed) && (0.0 != fPix1) )
		{
		  // if pixel in m_poImage is less < threshold,
		  // then check pixel in poImageOut is not "set"
		  // if not set, then set it and set nNewPixels to true
		    
		  if (fPix2 == fPix1)
		    {
		      // Not set, so set to 0.0

		      bSetToZero = TRUE;
		    }
		}
	      else if (0.0 != fPix2)
		{
		  // If 5 of 8 nearest neighbors are set to 0, then set this
		  // this one to 0 as well

		  fPix1 = fPix2;
		  int nCount = 0;
		  float fAvg = 0.0;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx+1, nyy-1);
		  if (0.0 == fPix2) 
		    nCount++;
		  else
		    fAvg += fPix2;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx+1, nyy);
		  if (0.0 == fPix2)
		    nCount++;
		  else
		    fAvg += fPix2;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx+1, nyy+1);
		  if (0.0 == fPix2)
		    nCount++;
		  else
		    fAvg += fPix2;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx, nyy-1);
		  if (0.0 == fPix2)
		    nCount++;
		  else
		    fAvg += fPix2;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx, nyy+1);
		  if (0.0 == fPix2)
		    nCount++;
		  else
		    fAvg += fPix2;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx-1, nyy-1);
		  if (0.0 == fPix2)
		    nCount++;
		  else
		    fAvg += fPix2;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx-1, nyy);
		  if (0.0 == fPix2)
		    nCount++;
		  else
		    fAvg += fPix2;
		  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx-1, nyy+1);
		  if (0.0 == fPix2)
		    nCount++;
		  else
		    fAvg += fPix2;

		  bSetToZero = (5 <= nCount);
		  if (bSetToZero) 
		    {
		      //cout << "5/8 set\n" << flush;
		    }
		  else
		    {
		      // Compute average of nearest non-zero neighbors
		      nCount = 8 - nCount;
		      if ( (0 < nCount) && (5 >= nCount) )
			{
			  // At least 3 neighboring pixels should be in the
			  // the shadow

			  fAvg = fAvg / float(nCount);
			  if (fPix1 < (0.85 * fAvg))
			    {
			      // If pixel value is less than 60% of average
			      // then it is in the shadow

			      bSetToZero = TRUE;
			      //cout << "0.85 avg set\n" << flush;
			    }
			}
		    }
		}
	      if (bSetToZero)
		{
		  fPix2 = 0.0;
		  (void) (poImageOut->*poImageOut->prnSetPixel)(nxx, nyy,
								fPix2);
		  if (1 == nLDebug)
		    (void) (m_poImage->*m_poImage->prnSetPixel)(nxx, nyy,
								fPix2);
		  nCountNewPixelsInShadowDir++;
		}
	    }
	}
      int nTest = 0;
      // Only call nRecurseShadow if a significant number of pixels in the
      // just tested region were also in the shadow.
      nTest = nEdgeSize * nEdgeSize * 10 / 25;
      if (nTest < nCountNewPixelsInShadowDir)
	{
	  m_nShadowPixCount = m_nShadowPixCount + nCountNewPixelsInShadowDir;
	  //cout << "ShadowPixCount: " << m_nShadowPixCount << endl << flush;
	  
	  // Do not call nRecurseShadow if you are near the edge of the 
	  // detector
	  int nTooClose;
	  nTooClose = max(5, nEdgeSize);
	  if (    (a2nCent[0] < nTooClose)
               || (a2nCent[1] < nTooClose)
               || ( m_poImage->nGetDimension(0) - a2nCent[0] < nTooClose)
               || ( m_poImage->nGetDimension(1) - a2nCent[1] < nTooClose) )
	    {
	      //cout << "Too close to edge\n" << flush;
	    }
	  else
	    {
	      
	      nStat = nRecurseShadow(a2nCent[0], a2nCent[1],
				     nEdgeSize, fAvg, fSD, 
				     fThreshold, fSlopePixVal,
				     poImageOut);
	    }
	  if (1 == nLDebug)
	    if (0 == (m_nRecurseCount % 500))
	      vRefresh();
	  if (0 != nStat) return (nStat);
	}
      else
	{
	  //cout << "No new pixels in shadow for "
	  //<< nCent0 << ", " << nCent1 << endl << flush;
	}
      nCountNewPixelsInShadow += nCountNewPixelsInShadowDir;
    } // end direction test

  return (0);
}

void
CXdisplay::vDrawResol(void)
{
  int i, j;
  int nStat1, nStat2;
  float a8fIceRings[8] = { 3.917, 3.684, 3.458, 2.683, 2.261, 2.081, 1.9584, 1.9272 }; 
  // Draw resolution arcs on the displayed image. 
  // For this to work, the source vector magnitude should be 1.

  if (   (NULL == m_poDetector)  // No detector info, so do not draw anything
      || (-1 == m_tImageProps.nResoCircles) ) // No circles, so do not draw 
    return;
//+2011-10-12 JWP
  if ("" != sGetEnv("DTDISPLAY_ICE_RINGS"))
    {
      Cstring sTemp;
      Cstring a8sTemp[8];
      sTemp = sGetEnv("DTDISPLAY_ICE_RINGS");
      j = split(sTemp, a8sTemp, 8, " ");
      for (i = 0; i < 8; i++)
	{
	  if (i < j) 
	    a8fIceRings[i] = atof(a8sTemp[i].string());
	  else
	    a8fIceRings[i] = 9999.0;
	  cout << "New ice reso: " << a8fIceRings[i] << endl;
	}
    }
//-2011-10-12 JWP

//  nStat = m_poDetector->poSpatial->nPixeltoMM(f1, f2, &f1MM, &f2MM);
//  if (0 == nStat)
  //printf("Draw Reso \n");
  float fResMin, fResMax;
  nStat1 = m_poDetector->nGetResolution(m_a3fS0, &fResMin, &fResMax);
  if (0 != nStat1)
    return;

//  cout << "Resmin: " << fResMin << ", Resmax: " << fResMax << '\n';

  // Change line color to red

  int nForeground = m_hGCValues.foreground;
  m_hGCValues.foreground = m_poColormap->ulGetColorIndex("red");
  XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCForeground, 
	    &m_hGCValues);

  float fPx0, fPx1, fMm0, fMm1;
  int   nDpyPx0, nDpyPx1, nDpyPx2, nDpyPx3;

/*
  // Quadrilateral drawing
  
  cout << "Q2: " << m_a5fQuad[0][0] << ", " << m_a5fQuad[0][1] << endl;
  if ( (0.0 <= m_a5fQuad[0][0]) && (0.0 <= m_a5fQuad[0][1]) )
    {
      m_a5fQuad[4][0] = m_a5fQuad[0][0];
      m_a5fQuad[4][1] = m_a5fQuad[0][1];

      nStat1 = nImgPixToDpyPix(m_a5fQuad[0][0], m_a5fQuad[0][1], 
			       &nDpyPx0, &nDpyPx1);
      if (0 == nStat1)
	{
	  for (i = 1; i < 5; i++)
	    {
	      nStat2 = nImgPixToDpyPix(m_a5fQuad[i][0], m_a5fQuad[i][1], 
				       &nDpyPx2, &nDpyPx3);
	      if ( (0 == nStat1) && (0 == nStat2) )
		{
		  // Both endpoints are on the displayed image, so draw the line

		  XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCtext,
			    nDpyPx0, nDpyPx1, nDpyPx2, nDpyPx3);
		}
	      nDpyPx0 = nDpyPx2;
	      nDpyPx1 = nDpyPx3;
	      nStat1  = nStat2;
	    }
	}
    }
*/
  // Mark current beam position (if on detector)

  float a3fX[3], a3fXR[3];
  float a3x3fMat[3][3];
  a3fX[0] = 0.0;
  a3fX[1] = 0.0;
  a3fX[2] = 0.0;

  nStat1 = m_poDetector->nCalcDetCoords(a3fX, m_a3fS0,
				       &fMm0, &fMm1, &fPx0, &fPx1);
  if (0 == nStat1)
    {
      nStat1 = nImgPixToDpyPix(fPx0, fPx1, &nDpyPx0, &nDpyPx1);
      if (0 == nStat1)
        {
	  // Draw + symbol

          XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCtext,
                    nDpyPx0-m_tImageProps.nIntegrateBox[0]/2, nDpyPx1,
		    nDpyPx0+m_tImageProps.nIntegrateBox[0]/2, nDpyPx1);
          XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCtext,
                    nDpyPx0, nDpyPx1-m_tImageProps.nIntegrateBox[1]/2, 
		    nDpyPx0, nDpyPx1+m_tImageProps.nIntegrateBox[1]/2);
        }
    }

  int nNumRings = m_tImageProps.nResoCircles;
  bool bDoIce = FALSE;
  if (8 == nNumRings) // If set in calling program, this is just a signal to
    {                 // use the a8fIceRings for the resolutions of the rings
      nNumRings = 8;
      bDoIce = TRUE;
      //cout << "ice rings: ";
      //for (i = 0; i < 8; i++)
      //    cout << a8fIceRings[i] << endl;
    }
  for (i = 1; i <= (nNumRings+1); i++)
    {
      // Compute initial reciprocal lattice vector a3fX at a given resolution
      // and get pixel coordinates thereof
      // 1. Get resolution that we are interested in.

      if (bDoIce && (9999.0 == a8fIceRings[i-1])) break; // This ring should NOT be drawn

      float fScale;
      
      fScale =   1. / (fResMax * fResMax)
	         * ((float) i / (float) (nNumRings+1));
      fScale = (float)sqrt(1.0/(double)fScale);

      // 2. Get a reciprocal lattice vector at this resolution that lies on
      //    the Ewald sphere

      // 2a. Find a vector perpendicular to S0

      a3fXR[0] = 1.0;
      a3fXR[1] = 1.0;
      a3fXR[2] = 1.0;
      vCross3D(a3fXR, m_a3fS0, a3fX);

      // 2b. Normalize its length

      (void) fNormVec3D(a3fX);

      // 3. Find 2theta of this resolution
      
      float fSinTheta = 0.5 / fScale;

      // Special Ice Ring code
      // nL = 2d sintheta; sintheta = nL / 2d
      
      if ( bDoIce && (a8fIceRings[i-1] >= fResMax) && (a8fIceRings[i-1] <= fResMin) )
	{
	  fSinTheta = m_fWavelength * 0.5 / a8fIceRings[i-1];
	}
      else if (bDoIce)
	continue;  // This means an ice ring was requested that would not occur on the image, so ...

      float fTheta = (float)asin((double)fSinTheta) / Gs_dRADIANS_PER_DEGREE;

      // 4. Rotate the negative of the source vector by the 2theta
      // angle around the vector just determined to give a
      // scattered beam wavevector S for this resolution
      
      vMulVec3DScalar(m_a3fS0, -1.0, a3fXR);
      vConvRotVec3DMat3D(2.0 * fTheta, a3fX, &a3x3fMat[0][0]);      
      vMulMat3DVec3D(a3x3fMat, &a3fXR[0], &a3fX[0]);
      
      // 5. Compute X = S0 + S
      //            -   --   -

      vAddVec3DVec3D(a3fX, m_a3fS0, a3fX);

      nStat1 = m_poDetector->nCalcDetCoords(a3fX, m_a3fS0,
					    &fMm0, &fMm1, &fPx0, &fPx1);
      if (0 == nStat1)
	nStat1 = nImgPixToDpyPix(fPx0, fPx1, &nDpyPx0, &nDpyPx1);

      int  nDoText = 0; // Draw resolution as text, but only at point 2
      for (j = 1; j <= 360; j = j+1)
	{
	  // Rotate a3fX by j degrees around source vector
	  // See if it impinges on detector

	  vConvRotVec3DMat3D((float) j, m_a3fS0, &a3x3fMat[0][0]);

	  vMulMat3DVec3D(a3x3fMat, &a3fX[0], &a3fXR[0]);
	  nStat2 = m_poDetector->nCalcDetCoords(a3fXR, m_a3fS0,
						&fMm0, &fMm1,
						&fPx0, &fPx1);
	  if (0 == nStat2)
	    nStat2 = nImgPixToDpyPix(fPx0, fPx1, &nDpyPx2, &nDpyPx3);
	  if ( (0 == nStat1) && (0 == nStat2) )
	    {
	      // Both endpoints are on the displayed image, so draw the line

	      nDoText++;
	      XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCtext,
			nDpyPx0, nDpyPx1, nDpyPx2, nDpyPx3);
	      if ( (2 == nDoText) && ("" == sGetEnv("DTDISPLAY_RINGS_NOLABEL") ) )
		{
		  if (bDoIce && (9999.0 == a8fIceRings[i-1]))
		    {
		      // No label
		    }
		  else
		    {
		  // Try to put text at end of 2nd segment
		  char cTemp[20];
		  float fValue = fLenVec3D(a3fX);
		  if (0.0 != fValue)
		    { 
		      fValue = m_fWavelength / fValue;
		      sprintf(cTemp, " %.2f", fValue);
		      int nTextLen = strlen(cTemp);
		      //printf("Reso calc: %.2f\n", fValue);
		      XDrawString(XtDisplay(m_hParent), m_hPixmap, m_hGCtext,
				  nDpyPx2, nDpyPx3, cTemp, nTextLen);
		    }
		    }
		}
	    }

	  // Set start of next line segment to end of last one

	  nDpyPx0 = nDpyPx2;
	  nDpyPx1 = nDpyPx3;
	  nStat1  = nStat2;
	}
    }

  // Return color back to what it was

  m_hGCValues.foreground = nForeground;
  XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCForeground,
	    &m_hGCValues);

  //vPlotDrawObjects();
}

void 
CXdisplay::vPlotDrawObjects(void)
{
  // This routine looks in the image header and draws lines based on
  // the numbers found in the DTDISPLAY_DRAW_LINES keyword.
  // The keyword has sequences of flag, Pix0, Pix1 that are [M|D] for move, 
  //    drawpixel numbers that are move/draw

  if (NULL == m_poImage)
    return;

  float fDegreeOffset = 0.0;
  int nStat;
  int nNumValues = 0;
  nStat = m_poImage->m_oHeader.nGetValue("ROTATION", &fDegreeOffset);
  //cout << "fDegreeOffset is: " << fDegreeOffset << endl << flush;
  nStat = m_poImage->m_oHeader.nGetValue("DTDISPLAY_DRAW_LINES",
					 &nNumValues);
  
  if (0 < nNumValues)
    {
      float *pfCoords;
      pfCoords = new float [nNumValues+1];
      nStat = m_poImage->m_oHeader.nGetValue("DTDISPLAY_DRAW_LINES",
					     nNumValues+1, pfCoords);
  
      if (0 != nStat)
	{
	  delete [] pfCoords;
	  // Error message here
	  return;
	}

      // Change line color to red

      int nForeground = m_hGCValues.foreground;
      m_hGCValues.foreground = m_poColormap->ulGetColorIndex("red");
      XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCForeground, 
		&m_hGCValues);

      float fPx0, fPx1, fMm0, fMm1;
      int   nDim0, nDim1;
      int   nDpyPx0, nDpyPx1, nDpyPx2, nDpyPx3;
      int   nStat0, nStat1, nStat2;

      m_poImage->nGetDimensions(&nDim0, &nDim1);

      // Loop through coordinates and draw the line

      int nx;
      int nNumLines = nNumValues / 4;
      float *pfTemp;
      pfTemp = &pfCoords[1];
      float fPx0Offset;

      for (nx = 0; nx < nNumLines; nx++)
	{
	  fPx0Offset = 20.0 - fDegreeOffset;
	  fPx0Offset *= 22.1667;

	  fPx0 = *pfTemp++;
	  fPx0 += fPx0Offset;

	  // RAPID pixels go from 0 to 4599 along Px0. 
	  // A full 360 degrees will be 7980 pixels.
	  // That is 22.16667 pixels per degree
	  // So an object must fit between 0 and 4599 pixels

	  if (7980.0 <= fPx0) fPx0 -= 7980.0;
	  if ((4600.0 - 7980.0) >= fPx0) fPx0 += 7980.0;
	  if (-1690.0 > fPx0) fPx0 = 4599.5;
	  else if (0.5 > fPx0) fPx0 = 0.5;
	  if (nDim0 <= fPx0) fPx0 = (float) nDim0 - 1.5;
	  fPx1 = *pfTemp++;
	  cout << "fPx0a,fPx1a " << fPx0 << ", " << fPx1 << endl;
	  nStat1 = nImgPixToDpyPix(fPx0, fPx1, &nDpyPx0, &nDpyPx1);
	  //if (0 == nStat1)
	    {
	      fPx0 = *pfTemp++;
	      fPx0 += fPx0Offset;
	      if (7980.0 <= fPx0) fPx0 -= 7980.0;
	      if ((4600.0 - 7980.0) >= fPx0) fPx0 += 7980.0;
	      if (-1690.0 > fPx0) fPx0 = 4599.5;
	      else if (0.5 > fPx0) fPx0 = 0.5;
	      if (nDim0 <= fPx0) fPx0 = nDim0 - 1.5;
	      cout << "fPx0b,fPx1b " << fPx0 << ", " << fPx1 << endl;
	      fPx1 = *pfTemp++;
	      nStat1 = nImgPixToDpyPix(fPx0, fPx1, &nDpyPx2, &nDpyPx3);
	    }
	    //if (0 == nStat1)
	    XDrawLine(XtDisplay(m_hParent), m_hPixmap, m_hGCtext,
		      nDpyPx0, nDpyPx1, nDpyPx2, nDpyPx3);
	}

      // Return color back to what it was

      m_hGCValues.foreground = nForeground;
      XChangeGC(XtDisplay(m_hParent), m_hGCtext, GCForeground,
		&m_hGCValues);

      delete [] pfCoords;
    }
  return;
}

void 
CXdisplay::vSetSpotInfo(const float fIntensity, const float fSigmaI, 
			const float fResolution, 
			const float fPx0, const float fPx1,
			const float fRotMid, const float fRotWidth)
{
  // Set spot info for use later in a call to 
  // m_poReflnlist->vAddSpot(...)

  m_tSpotInfo.nStat        = 0;  // Spot info is usable
  m_tSpotInfo.fIntensity   = fIntensity;
  m_tSpotInfo.fSigmaI      = fSigmaI;
  m_tSpotInfo.fResolution  = fResolution;
  m_tSpotInfo.fCentPx0     = fPx0;
  m_tSpotInfo.fCentPx1     = fPx1;
  m_tSpotInfo.fRotMid      = fRotMid;
  m_tSpotInfo.fRotWidth    = fRotWidth;
}

void
CXdisplay::vSetLineWidth(const float fLineWidth )
{
  // Set line width for thick rubberband in IMAGE pixels.  
  // Convert to display pixels first

  if (0.0 < fLineWidth)
    m_fLineWidth = fLineWidth;
  //  cout << "Setting rubberband line width to " << m_fLineWidth << endl;
  if (NULL != m_poRubberband)
    {
      int nStat;
      float fPX0, fPY0;
      int   nPX0, nPY0, nPX1, nPY1;
      nPX0 = m_nWidth / 2;
      nPY0 = m_nHeight  / 2;

      nStat = nDpyPixToImgPix(nPX0, nPY0, &fPX0, &fPY0);
      if (0 == nStat) 
	{
	  // Image coordinate near middle of display

	  nStat = nImgPixToDpyPix(fPX0, fPY0 + fLineWidth, &nPX1, &nPY1);
	  if (0 == nStat) 				  
	    {
	      fPX0 = (float) (nPX0 - nPX1);
	      fPY0 = (float) (nPY0 - nPY1);
	      fPX0 = fPX0 * fPX0  +  fPY0 * fPY0;
	      if (0 < fPX0)
		fPX0 = (float)sqrt((double)fPX0);
	      else
		fPX0 = 0.0;

	      m_nLineWidth = nint(fPX0);
	      //	      cout << "Setting linewidth to " << fLineWidth << endl;
	      m_poRubberband->vSetLineWidth(m_nLineWidth);
	    }
	}
    }
}

int
CXdisplay::nCalcEllipse(float fObsPx0, float fObsPx1, 
			float fSemiMajor, float fSemiMinor, 
			float fTilt, int nPoints, XPoint* p36hPoints)
{
  // Calculate the points (in display pixels) for drawing an ellipse
  // centered at fObsPx0, fObsPx1 (in image pixels) with semimajor and semiminor
  // axes of length fSemiMajor and fSemiMinor (in image pixels), 
  // tilted by fTilt wrt the fObsPx0 axis.
  //
  //  a = fSemiMajor, b = fSemiMinor
  //
  //  x**2/a**2 + y**2/b**2 = 1 = cos**2 + sin**2
  
  // fTilt will generate a 2x2 rotation matrix: |  cosT -sinT |
  //                                            |  sinT  cosT |
  // Note this appears to be the negative of the fTilt since in X windows
  // The origin is at the upper left, x across and y DOWN.

  // Remember XPoint contains (short x, short y)
  int nStat;
  int nx;
  
  int nDpyPx0, nDpyPx1;
  float fX, fY, fA, fB, fXp, fYp;
  float fDegree, fSinD, fCosD, fSinT, fCosT;
  fSinT   = sin((double)(fTilt * Gs_fRADIANS_PER_DEGREE));
  fCosT   = cos((double)(fTilt * Gs_fRADIANS_PER_DEGREE));
  fA      = fSemiMajor;
  fB      = fSemiMinor;

  float fDelta;
  fDelta = 360.0 / (float)(nPoints-1);
  fDegree = 0.0;
  for (nx = 0; nx < nPoints; nx++)
    {
      fSinD   = sin((double)(fDegree * Gs_fRADIANS_PER_DEGREE));
      fCosD   = cos((double)(fDegree * Gs_fRADIANS_PER_DEGREE));
      fX      =  fA * fCosD;
      fY      =  fB * fSinD;
      fXp     =  (fCosT * fX) + (-fSinT * fY) + fObsPx0;
      fYp     =  (fSinT * fX) + ( fCosT * fY) + fObsPx1;

      nStat = nImgPixToDpyPix( fXp, fYp, &nDpyPx0, &nDpyPx1);
      //      if (0 != nStat)
      // return (nStat);
      p36hPoints[nx].x = nDpyPx0;
      p36hPoints[nx].y = nDpyPx1;
      fDegree += fDelta;
    }
  return (nStat);
}

void 
CXdisplay::vSwapBytes(const int nNumChars, unsigned char *pucArray)
{
  // Swap pairs of bytes in the input char array

  int i;
  unsigned short int uiTemp;
  unsigned short int *puiTemp;
  puiTemp = (unsigned short int*) pucArray;
  i = nNumChars / 2;
  while (0 < i)
    {
      *puiTemp = ((*puiTemp & 0xff) << 8) | ((*puiTemp) >> 8);
      puiTemp++;
      i--;
    }
}
//
//+Private functions

