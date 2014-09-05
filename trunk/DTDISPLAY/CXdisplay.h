#ifndef DT_CXDISPLAY_H
#define DT_CXDISPLAY_H
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
// CXdisplay.h        Initial author: J.W. Pflugrath           18-Jul-1995
//    This file is the header file for class CXdisplay.
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

#include "UIComponent.h"

#include <math.h>
#include <Xm/Xm.h>
#include <Xm/Text.h>
#include <X11/keysym.h>
#include "Dtrek.h"
#include "Cstring.h"
//#include <iostream.h>
//#include <iomanip.h>
#include "Cimage.h"
#include "CXcolormap.h"
#include "CXrubberband.h"
#include "CXcursor.h"
#include "Creflnlist.h"
#include "Cdetector.h"
#include "dtrekvec.h"

//+Definitions and constants

#ifndef XK_KP_Up
#define XK_KP_Up XK_Up
#endif
#ifndef XK_KP_Right
#define XK_KP_Right XK_Right
#endif
#ifndef XK_KP_Left
#define XK_KP_Left XK_Left
#endif
#ifndef XK_KP_Down
#define XK_KP_Down XK_Down
#endif

//+Code begin

typedef struct _tagImagePropsDisplay
{
  float fZoom;
  float fScaleMin;
  float fScaleMax;
  float fScaleMinSD;
  float fScaleMaxSD;
  float fImageMin;
  float fImageMax;
  float fAvg;
  float fSD;
  float fAspectRatio;        // Pixel size aspect ratio (fast / slow)
  int   nDim[2];             // Image dimensions
  int   nOrig[2];            // Origin in image pixels of displayed image
  int   nExt[2];             // Extents in image pixels of displayed image
  int   nIntegrateBox[2];    // Integration box size in pixels
  int   nAutoRescale;        // Flag whether to auto re-scale new images or not
  int   nResoCircles;        // Number of resolution circles to display
  int   nOrient;             // Orient flag (0 -> 7)
  int   nAspectMode;         // Flag to preserve (=1) or not pixel size aspect
  int   nDisplayReflns;      // Flag to display (>0) reflections if available
  int   nDisplayPixValues;   // Flag to display (>0) pixel values on zoomed imgs
  int   nPlotRockingCurve;   // Plot a rocking curve in lower right corner
  float fSatPixValue;        // Saturated pixel value
  float fTooLowPixValue;     // Too low a pixel value. Default: -2
  Cstring sSatPixColor;      // Color to plot saturated pixels in, but must be
                             //   in the color table already.
  Cstring sTooLowPixColor;   // Color to plot too low pixels in, but must be
                             //   in the color table already. Default: yellow
} tagImagePropsDisplay;

typedef struct _tagReflnPropsDisplay 
{
  Cstring sReflnColor;
  int     nReflnPredObs;
  int     nReflnSymbol;
  float   a6fReflnSize[6];
  float   fImageRotStart;
  float   fImageRotEnd;
} tagReflnPropsDisplay;

typedef struct _tagSpotInfo
{
  int   nStat;
  float fIntensity;
  float fSigmaI;
  float fResolution;
  float fCentPx0;
  float fCentPx1;
  float fRotMid;
  float fRotWidth;
} tagSpotInfo;

//class CXGcontext;
class CXdisplay : public UIComponent {

public:

  Cimage        *m_poImage;
  CXcolormap    *m_poColormap;

//  CXGcontext  *m_poGcontext;
  GC             m_hGC;
  GC             m_hGCtext;
  GC             m_hGCrefln;
  XGCValues      m_hGCValues;
  Pixmap         m_hPixmap;
  XImage        *m_phXImage;
  XImage        *m_phXImageScale;
  XFontStruct   *m_pXFont;

  Widget         m_hParent;     // Should be a drawing area widget
  Dimension      m_nWidth;      // Width of drawing area
  Dimension      m_nHeight;     // Height of drawing area
  Dimension      m_nWidthUsed;  // Width used for image  (depends on zoom)
  Dimension      m_nHeightUsed; // Height used for image (depends on zoom)
  Visual        *m_pVisual;     // Visual of window in X Windows sense
  unsigned int   m_ulDepthImage;// Depth of image in X Windows sense
  int            m_nBytesPerPixel;  // Bytes per pixel for XImage's
    
  unsigned char *m_pucData;  // Pointer to image data converted to index in CLUT
  int            m_nDataSize;
  unsigned char *m_pucScale; // Pointer to image data converted to index in CLUT
  int            m_nScaleSize;
  int            m_nScaleHeight;
  int            m_nScaleWidth;
  int            m_nScaleMode;
  float          m_fScaleMin;
  float          m_fScaleMax;

  int     m_nRecurseCount;
  int     m_nShadowPixCount;


  bool           m_bLineInRefln;
  bool           m_bPixelSwapBytes;
  bool           m_bHasFloatValues;

  tagImagePropsDisplay  m_tImageProps;
  tagReflnPropsDisplay  m_tReflnProps;
  tagSpotInfo           m_tSpotInfo;

  int            m_nMin, m_nMax;
  int            m_nStart0, m_nStart1;
  float          m_fStep0,  m_fStep1;
  int            m_nDir0,   m_nDir1;
  int            m_nLineWidth;
  float          m_fLineWidth;
  CXrubberband  *m_poRubberband;
  CXcursor      *m_poCursor;
  Creflnlist    *m_poReflnlist;
  int            m_nXPrev, m_nYPrev; // Previous values of cursor position
  int            m_nCursorMode;
  int            m_nBeamstopShadowMode;
  unsigned long int m_ulSatColorIndex; // Color table index of a saturated pixel
  unsigned long int m_ulTooLowColorIndex; // Color table index of a too low pixel
  unsigned long int m_ulColorBlue;
  unsigned long int m_ulColorWhite;

  float          m_fPx0Prev, m_fPx1Prev;
  float          m_a5fQuad[5][2];

  Cdetector     *m_poDetector;
  float          m_fWavelength;
  float          m_a3fS0[3];
  Cstring        m_a10sReflnColors[10];

  XtPointer       m_pObj;          // For callback, pointer to an object
  void          (*m_prvRubberbandCallback)(XtPointer, const int nMode,
					   const float fPx0, const float fPx1,
					   const float fPx2, const float fPx3);

  UICallbackStruct *_clientDataStructs;
////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

  CXdisplay (char *, Widget, Cimage *, CXcolormap *, 
	     tagImagePropsDisplay *ptPropsIn = NULL,
	     tagReflnPropsDisplay *ptReflnPropsIn = NULL);
  ~CXdisplay ();
////////////////////////////////////////////////////////////////////////
//+Member fFunction prototypes
 public:

int
nAvgSD(Cimage& oImage, const int nPxArea[4], 
       const float fExcludeMin,
       const float fExcludeMax,
       float *pfMin, float *pfMax,
       float *pfAvg, float *pfSD);

void vImageToCLUT(Cimage& oImage, const float fMin, const float fMax,
		  const int nMin, const int nMax);
//

void vEventCallback(Widget, XtPointer client_data, XtPointer call_data);

int nDpyPixToImgPix(const int nX, const int nY, float *pf1, float *pf2);
int nImgPixToDpyPix(const float fPx0In, const float fPx1In, int *pnX, int *pnY);
int nGetImageProps(tagImagePropsDisplay *ptOutputProps);
int nSetImageProps(const tagImagePropsDisplay *ptInputProps);
int nGetReflnProps(tagReflnPropsDisplay *ptOutputProps);
int nSetReflnProps(const tagReflnPropsDisplay *ptInputProps, bool bRefresh=TRUE);
void vRefresh(void);
void vPan(const int nX, const int nY);
int nSetRegion(const int nOrigIn0, const int nOrigIn1, 
		const int nExtIn0,  const int nExtIn1);
void vPlotReflnlist(Creflnlist *poReflnlist = NULL);
void vPlotPixelValues(void);
void vPlotDrawObjects(void);
int  nCreatePS(const Cstring& sPSFile, const Cstring& sComment="", 
	       const Boolean& bPrintName=True);

inline void  vSetWavelength(const float fWavelength)
  {m_fWavelength = fWavelength;}
inline void  vSetSatPixValue(const float fSat)
  {m_tImageProps.fSatPixValue = fSat;}
inline void  vSetTooLowPixValue(const float fTooLow)
  {m_tImageProps.fTooLowPixValue = fTooLow;}
inline float fGetSatPixValue(void) {return (m_tImageProps.fSatPixValue); }
inline float fGetTooLowPixValue(void) {return (m_tImageProps.fTooLowPixValue); }
void vSetSatPixColor(const char *pcColor);
void vSetTooLowPixColor(const char *pcColor);

inline void vSetZoom(const float fZoom = 1.0) { m_tImageProps.fZoom = fZoom; }
void vSetResol(Cdetector *poDetector, const float *pfS0, 
	       const float fWavelength, const int nRings=5);

void vSetRefCursorColor(const char *pcColor);
void vSetSpotInfo(const float fIntensity, const float fSigmaI, 
                  const float fResolution, const float fPx0, const float fPX1,
                  const float fRotMid, const float fRotWidth);

 inline void vSetCursorMode(const int nMode) {m_nCursorMode = nMode; }
 void vSetLineWidth(const float fLineWidth );
 int  nCalcDrawBeamstopShadow(float fPixIn0, float fPixIn1, float fMaxValue = 0.0f);
 int  nCDBeamstopShadow(void);
 float fLineShadow(float fPx0, float fPx1, float fSlope, float fIntercept, 
                  float fWidth, int nNumPoints, int nSmooth = 0, int nDir = 0);
void vSetBeamstopShadowMode(const int nMode);
 inline int  nGetBeamstopShadowMode(void) {return (m_nBeamstopShadowMode); }

private:

 bool m_bSpacebarHeldDown;
 bool m_bControlKeyHeldDown;

static void vEventCallbackCallback(Widget, XtPointer, XtPointer);

void vDeleteXStuff(void);
int  nCreatePixmap(void);
void vPixelInfo(const int nX=-1, const int nY=-1, const int nInfoMode=0,
		const int nButton=1);

void vDrawProfile(const int nX=-1, const int nY=-1, const int nInfoMode=0);
void vDrawColorScale(const int nMode=0, const int nX=0, const int nY=0);
void vDrawResol(void);
int  nCalcEllipse(float, float, float, float, float, int, XPoint*);
int  nRecurseShadow(int nCent0, int nCent1, int nEdgeSize, float fAvg, float fSD,
		    float fThreshold, float fSlopePixVal, 
		    Cimage *poImageOut);
void vSwapBytes(const int nNumChars, unsigned char *pucArray);

};  // end of class CXdisplay

#endif   // DT_CXDISPLAY_H
