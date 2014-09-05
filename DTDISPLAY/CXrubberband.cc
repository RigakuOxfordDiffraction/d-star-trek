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
// CXrubberband.cc       Initial author: J.W. Pflugrath           18-Jul-1995
//    This file contains the member functions of class CXrubberband.
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

#include "CXrubberband.h"              // Class definition and prototypes

//+Code begin

//+Definitions, constants, and initialization of static member variables

#define XRUB_LINE_WIDTH 1
//+Public functions

// Constructors, destructors and assignments

CXrubberband::CXrubberband(Widget hParentIn, Pixmap *phPixmapIn, 
			   const int nStartX, const int nStartY,
			   const int nWidthMax, const int nHeightMax,
			   const int nModeIn)
{
  // Construct a CXrubberband object

  m_hParent         = hParentIn;
  m_phParentPixmap  = phPixmapIn;  // Use a pointer, since parent may change
  m_nColorAllocated = 0;           // pixmap out from under us

  // Initialize starting points of the rectangle, circle or line

  m_nMode    = nModeIn;
  m_nOrigX   = nStartX;
  m_nOrigY   = nStartY;
  m_nEndX    = m_nOrigX;
  m_nEndY    = m_nOrigY;
  m_nHeight  = 1;
  m_nWidth   = 1;
  m_nDrawX   = m_nOrigX;
  m_nDrawY   = m_nOrigY;
  m_nLineWidth = 0;

  // Get width and height of parent widget

  XtVaGetValues(m_hParent, XmNwidth, &m_nParentWidth,
		XmNheight, &m_nParentHeight, NULL);

  if ( (0 < nWidthMax) && (m_nParentWidth > nWidthMax) )
    m_nParentWidth = nWidthMax;
  if ( (0 < nHeightMax) && (m_nParentHeight > nHeightMax) )
    m_nParentHeight = nHeightMax;

  // Get a graphic context.  Replace this later with a Cgc class.

  XGCValues hXGCv;
  hXGCv.line_width = XRUB_LINE_WIDTH;
  if (m_nColorAllocated == 1)
    hXGCv.foreground = m_hColor.pixel;
  else
    hXGCv.foreground = BlackPixel(XtDisplay(m_hParent),
				  DefaultScreen(XtDisplay(m_hParent)));
  hXGCv.function = GXorReverse;
  m_hGC     = XCreateGC(XtDisplay(m_hParent),
			XtWindow(m_hParent), 
			GCLineWidth | GCFunction | GCForeground, &hXGCv);

  m_hGCValues.line_width = XRUB_LINE_WIDTH;
  // Change cursor?
}

CXrubberband::~CXrubberband() 
{
  vRestoreOutline();                

  XSync(XtDisplay(m_hParent), False); // Be sure to flush and process all events
  XmUpdateDisplay(m_hParent);         // that this object may have generated
                                      // before deleting the object

  if (NULL != m_hGC)
    XFreeGC(XtDisplay(m_hParent), m_hGC);
}

//+Public functions

void
CXrubberband::vDraw(const int nXIn, const int nYIn, const bool bRestore)
{

  // Restore previous rubberbanded outline 

  if (bRestore) vRestoreOutline();

  if ( (RB_TLIN_MODE != m_nMode)
       || (XRUB_LINE_WIDTH != m_hGCValues.line_width) )
    {
      m_hGCValues.line_width = XRUB_LINE_WIDTH;
      m_hGCValues.cap_style  = CapButt;
      XChangeGC(XtDisplay(m_hParent), m_hGC, GCCapStyle | GCLineWidth, &m_hGCValues);
      //      XChangeGC(XtDisplay(m_hParent), m_hGC, GCLineWidth, &m_hGCValues);
    }

  m_nEndX = nXIn;
  m_nEndY = nYIn;
  if (0 > m_nEndX) m_nEndX = 0;
  if (0 > m_nEndY) m_nEndY = 0;
  if (m_nEndX >= m_nParentWidth)  m_nEndX = m_nParentWidth-1;
  if (m_nEndY >= m_nParentHeight) m_nEndY = m_nParentHeight-1;
  
  if (RB_RECT_MODE == m_nMode)
    {
      // Draw a rectangle on the WINDOW, not on the PIXMAP!!

      m_nDrawX  = m_nOrigX;
      m_nDrawY  = m_nOrigY;
      m_nWidth  = m_nEndX - m_nOrigX;
      m_nHeight = m_nEndY - m_nOrigY;
      if (0 > m_nWidth)
	{
	  m_nWidth = -m_nWidth;
	  m_nDrawX = m_nEndX;
	}
      else if (0 == m_nWidth)
	{
	  m_nWidth = 1;
	}
      if (0 > m_nHeight)
	{
	  m_nHeight = -m_nHeight;
	  m_nDrawY  = m_nEndY;
	}
      else if (0 == m_nHeight)
	{
	  m_nHeight = 1;
	}

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nOrigX, m_nEndY);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nEndY, m_nEndX, m_nEndY);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nEndX, m_nEndY, m_nEndX, m_nOrigY);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nEndX, m_nOrigY, m_nOrigX, m_nOrigY);
    }
  else if ( (RB_TLIN_MODE == m_nMode) || (RB_LINE_MODE == m_nMode) )
    {
      // Draw a special THICK line

      m_nDrawX  = m_nOrigX;
      m_nDrawY  = m_nOrigY;
      m_nWidth  = m_nEndX - m_nOrigX;
      m_nHeight = m_nEndY - m_nOrigY;
      if (0 > m_nWidth)
	{
	  m_nWidth = -m_nWidth;
	  m_nDrawX = m_nEndX;
	}
      else if (0 == m_nWidth)
	{
	  m_nWidth = 1;
	}
      if (0 > m_nHeight)
	{
	  m_nHeight = -m_nHeight;
	  m_nDrawY  = m_nEndY;
	}
      else if (0 == m_nHeight)
	{
	  m_nHeight = 1;
	}

      // Set line thickness

      //std::cout << "Trying to set LINE WIDTH to: " << m_nLineWidth << std::endl;
      //std::cout << "m_nCircleDiam, m_hGCValues.line_width: " << m_nLineWidth << ", "
      // << m_hGCValues.line_width << std::endl;

      if (m_nLineWidth != m_hGCValues.line_width)
	{
	  //std::cout << "Setting LINE WIDTH to: " << m_nLineWidth << std::endl;
	  m_hGCValues.line_width  = m_nLineWidth;
	  m_hGCValues.cap_style   = CapRound;
	  //XChangeGC(XtDisplay(m_hParent), m_hGC, GCLineWidth, &m_hGCValues);
	  XChangeGC(XtDisplay(m_hParent), m_hGC, GCLineWidth | GCCapStyle, &m_hGCValues);
	}
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nEndX, m_nEndY);
    }
  else if (RB_LINE_MODE == m_nMode)
    {
      // Draw a single line

      m_nDrawX  = m_nOrigX;
      m_nDrawY  = m_nOrigY;
      m_nWidth  = m_nEndX - m_nOrigX;
      m_nHeight = m_nEndY - m_nOrigY;
      if (0 > m_nWidth)
	{
	  m_nWidth = -m_nWidth;
	  m_nDrawX = m_nEndX;
	}
      else if (0 == m_nWidth)
	{
	  m_nWidth = 1;
	}
      if (0 > m_nHeight)
	{
	  m_nHeight = -m_nHeight;
	  m_nDrawY  = m_nEndY;
	}
      else if (0 == m_nHeight)
	{
	  m_nHeight = 1;
	}

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nEndX, m_nEndY);
    }
  else if (RB_CIRC_MODE == m_nMode)
    {
//+2011-03-02 JWP
// Draw two circles, the original one and one at half the diameter
//-2011-03-02 JWP
      // Draw a circle on the WINDOW, not on the PIXMAP!!
      // Center of circle is original coordinate, 
      // perimeter is at current coordinate
      
      int nRadius;
      m_nDrawX  = m_nOrigX - m_nEndX;
      m_nDrawY  = m_nOrigY - m_nEndY;
      nRadius = (int) sqrt( (double) (   m_nDrawX * m_nDrawX 
				      +  m_nDrawY * m_nDrawY) );
      if (0 == nRadius)
	nRadius = 1;
      m_nHeight = 2 * nRadius;
      // m_nCircleDiam = m_nHeight;
      //cout << "Circle diam: " << m_nCircleDiam << endl;
      m_nWidth  = m_nHeight;
      m_nDrawX  = m_nOrigX - nRadius;
      m_nDrawY  = m_nOrigY - nRadius;
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX-m_nWidth/2, m_nOrigY,
		m_nOrigX+m_nWidth/2, m_nOrigY);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY-m_nHeight/2, 
		m_nOrigX, m_nOrigY+m_nHeight/2);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nEndX, m_nEndY);

      XDrawArc(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
	       m_nDrawX, m_nDrawY, m_nWidth, m_nHeight, 0, 23040);
//+2011-03-02 JWP
      m_nDrawX2  = m_nOrigX - nRadius/2;
      m_nDrawY2  = m_nOrigY - nRadius/2;
      XDrawArc(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
	       m_nDrawX2, m_nDrawY2, m_nWidth/2, m_nHeight/2, 0, 23040);
//-2011-03-02 JWP
    }
  else if (RB_CIR2_MODE == m_nMode)
    {
      // Draw a circle on the WINDOW, not on the PIXMAP!!
      // Center of circle is at current coordinate, 
      // perimeter is at original coordinate (possibly panned)
      
      int nRadius;

      m_nDrawX  = m_nOrigX - m_nEndX;
      m_nDrawY  = m_nOrigY - m_nEndY;
      nRadius = (int) sqrt( (double) (   m_nDrawX * m_nDrawX 
				      +  m_nDrawY * m_nDrawY) );
      if (0 == nRadius)
	nRadius = 1;
      m_nHeight = 2 * nRadius;
      m_nWidth  = m_nHeight;
      m_nDrawX  = m_nEndX - nRadius;
      m_nDrawY  = m_nEndY - nRadius;
//      m_nDrawX  = m_nOrigX;
//      m_nDrawY  = m_nOrigY;
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nEndX-m_nWidth/2, m_nEndY,
		m_nEndX+m_nWidth/2, m_nEndY);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nEndX, m_nEndY-m_nHeight/2, 
		m_nEndX, m_nEndY+m_nHeight/2);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nEndX, m_nEndY);

      XDrawArc(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
	       m_nDrawX, m_nDrawY, m_nWidth, m_nHeight, 0, 23040);
    }
};

void
CXrubberband::vGetCoords (int *pnStartX, int *pnStartY, int *pnEndX, 
			  int *pnEndY)
{
  // Return the box that was rubberbanded.  Start is guaranteed to
  // be less than the end;

  if (RB_RECT_MODE == m_nMode)
    {
      *pnStartX = min(m_nOrigX, m_nEndX);
      *pnStartY = min(m_nOrigY, m_nEndY);
      *pnEndX   = max(m_nOrigX, m_nEndX);
      *pnEndY   = max(m_nOrigY, m_nEndY);
    }
  else if (RB_LINE_MODE == m_nMode)
    {
      *pnStartX = m_nOrigX;
      *pnStartY = m_nOrigY;
      *pnEndX   = m_nEndX;
      *pnEndY   = m_nEndY;
    }
  else if (RB_TLIN_MODE == m_nMode)
    {
      *pnStartX = m_nOrigX;
      *pnStartY = m_nOrigY;
      *pnEndX   = m_nEndX;
      *pnEndY   = m_nEndY;
    }
  else if (RB_CIRC_MODE == m_nMode)
    {
      *pnStartX = m_nOrigX;
      *pnStartY = m_nOrigY;
      *pnEndX   = m_nEndX;
      *pnEndY   = m_nEndY;
    }
  else if (RB_CIR2_MODE == m_nMode)
    {
      *pnStartX = m_nOrigX;
      *pnStartY = m_nOrigY;
      *pnEndX   = m_nEndX;
      *pnEndY   = m_nEndY;
    }
};

int
CXrubberband::nSetLineColor(const unsigned int nRed, 
		      const unsigned int nGreen, const unsigned int nBlue)
{
  // First check of hColor has an allocated color.

  if (m_nColorAllocated != 0) 
    {
      m_hColor.red   = nRed;
      m_hColor.green = nGreen;
      m_hColor.blue  = nBlue;
//      XStoreColor(XtDisplay(m_hParent),
//		  m_hParentCmap,
//		  &m_hColor);
      return (0);
    }
  else
    return (-1);

}

int
CXrubberband::nSetLineColor(const char *cColorName)
{
  // Set the color of the rubberband outline to the input
  // color specified by char string in cColorName.  The string must match
  // a color in the X servers color database or an error is returned.

  XColor hXcolorExact, hXcolorScreen;
/*
  Status nStat = XLookupColor(XtDisplay(m_hParent),
			      m_hParentCmap,
			      cColorName, 
			      &hXcolorExact, &hXcolorScreen);
  if (nStat != 0) 
    {
      return( nSetLineColor(hXcolorExact.red,
			    hXcolorExact.green, hXcolorExact.blue));
    }
  else
    {
      return (-2);
    }
*/
  return (-1);
}

void
CXrubberband::vRestoreOutline(void)
{
  // Restore the pixmap displayed in the window by 
  // copying from the offscreen pixmap to the window the lines/area that were
  // previously drawn in vDraw;

  if (RB_RECT_MODE == m_nMode)
    {
      // For XCopyArea, the width and height must be positive.
      // Left column

      XCopyArea (XtDisplay(m_hParent), *m_phParentPixmap, XtWindow(m_hParent),
		 m_hGC, 
		 m_nDrawX, m_nDrawY,
		 1, m_nHeight,
		 m_nDrawX, m_nDrawY);

      // Right column

      XCopyArea (XtDisplay(m_hParent), *m_phParentPixmap, XtWindow(m_hParent),
		 m_hGC, 
		 m_nDrawX+m_nWidth, m_nDrawY,
		 1, m_nHeight,
		 m_nDrawX+m_nWidth, m_nDrawY);

      // Top row

      XCopyArea (XtDisplay(m_hParent), *m_phParentPixmap, XtWindow(m_hParent),
		 m_hGC, 
		 m_nDrawX, m_nDrawY,
		 m_nWidth, 1,
		 m_nDrawX, m_nDrawY);

      // Bottom row

      XCopyArea (XtDisplay(m_hParent), *m_phParentPixmap, XtWindow(m_hParent),
		 m_hGC, 
		 m_nDrawX, m_nDrawY+m_nHeight,
		 m_nWidth, 1,
		 m_nDrawX, m_nDrawY+m_nHeight);
    }
  else if (RB_LINE_MODE == m_nMode)
    {
      // For XCopyArea, the width and height must be positive.

      // Assumes m_hGC is XOrReverve
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nEndX, m_nEndY);
    }
  else if (RB_TLIN_MODE == m_nMode)
    {
      // Thick line mode

      // Assumes m_hGC is XOrReverve
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nEndX, m_nEndY);
    }
  else if (RB_CIRC_MODE == m_nMode)
    {
      // Assumes m_hGC is XOrReverve
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX-m_nWidth/2, m_nOrigY,
		m_nOrigX+m_nWidth/2, m_nOrigY);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY-m_nHeight/2, 
		m_nOrigX, m_nOrigY+m_nHeight/2);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nOrigX, m_nOrigY, m_nEndX, m_nEndY);

      XDrawArc(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
	       m_nDrawX, m_nDrawY, m_nWidth, m_nHeight, 0, 23040);
//+2011-03-02 JWP
      XDrawArc(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
	       m_nDrawX2, m_nDrawY2, m_nWidth/2, m_nHeight/2, 0, 23040);
//-2011-03-02 JWP
    }

  else if (RB_CIR2_MODE == m_nMode)
    {
      // Assumes m_hGC is XOrReverve

      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nEndX-m_nWidth/2, m_nEndY,
		m_nEndX+m_nWidth/2, m_nEndY);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nEndX, m_nEndY-m_nHeight/2, 
		m_nEndX, m_nEndY+m_nHeight/2);
      XDrawLine(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
		m_nEndX, m_nEndY, m_nOrigX, m_nOrigY);

      XDrawArc(XtDisplay(m_hParent), XtWindow(m_hParent), m_hGC,
	       m_nDrawX, m_nDrawY, m_nWidth, m_nHeight, 0, 23040);
    }

};

void
CXrubberband::vPan(const int nXIn, const int nYIn)
{
  // Move the rubberband object around without changing its size
 
  // Restore previous rubberbanded outline 

  vRestoreOutline();

  // Compute shift of old versus new position

  int nDeltaX = m_nEndX - nXIn;
  int nDeltaY = m_nEndY - nYIn;
  
  // Add shift to old coordinates and re-draw

  m_nEndX  = nXIn;
  m_nEndY  = nYIn;
  m_nOrigX = m_nOrigX - nDeltaX;
  m_nOrigY = m_nOrigY - nDeltaY;
  vDraw(m_nEndX, m_nEndY, False);
}
