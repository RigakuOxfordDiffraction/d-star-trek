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
// CXcolormap.cc            Initial author: J.W. Pflugrath           18-Jul-1995
//    This file contains the member functions of class CXcolormap.
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

#include "CXcolormap.h"              // Class definition and prototypes
#include "dtreksys.h"                // Has sGetEnv prototype
using namespace std;

//+Code begin

//+Definitions, constants, and initialization of static member variables

const unsigned long   CXcolormap::ms_ulCOLORERROR = 65536;

//+Public functions


// Constructors, destructors and assignments

CXcolormap::CXcolormap(const char *pcWidgetName, Widget hParentIn, 
		       const unsigned int unRequestedIn, 
		       const int nTolerance) : UIComponent(pcWidgetName)
{
  // Construct a CXcolormap object
  // Get a colormap of up to nRequested colors.
  // If the display supports less colors, get a colormap anyways and 
  // let the calling program think that it has all the requested colors.
  // See D. A. Young (1990) "The X Window System Programming and Applications
  // with Xt OSF/Motif Edition", Prentice-Hall pp.179-181.

  unsigned int unRequested;
  int i;
  int nStat;

  unRequested     = unRequestedIn;
  m_hParent       = hParentIn;
  m_hColormap     = (Colormap) NULL;
  m_pulColorIndex = NULL;
  m_pulColorAlloc = NULL;
  m_pnColorUsed   = NULL;
  m_phXcolor      = NULL;
  m_nNumRamp      = 1;
  m_nNumAllocated = 0;
  m_bNotDefault   = FALSE;
  m_nExtra        = 9;      // Try to allocate nine extra colors
  m_nTolerance    = nTolerance;
  m_nToleranceExtra = 1000;

  if ("" != sGetEnv("DTDISPLAY_COLOR_TOLEXTRA"))
    {
      m_nToleranceExtra = atoi(sGetEnv("DTDISPLAY_COLOR_TOLEXTRA").string());
      cout << "INFO: Color match extra tolerance set to " 
           << m_nToleranceExtra << endl;
    }

  // Get the default visual and if PseudoColor or TrueColor do different
  // things.  If neither of these two, then you are in trouble.

  m_phDefaultVisual = DefaultVisual(XtDisplay(m_hParent),
				    DefaultScreen(XtDisplay(m_hParent)));

  // Note: c_class is for C++ only!


  if (PseudoColor == m_phDefaultVisual->c_class)
    {
      //cout << "PseudoColor visual!\n";
    }
  else if (TrueColor == m_phDefaultVisual->c_class)
    {
      //cout << "TrueColor visual!\n";
    }
  else
    {
      cout << "FATAL ERROR: CXcolormap works only with "
              "PseudoColor and TrueColor default visuals!\n" << flush;
      exit (2);
    }

  // Start with the default colormap.

  m_hColormap = DefaultColormapOfScreen(XtScreen(m_hParent));

  // Get the number of available colors, this does not tell us the
  // true number of colors available, especially if visual class is TrueColor

  m_nAvail = XDisplayCells(XtDisplay(m_hParent),
			   DefaultScreen(XtDisplay(m_hParent)));

  // m_nAllocated is the number allocated in the arrays of this object.
  // The idea is to allocate the number of requested colors whether the
  // display supports them or not, then set them up with indices into the
  // display colortable.

  if ("" != sGetEnv("DTDISPLAY_NUMGRAYLEVELS"))
    {
      unRequested = (unsigned int) atoi(sGetEnv("DTDISPLAY_NUMGRAYLEVELS").string());
      unRequested = min(128, (int)unRequested);
      unRequested = max(1, (int)unRequested);      
      cout << "INFO: Num requested gray levels set to min(128," << unRequested 
           << ")\n"
           << "      by env var DTDISPLAY_NUMGRAYLEVELS.\n";
    }
  m_nAllocated    = unRequested + m_nExtra;

  if (TrueColor == m_phDefaultVisual->c_class)
    {
      // Don't worry about m_nAvail

      m_nAvail        = max(256, m_nAllocated);
    }

  m_pulColorIndex = new unsigned long [m_nAllocated];
  m_pulColorAlloc = new unsigned long [m_nAllocated];
  m_pnColorUsed   = new int [m_nAvail];
  m_phXcolor      = new XColor [m_nAllocated];

  for (i = 0; i < m_nAvail; i++)
    {
      m_pnColorUsed[i] = 1;        // Mark all color cells as used for now
    }

  if (m_nAvail < unRequested) 
    {
      //cout << "m_nAvail, unRequested: " << m_nAvail << ", " << unRequested
      //     << endl;
      
      cout << "WARNING: display does not support the requested number of colors!\n";

      // Alternate black and white pixel into our m_pulColorIndex array

      for (i = 0; i < m_nAllocated; i++)
	{
	  if (0 == (i % 2))
	    {
	      m_pulColorIndex[i] = BlackPixel(XtDisplay(m_hParent), 
					   DefaultScreen(XtDisplay(m_hParent)));
	    }
	  else
	    {
	      m_pulColorIndex[i] = WhitePixel(XtDisplay(m_hParent), 
					   DefaultScreen(XtDisplay(m_hParent)));
	    }
	  m_pulColorAlloc[i] = ms_ulCOLORERROR;
	}
//+jwp 26-Jun-1999 do not return
//      return;
//-jwp 26-Jun-1999
    }
  else
    {
      m_nNumAllocated = 0;
    }

//+jwp 30-Jul-1999  

  // Want TRUE for TrueColor mode and for first attempt with PseudoColor
  // unless over-ridden with DTDISPLAY_NOTREADONLY envvar

  m_bReadOnlyMode   = TRUE;
  if (   ("" != sGetEnv("DTDISPLAY_NOTREADONLY"))
      && (PseudoColor == m_phDefaultVisual->c_class) )
    m_bReadOnlyMode = FALSE;


  // Some stuff with read only color cells

  if (m_bReadOnlyMode)
    {
      // Try to allocate colors in read-only cells

      nStat = nLoadGrey(unRequested);
      //cout << "1st nLoadGrey nStat: " << nStat << endl;
      if (0 == nStat)
	{
	  // Getting all the read-only colors worked, so just quit while
	  // we are ahead

	  return;
	}
      else if (0 < m_nNumAllocated)
	{
	  //cout << "just before 1XFreeColors: " << endl << flush;
	  XFreeColors(XtDisplay(m_hParent), m_hColormap, m_pulColorAlloc,
		      m_nNumAllocated, 0);
	  m_nNumAllocated = 0;
	  XFlush(XtDisplay(m_hParent));

	  //cout << "just after 1XFreeColors: " << endl << flush;
	}
      else
	{
	  cout << "ERROR - programmer error in CXcolormap.\n" << flush;
	  exit (1);
	}
    }

  // OK, readonly failed, so try to create a colormap and allocate read/write
  // color cells.  
  // The following code only works with a PseudoColor visual,
  // so check if visual is OK, if not EXIT.

  if (PseudoColor != m_phDefaultVisual->c_class)
    {
      cout << "FATAL ERROR: CXcolormap not readonly works only with "
              "PseudoColor default visual!\n" << flush;
      exit (2);
    }

  m_bReadOnlyMode = FALSE;

//-jwp 30-Jul-1999

  // m_pulColorIndex will contain an index into the colormap, the range of the
  // index will be between 0 <= index < m_nAvail, so we can fake that we have 
  // all requested colors.
  
  // Try to allocate the required read/write color cells in the default colormap
  // The indices (.pixel values) of the allocated cells are returned in
  // m_pulColorAlloc.
  
  m_nNumAllocated = min(m_nAvail, m_nAllocated);
  nStat = XAllocColorCells(XtDisplay(m_hParent), m_hColormap, FALSE, NULL, 0,
			   m_pulColorAlloc, m_nNumAllocated);

  // Normally nStat is 1 (non-zero indicates SUCCESS)

  //cout << "nStat from XAllocColorCells in CXcolormap: " << nStat << endl;
//+jwp 11-Aug-1997
// Force using our own colormap...
// If you do not do this, then IRIX 6.3 does not do the correct colormap,
// BUT if you DO do this, then an xwd of dtdisplay does not have right
// colors either.
/*  There should be a #ifdef around this piece of code
    // Probably should use XCopyColormapAndFree (see p. 416)

  if (0 != nStat)
    {
      // Free the colors we got

      XFreeColors(XtDisplay(m_hParent), m_hColormap, m_pulColorAlloc,
		  m_nNumAllocated, 0);
    }
  // Set status to 0 to force creating a new colormap below...
  nStat = 0;  
*/
//-jwp 11-Aug-1997

  if (0 != nStat)
    {
      // XAllocColorCells was successful

      // Mark all allocated cells as usable
      for (i = 0; i < m_nAllocated; i++)
	{
	  m_pnColorUsed[m_pulColorAlloc[i]] = 0;  // Mark the allocated cells as
	}                                         // unused
    }
  else
    {
      // Error (XAllocColorCells failed), so get new colormap to work with
      // This can only work if the server supports more than one colormap.
      // Since many servers support only a single colormap, the following will
      // fail!

      m_bNotDefault = TRUE;

      // Get the colors of the default colormap and put them in the new colormap
      // This should reduce colormap flashing.

      delete [] m_phXcolor;
      delete [] m_pulColorIndex;
      delete [] m_pulColorAlloc;
      delete [] m_pnColorUsed;
      m_nAllocated    = m_nAvail;
      m_pulColorIndex = new unsigned long [m_nAllocated];
      m_pulColorAlloc = new unsigned long [m_nAllocated];
      m_pnColorUsed   = new int [m_nAllocated];
      m_phXcolor      = new XColor [m_nAllocated];
      for (i = 0; i < m_nAllocated; i++)
	{
	  m_pulColorIndex[i]  = m_nAllocated - i - 1;  // Load backwards
	  m_pulColorAlloc[i]  = m_pulColorIndex[i];
	  m_pnColorUsed[m_pulColorAlloc[i]]    = 0;    // This color cell is 
	  m_phXcolor[i].pixel = m_pulColorIndex[i];    // not used
	  m_phXcolor[i].flags = DoRed | DoGreen | DoBlue;
	}

      // Get all the colors in the current colormap 

      //cout << "m_nAvail, m_nAllocated: " << m_nAvail << ", " << m_nAllocated
      //     << endl << flush;
    
      XQueryColors(XtDisplay(m_hParent), m_hColormap, m_phXcolor, 
		   m_nNumAllocated);

      // Create the new colormap and allocate all the cells

      m_hColormap = XCreateColormap(XtDisplay(m_hParent), 
				    XtWindow(m_hParent),
				    DefaultVisual(XtDisplay(m_hParent),
		         		  DefaultScreen(XtDisplay(m_hParent))),
				    AllocAll);

      // Store the colors of the default colormap (which we just XQueried into
      // m_phXcolor) into the new colormap.

      m_nNumAllocated = min(m_nAllocated, m_nAvail);
      XStoreColors(XtDisplay(m_hParent), m_hColormap, m_phXcolor, 
		   m_nNumAllocated);
      
      // Let the window manager install the colormap
      // (This should NOT be called?)

      XInstallColormap(XtDisplay(m_hParent), m_hColormap);

      // Now set colormap property on a few windows
      // Climb the widget until we get to the top. From Vol. 6A, p. 211.

      Widget wTop, wLast;
      wTop = m_hParent;
      Cardinal nNumWidgets = 0;
      Widget a100WidgetList[102];
      a100WidgetList[0] = wTop;
      while (wTop  && !XtIsWMShell(wTop))
	{
//	  vSetWindow(wTop);       // Maybe move this outside of while
	  wLast = wTop;
	  wTop = XtParent(wTop);
	  if (100 > nNumWidgets)
	    {
	      nNumWidgets++;
	      a100WidgetList[nNumWidgets] = wTop;
	    }
	}
      vSetWindow(wLast);
      // See p. 417-420 of "Professional Graphics Programming in the
      // X Window System" by Eric Johnson & Kevin Richard 

      //cout << "numcolormap w: " << nNumWidgets << endl;
      XtSetWMColormapWindows(m_hParent, a100WidgetList, nNumWidgets);
     }
  nStat = nLoadGrey(unRequested);
}

CXcolormap::~CXcolormap() 
{

  if (m_bNotDefault)
    {
//      XFreeColors(XtDisplay(m_hParent), m_hColormap, m_pulColorAlloc,
//		  m_nAllocated, 0);
      XFreeColormap(XtDisplay(m_hParent), m_hColormap);
    }
  else
    {
      XFreeColors(XtDisplay(m_hParent), m_hColormap, m_pulColorAlloc,
		  m_nNumAllocated, 0);
    }
  if (NULL != m_pulColorIndex)
    {
      delete [] m_pulColorIndex;
      m_pulColorIndex = NULL;
    }

  if (NULL != m_pulColorAlloc)
    {
      delete [] m_pulColorAlloc;
      m_pulColorAlloc = NULL;
    }

  if (NULL != m_pnColorUsed)
    {
      delete [] m_pnColorUsed;
      m_pnColorUsed = NULL;
    }

  if (NULL != m_phXcolor)
    {
      delete [] m_phXcolor;
      m_phXcolor = NULL;
    }
}

//+Public functions

int
CXcolormap::nLoadGrey(const int nNumGreyIn)
{
  // Load greyscale into the colormap.  The problem is that the colormap
  // has at most m_nAvail cells.  We want the calling program to think that
  // there are at least m_nNumAllocated cells.  The calling program will 
  // use the array
  // m_pulColorIndex[0->(m_nNumAllocated-1)] to index the colorcells. 
  // For color displays,
  // there should be one-to-one correspondence, but for black&white displays
  // there will be a difference.  Also if there are enough extra colorcells
  // left over, go ahead and set some primary colors.

  int i;
  int nNumGrey;
  int nStat;
      
  m_nNumRamp = nNumGreyIn;

  nNumGrey       = nNumGreyIn;
  if (nNumGrey > m_nAvail)
    return (-1);                // We cannot do this, so return with error

  int nIncrement = -65535 / (nNumGrey-1);

  int nStart     = 65535;
  int nGrey      = nStart;

  nStat = 0;
  for (i = 0; (i < nNumGrey) && (0 == nStat); i++)
    {
      nStat = nSetColor(i, nGrey, nGrey, nGrey, m_nTolerance);
      nGrey = nGrey + nIncrement;
    }

  // Load extra colors if possible, note that nSetColor will not exceed the
  // limits of the color table

  if (0 == nStat) nStat = nSetColor(i++, "black", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "red", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "orange", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "yellow", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "blue", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "green", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "cyan", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "magenta", m_nTolerance);
  if (0 == nStat) nStat = nSetColor(i++, "white", m_nTolerance);

  return (0);
}

int
CXcolormap::nSetColor(const int nIndex, const unsigned int unRed, 
		      const unsigned int unGreen, const unsigned int unBlue,
		      const int nTolerance)
{
  // Set the color of an allocated cell in the colormap to the input
  // color specified by RGB values if the color does not already exist
  // in the colormap.  If it does exist in an allocated part of the colormap,
  // then change the value of m_pulColorIndex[nIndex] to the index of the cell
  // in the colormap.

  // Problems: we need a way to mark color index i as already used, so that
  // we do not store a color into that index.  This is the only routine that
  // should call the XStoreColor routine, so we have a chance here to test
  // and catch problems.

  int    i;
  XColor hXColor;
  int    nStat = 0;

  if ( (0 > nIndex) || (nIndex >= m_nAllocated) || (nIndex >= m_nAvail) ||
       (2 == m_nAvail)  )
    {
      return (-1);  // Index out of range of color table
    }
  else
    {
      hXColor.red   = (unsigned short) unRed;
      hXColor.green = (unsigned short) unGreen;
      hXColor.blue  = (unsigned short) unBlue;
      hXColor.flags = DoRed | DoGreen | DoBlue;
      hXColor.pixel = 0;
      if (m_bReadOnlyMode)
	{
	  // Here for ReadOnly colormap

	  // Danger, we do not want to re-allocate if it is already
	  // allocated!

	  nStat = XAllocColor(XtDisplay(m_hParent), m_hColormap, 
			      &hXColor);
	  if (0 != nStat)
	    { 
	      // Success!
	      // What to do here is not quite correct

	      m_pulColorAlloc[m_nNumAllocated] = hXColor.pixel;
	      m_pulColorIndex[nIndex]          = hXColor.pixel;
	      m_nNumAllocated++;
	      m_pnColorUsed[nIndex] = 1;
	      nStat = 0;
	    }
	  else
	    {
	      // Failure :(
	      cout << "ERROR in CXColormap:nSetColor() for nIndex: " 
                   << nIndex << endl;
	      cout << "INFO: consider reducing the number of gray levels by setting the\n"
                   << "      environment variable DTDISPLAY_NUMGRAYLEVELS to " << nIndex - 1 << endl;
	      nStat = nIndex;
	      m_pulColorIndex[nIndex] = 0;
	      m_pulColorAlloc[nIndex] = 0;
	    }
	}
      else
	{
	  hXColor.pixel = ulGetColorIndex(unRed, unGreen, unBlue, nTolerance);
	  if (ms_ulCOLORERROR != hXColor.pixel)
	    {
	      // Return index of existing color
	      m_pulColorIndex[nIndex] = hXColor.pixel;
	    }
	  else
	    {
	      // Here for Read/Write colormap
	      // Find an unused color cell to store the color in ...

	      for (i = 0; (i < m_nAllocated)
	                     && (0 != m_pnColorUsed[m_pulColorAlloc[i]]); i++) ;

	      if (i < m_nAllocated)
		{
		  // Found an unused cell

		  m_pulColorIndex[nIndex] = m_pulColorAlloc[i];
		  hXColor.pixel = m_pulColorIndex[nIndex];
		  XStoreColor(XtDisplay(m_hParent), m_hColormap, &hXColor);
		}
	      else
		{
		  // No unused colors, so set index to first color and return error
		  m_pulColorIndex[nIndex] = 0;
		  return (-1);
		}
	    }
	  // Mark this color cell as used...so we do not store a new color over it
	  m_pnColorUsed[m_pulColorIndex[nIndex]] = 1;
	}
    }
  return (nStat);
}

int
CXcolormap::nSetColor(const int nIndex, const char *cColorName,
		      const int nTolerance)
{
  // Set the color of an allocated cell in the colormap to the input
  // color specified by char string in cColorName if the color does not 
  // already exist in the colormap.  The string must match
  // a color in the X servers color database or an error is returned.
  // If it does exist in an allocated part of the colormap,
  // then change the value of m_pulColorIndex[nIndex] to the index of the cell
  // in the colormap.

  XColor hXcolorExact, hXcolorScreen;
  Status nStat = XLookupColor(XtDisplay(m_hParent), m_hColormap, cColorName, 
			      &hXcolorExact, &hXcolorScreen);
  if (0 != nStat) 
    {
      return (nSetColor(nIndex, hXcolorExact.red,
			hXcolorExact.green, hXcolorExact.blue, nTolerance));
    }
  else
    {
      return (nStat);
    }
}

unsigned long
CXcolormap::ulGetColorIndex(const char *pcColorName, const int nTolerance)
{
  // Return the color index in the allocated cells of the colormap of the
  // named color cColorName or ms_ulCOLORERROR if the color is not in the 
  // colormap

  XColor hXcolorExact, hXcolorScreen;
  Status nStat;

  nStat = XLookupColor(XtDisplay(m_hParent), m_hColormap, pcColorName, 
		      &hXcolorExact, &hXcolorScreen);
  if (1 == nStat)
    {
      // Get all the colors in the current colormap 

      return (ulGetColorIndex(hXcolorExact.red, 
			      hXcolorExact.green, hXcolorExact.blue,
			      nTolerance));
    }
  return (ms_ulCOLORERROR); // Need inquire color or something
}

//#define abs(f)			( f >=0 ? f : -(f))

unsigned long
CXcolormap::ulGetColorIndex(const unsigned int unRed, 
			    const unsigned int unGreen,
			    const unsigned int unBlue, const int nToleranceIn)
{
  // Return the color index in the allocated cells of the colormap 
  // of a color close to specified by the RGB inputs within a given tolerance
  // or ms_ulCOLORERROR if the color is not in the colormap

  int    i;
  int    nDiff, nMinDiff, nIndex;
  int    nTolerance;

  nTolerance = nToleranceIn;

  // For some reason, the SUNOS doesn't return the rgb values that
  // were allocated, use DTDISPLAY_COLOR_TOLEXTRA to set m_nToleranceExtra

  nTolerance += m_nToleranceExtra;

  // Get all the allocated colors in the current colormap 
  // Do this each time, since the colors may have changed since the last
  // call

  for (i = 0; i < m_nNumAllocated; i++)
    {
      m_phXcolor[i].pixel = m_pulColorAlloc[i];
    }
  XQueryColors(XtDisplay(m_hParent), m_hColormap, m_phXcolor, m_nNumAllocated);

  nMinDiff =  nTolerance + 1;
  nIndex  =  0;                             // For safety
  for (i = 0; i < m_nNumAllocated; i++)
    {
      nDiff = 0;
      nDiff =  max(nDiff, ABS((int) m_phXcolor[i].red   - (int) unRed));
      nDiff =  max(nDiff, ABS((int) m_phXcolor[i].green - (int) unGreen));
      nDiff =  max(nDiff, ABS((int) m_phXcolor[i].blue  - (int) unBlue));
      if (nDiff <= nMinDiff) 
	{
	  // Save index of the closest one

	  nMinDiff = nDiff;
	  nIndex  = i;
	}
    }

  if (nMinDiff <= nTolerance)
    {
      // If the closest one meets the tolerance, then return the index

      return (m_phXcolor[nIndex].pixel);
    }
  else
    {
      // Otherwise return error
      
      return (ms_ulCOLORERROR); // Need inquire color or something
    }
}

void
CXcolormap::vSetWindow(Widget w)
{
  XSetWindowColormap(XtDisplay(w), XtWindow(w), m_hColormap);
  XtVaSetValues(w, XmNcolormap, m_hColormap, NULL);
}

void
CXcolormap::vInstall(void)
{
  // Force the window manager install the colormap
  // Even though the X Window and Motif manuals suggest that the window
  // manager is responsible for automatically installing the colormap, I
  // I have found this is not always the case, so we have a way to install
  // it ourselves.

  XInstallColormap(XtDisplay(m_hParent), m_hColormap);
}
//+Private functions

