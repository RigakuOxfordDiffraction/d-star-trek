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
// CXprop.cc            Initial author: J.W. Pflugrath           18-Mar-1996
//    This file contains the member functions of class CXprop.
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

#include "CXprop.h"              // Class definition and prototypes
#include "dtreksys.h"            // For sGetEnv()
#include "X11/Xproto.h"          // For defs of X protocol request numbers

using std::cout;
using std::endl;
using std::flush;
using std::cerr;

//+Code begin


CXprop *CXprop::ms_pThis         = NULL;

CXprop::CXprop(const Cstring& sWindowName, const char *pcDisplayName)
{
  // This routine may be the FIRST X Window routine accessed, so the
  // DISPLAY environment variable must be legitimate, so call vModifyDISPLAY().

  vModifyDISPLAY();

  m_bAddedWindowID = False;
  m_hWindow        = 0;
  m_sWindowName    = sGetEnv("DTDISPLAY_PREFIX") + sWindowName;
  m_hDisplay       = XOpenDisplay(pcDisplayName);
  m_bUsingLocalErrorHandler = False;
  m_prnOrigErrorHandler = NULL;

  if (NULL == m_hDisplay)
    {
      // Some kind of error.

      return;
    }
  vIgnoreBadWindow(True);  // Trap BadWindow errors on X_ChangeProperty requests
}

CXprop::~CXprop()
{
  if (NULL != m_hDisplay)
    {
      if (m_bAddedWindowID)
	{
	  // If this instance object actually added the window ID property 
	  // to the root window, then
	  // Delete all properties placed on this window that begin with
	  // the same prefix.

	  // WARNING: be careful with this since we added DTDISPLAY_PREFIX

	  Cstring sPrefix = m_sWindowName.before("_");
	  m_hWindow       = hGetWindowID();
	  if (0 != m_hWindow)
	    {
	      char *pcTemp;
	      Atom *pAtoms;
	      int   nNumProps;
	      int   i;
	      pAtoms = XListProperties(m_hDisplay, m_hWindow, &nNumProps);
	      
	      if (NULL != pAtoms)
		{
		  for (i = 0; i < nNumProps; i++)
		    {
		      pcTemp = XGetAtomName(m_hDisplay, pAtoms[i]);
		      if (   (NULL != pcTemp) 
			  && (0 == Cstring(pcTemp).find(sPrefix)))
			XDeleteProperty(m_hDisplay, m_hWindow, pAtoms[i]);
		      else
			cout << "Property name: " << *pcTemp << endl;
		      XFree(pcTemp);
		    }
		  XFree((char *)pAtoms);
		}
	    }
	  
	  (void) nSetWindowID(0);  // Set window id to 0

	  // If this object added a window ID to the root window, then delete
	  // that property, too.
	  
	  Atom hAtom = XInternAtom(m_hDisplay, m_sWindowName.string(), True);
	  if (None != hAtom)
	    XDeleteProperty (m_hDisplay, DefaultRootWindow(m_hDisplay), hAtom);
	}
      XCloseDisplay(m_hDisplay);  // Be careful!
    }
}

Atom
CXprop::hSetProperty(const Cstring& sPropertyName, const Cstring& sPropText,
		     bool bIfOnly)
{
  // Set the text property with name sPropertyName to the value sPropText
  // for the window specified by m_sWindowName.
  // Return the Atom whose property was set, or None on failure.

  Atom   hAtom;

  // Get window id specified by m_sWindowName

  m_hWindow = hGetWindowID();

  if (0 == m_hWindow)
    {
      return (None);  // Error occurred
    }

  // Now set the named property on the window,
  // If bIfOnly = True,  do not create it if it does not exist (default)
  //            = False,        create it if it does not exist

  hAtom     = XInternAtom(m_hDisplay, sPropertyName.string(), bIfOnly);
  if (None == hAtom)
    {
      return (None);  // No atom with that name on the display
    }
  XTextProperty hTextProperty;
  char *pcTemp = sPropText.string();
  XStringListToTextProperty(&pcTemp, 1, &hTextProperty);
  XSetTextProperty(m_hDisplay, m_hWindow, &hTextProperty,hAtom);
  XFlush(m_hDisplay);
  if (NULL != hTextProperty.value)
    {
      XFree((char *)hTextProperty.value);
      hTextProperty.value = NULL;
    }
  return (hAtom);
}

int
CXprop::nGetProperty(const Cstring& sPropertyName, Cstring *psPropText)
{
  // Get the text property with name sPropertyName to the value *psPropText.
  // Return 0 on success, non-0 on failure.

  Atom   hAtom;
  Status nStat;

  // Get window id specified by m_sWindowName

  m_hWindow = hGetWindowID();

  if (0 == m_hWindow)
    {
      *psPropText = "Error finding Window ID.";
      return (-1);
    }

  // Now set the named property on the window, 
  // (create it if it does not exist)

  hAtom     = XInternAtom(m_hDisplay, sPropertyName.string(), False);

  if (None == hAtom)
    {
      *psPropText = "ERROR finding property Atom.";
      return (-2);  // No atom with that name on the display
    }

  XTextProperty hTextProperty;
  hTextProperty.value = NULL;
  nStat = XGetTextProperty(m_hDisplay, m_hWindow, &hTextProperty, hAtom);

  if ( (0 != nStat) && (NULL != hTextProperty.value) )
    {
//      cout << "Text property is: " << hTextProperty.value << endl;
      *psPropText = (char *)hTextProperty.value;
//      cout << "Text string is: " << *psPropText << endl;
      XFree((char *)hTextProperty.value);
      hTextProperty.value = NULL;   
      return (0);
    }
  else
    {
      *psPropText = "ERROR finding property text.";
      return (-3);
    }
}

int
CXprop::nGetProperty(const Atom& hAtomIn, Cstring *psPropText)
{
  // Get the text property with atom hAtomIn to the value *psPropText.
  // Return 0 on success, non-0 on failure.

  Status nStat;

  // Get window id specified by m_sWindowName

  m_hWindow = hGetWindowID();

  if (0 == m_hWindow)
    {
      *psPropText = "Error finding Window ID.";
      return (-1);
    }

  // Now set the named property on the window, 
  // (create it if it does not exist)

  XTextProperty hTextProperty;
  hTextProperty.value = NULL;
  nStat = XGetTextProperty(m_hDisplay, m_hWindow, &hTextProperty, hAtomIn);

  if ( (0 != nStat) && (NULL != hTextProperty.value) )
    {
//      cout << "Text property is: " << hTextProperty.value << endl;
      *psPropText = (char *)hTextProperty.value;
//      cout << "Text string is: " << *psPropText << endl;
      XFree((char *)hTextProperty.value);
      hTextProperty.value = NULL;   
      return (0);
    }
  else
    {
      *psPropText = "ERROR finding property text.";
      return (-3);
    }
}

Window
CXprop::hGetWindowID(void)
{

  // Get window id specified by m_sWindowName, but do not create an Atom if
  // it does not exist.

  Atom hAtom;

  if (NULL == m_hDisplay)
    return ( (Window) 0);  // No display available

  hAtom = XInternAtom(m_hDisplay, m_sWindowName.string(), True);
  if (None == hAtom)
    {
      return ((Window)0);  // No atom with that name on the display
    }

  Window       *phWindow;
  Window        hWindow;
  Atom          type;
  int           format;
  unsigned long nitems;
  unsigned long left;
		      
  // Get the current value of the property associated with the atom

  if (   (Success == XGetWindowProperty(m_hDisplay,
					DefaultRootWindow(m_hDisplay),
					hAtom, 0, 4, False, XA_WINDOW,
					&type, &format, &nitems, &left,
					(unsigned char **)&phWindow) )
      && (XA_WINDOW == type) )
    {
      hWindow = *phWindow;
      XFree((char *)phWindow);

      // Does hWindow still really exists?
//      if (XtIsWindow
      return (hWindow);
    }
  else
    {
      return (0);
    }
}

int
CXprop::nSetWindowID(Window hWindow)
{
  // Set the window ID in a property on the root window of the display.
  // Create the property if it does not exist.

  Atom hAtom;
  hAtom   = XInternAtom(m_hDisplay, m_sWindowName.string(), False);

  if (None == hAtom)
    {
      return (-1);  // Problem with getting the Atom
    }

  // Store the window ID on the root window as property named by m_sWindowName

  XChangeProperty (m_hDisplay, DefaultRootWindow(m_hDisplay),
		   hAtom, XA_WINDOW,
		   32, PropModeReplace, (unsigned char *)&hWindow, 1);
  XFlush(m_hDisplay);
  m_bAddedWindowID = True;
  return (0);
}

void
CXprop::vIgnoreBadWindow(bool bFlag)
{
  if (bFlag)
    {
      if (!m_bUsingLocalErrorHandler)
	{
	  m_bUsingLocalErrorHandler = True;
	  m_prnOrigErrorHandler = XSetErrorHandler(nCXpropErrorHandler);
	  ms_pThis = this;
	}
    }
  else
    {
      if (NULL != m_prnOrigErrorHandler)
	{
	  m_bUsingLocalErrorHandler = False;
	  XSetErrorHandler(m_prnOrigErrorHandler);
	}
    }
}

int
CXprop::nCXpropErrorHandler(Display *dpy, XErrorEvent *pEvent)
{
  // This is a static function, so we cannot use non-static member variables
  // in this routine.  Hence the use of the static ms_pThis member variable.
  // This kludge will only work flawlessly if there is a single instance
  // of the CXprop class within the program.  It should still work in almost
  // all other cases.

  if (   (BadWindow        == pEvent->error_code)
      && (X_ChangeProperty == pEvent->request_code) )
    {
//      cerr << "INFO: Bad Window error for X_ChangeProperty ignored!\n";
      return (0);
    }
  if (   (BadValue      == pEvent->error_code)
      && (X_QueryColors == pEvent->request_code) )
    {
      cerr << "INFO: BadValue error for X_QueryColors ignored!\n";
      return (0);
    }
  else if ( (NULL != ms_pThis) && (NULL != ms_pThis->m_prnOrigErrorHandler) )
    {
      return(ms_pThis->m_prnOrigErrorHandler(dpy, pEvent));
    }
  return (-1);
}
