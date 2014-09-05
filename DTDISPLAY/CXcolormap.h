#ifndef DT_CXCOLORMAP_H
#define DT_CXCOLORMAP_H
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
// CXcolormap.h        Initial author: J.W. Pflugrath           21-Jul-1995
//    This file is the header file for class CXcolormap.
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

#include "Xm/Xm.h"

//#include <iostream.h>
//#include <iomanip.h>
#include "Dtrek.h"
#include "Cstring.h"

//+Definitions and constants

//+Code begin

class CXGcontext;
class CXcolormap : public UIComponent {

public:

  Widget         m_hParent;     // Handle back to X11
  Colormap       m_hColormap;   // Handle to X11 colormap
  XColor        *m_phXcolor;
  int            m_nAvail;      // The maximum number of available colors
                                // in the colormap.  It could be 2.
  int            m_nNumRamp;    // Num colors in the color ramp
  int            m_nAllocated; // The number of colorcells to allocate in cmap
  int            m_nNumAllocated; // Number colorcells allocated in cmap
  int            m_nExtra;     // The number of extra colors to try to allocate.
  bool           m_bNotDefault;
  unsigned long *m_pulColorIndex;
  unsigned long *m_pulColorAlloc; // Pointer to allocated cells in a colormap
  int           *m_pnColorUsed;   // Pointer to used (i.e. set by this object)
                                  //        cells in a colormap
  int            m_nTolerance;
  int            m_nToleranceExtra;

  static const unsigned long   ms_ulCOLORERROR;
  int            m_nDefaultDepth;
  Visual        *m_phDefaultVisual;
  bool           m_bReadOnlyMode; 

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

   CXcolormap (const char *, Widget, const unsigned int unRequested, 
	       const int nTolerance=0);
  ~CXcolormap ();

////////////////////////////////////////////////////////////////////////
//+Member fFunction prototypes
 public:

//

inline int nNumColors(void) { return (m_nNumRamp); }
inline int nGetTolerance(void) { return (m_nTolerance); }
inline void vSetTolerance(const int nTolIn) { m_nTolerance = nTolIn; }
inline unsigned char ucColorIndex(int nValue) 
  { return ((unsigned char) m_pulColorIndex[nValue]); }
inline unsigned long ulColorIndex(int nValue) 
  { return (m_pulColorIndex[nValue]); }
inline unsigned short uiColorIndex(int nValue) 
  { return ( (unsigned short)m_pulColorIndex[nValue]); }
int nLoadGrey(const int nNumGreyIn);
int nSetColor(const int nIndex, const unsigned int nRed,
	      const unsigned int nGreen, const unsigned int nBlue, 
	      const int nTolerance = 0);

void vSetWindow(Widget w);
void vInstall(void);

int nSetColor(const int nIndex, const char *cColorName, 
	      const int nTolerance = 0);
unsigned long ulGetColorIndex(const unsigned int unRed, 
			      const unsigned int unGreen, 
			      const unsigned int unBlue,
			      const int nTolerance=0);
unsigned long ulGetColorIndex(const char *cColorName, const int nTolerance = 0);

private:


};  // end of class CXcolormap

#endif   // DT_CXCOLORMAP_H
