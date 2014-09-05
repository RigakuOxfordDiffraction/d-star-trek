#ifndef DT_CXRUBBERBAND_H
#define DT_CXRUBBERBAND_H
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
// CXrubberband.h        Initial author: J.W. Pflugrath           21-Jul-1995
//    This file is the header file for class CXrubberband.
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

#include "Xm/Xm.h"

//#include <iostream.h>
//#include <iomanip.h>
#include "Dtrek.h"

//+Definitions and constants

#define RB_NONE_MODE 0
#define RB_LINE_MODE 1
#define RB_CIRC_MODE 2
#define RB_POIN_MODE 3
#define RB_CIR2_MODE 4
#define RB_RECT_MODE 5
#define RB_TLIN_MODE 6

//+Code begin

class CXrubberband {

public:

  Widget         m_hParent;         // Handle back to X11
  Pixmap        *m_phParentPixmap;  // Pointer to a pixmap

  int            m_nMode;
  int            m_nOrigX, m_nOrigY;
  int            m_nEndX, m_nEndY;

  int            m_nWidth, m_nHeight;
  int            m_nDrawX, m_nDrawY;
  int            m_nDrawX2, m_nDrawY2;
  int            m_nLineWidth;    // This will be line thickness, too
  GC             m_hGC;
  XGCValues      m_hGCValues;

  Dimension      m_nParentWidth;
  Dimension      m_nParentHeight;
  XColor         m_hColor;        // Something to hold the color
  int            m_nColorAllocated;

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

   CXrubberband (Widget, Pixmap *, const int nStartX, const int nStartY, 
		 const int nWidthMax = 0, const int nHeightMax = 0,
		 const int nMode = 0);
  ~CXrubberband ();

////////////////////////////////////////////////////////////////////////
//+Member fFunction prototypes
 public:

//

void vDraw(const int nXIn, const int nYIn, const bool bRestore=True);
void vPan(const int nXIn, const int nYIn);
void vGetCoords (int *pnStartX, int *pnStartY, int *pnEndX, int *pnEndY);
int nSetLineColor(const char *cColorName);
int nSetLineColor(const unsigned int nRed, 
		  const unsigned int nGreen, const unsigned int nBlue);

inline int nGetMode(void) { return (m_nMode); }
inline void vSetLineWidth(const int nLineWidth ) { m_nLineWidth = nLineWidth; }

private:

void vSaveOutline(void);
void vRestoreOutline(void);

};  // end of class CXrubberband

#endif   // DT_CXRUBBERBAND_H

