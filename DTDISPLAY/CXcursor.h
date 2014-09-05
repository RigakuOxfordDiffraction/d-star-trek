#ifndef DT_CXCURSOR_H
#define DT_CXCURSOR_H
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
// CXcursor.h        Initial author: J.W. Pflugrath           21-Jul-1995
//    This file is the header file for class CXcursor.
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
#include <X11/cursorfont.h>
#include "UIComponent.h"
// For typdef of bool in some cases:
#include "Dtrek.h"
//+Definitions and constants

//+Code begin

class CXcursor {

public:

  Widget       m_wTop;             // Top widget;

  // Probably want the following to all be static since there is only 1 cursor
  // associated with the top widget.

  Cursor       m_hCursor;          // Handle to a cursor
  unsigned int m_unLast[10];       // Stack of 10 last used cursors numbers
  int          m_nNumStacked;      // Number on the stack

//+Constructors, destructors and assignments

   CXcursor (Widget wShellIn);
  ~CXcursor ();

////////////////////////////////////////////////////////////////////////
//+Member fFunction prototypes
 public:

  void vSetWait(void);
  void vSetHand(void);
  void vRestore(void);
  void vReset(void);
  void vSetFont(const unsigned int unNum, const bool bPush=FALSE);
//

private:

};  // end of class CXcursor

#endif   // DT_CXCURSOR_H

