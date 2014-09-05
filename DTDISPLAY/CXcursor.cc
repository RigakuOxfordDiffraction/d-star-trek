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
// CXcursor.cc       Initial author: J.W. Pflugrath           28-Aug-1995
//    This file contains the member functions of class CXcursor.
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
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "CXcursor.h"              // Class definition and prototypes

//+Code begin

//+Definitions, constants, and initialization of static member variables

//+Public functions

// Constructors, destructors and assignments

CXcursor::CXcursor(Widget wShellIn)
{
  // Construct a CXcursor object

  // Initialize member variables with defaults

  m_wTop = wShellIn;
  while (m_wTop && !XtIsWMShell(m_wTop))
    m_wTop = XtParent(m_wTop);
  m_hCursor = 0;
  for (int i = 0; i < 10; i++) m_unLast[i] = 0;
  m_nNumStacked = 0;
}

CXcursor::~CXcursor() 
{
  if (0 != m_hCursor)
    {
      vRestore();
      XFreeCursor(XtDisplay(m_wTop), m_hCursor);
    }
}

void
CXcursor::vSetFont(const unsigned int unNum, const bool bPush)
{
  m_hCursor = XCreateFontCursor(XtDisplay(m_wTop), unNum);
  XDefineCursor(XtDisplay(m_wTop), XtWindow(m_wTop), m_hCursor);
  XFlush(XtDisplay(m_wTop));
  if (bPush)
    {
      // Push this cursor "font" onto the stack

      m_nNumStacked++;
      if (9 < m_nNumStacked) m_nNumStacked = 9;   // Prevent stack overflow
      m_unLast[m_nNumStacked] = unNum;
    }
}

void
CXcursor::vSetWait(void)
{
  vSetFont(XC_watch, TRUE);
}

void
CXcursor::vReset(void)
{
  while (0 < m_nNumStacked)
    {
      vRestore();
    }
  vRestore();
}

void
CXcursor::vRestore(void)
{
  m_nNumStacked--;
  if (0 >= m_nNumStacked)
    {
      // Prevent stack underflow

      XUndefineCursor(XtDisplay(m_wTop), XtWindow(m_wTop));
      m_nNumStacked = 0;
    }
  else
    {
      vSetFont(m_unLast[m_nNumStacked], FALSE);
    }
}

void
CXcursor::vSetHand(void)
{
  vSetFont(XC_fleur, TRUE);
}

