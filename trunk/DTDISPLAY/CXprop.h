#ifndef DT_CXPROP_H
#define DT_CXPROP_H
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
// CXprop.h        Initial author: J.W. Pflugrath           18-Mar-1996
//    This file is the header file for class CXprop.
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

//#include <iostream.h>
#include "X11/Intrinsic.h"
#include "X11/Xatom.h"
// For typedef of bool in some cases:
#include "Dtrek.h"
#include "Cstring.h"

//+Definitions and constants

//+Code begin
class CXprop;  // forward declaration

class CXprop {

protected:

  Display       *m_hDisplay;
  Window         m_hWindow;
  Cstring        m_sWindowName;
  Bool           m_bAddedWindowID;
  Bool           m_bUsingLocalErrorHandler;
  XErrorHandler  m_prnOrigErrorHandler;
  static CXprop *ms_pThis;
public:

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

  CXprop (const Cstring& sWindowName, const char *pcDisplayName = NULL);
  ~CXprop ();

////////////////////////////////////////////////////////////////////////
//+Member function prototypes

  Atom  hSetProperty(const Cstring& sPropertyName, const Cstring& sPropText,
		     bool bIfOnly = True);
  int   nGetProperty(const Cstring& sPropertyName, Cstring *psPropText);
  int   nGetProperty(const Atom& hAtomIn, Cstring *psPropText);
  int   nSetWindowID(Window hWindow);
  void  vIgnoreBadWindow(bool bFlag);
  static int nCXpropErrorHandler(Display *, XErrorEvent *);

private:

  Window hGetWindowID(void);

};  // end of class CXprop

#endif   // DT_CXPROP_H

