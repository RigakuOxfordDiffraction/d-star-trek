#ifndef DT_CXHELP_H
#define DT_CXHELP_H
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
// CXhelp.h        Initial author: J.W. Pflugrath           21-Jul-1995
//    This file is the header file for class CXhelp.
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
#include <stdio.h>

#include "Dtrek.h"
#include "dtreksys.h"
#include "Cstring.h"
#include "UIComponent.h"

//+Definitions and constants

typedef struct _tagHelpResource 
{
  char *pcHelpLine;
  char *pcHelpHTML;
} tagHelpResource;

//+Code begin

class CXhelp {

public:


////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

   CXhelp (Widget wShellIn);
  ~CXhelp ();

////////////////////////////////////////////////////////////////////////
//+Member function prototypes
 public:

inline  void vDoNotPost(void) { m_nDoPost = 0; }
inline  void vDoPost(void)    { m_nDoPost = 1; }

virtual void vHelpCB(Widget, XtPointer, XtPointer);
virtual void vPostOneLineCB(Widget, XtPointer, XtPointer);
virtual void vUnpostOneLineCB(Widget, XtPointer, XtPointer);
virtual void vHTMLHelp(Widget, XtPointer, XtPointer);

static void vPostOneLineCBCallback(Widget, XtPointer, XtPointer);
static void vUnpostOneLineCBCallback(Widget, XtPointer, XtPointer);
static void vHTMLHelpCallback(Widget, XtPointer, XtPointer);

//

private:

  Widget  m_wParentShell;
  int     m_nDoPost;
  String  m_pcProgname;
  Cstring m_sSavedTitle;
  Cstring m_sHelpLine;
  Cstring m_sDelimiter;
  Cstring m_sPrefix;       // Prefix prepended to help filenames (e.g. "cal_")
  Cstring m_sSuffix;       // Suffix added to help filenames     (e.g. ".html")
  Cstring m_sDirectory;    // Directory for where files are      (e.g. "./")
  Cstring m_sHTMLViewer;   // Name of html viewer (e.g. "mosaic")
  Cstring m_sSeparator;    // Separator in incoming help lines   (e.g. "$$");
  Cstring m_sTmpDirectory; // Directory for temporary files (e.g. "/tmp")

  static XtResource m_hHelpResources[];
  tagHelpResource m_tHelpResource;

  static XtActionsRec m_hHelpActions[];

  static Cstring ms_sDelimiter;
  static Cstring ms_sPrefix;       // Prefix prepended to help filenames (e.g. "cal_")
  static Cstring ms_sSuffix;       // Suffix added to help filenames     (e.g. ".html")
  static Cstring ms_sDirectory;    // Directory for where files are      (e.g. "./")
  static Cstring ms_sHTMLViewer;   // Name of html viewer (e.g. "mosaic")
  static Cstring ms_sSeparator;    // Separator in incoming help lines   (e.g. "$$");
  static Cstring ms_sTmpDirectory; // Directory for temporary files (e.g. "/tmp")


virtual bool bHTMLViewerIsActive(int *pnPID);

};  // end of class CXhelp

/*
 * nonmember functions for hooking into Motif
 */
void vPostHelp(Widget, XEvent *, String *, Cardinal *);
void vUnPostHelp(Widget, XEvent *, String *, Cardinal *);

#endif   // DT_CXHELP_H

