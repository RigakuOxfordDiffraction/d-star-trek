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
// CXhelp.cc       Initial author: J.W. Pflugrath           28-Aug-1995
//    This file contains the member functions of class CXhelp.
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
//    This implements a help callback procedure for X-window widgets.  Two forms
//    of help are implemented.  The first form is for helpCallbacks.  If
//    the help callback ultimately calls CXhelp.vHTMLHelp(...), then 
//    help is displayed by using hypertext markup language (HTML)
//    files and HTML viewer such as Mosaic or NetScape.  The second form is
//    a one-line help message appended to a window title.  This is usually
//    activated by a armCallback or focusCallback and deactivated by a
//    disarmCallback or losingFocusCallback.
//
//    The following X resources may be used by this class:
//    ID                   Default        Description
//    
//    helpHTMLViewer       netscape       Program invoked to view html files
//    helpPrefix           ""             Prefix added to html file basenames
//    helpSuffix           .html          Suffix added to html file basenames
//    helpDelimiter        " -- Help: "   Text between window title and help
//    helpDirectory        ./             Directory for html files
//    helpTmpDirectory     /tmp/          Temporary directory (do not change)
//    helpSeparator        $$             Separator used by calling procedures
//                                        to separate a file basename or 
//                                        resource name from one-line of help 
//                                        text
//
//    Furthermore, when calling the vHTMLHelp member function, additional
//    resources may be used.  These are specified by the file basename passed
//    to the procedure with the string defined by the m_ssSuffix variable.  For
//    example, if the basename passed were "FileOpen", then the html file
//    would be helpDirectory + helpPrefix + "FileOpen" + helpSuffix or using
//    the defaults:
//                 ./FileOpen.html
//    The resource "htmlFileOpen" could be used to specify the file basename.
//    If the resource file contained the line
//     *htmlFileOpen:      openfile
//    Then after translation the html file would be:
//                 ./openfile.html
//    This has been tested with mosaic and netscape, but not with other html
//    viewers.
//    If the file passed to the HTML viewer does not exist, then the viewer
//    will report it, but this class does not do so yet.
//
//    If the basename for the html file is "", then the program name is used
//    for the basename.
//
//    For the one-line help mode, the resource name becomes 
//    *basename + *widget + ".helpLine".  For example:
//    *dtdisplay*tfRotStart.helpLine:  Rotation start angle in degrees
//
//    The resource prefixes "helpLine" and "html" are used to distinguish
//    resources between the two forms of help.
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "CXhelp.h"              // Class definition and prototypes
#include <string.h>

using namespace std;

//+Code begin

//+Definitions, constants, and initialization of static member variables

Cstring CXhelp::ms_sDelimiter    = " -- Help: ";
Cstring CXhelp::ms_sPrefix       = "";
Cstring CXhelp::ms_sSuffix       = ".html";
Cstring CXhelp::ms_sDirectory    = "./";
Cstring CXhelp::ms_sTmpDirectory = "/tmp/";
Cstring CXhelp::ms_sHTMLViewer   = "$(netscape)";
Cstring CXhelp::ms_sSeparator    = "$$";
XtResource CXhelp::m_hHelpResources[] =
{
  {
    "helpLine",                             // Name of resource
    "HelpLine",                             // Class of resource
    XmRString,                              // Required type
    sizeof(char *),                         // Size of expected type
    XtOffset(tagHelpResource*, pcHelpLine), // Destination address
    XmRString,                              // Type of default value
    (XtPointer) ""                          // Default value
  },
  {
    "helpHTML",                             // Name of resource
    "HelpHTML",                             // Class of resource
    XmRString,                              // Required type
    sizeof(char *),                         // Size of expected type
    XtOffset(tagHelpResource*, pcHelpHTML), // Destination address
    XmRString,                              // Type of default value
    (XtPointer) "No HTML help available"    // Default value

  }
};
XtActionsRec CXhelp::m_hHelpActions[] =
{
  {
    "vPostHelp",
    vPostHelp
  },
  {
    "vUnPostHelp",
    vUnPostHelp
  }
};

/*
 * Pointer to class, init'ed in creation operator.  This is used for hooks 
 * into MOTIF.  The hooks call the (public) member functions vPostOneLineCB()
 * and vUnpostOneLineCB() of the instance that this variable points to.  
 * The first instance of this class created is assigned to this variable.
 * If the instance is destroyed, then the help posts may not work correctly
 * until another instance of this class is created.
 */

CXhelp *g_poCXhelp = NULL;

//+Public functions

// Constructors, destructors and assignments

CXhelp::CXhelp(Widget wShellIn)
{
  // Construct a CXhelp object

  // Initialize member variables with defaults

  m_nDoPost       = 2;
  m_sHelpLine     = "NoHeLpIsGoOdHeLp";

  m_sDelimiter    = ms_sDelimiter;
  m_sPrefix       = ms_sPrefix;
  m_sSuffix       = ms_sSuffix;
  m_sDirectory    = ms_sDirectory;
  m_sHTMLViewer   = ms_sHTMLViewer;
  m_sSeparator    = ms_sSeparator;
  m_sTmpDirectory = ms_sTmpDirectory;

  // Get the parent shell so we know which window title gets the help lines

  m_wParentShell = wShellIn;
  while (m_wParentShell && !XtIsWMShell(m_wParentShell))
    m_wParentShell = XtParent(m_wParentShell);

  // Get a pointer to the current title 

  String pcTitle;
  XtVaGetValues(m_wParentShell, XmNtitle, &pcTitle, NULL);
  m_sSavedTitle = (const char *)pcTitle;
  
  // Get the name of the program so we can use XGetDefaults to find resources

  String pcClass;
  XtGetApplicationNameAndClass(XtDisplay(m_wParentShell), &m_pcProgname, &pcClass);

  // Examine resource database for any user preferences

  String  pcText;
  pcText = XGetDefault (XtDisplay(m_wParentShell), m_pcProgname, "helpHTMLViewer");
  if (pcText != NULL)
    {
      m_sHTMLViewer = pcText;
    }
  pcText = XGetDefault (XtDisplay(m_wParentShell), m_pcProgname, "helpPrefix");
  if (pcText != NULL)
    {
      m_sPrefix = pcText;
    }
  pcText = XGetDefault (XtDisplay(m_wParentShell), m_pcProgname, "helpSuffix");
  if (pcText != NULL)
    {
      m_sSuffix = pcText;
    }
  pcText = XGetDefault (XtDisplay(m_wParentShell), m_pcProgname, "helpDelimiter");
  if (pcText != NULL)
    {
      m_sDelimiter = pcText;
    }
  pcText = XGetDefault (XtDisplay(m_wParentShell), m_pcProgname, "helpDirectory");
  if (pcText != NULL)
    {
      m_sDirectory = pcText;
    }
  pcText = XGetDefault (XtDisplay(m_wParentShell), m_pcProgname, "helpTmpDirectory");
  if (pcText != NULL)
    {
      m_sTmpDirectory = pcText;
    }
  pcText = XGetDefault (XtDisplay(m_wParentShell), m_pcProgname, "helpSeparator");
  if (pcText != NULL)
    {
      m_sSeparator = pcText;
    }

  // Translate potential symbols

  m_sTmpDirectory = sTransSymbol(m_sTmpDirectory);
  m_sDirectory    = sTransSymbol(m_sDirectory);
  m_sHTMLViewer   = sTransSymbol(m_sHTMLViewer);

// Add MOTIF actions (via resource file for now)

  //if(NULL == g_poCXhelp){
     //g_poCXhelp = this;
  //}
  //XtAppAddActions(XtWidgetToApplicationContext(m_wParentShell),
                  //m_hHelpActions,XtNumber(m_hHelpActions));

}

CXhelp::~CXhelp() 
{
  // if(g_poCXhelp == this)
        //g_poCXhelp = NULL;
}

void
CXhelp::vHelpCB(Widget w, XtPointer clientData, XtPointer callData)
{
  char *pcHelp = (char *)callData;
  cout << "Help text is: " << pcHelp << endl;
//  XtVaSetValues(m_wParentShell, XmNtitle, pcHelp, NULL);
}

void
CXhelp::vPostOneLineCBCallback(Widget w, 
			       XtPointer clientData, XtPointer callData)
{
    UICallbackStruct *data = (UICallbackStruct *) clientData;
    CXhelp *obj = (CXhelp *)data->object;
    obj->vPostOneLineCB(w, data->client_data, callData);
}

void
CXhelp::vPostOneLineCB(Widget w, XtPointer clientData, XtPointer callData)
{
  if (0 < m_nDoPost)
    {
      Cstring sTemp;

      // Get a pointer to the current title of the parent shell
      // so we know which window title gets the help lines

      m_wParentShell = w;
      while (m_wParentShell && !XtIsWMShell(m_wParentShell))
	m_wParentShell = XtParent(m_wParentShell);

      String pcTitle;
      XtVaGetValues(m_wParentShell, XmNtitle, &pcTitle, NULL);
      m_sSavedTitle = (const char *)pcTitle;

      // Remove any previous delimiter and help line from title

      sTemp = m_sSavedTitle.before(m_sDelimiter);
      if ("" != sTemp)
	m_sSavedTitle = sTemp;

      sTemp       = "";
      if (2 == m_nDoPost)
	{
	  // Try to get helpline text from resource database first. 

	  if (NULL != clientData)
	    {
	      m_hHelpResources[0].default_addr = clientData;
	    }
	  else
	    {
	      m_hHelpResources[0].default_addr = (XtPointer)"No help available";
	    }
	  XtGetSubresources(XtParent(w),
			    (XtPointer) &m_tHelpResource,
			    XtName(w),
			    "",             // We do know the classname
			    m_hHelpResources,
			    1,
			    NULL,
			    0);

	  sTemp = m_tHelpResource.pcHelpLine;
	}

      // Now add the current help line to the title and set the title in widget

      if ("" != sTemp)            // This should ALWAYS be true now
	{
	  // A help line was found either in the database or in clientData
	  // or the default "No help available" is used

	  m_sHelpLine = sTemp;
	  XtVaSetValues(m_wParentShell, XmNtitle, 
			(m_sSavedTitle + m_sDelimiter + m_sHelpLine).string(), 
			NULL);
	}
    }
}

void
CXhelp::vUnpostOneLineCBCallback(Widget w, 
				 XtPointer clientData, XtPointer callData)
{
    UICallbackStruct *data = (UICallbackStruct *) clientData;
    CXhelp *obj = (CXhelp *)data->object;
    obj->vUnpostOneLineCB(w, data->client_data, callData);
}

void
CXhelp::vUnpostOneLineCB(Widget w, XtPointer clientData, XtPointer callData)
{
  // Unpost the help line ONLY if the title is the same as set before
  // That is the title was changed only by the vPostOneLineCB method.

  String pcTitle;
  Cstring sTemp;

  XtVaGetValues(m_wParentShell, XmNtitle, &pcTitle, NULL);
  sTemp = (const char *)pcTitle;
  if (sTemp.index(m_sSavedTitle) == 0)
      {
	XtVaSetValues(m_wParentShell, XmNtitle, m_sSavedTitle.string(), NULL);
      }
}

void
CXhelp::vHTMLHelpCallback(Widget w, 
				 XtPointer clientData, XtPointer callData)
{
    UICallbackStruct *data = (UICallbackStruct *) clientData;
    CXhelp *obj = (CXhelp *)data->object;
    obj->vHTMLHelp(w, data->client_data, callData);
}

void
CXhelp::vHTMLHelp(Widget w, XtPointer clientData, XtPointer callData)
{
  // Display with the HTML viewer the file coded in the clientData argument.
  // (If the file does not exist, then post an error message.)

  char   *pcHelp = (char *)clientData;
  Cstring sTemp;
  Cstring sFilename;
  int     nPID;
  FILE    *fp;

  sFilename = pcHelp;
  sTemp     = sFilename.before(m_sSeparator);
  if ("" == sTemp)
    {
      // No separator found, so use cProgname as the filename
      sTemp = m_pcProgname;
    }

  // See if this needs a translation via the resource database

  pcHelp = XGetDefault(XtDisplay(w), m_pcProgname,
		       ("html"+sTemp).string());
  if (NULL == pcHelp)
    {
      // Try the other way if the first is null.

      pcHelp = XGetDefault(XtDisplay(w), m_pcProgname,
			   (sTemp+"html").string());
    }
  if (NULL != pcHelp)
    {
      sTemp = pcHelp;
    }

  // See if a version of the HTML viewer is running in this process tree...

  sTemp = sTransSymbol(m_sDirectory + m_sPrefix + sTemp + m_sSuffix);

  if (!bHTMLViewerIsActive(&nPID))
    {
      // Viewer is not active, so start one up
      
      sTemp = m_sHTMLViewer + " " + sTemp + " &";
      nDoSystemCommand(sTemp);
    }
  else if (0 <= m_sHTMLViewer.find("mosaic"))
    {
      // Mosaic viewer is active, so work with it.
      // There is a feature of NCSA mosaic: If you create a file
      // that contains a html filename in the /tmp directory that is called
      // "Mosaic.pid" where pid is the process id of a running mosaic process
      // and then send a "kill -USR1" signal to that process, then it will
      // automatically load and view the file specified by the filename.

      char cCommand[255];
      sprintf(cCommand, "%sMosaic.%d", m_sTmpDirectory.string(), nPID);
      fp = fopen(cCommand, "w");
      if (!fp) return;
      fprintf( fp,"goto\n%s\n", sTemp.string());
      fclose( fp );
      sprintf(cCommand,"kill -USR1 %d", nPID);
      nDoSystemCommand((Cstring)cCommand);
    }
  else if (0 <= m_sHTMLViewer.find("netscape"))
    {
      // Netscape viewer is active.
      // Use the -remote feature of netscape.

      sTemp = m_sHTMLViewer + " -remote 'openFile(" + sTemp + ")'";
      nDoSystemCommand(sTemp);
    }
  else if (0 <= m_sHTMLViewer.find("mozilla"))
    {
      // Mozilla viewer is active.
      // Use the -remote feature of mozilla

      sTemp = m_sHTMLViewer + " -remote 'openFile(" + sTemp + ")'";
      //      printf("***INFO: help command: %s\n", sTemp.string());
      nDoSystemCommand(sTemp);
    }
}

bool
CXhelp::bHTMLViewerIsActive(int *pnPID)
{
  // See if an HTML viewer is running in this process tree and if so
  // return its process ID in *pnPID.
  // Returns True if it is running
  // Returns False if it is not running

  Cstring sTemp;
  FILE    *fp;
  char    cLine[85];

  sTemp = "ps >";
  sTemp += m_sTmpDirectory + "HELP.tmp";

  nDoSystemCommand(sTemp);

  sTemp = m_sTmpDirectory + "HELP.tmp";

  // The stuff below should be replaced with c++ i/o routines

  fp = fopen(sTemp.string(), "r");
  
  if (!fp) 
    return (False);

  while (!feof(fp) )
    {
      fgets(cLine, 80, fp);
      if (strstr(cLine, sFileGetBasename(m_sHTMLViewer).string()))
	{
	  fclose (fp);
	  (void) nFileDelete(sTemp);
	  sscanf (cLine, "%d", pnPID);
	  return (True);
	}
    }
  fclose (fp);
  (void) nFileDelete(sTemp);
  return (False);
}

/****************************************************************************
 * Non member functions which are used for MOTIF hooks.                     *
 ****************************************************************************/
void
vPostHelp(Widget w,
          XEvent *event,
          String *params,
          Cardinal *nparams)
{
  g_poCXhelp->vPostOneLineCB(w, NULL, NULL);
  return;
}
void
vUnPostHelp(Widget w,
            XEvent *event,
            String *params,
            Cardinal *nparams)
{
  g_poCXhelp->vUnpostOneLineCB(w, NULL, NULL);
  return;
}

