//
// README: Portions of this file are merged at file generation
// time. Edits can be made *only* in between specified code blocks, look
// for keywords <Begin user code> and <End user code>.
//
//////////////////////////////////////////////////////////////
//
// Header file for CUIscrolledlist
//
//    Created by Builder Xcessory Version 5.0.
//    Generated by Code Generator Xcessory 5.0 (05/22/98) .
//
//    This class is a user interface "component", as described
//    in "Object-Oriented Programming with C++ and OSF/Motif",
//    by Douglas Young, Prentice Hall, 1992. ISBN 0-13-630252-1
//
//////////////////////////////////////////////////////////////


// Begin user code block <file_comments>
//
// Copyright (c) 1997 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// CUIscrolledlist.h           Initial author: J.W. Pflugrath       9-Apr-1997
//    Class definition for a user interface Motif scrolledlist.
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

// End user code block <file_comments>

#ifndef CUIscrolledlist_H
#define CUIscrolledlist_H
#include "UIComponent.h"

//
// Globally included information (change thru Output File Names Dialog).
//


//
// Class Specific Includes (change thru the class in Resource Editor).
//



// Begin user code block <head>
#include "Dtrek.h"
#include "Cstring.h"

typedef void (*ScrolledListCallback) (XtPointer, Widget, XtPointer);

// End user code block <head>

class CUIscrolledlist : public UIComponent

// Begin user code block <superclass>
// End user code block <superclass>
{

// Begin user code block <friends>
// End user code block <friends>

  public:

    CUIscrolledlist(const char *, Widget);
    CUIscrolledlist(const char *);
    virtual ~CUIscrolledlist();
    virtual void create(Widget);
    const char *const  className();
    void set_cuiscrolledlist_width(XtPointer);
    void set_cuiscrolledlist_height(XtPointer);
    void set_lilist_visibleItemCount(XtPointer);
    void set_lilist_selectionPolicy(XtPointer);
    void set_lilist_listSizePolicy(XtPointer);
    void set_lbtitle_width(XtPointer);
    void set_lbtitle_height(XtPointer);
    void set_lbtitle_leftOffset(XtPointer);
    
    // Begin user code block <public>


      XtPointer             m_pObj;        // For callback, pointer to an object
      ScrolledListCallback  m_prvScrolledListCallback;

      void vAddItem(const Cstring &rsItem, const int nPosition=0);
      void vAddOtherList(const Widget hOtherList);
      void vAddOtherList(const Widget hOtherList, const Cstring& rsTemplateIn);
      void vDeleteItem(const Cstring &rsItem);
      void vDeleteItem(const int nPosition);
      void vDeleteAll(void);
      void vDeleteNonImages(Cstring& rsScanTemplateIn);
      void vGetSelection(int *pnNumItems, Cstring **ppsItems);
      void vSetSelection(int nPosition,     Boolean bNotify=False, Boolean bMultiple=False);
      void vSetSelection(XmString sxString, Boolean bNotify=False, Boolean bMultiple=False);
      void vSetSelection(Cstring &rsString, Boolean bNotify=False, Boolean bMultiple=False);
      void vSetSelectionAll(Boolean bNotify=False);
      void vShowButtons(const Boolean bYesNo);
      void vSetSelectPolicy(const unsigned char ucPolicy);
      virtual Widget wGetListWidget(void);

    // End user code block <public>
  protected:
    // Classes created by this class
    
    // Widgets created by this class
    Widget _CUIscrolledlist;
    Widget _pbDeselectAll;
    Widget _pbSelectAll;
    Widget _swList;
    Widget _liList;
    Widget _lbTitle;
    Widget _frFrame;
    
    // These virtual functions are called from the private callbacks 
    // or event handlers intended to be overridden in derived classes
    // to define actions
    
    virtual void vSelectCB(Widget, XtPointer, XtPointer);
    
    // Begin user code block <protected>
    // End user code block <protected>
  private: 
    
    //
    // Default application and class resources.
    //
    static String         _defaultCUIscrolledlistResources[];
    static UIAppDefault   _appDefaults[];
    static Boolean        _initAppDefaults;
    //
    // Callback client data.
    //
    UICallbackStruct  *_clientDataStructs;
    
    //
    // Callbacks to interface with Motif.
    //
    static void vSelectCBCallback(Widget, XtPointer, XtPointer);
    
    // Begin user code block <private>
    // End user code block <private>
};

// Begin user code block <tail>
// End user code block <tail>
#endif