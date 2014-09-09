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
//
// README: Portions of this file are merged at file generation
// time. Edits can be made *only* in between specified code blocks, look
// for keywords <Begin user code> and <End user code>.
//
//////////////////////////////////////////////////////////////
//
// Header file for CUIheaderEdit
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
// Copyright 1997 Molecular Structure Corporation
//                3200 Research Forest Drive
//                The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved
//
// CUIheaderEdit.h  Initial author: T.L. Hendrixson          10-Apr-1997
//    Part of this file was automatically generated by ICS Builder Xcessory.
//    It is the header file for the CUIheaderEdit class.
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

#ifndef CUIheaderEdit_H
#define CUIheaderEdit_H
#include "UIComponent.h"

//
// Globally included information (change thru Output File Names Dialog).
//


//
// Class Specific Includes (change thru the class in Resource Editor).
//



// Begin user code block <head>

#include "Cstring.h"
#include "Cimage_header.h"
#include "Cspacegroup.h"

#ifdef NEED_BOOL
#include "Bool.h"
#endif

enum eKeyword_Action {
   eKeywordAction_None,
   eKeywordAction_Add,
   eKeywordAction_Delete,
   eKeywordAction_Replace
};

// End user code block <head>

class CUIheaderEdit : public UIComponent

// Begin user code block <superclass>
// End user code block <superclass>
{

// Begin user code block <friends>
// End user code block <friends>

  public:

    CUIheaderEdit(const char *, Widget);
    CUIheaderEdit(const char *);
    virtual ~CUIheaderEdit();
    virtual void create(Widget);
    const char *const  className();
    void set_lcrystorient3_labelString(String);
    
    // Begin user code block <public>

    bool bHasHeaderChanged();
    void vGetHeader(Cimage_header *);
    void vSetHeader(Cimage_header *);
    
    XtPointer m_pObj;
    void (*m_prvHeaderModNotifyCallback)(XtPointer pClassPtr, 
                                         Widget wBaseWidgetOfClass,
                                         XtPointer pBoolValue,
                                         XtPointer pCallData);

    void vSetHeaderModNotifyCallback(XtPointer pClassPtr,
                                     void (*function)(XtPointer, Widget,
                                                      XtPointer, XtPointer));

   static Cstring ms_sCrystComment;
   static Cstring ms_sCrystMosaicity;
   static Cstring ms_sCrystOrient;
   static Cstring ms_sCrystSpaceGroupNumber;
   static Cstring ms_sCrystUnitCell;
   static Cstring ms_sDetDimensions;
   static Cstring ms_sDetDirectBeam;
   static Cstring ms_sDetGonio;
   static Cstring ms_sDetGonioNames;
   static Cstring ms_sDetNames;
   static Cstring ms_sDetSpatialType;
   static Cstring ms_sDetSpatialInfo;
   static Cstring ms_sSourceRot;
   static Cstring ms_sSourceWavelength;

    // End user code block <public>
  protected:
    // Classes created by this class
    
    // Widgets created by this class
    Widget _CUIheaderEdit;
    Widget _dsGetValueError;
    Widget _mdGetValueError;
    Widget _dsValueError;
    Widget _mdValueError;
    Widget _dsReplaceValueError;
    Widget _mdReplaceValueError;
    Widget _dsAddKeywordError;
    Widget _mdAddKeywordError;
    Widget _dsDeleteKeywordError;
    Widget _mdDeleteKeywordError;
    Widget _dsKeywordLengthError;
    Widget _mdKeywordLengthError;
    Widget _dsKeywordCharacterError;
    Widget _mdKeywordCharacterError;
    Widget _dsProtectedKeywordWarning;
    Widget _mdProtectedKeywordWarning;
    Widget _foEdit;
    Widget _lTab;
    Widget _foPageMisc;
    Widget _pbDeleteKeyword;
    Widget _pbReplaceValue;
    Widget _pbAddKeyword;
    Widget _lValue;
    Widget _swValue;
    Widget _tValue;
    Widget _lHeader;
    Widget _tfKeyword;
    Widget _seMain;
    Widget _lKeyword;
    Widget _swHeader;
    Widget _liHeader;
    Widget _foPageCryst;
    Widget _frCrystComment;
    Widget _lCrystComment;
    Widget _swCrystComment;
    Widget _tCrystComment;
    Widget _frCrystSpaceGroup;
    Widget _lCrystSpaceGroup;
    Widget _rcCrystSpaceGroup;
    Widget _lCrystSpaceGroupName;
    Widget _tfCrystSpaceGroupName;
    Widget _lCrystSpaceGroupNumber;
    Widget _tfCrystSpaceGroupNumber;
    Widget _frCrystMosaicity;
    Widget _lCrystMosaicity;
    Widget _tfCrystMosaicity;
    Widget _frCrystOrient;
    Widget _lCrystOrient;
    Widget _rcCrystOrient;
    Widget _lCrystOrient1;
    Widget _lCrystOrient2;
    Widget _lCrystOrient3;
    Widget _tfCrystOrient1;
    Widget _tfCrystOrient2;
    Widget _tfCrystOrient3;
    Widget _frUnitCell;
    Widget _lCrystUnitCell;
    Widget _rcCrystUnitCell;
    Widget _lCrystUnitCellA;
    Widget _lCrystUnitCellB;
    Widget _lCrystUnitCellC;
    Widget _tfCrystUnitCellA;
    Widget _tfCrystUnitCellB;
    Widget _tfCrystUnitCellC;
    Widget _lCrystUnitCellAlpha;
    Widget _lCrystUnitCellBeta;
    Widget _lCrystUnitCellGamma;
    Widget _tfCrystUnitCellAlpha;
    Widget _tfCrystUnitCellBeta;
    Widget _tfCrystUnitCellGamma;
    Widget _foPageDet;
    Widget _frDetDirectBeam;
    Widget _rcDetDirectBeam;
    Widget _lDetDirectBeamX;
    Widget _tfDetDirectBeamX;
    Widget _lDetDirectBeamY;
    Widget _tfDetDirectBeamY;
    Widget _lDetDirectBeam;
    Widget _frDetRot;
    Widget _rcDetRot;
    Widget _lDetRot1;
    Widget _lDetRot2;
    Widget _lDetRot3;
    Widget _tfDetRot1;
    Widget _tfDetRot2;
    Widget _tfDetRot3;
    Widget _lDetRot;
    Widget _frDetTrans;
    Widget _rcDetTrans;
    Widget _lDetTrans1;
    Widget _lDetTrans2;
    Widget _lDetTrans3;
    Widget _tfDetTrans1;
    Widget _tfDetTrans2;
    Widget _tfDetTrans3;
    Widget _lDetTrans;
    Widget _foPageSource;
    Widget _frSourceRot;
    Widget _lSourceRot;
    Widget _rcSourceRot;
    Widget _tfSourceRot1;
    Widget _tfSourceRot2;
    Widget _frSourceWavelength;
    Widget _lSourceWavelength;
    Widget _tfSourceWavelength;
    Widget _pbPageMisc;
    Widget _pbPageCryst;
    Widget _pbPageDet;
    Widget _pbPageSource;
    
    // These virtual functions are called from the private callbacks 
    // or event handlers intended to be overridden in derived classes
    // to define actions
    
    virtual void ClassHelpCB(Widget, XtPointer, XtPointer);
    virtual void ChangeHeaderCB(Widget, XtPointer, XtPointer);
    virtual void DeleteKeywordCB(Widget, XtPointer, XtPointer);
    virtual void ReplaceValueCB(Widget, XtPointer, XtPointer);
    virtual void AddKeywordCB(Widget, XtPointer, XtPointer);
    virtual void KeywordChangedCB(Widget, XtPointer, XtPointer);
    virtual void KeywordSelectedCB(Widget, XtPointer, XtPointer);
    virtual void ChangeCommentCB(Widget, XtPointer, XtPointer);
    virtual void ChangeCrystSpaceGroup(Widget, XtPointer, XtPointer);
    virtual void ChangeCrystMosaicityCB(Widget, XtPointer, XtPointer);
    virtual void ChangeCrystOrientCB(Widget, XtPointer, XtPointer);
    virtual void ChangeUnitCellCB(Widget, XtPointer, XtPointer);
    virtual void ChangeDirectBeamCB(Widget, XtPointer, XtPointer);
    virtual void ChangeDetGonioCB(Widget, XtPointer, XtPointer);
    virtual void ChangeSourceRotCB(Widget, XtPointer, XtPointer);
    virtual void ChangeSourceWavelengthCB(Widget, XtPointer, XtPointer);
    virtual void ChangeMenuCB(Widget, XtPointer, XtPointer);
    
    // Begin user code block <protected>
    // End user code block <protected>
  private: 
    
    //
    // Default application and class resources.
    //
    static String         _defaultCUIheaderEditResources[];
    static UIAppDefault   _appDefaults[];
    static Boolean        _initAppDefaults;
    //
    // Callback client data.
    //
    UICallbackStruct  *_clientDataStructs;
    
    //
    // Callbacks to interface with Motif.
    //
    static void ClassHelpCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeHeaderCBCallback(Widget, XtPointer, XtPointer);
    static void DeleteKeywordCBCallback(Widget, XtPointer, XtPointer);
    static void ReplaceValueCBCallback(Widget, XtPointer, XtPointer);
    static void AddKeywordCBCallback(Widget, XtPointer, XtPointer);
    static void KeywordChangedCBCallback(Widget, XtPointer, XtPointer);
    static void KeywordSelectedCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeCommentCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeCrystSpaceGroupCallback(Widget, XtPointer, XtPointer);
    static void ChangeCrystMosaicityCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeCrystOrientCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeUnitCellCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeDirectBeamCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeDetGonioCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeSourceRotCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeSourceWavelengthCBCallback(Widget, XtPointer, XtPointer);
    static void ChangeMenuCBCallback(Widget, XtPointer, XtPointer);
    
    // Begin user code block <private>

    void vAddKeywordError(const Cstring &);
    void vAddStringToList(const Cstring & , const Widget);
    void vDeleteKeywordError(const Cstring &);
    void vErrorWithString(Widget , XmString , const Cstring &);
    void vFillKeyword(const Cstring &);
    void vFillValue(const Cstring &);
    void vGetValueError(const Cstring &);
    void vHeaderHasChanged(bool , XtPointer);
    void vKeywordCharacterError(const Cstring &);
    void vKeywordError(const Cstring &);
    void vKeywordLengthError(const Cstring & , int);
    void vReplaceValueError(const Cstring &);
    void vSetKeywordButtons();
    void vValueError(const Cstring &);
    XmString m_axsAddKeywordError[2];
    XmString m_axsDeleteKeywordError[2];
    XmString m_axsGetValueError[2];
    XmString m_axsKeywordCharacterError[2];
    XmString m_axsKeywordLengthError[2];
    XmString m_axsReplaceValueError[2];
    XmString m_axsValueError[2];
    bool m_bHeaderHasChanged;
    Cimage_header m_oHeader;
    Cimage_header m_oWorkingHeader;
    
   Cspacegroup m_oSpaceGroup;

   Widget m_awCrystUnitCell[6];
   Widget m_awCrystOrient[3];
   Widget m_awDetDirectBeam[2];
   Widget m_awDetGonio[6];
   Widget m_awDetGonioNames[6];
   Widget m_awSourceRot[2];

   XmString m_xsCrystSpaceGroupNumber;
   XmString m_axsDefaultDetGonioNames[6];

   Cstring sGetDetectorPrefix(void);
   void vGetDirectBeamPosition(Cstring sPrefix, float af[2]);

   int nSearchAlphaList(Widget wList, const Cstring &sItem);

// Update routines

   void vUpdateDisplayedPage(void);
   void vUpdatePageCryst(void);
   void vUpdatePageDet(void);
   void vUpdatePageMisc(void);
   void vUpdatePageSource(void);

// Convenience routines ala XmWidgetSetSomething()

   void XmLabelSetString(Widget w, char *pc);
   void XmLabelSetString(Widget w, Cstring &s);
   Cstring XmStringGetCstring(XmString xsString);
   int XmStringSplit(XmString xsString, Cstring sToken, int nMaxParts,
                     XmString *xsParts);

// routines for changing the working header

   int nChangeWorkingHeader(const enum eKeyword_Action eAction, 
                            const Cstring &sKeyword, const int nValue);
   int nChangeWorkingHeader(const enum eKeyword_Action eAction,
                            const Cstring &sKeyword, const float fValue);
   int nChangeWorkingHeader(const enum eKeyword_Action eAction,
                            const Cstring &sKeyword, const int nNumValues,
                            const float *pfValues);
   int nChangeWorkingHeader(const enum eKeyword_Action eAction,
                            const Cstring &sKeyword, 
                            const Cstring &sValue);
   int nChangeWorkingHeader(const enum eKeyword_Action eAction,
                            const Cstring &sKeyword, 
                            const char *pcValue = "dummy")
   { Cstring s = pcValue; return nChangeWorkingHeader(eAction,sKeyword,s); }


    // End user code block <private>
};

// Begin user code block <tail>

bool bKeywordProtected(const Cstring &sKeyword);

// End user code block <tail>
#endif
