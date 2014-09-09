//*****************************************************************
//
// Copyright (c) 2004 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTCoreTclParam.cc   Initial author: RB           05-APR-2004
// This file is a DTREK version of SSI core file  CoreTclParam.h

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

//////////////////////////////////////////////////////////////////////////////////
//                       CoreTclParam.h
//////////////////////////////////////////////////////////////////////////////////
/*
The TclParam class is used for handling the SSI standard parameter format
for Tcl commands.  It can be used both to construct the Tcl command in the
proper format and to parse Tcl commands that have come through the Tcl interpreter.

To construct a Tcl command your code should look something like:

     CCoreTclParam pTclCom ("Command");
     pTclCom.SetParamString ("Param1", "String1");
     pTclCom.SetParamInt ("Param2", 14);
     ...
     CString command = pTclCom.ConstructCommand();

When a command has been parsed by the Tcl interpreter it will have an argc, argv
pair.  To use CCoreTclParam to get the parameters for this command your code should
do the following:

     CCoreTclParam pTclCom (argc, argv); //These are the argc argv from Tcl
     CString p1String = pTclCom.GetParamString ("Param1", found);
     if (found)
          ...
     int p2Int = pTclCom.GetParamInt ("Param2", found, valid);
     if (found && valid)
        ....
*/

#ifndef DT_TCLPARAM_H
#define DT_TCLPARAM_H

//    #define NO_MFC

#include <ctype.h>

#if !defined(WIN32) && !defined(__APPLE__)
    #include <iostream.h>
    #include <iomanip.h>
    #include <string.h>
#endif // !WIN32

//#ifndef NO_MFC
//#include <afx.h>

//#else  // NO_MFC
#include "Cstring.h"

#define CString MSCString

#ifndef TRUE
    #define TRUE 1
    #define FALSE 0
#endif


#define TRACE printf
//#endif    // NO_MFC

#ifndef WIN32
#define _T
#endif

#ifndef ON
#define ON 1
#define OFF 0
#endif

//////// NOTE: Any Class definitions or stuctures that are exported should be
// defined after the pragma pack (8)
#ifdef USE_CORE_DLL      // Using a core dll 
	#ifdef EXPORT_CORE_DLL	// exporting the core DLL
		#define CORE_EXPORT _declspec(dllexport)
	#else					// Importing the core dll
		#define CORE_EXPORT _declspec(dllimport)
	#endif
     #pragma pack (8)    // Default packing size
#else
     #define CORE_EXPORT // Core code linked statically
#endif

struct NameValue {
   NameValue() : next(NULL),name(""),value("")
   {}

     NameValue*     next;
     CString        name;
     CString        value;
};

// Functions for handling brace characters in strings
int TclHideSpecialChars( char* Strng );
int TclShowSpecialChars( char* Strng );
int TclHideChars( char* Strng );

class CORE_EXPORT CCoreTclParam {

public:

private:
    CString         tclCommand;
     NameValue*     paramListHead;
     NameValue*     paramListLast;


////////////////////////////////////////////////////////////////////////
//        Constructors and Destructors
public:
CCoreTclParam ();        // default constructor
CCoreTclParam (const CString command);
//This constructor is used to parse the command list from a TCL command 
//that follows SSI standards.  Since the type of the parameters is not known
//at this time, validity checking is done when the parameters are accessed.


CCoreTclParam (const int argc, char* argv[]);

~CCoreTclParam ();

////////////////////////////////////////////////////////////////////////
//        Member Function prototypes
public:
     CString GetCommand();
     bool SetCommand( const CString CommandName);

// These methods are used to define a parameter list
bool SetParamString (const CString sName, const CString sValue);
bool SetParamChar (const CString sName, const char cValue);
bool SetParamDouble (const CString sName, const double dValue, const int precision = 5);
bool SetParamFloat (const CString sName, const float fValue, const int precision = 5);
bool SetParamInt (const CString sName, const int iValue);
bool SetParamShort (const CString sName, const short sValue);
bool SetParamAddress (const CString sName, const unsigned int pObj);
bool SetParamBool (const CString sName, const CString sValue);
bool SetParamBool (const CString sName, const bool bValue);

// This method is used to generate the command
CString ConstructCommand (void);

//These methods are used to access individual parameters
CString GetParamString (const CString sName, bool& found);
char GetParamChar (const CString sName, bool& found);

//Note that these functions check the validity of the number.
double GetParamDouble (const CString sName, bool& found, bool &valid);
float GetParamFloat (const CString sName, bool& found, bool &valid);
int GetParamInt (const CString sName, bool& found, bool &valid);
short GetParamShort (const CString sName, bool& found, bool &valid);
unsigned int GetParamAddress (const CString sName, bool& found, bool &valid);
bool GetParamBool (const CString sName, bool& found, bool &valid);

//These routines can be used if you want to loop over the actual prameters rather
//than the expected parameters.  Note that index is 0 based.
int CountParams (void);
CString GetParamName(int index);

//This is for debugging use only
void PrintParamList(void);

private:
//char * copyString (const CString s);
NameValue* FindNameValue (const CString sName);
NameValue* NewNameValue (const CString sName);
bool IsValidDouble (const char * dString);
bool IsValidInt (const char * dString);
bool IsValidAddress (const char * dString);
bool IsHexDigit (const char x);

};  // End of class CCoreTclParam


#ifdef USE_CORE_DLL      // Using a core dll 
     #pragma pack ()     // restore value
#endif

#endif

