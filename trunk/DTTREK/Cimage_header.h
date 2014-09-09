#ifndef DT_CIMAGE_HEADER_H
#define DT_CIMAGE_HEADER_H
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
// Cimage_header.h     Initial author: J.W. Pflugrath           03-Mar-1995
//    This file is the header file for class Cimage_header.
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

#include "Dtrek.h"
#include "Cstring.h"
#include "dtreksys.h"
#include "dtarray.h"

//+External prototypes

//+Definitions and constants
const int    MAX_NUMBER_LEN = 50;
const int    LEN_MULT       = 512;
const char   PAD_CHAR       = ' ';
const int    SIZE_FIELD     =  5;
const int    DIG_PRECISION  =  6;
//const int    MAX_HEADER_LENGTH = 512*195;

enum eImage_header_states {
  eImage_header_unknown_state,
  eImage_header_available_state,
  eImage_header_empty_state
};

enum eImage_data_types {
  eImage_other_type,
  eImage_byte,
  eImage_ubyte,
  eImage_I2,
  eImage_I4,
  eImage_uI2,
  eImage_uI4,
  eImage_realIEEE,
  eImage_compressed
};

enum eImage_format_types {
  enImage_format_unknown=-1,
  enImage_format_RIGAKU_CCD,
  enImage_format_RAXIS,
  enImage_format_DIP,
  enImage_format_MEDOPTICS,
  enImage_format_BRUKER_SIEMENS,
  enImage_format_WINBMP,
  enImage_format_ADSC,
  enImage_format_BRANDEIS,
  enImage_format_MAR_CCD,
  enImage_format_MAR_IP,
  enImage_format_CBFDTREK,
  enImage_format_BAS2000,
  enImage_format_HIPIC
};

//+Code begin

class DTREK_EXPORT Cimage_header {
    
private:
    std::map<Cstring,Cstring,std::less<Cstring> > m_oKeywords;
    
  Cstring           m_sFilename;            // Original filename
  int               m_nLength;

  static Cstring  ms_sEqual;                // Equal sign character
  static Cstring  ms_sSemiNew;              // Semicolon + Newline characters
  static Cstring  ms_sStarta;               // Part a of sStart
  static Cstring  ms_sStartb;               // Part b of sStart
  static Cstring  ms_sStart;                // Starting text of header
  static Cstring  ms_sEnd;                  // Ending   text of header

  int m_nFindKeyword;

public:

  bool m_bAllowEnvironmentOverride;

  eImage_header_states m_eThe_State;

  static Cstring ms_sComment;
  static Cstring ms_sComment2;
  static Cstring ms_sDtdisplayOrientation;
  static Cstring ms_sHeaderBytes;
  static Cstring ms_sByteOrder;
  static Cstring ms_sByteOrderOriginal;
  static Cstring ms_sCompression;
  static Cstring ms_sDim;
  static Cstring ms_sDataType;
  static Cstring ms_sFilename;
  static Cstring ms_sNumReflns;
  static Cstring ms_sSaturatedValue;
  static Cstring ms_sMinRawPixOKValue;
  static Cstring ms_sRaxis;
  static Cstring ms_sRxPrefix;
  static Cstring ms_sRaxisCompressionRatio;
  static Cstring ms_sRaxisReadLines;
  static Cstring ms_sRaxisReadStart;
  static Cstring ms_sRaxisRecordSize;
  static Cstring ms_sSize1;
  static Cstring ms_sSize2;
  static Cstring ms_sType;

  static Cstring ms_sShoeboxSize;
  static Cstring ms_sShoeboxOrig;
  static Cstring ms_sShoeboxDir;
  static Cstring ms_sCog1;
  static Cstring ms_sCog2;
  static Cstring ms_sFint;

  static int     ms_nVerbose;

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Cimage_header ();
Cimage_header (const Cimage_header& header);
Cimage_header (const Cstring& sFilename);
Cimage_header (const char *pcHeader);
~Cimage_header ();

Cimage_header&  operator= (const Cimage_header& header);

////////////////////////////////////////////////////////////////////////
//+Member function prototypes


// Get length of header including padding

int nGetImageSize(void);

eImage_data_types eGetImageDataType(void);

eImage_format_types enGetOriginalFormat();

static Cstring sGetFormatTypeName(eImage_format_types enType);
static eImage_format_types enGetFormatTypeByName(Cstring& sFormatName);

// Get state of the header

inline eImage_header_states eGetState(void)
   { return (m_eThe_State); }

// Set state of the header

inline void vSetState (const eImage_header_states eState)
  { m_eThe_State = eState; }

inline void vSetEnvironmentOverride(bool bValue)
{  m_bAllowEnvironmentOverride = bValue; return; }


int nGetDimensions(int *pnDim0, int *pnDim1);

// Get value of a keyword as Cstring, integer, float, double,
//    or arrays of integer, float, double

int nGetValue(const Cstring& sKeyword, Cstring* sValue);
int nGetValue(const Cstring& sKeyword, bool *bValue);
int nGetValue(const Cstring& sKeyword, int *nValue);
int nGetValue(const Cstring& sKeyword, long *lValue);
int nGetValue(const Cstring& sKeyword, float *fValue);
int nGetValue(const Cstring& sKeyword, double *dValue);
int nGetValue(const Cstring& sKeyword, const int nDim, Cstring sValue[]);
int nGetValue(const Cstring& sKeyword, const int nDim, bool   bValue[]);
int nGetValue(const Cstring& sKeyword, const int nDim, int    nValue[]);
int nGetValue(const Cstring& sKeyword, const int nDim, long   lValue[]);
int nGetValue(const Cstring& sKeyword, const int nDim, float  fValue[]);
int nGetValue(const Cstring& sKeyword, const int nDim, double fValue[]);

int nGetValue(const Cstring& sKeyword, std::vector<int>&    anVec, char* pcListSeparators);
int nGetValue(const Cstring& sKeyword, std::vector<double>& afVec, char* pcListSeparators);

// Add KEYWORD=VALUE; to header as Cstring, integer, float, double,
//    or as arrays of integer float, double; Specify digits to right of
//    decimal point for floats and doubles

int nAddValue(const Cstring& sKeyword, const Cstring& sValue);
int nAddValue(const Cstring& sKeyword, const bool&   bValue);
int nAddValue(const Cstring& sKeyword, const char *pcValue)
   { Cstring s = pcValue; return nAddValue(sKeyword,s); }
int nAddValue(const Cstring& sKeyword, const int&    nValue);
int nAddValue(const Cstring& sKeyword, const long&   lValue);
int nAddValue(const Cstring& sKeyword, const float&  fValue);
int nAddValue(const Cstring& sKeyword, const double& dValue);
int nAddValue(const Cstring& sKeyword, const int&    nValue, const int nWid);
int nAddValue(const Cstring& sKeyword, const long&   lValue, const int nWid);
int nAddValue(const Cstring& sKeyword, const float&  fValue, const int nDec);
int nAddValue(const Cstring& sKeyword, const double& dValue, const int nDec);
int nAddValue(const Cstring& sKeyword, const int nDim, const bool   bValue[]);
int nAddValue(const Cstring& sKeyword, const int nDim, const int    nValue[]);
int nAddValue(const Cstring& sKeyword, const int nDim, const long   lValue[]);
int nAddValue(const Cstring& sKeyword, const int nDim, const float  fValue[]);
int nAddValue(const Cstring& sKeyword, const int nDim, const double dValue[]);
int nAddValue(const Cstring& sKeyword, const int nDim, const Cstring sValue[]);
int nAddValue(const Cstring& sKeyword, const int nDim, const int    nValue[],
              const int nDec);
int nAddValue(const Cstring& sKeyword, const int nDim, const float  fValue[],
              const int nDec);
int nAddValue(const Cstring& sKeyword, const int nDim, const double dValue[],
              const int nDec);

// Replace KEYWORD=VALUE; in header as Cstring, integer, float, double,
//    or as arrays of integer float, double; Specify digits to right of
//    decimal point for floats and doubles

int nReplaceValue(const Cstring& sKeyword, const Cstring& sValue);
int nReplaceValue(const Cstring& sKeyword, const bool&   bValue);
int nReplaceValue(const Cstring& sKeyword, const char *pcValue)
   { Cstring s = pcValue; return nReplaceValue(sKeyword,s); }
int nReplaceValue(const Cstring& sKeyword, const int&    nValue);
int nReplaceValue(const Cstring& sKeyword, const long&   lValue);
int nReplaceValue(const Cstring& sKeyword, const float&  fValue);
int nReplaceValue(const Cstring& sKeyword, const double& dValue);
int nReplaceValue(const Cstring& sKeyword, const int&    nValue, const int nWid);
int nReplaceValue(const Cstring& sKeyword, const long&   lValue, const int nWid);
int nReplaceValue(const Cstring& sKeyword, const float&  fValue, const int nDec);
int nReplaceValue(const Cstring& sKeyword, const double& dValue, const int nDec);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const bool   bValue[]);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const int    nValue[]);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const long   lValue[]);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const float  fValue[]);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const double dValue[]);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const Cstring sValue[]);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const int    nValue[],
              const int nDec);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const float  fValue[],
                  const int nDec);
int nReplaceValue(const Cstring& sKeyword, const int nDim, const double dValue[],
                  const int nDec);

int nReplaceValueDictionary(Cstring sDictionary,Cstring sName,double fValue,int nIndex = -1,Cstring sPrefix = "");
int nReplaceValueDictionary(Cstring sDictionary,Cstring sName,Cstring sValue,int nIndex = -1,Cstring sPrefix = "");
int nReplaceValueDictionary(Cstring sDictionary,Cstring sName,int nValue,int nIndex = -1,Cstring sPrefix = "");
int nGetValueDictionary(Cstring sDictionary,Cstring sName,double& fValue,int nIndex = -1,Cstring sPrefix = "");
int nGetValueDictionary(Cstring sDictionary,Cstring sName,Cstring& sValue,int nIndex = -1,Cstring sPrefix = "");
int nGetValueDictionary(Cstring sDictionary,Cstring sName,int& nValue,int nIndex = -1,Cstring sPrefix = "");
int nGetDictionaryNames(Cstring sDictionary, std::vector<Cstring>& asNames);
int nClearDictionary(Cstring sDictionary,Cstring sMask);

// Delete KEYWORD and value from image_header

int nDelete(const Cstring& sKeyword);

// Delete KEYWORDs and values that match Keyword mask

int nDeleteMask(const Cstring& sKeywordMask);

// Verbosity methods

static void vSetVerbose(const int nVerbose);
static int  nGetVerbose(void);

// Copy KEYWORDs and values from text of another header

int nCopyMask(Cimage_header& oHeaderIn, const Cstring& sKeywordMask);
int nCopyMask(const Cstring& sText, const Cstring& sKeywordMask);
int nCopyMask(Cimage_header& oHeaderIn,Cstring& sPrefixIn,Cstring& sPrefixOut);

// Empty the header

int nEmpty(void);

// List header on standard output

int nList(void);

// Read the header from an already open file

int nRead(int* pnFile, const Cstring& rsFilename = "", 
	  const Cstring& rsOrigFilename="");

// Write the header to an already open file

int nWrite(int* pnFile);

// Write the header to a file (open and close the file)

int nWrite(const Cstring& sFilename);

int nUnbinHeader(void);

inline void vSetFilename(const Cstring& rsFilename) { m_sFilename = rsFilename; }
inline bool bIsAvailable(void)
{  return (eImage_header_available_state == m_eThe_State); }


// function below returns 0 if keyword found, -1 if at end of header,
// and +ve number if something went wrong (sCurrentKeyword not in header, etc)

bool bKeywordExists(const Cstring &sKeyword);

void vHeapSort(int nObjects, Cstring **ppsObjects);
int  nParse(Cstring& sText,bool bBuild = false);
int  nLength() { return m_nLength; };
int  nFindKeywordMask(const Cstring& sMask, std::vector<Cstring>& asKeywords);

// These functions should be removed after the new changes have been tested.
int  nBuildText(Cstring& sTemp) { return nParse(sTemp,true); };
Cstring sGet() { Cstring sTemp; nParse(sTemp,true); return sTemp;};

/*
 * Some boolean checks for keywords and values.
 */
bool bIsKeywordValid(const Cstring &sKeyword);
bool bIsValueValid(const Cstring &sValue);
bool bTextHasEnd(const Cstring& sText);


private:

// Convert int, float, double to Cstring
// changed to use MSCString creation operators - tlh
Cstring    sIntToCstring(const int n){return (MSCString(n));}
Cstring    sIntToCstring(const int n,    const int nWid)
                        {return (MSCString(n,(int)nWid));}
Cstring  sFloatToCstring(const float f){return (MSCString(f));}
Cstring  sFloatToCstring(const float f,  const int nDec)
                        {return (MSCString(f,0,nDec));}
Cstring sDoubleToCstring(const double d){return (MSCString(d));}
Cstring sDoubleToCstring(const double d, const int nDec)
                        {return (MSCString(d,0,nDec));}

};  // end of class Cimage_header

#endif   // DT_CIMAGE_HEADER_H
