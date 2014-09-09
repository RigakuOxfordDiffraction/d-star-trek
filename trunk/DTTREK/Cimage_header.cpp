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
// Cimage_header.cc     Initial author: J.W. Pflugrath           03-Mar-1995
//    This file contains the member functions of class Cimage_header.
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
//  Cimage_header.cc implements the member functions of the Cimage_header
//  class.
//
//  The image header format consists of an ASCII string that has a length
//  that is a multiple of LEN_MULT (=512)characters with a minimum length
//  of LEN_MULT and maximum length of an integer whose ASCII representation
//  fits in SIZE_FIELD (=5) digits.  PAD_CHAR (Space ' ') characters will be
//  used to pad the image_header so its length is a multiple of LEN_MULT
//  characters.
//  The string consists of human readable text in the form:
//    KEYWORD=VALUE;^J
//  KEYWORD is a case-sensitive string bounded by whitespace on the left,
//             with no whitespace within.  KEYWORD must begin with a letter or
//             underscore and can contain only letters, digits and underscores.
//             The length of KEYWORD no longer has a limit.
//  =       is the equal sign; there is no whitespace on the left side, but
//             may or may not be on the right.
//  VALUE   is the value of the keyword, which may be a number (integer or
//             float), a string or arrays of numbers.  A string
//             value cannot contain the {, } and ; characters.
//             Elements in arrays are separated by whitespace.
//  ;       is the semicolon character; it may or may not be preceded by
//             whitespace.
//  ^J      is the newline character (ASCII 10).
// Whitespace is any length sequence of space, tab, and/or newline characters.
//
// Other image_header requirements:
// 1. The header begins with the 2-character sequence { ^J
//      where ^J is ASCII 10.
// 2. The third character is the position of the start of the required keyword=
//      pair:
//        HEADER_BYTES=value;
//      which gives the length in bytes of the entire header and value is
//      ALWAYS 5 bytes in length, padded on left with spaces.
// 3. The unpadded part of the header always ends with the 4 character sequence:
//      } ^J ^L ^J
//      where ^J is ASCII 10 and ^L is ASCII 12.
//      With this one can use the Unix 'more' command to view image
//      headers and avoid viewing binary data.  Padding will occur after the
//      ^L ^J if necessary.
// 4. The following keywords are also required for images and must be placed in
//      the image header when an image is constructed:
//        DIM   = value;       The number of dimensions in the image (usually 2)
//        SIZE1 = value;       The number of pixels along the 1st direction.
//        SIZE2 = value;       The number of pixels along the 2nd direction.
//        TYPE  = mad;         If you want to read these images with MADNES;
//        BYTE_ORDER = big_endian | little_endian  The byte order of the data
//        Data_type = signed char | unsigned char | short int | long int |
//                    unsigned short int | unsigned lont int | float IEEE |
//                    Compressed | Other
//                    signed char        is signed one byte integer value
//                    unsigned char      is unsigned one byte integer value
//                    short int          is signed 2-byte integer value
//                    long int           is signed 4-byte integer value
//                    unsigned short int is unsigned 2-byte integer value
//                    unsigned long int  is signed 4-byte integer value
//                    float IEEE         is IEEE floating point values
//                    Compressed means the data is compressed.  Then there must
//                               be a COMPRESSION keyword for the algorithm.
//                    Other      is some other representation which will make
//                               the data unusable by most programs.  There
//                               should be an DATA_OTHER keyword with more
//                               information.
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files
#include <string.h>
#include "Cimage_header.h"      // Class definition and prototypes
#include "dskio.h"
#include "dtrekdefs.h"
#include "dtreksys.h"
#include "winbmp.h"

#ifdef DTREK_CBF
#include "cbfdtrek.h"
#endif
#ifdef RAXIS
#include "raxis.h"
#endif
#ifdef MEDOPTICS
#include "medoptics.h"
#endif
#ifdef DIP2030
#include "dip2030.h"
#endif
#ifdef SIEMENS
#include "bs.h"
#endif
#ifdef NONIUS
#include "nonius.h"
#endif
#ifdef MARCCD
#include "marccd.h"
#endif
#ifdef MARIP
#include "marip.h"
#endif
#ifdef ADSCCCD
#include "adscccd.h"
#endif
#ifdef BRANDEISCCD
#include "brandeisccd.h"
#endif
#ifdef BAS2000
#include "bas2000.h"
#endif
#ifdef HIPIC
#include "hipic.h"
#endif

#ifdef WIN32
#include <time.h>
#endif

using namespace std;


//+Code begin

//+Definitions, constants, and initialization of static member variables

Cstring Cimage_header::ms_sComment               = D_K_Comment;
Cstring Cimage_header::ms_sComment2              = D_K_Comment2;
Cstring Cimage_header::ms_sDtdisplayOrientation  = D_K_DtdisplayOrientation;
Cstring Cimage_header::ms_sHeaderBytes           = D_K_HeaderBytes;
Cstring Cimage_header::ms_sByteOrder             = D_K_ByteOrder;
Cstring Cimage_header::ms_sByteOrderOriginal     = D_K_ByteOrderOriginal;
Cstring Cimage_header::ms_sCompression           = D_K_Compression;
Cstring Cimage_header::ms_sDim                   = D_K_Dim;
Cstring Cimage_header::ms_sDataType              = D_K_DataType;
Cstring Cimage_header::ms_sFilename              = D_K_Filename;
Cstring Cimage_header::ms_sNumReflns             = D_K_NumReflns;
Cstring Cimage_header::ms_sRaxis                 = D_K_Raxis;
Cstring Cimage_header::ms_sRxPrefix              = D_K_RxPrefix;
Cstring Cimage_header::ms_sRaxisCompressionRatio = D_K_RaxisCompressionRatio;
Cstring Cimage_header::ms_sRaxisReadLines        = D_K_RaxisReadLines;
Cstring Cimage_header::ms_sRaxisReadStart        = D_K_RaxisReadStart;
Cstring Cimage_header::ms_sRaxisRecordSize       = D_K_RaxisRecordSize;
Cstring Cimage_header::ms_sSaturatedValue        = D_K_SaturatedValue;
Cstring Cimage_header::ms_sMinRawPixOKValue      = D_K_MinRawPixOKValue;
Cstring Cimage_header::ms_sSize1                 = D_K_Size1;
Cstring Cimage_header::ms_sSize2                 = D_K_Size2;
Cstring Cimage_header::ms_sType                  = D_K_Type;

Cstring Cimage_header::ms_sShoeboxSize           = D_K_ShoeboxSize;
Cstring Cimage_header::ms_sShoeboxOrig           = D_K_ShoeboxOrig;
Cstring Cimage_header::ms_sShoeboxDir            = D_K_ShoeboxDir;
Cstring Cimage_header::ms_sCog1                  = D_K_Cog1;
Cstring Cimage_header::ms_sCog2                  = D_K_Cog2;
Cstring Cimage_header::ms_sFint                  = D_K_Fint;

Cstring Cimage_header::ms_sEqual   = "=";          // Equal sign character
Cstring Cimage_header::ms_sSemiNew = ";\n";        // Semicolon + Newline
                                                   //    characters
Cstring Cimage_header::ms_sStarta  = "{\n";        // Part a of ms_sStart
Cstring Cimage_header::ms_sStartb  = Cimage_header::ms_sHeaderBytes + "=";   // Part b of ms_sStart
Cstring Cimage_header::ms_sStart   = "{\n" + Cimage_header::ms_sHeaderBytes + "=";// Starting text of header
Cstring Cimage_header::ms_sEnd     = "}\n\f\n";    // Ending   text of header

int Cimage_header::ms_nVerbose = 1;





//+Public functions

// Constructors, destructors and assignments

Cimage_header::Cimage_header ()
{
    
    nEmpty();
    m_sFilename    = "";
    m_nLength = 0;
    m_nFindKeyword = 0;
    m_eThe_State   = eImage_header_available_state;
    m_bAllowEnvironmentOverride = TRUE;
}

Cimage_header::Cimage_header (const Cstring& sFilenameIn)
{
    // Read only the header of an image file.
    // Open, read, and close the file within this constructor.
    
    int nFile;
    int nLength;
    int nStat;
    Cstring sFilename;
    Cstring sTemp;
    Cstring sUncompress;
    
    m_nFindKeyword = 0;
    
    nEmpty();
   
   
    nFile       = 1;
    sTemp       = sTransSymbol(sFilenameIn);
    m_sFilename = sTemp;
    m_nLength   = 0;
    m_bAllowEnvironmentOverride = TRUE;


    // If filename ends in .Z, uncompress the file to a temporary file first,
    // then read it in.
    
    nLength   = sTemp.length();
    sFilename = sTemp;
    if ( ('Z' == sTemp.GetAt(nLength-1)) && ('.' == sTemp.GetAt(nLength-2)) )
    {
#ifdef WIN32
        
#ifdef SSI_PC
        //sUncompress = sFileBuildName(GetInstallDirectory(),"gzip");
        //sUncompress += " -d -c ";
        sUncompress = "gzip -d -c ";
#else
        sUncompress = "gzip -d -c ";
#endif  // SSI_PC
        
#else
        sUncompress = "zcat ";
#endif  // WIN32
        sTemp       = sTemp.before(nLength-2);
    }
    else if (   ('z' == sTemp.GetAt(nLength-1))
        && ('.' == sTemp.GetAt(nLength-3)) && ('g' == sTemp.GetAt(nLength-2)))
    {
        //cout << "gunzip!\n";
#ifdef WIN32
        
#ifdef SSI_PC
        //sUncompress = sFileBuildName(GetInstallDirectory(),"gzip");
        //sUncompress += " -d -c ";
        sUncompress = "gzip -d -c ";
#else
        sUncompress = "gzip -d -c ";
#endif  // SSI_PC
        
#else
        sUncompress = "gunzip -c ";
#endif  // WIN32
        sTemp       = sTemp.before(nLength-3);
    }
  else if (   ('2' == sTemp.GetAt(nLength-1)) && ('z' == sTemp.GetAt(nLength-2))
           && ('b' == sTemp.GetAt(nLength-3)) && ('.' == sTemp.GetAt(nLength-4)))
    {
#ifdef WIN32

#ifdef SSI_PC
      //sUncompress = sFileBuildName(GetInstallDirectory(),"bzip2");
      //sUncompress += " -d -c ";
      sUncompress = "";  // Not supported on Windows (TM) for now.
#else
      sUncompress = "bunzip2 -k -c ";
#endif  // SSI_PC

#else
      sUncompress = "bunzip2 -k -c ";
#endif  // WIN32
      sTemp       = sTemp.before(nLength-4);
    }
    else
    {
        sUncompress = "";
    }
    
    if ("" != sUncompress)
    {
        // File is most likely compressed, so uncompress it to a temporary file.
        // TODO: probably should put temp file in the /tmp directory
        
        if ("" != sGetEnv("TMPDIR"))
        {
            // If TMPDIR is defined, then use that as the directory
            
            sTemp = sTransSymbol("$(TMPDIR)") + '/' + sFileGetBasename(sTemp);
        }
        
        if (bFileExists(sTemp))
            sTemp = sFileGetTempName(".", "IMG");
        
        if (0 < ms_nVerbose)
            cout << "...uncompressing " << sFilename << " to " << sTemp << "...\n";
#ifdef SSI_PC
        nDoSystemCommand2(sUncompress + "\"" + sFilename + "\"", sTemp);
#else
        nDoSystemCommand(sUncompress + sFilename + " > " + sTemp);
#endif
        if (0 < ms_nVerbose)
            cout << "...uncompress done." << endl;
        sFilename = sTemp;
    }
    
    nLength   = sFilename.length();
    
    (void) dskbor(&nFile, sFilename.string(), &nLength, &nStat);
    if (0 != nStat)
    {
        cout << "     Error opening file " << sFilename << " Error is: "
            << nStat << endl;
    }
    else
    {
        //cout << "     File " << sFilename << " successfully opened.\n";
        
        //  Next read the header;  What if it is not a d*TREK image???
        nStat = nRead(&nFile);
        cout << "     Header of file " << sFilename << " successfully read." << endl;
        

      (void) dskbcr(&nFile, &nLength);
      if ( (0 == nStat) && (0 != nLength) ) nStat = nLength; // Report close err
    }
  if (0 == nStat) 
      m_eThe_State = eImage_header_available_state;

  if ("" != sUncompress)
    (void) nFileDelete(sFilename);
}

Cimage_header::Cimage_header (const char *pcHeader)
{
  Cstring sText;

  sText = pcHeader;

  if (nParse(sText))
      m_eThe_State = eImage_header_unknown_state;
  else
    m_eThe_State = eImage_header_available_state;

  m_bAllowEnvironmentOverride = TRUE;

  (void) nGetValue(ms_sFilename, &m_sFilename);

  m_nFindKeyword = 0;
}

Cimage_header::Cimage_header (const Cimage_header &roHeader)
{

  m_bAllowEnvironmentOverride = TRUE;
  m_oKeywords = roHeader.m_oKeywords;
  m_eThe_State = roHeader.m_eThe_State;
  m_sFilename  = roHeader.m_sFilename;
  m_nLength    = roHeader.m_nLength;

  // Do we want these to be copied from the incoming roHeader?

  m_nFindKeyword = 0;
}

Cimage_header::~Cimage_header () {}

Cimage_header&
Cimage_header::operator= (const Cimage_header& roHeader)
{
  // Simple copy constructor
  m_oKeywords = roHeader.m_oKeywords;
  m_eThe_State = roHeader.m_eThe_State;
  m_sFilename  = roHeader.m_sFilename;
  m_nLength    = roHeader.m_nLength;

  // Do we want these to be copied from the incoming roHeader?

  return *this;
}
// Member functions

int
Cimage_header::nGetDimensions(int *pnDim0, int *pnDim1)
{
  int nStat;
  nStat = nGetValue(ms_sSize1, pnDim0);
  nStat = nStat + nGetValue(ms_sSize2, pnDim1);
  return (nStat);
}
// Get value of a keyword as Cstring, integer, float, double,
//    or arrays of Cstring, integer, float, double

int
Cimage_header::nGetValue(const Cstring& rsKeyword, Cstring *psValue)
{
  Cstring sTemp1, sTemp2;

  // First look for a symbol (aka environment variable or logical name)
  // for the sKeyWord.  If it exists, use that as an override of what
  // is in the header.

  if (m_bAllowEnvironmentOverride)
     sTemp1 = sGetEnv(rsKeyword);
  else
     sTemp1 = "";

  if ("" != sTemp1)
    {
      //      cout << "Overriding keyword: " << rsKeyword << " found: " << sTemp1 << '\n';
      *psValue = sTransSymbol(sTemp1);
      (void) nReplaceValue(rsKeyword, sTemp1);  // Edit the header a little
      return (0);
     }
  else
    {

      // There was no overriding symbol , so look in header

      std::map<Cstring,Cstring,std::less<Cstring> >::iterator oIt1;
      oIt1 = m_oKeywords.find(rsKeyword);
      if (oIt1 != m_oKeywords.end())
        {
          *psValue = (*oIt1).second;

          // Added 25-Jul-2002 by jwp

          *psValue = sTransSymbol(*psValue);
	  //+JWP 2008-07-25
	  // Remove leading spaces from the value
	  
	  //-JWP 2008-07-25
          return 0;
        }
    }
  return (2);
}

int
Cimage_header::nGetValue(const Cstring& sKeyword, bool *bValue)
{
  int     nError;
  Cstring sTemp;
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      sTemp.downcase();
      // Strip off leading whitespace at beginning of sTemp
      while (   (' '  == sTemp.GetAt(0))
             || ('\n' == sTemp.GetAt(0))
             || ('\t' == sTemp.GetAt(0)))
        sTemp = sTemp.after(0);
      if ('t' == sTemp.GetAt(0))
         *bValue = TRUE;
      else if ('f' == sTemp.GetAt(0))
         *bValue = FALSE;
      else
         return 1;
    }
  return (nError);
}

int
Cimage_header::nGetValue(const Cstring& sKeyword, int *nValue)
{
  int     nError;
  Cstring sTemp;
//  extern const Regex RXwhite;      // = "[ \n\t]+"
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      // Strip off leading whitespace at beginning of sTemp
//    if (sTemp.index(RXwhite) == 0) sTemp = sTemp.after(RXwhite);
//    nError = sscanf((const char*)sTemp, "%d", nValue);
      while (   (' '  == sTemp.GetAt(0))
             || ('\n' == sTemp.GetAt(0))
             || ('\t' == sTemp.GetAt(0)))
        sTemp = sTemp.after(0);
      nError = sscanf(sTemp.string(), "%d", nValue);
      if (1 == nError) return 0;
    }
  return (nError);
}

int
Cimage_header::nGetValue(const Cstring& sKeyword, long *nValue)
{
  int     nError;
  Cstring sTemp;
//  extern const Regex RXwhite;      // = "[ \n\t]+"
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      // Strip off leading whitespace at beginning of sTemp
//    if (sTemp.index(RXwhite) == 0) sTemp = sTemp.after(RXwhite);
//    nError = sscanf((const char*)sTemp, "%d", nValue);
      while (   (' '  == sTemp.GetAt(0))
             || ('\n' == sTemp.GetAt(0))
             || ('\t' == sTemp.GetAt(0)))
        sTemp = sTemp.after(0);
      nError = sscanf(sTemp.string(), "%ld", nValue);
      if (1 == nError) return 0;
    }
  return (nError);
}

int
Cimage_header::nGetValue(const Cstring& sKeyword, float *fValue)
{
  int     nError;
  Cstring sTemp;
//  extern const Regex RXwhite;      // = "[ \n\t]+"
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      // Strip off leading whitespace at beginning of sTemp
//    if (sTemp.index(RXwhite) == 0) sTemp = sTemp.after(RXwhite);
//    nError = sscanf((const char*)sTemp, "%f", fValue);
      while (   (' '  == sTemp.GetAt(0))
             || ('\n' == sTemp.GetAt(0))
             || ('\t' == sTemp.GetAt(0)))
        sTemp = sTemp.after(0);
      nError = sscanf(sTemp.string(), "%f", fValue);
      if (1 == nError) return 0;
    }
  return nError;
}


int
Cimage_header::nGetValue(const Cstring& sKeyword, double *dValue)
{
  int nError;
  Cstring sTemp;
//  extern const Regex RXwhite;      // = "[ \n\t]+"
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      // Strip off leading whitespace at beginning of sTemp

//    if (sTemp.index(RXwhite) == 0) sTemp = sTemp.after(RXwhite);
      while (   (' '  == sTemp.GetAt(0))
             || ('\n' == sTemp.GetAt(0))
             || ('\t' == sTemp.GetAt(0)))
        sTemp = sTemp.after(0);
      //    nError = sscanf((const char*)sTemp, "%lf", dValue);
      nError = sscanf(sTemp.string(), "%lf", dValue);
      if (1 == nError) return 0;
  }
  return nError;
}

int
Cimage_header::nGetValue(const Cstring& sKeyword, const int nDim,
                         Cstring sValue[])
{
  int     i;
  int     nError;
  Cstring sTemp;
//  extern const Regex RXwhite;      // = "[ \n\t]+"
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      for (i = 0; (i < nDim) && (sTemp != ""); i++)
        {
          // Strip whitespace at beginning of sTemp if necessary

          while (   (' '  == sTemp.GetAt(0))
                 || ('\n' == sTemp.GetAt(0))
                 || ('\t' == sTemp.GetAt(0)))
            sTemp = sTemp.after(0);

          // Now look for a space

          if (sTemp.contains(' '))
            {
              sValue[i] = sTemp.before(' ');
              sTemp = sTemp.after(' '); // Could turn out to be a null string
            }
          else
            {
              sValue[i] = sTemp;
              sTemp = "";  // Set sTemp to null string, because nothing left
            }
        }
      return (i - nDim);
    }
  return nError;
}


int
Cimage_header::nGetValue(const Cstring& sKeyword,
                         const int nDim, bool bValue[])
{
  int     i;
  int     nError;
  Cstring sTemp;
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      sTemp.downcase();
      nError = 1;
      for (i = 0; i < nDim && 1 == nError; i++)
        {
          // Strip off leading whitespace at beginning of sTemp

          while (   (' '  == sTemp.GetAt(0))
                 || ('\n' == sTemp.GetAt(0))
                 || ('\t' == sTemp.GetAt(0)))
            sTemp = sTemp.after(0);
          if ('t' == sTemp.GetAt(0))
             bValue[i] = TRUE;
          else if ('f' == sTemp.GetAt(0))
             bValue[i] = FALSE;
          else
             return 1;
          sTemp = sTemp.after(' ');
        }
      if (1 == nError) return 0;
    }
  return nError;
}

int
Cimage_header::nGetValue(const Cstring& sKeyword,
                         const int nDim, int nValue[])
{
  int     i;
  int     nError;
  Cstring sTemp;
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      nError = 1;
      for (i = 0; i < nDim && 1 == nError; i++)
        {
          // Strip off leading whitespace at beginning of sTemp

          while (   (' '  == sTemp.GetAt(0))
                 || ('\n' == sTemp.GetAt(0))
                 || ('\t' == sTemp.GetAt(0)))
            sTemp = sTemp.after(0);

          nError = sscanf(sTemp.string(), "%d", &nValue[i]);
          sTemp = sTemp.after(' ');
        }
      if (1 == nError) return 0;
    }
  return nError;
}

int
Cimage_header::nGetValue(const Cstring& sKeyword,
                         const int nDim, long nValue[])
{
  int     i;
  int     nError;
  Cstring sTemp;
  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      nError = 1;
      for (i = 0; i < nDim && 1 == nError; i++)
        {
          // Strip off leading whitespace at beginning of sTemp

          while (   (' '  == sTemp.GetAt(0))
                 || ('\n' == sTemp.GetAt(0))
                 || ('\t' == sTemp.GetAt(0)))
            sTemp = sTemp.after(0);

          nError = sscanf(sTemp.string(), "%ld", &nValue[i]);
          sTemp = sTemp.after(' ');
        }
      if (1 == nError) return 0;
    }
  return nError;
}

int
Cimage_header::nGetValue(const Cstring& sKeyword,
                         const int nDim, float fValue[])
{
  int     i;
  int     nError;
  Cstring sTemp;

  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      nError = 1;
      for (i = 0; i < nDim && 1 == nError; i++)
        {
          while (   (' '  == sTemp.GetAt(0))
                 || ('\n' == sTemp.GetAt(0))
                 || ('\t' == sTemp.GetAt(0)))
            sTemp = sTemp.after(0);
          nError = sscanf(sTemp.string(), "%f", &fValue[i]);
          sTemp = sTemp.after(' ');
        }
      if (1 == nError) return 0;
    }
  return nError;
}

int Cimage_header::nGetValue(const Cstring& sKeyword,
                             const int nDim, double dValue[])
{
  int     i;
  int     nError;
  Cstring sTemp;

  nError = nGetValue(sKeyword, &sTemp);
  if (0 == nError)
    {
      nError = 1;
      for (i = 0; i < nDim && 1 == nError; i++)
        {
          // Strip off leading whitespace at beginning of sTemp

      while (   (' '  == sTemp.GetAt(0))
             || ('\n' == sTemp.GetAt(0))
             || ('\t' == sTemp.GetAt(0)))
          while (   (0 == sTemp.index(' '))
                 || (0 == sTemp.index('\n'))
                 || (0 == sTemp.index('\t')))
            sTemp = sTemp.after(0);

          nError = sscanf(sTemp.string(), "%lf", &dValue[i]);
          sTemp = sTemp.after(' ');
        }
      if (1 == nError) return 0;
    }
  return nError;
}

// Add KEYWORD=VALUE; to header as integer, float, double, string
//    or as arrays of integer float, double; Specify digits to right of
//    decimal point for floats and doubles


int Cimage_header::nAddValue(const Cstring& sKeyword, const bool &bValue)
{
   Cstring s;

   if (TRUE == bValue)
      s = 't';
   else
      s = 'f';
   return (nAddValue(sKeyword, s));
}


int Cimage_header::nAddValue(const Cstring& sKeyword, const int& nValue)
{
  return (nAddValue(sKeyword, sIntToCstring(nValue)));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const int& nValue,
                             const int nWid)
{
  return (nAddValue(sKeyword, sIntToCstring(nValue, nWid)));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const long& lValue)
{
  Cstring s;
  s = lValue;
  return (nAddValue(sKeyword, s));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const long& lValue,
                             const int nWid)
{
  Cstring s;
  s.assign(lValue,nWid);
  return (nAddValue(sKeyword, s));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const float& fValue)
{
  return (nAddValue(sKeyword, sFloatToCstring(fValue)));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const float& fValue,
                             const int nDec)
{
  return (nAddValue(sKeyword, sFloatToCstring(fValue, nDec)));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const double& dValue)
{
  return (nAddValue(sKeyword, sDoubleToCstring(dValue)));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const double& dValue,
                             const int nDec)
{
  return (nAddValue(sKeyword, sDoubleToCstring(dValue, nDec)));
}

int Cimage_header::nAddValue(const Cstring& sKeyword, const Cstring& sValue)
{
// Valid characters for the keyword are:
//     a-z, A-Z, 0-9, _
  char *pc,*pcEnd;
  Cstring sTemp;
  Cstring sCurrKey,sNextKey;
    //  extern const Regex RXwhite;      // = "[ \n\t]+"
    //  extern const Regex RXidentifier; // = "[A-Za-z_][A-Za-z0-9_]*"
    //  Regex  RXinvalidchars = ("[{};]");
    //
    //Test sKeyword first for invalid characters.
    //  if (sKeyword.matches(RXidentifier) == 0) return (1);
   pc    = sKeyword.string();
   pcEnd = pc + sKeyword.length();
   while (   ('a' <= *pc && 'z' >= *pc)
          || ('A' <= *pc && 'Z' >= *pc)
          || ('_' ==   *pc)
          || ('0' <= *pc && '9' >= *pc) )
     pc++;
   if (pc != pcEnd)
     return (1);

  // Test if sKeyword is too short

  if (1 > sKeyword.length())
    {
      cout << "ERROR in Cimage_header, sKeyword too short!\n";
      return (2);
    }

  
  //Test if sValue contains any invalid characters
  //  if (sValue.contains(RXinvalidchars) == 1) return (4);
  pc    = sValue.string();
  pcEnd = pc + sValue.length();
  while(   ('\0' != *pc)
        //&& ('{'  != *pc)
        //&& ('}'  != *pc)
        && (';'  != *pc) )
    pc++;
  if (pc != pcEnd)
    return(4);

  // Everything is OK, add it.

  m_oKeywords.insert(std::pair<Cstring,Cstring>(sKeyword,sValue));
  return 0;
};


int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const bool bValue[])
{ int i;
  Cstring sTemp;
  if (TRUE == bValue[0])
     sTemp = 't';
  else
     sTemp = 'f';
  for (i = 1; i < nDim; i++){
    sTemp += " ";
    if (TRUE == bValue[0])
       sTemp += 't';
    else
       sTemp += 'f';
  }
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const int nValue[])
{ int i;
  Cstring sTemp;
  sTemp = sIntToCstring(nValue[0]);
  for (i = 1; i < nDim; i++)
    sTemp = sTemp + " " + sIntToCstring(nValue[i]);
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const long lValue[])
{ int i;
  Cstring sTemp;
  sTemp = lValue[0];
  for (i = 1; i < nDim; i++){
    sTemp += ' ';
    sTemp += lValue[i];
  }
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const float fValue[])
{ int i;
  Cstring sTemp;
  sTemp = sFloatToCstring(fValue[0]);
  for (i = 1; i < nDim; i++)
    sTemp = sTemp + " " + sFloatToCstring(fValue[i]);
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const double dValue[])
{ int i;
  Cstring sTemp;
  sTemp = sDoubleToCstring(dValue[0]);
  for (i = 1; i < nDim; i++)
    sTemp = sTemp + " " + sDoubleToCstring(dValue[i]);
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const Cstring sValue[])
{ int i;
  Cstring sTemp;
  sTemp = sValue[0];
  for (i = 1; i < nDim; i++)
    sTemp = sTemp + " " + sValue[i];
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const int nValue[],    const int nWid)
{ int i;
  Cstring sTemp;
  sTemp = sIntToCstring(nValue[0], nWid);
  for (i = 1; i < nDim; i++)
    sTemp = sTemp + " " + sIntToCstring(nValue[i], nWid);
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const float fValue[],  const int nDec)
{ int i;
  Cstring sTemp;
  sTemp = sFloatToCstring(fValue[0], nDec);
  for (i = 1; i < nDim; i++)
    sTemp = sTemp + " " + sFloatToCstring(fValue[i], nDec);
  return (nAddValue(sKeyword, sTemp));
}

int
Cimage_header::nAddValue(const Cstring& sKeyword, const int nDim,
                         const double dValue[], const int nDec)
{ int i;
  Cstring sTemp;
  sTemp = sDoubleToCstring(dValue[0], nDec);
  for (i = 1; i < nDim; i++)
    sTemp = sTemp + " " + sDoubleToCstring(dValue[i], nDec);
  return (nAddValue(sKeyword, sTemp));
}

// Replace KEYWORD=VALUE; to header as integer, float, double, string
//    or as arrays of integer float, double; Specify digits to right of
//    decimal point for floats and doubles


int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const bool& bValue)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, bValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int& nValue)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int& nValue,
                             const int nWid)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nValue, nWid));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const long& lValue)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, lValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const long& lValue,
                             const int nWid)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, lValue, nWid));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const float& fValue)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, fValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const float& fValue,
                             const int nDec)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, fValue, nDec));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const double& dValue)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, dValue));
}

int Cimage_header::nReplaceValue(const Cstring& sKeyword, const double& dValue,
                                 const int nDec)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, dValue, nDec));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const Cstring& sValue)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, sValue));
}


int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const bool bValue[])
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, bValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const long lValue[])
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, lValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const int nValue[])
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, nValue));
}

int Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                                 const float fValue[])
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, fValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const double dValue[])
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, dValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const Cstring sValue[])
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, sValue));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const int nValue[],    const int nDec)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, nValue, nDec));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const float fValue[],  const int nDec)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, fValue, nDec));
}

int
Cimage_header::nReplaceValue(const Cstring& sKeyword, const int nDim,
                             const double dValue[], const int nDec)
{
  // Delete all occurences of the keyword

  while (0 == nDelete(sKeyword));
  return (nAddValue(sKeyword, nDim, dValue, nDec));
}

// Delete KEYWORD

int
Cimage_header::nDelete(const Cstring& sKeyword)
{
    if (m_oKeywords.erase(sKeyword) < 1)
        return 1;
    return 0;
}


int Cimage_header::nDeleteMask(const Cstring &sKeywordMask)
{
    std::vector<Cstring>          asKeywords;
    
    if( nFindKeywordMask(sKeywordMask, asKeywords) )
        return 1;

    for(int nx=0; nx < asKeywords.size(); nx++) 
        m_oKeywords.erase(asKeywords[nx]);
    
    return 0;
}
    

int
Cimage_header::nCopyMask(const Cstring &sHeaderIn, const Cstring &sKeywordMask)
{
   Cimage_header oHeader(sHeaderIn.string());
   return nCopyMask(oHeader,sKeywordMask);
}

int Cimage_header::nCopyMask(Cimage_header& oHeaderIn,Cstring& sPrefixIn,Cstring& sPrefixOut)
{
    Stringreg oReg;
    Cstring sKeywordMask;
    Cstring sKeywordMaskNoStar;
    Cstring sKeywordMaskNoStarEnd;
    Cstring sMask;
    Cstring sTemp,sTemp2;
    
    std::map<Cstring,Cstring,std::less<Cstring> >::iterator oIt1;
    std::map<Cstring,Cstring,std::less<Cstring> >::iterator oIt2;
        
    sMask = sPrefixIn;
    
    // Parse out next ';'
    while (sMask.length()) {
        sKeywordMask = sMask.before(';');
        
        sKeywordMaskNoStar = sKeywordMask;
        sKeywordMaskNoStar = sKeywordMaskNoStar.before('*');
        oIt1 = oHeaderIn.m_oKeywords.lower_bound(sKeywordMaskNoStar);
        sKeywordMaskNoStarEnd += '\xff';
        oIt2 = oHeaderIn.m_oKeywords.upper_bound(sKeywordMaskNoStarEnd);            
        while (oIt1 != oIt2) {
            if (!oReg.nParse((*oIt1).first,sKeywordMask)) {
                sTemp = sPrefixOut;
                sTemp2 = (*oIt1).first;
                if (sKeywordMaskNoStar.length())
                    sTemp += sTemp2.after(sKeywordMaskNoStar);
                else
                    sTemp += sTemp2;
                nReplaceValue(sTemp,(*oIt1).second);                
            };
            ++oIt1;
        };
        sMask = sMask.after(';');
    }; 
    return (0);

};


int Cimage_header::nFindKeywordMask(const Cstring& _sMask, std::vector<Cstring>& asKeywords)
{
    Stringreg oReg;
    Cstring sKeywordMask;
    Cstring sKeywordMaskNoStar;
    Cstring sMask;
    
    std::map<Cstring,Cstring,std::less<Cstring> >::iterator oIt1;
    std::map<Cstring,Cstring,std::less<Cstring> >::iterator oIt2;
    
    asKeywords.clear();
    sMask = _sMask;
    
    // Parse out next ';'
    while (sMask.length()) {
        sKeywordMask = sMask.before(';');
        
        sKeywordMaskNoStar = sKeywordMask;
        sKeywordMaskNoStar = sKeywordMaskNoStar.before('*');
        oIt1 = m_oKeywords.lower_bound(sKeywordMaskNoStar);
        sKeywordMaskNoStar += '\xff';
        oIt2 = m_oKeywords.upper_bound(sKeywordMaskNoStar);            
        while (oIt1 != oIt2) 
        {
            if (!oReg.nParse((*oIt1).first,sKeywordMask)) 
            {
                asKeywords.push_back((*oIt1).first);
            }
            
            ++oIt1;
        }
        sMask = sMask.after(';');
    } 

    return (0);
}


int Cimage_header::nCopyMask(Cimage_header &oHeaderIn, const Cstring &sKeywordMask)
{
    std::vector<Cstring>        asKeywords;
    Cstring                     sData;

    if( oHeaderIn.nFindKeywordMask(sKeywordMask, asKeywords) )
        return 1;
    
    for(int nx=0; nx < asKeywords.size(); nx++)
    {
        if (oHeaderIn.nGetValue(asKeywords[nx], &sData))
            return 1;
        
        nReplaceValue(asKeywords[nx], sData);
    }
    
    return 0;
}

// Empty the header
int
Cimage_header::nEmpty(void)
{
    m_oKeywords.clear();
  m_eThe_State    = eImage_header_empty_state;
  return (0);
}

// Read the header from an already open file

bool
Cimage_header::bTextHasEnd(const Cstring& sText) {
    int nPos;
    char cChar;
    const int nMaxNonWhiteCount = 2;
    int anNonWhite[nMaxNonWhiteCount];
    int nNonWhiteCount;
    nPos = sText.length() - 1;
    nNonWhiteCount = 0;
    while ((nNonWhiteCount < nMaxNonWhiteCount) && (nPos>=0)) {
        cChar = sText.GetAt(nPos--);
        if ((cChar != ' ') && (cChar != '\t') && (cChar!='\r') && (cChar!='\n'))
            anNonWhite[nNonWhiteCount++] = (int) cChar;
    };
    if ((anNonWhite[0] == '\f') && (anNonWhite[1] == /* { */ '}'))
        return true;
    return false;    
};

int
Cimage_header::nRead(int* pnFile, const Cstring& rsFilename,
                     const Cstring& rsOrigFilename)
{
  int   nStat;
  int   nLength;
  int   nHeaderBytes;
  char *pcBuffer;
  int i;
  int nLen_Mult; // Kludge for keeping NT compiler happy on dskbr() calls
                 // as cannot convert from const int* to int*

  int nDim0 = -1;
  int nDim1 = -1;
  long lFileSize;
  Cstring sTemp;
  Cstring sText;
  Cstring sNameForTemplate;
  bool bIsDtrek;
#ifdef BAS2000
  Cstring sNameForBAS2000Inf = "";//BAS2000
#endif
  if (100 >= LEN_MULT)
    {
      cout << "     Error in Cimage_header::nRead, bad LEN_MULT.\n";
      return (-2);
    }

  static int s_nTimes = 0;
  if (0 == (s_nTimes % 50))
    {
    }
  s_nTimes++;

  if ("" != rsFilename)
    {
      m_sFilename = rsFilename;
    }
  
  sNameForTemplate = m_sFilename;
  if ( ("" != rsOrigFilename) && (rsOrigFilename != rsFilename) )
    {
      sNameForTemplate = rsOrigFilename;
    }
#ifdef BAS2000
	if(-1 != sNameForTemplate.find(".img"))
	{
		sNameForBAS2000Inf = sNameForTemplate.before(".img");
		sNameForBAS2000Inf += ".inf";
	}
	else if(-1 != sNameForTemplate.find(".try"))
	{
		sNameForBAS2000Inf = sNameForTemplate.before(".try");
		sNameForBAS2000Inf += ".trf";
	}
#endif
  lFileSize = lFileGetSize(m_sFilename);

  pcBuffer = new char[LEN_MULT+1];
  pcBuffer[LEN_MULT] = '\0';
  nLen_Mult          = LEN_MULT;
  (void) dskbr(pnFile, pcBuffer, &nLen_Mult, &nStat);
  nEmpty();
  bIsDtrek = false;

  if (0 == nStat)
  {      
      sText = pcBuffer;
      if (0 == sText.find(ms_sStart)) {
          //  Read more of the header if there is some.
          bIsDtrek = true;
          

          while ((!nStat) && (!bTextHasEnd(sText)) /*&& (sText.length() <= MAX_HEADER_LENGTH)*/ ) {
              nLen_Mult          = LEN_MULT;
              memset(pcBuffer,0,nLen_Mult);
              (void) dskbr(pnFile, pcBuffer, &nLen_Mult, &nStat);

              // Change NULL characters into spaces.  
              for (i=0;i<nLen_Mult;i++) {
                  if (pcBuffer[i] == 0)
                      pcBuffer[i] = ' ';
              };
             
              sText += (sTemp = pcBuffer).remove('\r');
              if (0 != nStat)
              {
                  nStat = 0;
              }
          }
          if (!nStat)
            nStat = nParse(sText);
          if ((!nStat) && (!nGetValue((Cstring) "HEADER_BYTES", &nHeaderBytes))) {
              // Read HEADER_BYTES, and determine if we read all blocks.
              nHeaderBytes -= sText.length();
              nLen_Mult          = LEN_MULT;
              while ((nHeaderBytes>0) && (!nStat)) 
                {
                  (void) dskbr(pnFile, pcBuffer, &nLen_Mult, &nStat);
                  nHeaderBytes -= nLen_Mult;
                }
              nStat = 0;
          };

      }
#ifdef RAXIS
      else if (   (0 == sText.index("R-AXIS"))
          || (0 == sText.index("RAXIS")) )
      {
          // It is probably an R-AXIS header,
          // so read the R-AXIS header and create a d*TREK header from it.

          nStat = nReadRAXISheader(pnFile, sNameForTemplate, pcBuffer, this);
      }
#endif
#ifdef SIEMENS
      else if (0 == sText.index("FORMAT"))
      {
          // It is probably a Siemens header,
          // so read the Siemens header and create a d*TREK header from it.

          nStat = nReadSiemensHeader(pnFile, sNameForTemplate, pcBuffer, this);
      }
#endif
#ifdef MEDOPTICS
      else if (   ( (10 + 1024 * 1024 * 2) == lFileSize)
          || ( (10 +  512 *  512 * 2) == lFileSize)
          || ( 45785088               == lFileSize)
          )
      {

          nStat = nReadMedOpticsHeader(pnFile, m_sFilename, pcBuffer, this);
      }
#endif
#ifdef DIP2030
      else if ( (18001024 == lFileSize) || (32001024 == lFileSize) )
        {
          nStat = nReadDIP2030Header(pnFile, lFileSize, 
                                     sNameForTemplate, pcBuffer, this);
          if (!nStat) 
            nStat  = (-1L == fseek(pFdskfile(pnFile),lFileSize-1024,SEEK_SET));
          if (!nStat)
            nStat = nReadDIP2030Trailer(*pnFile, this);
          fseek(pFdskfile(pnFile),0,SEEK_SET);
          
      }
#endif
#ifdef MARIP
      else if (   (0 == strncmp("mar", (char *)(pcBuffer+207), 3))
          || (0 == strncmp("mar", (char *)(pcBuffer+128), 3))
          || (0 == strncmp("MAR", (char *)(pcBuffer+124), 3)) )
      {
          // It is a MAR IP header (we hope!)
          

          nStat = nReadMARIPHeader(pnFile, sNameForTemplate, pcBuffer, this);
      }
#endif
#ifdef NONIUS
      else if (    (0 == strncmp("No", (char *)(pcBuffer), 2))
		   && (0 == strncmp("ID =", (char *)(pcBuffer+80), 4) ) )
	{
	  // It is a Nonius KappaCCD header (we hope!)
          
	  nStat = nReadKappaCCDHeader(pnFile, sNameForTemplate, pcBuffer, this);
	}
#endif      
	  //else if (921654 == lFileSize)
	  else if (0 == strncmp("BM", (char *)pcBuffer, 2)) 
	{
          // A certain kind of Windows bitmap

          nStat = nReadWindowsBMPHeader(pnFile, sNameForTemplate, pcBuffer, this);
	}
#ifdef DTREK_CBF
      else if (0 == strncmp("###CBF:", (char *)pcBuffer, 7))
	{
	  nStat = nReadCBFHeader(pnFile, sNameForTemplate, pcBuffer, this);
	}
#endif
#ifdef BAS2000
      else if (bFileExists(sNameForBAS2000Inf))
      {
          nStat = nReadBAS2000Header(pnFile, sNameForTemplate, pcBuffer, this);
      }
#endif
#ifdef HIPIC
      else if (0 == strncmp("IM", (char *)pcBuffer, 2))
      {
          nStat = nReadHiPicHeader(pnFile, sNameForTemplate, pcBuffer, this);
      }
#endif
#ifdef MARCCD
      // WARNING!!! nReadMARCCDHeader MUST be the LAST one checked!
      else if (!nReadMARCCDHeader(pnFile, sNameForTemplate, pcBuffer, this))
      {
          nStat = 0;
      }
#endif
      else 
	nStat = 1;    
  } else
    nStat = 1;
      
  if ((nStat) && /*(lFileSize < MAX_HEADER_LENGTH) &&*/ (!bIsDtrek)) {
        int nExtFile = 2;
        int nExtLength;
        char* pcExtBuffer;
        
        
        nExtLength = m_sFilename.length();
        (void) dskbor(&nExtFile, m_sFilename.string(), &nExtLength, &nStat);
        nExtLength = lFileSize;
        if (!nStat) {
            
            pcExtBuffer = new char[nExtLength + 1];
            (void) dskbr(&nExtFile, pcExtBuffer, &nExtLength, &nStat);              
            (void) dskbcr(&nExtFile, &nLength);
            pcExtBuffer[nExtLength] = '\0';
            if (!nStat) {
                sText = pcExtBuffer;
                sText = sText.remove('\r');
                nStat = nParse(sText);
            };
            delete[] pcExtBuffer;
        };          
    };
    
    


  delete [] pcBuffer;

  if (0 == nStat) m_eThe_State = eImage_header_available_state;
  Cstring sDetName = "";
  (void) nGetValue(D_K_DetectorNames, &sDetName);
  if (   (0 == nGetValue(D_K_Rotation, &sTemp))
	 && ("" == sGetEnv("DTREK_IGNORE_ROTATION"))
	 && ("ADSC914_" != sDetName) )
    {
      // Do nothing, since header is probably OK
      // unless it is a SBC-CAT ADSC 914 detector
    }
#ifdef ADSCCCD
  // ADSC CCD images have a similar format to d*TREK images,
  // BUT they are missing lots of required d*TREK keywords in the
  // header, so edit the header here if we suspect it is an ADSC CCD image.
  // BUT be careful, do not edit the header if it has already been edited.

  else if (   (0 == nGetValue("OSC_RANGE",  &sTemp))
           && (0 == nGetValue("OSC_START",  &sTemp))  )
    //  else if (   ("ADSC914_" == nGetValue("PHIST",  &sTemp))
    {
      // Found keywords thought to be only in ADSC images,
      // so go generate missing keywords and info

      (void) nEditADSCCCDHeader(sNameForTemplate, this);
    }
#endif
#ifdef BRANDEISCCD
  else if (   (0 == nGetValue("PHIST",  &sTemp))
           && (0 == nGetValue("PHIINC",  &sTemp))  )
    {
      // Found keywords thought to be only in Brandeis B4 images,
      // so go generate missing keywords and info

      (void) nEditBrandeisCCDHeader(sNameForTemplate, this);
    }
#endif

  // If this header has an attached image,

  // Look for SCAN_TEMPLATE keyword and update the directory specification

  if ( (0 > nDim0) || (0 > nDim1) )
    {
      // If either of the dimensions is less than 0, get them from the header

      nDim0 = nDim1 = 0;

      (void) nGetValue(ms_sSize1, &nDim0);
      (void) nGetValue(ms_sSize2, &nDim1);
    }
  
  if ( (0 < nDim0) && (0 < nDim1) )
    {
      // This is an image, set scan_template based on the image filename,
      // NOT what the scan_template is in the header!

      // Try to decide if this is a legitimate image or a DTCOLLECT_HEADER

      if (0 == nGetValue(D_K_ScanTemplate, &sTemp) )
        {
          if ( lFileSize < (nDim0 * nDim1) )
            {
              // Cannot contain pixel values, so not a real image
	      // This could be dangerous with compressed images!
            }
          else
            {
              // See how many ?s there are in the template found in the
              // header (nPlaces).
              // If there are more than zero, feed that to sBuildScanTemplate,
              // otherwise default to 3.

              int nPlaces = 0, idx = -1;
              while ( (idx=sTemp.find('?', idx+1)) != -1)
                {
                  nPlaces++;
                }
              if (nPlaces <= 0) nPlaces = 3;
              sTemp = sBuildScanTemplate( sNameForTemplate, nPlaces);
              if ("" != sTemp)
                nReplaceValue(D_K_ScanTemplate, sTemp);
            }
        }
    }
  //+JWP 24-Mar-2014
  if ("" != sGetEnv("DTREK_UNBIN"))
    {
      // Try this to change things in the head if unbinning is asked for
      // note that nUnbinHeader() should not unbin if already unbinned
      // Only unbin if there is no image data (i.e. this is a true .head file)
      (void) nGetValue(ms_sSize1, &nDim0);
      (void) nGetValue(ms_sSize2, &nDim1);
      if ( (0 == nDim0) && (0 == nDim1) )
	nUnbinHeader();
    }
  //-JWP 24-Mar-2014
  return (nStat);
}

// Write the header to an already open file

int Cimage_header::nWrite(int* pnFile)
{
  int nLength;
  int nStat;
  Cstring sText;
  
  nReplaceValue("DTREK_DATE_TIME", sGetDate() + ' ' + sGetTime());
  nReplaceValue("DTREK_VERSION", (Cstring)D_K_DTREKVersion);
  nReplaceValue("DTREK_MODULE", sDtrekGetModuleName());

  //+JWP 24-Mar-2014
  if ("" != sGetEnv("DTREK_UNBIN"))
    {
      // Try this to change things in the head if unbinning is asked for
      // nUnbinHeader() should not unbin if already unbinned
      nUnbinHeader();
      
    }
  //-JWP 24-Mar-2014

  if (nParse(sText,true))
    return 1;
  

  nLength = sText.length();
  (void) dskbw(pnFile, sText.string(), &nLength, &nStat);
  return (nStat);
}

int Cimage_header::nWrite(const Cstring& sFilename)
{
  // Write header to a file all by itself

  int nStat;   // General status
  int nFile;   // File number for dsk* routines
  int nBytes;  // Size in chars of filename
  int nSize;   // Size in bytes of file
  int nDim0 = 0;
  int nDim1 = 0;

  // Get the SIZE1, SIZE2 parameters

  (void) nGetValue(ms_sSize1, &nDim0);
  (void) nGetValue(ms_sSize2, &nDim1);

  (void) nReplaceValue(ms_sSize1, (int) 0);
  (void) nReplaceValue(ms_sSize2, (int) 0);

  // If the overwrite environment variable is set, then make sure
  // any existing file gets a version number appended, so it is
  // not overwritten

  nStat = nFileAppendVersion(sFilename, TRUE);
  if (0 != nStat)
    {
      // Headers are so important, this is a FATAL error

      cout << "FATAL ERROR renaming header file: " << sFilename << "!\n";
#ifdef SSI_PC
          return nStat;
#else
      exit (nStat);
#endif
    }

  nFile  = 1;
  nSize  = 0;
  nBytes = sFilename.length();

  (void) dskbow(&nFile, sFilename.string(), &nBytes, &nSize, &nStat);
  if (0 == nStat)
    {
      nStat = nWrite(&nFile);
    }
  (void) dskbcw(&nFile, &nBytes);
  if ( (0 == nStat) && (0 != nBytes) ) 
      nStat = nBytes;
  (void) nReplaceValue(ms_sSize1, nDim0);
  (void) nReplaceValue(ms_sSize2, nDim1);
  return (nStat);
}

int Cimage_header::nList(void)
{
    Cstring sText;
    
    
    if (!nParse(sText,true))
    {
        cout << "     Header listing:: \n" << sText.through(ms_sEnd) << "\n:: end of header listing.\n";
        return 0;
    }
    else
        cout << "     Header cannot be listed: it has no characters!\n";
    return (1);
}

int
Cimage_header::nGetImageSize(void)
{
  // Return total size of image include header size and binary data

  int nDim0, nDim1, nBytesPerPixel;

  nDim0 = nDim1 = 0;

  (void) nGetValue(ms_sSize1, &nDim0);
  (void) nGetValue(ms_sSize2, &nDim1);
  nBytesPerPixel = 2;  // Assume 2-byte pixels for now
  return (nLength() + (nBytesPerPixel * nDim0 * nDim1));
}

bool
Cimage_header::bKeywordExists(const Cstring &sKeyword)
{
   std::map<Cstring,Cstring,std::less<Cstring> >::iterator oIt1;
   oIt1 = m_oKeywords.find(sKeyword);
   if (oIt1 != m_oKeywords.end()) 
       return TRUE;
   else
       return FALSE;
}

bool
Cimage_header::bIsValueValid(const Cstring &sValue)
{
   char *pc,*pcEnd;

   pc    = sValue.string();
   pcEnd = pc + sValue.length();
   while (   ('\0' != *pc)
          && ('{'  != *pc)
          && ('}'  != *pc)
          && (';'  != *pc) )
     pc++;
   if (pc != pcEnd)
     return FALSE;

   return TRUE;
}

bool
Cimage_header::bIsKeywordValid(const Cstring &sKeyword)
{
   char *pc,*pcEnd;

   if (32 < sKeyword.length() || 1 > sKeyword.length())
     return FALSE;

   pc    = sKeyword.string();
   pcEnd = pc + sKeyword.length();
   while (   ('a' <= *pc && 'z' >= *pc)
          || ('A' <= *pc && 'Z' >= *pc)
          || ('_' ==   *pc)
          || ('0' <= *pc && '9' >= *pc) )
      pc++;
   if (pc != pcEnd)
      return FALSE;

   return TRUE;

}

/****************************************************************************
 * This routine sorts an array of Cstring pointers in lexigraphical order   *
 * using the Heapsort algorithm.                                            *
 * Taken from: Numerical Recipes in C, p.247                                *
 ****************************************************************************/
void Cimage_header::vHeapSort(int       n,   /* the number of objects to sort */
                              Cstring **pps) /* array[0..n-1] of ptrs to sort */
{
   int i,j,l,ir;
   Cstring *rra;
   Cstring **ra;

   ra = pps-1;  // Kludge from NR for C arrays to start indexing at 1

   l = (n >> 1)+1;
   ir = n;

   if (n==1)
       return;
   for(;;){
      if(l > 1)
         rra = ra[--l];
      else{
         rra = ra[ir];
         ra[ir] = ra[1];
         if(--ir == 1){
            ra[1] = rra;
            return;
         }
      }
      i = l;
      j = l << 1;
      while(j <= ir){
         if(j < ir && *ra[j] < *ra[j+1])
            j++;
         if(*rra < *ra[j]){
            ra[i] = ra[j];
            j += (i=j);
         }
         else
            j = ir+1;
      }
      ra[i] = rra;
   }

}

void
Cimage_header::vSetVerbose(const int nVerbose)
{
  ms_nVerbose = nVerbose;
}

int
Cimage_header::nGetVerbose(void)
{
  return (ms_nVerbose);
}

eImage_data_types
Cimage_header::eGetImageDataType(void)
{
  eImage_data_types eImageDataType =   eImage_other_type;
  Cstring sTemp;

  if (0 == nGetValue(ms_sDataType, &sTemp))
    {
      if (D_K_ShortInt == sTemp)
        {
          eImageDataType = eImage_I2;
        }
      else if (D_K_UnsignedShortInt == sTemp)
        {
          eImageDataType  = eImage_uI2;
        }
      else if (D_K_LongInt == sTemp)
        {
          eImageDataType = eImage_I4;
        }
      else if (D_K_UnsignedLongInt == sTemp)
        {
          eImageDataType = eImage_uI4;
        }
      else if (D_K_SignedChar == sTemp)
        {
          eImageDataType = eImage_byte;
        }
      else if (D_K_UnsignedChar == sTemp)
        {
          eImageDataType = eImage_ubyte;
        }
      else if ( (D_K_FloatIEEE == sTemp) || ("float" == sTemp) )
        {
          eImageDataType = eImage_realIEEE;
        }
      else if (D_K_Compressed == sTemp)
        {
          eImageDataType = eImage_compressed;
        }
      else
        {
          eImageDataType = eImage_other_type;
        }
    }
  return (eImageDataType);
}


                    
int Cimage_header::nParse(Cstring& sText,bool bBuild) {
    int nStat;
    int nPass;
    int nLength;
    Cstring sTemp;
    Cstring sData;
    Cstring sValue;
    Cstring sKeyword;

    
    nStat = 0;
    
    if (bBuild)
      {
        sText = ms_sStart;
        nStat = 0;
        std::map<Cstring,Cstring,std::less<Cstring> >::iterator oMapi;
        for (nPass = 0; nPass < 2; nPass++)
          {
            if (nPass == 0)
              nLength = 0;
            else
              {
                // Make sure you use 5 digits and nothing more or less
                //cout << "AnLength is: " << nLength << endl << flush;
                Cstring sSize;

                // The 7 in the next line is for the #####;^n added to sText

                nLength = ((nLength + ms_sStart.length() + 7 + ms_sEnd.length())
                           / LEN_MULT) * LEN_MULT + LEN_MULT;
                //cout << "BnLength is: " << nLength << endl << flush;
                sSize = Cstring(nLength, 5);      
                //cout << "sSize >>>" << sSize << "<<<\n" << flush;
                sText += sSize;
                sText += ";\n";
              }
            for (oMapi = m_oKeywords.begin(); oMapi != m_oKeywords.end(); ++oMapi)
              {
                if (nPass == 0)
                  {
                    nLength += 3+(*oMapi).first.length() + (*oMapi).second.length();
                  }
                else if (!(*oMapi).first.contains("HEADER_BYTES"))
                  {
                    sText += (*oMapi).first;
                    sText += "=";
                    sText += (*oMapi).second;
                    sText += ";\n";
                  }
              }
          }
        sText += ms_sEnd;
        sText += replicate(' ',nLength - sText.length());
        if (0 != (sText.length() % LEN_MULT))
          {
            cout << "\n\n**************************************";
            cout << "\nWARNING!! image header length problem!"
                 << "\nPlease contact Rigaku/MSC, " 
                 << sText.length() % LEN_MULT;
            cout << "\n**************************************\n\n" << flush;
          }
      } 
    else
      {
        int nFirstChar,nLastChar;
        nEmpty();
        nStat = 0;
        for (nFirstChar = 0; nFirstChar < sText.length();)
          {
            nLastChar = nFirstChar;
            while ((nLastChar< sText.length()) && (sText.GetAt(nLastChar)!=';'))
              nLastChar++;
            sTemp = sText.substr(nFirstChar,nLastChar-nFirstChar);

            // Parse sTemp containing the string itself.
            if (sTemp.contains('='))
              {
                sKeyword = sTemp.before('=').remove(' ').remove('\t');
                sData = sTemp.after('=');
		//+JWP 2008-07-25
		// Remove leading spaces from the value
		while (' ' == sData.GetAt(0))
		  sData = sData.after(0);
		//-JWP 2008-07-25
                if (sKeyword.contains((Cstring) "HEADER_BYTES")) 
                  sKeyword = "HEADER_BYTES";
                    
                if (sKeyword.length())
                  {
                    m_oKeywords.erase(sKeyword);
                    m_oKeywords.insert(std::pair<Cstring,Cstring>(sKeyword,sData));
                  }
              } 
            else if (nLastChar!=sText.length())
              {
                nStat = 1;
                break;
              }
            nFirstChar = nLastChar;
            while (   (nFirstChar< sText.length()) 
                   && (sText.GetAt(nFirstChar)!='\n')
                   && (sText.GetAt(nFirstChar)!='\r'))
              nFirstChar++;
            while (   (nFirstChar< sText.length()) 
                   && (((sText.GetAt(nFirstChar)=='\n'))
                   || (sText.GetAt(nFirstChar)=='\r')))
              nFirstChar++;
          }

        if (0 == nStat) 
          m_eThe_State = eImage_header_available_state;
        else
          m_eThe_State = eImage_header_unknown_state;
      }
    return nStat;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cimage_header::nReplaceValueDictionary(Cstring sDictionary,Cstring sName,double fValue,int nIndex,Cstring sPrefix)
{
    Cstring         sValue(fValue);
    
    return nReplaceValueDictionary(sDictionary,sName,sValue,nIndex,sPrefix);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cimage_header::nReplaceValueDictionary(Cstring sDictionary,Cstring sName,int nValue,int nIndex,Cstring sPrefix)
{
    Cstring sValue(nValue);
    
    return nReplaceValueDictionary(sDictionary,sName,sValue,nIndex,sPrefix);
}

int Cimage_header::nReplaceValueDictionary(Cstring sDictionary,
                                           Cstring sName,
                                           Cstring sValue,
                                           int nIndex,
                                           Cstring sPrefix)
{
    Cstring sDictionaryNames;
    Cstring sNameArray;
    Cstring sTemp;
    int nx;

    if (sName == "NAMES")
        return 1;

    sDictionaryNames = sDictionary;
    sDictionaryNames += "_NAMES";

    nx = sName.find('*');
    if (nx==0)
        sName.replace("*",nx,sPrefix);
    
    nx = sName.find('*');
    if( nx > 0 ) 
    {
        sTemp = nIndex;
        sName.replace("*", nx, sTemp);
    }
    
    sTemp = sDictionary;
    sTemp += "_";
    sTemp += sName;
    sName = sTemp;

    if (nGetValue(sDictionaryNames,&sNameArray)) 
        sNameArray = "";
    
    nx = sNameArray.find(sName);
    
    if (nx < 0)
    {
        if (sNameArray.length())
            sNameArray += " ";
        sNameArray += sName;
        nReplaceValue(sDictionaryNames,sNameArray);
    }
    
    return nReplaceValue(sName,sValue);
}

int Cimage_header::nGetValue(const Cstring& sKeyword, std::vector<int>& anVec, char* pcListSeparators)
{
    Cstring     sTemp("");
    
    if( 0 != nGetValue(sKeyword, &sTemp) )
        return 1;
    
    if( 0 != sTemp.nListToVector(anVec, pcListSeparators) )
        return 1;

    return 0;
}

int Cimage_header::nGetValue(const Cstring& sKeyword, std::vector<double>& afVec, char* pcListSeparators)
{
    Cstring     sTemp("");
    
    if( 0 != nGetValue(sKeyword, &sTemp) )
        return 1;
    
    if( 0 != sTemp.nListToVector(afVec, pcListSeparators) )
        return 1;

    return 0;
}


int Cimage_header::nGetValueDictionary(Cstring sDictionary,Cstring sName,double& fValue,int nIndex,Cstring sPrefix) {
    Cstring sValue;
    if (nGetValueDictionary(sDictionary,sName,sValue,nIndex,sPrefix))
        return 1;
    if (1 == sscanf(sValue.string(),"%lf",&fValue))
        return 0;
    
    return 1;
};

int Cimage_header::nGetValueDictionary(Cstring sDictionary,Cstring sName,int& nValue,int nIndex,Cstring sPrefix) {
    Cstring sValue;
    if (nGetValueDictionary(sDictionary,sName,sValue,nIndex,sPrefix))
        return 1;
    if (1 == sscanf(sValue.string(),"%d",&nValue))
        return 0;
    return 1;
};

int Cimage_header::nGetValueDictionary(Cstring sDictionary,Cstring sName,Cstring& sValue,int nIndex,Cstring sPrefix) {
    Cstring sNameArray;
    Cstring sDictionaryNames;
    int nx;
    Cstring sTemp;

    sDictionaryNames = sDictionary;
    sDictionaryNames += "_NAMES";
    if (nGetValue(sDictionaryNames,&sNameArray)) 
        return 1;

    nx = sName.find('*');
    if (nx==0)
        sName.replace("*",nx,sPrefix);
    nx = sName.find('*');
    if (nx>0) {
        sTemp = nIndex;
        sName.replace("*",nx,sTemp);
    };
    
    sTemp = sDictionary;
    sTemp += "_";
    sTemp += sName;
    sName = sTemp;

    nx = sNameArray.find(sName);
    if (nx < 0)
        return 1;
    return nGetValue(sName,&sValue);
};


int Cimage_header::nGetDictionaryNames(Cstring sDictionary, std::vector<Cstring>& asNames)
{
    Cstring         sDictionaryNames = sDictionary;
    sDictionaryNames += "_NAMES";
    
    Cstring         sNameArray("");
    if( nGetValue(sDictionaryNames, &sNameArray) ) 
        return 1;
    
    asNames.clear();
    Cstring         sTemp("");
    while( sNameArray.length() ) 
    {
        sTemp = sNameArray.before(' ');
        
        if( sTemp.length() )
            asNames.push_back(sTemp);
        
        sNameArray = sNameArray.after(' ');
    }

    return 0;
}

int Cimage_header::nClearDictionary(Cstring sDictionary, Cstring sPattern)
{
    std::vector<Cstring>        asNames;
    Cstring                     sDictionaryNames("");
    Cstring                     sNameArray("");
    Stringreg                   oReg;
    int                         nx = 0;

    sDictionaryNames = sDictionary;
    sDictionaryNames += "_NAMES";
    if (nGetDictionaryNames(sDictionary,asNames))
      return 1;
    sPattern = (sDictionary + sPattern);
    for (nx = 0; nx < asNames.size(); nx++)
      {
        if ((sPattern.length() == 0) || (!oReg.nParse(asNames[nx],sPattern)))
          nDelete(asNames[nx]);
        else
          {
            sNameArray += asNames[nx];
            sNameArray += " ";
          }
      }
    nReplaceValue(sDictionaryNames,sNameArray);
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// Image format type stuff
typedef struct
{
    eImage_format_types eType;
    char*               sName;
}IMAGE_FORMAT_TYPE;

static const IMAGE_FORMAT_TYPE s_c_stImageFormatTypes[]=
{
    enImage_format_unknown,         "unknown",                          
    enImage_format_RIGAKU_CCD,      "RigakuCCD",
    enImage_format_RAXIS,           "RAXIS",
    enImage_format_DIP,             "DIP",
    enImage_format_MEDOPTICS,       "Medoptics",
    enImage_format_BRUKER_SIEMENS,  "Bruker_Siemens",
    enImage_format_WINBMP,          "Windows_Bitmap",
    enImage_format_ADSC,            "ADSC",
    enImage_format_BRANDEIS,        "Brandeis",
    enImage_format_CBFDTREK,        "CBF_DTREK",
    enImage_format_BAS2000,         "BAS2000",
    enImage_format_HIPIC,           "HiPic",

    enImage_format_unknown,         NULL                // by design NULL must be the last string value in this array           
};

eImage_format_types Cimage_header::enGetFormatTypeByName(Cstring& sFormat)
{
    int     ii = 0;
    while( NULL != s_c_stImageFormatTypes[ii].sName )
    {
        if( 0 == strcmp(sFormat.string(), s_c_stImageFormatTypes[ii].sName) )
            break;
        ii++;
    }
    
    return s_c_stImageFormatTypes[ii].eType;
}

Cstring Cimage_header::sGetFormatTypeName(eImage_format_types enType)
{
    int     ii = 0;
    while( NULL != s_c_stImageFormatTypes[ii].sName )
    {
        if( enType == s_c_stImageFormatTypes[ii].eType )
            break;
        ii++;
    }
    
    return Cstring(s_c_stImageFormatTypes[ii].sName);
}


eImage_format_types Cimage_header::enGetOriginalFormat()
{
    Cstring     sFormat("");
    
    if( 0 != nGetValue(Cstring(D_K_OriginalImageFormat), &sFormat) )
        return enImage_format_RIGAKU_CCD;  // If the keyword is absent, must assume the original format is d*TREK
    
    return enGetFormatTypeByName(sFormat);
}
int 
Cimage_header::nUnbinHeader(void)
{
  // Check if a header has already had changes applied for the Cimage::nUnbin() method,
  // and if not, then apply the changes.
  // 
  // This method was created so that Cimage_header operations were moved out of Cimage class.
  // This method is called by the Cimage::nUnbin(void) method.
  // The data in any Cimage object is not changed in this method,
  // but should be changed in the Cimage::nUnbin(() method.

  // This method does NOT know about Cimage::m_nDim[0], m_nDim[1];

  // Probably should do this only if (eImage_uI2 != m_eData_type)

  // Unbin the image data and adjust the image header if needed.
  //cout << "DTREK_UNBIN is set B\n";

  int nStat = 0;
  float a4fValues[4];
  Cstring  sDetPrefix = "";

  nStat = nGetValue("DTREK_UNBINNED", 2, &a4fValues[0]);
  if (0 == nStat)
    {
      // Header ALREADY has flag that it was UNBINNED, 
      // so DO NOT unbin again.
      
      //std::cout << "WARNING: In Cimage_header.nUnbinHeader(), already UNBINNED\n\n";
      return (-1);
    }
  // Otherwise, no hint that this header object was unbinned,
  // so apply the changes.

  // Probably should do this only if (eImage_uI2 != m_eData_type)

  // Modify this header
  // Make sure SIZE1, SIZE2
  // *_SPATIAL_DISTORTION_INFO  (pixel size, beam center)

  nGetValue(Cstring(D_K_DetectorNames), 1, &sDetPrefix);
  //  cout << "DetPrefix: " << sDetPrefix << endl;

  nGetValue(sDetPrefix + Cstring(D_K_SpatialDistortionInfo), 4, a4fValues);

  /*  cout << "SpatInf: " 
       << a4fValues[0] << ", "
       << a4fValues[1] << ", "
       << a4fValues[2] << ", "
       << a4fValues[3] << endl;
  */

  a4fValues[0] *= 2.0;
  a4fValues[1] *= 2.0;
  a4fValues[2] *= 0.5;
  a4fValues[3] *= 0.5;

  nReplaceValue(sDetPrefix + Cstring(D_K_SpatialDistortionInfo), 4, a4fValues);
  /*
  cout << "SpatInf: " 
       << a4fValues[0] << ", "
       << a4fValues[1] << ", "
       << a4fValues[2] << ", "
       << a4fValues[3] << endl;
  */

  Cstring sTemp = "";
  nStat = nGetValue(sDetPrefix + Cstring(D_K_SpatialBeamPosition), &sTemp);
  if (0 == nStat)
    {
      nReplaceValue(sDetPrefix + Cstring(D_K_SpatialBeamPosition), 2, a4fValues);      
    }

  nStat = nGetValue(sDetPrefix + Cstring(D_K_DetectorDimensions), 2, a4fValues);
  if (0 != nStat)
    {
      cout << "WARNING: No detector dimensions!\n" << endl;
    }
  a4fValues[0] *= 2.0;
  a4fValues[1] *= 2.0;
  nReplaceValue(sDetPrefix + Cstring(D_K_DetectorDimensions), 2, a4fValues);

  nStat = nGetValue(sDetPrefix + Cstring(D_K_UnbinnedDimensions), 2, a4fValues);
  if (0 != nStat)
    {
      if (2 < ms_nVerbose)
	cout << "WARNING: No detector unbinned dimensions!\n" << endl;
      nStat = nGetValue(Cstring(D_K_UnbinnedDimensions), 2, a4fValues);
      if ( (0 == nStat) && (0.0 >= a4fValues[0]) && (0.0 >= a4fValues[1]))
        nDelete(Cstring(D_K_UnbinnedDimensions));
    }
  else
    {
      a4fValues[0] *= 2.0;
      a4fValues[1] *= 2.0;
      nReplaceValue(sDetPrefix + Cstring(D_K_UnbinnedDimensions), 2, a4fValues);
    }

  // Add keyword to header of this image, so we know it is unbinned already

    a4fValues[0] = 0.5;
    a4fValues[1] = 0.5;

  nReplaceValue("DTREK_UNBINNED", 2, &a4fValues[0]);
  nStat = 0;
  return (nStat);
}

////////////////////////////////////////////////////////////////////////////////////
