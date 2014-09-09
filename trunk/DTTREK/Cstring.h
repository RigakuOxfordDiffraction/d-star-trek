#ifndef MSC_MSCSTRING_H
#define MSC_MSCSTRING_H

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

/****************************************************************************
 *                               Header Files                               *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

//#if !defined(VC8) && !defined(__APPLE__)
#if defined(VC6)
    #include <iostream.h>
    #include <iomanip.h>
	#include <fstream.h>
#else 
    #include <iostream>
    #include <iomanip>
	#include <fstream>
	
        // don't do this in a header file! see: http://www.parashift.com/c++-faq-lite/coding-standards.html#faq-27.5
        // using namespace std;
#endif

#include <vector>

#include "dtrekdefs.h"

/****************************************************************************
 *                           Global Variables                               *
 ****************************************************************************/
const char NULL_CHARACTER = '\0';
const int  gs_Cstring_Free_Bins = 14;               // 8==> Largest pointer bin at 1 << (8-1) == 128.  
const int  gs_Cstring_Max_Free = 512;
const int  gs_Cstring_Minimal_Alloc_Log2 = 3;       // 3==> Smallest alloc size of 8.
const int  gs_Cstring_Minimal_Of_Each_Type = 10;    // Minimal number of pointers in each bin.

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 ** Class definition for string class.  Used MSCString as opposed to       **
 ** Cstring as Cstring used by others.                                     **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/
class DTREK_EXPORT MSCString {

public:

// constructors

   MSCString(){vCreate();}
   MSCString(const MSCString &s, const int nPosition = 0, int nLength = -1)
      {vCreate(); assign(s,nPosition,nLength);}
   MSCString(const char c, const int nRepetitions = 1)
      {vCreate(); assign(c,nRepetitions);}
   MSCString(const char *pc, int nLength = -1)
      {vCreate(); assign(pc,nLength);}
   MSCString(const short i, int nFieldWidth = 0)
      {vCreate(); assign((long)i,nFieldWidth);}
   MSCString(const int n, int nFieldWidth = 0)
      {vCreate(); assign((long)n,nFieldWidth);}
   MSCString(const long l, int nFieldWidth = 0)
      {vCreate(); assign(l,nFieldWidth);}
   MSCString(const unsigned long ul, int nFieldWidth = 0)
      {vCreate(); assign(ul,nFieldWidth);}
   MSCString(const float f, int nFieldWidth = 0, const int nDecimal = 4)
      {vCreate(); assign((double)f,nFieldWidth,nDecimal);}
   MSCString(const double d, int nFieldWidth = 0, const int nDecimal = 4)
      {vCreate(); assign(d,nFieldWidth,nDecimal);}

// destructor

   ~MSCString();
//   ~MSCString() {vDestroy();}

// various operators (unary)

   MSCString &operator=(const MSCString &s){return (assign(s));        }
   MSCString &operator=(const char c)      {return (assign(c));        }
   MSCString &operator=(const char *pc)    {return (assign(pc));       }
   MSCString &operator=(const short i)     {return (assign((long)i));  }
   MSCString &operator=(const int n)       {return (assign((long)n));  }
   MSCString &operator=(const long l)      {return (assign(l));        }
   MSCString &operator=(const float f)     {return (assign((double)f));}
   MSCString &operator=(const double d)    {return (assign(d));        }

   MSCString &operator+=(const MSCString &s){return (append(s));        }
   MSCString &operator+=(const char c)      {return (append(c));        }
   MSCString &operator+=(const char *pc)    {return (append(pc));       }
   MSCString &operator+=(const short i)     {return (append((long)i));  }
   MSCString &operator+=(const int n)       {return (append((long)n));  }
   MSCString &operator+=(const long l)      {return (append(l));        }
   MSCString &operator+=(const float f)     {return (append((double)f));}
   MSCString &operator+=(const double d)    {return (append(d));        }

   
   char operator[](size_t n) const
      {return ((n < m_Length) ? m_pcString[n] : NULL_CHARACTER);}
   char operator[](int n) const
      {return (((size_t)n < m_Length) ? m_pcString[n] : NULL_CHARACTER);}
   char &operator[](size_t n){return m_pcString[n];}
   

// member functions

   MSCString &append(const MSCString &s, const int nPosition = 0,
                   int nLength = -1);
   MSCString &append(const char c, const int nRepetitions = 1);
   MSCString &append(const char *pc, int nLength = -1);
   MSCString &append(const short i, int nFieldWidth = 0)
      {return (append((long)i,nFieldWidth));}
   MSCString &append(const int n, int nFieldWidth = 0)
      {return (append((long)n,nFieldWidth));}
   MSCString &append(const long l, int nFieldWidth = 0);
   MSCString &append(const float f, int nFieldWidth = 0, const int nDecimal = 4)
      {return (append((double)f,nFieldWidth,nDecimal));}
   MSCString &append(const double d, int nFieldWidth = 0,
                     const int nDecimal = 4);

   MSCString &assign(const MSCString &s, const int nPosition = 0,
                   int nLength = -1);
   MSCString &assign(const char c, const int nRepetitions = 1);
   MSCString &assign(const char *pc, int nLength = -1);
   MSCString &assign(const short i, int nFieldWidth = 0)
      {return (assign((long)i,nFieldWidth));}
   MSCString &assign(const int n, int nFieldWidth = 0)
      {return (assign((long)n,nFieldWidth));}
   MSCString &assign(const long l, int nFieldWidth = 0);
   MSCString &assign(const unsigned long ul, int nFieldWidth = 0);
   MSCString &assign(const float f, int nFieldWidth = 0, const int nDecimal = 4)
      {return (assign((double)f,nFieldWidth,nDecimal));}
   MSCString &assign(const double d, int nFieldWidth = 0,
                     const int nDecimal = 4);


   int nExpand(const int nRequestedSize = 0);

   bool bParseFloat(float& fFloat);
   bool bMakeFloat(float fFloat,float fRoundError= 0.01,int nMaxNumer=10,int nMaxDenom=10);

   int compare(const MSCString &s) const;
   int compare(const char *pc) const;

   int compare_no_case(const MSCString &s) const;
   int compare_no_case(const char *pc) const;

   int find(const MSCString &s, const int nPosition = 0) const;
   int find(const char *pc, const int nPosition = 0, int nLength = -1) const;
   int find(const char c, const int nPosition = 0) const;
   
   int nRead(FILE* pFIn);

   int nReverseFind(const MSCString &s, int nPosition = -1) const;
   int nReverseFind(const char *pc, int nPosition = -1, int nLength = -1) const;
   int nReverseFind(const char c, int nPosition = -1) const;
   bool bRegularMatch(const MSCString& sPattern) const;

    // Useful methods for filenames
    MSCString sGetFilePath() const;
    MSCString sGetFileExt() const;
    MSCString sGetFileName() const;
    MSCString sGetFileTitle() const;



// member functions to be added later.
//   MSCString &insert(const int nPosition1, const MSCString &s,
//                     const int nPosition2 = 0, int nLength = -1);
//   MSCString &insert(const int nPosition1, const char *pc, int nLength - -1);
//   MSCString &insert(const int nPosition1, const char c, int nRepetitions = 1);
//
   MSCString &remove(const int nPosition = 0, int nLength = -1);
   MSCString &remove(const char c);
   MSCString &replace(const int nPosition,const int nLength,MSCString sNew);
   MSCString &replace(MSCString sToFind,const int nStart,MSCString sNew);

// These functions return the value of member variables which are private.
// See also Microsoft-isms (below) for cast to char *.

   size_t reserved(void) const {return m_ReservedLength;}
   size_t length(void) const   {return m_Length;}
   char *string(void) const    {return m_pcString;}

// substring function, used by functions below

   inline MSCString substr(const int nPosition, int nLength = -1)
      {return MSCString((*this),nPosition,nLength);}

// Return string with all characters to the left of the argument

   MSCString Left(const int n)        {return (substr(0,n))     ;  }
   MSCString left(const int n)        {return (substr(0,n))     ;  }
   MSCString left(const MSCString &s) {return (substr(0,find(s))); }
   MSCString left(const char c)       {return (substr(0,find(c))); }
   MSCString left(const char *pc)     {return (substr(0,find(pc)));}

// Return string with all characters to the right of the argument

   MSCString right(const MSCString &s);
   MSCString right(const char c);
   MSCString right(const char *pcText);
   MSCString Right(int n);
   MSCString right(int n) {return Right(n);}

// Return string with all characters from beginning up through the argument

   MSCString through(const MSCString &s);
   MSCString through(const char c);
   MSCString through(const char *pcText);

// Return string with all characters from argument to the end

   MSCString from(const int nFrom)    {return (substr(nFrom));       }
   MSCString from(const MSCString &s) {return (substr(find(s)));     }
   MSCString from(const char c)       {return (substr(find(c)));     }
   MSCString from(const char *pcText) {return (substr(find(pcText)));}

// gnu-isms

   MSCString after(const int n)       {return (substr(n+1));}
   MSCString after(const MSCString &s){return (right(s));   }
   MSCString after(const char *pc)    {return (right(pc));  }
   MSCString after(const char c)      {return (right(c));   }

   MSCString at(int nPosition, int nLength){return substr(nPosition,nLength);}
   MSCString before(const int  n)      {return (left(n)); }
   MSCString before(const MSCString &s){return (left(s)); }
   MSCString before(const char *pc)    {return (left(pc));}
   MSCString before(const char c)      {return (left(c)); }

   bool contains(const MSCString &s) const {return (-1 != find(s)); }
   bool contains(const char *pc) const     {return (-1 != find(pc));}
   bool contains(const char c) const       {return (-1 != find(c)); }
   void downcase(void);
   void upcase(void);
   int index(const MSCString &s) const {return find(s); }
   int index(const char *pc) const     {return find(pc);}
   int index(const char c) const       {return find(c); }

// Microsoft-isms

   
   inline operator const char *() const {return m_pcString;}
   

   size_t GetLength(void) const         {return m_Length;  }
   char GetAt(int n = 0) const
      {return ((0 <= n && (size_t)n < m_Length) ? m_pcString[n] : NULL_CHARACTER);}
   void SetAt(int nPosition = 0, char c = NULL_CHARACTER)
      {if(0 <= nPosition && (size_t)nPosition < m_Length) m_pcString[nPosition] = c; return;}
#ifdef WIN32
    void Format( char* sFormat, ... );
#endif


private:

// member variables

   size_t m_ReservedLength; // allocated length is 1 more (for the \0 char)
   size_t m_Length;
   char *m_pcString;

   static int    ms_anFreeBinsStart[gs_Cstring_Free_Bins];
   static int    ms_anFreeBinsCount[gs_Cstring_Free_Bins];
   static char*  ms_apcFreeBins[gs_Cstring_Max_Free];
   static int    ms_nNumThis;
   static int    ms_nFreeCount;
   

// member functions
   void vFreeDelete(void);
   int nGetBin(int nSize);
   char* pcFreeFromBin(int nBin);
   int nInsertInBin(char* pcPointer,int nBin);


   void vCreate(void);
   void vDestroy(void);

public:   
   int nListToVector(std::vector<double>& vec, char cListSeparator);

   int nListToVector(std::vector<int>&          vecInts, char* pcListSeparators);
   int nListToVector(std::vector<double>&    vecDoubles, char* pcListSeparators);
   int nListToVector(std::vector<MSCString>& vecStrings, char* pcListSeparators);

   int nListToVectorOfPairs(std::vector< std::pair<int,int> >& vecPairsOfInts, char* pcListSeparators, char* pcPairSeparators);

   bool bIsNumeric();

   void vAbbreviateFilePath(const int cnMaxNumberOfFilePathChars);
};
/****************************************************************************
 ****************************************************************************
 **                                                                        **
 ** End of MSCString class definition.                                     **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 *                    Overloaded operators (binary)                         *
 ****************************************************************************/
inline bool operator==(const MSCString &sLHS, const MSCString &sRHS)
   {return (0 == sLHS.compare(sRHS)); }
inline bool operator==(const MSCString &sLHS, const char *pcRHS)
   {return (0 == sLHS.compare(pcRHS));}
inline bool operator==(const char *pcLHS, const MSCString &sRHS)
   {return (0 == sRHS.compare(pcLHS));}

inline bool operator!=(const MSCString &sLHS, const MSCString &sRHS)
   {return !(sLHS == sRHS); }
inline bool operator!=(const MSCString &sLHS, const char *pcRHS)
   {return !(sLHS == pcRHS);}
inline bool operator!=(const char *pcLHS, const MSCString &sRHS)
   {return !(pcLHS == sRHS);}

inline bool operator>=(const MSCString &sLHS, const MSCString &sRHS)
   {return ( 0 <= sLHS.compare(sRHS)); }
inline bool operator>=(const MSCString &sLHS, const char *pcRHS)
   {return ( 0 <= sLHS.compare(pcRHS));}
inline bool operator>=(const char *pcLHS, const MSCString &sRHS)
   {return ( 0 >= sRHS.compare(pcLHS));}

inline bool operator<=(const MSCString &sLHS, const MSCString &sRHS)
   {return ( 0 >= sLHS.compare(sRHS)); }
inline bool operator<=(const MSCString &sLHS, const char *pcRHS)
   {return ( 0 >= sLHS.compare(pcRHS));}
inline bool operator<=(const char *pcLHS, const MSCString &sRHS)
   {return ( 0 <= sRHS.compare(pcLHS));}

inline bool operator>(const MSCString &sLHS, const MSCString &sRHS)
   {return ( 0 < sLHS.compare(sRHS)); }
inline bool operator>(const MSCString &sLHS, const char *pcRHS)
   {return ( 0 < sLHS.compare(pcRHS));}
inline bool operator>(const char *pcLHS, const MSCString &sRHS)
   {return ( 0 > sRHS.compare(pcLHS));}

inline bool operator<(const MSCString &sLHS, const MSCString &sRHS)
   {return ( 0 > sLHS.compare(sRHS)); }
inline bool operator<(const MSCString &sLHS, const char *pcRHS)
   {return ( 0 > sLHS.compare(pcRHS));}
inline bool operator<(const char *pcLHS, const MSCString &sRHS)
   {return ( 0 < sRHS.compare(pcLHS));}

inline MSCString operator+(const MSCString& sLHS, const char cRHS)
   {MSCString sTemP = sLHS; return (sTemP += cRHS); }
inline MSCString operator+(const MSCString& sLHS, const MSCString& sRHS)
   { MSCString sTemP = sLHS; return (sTemP += sRHS); }
//   {return (MSCString(sLHS) += sRHS); }
inline MSCString operator+(const MSCString& sLHS, const char *pcRHS)
   {MSCString sTemP = sLHS; return (sTemP += pcRHS);}
//   {return (MSCString(sLHS) += pcRHS);}
inline MSCString operator+(const char *pcLHS, const MSCString& sRHS)
   {MSCString sTemP = pcLHS; return (sTemP += sRHS);}
//   {return (MSCString(pcLHS) += sRHS);}

#if defined(VC6)
istream &operator>>(istream &oIstream, MSCString &s);
DTREK_EXPORT istream &getline(istream &oIstream, MSCString &sIn, char cDelimiter = '\n');
inline ostream &operator<<(ostream &o, const MSCString &s)
#else
std::istream &operator>>(std::istream &oIstream, MSCString &s);
DTREK_EXPORT std::istream &getline(std::istream &oIstream, MSCString &sIn, char cDelimiter = '\n');
inline std::ostream &operator<<(std::ostream &o, const MSCString &s)
#endif

#ifdef VC9
{return (o.write(s.string(),(std::streamsize)s.length()));}
#else
{return (o.write(s.string(),s.length()));}
#endif

/****************************************************************************
 *                           typedef for Cstring                            *
 ****************************************************************************/
typedef MSCString Cstring;


/****************************************************************************
 *                           Function prototypes                            *
 * These prototypes are for convenience functions that are not a part of    *
 * MSCString class per se, but are of general use.                          *
 ****************************************************************************/
DTREK_EXPORT int split(const MSCString &s, MSCString sPieces[], int nMaxPieces,
          const char *pcSeparators);
DTREK_EXPORT int split(const MSCString &s, MSCString sPieces[], int nPieces,
          const MSCString &sSeparators);

DTREK_EXPORT int splitcount(const MSCString &s, const char *pcSeparators);
DTREK_EXPORT int splitcount(const MSCString &s, const Cstring &sSeparators);

DTREK_EXPORT MSCString replicate(const char c, const int nRepetitions);
DTREK_EXPORT MSCString replicate(const char *pc, const int nRepetitions);


struct strreg_num {
    enum {DOUBLE,LONG};
    int type;
    union { double d; long l; };
};

class DTREK_EXPORT Stringreg {
    int strreg_locs[10];	// Locations of start of * subsections in the target string.
    int strreg_lens[10];	// The length of these.
    int strwht_locs[10];    // The length of each white space subsection in the target string.
    int strwht_lens[10];
    strreg_num strreg_nums[10];	// Numbers read in.
    int strreg_nums_ct;
    int strreg_loc;
    int strwht_loc;

    int first_pat_loc;      // First character w.r.t. strreg_str

    Cstring strreg_str;

public:

    int nGetPatternStart(); 
    int nGetVarLength(int nVar) { return strreg_lens[nVar]; };
    int nGetVarStart(int nVar)  { return strreg_locs[nVar]; };
    int nParse(const Cstring& sString,const Cstring& sPat);
    int nLoadInt(int nVarNum,int& nValue);
    int nArgCt() { return strreg_loc; };
    int nLoadFloat(int nVarNum,float& fValue);
    int nLoadString(int varnum,Cstring& sBuf);
    int nPrint(Cstring& sFormat,Cstring& sBuf);
};


#endif   /* MSC_MSCSTRING_H */
