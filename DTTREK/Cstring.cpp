// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cstring.cc            Initial author: T.L. Hendrixson          29-Feb-1996
//    This file contains the member functions of class MSCString.
//    This class is intended to be a replacement for the Cstring class from the
//    gnu library and is loosely modelled after the proposed ANSI standard
//    string class.  In order to reduce the amount of recoding, some gnu-isms
//    have been retained (see the header file).  Note that Cstring has been
//    typedef'd as MSCString in the header file Cstring.h, so the two
//    designations can be used interchangably.
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
#include <assert.h>
#include <ctype.h>

#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

#include "Cstring.h"

#include "minmax.h"
#ifdef NOANSIMATH
#include "ansimath.h"
#endif

#if !defined(VC6)
using namespace std;
#endif

/****************************************************************************
 *                           Macro Definitions                              *
 ****************************************************************************/
#if !defined(TRUE) || (1 != (TRUE))
#define TRUE (1)
#endif
#if !defined(FALSE) || (0 != (FALSE))
#define FALSE (0)
#endif

// A little more than 512

int    MSCString::ms_anFreeBinsStart[gs_Cstring_Free_Bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int    MSCString::ms_anFreeBinsCount[gs_Cstring_Free_Bins] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
char*  MSCString::ms_apcFreeBins[gs_Cstring_Max_Free];
int    MSCString::ms_nNumThis = 0;
int    MSCString::ms_nFreeCount = 0;


/****************************************************************************
 *                           Global Variables                               *
 ****************************************************************************/
const int FAILURE = 0;
const int SUCCESS = 1;
MSCString::~MSCString()
{
  vDestroy();
}
/****************************************************************************
 * This routine is the creator for the string class.                        *
 ****************************************************************************/
void MSCString::vCreate(void)
{
#ifndef _MT
  ms_nNumThis++;
#endif
  m_pcString = NULL;
  m_ReservedLength = m_Length = 0;
  if (SUCCESS == nExpand())
    m_pcString[m_Length] = NULL_CHARACTER;
  return;
}

/****************************************************************************
 * Local Heap management functions */
/****************************************************************************/

int MSCString::nGetBin(int nSize) {
    int nx;
    // We don't use first few bins, since our minimal size is gs_Cstring_Minimal_Alloc_Log2
    if (nSize<(1 << gs_Cstring_Minimal_Alloc_Log2))
        return -1;
    for (nx=gs_Cstring_Minimal_Alloc_Log2 ;(nSize>(1 << nx)) && (nx<gs_Cstring_Free_Bins);nx++);
    // If the pointer is larger than the maximum, don't save in a bin.
    if (nx==gs_Cstring_Free_Bins)
        return -1;
    // Otherwise, if the pointer is not a power of 2, then we store it in the next lower bin.
    // This is why allocations should be in powers of 2.
    else if (nSize<(1 << nx))
        nx--;
    return nx;
};

char* MSCString::pcFreeFromBin(int nBin) {
    char* pcRet;
    int nx;
    if (ms_anFreeBinsCount[nBin]==0)
        return NULL;
    pcRet = ms_apcFreeBins[ms_anFreeBinsStart[nBin]];
    nx = gs_Cstring_Max_Free - ms_anFreeBinsStart[nBin]-1;
    if (nx)
        memmove(&ms_apcFreeBins[ms_anFreeBinsStart[nBin]],&ms_apcFreeBins[ms_anFreeBinsStart[nBin]+1],nx*sizeof(char*));
    ms_anFreeBinsCount[nBin]--;
    for (nx=nBin+1;nx<gs_Cstring_Free_Bins;nx++)
        ms_anFreeBinsStart[nx]--;
    ms_nFreeCount--;
    return pcRet;
};

int MSCString::nInsertInBin(char* pcPointer,int nBin) {
    int nx;

    if (ms_nFreeCount==gs_Cstring_Max_Free) {
        // Perhaps we should make room for this pointer.
        if (ms_anFreeBinsCount[nBin]<gs_Cstring_Minimal_Of_Each_Type) {
            int nFullestBin = 0;
            for (nx=0;nx<gs_Cstring_Free_Bins;nx++) {
                if (ms_anFreeBinsCount[nx]>ms_anFreeBinsCount[nFullestBin]) {
                    nFullestBin = nx;
                };
            };
            if (ms_anFreeBinsCount[nFullestBin]>gs_Cstring_Minimal_Of_Each_Type) {
                // Make room by deleteting one pointer.
                delete[] pcFreeFromBin(nFullestBin);
            } else
                return 1;
        } else
            return 1;
    };
    nx = gs_Cstring_Max_Free - ms_anFreeBinsStart[nBin]-1;
    if (nx)
        memmove(&ms_apcFreeBins[ms_anFreeBinsStart[nBin]+1],&ms_apcFreeBins[ms_anFreeBinsStart[nBin]],nx*sizeof(char*));
    ms_apcFreeBins[ms_anFreeBinsStart[nBin]] = pcPointer;
    ms_anFreeBinsCount[nBin]++;
    for (nx=nBin+1;nx<gs_Cstring_Free_Bins;nx++)
        ms_anFreeBinsStart[nx]++;    
    ms_nFreeCount++;
    return 0;
};

void MSCString::vDestroy(void)
{
#ifndef _MT
    ms_nNumThis--;
    if (0 == ms_nNumThis)
    {
#endif
        if (NULL != m_pcString)
        {
            delete [] m_pcString;
            m_pcString = NULL;
        }
#ifndef _MT
        vFreeDelete();      
    } 
    else 
    {
        m_pcString[0] = NULL_CHARACTER;
        
        int     nBin = nGetBin(m_ReservedLength+1);
        
        if ((nBin==-1) || (nInsertInBin(m_pcString,nBin)))
        {
            // No room for it on the free list or it was too big, so delete it
            delete [] m_pcString;
        } 
        m_pcString = NULL;
    }
#endif
    m_ReservedLength = m_Length = 0;
    return;
}

int MSCString::nExpand(const int nRequestedSize)
{
   char *pc = NULL;
   char *pcOld = NULL;
   int nSize;

   if ((0 < nRequestedSize) && (m_ReservedLength > (size_t)nRequestedSize) )
       return SUCCESS;
   
#ifdef _MT

   pc = m_pcString;

   // Make the new size a multiple of 32, which it will be when we add 1
   // when we do the allocation.

   nSize = nRequestedSize | 31;
   m_pcString = new char[nSize+1];

   // If we had an old string, copy it over.

   if(NULL != pc){
      memcpy(m_pcString,pc,m_Length);
      m_pcString[m_Length] = NULL_CHARACTER;
      delete [] pc;
   }

#else

   // Save old pointer for now.  We will delete it later.
   pcOld = m_pcString;
   
   // We need a string of larger size
   nSize = max((1 <<gs_Cstring_Minimal_Alloc_Log2),nRequestedSize+1);

   int      nBin = nGetBin(nSize);
   if ((nBin==-1) || (nSize>(1 << (gs_Cstring_Free_Bins-1)))) {
       // The allocation is just too big.  Bypass...
#ifdef DTREK_MEMORY_DEBUG
       g_nDTMAlloc = 1;
#endif
       m_pcString = new char [nSize];
   } else {
       // For purposes of memory management, we must make the new size a power of 2.
       if (nSize!=(1 << nBin)) {
           nBin++;
           nSize = (1 << nBin);
       };

       // Check this bin, and possibly next bin for a fit.  If neither work, then we must allocate.
       pc = pcFreeFromBin(nBin);
       if ((!pc) && (nBin+1<gs_Cstring_Free_Bins)) {
           nBin++;
           pc = pcFreeFromBin(nBin);
       };
       if (pc) {
           nSize = (1 << nBin);
           m_pcString = pc;
       } else {
           // Allocate.
#ifdef DTREK_MEMORY_DEBUG
           g_nDTMAlloc = 1;
#endif
           m_pcString = new char[nSize];
       };
   };

       // Copy previous data over.
       if (pcOld) 
           memcpy(m_pcString, pcOld, m_Length);
       m_pcString[m_Length] = NULL_CHARACTER;

   // Delete old pointer.
   nBin = nGetBin(m_ReservedLength+1);
   if ((nBin==-1) || (nInsertInBin(pcOld,nBin))) {
       // No room for it on the free list or it was too big, so delete it
       if (pcOld)
           delete [] pcOld;
   } 

#endif  /* _MT */

   m_ReservedLength  = nSize-1;
   return (SUCCESS);
}

void
MSCString::vFreeDelete(void) {
    int nx,ny;
    for (nx=0;nx<gs_Cstring_Free_Bins;nx++) {
        for (ny=0;ny<ms_anFreeBinsCount[nx];ny++) {
            delete[] ms_apcFreeBins[ms_anFreeBinsStart[nx]+ny];
            ms_apcFreeBins[ms_anFreeBinsStart[nx]+ny] = NULL;
        };
        ms_anFreeBinsStart[nx] =0;
        ms_anFreeBinsCount[nx] = 0;
    };
    ms_nFreeCount = 0;
}



/****************************************************************************
 * This routine concatenates one string class variable to the end of        *
 * another string class variable.                                           *
 ****************************************************************************/
MSCString &MSCString::append(const MSCString &s,
                             const int nPosition,
                             int nLength)
{
   int nTotalLength;
   if(0 == s.length())
      return (*this);
   if((size_t)nPosition > s.length() || nPosition < 0){
      cout << "Position out of range when trying to append:\n";
      cout << s << endl;
      exit(1);
   }
   if(nLength < 0)
      nLength = s.length();
   else
      nLength = min((int)s.m_Length-nPosition,nLength);
   nTotalLength = m_Length+nLength;
   if(SUCCESS == nExpand(nTotalLength)){
      memcpy(m_pcString+m_Length,s.m_pcString+nPosition,nLength);
      m_Length = nTotalLength;
      m_pcString[m_Length] = NULL_CHARACTER;
   }
   return (*this);
}
/****************************************************************************
 * This routine concatenates one or more characters to the end of a string  *
 * class variable.                                                          *
 ****************************************************************************/
MSCString &MSCString::append(const char c,
                             const int nRepetitions)
{
   int nTotalLength;
   if(nRepetitions <= 0)
      return (*this);
   nTotalLength = m_Length+nRepetitions;
   if(SUCCESS == nExpand(nTotalLength)){
      memset(m_pcString+m_Length,c,nRepetitions);
      m_Length = nTotalLength;
      m_pcString[m_Length] = NULL_CHARACTER;
   }
   return (*this);
}
/****************************************************************************
 * This routine concatenates a character string to the end of a string      *
 * class variable.                                                          *
 ****************************************************************************/
MSCString &MSCString::append(const char *pc,
                             int nLength)
{
   int n;
   int nTotalLength;
   n = strlen(pc);
   if(nLength < 0 || nLength > n)
      nLength = n;
   if(0 == (nTotalLength=nLength))
      return (*this);
   nTotalLength += m_Length;
   if(SUCCESS == nExpand(nTotalLength)){
      memcpy(m_pcString+m_Length,pc,nLength);
      m_Length = nTotalLength;
      m_pcString[m_Length] = NULL_CHARACTER;
   }
   return (*this);
}
/****************************************************************************
 * This routine appends the string representation of a long integer to a    *
 * string.                                                                  *
 ****************************************************************************/
MSCString &MSCString::append(const long l,
                             int nFieldWidth)
{
   int n;
   char *format = NULL;
// calculate the minimum fieldwidth required to store this number
   //n = (0 == l) ? 1 : (1+(int)log10(f   (double)l)));
   //if(0 > l)  // negative sign
      //n++;
// Do it this way because g++ compiler breaks code if compiled -O3 using
// above code
   if(0 == l)
      n = 1;
   else if(0 > l)
      n = 2+(int)log10(fabs((double)l));
   else
      n = 1+(int)log10(fabs((double)l));
// Was supplied fieldwidth big enough?
   if(nFieldWidth < n)
      nFieldWidth = n;
// generate format string
   n = (int)log10((double)nFieldWidth)+5; // %#ld\0
   format = new char[n];
   sprintf(format,"%%%dld",nFieldWidth);
// expand MSCString and print number to MSCString
   if(SUCCESS == nExpand(m_Length+nFieldWidth)){
      sprintf(m_pcString+m_Length,format,l);
      m_Length += nFieldWidth;
   }
   delete []format;
   return (*this);
}
/****************************************************************************
 * This routine appends the string representation of a double floating      *
 * point number to a string.                                                *
 ****************************************************************************/
MSCString &MSCString::append(const double d,
                             int nFieldWidth,
                             const int nDecimal)
{
   int n,ndec;
   char *format;
// calculate the minimum fieldwidth required to store this number
// nDecimal (behind the decimal) + decimal point + integer part
   ndec = (nDecimal < 0) ? 0 : nDecimal;
   n = ndec+1;
   if(1. > fabs(d))
      n++;  // one space for the zero in front of the decimal point
   else
      n += 1+(int)log10(fabs(d));
   if(0 > d)  // negative sign
      n++;
// Was supplied fieldwidth big enough?
   if(nFieldWidth < n)
      nFieldWidth = n;
// generate format string
   n = (int)log10((double)nFieldWidth)+(int)log10((double)ndec)+7; // %#.#lf\0
   format = new char[n];
   sprintf(format,"%%%d.%dlf",nFieldWidth,ndec);
// expand MSCString and print number to MSCString
   if(SUCCESS == nExpand(m_Length+nFieldWidth)){
      sprintf(m_pcString+m_Length,format,d);
      m_Length += nFieldWidth;
   }
   delete []format;
   return (*this);
}
/****************************************************************************
 * This routine copies the specified character to the original string,      *
 * copying the character nRepetitions times.                                *
 ****************************************************************************/
MSCString &MSCString::assign(const char c,
                             const int nRepetitions)
{
   if ( (0 >= nRepetitions) || (NULL_CHARACTER == c) )
     {
       m_Length = 0;
       m_pcString[m_Length] = NULL_CHARACTER;
     }
   else if (SUCCESS == nExpand(nRepetitions))
     {
       memset(m_pcString, c, nRepetitions);
       m_Length = nRepetitions;
       m_pcString[m_Length] = NULL_CHARACTER;
     }
   return (*this);
}
/****************************************************************************
 * This routine copies the specified character sequence to the original     *
 * string, copying nLength characters, or all of the specified character    *
 * sequence if nLength is less than zero.                                   *
 ****************************************************************************/
MSCString &MSCString::assign(const char *pc,
                             int nLength)
{
   int n;
   if (NULL != pc)
     n = strlen(pc);
   else
     n = 0;
   if(nLength < 0 || nLength > n)
      nLength = n;
   if(0 == nLength){
      m_Length = 0;
      m_pcString[m_Length] = NULL_CHARACTER;
   }
   else{
      if(SUCCESS == nExpand(nLength)){
         memcpy(m_pcString,pc,nLength);
         m_Length = nLength;
         m_pcString[m_Length] = NULL_CHARACTER;
      }
   }
   return (*this);
}
/****************************************************************************
 * This routine copies the specified string to the original string,         *
 * starting with the nPosition offset of the specified string, copying      *
 * nLength characters, or all of the remaining specified string if nLength  *
 * is less than zero.                                                       *
 ****************************************************************************/
MSCString &MSCString::assign(const MSCString &s,
                             const int nPosition,
                             int nLength)
{
   if (0 == s.length()){
      if (m_ReservedLength == 0)
          nExpand();
      m_Length = 0;
      m_pcString[m_Length] = NULL_CHARACTER;
   }
   else if((size_t)nPosition > s.length() || nPosition < 0)
     {
      cout << "Position out of range when trying to assign:\n";
      cout << s << endl;
       m_Length = 0;
       m_pcString[m_Length] = NULL_CHARACTER;
       exit(1);
     }
   else
     {
       if(nLength < 0)
         nLength = s.length()-nPosition;
       else
         nLength = min((int)s.m_Length-nPosition,nLength);
       if(SUCCESS == nExpand(nLength)){
         memcpy(m_pcString,s.m_pcString+nPosition,nLength);
         m_Length = nLength;
         m_pcString[m_Length] = NULL_CHARACTER;
       }
     }
   return (*this);
}
/****************************************************************************
 * This routine generates the string representation of a long integer.      *
 ****************************************************************************/
MSCString &MSCString::assign(const long l,
                             int nFieldWidth)
{
   int n;
   char *format = NULL;
// calculate the minimum fieldwidth required to store this number
   //n = (0 == l) ? 1 : (1+(int)log10(fabs((double)l)));
   //if(0 > l)  // negative sign
      //n++;
// g++ compiler breaks on above code using -O3 compilation flag
   if(0 == l)
      n = 1;
   else if(0 > l)
      n = 2+(int)log10(fabs((double)l));
   else
      n = 1+(int)log10(fabs((double)l));
// Was supplied fieldwidth big enough?
   if(nFieldWidth < n)
      nFieldWidth = n;
// generate format string
   n = (int)log10((double)nFieldWidth)+5; // %#ld\0
   format = new char[n];
   sprintf(format,"%%%dld",nFieldWidth);
// expand MSCString and print number to MSCString
   if(SUCCESS == nExpand(nFieldWidth)){
      sprintf(m_pcString,format,l);
      m_Length = nFieldWidth;
   }
   delete []format;
   return (*this);
}
MSCString &MSCString::assign(const unsigned long ul,
                             int nFieldWidth)
{
   int n;
   char *format = NULL;
// calculate the minimum fieldwidth required to store this number
   //n = (0 == l) ? 1 : (1+(int)log10(fabs((double)l)));
   //if(0 > l)  // negative sign
      //n++;
// g++ compiler breaks on above code using -O3 compilation flag
   if(0 == ul)
      n = 1;
//   else if(0 > ul)
//      n = 2+(int)log10(fabs((double)ul));
   else
      n = 1+(int)log10(fabs((double)ul));
// Was supplied fieldwidth big enough?
   if(nFieldWidth < n)
      nFieldWidth = n;
// generate format string
   n = (int)log10((double)nFieldWidth)+5; // %#lu\0
   format = new char[n];
   sprintf(format,"%%%dlu",nFieldWidth);
// expand MSCString and print number to MSCString
   if(SUCCESS == nExpand(nFieldWidth)){
      sprintf(m_pcString,format,ul);
      m_Length = nFieldWidth;
   }
   delete []format;
   return (*this);
}
/****************************************************************************
 * This routine generates the string representation of a double floating    *
 * point number.                                                            *
 ****************************************************************************/
MSCString &MSCString::assign(const double d,
                             int nFieldWidth,
                             const int nDecimal)
{
   int n,ndec;
   char *format;
// calculate the minimum fieldwidth required to store this number
// nDecimal (behind the decimal) + decimal point + integer part
   ndec = (nDecimal < 0) ? 0 : nDecimal;
   n = ndec+1;
   if(1. > fabs(d))
      n++;  // one space for the zero in front of the decimal point
   else
      n += 2+(int)log10(fabs(d));
   if(0 > d)  // negative sign
      n++;
// Was supplied fieldwidth big enough?
   if(nFieldWidth < n)
      nFieldWidth = n;
// generate format string
   n = (int)log10((double)nFieldWidth)+(int)log10((double)ndec)+7; // %#.#lf\0
//+29-Mar-2000 jwp: prevent bogus values of n
   if (1 > n) n = 1;
   if (3000 < n) n = 3000;
//-jwp
   format = new char[n];
   sprintf(format,"%%%d.%dlf",nFieldWidth,ndec);
// expand MSCString and print number to MSCString
   if(SUCCESS == nExpand(nFieldWidth)){
      sprintf(m_pcString,format,d);
      m_Length = nFieldWidth;
   }
   delete []format;
   return (*this);
}

/****************************************************************************
 * This routine parses an input string that might be in any of various      *
 * formats, including the fractional format                                 *
 ****************************************************************************/

bool MSCString::bParseFloat(float& fFloat) {
    char* pcData;


    pcData = string();

    // Check to see if we are just dealing with a vanilla floating point or integer.
    if (strlen(pcData)==strspn(pcData,"\t\n\r .0123456789e+-")) {
        if (1!=sscanf(pcData,"%f",&fFloat))
            return FALSE;
        else
            return TRUE;
    };
    // Check to see if we are dealing with a fractional number.
    if (strlen(pcData)==strspn(pcData,"\t\n\r0123456789.+-/")) {
        float fNumer;
        float fDenom;
        char* pcSlash;
        // We want to parse the numerator and denomenator.  To do this, we find the
        // '/' character, and set it to '\x0'.  Then, we parse the numerator and denomenator.
        // If all goes well, the '/' character is reinstated.
        pcSlash = strchr(pcData,'/');
        if (!pcSlash)
            return FALSE;
        *pcSlash = '\x0';
        if ((1!=sscanf(pcData,"%f",&fNumer)) || (1!=sscanf(pcSlash+1,"%f",&fDenom))) {
            *pcSlash = '/';
            return FALSE;
        };
        *pcSlash = '/';
        fFloat = ((float) fNumer)/((float) fDenom);
        return TRUE;
    };
    // Other formats can be added here...

    return FALSE;
};

/****************************************************************************
 * This routine tries to find the best fractional representation of a number*
 * within tollerance.  If it cannot find a reasonable number, then it simply*
 * prints out the floating point representation in (*this)                  *
 * Returns TRUE if it works                                                 *
 ****************************************************************************/

bool MSCString::bMakeFloat(float fFloat,float fRoundError,int nMaxNumer,int nMaxDenom)
{
    int nNumer,nDenom;
    int nChooseNumer,nChooseDenom;
    bool bSolutionFound;
    float fTestFloat;
    char pcBuf[20];

    bSolutionFound = FALSE;


    // See if the number is near an integer.
    if (fabs(fFloat - floor(fFloat))<fRoundError) {
        sprintf(pcBuf,"%d",(int) floor(fFloat));
        (*this) = pcBuf;
        return TRUE;
    };
    if (fabs(fFloat - (1+floor(fFloat))) < fRoundError) {
        sprintf(pcBuf,"%d",(int) (1+floor(fFloat)));
        (*this) = pcBuf;
        return TRUE;
    };
    
    for (nNumer=-nMaxNumer;nNumer<=nMaxNumer;nNumer++) {
        for (nDenom=2;nDenom<=nMaxDenom;nDenom++) {
            fTestFloat = ((float) nNumer)/((float) nDenom);
            if (fabs(fTestFloat - fFloat)< fRoundError) {
                if (!bSolutionFound) {
                    nChooseNumer = nNumer;
                    nChooseDenom = nDenom;
                    bSolutionFound = TRUE;
                } else if (fabs((double)nNumer)+fabs((double)nDenom)<fabs((double)nChooseNumer)+fabs((double)nChooseDenom)) {
                    nChooseNumer = nNumer;
                    nChooseDenom = nDenom;
                };
            };                
        };
    };

    if (bSolutionFound) {
        sprintf(pcBuf,"%d/%d",nChooseNumer,nChooseDenom);
        (*this) = pcBuf;
        return TRUE;
    } else {
        sprintf(pcBuf,"%6.3f",fFloat);
        (*this) = pcBuf;
        return FALSE;
    };
};





/****************************************************************************
 * This routine returns the lexigraphical comparison of the original string *
 * with the indicated string.                                               *
 ****************************************************************************/
int MSCString::compare(const MSCString &s) const
{
   return (strcmp(m_pcString,s.string()));
}
/****************************************************************************
 * This routine returns the lexigraphical comparison of the original string *
 * with the indicated character sequence.                                   *
 ****************************************************************************/
int MSCString::compare(const char *pc) const
{
   return (strcmp(m_pcString,pc));
}

/****************************************************************************
 * This routine returns the lexigraphical comparison of the original string *
 * with the indicated string.                                               *
 ****************************************************************************/
int MSCString::compare_no_case(const MSCString &s) const
{
#ifdef WIN32
   return _strcmpi(m_pcString,s.string());
#else
   return strcasecmp(m_pcString,s.string());
#endif
}
/****************************************************************************
 * This routine returns the lexigraphical comparison of the original string *
 * with the indicated character sequence.                                   *
 ****************************************************************************/
int MSCString::compare_no_case(const char *pc) const
{
#ifdef WIN32
   return _strcmpi(m_pcString,pc);
#else
   return strcasecmp(m_pcString,pc);
#endif
}

/****************************************************************************
 * This routine finds the first offset of the indicated character, starting *
 * at the offset position in the original string.                           *
 ****************************************************************************/
int MSCString::find(const char c,
                    const int nPosition) const
{
   char *pc;
   if(0 > nPosition || (size_t)nPosition >= m_Length)
      return -1;
   pc = strchr(m_pcString+nPosition,(int)c);
   return ((NULL == pc) ? -1 : pc-m_pcString);
}
/****************************************************************************
 * This routine finds the first offset of the indicated character sequence, *
 * starting at the offset position in the original string.                  *
 ****************************************************************************/
int MSCString::find(const char *pc,
                    const int nPosition,
                    int nLength) const
{
   char *pcAt,*pcEnd;
   int n;
/*
 * If no comparison length specified, use length of indicated sequence.
 */
   if(0 > nLength)
      nLength = strlen(pc);
/*
 * If meaningless start postion in original string or indicated string has no
 * length, return nothing found.
 */
   if(0 == nLength || 0 > nPosition || (size_t)nPosition >= m_Length)
      return -1;
/*
 * Figure out the last position in the original string at which a comparison
 * has any meaning.  If no comparisons can be meaningful, then return nothing
 * found.
 */
   n = m_Length-nPosition;
   if(n < nLength)
      return -1;
/*
 * If only one position in the original string can be compared, then no need
 * to do the rest of this routine.
 */
   pcAt = m_pcString+nPosition;
   if(n == nLength){
      return ((0 == memcmp(pcAt,pc,nLength)) ? (pcAt - m_pcString) : -1);
   }
/*
 * Keep looking for the string until we find it or we run out of things to
 * compare.  Note that you can't use strstr() since you may be looking only
 * for a part of pc.
 */
   pcEnd = pcAt+(n-nLength)+1;
   for(NULL; pcAt < pcEnd; pcAt++){
      if(NULL == (pcAt=strchr(pcAt,(int)(*pc))))
         break;
      if (pcAt >= pcEnd)
          break;
      if(0 == memcmp(pcAt,pc,nLength))
         return (pcAt-m_pcString);
   }
   return -1;
}
/****************************************************************************
 * This routine finds the first offset of the indicated string, starting at *
 * the offset position in the original string.                              *
 ****************************************************************************/
int MSCString::find(const MSCString &s,
                    const int nPosition) const
{
   int n;
   n = s.length();
   return (find(s.string(),nPosition,n));
}


/****************************************************************************
 * This routine finds the last offset of the indicated character, starting  *
 * at the offset position in the original string.                           *
 ****************************************************************************/
int MSCString::nReverseFind(const char c,
                    int nPosition) const
{
   char *pc, *pcSearchStr;
   int idx;

   /*
    * If no position specified, use last position in string.
    */
   if(nPosition == -1) {
       nPosition = m_Length-1;
   }

   if(0 > nPosition || (size_t)nPosition >= m_Length)
      return -1;

   // Make a copy of the search string terminating one character after
   // index nPosition.
   pcSearchStr = new char[nPosition+2]; // +1 for NULL terminator
   strncpy(pcSearchStr,m_pcString,nPosition+1);
   pcSearchStr[nPosition+1] = '\0';

   // Do the reverse search
   pc = strrchr(pcSearchStr,(int)c);
   idx = (pc) ? pc-pcSearchStr : -1;
   
   delete pcSearchStr;
   return idx;
}

/****************************************************************************
 * This routine finds the last offset of the indicated character sequence,  *
 * starting at the offset position in the original string.                  *
 ****************************************************************************/
int MSCString::nReverseFind(const char *pc,
                    int nPosition,
                    int nLength) const
{
   char *pcAt,*pcEnd;
   //int n;

   /*
    * If no position specified, use last position in string.
    */
   if(nPosition == -1) {
       nPosition = m_Length-1;
   }

/*
 * If no comparison length specified, use length of indicated sequence.
 */
   if(0 > nLength)
      nLength = strlen(pc);
/*
 * If meaningless start postion in original string or indicated string has no
 * length, return nothing found.
 */
   if(0 == nLength || 0 > nPosition || (size_t)nPosition >= m_Length)
      return -1;
/*
 * Figure out the last position in the original string at which a comparison
 * has any meaning.  If no comparisons can be meaningful, then return nothing
 * found.
 */
   //if(nPosition+1 < nLength)
   //   return -1;
/*
 * If only one position in the original string can be compared, then no need
 * to do the rest of this routine.
 */
   pcAt = m_pcString+nPosition;//-nLength;
   //if(nPosition+1 == nLength){
   //   return ((0 == memcmp(pcAt,pc,nLength)) ? pcAt-m_pcString : -1);
   //}
/*
 * Keep looking for the string until we find it or we run out of things to
 * compare.  Note that you can't use strstr() since you may be looking only
 * for a part of pc.
 */
   pcEnd = m_pcString-1;
   for(NULL; pcAt > pcEnd; pcAt--){
       // Find last occurrence of pc[0] beginning the search at pcAt.
       if(NULL == (pcAt=m_pcString+nReverseFind((int)(*pc),pcAt-m_pcString)))
           break;
       if(0 == memcmp(pcAt,pc,nLength))
           return (pcAt-m_pcString);
   }
   return -1;
}

/****************************************************************************
 * This routine finds the last offset of the indicated string, starting at  *
 * the offset position in the original string.                              *
 ****************************************************************************/
int MSCString::nReverseFind(const MSCString &s,
                    int nPosition) const
{
   int n;
   n = s.length();

   return (nReverseFind(s.string(),nPosition,n));
}

/****************************************************************************
 * This routine provides an alternate pathway into Stringreg
 ****************************************************************************/

bool MSCString::bRegularMatch(const Cstring& sPattern) const  {
    static Stringreg oReg;

    return (0==oReg.nParse(*this,sPattern));
};

/****************************************************************************
 * These functions remove characters from a string.
 ****************************************************************************/

MSCString& MSCString::remove(const int nPosition, int nLength) {
    Cstring sTemp;
    if ((nPosition<0) || (nPosition>=m_Length))
        return (*this);
    sTemp = before(nPosition);
    if (nLength!=-1)
        sTemp += after(min((int) length()-1,nPosition+nLength-1));
    (*this) = sTemp;
    return (*this);
};

MSCString& MSCString::remove(const char c) {
    Cstring sResidual;
    int nPosition;
    
    sResidual = (*this);
    (*this) = "";
    while (-1 !=(nPosition = sResidual.find(c))) {
        (*this) += sResidual.before(nPosition);
        sResidual = sResidual.after(nPosition);
    }
    (*this) += sResidual;   
    return (*this);
};

/****************************************************************************
 * These routines will replace string values with new text, or do nothing   *
 * if a string is not found                                                 *
 ****************************************************************************/



MSCString& MSCString::replace(const int nPosition,const int nLength,MSCString sNew) {
    Cstring sTemp;
    if ((nPosition < 0) || (nPosition + nLength > length()) || (!nLength))
        return *this;
    sTemp = this->before(nPosition);
    sTemp += sNew;
    sTemp += this->after(nPosition + nLength - 1);
    (*this) = sTemp;
    return *this;

};

MSCString& MSCString::replace(MSCString sToFind,const int nStart,MSCString sNew) {
    return replace(this->find(sToFind,nStart),sToFind.length(),sNew);
};



/****************************************************************************
 * This routine returns a substring containing all of the string from the   *
 * start up to and including the indicated string.                          *
 ****************************************************************************/
MSCString MSCString::through(const MSCString &s)
{
   int n;
   n = find(s.string());
   return ((0 <= n) ? (substr(0,n+s.length())) : MSCString(""));
}
/****************************************************************************
 * This routine returns a substring containing the rightmost n characters.  *
 ****************************************************************************/
MSCString MSCString::Right(int n)
{
   size_t index;
   if(0 >= n)
      return MSCString( "" );
   if((size_t)n <= m_Length){
      index = m_Length-(size_t)n;
      return MSCString( &(m_pcString[index]) );
   }
   else
      return MSCString( "" );
}
/****************************************************************************
 * This routine returns a substring containing all of the string from the   *
 * start up to and including the indicated character sequence.              *
 ****************************************************************************/
MSCString MSCString::through(const char *pc)
{
   int n;
   n = find(pc);
   return ((0 <= n) ? substr(0,n+strlen(pc)) : MSCString(""));
}
/****************************************************************************
 * This routine returns a substring containing all of the string from the   *
 * start up to and including the indicated character.                       *
 ****************************************************************************/
MSCString MSCString::through(const char c)
{
   int n;
   n = find(c);
   return ((0 <= n) ? substr(0,n+1) : MSCString(""));
}
/****************************************************************************
 * This routine returns a substring containing all of the string after the  *
 * indicated string.                                                        *
 ****************************************************************************/
MSCString MSCString::right(const MSCString &s)
{
   int n;
   n = find(s.string());
   return ((0 <= n) ? substr(n+s.length()) : MSCString(""));
}
/****************************************************************************
 * This routine returns a substring containing all of the string after the  *
 * indicated character sequence.                                            *
 ****************************************************************************/
MSCString MSCString::right(const char *pc)
{
   int n;
   n = find(pc);
   return ((0 <= n) ? substr(n+strlen(pc)) : MSCString(""));
}
/****************************************************************************
 * This routine returns a substring containing all of the string after the  *
 * indicated character.                                                     *
 ****************************************************************************/
MSCString MSCString::right(const char c)
{
   int n;
   n = find(c);
   return ((0 <= n) ? substr(n+1) : MSCString(""));
}
/****************************************************************************
 * This routine defines the string version of the i/o operator >>.          *
 ****************************************************************************/
istream &operator>>(istream &oIstream,
                    Cstring &sIn)
{
   sIn.assign(NULL_CHARACTER);   // Set string to null character

   // Eat any initial whitespace
   char c;
   do{
      oIstream.get(c);
   }while(0 != isspace(static_cast<int>(c)) && oIstream.good());

   // Read until we hit whitespace
   if(oIstream.good()){
      do{
         sIn.append(c);
         oIstream.get(c);
      }while(0 == isspace(static_cast<int>(c)) && oIstream.good());

      // If we read a newline character, put it back
      if(oIstream.good()){
         if('\n' == c)
            oIstream.putback(c);
      }
   }

   return (oIstream);
}
/****************************************************************************
 * This routine defines the string version of the i/o routine getline().    *
 ****************************************************************************/
istream &getline(istream &oIstream,
                 MSCString &sIn,
                 char cDelimiter)
{
   sIn.assign(NULL_CHARACTER);

   char c;
   oIstream.get(c);
   while(oIstream.good() && cDelimiter != c){
      sIn.append(c);
      oIstream.get(c);
   }

   return (oIstream);
}
/****************************************************************************
 * This routine converts all of the characters in a string to lowercase.    *
 ****************************************************************************/
void MSCString::downcase(void)
{
   size_t n;
   for(n = 0; n < m_Length; n++){
      m_pcString[n] = tolower(m_pcString[n]);
      //if(0 != isupper(m_pcString[n]))
         //m_pcString[n] += 'a'-'A';
   }
   return;
}
/****************************************************************************
 * This routine converts all of the characters in a string to uppercase.    *
 ****************************************************************************/
void MSCString::upcase(void)
{
   size_t n;
   for(n = 0; n < m_Length; n++){
      m_pcString[n] = toupper(m_pcString[n]);
      //if(0 != islower(m_pcString[n]))
         //m_pcString[n] += 'A'-'a';
   }
   return;
}
/****************************************************************************
 ****************************************************************************
 **                                                                        **
 ** The following are convenience functions and are not directly           **
 ** a part of the string class.                                            **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************
 * This routine splits a MSCstring up into pieces according to the separa-  *
 * tors passed.  Pieces are formed when any character in the separator      *
 * sequence is found in the MSCstring.  nMaxPieces is the maximum number of *
 * pieces that will be generated.                                           *
 ****************************************************************************/
int split(const MSCString &sString,
          MSCString sPieces[],
          int nMaxPieces,
          const MSCString &sSeparators)
{
   int i,j,k,kk,nNumSeparators;
   MSCString sTemp;
   nNumSeparators = sSeparators.length();
// remove leading separators from the string
   sTemp = sString;
   while(-1 != sSeparators.find(sTemp.GetAt(0)) && "" != sTemp)
      sTemp = sTemp.substr(1);
   if("" == sTemp)
      return 0;
// break the string into pieces
   for(i = 0; i < nMaxPieces; i++)
     {
     // remove leading separators
       while(-1 != sSeparators.find(sTemp.GetAt(0)) && "" != sTemp)
          sTemp = sTemp.substr(1);
       if("" == sTemp)
          return (i);
     // find earliest position of any separator in the string
       k = sTemp.length()+1;
       for(j = 0; j < nNumSeparators; j++)
         {
           kk = sTemp.find(sSeparators.GetAt(j));
           if(0 < kk)
             k = min(k,kk);
         }
     // If no more separators in string, put string in pieces array and leave
       if(sTemp.length() < (size_t)k)
         {
           sPieces[i++] = sTemp;
           break;
         }
     // Put the segment of the string that is in front of the separator in
     // the pieces array, and remove this piece and separator from the string
       sPieces[i] = sTemp.before(k);
       sTemp = sTemp.after(k);
     }
// return the number of pieces found
   return (i);
}
/****************************************************************************
 * This version of split() uses char * for the separators.                  *
 ****************************************************************************/
int split(const MSCString &s,
          MSCString sPieces[],
          int nMaxPieces,
          const char *pcSeparators)
{
   MSCString sSep = pcSeparators;
   return (split(s,sPieces,nMaxPieces,sSep));
}
/****************************************************************************
 * This routine counts the number of pieces that would be formed if a       *
 * MSCString was split up into pieces according to the separators passed.   *
 * Pieces are formed when any character in the separator sequence is found  *
 * in the MSCstring.                                                        *
 ****************************************************************************/
int splitcount(const MSCString &sString,
               const MSCString &sSeparators)
{
   int i,j,k,kk,nNumSeparators;
   MSCString sTemp;
   nNumSeparators = sSeparators.length();
// remove leading separators from the string
   sTemp = sString;
   while(-1 != sSeparators.find(sTemp.GetAt(0)) && "" != sTemp)
      sTemp = sTemp.substr(1);
   if("" == sTemp)
      return 0;
// break the string into pieces
   for(i = 0; i < INT_MAX; i++)
     {
     // remove leading separators
       while(-1 != sSeparators.find(sTemp.GetAt(0)) && "" != sTemp)
          sTemp = sTemp.substr(1);
       if("" == sTemp)
          return (i);
     // find earliest position of any separator in the string
       k = sTemp.length()+1;
       for(j = 0; j < nNumSeparators; j++)
         {
           kk = sTemp.find(sSeparators.GetAt(j));
           if(0 < kk)
             k = min(k,kk);
         }
     // If no more separators in string, increment counter and leave
       if(sTemp.length() < (size_t)k)
         {
           i++;
           break;
         }
     // Remove this piece and separator from the string
       sTemp = sTemp.after(k);
     }
// return the number of pieces found
   return (i);
}
/****************************************************************************
 * This version of splitcount() uses char * for the separators.             *
 ****************************************************************************/
int splitcount(const MSCString &s,
               const char *pcSeparators)
{
   MSCString sSep = pcSeparators;
   return (splitcount(s,sSep));
}
/****************************************************************************
 * This routine generates a string composed of a repeated character         *
 * sequence.                                                                *
 ****************************************************************************/
MSCString replicate(const char *pc,
                    const int nRepetitions)
{
   MSCString replica;
   int i;
   replica = pc;
   for(i = 1; i < nRepetitions; i++)
      replica += pc;
   return replica;
}
/****************************************************************************
 * This routine generates a string composed of one character repeated.      *
 ****************************************************************************/
MSCString replicate(const char c,
                    const int nRepetitions)
{
   return (MSCString(c,nRepetitions));
}


/*****************************************************************************
 *
 * Returns the path or directory of a filename.  For example:
 * C:\MyDir\MyFile.txt ---> C:\MyDir
 *****************************************************************************/
MSCString MSCString::sGetFilePath() const
{
    int idx, i, j;
    MSCString str(m_pcString);

    i = nReverseFind('\\');
    j = nReverseFind('/');
    idx = i>j?i:j;

    if( idx > -1 ) {
        return str.left(idx);    
    }
    str = "";
    return str;
}

/*****************************************************************************
 * Returns the extension of a filename.  For example:
 * C:\MyDir\MyFile.txt ---> txt
 *****************************************************************************/

MSCString MSCString::sGetFileExt() const
{
    MSCString str(m_pcString);
    int idx = str.nReverseFind('.');

    if( idx > -1 ) {
        return str.from(idx+1);    
    }
    str = "";
    return MSCString();
}

/*****************************************************************************
 * Returns the filename without directory.  For example:
 * C:\MyDir\MyFile.txt ---> MyFile.txt
 *****************************************************************************/

MSCString MSCString::sGetFileName() const
{
    MSCString str(m_pcString);
    int idx, i, j;

    i = nReverseFind('\\');
    j = nReverseFind('/');
    idx = i>j?i:j;

    if( idx > -1 ) {
        return str.from(idx+1); 
    }
    return str;
}

/*****************************************************************************
 * Returns the title of a filename.  For example:
 * C:\MyDir\MyFile.txt ---> MyFile
 *****************************************************************************/

MSCString MSCString::sGetFileTitle() const
{
    MSCString title = sGetFileName();
    int idx = title.nReverseFind('.');
    if( idx > -1 ) {
        return title.left(idx);
    }
    return title;
}

// Similar to the MFC CString Format method.
// This is _NOT_ UNIX compatible!!!

#ifdef WIN32
void MSCString::Format( char* sFormat, ... )
{
    size_t size = 64;
    char* buffer = new char[size+1];

    va_list marker;
    // Initialize variable arguments.
    va_start( marker, sFormat );
    // Loop until the buffer is big enough.
    while( _vsnprintf( buffer, size, sFormat, marker ) == -1 ) {
        buffer[size] = '\0';
        delete [] buffer;
        size*=2;
        buffer = new char[size+1];
    }
    buffer[size] = '\0';
    assign( buffer );

    delete [] buffer;
}
#endif


int MSCString::nRead(FILE* pFIn) {
    char pcBuf[50];
    char* pcStat;
    int nx;
    (*this) = "";
    if ((feof(pFIn)) || (ferror(pFIn)))
        return 1;
    do {
    if (!(pcStat = fgets(&pcBuf[0],50-1,pFIn))) {
        if (ferror(pFIn))
            return 1;
    };
    nx = strlen(pcBuf);
    if ((nx) && (pcStat))
        (*this) += pcStat;
    } while ((pcStat) && (nx) && (pcBuf[nx-1]!='\n'));
    return 0;
};

// Parse a string as a list of numbers and save them in a vector
int MSCString::nListToVector(std::vector<MSCString>& saArray, char* pcListSeparators)
{
    MSCString       sCopy(m_pcString);

    saArray.clear(); 
    MSCString      ss;

    char*       pc = strtok(sCopy.m_pcString, pcListSeparators);
    if( !pc )
        return -1;
    
    ss = pc;
    saArray.push_back(ss);
    
    while(NULL != (pc = strtok(NULL, pcListSeparators)) )
    {
        ss = pc;
        saArray.push_back(ss);
    }

    return 0;
}
//////////////////////////////////////////////////////////////////////////////
int MSCString::nListToVector(std::vector<int>& naArray, char* pcListSeparators)
{
    MSCString       sCopy(m_pcString);

    naArray.clear(); 
    int      nn=0;

    char*       pc = strtok(sCopy.m_pcString, pcListSeparators);
    if( !pc )
        return -1;
    
    nn = atoi(pc);
    naArray.push_back(nn);
    
    while(NULL != (pc = strtok(NULL, pcListSeparators)) )
    {
        nn = atoi(pc);
        naArray.push_back(nn);
    }

    return 0;
}
//////////////////////////////////////////////////////////////////////////////
int MSCString::nListToVector(std::vector<double>& daArray, char* pcListSeparators)
{
    MSCString       sCopy(m_pcString);
    daArray.clear(); 
    double      dd=0;

    char*       pc = strtok(sCopy.m_pcString, pcListSeparators);
    if( !pc )
        return -1;
    
    dd = atof(pc);
    daArray.push_back(dd);
    
    while(NULL != (pc = strtok(NULL, pcListSeparators)) )
    {
        dd = atof(pc);
        daArray.push_back(dd);
    }

    return 0;
}
//////////////////////////////////////////////////////////////////////////////
int MSCString::nListToVector(std::vector<double>& daArray, char cListSeparator)
{
    MSCString       sCopy(m_pcString);
    daArray.clear(); 
    double      dd=0;

    char*       pc = strtok(sCopy.m_pcString, &cListSeparator);
    if( !pc )
        return -1;
    
    dd = atof(pc);
    daArray.push_back(dd);
    
    while(NULL != (pc = strtok(NULL, &cListSeparator)) )
    {
        dd = atof(pc);
        daArray.push_back(dd);
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int MSCString::nListToVectorOfPairs(std::vector< std::pair<int,int> >& vecPairsOfInts, char* pcListSeparators, char* pcPairSeparators)
{
    vecPairsOfInts.clear();
    
    std::vector<Cstring>   vecStringsOfIntPairs;
    
    if( 0 != nListToVector(vecStringsOfIntPairs, pcListSeparators) )
        return 1;

    std::vector<int>       vecInts;    
    int     n0 = -1;
    int     n1 = -1;
    for(int ii=0; ii < vecStringsOfIntPairs.size(); ii++)
    {
        vecStringsOfIntPairs[ii].nListToVector(vecInts, pcPairSeparators);

        if( vecInts.size() == 1 )
        {
            n0 = vecInts[0];
            n1 = n0;
        }
        else if (vecInts.size() == 2 )
        {
            n0 = vecInts[0];
            n1 = vecInts[1];
        }

        vecPairsOfInts.push_back(std::pair<int,int>(n0, n1));
    }
    
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool MSCString::bIsNumeric()
{
    if( m_Length == 0 )
        return false;
    
    if( m_pcString[m_Length] != '\0' )
        return false; // should not happen

    char    cc;
    for(int ii=0; ii < m_Length; ii++) 
    {
        cc = m_pcString[ii];

        if( cc < '0' || cc > '9' )
        {
            if( cc == '.' || (cc == '-' && ii == 0) )
                continue;
            else
                return false;
        }
    }
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////
// If the original filepath is longer than cnMaxNumberOfFilePathChars
// it will be abbreviated. Currently supported scheme of abbreviation is as follows:
//
// 1) Create an empty string S, N chars long, where N is the max allowed number of file path chars.
// 2) Copy the first (N/2 - 3) characters from the beginning of the file path to the first (N/2 - 3)
//    characters of the string S.
// 3) Append three dots to the string S.
// 4) Copy N/2 chars from the end of the file path to the right N/2 chars of string S. 
void MSCString::vAbbreviateFilePath(const int cnMaxNumberOfFilePathChars)
{
    if( m_Length <= cnMaxNumberOfFilePathChars )
        return;  // nothing to do
    
    char*       pcBuffer = new char[cnMaxNumberOfFilePathChars+1];  // +1 for end NULL character
    
    int         ii = 0;
    
    ////////////////////////////////////////////////////////
    // Initialize to spaces
    for(ii = 0; ii < cnMaxNumberOfFilePathChars; ii++)
    {
        pcBuffer[ii] = ' ';
    }
    pcBuffer[cnMaxNumberOfFilePathChars] = '\0';
    ////////////////////////////////////////////////////////

    const int cnNumberOfDots = 3; 
    const int cnNumberOfCharsBeforeBreak = cnMaxNumberOfFilePathChars/2;

    // A safety check
    if( cnNumberOfDots > cnNumberOfCharsBeforeBreak )
    {
        delete [] pcBuffer;
        return; 
    }
    
    for(ii=0; ii < cnNumberOfCharsBeforeBreak - cnNumberOfDots; ii++)
    {
        pcBuffer[ii] = GetAt(ii);
    }

    for(;ii < cnNumberOfCharsBeforeBreak; ii++)
    {
        pcBuffer[ii] = '.';
    }

    for(ii=0; ii < cnMaxNumberOfFilePathChars - cnNumberOfCharsBeforeBreak; ii++)
    {
       pcBuffer[cnMaxNumberOfFilePathChars - 1 - ii] = GetAt(m_Length - 1 - ii); 
    }
    
    assign(pcBuffer);

    delete [] pcBuffer;
}
//////////////////////////////////////////////////////////////////////////////////
int Stringreg::nGetPatternStart() { 
    return ((strreg_locs[0]==-1)?first_pat_loc:min(first_pat_loc,strreg_locs[0])); 
};

int Stringreg::nParse(const Cstring& sString,const Cstring& sPat) {
    int px,sx;
    int posstack[400];
    int posptr = -1;
    char* pat = sPat.string();
    char* str = sString.string();
    
    strreg_locs[0]=-1;  // This flag indicates that a * or ? has not yet been initialized.
    strreg_loc=0;
    strreg_str=sString; // Saved for printing purposes.
    strreg_nums_ct=0;
    strwht_locs[0] = -1;
    strwht_loc = 0;
    
    // RB increasing the buffer size
    char        patbuf[16];

    memset(patbuf, 0, 16);
    
    patbuf[0] = ' ';
    
    for(px=0,sx=0; 1;)
    {
        if (strreg_loc>=10) 
            return 1;
        if (pat[px]==0)
        { 
            if (str[sx]==0) 
                return 0; 
            else 
                goto retry;
        }
        if (pat[px] != '*')
        {
            if (str[sx]==0)
            {
                if (pat[px]==0) 
                    return 0; 
                else 
                    goto retry;
            }
        }
        if (pat[px] == '?')
        {
            strreg_locs[strreg_loc]=sx;
            strreg_lens[strreg_loc]=1;
            strreg_locs[strreg_loc+1]=-1;
            strreg_loc++;px++;  sx++;
            continue; 
        }
        if (pat[px] == '*')
        {
            if (strreg_locs[strreg_loc] == -1)
                strreg_locs[strreg_loc]=sx;
            if ((pat[px+1]) && (strchr("*?",pat[px+1]))) 
                return 1;
            
            // Expedite a '*'
            
            if (!pat[px+1])
            { 
                strreg_lens[strreg_loc] = strlen(&str[sx]);
                strreg_loc++;
                return 0; 
            }
            
            //      char* patbuf=" ";
            if (   (pat[px+1]=='.') 
                && (pat[px+2]) 
                && (pat[px+2]!='?') 
                && (pat[px+2]!='*'))
            {
                // Since whitespace is optional, we need to look-ahead beyond it possibly.
                //patbuf="\t\n ";
                strcpy(patbuf, "\t\n ");
                if (pat[px+2] == '#')
                    //patbuf = "\t\n 0123456789";
                    strcpy(patbuf, "\t\n 0123456789-+");
                else
                    if (pat[px+2] == '\\') 
                        patbuf[3]=pat[px+3];
                    else
                        patbuf[3]=pat[px+2];
            } 
            else
                if (pat[px+1] == '.')
                    //patbuf = "\t\n "; 
                    strcpy(patbuf, "\t\n "); 
                else
                    if (pat[px+1] == '\\')
                        patbuf[0]=pat[px+2];
                    else
                        if (pat[px+1] == '#')
                            //patbuf="0123456789.";
                            strcpy(patbuf, "0123456789.-+");
                        else 
                            patbuf[0] = pat[px+1];
                        if ((str[sx]) && strchr(patbuf,str[sx]))
                        {
                            posstack[++posptr] = px;
                            posstack[++posptr] = sx+1;
                            posstack[++posptr] = strreg_loc;
                            strreg_lens[strreg_loc]=sx-strreg_locs[strreg_loc];
                            strreg_locs[strreg_loc+1]=-1;
                            strreg_loc++;px++;
                            continue;
                        };
                        if (str[sx]==0) goto retry;
                        sx++;
                        continue;
        }
        
        if (pat[px]=='#') {
            char numbuf[20]; int x;
            for (x=0;(str[sx+x]) && (strchr("0123456789.-+",str[sx+x]));x++);
            if ((x==0) || ((x==1) && (str[sx]=='.'))) goto retry;
            if (x>=20) goto retry;
            if (strreg_nums_ct==10) 
                return 1;
            strncpy(numbuf,&str[sx],x); numbuf[x]=0;
            if (strchr(numbuf,'.')) {
                if (1!=sscanf(numbuf,"%lf",&strreg_nums[strreg_nums_ct].d)) goto retry;
                strreg_nums[strreg_nums_ct].type=strreg_num::DOUBLE;
            } else {
                if (1!=sscanf(numbuf,"%ld",&strreg_nums[strreg_nums_ct].l)) goto retry;
                strreg_nums[strreg_nums_ct].type=strreg_num::LONG;
            };
            strreg_nums_ct++;
            sx+=x;
            px++;
            continue;
        };
        if (pat[px]=='.') {
            int x;
            for (x=0;(str[sx+x]) && (strchr("\t\n\r ",str[sx+x]));x++);
            strwht_locs[strwht_loc] = sx;
            strwht_lens[strwht_loc] = x;
            strwht_loc++;
            sx+=x;
            px++;
            continue;
        };
        if (pat[px]=='\\') 
            px++;
        if (pat[px]==str[sx]) { 
            if (px==0)
                first_pat_loc = sx;
            px++; 
            sx++; 
            continue; 
        };
retry:
        if (posptr == -1) 
        {
            return 1;
        }
        strreg_loc=posstack[posptr--];
        sx=posstack[posptr--];
        px=posstack[posptr--];
    }; // for   
    
    // Look like it never reaches here;
    
    return 0;
};



int Stringreg::nPrint(Cstring& sFormat,Cstring& sBuf) {
    int x,y,z,s;
    char numbuf[50];
    char* fmt;
    int nDefaultArg = 0;
    int nDefaultWht = 0;
    
    sBuf = "";
    fmt = sFormat.string();
    s=strlen(fmt);
    for (x=0,y=0;x<s;x++) {
        if ((fmt[x]=='*') || (fmt[x]=='?')) {
            for (z=0;z<strreg_lens[nDefaultArg];z++) 
                sBuf+=strreg_str.GetAt(z+strreg_locs[nDefaultArg]);
            nDefaultArg++;
        } else if (fmt[x]=='%') {
            x++;
            if (fmt[x]=='.') {
                for (z=0;z<strwht_lens[nDefaultWht];z++) 
                    sBuf+=strreg_str.GetAt(z+strwht_locs[nDefaultWht]);
                nDefaultWht++;
            } else if (fmt[x]=='#') {
                x++;
                if ((!isdigit(fmt[x])) || (strreg_nums_ct<=fmt[x]-'0')) 
                    return 1;
                z=fmt[x]-'0';
                if (strreg_nums[z].type==strreg_num::LONG) sprintf(numbuf,"%ld",strreg_nums[z].l); else
                    sprintf(numbuf,"%lf",strreg_nums[z].d);
                for (z=0;z<strlen(numbuf);) 
                    sBuf+=numbuf[z++];
                continue;
            } else {
                if (!isdigit(fmt[x])) 
                    return 1;
                if (fmt[x]=='0')
                    sBuf+=strreg_str;
                else if (fmt[x]-'1'>=strreg_loc) 
                    return 1;
                else for (z=0;z<strreg_lens[fmt[x]-'1'];z++) 
                    sBuf+=strreg_str.GetAt(z+strreg_locs[fmt[x]-'1']);
            };
        } else sBuf+=fmt[x];
    };
    return 0;
};

int Stringreg::nLoadString(int varnum,Cstring& sBuf) {
    int y,z;
    if (varnum>=strreg_loc) 
        return 1;
    for (z=0,y=0;z<strreg_lens[varnum];z++) sBuf+=strreg_str.GetAt(z+strreg_locs[varnum]);
    return 0;
};

int Stringreg::nLoadInt(int varnum,int& nValue) {
    if (strreg_nums_ct<=varnum) 
        return 1;
    if (strreg_nums[varnum].type==strreg_num::LONG)
        nValue =(strreg_nums[varnum].l);
    else if (strreg_nums[varnum].type==strreg_num::DOUBLE)
        nValue = (int) (strreg_nums[varnum].d);
    return 0;
};
int Stringreg::nLoadFloat(int varnum,float& fValue) {
    if (strreg_nums_ct<=varnum) 
        return 1;
    if (strreg_nums[varnum].type==strreg_num::LONG)
        fValue = (strreg_nums[varnum].l);
    else if (strreg_nums[varnum].type==strreg_num::DOUBLE)
        fValue = (strreg_nums[varnum].d);
    return 0;
};

