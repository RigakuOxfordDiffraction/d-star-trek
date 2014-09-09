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
// Crefln.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Crefln which implements
//    the reflection encapsulation of d*TREK.
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
//   Reflections consist of an array of at least 3 integer values, an array of
//   of least 2 float values, possibly an array of Cstring values and a pointer
//   to a reflection list which holds group properties of the reflection.
//   All reflections in a list must have arrays with the same number of
//   elements. The lengths of the arrays are dynamically controlled.  The first
//   the elements of the integer array contain the Miller indices H, K, L.
//   The subsequent elements, if any, are defined by the Creflnlist class.
//   The first 2 elements of the float array contain the intensity and 
//   standard deviation of the intensity.
//
//   This class works closely with the Creflnlist class.
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Crefln.h"         // Class definition and prototypes
#include "string.h"
#include "memory.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

int Crefln::ms_nErrorSpecial      =  0;   /*A*/
int Crefln::ms_nErrorOnEdge0      =  1;   /*B*/
int Crefln::ms_nErrorOffEdge0     =  2;   /*C*/
int Crefln::ms_nErrorOnEdge1      =  3;   /*D*/
int Crefln::ms_nErrorOffEdge1     =  4;   /*E*/
int Crefln::ms_nErrorOnEdge2      =  5;   /*F*/
int Crefln::ms_nErrorOffEdge2     =  6;   /*G*/
int Crefln::ms_nErrorHKLBound     =  7;   /*H*/
int Crefln::ms_nErrorRing         =  8;   /*I*/
int Crefln::ms_nErrorTooDark      =  9;   /*J*/
int Crefln::ms_nErrorBackground   = 10;   /*K*/
int Crefln::ms_nErrorBackgroundSD = 11;   /*L*/
int Crefln::ms_nErrorNonunfA      = 12;   /*M*/
int Crefln::ms_nErrorNonunfB      = 13;   /*N*/
int Crefln::ms_nErrorPartialStart = 14;   /*O*/
int Crefln::ms_nErrorPartialEnd   = 15;   /*P*/
int Crefln::ms_nErrorRotTooWide   = 16;   /*Q*/
int Crefln::ms_nErrorIntProblem   = 17;   /*R*/
int Crefln::ms_nErrorTooFewPix    = 18;   /*S*/
int Crefln::ms_nErrorUnaccountedPix = 19; /*T*/
int Crefln::ms_nErrorOverlap        = 20; /*U*/


Cstring Crefln::ms_asErrorsNames[] =  {
"Special",
"OnEdge0",
"OffEdge0",
"OnEdge1",  
"OffEdge1",  
"OnEdge2",  
"OffEdge2",     
"HKLBound",    
"Ring",     
"TooDark",
"Background",
"BackgroundSD",
"NonunfA",
"NonunfB",     
"PartialStart",
"PartialEnd",
"RotTooWide",  
"IntProblem",  
"TooFewPix",  
"UnaccountedPix",
"Overlap",
"",
"",
"",
"",
"",
"",
"",
"",
"",
"",
""
};

Cstring Crefln::ms_asErrorStrings[] = {
  "Processed special",
  "On edge 0",
  "Off edge 0",
  "On edge 1",
  "Off edge 1",
  "On edge 2",
  "Off edge 2",
  "Outside HKL bounds",
  "Ring",
  "Too dark",                 
  "Bad background",           
  "Bad background sd",
  "Bad non-uniformity A",
  "Bad non-uniformity B",
  "Partial at scan start",
  "Partial at scan end",
  "Rotation too wide",
  "Integration problem",
  "Too few pixels",
  "Fuzzy spot shape",
  "Overlap",
  "?",
  "?",
  "?",
  "?",
  "?",
  "?"
};

Cstring Crefln::ms_asErrorPhrases[] = {
  "Not processed as strong, but otherwise OK",
  "Some peak on box edge in 1st dim",
  "Some peak outside box edge in 1st dim",
  "Some peak on box edge in 2nd dim",
  "Some peak outside box edge in 2nd dim",
  "Some peak on box edge in 3rd dim",
  "Some peak outside box edge in 3rd dim",
  "Mask lies outside HKL bounds",
  "Intersection with excluded resolution ring",
  "Some peak pixels saturated",   
  "Error in background determination", 
  "Error in background sigma determination",
  "Any shoebox pixels flagged as bad",
  "Peak pixels flagged as bad or in shadow",
  "Reflns incomplete at start of the rotation",
  "Reflns incomplete at end of the rotation",
  "Reflns predicted to be on too many images",
  "Bad mask computation",
  "Too few pixels in peak mask",
  "Fuzzy spot profile",
  "Spatial overlaps (dtprofit)",
  "?",
  "?",
  "?",
  "?",
  "?",
  "?"
};



// Constructors, destructors and assignments


Crefln::Crefln(Creflnlist* poReflnlistIn, 
               const int nHval, const int nKval, const int nLval)
{
  m_poReflnlist = poReflnlistIn;
  m_poData = NULL;
  m_nSelect = 0;
  (void) nInitValues();
  vSetH(nHval);
  vSetK(nKval);
  vSetL(nLval);
}

Crefln::Crefln(const Crefln& oOther)
{
  m_poData = NULL;
  m_a2cOffsets[0] = oOther.m_a2cOffsets[0];
  m_a2cOffsets[1] = oOther.m_a2cOffsets[1];
  m_nSelect = 0;

  // Now copy the different Crefln fields to the new Crefln

  m_eThe_State  = oOther.m_eThe_State;
  m_poReflnlist = oOther.m_poReflnlist;
  (void) nInitValues();
  memcpy(&m_poData->m_pnIntField[0],
         &oOther.m_poData->m_pnIntField[0],
         m_nTotalDataSize);
  vCheckStrings();

  m_pvUserData = oOther.m_pvUserData;
}

Crefln::Crefln(Creflnlist& oReflnlist, const Crefln& oOther)
{
  if (&oReflnlist == oOther.m_poReflnlist)
    cout << "Warning in Crefln:: FROM and TO Creflnlists are the same!\n";

  m_poData = NULL;
  m_a2cOffsets[0] = oOther.m_a2cOffsets[0];
  m_a2cOffsets[1] = oOther.m_a2cOffsets[1];
  m_nSelect = 0;

  // Now copy the different Crefln fields to the new Crefln

  m_eThe_State  = oOther.m_eThe_State;
  m_poReflnlist = &oReflnlist;
  m_pvUserData = oOther.m_pvUserData;
  (void) nInitValues();
  memcpy(&m_poData->m_pnIntField[0],
         &oOther.m_poData->m_pnIntField[0],
         m_nTotalDataSize);
  vCheckStrings();
}

inline void Crefln::vCheckStrings() { 
    int i;   
    if (m_poReflnlist) {
        for(i=0;i<m_poReflnlist->m_nCstringReflnFields;i++) {
            if (pcGetField(i)[0]==0)
                vSetField(i,m_poReflnlist->ms_ssIsEmpty);
        };
    };
};


Crefln& Crefln::operator=(const Crefln& oOther)
{


  if (&oOther == this)     // these 2 statements prevent ...
    return *this;          //   ... assigning an object to itself

  int i;

  // Copy the other Crefln state, int, float and Cstring Fields to this Crefln
  // Do not copy the pointer poReflnlist!
  // Set m_pcSelect to all 0's!
  // Only copy the Fields we have space for  (do we need warning message
  //                                         if Crefln is truncated?)
  // There is a potential bug here, because the Field categories may not match.
  // to avoid this bug, check that the Creflns belong to the same list (thus,
  // the fields must match.
  // 

  if (m_poReflnlist == oOther.m_poReflnlist) 
    {
      m_eThe_State  = oOther.m_eThe_State;
      m_nSelect = 0;
      memset(&m_poData->m_pnIntField[0], 
             0,
             m_nTotalDataSize);
      if ((oOther.m_a2cOffsets[0]==m_a2cOffsets[0]) && (oOther.m_a2cOffsets[1]==m_a2cOffsets[1]) && (m_nTotalDataSize == oOther.m_nTotalDataSize)) {
        memcpy(&m_poData->m_pnIntField[0],
               &oOther.m_poData->m_pnIntField[0],
               m_nTotalDataSize);

      } else {
          // Oh no! we are trying to fool the reflection as to how the data is stored.
          // This requires three distinct copies.
          memcpy(&m_poData->m_pnIntField[0],
                 &oOther.m_poData->m_pnIntField[0],
                 sizeof(int)* m_poReflnlist->m_nIntReflnFields);
          memcpy(&m_poData->m_pfFloatField[m_a2cOffsets[0]],
                 &oOther.m_poData->m_pfFloatField[oOther.m_a2cOffsets[0]],
                 sizeof(float)* m_poReflnlist->m_nFloatReflnFields);
          if (m_poReflnlist->m_nCstringReflnFields)
              memcpy(&m_poData->m_pcStringField[m_a2cOffsets[1]],
                     &oOther.m_poData->m_pcStringField[oOther.m_a2cOffsets[1]],
                     nCREFLN_STRING_SIZE* m_poReflnlist->m_nCstringReflnFields);
      };
      vCheckStrings();
  } else  {
      int     nTempOther;
      Cstring sTemp;
      
      for (i = 0; i < m_poReflnlist->m_nIntReflnFields; i++) 
      {
          
          // Get name of each THIS integer field and see if it exists in OTHER
          // list
          
          sTemp      = m_poReflnlist->sGetFieldName(i, eReflnField_int_type);
          nTempOther = oOther.m_poReflnlist->nGetFieldIndex(sTemp, 
              eReflnField_int_type);
          if (0 <= nTempOther)
          {
              // There is a match
              
              vSetField(i, oOther.nGetField(nTempOther));
          }
      }
      //
      for (i = 0; i < m_poReflnlist->m_nFloatReflnFields; i++) 
      {
          
          // Get name of each THIS float field and see if it exists in OTHER
          // list
          
          sTemp      = m_poReflnlist->sGetFieldName(i, eReflnField_float_type);
          nTempOther = oOther.m_poReflnlist->nGetFieldIndex(sTemp, 
              eReflnField_float_type);
          if (0 <= nTempOther)
          {
              // There is a match
              
              vSetField(i, oOther.fGetField(nTempOther));
          }
      }
      //
      for (i = 0; i < m_poReflnlist->m_nCstringReflnFields; i++) 
      {
          
          // Get name of each THIS Cstring field and see if it exists in OTHER
          // list
          
          sTemp      = m_poReflnlist->sGetFieldName(i, eReflnField_Cstring_type);
          nTempOther = oOther.m_poReflnlist->nGetFieldIndex(sTemp, 
              eReflnField_Cstring_type);
          if (0 <= nTempOther)
          {
              // There is a match
              
              vSetField(i, oOther.sGetField(nTempOther));
          }
          else if ("" == sGetField(i))
          {  // Do not allow empty string fields
              vSetField(i, m_poReflnlist->ms_ssIsEmpty);
          }
      }
  }

  if (NULL != oOther.m_pvUserData)
    {
      m_pvUserData = oOther.m_pvUserData;
    }

  return (*this);
}

Crefln::~Crefln()
{
  // Higher level procedure is responsible for any memory leak
  // with m_pvUserData

  m_pvUserData  = NULL;
  if (m_poData)
      delete[] ((float*) m_poData);
  m_poData = NULL;
}

int Crefln::nInitValues(void)
{
//  cout << "Crefln::nInitValues called\n";

  m_eThe_State  = eRefln_unknown_state;
  if (NULL == m_poReflnlist)
    {

      nComputeOffsets(3,2,0);
      m_nSelect = 1 + (m_nTotalDataSize /sizeof(float));
      m_poData = (CreflnData*) new float [m_nSelect];
      //      m_poData = (CreflnData*) new char[m_nTotalDataSize];
      //+jwp 19-Jun-2001
      memset(m_poData, 0, m_nSelect * sizeof(float));
      //-jwp 19-Jun-2001
      m_nSelect = 0;
      vSetField(0,0);
      vSetField(1,0);
      vSetField(2,0);
      vSetField(0,(float) -999.0);
      vSetField(1,(float) -999.0);      
    }
  else
    {
      int i;

      // Allocate memory for the different Crefln fields
      
      nComputeOffsets(m_poReflnlist->m_nIntReflnFields,m_poReflnlist->m_nFloatReflnFields,m_poReflnlist->m_nCstringReflnFields);
      i = 1 + (m_nTotalDataSize / sizeof(float));
      m_poData = (CreflnData*) new float [i];
      //      m_poData = (CreflnData*) new char[m_nTotalDataSize];
      //+jwp 19-Jun-2001
      memset(m_poData, 0, i * sizeof(float));
      //-jwp 19-Jun-2001

      for (i = 0; i < m_poReflnlist->m_nIntReflnFields; i++) 
          vSetField(i,0);
      for (i = 0; i < m_poReflnlist->m_nFloatReflnFields; i++) 
          vSetField(i,(float) -999.0);
      for (i = 0; i < m_poReflnlist->m_nCstringReflnFields; i++) 
          vSetField(i,"");
      m_nSelect = 0;
    }
  m_pvUserData      = NULL;
  return (0);
}

int Crefln::nList(const int nType)
{
  int i;
  if (0 == nType)
    (void) m_poReflnlist->nListFields();
  else if (1 == nType) 
    {
      cout << "HKL: " << nGetH() << " " << nGetK() << " " << nGetL() << endl;
    }
  else if (2 == nType)
    {
      for (i = 0; i < m_poReflnlist->m_nIntReflnFields; i++)
        {
          cout << m_poReflnlist->sGetFieldName(i, eReflnField_int_type)
               << ": " << nGetField(i) << endl;
        }
      for (i = 0; i < m_poReflnlist->m_nFloatReflnFields; i++)
        {
          cout << m_poReflnlist->sGetFieldName(i, eReflnField_float_type)
               << ": " << fGetField(i) << endl;
        }
      for (i = 0; i < m_poReflnlist->m_nCstringReflnFields; i++)
        {
          cout << m_poReflnlist->sGetFieldName(i, eReflnField_Cstring_type)
               << ": " << sGetField(i) << endl;
        }
    }
  return (0);
}

void Crefln::vSetField(const int nFieldIndex, const Cstring& sVal) {
    strncpy(&m_poData->m_pcStringField[m_a2cOffsets[1]+nFieldIndex][0],
            sVal.string(),
            nCREFLN_STRING_SIZE);
    m_poData->m_pcStringField[m_a2cOffsets[1]+nFieldIndex][nCREFLN_STRING_SIZE-1] = 0;
}

void Crefln::vSetField(const int nFieldIndex, const char* pcVal) {
    strncpy(&m_poData->m_pcStringField[m_a2cOffsets[1]+nFieldIndex][0],
            pcVal,
            nCREFLN_STRING_SIZE);
    m_poData->m_pcStringField[m_a2cOffsets[1]+nFieldIndex][nCREFLN_STRING_SIZE-1] = 0;
}


int Crefln::nComputeOffsets(int nInts,int nFloats,int nStrings) {
    int nx;
    nx = sizeof(int)*nInts;
    m_a2cOffsets[0] = (nx)/sizeof(float) + ((nx % sizeof(float))!=0);
    nx = sizeof(int)*nInts + sizeof(float)*nFloats;
    m_a2cOffsets[1] = (nx)/nCREFLN_STRING_SIZE + ((nx % nCREFLN_STRING_SIZE)!=0);

    m_nTotalDataSize = m_a2cOffsets[1]*nCREFLN_STRING_SIZE + nStrings*nCREFLN_STRING_SIZE;
    return 0;
};


int
Crefln::nUpdateHeader(Cimage_header *poHeader, const Cstring& sPre)
{
  // Place the reflection information in an image header.
  // This may be useful for debugging or for placing standard
  // reflection measurements in the header.
  // The header contains keywords and values.  The keywords of
  // a reflection are the reflnlist field names with the prefix
  // sPre prepended to them.
  //
  // *poHeader   -- the header object to update.
  // sPre        -- the prefix string before each keyword, may be blank.
  // Returns     -- 0 success, not 0 for failure

  int i;         // Loop counter
  int nStat;     // Local status

  nStat = 0;
  for (i = 0; i < m_poReflnlist->m_nIntReflnFields && (0 == nStat); i++)
    {
      nStat = poHeader->nReplaceValue(sPre
                      + m_poReflnlist->sGetFieldName(i, eReflnField_int_type), 
                                      nGetField(i));
    }
  for (i = 0; i < m_poReflnlist->m_nFloatReflnFields && (0 == nStat); i++)
    {
      nStat = poHeader->nReplaceValue(sPre
                      + m_poReflnlist->sGetFieldName(i, eReflnField_float_type),
                                      fGetField(i));
    }
  for (i = 0; i < m_poReflnlist->m_nCstringReflnFields && (0 == nStat); i++)
    {
      nStat = poHeader->nReplaceValue(sPre
               + m_poReflnlist->sGetFieldName(i, eReflnField_Cstring_type),
                                      sGetField(i));
    }
  return (nStat);
}

void
Crefln::vWrite(FILE* pFOut, const char *pcSelectIn)
{
    // Write a reflection to the specified output stream
    // Status of the write is returned in *pFOut
    
    int j;
    float f0;
    int nTemp;
    char* pcTemp;
    int nChars = 0;
    int nWarn  = 0;
    int nTotal = 0;
    int nPass  = 0;
    static char a2048cTemp[2048];
    static char a40cTemp1[40];
    static char a40cTemp2[40];
    
    // If we are writing binary reflection lists, then use a sepearate routine.
    if (m_poReflnlist->m_poBinaryHeader != NULL)
    {
        m_poReflnlist->m_poBinaryHeader->nWrite(pFOut,*this,pcSelectIn);
        return;
    }
    
    // Write out only the selected fields, but only if the reflection
    // is 'selected' or 'not ignored'
    
    const char *pcIn = pcSelectIn;
    const int  *pnIn = m_poReflnlist->m_pnFortranFormat;
    const int  *pnIn2 = m_poReflnlist->m_pnFortranFormat2;
    nPass = 0;
    while (nPass < 2)
    {
        //          nPass++;
        nPass = 2;
        nTotal = 0;
        if ((!pcIn) || (m_nSelect == pcIn[m_poReflnlist->m_nTotalFieldsPlus1-1]))
        {
            for (j = 0; j < m_poReflnlist->m_nIntReflnFields; j++)
            {

                if ((!pcIn) || (0 == *pcIn++))
                {
                    if (pnIn) {

                        sprintf(&a40cTemp1[0],"%%%dd",(*pnIn));

                        if ((*pnIn) != sprintf(&a2048cTemp[nTotal],a40cTemp1,nGetField(j))) {
                            cout << "WARNING in Crefln:  could not format integer " << nGetField(j) << "\n" 
                                << " for formated output.\n";

                        };
                        nChars = *(pnIn);
                    } else if (j < 3)
                        nChars = sprintf(&a2048cTemp[nTotal]," %5d",nGetField(j));
                    else
                        nChars = sprintf(&a2048cTemp[nTotal]," %d",nGetField(j));
                    if (0 >= nChars)
                    {
                        cout << "WARNING in Crefln: error transferring int " << j << ": "
                            << nGetField(j) << endl;
                        nWarn++;
                        nChars = 0;
                    }
                    nTotal += nChars;
                    a2048cTemp[nTotal] = 0;
                }
                if (pnIn) 
                    pnIn++;
                if (pnIn2)
                    pnIn2++;
            }
            for (j = 0; j < m_poReflnlist->m_nFloatReflnFields; j++) 
            {

                if ((!pcIn) || (0 == *pcIn++))
                {
                    f0 = fGetField(j);
                    if (pnIn) {
                        int nDecimal;
                        nTemp = sprintf(&a40cTemp1[0],"%-38.8f",f0);
                        for (nDecimal = 0; nDecimal < nTemp; nDecimal++) {
                            if (a40cTemp1[nDecimal] == '.')
                                break;
                        };
                        for (nTemp = 0; nTemp < (*pnIn); nTemp++)
                            a2048cTemp[nTotal + nTemp] = ' ';
                        
                        if (nDecimal < (*pnIn)) {
                            
                            if ((*pnIn2)==0)
                                nDecimal--;
                            else {
                                int nx;
                                for (nx = 0; nx < (*pnIn2); nx++) {
                                    if (nDecimal + 1 < (*pnIn))
                                        nDecimal++;
                                };
                            };

                            for (nTemp = (*pnIn) - 1;(nDecimal >= 0); nDecimal--,nTemp--)
                                a2048cTemp[nTotal + nTemp] = a40cTemp1[nDecimal];
                            nTotal += (*pnIn);
                            a2048cTemp[nTotal] = 0;
                        } else {
                            cout << "WARNING in Crefln:  could not format float " << fGetField(j) << "\n" 
                                << " for formated output.\n";
                        };
                    } else {
                        nTemp = (int) ABS(f0);
                        if (ABS(f0) >= 1000000)
                            nTemp = 1;
                        else if (nTemp < 1)
                            nTemp = 6;
                        else if (nTemp < 10)
                            nTemp = 5;
                        else if (nTemp < 100)
                            nTemp = 4;
                        else if (nTemp < 1000)
                            nTemp = 3;
                        else if (nTemp < 10000)
                            nTemp = 2;
                        else
                            nTemp = 1;
                        nChars = sprintf(&a2048cTemp[nTotal]," %.*f", nTemp, f0);
                        if (0 >= nChars)
                        {
                            cout << "WARNING in Crefln: error transferring float " << f0 << ": "
                                << nGetField(j) << endl;
                            nWarn++;
                            nChars = 0;
                        }
                        nTotal += nChars;
                    };
                }
                if (pnIn) 
                    pnIn++;
                if (pnIn2)
                    pnIn2++;
            }
            // We should have ALL numbers up to this point in a2048cTemp, 
            // so anything but [ -+.0-9] indicates an error
            
            j = 0;
            pcTemp = a2048cTemp;
            while (j < nTotal)
            {
                // Assume ASCII collating sequence, ignore /,
                if ( (*pcTemp != ' ') && ((*pcTemp < '+') || (*pcTemp > '9')) )
                {
                    cout << "WARNING Removing bogus characters in numeric fields: " << nPass << "\n"
                            << a2048cTemp << endl;
                    // This will replace it with -9[9..], unless
                    // only a single character is wrong.
                    
                    char cSubs = '-';
                    do
                    {
                        a2048cTemp[j] = cSubs;
                        if ('-' == cSubs)
                            cSubs = '9';
                        else
                            cSubs = ' ';
                        j++;
                        pcTemp++;
                    } while ( (j < nTotal) && (a2048cTemp[j] != ' ') );
                    
                    if (2 >= nPass)
                    {
                        cout << "REPAIRED: " << nPass << "\n" << a2048cTemp << endl;
                        
                        // Try to make fSigmaI < 0, so it will be excluded later;
                        
                        vSetSigmaI(-ABS(fGetSigmaI()));
                    }
                }
                j++;
                pcTemp++;
            }
            
            for (j = 0; j < m_poReflnlist->m_nCstringReflnFields; j++)
            {

                if ((!pcIn) || (0 == *pcIn++))
                {
                    if (pnIn) {
                        
                        sprintf(&a40cTemp1[0],"%%%ds",(*pnIn));
                        
                        if ((*pnIn) != sprintf(&a2048cTemp[nTotal],a40cTemp1,sGetField(j).string())) {
                            cout << "WARNING in Crefln:  could not format string " << sGetField(j) << "\n" 
                                << " for formated output.\n";
                            
                        };
                        nChars = *(pnIn);
                        nTotal += nChars;
                        a2048cTemp[nTotal] = 0;
                    } else {
                        
                        nChars = sprintf(&a2048cTemp[nTotal]," %s", sGetField(j).string());
                        if (0 >= nChars)
                        {
                            cout << "WARNING in Crefln: error transferring string " << j << ": "
                                << sGetField(j) << endl;
                            nWarn++;
                            nChars = 0;
                        }
                        nTotal += nChars;
                        a2048cTemp[nTotal] = 0;
                    };
                }
                if (pnIn)
                    pnIn++;
                if (pnIn2)
                    pnIn2++;
            }
            
            
            if ( (0 == nWarn) && (nTotal < 2048) )
            {
                nTemp = fprintf(pFOut,"%s\n", a2048cTemp);
                if (nTemp != nTotal+1)
                    cout << "WARNING in Crefln: error writing refln:\n" 
                    << a2048cTemp << endl << flush;
            }   
            else
                cout << "WARNING in Crefln: refln not written because of previous warning(s).\n" << flush;            
        }
   }
   
}

int
Crefln::nRead(FILE* pFIn)
{
  // Read a line that contains a reflection from the input stream
  // specified by *poIn
  // Returns an error status, 0 if success otherwise
  // the tens digit says error occurred reading an integer field
  // the thousands digits says an error occurred reading a float field
  // the hundred thousands digits says an error occurred reading a
  //  string field

  int     i;      // Loop counter
  int     nTemp;
  float   fTemp;
  int     nStat = 0;

  // Read all the integer fields

  for (i = 0; (i < m_poReflnlist->m_nIntReflnFields) && (0 == nStat); i++) 
  {
      if (1==fscanf(pFIn,"%d",&nTemp))
          vSetField(i, nTemp);
      else
          nStat = (i+1) * 10;
  }

  // Read all the float fields

  for (i = 0; (i < m_poReflnlist->m_nFloatReflnFields) && (0 == nStat); i++) 
  {
      if (1==fscanf(pFIn,"%f",&fTemp))
          vSetField(i, fTemp);
      else
          nStat = (i+1) * 1000;
  }

  // Read all the Cstring fields

  for (i = 0; (i < m_poReflnlist->m_nCstringReflnFields) && (0 == nStat); i++) 
    {
      char pcBuffer[200];

      if (1==fscanf(pFIn,"%s",pcBuffer))
          vSetField(i, pcBuffer);
      else
          nStat = (i+1) * 100000;
    }

  return (nStat);
}

int
Crefln::nPackHKL(const int *pnHKL) const
{
  // Pack Miller indices into a single 32-bit signed integer for sorting:
  // Bit   0    is unused in this routine, but holds +- anom flag: 0 = +, 1 = -
  //                                       in other routines.
  // Bits  1-10 hold L + 512;     -511 <= L <= 511;
  // Bits 11-20 hold K + 512;     -511 <= K <= 511;
  // Bits 21-30 hold H + 512;     -511 <= H <= 511;
  // All bits are set if any of the input indices are out of range, that is
  //             the result is -1.
  // Although the packing algorithm will accommodate an index one less than
  // the above limits, do not allow that since the nReduceHKL routine will 
  // try to change the sign anyways.

  const int *pnTemp;
  if (NULL == pnHKL)
    {
      pnTemp = m_poData->m_pnIntField;
    }
  else
    {
      pnTemp = pnHKL;
    }
  if ( (pnTemp[0] < -511) || (pnTemp[0] > 511) ||
       (pnTemp[1] < -511) || (pnTemp[1] > 511) ||
       (pnTemp[2] < -511) || (pnTemp[2] > 511) ) return (-1);
  return (  ((pnTemp[0] + 512) << 21)
          | ((pnTemp[1] + 512) << 11)
          |  (pnTemp[2] + 512) <<  1);
}

int
Crefln::nUnPackHKL(const int nPacked, int *pnHKL)
{
  // Unpack 3 Miller indices that were packed by nPack

  int *pnTemp;
  if (NULL == pnHKL)
    {
      pnTemp = m_poData->m_pnIntField;
    }
  else
    {
      pnTemp = pnHKL;
    }
  pnTemp[0] = ((nPacked & 0x7fe00000) >> 21) -  512;  //(1023L << 21) = 0x7fe00000
  pnTemp[1] = ((nPacked &   0x1ff800) >> 11) -  512;  //(1023L << 11) =   0x1ff800
  pnTemp[2] = ((nPacked &      0x7fe) >>  1) -  512;  //(1023L <<  1) =    0x7fe

  if (nPacked < 0) return (-1);
  return (0);
}

void Crefln::vReindex(const float *pfMatIn)
{
  // Reindex the hkl by multiplying the 3x3 matrix *pfMatIn.
  // No check is made to ensure that *pfMatIn is a valid reindexing matrix
  // (right-handedness, determinant==1, etc)

  int   i;
  float a3fHKLOld[3];
  float a3fHKLNew[3];

  for (i = 0; i < 3; i++)
    {
      a3fHKLOld[i] = (float) m_poData->m_pnIntField[i];
    }

  vMulMatNDVecND(3, pfMatIn, a3fHKLOld, a3fHKLNew);

  // Round-off indices to nearest integer hkl

  for (i = 0; i < 3; i++)
    {
      if (0 > a3fHKLNew[i])
        {
          m_poData->m_pnIntField[i] = (int) (a3fHKLNew[i] - 0.5);
        }
      else
        {
          m_poData->m_pnIntField[i] = (int) (a3fHKLNew[i] + 0.5);
        }
    }
}



bool Crefln::bEllipsoidAvailable(bool bCheckParamFieldsAlso) {
    return ((!bCheckParamFieldsAlso || ((m_poReflnlist->m_nFI_fEllipseMajorAxis>=0) && 
        (m_poReflnlist->m_nFI_fEllipseMinorAxis>=0) && 
        (m_poReflnlist->m_nFI_fEllipseAxisMajorOffset>=0)))  &&
        (m_poReflnlist->m_nFI_fEllipsoidA00>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidA01>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidA11>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidb0>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidb1>=0) &&        
        (m_poReflnlist->m_nFI_fEllipsoidc>=0) &&        
        ((m_poReflnlist->m_nFI_fObsPx0>=0) || (m_poReflnlist->m_nFI_fCalcPx0>=0)) &&
        ((m_poReflnlist->m_nFI_fObsPx1>=0) || (m_poReflnlist->m_nFI_fCalcPx1>=0)) );
};

int
Crefln::nConvertToAxialEllipse(double* pfMaxSize,double* pfMinSize) {
    bool bSetParamFields;

    bSetParamFields = (!pfMaxSize && !pfMinSize);
  
    if (bEllipsoidAvailable(bSetParamFields)) {
            
            
        double a7fEllipsoid[7];
        double a2fRelOffset[2];
        double a2fMajorMinor[2];
        double fRotOffset;
        double a2fObsPx[2];
        a7fEllipsoid[0] = fGetField(m_poReflnlist->m_nFI_fEllipsoidA00);
        a7fEllipsoid[1] = fGetField(m_poReflnlist->m_nFI_fEllipsoidA01);
        a7fEllipsoid[2] = fGetField(m_poReflnlist->m_nFI_fEllipsoidA01);
        a7fEllipsoid[3] = fGetField(m_poReflnlist->m_nFI_fEllipsoidA11);
        a7fEllipsoid[4] = fGetField(m_poReflnlist->m_nFI_fEllipsoidb0);
        a7fEllipsoid[5] = fGetField(m_poReflnlist->m_nFI_fEllipsoidb1);
        a7fEllipsoid[6] = fGetField(m_poReflnlist->m_nFI_fEllipsoidc);

        a2fObsPx[0] = fGetField(m_poReflnlist->m_nFI_fObsPx0);
        a2fObsPx[1] = fGetField(m_poReflnlist->m_nFI_fObsPx1);

        if (!nGetEllipsoidTiltAngleCenter(&a7fEllipsoid[0],&a7fEllipsoid[4],a7fEllipsoid[6],&a2fRelOffset[0],&a2fMajorMinor[0],&fRotOffset)) {
            
            if (bSetParamFields) {
                vSetField(m_poReflnlist->m_nFI_fEllipseMajorAxis,(float) a2fMajorMinor[0]);
                vSetField(m_poReflnlist->m_nFI_fEllipseMinorAxis,(float) a2fMajorMinor[1]);
                vSetField(m_poReflnlist->m_nFI_fEllipseAxisMajorOffset,(float) fRotOffset);
                vSetField(m_poReflnlist->m_nFI_fObsPx0,(float) (0.5 + a2fObsPx[0] + a2fRelOffset[0]));
                vSetField(m_poReflnlist->m_nFI_fObsPx1,(float) (0.5 + a2fObsPx[1] + a2fRelOffset[1]));
            } else {
                if (pfMaxSize)
                    *pfMaxSize = a2fMajorMinor[0]*2.0;
                if (pfMinSize)
                    *pfMinSize = a2fMajorMinor[1]*2.0;
            };
        } else
            return 1;
        return 0;
    };
    return 1;   
};

int 
Crefln::nPutGetEllipse(bool bPutToRefln,float* a2x2fEllipseA,float* a2fEllipseb,float* pfEllipsec,float* a2fOffset) {
    double a2x2fEllipseA_[2][2];
    double a2fEllipseb_[2];
    double fEllipsec_;
    double a2fOffset_[2];

    if (bPutToRefln) {
        vCopyVecND(4,a2x2fEllipseA,&a2x2fEllipseA_[0][0]);
        vCopyVecND(2,a2fEllipseb,&a2fEllipseb_[0]);
        fEllipsec_ = *pfEllipsec;
        if (a2fOffset)
            vCopyVecND(2,a2fOffset,&a2fOffset_[0]);
    };
    if (nPutGetEllipse(bPutToRefln,&a2x2fEllipseA_[0][0],&a2fEllipseb_[0],&fEllipsec_,a2fOffset?(&a2fOffset_[0]):NULL))
        return 1;
    if (!bPutToRefln) {
        vCopyVecND(4,&a2x2fEllipseA_[0][0],a2x2fEllipseA);
        vCopyVecND(2,&a2fEllipseb_[0],a2fEllipseb);
        *pfEllipsec = fEllipsec_;
        if (a2fOffset)
            vCopyVecND(2,&a2fOffset_[0],a2fOffset);
    };
    return 0;
};

int
Crefln::nPutGetEllipse(bool bPutToRefln,double* a2x2fEllipseA,double* a2fEllipseb,double* pfEllipsec,double* a2fOffset) {
    if ((m_poReflnlist->m_nFI_fEllipsoidA00>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidA01>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidA11>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidb0>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidb1>=0) &&
        (m_poReflnlist->m_nFI_fEllipsoidc>=0)) {
        if (bPutToRefln) {
            vSetField(m_poReflnlist->m_nFI_fEllipsoidA00,(float) a2x2fEllipseA[0]);
            vSetField(m_poReflnlist->m_nFI_fEllipsoidA01,(float) ((a2x2fEllipseA[1] + a2x2fEllipseA[2])*0.5));
            vSetField(m_poReflnlist->m_nFI_fEllipsoidA11,(float) a2x2fEllipseA[3]);
            vSetField(m_poReflnlist->m_nFI_fEllipsoidb0,(float) a2fEllipseb[0]);
            vSetField(m_poReflnlist->m_nFI_fEllipsoidb1,(float) a2fEllipseb[1]);
            vSetField(m_poReflnlist->m_nFI_fEllipsoidc,(float) *pfEllipsec);
            if ((a2fOffset) && (m_poReflnlist->m_nFI_fObsPx0>=0) && (m_poReflnlist->m_nFI_fObsPx1>=0)) {
                vSetField(m_poReflnlist->m_nFI_fObsPx0,(float) a2fOffset[0]);
                vSetField(m_poReflnlist->m_nFI_fObsPx1,(float) a2fOffset[1]);
            } else if ((a2fOffset) && (m_poReflnlist->m_nFI_fCalcPx0>=0) && (m_poReflnlist->m_nFI_fCalcPx1>=0)) {
                vSetField(m_poReflnlist->m_nFI_fCalcPx0,(float) a2fOffset[0]);
                vSetField(m_poReflnlist->m_nFI_fCalcPx1,(float) a2fOffset[1]);
            };
        } else {
            a2x2fEllipseA[0] = fGetField(m_poReflnlist->m_nFI_fEllipsoidA00);
            a2x2fEllipseA[1] = a2x2fEllipseA[2] = fGetField(m_poReflnlist->m_nFI_fEllipsoidA01);
            a2x2fEllipseA[3] = fGetField(m_poReflnlist->m_nFI_fEllipsoidA11);
            a2fEllipseb[0] = fGetField(m_poReflnlist->m_nFI_fEllipsoidb0);
            a2fEllipseb[1] = fGetField(m_poReflnlist->m_nFI_fEllipsoidb1);
            pfEllipsec[0] = fGetField(m_poReflnlist->m_nFI_fEllipsoidc);
            if ((a2fOffset) && (m_poReflnlist->m_nFI_fObsPx0>=0) && (m_poReflnlist->m_nFI_fObsPx1>=0)) {
                a2fOffset[0] = fGetField(m_poReflnlist->m_nFI_fObsPx0);
                a2fOffset[1] = fGetField(m_poReflnlist->m_nFI_fObsPx1);
            } else if ((a2fOffset) && (m_poReflnlist->m_nFI_fCalcPx0>=0) && (m_poReflnlist->m_nFI_fCalcPx1>=0)) {
                a2fOffset[0] = fGetField(m_poReflnlist->m_nFI_fCalcPx0);
                a2fOffset[1] = fGetField(m_poReflnlist->m_nFI_fCalcPx1);
            };
        };
        return 0;
    } else
        return 1;            

};

