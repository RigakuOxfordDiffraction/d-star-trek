//
// Copyright (c) 1998 Molecular Structure Corporation
// All rights reserved.
//
// dtreflnoverlap.cc        Initial author: J.W. Pflugrath       ??-Jul-1998
//    Determines spatial overlaps within a Creflnlist.
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
//  
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <iostream.h>
#include "Creflnlist.h"   // includes Crefln.h
#include "Cimage.h"
#include "Ccrystal.h"

//+Code begin

int main()
{
  Creflnlist *poReflnlist;
  Cstring     sFilename;

  int         nStat;
  int         i, j;
  int   a3nFields[3];
  float a3fOverlapCloseness[3];
  int nNumOverlaps;
  int nResultField;
  float a3fWidth[3];

  //////////////////////////////////////////////////////////////  
  // Begin of testing of overlap method
  //////////////////////////////////////////////////////////////

  cout << "Enter the file name of the reflection list (e.g., dtpredict.ref): ";
  cin >> sFilename;

  poReflnlist = new Creflnlist(sFilename);

  // Try out the ovlp function

  cout << "Enter 0 for observed, 1 for calculated fields: ";
  cin >> i;
  if (0 == i)
    {
      a3nFields[0] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfObsPx0);
      a3nFields[1] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfObsPx1);
      a3nFields[2] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfObsRotMid);
    }
  else
    {
      a3nFields[0] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcPx0);
      a3nFields[1] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcPx1);
      a3nFields[2] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcRotMid);
    }

//  cout << "a3nFields" << ' ' << a3nFields[0]
//       << ' ' << a3nFields[1]  << ' ' << a3nFields[2] << endl;

  cout << "Enter 3 closeness numbers for Px0 (px), Px1 (px), Rot (deg): ";
  cin >> a3fOverlapCloseness[0] >> a3fOverlapCloseness[1] 
      >> a3fOverlapCloseness[2];

  cout << "Enter refln size in pixels in 2 pixel directions: ";
  cin >> a3fWidth[0] >> a3fWidth[1];

  if (0 == i)
    a3fWidth[2] = -poReflnlist->nGetFieldIndex(poReflnlist->ms_sfObsRotWidth);
  else
    a3fWidth[2] = -poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcRotWidth);

  if (0.0 < a3fWidth[2])
    cout << "WARNING no rotation width field!" << endl;

  cout << "...working..." << endl;

  nResultField = poReflnlist->nGetFieldIndex(poReflnlist->ms_snNonunfFlag,
					     eReflnField_int_type);
//  cout << "fieldname: " << poReflnlist->ms_snNonunfFlag << endl;
//  cout << "nResultField: " << nResultField << endl;
//  if (0 > nResultField) nResultField = 4;
  nNumOverlaps = poReflnlist->nOverlapCheck(3, (const int*)a3nFields,
					    (const double*)a3fOverlapCloseness, 
					    (const double*)a3fWidth, 
					    nResultField, 0x3);

  cout << "Number of overlaps: " << nNumOverlaps << endl;

  cout << "Write out overlaps, non-overlaps, all reflns [0,1,2]: ";
  cin >> nNumOverlaps;

  if (0 == nNumOverlaps)
    {
      // Write out only overlapped reflns

      poReflnlist->nSelect("-nNonunf_flag==0");
      poReflnlist->nDelete((const char*)NULL, (const bool)FALSE);
    }
  else if (1 == nNumOverlaps)
    {
      // Write out only non-overlapped reflns

      poReflnlist->nSelect("-nNonunf_flag!=0");
      poReflnlist->nDelete((const char*)NULL, (const bool)FALSE);
    }

  nStat = poReflnlist->nWrite("overlap.ref");
  if (0 == nStat)
    {
      cout << "File overlap.ref written." << endl;
    }
  else
    {
      cout << "ERROR writing file overlap.ref." << endl;
    }

  //////////////////////////////////////////////////////////////
  // End of testing of ovlp method
  //////////////////////////////////////////////////////////////

  delete poReflnlist;

  return (nStat);
}
