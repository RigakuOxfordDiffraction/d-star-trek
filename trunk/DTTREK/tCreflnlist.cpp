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
// tCreflnlist.cc        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file tests the classes Creflnlist and Crefln.
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
  Creflnlist  oReflnlist;
  Creflnlist *poReflnlist;
  Cstring     sFilename;

  int         nStat, nH, nK, nL;
  int         i, j;
  int   a3nFields[3];
  float a3fOverlapCloseness[3];
  float a3fWidth[3];

  //////////////////////////////////////////////////////////////  
  // begin of testing of ovlp method
  //////////////////////////////////////////////////////////////

  cout << "Enter the file name of the reflection list (e.g., dtpredict.ref): ";
  cin >> sFilename;

  poReflnlist = new Creflnlist(sFilename);

  // try out the ovlp function
  cerr << "testing nOverlapCheck...\n";
  a3nFields[0] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcXmm);
  a3nFields[1] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcYmm);
  a3nFields[2] = poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcRotMid);
  cout << "Enter 3 closeness numbers for Xmm, Ymm, Rot: ";
  cin >> a3fOverlapCloseness[0] >> a3fOverlapCloseness[1] 
      >> a3fOverlapCloseness[2];
  a3fWidth[0] = 0;
  a3fWidth[1] = 0;
  a3fWidth[2] = -poReflnlist->nGetFieldIndex(poReflnlist->ms_sfCalcRotWidth);
  int nNumOverlaps;
  int nResultField;
  nResultField = poReflnlist->nGetFieldIndex(poReflnlist->ms_snNonunfFlag);
  nNumOverlaps = poReflnlist->nOverlapCheck(3, a3nFields,
					    a3fOverlapCloseness, 
					    a3fWidth, nResultField, 0x3);

  cout << "Number of overlaps: " << nNumOverlaps << endl;
  // Write out only overlapped reflns
  poReflnlist->nSelect("-nNonunf_flag==0");
  poReflnlist->nDelete((const char*)NULL, (const bool)FALSE);
  poReflnlist->nWrite("overlap.ref");
  //////////////////////////////////////////////////////////////
  // end of testing of ovlp method
  //////////////////////////////////////////////////////////////

  Crefln* poRefln1 = new Crefln (&oReflnlist);

  Crefln* poRefln2 = new Crefln(*poRefln1);

  Crefln* poRefln3 = new Crefln(*poRefln1);

  *poRefln3 = *poRefln2;
/*
  poRefln1->nList();
  poRefln2->nList();
  poRefln3->nList();
*/
  Cimage_header   *poHeader;
  Ccrystal    *poCrystal;
  Cstring      sName;
  
  
  cout << "Enter name of image to get crystal info from: ";
  cin >> sName;

  poHeader  = new Cimage_header(sName);

  poCrystal = new Ccrystal(*poHeader);

  poReflnlist = new Creflnlist(poCrystal, 10., 5.0);
  nStat = poReflnlist->nCountUnique();

  cout << "There were " << nStat << " unique reflections"
       << " in the list of " << poReflnlist->nGetNumReflns() << "." << endl;

  cout << "Enter name of reflection file to write: ";
  cin >> sFilename;
  poReflnlist->nWrite(sFilename, poReflnlist->pnGetSortIndex());
  

  cout << "Enter name of reflection file: ";
  cin >> sFilename;
  poReflnlist = new Creflnlist(sFilename);
  if (!poReflnlist->bIsAvailable())
    {
      cerr << "Error reading reflection list file: " << sFilename << endl;
      delete poReflnlist;
    }

  int nNumRef = poReflnlist->nGetNumReflns();
  float fSum = 0.0;
  float fSumI = 0;
  Crefln *poRefln;

  j = 0;
  for (i = 0; i < nNumRef; i++)
    {
      poRefln = poReflnlist->poGetRefln(i);
      if (0.0 < poRefln->fGetSigmaI())
	{
	  fSum = fSum   + poRefln->fGetIntensity() / poRefln->fGetSigmaI();
	  fSumI = fSumI + poRefln->fGetIntensity();
	  j++;
	}
    }
  if (0 < j)
    {
      cout << "Average I/SigmaI for " << j << " reflns is: " << fSum / (float)j
           << "\nAverage I        for " << j << " reflns is: " << fSumI / (float)j
           << endl;
    }

  poReflnlist->nReduce(*poCrystal);
  nStat = poReflnlist->nCountUnique();

  cout << "There were " << nStat << " unique reflections"
       << " in the list of " << poReflnlist->nGetNumReflns() << "." << endl;

  cout << "Enter name of reflection file to write: ";
  cin >> sFilename;
  poReflnlist->nWrite(sFilename, poReflnlist->pnGetSortIndex());

  (void) poReflnlist->nList((const int)10);
  delete poReflnlist;

  poReflnlist = new Creflnlist((const int)2);

  for (i = 0; i < 10; i++) {
    cout << "Enter an hkl value: ";
    cin >> nH >> nK >> nL;
    nStat = poReflnlist->nInsert(nH, nK, nL);
  }

  cout << "0\n";
  (void) poReflnlist->nList((const int)0);
  cout << "\n\n5\n";
  (void) poReflnlist->nList((const int)5);
  cout << "\n\n20\n";
  (void) poReflnlist->nList((const int)20);


  cout << "Enter name of reflection file to append to the list: ";
  cin >> sFilename;
  poReflnlist->nRead(sFilename);


  cout << "Enter name of reflection file to write: ";
  cin >> sFilename;
  poReflnlist->nWrite(sFilename);
  (void) poReflnlist->nList((const int)20);

  delete poReflnlist;

  return 0;
}
