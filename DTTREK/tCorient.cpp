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
// tCpredict.cc     Initial author: J.W. Pflugrath        2-May-1995
//    This file tests the class Cpredict.
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
//    Example:  % tCorient
//                 Enter image to get crystal and goniometer info from: 1.img
//                  
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include <iostream.h>
#include "Cstring.h"
#include "Dtrek.h"
#include "Cimage_header.h"
#include "Corient.h"
#include "dtrekvec.h"

int main()
{
  int            nStat;
  Cstring        sName;
  Cimage_header *poHeader;
  Corient       *poOrient;
  Cstring        sOrientWish;

  cout << "Enter name of file to get crystal and goniometer info from: ";
  cin >> sName;
  poHeader = new Cimage_header(sName);
  
  if (!poHeader->bIsAvailable())
    {
      cerr << "Problem reading header of file " << sName;
      delete poHeader;
      return (-1);
    }

  poOrient  = new Corient(*poHeader);

  nStat = poOrient->nCalcAngles();
  nStat = poOrient->nList();

  cout << "\nEnter orient request: ";
  getline(cin, sOrientWish);

  nStat = poOrient->nSetOrient(sOrientWish);
  nStat = poOrient->nCalcOrient();
  nStat = poOrient->nPickSoln();
  if (0 <= nStat)
    {
      nStat = poOrient->nUpdateHeader(poHeader);
      nStat = poOrient->nList();
      cout << "\nEnter name of image output file to write the results to: ";
      cin >> sName;
      if (cin)
	{
	  nStat = poHeader->nWrite(sName);
	}
    }
  delete poOrient;
  delete poHeader;
  return (nStat);
}



