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
// tCFind.cc     Initial author: J.W. Pflugrath           17-Apr-1995
//    This file tests the class Cfind.
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
//    Example:  % tCfind
//                 Enter image to get scan info from: 1.img
//                  
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include <iostream.h>
//#include "String.h"
#include "Cstring.h"
#include "Dtrek.h"
#include "Cfind.h"       // Includes many other include files
#include "Cimage.h"
#include "Cscan.h"
#include "Creflnlist.h"

int main()
{
  int          nStat;
//  String       sName;
  Cstring       sName;
  Cimage      *poImage;
  Cscan       *poScan;
  Creflnlist  *poReflnlist;
  Cfind       *poFind;

  Cimage_header *poImageHeader;

  cout << "Enter the name of image with info to construct SCAN from: ";
  cin >> sName;

  poImage     = new Cimage(sName);           // Create (read) image

//  poScan      = new Cscan(poImage->oHeader); // Create scan

  poReflnlist = new Creflnlist();            // Create empty reflection list

  poReflnlist->nListFields();
  poReflnlist->nList();

//  poFind      = new Cfind(*poImage, poReflnlist); // Create find object.

  poFind = new Cfind(poImage->oHeader, poReflnlist);

  poReflnlist->nListFields();
  poReflnlist->nList();

  delete poImage;
  cout << "Enter name of image to search for peaks:";
  cin >> sName;

  poImage = new Cimage(sName);

  int   nArea[4];
  float fExcludeMin, fExcludeMax;
  float fAvg, fSD;

/*
  cout << "Enter origin, extents of area to calculate background in:";
  cin >>  nArea[0] >> nArea[1] >> nArea[2] >> nArea[3];
  cout << "Enter min, max exclusion limits:";
  cin >> fExcludeMin >> fExcludeMax;
  
  nStat = poFind->nBackgroundAvgSD(*poImage, nArea, fExcludeMin, fExcludeMax,
                           &fAvg, &fSD);                          
*/

  nStat = poFind->nFind2D(*poImage, (const float)0.0);

//  cout << "Average is: " << fAvg << "\nSD is: " << fSD << endl;

  poReflnlist->nListFields();
  poReflnlist->nList();

  cout << "Enter reflection list filename:";
  cin >> sName;
  poReflnlist->nWrite(sName);


  return (0);
}



