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
//    Example:  % tCpredict
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
#include "Cpredict.h"       // Includes many other include files
#include "Cspacegroup.h"
#include "Cimage.h"
#include "Csource.h"
#include "Cdetector.h"
#include "Ccrystal.h"
#include "Cgoniometer.h"
#include "Corient.h"
#include "dtrekvec.h"

int main()

{
  int          nStat;
//  String       sName;
  Cstring       sName;
  Cspacegroup  *poSp;
  int          nHKL[3];
  int          nHKLr[3];
  Cimage      *poImage;
  Ccrystal    *poCrystal;
  Corient     *poOrient;
  float       fVec[3];
  float       fAngles[6];
//  String      sOrientWish;
  Cstring      sOrientWish;
 
  Creflnlist  oReflnList;
  Crefln      oRefln(&oReflnList);

  cout << "Enter name of image to get crystal info from: ";
  cin >> sName;

  poImage = new Cimage(sName);

  poOrient  = new Corient(poImage->oHeader);

  cout << "Enter orient string: ";
  cin  >> sOrientWish;
  nStat = poOrient->nSetOrient(sOrientWish);
  nStat = poOrient->nCalcOrient();
  (void) poOrient->nListSoln();

  do {
    cout << "Enter a 3 goniometer angles: ";
    cin >> fVec[0] >> fVec[1] >> fVec[2];
    nStat = poOrient->nCalcAngles((float *)fVec);
    nStat = poOrient->nList();
  } while ( (fVec[0] != 0.0) || (fVec[1] != 0.0) || (fVec[2] != 0.0));

  poCrystal = new Ccrystal(poImage->oHeader);
  poCrystal->nList();
  float fTemp;
    do {

    cout << "Enter max resol: ";
    cin >> fTemp;

    nStat = poCrystal->nCalcGetNumUniqueRefs(fTemp, 0.0);
    cout << "Num unique " << nStat << endl;
  } while ( (fTemp > 1.5) && (nStat >= 0) );

  do {
    cout << "Enter the number of spacegroup: ";
    cin >> nStat;
    poSp = new Cspacegroup (nStat);
//    cout << "Spacegroup name: " << poSp->sName << endl;
    cout << "Spacegroup name: " << poSp->sGetName() << endl;
    nStat = poSp->nReadEquivPos();
    
    // 
    do {
      cout << "Enter an HKL: ";
      cin >> nHKL[0] >> nHKL[1] >> nHKL[2];
      oRefln.vSetH(nHKL[0]);
      oRefln.vSetK(nHKL[1]);
      oRefln.vSetL(nHKL[2]);
//      nStat = poSp->nReduceHKL(nHKL, nHKLr);
      nStat = poSp->nReduceHKL(&oRefln, nHKLr);
      cout << "Reduced hkl: " << nHKLr[0] << ", "
           << nHKLr[1] << ", "
	   << nHKLr[2] << "   nStat: " << nStat << endl;
    } while ( (nHKL[0] != 0) || (nHKL[1] != 0) || (nHKL[2] != 0) );
  } while (sName != "q");

  return (0);
}



