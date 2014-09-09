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
// tCrefine.cc     Initial author: J.W. Pflugrath       12-May-1995
//    This file tests the class Crefine.
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
//    Example:  % tCrefine
//                 Enter image to get scan info from: 1.img
//                  
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include <iostream.h>
#include "Dtrek.h"
#include "dtrekvec.h"
#include "Cstring.h"
#include "Crefine.h"       // Includes many other include files
#include "Creflnlist.h"
#include "Cimage.h"

int main()
{
  int            i;
  int            nStat;
  Cstring        sName;
  Cimage_header *poHeader;
  Creflnlist    *poReflnlistIn;
  Crefine       *poRefine;

  // Read the name of an image and read the header only from that image file

  cout << "Enter the name of image with info to construct things from: ";
  cin >> sName;

  poHeader     = new Cimage_header(sName);   // Create (read) image header

  // Read a reflection list from a file

  cout << "Enter the name of the input reflection list file: ";
  cin >> sName;
  poReflnlistIn = new Creflnlist(sName);
  poReflnlistIn->nList();

  poRefine   = new Crefine(*poHeader, poReflnlistIn);

//  poRefine->poReflnlist->nListFields();
//  poRefine->poReflnlist->nList();

  nStat = poRefine->nCopyReflns();

  poRefine->nList(1);

  poRefine->nSetup();

  cout << "tCrefine can refine the following properties "
       << "(T=translation, R=rotation:\n" 
       << "Detector            Source      Crystal\n"
       << "T1 T2 T3 R1 R2 R3   R1 R2 Wvl   R1 R2 R3 a b c alp bet gam "
       << "mosaicity."
       << endl;
  cout << "Refinement flags settings (0 = fixed; 1 = refined): " << endl;
  printf("%2d %2d %2d %2d %2d %2d   %2d %2d %3d   %2d %2d %2d %1d %1d %1d %3d %3d %3d %2d",
	 poRefine->m_nRefineFlags[0],
	 poRefine->m_nRefineFlags[1],
	 poRefine->m_nRefineFlags[2],
	 poRefine->m_nRefineFlags[3],
	 poRefine->m_nRefineFlags[4],
	 poRefine->m_nRefineFlags[5],
	 poRefine->m_nRefineFlags[6],
	 poRefine->m_nRefineFlags[7],
	 poRefine->m_nRefineFlags[8],
	 poRefine->m_nRefineFlags[9],
	 poRefine->m_nRefineFlags[10],
	 poRefine->m_nRefineFlags[11],
	 poRefine->m_nRefineFlags[12],
	 poRefine->m_nRefineFlags[13],
	 poRefine->m_nRefineFlags[14],
	 poRefine->m_nRefineFlags[15],
	 poRefine->m_nRefineFlags[16],
	 poRefine->m_nRefineFlags[17],
	 poRefine->m_nRefineFlags[18]);
/*
  for (i = 0; (i < nMAXPARAMS+1); i++)
    {
      cout << poRefine->m_nRefineFlags[i] << " ";
    }
*/
  cout << "\nEnter new flag or -1 signal the end of input:" << endl;

 for (i = 0; (i < nMAXPARAMS+1) && (nStat != -1); i++) 
   {
     cin >> nStat;
     if ( (0 == nStat) || (1 == nStat) ) poRefine->m_nRefineFlags[i] = nStat;
   }

  cout << "Enter verbose level: ";
  cin >> i;
  poRefine->vSetVerboseLevel(i);
  cout << "Enter number of refinement loops: ";
  cin >> i;
  nStat = poRefine->nRefineLoop(i);
  
  if (0 == nStat)
    {
      // Write results to output

      poRefine->nUpdateHeader(poHeader);

      //  cout << "Enter new header file name to place results in. ";

      Cimage oImage(0,0);
      oImage.oHeader = *poHeader;
      nStat = oImage.nWrite();
    }

  delete poHeader;
  delete poRefine;
  delete poReflnlistIn;
  return (nStat);
}

