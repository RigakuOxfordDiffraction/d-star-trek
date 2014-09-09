//
// Copyright (c) 2006 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtdezinger.cc     Initial author: J.W. Pflugrath           10-Mar-1995
//    This tries to remove zingers from an input images.
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
//    dtdezinger reads from stdin the filename of an image that should have
//    zingers removed.  It tries to remove zingers by replacing pixels that
//    have a value above 60,000 with an average of the none 60K pixels around
//    it. Then it reads from stdin the name of the image to write out and 
//    does so.
//   
//
//    Example:  % dtdezinger
//                  zinger.img  nozinger.img
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream

#include "Cstring.h"
#include "Cimage.h"
#include "minmax.h"

using std::cout;
using std::endl;
using std::flush;
using std::cin;
using std::cerr;

int main()

{
  Cimage  *poInput_image, *poOutput_image;
  int      nDim0, nDim1;
  int      nStat;
  int      i, j, k, l, lmin, lmax, kmin, kmax;
  
  Cstring   sOut,sComment;
  unsigned short int  uiFrom;
  unsigned short int  uiMax;
  unsigned short int  uiLast;
  float               fSum;
  int                 nSum;

  cerr << "\n\ndtdezinger: Copyright (c) 2006 Rigaku\n"
       << endl;

  cout << "dtdezinger: Enter input image name: ";
  if (!(cin >> sOut) )
    {
      cerr << "dtdezinger: Error reading input image filename.\n";
      return 1;
    }

  cout << "dtdezinger: Enter max legitimate value in input image: ";
  if (!(cin >> uiMax) )
    {
      cerr << "dtdezinger: Error reading max value: \n";
      return 1;
    }

  poInput_image = new Cimage(sOut);
  
  (void) poInput_image->nGetDimensions(&nDim0, &nDim1);

  poOutput_image = new Cimage(nDim0, nDim1, poInput_image->m_eData_type);

  sComment = "Made by dtdezinger from image: " + sOut;

  uiLast = 0;
  int nZingers = 0;

  for (i = 0; i < nDim1; i++) {
    for (j = 0; j < nDim0; j++) {
      nStat = poInput_image->nGetPixel(j, i, &uiFrom);
      if (uiFrom <= uiMax) {
	nStat = poOutput_image->nSetPixel(j, i, uiFrom);
      }
      else {
	nZingers++;
	// Compute average around bogus pixel and replace it with the average
	kmin = max(0, j - 1);
	kmax = min(nDim1, j + 1);
        lmin = max (0, l - 1);
        lmax = min (nDim1, l + 1);
	fSum = 0.0;
        nSum = 0;
	for (k = kmin; k < kmax; k++) {
	  for (l = lmin; l < lmax; l++) {
	    nStat = poInput_image->nGetPixel(l, k, &uiFrom);
	    if ( (nStat == 0) && (uiFrom <= uiMax) ) {
	      nSum++;
	      fSum = fSum + (float) uiFrom;
	    }
	  }
        }
	if (nSum > 0) {
	  uiFrom = (unsigned short int) (fSum / (float) nSum);
	  uiLast = uiFrom;
	  nStat = poOutput_image->nSetPixel(j, i, uiFrom);
	}
	else
	  nStat = poOutput_image->nSetPixel(j, i, uiLast);
      } 
    }
  }

  // Write out image

  (void) poOutput_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment,
						 sComment);

  cout << "dtdezinger: There were " << nZingers << " zinger pixels found.\n\n";
  cout << "dtdezinger: Enter output image name: ";

  if ( !(cin >> sOut) ) {
    cerr << "\ndtdezinger: Error reading output filename!\n";
  }
  nStat = poOutput_image->nWrite(sOut);

  if (nStat != 0)
    cerr << "dtdezinger: Error writing output image: " << nStat << "!\n"; 
  else
    cerr << "dtdezinger: Done" << endl;
  delete poInput_image;
  delete poOutput_image;

  return (nStat);
}
