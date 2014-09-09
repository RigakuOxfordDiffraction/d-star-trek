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
// dtdarksubtract.cc     Initial author: J.W. Pflugrath           10-Mar-1995
//    Subtract a dark image multiplied by a scale factor from another image.
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
//    dtdarksubtract reads from stdin the name of dark image and a scale 
//    factor to multiply it by.  Next it reads in the name of another image.
//    A new image is created such that new = old - (scale * dark). Any negative
//    pixel that results is given a value of 0.
//
//    Example:  % dtdarksubtract
//                  dark.img  1.10
//                  input.img
//                  output.img
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include "Cstring.h"
#include "Cimage.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

int main()

{
  Cimage  *poDark_image;
  Cimage  *poInput_image;
  int     nStat, nSub;
  int     i, j;
  int     nTruncate;
  float   fScale    = 1.0;
  float   fScaleOut = 1.0;
  Cstring sName, sComment;

  cout << "\n\ndtdarksubtract: Copyright (c) 2006 Rigaku\n"
       << endl;



  nStat = 0;
  cout << "dtdarksubtract: Enter dark image name: ";
   if (cin >> sName) {
    poDark_image = new Cimage (sName);
    if (poDark_image->bIsAvailable())
      {
	sComment = "Made by dtdarksubtract from dark image: " + sName;
	cout << "dtdarksubtract: Enter a scale factor to multiply dark image by: ";
	cin >> fScale;
	cout << "dtdarksubtract: Enter input image name: ";
	cin >> sName;
	poInput_image = new Cimage(sName);

	cout << "dtdarksubtract: Enter a scale factor to multiply output image by: ";
	cin >> fScaleOut;

	if ( (poInput_image->bIsAvailable())
	    && (poDark_image->nGetDimension(0) == poInput_image->nGetDimension(0))
	    && (poDark_image->nGetDimension(1) == poInput_image->nGetDimension(1)) )
	  {
	    j = poDark_image->nGetDimension(0) * poDark_image->nGetDimension(1);
	    
	    cout << "INFO: Scale factor for dark image is   " << fScale << endl;
	    cout << "INFO: Scale factor for output image is " << fScaleOut << endl;

	    nTruncate = 0;
	    float fValIn, fValDark;
	    float fValOut;
	    int nDim0, nDim1;
	    poInput_image->nGetDimensions(&nDim0, &nDim1);
	    for (j = 0; j < nDim1; j++) 
	      {
		for (i = 0; i < nDim0; i++) 
		  {
		    fValIn   = (poInput_image->*poInput_image->prfGetPixel)(i, j);
		    fValDark = (poDark_image->*poDark_image->prfGetPixel)(i, j);
		    fValOut  = fValIn -  (fScale * fValDark);
		    if (fValOut < 0.0) 
		      {
			fValOut = 0.0;
			nTruncate++;
		      }
		    fValOut *= fScaleOut;
		    (poInput_image->*poInput_image->prnSetPixel)(i, j, fValOut);
		  }
	      }
	  }
      }
  }

  cout << "Number of pixels truncated to 0: " << nTruncate << endl;
  // Write out image
  
  sComment = sComment + " and input image: " + sName;
  (void) poInput_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sComment);
  cout << "dtdarksubtract: Enter output image name: ";
  cin >> sName;
  nStat = poInput_image->nWrite(sName);

  if (nStat != 0)
    cout << "dtdarksubtract: Error writing output image: " << nStat << "!\n"; 
  else
    cout << "dtdarksubtract: Done" << endl;
  delete poDark_image;
  delete poInput_image;
  return (nStat);
}
