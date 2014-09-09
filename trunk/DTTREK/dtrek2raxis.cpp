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
// dtrek2raxis.cc     Initial author: J.W. Pflugrath           05-Jul-1995
//    This reads a D*TREK style image and writes an R-AXIS style image.
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
//    dtrek2raxis reads from stdin the filename of a D*TREK
//    to read and the name of a R-AXIS image to write.  It reads in the input
//    input and writes out the output image.
//
//    Example:  % dtrek2raxis
//                  dtrek.img raxis.img
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
//#include <iostream.h>
#include "Cstring.h"
#include "Cimage.h"
#include "raxis.h"

using std::cin;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;

//+Function prototypes

int main()

{
  Cimage  *poImage;
  Cstring  sIn, sOut;
  int      nStat;

  cerr << "\n\ndtrek2raxis: Copyright (c)2006 Rigaku\n"
       << endl;

  cout << "dtrek2raxis: Enter input d*TREK image name: ";
  if (!(cin >> sIn)) {
    cerr << "dtrek2raxis: ABORTING: invalid image name." << endl;
    return (1);
  }

  poImage = new Cimage(sIn);
//  if (!poImage->bIsAvailable) return (-1);
  if (!poImage->bIsAvailable()) return (-1);

  float fBad;
  int   nDim0, nDim1;
  int   i, j;
  (void) poImage->nGetDimensions(&nDim0, &nDim1);

  cout << "dtrek2raxis: Enter pixel value for bad pixels: ";
  cin >> fBad;

  
  if (fBad != 65535.0)
    {
      for (j = 0; j < nDim1; j++)
	{
	  for (i = 0; i < nDim0; i++)
	    {
	      if ( (poImage->*poImage->prfGetPixel)(i, j) == 65535.0)
		{ 
		  (void) poImage->nSetPixel(i, j, fBad);
		}
	    }
	}
    }
   
  cout << "dtrek2raxis: Enter output R-AXIS image name: ";
  if (!(cin >> sOut)) {
    cerr << "dtrek2raxis: ABORTING: invalid image name." << endl;
    return (1);
  }

  // Write out image

  nStat = nWriteRAXIS(*poImage, sOut);

  if (nStat != 0)
    cerr << "dtrek2raxis: Error writing output image: " << nStat << "!\n"; 
  else
    cerr << "dtrek2raxis: Done" << endl;
  delete poImage;
  return (nStat);
}
