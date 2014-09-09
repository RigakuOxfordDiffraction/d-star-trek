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
// dtsmooth.cc     Initial author: J.W. Pflugrath           25-Jun-1996
//    This smooths an image 
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
//    dtsmooth reads from stdin the filename of a d*TREK image file
//    to read, min pixel value, max pixel value, smooth number (nS),  and 
//    the name of an image file to write.  It reads in the input
//    image, smooths it and writes out the result to the output image.
//
//    Example:  % dtsmooth
//                  unsmoothed.img
//                  0 65535 5
//                  smoothed.img
//   An output pixel at (i,j) is the average of the +-nS pixels centered on
//   (i,j) with pixel (i,j) weighted by 8 in the sum.

//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Cstring.h"
#include "Cimage.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

int main()

{
  Cimage  *poImage;
  Cstring  sIn, sOut;
  int      nStat;
  int      nDim0, nDim1;
  int      nPixel;
  float    fMin, fMax;

  cout << "\n\ndtsmooth: Copyright (c) 2006 Rigaku\n";

  cout << "dtsmooth: Enter input unsmoothed image name: ";
  if (!(cin >> sIn)) 
    {
      cout << "dtsmooth: ABORTING: invalid image name.\n";
      return (1);
    }

  poImage = new Cimage(sIn);
  if (!poImage->bIsAvailable()) return (-1);

  cout << "\ndtsmooth: Enter min, max allowed pixel values: ";
  if (!(cin >> fMin))
    {
      cout << "dtsmooth: ABORTING: invalid min allowed pixel value.\n";
      return (1);
    }
  if (!(cin >> fMax))
    {
      cout << "dtsmooth: ABORTING: invalid min allowed pixel value.\n";
      return (1);
    }
  cout << "\ndtsmooth: Enter smooth number: ";
  if (!(cin >> nPixel))
    {
      cout << "dtsmooth: ABORTING: invalid smooth number:\n";
      return (1);
    }

  cout << "\ndtsmooth: Enter output smoothed image name: ";
  if (!(cin >> sOut))
    {
      cout << "dtsmooth: ABORTING: invalid image name.\n";
      return (1);
    } 

  cout << "...working...\n" << flush;

  nStat = poImage->nGetDimensions(&nDim0, &nDim1);

  Cimage oOutImage (*poImage, NULL);

  int i, j, ii, jj;    // Loop counters
  int jlo, jhi, ilo, ihi;

  float fSum, fTemp;
  int   nSum;
  for (j = 0; j < nDim1; j++)
    {
      jlo = max(0, j - nPixel);
      jhi = min(nDim1-1, j + nPixel);
      for (i = 0; i < nDim0; i++)
        {
	  // This is not the most efficient way, but what the hey...

	  ilo = max(0, i - nPixel);
	  ihi = min(nDim0-1, i + nPixel);
	  fSum = 0.0;
	  nSum = 0;

	  for (jj = jlo; jj <= jhi; jj++)
	    {
	      for (ii = ilo; ii <= ihi; ii++)
		{
		  fTemp = (poImage->*poImage->prfGetPixel)(ii, jj);
		  if ( (fMax >= fTemp) && (fMin <= fTemp) )
		    {
		      fSum += fTemp;
		      nSum++;
		      if (  (ii == i)  && (jj == j) )
			{
			  // Overweight central pixel by 8 more

			  fSum = fSum + (8.0f * fTemp);
			  nSum = nSum + 8;
			}
		    }
		}
	    }
	  if (0 < nSum)
	    {
	      fSum = fSum / (float) nSum;
	      (oOutImage.*oOutImage.prnSetPixel)(i, j, fSum);
	    }
	  else
	    {
	      (oOutImage.*oOutImage.prnSetPixel)(i, j, 
				      (poImage->*poImage->prfGetPixel)(i, j));
	    }
        }
    }

  // Write out image, but first adjust header

  oOutImage.m_oHeader = poImage->m_oHeader;
  nStat = oOutImage.nWrite(sOut);

  if (0 != nStat)
    cout << "\ndtsmooth: Error writing output image: " << nStat << "!\n"; 
  else
    cout << "\ndtsmooth: Done" << endl;
  delete poImage;
  return (nStat);
}
