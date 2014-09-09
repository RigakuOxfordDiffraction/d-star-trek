//
// Copyright (c) 2006 Rigaku
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtmakemarks.cc     Initial author: J.W. Pflugrath           10-Mar-1995
//    This reads an input image and labels pixels.
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
//    dtmakemarks reads from stdin the filename of an image and the extents
//    of a submodule in the image.  It also reads in a minimum pixel value.
//    Then it creates a new image with 1000*module number for pixel values
//    except for those pixels which fall below the minimum pixel value which
//    are set to 0.
//
//    Example:  % dtmakemarks
//                  input.img  1024 1024 10000 mark.img
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include "Cstring.h"
#include "Cimage.h"

int main()

{
  Cimage  *poBig_image, *poMark_image;
  int      nDim0, nDim1;
  int      nBin0, nBin1;
  int      nStat;
  int      i, j, k, l;
  
  Cstring             sOut,sComment;
  unsigned short int  iFrom;
  unsigned short int* piTo;
  unsigned short int  iMark;
  unsigned short int  iMin;
  int                 nBad;
  int                 nGood;


  cerr << "\n\ndtmakemarks: Copyright (c) 2006 Rigaku\n"
       << endl;

  cout << "dtmakemarks: Enter input image name,\n" <<
          "             the extents of a submodule, and min pix value: ";
  if (!(cin >> sOut >> nBin0 >> nBin1 >> iMin) || (nBin0 < 1 || nBin1 < 1) ) {
    cerr << "dtmakemarks: ABORTING: invalid extents." << endl;
    return 1;
  }

  poBig_image = new Cimage(sOut);
  
  nDim0 = poBig_image->m_nDim[0];
  nDim1 = poBig_image->m_nDim[1];

  cout << "    Input image size: " << nDim0 << ", " << nDim1 << ".\n";
  nDim0 = nDim0 / nBin0;
  nDim1 = nDim1 / nBin1;
  if ( ( (nBin0 * nDim0) != poBig_image->m_nDim[0]) ||
       ( (nBin1 * nDim1) != poBig_image->m_nDim[1]) ) {
    cerr << "\ndtmakemarks: Error in submodule extents!\n"; 
    return (3);
  }

  poMark_image = new Cimage(poBig_image->m_nDim[0], poBig_image->m_nDim[1],
			    eImage_I2);

  piTo = poMark_image->m_The_Data.pui;
  nBad  = 0;
  nGood = 0;
  iMark = 0;
  unsigned short int iZero = 0;

  for (i = 0; i < poBig_image->m_nDim[1]; i = i + nBin1)
    {
      for (j = 0; j < poBig_image->m_nDim[0]; j = j + nBin0)
	{
	  iMark = iMark + 1000;
	  cout << "    New submodule starts at " << j << ", " << i << ".\n"; 
	  for (k = 0; k < nBin1; k++)
	    {
	      for (l = 0; l < nBin0; l++)
		{
		  nStat = poBig_image->nGetPixel(j+l, i+k, &iFrom);
		  if (nStat != 0)
		    {
		      cerr << "\ndtmakemarks: Error getting pixel from image!\n";
		      return (4);
		    }
		  if (iFrom >= iMin)
		    {
		      nGood++;
		      nStat = poMark_image->nSetPixel(j+l, i+k, iMark);
		    }
		  else
		    {
		      nBad++;
		      nStat = poMark_image->nSetPixel(j+l, i+k, iZero);
		    }
		}
	    }
	}
    }
  // Write out image

  sComment = "Made by dtmakemarks from image " + sOut;
  (void) poMark_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sComment);

  cout << "dtmakemarks: Enter output image name: ";

  if ( !(cin >> sOut) ) {
    cerr << "\ndtmakemarks: Error reading output filename!\n";
  }
  cout << "dtmakemarks:  Number of good pixels " << nGood << endl;
  cout << "              Number of  bad pixels " << nBad << endl;

  nStat = poMark_image->nWrite(sOut);

  if (nStat != 0)
    cerr << "dtmakemarks: Error writing output image: " << nStat << "!\n"; 
  else
    cerr << "dtmakemarks: Done" << endl;
  delete poBig_image;
  delete poMark_image;

  return (nStat);
}
