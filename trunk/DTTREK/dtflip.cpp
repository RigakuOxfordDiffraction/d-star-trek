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
// dtflip.cc     Initial author: J.W. Pflugrath           13-Jun-1996
//    This flips an image 
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
//    dtflip reads from stdin the filename of a D*TREK image file
//    to read and the name of an image file to write.  It reads in the input
//    image, flips it and writes out the result to the output image.
//
//    Example:  % dtflip
//                  unflipped.img
//                  flipped.img
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Cstring.h"
#include "Cimage.h"

int main()

{
  Cimage  *poImage;
  Cstring  sIn, sOut;
  int      nStat;
  int      nDim0, nDim1;

  cerr << "\n\ndtflip: Copyright (c) 2006 Rigaku\n"
       << endl;

  cout << "dtflip: Enter input unflipped image name: ";
  if (!(cin >> sIn)) 
    {
      cerr << "dtflip: ABORTING: invalid image name." << endl;
      return (1);
    }

  poImage = new Cimage(sIn);
//  if (!poImage->bIsAvailable) return (-1);
  if (!poImage->bIsAvailable()) return (-1);

  cout << "\ndtflip: Enter output flipped image name: ";
  if (!(cin >> sOut))
    {
      cerr << "dtflip: ABORTING: invalid image name." << endl;
      return (1);
    } 

  nStat = poImage->nGetDimensions(&nDim0, &nDim1);

  Cimage oOutImage (*poImage, NULL);
//  Cimage oOutImage(nDim1, nDim0, eImage_uI2);

  int i, j, ii, jj, iii;    // Loop counters
  unsigned short int uiTemp1, uiTemp2, uiTemp3, uiTemp4;

  poImage->nSetNextPixel(0, 0);

  for (j = 0; j < nDim1; j++)
    {
//      jj = nDim0 - j - 1;
      jj = nDim1 - j - 1;
      for (i = 0; i < nDim0; i++)
        {
          uiTemp1 = poImage->uiGetNextPixel();   
//	  ii = nDim1 - i - 1;
	  ii = nDim0 - i - 1;
          oOutImage.nSetPixel(     i, jj, uiTemp1);
        }
    }

  // Write out image, but first adjust header

  oOutImage.oHeader = poImage->oHeader;
//  oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sSize1, nDim1);
//  oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sSize2, nDim0);

  nStat = oOutImage.nWrite(sOut);

  if (0 != nStat)
    cerr << "\ndtflip: Error writing output image: " << nStat << "!\n"; 
  else
    cerr << "\ndtflip: Done" << endl;
  delete poImage;
  return (nStat);
}
