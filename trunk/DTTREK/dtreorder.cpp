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
// dtreorder.cc     Initial author: J.W. Pflugrath           13-Sep-1995
//    This re-orders a multiple readout CCD image so that it makes sense.
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
//    dtreorder reads from stdin the filename of a D*TREK image file
//    to read and the name of an image file to write.  It reads in the input
//    image, reorders it and writes out the result to the output image.
//
//    Example:  % dtreorder
//                  scrambled.img 
//                  unscrambled.img
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <iostream.h>
#include "Cstring.h"
#include "Cimage.h"

int main()

{
  Cimage  *poImage;
  Cstring  sIn, sOut;
  int      nStat;
  int      nDim0, nDim1;

  cerr << "\n\ndtreorder: Copyright (c) 2006 Rigaku\n"
       << endl;

  cout << "dtreorder: Enter input scrambled image name: ";
  if (!(cin >> sIn)) {
    cerr << "dtreorder: ABORTING: invalid image name." << endl;
    return (1);
  }

  poImage = new Cimage(sIn);
//  if (!poImage->bIsAvailable) return (-1);
  if (!poImage->bIsAvailable()) return (-1);

  cout << "\ndtreorder: Enter output unscrambled image name: ";
  if (!(cin >> sOut)) {
    cerr << "dtreorder: ABORTING: invalid image name." << endl;
    return (1);
  }

  // Do some descrambling here


  Cimage oOutImage (*poImage, NULL);

  int i, j, ii, jj, iii;    // Loop counters
  unsigned short int uiTemp1, uiTemp2, uiTemp3, uiTemp4;

  nStat = poImage->nGetDimensions(&nDim0, &nDim1);
  poImage->nSetNextPixel(0, 0);

//  printf("\a \n nDim0 = %d, nDim1 = %d \n\n", nDim0,nDim1);

  // nDim0 is usually 2086  fast pixel direction use variables i, ii
  // nDim1 is usually 2046  slow pixel direction use variables j, jj
  // So here are some possibilities:
  //
  //   i :                  0 --- 1042   left half normal
  //  ii :               1043 --- 2085   right half normal  (i + nDim0 / 2)
  // nDim0 - 1 - i:      2085 --- 1043   right half reverse
  // nDim0 / 2 - 1 - i:  1042 --- 0      left half reverse

  //   j :                  0 --- 1022   top half normal
  //  jj :               2045 --- 1023   bottom half reverse  (nDim1 - 1 - j)
  // j + nDim1 / 2:      1023 --- 2045   bottom half normal
  // nDim1 / 2 - 1 - j:  1022 --- 0      top half reverse

  
  for (j = 0; j < 1023; j++)
    {
      jj = 1023+j;

      ii = 4172/4;
      for (i = 0; i < ii; i++)                    
        {
          
	  //Quadrant A
          uiTemp1 = poImage->uiGetNextPixel();   
          uiTemp1 = (uiTemp1 & 0xFFFF); 

	  //Quadrant B 
          uiTemp2 = poImage->uiGetNextPixel();
          uiTemp2 = (uiTemp2 & 0xFFFF); 

	  //Quadrant C
          uiTemp3 = poImage->uiGetNextPixel();
          uiTemp3 = (uiTemp3 & 0xFFFF); 

	  //Quadrant D
          uiTemp4 = poImage->uiGetNextPixel();
          uiTemp4 = (uiTemp4 & 0xFFFF); 

          oOutImage.nSetPixel(     i,  j, uiTemp3);
          oOutImage.nSetPixel(  ii+i,  j, uiTemp4);  
          oOutImage.nSetPixel(1043-i-1, 3069-jj-1, uiTemp1);
          oOutImage.nSetPixel(2086-i-1, 3069-jj-1, uiTemp2);



//          oOutImage.nSetPixel(0*ii+i,  j, (unsigned short int) 100);
//          oOutImage.nSetPixel(1*ii+i,  j, (unsigned short int) 200);
//          oOutImage.nSetPixel(0*ii+i, jj, (unsigned short int) 300);
//          oOutImage.nSetPixel(1*ii+i, jj, (unsigned short int) 400);

        }
    }

  // Write out image

  nStat = oOutImage.nWrite(sOut);

  if (nStat != 0)
    cerr << "\ndtreorder: Error writing output image: " << nStat << "!\n"; 
  else
    cerr << "\ndtreorder: Done" << endl;
  delete poImage;
  return (nStat);
}


