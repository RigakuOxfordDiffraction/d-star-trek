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
// dtcompose.cc     Initial author: J.W. Pflugrath           10-Mar-1995
//    This composes a larger image from several smaller images.
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
//    dtcompose reads from stdin the filename and dimensions of an image 
//    to create.
//    Next it reads additional image filenames and their position (in pixels)
//    in the new image.  The pixel values of these additional images are 
//    placed in the composed image.
//    The new image has the same Data_type as the first image placed in it.
//    The Data_types of the subsequent images are converted to this Data_type.
//    When there is an end-of-file on stdin, the composed image is written
//    to the filename given at the beginning.
//    Since each subsequent is placed in the composed image in the order given,
//    it is possible to replace values in the composed image that came from an
//    input image file.
//
//    Information about the progress of dtcompose is written to stderr.
//
//    Example:  % dtcompose 
//                  big.img  3072 3072
//                  sub1.img    0    0
//                  sub2.img 1024    0
//                  sub3.img 2048    0
//                  sub4.img    0 1024
//                  sub5.img 1024 1024
//                  sub6.img 2048 1024
//                  sub7.img    0 2048
//                  sub8.img 1024 2048
//                  sub9.img 2048 2048
//                  ^D
//
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
  Cimage  *poSmall_image;
  Cimage  *poBig_image;
  int      nDim0, nDim1;
  int     nStat;
  int     nTo[2] = {0, 0};
  int     nFrom[4] = {0, 0, 1024, 1024};
  int     i;
  Cstring  sOut, sName, sComment;

  vDtrekSetModuleName("dtcompose");
  vPrintCopyrightInfo();
  
  sComment = "Made by dtcompose from images: ";

  cout << "dtcompose: Enter output image name and its dimensions: ";
  if (!(cin >> sOut >> nDim0 >> nDim1) || (nDim0 < 0 || nDim1 < 0) ) {
    cout << "dtcompose: ABORTING: invalid image dimensions." << endl;
    return 1;
  }


  int nFirst = 0;
  nStat = 0;
  cout << "dtcompose: Enter input image name and\n";
  cout << "           its origin in output image: ";
  while ( (nStat == 0) && (cin >> sName >> nTo[0] >> nTo[1]) ) 
    {
      poSmall_image = new Cimage (sName);
      if (poSmall_image->bIsAvailable()) nFirst++;
      if (1 == nFirst) 
	{
	  // Very first subimage read in, so create output image and initialize
	  // to all 0s.

	  poBig_image = new Cimage(nDim0, nDim1, poSmall_image->m_eData_type);
      
	  //      (void) poBig_image->nSetNextPixel(0, 0);
	  //      for (i = 0; i < poBig_image->nData_size; i++)
	  //	(void) nPutNextPixel(0);

	  char* data;
	  data = poBig_image->m_The_Data.pc;
	  for (i = 0; i < poBig_image->m_nData_size; i++)
	    *data++ = 0;

	  // Copy header from first image to output image.

	  poBig_image->m_oHeader.nCopyMask(poSmall_image->m_oHeader,
					   Cstring("*"));
	}

      //  Small_image has been read in; Big_image is ready to receive it.

      nFrom[2] = poSmall_image->m_nDim[0];
      nFrom[3] = poSmall_image->m_nDim[1];
      nStat = poBig_image->nCopyRect(*poSmall_image, nFrom, nTo);
      if (nStat == 0)
	{
	  if (sComment.length() < 240) 
	    {
	      sComment = sComment + sName + ", ";
	    }
	  cout << "dtcompose: Enter input image name and\n";
	  cout << "           its origin in output image (^D to finish): ";
	}
      else
	{
	  cout << "dtcompose: Problem adding image " << sName << endl;
	}
      delete poSmall_image;
    }
  if ( (nStat == 0) && (nFirst > 0) ) 
    {

    // Write out image

      (void) poBig_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment,
						  sComment);
      nStat = poBig_image->nWrite(sOut);
    }
  if ( (nStat != 0) && (nFirst > 0) ) 
    cout << "dtcompose: Error writing output image: " << nStat << "!\n"; 
  else
    cout << "dtcompose: Done" << endl;
  delete poBig_image;
  return (nStat);
}
