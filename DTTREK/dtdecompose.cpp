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
// dtdecompose.cc     Initial author: J.W. Pflugrath           10-Mar-1995
//    This decomposes a larger image into several images.
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
//    dtdecompose reads from stdin the filename of the input image 
//    to decompose.
//    Next it reads image filenames and the origin and extent
//    (in pixels) of the rectangle of pixels in the input image to be copied
//    into the output image.  It creates the new image and writes it out.
//    The new image has the same Data_type as the input image.
//    When there is an end-of-file on stdin, dtdecompose exits.
//
//    Example:  % dtdecompose 
//                  big.img
//                     0    0 1024 1024 sub1.img
//                  1024    0 1024 1024 sub2.img
//                  2048    0 1024 1024 sub3.img
//                     0 1024 1024 1024 sub4.img
//                  1024 1024 1024 1024 sub5.img
//                  2048 1024 1024 1024 sub6.img
//                     0 2048 1024 1024 sub7.img
//                  1024 2048 1024 1024 sub8.img
//                  2048 2048 1024 1024 sub9.img
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
  Cimage  oBig_image;
  Cimage *poSmall_image;
  int     nStat;
  int     nTo[2] = {0, 0};
  int     nFrom[4];
  Cstring sOut, sName, sComment;
  char    cComment[100];

  vDtrekSetModuleName("dtdecompose");
  vPrintCopyrightInfo();

  sComment = "Made by dtdecompose from image: ";

  cout << "dtdecompose: Enter input image name: ";
  if (!(cin >> sName)) {
    cout << "dtdecompose: ABORTING: invalid name." << endl;
    return 1;
  }

  nStat = oBig_image.nRead(sName);
  if ( (nStat != 0) || (!oBig_image.bIsAvailable()) ) {
    cout << "Error reading input image.\n";
    return (1);
  }
  cout << "INFO: dimensions of the input image are: " << oBig_image.nGetDimension(0) 
       << " by " << oBig_image.nGetDimension(1) << " pixels.\n" << flush; 
//+JWP
  float fSatVal;
  fSatVal = oBig_image.fGetSatValue();
  cout << "INFO: saturated value in the input image is: " << fSatVal << endl;
//-JWP

  while (nStat == 0) {      
    cout << "\ndtdecompose: Enter origin and \n";
    cout << "             extents to extract and output filename: ";
    cin >> nFrom[0] >> nFrom[1] >> nFrom[2] >> nFrom[3] >> sOut;
    if (!cin) {
      cout << "\ndtdecompose: Done" << endl;
      return 0;
    }
    else {
       
      // Construct the "To" image

      poSmall_image = new Cimage(nFrom[2], nFrom[3], oBig_image.m_eData_type);
      nStat = poSmall_image->nCopyRect(oBig_image, nFrom, nTo);
      if (nStat != 0) {
	cout << "dtdecompose: Error copying to small image: " << nStat << endl; 
	return (1);
      }
      else {

	// Copy the header, but reset SIZE1 and SIZE2

	poSmall_image->m_oHeader = oBig_image.m_oHeader;
	poSmall_image->m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, nFrom[2]);
	poSmall_image->m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, nFrom[3]);
	(void) sprintf(cComment, " with origin at %d, %d and extents %d, %d",
                       nFrom[0], nFrom[1], nFrom[2], nFrom[3]);
	sComment = "Subimage made from " + sName + cComment;
	poSmall_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sComment);
	//+JWP 2008-09-03
	// Make sure saturated value is the same internally
	poSmall_image->vSetSatValue(fSatVal);
	//-JWP 2008-09-03

	nStat = poSmall_image->nWrite(sOut);
	if (nStat == 0) 
	  {
	    cout << "dtdecompose: Wrote image " << sOut << cComment << endl;
	  }
	else
	  {
	    cout << "dtdecompose: Error writing output image " << sOut << endl; 
	    return (nStat);
	  }
       delete poSmall_image;
      }
    }
  }
  return (nStat);
}
