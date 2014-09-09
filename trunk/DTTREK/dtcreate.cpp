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
// dtcreate.cc     Initial author: J.W. Pflugrath           10-Mar-1995
//    This creates an image with a single value for all pixels.
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
//    dtcreate reads from stdin the filename and dimensions of an image 
//    to create and a constant value for all pixels in the image.
//
//    Example:  % dtcreate
//                  big.img  3072 3072 100
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
  Cimage  *poBig_image;
  int      nDim0, nDim1;
  int     nStat;
  int     i, iValue;
  Cstring  sOut,sComment;

  cout << "\n\ndtcreate: Copyright (c) 2006 Rigaku\n"
       << endl;

  sComment = "Made by dtcreate with value: ";

  cout << "dtcreate: Enter output image name,\n" <<
          "          its dimensions and constant value: ";
  if (!(cin >> sOut >> nDim0 >> nDim1) || (nDim0 < 0 || nDim1 < 0) )
    {
      cout << "dtcreate: ABORTING: invalid image dimensions." << endl;
      return 1;
    }
  if (!(cin >> iValue))
    {
      cout << "dtcreate: Error reading constant value.\n";
      return 3;
    }

  poBig_image = new Cimage(nDim0, nDim1, eImage_I2);
  
  short int* piData;
  piData = poBig_image->m_The_Data.pi;        // This only works with short int!
  for (i = 0; i < (nDim0 * nDim1); i++)
    *piData++ = iValue;

  // Write out image

  (void) poBig_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment,
					      sComment);
  nStat = poBig_image->nWrite(sOut);

  if (nStat != 0)
    cout << "dtcreate: Error writing output image: " << nStat << "!\n"; 
  else
    cout << "dtcreate: Done" << endl;
  delete poBig_image;
  return (nStat);
}
