//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// tCimage.cc            Initial author: J.W. Pflugrath           03-Mar-1995
//    This file tests the class Cimage.
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
//  
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <iostream.h>
#include "Cimage.h"

//+Code begin

int main()
{
  Cimage array;
  Cimage darray[10];
  Cimage array2("test.img");
  int nDim1  = 1024;
  short int *piInt;
  int i;
  int nWhere[2];

  array.oHeader.nList();
  array2.oHeader.nList();

// Test nGetPixel, nPutPixel;
  (void) array2.nGetDimensions(&nWhere[0], &nWhere[1]);
  nWhere[0] = nWhere[0] / 2;
  nWhere[1] = nWhere[1] / 2;
  short int iPixel;
  (void) array2.nGetPixel(nWhere, &iPixel);
  cout << "Pixel at " << nWhere[0] << ", " << nWhere[1] << " is " 
       <<  iPixel << ".\n";
  
  iPixel = 9999;
  (void) array2.nSetPixel(nWhere, iPixel);
  (void) array2.nGetPixel(nWhere, &iPixel);
  cout << "Pixel at " << nWhere[0] << ", " << nWhere[1] << " is now " 
       <<  iPixel << ".\n\n";

  (void) array2.nSetNextPixel(nWhere);
  for (i=0; i < 10; i++)
      cout << "Next pixel is : " << array2.iGetNextPixel() << endl;

  cout << endl;
  (void) array2.nSetNextPixel(nWhere);
  for (i=0; i < 10; i++) array2.vSetNextPixel(iPixel);

  (void) array2.nSetNextPixel(nWhere);
  for (i=0; i < 10; i++)
    cout << "Next pixel is : " << array2.iGetNextPixel() << endl;


//  for (i=0; i< 10; i++) {
//    nDim1 = darray[i].nRead("test.img");
//    cout << "Status from darray[i].nRead: " << nDim1 << ".\n";
//  }

  nDim1 = array.nRead("test.img");
  cout << "Status from array.nRead: " << nDim1 << ".\n";
  array.oHeader.nList();

  cout << "Try a read without a file name!\n";
  nDim1 = array2.nRead();

  cout << "Status from array2.nRead: " << nDim1 << ".\n";

// Test Cimage::nList

  int nRect[4];// = {200, 300, 5, 10};   // fast, slow start; fast, slow extents
  do {
    cout << "Input center, widths of rectangle to list: ";
    cin >> nRect[0] >> nRect[1] >> nRect[2] >> nRect[3];
    array2.nList(nRect);
  } while (nRect[3] != 0);

// Test prfGetPixel

  (void) array2.nGetDimensions(&nWhere[0], &nWhere[1]);
  nWhere[0] = nWhere[0] / 2;
  nWhere[1] = nWhere[1] / 2;
  float fTemp = (array2.*array2.prfGetPixel)(nWhere[0]+i, nWhere[1]+i);
  for (i = 0; i < 10; i++)
    {
      cout << "Pixel at " << nWhere[0]+i << ", " << nWhere[1]+i << " is " 
       <<  (array2.*array2.prfGetPixel)(nWhere[0]+i, nWhere[1]+i) << ".\n\n";
    }
  
// Test Cimage::nGetRect, 3rd ctor and nWrite

  do {
    cout << "Input origin and extents of subimage to copy to new image: ";
    cin >> nRect[0] >> nRect[1] >> nRect[2] >> nRect[3];
//    short int *piShort = new (short int) [nRect[2]*nRect[3]];
    short *piShort = new short[nRect[2]*nRect[3]];
    nDim1 = array2.nGetRect(nRect, (char **)&piShort);

    if (nDim1 == 0) {
      Cimage *image1 = new Cimage(nRect[2], nRect[3], eImage_I2, piShort);
      cout << "\n\nListing of image1.oHeader\n";  image1->oHeader.nList();
      nDim1 = image1->nWrite();
      delete image1;
    }
    else
      cout << "Error in nGetRect\n";
    delete [] piShort;
  } while (nRect[3] != 0);

  {
    Cimage *image2 = new Cimage("test.img");
    cout << "\n\nListing of image2.oHeader\n";  image2->oHeader.nList();
    delete image2;
  }
  { 
    Cimage *image3 = new Cimage("");
    cout << "\n\nListing of image3.oHeader\n";  image3->oHeader.nList();
    delete image3;
  }
  array2.oHeader.nList();

  array.nWrite("writetest.img");
  delete piInt;

  short int* piData;
  piData = array.The_Data.pi;
  for (i = 0; i < array.nGetDimension(); i++)
    *piData++ = 100;
  array.nWrite();

  return 0;
}

