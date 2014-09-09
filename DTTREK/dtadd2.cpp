// Jiffy that adds two images together and writes a third image of the sum
//
// Copyright (c) 2006 Rigaku
//
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
 
//+Include files

#include  "Cstring.h"

#include "Cimage.h"
#include "Cscan.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

//+Code begin

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  int     i, j;
  int     nStat;
  Cstring sImage1, sImage2, sImage3;

  Cimage *poImage1;
  Cimage *poImage2;
  Cimage *poImage3;

  vDtrekSetModuleName("dtadd2");
  cout << "\ndtadd2:  Copyright (c) 2006 Rigaku\n";
  cout << D_K_DTREKVersion << endl << flush;
  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;
  // Parse command line arguments
  
  argc--; argv++;

  nStat = 1;

  if (3 > argc)
    {
      cout << "ERROR - not enough input image filenames!\n";
      cout << "Usage: dtadd in_image1 in_image2 out_sum_image\n";
      return (1);
    }

  sImage1 = (const char*) argv[0];
  sImage2 = (const char*) argv[1];
  sImage3 = (const char*) argv[2];

  poImage1 = new Cimage(sImage1);
  if (!poImage1->bIsAvailable())
    {
      cout << "ERROR - image " << sImage1 << " not available!\n";
      delete poImage1;
      return (2);
    }
  poImage2 = new Cimage(sImage2);
  if (!poImage2->bIsAvailable())
    {
      cout << "ERROR - image " << sImage2 << " not available!\n";
      delete poImage1;
      delete poImage2;
      return (2);
    }

  float f1, f2, fValue, fTime;
  float fSatValue;

  int   nNumPixels;
  int   nDim0, nDim1;
  int   nSat1, nSat2, nSat3;

  nSat1 = nSat2 = nSat3 = 0;
  nDim0 = poImage1->nGetDimension(0);
  nDim1 = poImage1->nGetDimension(1);

  fSatValue = poImage1->fGetSatValue();

  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
	{
	  f1 =  (poImage1->*poImage1->prfGetPixel)(i, j);
	  if (f1 > fSatValue) 
	    {
	      nSat1++;
	    }
	  f2 =  (poImage2->*poImage2->prfGetPixel)(i, j);
	  if (f2 > fSatValue) 
	    {
	      nSat2++;
	    }

	  // Add the two values together

	  f1 = f1 + f2;
	  if (f1 > fSatValue) 
	    {
	      f1 = fSatValue;
	      nSat3++;
	    }
	  (poImage1->*poImage1->prnSetPixel)(i, j, f1);
	}
    }

  cout << "Num saturated in image 1, 2 and output: " << nSat1 << ", "
       << nSat2 << ", " << nSat3 << endl;

  // Change scan end, change time, change increment

  Crotation *poRotation;
  poRotation = new Crotation (poImage2->m_oHeader);
  if (poRotation->bIsAvailable())
    {
      f1    = poRotation->fGetRotEnd();
      fTime = poRotation->fGetExposureTime();
    }
  delete poRotation;
  poRotation = new Crotation (poImage1->m_oHeader);
  if (poRotation->bIsAvailable())
    {
      // Only set the RotEnd value if it is greater than the RotStart value
      float fRotStart;
      fRotStart = poRotation->fGetRotStart();
      if (f1 > fRotStart)
	{
	  poRotation->vSetRotEnd(f1); 
	  poRotation->vSetIncrement(f1 - fRotStart);
	}
      // Adjust the exposure time, too
      fTime = fTime + poRotation->fGetExposureTime();
      poRotation->vSetExposureTime(fTime);

      poRotation->nUpdateHeader(&poImage1->m_oHeader, (Cstring)"");      
    }
  delete poRotation;

  nStat = poImage1->nWrite(sImage3);

  delete poImage1;
  delete poImage2;
  return (nStat);
}

