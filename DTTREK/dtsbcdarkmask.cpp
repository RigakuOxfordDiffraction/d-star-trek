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
//   dtsbcdarkmask  takes [averaged] dark and flood images and constructs
//                  an image designating which pixels in data image have
//                  no X-ray contributions.  Such pixels can be used to
//                  create a dark image from a set of data images with the
//                  dtsbcdark program.
//   Usage: dtsbcdarkmask dark.dtav flood.dtav min_value max_value
//
//+Include files

#include "Cimage.h"
#include "Cstring.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

//+Code begin

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  int   i, j;
  int nDim0, nDim1;
  Cstring  sFlood, sDark;

  Cimage        *poImgFlood;
  Cimage        *poImgDark;
  Cstring        sTemp;
  int            nStat;
  unsigned short int uiMaxValue;
  unsigned short int uiMinValue;

  vDtrekSetModuleName("dtsbcdarkmask");
  vPrintCopyrightInfo();

  // Copy command line to output log

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;
  cout << "\ndtsbcdarkmask: Select dark pixels in a flood or data image and\n"
       << "               output their locations in a darkmask image.\n\n";

  argc--;
  argv++;

  if (4 > argc)
    {
      cout << "ERROR: not enough input arguments!\n"
           << "Usage:  dtsbcdarkmask flood_image_name dark_image_name min_value max_value\n" << endl;
      exit (1);
    }
  sFlood     = argv[0];
  sDark      = argv[1];
  uiMinValue = (unsigned short int) atoi(argv[2]);
  uiMaxValue = (unsigned short int) atoi(argv[3]);
  poImgFlood = new Cimage(sFlood);
  poImgDark  = new Cimage(sDark);

  if (20000 < uiMaxValue)
    {
      cout << "ERROR, bogus max value: " << uiMaxValue << "!\n" << endl;
      exit (2);
    }
  if (uiMinValue < uiMinValue)
    {
      cout << "ERROR, bogus min value: " << uiMinValue << "!\n" << endl;
      exit (2);
    }

  if (!poImgFlood->bIsAvailable())
    {
      cout << "ERROR, flood image: " << sFlood << ",  is unavailable!\n" << endl;
      exit (2);
    }
  if (!poImgDark->bIsAvailable())
    {
      cout << "ERROR, dark image: " << sDark << ", is unavailable!\n" << endl;
      exit (2);
    }

  unsigned short int uiDark, uiFlood;

  poImgFlood->nGetDimensions(&nDim0, &nDim1);

  // Mark the "fiducial lines" as bad

  int a6BadLines[6];
  if (3072 == nDim0)
    {
      a6BadLines[0] =    0;
      a6BadLines[1] = 1023;
      a6BadLines[2] = 1024;
      a6BadLines[3] = 2047;
      a6BadLines[4] = 2048;
      a6BadLines[5] = 3071;
    }
  else
    {
      a6BadLines[0] =    0;
      a6BadLines[1] =  511;
      a6BadLines[2] =  512;
      a6BadLines[3] = 1023;
      a6BadLines[4] = 1024;
      a6BadLines[5] = 1535;
    }
  uiFlood = 0;
  for (j = 0; j < 6; j++)
    {
      for (i = 0; i < nDim0; i++)
	{
	  poImgFlood->nSetPixel(i, a6BadLines[j], uiFlood);
	}
      for (i = 0; i < nDim1; i++)
	{
	  poImgFlood->nSetPixel(a6BadLines[j], i, uiFlood);
	}
    }

  // Compare the dark and flood images.  Only pixels that are dark in both
  // the flood and dark images are valid and go in the dark mask.

  poImgFlood->nSetNextPixel(0, 0);
  poImgDark->nSetNextPixel(0, 0);
  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
	{
	  uiDark  = poImgDark->uiGetNextPixelNoInc();
	  uiFlood = poImgFlood->uiGetNextPixel();
	  if ( ( (uiDark + 10) < uiFlood) || (0 == uiFlood) 
	       || (uiDark < uiMinValue)
	       || (uiDark > uiMaxValue) )
	    {
	      // Bad pixel, so set to 0

	      uiDark = 0;
	    }
	  poImgDark->vSetNextPixel(uiDark);
	}
    }

  // Count the number of dark pixels in each readout area

  double a18dSum[6][3];
  double a18dNum[6][3];
  int ii, jj;
  for (j = 0; j < 3; j++)
    {
      for (i = 0; i < 6; i++)
	{
	  a18dSum[i][j] = 0.0;
	  a18dNum[i][j] = 0.0;
	}
    }

  poImgDark->nSetNextPixel(0, 0);
  for (j = 0; j < nDim1; j++)
    {
      jj = j / (nDim1 / 3);
      for (i = 0; i < nDim0; i++)
	{
	  ii = i / (nDim0 / 6);
	  uiDark  = poImgDark->uiGetNextPixel();
	  if (0 < uiDark)
	    {
	      a18dSum[ii][jj] += (double)uiDark;
	      a18dNum[ii][jj] += 1.0;
	    }
	}
    }

  printf("\nArea[i][j]  # dark pxls  Average\n"
         "--------------------------------\n");
  for (j = 0; j < 3; j++)
    {
      for (i = 0; i < 6; i++)
	{
	  printf("    %2d %2d   %8d", 
		 i, j, (int)a18dNum[i][j]);
	  if (0.0 >= a18dNum[i][j])
	    {
	      nStat = 1;
	      printf("    ERROR, no dark pixels!\n");
	    }
	  else
	    {
	      a18dSum[i][j] /= a18dNum[i][j];
	      a18dSum[i][j] += 0.5;   // Round-up!
	      printf("   %9d\n", (int)a18dSum[i][j]);
	    }
	}
    }
  printf("--------------------------------\n\n");


  sDark = (Cstring)"darkmask" + Cstring(nDim0) + ".img";
  nStat = poImgDark->nWrite(sDark);

  delete poImgDark;
  delete poImgFlood;

  return (nStat);
}






