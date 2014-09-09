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
//   dtsbcdark  takes a dark mask image, a data image and number of data images
//              to use and constructs a dark image for use in correcting
//              data images.
//   Usage: dtsbcdark  dtsbcdark.mask  data001.img 10
//
//+Include files

//#include <iostream.h>
//#include <iomanip.h>
#include "Cimage.h"
#include "Cstring.h"
#include "Cscan.h"

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
  Cstring  sData, sMask;

  Cscan         *poScan;
  Cimage        *poImgData;
  Cimage        *poImgMask;
  Cstring        sTemp;
  int            nStat;
  int            nNumImg;
  unsigned short int uiMaxValue;
  unsigned short int uiMinValue;

  vDtrekSetModuleName("dtsbcdark");
  vPrintCopyrightInfo();

  // Copy command line to output log

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;

  cout << "\ndtsbcdark:  Create a dark image from the average of\n"
       << "            a subset of dark pixels in a set of data images.\n"
       << "Assumption: Each CCD readout area has a single CONSTANT dark value"
       << "\n            for ALL pixels in that area.  There are 18 areas.\n"
       << "            The darkmask image was created by dtsbcdarkmask when\n"
       << "            the detector was calibrated.\n";

  argc--;
  argv++;
  if (3 > argc)
    {
      cout << "\nERROR: not enough input arguments!\n"
           << "Usage: dtsbcdark darkmask_image data_image num_data_imgs_to_average [min_value [max_value]]\n"
           << "       (Use dtsbcdarkmask to create the darkmask_image.)\n" << endl;
      exit (1);
    }
  sMask      = argv[0];
  sData      = argv[1];
  nNumImg    = atoi(argv[2]);
  uiMinValue = 1;
  uiMaxValue = 20000;
  if (3 < argc)
    uiMinValue = (unsigned short int) atoi(argv[3]);
  if (4 < argc)
    uiMaxValue = (unsigned short int) atoi(argv[4]);
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

  poImgMask  = new Cimage(sMask);
  poImgData  = new Cimage(sData);

  if (!poImgMask->bIsAvailable())
    {
      cout << "ERROR, Mask image: " << sMask << ", is unavailable!\n" << endl;
      exit (2);
    }
  if (!poImgData->bIsAvailable())
    {
      cout << "ERROR, Data image: " << sData << ",  is unavailable!\n" << endl;
      exit (2);
    }
  if (0 >= nNumImg)
    {
      cout << "ERROR, num images to use: " << nNumImg << " is bogus!\n" << endl;
      exit (2);
    }

  unsigned short int uiMask, uiData;

  poImgData->nGetDimensions(&nDim0, &nDim1);

  if (    (poImgMask->nGetDimension(0) != nDim0)
       || (poImgMask->nGetDimension(1) != nDim1) )
    {
      cout << "ERROR, mask and data image do not have the same dimensions!\n" 
	   << endl;
      exit (3);
    }

  // Average nNumImg data images to get some averaged dark pixels to use
  // TODO

  cout << "Averaging " << nNumImg << " data images..." << flush;

  poScan = new Cscan(poImgData->m_oHeader);
  poScan->vSetSeqStart(poScan->nGetSeqNum(sData));

  nStat = poScan->nAvgSD(poImgData, nNumImg, NULL, NULL);
  if (bFileExists(sData+"sd"))  // We don't need the sd image that was written
    {
      nFileDelete(sData+"sd");
      cout << "...temporary file  " << sData << "sd  deleted.\n";
    }
      
  if (0 == nStat)
    {
      cout << "...done averaging.\n";
      cout << "Min, Max allowed dark pixel value: " << uiMinValue << ", " << uiMaxValue
           << endl << flush;
    }
  else
    {
      cout << "...ERROR problem with averaging.\n" << flush;
      exit (4);
    }

  // Now build the output dark image

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

  poImgMask->nSetNextPixel(0, 0);
  poImgData->nSetNextPixel(0, 0);
  for (j = 0; j < nDim1; j++)
    {
      jj = j / (nDim1 / 3);
      for (i = 0; i < nDim0; i++)
	{
	  ii = i / (nDim0 / 6);
	  uiMask  = poImgMask->uiGetNextPixel();
	  uiData  = poImgData->uiGetNextPixel();
	  if (   (0 < uiMask)
              && (uiData >= uiMinValue)
	      && (uiData <= uiMaxValue) )
	    {
              // Use only those data_image pixels 
	      //    flagged as dark by the dark mask.

	      a18dSum[ii][jj] += (double)uiData;
	      a18dNum[ii][jj] += 1.0;
	    }
	}
    }

  // Calculate the average dark pixel value for this image readout area
  
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
	      printf("   ERROR, no dark pixels!\n");
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

  // Create the dark image

  poImgData->nSetNextPixel(0, 0);
  for (j = 0; j < nDim1; j++)
    {
      jj = j / (nDim1 / 3);
      for (i = 0; i < nDim0; i++)
	{
	  ii = i / (nDim0 / 6);

	  uiData  = (unsigned short int) (a18dSum[ii][jj]);
	  //!!!uiData +=1;   // I find that we underestimate the dark by 1 ADU
	  poImgData->vSetNextPixel(uiData);
	}
    }

  // The name should be the input data image name with _dark appended!

  nStat = poImgData->nWrite(sData + "_dark");

  return (nStat);
}


