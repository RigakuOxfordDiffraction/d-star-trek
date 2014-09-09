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
// dtreference.cc     Initial author: J.W. Pflugrath           09-Sep-1995
//    This make a non-uniformity of response reference image from an
//    experimentally obtained image.
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
//    dtreference reads from stdin the filename of the input image 
//    to reference.
//
//    Example:  % dtreference 
//                  input.img
//                  550.2 549.4 0.100 0.100 output.img
//                  ^D
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include <iostream.h>
#include "Cstring.h"
#include "Cimage.h"
#include "Cspatial.h"

// Function prototypes in this source module

float fInterpValue(const float fXmm, const float fYmm, Cimage *poImage);

int main()

{
  Cimage   *poInput_image;
  Cimage   *poTemp_image;
  Cimage   *poOutput_image;
  Cspatial *poSpatialIn;
  Cspatial *poSpatialOut;
  int       nStat;
  float     fSpatial[4];
  int       nDim0, nDim1;
  int i, j;   // Loop counters

  Cstring   sOut, sName, sComment;
  char      cComment[100];

  cerr << "\n\ndtreference: Copyright (c) 2006 Rigaku\n";

  sComment = "Made by dtreference from image: ";

  cout << "dtreference: Enter input image name: ";
  if (!(cin >> sName))
    {
      cerr << "dtreference: ABORTING: invalid name.\n";
      return (1);
    }

  poInput_image = new Cimage (sName);
  if ( !poInput_image->bIsAvailable() )
    {
      cerr << "Error reading input image.\n";
      return (1);
    }

  poSpatialIn = new Cspatial (poInput_image->m_oHeader);
  if (!poSpatialIn->bIsAvailable())
    {
      cerr << "dtreference: error with input image spatial distortion!\n";
      return (1);
    }

  cout << "dtreference: Enter output image dimensions, beam center and pixel sizes: ";
  cin >> nDim0 >> nDim1 >> fSpatial[0] >> fSpatial[1] >> fSpatial[2] >> fSpatial[3];

  poSpatialOut = new Cspatial(fSpatial[0], fSpatial[1], fSpatial[2], fSpatial[3],
			      nDim0, nDim1, 0, -1, 1, 0);

  nStat        = poSpatialOut->nSetMaxPixel(nDim0-1, nDim1-1);

  // Find mm limits of output image and if they exist in the input image.
  float fXYmm[4][2];

  nStat = nStat + poSpatialOut->nPixeltoMM((float) 0.0, (float) 0.0, &fXYmm[0][0], &fXYmm[0][1]);
  nStat = nStat + poSpatialOut->nPixeltoMM((float) 0.0, (float) nDim1-1.0, &fXYmm[1][0], &fXYmm[1][1]);
  nStat = nStat + poSpatialOut->nPixeltoMM((float) nDim0-1.0, (float) nDim1-1.0, &fXYmm[2][0], &fXYmm[2][1]);
  nStat = nStat + poSpatialOut->nPixeltoMM((float) nDim0-1.0, (float) 0.0, &fXYmm[3][0], &fXYmm[3][1]);

  if (0 != nStat)
    {
      cerr << "dtreference: error with output image mm limits!\n";
      return (1);
    }
  
  float f12pix[4][2];
  nStat = nStat + poSpatialIn->nMMtoPixel(fXYmm[0][0], fXYmm[0][1], &f12pix[0][0], &f12pix[0][1]);
  nStat = nStat + poSpatialIn->nMMtoPixel(fXYmm[1][0], fXYmm[1][1], &f12pix[1][0], &f12pix[1][1]);
  nStat = nStat + poSpatialIn->nMMtoPixel(fXYmm[2][0], fXYmm[2][1], &f12pix[2][0], &f12pix[2][1]);
  nStat = nStat + poSpatialIn->nMMtoPixel(fXYmm[3][0], fXYmm[3][1], &f12pix[3][0], &f12pix[3][1]);

  if (0 != nStat)
    {
      cerr << "dtreference: error with input image mm limits!\n";
      return (2);
    }
  
  // Compute average, min and max pixel value in input reference image

  int   nPxArea[4];
  float fAvg, fSD, fMin, fMax;

  // Find the pixel corners

  float fMin0, fMin1, fMax0, fMax1;
  
  fMin0 = f12pix[0][0];
  fMin1 = f12pix[0][1];
  fMax0 = fMin0;
  fMax1 = fMin1;

  for (i = 0; i < 4; i++)
    {
      if (f12pix[i][0] < fMin0)
	{
	  fMin0 = f12pix[i][0];
	}
      if (f12pix[i][0] > fMax0)
	{
	  fMax0 = f12pix[i][0];
	}
      if (f12pix[i][1] < fMin1)
	{
	  fMin1 = f12pix[i][1];
	}
      if (f12pix[i][1] > fMax1)
	{
	  fMax1 = f12pix[i][1];
	}
    }
  nPxArea[0] = (int) (fMin0 - (float) 0.5);
  nPxArea[1] = (int) (fMin1 - (float) 0.5);
  nPxArea[2] = (int) (fMax0 - fMin0 + (float) 0.5);
  nPxArea[3] = (int) (fMax1 - fMin1 + (float) 0.5);

  poInput_image->nAvgSD(nPxArea, (float) 0.0, (float)1.0E20, &fAvg, &fSD, &fMin, &fMax);
    
  // Create output image of proper size and data type

  poTemp_image = new Cimage(nDim0, nDim1, eImage_uI2);

  // Now loop over all pixels in the output image and for each compute:
  // mm coordinate
  // pixel value in input image
  // pixel value scaled between 0 and 32767 (could be 65535)

  unsigned short int uiTemp;
  float fValue;
  float fScale;

  if (fMin == fMax)
    {
      fMin = fMax - 1.0;   // This prevents potential divide by 0.0
    }
  fScale = (fMax - fMin) / 32767.0;
  fScale = 1.0 / fScale;
    
  poTemp_image->nSetNextPixel(0, 0);

  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
	{
	  // Get mm coordinates of the output pixel

	  nStat = poSpatialOut->nPixeltoMM((float) i, (float) j,
					   &fXYmm[0][0], &fXYmm[0][1]);

	  if (0 != nStat)
	    {
	      cerr << "Error in output image mm.\n";
	      return (3);
	    }

	  // Get pixel coordinates of the mm coords in the input image
	  
	  nStat = poSpatialIn->nMMtoPixel(fXYmm[0][0], fXYmm[0][1],
					  &f12pix[0][0], &f12pix[0][1]);

	  if (0 != nStat)
	    {
	      cerr << "Error in input image mm.\n";
	      return (4);
	    }

	  // Next few lines should be replaced with some interpolation
	  // or other function

	  fValue = fInterpValue(f12pix[0][0], f12pix[0][1], poInput_image);
	  if (0.0 >= fValue)
	    {
	      cerr << "Error in input image mm.\n";
	      return (4);
	    }

  //	  fValue = (fValue - fMin) * fScale;
	  uiTemp = (unsigned short int) min(32767.0,fValue);
	  poTemp_image->vSetNextPixel(&uiTemp);
	}
    }

  // Free up some memory since we are not using the input image anymore.

  delete poInput_image;

  // Create the output image by smoothing the temporary image

  poOutput_image = new Cimage(nDim0, nDim1, eImage_uI2);
  nPxArea[2] = 3;
  nPxArea[3] = 3;
  poOutput_image->nSetNextPixel(0, 0);
  for (j = 0; j < nDim1; j++)
    {
      if ( (j == 0) || (j == nDim1-1))
      {
	// First and last lines are special

	for (i = 0; i < nDim0; i++)
	  {
	    uiTemp = (unsigned short int) 
	               (poTemp_image->*poTemp_image->prfGetPixel)(i, j);
	    poOutput_image->vSetNextPixel(&uiTemp);
	  } // end i loop
      }
      else
	{
	  // First pixel in line is transferred
	  uiTemp = (unsigned short int) 
	    (poTemp_image->*poTemp_image->prfGetPixel)(0, j);
	  poOutput_image->vSetNextPixel(&uiTemp);

	  nPxArea[1] = j-1;
	  nPxArea[0] = 0;
	  for (i = 1; i < nDim0-1; i++)
	    {
	      nStat      = poTemp_image->nAvgSD(nPxArea,
						(float) 0.0, (float)1.0E20, 
						&fAvg, &fSD, &fMin, &fMax);
	      // What if nStat is not 0?

	      uiTemp = (unsigned short int) fAvg;
	      poOutput_image->vSetNextPixel(&uiTemp);
	      nPxArea[0]++;
	    }

	  // Last pixel in line is transferred
	  uiTemp = (unsigned short int) 
	    (poTemp_image->*poTemp_image->prfGetPixel)(nDim0-1, j);
	  poOutput_image->vSetNextPixel(&uiTemp);
	}
    }

  if (0 == nStat)
    {
      // Write out the smoothed image

      cout << "dtreference: Enter output image name: ";
      if (!(cin >> sOut))
	{
	  cerr << "dtreference: ABORTING: invalid name.\n";
	  return (5);
	}
      sComment = sComment + sName;
      poOutput_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sComment);
      poSpatialOut->nUpdateHeader(&poInput_image->m_oHeader);
      poOutput_image->m_oHeader.nReplaceValue(Cimage_header::ms_sType, "mad");
      nStat = poOutput_image->nWrite(sOut);
    }
  delete poTemp_image;
  delete poOutput_image;
  return (nStat);
}

float fInterpValue(const float f0pix, const float f1pix, Cimage *poImage)
{
  int nPxArea[4];
  int nStat;
  float fAvg, fSD, fMin, fMax;

  // This should probably be done with mm coordinates and not pixel coordinates

  nPxArea[0] = (int) (f0pix - (float) 1.5);
  nPxArea[1] = (int) (f1pix - (float) 1.5);
  nPxArea[2] = 4;
  nPxArea[3] = 4;
  nStat      = poImage->nAvgSD(nPxArea,
			  (float) 0.0, (float)1.0E20, 
			  &fAvg, &fSD, &fMin, &fMax);
    
  if (0 == nStat)
    {
      return (fAvg);
    }
  return (-999.0);
}
