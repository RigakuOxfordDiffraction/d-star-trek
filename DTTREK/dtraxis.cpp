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
// dtraxis.cc     Initial author: J.W. Pflugrath           13-Sep-1995
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
//    dtraxis reads from stdin the filename of a D*TREK image file
//    to read and the name of an image file to write.  It reads in the input
//    image, reorders it and writes out the result to the output image.
//
//    Example:  % dtraxis
//                  scrambled.img 
//                  unscrambled.img
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

//#include <iostream.h>
#include "Cstring.h"
#include "Cimage.h"
#include "Cspatial.h"
#include "dtrekdefs.h"
#include "raxis.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;
using std::cerr;

int main()

{
  Cimage  *poImage;
  Cstring  sIn, sOut;
  int      nStat;
  int      nDim0, nDim1;
  int      nOutDim0, nOutDim1;
  int      nOffSet0, nOffSet1;
  float    fSatVal;

  cerr << "\n\ndtraxis: Copyright (c) 2006 Rigaku\n"
       << endl;

  cout << "dtraxis: Enter input corrected CCD image name: ";
  if (!(cin >> sIn)) 
    {
      cerr << "dtraxis: ABORTING: invalid image name.\n";
      return (1);
    }

  poImage = new Cimage(sIn);
//  if (!poImage->bIsAvailable) return (-1);
  if (!poImage->bIsAvailable()) 
    {
      cerr << "dtraxis: ABORTING: image " << sIn << " not available.\n";
      return (-1);
    }

  cout << "\ndtraxis: Enter output R-AXIS image name: ";
  if (!(cin >> sOut))
    {
      cerr << "dtraxis: ABORTING: invalid image name." << endl;
      return (1);
    } 

//  cout << "dtrek2raxis: Enter pixel value for bad pixels: ";
//  cin >> fBad;

  // Do some descrambling here

//  Cimage oOutImage (*poImage, NULL);

  nStat = poImage->nGetDimensions(&nDim0, &nDim1);
  if ( (1900 >= nDim0) && (1900 >= nDim1) )
    {
      // Assume a 1536 by 1536 input image
      nOutDim0 = 1900;
      nOutDim1 = 1900;
      // nOffSet0 = (1900 - 1536) / 2;
      // nOffSet1 = (1900 - 1536) / 2;
      nOffSet0 = 0;
      nOffSet1 = 0;
    }
  else
    {
      // Assume a 3072 by 3072 input image

      nOutDim0 = 3000;
      nOutDim1 = 3000;
      nOffSet0 = 0;
      nOffSet1 = 0;
    }

  int i, j, jj;    // Loop counters
  unsigned short int uiTemp1;

  Cimage oOutImage(nOutDim1, nOutDim0, eImage_uI2);
  oOutImage.nSetNextPixel(0, 0);
  for (j = 0; j < nOutDim1; j++)
    {
      for (i = 0; i < nOutDim0; i++)
	{
	  // Zero entire output input
	  oOutImage.vSetNextPixel((unsigned short int) 0);
	}
    }

  // Now transfer input image to output image

  int   nExt0, nExt1;
  float fPixel;
  float fCompressInfo[2];

  // The RAXIS images have compressed pixels if the pixel value is above 32767.
  // So use the compression factor for the specific R-AXIS image (whether R-AXIS II or R-AXIS IV)

  if (1900 < nOutDim0)
    {
      // RAXIS 4
      fCompressInfo[0] = 32.0;
      fCompressInfo[1] = 32767.0;
      }
  else
    {
      // RAXIS 2
      fCompressInfo[0] = 8.0;
      fCompressInfo[1] = 32767.0;
    }

  nExt0 = nDim0;
  nExt1 = nDim1;

  fSatVal = poImage->fGetSatValue();
  cout << "Saturated value: " << fSatVal << '\n';

  poImage->nSetNextPixel(0, 0);
  for (j = 0; j < nDim0; j++)
    {
      jj = nDim1 - j - 1;                           // Reverse slow
      for (i = 0; i < nDim1; i++)
        {
          uiTemp1 = poImage->uiGetNextPixel();   
	  fPixel  = (float) uiTemp1;

	  // Convert to RAXISIV or RAXIS II leave saturated pixels as saturated

	  if (fSatVal <= fPixel)
	    {
	      // Pixel is saturated, make sure stays saturated

	      uiTemp1 = 65535;
	    }
	  else if (32768.0 <= fPixel)
	    {
	      uiTemp1 = (unsigned short int) ((fPixel / fCompressInfo[0]) +
					fCompressInfo[1]);
	    }
	  //	  ii = nOutDim0 - i - 1;                       // Reverse fast
	  //	  ii = nOut0 - i - 1;                       // Reverse fast

	  // Switch directions and apply offsets to center the output image
	  // Note that nSetPixel checks bounds so we cannot go out of bounds

          oOutImage.nSetPixel(    j+nOffSet1, i+nOffSet0, uiTemp1);
        }
    }

  // Write out image, but first adjust header

  oOutImage.m_oHeader = poImage->m_oHeader;
  oOutImage.m_oHeader.nReplaceValue(Cimage_header::ms_sSize1, nOutDim1);
  oOutImage.m_oHeader.nReplaceValue(Cimage_header::ms_sSize2, nOutDim0);
  oOutImage.m_oHeader.nReplaceValue(Cimage_header::ms_sRaxisCompressionRatio,
				    fCompressInfo[0], 3);


  Cstring sPrefix;
  float a4fTemp[4];
  nStat = oOutImage.m_oHeader.nGetValue(Cdetector::ms_sDetectorNames, &sPrefix);
  if (0 != nStat) sPrefix = "";
  
  nStat = oOutImage.m_oHeader.nGetValue(sPrefix + 
					Cspatial::ms_sSpatialDistortionInfo,
					4, a4fTemp);
  // Because the correct program may have written the image, you must
  // delete ALL the SPATIAL_DISTORTION_INFO keywords

  while (0 == oOutImage.m_oHeader.nDelete(sPrefix + 
					  Cspatial::ms_sSpatialDistortionInfo))
    ;
  float fTemp;
  a4fTemp[0] = a4fTemp[0] + nOffSet0;
  a4fTemp[1] = a4fTemp[1] + nOffSet1;
  if (3000 == nOutDim0)
    {
      fTemp      = a4fTemp[0];
      a4fTemp[0] = a4fTemp[1];
      a4fTemp[1] = fTemp;
      a4fTemp[2] = 0.100;
      a4fTemp[3] = 0.100;
    }
  else
    {
      fTemp      = a4fTemp[0];
      a4fTemp[0] = a4fTemp[1];
      a4fTemp[1] = fTemp;
      a4fTemp[2] = 0.101;
      a4fTemp[3] = 0.101;

      // denzo expects R-AXIS II images to be padded to 4096 byte records

      oOutImage.m_oHeader.nReplaceValue(sPrefix + 
					Cimage_header::ms_sRaxisRecordSize,
					(int) 4096);
    }
  oOutImage.m_oHeader.nReplaceValue(sPrefix + 
				    Cspatial::ms_sSpatialDistortionInfo,
				    4, a4fTemp, 4);

  nStat = nWriteRAXIS(oOutImage, sOut);

  if (0 != nStat)
    cerr << "\ndtraxis: Error writing output R-AXIS image: " << nStat << "!\n"; 
  else
    cerr << "\ndtraxis: Done" << endl;
  delete poImage;
  return (nStat);
}
