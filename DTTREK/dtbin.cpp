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
// dtbin.cc     Initial author: J.W. Pflugrath           10-Mar-1995
//    This bins pixels in an input image to create an output image.
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
//    dtbin reads from stdin the filename of an image and the number of
//    pixels to bin in the fast and the slow direction. It reads the image
//    then creates a new image by binning the pixels, then reads from stdin
//    the name of the image to write out and does so.
//   
//
//    Example:  % dtbin
//                  unbinned.img  2 2 binned.img
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include "Cstring.h"
#include "Cimage.h"
#include "Cdetector.h"

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
#endif

int main()
{
  Cimage  *poBig_image, *poBinned_image;
  int      nDim0, nDim1;
  int      nBigDim0, nBigDim1;
  int      nBin0, nBin1;
  int      nStat;
  int      i, j, k, l;
  bool     bSaturated;
  Cstring  sOut,sComment;
  unsigned short int  uiFrom;
  unsigned short int  uiBinned;
  unsigned short int  uiSatValue;

  cout << "\n\ndtbin: Copyright (c) 2006 Rigaku\n"
       << endl;

  cout << "dtbin: Enter input image name,\n" <<
          "          and its binning factors: ";
  if (!(cin >> sOut >> nBin0 >> nBin1) || (1 > nBin0 || 1 > nBin1) )
    {
      cout << "dtbin: ABORTING: invalid binning factors." << endl;
      return 1;
    }

  poBig_image = new Cimage(sOut);
  
  if (!poBig_image->bIsAvailable())
    {
      cout << "dtbin: ABORTING: invalid input image.\n";
      return 1;
    }

  poBig_image->nGetDimensions(&nBigDim0, &nBigDim1);

  nDim0 = nBigDim0 / nBin0;
  nDim1 = nBigDim1 / nBin1;
  if (   ( (nBin0 * nDim0) != nBigDim0)
      || ( (nBin1 * nDim1) != nBigDim1) )
    {
      nBigDim0 = nBin0 * nDim0;
      nBigDim1 = nBin1 * nDim1;
      cout << "\ndtbin: WARNING binning factors preclude using entire input image!\n"
           << "         Output image will be " << nDim0 << " by " << nDim1
           << " pixels in size.\n" << flush;
    }

  poBinned_image = new Cimage(nDim0, nDim1, eImage_uI2);

  sComment = "Made by dtbin by binning " + sOut;

  uiSatValue = (unsigned short int) poBig_image->fGetSatValue();

  poBinned_image->nSetNextPixel(0, 0);

  float fSum;
  float fNum_pixels = nBin0 * nBin1;
  
  for (i = 0; i < nBigDim1; i = i + nBin1)
    {
      for (j = 0; j < nBigDim0; j = j + nBin0)
	{
	  fSum = 0.0;
	  bSaturated = FALSE;
	  for (k = 0; k < nBin1; k++)
	    {
	      for (l = 0; l < nBin0; l++)
		{
		  // cout << "i, j, k, l, i+k, j+l : " << i << j << k << l << i+k << j+l << endl;
		  nStat = poBig_image->nGetPixel(j+l, i+k, &uiFrom);
		  if (0 != nStat)
		    {
		      cout << "\ndtbin: Error getting pixel from image!\n";
		      return (4);
		    }
		  if (uiFrom < uiSatValue)
		    fSum = fSum + (float) uiFrom;
		  else
		    bSaturated = TRUE;
		}
	    }
	  if (!bSaturated)
	    uiBinned = (unsigned short int) (fSum / fNum_pixels);
	  else
	    uiBinned = uiSatValue;
	  poBinned_image->vSetNextPixel(uiBinned);
	}
    }

  // Write out image

  int	   a2nDetDim[2];
  float	   a2fBeamPos[4];
  Cstring  sTemp;

  // Copy input header to output header

  poBinned_image->m_oHeader = poBig_image->m_oHeader; 

  // Modify the header

  Cstring sDetName = "";
  (void) poBinned_image->m_oHeader.nGetValue(Cdetector::ms_sDetectorNames, 
					   1, &sDetName);

  nStat = poBinned_image->m_oHeader.nGetValue(sDetName + 
					    Cdetector::ms_sDetectorDimensions,
					    2, a2nDetDim);
  if (0 == nStat)
    {
      // Found the keyword, so adjust dimensions according to binning factors

      a2nDetDim[0] = a2nDetDim[0] / nBin0;
      a2nDetDim[1] = a2nDetDim[1] / nBin1;
      nStat = poBinned_image->m_oHeader.nReplaceValue(sDetName + 
					    Cdetector::ms_sDetectorDimensions, 
						    2, a2nDetDim);
    }

  // Delete active map bitmap keywords

  (void) poBinned_image->m_oHeader.nDelete(Cstring(D_K_BitmapSize));
  (void) poBinned_image->m_oHeader.nDelete(Cstring(D_K_BitmapType));

/******
  (void) poBinned_image->m_oHeader.nGetValue(sDetName +
						 Cspatial::ms_sSpatialBeamPosn,
						 &sTemp);
  cout << "Beamposn>>" << sTemp << "<<" << endl;
  nStat = poBinned_image->m_oHeader.nGetValue(sDetName +
					      Cspatial::ms_sSpatialDistortionInfo,
					      &sTemp);
  cout << "DistInfo>>" << sTemp << "<<" << endl;
      
****/
  nStat = poBinned_image->m_oHeader.nGetValue(sDetName +
					      Cspatial::ms_sSpatialBeamPosn,
					      2, a2fBeamPos);
  if (0 == nStat)
    {
      // Found the keyword, so adjust direct beam position

      a2fBeamPos[0] = a2fBeamPos[0] / (float)nBin0;
      a2fBeamPos[1] = a2fBeamPos[1] / (float)nBin1;

      (void) poBinned_image->m_oHeader.nReplaceValue(sDetName +
						     Cspatial::ms_sSpatialBeamPosn,
						     2, a2fBeamPos, 1);
    }

  // Look for spatial distortion info and 
  // change the pixel size if required as well.

  sTemp = "";
  nStat = poBinned_image->m_oHeader.nGetValue(sDetName +
		Cspatial::ms_sSpatialDistortionType, &sTemp);
  if (0 == nStat)
    {
      if (Cspatial::ms_sSpatialTypeSimple == sTemp)
	{
	  nStat = poBinned_image->m_oHeader.nGetValue(sDetName +
					      Cspatial::ms_sSpatialDistortionInfo,
					      4, a2fBeamPos);
	  if (0 == nStat)
	    {
	      // Found the keyword, so adjust direct beam position

	      a2fBeamPos[0] = a2fBeamPos[0] / (float)nBin0;
	      a2fBeamPos[1] = a2fBeamPos[1] / (float)nBin1;
	      a2fBeamPos[2] = a2fBeamPos[2] * (float)nBin0;
	      a2fBeamPos[3] = a2fBeamPos[3] * (float)nBin1;

	      (void) poBinned_image->m_oHeader.nReplaceValue(sDetName +
						     Cspatial::ms_sSpatialDistortionInfo,
						     4, a2fBeamPos, 5);
	    }
	}
    }

  (void) poBinned_image->m_oHeader.nReplaceValue(Cimage_header::ms_sComment,sComment);

  cout << "dtbin: Enter output image name: ";

  if ( !(cin >> sOut) )
    {
      cout << "\ndtbin: Error reading output filename!\n";
    }
  nStat = poBinned_image->nWrite(sOut);

  if (0 != nStat)
    cout << "dtbin: Error writing output image: " << nStat << "!\n"; 
  else
    cout << "dtbin: Done.\n";
  delete poBig_image;
  delete poBinned_image;

  return (nStat);
}
