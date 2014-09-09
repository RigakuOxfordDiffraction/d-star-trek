//
// Copyright (c) 2007 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainBeamMask.cpp       Initial author: RB          14-Nov-2007
// This file contains the member functions of class CDTMainBeamMask

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
//    dtbeammask reads an image file, determines where the direct beam
//    would impinge on the image and then uses an edge detection and 
//    thresholding algorithm to set pixel values contiguous with the beam
//    position to a value of 0.  It writes out the modified image.
//    Usage: dtbeammask in_image_filename [-nobeam] [-opt1 ... [-optN ...]] out_image_filename
//
//    WARNING this source code is lifted from CXdisplay.cc, so it is
//    duplicated.  The code here and the code in CXdisplay IS DIFFERENT.
//    Yes, I know this is not good programming practice. -jwp
//
//+ToDo
//
//   Error messages need implementing
//

#include "DTMainBeamMask.h"
#include "Cstring.h"
#include "Cimage.h"
#include "Csource.h"
#include "Cdetector.h"

using std::cout;
using std::endl;
using std::flush;

CDTMainBeamMask::CDTMainBeamMask()
{
    m_nRecurseCount = 0;
    m_nShadowPixCount = 0;
    m_fEraseValue = 0.0;
}

int CDTMainBeamMask::nExecute(unsigned int argc,     // Command line argument count
                             char *argv[]) // Pointers to command line args
{
  Cimage  *poImage = NULL;
  int     nStat = 0;
  Cstring sIn, sOut, sComment;
  char    cComment[100];

  cout << "\n\ndtbeammask: Copyright (c) 2006 Rigaku\n";
  cout << D_K_DTREKVersion << endl << flush;
  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;
  bool bNeedHelp = FALSE;
  if (1 < argc)
    {
      sOut = Cstring(argv[1]);
      bNeedHelp = ( (sOut == Cstring("-h")) || (sOut == Cstring("-help")) );
    }

    if ((3 > argc) || bNeedHelp)
    {
      cout << "dtbeammask creates a mask file for use with d*TREK by setting pixel\n"
           << "           values in the beamstop shadow to ERASEVALUE either automatically\n"
           << "           or explicitly.  It reads an image file, determines where the direct\n"
           << "           beam would impinge on the image and then uses an edge detection\n"
           << "           and thresholding algorithm to set pixel values contiguous with the\n"
           << "           position to ERASEVALUE.  Command line options allow one to set to \n"
           << "           ERASEVALUE explicitly pixels, lines, rows, circles, and/or regions.\n"
           << "           It writes out a modified image.  dtbeammask is useful in scripts.\n\n"

           << "Usage: dtbeammask in_image_file [-nobeam] [-opt1 ... [-optN ...]] out_image_file\n" 
           << "  where -opt1, -opt2, ..., -optN is one of\n"
           << "   -erasevalue ERASEVALUE (default is 0 or -1 depending on image pixel format)\n"
           << "   -seedpoint, -pixel, -line, -row, -column, -circle,\n"
           << "   -outsidecircle, -insidecircle and ... are the arguments to the option.\n"
           << "  Examples:\n"
           << "   dtbeammask in01.img -erasevalue 0 -pixel 100 100 -erasevalue 200 \\\n"
           << "              -row 200 0 512 out001.img\n"
           << "   dtbeammask in01.img -col 200 512 1023 -circle 512 512 45 out001.img\n"
           << "   dtbeammask in01.img -rect 490 0 30 512 -outsidecircle 512 512 500 out001.img\n\n";
      return (-1);
    }

  sOut = Cstring(argv[argc-1]);
  cout << "INFO output image is: " << sOut << endl;
  sIn  = Cstring(argv[1]);
  sComment = "Made by dtbeammask from image: ";

  poImage = new Cimage(sIn);
  if (NULL == poImage)
    {
      cout << "Error reading input image.\n";
      return (1);
    }
  else if (!poImage->bIsAvailable())
    {
      cout << "Error reading input image.\n";
      return (1);
    }
  float fBeamPx0, fBeamPx1;
  Cstring sArg, sLine;  
  int nArg;
  nArg = 2;
  argc = argc - 1; // Last arg is filename and already parsed.

  if (eImage_I4 == poImage->nGetDataType())
    {
      m_fEraseValue = -1.0;
    }

  sArg = argv[nArg];
  if ("-nobeam" == sArg) 
    {
      // Do nothing
      nArg++;
    }
  else
    {
      // Start with the predicted beam center
      // ONLY IF -nobeam not the FIRST option

      // Look ahead and see if -erasevalue is the option right after the image name

      if (2 < argc)
	{
	  // Do not increment argc, so this will be parsed again below

	  sArg = Cstring(argv[2]);
	  if ("-erasevalue" == sArg) 
	    {
	      m_fEraseValue = atof(argv[3]);
	      cout << "INFO: ERASEVLAUE set to " << m_fEraseValue << endl;
	    }
        }

      // Get beam center
      float fCutoffValue = 0.0;

      m_nRecurseCount   = 0;
      m_nShadowPixCount = 0;
      nStat = nCDBeamstopShadow(poImage, &fBeamPx0, &fBeamPx1);
      cout << "Direct beam is at pixel " << fBeamPx0 << ", " << fBeamPx1 << endl;
      if (0 == nStat)
        nStat = nCalcDrawBeamstopShadow(fBeamPx0, fBeamPx1,
                                        fCutoffValue, poImage);

    }

  // Now look at other command line options to edit the
  // nEditPixels(Cstring sLine, float fDefaultPixVal, Cimage *poImage)

  // NOTE: -seedpoint and -pixel ARE DIFFERENT!

  // TODO: This should understand -circle, -pixel, -rect options

  while ( (0 == nStat) && (nArg < argc) )
    {
      sArg = Cstring(argv[nArg]);
      if ("-seedpoint" == sArg) 
        {
          nArg++;
          fBeamPx0 = atof(argv[nArg++]);
          fBeamPx1 = atof(argv[nArg++]);
          cout << "Seed point: " 
               << fBeamPx0 << " " << fBeamPx1 << endl;

          //nStat = nCalcDrawBeamstopShadow(fBeamPx0, fBeamPx1,
          // fCutoffValue, poImage);
        }
      else if ("-erasevalue" == sArg) 
        {
          nArg++;
          m_fEraseValue = atof(argv[nArg++]);

          //nStat = nCalcDrawBeamstopShadow(fBeamPx0, fBeamPx1,
          // fCutoffValue, poImage);
        }
      else if ('-' == sArg.GetAt(0))
        {
          // -pixel, -rect, -line, -row, -circle, ....
          //
          // Construct the sLine text from the argv

          sLine = sArg.after(0);
          nArg++;
          while (nArg < argc)
            {
              sArg = Cstring(argv[nArg]);
              if (sArg.bIsNumeric())
                {
                  sLine = sLine + ' ' + sArg;
                  nArg++;
                }
              else
                {
                  break;  // Leave the while loop
                }
            }

          cout << "Edit string: " << sLine << endl;
          cout << "Erase value used: " << m_fEraseValue << endl;

          nStat = nEditPixels(sLine, m_fEraseValue, poImage);
          if (0 != nStat) cout << "ERROR in dtbeammask, problem with edit.\n" << endl;
        }
      else 
        {
          cout << "ERROR in dtbeammask, UNKNOWN command line argument: " 
               << sArg << endl;
          nStat = 1;
        }
    }

  if (0 == nStat)
    {
      if (eImage_I4 == poImage->nGetDataType())
	{
	  cout << "      Attempting conversion to unsigned short int ...\n";
	  
	  //+JWP 2008-02-07
	  nStat = poImage->nConvertIntToUIShort();
	  if (0 == nStat)
	    cout << "      ... conversion successful.\n";
	  else
	    cout << "ERROR: Conversion UNSUCCESSFUL.\n";
	  //-JWP 2008-02-07
	}
      nStat = poImage->nWrite(sOut);
    }
  delete poImage;
  poImage = NULL;
  return (nStat);
}


int CDTMainBeamMask::nCalcDrawBeamstopShadow(float fPixIn0,
                                            float fPixIn1,
                                            float fCutoffValue,
                                            Cimage *poImage)
{
  // Calculate and Draw (i.e. mark) the pixels 
  // that are in the beam stop shadow

  int nStat1;

  float fBeamPx0, fBeamPx1;

  //m_nBeamstopShadowMode = 0;
  fBeamPx0 = fPixIn0;
  fBeamPx1 = fPixIn1;

  int nx, ny, nStat;
  int nWidth = 6;

  // Algorithm
  // Find low pixel values that are connected through other low pixel values
  // to the input pixel position.

  // 1. Start at beam center pixel.
  // 2. Average NxN pixel area and get average and SD.
  // 3. Determine a cutoff value; this cutoff value will change as 
  //    resolution changes.
  // 4. Set pixel values less than or equal to this cutoff to blue.

  if (-1.0 < fBeamPx0)
    {
      // Beam center on the image
      int nEdgeSize = 5;
      C3Ddata *po3Ddata;
      po3Ddata = new C3Ddata(nEdgeSize, nEdgeSize, 1, nint(fBeamPx0)-3, nint(fBeamPx1)-3, 0);
      nStat1 = po3Ddata->nFill2D(poImage, -1);

      float a3fBkg[3], fAvg, fSD, fMaxValue, a3fPeakCent[3], a3fTotCent[3];
      nStat1 = po3Ddata->nBackgroundAvgSD(0, &a3fBkg[0], &fAvg, &fSD, 
                                          &fMaxValue,
                                          &a3fPeakCent[0], &a3fTotCent[0]);

      delete po3Ddata;
      //cout << "Average is: " << fAvg << endl;

      float fThreshold = max(fAvg + 3.0 * fSD, 3.0 * fAvg);
      //cout << "Enter new threshold: ";
      //cin >> fThreshold;
      //fThreshold = m_tImageProps.fScaleMax;
      if (0.0f >= fCutoffValue)
        {
          cout << "Threshold is " << fThreshold << endl << flush;
        }

      cout << "INFO: Automatically trying to flag beamstop shadow pixels...\n"
           << flush;
      // Now 
      int nCent0, nCent1;
      nCent0 = nint(fBeamPx0);
      nCent1 = nint(fBeamPx1);
      
      m_nRecurseCount   = 0;
      m_nShadowPixCount = 0;
      float fSlope = 50.0;
      Cimage oImage(*poImage, poImage->pcTheData());
      nEdgeSize = 5;
      nStat1 = nRecurseShadow(nCent0, nCent1, nEdgeSize, fAvg, fSD, fThreshold,
                              fSlope, &oImage, poImage);

      cout << "INFO: Number of pixels added to the beamstop shadow: "
           <<  m_nShadowPixCount << endl << flush;
    }
  //vRefresh();
  return (0);
}

int CDTMainBeamMask::nCDBeamstopShadow(Cimage *poImage, float *pfBeamPx0, float *pfBeamPx1)
{
  // This should probably go in the Cimage class.

  // Algorithm:
  // 1. Predict pixel position of beam center.
  // 2. Computer a radius of 1/6 the edge of the image.
  // 3. At that radius, for angle = 0 to 360, step 10, check an area of 
  //    5 x 5 pixels and get an average background value.
  // 4. Find the angle of the lowest average pixel value.
  // 5. Use a simplex refinement to go to lowest average angle.
  // 6. Repeat steps 3-5 for a radius of 1/4 the edge of the image.
  // 7. Connect the two points for a line.
  // 8. Move out from beam center along the line and do a perpendicular
  //    line at each pixel, determine
  // 9. Set pixel values along the line below a threshold to 0.

  int i, j;
  int nStat1, nStat2;
  *pfBeamPx0 = 0.0;
  *pfBeamPx1 = 0.0;

  Cimage_header oHeader;
  oHeader = poImage->oGetHeader();

  Csource   *poSource   = NULL;
  Cdetector *poDetector = NULL;
  poDetector = new Cdetector(oHeader);
  poSource   = new Csource(oHeader);

  if (   !poSource->bIsAvailable() 
      || !poDetector->bIsAvailable() )
    {
      cout << "ERROR something wrong with image header.\n";
      delete poSource;
      delete poDetector;
      return (1);
    }

  float fPx0, fPx1, fMm0, fMm1;

  // Get current beam position (if on detector)

  float a3fS0[3];
  float a3fX[3], a3fXR[3];
  float a3x3fMat[3][3];
  a3fX[0] = 0.0;
  a3fX[1] = 0.0;
  a3fX[2] = 0.0;
  float fBeamPx0, fBeamPx1;
  poSource->vCalcGetS0(a3fS0);
  nStat1 = poDetector->nCalcDetCoords(a3fX, a3fS0,
                                      &fMm0, &fMm1, &fBeamPx0, &fBeamPx1);
  if (0 != nStat1)
    return (1);

  *pfBeamPx0 = fBeamPx0;
  *pfBeamPx1 = fBeamPx1;
  
  delete poSource;
  delete poDetector;
  return (0);
}

float CDTMainBeamMask::fLineShadow(float fPx0, 
                                  float fPx1, 
                                  float fSlope, 
                                  float fIntercept, 
                                      float fWidth,
                                  int nNumPoints,
                                  int nSmooth, 
                                  int nDir,
                                      Cimage *poImage)
{
  // TODO: If Abs(fSlope) > 1.0, then loop over pixels a different way
  //       to find the diffuse edge

  // JWP:  Maybe this routine can be used to get the threshold value
  //       in a local area for the other routine?
  //
  // Given a line: Y = fSlope * X + fIntercept, 
  // a point on/near the line: fPx0, fPx1
  // and a direction to go in: nDir
  // Do the following:
  // 1. Copy pixel values in an image to a temporary variable
  // 2. Smooth the values of the points
  // 3. Compute 1st and 2nd derivatives of the values
  // 4. Determine a threshold for the values, where point at fPx0, fPx1
  //    is considered low, and something at the end is considered high.
  // 5. Set pixel values in the original image to 0 from fPx0, fPx1 to
  //    to cutoff point.

  // Cimage is poImage

  int nx, ny;
  int np, nxx;
  int nDim0, nDim1;
  float *pfValues;
  float *pfSmoothValues;
  float *pfSmSortedValues;
  float *pfFirstDeriv;
  float *pfSecondDeriv;

  pfValues       = new float [nNumPoints+1];
  pfSmoothValues = new float [nNumPoints+1];
  pfSmSortedValues = new float [nNumPoints+1];
  pfFirstDeriv   = new float [nNumPoints+1];
  pfSecondDeriv  = new float [nNumPoints+1];

  // 1. Copy pixel values in an image to pfValues

  nDim0 = poImage->nGetDimension(0);
  nDim1 = poImage->nGetDimension(1);

  if (0 == nDir) nDir  = 1;

  // 1. Copy nNumPoints image values into pfValues;

  // Find limits within the image for 

  nx    = nint(fPx0);
  ny    = nint(fPx1);
  int nw;
  float fSum;
  int nCount;
  int nWidth = nint(fWidth);
  if (0.0 == fSlope)
    {
      for (np = 0; np < nNumPoints; np++)
        {
          fSum = 0.0;
          nCount = 0;
          pfValues[np] = 0.0;
          // ny does not change
          ny = nint(fPx1);
          if ( (0 <= nx) && (nDim0 > nx) )
            {
              for (nw = 0; nw < nWidth; nw++)
                {
                  if ( (0 <= ny) && (nDim1 > ny) )
                    {
                      fSum += (poImage->*poImage->prfGetPixel)(nx, ny);
                      nCount++;
                    }
                  ny++; // Iterate over width which changes the Y intercept
                }
            }
          if (0 < nCount) pfValues[np] = fSum / (float) nCount;
          nx += nDir;
        }
    }
  else // Look for an edge along a vertical line, so look along ny and not nx
    {
      for (np = 0; np < nNumPoints; np++)
        {
          fSum = 0.0;
          nCount = 0;
          pfValues[np] = 0.0;
          nx = nint(fPx0);
          if ( (0 <= ny) && (nDim1 > ny) )
            {
              for (nw = 0; nw < nWidth; nw++)
                {
                  if ( (0 <= nx) && (nDim0 > nx) )
                    {
                      fSum += (poImage->*poImage->prfGetPixel)(nx, ny);
                      nCount++;
                    }
                nx++; // Iterate over width which changes the X intercept
                }
            }
          if (0 < nCount) pfValues[np] = fSum / (float) nCount;
          ny += nDir;
        }
    }

  // 2. Smooth the pfValues array into pfSmoothValues;
  //    Use 1/20th of the total points, but 5 points minimum for smoothing
  
  int nSmPoints = 0;
  nCount = 0;
  fSum = 0.0f;
  
  nSmPoints = max(5, nNumPoints / 20);
  for (np = 0; np < nNumPoints; np++)
    {
      fSum   = 0.0f;
      nCount = 0;
      pfSmoothValues[np] = 0.0;
      for (nx = max(0, (np - nSmPoints /2)); 
           nx < min(nNumPoints, np + (nSmPoints/2)); nx++)
        {
          if (0.0f < pfValues[nx])
            {
              nCount++;
              fSum += pfValues[nx];
            }
        }
      if (nCount > 0) pfSmoothValues[np] = fSum / (float) nCount;
      pfSmSortedValues[np] = pfSmoothValues[np];
    }

  // Find the median smoothed pixel value 

  qsort(&pfSmSortedValues[0], nNumPoints, sizeof(float), float_cmp);
  float fMedianValue;
  fMedianValue = pfSmSortedValues[nNumPoints/2];
  
  // Remove the top 15% of pixel values

  float fTopValueLimit;

  fTopValueLimit = pfSmSortedValues[nNumPoints*85/100];
  
  // Set all pixel values above the fTopValueLimit to fTopValueLimit
  for (np = 0; np < nNumPoints; np++)  
    {
      if (fTopValueLimit < pfSmoothValues[np])
        pfSmoothValues[np] = fTopValueLimit;
    }

  // Find the location of the first pixel above the median value

  int nMedi = 0;
  while (pfSmoothValues[nMedi] < fMedianValue)
    {
      nMedi++;
    }
  //cout << "nMedi = " << nMedi << endl;

  // 3. Compute 1st and 2nd derivatives

  pfFirstDeriv[0]  = 0.0f;
  pfSecondDeriv[0] = 0.0f;
  for (np = 1; np < nNumPoints; np++)
    {
      if (0.0 < pfValues[np])
        {
          pfFirstDeriv[np]  = pfSmoothValues[np] - pfSmoothValues[np - 1];
          pfSecondDeriv[np] = pfFirstDeriv[np]   - pfFirstDeriv[np - 1];
        }
      else
        {
          pfFirstDeriv[np]  = -9999.0f;
          pfSecondDeriv[np] = 0.0f;
        }
    }

  //  4. Determine a threshold for the values.
  //    a. Find point maximum 1st derivative starting at the expected high value edge
  //       but only look up to the median pixel value
  //    b. Find the average of the 1st derivative at the nSmPoints high end (should be low)
  //    c. Find the point between (a) and then that is at least 10% heigher than value (b).

  int nMax;
  int nMin;
  nMax  = 0;
  nMin  = 0;
  for (np = 1; np < nMedi; np++)
    {
      if (pfFirstDeriv[nMax] < pfFirstDeriv[np])
        {
          nMax = np;
        }
      if (pfFirstDeriv[nMin] > pfFirstDeriv[np])
        {
          nMin = np;
        }
    }

  // If the minimum first derivative is hugely negative, this suggests a peak
  // in the array that we do not want to include in the mask.

  //if (nMin > nMax) cout << " WARNING peak for fpixin1: " << fPx1 << endl;

  np = nMedi;
  float a2fAvg[2];
  for (nx = 0; nx < 2; nx++)
    {
      fSum   = 0.0f;
      nCount = 0;
      while ( (0 <= np) && (nCount < nSmPoints) )
        {
          if (-9999.0f != pfFirstDeriv[np])
            {
              fSum += pfFirstDeriv[np];
              nCount++;
            }
          np--;
        }
      a2fAvg[nx] = fSum / (float) nCount;
    }

  // Of these two averaged derivatives, select the slope closest to 0;

  if (ABS(a2fAvg[0]) > ABS(a2fAvg[1]))
    {
      a2fAvg[0] = a2fAvg[1];
    }
  fSum = a2fAvg[0];

  // fSum contains the average 1st derivative closest to "flat"

  int nFlat;
  np = nMax;
  float fThreshold = -999.0;
  while (   (pfFirstDeriv[np] != -9999.0)
         && (pfFirstDeriv[np] > fSum)
         && (np < nNumPoints) )
    {
      fThreshold = pfSmoothValues[np];
      np++;
    }
  nFlat = np;
  nFlat = min(nMax + 2, nNumPoints);

  // 4.5
  // Report the max value to use as a threshold or cutoff
  // elsewhere:

  if (pfSmoothValues[nFlat] > fThreshold)
    fThreshold = pfSmoothValues[nFlat];

  delete [] pfFirstDeriv;
  delete [] pfSecondDeriv;
  delete [] pfValues;
  delete [] pfSmoothValues;
  delete [] pfSmSortedValues;
  return (fThreshold);
}

int CDTMainBeamMask::nRecurseShadow(int nCent0, 
                                   int nCent1, 
                                       int nEdgeSize,
                                   float fAvg,
                                   float fSD,
                                       float fThreshold,
                                   float fSlopePixVal, 
                                       Cimage *poImage,
                                   Cimage *poImageOut)
{
  // fSlopePixVal is not used just yet.
  
  // Look for pixels in poImage which have a value below fThreshold
  // in a region centered on pixel coordinate (nCent0, nCent1) 
  // of area nEdgeSize by nEdgeSize.  Any pixel values below fThreshold
  // will get set to a value of 0 in *poImageOut.
  // If there are no new pixels below the threshold, then exit, otherwise
  // if there are new pixels below the threshold, call this routine again
  // with an adjacent value of (nCent0, nCent1).  For the adjacent area 
  // go out in 4 different directions.

  int nLDebug = 1;
  int nStat;
  int  nCountNewPixelsInShadow;
  int  nCountNewPixelsInShadowDir;
  int nDir;
  int a4x2nEdge[4][2];
  int a2nCent[2];
  int nHalfEdge;
  a2nCent[0] = nCent0;
  a2nCent[1] = nCent1;
  int a5x2nDir[5][2] = { {0, 0}, {1, 0}, {-2, 0}, {0, 1}, {0, -2} };
      
  nHalfEdge = nEdgeSize / 2;

  m_nRecurseCount++;
  //cout << "nRecurseShadow " << nCent0 << ", " << nCent1 
  // << "   count: " << m_nRecurseCount << endl << flush;
  int nTooMany;
  nTooMany = poImageOut->nGetDimension() / 50;  // About 2% of pixels

  // Do not let recursion go forever
  if (   (nTooMany < m_nRecurseCount) 
         || (m_nShadowPixCount > (nTooMany * 5) ))
    {
      cout << "WARNING: Too many calls to CXdisplay::nRecurseShadow(...)\n"
           << "         Perhap the threshold is too high.\n"
           << endl;
      return (1);
    }
  // Try to calculate new average threshold here?
  C3Ddata *po3Ddata;
  po3Ddata = new C3Ddata(nEdgeSize, nEdgeSize, 1, 
                         a2nCent[0] - nHalfEdge,
                         a2nCent[1] - nHalfEdge,
                         0);
  nStat = po3Ddata->nFill2D(poImage, -1);
  float a3fBkg[3], fAvg2, fAvg3, fSD2, fMaxValue, a3fPeakCent[3], a3fTotCent[3];
  nStat = po3Ddata->nBackgroundAvgSD(0, &a3fBkg[0], &fAvg2, &fSD2, 
                                     &fMaxValue,
                                     &a3fPeakCent[0], &a3fTotCent[0]);
  delete po3Ddata;
  fAvg3 = fAvg2;
  po3Ddata = new C3Ddata(nEdgeSize, nEdgeSize, 1, 
                         a2nCent[0] - 20 * nHalfEdge,
                         a2nCent[1] - 20 * nHalfEdge,
                         0);
  nStat = po3Ddata->nFill2D(poImage, -1);
  if (0 == nStat)
    nStat = po3Ddata->nBackgroundAvgSD(0, &a3fBkg[0], &fAvg3, &fSD2, 
                                       &fMaxValue,
                                       &a3fPeakCent[0], &a3fTotCent[0]);
  delete po3Ddata;

  // Proceed only if the AVERAGE is less than the threshold

  if (fAvg2 > fThreshold)
    {
      //      cout << "Avg > Thresh: " << fAvg2 << " > " << fThreshold << endl;
      return (0);
    }

  // Try to find a suitable threshold value in this region


  float fThresholdUsed;  // Maybe threshold should be lowered as one goes to 
                         // to higher resolution?

 //                    pix0      pix1   slope  intercept
  float fPx0, fPx1;
  fPx0 = (float)nint((float)nCent0);
  fPx1 = (float)nint((float)nCent1);
  int nWidth = 12;
  float a4fThresh[4] = {0.0, 0.0, 0.0, 0.0};
  // y = slope * x   + y-intercept
  // 2 horizontal lines
  //                           x    y   slope intercept
  a4fThresh[0] = fLineShadow(fPx0, fPx1, 0.0, fPx1, 
                      (float)nWidth, 151, 7, 1, poImage);
  a4fThresh[1] = fLineShadow(fPx0, fPx1, 0.0, fPx1, 
                      (float)nWidth, 151, 7, -1, poImage);

  // 2 more vertical lines

  a4fThresh[2] = fLineShadow(fPx0, fPx1, 1000.0, fPx0,
                      (float)nWidth, 151, 7, 1, poImage);
  a4fThresh[3] = fLineShadow(fPx0, fPx1, 1000.0, fPx0, 
                      (float)nWidth, 151, 7, -1, poImage);
  int np;
  int nCount = 0;
  float fSum = 0.0;

  // Average the top two thresholds

  qsort(&a4fThresh[0], 4, sizeof(float), float_cmp);

  for (np = 2; np < 4; np++)
    {
      if (0.0 < a4fThresh[np])
        {
          fSum += a4fThresh[np];
          nCount++;
        }
    }
  if (0 < nCount)
    fSum = fSum / (float)nCount;

  fThresholdUsed = 0.6 * fSum;

  if (0.0 >= fThresholdUsed)
    {
      // fLineShadow did not do a good job;
      fThresholdUsed = fThreshold;
    }
  //cout << "nCent0,1: " << nCent0 << ", " <<  nCent1 << endl;
  //cout << "tholdused, hold: " << fThresholdUsed << "   " << fThreshold << endl;
  nCountNewPixelsInShadow = 0;

  for (nDir = 0; nDir < 5; nDir++)
    {
      nCountNewPixelsInShadowDir = 0;

      // Go look in 5 directions (really 4, but include the center)

      // Shift to a new center to look at

      a2nCent[0] += (18 * nEdgeSize * a5x2nDir[nDir][0] / 20);
      a2nCent[1] += (18 * nEdgeSize * a5x2nDir[nDir][1] / 20);

      int nx, ny;
      int nxx, nyy;
      bool bSetToZero;
      float fPix1;
      float fPix2;
      float fPixPrev = 0.0f;
      nyy = a2nCent[1] - nHalfEdge;
      for (ny = 0; ny < nEdgeSize; ny++)
        {
          nyy++;
          nxx = a2nCent[0] - nHalfEdge;
          for (nx = 0; nx < nEdgeSize; nx++)
            {
              nxx++;
              fPix1 = (poImage->*poImage->prfGetPixel)(nxx, nyy);
              fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx, nyy);       
              bSetToZero = FALSE;
              if ( (fPix1 < fThresholdUsed) && (m_fEraseValue < fPix1) )
                {
                  // if pixel in poImage is less < threshold,
                  // then check pixel in poImageOut is not "set"
                  // if not set, then set it and set nNewPixels to true
                    
                  if (fPix2 == fPix1)
                    {
                      // Not set, so set to m_fEraseValue

                      bSetToZero = TRUE;
                    }
                }
              else if (m_fEraseValue < fPix2)
                {
                  // If 5 of 8 nearest neighbors are less than or equal to m_fEraseValue, then set this
                  // this one to 0 as well

                  fPix1 = fPix2;
                  int nCount = 0;
                  float fAvg = 0.0;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx+1, nyy-1);
                  if (m_fEraseValue >= fPix2) 
                    nCount++;
                  else
                    fAvg += fPix2;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx+1, nyy);
                  if (m_fEraseValue >= fPix2)
                    nCount++;
                  else
                    fAvg += fPix2;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx+1, nyy+1);
                  if (m_fEraseValue >= fPix2)
                    nCount++;
                  else
                    fAvg += fPix2;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx, nyy-1);
                  if (m_fEraseValue >= fPix2)
                    nCount++;
                  else
                    fAvg += fPix2;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx, nyy+1);
                  if (m_fEraseValue >= fPix2)
                    nCount++;
                  else
                    fAvg += fPix2;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx-1, nyy-1);
                  if (m_fEraseValue >= fPix2)
                    nCount++;
                  else
                    fAvg += fPix2;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx-1, nyy);
                  if (m_fEraseValue >= fPix2)
                    nCount++;
                  else
                    fAvg += fPix2;
                  fPix2 = (poImageOut->*poImageOut->prfGetPixel)(nxx-1, nyy+1);
                  if (m_fEraseValue >= fPix2)
                    nCount++;
                  else
                    fAvg += fPix2;

                  bSetToZero = (5 <= nCount);
                  if (bSetToZero) 
                    {
                      //cout << "5/8 set\n" << flush;
                    }
                  else
                    {
                      // Compute average of nearest non-zero neighbors
                      nCount = 8 - nCount;
                      if ( (0 < nCount) && (5 >= nCount) )
                        {
                          // At least 3 neighboring pixels should be in the
                          // the shadow

                          fAvg = fAvg / float(nCount);
                          if (fPix1 < (0.85 * fAvg))
                            {
                              // If pixel value is less than 60% of average
                              // then it is in the shadow

                              bSetToZero = TRUE;
                              //cout << "0.85 avg set\n" << flush;
                            }
                        }
                    }
                }
              if (bSetToZero)
                {
                  fPix2 = m_fEraseValue;
                  (void) (poImageOut->*poImageOut->prnSetPixel)(nxx, nyy,
                                                                fPix2);
                  if (1 == nLDebug)
                    (void) (poImage->*poImage->prnSetPixel)(nxx, nyy,
                                                                fPix2);
                  nCountNewPixelsInShadowDir++;
                }
            }
        }
      int nTest = 0;
      // Only call nRecurseShadow if a significant number of pixels in the
      // just tested region were also in the shadow.
      nTest = nEdgeSize * nEdgeSize * 10 / 25;
      if (nTest < nCountNewPixelsInShadowDir)
        {
          m_nShadowPixCount = m_nShadowPixCount + nCountNewPixelsInShadowDir;
          //cout << "ShadowPixCount: " << m_nShadowPixCount << endl << flush;
          
          // Do not call nRecurseShadow if you are near the edge of the 
          // detector
          int nTooClose;
          nTooClose = max(5, nEdgeSize);
          if (    (a2nCent[0] < nTooClose)
               || (a2nCent[1] < nTooClose)
               || ( poImage->nGetDimension(0) - a2nCent[0] < nTooClose)
               || ( poImage->nGetDimension(1) - a2nCent[1] < nTooClose) )
            {
              //cout << "Too close to edge\n" << flush;
            }
          else
            {
              
              nStat = nRecurseShadow(a2nCent[0], a2nCent[1],
                                     nEdgeSize, fAvg, fSD, 
                                     fThreshold, fSlopePixVal,
                                     poImage, poImageOut);
            }
          if (1 == nLDebug)
            if (0 == (m_nRecurseCount % 500))
              {
                //vRefresh();
              }
          if (0 != nStat) return (nStat);
        }
      else
        {
          //cout << "No new pixels in shadow for "
          //<< nCent0 << ", " << nCent1 << endl << flush;
        }
      nCountNewPixelsInShadow += nCountNewPixelsInShadowDir;
    } // end direction test

  return (0);
}

int CDTMainBeamMask::nEditPixels(Cstring sLine, float fDefaultPixVal, Cimage *poImage)
{
  // It would probably be best if this routine was in the Cimage:: class.

  // This routine was lifted from the Ccalibrate:: routine of a similar name
  // and modified for use here.

  // Read the file rsFilename for a list of bad pixel commands
  // Valid commands are ROW, COLUMN, LINE, RECT, PIXEL, CIRCLE, INSIDECIRCLE,
  // OUTSIDECIRCLE


  int     nStat = 0;
  Cstring asTokens[10];
  int     nNumToks;
  int     i, j, ii, jj;
  float   fPixVal;

  // asTokens will be an array of strings

  nNumToks = split(sLine, asTokens, 10, " ,\t");
  if (0 == nNumToks)
    asTokens[0] = "";
  asTokens[0].downcase();
  if (   (0 == nNumToks)
      || (0 == asTokens[0].index('!')) )
    {
      // Skip the commented line or empty line
    }
  else if (   (asTokens[0] == "col")
           || (asTokens[0] == "column") )
    {
      // column i, start, end, flag
      if (3 < nNumToks)
        {
          if (4 < nNumToks)
            fPixVal = atof(asTokens[4].string());
          else
            fPixVal = fDefaultPixVal;
          i = atoi(asTokens[1].string());
          for (j =  atoi(asTokens[2].string());
               j <= atoi(asTokens[3].string()); j++)
            {
              (void) (poImage->*poImage->prnSetPixel)(i, j, fPixVal);
            }
        }
    }
  else if (   (asTokens[0] == "line")
           || (asTokens[0] == "row") )
    {
      // row j, start, end | flag
      if (3 < nNumToks)
        {
          if (4 < nNumToks)
            fPixVal = atof(asTokens[4].string());
          else
            fPixVal = fDefaultPixVal;
          j = atoi(asTokens[1].string());
          for (i = atoi(asTokens[2].string()); 
               i <= atoi(asTokens[3].string()); i++)
            {
              (void) (poImage->*poImage->prnSetPixel)(i, j, fPixVal);
            }
        }
    }
  else if (   (asTokens[0] == "rectangle")
              || (asTokens[0] == "rect") )
    {
      // Rect orig0, orig1, ext0, ext1 | flag
      
      if (4 < nNumToks)
        {
          if (5 < nNumToks)
            fPixVal = atof(asTokens[5].string());
          else
            fPixVal = fDefaultPixVal;
          int ii, jj;
          jj = atoi(asTokens[2].string());
          for (j = 0; j < atoi(asTokens[4].string()); j++, jj++)
            {
              ii = atoi(asTokens[1].string());
              for (i = 0; i < atoi(asTokens[3].string()); i++, ii++)
                {
                  (void) (poImage->*poImage->prnSetPixel)(ii, jj, fPixVal);
                }
            }
        }
      
    }
  else if (   (asTokens[0] == "insidecircle")
           || (asTokens[0] == "circle") )
    {
      // Circle center0, center1, radius0 | flag
      
      if (3 < nNumToks)
        {
          if ( 4 < nNumToks)
            fPixVal = atof(asTokens[4].string());
          else
            fPixVal = fDefaultPixVal;
          int ii, jj;
          int icen, jcen, iradsq;
          icen   = atoi(asTokens[1].string());
          jcen   = atoi(asTokens[2].string());
          iradsq = atoi(asTokens[3].string());
          cout << "center, radius: " << icen << ", " << jcen << ", " 
               << iradsq << endl;
          
          iradsq = iradsq * iradsq;
          int nDim0, nDim1;
          (void)poImage->nGetDimensions(&nDim0, &nDim1);
          for (j = 0; j < nDim1; j++)
            {
              jj = (j-jcen) * (j-jcen);
              for (i = 0; i < nDim0; i++)
                {
                  ii = (i-icen) * (i-icen);
                  if (iradsq > (ii + jj) )
                    {
                      (void) (poImage->*poImage->prnSetPixel)(i, j, fPixVal);
                    }
                }
            }
        }
    }
  else if (asTokens[0] == "outsidecircle")
    {
      // Circle center0, center1, radius0 | flag
      
      if (3 < nNumToks)
        {
          if ( 4 < nNumToks)
            fPixVal = atof(asTokens[4].string());
          else
            fPixVal = fDefaultPixVal;
          int ii, jj;
          int icen, jcen, iradsq;
          icen   = atoi(asTokens[1].string());
          jcen   = atoi(asTokens[2].string());
          iradsq = atoi(asTokens[3].string());
          cout << "center, radius: " << icen << ", " << jcen << ", " 
               << iradsq << endl;
          
          iradsq = iradsq * iradsq;
          
          int nDim0, nDim1;
          (void)poImage->nGetDimensions(&nDim0, &nDim1);
          for (j = 0; j < nDim1; j++)
            {
              jj = (j-jcen) * (j-jcen);
              for (i = 0; i < nDim0; i++)
                {
                  ii = (i-icen) * (i-icen);
                  if (iradsq < (ii + jj) )
                    {
                      (void) (poImage->*poImage->prnSetPixel)(i, j, fPixVal);
                    }
                }
            }
        }
    }
  else if (asTokens[0] == "pixel")
    {
      // pixel i j | flag
      if (2 < nNumToks)
        {
          if (3 < nNumToks)
            fPixVal = atof(asTokens[3].string());
          else
            fPixVal = fDefaultPixVal;
          i = atoi(asTokens[1].string());
          j = atoi(asTokens[2].string());
          (void) (poImage->*poImage->prnSetPixel)(i, j, fPixVal);
        }
    }
  else
    {
      cout << "WARNING: unknown bad pixel keyword: " << asTokens[0] << endl;
      nStat = 1;
    }
  return (nStat);
}
void CDTMainBeamMask::vError(const int nErrorNum, const Cstring& sMessage)
{
}
