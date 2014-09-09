//#
//#Copyright 1998 Molecular Structure Corporation
//#               9009 New Trails Drive
//#               The Woodlands, TX, USA  77381
//#
//#The contents are unpublished proprietary source
//#code of Molecular Structure Corporation
//#
//#All rights reserved
//#
//#Ctransform.cpp      Initial author:            Jan 1999
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

// Description:
//
//
// ToDo:
//
//

/****************************************************************************
 *                              Include Files                               *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "DevErrCode.h"
#include "Ctransform.h"

#ifdef ANL_TIMER
#include "anlTimer.h"
#endif


#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

/****************************************************************************
 *                         Static Member Variables                          *
 ****************************************************************************/

/****************************************************************************
 *                             Global variables                             *
 ****************************************************************************/

/****************************************************************************
 *                            External variables                            *
 ****************************************************************************/

/****************************************************************************
 *                           Function prototypes                            *
 ****************************************************************************/

/****************************************************************************
 *                           Non-class functions                            *
 ****************************************************************************/

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                         Public Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 * This is the constructor for the class.                                   *
 ****************************************************************************/
Ctransform::Ctransform(void)
{
  m_a2nTOutDim[0]    = 0;
  m_a2nTOutDim[1]    = 0;
  m_poImgTransform   = NULL;
  m_sTransformName   = sTransSymbol("$(TRANSFORM)");

  m_lNumTime = 0;
  m_dSumTime = m_dSumTime2 = 0;
  
  m_poImgBadInterp   = NULL;
  m_nNumBadInterpOut = 0;
}

/****************************************************************************
 * This is the destructor for the class.                                    *
 ****************************************************************************/
Ctransform::~Ctransform()
{
  if (NULL != m_poImgTransform)
    {
      delete m_poImgTransform;
      m_poImgTransform = NULL;
    }
  if (NULL != m_poImgBadInterp)
    {
      delete m_poImgBadInterp;
      m_poImgBadInterp = NULL;
      m_nNumBadInterpOut = 0;
    }
}

/****************************************************************************
 ****************************************************************************/
int
Ctransform::nTransform(Cimage *poImageIn, Cimage *poImageOut,
                                 Cimage *poImageDark)
{
  // Transform the input image into the output image by using
  // the m_poImgTransform information which applies the non-uniformity
  // and spatial distortion corrections simultaneously.
  //
  // The transform image m_poImgTransform has the following properties:
  // There are pairs of values:  (input_pixel_offset, 65536 * scale_factor)
  // sorted in the order of which output_pixel they contribute to.
  // A scale_factor=0 means the output_pixel is bad and/or has no input_pixels
  // that contribute to it.  If an output_pixel has contributions from more
  // than one input_pixel, then the contributions appear in the image from
  // lowest scale_factor to highest scale_factor.  Thus, if a scale_factor is
  // lower than the previous scale_factor, that means a new output_pixel is
  // being worked on.  A scale_factor of 1 is an indication to go to the
  // next output_pixel and the associated input_pixel does not contribute
  // to the output pixel.  In this way, the m_poImgTransform does not hold
  // any explicit information about which output_pixel the input_pixel(s)
  // contribute to.

  // If NULL != poImageDark, then *poImageIn is modified during the transform.

  // CAVEAT:  This ONLY works for (unsigned short int) images!

  int nStat = 0;
  if (NULL == m_poImgTransform)
    {
      poImageIn->nGetDimensions(&m_a2nTOutDim[0], &m_a2nTOutDim[1]);
      cout << "Reading transform image... " << m_sTransformName << endl;
      m_poImgTransform = new Cimage(m_sTransformName);
      if (m_poImgTransform->bIsAvailable())
        {
          (void) m_poImgTransform->m_oHeader.nGetValue("OUTPUT_SIZE", 2,
                                                       m_a2nTOutDim);
        }
    }
  if (!m_poImgTransform->bIsAvailable())
    {
        cout << "ERROR, transform image not available!" << endl;
      nStat = 3;
    }
  else
    {
      // Do the transformation

      int nDim0, nDim1;
      poImageIn->nGetDimensions(&nDim0, &nDim1);

      float fCalibScale = 1.0;
      int nScale;
      int nScalePrev;
      int nSaturated;
      int nPedestal      = 20;
      int nInPedestal    = 0;
      int nBadValue      = 0;
      int nInOffsetPixel;

      int nNumTruncUp;
      int nNumTruncDown;
      int nNumBad;
      int nNumAvail;
      int nOutput;
      int nInput;
      int nDark;
      int nOutDim0, nOutDim1;
      short int iInOffsetDelta;
      bool bSaturated = FALSE;
      double dTime;

      float fSumPix, fPix;
      float fBadValue;
      int   nNumPix, nOffset;
      unsigned short int uiPix;
      int i;
      int nJ0, nJ1, nJV, nJF;

#ifdef ANL_TIMER
      anl_reset_timer(0, "Ccal::nTransform, Tran");
      anl_start_timer(0);
#endif

      // Get saturated value and number of input values for the transform

      nSaturated = (int) poImageIn->fGetSatValue();
      nSaturated = min(65535, nSaturated);

      if (0 != m_poImgTransform->m_oHeader.nGetValue("CALIB_PEDESTAL", &nPedestal))
        (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_PEDESTAL", &nPedestal);
      (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_BADVALUE", &nBadValue);
      (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_SCALE", &fCalibScale);
      int nNumBadOut = 0;
      (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_BI_OUTPUT", &nNumBadOut);

      poImageOut->nGetDimensions(&nOutDim0, &nOutDim1);
      cout << "Expected output size from Transform header: "
           << m_a2nTOutDim[0] << " by " << m_a2nTOutDim[1]
           << "\nActual output size from output header: "
           << nOutDim0 << " by " << nOutDim1 << endl;
      if (   (nOutDim0 != m_a2nTOutDim[0])
          || (nOutDim1 != m_a2nTOutDim[1]) )
          cout << "ERROR in nTransform, mismatch of output image dimensions!\n"
        << flush;

      cout << "Transform scale factor: " << fCalibScale << endl;

      if (NULL == poImageDark)
        {
          if (0 != poImageIn->m_oHeader.nGetValue("DARK_PEDESTAL", &nInPedestal))
            nInPedestal = 0;
        }
      else
        {
          // Subtract dark image, keeping saturated pixels saturated
          // TODO: Special problem with pedestal!  Should it be added
          //       back in during dark subtraction?  Then removed later?

#ifdef ANL_TIMER
          anl_reset_timer(1, "Ccal::nTransform, Dark");
          anl_start_timer(1);
#else
          dTime = dGetTimeOfDay();
#endif
          nNumAvail = poImageDark->nGetDimension();
          poImageIn->nSetNextPixel(0, 0);
          poImageDark->nSetNextPixel(0, 0);
          nNumTruncUp   = 0;
          nNumTruncDown = 0;
          nNumBad         = 0;
          while (0 < nNumAvail)
            {
              nOutput =  (int)poImageIn->uiGetNextPixelNoInc();
              nDark   =  (int)poImageDark->uiGetNextPixel();
              if (nSaturated > nOutput)
                {
                  nOutput = nOutput - nDark + nInPedestal;
                  if (0 > nOutput)
                    {
                      nOutput = 0;
                      nNumTruncUp++;
                    }
                  else if (65535 < nOutput)
                    {
                      nOutput = 65535;
                      nNumTruncDown++;
                    }
                }
              poImageIn->vSetNextPixel( (unsigned short int) nOutput);
              nNumAvail--;
            }
#ifdef ANL_TIMER
          anl_stop_timer(1);
#else
          dTime = dGetTimeOfDay()-dTime;
          cout << "\nTime for CTransform dark subtraction: "
               << dTime << " seconds" << endl;
#endif
          cout << "For dark subtraction: "
               << "\nNumber truncated down to 65535: " << nNumTruncDown
               << "\nNumber truncated up   to     0: " << nNumTruncUp
               << endl;
        }

      int nNumWritten = 0;
      nNumAvail       = 0;
      nNumTruncUp     = 0;
      nNumTruncDown   = 0;
      nNumBad         = 0;

      nStat = m_poImgTransform->m_oHeader.nGetValue("CALIB_USED", &nNumAvail);

      // If interpolation of bad output pixels is to be done, 
      //  then prepare for that

      if (0 < nNumBadOut)
        {
          if (0 == m_nNumBadInterpOut)
            {
              m_nNumBadInterpOut = nNumBadOut;
              m_poImgBadInterp   = new Cimage(4, m_nNumBadInterpOut, eImage_I4);
              (void) m_poImgBadInterp->nSetNextPixel(0, 0);
            }
          else 
            {
              nNumBadOut = 0;
            }
        }

      // Prime the algorithm

      fCalibScale = fCalibScale / 65536.0f;

#ifndef ANL_TIMER
      dTime = dGetTimeOfDay();
#endif

      poImageOut->nSetNextPixel(0, 0);
      m_poImgTransform->nSetNextPixel(0, 0);
      nOutput         = 0;
      nInOffsetPixel  = 0;
      nScalePrev      = -1;

      // Loop through all available input_pixel contributions

      while (0 < nNumAvail)
        {
          iInOffsetDelta  = m_poImgTransform->iGetNextPixel();
          nScale          = (int)m_poImgTransform->uiGetNextPixel();
          nInOffsetPixel += iInOffsetDelta;
          nNumAvail--;

          // A value of -32768 for iInOffsetDelta (if nScale not 0 nor 1)
          // indicates next pixel pair has the actual nInOffsetPixel

//          if (   (0 == iInOffsetDelta)
//              && (1 < nScale)
//              && (0 < nNumAvail) )
          if (-32768 == iInOffsetDelta) // fewest
            if (1 < nScale)        // next fewest
              if (0 < nNumAvail)   // really never happens
                {
                  // This means the next "pixels" contain the actual nInOffsetPixel
                  // value

                  nInOffsetPixel  = 65536 * m_poImgTransform->iGetNextPixel();
                  nInOffsetPixel += (int)m_poImgTransform->uiGetNextPixel();
                  nNumAvail--;
                }


          if (0 == nScalePrev)
            {
              // Output pixel is flagged as bad

              nNumBad++;
              nOutput = nBadValue;
              poImageOut->vSetNextPixel( (unsigned short int) nOutput);
              if (0 != nNumBadOut)
                {
                  // Do this only the 1st image through here 
                  // AND if interp required

                  // Store list of bad pixels in the output image

                  // debug line:
//                  if (nNumBadOut == m_nNumBadInterpOut)
//                    cout << "FIRST bad interpout!\n" << flush;

                  nNumBadOut--;
                  m_poImgBadInterp->vSetNextPixel(
                          (INT4) (nNumWritten % nOutDim0)); // px0
                  m_poImgBadInterp->vSetNextPixel(
                          (INT4) (nNumWritten / nOutDim0)); // px1
                  m_poImgBadInterp->vSetNextPixel(
                          (INT4)0); // value
                  m_poImgBadInterp->vSetNextPixel(
                          (INT4)0); // flag?
                }
              nNumWritten++;
              nOutput    = 0;
              bSaturated = FALSE;
            }
          else if (nScale < nScalePrev)
            {
              // Write out previous output pixel, but round-up,
              // scale properly and add in the pedestal.

              //          nOutput = ((nOutput + 32768) >> 16) + nPedestal;
              //          nOutput = ((nOutput+32768) / 65536) + nPedestal;

              nOutput = int(((float) (unsigned int)nOutput * fCalibScale) + 0.5) + nPedestal;

              if ( bSaturated || (nOutput >= nSaturated) )
                {
                  nNumTruncDown++;
                  nOutput = nSaturated;
                }
              else if (0 <= nOutput)
                {
                  // No nothing, this should be the vast majority of the pixels
                }
              else if (-16384 > nOutput)
                {
                  // Occurs on overflows into the sign bit of nOutput

                  nNumTruncDown++;
                  nOutput = nSaturated;
                }
              else if (0 > nOutput)
                {
                  nNumTruncUp++;
                  nOutput = 0;
                }

              poImageOut->vSetNextPixel( (unsigned short int) nOutput);
              nNumWritten++;
              nOutput    = 0;
              bSaturated = FALSE;
            }
          if (1 < nScale)
            {
              // This pixel contributes to currently being summed output pixel

              nInput = (int) poImageIn->uiGetPixel(nInOffsetPixel);
              if (nInput >= nSaturated)
                {
                  bSaturated = TRUE;
                }
              nInput  -= nInPedestal;
	      if (0 > nInput) nInput = 0;  // nInput cannot be negative!
              nOutput += (nScale * nInput);
            }
          nScalePrev = nScale;
        } // end while

      // If requested, apply any interpolation of bad output pixels
      // after all output pixels are available


      if (NULL != m_poImgBadInterp)
        {
//          cout << "::nTransform, doing interpolation!\n" << flush;
//          m_poImgBadInterp->nWrite("badinterp.img");

          nInterpolationFunction(nOutDim0, nOutDim1, nBadValue,
				 m_nNumBadInterpOut, m_poImgBadInterp, poImageOut);

          // Now that all interpolated values are calculated,
          // copy them into the output image
          
          m_poImgBadInterp->nSetNextPixel(0, 0);          
          for (i = 0; i < m_nNumBadInterpOut; i++)
            {
              nJ0     = (int) m_poImgBadInterp->lGetNextPixel(); // px0
              nJ1     = (int) m_poImgBadInterp->lGetNextPixel(); // px1
              uiPix   = (unsigned short int) m_poImgBadInterp->lGetNextPixel(); // value
              nJF     = (int) m_poImgBadInterp->lGetNextPixel(); // flag
              if (1 < uiPix)
                {
                  nOffset = nJ0 + (nOutDim0 * nJ1);
                  (void)poImageOut->vSetPixel(nOffset, uiPix);
/*
                  (void)poImageOut->nSetPixel(nJ0, nJ1, uiPix);
                  cout << "OKpix nOffset, px0, px1: " << nOffset << ", "
                       << nOffset % nOutDim0 << ", " << nOffset / nOutDim0
                       << "   :   " << nJ0 << ", " << nJ1 << ", " << uiPix
                       << endl << flush;
                }
              else
                {
                  cout << "baduipix nOffset, px0, px1: " << nOffset << ", "
                       << nOffset % nOutDim0 << ", " << nOffset / nOutDim0
                       << "   :   " << uiPix
                       << endl << flush;
*/
                }
            };
        } // end do interpolation if


#ifndef ANL_TIMER
      dTime = dGetTimeOfDay()-dTime;
      m_lNumTime++;
      m_dSumTime += dTime;
      m_dSumTime2 += dTime*dTime;
      cout << "\nTime for CTransform spatial and nonunf correction: "
           << dTime << " seconds" << endl;
#endif

      if (nNumWritten < poImageOut->nGetDimension())
        {
          // One should NOT get here!

          cout << "\nWARNING in nTransform, output is lite!\n";
        }

//      std::cout << "nNumWritten: " << nNumWritten << endl;
      cout << "For spatial and nonunf correction: "
           << "\nNumber truncated down to 65535: " << nNumTruncDown
           << "\nNumber truncated up   to     0: " << nNumTruncUp
           << "\nNumber flagged as bad:          " << nNumBad
           << "  (" << nNumBad / (nOutDim0 * nOutDim1) * 100 << "%)"
           << endl;
#ifdef ANL_TIMER
      anl_stop_timer(0);
      if (NULL != poImageDark) anl_print_timer(1);
      anl_print_timer(0);
#endif

/*
 * Add some nonunf and spatial distortion info to the output header.
 */

   // Copy the input header to the output header to make sure that things
   // get kept correct.

      poImageOut->m_oHeader = poImageIn->m_oHeader;

   // Add the nonunf and spatial distortion info to the header
   // Ignore any errors that might occur in the addition, since it
   // doesn't affect the transform.

      nAddSpatialInfo(poImageOut->m_oHeader);
      nAddNonunfInfo(poImageOut->m_oHeader);

    }

  return (nStat);
}

/****************************************************************************
 * This routine adds the nonunf info to the header.                         *
 ****************************************************************************/
int
Ctransform::nAddNonunfInfo(Cimage_header &oHeader)
{
   int nStatus;
   Cstring sKeyword,sDetName;


   sDetName = sGetDetectorPrefix(oHeader);
   sKeyword  = sDetName;
   sKeyword += D_K_NonunfInfo;
   nStatus = oHeader.nReplaceValue(sKeyword,"None");
   if (0 == nStatus)
     {
       sKeyword = sDetName;
       sKeyword += D_K_NonunfType;
       nStatus = oHeader.nReplaceValue(sKeyword,"None");
     }

   return nStatus;
}

/****************************************************************************
 * This routine adds the spatial distortion info to the header.             *
 ****************************************************************************/
int
Ctransform::nAddSpatialInfo(Cimage_header &oHeader)
{
   int nStatus;
   double d,adInfo[4],mmDimensions[2],PixelDimensions[2],UnbinnedPixelDimensions[2];
   double DirectBeamPosition[2];
   Cstring sKeyword,sDetName;

   sDetName = sGetDetectorPrefix(oHeader);

   // Get the dimensions of the detector in pixels.
   // These may be binned pixels.
   sKeyword = sDetName;
   sKeyword += D_K_DetectorDimensions;
   oHeader.nGetValue(sKeyword,2,PixelDimensions);
   // Get # of unbinned pixels and size of detector
   sKeyword = sDetName;
   sKeyword += D_K_UnbinnedDimensions;
   if (0 != oHeader.nGetValue(sKeyword,2,UnbinnedPixelDimensions))
   {
     UnbinnedPixelDimensions[0] = PixelDimensions[0];
     UnbinnedPixelDimensions[1] = PixelDimensions[1];
   }
   sKeyword  = sDetName;
   sKeyword += D_K_DetectorSize;
   oHeader.nGetValue(sKeyword,2,mmDimensions);

   // See if we know where the unbinned direct beam is at.
   // If not, then just use the midpoint of the unbinned detector
   // when building the spatial distortion info.

   sKeyword = sDetName;
   sKeyword += D_K_UnbinnedBeamPosition;
   if (0 != oHeader.nGetValue(sKeyword,2,DirectBeamPosition))
   {
     DirectBeamPosition[0] = UnbinnedPixelDimensions[0]/2;
     DirectBeamPosition[1] = UnbinnedPixelDimensions[1]/2;
   }

   d         = UnbinnedPixelDimensions[0] / PixelDimensions[0];
   adInfo[0] = DirectBeamPosition[0] / d;
   d         = UnbinnedPixelDimensions[1] / PixelDimensions[1];
   adInfo[1] = DirectBeamPosition[1] / d;
   adInfo[2] = mmDimensions[0] / PixelDimensions[0];
   adInfo[3] = mmDimensions[1] / PixelDimensions[1];

   // Add what we think is the spatial distortion info

   sKeyword  = sDetName;
   sKeyword += D_K_SpatialDistortionInfo;
   nStatus   = oHeader.nReplaceValue(sKeyword,4,adInfo);

   // Spatial distortion beam position

   if (0 == nStatus)
     {
       sKeyword  = sDetName;
       sKeyword += D_K_SpatialBeamPosition;
       nStatus   = oHeader.nReplaceValue(sKeyword,2,adInfo);
     }

   // Spatial distortion type

   if (0 == nStatus)
     {
       sKeyword  = sDetName;
       sKeyword += D_K_SpatialDistortionType;
       nStatus   = oHeader.nReplaceValue(sKeyword,D_K_SpatialTypeSimple);
     }

   // Spatial distortion vectors

   if (0 == nStatus)
     {
       sKeyword  = sDetName;
       sKeyword += D_K_SpatialDistortionVectors;
       nStatus   = oHeader.nReplaceValue(sKeyword,"0 -1 1 0");
     }

   return nStatus;
}

/****************************************************************************
 * Returns timing information concerning the amount of time that it takes   *
 * to actually transform an image.                                          *
 ****************************************************************************/
int
Ctransform::nGetTimingInfo(double *pdValue,
                           double *pdSigma)
{
  if( (NULL == pdValue) || (NULL == pdSigma) )
    {
      return DEV_FAILED;
    }
  else if (0 == m_lNumTime)
    {
      *pdValue = *pdSigma = 0.0;
      return DEV_FAILED;
    }

   if (1 == m_lNumTime)
     {
       *pdValue = m_dSumTime;
       *pdSigma = 0.0;
     }
   else
     {
       *pdValue = m_dSumTime / (double)m_lNumTime;
       *pdSigma = sqrt((m_dSumTime2*m_lNumTime-m_dSumTime*m_dSumTime)
                       / (double)m_lNumTime / (double)(m_lNumTime-1));
     }

  return DEV_SUCCESS;
}

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                       Protected Member Functions                       **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 ****************************************************************************
 **                                                                        **
 **                        Private Member Functions                        **
 **                                                                        **
 ****************************************************************************
 ****************************************************************************/

/****************************************************************************
 * This routine returns the detector prefix listed in the header.  It is    *
 * assumed that the last detector name listed is the prefix to use.  If no  *
 * detector prefix can be located, then the prefix "CCD_" is added to the   *
 * header and is used.                                                      *
 ****************************************************************************/
Cstring
Ctransform::sGetDetectorPrefix(Cimage_header &oHeader)
{
   int nDetectors,nStatus;
   Cstring sKeyword,sDetName,*psDetNames;

   nDetectors = 0;
   sKeyword   = D_K_DetectorNumber;
   nStatus    = oHeader.nGetValue(sKeyword,&nDetectors);
   if (0 == nStatus)
     {
       if (0 < nDetectors)
         psDetNames = new Cstring [nDetectors];
       else
         nStatus = -1;
     }
   if (0 == nStatus)
     {
       sKeyword = D_K_DetectorNames;
       nStatus  = oHeader.nGetValue(sKeyword,nDetectors,psDetNames);
     }
   if (0 != nStatus)
     {
       sDetName = "CCD_";
       //cerr << "\nWARNING Could not get values for\n"
            //<< sKeyword
            //<< "\nin detector header.\nAdding \"" << sDetName
            //<< "\" as detector name prefix."
            //<< endl << endl;
       nStatus = oHeader.nReplaceValue(D_K_DetectorNumber,1);
       nStatus = oHeader.nReplaceValue(D_K_DetectorNames,sDetName);
     }
   else
     sDetName = psDetNames[nDetectors-1];
   if (0 < nDetectors)
     delete [] psDetNames;

   return sDetName;
}

