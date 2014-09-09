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
// Cnonunf.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Cnonunf which implements
//    the non-uniformity correction  encapsulation of d*TREK.
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

#include "Dtrek.h"
#include "Cnonunf.h"         // Class definition and prototypes
#include "dtrekdefs.h"
#include "Cscan.h"
#include "ImagePool.h"
#include "CImageWait.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

//+Definitions, constants, and initialization of static member variables

int Cnonunf::ms_nBadPixFlag = 2048;
int Cnonunf::ms_nSatPixFlag = 4096;
Cstring Cnonunf::ms_sNonunfStateUnknown       = D_K_NonunfStateUnknown;
Cstring Cnonunf::ms_sNonunfStateNone          = D_K_NonunfStateNone;
Cstring Cnonunf::ms_sNonunfStateSimpleScale   = D_K_NonunfStateSimpleScale;
Cstring Cnonunf::ms_sNonunfStateSimpleMask    = D_K_NonunfStateSimpleMask;
Cstring Cnonunf::ms_sNonunfStateDarkNonunf    = D_K_NonunfStateDarkNonunf;
Cstring Cnonunf::ms_sNonunfStateDarkDCoffsetNonunf
                                              = D_K_NonunfStateDarkDCoffsetNonunf;
Cstring Cnonunf::ms_sNonunfType               = D_K_NonunfType;
Cstring Cnonunf::ms_sNonunfDenominator        = D_K_NonunfDenominator;
Cstring Cnonunf::ms_sNonunfNumerator          = D_K_NonunfNumerator;
Cstring Cnonunf::ms_sNonunfFlag1              = D_K_NonunfFlag1;
Cstring Cnonunf::ms_sNonunfFlag2              = D_K_NonunfFlag2;
Cstring Cnonunf::ms_sNonunfFlag3              = D_K_NonunfFlag3;
Cstring Cnonunf::ms_sNonunfFlag4              = D_K_NonunfFlag4;
Cstring Cnonunf::ms_sNonunfInfo               = D_K_NonunfInfo;
Cstring Cnonunf::ms_sNonunfKey                = D_K_NonunfKey;

const bool c_bDeleteImageIfNotReferenced = false;

//+Code begin

//+Public functions


// Constructors, destructors and assignments

Cnonunf::Cnonunf()
{
 (void) nInitValues();
}

Cnonunf::Cnonunf(const Cstring& sNonunfName, const Cstring& sDarkName)
{
  (void) nInitValues(sNonunfName, sDarkName);
}

int Cnonunf::nInitValues(const Cstring& sNonunfName, const Cstring& sDarkName)
{
  (void)nInitValues();
  cout << "...reading nonunf and dark image files..." << endl;
  
  m_poNonunfImage = m_pImagePool->poGetImage(sNonunfName, !c_bDeleteImageIfNotReferenced);

  m_poDarkImage = m_pImagePool->poGetImage(sDarkName, !c_bDeleteImageIfNotReferenced);

  m_poDCoffsetImage = NULL;

  m_poNonunfImage->m_bReadApplyGonioMask = TRUE;
  m_poNonunfImage->nComputeApplyActiveMask();

  if ( (m_poNonunfImage->bIsAvailable()) && (m_poDarkImage->bIsAvailable()) )
    {
      // Look at nonunf image and make some adjustments

      int i, j, k;
      int nStat;
      unsigned short int uiPix;

      if (eImage_uI2 != m_poNonunfImage->nGetDataType())
        {
          cout << "INFO: non-uniformity image is not unsigned short int!\n";
        }

      nStat = 0;
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag1,
                                                         &m_a4nBadFlag[0]);
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag2,
                                                         &m_a4nBadFlag[1]);
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag3,
                                                         &m_a4nBadFlag[2]);
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag4,
                                                         &m_a4nBadFlag[3]);

      //  Look for presence of ms_sNonunfDenominator and ms_sNonunfNumerator keywords.
      //  The header should have only ONE of these, not both.

      i = m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfDenominator,
                                             &m_fDenominator);
      j = m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfNumerator,
                                             &m_fNumerator);

      //  Check for errors in the nonunf image header info

      if (0 == nStat)
        {
          if ( (0 == i) && (0 == j) )
            // Both keywords present, so invalid header
            nStat = 2;
          else if ( (0 == i) && (0.0 >= m_fDenominator) )
            // Invalid denominator,  so invalid header
            nStat = 3;
          else if ( (0 == j) && (0.0 >= m_fNumerator) )
            {
              // Invalid numerator,  so invalid header
              nStat = 4;
            }
        }
    if (0 != nStat)
      {
        // Something wrong with header keywords

        cout << "Error with mask/non-uniformity image header info!\n";
        m_eThe_State    = eNonunf_unknown_state;
        m_prfNonunfFactor = &Cnonunf::fUnknownNonunfFactor; // Pointer to Function
      }
    else
      {
        // All is OK with the header keywords

        if (0 == i)
          m_prfNonunfFactor = &Cnonunf::fDarkDenomNonunfFactor;
        else
          m_prfNonunfFactor = &Cnonunf::fDarkNumerNonunfFactor;

        // Loop through all pixels in nonunf image and set bad pixels to
        // negative values.  Bad pixels have values that are negative to
        // start with or are equal one of the NONUNF_BAD* flags

        // Assumption: nonunf image is SIGNED!!!

        (void) m_poNonunfImage->nSetNextPixel(0, 0);
        int nDim[2];
        int nBad   = 0;
        int nTotal = m_poNonunfImage->nGetDimensions(&nDim[0], &nDim[1]);
        m_poNonunfImage->nSetNextPixel(0, 0);
        for (j = 0; j < nDim[1]; j++)
          {
            for (i = 0; i < nDim[0]; i++)
              {
                uiPix = m_poNonunfImage->uiGetNextPixel();
                for (k = 0; (k < 4) && (0 < uiPix); k++)
                  {
                    if (   (uiPix == m_a4nBadFlag[k])
                        && (uiPix > 0) )
                      {
                        uiPix = 0;
                        // Modify the image
                        (void) m_poNonunfImage->nSetPixel(i, j, uiPix);
                      }
                  } // end of k loop
                if (0 == uiPix) nBad++;
              } // end of i loop
          } // end of j loop
        cout << "There were " << nTotal << " pixels, with " << nBad
             << " bad pixels in the nonunf file." << endl;
        m_eThe_State      = eNonunf_dark_nonunf_state;
      }
    }
  else
    {
      cout << "Error reading one or more mask/non-uniformity files\n";
      m_eThe_State      = eNonunf_unknown_state;
      m_prfNonunfFactor = &Cnonunf::fUnknownNonunfFactor; // Pointer to Function
      
      return 1;
    }

  // Probably should check that images are of the right size and also
  // really are nonunf and dark images!

  return (0);
}

Cnonunf::Cnonunf(const Cstring& sNonunfName, const Cstring& sDarkName,
                 const Cstring& sDCoffsetName)
{
  (void) nInitValues(sNonunfName, sDarkName, sDCoffsetName);
}

Cnonunf::Cnonunf(const Cstring& sNonunfName)
{
  (void) nInitValues(sNonunfName);
}

int
Cnonunf::nInitValues(const Cstring& sNonunfName, const Cstring& sDarkName,
                 const Cstring& sDCoffsetName)
{
  (void) nInitValues();
  cout << "...reading nonunf, dark and DCoffset image files..." << endl;
  
  m_poNonunfImage = m_pImagePool->poGetImage(sNonunfName, !c_bDeleteImageIfNotReferenced);

  m_poDarkImage = m_pImagePool->poGetImage(sDarkName, !c_bDeleteImageIfNotReferenced);

  m_poDCoffsetImage = m_pImagePool->poGetImage(sDCoffsetName, !c_bDeleteImageIfNotReferenced);
  
  if (   (m_poNonunfImage->bIsAvailable())
      && (m_poDarkImage->bIsAvailable())
      && (m_poDCoffsetImage->bIsAvailable()) )
    {
      // Look at nonunf image and make some adjustments

      int i, j, k;
      int nStat;
      unsigned short int uiPix;

      if (eImage_uI2 != m_poNonunfImage->nGetDataType())
        {
          cout << "INFO: non-uniformity image is not unsigned short int!\n";
        }

      nStat = 0;
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag1,
                                                         &m_a4nBadFlag[0]);
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag2,
                                                         &m_a4nBadFlag[1]);
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag3,
                                                         &m_a4nBadFlag[2]);
      nStat = nStat + m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfFlag4,
                                                         &m_a4nBadFlag[3]);
      //  Look for presence of NONUNF_DENOMINATOR and NONUNF_NUMERATOR keywords.
      //  The header should have only ONE of these, not both.

      i = m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfDenominator,
                                             &m_fDenominator);
      j = m_poNonunfImage->m_oHeader.nGetValue(ms_sNonunfNumerator,
                                             &m_fNumerator);

      //  Check for errors in the nonunf image header info

      if (0 == nStat)
        {
          if ( (0 == i) && (0 == j) )
            // Both keywords present, so invalid header
            nStat = 2;
          else if ( (0 == i) && (0.0 >= m_fDenominator) )
            // Invalid denominator,  so invalid header
            nStat = 3;
          else if ( (j == 0) && (0.0 >= m_fNumerator) )
            {
              // Invalid numerator,  so invalid header
              nStat = 4;
            }
        }
    if (0 != nStat)
      {
        // Something wrong with header keywords
        cout << "Error with mask/non-uniformity image header info!\n";
        m_eThe_State      = eNonunf_unknown_state;
        m_prfNonunfFactor = &Cnonunf::fUnknownNonunfFactor; // Pointer to Function
      }
    else
      {
        // All is OK with the header keywords

        if (0 == i)
          m_prfNonunfFactor = &Cnonunf::fDarkDenomNonunfFactor;
        else
          m_prfNonunfFactor = &Cnonunf::fDarkNumerNonunfFactor;

        // Loop through all pixels in nonunf image and set bad pixels to
        // negative values.  Bad pixels have values that are negative to
        // start with or are equal one of the NONUNF_BAD* flags

        (void) m_poNonunfImage->nSetNextPixel(0, 0);
        int nDim[2];
        int nBad   = 0;
        int nTotal = m_poNonunfImage->nGetDimensions(&nDim[0], &nDim[1]);
        m_poNonunfImage->nSetNextPixel(0, 0);
        for (j = 0; j < nDim[1]; j++)
          {
            for (i = 0; i < nDim[0]; i++)
              {
                uiPix = m_poNonunfImage->uiGetNextPixel();
                for (k = 0; (k < 4) && (0 < uiPix); k++)
                  {
                    if (   (uiPix == m_a4nBadFlag[k])
                        && (uiPix > 0) )
                      {
                        uiPix = 0;
                        // Modify the image
                        (void) m_poNonunfImage->nSetPixel(i, j, uiPix);
                      }
                  } // end of k loop
                if (0 == uiPix) nBad++;
              } // end of i loop
          } // end of j loop

        cout << "There were " << nTotal << " pixels, with " << nBad
             << " bad pixels in the nonunf file." << endl;
        m_eThe_State      = eNonunf_dark_dcoffset_state;
      }
    }
  else
    {
      cout << "Error reading one or more mask/non-uniformity files\n";
      m_eThe_State      = eNonunf_unknown_state;
      m_prfNonunfFactor = &Cnonunf::fUnknownNonunfFactor; // Pointer to Function
    
      return 1;
    }

  // Probably should check that images are of the right size and also
  // really are nonunf, dark and dcoffset images!

  return (0);
}

int Cnonunf::nInitValues(const Cstring& sNonunfName, Cimage_header *poHeader)
{
  // Called when simple_mask is the type
  (void) nInitValues();
  
  if (NULL != m_poNonunfImage)
  {
    delete      m_poNonunfImage;
    m_poNonunfImage = NULL;	
  }
  
  int           nStat = 1;
  
  Cstring       sTemp(sNonunfName);
  
  //+JWP 2008-02-07

  // This will only help if the image is 
  // in the original CBF format when read in below!
  // That may not be the case!
  // Since it may be a Pilates long int image, make sure that
  // image is read as unsigned short int

  Cstring sCBF = "";
  sCBF = sGetEnv(D_K_DtrekCBFPixDataType);

  // Try to be sure it is 2-bytes per pixel even if a CBF image

//+JWP 2011-01-26  
//  No!  We can now have signed long int because it will be converted below anyways, 
// so comment out this bit.
  /*****
	//cout << "Setting CBF to 2\n" << flush;
	//nPutEnv(D_K_DtrekCBFPixDataType, "2");
  *****/
//-JWP 2011-01-26  
  //-JWP 2008-02-07  

  if ( ("None" == sTemp) || ("FirstScanImage" == sTemp) )
    {
      // Simple_mask, but name is "None", so try to use the FIRST image
      // defined in the scan defined in the input header

      if ( (NULL != poHeader) && poHeader->bIsAvailable())
        {
          Cscan *poScan;
          poScan = new Cscan(*poHeader);
          if (poScan->bIsAvailable())
          {
              sTemp = poScan->sGetImageName();
              cout << "TINFO in Cnonunf: using\n     " << sTemp
                 << "\n     as the simple mask/nonunf file, was "
                   << sNonunfName << ".\n";
              
              if( THE_IMAGE_WAIT_PTR->bWaitRequired() && !THE_IMAGE_WAIT_PTR->bDoWaitForFirstImage(sTemp) )
		{
                    nStat = 1;    // couldn't open the image
		}
              else
		{
		  // This is where a mask image is read in.
		  m_poNonunfImage = m_pImagePool->poGetImage(sTemp, !c_bDeleteImageIfNotReferenced);
		}

              if( m_poNonunfImage && m_poNonunfImage->bIsAvailable() )
                  nStat = 0;
              else
                  nStat = 1;
          }
          
          delete poScan;
          
          poScan = NULL;
        }
      
      if (0 != nStat)
        {
          cout << "WARNING in Cnonunf: mask file " << sTemp
               << " not available.\n";
        }
    }
  else
    {
//+2010-02-03 JWP
      // Look for the text %2T and %DX in the mask filename.  
      // If found, replace with the integral 2theta value or the distance value.
      // Examples:
      //   beam.mask -> beam.mask
      //   beam%2T%DX.mask -> beam2T0DX70.mask   (2theta = 0.1, Dist=69.8)
      //   beam%2T_%DX.mask -> beam2T-10_DX100.mask (2theta = -0.3, Dist = 100.1) 

      Cstring s2T = "%2T";
      Cstring sDX = "%DX";
      Cstring sBefore, sAfter;
      Cstring sFallback = "";
      const int n2TIndex = 1;
      const int nDXIndex = 5;
      int nLocalStat = 1;  // 1 marks a failure
      float a6fDetGonio[6];
      sTemp = sNonunfName;
      if (sTemp.contains(s2T) || sTemp.contains(sDX) )
        {
          // Get the detector goniometer values directly from the header string if it is available

          if ( (NULL != poHeader) && poHeader->bIsAvailable())
            {
              // JWP: Unfortunately, this is a direct look into the header 
              //      and does not use encapsulation with the Cdetector and Cgoniometer classes.

              // Find 2THETA in the header and round-off to integer
              Cstring sDet = "";
              nLocalStat = poHeader->nGetValue(Cstring(D_K_DetectorNames), 1, &sDet);
              a6fDetGonio[n2TIndex] = 0.0;  // 2theta
              a6fDetGonio[nDXIndex] = 0.0;  // Distance 
              if (0 == nLocalStat)
                nLocalStat = poHeader->nGetValue(sDet+Cstring(D_K_GonioValues), 6, &a6fDetGonio[0]);
            }

          // A problem is that the refined values of the detector goniometer may make the nonunf filename
          // different, so we will want to allow for +-1 degree or +-1 mm in the file name, just in case.  
          // At this point if (0 == nLocalStat) we probably have a value for 2theta and a value for distance.
          // So want to create the filename, then see if the file exists.  
          // If the file exists, we are done.
          // If the file does not exist, we want to try slightly different 2theta and distance's to see if it
          // exists.
          int nTry;
          float af2x6[10][2] ={ { 0.0,  0.0},
                                { 0.0,  1.0},
                                { 0.0, -1.0},
                                { 1.0,  0.0},
                                { 1.0,  1.0},
                                { 1.0, -1.0},
                                {-1.0,  0.0},
                                {-1.0,  1.0},
                                {-1.0, -1.0},
                                { 0.0,  0.0} }; // Try 0 0 one more time because if we get to use this one,
                                                //    the error message will have the desired filename in it.
          float f2Theta = a6fDetGonio[n2TIndex];
          float fDist   = a6fDetGonio[nDXIndex];
          for (nTry = 0; nTry < 10; nTry++)
            {
              f2Theta = a6fDetGonio[n2TIndex] + af2x6[nTry][0];
              fDist   = a6fDetGonio[nDXIndex] + af2x6[nTry][1];

              if (sTemp.contains(s2T))
                {
                  sBefore = sTemp.before(s2T);
                  sAfter  = sTemp.after(s2T);
                  sFallback = sBefore + sAfter;
                  s2T = "";
                  if (0 == nLocalStat)
                    {
                      s2T = Cstring(f2Theta, 0, 0); // Convert to integer
                      s2T = s2T.remove(' ');        // Remove spaces
                      if ("-0" == s2T) s2T = "0";   // Remove minus sign on a zero.
                      s2T = Cstring("2T") + s2T;

                    }
                  sTemp = sBefore + s2T + sAfter;
                }
      
              if (sTemp.contains(sDX))
                {
                  sBefore = sTemp.before(sDX);
                  sAfter  = sTemp.after(sDX);
                  sFallback = sFallback.before(sDX) + sFallback.after(sDX);
                  sDX = "";
                  if (0 == nLocalStat)
                    {
                      sDX = Cstring(fDist, 0, 0); // Convert to integer
                      sDX = sDX.remove(' ');      // Remove spaces
                      sDX = Cstring("DX") + sDX;
                    }
                  sTemp = sBefore + sDX + sAfter;
                }
              if (bFileExists(sTemp))
                break;  // File exists break out of loop
            } // end of nTry
          if (10 <= nTry)
            {
              // File did not exist, but name has %2T or %DX in it.  
              // Remove these and see if the file exists.
              if (bFileExists(sFallback))
                {
                  cout << "WARNING, file " << sTemp << " does not exist,\n     but file " << sFallback
                       << " does, so use it.\n\n";
                  sTemp = sFallback;
                }
            }
        }  // end of block with %2T || %DX

      cout << "...reading mask file: " << sTemp << " ..." << endl;
      
      m_poNonunfImage = m_pImagePool->poGetImage(sTemp, !c_bDeleteImageIfNotReferenced);

//-2010-02-03 JWP
      
      if( !m_poNonunfImage->bIsAvailable() )  
        nStat = 1;
    }

  nStat = 1;
  if (NULL == m_poNonunfImage)
    {
      // Did not attemp to instantiate a nonunf image above
    }
  else if (m_poNonunfImage->bIsAvailable())
    {
      nStat = 0;
      if (eImage_uI2 != m_poNonunfImage->nGetDataType())
        {
          cout << "INFO: non-uniformity image is not unsigned short int!\n";
          if (eImage_I4 == m_poNonunfImage->nGetDataType())
            {
              cout << "      Attempting conversion ...\n";

              //+JWP 2008-02-07
              nStat = m_poNonunfImage->nConvertIntToUIShort();
              if (0 == nStat)
                cout << "      ... conversion successful.\n";
              else
                cout << "ERROR: Conversion UNSUCCESSFUL.\n";
              //-JWP 2008-02-07
            }
          else
            {
              cout << "SEVERE ERROR.\n";
              exit (-1);
            }
        }

      // While we are here, check and fix up some problems
      // with MAR mask images

      sTemp = "";
      m_poNonunfImage->m_oHeader.nGetValue(Cstring(D_K_DetectorNames), &sTemp);
      if ( ("MARIP_" == sTemp) || ("MARCCD_" == sTemp) )
        {
          // It is a MAR detector, so we need to flag inactive pixels
          // outside the circle, but they may not have a value of 0
          
          float fTemp = 0.0;
          m_poNonunfImage->m_oHeader.nGetValue(Cstring("MARCCD_DETECTOR_SIZE"),
                                               &fTemp);
          if (   ( (220.0 <= fTemp) && (230.0 >= fTemp) )
              || ( (290.0 <= fTemp) && (310.0 >= fTemp) ) )
            {
              // It is a MAR 225 CCD which is square
              // It is a MAR 300 CCD which is square
            }
          else
            {
              //cout << "MAR circle inscribed in square\n";
              unsigned short int uiTemp0, uiTemp1, uiTemp2;
              m_poNonunfImage->nGetPixel(0, 0, &uiTemp0);
              m_poNonunfImage->nGetPixel(0, 1, &uiTemp1);
              m_poNonunfImage->nGetPixel(0, 2, &uiTemp2);
              if (    (0 != uiTemp0) && (200 > uiTemp0)
                   && (uiTemp1 == uiTemp0)
                   && (uiTemp2 == uiTemp0) )
                {
                  cout << "INFO: MAR image with non-zero inactive pixels.\n"
                       << "      Setting values <= " << uiTemp2 
                       << " to 0 in the mask file.\n";
                  m_poNonunfImage->nSetNextPixel(0,0);
                  int nx, ny;
                  ny = m_poNonunfImage->nGetDimension();
                  uiTemp0 = 0;
                  for (nx = 0; nx < ny; nx++)
                    {
                      uiTemp1 = m_poNonunfImage->uiGetNextPixelNoInc();
                      if (uiTemp1 <= uiTemp2)
                        {
                          // Reset to 0
                          m_poNonunfImage->vSetNextPixel(uiTemp0);
                        }
                      else
                        {
                          m_poNonunfImage->vSetNextPixel(uiTemp1);
                        }
                    }
                }
            }
        }
  }
  
  if (0 == nStat)
  {
      m_eThe_State      = eNonunf_simple_mask_state;
      m_fScaleFactor    = 1.0;
      m_prfNonunfFactor = &Cnonunf::fSimpleNonunfFactor;  // Pointer to Function
  }
  else
  {
      m_eThe_State      = eNonunf_unknown_state;
      m_prfNonunfFactor = &Cnonunf::fUnknownNonunfFactor; // Pointer to Function
      cout << "ERROR with mask/non-uniformity image header info!\n";
      
      nPutEnv(D_K_DtrekCBFPixDataType, sCBF);
      return nStat;
  }
  // Restore whatever it was there (or if nothing, this unsetenv's it)
  nPutEnv(D_K_DtrekCBFPixDataType, sCBF);
  return 0;
}

Cnonunf::Cnonunf(Cimage_header& oHeader, const Cstring& sPrefix)
{
  (void) nInitValues(oHeader, sPrefix);
}

int
Cnonunf::nInitValues(Cimage_header& oHeader, const Cstring& sPre)
{
  Cstring sTemp[3];
  Cstring sPrefix;
  int nStat;
  int nTemp;

  nStat = 0;
  (void) nInitValues();

  sPrefix = sPre;
  if ("" == sPrefix)
    {
      // Null prefix string, see if the Cdetector::ms_sDetectorNames keyword
      // exists in the header.  If so, change the prefix to the first name
      // in the list of names.  This call is OK, since if there is no
      // keyword, then sPrefix is set to "" anyways.

      (void) oHeader.nGetValue(D_K_DetectorNames, 1, &sPrefix);
    }
  if (!oHeader.bIsAvailable())
    {
      cout << "Cnonunf ERROR: image header is not valid.  Cannot construct mask/nonunf object!\n";
      nStat = 1;
    }
  else if (0 == oHeader.nGetValue(sPrefix + ms_sNonunfKey, &sTemp[0]))
    {
      cout << "Cnonunf ERROR: cannot use database key yet!\n";
      nStat = 2;
    }
        
  // Try to get required information from the header

  nTemp = oHeader.nGetValue(sPrefix + ms_sNonunfType, &sTemp[0]);
  //debug
  cout << sPrefix << ms_sNonunfType << ": >>" << sTemp[0] << "<<\n";
  sTemp[0] = sTransSymbol(sTemp[0]);

  if (0 == nTemp)
    {
      if (ms_sNonunfStateUnknown == sTemp[0])
        {
          m_eThe_State = eNonunf_unknown_state;
        }
      else if (ms_sNonunfStateNone == sTemp[0])
        {
          m_eThe_State = eNonunf_none_state;
          m_prfNonunfFactor = &Cnonunf::fNoneNonunfFactor;// Pointer to Function
        }
      else if (ms_sNonunfStateSimpleScale == sTemp[0])
        {
          m_fScaleFactor = 1.0;
          nTemp = oHeader.nGetValue(sPrefix + ms_sNonunfInfo, &m_fScaleFactor);
          if (0 == nTemp)
            {
              m_eThe_State = eNonunf_simple_scale_state;

              // Pointer to Function

              m_prfNonunfFactor = &Cnonunf::fSimpleNonunfFactor;
            }
          else
            {
              cout << "Cnonunf: error reading " << sPrefix << ms_sNonunfInfo
                   << " from header!\n";
              nStat++;
            }   
        }
      else if( ms_sNonunfStateSimpleMask == sTemp[0] )
      {
            nTemp = oHeader.nGetValue(sPrefix + ms_sNonunfInfo, &sTemp[0]);
      
            if( 0 == nTemp && sTemp[0].length() > 0 )
            {
                m_eThe_State = eNonunf_simple_mask_state;
                
                if( 0 != nInitValues(sTemp[0], &oHeader) )
                {
                    cout << "Cnonunf: error reading simple Mask/Nonunf file name from header!\n";
                    nStat++;
                }   
            }
            else
            {
                cout << "Cnonunf: error reading simple Mask/Nonunf file name from header!\n";
                nStat++;
            }   
      }
      else if (ms_sNonunfStateDarkNonunf == sTemp[0])
        {
          int nFiles = 2;
          nTemp = oHeader.nGetValue(sPrefix + ms_sNonunfInfo, nFiles, sTemp);
          if (0 == nTemp && sTemp[0].length() > 0 && sTemp[1].length() > 0 )
            {
              (void) nInitValues(sTemp[0], sTemp[1]);
            }
          else
            {
              cout << "Cnonunf: error reading Nonunf and Dark file names from header!\n";
              nStat++;
            }   
        }
      else if (ms_sNonunfStateDarkDCoffsetNonunf == sTemp[0])
        {
          int nFiles = 3;
          nTemp = oHeader.nGetValue(sPrefix + ms_sNonunfInfo, nFiles, sTemp);
          if (0 == nTemp && sTemp[0].length() > 0 && sTemp[1].length() > 0 && sTemp[2].length() > 0 )
            {
              (void) nInitValues(sTemp[0], sTemp[1], sTemp[2]);
            }
          else
            {
              cout << "Cnonunf: error reading Nonunf, Dark and DCoffset file names from header!\n";
              nStat++;
            }   
        }
      else
        {
          cout << "Cnonunf: invalid " << sPrefix << ms_sNonunfInfo << " in header!\n";
          nStat++;
        }
    }
  else
    {
      cout << "CNonunf ERROR: reading " << sPrefix << ms_sNonunfType << " keyword!\n";
      nStat++;
    }
  m_fMinRawPixOKValue = 1.0;
  nTemp = oHeader.nGetValue("DTREK_NONUNF_OKVALUE", &m_fMinRawPixOKValue);
  cout << "Min raw image pixel OK value in mask/nonunf/image file: " << m_fMinRawPixOKValue << endl;

  if (0 != nStat) (void) nInitValues(); // Use default constructor on all errors
  return (0);
}

Cnonunf::~Cnonunf()
{
  //  cout << "Cnonunf::destructor called!\n";

  // Be sure to 'delete' any members allocated with 'new' here

  if (NULL != m_poNonunfImage)
    {
      m_pImagePool->bFreeImage(m_poNonunfImage->sGetName(), c_bDeleteImageIfNotReferenced);
      m_poNonunfImage = NULL;
    }

  if (NULL != m_poDarkImage)
    {
      m_pImagePool->bFreeImage(m_poDarkImage->sGetName(), c_bDeleteImageIfNotReferenced);
      m_poDarkImage = NULL;
    }

  if (NULL != m_poDCoffsetImage)
    {
      m_pImagePool->bFreeImage(m_poDCoffsetImage->sGetName(), c_bDeleteImageIfNotReferenced);
      m_poDCoffsetImage = NULL;
    }
}

int Cnonunf::nInitValues(void)
{
  //  cout << "Cnonunf::nInitValues called\n";
  m_eThe_State      = eNonunf_unknown_state;
  m_sNonunf_key     = "Unknown";
  m_nMethod         = 0;
  m_fScaleFactor    = 1.0;
  m_fMinFactorLimit = 0.0;
  m_fMaxFactorLimit = 32767.0;
  m_fNumerator      = 0.0;
  m_fDenominator    = 0.0;
  m_a4nBadFlag[0]   = 0;
  m_a4nBadFlag[1]   = 0;
  m_a4nBadFlag[2]   = 0;
  m_a4nBadFlag[3]   = 0;
  m_poNonunfImage   = NULL;
  m_poDarkImage     = NULL;
  m_poDCoffsetImage = NULL;
  m_prfNonunfFactor = &Cnonunf::fUnknownNonunfFactor;  // Pointer to Function
  m_fBadFlag        = -99999.00;  // Needs to be a very large NEGATIVE value

  // m_fMinRawPixOKValue sets the lowest OK value for a raw pixel in an image.
  // Note that this value is for the RAW image and NOT the NONUNF image or mask.
  // The minimum OK value for a valid pixel in the NONUNF image or mask is still 0.

  m_fMinRawPixOKValue  =  1.0;  // Could be changed to 0.0 for some CBF-style images later

  if ("" != sGetEnv("DTREK_NONUNF_OKVALUE"))
    {
      m_fMinRawPixOKValue  =  atof(sGetEnv("DTREK_NONUNF_OKVALUE").string());
      if (1.0 != m_fMinRawPixOKValue)
	cout << "INFO: Min raw image pixel OK value in image file: " 
             << m_fMinRawPixOKValue 
	     << "\nNOTE: nonunf/mask min value is kept at 1.\n";
    }
  
  m_pImagePool =  CImagePool::poGetInstance();
  
  return (0);
}

int Cnonunf::nList(void)
{
  cout << "Nonunf key:  " << m_sNonunf_key
       << "\nNonunf state: ";
  cout << "Minimum raw (before correction) pixel value that is still OK: "
       << m_fMinRawPixOKValue << endl;
  if (eNonunf_unknown_state == m_eThe_State)
    cout << ms_sNonunfStateUnknown;
  else if (eNonunf_none_state == m_eThe_State)
    cout << ms_sNonunfStateNone;
  else if (eNonunf_simple_scale_state == m_eThe_State)
    cout << ms_sNonunfStateSimpleScale << ": " << m_fScaleFactor;
  else if (eNonunf_dark_nonunf_state == m_eThe_State)
    {
      cout << ms_sNonunfStateDarkNonunf;
      cout << "\nNonunf image and dark image will be used"
           << "\nNon-uniformity image is: ";
      if (NULL == m_poNonunfImage)
        cout << "undefined";
      else
        cout << m_poNonunfImage->sGetName();

      cout << "\nDark image is: ";
      if (NULL == m_poDarkImage)
        cout << "undefined";
      else
        cout << m_poDarkImage->sGetName();
    }
  else if (eNonunf_dark_dcoffset_state == m_eThe_State)
    {
      cout << ms_sNonunfStateDarkDCoffsetNonunf;
      cout << "\nDCoffset and dark images will be used";
      cout << "\nNon-uniformity image is: ";
      if (NULL == m_poNonunfImage)
        cout << "undefined";
      else
        cout << m_poNonunfImage->sGetName();

      cout << "\nDark image is: ";
      if (NULL == m_poDarkImage)
        cout << "undefined";
      else
        cout << m_poDarkImage->sGetName();

      cout << "\nDC offset image is: ";
      if (NULL == m_poDCoffsetImage)
        cout << "undefined";
      else
        cout << m_poDCoffsetImage->sGetName();
    }
  cout << endl;
  return (0);
}

bool Cnonunf::bPixelIsBad(const int nPx0, const int nPx1)
{
  unsigned short int uiPix;
  int nStat;
  if (eNonunf_unknown_state == m_eThe_State)
    return (TRUE);   // It is BAD!
  else if (eNonunf_none_state == m_eThe_State)
    return (FALSE);  // It is always GOOD!

  // DANGER: Next line assumes NonunfImage is UNSIGNED SHORT INT.

  nStat = m_poNonunfImage->nGetPixel(nPx0, nPx1, &uiPix);
  return ( (0 != nStat) || (0 == uiPix) );
}

float Cnonunf::fDarkDenomNonunfFactor(const int nPx0, const int nPx1)
{
  float fPix;
  fPix = (m_poNonunfImage->*m_poNonunfImage->prfGetPixel)(nPx0, nPx1);

  if (0.0 < fPix)
    return ( fPix / m_fDenominator);
  else
    {
      // 0 or negative is a bad nonunf pixel

      return (m_fBadFlag);  // Should RETURN badpix flag perhaps
    }
}

float Cnonunf::fDarkNumerNonunfFactor(const int nPx0, const int nPx1)
{
  float fPix;
  fPix = (m_poNonunfImage->*m_poNonunfImage->prfGetPixel)(nPx0, nPx1);

  if (0.0 < fPix)
    return ( m_fNumerator / fPix);
  else
    {
      // 0 or negative is a bad nonunf pixel

      return (m_fBadFlag);  // Should RETURN badpix flag perhaps
    }
}

float Cnonunf::fGetMaskValue(const int nPx0, const int nPx1)
{
  return ((m_poNonunfImage->*m_poNonunfImage->prfGetPixel)
                              (nPx0, nPx1));
}

float Cnonunf::fSimpleNonunfFactor(const int nPx0, const int nPx1)
{
  return (m_fScaleFactor);
}

float Cnonunf::fNoneNonunfFactor(const int nPx0, const int nPx1)
{
  return (1.0);
}

float Cnonunf::fUnknownNonunfFactor(const int nPx0, const int nPx1)
{
  return (0.0);  // Really make sure if unknown that any attempted use yields
                 //  garbage!
}

int
Cnonunf::nCorrectImage(Cimage *poImage)
{
  int i, nNumPixels;
  unsigned short int uiNonunf;

  if (eNonunf_none_state == m_eThe_State)
    {
      // There is no non-uniformity nor bad pixels so the image needs
      //   no correction!

      // Presumably bad pixels were set to 0 by any 
      // external embedded mask routine (see the Cimage class)

      return (0);
    }
  else if (eNonunf_simple_mask_state == m_eThe_State)
    {
      if (eImage_uI2 == poImage->nGetDataType())
	{
	  // For speed treat the unsigned short int as a special case
	  unsigned short int uiImage;

	  nNumPixels = poImage->nGetDimension();
	  m_poNonunfImage->nSetNextPixel(0,0);
	  poImage->nSetNextPixel(0,0);
	  for (i = 0; i < nNumPixels; i++)
	    {
	      uiNonunf = m_poNonunfImage->uiGetNextPixel();
	      uiImage  = poImage->uiGetNextPixelNoInc();
	      if (0 == uiNonunf)
		uiImage = 0;
	      poImage->vSetNextPixel(uiImage);
	    }
	}
      else
	{
	  // All other will use the slow method
	  float fPixVal;
	  int nDim0, nDim1;
	  int j;

	  nNumPixels = poImage->nGetDimension();
	  poImage->nGetDimensions(&nDim0, &nDim1);
	  m_poNonunfImage->nSetNextPixel(0,0);
	  for (j = 0; j < nDim1; j++)
	    {
	      for (i = 0; i < nDim0; i++)
		{
		  uiNonunf = m_poNonunfImage->uiGetNextPixel();
		  fPixVal = (poImage->*poImage->prfGetPixel)(i, j);
		  if (0 == uiNonunf)
		    fPixVal = 0;
		  (poImage->*poImage->prnSetPixel)(i, j, fPixVal);
		}
	    }
	}
    }
  else if (eNonunf_dark_nonunf_state == m_eThe_State)
    {
      // Subtract off Dark current image, apply nonunf scale factor and
      // set any bad pixels to 0.  This does not check for
      // saturated pixels.  It also assumes the images are unsigned short
      // int.

      int   nDim0, nDim1;
      int i, nNumPixels;
      float fPixel, fNonunf;
      float fSatVal    = poImage->fGetSatValue();
      unsigned short int uiDark, uiNonunf, uiImage;

      if (eImage_uI2 != poImage->nGetDataType())
        {
          cout << "WARNING in Cnonunf::nCorrectImage(), unsupported image data type!\n";
          return (-1);
        }
      nNumPixels = poImage->nGetDimension();
      poImage->nGetDimensions(&nDim0, &nDim1);
      m_poNonunfImage->nSetNextPixel(0,0);
      m_poDarkImage->nSetNextPixel(0,0);
      poImage->nSetNextPixel(0,0);

      if (0.0 == m_fDenominator)
        {
          for (i = 0; i < nNumPixels; i++)
            {
              uiNonunf = m_poNonunfImage->uiGetNextPixel();
              uiDark   = m_poDarkImage->uiGetNextPixel();
              uiImage  = poImage->uiGetNextPixelNoInc();
              if (0 >= uiNonunf)
                {
                  uiImage = 0;
                }
              else
                {
                  fNonunf = m_fNumerator / (float)uiNonunf;
                  fPixel = (float)(uiImage - uiDark) * fNonunf;
                  if (0.0 > fPixel) fPixel = 0.0;
                  uiImage = (unsigned short int) fPixel;
                }
              poImage->vSetNextPixel(uiImage);
            }

        }
      else
        {
          for (i = 0; i < nNumPixels; i++)
            {
              uiNonunf = m_poNonunfImage->uiGetNextPixel();
              uiDark   = m_poDarkImage->uiGetNextPixel();
              uiImage  = poImage->uiGetNextPixelNoInc();
              if (0 >= uiNonunf)
                {
                  uiImage = 0;
                }
              else
                {
                  fNonunf = (float)uiNonunf / m_fDenominator;
                  fPixel = (float)(uiImage - uiDark) * fNonunf;
                  if (0.0 > fPixel) fPixel = 0.0;
                  uiImage = (unsigned short int) fPixel;
                }
              poImage->vSetNextPixel(uiImage);
            }
        }
      return (0);
    }
  else if (eNonunf_dark_dcoffset_state == m_eThe_State)
    {
      // Subtract off DCoffset and Dark current image and
      //    flag any bad pixels as bad

      return (1); // Not implemented
    }
  else if (eNonunf_simple_scale_state == m_eThe_State)
    {
      // Scale all the pixels and flag any bad pixels as bad

      return (1); // Not implemented
    }
  else
    {
      // Non-uniformity information is not available, so no correction is
      //   possible, so return error.

      cout << "Cnonunf::nCorrectImage - unknown nonunf state!\n";
      return (1);
    }
  return (0);  // This keeps VC++5.0 from throwing error
}

int
Cnonunf::nCorrect3Ddata(C3Ddata *po3Ddata, float *pfPerLayerFactors)
{
  int i, j, k; // Loop counters

  float fDark;
  float fNonunf;
  int   nPx0;
  int   nPx1;
  int   nPixPerLayer;
  float *pfTemp1, *pfTemp2;
  float fSatVal;
  float *pfFactor;
  float *pfLayerFactor;
  int   nNumBadPix =  0;
  int   nNumSatPix =  0;
  int   nStat = 0;

  pfLayerFactor = new float [po3Ddata->nGetExt(2)];

  for (i = 0; i < po3Ddata->nGetExt(2); i++)
    {
      if (NULL == pfPerLayerFactors)
        pfLayerFactor[i] = 1.0;      // Layer factors not supplied, use 1.0
      else
        pfLayerFactor[i] = pfPerLayerFactors[i];
    }

  if (eNonunf_none_state == m_eThe_State)
    {
      // There is no non-uniformity nor bad pixels so the image needs
      //   no correction!

      if (e3Ddata_float != po3Ddata->eGetType())
	nStat = po3Ddata->nConvertToFloat();

      //TJN: This loop was commented out, but has now been reinstated.
      //     We now regard a pixel value below m_fMinRawPixOKValue in the 
      //        uncorrected shoebox as a bad pixel.
      //     These pixels might have been masked out by the active mask.

      if (0 == nStat)
        {
          fSatVal = po3Ddata->m_poImage->fGetSatValue();
          pfTemp1 = po3Ddata->m_pfData;
          i       = po3Ddata->nGetExt(0)
                    * po3Ddata->nGetExt(1) * po3Ddata->nGetExt(2);
          for (j = 0; j < i; j++)
            {
              if (*pfTemp1 < m_fMinRawPixOKValue)
                {
                  *pfTemp1 = m_fBadFlag;
                } 
              else if (fSatVal <= *pfTemp1)
                {
                  nNumSatPix++;
                }
              pfTemp1++;
            }
        }
    }
  else if (eNonunf_dark_nonunf_state == m_eThe_State)
    {
      if (e3Ddata_float != po3Ddata->eGetType())
	nStat = po3Ddata->nConvertToFloat();
      if (0 == nStat)
        {
          fSatVal      = po3Ddata->m_poImage->fGetSatValue();
          nPixPerLayer = po3Ddata->nGetExt(0) * po3Ddata->nGetExt(1);
          pfTemp1      = po3Ddata->m_pfData;

          nPx1         = po3Ddata->nGetOffset(1);
          for (j = 0; j < po3Ddata->nGetExt(1); j++)
            {
              nPx0 = po3Ddata->nGetOffset(0);
              for (i = 0; i < po3Ddata->nGetExt(0); i++)
                {
                  fNonunf = (this->*m_prfNonunfFactor)(nPx0, nPx1);
                  pfTemp2 = pfTemp1;
                
                  if ( (fNonunf == m_fBadFlag) || (0 >= fNonunf) )
                    {
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
                          nNumBadPix++;
                          *pfTemp2 = m_fBadFlag;
                          pfTemp2  = pfTemp2 + nPixPerLayer;
                        }
                    }
                  else
                    {
                      fDark  = (m_poDarkImage->*m_poDarkImage->prfGetPixel)
                               (nPx0, nPx1);
                      pfFactor = pfLayerFactor;
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
			  if (*pfTemp2 < m_fMinRawPixOKValue)
                            {
                              *pfTemp2 = m_fBadFlag;
                              nNumBadPix++;
                            } 
                          else if (*pfTemp2 < fSatVal)
                            {
                              *pfTemp2 = *pfFactor
                                         * (*pfTemp2 - fDark) * fNonunf;
                            }
                          else
                            {
                              nNumSatPix++;
                            }
                          pfFactor++;
                          pfTemp2 = pfTemp2 + nPixPerLayer;
                        }
                    }
                  pfTemp1++;  // Next pixel in layer
                  nPx0++;     // Next pixel in nonunf and dark arrays
                }  // end of i loop
              nPx1++;
            } // end of j loop
        }
    }
  else if (eNonunf_dark_dcoffset_state == m_eThe_State)
    {
      if (e3Ddata_float != po3Ddata->eGetType())
	nStat = po3Ddata->nConvertToFloat();
      if (0 == nStat)
        {
          // Subtract off DCoffset and Dark current image and
          //    flag any bad pixels as bad

          fSatVal      = po3Ddata->m_poImage->fGetSatValue();
          nPixPerLayer = po3Ddata->nGetExt(0) * po3Ddata->nGetExt(1);
          pfTemp1      = po3Ddata->m_pfData;

          nPx1         = po3Ddata->nGetOffset(1);
          for (j = 0; j < po3Ddata->nGetExt(1); j++)
            {
              nPx0 = po3Ddata->nGetOffset(0);
              for (i = 0; i < po3Ddata->nGetExt(0); i++)
                {
                  fNonunf = (this->*m_prfNonunfFactor)(nPx0, nPx1);
                  pfTemp2 = pfTemp1;
        
                  if ( (fNonunf == m_fBadFlag) || (0 >= fNonunf) )
                    {
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
                          *pfTemp2 = m_fBadFlag;
                          pfTemp2  = pfTemp2 + nPixPerLayer;
                        }
                    }
                  else
                    {
                      fDark  = (m_poDarkImage->*m_poDarkImage->prfGetPixel)
                                                               (nPx0, nPx1);
                      pfFactor = pfLayerFactor;
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
			  if (*pfTemp2 < m_fMinRawPixOKValue)
                            {
                              *pfTemp2 = m_fBadFlag;
                              
                            } 
                          else if (*pfTemp2 < fSatVal)
                            {
                              // What about DC offset??
                              
                              *pfTemp2 = *pfFactor
                                         * (*pfTemp2 - fDark) * fNonunf;
                            }
                          pfFactor++;
                          pfTemp2 = pfTemp2 + nPixPerLayer;
                        }
                    }
                  pfTemp1++;  // Next pixel in layer
                  nPx0++;     // Next pixel in nonunf array
                }  // end of i loop
              nPx1++;
            } // end of j loop
        }
    }
  else if (eNonunf_simple_scale_state == m_eThe_State)
    {
      if (e3Ddata_float != po3Ddata->eGetType())
	nStat = po3Ddata->nConvertToFloat();
      if (0 == nStat)
        {
          // Scale all the pixels and flag any bad pixels as bad

          fSatVal      = po3Ddata->m_poImage->fGetSatValue();
          nPixPerLayer = po3Ddata->nGetExt(0) * po3Ddata->nGetExt(1);
          pfTemp1      = po3Ddata->m_pfData;

          nPx1         = po3Ddata->nGetOffset(1);
          for (j = 0; j < po3Ddata->nGetExt(1); j++)
            {
              nPx0 = po3Ddata->nGetOffset(0);
              for (i = 0; i < po3Ddata->nGetExt(0); i++)
                {
                  fNonunf = (this->*m_prfNonunfFactor)(nPx0, nPx1);
                  pfTemp2 = pfTemp1;
        
                  if ( (fNonunf == m_fBadFlag) || (0 >= fNonunf) )
                    {
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
                          *pfTemp2 = m_fBadFlag;
                          pfTemp2  = pfTemp2 + nPixPerLayer;
                        }
                    }
                  else
                    {
                      fDark  = (m_poDarkImage->*m_poDarkImage->prfGetPixel)
                                                               (nPx0, nPx1);
                      pfFactor = pfLayerFactor;
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
			  if (*pfTemp2 < m_fMinRawPixOKValue)
                            {
                              *pfTemp2 = m_fBadFlag;
                            } 
                          else if (*pfTemp2 < fSatVal)
                            {
                              // What about DC offset??
                              *pfTemp2 = *pfFactor
                                         * (*pfTemp2 - fDark) * fNonunf;
                            }
                          pfFactor++;
                          pfTemp2 = pfTemp2 + nPixPerLayer;
                        }
                    }
                  pfTemp1++;  // Next pixel in layer
                  nPx0++;     // Next pixel in nonunf array
                }  // end of i loop
              nPx1++;
            } // end of j loop
        }
    }
  else if (eNonunf_simple_mask_state == m_eThe_State)
    {
      if (e3Ddata_float != po3Ddata->eGetType())
	nStat = po3Ddata->nConvertToFloat();
      if (0 == nStat)
        {
          fSatVal      = po3Ddata->m_poImage->fGetSatValue();
          nPixPerLayer = po3Ddata->nGetExt(0) * po3Ddata->nGetExt(1);
          pfTemp1      = po3Ddata->m_pfData;

          nPx1         = po3Ddata->nGetOffset(1);
          for (j = 0; j < po3Ddata->nGetExt(1); j++)
            {
              nPx0 = po3Ddata->nGetOffset(0);
              for (i = 0; i < po3Ddata->nGetExt(0); i++)
                {
                  if( !m_poNonunfImage )
                  {
                     nStat = 1;
                     goto leave;
                  }

                  fNonunf = (m_poNonunfImage->*m_poNonunfImage->prfGetPixel)
                              (nPx0, nPx1);
                  pfTemp2 = pfTemp1;
                
                  if ( (fNonunf == m_fBadFlag) || (0 >= fNonunf) )
                    {
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
                          nNumBadPix++;
                          *pfTemp2 = m_fBadFlag;
                          pfTemp2  = pfTemp2 + nPixPerLayer;
                        }
                    }
                  else
                    {
                      pfFactor = pfLayerFactor;
                      for (k = 0; k < po3Ddata->nGetExt(2); k++)
                        {
			  if (*pfTemp2 < m_fMinRawPixOKValue)
                            {
                              *pfTemp2 = m_fBadFlag;
                              nNumBadPix++;
                            } 
                          else if (*pfTemp2 < fSatVal)
                            {
                              *pfTemp2 = *pfFactor * (*pfTemp2);
                            }
                          else
                            {
                              nNumSatPix++;
                            }
                          pfFactor++;
                          pfTemp2 = pfTemp2 + nPixPerLayer;
                        }
                    }
                  pfTemp1++;  // Next pixel in layer
                  nPx0++;     // Next pixel in nonunf and dark arrays
                }  // end of i loop
              nPx1++;
            } // end of j loop
        }
    }
  else
    {
      // Non-uniformity information is not available, so no correction is
      //   possible, so return error.

      cout << "Cnonunf::nCorrect3Ddata - unknown nonunf state!\n";
      nStat = 1;
    }

  if (0 < nNumBadPix)
    nStat = nStat + ms_nBadPixFlag;
  if (0 < nNumSatPix)
    nStat = nStat + ms_nSatPixFlag;

leave:
  delete [] pfLayerFactor;
  
  return (nStat);
}

int
Cnonunf::nCorrectPixel(const int nPx0, const int nPx1, Cimage *poImage,
                       float *pfValue, const bool bCheckBounds)
{
  // Get pixel (nPx0, nPx1) out of poImage and correct for non-uniformity
  // of response.  Return corrected value in *pfValue.  Return 0 on success,
  // otherwise on bad-nonuniformity or out-of-bounds

  int   nDim0, nDim1;
  float fPixel;
  float fNonunf;
  float fDark;

  *pfValue = 0.0;

  if (bCheckBounds)
    {
      // Check bounds

      (void) poImage->nGetDimensions(&nDim0, &nDim1);
      if (   (0 > nPx0)
          || (0 > nPx1)
          || (nPx0 >= nDim0)
          || (nPx1 >= nDim1))
        return (3);                 // Out-of-bounds
    }

  // Get pixel value as floating point

  fPixel = (poImage->*poImage->prfGetPixel)(nPx0, nPx1);

  if (eNonunf_none_state == m_eThe_State)
    {
      // There is no non-uniformity nor bad pixels so the pixel needs
      //   no correction!

      *pfValue = fPixel;
      return (0);
    }
  else if (eNonunf_dark_nonunf_state == m_eThe_State)
    {
      fNonunf = (this->*m_prfNonunfFactor)(nPx0, nPx1);
      if ( (fNonunf == m_fBadFlag) || (0 >= fNonunf) )
        {
          *pfValue = m_fBadFlag;
          return (1);
        }
      else
        {
          fDark   = (m_poDarkImage->*m_poDarkImage->prfGetPixel)(nPx0, nPx1);
          if (fPixel < poImage->fGetSatValue())
            *pfValue = (fPixel - fDark) * fNonunf;
          return (0);
        }
    }
  else if (eNonunf_dark_dcoffset_state == m_eThe_State)
    {
      // Subtract off DCoffset and Dark current image

      fNonunf = (this->*m_prfNonunfFactor)(nPx0, nPx1);
      if ( (fNonunf == m_fBadFlag) || (0 >= fNonunf) )
        {
          *pfValue = m_fBadFlag;
          return (1);
        }
      else
        {
          fDark   = (m_poDarkImage->*m_poDarkImage->prfGetPixel)(nPx0, nPx1);
          if (fPixel < poImage->fGetSatValue())
            *pfValue = (fPixel - fDark) * fNonunf;
          return (0);
        }
    }
  else if (eNonunf_simple_mask_state == m_eThe_State)
    {
      fNonunf = (this->*m_prfNonunfFactor)(nPx0, nPx1);
      if ( (fNonunf == m_fBadFlag) || (0 >= fNonunf) )
        {
          *pfValue = m_fBadFlag;
          return (1);
        }
      else
        {
          *pfValue = fPixel;
          return (0);
        }
    }
  else if (eNonunf_simple_scale_state == m_eThe_State)
    {
      // Scale the pixels and flag any bad pixels as bad

      cout << "Cnonunf::correct simple not implemented yet!\n";
      *pfValue = fPixel;
      return (2);
    }
  else if (eNonunf_unknown_state == m_eThe_State)
    {
      // Non-uniformity information is not available, so no correction is
      //   possible, so return error.

      cout << "Cnonunf::nCorrectPixel - unknown nonunf state!\n";
      *pfValue = fPixel;
      return (1);
    }
  return (0);
}
