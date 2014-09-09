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
// Cscan.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class Cscan which implements
//    the scan encapsulation of d*TREK.
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
#include "Cscan.h"         // Class definition and prototypes
#include "dtrekdefs.h"
//+Definitions, constants, and initialization of static member variables

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


Cstring Cscan::ms_sScanPrefix     = D_K_ScanPrefix;
Cstring Cscan::ms_sScanTemplate   = D_K_ScanTemplate;
Cstring Cscan::ms_sScanKey        = D_K_ScanKey;
Cstring Cscan::ms_sScanSeqInfo    = D_K_ScanSeqInfo;
Cstring Cscan::ms_sScanWavelength = D_K_ScanWavelength;
Cstring Cscan::ms_sScanWavelengthOpts = D_K_ScanWavelengthOpts;
Cstring Cscan::ms_sScanCrysDatum  = D_K_ScanCrysDatum;
Cstring Cscan::ms_sScanDetDatum   = D_K_ScanDetDatum;
Cstring Cscan::ms_sScanDezinger   = D_K_ScanDezinger;
Cstring Cscan::ms_sScanMode       = D_K_ScanMode;
Cstring Cscan::ms_sScanModeUnknown= D_K_ScanModeUnknown;
Cstring Cscan::ms_sScanModeStillO = D_K_ScanModeStillO;
Cstring Cscan::ms_sScanModeStillC = D_K_ScanModeStillC;
Cstring Cscan::ms_sScanModeScanO  = D_K_ScanModeScanO;
Cstring Cscan::ms_sScanModeScanC  = D_K_ScanModeScanC;
Cstring Cscan::ms_sScanDetBinMode = D_K_ScanDetBinMode;
Cstring Cscan::ms_sScanDetectorOptions = D_K_ScanDetectorOptions;


//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cscan::Cscan()
{
 (void) nInitValues();
}

Cscan::Cscan(const Cstring sTemplate, const int nSeqNum = 0,
	     const int nStepNum = 1)
{
  (void) nInitValues();
  m_sFileTemplate    = sTemplate;
  m_nFileStartSeqNum = nSeqNum;
  m_nFileStepNum     = nStepNum;
  m_nFileSeqNum      = m_nFileStartSeqNum;
}

Cscan::Cscan(Cimage_header& oHeader, const Cstring &sPre)
                                     // Construct a scan by looking at keywords
{                                    // in an image header
  //bool bValue;
  int nStat, nTemp, i;
  int a2nOpts[2];
  Cstring sTemp;
  Cstring sPrefix;

  sPrefix = sPre;

  nStat = nInitValues();
  if (!oHeader.bIsAvailable())
    {
      cout << "Cscan ERROR: image header is not valid.  Cannot construct scan!\n";
      nStat = 1;
    }
  else if (oHeader.nGetValue(sPrefix + ms_sScanKey, &sTemp) == 0)
    {
      cout << "Cscan ERROR: cannot use database key yet!\n";
      nStat = 2;
    }
  else
    {
      // Try to get required information from the header

      float fValues[10];
      delete m_poRotation;   // Delete m_poRotation made by nInitValues
      m_poRotation = new Crotation(oHeader, sPrefix + ms_sScanPrefix);

      // Now try to read some other scan information

      nTemp = oHeader.nGetValue(sPrefix + ms_sScanTemplate, &sTemp);
      if (0 == nTemp)
	{
	  m_sFileTemplate = sTemp;
	}
      else
	{
	  nTemp = oHeader.nGetValue(Cimage_header::ms_sFilename, &sTemp);
	  if (0 == nTemp)
	    {
	      m_sFileTemplate = sBuildScanTemplate(sTemp);

//	      cout << "Cscan ScanTemplate is: " << m_sFileTemplate << endl;
	    }
	  else
	    {
	      cout << "Cscan WARNING: no " << sPrefix + ms_sScanTemplate << " keyword!\n";
	      cout << "Cscan WARNING: no FILENAME keyword!\n";
	      nStat++;
	    }
	}
      nTemp = oHeader.nGetValue(sPrefix + ms_sScanSeqInfo, 3, fValues);
      if (0 == nTemp)
	{
	  m_nFileStartSeqNum = (int) fValues[0];
	  m_nFileStepNum     = (int) fValues[1];
	  m_nNumImages       = (int) fValues[2];
	  m_nFileSeqNum      = m_nFileStartSeqNum;
	}
      else
	{
	  cout << "Cscan WARNING: no " << sPrefix + ms_sScanSeqInfo << " keyword!\n";
	  nStat++;
	}

      nTemp = oHeader.nGetValue(sPrefix + ms_sScanWavelength, 1, fValues);
      if (0 == nTemp)
	{
	  m_fWavelength = fValues[0];
	}
      else
	{
	  // Remain quiet about no wavelength
	}
      nTemp = oHeader.nGetValue(sPrefix + ms_sScanWavelengthOpts, 2, a2nOpts);
      if (0 == nTemp)
	{
	  m_nWavelengthOption   = a2nOpts[0];
	  m_nWavelengthOptimize = a2nOpts[1];
	}
      else
	{
	  // Remain quiet about no wavelength options
	}

      nTemp = oHeader.nGetValue(sPrefix + ms_sScanCrysDatum, &m_nNumCrysDatum);
      if ( (0 == nTemp) && (0 < m_nNumCrysDatum) )
	{
	  nTemp = oHeader.nGetValue(sPrefix + ms_sScanCrysDatum,
				    m_nNumCrysDatum+1, fValues);
	  if (0 == nTemp)
	    {
	      if (NULL != m_pfCrysDatum) // This should always be false
		{                        // because of call to ::nInitValues()
		  delete [] m_pfCrysDatum;
		}
	      m_pfCrysDatum = new float [m_nNumCrysDatum];
	      for (i = 0; i < m_nNumCrysDatum; i++)
		{
		  m_pfCrysDatum[i] = fValues[i+1];
		}
	    }
	}
      else
	{
	  // Remain quiet about no crystal datum, but set them to 0.0
	}

      nTemp = oHeader.nGetValue(sPrefix + ms_sScanDetDatum, &m_nNumDetDatum);
      if ( (0 == nTemp) && (m_nNumDetDatum > 0) )
	{
	  nTemp = oHeader.nGetValue(sPrefix + ms_sScanDetDatum,
				    m_nNumDetDatum+1, fValues);
	  if (0 == nTemp)
	    {
	      if (m_pfDetDatum != NULL) // This should always be false
		{                       // because of call to ::nInitValues()
		  delete [] m_pfDetDatum;
		}
	      m_pfDetDatum = new float [m_nNumDetDatum];
	      for (i = 0; i < m_nNumDetDatum; i++)
		{
		  m_pfDetDatum[i] = fValues[i+1];
		}
	    }
	}
      else
	{
	  // Remain quiet about no detector datum, but set them to 0.0
	}

      nTemp = oHeader.nGetValue(sPrefix + ms_sScanMode, &sTemp);
      if (0 == nTemp)
	{
	  if (sTemp == ms_sScanModeStillO)
	    {
	      m_eThe_Mode = eScanMode_StillOpen;
	    }
	  else if (sTemp == ms_sScanModeStillC)
	    {
	      m_eThe_Mode = eScanMode_StillClosed;
	    }
	  else if (sTemp == ms_sScanModeScanO)
	    {
	      m_eThe_Mode = eScanMode_ScanOpen;
	    }
	  else if (sTemp == ms_sScanModeScanC)
	    {
	      m_eThe_Mode = eScanMode_ScanClosed;
	    }
	}
      else
	{
	  m_eThe_Mode = eScanMode_Unknown;
	}

/*
 * Detector specific options.
 * Maybe these should be an option string?
    nTemp = oHeader.nGetValue(sPrefix + ms_sScanDezinger, &bValue);
    if ( (0 == nTemp) )
       m_bDezingerImage = bValue;
    nTemp = oHeader.nGetValue(sPrefix + ms_sScanDetBinMode, &sTemp);
    if(0 == nTemp){
       if("2x2" == sTemp.substr(0,3))
          m_nDetBinMode = 2;
       else
          m_nDetBinMode = 1;
    }
 */
      nTemp = oHeader.nGetValue(sPrefix + ms_sScanDetectorOptions,&sTemp);
      if (0 == nTemp)
	m_sDetectorOptions = sTemp;

    }

  if (0 == nStat)
    {
      m_eThe_State = eScan_available_state; // Though m_poSpatial not available
    }                                     //    nor m_poNonunf
}

Cscan::~Cscan()
{
  if (NULL != m_poRotation)
    {
      delete m_poRotation;
    }
  if (NULL != m_pfCrysDatum)
    {
      delete [] m_pfCrysDatum;
      m_pfCrysDatum = NULL;
      m_nNumCrysDatum = 0;
    }
  if (NULL != m_pfDetDatum)
    {
      delete [] m_pfDetDatum;
      m_pfDetDatum = NULL;
      m_nNumDetDatum = 0;
    }
  if (m_bNewNonunf && (NULL != m_poNonunf) )
    {
      delete m_poNonunf;
      m_poNonunf = NULL;
      m_bNewNonunf = FALSE;
    }
}

int Cscan::nInitValues(void)
{
  m_eThe_State              = eScan_unknown_state;
  m_sScan_key               = "unknown";
  m_sFileTemplate           = "";
  m_nFileStartSeqNum        = 0;
  m_nFileStepNum            = 1;
  m_nFileSeqNum             = m_nFileStartSeqNum;
  m_nNumImages              = 0;
  m_nNumImagesAlloc         = 0;
  m_nNumImagesAvail         = 0;
  m_poTheImages             = NULL;
  m_nNumDarkImagesAlloc     = 0;
  m_nNumDarkImagesAvail     = 0;
  m_poTheDarkImages         = NULL;
  m_nNumDCoffsetImagesAlloc = 0;
  m_nNumDCoffsetImagesAvail = 0;
  m_poTheDCoffsetImages     = NULL;
  m_poRotation              = new Crotation (); // Default rotation
  m_poSpatial               = NULL;
  m_poNonunf                = NULL;
  m_bNewNonunf              = FALSE;
  m_nNumDetectors           = 0;
  m_nNumCounters            = 0;
  m_nNumDetDatum            = 0;
  m_pfDetDatum              = NULL;
  m_nNumCrysDatum           = 0;
  m_pfCrysDatum             = NULL;
  m_fWavelength             = 0.0;   // Unknown wavelength!
  m_nWavelengthOption       = 0;
  m_nWavelengthOptimize     = 0;
  //m_bDezingerImage          = FALSE;
  //m_nDetBinMode             = 1;
  m_sDetectorOptions        = "None";

  return (0);
}

int Cscan::nList(void)
{
  cout << "Scan key:  " << m_sScan_key << endl;
  if (NULL != m_poRotation)
    (void) m_poRotation->nList();
  else
    cout << "Rotation value unknown.\n";
  cout << "Template: " << m_sFileTemplate
       << "\n   Start: " << m_nFileStartSeqNum
       << "\n    Step: " << m_nFileStepNum << endl;
  return (0);
}

int Cscan::nSetRotation(const Crotation& oRotation)
{
  if (NULL != m_poRotation)
    {
      delete m_poRotation;  // Delete current Crotation object and
    }
  m_poRotation = new Crotation(oRotation);
  return (0);
}

int Cscan::nSetNonunf(const Cnonunf& oNonunf)
{
  if (NULL != m_poNonunf)
    {
      delete m_poNonunf;  // Delete current Cnonunf object and
    }
  m_poNonunf = new Cnonunf(oNonunf);   // Hmmm, this constructor does not exist!
  m_bNewNonunf = TRUE;
  if (!m_poNonunf->bIsAvailable())
    {
      cout << "Cscan::nSetNonunf error!" << endl;
      return (1);
    }
  return (0);
}

int Cscan::nGetImage(Cimage* poImage)
{
  int nStat;
  Cstring sName;
  nStat  = nGetImageName(&sName);
  if (0 == nStat)
    {
      if (NULL == poImage)
	{
// This does not work since poImage is deleted when it goes out of scope
//	  poImage = new Cimage (sName);
// Instead return error
	  return (-1);
	}
      else
	{
	  return (poImage->nRead(sName));
	}
    }
  else
    return (nStat);
}
Cstring 
Cscan::sGetImageName(const int nSeqNum)
{
  int nCurrSeqNum;
  nCurrSeqNum = m_nFileSeqNum; // Save
  if( nSeqNum != -999999 )
    {
      vSetSeqNum( nSeqNum );
    }
  Cstring sName;
  if (0 != nGetImageName( &sName ))
    {
      sName = "INVALID_IMAGE.001";
    }
  m_nFileSeqNum = nCurrSeqNum;  // Restore
  return (sName);
}


int Cscan::nGetImageName(Cstring* psName)
{

  // Create a filename from the m_sFileTemplate and m_nFileSeqNum
  // Perhaps we need a ::sGetImageName(const int nSeqNum);

  int  i, nLen, nFirstHash;
  int  nStat;
  bool bNegative;

  char *cTemp;
  cTemp = new char [m_sFileTemplate.length()+10];  // This is max length we'll need

  *psName = m_sFileTemplate;

  bNegative = (0 > m_nFileSeqNum);        // Take care of negative sign later

  // Count number of available places 
  //   (should be a member variable and done when m_sFileTemplate assigned)

  int nPlaces = 0;
  for (i = 0; i < psName->length(); i++)
    {
      if (   ('?' == psName->GetAt(i))
	  || ('#' == psName->GetAt(i)))
	nPlaces++;
    }

  nLen = sprintf(cTemp, "%d", abs(m_nFileSeqNum));

  if (0 >= nLen)
    {
      cout << "ERROR: Cscan: writing m_nFileSeqNum: " << m_nFileSeqNum << endl;
      nStat = nLen;
    }
  else
    {
// mrp - eliminate troublesome code that was apparently designed to handle old Siemens style image numbers, eg. *.a00
// mrp This code causes bug 5851
// mrp This format was discontinued in 1999, so supporting it is not worth the trouble
//      if (nLen == (nPlaces+1))
//	{
//	  // If nLen > nPlaces, then modify cTemp so first digit is 10=a, 11=b

//	  int j;
//          char acTwo[3];
//	  acTwo[0] = cTemp[0];
//	  acTwo[1] = cTemp[1];
//	  acTwo[2] = '\0';
//	  j = atoi(acTwo) - 10;
//	  cTemp[0] = char (j + 'a');
//	  for (i = 1; i < nLen; i++)  // This copies the '\0', too.
//	    {
//	      cTemp[i] = cTemp[i+1];
//	    }
//	  nLen--;
//	}

      nFirstHash = -1;
      for (i = psName->length()-1; i >= 0; i--)
	{
	  // Work backwards through psName
	  if (   ('?' == psName->GetAt(i))
	      || ('#' == psName->GetAt(i)))
	    {
	      nFirstHash = i;
	      if (0 < nLen)
		psName->SetAt(i, cTemp[--nLen]);
	      else if (0 == nLen)                 // No more digits in cTemp
		psName->SetAt(i, '0');            //   so insert leading zeroes
	    }
	}
      if (bNegative)
	{
	  if ( (0 <= nFirstHash) && ('0' == psName->GetAt(nFirstHash)) )
	    psName->SetAt(nFirstHash, '-'); // Have negative sign and have space
	                                    // for it in psName
	  else
	    nLen = 1;
	}
      if ( (1 <= nLen) || (-1 == nFirstHash) )
	{
	  // All the ? marks were used up or there were no ? marks
//mrp	  cout << "ERROR: Cscan not enough ?'s in m_sFileTemplate!\n";
	  cout << "WARNING: Cscan not enough ?'s in m_sFileTemplate!\n"; //mrp
//mrp	  nStat = 1;
	  nStat = 0; //mrp
	}
      else
	{
	  nStat = 0;
	}
    }

//  cout << "TEST: in: " << m_nFileSeqNum << " out: "
//       << nGetSeqNum(*psName) << endl << flush;

  delete [] cTemp;
  cTemp = NULL;

  return (nStat);
}

int
Cscan::nOverlay(Cimage *poImage, const int nNum,
		void (*prvProgressCallbackProc)(void *pObj,  int *pnImg),
		void *pObj)
{
  // Overlay a sequence of images: find the highest per pixel value in a scan
  // Return 0 - success
  //       >1 - error occurred at image

  int nStat;
  int i, j;

  nStat = 1;
  vInitSeqNum();                        // Initialize to first sequence
  nGetImage(poImage);                   // Get the first image
  if (!poImage->bIsAvailable())         // If error getting image, return
    return (nStat);

  // Get scan and rotation info of the modified image before it is modified
  
  Cscan     *poScanTmp = new Cscan(poImage->m_oHeader);
  Crotation *poRotnTmp = new Crotation(poImage->m_oHeader);

  //Debug poScanTmp->nList();
  float fIncrement = 0.0;
  fIncrement = poRotnTmp->fGetIncrement();

  Cimage oImage;
  int nPixels;
  nPixels = poImage->nGetDimension();   // Get dimensions
  unsigned short uiIn, uiOut;
  for (i = 1; i < nNum; i++)
    {
      if (NULL != prvProgressCallbackProc)
	{
	  j = i;
	  (*prvProgressCallbackProc)(pObj, &j);
	  if (j != i)
	    {
	      return (-1);
	    }
	}
      vNextSeqNum();
      nGetImage(&oImage);               // Get next image
      nStat++;
      if (!oImage.bIsAvailable())
	break;                          // Break on error getting image

      if (2 == oImage.nBytesPerPixel())
	{
	  oImage.nSetNextPixel(0, 0);       // Start at first pixel in image
	  poImage->nSetNextPixel(0, 0);
	  for (j = 0; j < nPixels; j++)     // Loop over all pixels
	    {
	      uiIn = oImage.uiGetNextPixel();        // Get next pixel in images
	      uiOut = poImage->uiGetNextPixelNoInc();
	      if (uiIn < uiOut)                      // Compare
		{
		  uiIn = uiOut;
		}
	      poImage->vSetNextPixel(uiIn);          // Keep largest
	    }
	}
      else
	{
	  int i;
	  int nDim0, nDim1;
	  float fPixIn, fPixOut;
	  (void) poImage->nGetDimensions(&nDim0, &nDim1);
	  for (j = 0; j < nDim1; j++)     // Loop over all pixels
	    {
	      for (i = 0; i < nDim0; i++)     // Loop over all pixels
		{
		  fPixIn  = (oImage.*oImage.prfGetPixel)(i,j);
		  fPixOut = (poImage->*poImage->prfGetPixel)(i,j);
		  if (fPixIn < fPixOut)                      // Compare
		    {
		      fPixIn = fPixOut;
		    }
		  (poImage->*poImage->prnSetPixel)(i, j, fPixIn); // Keep largest
		}
	    }
	}

      // Change scan info in the modified image as it has been modified already

      poRotnTmp->vSetIncrement( nNum * fIncrement);
      poRotnTmp->vSetRotEnd(poRotnTmp->fGetRotEnd() + fIncrement);
      poRotnTmp->nUpdateHeader((Cimage_header*)&poImage->m_oHeader);

      poScanTmp->m_poRotation->vSetRotEnd(poRotnTmp->fGetRotEnd());
      poScanTmp->m_poRotation->vSetIncrement( nNum * fIncrement);

      poScanTmp->nUpdateHeader((Cimage_header*)&poImage->m_oHeader);
      //Debug poScanTmp->nList();
      
    }

  delete poScanTmp;
  delete poRotnTmp;

  if (nStat >= nNum)
    nStat = 0;
  return (nStat);
}

int
Cscan::nAdd(Cimage *poImage, const int nNum,
                void (*prvProgressCallbackProc)(void *pObj,  int *pnImg),
                void *pObj)
{
  // Add a sequence of images to make a new image, watch out for bad pixels 
  //   and for saturated pixels.
  // Return 0 - success
  //       >1 - error occurred at image

  int nStat;
  int i, j;

  nStat = 1;
  vInitSeqNum();                      // Initialize to first sequence
  nGetImage(poImage);                 // Get the first image, this is modified
  if (!poImage->bIsAvailable())       // If error getting image, return
    return (nStat);

  // Get scan and rotation info of the modified image before it is modified
  
  Cscan     *poScanTmp = new Cscan(poImage->m_oHeader);
  Crotation *poRotnTmp = new Crotation(poImage->m_oHeader);

  //Debug poScanTmp->nList();
  float fIncrement = 0.0;
  fIncrement = poRotnTmp->fGetIncrement();

  Cimage oImage;
  int nPixels;
  nPixels = poImage->nGetDimension();   // Get dimensions
  unsigned short uiIn, uiOut;
  for (i = 1; i < nNum; i++)
    {
      if (NULL != prvProgressCallbackProc)
        {
          j = i;
          (*prvProgressCallbackProc)(pObj, &j);
          if (j != i)
            {
              return (-1);
            }
        }
      vNextSeqNum();
      nGetImage(&oImage);               // Get next image
      nStat++;
      if (!oImage.bIsAvailable())
        break;                          // Break on error getting image

      // Always use this method when adding images because real math is used and it is safer

      int i;
      int nDim0, nDim1;
      float fPixIn, fPixOut;
      float fSatVal;
      fSatVal = oImage.fGetSatValue();
      (void) poImage->nGetDimensions(&nDim0, &nDim1);
      for (j = 0; j < nDim1; j++)     // Loop over all pixels
        {
          for (i = 0; i < nDim0; i++)     // Loop over all pixels
            {
              fPixIn  = (oImage.*oImage.prfGetPixel)(i,j);
              fPixOut = (poImage->*poImage->prfGetPixel)(i,j);

              if (fPixIn >= 0)
                fPixIn = fPixOut + fPixIn;  // Add the pixel values if input greater than 0
              if (fPixIn > fSatVal) fPixIn = fSatVal;

              (poImage->*poImage->prnSetPixel)(i, j, fPixIn); // Put sum in the result image
            }
        }

      // Change scan info in the modified image as it has been modified already

      poRotnTmp->vSetIncrement( nNum * fIncrement);
      poRotnTmp->vSetRotEnd(poRotnTmp->fGetRotEnd() + fIncrement);
      poRotnTmp->nUpdateHeader((Cimage_header*)&poImage->m_oHeader);

      poScanTmp->m_poRotation->vSetRotEnd(poRotnTmp->fGetRotEnd());
      poScanTmp->m_poRotation->vSetIncrement( nNum * fIncrement);

      poScanTmp->nUpdateHeader((Cimage_header*)&poImage->m_oHeader);
      //Debug poScanTmp->nList();
      
    } // Next image

  delete poScanTmp;
  delete poRotnTmp;

  if (nStat >= nNum)
    nStat = 0;
  return (nStat);
}

int
Cscan::nTile(Cimage* poImage, const int nNum,
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg),
	     void *pObj)
{
  // Tile a sequence of images
  // Return 0 - success
  //       >1 - error occurred at image

  int nStat;
  int i, j, k, l, m;   // Loop counters
  int nRows, nCols;
  int nRowOff, nColOff;
  int n0Off, n1Off;
  int nDim0, nDim1;
  float fTemp;
  float fStep0, fStep1;

  nStat = 1;
  if (1 > nNum) return (nStat);         // No images to tile

  vInitSeqNum();                        // Initialize to first sequence
  nGetImage(poImage);                   // Get the first image
  if (!poImage->bIsAvailable())         // If error getting image, return
    {
      return (nStat);
    }

  fTemp  = sqrt((double)nNum);
  nCols  = (int) (fTemp + 0.5);
  nRows  = nNum / nCols;
  if ( (nRows * nCols) < nNum) nRows++;
  (void) poImage->nGetDimensions(&nDim0, &nDim1);   // Get dimensions

  fStep0 = (float) nCols;
  fStep1 = max(fStep0,(float) nRows);
  fStep0 = fStep1;                        // Keep aspect ratio the same

  // Shrink first image into its corner

  for (i = 0; i < (nDim1 / nRows); i++)
    {
      k = (int) ((float) i * fStep1);
      for (j = 0; j < (nDim0 / nCols); j++)
	{
	  l = (int) ((float) j * fStep0);
	  fTemp = (poImage->*poImage->prfGetPixel)(l, k);
	  //fTemp = poImage->fGetPixel(l, k);
	  (poImage->*poImage->prnSetPixel)(j, i, fTemp);
	  //poImage->nSetPixel(j, i, fTemp);
	}
    }

  // Now work with the subsequent images

  nRowOff = 0;
  nColOff = 0;

  Cimage oImage;
  for (m = 1; m < nNum; m++)
    {
      if (NULL != prvProgressCallbackProc)
	{
	  i = m;
	  (*prvProgressCallbackProc)(pObj, &i);
	  if (i != m)
	    {
	      return (-1);
	    }
	}
      vNextSeqNum();
      nGetImage(&oImage);               // Get next image
      nStat++;
      if (!oImage.bIsAvailable())
	break;                          // Break on error getting image

      nColOff++;
      if (nColOff >= nCols)
	{
	  nColOff = 0;
	  nRowOff++;
	}

      n0Off = nColOff * (nDim0 / nCols);
      n1Off = nRowOff * (nDim1 / nRows);
      for (i = 0; i < (nDim1 / nRows); i++)
	{
	  k = (int) ((float) i * fStep1);
	  for (j = 0; j < (nDim0 / nCols); j++)
	    {
	      l = (int) ((float) j * fStep0);
	      fTemp = (oImage.*oImage.prfGetPixel)(l, k);
	      //fTemp = oImage.fGetPixel(l, k);
	      (poImage->*poImage->prnSetPixel)(n0Off+j, n1Off+i, fTemp);
	      //poImage->nSetPixel(n0Off+j, n1Off+i, fTemp);
	    }
	}
    }

  if (nStat >= nNum)
    {
      nStat = 0;

      // Indicate in header that image is tiled

      int a2nTemp[2];
      a2nTemp[0] = nRows;
      a2nTemp[1] = nCols;
      poImage->m_oHeader.nReplaceValue(D_K_DtdisplayTile, 2, a2nTemp);
    }

  return (nStat);
}

int
Cscan::nUnderlay(Cimage* poImage, const int nNum,
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg),
	     void *pObj)
{
  // Underlay a sequence of images: find the lowest per pixel value in a scan
  // Return 0 - success
  //       >1 - error occurred at image

  int nStat;
  int i, j;

  nStat = 1;
  vInitSeqNum();                        // Initialize to first sequence
  nGetImage(poImage);                   // Get the first image
  if (!poImage->bIsAvailable())         // If error getting image, return
    return (nStat);
  Cimage oImage;
  int nPixels;
  nPixels = poImage->nGetDimension();   // Get dimensions
  unsigned short uiIn, uiOut;
  for (i = 1; i < nNum; i++)
    {
      if (NULL != prvProgressCallbackProc)
	{
	  j = i;
	  (*prvProgressCallbackProc)(pObj, &j);
	  if (j != i)
	    {
	      return (-1);
	    }
	}
      vNextSeqNum();
      nGetImage(&oImage);               // Get next image
      nStat++;
      if (!oImage.bIsAvailable())
	break;                          // Break on error getting image
      if (2 == oImage.nBytesPerPixel())
	{
	  oImage.nSetNextPixel(0, 0);       // Start at first pixel in image
	  poImage->nSetNextPixel(0, 0);
	  for (j = 0; j < nPixels; j++)     // Loop over all pixels
	    {
	      uiIn = oImage.uiGetNextPixel();        // Get next pixel in images
	      uiOut = poImage->uiGetNextPixelNoInc();
	      if ( ((uiIn > uiOut) || (0 == uiIn)) && (0 != uiOut) )    // Compare
		{
		  uiIn = uiOut;
		}
	      poImage->vSetNextPixel(uiIn);          // Keep smallest
	    }
	}
      else
	{
	  int i;
	  int nDim0, nDim1;
	  float fPixIn, fPixOut;
	  (void) poImage->nGetDimensions(&nDim0, &nDim1);
	  for (j = 0; j < nDim1; j++)     // Loop over all pixels
	    {
	      for (i = 0; i < nDim0; i++)     // Loop over all pixels
		{
		  fPixIn  = (oImage.*oImage.prfGetPixel)(i,j);
		  fPixOut = (poImage->*poImage->prfGetPixel)(i,j);
		  //if ( ((fPixIn > fPixOut) || (-1 == fPixIn)) && (-1 != fPixOut) )    // Compare
		  if ( ((fPixIn > fPixOut) || (0 == fPixIn)) && (0 != fPixOut) )    // Compare
		    {
		      fPixIn = fPixOut;
		    }
		  (poImage->*poImage->prnSetPixel)(i, j, fPixIn); // Keep smallest
		}
	    }
	}
    }
  if (nStat >= nNum)
    nStat = 0;
  return (nStat);
}

int
Cscan::nSetDatum(const int nCrys, float *pfCrys, const int nDet, float *pfDet)
{
  int i;
  int nStat = 0;
  //  cout << "Cscan:nSetDatum, nCrys, nDet: " << nCrys << ", " << nDet << endl << flush;
  if (nCrys > 0)
    {
      if (NULL == m_pfCrysDatum)
	{
	  m_pfCrysDatum = new float [nCrys];
	}
      else if (nCrys > m_nNumCrysDatum)
	{
	  delete [] m_pfCrysDatum;
	  m_pfCrysDatum = new float [nCrys];
	}
      m_nNumCrysDatum = nCrys;
      for (i = 0; i < m_nNumCrysDatum; i++)
	{
	  m_pfCrysDatum[i] = pfCrys[i];
	}
    }
  else
    {
      nStat = 1;
    }

  if (nDet > 0)
    {
      if (NULL == m_pfDetDatum)
	{
	  m_pfDetDatum = new float [nDet];
	}
      else if (nDet > m_nNumDetDatum)
	{
	  delete [] m_pfDetDatum;
	  m_pfDetDatum = new float [nDet];
	}
      m_nNumDetDatum = nDet;
      for (i = 0; i < m_nNumDetDatum; i++)
	{
	  m_pfDetDatum[i] = pfDet[i];
	}
    }
  else
    {
      nStat = nStat + 1;
    }
  return (nStat);
}

int
Cscan::nGetDatum(int *pnCrys, float *pfCrys, int *pnDet, float *pfDet)
{
  // Get the datum values in the scan object and return them in the
  // passed arguments.
  // On entry, *pnCrys has the dimensions of the *pfCrys array
  //           *pnDet  has the dimensions of the *pfDet  array
  // On exit   *pnCrys has the number of datum values transferred to *pfCrys
  //           *pnDet  has the number of datum values transferred to *pfDet
  // If there was enough space in the input arguments, return 0
  // otherwise return 1 if not enough dimensions in *pfCrys
  //                  2 if not enough dimensions in *pfDet
  //                  3 if not enough dimensions in both

  int i;
  int nStat = 0;
  int nDim;
  nDim = m_nNumCrysDatum;
  //  cout << "Cscan::nGetDatum " << m_nNumCrysDatum << ", " << m_nNumDetDatum
  //  << endl << flush;
  if (*pnCrys < nDim)
    {
      nStat = 1; // Not enough space
      nDim = *pnCrys;
    }
  *pnCrys = nDim;
  for (i = 0; i < nDim; i++)
    {
      pfCrys[i] = m_pfCrysDatum[i];
    }

  nDim = m_nNumDetDatum;
  if (*pnDet < nDim)
    {
      nStat = nStat + 2;
      nDim = *pnDet;
    }
  *pnDet = nDim;
  for (i = 0; i < m_nNumDetDatum; i++)
    {
      pfDet[i] = m_pfDetDatum[i];
    }
  return (nStat);
}

int
Cscan::nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre)
{
  // Place value of scan into the image header

  //bool bTemp;
  int     i;
  int     a3nTemp[3];
  float  *pfTemp;
  Cstring sTemp;

  // Add path to scan template if it has no path

  sTemp = sFileGetDirectory(m_sFileTemplate)   // Return cwd if no directory
          + sFileGetBasename(m_sFileTemplate);
  (void) poHeader->nReplaceValue(sPre + ms_sScanTemplate, sTemp);

  a3nTemp[0] = m_nFileStartSeqNum;
  a3nTemp[1] = m_nFileStepNum;
  a3nTemp[2] = m_nNumImages;
//  a3nTemp[2] = m_nFileSeqNum;

  (void) poHeader->nReplaceValue(sPre + ms_sScanSeqInfo, 3, a3nTemp);
  (void) poHeader->nReplaceValue(sPre + ms_sScanWavelength, m_fWavelength, 6);
  a3nTemp[0] = m_nWavelengthOption;
  a3nTemp[1] = m_nWavelengthOptimize;
  (void) poHeader->nReplaceValue(sPre + ms_sScanWavelengthOpts, 2, a3nTemp);

  pfTemp = new float [m_nNumCrysDatum+1];
  pfTemp[0] = (float) m_nNumCrysDatum;
  for (i = 0; i < m_nNumCrysDatum; i++)
    {
      pfTemp[i+1] = m_pfCrysDatum[i];
    }
  (void) poHeader->nReplaceValue(sPre + ms_sScanCrysDatum,
				 m_nNumCrysDatum+1, pfTemp, 3);
  delete [] pfTemp;

  pfTemp = new float [m_nNumDetDatum+1];
  pfTemp[0] = (float) m_nNumDetDatum;
  for (i = 0; i < m_nNumDetDatum; i++)
    {
      pfTemp[i+1] = m_pfDetDatum[i];
    }
  (void) poHeader->nReplaceValue(sPre + ms_sScanDetDatum,
				 m_nNumDetDatum+1, pfTemp, 3);
  delete [] pfTemp;


  if (m_eThe_Mode == eScanMode_StillClosed)
    {
      sTemp = ms_sScanModeStillC;
    }
  else if (m_eThe_Mode == eScanMode_StillOpen)
    {
      sTemp = ms_sScanModeStillO;
    }
  else if (m_eThe_Mode == eScanMode_ScanClosed)
    {
      sTemp = ms_sScanModeScanC;
    }
  else if (m_eThe_Mode == eScanMode_ScanOpen)
    {
      sTemp = ms_sScanModeScanO;
    }
  else
    {
      sTemp = ms_sScanModeUnknown;
    }

  (void) poHeader->nReplaceValue(sPre + ms_sScanMode, sTemp);
				
  (void) m_poRotation->nUpdateHeader(poHeader, sPre + ms_sScanPrefix);

  // These are no longer supported, replaced by the detector option string.
  //bTemp = m_bDezingerImage;
  //(void) poHeader->nReplaceValue(sPre + ms_sScanDezinger,bTemp);
  //sTemp = (2 == m_nDetBinMode) ? "2x2" : "1x1";
  //(void) poHeader->nReplaceValue(sPre + ms_sScanDetBinMode, sTemp);
  // Since they are no longer supported, remove them from the header if they
  // are there.
  (void) poHeader->nDelete(sPre + ms_sScanDezinger);
  (void) poHeader->nDelete(sPre + ms_sScanDetBinMode);

  sTemp = m_sDetectorOptions;
  (void) poHeader->nReplaceValue(sPre + ms_sScanDetectorOptions,sTemp);

  return (0);
}


void
Cscan::vDeleteFromHeader(Cimage_header *poHeader, const Cstring &sPre)
{
  // Delete this scan from the header

  (void) poHeader->nDelete(sPre + ms_sScanTemplate);
  (void) poHeader->nDelete(sPre + ms_sScanSeqInfo);
  (void) poHeader->nDelete(sPre + ms_sScanWavelength);
  (void) poHeader->nDelete(sPre + ms_sScanCrysDatum);
  (void) poHeader->nDelete(sPre + ms_sScanDetDatum);
  (void) poHeader->nDelete(sPre + ms_sScanMode);
  (void) poHeader->nDelete(sPre + ms_sScanDezinger);
  (void) poHeader->nDelete(sPre + ms_sScanDetBinMode);
  (void) poHeader->nDelete(sPre + ms_sScanDetectorOptions);
  m_poRotation->vDeleteFromHeader(poHeader, sPre + ms_sScanPrefix);
}

int
Cscan::nGetNumImages(void)
{
  if (0.0 != m_poRotation->fGetIncrement())
    {
      return ((int)(  (m_poRotation->fGetRotEnd() - m_poRotation->fGetRotStart())
		    / m_poRotation->fGetIncrement()  + 0.5)); // Round-up
    }
  else
    {
      if (0 > m_nNumImages)
	return (1);
      return (m_nNumImages);
    }
}

int
Cscan::nAvgSD(Cimage* poImage, const int nNum,
	     void (*prvProgressCallbackProc)(void *pObj, int *pnImg),
	     void *pObj)
{
  // Compute average and standard deviation of a series of images.
  // Return 0 - success
  //       >1 - error occurred at image

  int   nStat;
  int   i, j, k;
  int   nDim0, nDim1;
  float fPixel;

  nStat = 1;
  vInitSeqNum();                        // Initialize to first sequence
  nGetImage(poImage);                   // Get the first image
  if (!poImage->bIsAvailable())         // If error getting image, return
    return (nStat);

  (void) poImage->nGetDimensions(&nDim0, &nDim1);

  // Two float images: one for sum of pixels for average
  //                   one for sum of squares for sd

  Cimage oSumImage(nDim0, nDim1, eImage_realIEEE);
  Cimage oSum2Image(nDim0, nDim1, eImage_realIEEE);

  // Initialize the sums (do it this way to avoid overwriting *poImage in case
  // there is an error down below)

  (void) oSumImage.nSetNextPixel(0, 0);
  (void) oSum2Image.nSetNextPixel(0, 0);
  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
	{
	  fPixel = (poImage->*poImage->prfGetPixel)(i,j);
	  oSumImage.vSetNextPixel(fPixel);
	  oSum2Image.vSetNextPixel(fPixel * fPixel);
	}
    }

  Cimage oImage;
  for (k = 1; k < nNum; k++)
    {
      if (NULL != prvProgressCallbackProc)
	{
	  j = k;
	  (*prvProgressCallbackProc)(pObj, &j);
	  if (j != k)
	    {
	      return (-1);
	    }
	}
      vNextSeqNum();
      nGetImage(&oImage);               // Get next image
      if (!oImage.bIsAvailable())
	break;                          // Break on error getting image
      nStat++;                          // Keep track of images read

      (void) oSumImage.nSetNextPixel(0, 0);  // Start at first pixel
      (void) oSum2Image.nSetNextPixel(0, 0); // Start at first pixel

      for (j = 0; j < nDim1; j++)       // Loop over slow dim
	{
	  for (i = 0; i < nDim0; i++)   // Loop over fast dim
	    {
	      // Get pixel in image

	      fPixel = (oImage.*oImage.prfGetPixel)(i,j);

	      // Add to sum and sum of squares

	      oSumImage.vSetNextPixel(oSumImage.fGetNextPixelNoInc() + fPixel);
	      oSum2Image.vSetNextPixel(oSum2Image.fGetNextPixelNoInc()
				       + fPixel * fPixel);
	    }
	}
    }

  // Compute average and sd on a pixel-by-pixel basis

  if (1 < nStat)
    {
      // Do this only if more than 1 image

      (void) oSumImage.nSetNextPixel(0, 0);
      (void) oSum2Image.nSetNextPixel(0, 0);

      float fNumPixels = (float) nStat;         // nStat holds actual num images
      float fNumNumPixM1  = fNumPixels * (fNumPixels - 1.0f);
      float fTemp1, fTemp2;

      for (j = 0; j < nDim1; j++)       // Loop over slow dim
	{
	  for (i = 0; i < nDim0; i++)   // Loop over fast dim
	    {
	      // For each pixel, compute sample standard deviation =
	      //
	      //        sqrt([sum (x**2) - (sum x)**2 / n] / (n - 1))
	      //  -or-
	      //        sqrt([n * sum (x**2) - (sum x)**2] / (n * (n - 1))
	      // Compute average =  (sum x) / n

	      fTemp1 =  oSumImage.fGetNextPixel();      // get and increment
	      fTemp2 = oSum2Image.fGetNextPixelNoInc(); // get and do not inc
	
	      fTemp2 = (fNumPixels * fTemp2 - (fTemp1 * fTemp1)) / fNumNumPixM1;
	      if (0.0 < fTemp2)
		{
		  fTemp2 = sqrt((double)fTemp2);
		}
	      else
		{
		  fTemp2 = sqrt((double)-fTemp2);
		}
	      fTemp1 = fTemp1 / fNumPixels;
//2013-05-10 +jwp
//	      poImage->nSetPixel(i, j, fTemp1);  // Original image has average
             (poImage->*poImage->prnSetPixel)(i, j, fTemp1);
//2013-05-10 -jwp
	      oSum2Image.vSetNextPixel(fTemp2);  // this one has std.dev
                                                 //  (now increment!)
	    }
	}
      Cstring sComment;
      char    a255cTemp[255];
      sprintf(a255cTemp, "Standard deviation of %d images with template: %s,"
	                 " start: %d and incr: %d",
	      nStat, m_sFileTemplate.string(),
	      nGetSeqNum(0), nGetSeqInc());
      sComment = a255cTemp;
      oSum2Image.m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sComment);
      nUpdateHeader(&oSum2Image.m_oHeader);  // Place scan info into header
      vInitSeqNum();                       // Reset to first image number
      (void) nGetImageName(&sComment);     //  so this will be first image name
      sComment = sComment + "sd";   // Add sd to image name

      // Strip off possible directory info for Unix, DOS and VMS systems

      while (0 <= sComment.index("/")) sComment = sComment.after("/"); // Unix
      while (0 <= sComment.index("]")) sComment = sComment.after("]"); // VMS
      while (0 <= sComment.index(":")) sComment = sComment.after(":"); // more VMS
      while (0 <= sComment.index("\\")) sComment = sComment.after("\\"); // DOS
      oSum2Image.nWrite(sComment);  // Write out standard deviation

      sprintf(a255cTemp, "Average of %d images with template: %s,"
	                 " start: %d and incr: %d",
	      nStat, m_sFileTemplate.string(),
	      nGetSeqNum(0), nGetSeqInc());
      sComment = a255cTemp;
      poImage->m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sComment);

    }
  if (nStat >= nNum)
    nStat = 0;
  return (nStat);
}

float
Cscan::fCalcGetRotStart(const int nSeqNum)
{
  // Return the rotation angle start for image with sequence number nSeqNum.

  float fTemp;

  fTemp = m_poRotation->fGetRotStart();
  if (0 != m_nFileStepNum)
    {
      fTemp = (float) ((nSeqNum - m_nFileStartSeqNum) / m_nFileStepNum)
	              * m_poRotation->fGetIncrement()  +  fTemp;
    }
  return (fTemp);
}

float
Cscan::fCalcGetRotEnd(const int nSeqNum)
{
  // Return the rotation angle end for image with sequence number nSeqNum.

  return (m_poRotation->fGetIncrement() + fCalcGetRotStart(nSeqNum));
}

int
Cscan::nGetSeqNum(const Cstring& rsImageName)
{
  // Get the sequence number of an image from its filename and the
  // current template.  Return -99999 if unable to determine.

  Cstring sScanBase;
  Cstring sFileBase;
  sFileBase = sFileGetBasename(rsImageName);
  sScanBase = sFileGetBasename(sGetTemplate());

  int nSeq     = 0;
  int nDecimal = 1;
  int i;
  char cFile;
  char cScan;
  for (i = sScanBase.length(); i >= 0; i--)
    {
      cFile = sFileBase.GetAt(i);
      cScan = sScanBase.GetAt(i);
      if (cFile != cScan)
	{
	  if ( ('?' != cScan) && ('#' != cScan) )
	    {
	      // Different from scan template,
	      // but not a legitimate wildcard char in the scan template

	      return (-99999);
	    }
	  if ('-' == cFile)
	    {
	      // Only 1 minus sign allowed, so we are done

	      nSeq = -nSeq;
	      return (nSeq);
	    }
	  else if ( ('0' <= cFile) && ('9' >= cFile) )
	    {
	      // A digit between 0 and 9, inclusive

	      nSeq      = nSeq + nDecimal * ((int) cFile - (int) '0');
	      nDecimal *= 10;
	    }
	  else if ( ('a' <= cFile) && ('z' >= cFile) )
	    {
	      // A letter between a and z, inclusive

	      nSeq      = nSeq + nDecimal * ((int) cFile - (int) 'a' + 10);
	      nDecimal *= 10;
	    }
	  else
	    {
	      return (-99999);
	    }
	}
    }
  return (nSeq);
}
