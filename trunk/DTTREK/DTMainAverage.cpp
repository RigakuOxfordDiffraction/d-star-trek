//
// Copyright (c) 1998-2006 Rigaku
//
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainAverage.cpp   Initial author: RB           12-Sep-2006
//
// Transferred from dtaverage.cc.
// dtaverage.cc     Initial author: J.W. Pflugrath           07-May-1996
// Average a series of images on a pixel by pixel basis.  Exclude
// outliers from the averaging process.  Outliers are specified by the
// naming of an average image, a standard deviation image and a sigma
// value.

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
 
#include "Dtrek.h"
#include "Cscan.h"

#include "DTMainAverage.h"

using std::cout;
using std::endl;
using std::flush;

int CDTMainAverage::nExecute(unsigned int argc,      // Command line argument count
                               char        *argv[])  // Pointers to command line args
{
    int          nStat;
    Cstring      sInTemplate  = "image.???";
    Cstring      sTemplate;
    Cstring      sOutFile;
    Cstring      sAvgFile;
    Cstring      sSDFile;
    Cstring      sOutSDFile;
    Cstring      sRejectFile;
    int          nSeqStart     = 1;
    int          nSeqIncr      = 1;
    int          nNumImages    = 1000;
    int          nNumSkipped   = 0;

    Cimage      *poImage1;
    Cimage      *poImage2;
    Cimage      *poImageAvg;
    Cimage      *poImageSD;
    Cscan       *poScanIn;

    int          i, j, k;  // Loop counters
    float        fSigma = 3.0;
    float        fMaxSD = 1000000.0;
    float        fMinSD = 1.0;
    float        fMax   = 1000000.0;
    float        fMin   = -1000000.0;
    float        fWarnLimit = 0.5;

    int          nRejects, nSumRejects;
    int          nRejects0SD;
    int          nTemp;
    float        fTemp;
    Cstring      sTemp;
    bool         bSkip = TRUE;
    bool         bDoUnderlay = FALSE;
    bool         bDoOverlay  = FALSE;
    bool         bDoAdd      = FALSE;
    char         cReject;

    float        fPixel, fAvgPixel, fSDPixel, fTest, fDelta;

    poImage1   = NULL;
    poImage2   = NULL;
    poImageAvg = NULL;
    poImageSD  = NULL;
    poScanIn   = NULL;

    vDtrekSetModuleName("dtaverage");
    vPrintCopyrightInfo();

    // Copy command line to output log
    cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

    // Parse command line arguments

    argc--; argv++;

    nStat = 1;

    for (i = 0; i < argc; i++) 
    {
      sTemp = (const char*) argv[i];
    //      cout << "Command line string: >>" << sTemp << "<<" << endl;
  
      if ("-scan" == sTemp) 
    {
      i++;
      if (i < argc)
	  {
	      sInTemplate = argv[i];
	      
          if( bCreateHeader(sInTemplate) )
          {
	          poScanIn = new Cscan(*m_poHeader);
	          if (poScanIn->bIsAvailable())
		        {
		          nSeqStart  = poScanIn->nGetSeqNum(0);
		          nSeqIncr   = poScanIn->nGetSeqInc();
		          nNumImages = poScanIn->nGetNumImages();
		        }
	          else
		        {
		          DTREK_ERROR(2, "ERROR: dtaverage - invalid scan in header!");
		        }
	          delete poScanIn;
	          poScanIn = NULL;
          }
	          
          delete m_poHeader;
	      m_poHeader = NULL;
	    }
      else
	    {
	      DTREK_ERROR(2, "ERROR: dtaverage - missing scan file name!");
	    }
    }
      else if ("-avg" == sTemp) 
    {
      i++;
      if (i < argc)
	    {
	      sAvgFile = argv[i];
	    }
      else
	    {
	      DTREK_ERROR(2, "ERROR: dtaverage - missing average image file name!");
	    }
    }
      else if ("-sd" == sTemp) 
    {
      i++;
      if (i < argc)
	    {
	      sSDFile = argv[i];
	    }
      else
	    {
	      DTREK_ERROR(2, "ERROR: dtaverage - missing sd image file name!");
	    }
    }
      else if ("-outsd" == sTemp) 
    {
      i++;
      if (i < argc)
	    {
	      sOutSDFile = argv[i];
	    }
      else
	    {
	      DTREK_ERROR(2, "ERROR: dtaverage - missing outsd image file name!");
	    }
    }
      else if ("-outrej" == sTemp) 
    {
      i++;
      if (i < argc)
	    {
	      sRejectFile = argv[i];
	    }
      else
	    {
	      DTREK_ERROR(2, "ERROR: dtaverage - missing outrej image file name!");
	    }
    }
      else if ("-start" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%d", &nTemp);
	      if (1 == nStat)
	    {
	      nSeqStart = nTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-inc" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%d", &nTemp);
	      if (1 == nStat)
	    {
	      nSeqIncr = nTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-num" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%d", &nTemp);
	      if (1 == nStat)
	    {
	      nNumImages = nTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-template" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      sInTemplate = argv[i];
    //	      cout << "sInTemplate is " << sInTemplate << endl;
	    }
      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
    }
      else if ("-out" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      sOutFile = argv[i];
	    }
      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
    }
      else if ( ("-add" == sTemp) || ("-sum" == sTemp) )
    {
      bDoOverlay  = FALSE;
      bDoUnderlay = FALSE;
      bDoAdd      = TRUE;
    }
      else if ("-overlay" == sTemp) 
    {
      bDoOverlay  = TRUE;
      bDoUnderlay = FALSE;
      bDoAdd      = FALSE;
    }
      else if ("-underlay" == sTemp) 
    {
      bDoOverlay  = FALSE;
      bDoUnderlay = TRUE;
      bDoAdd      = FALSE;
    }
      else if ("-sigma" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nStat)
	    {
	      fSigma = fTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-maxsd" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nStat)
	    {
	      fMaxSD = fTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-minsd" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nStat)
	    {
	      fMinSD = fTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-max" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nStat)
	    {
	      fMax = fTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-min" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nStat)
	    {
	      fMin = fTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-warnlimit" == sTemp) 
    {
      i++;
      if (i < argc) 
	    {
	      nStat = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nStat)
	    {
	      fWarnLimit = fTemp;
	    }
	      else
	    DTREK_ERROR(6, "ERROR: dtaverage - invalid argument!");
	    }
    }
      else if ("-noskip" == sTemp) 
    {
      bSkip = FALSE;
    }
      else if ("-help" == sTemp)
    {
      DTREK_ERROR(0, "");
    }
      else
    {
      DTREK_ERROR(0,"" );
    }
    }

    // Create default scan objects;

    poScanIn  = new Cscan();
    poScanIn->vSetTemplate(sInTemplate);
    poScanIn->m_poRotation->vSetIncrement(0.0);
    poScanIn->vSetSeqStart(nSeqStart);
    poScanIn->vSetSeqInc(nSeqIncr);
    poScanIn->vSetNumImgs(nNumImages);

    ////////////////////////////////////////////////////////////////////////
    // Done with command line arguments so now go do the deed
    ////////////////////////////////////////////////////////////////////////

    poScanIn->vInitSeqNum();

    // Get first input image name and read in image

    nStat    = poScanIn->nGetImageName(&sInTemplate);

    if (0 != nStat)
    DTREK_ERROR(0, "First image name invalid!");

    // Set output file names if not set already
    // Strip off directory from sInTemplate for potential use in output files.

    sTemplate = sInTemplate;
    while (0 <= sTemplate.find("/")) sTemplate = sTemplate.after("/");

    if ("" == sSDFile)
    {
      sSDFile = sTemplate + "sd";
    }
    if ("" == sAvgFile)
    {
      sAvgFile = sTemplate + "av";
    }
    if ("" == sOutFile)
    {
      if (bDoUnderlay)
    sOutFile = sTemplate + "dtun";
      else if (bDoOverlay)
    sOutFile = sTemplate + "dtov";
      else if (bDoAdd)
    sOutFile = sTemplate + "dtad";
      else
    sOutFile = sTemplate + "dtav";
    }

    if ("" == sOutSDFile)
    {
      sOutSDFile = sTemplate + "dtsd";
    }

    if ("" == sRejectFile)
    {
      sRejectFile = sTemplate + "rejects";
    }

    cout << "\nTemplate:    " << poScanIn->sGetTemplate()
       << "\nSeq start:   " << nSeqStart
       << "\nSeq incr:    " << nSeqIncr
       << "\nNum imgs:    " << nNumImages;
    if (bDoUnderlay)
    {
      cout << "\nINFO: Performing underlay operation"
       << "\nOutput overlayed image file: " << sOutFile << "\n\n";
    }
    else if (bDoOverlay)
    {
      cout << "\nINFO: Performing overlay operation"
       << "\nOutput overlayed image file: " << sOutFile << "\n\n";
    }
    else if (bDoAdd)
    {
      cout << "\nINFO: Performing add operation"
       << "\nOutput added (summed) image file: " << sOutFile << "\n\n";
    }
    else
    {
      cout
       << "\nIn Avg img:  " << sAvgFile 
       << "\nIn SD  img:  " << sSDFile
       << "\nSigma:       " << fSigma
       << "\nMin:         " << fMin
       << "\nMax:         " << fMax
       << "\nMaxSD:       " << fMaxSD
       << "\nMinSD:       " << fMinSD
       << "\nOut rej img: " << sRejectFile
       << "\nOut SD img:  " << sOutSDFile
       << "\nOut Avg img: " << sOutFile << "\n\n";
      poImageAvg = new Cimage(sAvgFile);
      if (!poImageAvg->bIsAvailable())
    DTREK_ERROR(2, "ERROR: dtaverage - invalid input average image!");

      poImageSD = new Cimage(sSDFile);
      if (!poImageSD->bIsAvailable())
    DTREK_ERROR(2, "ERROR: dtaverage - invalid input sd image!");
    }

    poImage1 = new Cimage(sInTemplate);
    if (!poImage1->bIsAvailable())
    DTREK_ERROR(0, "First image in scan/series not available!");

    int nDim0, nDim1;

    (void) poImage1->nGetDimensions(&nDim0, &nDim1);

    // Now loop through images

    if (bDoUnderlay)
      {
        nStat = poScanIn->nUnderlay(poImage1, nNumImages, NULL, NULL);
        if (0 == nStat)
          {
            nStat =poImage1->nWrite(sOutFile);
            if (0 == nStat)
              cout << "INFO: underlay successful. Done.\n";
          }
        else
          {
            cout << "WARNING: underlay failed. Done.\n";
          }
        if (NULL != poImage1)
          delete poImage1;
        if (NULL != poScanIn)
          delete poScanIn;
        return (nStat);
      }
    else if (bDoOverlay)
      {
        nStat = poScanIn->nOverlay(poImage1, nNumImages, NULL, NULL);
        if (0 == nStat)
          {
            nStat =poImage1->nWrite(sOutFile);
            if (0 == nStat)
              cout << "INFO: overlay successful. Done.\n";
          }
        else
          {
            cout << "WARNING: overlay failed. Done.\n";
          }
        if (NULL != poImage1)
          delete poImage1;
        if (NULL != poScanIn)
          delete poScanIn;
        return (nStat);
      }

    else if (bDoAdd)
      {
        nStat = poScanIn->nAdd(poImage1, nNumImages, NULL, NULL);
        if (0 == nStat)
          {
//+JWP 2010-04-28
	    // Somehow figure out the SCAN_SEQ_INFO that needs to be updated.
	    // Get the sequence start for this starting image 

	    int nSeqNumNew = 0;
	    Cstring sTemplate;
	    int nPlaces = 4;

	    // Count number of ? in the template on the command line, but if less than 2, use 4 anyways
	    nPlaces   = poScanIn->nGetNumWildcards();
	    if (nPlaces < 2) nPlaces = 4;
	    sTemplate = sBuildScanTemplate(sOutFile, nPlaces, &nSeqNumNew);

	    // Set the number of the first image in the scan
	    // We may wish to change the other numbers of the scan sequence info later

	    Cscan *poScanTmp = new Cscan(poImage1->m_oHeader);
	    poScanTmp->vSetSeqStart(nSeqNumNew);

	    // Try to reduce the number of images in the scan object to what we think it might be

	    nNumImages = poScanTmp->nGetNumImages() / nNumImages;
	    poScanTmp->vSetNumImgs(nNumImages);

	    poScanTmp->nUpdateHeader((Cimage_header*)&poImage1->m_oHeader);
	    delete poScanTmp;
//-JWP 2010-04-28
            nStat =poImage1->nWrite(sOutFile);
            if (0 == nStat)
              cout << "INFO: overlay successful. Done.\n";
          }
        else
          {
            cout << "WARNING: overlay failed. Done.\n";
          }
        if (NULL != poImage1)
          delete poImage1;
        if (NULL != poScanIn)
          delete poScanIn;
        return (nStat);
      }

    // Stop here for underlay and overlay and add operations

    // Normal operation continue on and do averaging

    nStat = 1;

    // Two float images: one for sum of pixels for average
    //                   one for sum of squares for sd

    Cimage oSumImage(nDim0, nDim1, eImage_realIEEE);
    Cimage oSum2Image(nDim0, nDim1, eImage_realIEEE);

    // One integer image, for counting number of contributing pixels, set to 0s

    Cimage oCountImage(nDim0, nDim1, eImage_uI2);

    // Initialize the sums to 0's

    (void) oSumImage.nSetNextPixel(0, 0);
    (void) oSum2Image.nSetNextPixel(0, 0);
    (void) oCountImage.nSetNextPixel(0, 0);
    for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
    {
      oSumImage.vSetNextPixel((float)0.0);
      oSum2Image.vSetNextPixel((float)0.0);
      oCountImage.vSetNextPixel((unsigned short int)0);
    }
    }

    nSumRejects = 0;
    nStat       = 0;
    for (k = 1; (k <= nNumImages) && (0 == nStat); k++)
    {
      (void) oSumImage.nSetNextPixel(0, 0);  // Start at first pixel
      (void) oSum2Image.nSetNextPixel(0, 0); // Start at first pixel
      (void) oCountImage.nSetNextPixel(0, 0); // Start at first pixel

      nRejects    = 0;
      nRejects0SD = 0;

      printf("Image num: %d\nPixel         Value   Avg     SD    Test   Delta  R\n", 
	     k-1);
      for (j = 0; j < nDim1; j++)       // Loop over slow dim
    {
      for (i = 0; i < nDim0; i++)   // Loop over fast dim
	    {
	      // Get pixel values in the input images

	      fPixel    = (poImage1->*poImage1->prfGetPixel)(i,j);
	      fAvgPixel = (poImageAvg->*poImageAvg->prfGetPixel)(i,j);
	      fSDPixel  = (poImageSD->*poImageSD->prfGetPixel)(i,j);

	      // Test to see if pixel in poImage1 should be excluded...

	      if (fSDPixel > fMaxSD) fSDPixel = fMaxSD;
	      if (fSDPixel < fMinSD) fSDPixel = fMinSD;

	      fTest = (float)fabs((double)(fSigma * fSDPixel));
	      fDelta = (float)fabs((double)(fPixel - fAvgPixel));
	      if (1.0 > fTest) fTest = 1.0;
	      if ( (   (fDelta > fTest)
		    || (fPixel < fMin)
		    || (fPixel > fMax) )
	      &&
	      (0.0 != fDelta) )
	    {
	      // Fails tests, increment pixel pointers into images

	      (void) oSumImage.fGetNextPixel();
	      (void) oSum2Image.fGetNextPixel();
	      (void) oCountImage.uiGetNextPixel();		  

	      nRejects++;
	      if (1000 > nRejects)
		    {
		      // List rejected pixel info, but only for first 1000

		      if (fPixel > fMax)
		    cReject = '>';
		      else if (fPixel < fMin)
		    cReject = '<';
		      else if (fDelta > fTest)
		    cReject = 'D';
		      else
		    cReject = '?';

		      printf("%5d %5d %6.0f %6.0f %7.1f %7.1f %6.0f %c\n",
			     i, j, fPixel, fAvgPixel, fSDPixel,
			     fTest, fDelta, cReject);
		    }
	    }
	      else
	    {
	      // Passes sigma exclusion test, so sum it in
	      // Add to sum and sum of squares

	      oSumImage.vSetNextPixel(oSumImage.fGetNextPixelNoInc()
				      + fPixel);
	      oSum2Image.vSetNextPixel(oSum2Image.fGetNextPixelNoInc()
				       + fPixel * fPixel);
	      oCountImage.vSetNextPixel((unsigned short int)(1
			     +oCountImage.uiGetNextPixelNoInc()));
	    }
	    }
    }
      nSumRejects += nRejects;
      if (fWarnLimit <= (100.0 * (float)nRejects/(float(nDim0 * nDim1))) )
    {
      cout << "WARNING, rejected " << nRejects << " pixels, which"
               << " is more than the warning limit of " 
               << fWarnLimit << "% pixels in a single image.\n";
    }
      cout << nRejects << " pixels rejected.\n";

      poScanIn->vNextSeqNum();
      if (nNumImages > k)
    {
      // Get next image only if needed

      poScanIn->nGetImage(poImage1);      // Get next image
      while (bSkip && !poImage1->bIsAvailable() && (nNumImages > k))
	    {
	      // In case of ERROR reading an image, just skip it if there are more images

	      cout << "\nWARNING, image not available, so skipping it!\n" << endl;
	      k++;
	      nNumSkipped++;
	      poScanIn->vNextSeqNum();
	      poScanIn->nGetImage(poImage1);      // Get next image
	    }
    }

      if (!poImage1->bIsAvailable())
    {
      // Error never got an image, change nStat to be number of images read

      cout << "\nWARNING, image not available!\n" << endl;
      nStat = k;
    }
    }

    cout << nSumRejects << " pixels rejected from all processed images.\n";

    // Compute average and sd on a pixel-by-pixel basis

    if (0 == nStat)
    nStat = nNumImages - nNumSkipped;
    else
    nStat = nStat - nNumSkipped;

    if (2 > nStat)
    {
      cout << "Less than 2 input images: " << nStat << ", cannot average!\n";
    }
    else
    {
      cout << "\nNumber of input images used: " << nStat << "\n\n";

      // Do this only if more than 1 image

      (void) oSumImage.nSetNextPixel(0, 0);
      (void) oSum2Image.nSetNextPixel(0, 0);
      (void) oCountImage.nSetNextPixel(0, 0);

      float fNumPixels;         // nStat holds actual num images
      float fNumNumPixM1;
      float fTemp1, fTemp2;
    //      int nPixels = nDim0 * nDim1;

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

	      fTemp1       =  oSumImage.fGetNextPixel();      // get and inc
	      fTemp2       = oSum2Image.fGetNextPixelNoInc(); // get, do not inc
	      fNumPixels   = (float) oCountImage.uiGetNextPixelNoInc();
	      fNumNumPixM1 = fNumPixels * (fNumPixels - 1.0);	  

	      // Set the pixel value in the count image to the number of
	      // images with this pixel rejected

	      oCountImage.vSetNextPixel((unsigned short int) (nStat 
				    - (int) fNumPixels));

	      if (1.0 < fNumPixels)
	    {
	      fTemp2       = (fNumPixels * fTemp2 - (fTemp1 * fTemp1))
		                 / fNumNumPixM1;
	      if (0.0 < fTemp2)
		    {
		      fTemp2 = (float)sqrt((double)fTemp2);
		    }
	      else
		    {
		      fTemp2 = (float)sqrt((double)-fTemp2);
		    }
	      fTemp1 = fTemp1 / fNumPixels;
	      poImage1->nSetPixel(i, j, fTemp1);  // Orig. image has average
	      oSum2Image.vSetNextPixel(fTemp2);   // this one has std.dev
		                                      //  (now increment!)
	    }
	      else
	    {
	      oSum2Image.vSetNextPixel((float)0.0);
	    }
	    }
    }
      char    a255cTemp[255];
      sprintf(a255cTemp, "Standard deviation of %d images with template: %s"
	                 ", start: %d, incr: %d and sigma: %f",
	      nStat, (poScanIn->sGetTemplate()).string(), 
	      poScanIn->nGetSeqNum(0), poScanIn->nGetSeqInc(), fSigma);
      sTemp = a255cTemp;
      oSum2Image.m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sTemp);

      // Place scan info into header

      poScanIn->nUpdateHeader(&oSum2Image.m_oHeader);

      // Strip off possible directory info for Unix, DOS and VMS systems
  
      cout << "\n";
      oSum2Image.nWrite(sOutSDFile);  // Write out standard deviation
  
      sprintf(a255cTemp, "Average of %d images with template: %s,"
	                 " start: %d, incr: %d; sigma: %f",
	      nStat, poScanIn->sGetTemplate().string(), 
	      poScanIn->nGetSeqNum(0), poScanIn->nGetSeqInc(), fSigma);
      sTemp = a255cTemp;
      poImage1->m_oHeader.nReplaceValue(Cimage_header::ms_sComment, sTemp);
      poImage1->nWrite(sOutFile);

      // Write out image with rejected pixels marked

      oCountImage.nWrite(sRejectFile);
    }      

    /*
      // Update header of output image, use first image for most of header info
  
      oOutImage.oHeader = poImage1->oHeader;
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sDataType, 
                                      Cstring("short int"));
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sCompression,
                                      Cstring("None"));
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sDataType, 
                                      Cstring("short int"));
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sComment2, 
				      Cstring("dtaverage: ") + sInTemplate
                                      + " - " + sOutFile);
    */

    if (NULL != poImage1)
    delete poImage1;
    if (NULL != poScanIn)
    delete poScanIn;
    if (nStat == nNumImages)
    nStat = 0;
    
    return nStat;
}

void CDTMainAverage::vError(const int nErrorNum, const Cstring& sMessage)
{
    if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
    
    cout << "\ndtaverage re-averages a scan or series of images.\n"
       << "          It rejects outliers in the averaging process and saves\n"
       << "          rejected pixel positions in a separate image file.\n"
       << "dtaverage - Usage:\n"
       << "dtaverage [-scan scanfile] [ options ...]\n\n"
       << "Command line options:     Description:\n\n"
       << " -scan   sName            File sName (no default) has the scan\n"
       << "                          definition.\n\n"
       << " -template sInTemplate    Input image file template. If sInTemplate\n"
       << "                          contains ? be sure to enclose it in quotes.\n"
       << "                          Default: image.???\n\n"
       << " -start  nSeqStart        Sequence start number. Default: 1.\n\n"
       << " -inc    nSeqIncr         Sequence increment. Default: 1.\n\n"
       << " -num    nNumImages       Number of images to process.  Default: 1000.\n\n"
       << " -underlay\n"
       << " -overlay                 The default operation is to average images,\n"
       << " -add                     but these operations force an underlay, overlay,\n"
       << "                          or add  operation similar to that in dtdisplay.\n\n"
       << " -avg    sAvgFile         Input averaged image file name. Default is\n"
       << "                          formed from the input file template by\n"
       << "                          stripping off any leading directory and\n"
       << "                          appending av to the name.\n\n"
       << " -sd     sSDFile          Input standard deviation image file name. Default is\n"
       << "                          formed from the input file template by\n"
       << "                          stripping off any leading directory and\n"
       << "                          appending sd to the name.\n\n"
       << " -sigma  fSigma           Pixels which differ by more than fSigma from\n"
       << "                          the average value will be excluded from the\n"
       << "                          calculation of the new average and standard\n"
       << "                          deviation images.  Default 3.\n\n"
       << " -maxsd  fMaxSD           Any pixel input standard deviation larger than\n"
       << "                          fMaxSD will be set to fMaxSD. Default 1000000.\n\n"
       << " -minsd  fMinSD           Any pixel input standard deviation smaller than\n"
       << "                          fMinSD will be set to fMinSD. Default 1.\n\n"
       << " -min    fMin             Pixels with values below fMin will be excluded\n"
       << "                          from the calculation of the new average and\n"
       << "                          standard deviation images.  Default -1000000.\n\n"
       << " -max    fMax             Pixels with values above fMax will be excluded\n"
       << "                          from the calculation of the new average and\n"
       << "                          standard deviation images.  Default 1000000.\n\n"

       << " -out    sOutFile         Output image file name. Default is\n"
       << "                          formed from the input file template by\n"
       << "                          stripping off any leading directory and\n"
       << "                          appending dtav|dtun|dtov|dtad to the name of\n"
       << "                          averaged|underlayed|overlayed|added image.\n\n"
       << " -outsd  sOutSDFile       Output standard deviation image file name. Default is\n"
       << "                          formed from the input file template by\n"
       << "                          stripping off any leading directory and\n"
       << "                          appending dtsd to the name.\n\n"
       << " -outrej sRejectFile      Output rejectfile name. Default is\n"
       << "                          formed from the input file template by\n"
       << "                          stripping off any leading directory and\n"
       << "                          appending rejects to the name.\n\n"
       << " -noskip                  Do not skip an image if not found.\n"
       << "                          Default: ignore missing images (skip).\n\n"
       << " -help                    Print this help text.\n\n";
    cout << "Examples:\n"
       << "dtaverage -template \"lys???.osc\" -start 1 -num 10\n"
       << "dtaverage -template \"lys???.osc\" -start 1 -num 10 -add -out lys0001add.img\n"
       << "dtaverage -template \"lys???.osc\" -start 1 -num 10 -underlay\n\n";
}




