//
// Copyright (c) 1996 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtimagederiv.cc     Initial author: J.W. Pflugrath           02-May-1996
//    ...
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
//    Example:  % dtimagederiv [ -scan headerfile] [options ...]
//                  
//    Command line options:               Default:
//                -template sTemplate     image.???
//                -start    nSeqStart     1
//                -inc      nSeqIncr      1
//                -num      nNumImages    1000
//                -out      sOutTemplate  (current dir) sTemplate + "sub"
//                -help
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Cscan.h"

//+Code begin

//+Definitions, constants, and initialization of static member variables

//+Public functions

void vError(const int nErrorNum=0, const Cstring& sError=NULL);

//

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
     
{
  int          nStat;
  Cstring      sInTemplate  = "image.???";
  Cstring      sOutTemplate;
  int          nSeqStart     = 1;
  int          nSeqIncr      = 1;
  int          nNumImages    = 1000;

  Cimage_header *poHeader;
  Cimage      *poImage1;
  Cimage      *poImage2;
  Cimage      *poImageTemp;
  Cscan       *poScanIn;
  Cscan       *poScanOut;

  int          i, j, k;  // Loop counters
  float        fPixel1, fPixel2;
  int          nPixel1, nPixel2;
  int          nTemp;
  float        fTemp;
  Cstring      sTemp;

  poImage1  = NULL;
  poImage2  = NULL;
  poScanIn  = NULL;
  poScanOut = NULL;

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
	      poHeader = new Cimage_header(sInTemplate);
	      if (poHeader->bIsAvailable())
		{
		  poScanIn = new Cscan(*poHeader);
		  if (poScanIn->bIsAvailable())
		    {
		      nSeqStart  = poScanIn->nGetSeqNum(0);
		      nSeqIncr   = poScanIn->nGetSeqInc();
		      nNumImages = poScanIn->nGetNumImages();
		    }
		  else
		    {
		      vError(2, "ERROR: dtimagederiv - invalid scan in header!");
		    }
		  delete poScanIn;
		  poScanIn = NULL;
		}
	      delete poHeader;
	      poHeader = NULL;
	    }
	  else
	    {
	      vError(2, "ERROR: dtimagederiv - missing scan file name!");
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
		vError(6, "ERROR: dtimagederiv - invalid argument!");
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
		vError(6, "ERROR: dtimagederiv - invalid argument!");
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
		vError(6, "ERROR: dtimagederiv - invalid argument!");
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
	    vError(6, "ERROR: dtimagederiv - invalid argument!");
	}
      else if ("-out" == sTemp) 
	{
	  i++;
	  if (i < argc) 
	    {
	      sOutTemplate = argv[i];
	    }
	  else
	    vError(6, "ERROR: dtimagederiv - invalid argument!");
	}
      else if ("-help" == sTemp)
	{
	  vError(0, "");
	}
      else
	{
	  vError(0,"" );
	}
    }

  // Create default scan objects;

  poScanIn  = new Cscan();
  poScanOut = new Cscan();
  poScanIn->vSetTemplate(sInTemplate);
  poScanIn->poRotation->vSetIncrement(0.0);
  poScanIn->vSetSeqStart(nSeqStart);
  poScanIn->vSetSeqInc(nSeqIncr);
  poScanIn->vSetNumImgs(nNumImages);

  if ("" == sOutTemplate) 
    sOutTemplate = sInTemplate + "sub";
  // Should have a copy constructor!
  // Strip off directory
  while (0 <= sOutTemplate.find("/")) sOutTemplate = sOutTemplate.after("/");
  poScanOut->vSetTemplate(sOutTemplate);
  poScanOut->vSetSeqStart(poScanIn->nGetSeqNum(0));
  poScanOut->vSetSeqInc(poScanIn->nGetSeqInc());

//  poScanIn->nList();

  cout <<   "Template:  " << sInTemplate
       << "\nSeq start: " << nSeqStart
       << "\nSeq incr:  " << nSeqIncr
       << "\nNum imgs:  " << nNumImages
       << "\nOut templ: " << sOutTemplate << "\n\n";

  ////////////////////////////////////////////////////////////////////////
  // Done with command line arguments so now go do the deed
  ////////////////////////////////////////////////////////////////////////

  poScanIn->vInitSeqNum();
  poScanOut->vInitSeqNum();

  // Get first input image name and read in image

  nStat    = poScanIn->nGetImageName(&sInTemplate);
  if (0 != nStat)
    vError(0, "First image name invalid!");

  poImage1 = new Cimage(sInTemplate);

  if (!poImage1->bIsAvailable())
    vError(0, "First image not available!");

  // Get next input image name and read in image

  poScanIn->vNextSeqNum();
  nStat    = poScanIn->nGetImageName(&sOutTemplate);
  if (0 != nStat)
    vError(0, "Next image name invalid!");

  poImage2 = new Cimage(sOutTemplate);

  if (!poImage2->bIsAvailable())
    vError(0, "Next image not available!");

  int nDim0, nDim1;

  (void) poImage1->nGetDimensions(&nDim0, &nDim1);

  // Create empty output image

  Cimage oOutImage(nDim0, nDim1, eImage_I2);  

  nStat = 0;

  for (k = 1; (k < poScanIn->nGetNumImages()) && (0 == nStat); k++)
    {
      // Set pointer to first pixel in all images

      (void) poImage1->nSetNextPixel(0, 0);
      (void) poImage2->nSetNextPixel(0, 0);
      (void) oOutImage.nSetNextPixel(0, 0);
      
      if (   (eImage_uI2            == poImage1->eData_type)
	  && (eImage_no_compression == poImage1->eCompression) )
	{
	  // Try to speed things up for this specific case

	  for (i = 0; i < nDim0 * nDim1; i++)
	    {
	      nPixel1 = (int) poImage1->uiGetNextPixel();
	      nPixel2 = (int) poImage2->uiGetNextPixel();
	      oOutImage.vSetNextPixel((short int) (nPixel2-nPixel1));
	    }
	}
      else
	{
	  // This is the general case for all image data types

	  for (j = 0; j < nDim1; j++)
	    {
	      for (i = 0; i < nDim0; i++)
		{
		  fPixel1 = (poImage1->*poImage1->prfGetPixel)(i,j);
		  fPixel2 = (poImage2->*poImage1->prfGetPixel)(i,j);
		  oOutImage.vSetNextPixel((short int) (fPixel2-fPixel1));
		}
	    }
	}

      // Update header of output image, use first image for most of header info
      
      oOutImage.oHeader = poImage1->oHeader;
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sDataType, 
				      Cstring("short int"));
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sCompression, 
				      Cstring("None"));
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sDataType, 
				      Cstring("short int"));
      oOutImage.oHeader.nReplaceValue(Cimage_header::ms_sComment2, 
				      Cstring("dtimagederiv: ") + sInTemplate
                                      + " - " + sOutTemplate);

      // Get name of output image and write it out...

      nStat = poScanOut->nGetImageName(&sTemp);
      poScanOut->vNextSeqNum();
      if (0 == nStat) nStat = oOutImage.nWrite(sTemp);
      if (0 != nStat)
	vError(8, "ERROR - dtimagederiv writing output image!");

      // Now switch image pointers and file names

      poImageTemp = poImage1;
      poImage1    = poImage2;
      poImage2    = poImageTemp;
      sInTemplate = sOutTemplate;

      // Get next input image name and read in the image

      poScanIn->vNextSeqNum();
      (void) poScanIn->nGetImageName(&sOutTemplate);
      nStat = poScanIn->nGetImage(poImage2);
      if ( (0 != nStat) || !poImage2->bIsAvailable() )
	vError(0, "Image i not available!");
    }

  if (NULL != poImage1)
    delete poImage1;
  if (NULL != poImage2)
    delete poImage2;
  if (NULL != poScanIn)
    delete poScanIn;
  if (NULL != poScanOut)
    delete poScanOut;

  return (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cerr << sMessage << endl;
    }
  cout << "\ndtimagederiv - Usage:\n"
       << "dtimagederiv [-scan scanfile] [ options ...]\n\n"
       << "Command line options:     Description:\n\n"
       << " -scan   sName            File sName (no default) has the scan\n"
       << "                          definition.\n\n"
       << " -start  nSeqStart        Sequence start number. Default: 1.\n\n"
       << " -inc    nSeqIncr         Sequence increment. Default: 1.\n\n"
       << " -num    nNumImages       Number of images to process.  Default: 1000.\n\n"
       << " -template sInTemplate    Input image file template. If sInTemplate\n"
       << "                          contains ? be sure to enclose it in quotes.\n"
       << "                          Default: image.???\n\n"
       << " -out    sOutTemplate     Output image file template. Default is\n"
       << "                          formed from the input file template by\n"
       << "                          stripping off any leading directory and\n"
       << "                          appending sub to the name.\n\n"
       << " -help                    Print this help text.\n\n";

  exit (nErrorNum);
}



