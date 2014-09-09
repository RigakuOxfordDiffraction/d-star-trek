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
// dtnonunfedit.cc     Initial author: J.W. Pflugrath           15-Nov-1996
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
//    dtnonunfedit creates a new nonunf file from an old nonunf file,
//    an image and a minimum.  Pixels with values below the minimum in
//    in the image are flagged as bad in the new nonunf file.
//    Example:  % dtnonunfedit nonunf.old image min nonunf.new
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

//#include <iostream.h>
#include "Cstring.h"
#include "Cimage.h"
#include "Cnonunf.h"

using std::cout;
using std::cerr;
using std::endl;
using std::flush;

//+Public functions

void vError(const int nErrorNum=0, const Cstring& sError=NULL);

//+ Code begin

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  Cstring sNonunfOld;
  Cstring sImage;
  Cstring sNonunfNew;
  float   fMin, fMax;
  int     nMin;
  int     nMax;

  int      nStat;
  Cimage  *poNonunf;
  Cimage  *poImage;
  int     nArg2;
  int     nBad;
  int     nNumBad;

  // Parse command line arguments

  argc--; argv++;

  vDtrekSetModuleName("dtnonunfedit");
  vPrintCopyrightInfo();

  if (4 > argc)
    {
      // Not enough input arguments and thats an error!

      vError(1, "ERROR - not enough input arguments!");
    }

  sNonunfOld = (const char*) argv[0];
  sNonunfNew = (const char*) argv[3];
  sImage     = (const char*) argv[1];
  nStat      = sscanf(argv[2], "%d", &nArg2);
  fMin      = (float)nArg2;
  fMax      = (float)abs(nArg2);

  if (0 <= nArg2)
    {
      cout << "\nINFO: Pixels with values less than " << fMin
           << "\n      in image "
           << sImage << " will flag the\n      corresponding pixel in "
           << sNonunfNew << " as bad.\n";
    }
  else
    {
      cout << "\nINFO: Pixels with values greater than " << fMax 
           << "\n      in image "
           << sImage << " will flag the\n      corresponding pixel in "
           << sNonunfNew << " as bad.\n";
    }

  Cstring sCBF;
  sCBF = sGetEnv(D_K_DtrekCBFPixDataType);
  if (0 < sCBF.length())
    {
      cout << "WARNING, possible CBF images will not be used as long int.\n"
           << flush;
      nPutEnv(D_K_DtrekCBFPixDataType, "2");
    }

  poNonunf = new Cimage(sNonunfOld);
  if (!poNonunf->bIsAvailable())
    {
      vError(1, "ERROR reading old nonunf image.");
    }

  poImage = new Cimage(sImage);
  if (!poImage->bIsAvailable())
    {
      vError(1, "ERROR reading test image.");
    }

  if (1 != nStat)
    {
      vError(1, "ERROR reading minimum value.");
    }

  poImage->nSetNextPixel(0,0);
  poNonunf->nSetNextPixel(0,0);

  int nDim0, nDim1;
  int nNumPixels;
  int i, j;
  float fBad;
  float fTest;

  nBad = 0;  // In case next line does not change nBad

  if (eImage_I4 == poNonunf->nGetDataType())
    {
      // If the data type of the image is long int, then we want
      // to use -1 as the default value for bad pixels.  This can
      // be overridden later.

      nBad = -1;
    }
  poNonunf->m_oHeader.nGetValue(Cnonunf::ms_sNonunfFlag1, &nBad);
  fBad = (float) nBad;

  (void) poNonunf->nGetDimensions(&nDim0, &nDim1);
  nNumPixels = nDim0 * nDim1;
  nNumBad = 0;
  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
	{
	  fTest = (poImage->*poImage->prfGetPixel)(i, j);
	  if (0 <= nArg2)
	    {
	      if (fTest < fMin)
		{
		  (poNonunf->*poNonunf->prnSetPixel)(i, j, fBad);
		  nNumBad++;  // This could have been bad to start with
		}
	    }
	  else
	    {
	      // nArg2 was negative!  So flag HOT pixels
	      if (fTest > fMax)
		{
		  (poNonunf->*poNonunf->prnSetPixel)(i, j, fBad);
		  nNumBad++;  // This could have been bad to start with
		}
	    }
	}
    }
  cout << "\nImage is " << nDim0 << " by " << nDim1 
       << " pixels.  Total number of pixels examined: " << nNumPixels << endl;
  Cstring sPercent;
  fBad = (float)nNumBad / (float)nNumPixels * 100.0;
  sPercent = Cstring(fBad, 5, 2);

  if (0 <= nArg2)
    {
      cout << "\nINFO: There were " << nNumBad << " pixels less than " << fMin
           << "\n      in " << sImage << " (" 
           << sPercent
           << "%) that were flagged as bad"
           << " in " << sNonunfNew << '\n' << endl;
    }
  else
    {
      cout << "\nINFO: There were " << nNumBad << " pixels greater than " << fMax
           << "\n      in " << sImage << " (" 
           << sPercent
           << "%) that were flagged as bad"
           << " in " << sNonunfNew << '\n' << endl;
    }

  // Write new nonunf image

  nStat = poNonunf->nWrite(sNonunfNew);

  if (0 != nStat)
    cerr << "\ndtnonunfedit: Error writing new nonunf image: " << nStat << "!\n"; 
  else
    cerr << "\ndtnonunfedit: Done.\n";

  delete poImage;
  delete poNonunf;
  if (0 < sCBF.length())
    {
      // Restore this

      nPutEnv(D_K_DtrekCBFPixDataType, sCBF);
    }
  exit (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum) 
    {
      cerr << sMessage << endl;
    }
  cout << "\n Usage:\n"
       << " dtnonunfedit sNonunfOld sNonunfEdit nMinGood sNonunfNew\n\n"
       << " Options:     Description:\n\n"
       << " sNonunfOld   An existing non-uniformity file.\n\n"
       << " sNonunfEdit  A file used to edit the existing non-uniformity file.\n"
       << "              Pixels with values less than nMinGood in this file will\n"
       << "              be flagged as bad in the new non-uniformity file.\n\n"
       << " nMinGood     Minimum allowed value for a good pixel in sNonunfEdit.\n"
       << "              If nMinGood < 0, then pixels in sNonunfEdit greater\n"
       << "              than the absolute value of nMinGood are flagged as bad.\n\n"
       << " sNonunfNew   The new non-uniformity file created by dtnonunfedit.\n\n";
  exit (nErrorNum);
}
