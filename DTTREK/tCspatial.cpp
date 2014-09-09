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
// tCspatial.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//    This file tests the class Cspatial.
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

#include "Cspatial.h"

//+Code begin

int main()
{
  Cspatial *poSPD;
  int      nStat;
  float fP1, fP2, fX, fY;
  Cstring sFileset;

  Cimage_header *poHeader;

  // Zeroeth check, initialization from a image header

  cout << "Enter name of header file to get spatial distortion info from: ";
  cin >> sFileset;
  poHeader = new Cimage_header(sFileset);
  
  if (!poHeader->bIsAvailable())
    {
      cerr << "Header not available!\n";
      exit (-1);
    }

  // Get the prefix for the 1st detector, if there isn't one, then
  //  the null string will be returned, which is OK 

  Cstring sDetName;
  (void) poHeader->nGetValue(Cdetector::ms_sDetectorNames, 1, &sDetName);

  poSPD = new Cspatial(*poHeader, sDetName);
  delete poHeader;
  if (!poSPD->bIsAvailable())
    {
      cerr << "Error reading spatial distortion info!\n";
      delete poSPD;
    }
  else
    {
      (void) poSPD->nList();
      delete poSPD;
    }
  cout << "\n\n";


  cout << "Enter name of spatial distortion fileset: ";
  cin >> sFileset;
  poSPD = new Cspatial(sFileset);
  if (!poSPD->bIsAvailable())
    {
      cerr << "Error reading spatial distortion info!\n";
      delete poSPD;
    }
  (void) poSPD->nList();

  poSPD->vSetVerbose(3);
  nStat = 0;
  while (0 == nStat)
    {
      cout << "Enter pixel coordinates: ";
      cin >> fP1 >> fP2;
      nStat = poSPD->nPixeltoMM(fP1, fP2, &fX, &fY);
      if (0 == nStat)
	{
	  nStat = poSPD->nMMtoPixel(fX, fY, &fP1, &fP2);
	  cout << "MM coordinates: " << fX << ", " << fY << endl;
	  cout << " back to pixels: " << fP1 << ", " << fP2 << endl;
	  if (0 != nStat) cout << " Error in MMtoPixel.\n";
	}
      else
	cout << "  Error in PixeltoMM.\n";
    }

  nStat = 0;
  while (0 == nStat)
    {
      cout << "Enter MM coordinates: ";
      cin >> fX >> fY;
      nStat = poSPD->nMMtoPixel(fX, fY, &fP1, &fP2);
      if (0 == nStat) 
	{
	  nStat = poSPD->nPixeltoMM(fP1, fP2, &fX, &fY);
	  cout << "Pixel coordinates: " << fP1 << ", " << fP2 << endl;
	  cout << "  back to mm:      " << fX << ", " << fY << endl;
	  if (0 != nStat) cout << " Error in PixeltoMM.\n";
	}
      else
	cout << " Error in MMtoPixel.\n";
    }

  delete poSPD;

  cout << " Done with Stanton interpolation tables, now use simple spd.\n";

  float Beam0, Beam1, PxSize0, PxSize1;
  Beam0   = 260.;
  Beam1   = 250.;
  PxSize0 = 0.2;
  PxSize1 = 0.2;

  poSPD = new Cspatial(Beam0, Beam1, PxSize0, PxSize1, 512, 512,
		       0.0, -1.0, 1.0, 0.0);

  (void) poSPD->nList();

  nStat = 0;
  while (0 == nStat)
    {
      cout << "Enter pixel coordinates: ";
      cin >> fP1 >> fP2;
      nStat = poSPD->nPixeltoMM(fP1, fP2, &fX, &fY);
      if (0 == nStat)
	{
	  nStat = poSPD->nMMtoPixel(fX, fY, &fP1, &fP2);
	  cout << "MM coordinates: " << fX << ", " << fY << endl;
	  cout << " back to pixels: " << fP1 << ", " << fP2 << endl;
	  if (0 != nStat) cout << " Error in MMtoPixel.\n";
	}
      else
	cout << "  Error in PixeltoMM.\n";
    }
  nStat = 0;
  while (0 == nStat) 
    {
      cout << "Enter MM coordinates: ";
      cin >> fX >> fY;
      nStat = poSPD->nMMtoPixel(fX, fY, &fP1, &fP2);
      if (0 == nStat)
	{
	  nStat = poSPD->nPixeltoMM(fP1, fP2, &fX, &fY);
	  cout << "Pixel coordinates: " << fP1 << ", " << fP2 << endl;
	  cout << "  back to mm:      " << fX << ", " << fY << endl;
	  if (nStat != 0) cout << " Error in PixeltoMM.\n";
	}
      else
	cout << " Error in MMtoPixel.\n";
    }
  
  delete poSPD;

  return 0;
}

