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
// dtrefln2raxis.cc     Initial author: J.W. Pflugrath           11-Jan-1996
//    This reads a D*TREK style reflection file and writes a R-AXIS style 
//    .int file.
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
//    dtrefln2raxis reads from the command line the filename of a D*TREK 
//    reflection file to read (default dtrefln.ref) and the name of a R-AXIS 
//    .int file to write (default dtrefln001.int).
//    It reads in the input file and writes out the output file in the proper
//    format.
//
//    Example:  % dtrefln2raxis  dtrefln.ref dtrefln001.int
//
//+ToDo
//
//   Better error messages need implementing.
//   It would be better to get some info from an image file header.
//
//+Include files

#include <iostream.h>
#include "Cstring.h"
#include "Crefln.h"
#include "Creflnlist.h"

#include "dskio.h"

//+Function prototypes

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
     
{
  Cstring  sIn  = "dtrefln.ref";
  Cstring  sOut = "dtrefln001.int";
  int      nStat;

  cerr << "\n\ndtrefln2raxis: Copyright (c) 2006 Rigaku\n\n";

  // Parse command line arguments, if any
  
  argc--; argv++; // Skip program name

  if (argc >= 1) 
    {
      sIn = (const char*) argv[0];
      argc--;
      argv++;
      if (argc >= 1)
	{
	  sOut = (const char*) argv[0];
	}
    }

  cout << "dtrefln2raxis: convert "
       << sIn  << " to "
       << sOut << '\n';

  Creflnlist oReflnlist(sIn);  // Create and read in the input reflection list.

  int nNumReflns = oReflnlist.nGetNumReflns();

  // Check reflection list for certain items, assume at least h,k,l are there

  nStat = 0;
  if (0 > oReflnlist.m_nFI_fIntensity)
    {
      cerr << "ERROR: Input reflection list is missing "
           << oReflnlist.ms_sfIntensity
           << " field!\n";
      nStat++;
    }
  if (0 > oReflnlist.m_nFI_fSigmaI)
    {
      cerr << "ERROR: Input reflection list is missing "
           << oReflnlist.ms_sfSigmaI
           << " field!\n";
      nStat++;
    }
  if (0 > oReflnlist.m_nFI_fObsPx0)
    {
      cerr << "ERROR: Input reflection list is missing "
           << oReflnlist.ms_sfObsPx0
           << " field!\n";
      nStat++;
    }
  if (0 > oReflnlist.m_nFI_fObsPx1)
    {
      cerr << "ERROR: Input reflection list is missing "
           << oReflnlist.ms_sfObsPx1
           << " field!\n";
      nStat++;
    }
  if (0 > oReflnlist.m_nFI_fObsRotMid)
    {
      cerr << "ERROR: Input reflection list is missing "
           << oReflnlist.ms_sfObsRotMid
           << " field!\n";
      nStat++;
    }
  if (0 != nStat)
    {
      cerr << "ERROR with input file " << sIn << "!\n";
      return (nStat);
    }

  int i;

  // Find min and max rotation range in the file
  // Find average value of pixel coordinates, too and use as center

  Crefln *poRefln = oReflnlist.poGetRefln(0);

  float fAvgPx0 = 0.0;
  float fAvgPx1 = 0.0;
  float fMinRot, fMaxRot, fRot;
  fMinRot = poRefln->fGetField(oReflnlist.m_nFI_fObsRotMid);
  fMaxRot = fMinRot;
  for (i = 1; i < nNumReflns; i++)
    {
      poRefln = oReflnlist.poGetRefln(i);
      fRot    = poRefln->fGetField(oReflnlist.m_nFI_fObsRotMid);
      if (fRot < fMinRot) fMinRot = fRot;
      if (fRot > fMaxRot) fMaxRot = fRot;
      fAvgPx0 = fAvgPx0 + poRefln->fGetField(oReflnlist.m_nFI_fObsPx0);
      fAvgPx1 = fAvgPx1 + poRefln->fGetField(oReflnlist.m_nFI_fObsPx1);
    }

  fAvgPx0 = fAvgPx0 / (float) nNumReflns;
  fAvgPx1 = fAvgPx1 / (float) nNumReflns;

  // Open output file

  int nSize     = sOut.length();
  int nFile     = 1;
  int nInitsize = 4096;     // Careful about this.  not used on unix, but vms...
  dskbow(&nFile, sOut.string(), &nSize, &nInitsize, &nStat);

  if (0 != nStat)
    {
      cerr << "ERROR opening output file " << sOut << "!\n";
      return (nStat);
    }

  // Write raxis.int style header

  nSize = sizeof(int);
  float fTemp;
  int   nTemp;

//Record 1
  nTemp = 3;                         // Code for X,Y, LP corrected Int, SigI
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&nNumReflns, &nSize, &nStat); // Number of reflns
  if (0 != nStat) return (nStat);

  nSize = sizeof(float);
  dskbw(&nFile, (char *)&fMinRot, &nSize, &nStat);    // Min rotation angle
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&fMaxRot, &nSize, &nStat);    // Max rotation anglep
  if (0 != nStat) return (nStat);

  nSize = sizeof(int);
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);

//Record 2
  fTemp = 2.0;                                       // Dummy overlap
  nSize = sizeof(float);
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  nTemp = 65535;                                     // Dummy overload
  nSize = sizeof(int);
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  nTemp = 20;
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy pixel overloads
  if (0 != nStat) return (nStat);

  fTemp = 20.0;
  nSize = sizeof(float);
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy central region
  if (0 != nStat) return (nStat);

//Record 3
  dskbw(&nFile, (char *)&fAvgPx0, &nSize, &nStat);    // Dummy primary beam
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&fAvgPx1, &nSize, &nStat);    // Dummy primary beam
  if (0 != nStat) return (nStat);
    
  fTemp = 0.0;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy omega
  if (0 != nStat) return (nStat);

  fTemp = 100.0;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy distance
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);

//Record 4
  fTemp = 1.0;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy YScale
  if (0 != nStat) return (nStat);

  fTemp = 0.0;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy tilt
  if (0 != nStat) return (nStat);

  fTemp = 0.0;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy twist
  if (0 != nStat) return (nStat);

  fTemp = 0.0;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy bulge
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);

//Record 5
  nTemp = 10;
  nSize = sizeof(int);
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nxs
  if (0 != nStat) return (nStat);

  nTemp = 10;
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nys
  if (0 != nStat) return (nStat);

  nTemp = 3;
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy ncr
  if (0 != nStat) return (nStat);

  nTemp = 1;
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nrx
  if (0 != nStat) return (nStat);

  nTemp = 1;
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nry
  if (0 != nStat) return (nStat);

//Record 6
  fTemp = 1.0;
  nSize = sizeof(float);
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // rms limit
  if (0 != nStat) return (nStat);

  fTemp = 1.0;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // rms limit
  if (0 != nStat) return (nStat);

  nTemp = 1;
  nSize = sizeof(int);
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // dummy flag 
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // padding
  if (0 != nStat) return (nStat);
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // padding
  if (0 != nStat) return (nStat);

//Record 7
  fTemp = 0.5;
  nSize = sizeof(float);
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // dummy rms
  if (0 != nStat) return (nStat);

  fTemp = 0.5;
  dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);      // dummy rms
  if (0 != nStat) return (nStat);

  nTemp = nNumReflns;
  nSize = sizeof(int);
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Accepted full
  if (0 != nStat) return (nStat);

  nTemp = 0;
  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Accepted partials
  if (0 != nStat) return (nStat);

  dskbw(&nFile, (char *)&nTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);
//
  // Now write reflection records

  int   nH, nK, nL, nPack;
  float fX, fY;
  float fInt, fSigI;

  for (i = 0; i < nNumReflns; i++)
    {
      poRefln = oReflnlist.poGetRefln(i);
      nPack   =   ((poRefln->nGetH()+128) << 20)
	        + ((poRefln->nGetK()+128) << 12)
		+ ((poRefln->nGetL()+128) << 4);
      nSize = sizeof(int);
      dskbw(&nFile, (char *)&nPack, &nSize, &nStat);
      if (0 != nStat) return (nStat);
      fTemp = poRefln->fGetField(oReflnlist.m_nFI_fObsPx0);
      nSize = sizeof(float);
      dskbw(&nFile, (char *)&fTemp,&nSize, &nStat);
      if (0 != nStat) return (nStat);
      fTemp = poRefln->fGetField(oReflnlist.m_nFI_fObsPx1);
      dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);
      if (0 != nStat) return (nStat);
      fTemp = poRefln->fGetIntensity();
      dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);
      if (0 != nStat) return (nStat);
      fTemp = poRefln->fGetSigmaI();
      dskbw(&nFile, (char *)&fTemp, &nSize, &nStat);
      if (0 != nStat) return (nStat);
    }

  // Close the output file
  
  (void) dskbcw(&nFile, &nStat);

  return (nStat);
}
