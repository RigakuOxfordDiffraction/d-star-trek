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
// dtraxis2refln.cc     Initial author: J.W. Pflugrath           08-Feb-1996
//    This reads a R-AXIS style .int file and writes a D*TREK style 
//    reflection file.
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
//    dtraxis2refln reads from the command line the filename of a R-AXIS
//    .int file to read (default dtrefln001.int) and the name of a d*TREK
//    reflection file to write (default dtrefln001.ref).
//    It reads in the input file and writes out the output file in the proper
//    format.
//
//    Example:  % dtraxis2refln  dtrefln001.int dtrefln001.ref 
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
  Cstring  sIn  = "dtrefln001.int";
  Cstring  sOut = "dtrefln001.ref";
  int      nStat;

  cerr << "\n\ndtraxis2refln: Copyright (c) 2006 Rigaku\n\n";

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

  Creflnlist oReflnlist;  // Create the output reflection list.

  // Add the fields we need ...

  (void) oReflnlist.nExpandGetField(oReflnlist.ms_sfObsXmm);
  (void) oReflnlist.nExpandGetField(oReflnlist.ms_sfObsYmm);
  (void) oReflnlist.nExpandGetField(oReflnlist.ms_sfObsRotMid);
  int nFI_nPartial = oReflnlist.nExpandGetField("nPartial");

  // Open input file

  int nSize     = sIn.length();
  int nFile     = 1;
  (void) dskbor(&nFile, sIn.string(), &nSize, &nStat);

  if (0 != nStat)
    {
      cerr << "ERROR opening input file " << sIn << "!\n";
      return (nStat);
    }

  // Read raxis.int style header

  nSize = 4;
  float fMaxRot, fMinRot;
  float fTemp;
  int   nTemp;
  int   nNumReflns;

//Record 1

// Code for X,Y, LP corrected Int, SigI
  nSize = sizeof(int);
  dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  if (3 != nTemp) 
    {
      cerr << "INPUT file" << sIn << " is not type 3!\n";
      cerr << "It has type " << nTemp << '\n';
      exit (1);
    }

  dskbr(&nFile, (char *)&nNumReflns, &nSize, &nStat); // Number of reflns
  if (0 != nStat) return (nStat);

  nSize = sizeof(float);
  dskbr(&nFile, (char *)&fMinRot, &nSize, &nStat);    // Min rotation angle
  if (0 != nStat) return (nStat);

  dskbr(&nFile, (char *)&fMaxRot, &nSize, &nStat);    // Max rotation angle
  if (0 != nStat) return (nStat);

  fMaxRot = (fMinRot + fMaxRot) * 0.5;

  nSize = sizeof(int);
  dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);

//Record 2
  fTemp = 2.0;                                       // Dummy overlap
  nSize = sizeof(float);
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  nTemp = 65535;                                     // Dummy overload
  nSize = sizeof(int);
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);
  if (0 != nStat) return (nStat);

  nTemp = 20;
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy pixel overloads
  if (0 != nStat) return (nStat);

  fTemp = 20.0;
  nSize = sizeof(float);
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy central region
  if (0 != nStat) return (nStat);

//Record 3
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);    // Dummy primary beam
  if (0 != nStat) return (nStat);

  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);    // Dummy primary beam
  if (0 != nStat) return (nStat);
    
  fTemp = 0.0;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy omega
  if (0 != nStat) return (nStat);

  fTemp = 100.0;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy distance
  if (0 != nStat) return (nStat);

  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);

//Record 4
  fTemp = 1.0;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy YScale
  if (0 != nStat) return (nStat);

  fTemp = 0.0;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy tilt
  if (0 != nStat) return (nStat);

  fTemp = 0.0;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy twist
  if (0 != nStat) return (nStat);

  fTemp = 0.0;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Dummy bulge
  if (0 != nStat) return (nStat);

  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);

//Record 5
  nTemp = 10;
  nSize = sizeof(int);
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nxs
  if (0 != nStat) return (nStat);

  nTemp = 10;
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nys
  if (0 != nStat) return (nStat);

  nTemp = 3;
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy ncr
  if (0 != nStat) return (nStat);

  nTemp = 1;
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nrx
  if (0 != nStat) return (nStat);

  nTemp = 1;
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Dummy nry
  if (0 != nStat) return (nStat);

//Record 6
  fTemp = 1.0;
  nSize = sizeof(float);
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // rms limit
  if (0 != nStat) return (nStat);

  fTemp = 1.0;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // rms limit
  if (0 != nStat) return (nStat);

  nTemp = 1;
  nSize = sizeof(int);
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // dummy flag 
  if (0 != nStat) return (nStat);

  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // padding
  if (0 != nStat) return (nStat);
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // padding
  if (0 != nStat) return (nStat);

//Record 7
  fTemp = 0.5;
  nSize = sizeof(float);
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // dummy rms
  if (0 != nStat) return (nStat);

  fTemp = 0.5;
  (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);      // dummy rms
  if (0 != nStat) return (nStat);


  nSize = sizeof(int);
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);     // Accepted full
  if (0 != nStat) return (nStat);

  nTemp = 0;
  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Accepted partials
  if (0 != nStat) return (nStat);

  (void) dskbr(&nFile, (char *)&nTemp, &nSize, &nStat);      // Padding
  if (0 != nStat) return (nStat);
//

  // Now read reflection records

  int   nH, nK, nL, nPack;
  int   nP;
  float fX, fY;
  float fInt, fSigI;

  Crefln oRefln(&oReflnlist);

  for (int i = 0; i < nNumReflns; i++)
    {
      (void) dskbr(&nFile, (char *)&nPack, &nSize, &nStat);
      if (0 != nStat) return (nStat);

      nH = (255 & (nPack >> 20)) - 128;
      nK = (255 & (nPack >> 12)) - 128;
      nL = (255 & (nPack >>  4)) - 128;
      nP = ( 15 & (nPack >> 28));

      oRefln.vSetH(nH);
      oRefln.vSetK(nK);
      oRefln.vSetL(nL);

      oRefln.vSetField(nFI_nPartial, nP);

      nSize = sizeof(float);
      (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);
      if (0 != nStat) return (nStat);
      oRefln.vSetField(oReflnlist.m_nFI_fObsXmm, fTemp);

      (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);
      if (0 != nStat) return (nStat);
      oRefln.vSetField(oReflnlist.m_nFI_fObsYmm, fTemp);
      oRefln.vSetField(oReflnlist.m_nFI_fObsRotMid, fMaxRot);

      (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);
      if (0 != nStat) return (nStat);
      oRefln.vSetIntensity(fTemp);

      (void) dskbr(&nFile, (char *)&fTemp, &nSize, &nStat);
      if (0 != nStat) return (nStat);
      oRefln.vSetSigmaI(fTemp);

      oReflnlist.nInsert(&oRefln);
    }

  // Close the input file
  
  (void) dskbcr(&nFile, &nStat);
  
  // Write out the reflnlist

  nStat = oReflnlist.nWrite(sOut);

  return (nStat);
}
