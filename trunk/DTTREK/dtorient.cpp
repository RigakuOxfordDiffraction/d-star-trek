//
// Copyright (c) 1997 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtorient.cc     Initial author: J.W. Pflugrath        27-May-1997
//    Orient the crystal in specific ways once an orientation matrix is known.
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
//    Example:  % dtorient header_file request
//                where request is: "a||x" "b||y"
//                                  "crys_vec1||lab_vec1" "crys_vec2||lab_vec2"
//                              or: rot 20 x
//                                  rot angle lab_vec 
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream
#include "Cstring.h"
#include "Dtrek.h"
#include "Cimage_header.h"
#include "Ccrystal.h"
#include "Corient.h"
#include "dtrekvec.h"

//+Public functions

#ifdef SSI_PC
#include "CrclIncludes.h"
#define main   dtorient_main
#define vError dtorient_vError
#endif

using std::cout;
using std::endl;
using std::flush;

void vError(const int nErrorNum=0, const Cstring& sError=NULL);
int  nNTHUfixup(Cimage_header&);
int  nNTHUfixup(Cstring&, Cimage_header&, Cimage_header&, Ccrystal*);

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  int            i;
  int            nStat;
  int            nStatNTHU;
  Cstring        sTemp;
  Cimage_header *poHeader     = NULL;
  Cimage_header  oHeaderNTHU;
  Ccrystal      *poCrystal    = NULL;
  Corient       *poOrient     = NULL;
  Cstring        sOrientWish;

  vDtrekSetModuleName("dtorient");
  vPrintCopyrightInfo();
  //cout << "\ndtorient:  Copyright (c) 2002 Rigaku/MSC, Inc.\n";
  //cout << D_K_DTREKVersion << endl;

  // Copy command line to output log

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;

  // Get header file to start with
  
  argc--;
  argv++;

  nStat = 1;

  if (1 > argc)
    {
      // No header file argument and thats an error!

      vError(1, "ERROR - no header filename!\n");
    }

  sTemp = (const char*) argv[0];
  if ( ("-help" == sTemp) || ("-h" == sTemp) ) 
    {
      vError(0, "");
    }
  else
    {
      // Read the name of a file and read the header only from that image file

      poHeader     = new Cimage_header(sTemp);   // Create (read) image header
    }

  if (!poHeader->bIsAvailable())
    {
      cout << "ERROR in dtorient: Header not available!\n";
      delete poHeader;
      poHeader = NULL;
      nStat    = 1;
    }
  else
    {
      // nNTHUfixup() returns 0 if it is an NTHU 4-axis header was input
      //                        and a modified poHeader
      //              returns 1 if it is NOT an NTHU 4-axis
      //                        and an unchanged poHeader

      oHeaderNTHU = *poHeader;          // Save copy of original
      nStatNTHU = nNTHUfixup(*poHeader);
      nStat = 0;
      poCrystal = new Ccrystal(*poHeader);
      if (!poCrystal->bIsAvailable())
	{
	  cout << "ERROR in dtorient: Crystal info not available!\n";
	  nStat = 1;
	}
    }
  if (0 != nStat)
    {
      if (NULL != poCrystal)
	{
	  delete poCrystal;
	  poCrystal = NULL;
	}
      vError(3, "");      // This exits or returns
    }

  poOrient  = new Corient(*poHeader);

  // Parse command line options

  int nMode = 0;
  sOrientWish = "";
  for (i = 1; i < argc; i++)
    {
      sTemp = (const char*) argv[i];
      cout << "Command line argument is: >>>" << sTemp << "<<<\n";
      if ( ("-rot" == sTemp) || ("rot" == sTemp) )
	{
	  sOrientWish = " rot";
	  nMode       = -1;
	}
      else if (0 == nMode)
	{
	  sOrientWish = sTemp + ';';
	  nMode = 1;
	}
      else
	{
	  sOrientWish = sOrientWish + ' ' + sTemp;
	}
    }

  nStat = poOrient->nCalcAngles();
  nStat = poOrient->nList();

//  cout << "\nEnter orient request: ";
//  getline(cin, sOrientWish);

  nStat = poOrient->nSetOrient(sOrientWish);
  nStat = poOrient->nCalcOrient();
  nStat = poOrient->nPickSoln();
  if (0 <= nStat)
    {
      nStat = poOrient->nUpdateHeader(poHeader);
      if ( (0 == nStatNTHU) && (sOrientWish.contains(')')) )
	{
	  // We need to modify poHeader to put a omega at the 2theta
	  // value of a refln given on the command line.

	  nStatNTHU = nNTHUfixup(sOrientWish, oHeaderNTHU, *poHeader, poCrystal);
	  if (0 != nStatNTHU)
	    cout << "ERROR orienting reflection!\n";
	}
      nStat = poOrient->nList();
      nStat = poHeader->nWrite(sDtrekGetPrefix() + "dtorient.head");
    }
  delete poOrient;
  delete poHeader;
  if (NULL != poCrystal)
    {
      delete poCrystal;
      poCrystal = NULL;
    }
  return (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << '\n';
    }
  cout << "\ndtorient - Usage:\n"
       << "dtorient sHeaderfile [ options ...]\n\n"
       << "Command line options: Description:\n\n"
       << "sHeaderfile         Name of a file with a header that has the crystal\n"
       << "                    properties.\n\n"
       << "-help               Print this help text.\n\n"
       << "Information about the crystal unit cell and orientation is\n"
       << "determined from the header.  If the header is incorrect,\n"
       << "it may be edited with dtheaderedit or dtprocess.  A new header,\n"
       << "*dtorient.head, is written with the results.\n\n"
       << "Examples:\n\n"
       << "   dtorient dtrefine.head \"a||omega\" \"b||source\"\n"
       << "   dtorient dtrefine.head rot 10 source\n"
       << "   dtorient dtrefine.head \"a||x\" \"b||z\"\n"
       << "   dtorient dtrefine.head \"( 1 4 2 ) || omega\" \"b||z\"\n"
       << "   dtorient dtrefine.head \"[ 1 1 0 ] ||omega\" \"b||z\"\n";


#ifndef SSI_PC
  exit (nErrorNum);
#endif
}


