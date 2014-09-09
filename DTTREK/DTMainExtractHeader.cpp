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
// DTMainExtractHeader.cpp  Initial author: RB           15-Mar-2004
//
//  This file contains the member functions of class CDTMainExtractHeader
//

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

#include "DTMainExtractHeader.h"
#include "Cimage_header.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

int CDTMainExtractHeader::nExecute(unsigned int argc, char* argv[])
{
  int          nStat;
  Cstring      sImageFile  = "scan.img";
  Cstring      sHeaderFile = "scan.head";

  Cimage_header *poHeader;

  poHeader    = NULL;

  vDtrekSetModuleName("dtextractheader");
  vPrintCopyrightInfo();

  // Parse command line arguments
  
  argc--; argv++;

  nStat = 0;

  if (2 <= argc)
    {
      sImageFile  = (const char*) argv[0];
      sHeaderFile = (const char*) argv[1];
      if ( ("-help" == sImageFile) || ("-h" == sImageFile) ) 
	{
	  DTREK_ERROR(0, "");
	}
    }

  // Create default scan and find objects;

  poHeader = new Cimage_header(sImageFile);
  if (poHeader->bIsAvailable())
    {
      (void) poHeader->nReplaceValue(Cimage_header::ms_sSize1, 0);
      (void) poHeader->nReplaceValue(Cimage_header::ms_sSize2, 0);
      nStat = poHeader->nWrite(sHeaderFile);
    }
  else
    {
      DTREK_ERROR(1, "ERROR: header not available!");
    }
  delete poHeader;
  poHeader = NULL;
  return (nStat);
}

void CDTMainExtractHeader::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
  cout << "\ndtextractheader - Usage:\n"
       << "      dtextractheader imagefile newheaderfile\n\n";
}
