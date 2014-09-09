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
 
//+Include files

#include "dtrekvec.h"
#include "Cscan.h"
#include "Cimage.h"
#include "Creflnlist.h"
#include "Crefln.h"
#include "Ccrystal.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

typedef struct _tagScalepackStringRec
  {
    char    ac_nH[5];
    char    ac_nK[5];
    char    ac_nL[5];
    char    ac_fInt[9];
    char    ac_fSigmaI[9];
    char    ac_fIntPlus[9];
    char    ac_fSigmaIPlus[9];
    char    ac_fIntMinus[9];
    char    ac_fSigmaIMinus[9];
  } tagScalepackStringRec;

//+Code begin

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  
  // Convert a scalepack .sca file to a Creflnlist

  int   i;
  Creflnlist *poReflnlist;
  tagScalepackStringRec tSCAline;
  Cstring sTemp;
  Cstring a15sTokens[15];

  int     nH, nK, nL;
  float   fInt, fSigmaI, fIntPlus, fIntMinus, fSigmaIPlus, fSigmaIMinus;
  int nReadIn = 0;
  int nStat;

  vDtrekSetModuleName("dtscalepack2refln");
  cout << "\ndtscalepack2refln:  Copyright (c) 2006 Rigaku\n";
  cout << D_K_DTREKVersion << endl << flush;
  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;
  // Parse command line arguments

  if (3 > argc)
    {
      cout << "Usage: dtscalepack2refln .sca_file .ref_file\n"
           << flush;
      return (1);
    }

  poReflnlist = new Creflnlist();
  (void) poReflnlist->nExpandGetField(Creflnlist::ms_sfIntensityPlus);
  (void) poReflnlist->nExpandGetField(Creflnlist::ms_sfSigmaIPlus);
  (void) poReflnlist->nExpandGetField(Creflnlist::ms_sfIntensityMinus);
  (void) poReflnlist->nExpandGetField(Creflnlist::ms_sfSigmaIMinus);

  Crefln oRefln(poReflnlist);

  ifstream oIn( sTransSymbol(argv[1]).string());
  if ((oIn.rdbuf()->is_open()))
    {
      if (!oIn)
	{
	  nStat = 1;
	  return (nStat);
	}
      else
	{
	  nStat = 0;

	  // Read integer on first line

	  getline(oIn, sTemp);

	  int nSkip = 1;
	  nSkip = atoi(sTemp.string());

	  // Skip the required number of lines

	  for (i = 0; i < nSkip; i++)
	    getline(oIn, sTemp);

	  // Next line has the cell and spacegroup.  Don't do anything with
	  // them for now.

	  getline(oIn, sTemp);
	  split(sTemp, a15sTokens, 6, ' ');

	  // Read the reflns

	  bool bNotDone = TRUE;
	  while (bNotDone)
	    {
	      getline(oIn, sTemp);

	      // We have to do some special parsing of the refln line because
	      // the numbers sometimes run together without spaces between them

	      memset(&tSCAline, char(0), sizeof(tSCAline));
	      nStat = sscanf(sTemp.string(), "%4c%4c%4c%8c%8c%8c%8c",
			     &tSCAline.ac_nH, &tSCAline.ac_nK, &tSCAline.ac_nL, 
			     &tSCAline.ac_fIntPlus, &tSCAline.ac_fSigmaIPlus,
			     &tSCAline.ac_fIntMinus, &tSCAline.ac_fSigmaIMinus);
	      if (5 <= nStat)
		{
		  sscanf(tSCAline.ac_nH,            "%4d", &nH);
		  sscanf(tSCAline.ac_nK,            "%4d", &nK);
		  sscanf(tSCAline.ac_nL,            "%4d", &nL);
		  sscanf(tSCAline.ac_fIntPlus,      "%8f", &fIntPlus);
		  sscanf(tSCAline.ac_fSigmaIPlus,   "%8f", &fSigmaIPlus);
		}
	      if (7 <= nStat)
		{
		  sscanf(tSCAline.ac_fIntMinus,     "%8f", &fIntMinus);
		  sscanf(tSCAline.ac_fSigmaIMinus,  "%8f", &fSigmaIMinus);
		}
	      else
		{
		  fIntMinus    =  0.0;
		  fSigmaIMinus = -1.0;
		}

	      bNotDone = (999 != nH) && (5 <= nStat);
	      if (bNotDone)
		{
		  nReadIn++;
		  oRefln.vSetH(nH);
		  oRefln.vSetK(nK);
		  oRefln.vSetL(nL);

		  // Build fInt and fSigmaI

		  if (5 == nStat)
		    {
		      fInt    = fIntPlus;
		      fSigmaI = fSigmaIPlus;
		    }
		  else if (7 <= nStat)
		    {
		      if ( (0.0 == fIntPlus) && (-1.0 == fSigmaIPlus) )
			{
			  fInt    = fIntMinus;
			  fSigmaI = fSigmaIMinus;
			}
		      else
			{
			  fInt    = 0.5 * (fIntPlus + fIntMinus);
			  fSigmaI = sqrt( (double)( fSigmaIPlus  * fSigmaIPlus
					          + fSigmaIMinus * fSigmaIMinus));
			}
		    }

		  oRefln.vSetIntensity(fInt);
		  oRefln.vSetSigmaI(fSigmaI);

		  oRefln.vSetField(poReflnlist->m_nFI_fIntensityPlus,
				   fIntPlus);
		  oRefln.vSetField(poReflnlist->m_nFI_fSigmaIPlus,
				   fSigmaIPlus);
		  oRefln.vSetField(poReflnlist->m_nFI_fIntensityMinus,
				   fIntMinus);
		  oRefln.vSetField(poReflnlist->m_nFI_fSigmaIMinus,
				   fSigmaIMinus);
		  poReflnlist->nInsert(&oRefln);
		}
	    }
	  cout << "Number of reflections read in: " << nReadIn << endl;
	  poReflnlist->nWrite(Cstring(argv[2]));
	}
      oIn.close();                 // Make sure we always close the stream!
    }
  delete poReflnlist;
  return (0);
}

