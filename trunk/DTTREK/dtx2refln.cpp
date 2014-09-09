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

using namespace std;

typedef struct _tagDenzoXStringRec
  {
    char    ac_nH[5];
    char    ac_nK[5];
    char    ac_nL[5];
    char    ac_nPartial[3];
    char    ac_fInt[9];
    char    ac_fOtherInt[9];
    char    ac_fProfileChiSq[8];
    char    ac_fSigmaI[7];
    char    ac_fOblique[7];
    char    ac_fObs_pixel1[8];
    char    ac_fObs_pixel0[8];
    char    ac_fLPOfactor[7];
    char    ac_fStrength[9];
  } tagDenzoXStringRec;

//+Code begin

int main(int   argc,      // Command line argument count
         char *argv[])    // Pointers to command line args
{
  
  // Convert a denzo .x file to a Creflnlist

  int   i;
 
  int nFI_nPartialFlag, nFI_fOtherInt, nFI_fProfileChiSq, nFI_fLPOfactor;
  int nFI_fStrength;
  Creflnlist *poReflnlist;

  tagDenzoXStringRec tDXline;
  memset(&tDXline, char(0), sizeof(tDXline));

  if (3 > argc)
    {
      cout << "Usage: dtx2refln .x_file .ref_file [color_num_for_partials]\n"
           << flush;
      return (1);
    }
  int   nColor = 0;  // Partials will have color 0=blue, choose 5 for yellow
  if (3 < argc)
    nColor = atoi(argv[3]);

  poReflnlist = new Creflnlist();
  nFI_nPartialFlag  = poReflnlist->nExpandGetField("nPartialFlag");
  (void)poReflnlist->nExpandGetField(poReflnlist->ms_snNonunfFlag);
  nFI_fOtherInt     = poReflnlist->nExpandGetField("fOtherInt");
  nFI_fProfileChiSq = poReflnlist->nExpandGetField("fProfileChiSq");
  (void)poReflnlist->nExpandGetField(poReflnlist->ms_sfOblique);
  (void)poReflnlist->nExpandGetField(poReflnlist->ms_sfObsPx0);
  (void)poReflnlist->nExpandGetField(poReflnlist->ms_sfObsPx1);
  (void)poReflnlist->nExpandGetField(poReflnlist->ms_sfObsRotMid);
  nFI_fLPOfactor    = poReflnlist->nExpandGetField("fLPOfactor");
  nFI_fStrength     = poReflnlist->nExpandGetField("fStrength");

  Cstring sTemp, sTemp2;
  Cstring a15sTokens[15];
  float   fRotStart, fRotEnd, fRotMid;

  int     nH, nK, nL, nPartial;
  float   fInt, fSigmaI, fOtherInt, fProfileChiSq, fOblique;
  float   fObs_pixel0, fObs_pixel1, fObs_rotmid;
  float   fLPOfactor, fStrength;

  int nStat;
  int nReadIn = 0;
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

	  // Skip first 4 lines

	  for (i = 0; i < 4; i++)
	    getline(oIn, sTemp);

	  // 5th line has rot start and end

	  getline(oIn, sTemp);
	  split(sTemp, a15sTokens, 2, ' ');
	  fRotStart = (float) atof(a15sTokens[0].string());
	  fRotEnd   = (float) atof(a15sTokens[1].string());

	  fRotMid   = 0.5f * (fRotStart + fRotEnd);

	  // Skip 6th line

	  getline(oIn, sTemp);

	  bool bNotDone = TRUE;
	  while (bNotDone)
	    {
	      getline(oIn, sTemp);

	      // We have to do some special parsing of the refln line because
	      // the numbers sometimes run together without spaces between them

	      nStat = sscanf(sTemp.string(), "%4c%4c%4c%2c%8c%8c%7c"
                    	                     "%6c%6c%7c%7c%6c%8c",
			     &tDXline.ac_nH, &tDXline.ac_nK, &tDXline.ac_nL, 
			     &tDXline.ac_nPartial,
			     &tDXline.ac_fInt, &tDXline.ac_fOtherInt, 
			     &tDXline.ac_fProfileChiSq, &tDXline.ac_fSigmaI,
			     &tDXline.ac_fOblique,
			     &tDXline.ac_fObs_pixel1, &tDXline.ac_fObs_pixel0, 
			     &tDXline.ac_fLPOfactor, &tDXline.ac_fStrength);
	      if (13 == nStat)
		{
		  sscanf(tDXline.ac_nH,            "%4d", &nH);
		  sscanf(tDXline.ac_nK,            "%4d", &nK);
		  sscanf(tDXline.ac_nL,            "%4d", &nL);
		  sscanf(tDXline.ac_nPartial,      "%2d", &nPartial);
		  sscanf(tDXline.ac_fInt,          "%8f", &fInt);
		  sscanf(tDXline.ac_fOtherInt,     "%8f", &fOtherInt); 
		  sscanf(tDXline.ac_fProfileChiSq, "%7f", &fProfileChiSq);
		  sscanf(tDXline.ac_fSigmaI,       "%6f", &fSigmaI);
		  sscanf(tDXline.ac_fOblique,      "%6f", &fOblique);
		  sscanf(tDXline.ac_fObs_pixel1,   "%7f", &fObs_pixel1);
		  sscanf(tDXline.ac_fObs_pixel0,   "%7f", &fObs_pixel0);
		  sscanf(tDXline.ac_fLPOfactor,    "%6f", &fLPOfactor);
		  sscanf(tDXline.ac_fStrength,     "%8f", &fStrength);
		}
	      bNotDone = (999 != nH) && (13 == nStat);
	      if (bNotDone)
		{
		  nReadIn++;
		  oRefln.vSetH(nH);
		  oRefln.vSetK(nK);
		  oRefln.vSetL(nL);
		  oRefln.vSetField(nFI_nPartialFlag, nPartial);

		  // Try to affect color of refln when plotting with dtdisplay

		  if (1 == nPartial)
		    nPartial = nColor;  // It is a partial refln
		  oRefln.vSetField(poReflnlist->m_nFI_nNonunfFlag, nPartial);

		  oRefln.vSetIntensity(fInt);
		  oRefln.vSetSigmaI(fSigmaI);
		  oRefln.vSetField(nFI_fOtherInt,    fOtherInt);
		  oRefln.vSetField(nFI_fProfileChiSq, fProfileChiSq);
		  oRefln.vSetField(nFI_fLPOfactor, fLPOfactor);
		  oRefln.vSetField(nFI_fStrength, fStrength);
		  oRefln.vSetField(poReflnlist->m_nFI_fObsPx0, fObs_pixel0);
		  oRefln.vSetField(poReflnlist->m_nFI_fObsPx1, fObs_pixel1);
		  oRefln.vSetField(poReflnlist->m_nFI_fObsRotMid, fRotMid);
		  oRefln.vSetField(poReflnlist->m_nFI_fOblique, fOblique);
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

