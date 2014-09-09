//
// Copyright (c) 2006 Rigaku
// All rights reserved.

// dtrebatch.cc     Initial author: J.W. Pflugrath           21-Apr-2000
//    This reads a d*TREK dtprofit.ref reflnlist file and changes the batch
//    id's based on the observed rotation angle centroid
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
//    dtrebatch reads from the command line the reflnlist filename
//    if present (default: *dtprofit.ref) and options
//
//    Syntax:   % dtrebatch dtprofit.ref -dbatch 2
//    Example:  % dtrebatch dtprofit.ref -dbatch 2
//
//+ToDo
//
//   Better error messages need implementing.
//
//+Include files

//#include <iostream.h>
#include "Dtrek.h"
#include "dtrekdefs.h"
#include "Cstring.h"
#include "Crefln.h"
#include "Creflnlist.h"
#include "dskio.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

#ifdef SSI_PC
#define main   dtprofit_main
#define vError dtprofit_vError
#endif

void vError(const int nErrorNum=0, const Cstring& sError=NULL);

#ifdef SSI_PC
#define main dtrebatch_main
#endif

int main(int   argc,    // Command line argument count
         char *argv[])  // Pointers to command line args

{
  int      nStat;
  int      i, j;
  Cstring  sIn            = sDtrekGetPrefix() + "dtprofit.ref";
  Cstring  sOut           = "rb_" + sIn;
  Creflnlist *poReflnlist = NULL;
  Cstring        sTemp;
  Cstring        sMissingOption
                    = "ERROR: dtrebatch - missing option argument!";
  Cstring        sInvalidArgument
                    = "ERROR: dtrebatch - invalid argument!";

  float          m_fOverallRotStart;
  float          m_fOverallRotEnd;
  float          m_fRotIncrement;
  float          m_fDegreesPerBatch = 1.0f;
  float          fTemp;
  int            m_nImagesPerBatch;
  int            nTemp;
  int            m_nVerbose = 0;
  Cstring        m_sBatchPrefix = "";
  bool           m_bRebatch = TRUE;
  char           a255cNum[256];

  cout << "\ndtrebatch:  Copyright (c) 2006 Rigaku\n";
  cout << D_K_DTREKVersion << endl << flush;

  // Copy command line to output log

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;

  // Parse command line arguments

  argc--; argv++;

  nStat = 1;

  if (1 > argc)
    {
      // No reflection list arguments and thats an error!

      vError(3, "ERROR - no reflection list filename!\n");
    }

  sIn = (const char*) argv[0];

  // Read a reflection list from a file

  cout << "Reflection list: " << sIn << endl << flush;
  argc--;  argv++;

  if (1 <= argc)
    {
      sTemp = (const char*) argv[0];
      if ( ("-help" == sTemp) || ("-h" == sTemp) )
	{
	  vError(0, "");
	}
    }
  // Parse any additional command line options

  for (i = 0; i < argc; i++)
    {
      sTemp = (const char*) argv[i];
      cout << "Command line string: >>" << sTemp << "<<" << endl << flush;

      if ("-verbose" == sTemp)
	{
	  i++;
	if (i < argc)
	  {
	    nStat = sscanf(argv[i], "%d", &nTemp);
	    if (1 == nStat)
	      m_nVerbose = nTemp;
	    else
	      vError(6, sInvalidArgument);
	  }
	else
	  {
	    vError(5, sMissingOption);
	  }
	}
      else if ("-norebatch" == sTemp)
	{
	  m_bRebatch = FALSE;
	}
      else if ("-dbatch" == sTemp)
	{
	  i++;
	if (i < argc)
	  {
	    nStat = sscanf(argv[i], "%f", &fTemp);
	    if (1 == nStat)
	      m_fDegreesPerBatch = fTemp;
	    else
	      vError(6, sInvalidArgument);
	  }
	else
	  {
	    vError(5, sMissingOption);
	  }
	}
      else
	{
	  vError(2, "ERROR - Unknown command line option");
	}
    }

  // Read the reflnlist

  cout << "Reading " << sIn << " ..." << endl << flush;
  poReflnlist = new Creflnlist (sIn);

  if (!poReflnlist->bIsAvailable())
    {
      cout << "dtrebatch: reflnlist " << sIn << " unavailable!\n" << flush;
      delete poReflnlist;
#ifdef SSI_PC
      return 1;
#else
      exit (1);
#endif
    }

  if (1 > poReflnlist->nGetNumReflns())
    {
      cout << "dtrebatch: reflnlist " << sIn << " has no reflections!\n" << flush;
      delete poReflnlist;
#ifdef SSI_PC
      return 1;
#else
      exit (1);
#endif
    }

  if (0 > poReflnlist->m_nFI_fObsRotMid)
    {
      cout << "dtrebatch: reflnlist " << sIn << " has no fObsRotMid field!\n" << flush;
      delete poReflnlist;
#ifdef SSI_PC
      return 1;
#else
      exit (1);
#endif
    }

  if (0 > poReflnlist->m_nFI_sBatch)
    {
      cout << "dtrebatch: reflnlist " << sIn << " has no sBatch field!\n" << flush;
      delete poReflnlist;
#ifdef SSI_PC
      return 1;
#else
      exit (1);
#endif
    }

  sOut = sFileGetDirectory(sIn) + "rb_" + sFileGetBasename(sIn);

  // Sort reflnlist on ObsRotMid angles and re-label sBatch's

  poReflnlist->vSort(eReflnField_float_type, poReflnlist->m_nFI_fObsRotMid, NULL);

  float   fRotMin, fRotMax;
  Crefln *poRefln;
  int    *pnIndex;
  Cstring sBatchMin, sBatchMax;
  pnIndex   = poReflnlist->pnGetSortIndex();
  poRefln   = poReflnlist->poGetRefln(pnIndex[0]);
  fRotMin   = poRefln->fGetField(poReflnlist->m_nFI_fObsRotMid);
  sBatchMin = poRefln->sGetField(poReflnlist->m_nFI_sBatch);
  poRefln   = poReflnlist->poGetRefln(pnIndex[poReflnlist->nGetNumReflns()-1]);
  fRotMax   = poRefln->fGetField(poReflnlist->m_nFI_fObsRotMid);
  sBatchMax = poRefln->sGetField(poReflnlist->m_nFI_sBatch);


  nTemp   = (int)fRotMin;
  fRotMin = (float) nTemp;
  nTemp   = nint(fRotMax);
  fRotMax = (float) nTemp;

  if (-1 < m_nVerbose)
    {
      cout << "RotMin,Max: " << fRotMin << ", " << fRotMax << endl;
      cout << "sBatchMin,Max: " << sBatchMin << ", " << sBatchMax << endl;
    }

  // Divide reflections into scaling batches.  One scaling batch per image

  int nBatchMin;
  float fDegreesPerBatch = 1.0;

  // Default is 1 degree per scaling batch unless set

  if (0.0 < m_fDegreesPerBatch)
    {
      // User input the desired degrees per batch 

      fDegreesPerBatch = m_fDegreesPerBatch;
    }

  cout << "\nScaling batch ids will be every " << fDegreesPerBatch
       << " degree";
  if (1.0 < fDegreesPerBatch)
    cout << 's';
  cout << ".\n" << endl;

  m_fOverallRotStart = fRotMin;
  m_fOverallRotEnd   = fRotMax;
  fRotMax = m_fOverallRotStart + fDegreesPerBatch;

  int *pnNumRefsPerBatch;
  pnNumRefsPerBatch = new int [10000];

  sTemp = sBatchMin;
  if ("" != m_sBatchPrefix)
    sTemp = sBatchMin.after(m_sBatchPrefix);
  nBatchMin = atoi(sTemp.string());
  int nBatchCount      = 1;
  pnNumRefsPerBatch[0] = 0;              // Not pnNumRefsPerBatch[0] is unused
  pnNumRefsPerBatch[nBatchCount] = 0;
  for (i = 0; i < poReflnlist->nGetNumReflns(); i++)
    {
      poRefln   = poReflnlist->poGetRefln(pnIndex[i]);
      if (poRefln->fGetField(poReflnlist->m_nFI_fObsRotMid) > fRotMax)
	{
	  // Increment sBatchMin somehow

	  if (1000 > nBatchMin)
	    sprintf(a255cNum, "%s%03d", m_sBatchPrefix.string(), nBatchMin+nBatchCount);
	  else
	    sprintf(a255cNum, "%s%d", m_sBatchPrefix.string(), nBatchMin+nBatchCount);
	  if (3 < m_nVerbose)
	    {
	      cout << "Numrefs in batch: " << nBatchCount << ", " 
		   << pnNumRefsPerBatch[nBatchCount] << endl;
	    }
	  nBatchCount++;
	  if (10000 <= nBatchCount)
	    {
	      cout << "ERROR too many batches!\n" << flush;
#ifdef SSI_PC
	      return (3);
#else
	      exit (3);
#endif
	    }

	  sBatchMin = (Cstring)a255cNum;
	  fRotMax = m_fOverallRotStart + (float)nBatchCount * fDegreesPerBatch;
	  //	  cout << "New batchname: " << sBatchMin << " for rot < " 
	  //	       << fRotMax << endl;
	  pnNumRefsPerBatch[nBatchCount] = 0;
	}
      poRefln->vSetField(poReflnlist->m_nFI_sBatch, sBatchMin);
      pnNumRefsPerBatch[nBatchCount]++;
    }

  // There are nBatch's, so look if the first and last batch have 
  // too few reflns.  If so, rename them

  // Average the number of reflns in batches 2 through nBatchCount

  fTemp = 0.0;
  for (i = 2; i < nBatchCount; i++)
    {
      fTemp = fTemp + pnNumRefsPerBatch[i];
    }
  fTemp = fTemp / (nBatchCount - 2);
  cout // << "...done.\n\n"
       << "Average number of reflns per batch: " << nint(fTemp) << "\n\n";

  if ( (fTemp / 3.0) > pnNumRefsPerBatch[1])
    {
      poRefln  = poReflnlist->poGetRefln(pnIndex[0]);
      sTemp    = poRefln->sGetField(poReflnlist->m_nFI_sBatch);
      cout << "WARNING first batch (" << sTemp
	   << ") has too few reflections ( " << pnNumRefsPerBatch[1]
	   << " ),\n";
      if (m_bRebatch)
	{
	  // First batch is lite, so re-batch it.

	  cout << "  so first batch combined with next batch (";
	  poRefln   = poReflnlist->poGetRefln(pnIndex[pnNumRefsPerBatch[1]]);
	  sTemp     = poRefln->sGetField(poReflnlist->m_nFI_sBatch);
	  cout << sTemp << ")!" << endl;
	  for (i = 0; i < pnNumRefsPerBatch[1]; i++)
	    {
	      poRefln   = poReflnlist->poGetRefln(pnIndex[i]);
	      poRefln->vSetField(poReflnlist->m_nFI_sBatch, sTemp);
	    }
	}
      else
	{
	  cout << "  BUT NOT COMBINED with next batch!\n" << endl;
	}
    }

  if ( (fTemp / 3.0) > pnNumRefsPerBatch[nBatchCount])
    {
      // First batch is lite, so re-batch it.

      poRefln  = poReflnlist->poGetRefln(pnIndex[poReflnlist->nGetNumReflns()-1]);
      sTemp    = poRefln->sGetField(poReflnlist->m_nFI_sBatch);
      cout << "WARNING last batch (" << sTemp
           << ") has too few reflections ( " << pnNumRefsPerBatch[nBatchCount]
	   << " ),\n";
      if (m_bRebatch)
	{
	  cout << "  so last batch combined with previous batch (";
	  j = poReflnlist->nGetNumReflns() - 1 - pnNumRefsPerBatch[nBatchCount];
	  poRefln   = poReflnlist->poGetRefln(pnIndex[j]);
	  sTemp     = poRefln->sGetField(poReflnlist->m_nFI_sBatch);
	  cout << sTemp << ")!" << endl;
	  for (i = 0; i < pnNumRefsPerBatch[nBatchCount]; i++)
	    {
	      j = poReflnlist->nGetNumReflns() - 1 - i;
	      poRefln   = poReflnlist->poGetRefln(pnIndex[j]);
	      poRefln->vSetField(poReflnlist->m_nFI_sBatch, sTemp);
	    }
	}
      else
	{
	  cout << "  BUT NOT COMBINED with previous batch!\n" << endl;
	}
    }
  poRefln  = poReflnlist->poGetRefln(pnIndex[0]);
  sBatchMin= poRefln->sGetField(poReflnlist->m_nFI_sBatch);
  
  poRefln  = poReflnlist->poGetRefln(pnIndex[poReflnlist->nGetNumReflns()-1]);
  sBatchMax= poRefln->sGetField(poReflnlist->m_nFI_sBatch);

  cout << "First, last batch ids: " << sBatchMin << "   " << sBatchMax << endl;
  delete [] pnNumRefsPerBatch;

  nStat = poReflnlist->nWrite(sTransSymbol(sOut));
  if (0 == nStat)
    {
      cout << "...done ";
    }
  else
    {
      cout << "ERROR! ";
    }
  cout << "writing reflnlist file: " << sOut << '\n' << flush;
  delete poReflnlist;
  return (nStat);
}

void vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << '\n';
    }
  cout << "\ndtrebatch - Usage:\n"
       << "dtrebatch  ref_file [options]\n\n"
       << "Options:\n\n"
       << "-dbatch fDegPerBatch    Specify degrees per batch\n\n" << flush;
#ifndef SSI_PC
  exit (nErrorNum);
#endif

}
