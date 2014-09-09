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
// Cscalemerge.cc        Initial author: J.W. Pflugrath           14-Jan-1996
//  This file contains the member functions of class Cscalemerge
//    which implements the scaling and merging encapsulation of d*TREK.
//    The methods are directly from either the paper
//    W. Kabsch (1988) J. Appl. Cryst. 21, 916-924.  "Evaluation of single-
//       Crystal X-Ray Diffraction Data from a Position-Sensitive Detector"
//    -or-
//    Fox and Holmes (1966) Acta Cryst. 20, 886-891.
//    Hamilton, Rollett, Sparks (1965) Acta Cryst. 18, 129-?.
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
//+Description
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Cscalemerge.h"         // Class definition and prototypes

//+Definitions, constants, and initialization of static member variables

Cstring Cscalemerge::ms_sf2STLsq    = "f2STLsq";
Cstring Cscalemerge::ms_snBatchIdx  = "nBatchIndex";

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cscalemerge::Cscalemerge(Creflnlist    *poReflnlistIn,
			 Cimage_header *poHeaderIn)
{
 (void) nInitValues();

 m_poCrystal     = new Ccrystal(*poHeaderIn);
 if (poReflnlistIn->bIsAvailable() && m_poCrystal->bIsAvailable())
   {
     m_eThe_State = eScalemerge_available_state;

     // Create the internal reflection list with the same fields
     // as the input reflection list

     m_poReflnlist = new Creflnlist(*poReflnlistIn);

     // Expand the working reflection list to include required fields


//     m_nFI_fOtherInt = m_poReflnlist->nExpandGetField(D_K_ProfitfOtherInt);
//     m_nFI_fOtherSig = m_poReflnlist->nExpandGetField(D_K_ProfitfOtherSig);

     (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_ssBatch);

     m_nFI_f2STLsq   = m_poReflnlist->nExpandGetField(ms_sf2STLsq);
     m_nFI_nBatchIdx = m_poReflnlist->nExpandGetField(ms_snBatchIdx);
     m_poReflnlist->nReduce(*m_poCrystal);  // This adds lots of fields, too.

     // Now our internal list has all the fields it needs, insert the input
     // list into it
     // (Optimize strategy: if *poReflnlistIn and m_poReflnlist have the
     //  same fields, then just change pointer.  In this case, make sure
     // destructor does not delete what m_poReflnlist points to.)

     m_poReflnlist->nInsertListFrom(*poReflnlistIn);
     m_nFI_nAnomFlag
       = m_poReflnlist->nGetFieldIndex(m_poReflnlist->ms_snFplusminusFlag);
   }
}

Cscalemerge::~Cscalemerge()
{
//  cout << "Cscalemerge destructor called!\n";
  if (NULL != m_poCrystal)
    {
      delete m_poCrystal;
      m_poCrystal = NULL;
    }

  if (NULL != m_poReflnlist)
    {
      delete m_poReflnlist;
      m_poReflnlist = NULL;
    }

  if (NULL != m_poReflnlistAve)
    {
      delete m_poReflnlistAve;
      m_poReflnlistAve = NULL;
    }

  if (NULL != m_poReflnlistComp)
    {
      delete m_poReflnlistComp;
      m_poReflnlistComp = NULL;
    }

  if (NULL != m_psBatchNames)
    {
      delete [] m_psBatchNames;
      m_psBatchNames = NULL;
    }
  vDeleteArrays();
}

int
Cscalemerge::nInitValues(void)
{
  m_eThe_State          = eScalemerge_unknown_state;
  m_sScalemerge_key     = "unknown";
  m_sBatchFixed         = Creflnlist::ms_ssIsEmpty;
  m_psBatchNames        = NULL;
  m_fScaleFixed         = 1.0;
  m_fBValueFixed        = 0.0;
  m_nBatchFixed         = 0;
  m_poCrystal           = NULL;
  m_poReflnlist         = NULL;
  m_poReflnlistAve      = NULL;
  m_poReflnlistComp     = NULL;
  m_bNewList            = FALSE;
  m_bLastCycle          = FALSE;
  m_bSigmaOverEstimate  = FALSE;
  m_nAnomFlag           = 0;
  m_nScaleAnom          = 0;
  m_nFixBFlag           = 0;
  m_nMethod             = 2;
  m_nVerbose            = 3;
  m_nCountOverlap       = 0;
  m_nNumBatches         = 0;
  m_nFI_nAnomFlag       = -1;
  m_nFI_fIntensityPlus  = -1;
  m_nFI_fSigmaIPlus     = -1;
  m_nFI_fIntensityMinus = -1;
  m_nFI_fSigmaIMinus    = -1;
  m_pfAab               = NULL;
  m_pfBa                = NULL;
  m_pfIh                = NULL;
  m_pfRha               = NULL;
  m_pfUha               = NULL;
  m_pfVha               = NULL;
  m_pfUh                = NULL;
  m_pfVh                = NULL;
  m_pfGa                = NULL;
  m_pnGaFixed           = NULL;
  m_pfWhla              = NULL;
  m_pfW1                = NULL;
  m_pnW2                = NULL;
  m_pfW3                = NULL;
  m_pfW4                = NULL;
  m_pnW5                = NULL;
  m_pnIndex             = NULL;
  m_pnNumContrib        = NULL;
  m_pnNumInBatch        = NULL;
  m_pnNumAccepted       = NULL;
  m_pnNumRejected       = NULL;
  m_pnNumDeselected     = NULL;
  m_pnNumOverlap        = NULL;
  m_nMaxSymRelated      = 5000;
  m_pnNumMultsBatch     = NULL;

  m_ptStatsReso         = NULL;
  m_ptStatsInt          = NULL;
  m_ptStatsBatch        = NULL;
  m_nNumBinsReso        = 10;
  m_nNumBinsInt         = 10;
  m_a2fRangeReso[0]     = 0.0;
  m_a2fRangeReso[1]     = 0.0025f;
  m_a2fRangeInt[0]      = 0.0;
  m_a2fRangeInt[1]      = 20.0;

  m_fErrorMul           = 1.0;
  m_fErrorAdd           = 0.00;
  m_fEigRho             = 0.00004f;
  m_fRejCriteria        = (float)10.0E10;       // Should not have any rejected
  m_nMaxCycles          = 30;
  m_fResolutionMin      = 999999.0;
  m_fResolutionMax      = 0.0001f;

  return (0);
}

int
Cscalemerge::nList(const int nNum)
{
  cout << "Cscalemerge listing:"
       << "\nBatch  fixed:      " << m_sBatchFixed
       << "\nScale  fixed:      " << m_fScaleFixed
       << "\nBvalue fixed:      " << m_fBValueFixed
       << "\nRej criteria:      " << m_fRejCriteria
       << "\nMul sig fact:      " << m_fErrorMul
       << "\nAdd sig fact:      " << m_fErrorAdd
       << "\n  Max cycles:      " << m_nMaxCycles
       << "\n   FixB Flag:      " << m_nFixBFlag
       << "\nOutput Anom Flag:  " << m_nAnomFlag
       << "\nScale  Anom Flag:  " << m_nScaleAnom
       << "\nResolution:        " << m_fResolutionMin << " to "
       <<                            m_fResolutionMax << " Angstroms"
       << endl;
  if (NULL != m_poReflnlist)
    {
      cout << "\nNumber input reflns: " << m_poReflnlist->nGetNumReflns() << endl;
    }
  if (NULL != m_poCrystal)
    m_poCrystal->nList(1);
  return (0);
}


int
Cscalemerge::nScaleSetup(void)
{
  // Should have reflection list and a few other items, so ready
  // to compute scale factors for batches in the reflnlist

  int      i, j, k;              // Loop variables
  int      nStat;                // Local status
  int      nNumReflns;           // Temp variable holding number of reflns in list
  float    fTemp;                // Temp variable

  Crefln  *poRefln;              // A pointer used for readability of code
  Cstring  sBatchPrev, sBatchCurr; // Some strings to hold batch names
  int     *pnNumReflns;
  float   *pfMin[5], *pfMax[5];  // Min, max of 5 batch variables
  int      nFI[5];               // Indexes of 5 batch variables
  int      nNumBatches = 0;
  bool     bNoMatch;
  int      nBatchIndex;
  float    a3x3fBMat[3][3];      // Local var for crystal B matrix
  float    a3fX[3];              // Local var for recip lattice vector Bh = x
  float    a3fH[3];              // hkl of refln
  int      nAlloc;
  int      nAllocBunch = 100;

  if (!bIsAvailable())
    {
      cerr << "Cscalemerge object not available!\n";
      return (-1);  // Do not do anything is assertion fails
    }

  if (NULL == m_poCrystal)
    {
      cerr << "Crystal object is NULL!\n";
      return (-2);
    }
  if (!m_poCrystal->bIsAvailable())
    {
      cerr << "Crystal object not available!\n";
      return (-2);
    }

  // Initialize some variables

  nAlloc = nAllocBunch;
  for (i = 0; i < 5; i++)
    {
      pfMin[i]    = new float [nAlloc];
      pfMax[i]    = new float [nAlloc];
    }

  if (NULL != m_psBatchNames) delete [] m_psBatchNames;
  m_psBatchNames  = new Cstring [nAlloc];
  pnNumReflns     = new int [nAlloc];

  // Figure out which batch to scale against

  // Get the crystal B matrix into local variable

  m_poCrystal->nCalcBMatrix();
  m_poCrystal->vGetBMatrix(&a3x3fBMat[0][0]);

  // The following fields will be used to divide batches into subbatches
  // Maybe use a member variable and a method to specify these

  nFI[0] = m_poReflnlist->m_nFI_fObsPx0;
  if (0 > nFI[0])
    {
      cerr << "No observed pixel 0 field!\n";
      nFI[0] = 0;
//      return (-2);
    }
  nFI[1] = m_poReflnlist->m_nFI_fObsPx1;
  if (0 > nFI[1])
    {
      cerr << "No observed pixel 1 field!\n";
      nFI[1] = 0;
//      return (-2);
    }
  nFI[2] = m_poReflnlist->m_nFI_fObsRotMid;
  if (0 > nFI[2])
    {
      cerr << "No observed rot mid field!\n";
      nFI[2] = 0;
//      return (-2);
    }
  nFI[3]          = m_nFI_f2STLsq;
  nFI[4]          = m_poReflnlist->m_nFI_fIntensity;

  // Find number of different batches in the reflection list

  sBatchPrev   = "";
  nBatchIndex  = -1;
  nNumReflns   = m_poReflnlist->nGetNumReflns();

  float f2STLsq, f2STLsqMin, f2STLsqMax;
  f2STLsqMin = 0.0;
  if (0.0 < m_fResolutionMin)
    f2STLsqMin = 1.f/(m_fResolutionMin * m_fResolutionMin);
  f2STLsqMax = 999999.0;
  if (0.0 < m_fResolutionMin)
    f2STLsqMax = 1.f/(m_fResolutionMax * m_fResolutionMax);

  int *pnDeleteArray;
  pnDeleteArray = new int [nNumReflns];

  for (i = 0; i < nNumReflns; i++)
    {
      poRefln  = m_poReflnlist->poGetRefln(i);

      // Compute reciprocal lattice vector for potential scaling
      // May be a factor of 4 or 0.25 here?

      a3fH[0] = (float) poRefln->nGetH();
      a3fH[1] = (float) poRefln->nGetK();
      a3fH[2] = (float) poRefln->nGetL();
      vMulMat3DVec3D(a3x3fBMat, a3fH, a3fX);

      f2STLsq = fDot3D(a3fX, a3fX);
      poRefln->vSetField(m_nFI_f2STLsq, f2STLsq);
      if ( (f2STLsq > f2STLsqMax) || (f2STLsq < f2STLsqMin) )
	{
	  // Mark reflection for deletion?

	  pnDeleteArray[i] = 1;
	}
      else
	{
	  // Passes resolution limits, so work with it
	
	  pnDeleteArray[i] = 0;

	  sBatchCurr = poRefln->sGetField(m_poReflnlist->m_nFI_sBatch);
	  bNoMatch   = (sBatchPrev != sBatchCurr);
	  for (j = 0; (j < nNumBatches) && bNoMatch; j++)
	    {
	      // Check that batchname is in not list of previously found
	      // batch names
	      // We need to do this faster, maybe with some kind of hashing
	      // or other direct searching

	      bNoMatch = (m_psBatchNames[j] != sBatchCurr);
	      if (!bNoMatch)
		nBatchIndex = j;        // Save batch index
	    }
	  if (bNoMatch)
	    {
	      // Still no match, so add this batch name to the list of batches

	      nBatchIndex = nNumBatches;
	      nNumBatches++;
	      if (nNumBatches > nAlloc)
		{
		  // Not enough room in our arrays, so alloc new space

		  nAlloc = nAlloc + nAllocBunch;
		  float *pfTempMin[5], *pfTempMax[5];
		  for (k = 0; k < 5; k++)
		    {
		      pfTempMin[k]    = new float [nAlloc];
		      pfTempMax[k]    = new float [nAlloc];
		    }
		  int     *pnTemp;
		  Cstring *psTemp;
		  pnTemp  = new int     [nAlloc];
		  psTemp  = new Cstring [nAlloc];
		
		  // And Copy old info to new location

		  for (j = 0; j < nBatchIndex; j++)
		    {
		      pnTemp[j]    = pnNumReflns[j];
		      psTemp[j]    = m_psBatchNames[j];
		      for (k = 0; k < 5; k++)
			{
			  pfTempMin[k][j] =  pfMin[k][j];
			  pfTempMax[k][j] =  pfMax[k][j];
			}
		    }

		  // Delete old location

		  delete [] m_psBatchNames;
		  m_psBatchNames = psTemp;
		  delete [] pnNumReflns;
		  pnNumReflns    = pnTemp;
		  for (k = 0; k < 5; k++)
		    {
		      delete [] pfMin[k];
		      pfMin[k] = pfTempMin[k];
		      delete [] pfMax[k];
		      pfMax[k] = pfTempMax[k];
		    }
		}

	      // Initialize this batches info

	      for (k = 0; k < 5; k++)
		{
		  pfMin[k][nBatchIndex] = poRefln->fGetField(nFI[k]);
		  pfMax[k][nBatchIndex] = poRefln->fGetField(nFI[k]);
		}
	      m_psBatchNames[nBatchIndex] = sBatchCurr;
	      pnNumReflns[nBatchIndex]    = 0;
	    } // endif bNoMatch

	  sBatchPrev = sBatchCurr;

	  // Save batch index in refln

	  poRefln->vSetField(m_nFI_nBatchIdx, nBatchIndex);

	  // Save min and max of 5 refln fields for possible grouping

	  for (k = 0; k < 5; k++)
	    {
	      fTemp = poRefln->fGetField(nFI[k]);
	      pfMin[k][nBatchIndex] = min(pfMin[k][nBatchIndex], fTemp);
	      pfMax[k][nBatchIndex] = max(pfMax[k][nBatchIndex], fTemp);
	    }
	  pnNumReflns[nBatchIndex]++;
	}
    }  // i, end of looping through reflns

  int nNumDel =m_poReflnlist->nDelete(1, pnDeleteArray);
  cout << "\nReflections exceeding resolution limits: " << nNumDel
       << "\n                                   Kept: " << m_poReflnlist->nGetNumReflns()
       << " for scaling.\n" << flush;

  delete [] pnDeleteArray;

  if (0 >= m_poReflnlist->nGetNumReflns())
    {
      cout << "ERROR, no reflections pass resolution limits!\n";
#ifdef SSI_PC
      return (-1);
#else
      exit (-1);
#endif
    }
  if (1 < m_nVerbose)
    {
      if (1 < nNumBatches)
	{
	  printf("\nThere are %d different batches ", nNumBatches);
	}
      else
	{
	  printf("\nThere is only %d batch ", nNumBatches);
	}
      printf("in the input reflection list.\n");
      printf("\nObserved position limits of the Batches\n");
      char *pcShortLine = "--------------------------------------------------------------------------------\n";
      if ( (0 < nFI[0]) || (0 < nFI[1]) || (0 < nFI[2]) )
	{
	  // Print these out only if the fields exist in the input reflnlist

	  printf(pcShortLine);
	  printf("%6s%7s %22s%22s%22s\n", "Batch", "Num",
	     (m_poReflnlist->sGetFieldName(nFI[0],eReflnField_float_type)).string(),
	     (m_poReflnlist->sGetFieldName(nFI[1],eReflnField_float_type)).string(),
	     (m_poReflnlist->sGetFieldName(nFI[2],eReflnField_float_type)).string());
	  printf("%6s%7s %11s%11s%11s%11s%11s%11s\n", "name", "refs",
		 "Min", "Max", "Min", "Max", "Min", "Max");
	  printf(pcShortLine);
	  for (j = 0; j < nNumBatches; j++)
	    {
	      printf("%6s%7d ", m_psBatchNames[j].string(), pnNumReflns[j]);
	      for (k = 0; k < 3; k++)
		{
		  printf(" %10.1f %10.1f", pfMin[k][j], pfMax[k][j]);
		}
	      printf("\n");
	    }
	  printf(pcShortLine);
	}

      printf("\n\nIntensity and Resolution limits of the Batches\n");
      printf(pcShortLine);
      printf("%6s%7s %22s%22s%22s\n", "Batch", "Num", "Intensity ",
	     "Resolution ", "[2sinT/lam]^2 ");
      printf("%6s%7s %11s%11s%11s%11s%11s%11s\n", "name", "refs",
	     "Min", "Max", "Min", "Max", "Min", "Max");
      printf(pcShortLine);
      for (j = 0; j < nNumBatches; j++)
	{
//	  printf("%6s%7d ", (const char *)m_psBatchNames[j], pnNumReflns[j]);
	  printf("%6s%7d ", m_psBatchNames[j].string(), pnNumReflns[j]);
	  k = 4;
	  printf(" %10.1f %10.1f", pfMin[k][j], pfMax[k][j]);
	  k = 3;
	  printf(" %10.2f %10.2f", 1./sqrtf(pfMin[k][j]), 1./sqrtf(pfMax[k][j]));
	  printf(" %10.4f %10.4f", pfMin[k][j], pfMax[k][j]);
	  printf("\n");
	}
      printf(pcShortLine);
      printf("\n");
    }

  // Save min, max of reso and intensity for all batches

  m_a2fRangeReso[0] = pfMin[3][0];
  m_a2fRangeReso[1] = pfMax[3][0];
  m_a2fRangeInt[0]  = pfMin[4][0];
  m_a2fRangeInt[1]  = pfMax[4][0];
  for (j = 1; j < nNumBatches; j++)
    {
      m_a2fRangeReso[0] = min(pfMin[3][j], m_a2fRangeReso[0]);
      m_a2fRangeReso[1] = max(pfMax[3][j], m_a2fRangeReso[1]);
      m_a2fRangeInt[0]  = min(pfMin[4][j], m_a2fRangeReso[0]);
      m_a2fRangeInt[1]  = max(pfMax[4][j], m_a2fRangeReso[1]);
    }

  m_a2fRangeReso[0] = m_a2fRangeReso[0] * sqrtf(m_a2fRangeReso[0]);
  m_a2fRangeReso[1] = m_a2fRangeReso[1] * sqrtf(m_a2fRangeReso[1]);

  m_fSlopeReso = (m_a2fRangeReso[1] - m_a2fRangeReso[0])
                  / (float) m_nNumBinsReso;

  // Make intensity really Intensity / SigmaI
  m_a2fRangeInt[1] = 20.0;
  m_a2fRangeInt[0] = 0.0;
  m_fSlopeInt  = (m_a2fRangeInt[1] - m_a2fRangeInt[0])
                  / (float) m_nNumBinsInt;

  ////////////////////////////////////////////////////////////////////////
  // Reduce reflections to asymmetric unit and sort on nPackedHKL       //
  // so that symmetry-related reflections appear next to each other     //
  //                                                                    //
  // There should be a way to avoid doing this if the input list is     //
  // known to have been reduced and sorted already.                     //
  ////////////////////////////////////////////////////////////////////////

  // If keeping I+, I- separate for scaling, then place note in output

  if (0 == m_nScaleAnom)
    {
      cout << "Scaling assumes I+ equals I-.\n";
    }
  else
    {
      cout << "Scaling assumes I+ DOES NOT EQUAL I-,\n"
           << "   so anomalous information kept separated for scaling!\n";
    }

  cout << "Sorting and reducing reflnlist to asymmetric unit ..." << endl << flush;
  nStat = m_poReflnlist->nReduce(*m_poCrystal, m_nScaleAnom);
  if (0 != nStat)
    {
      cerr << "ERROR Cscalemerge: failure to reduce dataset!\n";
#ifdef SSI_PC
      return (-1);
#else
      exit (-1);
#endif
    }
  cout << "... done sorting" << endl << flush;

  // Allocate space for all arrays based on number of batches
  // The array names are from p. 923 of Kabsch's paper with a for alpha,

  // If we use the Kabsch method, then we need a better of way of
  // determining m_nNumSubBatchesPerBatch...(use pfMin, pfMax)

  m_nNumBatches = nNumBatches;
  m_nNumSubBatchesPerBatch = 1;
  if ( (1 == m_nMethod) || (2 == m_nMethod) )
    m_nNumSubBatchesPerBatch = 2;
  m_nNumParams = nNumBatches * m_nNumSubBatchesPerBatch;

  // Determine which batch to fix ...

  if (Creflnlist::ms_ssIsEmpty == m_sBatchFixed)
    {
      // No batch specified to fix, so use first one in the reflnlist

      poRefln = m_poReflnlist->poGetRefln(0);
      m_sBatchFixed = poRefln->sGetField(m_poReflnlist->m_nFI_sBatch);
      m_nBatchFixed = 0;
    }
  else
    {
      m_nBatchFixed = -1;
      for (i = 0; (i < m_nNumBatches) && (0 > m_nBatchFixed); i++)
	{
	  if (m_sBatchFixed == m_psBatchNames[i])
	    {
	      m_nBatchFixed = i;
	    }
	}
      if (0 > m_nBatchFixed)
	{
	  cerr << "Batch name " << m_sBatchFixed << " not found.\n";
	  m_nBatchFixed = 0;
	  m_sBatchFixed = m_psBatchNames[0];
	}
    }

  // Clean-up:  Delete new'd variables

  delete [] pnNumReflns;  // Maybe all these should be member variables?

  for (j = 0; j < 5; j++)
    {
      delete [] pfMin[j];
      delete [] pfMax[j];
    }

  return (nStat);
}

void
Cscalemerge::vDeleteArrays(void)
{
  if (NULL != m_pfAab)
    {
      delete [] m_pfAab;
      m_pfAab = NULL;
    }
  if (NULL != m_pfBa)
    {
      delete [] m_pfBa;
      m_pfBa = NULL;
    }
  if (NULL != m_pfIh)
    {
      delete [] m_pfIh;
      m_pfIh = NULL;
    }
  if (NULL != m_pfRha)
    {
      delete [] m_pfRha;
      m_pfRha = NULL;
    }
  if (NULL != m_pfUha)
    {
      delete [] m_pfUha;
      m_pfUha = NULL;
    }
  if (NULL != m_pfVha)
    {
      delete [] m_pfVha;
      m_pfVha = NULL;
    }
  if (NULL != m_pfUh)
    {
      delete [] m_pfUh;
      m_pfUh = NULL;
    }
  if (NULL != m_pfVh)
    {
      delete [] m_pfVh;
      m_pfVh = NULL;
    }
  if (NULL != m_pfGa)
    {
      delete [] m_pfGa;
      m_pfGa = NULL;
    }
  if (NULL != m_pnGaFixed)
    {
      delete [] m_pnGaFixed;
      m_pnGaFixed = NULL;
    }
  if (NULL != m_pfW1)
    {
      delete [] m_pfW1;
      m_pfW1 = NULL;
    }
  if (NULL != m_pnW2)
    {
      delete [] m_pnW2;
      m_pnW2 = NULL;
    }
  if (NULL != m_pfW3)
    {
      delete [] m_pfW3;
      m_pfW3 = NULL;
    }
  if (NULL != m_pfW4)
    {
      delete [] m_pfW4;
      m_pfW4 = NULL;
    }
  if (NULL != m_pnW5)
    {
      delete [] m_pnW5;
      m_pnW5 = NULL;
    }
  if (NULL != m_pfWhla)
    {
      delete [] m_pfWhla;
      m_pfWhla = NULL;
    }
  if (NULL != m_pnIndex)
    {
      delete [] m_pnIndex;
      m_pnIndex = NULL;
    }

  if (NULL != m_ptStatsReso)
    {
      delete [] m_ptStatsReso;
      m_ptStatsReso = NULL;
    }
  if (NULL != m_ptStatsInt)
    {
      delete [] m_ptStatsInt;
      m_ptStatsInt = NULL;
    }
  if (NULL != m_ptStatsBatch)
    {
      delete [] m_ptStatsBatch;
      m_ptStatsBatch = NULL;
    }
  if (NULL != m_pnNumMultsBatch)
    {
      delete [] m_pnNumMultsBatch;
      m_pnNumMultsBatch = NULL;
    }
}

int
Cscalemerge::nScale1(const int nMaxCyclesIn, Creflnlist *poReflnlistAve)
{
  // Perform scaling procedure according to the method specified by m_nMethod

  int     i, j, k, l;      // Loop counters
  int     nCyc;            // Iterative cycle number
  Crefln *poRefln;         // Handy pointer to a refln in reflnlist
  int     nNumRelated;     // Number of sym-related reflns to a given unique hkl
  int     nNumReflns;      // Number of reflns in reflnlist
  int     nPack;           // Packed reduced HKL value for a refln
                           //   all sym-related reflns have the same nPack
  int     nMaxCycles;
  int     nPackPrev;       // Packed reduced HKL value for prev refln
  int     nFirst;          // Index of index of first refln in a sym-rel group
  int     nLast;           // Index of index of last  refln in a sym-rel group
  bool    bSymRelated;     // Flag whether refln is sym-related to prev
  bool    bLastRefln;      // Flag for last refln
  int     nAccepted;       // Number of accepted refln in group of sym-related
  int     nRejected;       // Number of rejected refln in group of sym-related
  int     nTotAccepted;    // Total number of accepted refln
  int     nTotRejected;    // Total number of rejected refln
  int     nUnique;         // Total number of unique reflns going into refnement
  int     nIdx, nIdxPrev;  // Batch index
  int     a10nMult[10];
  float   fRejCriteria;    // Rejection criteria when averaging sym-related refs
  float   fSumResidSq;     // Sum of         Whl *(Ihl - <Ih>)**2      for all h
  float   fSumResid;       // Sum of              |Ihl - <Ih>|         for all h
  float   fResidSq;        //                Whl *(Ihl - <Ih>)**2  for a given h
  float   fResid;          // Sum of              |Ihl - <Ih>|     for a given h
  float   fSumAvgRootInt;
  float   fSumAvgInt;
  float   fSumAvgSig;

  float  *pfAab;           // Convenience pointer in m_pfAab

  // Along with certain member variables, the following are used in the
  // sums in normal equations.  The "Eqn #" in the comments reference equations
  // on p. 923 of Kabsch's paper.  The first equation on p. 923 is 1 with the
  // others numbered consecutively.

  float **ppfMat, *pfEigval, *pfOffDiag;  // Matrices, vectors for eigenvalue
  float  *pfAabDiag;                      // work
  float   fEigvalMax;
  float   fTemp;                          // A temp variable
  int    *pnTemp;
  float   fSigmaI;

  vDeleteArrays();  // Prevent memory leaks

  // Along with certain local variables, the following are used in the
  // sums in normal equations.
  // If you use the Kabsch method, then
  // the "Eqn #" in the comments reference equations
  // on p. 923 of Kabsch's paper.  The first equation on p. 923 is 1 with the
  // others numbered consecutively.  See vKabsch()

  m_pfAab     = new float [m_nNumParams * m_nNumParams]; // Square matrix
  m_pfBa      = new float [m_nNumParams];  // Right-hand side of norm eqn
  m_pfRha     = new float [m_nNumParams];  // Vha - (Ga * Uha * Ih)
  m_pfUha     = new float [m_nNumParams];  // Sum(Whla / sigmaIhl**2)
  m_pfVha     = new float [m_nNumParams];  // Sum(Whla * Ihl / sigmaIhl**2)
  m_pfUh      = new float [m_nNumParams];  // Sum(Ga * Ga * Uha)
  m_pfVh      = new float [m_nNumParams];  // Sum(Ga * Vha)
  m_pfGa      = new float [m_nNumParams];  // Scale factors at position a
  m_pnGaFixed = new int   [m_nNumParams];  // Flag for whether Ga is fixed
  m_pfWhla    = new float [m_nNumParams];  // Weight applied to observation l
                                              //   for unique h and position a
  m_pfW1      = new float [m_nMaxSymRelated];
  m_pnW2      = new int   [m_nMaxSymRelated];
  m_pfW3      = new float [m_nMaxSymRelated];
  m_pfW4      = new float [m_nMaxSymRelated];
  m_pnW5      = new int   [m_nMaxSymRelated];

  m_pnNumContrib    = new int [m_nNumParams];
  m_pnNumInBatch    = new int [m_nNumBatches+1];
  m_pnNumAccepted   = new int [m_nNumBatches+1];
  m_pnNumRejected   = new int [m_nNumBatches+1];
  m_pnNumDeselected = new int [m_nNumBatches+1];
  m_pnNumOverlap    = new int [m_nNumBatches * (m_nNumBatches+1)];

  m_pnNumMultsBatch = new int [m_nNumBatches+1];

  m_ptStatsReso     = new tagStats [m_nNumBinsReso+1];
  m_ptStatsInt      = new tagStats [m_nNumBinsInt+1];
  m_ptStatsBatch    = new tagStats [m_nNumBatches+1];

  ppfMat            = new float* [m_nNumParams+1];
  pfEigval          = new float [m_nNumParams];
  pfOffDiag         = new float [m_nNumParams];
  pfAabDiag         = new float [m_nNumParams];
  ppfMat[0] = NULL;
  for (i = 0; i < m_nNumParams; i++)
    {
      ppfMat[i+1] = new float [m_nNumParams+1];
    }

  // Preliminaries before cycling

  Crefln  *poReflnAve;
  poReflnAve    = NULL;

  bool bNewList = FALSE;
  if (NULL == poReflnlistAve)
    {
      poReflnlistAve = new Creflnlist;
      bNewList       = TRUE;
    }

  if (0 != poReflnlistAve->nGetNumReflns())
    {
      cerr << "ERROR, output averaged reflnlist must start off empty!\n";
#ifdef SSI_PC
      return (1);
#else
      exit (1);
#endif
    }

  if (0 != m_nAnomFlag)
    {
      // Add required fields to poReflnlistAve for anomalous dispersion info

      m_nFI_fIntensityPlus  = poReflnlistAve->nExpandGetField("fIntensity+");
      m_nFI_fSigmaIPlus     = poReflnlistAve->nExpandGetField("fSigmaI+");
      m_nFI_fIntensityMinus = poReflnlistAve->nExpandGetField("fIntensity-");
      m_nFI_fSigmaIMinus    = poReflnlistAve->nExpandGetField("fSigmaI-");
    }

  poReflnAve = new Crefln(poReflnlistAve);

  for (i = 0; i < m_nNumParams; i++)
    {
      m_pfGa[i]      = 1.0;      // Start with all parameters = 1.0
      m_pnGaFixed[i] = 0;        // and all not fixed
    }

  if ( (1 == m_nMethod) || (2 == m_nMethod) )
    {
      // Set all B factors to 0.0
      for (i = m_nNumBatches; i < m_nNumParams; i++)
        {
          m_pfGa[i] = 0.0;

	  // Fix all B factors if requested
	  if (1 == m_nFixBFlag) m_pnGaFixed[i] = 1;
        }

      // Anything else fixed?

      for (i = 0; i < m_nNumBatches; i++)
	{
	  if (m_psBatchNames[i] == m_sBatchFixed)
	    {
//	      m_pnGaFixed[i]                = 1;  // Fix ki, no, let eigval filt
//	      m_pnGaFixed[i+m_nNumBatches]  = 1;  // Fix Bi   take care of it
//	      m_pfGa[i]                     = m_fScaleFixed;  // May want these
//	      m_pfGa[i+m_nNumBatches]       = m_fBValueFixed; // applied later!
	      m_pfGa[i]                     = 1.;              //??
	      m_pfGa[i+m_nNumBatches]       = 0.;              //??
	    }
	}
    }

  nNumReflns = m_poReflnlist->nGetNumReflns();

  m_pnIndex = new int [nNumReflns];

  m_poReflnlist->vSort(eReflnField_int_type, m_poReflnlist->m_nFI_nPackedHKL,
		       m_pnIndex);

  bool bConverged = FALSE;

  nMaxCycles   = nMaxCyclesIn;
  m_bLastCycle = FALSE;

  if (0 > nMaxCycles) nMaxCycles = m_nMaxCycles;

  for (nCyc = 1; (nCyc <= nMaxCycles) && (!bConverged); nCyc++)
    {
      //
      if (3 < m_nVerbose)
	cout << "\nStarting cycle: " << nCyc << "\n" << flush;

      if (nCyc == nMaxCycles)
	{
	  // Re-Select reflections with I/sigmaI that were deselected
	  // so that the summary tables use all data!

	  (void) m_poReflnlist->nSelect((Cstring)"+fIntensity/fSigmaI>-1000.0");
	
	  // Set last cycle flag to TRUE

	  m_bLastCycle = TRUE;
	}

      if ( (4 < nMaxCycles) && (nCyc <= 2) )
	{
	  // On first two cycles of extended cycle loops, do not reject any reflns
	  fRejCriteria = 1.0E10;
	}
      else
	{
	  fRejCriteria = m_fRejCriteria;
	}

      // Zero all sums

      pnTemp = m_pnNumOverlap;
      for (i = 0; i < m_nNumBatches+1; i++)
	{
	  m_pnNumInBatch[i]    = 0;
	  m_pnNumAccepted[i]   = 0;
	  m_pnNumRejected[i]   = 0;
	  m_pnNumDeselected[i] = 0;
	  for (j = 0; j < m_nNumBatches; j++)
	    {
	      *pnTemp++ = 0;
	    }
	}

      pfAab = m_pfAab;
      for (i = 0; i < m_nNumParams; i++)
	{
	  m_pnNumContrib[i]    = 0;
	  m_pfBa[i]            = 0.0;
	  m_pfRha[i]           = 0.0;
	  for (j = 0; j < m_nNumParams; j++)
	    {
	      *pfAab++ = 0.0;       // Zero normal matrix via temporary ptr
	    }
	}

      for (i = 0; i < m_nNumBinsInt+1; i++)
	{
	  vInitStats(&m_ptStatsInt[i]);
	}
      for (i = 0; i < m_nNumBinsReso+1; i++)
	{
	  vInitStats(&m_ptStatsReso[i]);
	}
      for (i = 0; i < m_nNumBatches+1; i++)
	{
	  vInitStats(&m_ptStatsBatch[i]);
	}

      nTotAccepted   = 0;
      nTotRejected   = 0;
      nUnique        = 0;
      fSumResid      = 0.0;
      fSumResidSq    = 0.0;
      fSumAvgInt     = 0.0;
      fSumAvgSig     = 0.0;
      fSumAvgRootInt = 0.0;

      for (i =0; i < 10; i++)
	a10nMult[i] = 0;

      // Make sure last refln is selected since required by loop logic below.

      j = nNumReflns - 1;
      while ( (j >= 0) &&
	      !m_poReflnlist->poGetRefln(m_pnIndex[j])->bIsSelected())
	{
	  poRefln = m_poReflnlist->poGetRefln(m_pnIndex[j]);
	  nIdx    = poRefln->nGetField(m_nFI_nBatchIdx);
	  m_pnNumDeselected[nIdx]++;
	  m_pnNumInBatch[nIdx]++;
	  j--;
	}
      nNumReflns = j + 1;

      // Loop through all reflections, but skip those not selected

      nPackPrev   = -1;
      nIdxPrev    = -1;
      nNumRelated = 0;
      nFirst      = -1;
      for (j = 0; (j < nNumReflns); j++)
	{
	  poRefln = m_poReflnlist->poGetRefln(m_pnIndex[j]);
	  nIdx    = poRefln->nGetField(m_nFI_nBatchIdx);
	  m_pnNumInBatch[nIdx]++;
	  nAccepted = 0;
	  if (!poRefln->bIsSelected())
	    {
	      m_pnNumDeselected[nIdx]++;      // Count deselected reflns
	    }
	  else
	    {
	      // Look at only selected reflections.
	      // First check if symmetry-related to previous reflection
	      // by looking at the nPackedHKL field
	
	      nPack   = poRefln->nGetField(m_poReflnlist->m_nFI_nPackedHKL);
	      fSigmaI = poRefln->fGetSigmaI();
	      if (0.0 > fSigmaI)
		{
		  // Restore reflns rejected by prev cycles
		
		  fSigmaI  = -fSigmaI;
		  poRefln->vSetSigmaI(fSigmaI);
		}

	      bSymRelated = (nPack == nPackPrev);
	      bLastRefln  = (j == nNumReflns-1);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
	      if (bSymRelated)
		{
		  nNumRelated++;     // Increment num of sym-related in group
		  nLast = j;         // Keep index of last refln in group
		}
	      if (!bSymRelated || (bSymRelated && bLastRefln) )
		{
		  // This reflection is not symmetry-related to previous
		  // (or is the last refln and is symmetry-related!)
		  // So if there is a previous group of symmetry-related
		  // reflections, then process them
		  // BUT we need to have at least 2 sym-related reflns or
		  // what is the point?

		  if (1 == nNumRelated)
		    {
		      // Single previous refln needs processing

		      fResidSq = fAverageSymRelated(0, nLast - nFirst + 1,
						    &m_pnIndex[nFirst],
						    fRejCriteria,
						    poReflnAve,
						    &nAccepted, &nRejected,
						    &fResid);
		      if ( (m_bLastCycle) && (0 < nAccepted) )
			{
			  poReflnlistAve->nInsert(poReflnAve);
			  a10nMult[min(9, nAccepted)]++;
			}
		    }
		  else if (1 < nNumRelated)
		    {
		      // Ok, process this group of sym-related reflns
		      // with indices in m_pnIndex from nFirst thru nLast
		      // Remember: there may be some intervening deselected refs
		      // so we cannot use nNumRelated

		      // We could make the first few cycles faster by
		      // not calling fAverageSymRelated and setting
		      // nAccepted to nNumRelated.  fResid would be wrong,
		      // but ...

		      fResidSq = fAverageSymRelated(0, nLast - nFirst + 1,
						    &m_pnIndex[nFirst],
						    fRejCriteria,
						    poReflnAve,
						    &nAccepted, &nRejected,
						    &fResid);
		      if ( (m_bLastCycle) && (0 < nAccepted) )
			{
			  poReflnlistAve->nInsert(poReflnAve);
			  a10nMult[min(9, nAccepted)]++;
			}

		      nTotRejected = nTotRejected + nRejected;

		      // Danger here, if nAccepted > 1, but the reflns
		      // belong to the same scaling batch, then there
		      // is no new information for the vKandB routine.

		      if (1 < nAccepted)
			{
			  nTotAccepted = nTotAccepted + nAccepted;
			  nUnique++;
			  fSumResidSq = fSumResidSq + fResidSq;
			  fSumResid   = fSumResid   + fResid;
			  fSumAvgInt  = fSumAvgInt
			    + poReflnAve->fGetIntensity() * (float) nAccepted;
			  fSumAvgSig  = fSumAvgSig
			    + poReflnAve->fGetSigmaI();
			  fSumAvgRootInt = fSumAvgRootInt
			    + poReflnAve->fGetIntensity() * (float) nAccepted;
			  a10nMult[min(9, nAccepted)]++;

			  // Add contribution of this group to normal matrix

			  if (!m_bLastCycle)
			    {
			      if (0 == m_nMethod)
				vKabsch(nLast-nFirst+1, &m_pnIndex[nFirst]);
			      else if (1 == m_nMethod)
				vKandB(nLast-nFirst+1, &m_pnIndex[nFirst]);
			      else if (2 == m_nMethod)
				vKandB2(nLast-nFirst+1, &m_pnIndex[nFirst]);
			    }
			}
		    }

		  // Initialize items for this new group of sym-related reflns
	
		  nPackPrev   = nPack;    // Keep pack value for this group
		  nIdxPrev    = nIdx;
		  nNumRelated = 1;        // We have 1 reflection in the group
		  nFirst      = j;        // Keep starting index of this group
		  nLast       = j;        // Keep ending   index of this group
		}

	      if (!bSymRelated && bLastRefln)
		{
		  // Process last singlet refln

		  fResidSq = fAverageSymRelated(0, nLast - nFirst + 1,
						&m_pnIndex[nFirst],
						fRejCriteria,
						poReflnAve,
						&nAccepted, &nRejected,
						&fResid);
		  if ( (m_bLastCycle) && (0 < nAccepted) )
		    {
		      poReflnlistAve->nInsert(poReflnAve);
		      a10nMult[min(9, nAccepted)]++;
		    }
		}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
	    } // endif refln selected
	} // end j loop over reflns

      int nNumRefined;
      if (m_bLastCycle)
	{
	  nNumRefined = 0;
	}
      else
	{
	  // Remove all rows and columns of fixed parameters

	  for (i = 0; i < m_nNumParams; i++)
	    {
	      // Fix any parameter which has 0 on the diagonal of the normal
	      // equations so that there is no divide by 0 below.
	
	      k = (m_nNumParams + 1) * i;
	      if ( (0.0 == m_pfAab[k]) && (0 == m_pnGaFixed[i]) )
		{
		  m_pnGaFixed[i] = 1;
		  cout << "WARNING: parameter " << i << " was changed to FIXED!\n";
		}
	    }

#ifdef SSI_PC
     cout << flush;
#endif
	
	  if (6 < m_nVerbose)
	    {
	      // Print out the matrix
	
	      cout << "m_pfAab before compress: \n";
	      vListMatMN(m_nNumParams, m_nNumParams, m_pfAab);
	      cout << "m_pfBa: \n" << flush;
	      vListMatMN(m_nNumParams, 1, m_pfBa);
	    }

	  // Remove rows and columns with fixed parameters
	
	  nNumRefined = nMatrixCompress(m_nNumParams, m_pnGaFixed,
					    m_pfAab, m_pfBa);
	  if (5 < m_nVerbose)
	    {
	      // Print out the matrix
	
	      cout << "m_pfAab after compress: \n";
	      vListMatMN(nNumRefined, nNumRefined, m_pfAab);
	      cout << "m_pfBa: \n";
	      vListMatMN(nNumRefined, 1, m_pfBa);
#ifdef SSI_PC
         cout << flush;
#endif
	    }
	}

      // We probably need a warning about number of parameters versus
      // number of observations here

      if (m_bLastCycle)
	{
	  // Do not refine anything on the last cycle

	  cout << "Last cycle: no shifts.\n" << flush;
	  nNumRefined = 0;
	}
      else if (0 == nNumRefined)
	{
	  // Nothing to refine, that's OK, we just calculate ChiSq and
	  // report
	    cout << "WARNING: Nothing to refine!\n" << flush;
	    nMaxCycles = nCyc + 1;  // Do just one more cycle for statistics.
	}
      else if (1 == nNumRefined)
	{
	  // Mmmm, only 1 refined parameter, this means we do not have
	  // a matrix, so this can be made easy
	
	  if (0.0 != m_pfAab[0])
	    {
	      // Compute the solution

	      pfOffDiag[0] = m_pfBa[0] / m_pfAab[0];
	      pfAabDiag[0] = 1.0;                      // Shift scale factor
	    }
	}
      else
	{

	  // Now solve normal equations for arrays m_pfAab, m_pfGa and m_pfBa:
	  //
	  //     Sum Aab delGb   = Ba         (Eqn 2)
	  //         === -----     --
	
	  // First apply some Eigenvalue filtering as describe by G. Bricogne
	  // in the Proceedings ... Vol III, pp. 67-68.
	  // See also Diamond, R. (1966) Acta Cryst. 21, 253-266.
	
	  // Scale the vector m_pfBa, by 1/sqrtf(A[a][a])
	  //                               (diagonal of Aab)**(-1/2)
	  // Scale the matrix m_pfAab, so that there are 1's along diagonal
	  // And place result in ppfMat for call to vTRED2.

	  for (i = 0; i < nNumRefined; i++)
	    {
	      k = (nNumRefined + 1) * i;
	      for (j = i; j < nNumRefined; j++)
		{
		  l = i * nNumRefined  +  j;     // l is element [i][j]

		  // Note i+1, j+1, not i,j as per kludge of Numerical Recipes:
		
                 ppfMat[i+1][j+1] = m_pfAab[l];  // place m_pfAab into ppfMat
                 ppfMat[j+1][i+1] = m_pfAab[l];  // place m_pfAab into ppfMat
		} // end j
	    } // end i

	  if (4 < m_nVerbose)
	    {
	      cout << "ppfMat before vTRED: \n" << flush;
	      for (i = 1; i <= nNumRefined; i++)
		vListMatMN(nNumRefined, 1, &ppfMat[i][1]);
	    }

	  // When passing these arrays, use kludge in Numerical Recipes...

	  vTRED2( nNumRefined, ppfMat, pfEigval-1, pfOffDiag-1);
	  vTQLI(  nNumRefined, pfEigval-1, pfOffDiag-1, ppfMat);
	  vEigsrt(nNumRefined, pfEigval-1, ppfMat);  //Sort in DESCENDING order!

	  if (4 < m_nVerbose)
	    {
	      cout << "ppfMat after vTRED: \n" << flush;
	      for (i = 1; i <= nNumRefined; i++)
		vListMatMN(nNumRefined, 1, &ppfMat[i][1]);
	    }

	  if (3 < m_nVerbose)
	    {
	      cout << "Eigenvalues are: \n" << flush;
	      vListMatMN(nNumRefined, 1, pfEigval);
	    }
	  // Save maximum eigenvalue for filtering

	  fEigvalMax = pfEigval[0];
	
	  // Form "filtered inverse" of Aab into m_pfAab
	  // Put shift vector in pfOffDiag
	  // Watch out here for difference between m_nNumParams and nNumRefined!

	  vZeroMat(1, m_nNumParams, pfOffDiag);          // Zero shift vector
	  vZeroMat(m_nNumParams, m_nNumParams, m_pfAab); // Zero filtered
                                                         //  inverse
	  pfAab = m_pfAab;
	  for (k = 0; k < nNumRefined; k++)          // Loop over columns
	    {
	      for (j = 0; j < nNumRefined; j++)      // loop over rows
		{
		  for (i = 0; i < nNumRefined; i++)  // loop over eigenvalues
		    {
		      if (pfEigval[i] >= fEigvalMax * m_fEigRho)
			{
			  // Include these eigenvectors in superposition on way
			  // to filtered inverse                T
			  // Sum over i ( 1/eigval[i]) * Ei * Ei )
			  //                             --   --  (Nx1 * 1xN)
			  //                                             = NxN
	  // m_pfAab[j][k] = SUM(i)[eigvec[j][i]*eigvec[k][i]/eigval[i]]

 			  *pfAab = *pfAab + ppfMat[k+1][i+1] * ppfMat[j+1][i+1]
			               / pfEigval[i];
			}
		    }      // end i
		  pfAab++;
		}  // end j
	    }  // end k
	
	  // Now have the matrix in m_pfAab to multiply m_pfBa by to get
	  // shifts
	
	  if (4 < m_nVerbose)
	    {
	      cout << "after superposition: m_pfAab: \n" << flush;
	      vListMatMN(nNumRefined, nNumRefined, m_pfAab);
	    }

	  vMulMatNDVecND(nNumRefined, m_pfAab, m_pfBa, pfOffDiag);
	}

      // OK, adjust pfOffDiag (which holds shifts) in case there were some
      // fixed parameters that were removed by the matrix compression above

      if (0 < nNumRefined)
	{
	  j = nNumRefined-1;
	  for (i = m_nNumParams-1; i >= 0; i--)
	    {
	      if (0 == m_pnGaFixed[i])
		{
		  // Was refined, so shift available
		
		  pfOffDiag[i] = pfOffDiag[j];
		  j--;
		}
	      else
		{
		  pfOffDiag[i] = 0.0; 	      // Was fixed, so shift is 0.0
		}
	    }
	}
      else
	{
	  for (i = 0; i < m_nNumParams; i++)
	    {
	      pfOffDiag[i] = 0.0;
	    }
	}

      if (3 < m_nVerbose)
	{
	  cout << "Shifts: pfOffDiag: \n" << flush;
	  vListMatMN(m_nNumParams, 1, pfOffDiag);
	}

      // Fix up number of contributors to the B factors
      // This m_pnNumContrib treatment throughout the program appears
      // kludgy or kludgey.
      // Fix shifts in Bfactors so that the fixed batch does not shift.

      if ( (1 == m_nMethod) || (2 == m_nMethod) )
	{
	  fTemp = pfOffDiag[m_nBatchFixed+m_nNumBatches];
	  for (i = m_nNumBatches; i < m_nNumParams; i++)
	    {
	      m_pnNumContrib[i] = m_pnNumContrib[i-m_nNumBatches];
	      pfOffDiag[i] = pfOffDiag[i] - fTemp;
	    }
	}

      fTemp = 0.0;
      for (i = 0; i < m_nNumParams; i++)
	{
	  // Apply shift to scale factors,
	  //   and sum new scale factors so they can be normalized
	
	  m_pfVha[i] = m_pfGa[i];                // Save old scale factor
	  if (0 < m_pnNumContrib[i])
	    {
	      // Add shifts, but only if there were contributors

	      m_pfGa[i]  = m_pfGa[i] + pfOffDiag[i];
	    }
	  fTemp      = fTemp + m_pfGa[i];        // Sum(Ga[i]) for normalizing
	}
      if (3 < m_nVerbose)
	{
	  cout << "m_pfGa: \n" << flush;
	  vListMatMN(m_nNumParams, 1, m_pfGa);
	}

      if ( (1 == m_nMethod) || (2 == m_nMethod) )
	{
	  // Rescale, so that fixed scale factor is constant = m_fScaleFixed

	  fTemp = m_pfGa[m_nBatchFixed];   //  * m_fScaleFixed;
	  for (i = 0; i < m_nNumBatches; i++)
	    {
	      m_pfGa[i] = m_pfGa[i] / fTemp;
	    }
	}

      if (3 < m_nVerbose)
	{
	  cout << "after rescale: m_pfGa: \n" << flush;
	  vListMatMN(m_nNumParams, 1, m_pfGa);
	}

      bConverged = TRUE;                           // Be optimistic!
//      fTemp = fTemp / m_nNumParams;
      for (i = 0; i < m_nNumParams; i++)
        {
//          m_pfGa[i] = m_pfGa[i] / fTemp;           // Normalize
          if (fabs(m_pfGa[i] - m_pfVha[i]) > m_fEigRho)
            {
              bConverged = FALSE;                  // We were too optimistic
            }
        }

      if (nCyc <= nMaxCycles/2) bConverged = FALSE;   // Force a few cycles

      if ( (3 < m_nVerbose) || (m_bLastCycle) )
	{
	  // Print some end-of-cycle statistics

	  printf("For cycle number %d\n\n", nCyc);
//	  printf("         Normalized  ChiSquared: %10.3f\n",
//		 fSumResidSq / (float)nUnique
//		             * sqrtf((float) nTotAccepted / (float) nUnique)
// **		    fSumResidSq / (float)nUnique * (nTotAccepted / nUnique)
// **               fSumResidSq / (float)nTotAccepted
// **               fSumResidSq / (float)nUnique
//		 );
	  printf("   Method 1     Expected Rmerge: %10.3f\n",
		 1./sqrtf(fSumAvgInt/nTotAccepted)
		 );
	  printf("   Method 2     Expected Rmerge: %10.3f\n",
		 (fSumAvgSig/nUnique) / (fSumAvgInt/nTotAccepted)
//		     * sqrtf(nUnique / nTotAccepted)
		 );
	  printf("                Actual   Rmerge: %10.3f\n",
		 fSumResid / fSumAvgInt
		 );
	
	  char     *pcLine =
   "--------------------------------------------------------------------------------\n";

	  // Print out multiplicity statistics

	  printf("\n\nMultiplicity of observed reflections\n");
	  printf(pcLine);
	  printf("%8s %7d %7d %7d %7d %7d %7d %7d %7d %7s\n", "Mult | ",
		 (int)1, (int)2, (int)3, (int)4, (int)5, (int)6,
		 (int)7, (int)8, ">8");
	  printf(pcLine);
	  printf("%8s", " Refs | ");
	  for (i = 1; i < 10; i++)
	    {
	      printf(" %7d", a10nMult[i]);
	    }
	  printf("\n%s", pcLine);
	  printf("*Reflections with a multiplicity of 1 are not used in\n"
		 " scale factor refinement nor in Rmerge calculations.\n\n");

	  // Print out number of overlapping reflns in each batch

	  if (0 != m_nCountOverlap)
	    vListOverlap();

	  // Print out number of reflns in each category for each batch

	  printf("\n\nReflections in input file\n");
	  printf(pcLine);
	  printf("%14s %9s %9s %9s %9s %9s\n",
		 "Batch", " Num", "  Num", " Num", "  Num", "    Num");
	  printf("%14s %9s %9s %9s %9s %9s\n",
		 " name", "refs", "excluded", "rejs", "ovlps", "singles");
	  printf(pcLine);

	  int nSumSingles = 0;
	  for (i = 0; i < m_nNumBatches; i++)
	    {
	      m_pnNumInBatch[m_nNumBatches]
		= m_pnNumInBatch[m_nNumBatches] + m_pnNumInBatch[i];
	      m_pnNumDeselected[m_nNumBatches]
		 = m_pnNumDeselected[m_nNumBatches] + m_pnNumDeselected[i];
	      m_pnNumRejected[m_nNumBatches]
		 = m_pnNumRejected[m_nNumBatches] + m_pnNumRejected[i];
	      m_pnNumAccepted[m_nNumBatches]
		 = m_pnNumAccepted[m_nNumBatches] + m_pnNumAccepted[i];
	      nSumSingles = nSumSingles + m_ptStatsBatch[i].nNumSingles;
	      printf("%14s %9d %9d %9d %9d %9d\n",
//		     (const char *)m_psBatchNames[i],
		     m_psBatchNames[i].string(),
		     m_pnNumInBatch[i], m_pnNumDeselected[i],
		     m_ptStatsBatch[i].nNumRejects, m_pnNumAccepted[i],
		     m_ptStatsBatch[i].nNumSingles);
	    }
	  i = m_nNumBatches;
	  printf(pcLine);
	  printf("%14s %9d %9d %9d %9d %9d\n",
		     "All batches",
		     m_pnNumInBatch[i], m_pnNumDeselected[i],
		     m_pnNumRejected[i], m_pnNumAccepted[i],
		     nSumSingles);                         // Kludge watch out!


	  if ( (3 < m_nVerbose) || ((1 < nMaxCycles) && (1 < m_nNumBatches)) )
	    {
	      // Print out scale factors and number of reflns contributing
	      // to those scale factors, but NOT if not refined (i.e. cycles<2)

	      printf("\n\nRefined scale factors\n");
	      printf(pcLine);
	      printf("%14s %9s %9s\n",
		     "Batch", " Num", "Scale");
	      printf("%14s %9s %9s %9s ",
		     " name", "ovlps", "  K  ", "*Shifts");
	      if ( (1 == m_nMethod) || (2 == m_nMethod) )
		{
		  printf ("%9s %9s\n", "B  ", "*Shifts");
		}
	      else
		{
		  printf("\n");
		}
	      printf(pcLine);
	      for (i = 0; i < m_nNumBatches; i++)
		{
		  printf("%14s %9d %9.4f %9.4f ",
//		     (const char *)m_psBatchNames[i], m_pnNumContrib[i],
		     m_psBatchNames[i].string(), m_pnNumContrib[i],
		     m_pfGa[i], pfOffDiag[i]);

		  // If KandB print out B factor and B factor shifts

		  if ( (1 == m_nMethod) || (2 == m_nMethod) )
		    {
		      printf("%9.4f %9.4f\n",
			     m_pfGa[i+m_nNumBatches], pfOffDiag[i+m_nNumBatches]);
		    }
		  else
		    {
		      printf("\n");
		    }
		}
	      printf(pcLine);
	      printf("*Shifts are for previous cycle only!\n");
	    }
	  else
	    {
	      printf("\n***Scale factors not refined by this program"
                     "\n   since number of cycles <= 1 OR number of batches <= 1.***\n");
	    }

	  vListStats();

//	  poReflnlistAve->nAddResol(*m_poCrystal);
//	  poReflnlistAve->nReduce(*m_poCrystal);
//	  nListCompleteness(poReflnlistAve);
	  nListCompleteness();
	}

      if (bConverged && (nCyc < nMaxCycles))
	{
	  nMaxCycles = nCyc + 1;  // Do just 1 more cycle to calculate results
	  bConverged = FALSE;     //  without any more shifts
	  cout << "Converged before max cycles, so one more for final statistics.\n" << flush;
	}

    } // end nCyc loop over refinement cycles

  // Done with cycling, apply scale factors

  nNumReflns = m_poReflnlist->nGetNumReflns();
  for (j = 0; j < nNumReflns; j++)
    {
      poRefln = m_poReflnlist->poGetRefln(m_pnIndex[j]);

      // Calculate the scale factor for this refln
      // Apply factor even if reflection is not selected as long as it is
      // available.

      fTemp = fCalcGetScale(poRefln);

      // And apply it if the scale factor is legitimate
      if (0.0 < fTemp)
	{
	  poRefln->vSetIntensity(fTemp * poRefln->fGetIntensity());
	  poRefln->vSetSigmaI(fTemp * poRefln->fGetSigmaI());
	}
    }

  for (i = 1; i < m_nNumParams+1; i++)
    {
      delete [] ppfMat[i];
    }

  if (NULL != poReflnAve)
    {
      delete poReflnAve;
      poReflnAve = NULL;
    }
  if (bNewList)
    {
      delete poReflnlistAve;
      poReflnlistAve = NULL;
    }

  delete [] pfAabDiag;
  delete [] pfEigval;
  delete [] pfOffDiag;
  delete [] ppfMat;
  delete [] m_pnNumContrib;
  delete [] m_pnNumInBatch;
  delete [] m_pnNumAccepted;
  delete [] m_pnNumRejected;
  delete [] m_pnNumDeselected;
  delete [] m_pnNumOverlap;
  return (0);
}


float
Cscalemerge::fAverageSymRelated(const int nMethod, const int nNumRef,
                                const int *pnIndex, const float fRejCriteria,
				Crefln *poReflnAve,
				int *pnAccepted, int *pnNewRejects,
				float *pfSumResid)
{
  // Average a group of symmetry-related reflections to get intensity
  // and sigmaI for this reduced hkl
  // Reflections may or may not have already had scale factors applied to them:
  //  scale factors come from the method fCalcGetScale which may return 1.0.

  // Average Intensity is SUM(Wi*Ii)/ SUM(Wi), where Ii is observed intensity
  // for measurement i, and Wi is the weight of the measurement i, defined in
  // fCalcGetWeight() usually as
  // Wi =   1 / [ (SigmaI*m_fErrorMul)**2 + (Ii*m_fErrorAdd)**2 ]
  //
  // Ignore reflections that are deselected and which have SigmaI <= 0.
  //
  // SigmaI = sqrt(1 / SUM(Wi)), so returned sigma is adjusted by Error factors

  int    i, j, k, m;
  Crefln *poRefln;
  float   fSumIntWgt  = 0.0;
  float   fSumWgt     = 0.0;
  int     nAccepted   = 0;
  bool    bNewRejects = TRUE;
  float   fInt, fSigma, fWeight;
  float   fTestDev, fTestVar, fResidSq;
  float   fWorst;
  int     nWorst;
  float   fSumResidSq;
  float   fSumResid;
  int     nH, nK, nL;

  // Check if scratch arrays have enough space

  if (nNumRef > m_nMaxSymRelated)
    {
      // Work vectors are too small, so increase their size

      m_nMaxSymRelated = nNumRef;
      delete [] m_pfW1;
      delete [] m_pnW2;
      delete [] m_pfW3;
      delete [] m_pfW4;
      delete [] m_pnW5;
      m_pfW1 = new float [m_nMaxSymRelated];  // Holds intensity for ref i
      m_pnW2 = new int   [m_nMaxSymRelated];  // Not used here
      m_pfW3 = new float [m_nMaxSymRelated];  // Holds weight for ref i
      m_pfW4 = new float [m_nMaxSymRelated];  // Holds scale factor for ref i
      m_pnW5 = new int   [m_nMaxSymRelated];  // Holds anomalous designation
    }

  // Compute initial sums and number of accepted reflections

  float *pfInt;
  float *pfWeight;
  float *pfScale;
  int   *pnAnomFlag;
  int   *pnBatchIdx;
  int    nPack, nPackPrev;
  float fAverageSigma;

  pfInt         = m_pfW1;
  pfWeight      = m_pfW3;
  pnBatchIdx    = m_pnW2;
  pfScale       = m_pfW4;
  pnAnomFlag    = m_pnW5;
  poRefln       = m_poReflnlist->poGetRefln(pnIndex[0]);
  nPackPrev     = poRefln->nGetField(m_poReflnlist->m_nFI_nPackedHKL);

  j = 0;
  fAverageSigma = 0.0;
  for (i = 0; i < nNumRef;
       i++, pfInt++, pfWeight++, pnBatchIdx++, pfScale++, pnAnomFlag++)
    {
      *pfWeight   = 0.0;
      *pnBatchIdx = -1;
      *pfScale    = -1.0;

      poRefln     = m_poReflnlist->poGetRefln(pnIndex[i]);
      fSigma      = poRefln->fGetSigmaI();
      *pnAnomFlag = poRefln->nGetField(m_nFI_nAnomFlag);
//	  fAverageSigma += fSigma;

      if ( (poRefln->bIsSelected()) && (0.0 < fSigma) )
        {
	  if (0 == j)
	    {
	      // Get reduced hkl only first time

	      nH = poRefln->nGetField(m_poReflnlist->m_nFI_nReducedH);
	      nK = poRefln->nGetField(m_poReflnlist->m_nFI_nReducedK);
	      nL = poRefln->nGetField(m_poReflnlist->m_nFI_nReducedL);
	      poReflnAve->vSetH(nH);
	      poReflnAve->vSetK(nK);
	      poReflnAve->vSetL(nL);
	      j = 1;
	    }

	  nPack      = poRefln->nGetField(m_poReflnlist->m_nFI_nPackedHKL);
	  if (nPack != nPackPrev)
	    {
	      cerr << "ERROR in fAverageSymRelated: NOT sym-related!\n";
	      cerr << "nPackPrev: " << nPackPrev << "  nPack: " << nPack << '\n';
	    }
	  *pfScale    = fCalcGetScale(poRefln);
          *pfInt      = *pfScale * poRefln->fGetIntensity();
          *pfWeight   = fCalcGetWeight(*pfInt, fSigma * *pfScale, poRefln);
	  *pnBatchIdx = poRefln->nGetField(m_nFI_nBatchIdx);
          fSumIntWgt  = fSumIntWgt + *pfWeight * *pfInt;
          fSumWgt     = fSumWgt    + *pfWeight;
	  fAverageSigma += fSigma;
          nAccepted++;
        }
    } // end of i refln loop

//	fAverageSigma/=nNumRef;
	fAverageSigma/=nAccepted;

  // Now test each measurement against the average* and reject the measurement
  // with the largest deviation until there are no more rejects or until
  // there are just 2 measurements left.  The average* is calculated from the
  // OTHER reflns in the distribution (i.e. the test refln is excluded from
  // the calculation of the average to be tested against).

  fSumResidSq = 0.0;
  *pnNewRejects = 0;
  while ( (bNewRejects) && (2 <= nAccepted) )
    {
      // Keep doing this while there are new rejects
      // Find the measurement that deviates the most from the result

      bNewRejects = FALSE;                // Be optimistic!
      fWorst      = -1.0;
      nWorst      = -1;
      fSumResidSq = 0.0;
      fSumResid   = 0.0;
      pfInt       = m_pfW1;
      pfWeight    = m_pfW3;
      for (i = 0; i < nNumRef; i++, pfInt++, pfWeight++)
        {
          if (0.0 < *pfWeight)
            {
              // Subtract out contribution of this refln from sums
              //  and add back in variance of this refln

              fTestVar = 1.0f / (fSumWgt - *pfWeight);

              // Compute  <I> - Ii, where <I> is mean I without this refln

              fTestDev = (fSumIntWgt - (*pfInt * *pfWeight)) * fTestVar
		         - *pfInt;

              fTestVar = fTestVar  + 1.0f / *pfWeight;

              // For testing, compute  (delI/sigma)^2 and save largest one

	      fSumResid   = fSumResid + fabs(fTestDev);
              fTestDev    = fTestDev * fTestDev;
              fResidSq    = fTestDev / fTestVar;
              fSumResidSq = fSumResidSq + fResidSq;
              if (fResidSq >= fWorst)
                {
                  fWorst = fResidSq;
                  nWorst = i;
                }
            }
        }  // end of i refln loop

      // Rejection criteria of worst reflection.
      // Later on we'll add a real chisq test and Bayesian statistics.
      // Do not delete the last 2 though!

      if (fWorst > fRejCriteria)
        {
	  if (2 < nAccepted)
	    {
	      // Reject the worst refln by setting its SigmaI to be negative
	      // and its weight to be 0.0

	      fSumIntWgt     = fSumIntWgt - (m_pfW1[nWorst] * m_pfW3[nWorst]);
	      fSumWgt        = fSumWgt    - m_pfW3[nWorst];
	      m_pfW3[nWorst] = 0.0;
	      fAverageSigma = ((fAverageSigma)*nAccepted -
			       poRefln->fGetSigmaI())/(nAccepted - 1.0);

	      poRefln        = m_poReflnlist->poGetRefln(pnIndex[nWorst]);
	      poRefln->vSetSigmaI(-poRefln->fGetSigmaI());
	      nAccepted--;
	      (*pnNewRejects)++;
	      bNewRejects = TRUE;
	    }
	  else if (2 == nAccepted)
	    {
	      // Reject both remaining reflections

	      *pnNewRejects = *pnNewRejects + 2;
	      nAccepted     = 0;
	    }
	}
    } // end while

  // Special treatment of just 2 observations
  //  1. They match, average as above
  //  2. They do not match, reject lowest observation (if low enough) because
  //         it is behind beamstop or other shadow
  //  3. They do not match, reject highest observation (if high enough) because
  //          it is a zinger, or contaminated by ice or salt diffraction
  //  4. They do not match, reject both because you cannot decide which one
  //         is a better estimate of the true result
  //  5. They do not match, average both as above anyways.
  // The treatment to use depends on what you know about your detector and
  // the experiment performed.  For example, if you have taken care to not
  // predict reflections in shadowed regions of the detector, then you should
  // not do #2.  Or if you have no ice diffraction and an imaging plate, then
  // you should not do #3.

  if (2 == nAccepted)
    {
      // Do something based on the mode.  Be sure to adjust fSumResidSq
      // if necessary
    }

  if ( (0 != m_nAnomFlag) && (m_bLastCycle) )
    {
      // Calculate Int+ and Int-
      // Even thought centric reflections should not have anomalous scattering
      // it is useful to still designate + and - for them.

      float   fIntPlus, fIntMinus, fSigmaIPlus, fSigmaIMinus;
      float   *pfSumInt, *pfSumWeight;

      fIntPlus = fIntMinus = fSigmaIPlus = fSigmaIMinus = 0.0;

      pfInt       = m_pfW1;
      pfWeight    = m_pfW3;
      pnAnomFlag  = m_pnW5;
      for (i = 0; i < nNumRef; i++, pfInt++, pfWeight++, pnAnomFlag++)
        {
          if (0.0 < *pfWeight)
	    {
	      if (0 < *pnAnomFlag)
		{
		  pfSumInt    = &fIntPlus;
		  pfSumWeight = &fSigmaIPlus;
		}
	      else
		{
		  pfSumInt    = &fIntMinus;
		  pfSumWeight = &fSigmaIMinus;
		}
	      *pfSumInt    += *pfInt * *pfWeight;
	      *pfSumWeight += *pfWeight;
	    }
	}
      if (0.0 < fSigmaIPlus)
	{
	  fIntPlus     = fIntPlus / fSigmaIPlus;
	  fSigmaIPlus  = 1.0f / sqrtf(fSigmaIPlus);
	}
      else
	{
	  fIntPlus    = -1.0;
	  fSigmaIPlus = -1.0;
	}
      poReflnAve->vSetField(m_nFI_fIntensityPlus, fIntPlus);
      poReflnAve->vSetField(m_nFI_fSigmaIPlus, fSigmaIPlus);
      if (0.0 < fSigmaIMinus)
	{
	  fIntMinus    = fIntMinus / fSigmaIMinus;
	  fSigmaIMinus = 1.0f / sqrtf(fSigmaIMinus);
	}
      else
	{
	  fIntMinus    = -1.0;
	  fSigmaIMinus = -1.0;
	}
      poReflnAve->vSetField(m_nFI_fIntensityMinus, fIntMinus);
      poReflnAve->vSetField(m_nFI_fSigmaIMinus, fSigmaIMinus);
    }

  if (0.0 < fSumWgt)
    {
      fInt   = fSumIntWgt / fSumWgt;
      fSigma = 1.0f / sqrtf(fSumWgt);
    }
  else
    {
      fInt   = -1.0;
      fSigma = -1.0;
    }

/*
//+jwp  31-Oct-1998
// Put average Int and Sigma in the scaled, but unaveraged reflnlist
  for (i = 0; i < nNumRef; i++)
    {
      poRefln     = m_poReflnlist->poGetRefln(pnIndex[i]);
      poRefln->vSetField(m_nFI_fOtherInt, fInt);
      poRefln->vSetField(m_nFI_fOtherSig, fSigma);
    }
//-jwp
*/
  // Accumulate statistics

  tagStats *ptBatch, *ptInt, *ptReso;
  float fReso;

  // Get resolution bin

  int nBinReso;
  poRefln     = m_poReflnlist->poGetRefln(pnIndex[0]);
  fReso       = poRefln->fGetField(m_nFI_f2STLsq);
  fReso       = fReso * sqrtf(fReso);
  nBinReso    = (int) ((fReso - m_a2fRangeReso[0]) / m_fSlopeReso);
  if (0 > nBinReso)
    {
      nBinReso = 0;
    }
  else if (m_nNumBinsReso <= nBinReso)
    {
      nBinReso = m_nNumBinsReso-1;
    }

  // Get intensity / sigma bin

  int nBinInt;
  fReso   = fInt / fSigma;
  nBinInt = (int) (fReso / m_fSlopeInt);
  if (0 > nBinInt)
    {
      nBinInt = 0;
    }
  else if (m_nNumBinsInt <= nBinInt)
    {
      nBinInt = m_nNumBinsInt-1;
    }

  ptReso  = &m_ptStatsReso[nBinReso];
  ptInt   = &m_ptStatsInt[nBinInt];

  // Get multiplicity bin.  0, 1, 2, 3, 4, 5-6, 7-8, 9-12, 13-19, >19
  // Get multiplicity bin.  0, 1, 2, 3, 4, 5-8, 9-12, >12

  int nBinMult;
  if (5 > nAccepted)
    {
      nBinMult = nAccepted;
    }
  else if (12 <  nAccepted)
    {
      nBinMult = 7;
    }
  else if (8 < nAccepted)
    {
      nBinMult = 6;
    }
  else if (4 < nAccepted)
    {
      nBinMult = 5;
    }

  if (1 > nAccepted)
    {
      // Do nothing but return bogus values for the averaged intensity and its
      // standard deviation (we need to count rejected reflections in the
      // correct statistical slots!)

      poReflnAve->vSetIntensity((float)-1.0);
      poReflnAve->vSetSigmaI((float)-1.0);

      pnBatchIdx  = m_pnW2;
      for (i = 0; i < nNumRef; i++, pnBatchIdx++)
	{
	  if (0 <= *pnBatchIdx)
	    {
	      // Only those reflections with non-negative batch indices are
	      // part of this group

	      ptBatch = &m_ptStatsBatch[*pnBatchIdx];
	      m_pnNumRejected[*pnBatchIdx]++;   // Count rejected refs

	      ptReso->nNumRejects++;
	      ptInt->nNumRejects++;
	      ptBatch->nNumRejects++;

	      ptBatch->nNumRefs++;
	      ptReso->nNumRefs++;
	      ptInt->nNumRefs++;
	    }
	}
    }
  else if (1 == nAccepted)
    {
      // Special case of just 1 refln.  Do not add in contributions to
      // statistics that require more than 1 refln for legitimacy.

      poReflnAve->vSetIntensity(fInt);
      if (m_bSigmaOverEstimate)
	poReflnAve->vSetSigmaI(max(fAverageSigma,fSigma));
      else
	poReflnAve->vSetSigmaI(fSigma);

      // Even though nAccepted==1, we could have nNumRef > 1, so we
      // must find which batch we are dealing with by scanning m_pnW2.

      i = 0;
      while ( (i < nNumRef) && (-1 == m_pnW2[i]) ) i++;

      if (i == nNumRef)
	{
	  // Did not find a batch number
	  // (This is defensive, we never expect to get here)

	  cerr << "fAvgSymRelated:: PROGRAMMER ERROR!\n";
	}
      else
	{                                           //  refln is the refln
	  ptBatch = &m_ptStatsBatch[m_pnW2[i]];
	  ptBatch->nNumRefs++;                      // Count all accepted refs
	  ptBatch->nNumSingles++;                   // Count singles
	  ptBatch->anNumInBinMults[nBinMult]++;     // Count singles again
	  ptBatch->anNumInBin[nBinInt]++;           // Count in I/sig ranges

	  ptReso->nNumRefs++;                       // Count all refs
	  ptReso->nNumSingles++;                    // Count singles
	  ptReso->anNumInBinMults[nBinMult]++;      // Count multiplicity
	  ptReso->anNumInBin[nBinInt]++;            // Count I/sig in reso

	  ptInt->nNumRefs++;                        // Count all refs
	  ptInt->nNumSingles++;                     // Count singles
	  ptInt->anNumInBinMults[nBinMult]++;       // Count multiplicity
	  ptInt->anNumInBin[nBinInt]++;


	  ptBatch->dSumIoverSig = ptBatch->dSumIoverSig + fReso;
	  ptReso->dSumIoverSig  = ptReso->dSumIoverSig  + fReso;
	  ptInt->dSumIoverSig   = ptInt->dSumIoverSig   + fReso;

	  // Singles DO NOT contribute to Rmerge and ChiSq

	  //      ptBatch->nNumUsed++;
	  //      ptBatch->dNumer      = ptBatch->dNumer       + fTestDev;
	  //      ptBatch->dDenom      = ptBatch->dDenom       + fInt;
	  //      ptBatch->dSumChiSq   = ptBatch->dSumChiSq    + fSumResidSq;
	}
    }
  else if (1 < nAccepted)
    {

      // Accumulate some statistics

      poReflnAve->vSetIntensity(fInt);
      if (m_bSigmaOverEstimate)
	poReflnAve->vSetSigmaI(max(fAverageSigma,fSigma));
      else
	poReflnAve->vSetSigmaI(fSigma);

      fSumResidSq = fSumResidSq / sqrtf((float)nAccepted);

      ptReso->nNumMults++;
      ptInt->nNumMults++;
      ptInt->dSumIoverSig  = ptInt->dSumIoverSig   + fReso;
      ptReso->dSumIoverSig = ptReso->dSumIoverSig  + fReso;

//      ptBatch->anNumInBinMults[nBinMult]++;     // Count singles again
//      ptBatch->anNumInBin[nBinInt]++;           // Count in I/sig ranges
      ptReso->anNumInBinMults[nBinMult]++;      // Count multiplicity
      ptReso->anNumInBin[nBinInt]++;            // Count I/sig in reso
      ptInt->anNumInBinMults[nBinMult]++;       // Count multiplicity
      ptInt->anNumInBin[nBinInt]++;

      // Find out how many different batches this contributes to

      ptBatch = &m_ptStatsBatch[0];
      for (i = 0; i < m_nNumBatches; i++, ptBatch++)
	{
	  j = 0;
	  m = 0;
	  pnBatchIdx = m_pnW2;
	  pfWeight    = m_pfW3;
	  for (k = 0; (k < nNumRef) && (0 == j); k++, pnBatchIdx++)
	    {
	      if ( (0 < *pfWeight) && (*pnBatchIdx == i) )
		{
		  // Found once that belongs to batch i

		  ptBatch->nNumMults++; // Count batch mults just once
		  j = 1;
		}
	    }
	}  // end i loop


      *pfSumResid = 0.0;
      pfInt       = m_pfW1;
      pnBatchIdx  = m_pnW2;
      pfWeight    = m_pfW3;
      pfScale     = m_pfW4;
      fWeight     = (float)(nAccepted-1) / (float) nAccepted;
      for (i = 0; i < nNumRef; i++, pfInt++, pfWeight++, pnBatchIdx++,
	                       pfScale++)
	{
	  ptBatch = &m_ptStatsBatch[*pnBatchIdx];
	  if (0.0 < *pfWeight)
	    {
	      fTestDev = fabs(*pfInt - fInt);             // |Icorr - <Ih>|
	      fResidSq = fTestDev * fTestDev;            // (Icorr - <Ih>)^2
	      fTestVar = *pfInt * *pfInt;

//+1-May	      fWorst   = *pfWeight * fResidSq * fWeight; // ChiSq
//	      fWorst   = *pfWeight * fResidSq / (float) nAccepted; // ChiSq
//	      fWorst   = *pfWeight * fResidSq * fWeight; // ChiSq
	      fWorst   = *pfWeight * fResidSq;
//-1-May
	      fReso    = *pfInt * sqrtf(*pfWeight);      // I/SigmaI,

	      *pfSumResid = *pfSumResid + fTestDev;
	      m_pnNumContrib[*pnBatchIdx]++;
	      m_pnNumAccepted[*pnBatchIdx]++;            // Count accepted refs

	      ptReso->nNumRefs++;
	      ptReso->nNumUsed++;
	      ptReso->dNumer       = ptReso->dNumer        + fTestDev;
	      ptReso->dDenom       = ptReso->dDenom        + fabs(fInt);
	      ptReso->dWtNumer     = ptReso->dWtNumer      + fResidSq;
	      ptReso->dWtDenom     = ptReso->dWtDenom      + fTestVar;
	      ptReso->dSumChiSq    = ptReso->dSumChiSq     + fWorst;

	      ptInt->nNumRefs++;
	      ptInt->nNumUsed++;
	      ptInt->dNumer        = ptInt->dNumer         + fTestDev;
	      ptInt->dDenom        = ptInt->dDenom         + fabs(fInt);
	      ptInt->dWtNumer      = ptInt->dWtNumer       + fResidSq;
	      ptInt->dWtDenom      = ptInt->dWtDenom       + fTestVar;
	      ptInt->dSumChiSq     = ptInt->dSumChiSq      + fWorst;

	      ptBatch->nNumRefs++;
	      ptBatch->nNumUsed++;
	      ptBatch->dNumer      = ptBatch->dNumer       + fTestDev;
	      ptBatch->dDenom      = ptBatch->dDenom       + fabs(fInt);
	      ptBatch->dWtNumer    = ptBatch->dWtNumer     + fResidSq;
	      ptBatch->dWtDenom    = ptBatch->dWtDenom     + fTestVar;
	      ptBatch->dSumIoverSig= ptBatch->dSumIoverSig + fReso;
	      ptBatch->dSumChiSq   = ptBatch->dSumChiSq    + fWorst;

	    }
	  else if (0 <= *pnBatchIdx)
	    {
	      m_pnNumRejected[*pnBatchIdx]++;   // Count rejected refs

	      ptReso->nNumRejects++;
	      ptInt->nNumRejects++;
	      ptBatch->nNumRejects++;

	      ptBatch->nNumRefs++;
	      ptReso->nNumRefs++;
	      ptInt->nNumRefs++;
	    }
	}

      if (0 != m_nCountOverlap)
	{
	  // Count overlaps between batches, note that the diagonal of
	  // m_pnNumOverlap counts symmetry-related reflns in the same batch,
	  // but does not count reflns that overlap just with themselves!
	  // Note that the UPPER RIGHT and LOWER LEFT of the matrix are not the
	  // same!  See documentation for more info.

	  for (i = 0; i < nNumRef; i++)
	    {
	      if (0 <= m_pnW2[i])
		{
		  for (j = 0; j < m_nNumBatches; j++)
		    {
		      m = 0;
		      for (k = 0; (k < nNumRef) && (0 == m); k++)
			{
			  if ( (k != i) && (j == m_pnW2[k]) )
			    {
			      m = 1;
			    }
			}
		      if (1 == m)
			{
			  // An overlap between batch m_pnW2[i} and batch j
			  // was found
			
			  k = m_pnW2[i] * m_nNumBatches + j;
			  m_pnNumOverlap[k]++;
			}
		    }
		}
	    } // end i loop
	}  // endif m_nCountOverap
    }

  *pnAccepted = nAccepted;
  return (fSumResidSq);
}

float
Cscalemerge::fCalcGetWeight(const float fIntensity, const float fSigma,
			    Crefln *poRefln)
{
  // Calculate a weight for a reflection to be used in various other
  // places
  // The weight W of the measurement is defined as
  // W =   1 / [ (SigmaI*m_fErrorMul)**2 + (I*m_fErrorAdd)**2 ]
  // Note:  This returns the weight even if fSigma is negative!
  // Also W >= 0.

  float fWeight;
  float fTemp;
  fWeight = fSigma * m_fErrorMul;
  fTemp   = fIntensity * m_fErrorAdd;
  fWeight = fWeight * fWeight  +  fTemp * fTemp;
  if (0.0 < fWeight)
    {
      fWeight = 1.0f / fWeight;
    }
  else
    {
      fWeight = 0.0;
    }
  return (fWeight);
//  return (1.0);         // Unit weights!


}

void
Cscalemerge::vCalcGetWeights(Crefln *poRefln)
{
  // Calculate the weights of a refln to use to apply the refln int and sigmao
  // to a scaling batch or subbatch for the Kabsch method!
  // The weight will be based strictly on position, either on the detector,
  // or in reciprocal space or ...
  // Normalize weights so Sum(Wgt) = 1.

  // Subbatches are divided into resolution and rotation angle range
  // weight = exp(-0.5 * delta**2), if weight < 0.018, weight = 0.0
  // if the refln does not belong to the batch in question, then the
  // weight is 0.0

  int   i;
  float fTest;
  fTest = 1.0f / (float) m_nNumParams;
  for (i = 0; i < m_nNumParams; i++)
    {
      m_pfWhla[i] = (float) 1.0 / fTest;   // For now equal weights
    }
}

float
Cscalemerge::fCalcGetScale(Crefln *poRefln)
{
  // Calculate the scale factor to apply to a given observed reflection.
  // The scale factor is a function of its "position" and the
  // least-squares determined scale factors at various "positions"
  // Maybe in the future also return d(scale)/d(Pi) (see routine vKandB()).
  //

  float fScale;
  if ( (1 == m_nMethod) || (2 == m_nMethod) )
    {
      // Fox and Holmes method
      //   Ki    =     ki    * exp(-2.0 *       Bi             * s^2)
      //  <Ih>   ~=    Ki * Ihi  (i.e. scale factor is applied to observed int)

      int nIdx    = poRefln->nGetField(m_nFI_nBatchIdx);
      if ( (nIdx < 0) || (nIdx >= m_nNumBatches) )
	return (-1.0);

      fScale = 0.25f * poRefln->fGetField(m_nFI_f2STLsq);
      if (0.0 > fScale)
	return (-1.0);

      fScale = m_pfGa[nIdx] * (float)exp(-2.0 * m_pfGa[m_nNumBatches + nIdx] * fScale);
      return (fScale);
    }

  // For other methods return 1.0 for now
  return (1.0);
}

void
Cscalemerge::vKabsch(const int nNumRefs, const int *pnIndex)
{
  // Place sums of into the normal matrix using the Kabsch method
  // for a group of symmetry-related reflns

  int i, j;
  Crefln *poRefln;
  float   fWeight;
  float   fInt, fSigmaI;

  float   fIh;             // Vh / Uh                   Kabsch Eqn  5
  float   fUh;             // Sum( Ga * Ga * Uha)              Eqn  9
  float   fVh;             // Sum( Ga * Vha)                   Eqn 10
  float   fVhaOvUh;        // Vha / Uh                     for Eqn  3
  float   fRha_VhaOvUh;    // (Rha - Vha) / Uh             for Eqn  3
  float  *pfAab;           // Convenience pointer in m_pfAab

  if (2 > nNumRefs)
    return;             // Nothing to do!

  // Zero some sums

  for (i = 0; i < m_nNumParams; i++)
    {
      m_pfUha[i] = 0.0;
      m_pfVha[i] = 0.0;
    }

  // 1. Determine which scalebatches these reflections belong to
  //             (maybe do this outside?)
  // 2. Determine which reflns should be rejected since they are outliers
  // 3. Determine weights of these reflns vs each scalebatch
  // 4. Normalize weights so sum(weights) = 1 (or nNumRef?)
  // 5. Compute sums

  for (i = 0; i < nNumRefs; i++)
    {
      poRefln = m_poReflnlist->poGetRefln(pnIndex[i]);
      fSigmaI = poRefln->fGetSigmaI();
      if ( (poRefln->bIsSelected()) && (0.0 < fSigmaI) )
	{
	  // Compute weights of this reflection into each scale batch group
          // subject to condition that Sum(Wi) = 1.0.
	  // Weights are returned m_pfWhla

	  vCalcGetWeights(poRefln);

          // Add weighted contribution of this refln to the sums

          for (j = 0; j < m_nNumParams; j++)
            {
              if (0.0 < m_pfWhla[j])
                {
		  fInt       = poRefln->fGetIntensity();
                  fWeight    = m_pfWhla[j] * fCalcGetWeight(fInt, fSigmaI, poRefln);
                  m_pfUha[j] = m_pfUha[j] + fWeight;                 // Eqn 7
                  m_pfVha[j] = m_pfVha[j] + fWeight * fInt;          // Eqn 8
                }
            } // end j loop over batches and subbatches
	}  // endif selection
    } // end of i refln loop

  fUh = 0.0;
  fVh = 0.0;

  for (i = 0; i < m_nNumParams; i++)
    {
      if (0.0 < m_pfUha[i])
	{
	  fUh = fUh + m_pfUha[i] * m_pfGa[i] * m_pfGa[i];            // Eqn 9
	  fVh = fVh + m_pfVha[i] * m_pfGa[i];                        // Eqn 10
	}
    } // end i loop

  // Compute normalizer for this group of sym-related reflections contribution
  // to the sums ...

  if (0.0 < fUh)
    {
      fIh = fVh / fUh;                                                // Eqn 5
      for (i = 0; i < m_nNumParams; i++)
	{
	  m_pfRha[i] = m_pfVha[i] -  m_pfGa[i] * m_pfUha[i] * fIh;    // Eqn 6
	}

      pfAab = m_pfAab;

      for (i = 0; i < m_nNumParams; i++)
	{
	  m_pfBa[i]  = m_pfBa[i] + fIh * m_pfRha[i];                  // Eqn 4
	
	  // Eqn 3 is complicated, so factorize a bit
	
	  fRha_VhaOvUh = (m_pfRha[i] - m_pfVha[i])                    // Eqn 3a
                         / fUh;
	  fVhaOvUh     = m_pfVha[i] / fUh;                            // Eqn 3b
	
	  // Even though matrix is symmetric, sum
	  // both upper and lower triangular elements
	  // for clarity, maybe change this later
	
	  for (j = 0; j < m_nNumParams; j++)
	    {
	      *pfAab = *pfAab + fRha_VhaOvUh * m_pfVha[j]             // Eqn 3c
                              + fVhaOvUh * m_pfRha[j];
	      pfAab++;
	    }
	
	  // Now add fIh * fIh * m_pfUha to diagnl elem.
	
	  j = (m_nNumParams + 1) * i;
	  m_pfAab[j] = m_pfAab[j] + fIh * fIh * m_pfUha[i];           // Eqn 3d
	}
    }

  // All done, so return to our regularly scheduled normal equations matrix

  return;
}

void
Cscalemerge::vKandB(const int nNumRefs, const int *pnIndex)
{
  // Place sums into the normal matrix using the Fox and Holmes method
  // for a group of symmetry-related reflns
  //
  //  Let Ihi be the ith observation for a given unique hkl h with
  //  Let Kj  be the scale factor such that Kj*Ihi ~= <Ih> the weighted average
  //  for all reflns with hkl symmetry equivalent to h.
  //  Whi = 1 / [(sighi * E1)^2 + (<Ih> * E2)^2], a weight for observation Ihi
  //  Let Kj = kj * exp(-2Bjs^2) where s = sinTheta/lambda
  //  Let Gj = 1 / Kj.
  //  <Ih> = (Sum Whi*Gj*Ihi) / (Sum Whi*Gj*Gj)   (ith refln uses scale fact Gj)
  //
  //  Minimize    ChiSq = Sum Sum [Whi(Ihi - Gj<Ih>)^2]
  //                       h   i
  //  in order to find kj, Bj for all batches of reflns.


  int i, j, k, mi, mj;
  int     nIdx;
  Crefln *poRefln;
  float   fSigmaI;
  float   fWeight;
  float   fInt;
  float   fSinTLSq;

  float   fIh;             // Vh / Uh = <Ih>
  float   fUh;             // Sum( Gj * Gj * whi)
  float   fVh;             // Sum( Gj * whi * Ihi)
  float   fTemp;           // Temp var
  float  *pfAab;           // Convenience pointer in m_pfAab

  if (2 > nNumRefs)
    return;             // Nothing to do!

  int nBoffset = nNumRefs * m_nNumBatches;

  if ((nNumRefs * m_nNumParams) > m_nMaxSymRelated)
    {
      // Work vectors are too small, so increase their size

      m_nMaxSymRelated = nNumRefs * m_nNumParams;
      // KandB
      delete [] m_pfW1;
      delete [] m_pnW2;
      delete [] m_pfW3;
      delete [] m_pfW4;
      delete [] m_pnW5;
      m_pfW1 = new float [m_nMaxSymRelated];
      m_pnW2 = new int   [m_nMaxSymRelated];
      m_pfW3 = new float [m_nMaxSymRelated];
      m_pfW4 = new float [m_nMaxSymRelated];
      m_pnW5 = new int   [m_nMaxSymRelated];
    }

  // Zero some sums and arrays

  for (i = 0; i < m_nNumParams; i++)
    {
      m_pfRha[i] = 0.0;
    }
  for (i = 0; i < nNumRefs * m_nNumParams; i++)
    {
      m_pfW1[i] = 0.0;
    }

  // 1. Compute scale factors Gi for this group of reflns
  //           Gi = 1 / Ki;  Ki = ki * exp(-2Bis^2); s = sin^2theta/lambda^2
  //    Place Gi in m_pfUha[i]

  //    ki is found in m_pfGa[i]            (i.e. 0, 1, 2, ...m_nNumBatches-1)
  //    Bi is found in m_pfGa[i + nBatches] (i.e. m_nNumBatches, nB+1, nB+2,...)

  //	Also compute -Gj/kj and Gj*2s^2,  [  del(Gj)/del(kj), del(Gj)/del(Bj) ]
  //     place these in m_pfVha[2*i], m_pfVha[2*i+1] respectively
	
  poRefln  = m_poReflnlist->poGetRefln(pnIndex[0]);
  fSinTLSq = poRefln->fGetField(m_nFI_f2STLsq);  // This same for all refs here
                                                 //  in sym-related group
  for (i = 0; i < m_nNumBatches; i++)
    {
      //   Ki    =     ki    *        exp(-2.0 *       Bi                  * s^2)
      m_pfUha[i] = m_pfGa[i] * (float)exp(-2.0 * m_pfGa[m_nNumBatches + i] * fSinTLSq);
      if (0.0 != m_pfUha[i])
	{
	  m_pfUha[i]                = 1.0f / m_pfUha[i];        // Gi = 1 / Ki
	  m_pfVha[i]                = -m_pfUha[i] / m_pfGa[i]; // Vha[0,1,2,...]
	  m_pfVha[i+ m_nNumBatches] =  m_pfUha[i] * 2.0f * fSinTLSq;
	}                                                      // Vha[nB,nB+1,..
    } // end i loop

  // 2. Compute <Ih> = Sum(Whi * Ihi * Gj) / Sum(Whi * Gj * Gj)

  fUh = 0.0;
  fVh = 0.0;

  for (i = 0; i < nNumRefs; i++)
    {
      poRefln   = m_poReflnlist->poGetRefln(pnIndex[i]);
      fSigmaI   = poRefln->fGetSigmaI();
      nIdx      = poRefln->nGetField(m_nFI_nBatchIdx);
      m_pnW2[i] = -1;
      if (nIdx >= m_nNumBatches)
	{
	  cerr << "BATCH NUMBER NOT VALID!\n";
	}
      else if (poRefln->bIsSelected())
	{
	  if (0.0 < fSigmaI)
	    {
//	      m_pnW2[i] = nIdx;
//	      m_pnNumContrib[nIdx]++;
//	      m_pnNumAccepted[nIdx]++;   // Count accepted refs
	      fInt    = poRefln->fGetIntensity();
	      fWeight = fCalcGetWeight(fInt, fSigmaI, poRefln);
	      fVh     = fVh + fWeight * m_pfUha[nIdx] * fInt;
	      fUh     = fUh + fWeight * m_pfUha[nIdx] * m_pfUha[nIdx];
	    }  // endif selection
	  else
	    {
//	      m_pnNumRejected[nIdx]++;   // Count rejected refs
	    }
	}
    } // end of i refln loop

  if (0.0 < fUh)
    {
      fIh = fVh / fUh;                            // <Ih>

      // Compute derivatives pairwise, for each kj and Bj.

      for (i = 0; i < m_nNumBatches; i++)
	{
	  // Compute del(Gi<Ih>)/del(Pk) for Pk = (kk,Bk),
	  //                                     k = 1, m_nNumBatches
	  // Place in m_pfW1(k) and m_pfW1(k+m_nNumBatches)

          // dGj/dkj in [i], dGj/dBj in [i+m_nNumBatches]

	  for (j = 0; j < nNumRefs; j++)
	    {
	      poRefln = m_poReflnlist->poGetRefln(pnIndex[j]);
	      fSigmaI = poRefln->fGetSigmaI();
	      if ( (poRefln->bIsSelected()) && (0.0 < fSigmaI) )
		{
		  nIdx    = poRefln->nGetField(m_nFI_nBatchIdx);
		  fInt    = poRefln->fGetIntensity();
		  fWeight = fCalcGetWeight(fInt, fSigmaI, poRefln);
                                                            // Could use fIh
                                                            // instead of fInt
		  fTemp   = fWeight * m_pfUha[i]            // Whj * Gi *
                            * (fInt - 2.0f * m_pfUha[nIdx] * fIh);//(Ihi-2Gj<Ih>)
		  fTemp   = fTemp / fUh;                    //   / Sum(Whi*Gi^2)
		
		  // Remember that m_pfVha[j]               =  -Gj/kj
		  //  and          m_pfVha[j+m_nNumBatches] = 2*Gj*s^2
		
		  k = i * m_nNumBatches + j;
		  m_pfW1[k]           = fTemp * m_pfVha[nIdx];
 		  m_pfW1[k+nBoffset]  = fTemp * m_pfVha[nIdx + m_nNumBatches];

		  // Add in extra terms if this refln belongs to batch i

		  if (i == nIdx)
		    {
		      // We are dealing with something on the diagonal

		      m_pfW1[k]          = m_pfW1[k] + fIh * m_pfVha[nIdx];
		      m_pfW1[k+nBoffset] = m_pfW1[k+nBoffset]
			                  + fIh * m_pfVha[nIdx + m_nNumBatches];
		    }
		} // end refln loop
	    }
	} // Ok, all derivatives are calculated

      if (5 < m_nVerbose)
	{
	  cout << "Deriv pfW1: \n" << flush;
	  vListMatMN(nNumRefs*m_nNumParams, 1, m_pfW1);
	}

      // Now add contributions to normal matrix
      // Even though matrix is symmetric, sum
      // both upper and lower triangular elements
      // for clarity, maybe change this later

      int ij, ji;
      for (k = 0; k < nNumRefs; k++)
	{
	  poRefln = m_poReflnlist->poGetRefln(pnIndex[k]);
	  fSigmaI = poRefln->fGetSigmaI();
	  if ( (poRefln->bIsSelected()) && (0.0 < fSigmaI) )
	    {
	      nIdx    = poRefln->nGetField(m_nFI_nBatchIdx);
	      fInt    = poRefln->fGetIntensity();
	      fWeight = fCalcGetWeight(fInt, fSigmaI, poRefln);

	      pfAab   = m_pfAab;
	      for (i = 0; i < m_nNumParams; i++)
		{
		  mi        = i * m_nNumBatches + k;
		  m_pfBa[i] = m_pfBa[i]  + (fInt - m_pfUha[nIdx] * fIh)
//		                            * m_pfW1[mi] * fWeight;
		                            * m_pfW1[mi];

		  ij = i * m_nNumParams - 1;
		  ji = i - m_nNumParams;
		  for (j = 0; j <= i; j++)
//		  for (j = 0; j < m_nNumParams; j++)
		    {
		      // A lot of these derivatives might be 0s, but go ahead
	
		      ij++;
		      ji = ji + m_nNumParams;

		      mj =  j * m_nNumBatches + k;
		      m_pfAab[ij] = m_pfAab[ij] + (m_pfW1[mi] * m_pfW1[mj]);
//		      m_pfAab[ij] = m_pfAab[ij] + (m_pfW1[mi] * m_pfW1[mj]) * fWeight;
		      m_pfAab[ji] = m_pfAab[ij];
//		      *pfAab = *pfAab + (m_pfW1[mi] * m_pfW1[mj]) * fWeight;
		      pfAab++;  // Next element
		    } // end j
		} // end i
	    } // endif
	} // end k refln loop
    } // end fUh > 0

  if (6 < m_nVerbose)
    {
      // Print out the matrix

      cout << "m_pfAab: \n";
      vListMatMN(m_nNumParams, m_nNumParams, m_pfAab);
      cout << "m_pfBa: \n";
      vListMatMN(m_nNumParams, 1, m_pfBa);
#ifdef SSI_PC
      cout << flush;
#endif
    }

  // All done, so return to our regularly scheduled normal equations matrix

  return;
}

void
Cscalemerge::vKandB2(const int nNumRefs, const int *pnIndex)
{
  // Place sums into the normal matrix using the Fox and Holmes method
  // for a group of symmetry-related reflns
  //
  //  Let Ihi be the ith observation for a given unique hkl h with
  //  Let Kj  be the scale factor such that Kj*Ihi ~= <Ih> the weighted average
  //  for all reflns with hkl symmetry equivalent to h.
  //  Whi = 1 / [(sighi * E1)^2 + (<Ih> * E2)^2], a weight for observation Ihi
  //  Let Kj = kj * exp(-2Bjs^2) where s = sinTheta/lambda
  //  Let Gj = 1 / Kj.
  //  <Ih> = (Sum Whi*Gj*Ihi) / (Sum Whi*Gj*Gj)   (ith refln uses scale fact Gj)
  //
  //  Minimize    ChiSq = Sum Sum [Whi(Ihi - Gj<Ih>)^2]
  //                       h   i
  //  in order to find kj, Bj for all batches of reflns.


  int i, j, k, m, ij, ji;
  int     nIdx;
  Crefln *poRefln;
  float   fSigmaI;
  float   fWeight;
  float   fInt;
  float   fSinTLSq;

  float   fIh;             // Vh / Uh = <Ih>
  float   fUh;             // Sum( Gj * Gj * whi)
  float   fVh;             // Sum( Gj * whi * Ihi)
  float   fTemp;           // Temp var

  if (2 > nNumRefs)
    return;             // Nothing to do!

//  int nBoffset = nNumRefs * m_nNumBatches;

  if ((nNumRefs * m_nNumParams) > m_nMaxSymRelated)
    {
      // Work vectors are too small, so increase their size

      m_nMaxSymRelated = nNumRefs * m_nNumParams;
      // KandB2
      delete [] m_pfW1;
      delete [] m_pnW2;
      delete [] m_pfW3;
      delete [] m_pfW4;
      delete [] m_pnW5;
      m_pfW1 = new float [m_nMaxSymRelated];
      m_pnW2 = new int   [m_nMaxSymRelated];
      m_pfW3 = new float [m_nMaxSymRelated];
      m_pfW4 = new float [m_nMaxSymRelated];
      m_pnW5 = new int   [m_nMaxSymRelated];
    }

  // Zero some sums and arrays

  for (i = 0; i < m_nNumParams; i++)
    {
      m_pfRha[i] = 0.0;
    }
  for (i = 0; i < nNumRefs * m_nNumParams; i++)
    {
      m_pfW1[i] = 0.0;
    }

  // 1. Compute scale factors Gi for this group of reflns
  //           Gi = 1 / Ki;  Ki = ki * exp(-2Bis^2); s = sin^2theta/lambda^2
  //    Place Gi in m_pfUha[i]

  //    ki is found in m_pfGa[i]            (i.e. 0, 1, 2, ...m_nNumBatches-1)
  //    Bi is found in m_pfGa[i + nBatches] (i.e. m_nNumBatches, nB+1, nB+2,...)

  //	Also compute -Gj/kj and Gj*2s^2,  [  del(Gj)/del(kj), del(Gj)/del(Bj) ]
  //     place these in m_pfVha[2*i], m_pfVha[2*i+1] respectively
	
  poRefln  = m_poReflnlist->poGetRefln(pnIndex[0]);
  fSinTLSq = poRefln->fGetField(m_nFI_f2STLsq);  // This same for all refs here
                                                 //  in sym-related group
  for (i = 0; i < m_nNumBatches; i++)
    {
      //   Ki    =     ki    *        exp(-2.0 *       Bi                  * s^2)
      m_pfUha[i] = m_pfGa[i] * (float)exp(-2.0 * m_pfGa[m_nNumBatches + i] * fSinTLSq);
      if (0.0 != m_pfUha[i])
	{
	  m_pfUha[i]                = 1.0f / m_pfUha[i];        // Gi = 1 / Ki
	  m_pfVha[i]                = -m_pfUha[i] / m_pfGa[i]; // Vha[0,1,2,...]
	  m_pfVha[i+ m_nNumBatches] =  m_pfUha[i] * 2.0f * fSinTLSq;
	}                                                      // Vha[nB,nB+1,..
    } // end i loop

  // 2. Compute <Ih> = Sum(Whi * Ihi * Gj) / Sum(Whi * Gj * Gj)

  fUh = 0.0;
  fVh = 0.0;

  for (i = 0; i < nNumRefs; i++)
    {
      poRefln   = m_poReflnlist->poGetRefln(pnIndex[i]);
      fSigmaI   = poRefln->fGetSigmaI();
      nIdx      = poRefln->nGetField(m_nFI_nBatchIdx);
      if (nIdx >= m_nNumBatches)
	{
	  cerr << "BATCH NUMBER NOT VALID!\n";
	}
      else if (poRefln->bIsSelected())
	{
	  if (0.0 < fSigmaI)
	    {
	      fInt    = poRefln->fGetIntensity();
	      fWeight = fCalcGetWeight(fInt, fSigmaI, poRefln);
	      fVh     = fVh + fWeight * m_pfUha[nIdx] * fInt;
	      fUh     = fUh + fWeight * m_pfUha[nIdx] * m_pfUha[nIdx];
	    }  // endif selection
	}
    } // end of i refln loop

  if (0.0 < fUh)
    {
      fIh = fVh / fUh;                            // <Ih>

      // Now add contributions to normal matrix
      // Matrix is symmetric, so compute a triangle and set other half to it

      int    nBi, nBj;
      float  fIntI, fSigmaII, fWeightI, fIntJ, fSigmaIJ, fWeightJ;
      float  fDelGi, fDelGj;
      float  fGi, fGj;
      Crefln *poReflnI;
      Crefln *poReflnJ;

      // Fill in only half the matrix (Lower Left Triangle),
      // since it is symmetric, then set  Mat[ji] = Mat[ij]

      for (i = 0; i < m_nNumParams; i++)
	{
	  nBi = i;                                             // Batch index
	  if (nBi >= m_nNumBatches) nBi = nBi - m_nNumBatches;
	  fGi    = m_pfUha[nBi];
	  fDelGi = m_pfVha[i];
	  ij = i * m_nNumParams - 1;
	  ji = i - m_nNumParams;
	  for (j = 0; j <= i; j++)
	    {
	      ij++;
	      ji = ji + m_nNumParams;
	      nBj   = j;                                        // Batch index
	      if (nBj >= m_nNumBatches) nBj = nBj - m_nNumBatches;
	      fGj    = m_pfUha[nBj];
	      fDelGj = m_pfVha[j];
	      for (k = 0; k < nNumRefs; k++)
		{
		  poReflnI = m_poReflnlist->poGetRefln(pnIndex[k]);
		  fSigmaII = poReflnI->fGetSigmaI();
		  if (   (poReflnI->bIsSelected()) && (0.0 < fSigmaII)
		      && (nBi == poReflnI->nGetField(m_nFI_nBatchIdx)) )
		    {
		      fIntI    = poReflnI->fGetIntensity();
		      fWeightI = fCalcGetWeight(fIntI, fSigmaII, poReflnI);

		      // Reflection contributes to parameter i

		      if (0 == j)
			{
			  // Do this only once

			  m_pfBa[i] = m_pfBa[i] + fDelGi * fWeightI
			                          * fIh * (fIntI - fGi * fIh);
			}
		      for (m = 0; m < nNumRefs; m++)
			{
			  poReflnJ = m_poReflnlist->poGetRefln(pnIndex[m]);
			  fSigmaIJ = poReflnJ->fGetSigmaI();
			  if (   (poReflnJ->bIsSelected()) && (0.0 < fSigmaIJ)
			      && (nBj == poReflnJ->nGetField(m_nFI_nBatchIdx)) )
			    {
			      // Reflection contributes to param j

			      fIntJ    = poReflnJ->fGetIntensity();
			      fWeightJ = fCalcGetWeight(fIntJ, fSigmaIJ,
							poReflnJ);
			      fTemp = 0.0;
			      if (i == j)
				{
				  // On diagonal
				
				  fTemp  = fIh * fIh * fWeightI;
				}
			      fTemp = fTemp + fWeightI * fWeightJ
				              * (fIntI * fIntJ
						 - fIh * (fGi * fIntJ +
							  fGj * fIntI));
			      m_pfAab[ij] = m_pfAab[ij]
				               + fDelGi * fDelGj * fTemp;
			      m_pfAab[ji] = m_pfAab[ij];
			    }
			} // end m refln loop
		    }
		} // end k refln loop
	    } // end j
	} // end i
    } // end fUh > 0

  if (6 < m_nVerbose)
    {
      // Print out the matrix

      cout << "m_pfAab: \n";
      vListMatMN(m_nNumParams, m_nNumParams, m_pfAab);
      cout << "m_pfBa: \n";
      vListMatMN(m_nNumParams, 1, m_pfBa);
#ifdef SSI_PC
      cout << flush;
#endif
    }

  // All done, so return to our regularly scheduled normal equations matrix

  return;
}

void
Cscalemerge::vRoughScale(void)
{
  // Compute and apply rough scaling factors between the different
  // batches of reflections, so that the starting shifts in the non-linear
  // least squares might be reasonable

  // For rough scaling, just compute average F in each data set

  int     i, j;
  float   fSigmaI, fInt;
  int     nNumRefs = m_poReflnlist->nGetNumReflns();
  Crefln *poRefln;

  // Zero the sums
  for (i = 0; i < m_nNumBatches; i++)
    {
      m_pfUha[i] = 0.0;
      m_pfVha[i] = 0.0;
    }

  for (i = 0; i < nNumRefs; i++)
    {
      poRefln = m_poReflnlist->poGetRefln(i);
      fSigmaI = poRefln->fGetSigmaI();
      if ( (poRefln->bIsSelected()) && (0.0 < fSigmaI) )
	{
	  fInt = poRefln->fGetIntensity();
	  j = poRefln->nGetField(m_nFI_nBatchIdx);
	  m_pfUha[j] = m_pfUha[j] + fInt;
          m_pfVha[j] = m_pfVha[j] + 1.0f;
	}
    }
}

void
Cscalemerge::vListOverlap(void)
{
  // List matrix of overlaps between scaling batches
  // The matrix may be large, so list out in columns of 10

  int i, j, k;

  int  *pnTemp;
  char *pcDash6 = "------";

  printf("\n\nOverlaps among scaling batches");
  for (k = 0; k < m_nNumBatches; k = k + 10)
    {
      printf("\n");

      // List batchnames across top

      printf("%6s--", pcDash6);
      for (j = k; (j < m_nNumBatches) && (j < k + 10); j++)
	{
	  printf(pcDash6);
	}
      printf("--%6s\n", pcDash6);

      printf("%5s | ", "Batch");
      for (j = k; (j < m_nNumBatches) && (j < k + 10); j++)
	{
//	  printf("%6s", (const char *)m_psBatchNames[j]);
	  printf("%6s", m_psBatchNames[j].string());
	}
      printf ("\n");
      printf("%6s--", pcDash6);
      for (j = k; (j < m_nNumBatches) && (j < k + 10); j++)
	{
	  printf(pcDash6);
	}
      printf("--%6s\n", pcDash6);

      for (i = 0; i < m_nNumBatches; i++)
	{
	  // List batchnames on left side

//	  printf("%5s | ", (const char *) m_psBatchNames[i]);
	  printf("%5s | ", m_psBatchNames[i].string());

	  pnTemp = &m_pnNumOverlap[i * m_nNumBatches] + k;

	  for (j = k; (j < m_nNumBatches) && (j < k + 10); j++)
	    {
	      printf("%6d", *pnTemp++);
	    } // end j (column)

	  // List batchnames on right side

//	  printf(" |%5s\n", (const char*) m_psBatchNames[i]);
	  printf(" |%5s\n", m_psBatchNames[i].string());
	} // end i (row)

      // List batchnames across bottom

      printf("%6s--", pcDash6);
      for (j = k; (j < m_nNumBatches) && (j < k + 10); j++)
	{
	  printf(pcDash6);
	}
      printf("--%6s\n", pcDash6);

      printf("%5s | ", "Batch");
      for (j = k; (j < m_nNumBatches) && (j < k + 10); j++)
	{
//	  printf("%6s", (const char *)m_psBatchNames[j]);
	  printf("%6s", m_psBatchNames[j].string());
	}
      printf ("\n");
      printf("%6s--", pcDash6);
      for (j = k; (j < m_nNumBatches) && (j < k + 10); j++)
	{
	  printf(pcDash6);
	}
      printf("--%6s\n", pcDash6);
    } // end k (group of 10)
  printf("\n");
}

void
Cscalemerge::vInitStats(tagStats *ptStats)
{
  // Initialize a tagStats structure to all zeroes

  int i;
  ptStats->nNumRefs    = 0;
  ptStats->nNumRejects = 0;
  ptStats->nNumUsed    = 0;
  ptStats->nNumMults   = 0;
  for (i = 0; i < 10; i++)
    {
      ptStats->anNumInBinMults[i] = 0;
      ptStats->anNumInBin[i] = 0;
    }
  ptStats->nNumSingles = 0;
  ptStats->dNumer      = 0.0;
  ptStats->dDenom      = 0.0;
  ptStats->dSumChiSq   = 0.0;
  ptStats->dWtNumer    = 0.0;
  ptStats->dWtDenom    = 0.0;
  ptStats->dSumIoverSig= 0.0;
  ptStats->dSumChiSq   = 0.0;
  ptStats->dRmerge     = 0.0;	
}

void
Cscalemerge::vListStats(void)
{
  // List all statistics?

  int i;
  tagStats *ptStats;
  tagStats *ptSCumul;
  float    *pfSlope;
  float     fLow;
  float     fResoLow, fResoHigh;
  int       nTemp;
  char      cTemp1, cTemp2;
  char     *pcLine =
    "--------------------------------------------------------------------------------\n";


  printf("\n\nIn the tables below Rmerge is defined as:\n");
  printf("    Rmerge = Sum Sum |Ihi - <Ih>| / Sum Sum <Ih>\n");
  printf("              h   i                  h   i\n");
  printf("    where Ihi is the ith used observation for unique hkl h,\n");
  printf("    and  <Ih> is the mean intensity for unique hkl h.\n");

  // List batch statistics
  
  // For batch reduced chi-squares, we need to get the total number of
  // multiply-measured unique reflections as well as the total number of
  // contributors used in the multiply-measured unique reflections.
  // To make things work out, we will apportion the number of multiply-measured
  // unique reflections to each batch according to the number of contributing
  // reflections in the batch.

  int nNumMults = 0;
  int nNumUsed  = 0;
  ptStats  =  m_ptStatsReso;
  for (i = 0; i < m_nNumBinsReso; i++, ptStats++)
    {
      nNumMults = nNumMults + ptStats->nNumMults;
      nNumUsed  = nNumUsed  + ptStats->nNumUsed;
    }
  //  cout << "nNumMults: " << nNumMults << endl << flush;
  float fTempBatch;
  fTempBatch = (float)nNumMults / (float)nNumUsed;
  //  cout << "fTempBatch: " << fTempBatch << endl;  

  printf("\n\nRmerge vs Batch\n");
  printf(pcLine);
  printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
	 "Batch",  "Average", "Num", " Num", " Num", "    Num", "<I/ ",
	 "ChiSq", "Rmerge", "Rmerge");
  printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
	 " name", " counts", "obs", "rejs", "ovlps", "single", "sig>",
	 "    ", " batch", " cumul");

  printf(pcLine);
  ptStats  =  m_ptStatsBatch;
  ptSCumul = &m_ptStatsBatch[m_nNumBatches];
  vInitStats(ptSCumul);                          // Initialize cumulative stats
  for (i = 0; i < m_nNumBatches; i++, ptStats++)
    {
      ptSCumul->nNumRefs    = ptSCumul->nNumRefs    + ptStats->nNumRefs;
      ptSCumul->nNumRejects = ptSCumul->nNumRejects + ptStats->nNumRejects;
      ptSCumul->nNumUsed    = ptSCumul->nNumUsed    + ptStats->nNumUsed;
      ptSCumul->nNumMults   = ptSCumul->nNumMults   + ptStats->nNumMults;
      ptSCumul->nNumSingles = ptSCumul->nNumSingles + ptStats->nNumSingles;
      ptSCumul->dNumer      = ptSCumul->dNumer      + ptStats->dNumer;
      ptSCumul->dDenom      = ptSCumul->dDenom      + ptStats->dDenom;
      ptSCumul->dWtNumer    = ptSCumul->dWtNumer    + ptStats->dWtNumer;
      ptSCumul->dWtDenom    = ptSCumul->dWtDenom    + ptStats->dWtDenom;
      ptSCumul->dSumIoverSig= ptSCumul->dSumIoverSig+ ptStats->dSumIoverSig;
      ptSCumul->dSumChiSq   = ptSCumul->dSumChiSq   + ptStats->dSumChiSq;
      nTemp                 = ptStats->nNumSingles  + ptStats->nNumMults;
      printf("%14s %7.0f %7d %5d %7d %7d ",
//	     (const char *) m_psBatchNames[i],
	     m_psBatchNames[i].string(),
	     ptStats->dDenom / ptStats->nNumUsed,
	     ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
	     ptStats->nNumSingles);
      if (0.0 != ptStats->dDenom)
	{
	  printf("%6.1f %6.2f %6.3f ",
		 ptStats->dSumIoverSig / (float)nTemp,
		 ptStats->dSumChiSq
//+1-May
//		     / (float)ptStats->nNumMults,
//		     / ((float)ptStats->nNumUsed - (float)ptStats->nNumMults),
		 / ((float)ptStats->nNumUsed - (fTempBatch * (float)ptStats->nNumUsed)),
		 ptStats->dNumer / ptStats->dDenom);
/*
      printf("\n%f  %d  %d  : SumChiSq, nNumUsed, m_nNumBatches\n", ptStats->dSumChiSq,
	     ptStats->nNumUsed,
	     m_nNumBatches);
*/
//+1-May

	}
      else
	printf("%6s %6s %6s ", "---", "---", "---");
      if (0.0 != ptSCumul->dDenom)
	printf("%6.3f\n", ptSCumul->dNumer / ptSCumul->dDenom);
      else
	printf("%6s\n", "---");
    }
  printf(pcLine);
  printf("%14s %7.0f %7d %5d %7d %7d ",
	 "All batches",
	 ptStats->dDenom / ptStats->nNumUsed,
	 ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
	 ptStats->nNumSingles);
  if (0.0 != ptStats->dDenom)
    {
      nTemp = ptStats->nNumSingles + ptStats->nNumMults;
      printf("%6.1f ", ptStats->dSumIoverSig / (float)nTemp);


      printf("%6.2f ", ptStats->dSumChiSq
//+1-May            / (float) ptStats->nNumMults);
//	     / ((float)ptStats->nNumUsed - (float)ptStats->nNumMults));
	     / ((float)(ptStats->nNumUsed - (fTempBatch * (float)ptStats->nNumUsed))));
      /*
      printf("\n%f  %d  %d  %f: SumChiSq, nNumUsed, m_nNumBatches, fTempBatch\n", 
	     ptStats->dSumChiSq,
	     ptStats->nNumUsed,
	     m_nNumBatches, fTempBatch);
      */
//-1-May

      printf("%6.3f %6.3f\n", ptStats->dNumer / ptStats->dDenom,
	     ptSCumul->dNumer / ptSCumul->dDenom);
    }
  else
    {
      printf("%6s %6s %6s %6s\n ", "---", "---", "---", "---");
    }

  // List intensity statistics

  printf("\n\nRmerge vs Intensity/SigmaI\n");
  printf(pcLine);
  printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
	 "Int/sigmaI ",  "Average", "Num", " Num", " Num", "   Num", "<I/ ",
	 "ChiSq", "Rmerge", "Rmerge");
  printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
	 "   range   ", " counts", "obs", "rejs", "ovlps", "mults", "sig>",
	 "norm", " shell", " cumul");

  printf(pcLine);
  ptStats  = &m_ptStatsInt[m_nNumBinsInt-1];
  ptSCumul = &m_ptStatsInt[m_nNumBinsInt];
  vInitStats(ptSCumul);                          // Initialize cumulative stats
  fLow     =  m_a2fRangeInt[1];
  pfSlope  = &m_fSlopeInt;
  for (i = m_nNumBinsInt-1; i >= 0; i--, ptStats--)
    {
      ptSCumul->nNumRefs    = ptSCumul->nNumRefs    + ptStats->nNumRefs;
      ptSCumul->nNumRejects = ptSCumul->nNumRejects + ptStats->nNumRejects;
      ptSCumul->nNumUsed    = ptSCumul->nNumUsed    + ptStats->nNumUsed;
      ptSCumul->nNumMults   = ptSCumul->nNumMults   + ptStats->nNumMults;
      ptSCumul->nNumSingles = ptSCumul->nNumSingles + ptStats->nNumSingles;
      ptSCumul->dNumer      = ptSCumul->dNumer      + ptStats->dNumer;
      ptSCumul->dDenom      = ptSCumul->dDenom      + ptStats->dDenom;
      ptSCumul->dWtNumer    = ptSCumul->dWtNumer    + ptStats->dWtNumer;
      ptSCumul->dWtDenom    = ptSCumul->dWtDenom    + ptStats->dWtDenom;
      ptSCumul->dSumIoverSig= ptSCumul->dSumIoverSig+ ptStats->dSumIoverSig;
      ptSCumul->dSumChiSq   = ptSCumul->dSumChiSq   + ptStats->dSumChiSq;
      nTemp                 = ptStats->nNumSingles  + ptStats->nNumMults;
      cTemp1 = ' ';
      cTemp2 = ' ';
      if ((m_nNumBinsInt-1) == i) cTemp1 = '>';
      if (0 == i)                 cTemp2 = '<';
      if ( (0.0 != ptStats->dDenom) && (0 != ptStats->nNumUsed) )
	{
	  printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d %7d %6.1f %6.2f %6.3f ",
		 cTemp2, fLow - *pfSlope, cTemp1, fLow,
		 ptStats->dDenom / (float) ptStats->nNumUsed,
		 ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
		 ptStats->nNumMults,
		 ptStats->dSumIoverSig / (float)nTemp,
//		 ptStats->dSumIoverSig / (float)ptStats->nNumMults,

		 ptStats->dSumChiSq
//                     / (float) ptStats->nNumMults,
                       / ((float)ptStats->nNumUsed - (float)ptStats->nNumMults),
//		 ptStats->dSumChiSq,
//		 ptStats->dWtNumer / ptStats->dWtDenom,
		 ptStats->dNumer / ptStats->dDenom);
	}
      else
	{
	  printf("   %c%2.0f -  %c%2.0f%8s%8d%7d %7d %7d %6s %6s %6s ",
		 cTemp2, fLow - *pfSlope, cTemp1, fLow,
		 "---",
		 ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
		 ptStats->nNumMults, "---", "---", "---");
	}
      if (0.0 != ptSCumul->dDenom)
	{
	  printf("%6.3f\n", ptSCumul->dNumer / ptSCumul->dDenom);
	}
      else
	{
	  printf("%6s\n", "---");
	}
      fLow = fLow - *pfSlope;
    }

  cTemp1 = '>';
  cTemp2 = '<';
  ptStats = ptSCumul;
  printf(pcLine);
  printf("   %c%2.0f -  %c%2.0f%8.0f%8d%7d %7d %7d ",
	 cTemp2, m_a2fRangeInt[0], cTemp1, m_a2fRangeInt[1],
	 ptStats->dDenom / ptStats->nNumUsed,
	 ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
	 ptStats->nNumMults);
  if (0.0 != ptStats->dDenom)
    {
      nTemp = ptStats->nNumSingles + ptStats->nNumMults;
      printf("%6.1f ", ptStats->dSumIoverSig / (float)nTemp);
//      printf("%6.1f ", ptStats->dSumIoverSig / (float)ptStats->nNumMults);

      printf("%6.2f ", ptStats->dSumChiSq
//                     / (float) ptStats->nNumMults);
                       / ((float)ptStats->nNumUsed - (float)ptStats->nNumMults));
//      printf("%6.2f ", ptStats->dSumChiSq);
//      printf("%6.2f ", ptStats->dWtNumer / ptStats->dWtDenom);
      printf("%6.3f %6.3f\n", ptStats->dNumer / ptStats->dDenom,
	     ptSCumul->dNumer / ptSCumul->dDenom);
      /*
      printf("\n%f  %d  %d  : SumChiSq, nNumUsed, nNumMults\n", 
	     ptStats->dSumChiSq,
	     ptStats->nNumUsed,
	     ptStats->nNumMults);
      */
    }
  else
    {
      printf("%6s %6s %6s %6s\n", "---", "---", "---", "---");
    }

  // List resolution statistics

  printf("\n\nRmerge vs Resolution\n");
  printf(pcLine);
  printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
	 "Resolution  ", "Average", "Num", " Num", " Num", "   Num", "<I/ ",
	 "ChiSq", "Rmerge", "Rmerge");
  printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
	 "   range    ", " counts", "obs", "rejs", "ovlps", "mults", "sig>",
	 "norm", " shell", " cumul");

  printf(pcLine);
  ptStats  =  m_ptStatsReso;
  ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
  vInitStats(ptSCumul);                          // Initialize cumulative stats
  fLow     =  m_a2fRangeReso[0];
  pfSlope  = &m_fSlopeReso;
  for (i = 0; i < m_nNumBinsReso; i++, ptStats++)
    {
      fResoLow              = 1.0f / powf(fLow, 0.333333f);
      fResoHigh             = 1.0f / powf(fLow + *pfSlope, 0.333333f);
      ptSCumul->nNumRefs    = ptSCumul->nNumRefs    + ptStats->nNumRefs;
      ptSCumul->nNumRejects = ptSCumul->nNumRejects + ptStats->nNumRejects;
      ptSCumul->nNumUsed    = ptSCumul->nNumUsed    + ptStats->nNumUsed;
      ptSCumul->nNumMults   = ptSCumul->nNumMults   + ptStats->nNumMults;
      ptSCumul->nNumSingles = ptSCumul->nNumSingles + ptStats->nNumSingles;
      ptSCumul->dNumer      = ptSCumul->dNumer      + ptStats->dNumer;
      ptSCumul->dDenom      = ptSCumul->dDenom      + ptStats->dDenom;
      ptSCumul->dWtNumer    = ptSCumul->dWtNumer    + ptStats->dWtNumer;
      ptSCumul->dWtDenom    = ptSCumul->dWtDenom    + ptStats->dWtDenom;
      ptSCumul->dSumIoverSig= ptSCumul->dSumIoverSig+ ptStats->dSumIoverSig;
      ptSCumul->dSumChiSq   = ptSCumul->dSumChiSq   + ptStats->dSumChiSq;
      nTemp                 = ptStats->nNumSingles  + ptStats->nNumMults;
      if ( (0.0 != ptStats->dDenom) && (0 != ptStats->nNumUsed) )
	{
	  printf("%6.2f -%5.2f%8.0f%8d%7d %7d %7d %6.1f %6.2f %6.3f ",
		 fResoLow, fResoHigh,
		 //		 1./sqrtf(fLow), 1./sqrtf(fLow + *pfSlope),
		 ptStats->dDenom / (float) ptStats->nNumUsed,
		 ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
		 ptStats->nNumMults,
		 ptStats->dSumIoverSig / (float)nTemp,
//		 ptStats->dSumIoverSig / (float)ptStats->nNumMults,

		 ptStats->dSumChiSq
//+1-May
//                     / (float) ptStats->nNumMults,
                     / ((float)ptStats->nNumUsed - (float)ptStats->nNumMults)
//		 / (float)ptStats->nNumMults
//		 / ((float)ptStats->nNumUsed / ((float)ptStats->nNumMults-1.0))
		 ,
//-1-May

//		 ptStats->dSumChiSq,
//		 ptStats->dWtNumer / ptStats->dWtDenom,
		 ptStats->dNumer / ptStats->dDenom);
	}
      else
	{
	  printf("%6.2f -%5.2f%8s%8d%7d %7d %7d %6s %6s %6s ",
		 fResoLow, fResoHigh,
		 "---",
		 ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
		 ptStats->nNumMults, "---", "---", "---");
	}
      if (0.0 != ptSCumul->dDenom)
	printf("%6.3f\n", ptSCumul->dNumer / ptSCumul->dDenom);
      else
	printf("%6s\n", "---");
      fLow = fLow + *pfSlope;
    }
  printf(pcLine);
  printf("%6.2f -%5.2f%8.0f%8d%7d %7d %7d ",
	 1.0f / powf(m_a2fRangeReso[0], 0.333333f),
	 1.0f / powf(m_a2fRangeReso[1], 0.333333f),
	 ptStats->dDenom / ptStats->nNumUsed,
	 ptStats->nNumRefs, ptStats->nNumRejects, ptStats->nNumUsed,
	 ptStats->nNumMults);
  if (0.0 != ptStats->dDenom)
    {
      nTemp  = ptStats->nNumSingles  + ptStats->nNumMults;
      printf("%6.1f ", ptStats->dSumIoverSig / (float)nTemp);
//      printf("%6.1f ", ptStats->dSumIoverSig / (float)ptStats->nNumMults);

      printf("%6.2f ", ptStats->dSumChiSq
//                     / (float) ptStats->nNumMults);
                       / ((float)ptStats->nNumUsed - (float)ptStats->nNumMults));
//      printf("%6.2f ", ptStats->dSumChiSq);

//      printf("%6.2f ", ptStats->dWtNumer / ptStats->dWtDenom);
      printf("%6.3f ", ptStats->dNumer / ptStats->dDenom);
      printf("%6.3f\n", ptSCumul->dNumer / ptSCumul->dDenom);

      /*
      printf("\n%f  %d  %d  : SumChiSq, nNumUsed, nNumMults\n", 
	     ptStats->dSumChiSq,
	     ptStats->nNumUsed,
	     ptStats->nNumMults);
      */
    }
  else
    {
      printf("%6s %6s %6s %6s\n", "---", "---", "---", "---");
    }
  fflush(stdout);
}

int
Cscalemerge::nListCompleteness(void)
{
  // Calculate and list completeness table

  int i;
  tagStats *ptStats;
  tagStats *ptSCumul;
  float    *pfSlope;
  float     fLow;
  float     fResoLow, fResoHigh;
  char     *pcLine =
    "--------------------------------------------------------------------------------\n";

  if (NULL == m_poReflnlistComp)
    {
      // Form a list of all possible hkl's for comparison to
      // You are stuck with the resolution range first given.

      m_poReflnlistComp = new Creflnlist(m_poCrystal,
					 1.0f/powf(m_a2fRangeReso[1], 0.33333f),
					 1.0f/powf(m_a2fRangeReso[0], 0.33333f),
					 m_nScaleAnom);
    }

  if (!m_poReflnlistComp->bIsAvailable())
    {
      return (-1);
    }

//  (void) m_poReflnlistComp->nWrite("complete1.ref");

  // List completeness statistics

  if (0 != m_nScaleAnom)
    {
      printf("\n\nIn the following completeness and redundancy tables, I+ and I-");
      printf("\nare treated as non-equivalent reflections.\n");
    }

  printf("\n\nCompleteness vs Resolution\n");
  printf(pcLine);
  printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
	 "Resolution  ", "  Calc", "Num", " Num", " Num", "   Num", "Num",
	 " Avg", "%Comp", "%Comp");
  printf("%13s %7s %7s %6s %7s %7s %6s %6s %6s %6s\n",
	 "   range    ", "unique", "obs", "rejs", "mults", "single",
	 "unique", "mult", "shell", "cumul");

  printf(pcLine);
  ptStats  =  m_ptStatsReso;
  ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
  vInitStats(ptSCumul);                          // Initialize cumulative stats
  fLow     =  m_a2fRangeReso[0];
  pfSlope  = &m_fSlopeReso;
  int nNumUniqueCalc, nNumUniqueFound;

  for (i = 0; i < m_nNumBinsReso; i++, ptStats++)
    {
      fResoLow              = 1.0f / powf(fLow, 0.333333f);
      fResoHigh             = 1.0f / powf(fLow + *pfSlope, 0.333333f);
      nNumUniqueCalc = m_poReflnlistComp->nCountUnique(fResoHigh, fResoLow,
						       m_nScaleAnom);
      nNumUniqueFound       = ptStats->nNumMults    + ptStats->nNumSingles;
      ptSCumul->nNumRefs    = ptSCumul->nNumRefs    + ptStats->nNumRefs;
      ptSCumul->nNumRejects = ptSCumul->nNumRejects + ptStats->nNumRejects;
      ptSCumul->nNumUsed    = ptSCumul->nNumUsed    + ptStats->nNumUsed;
      ptSCumul->nNumMults   = ptSCumul->nNumMults   + ptStats->nNumMults;
      ptSCumul->nNumSingles = ptSCumul->nNumSingles + ptStats->nNumSingles;
      ptSCumul->dNumer      = ptSCumul->dNumer      + ptStats->dNumer;
      ptSCumul->dDenom      = ptSCumul->dDenom      + (float)nNumUniqueCalc;
      ptSCumul->dWtNumer    = ptSCumul->dWtNumer    + ptStats->dWtNumer;
      ptSCumul->dWtDenom    = ptSCumul->dWtDenom    + ptStats->dWtDenom;
      ptSCumul->dSumIoverSig= ptSCumul->dSumIoverSig+ ptStats->dSumIoverSig;
      ptSCumul->dSumChiSq   = ptSCumul->dSumChiSq   + ptStats->dSumChiSq;

      if ( (0 != nNumUniqueCalc) && (0 != ptStats->nNumUsed)
	                         && (0 != nNumUniqueFound) )
	{
	  printf("%6.2f -%5.2f%8d%8d%7d%8d%8d%7d%7.2f%7.1f ",
		 fResoLow, fResoHigh,
		 nNumUniqueCalc,
		 ptStats->nNumRefs, ptStats->nNumRejects,
		 ptStats->nNumMults,
		 ptStats->nNumSingles,
		 nNumUniqueFound,
		 (float) (ptStats->nNumUsed + ptStats->nNumSingles)
		       / (float) nNumUniqueFound,
		 min(100., 100. * (float) nNumUniqueFound / (float) nNumUniqueCalc));
	}
      else if (0 < nNumUniqueCalc)
	{
	  printf("%6.2f -%5.2f%8d%8d%7d%8d%8d%7d%7s%7.1f ",
		 fResoLow, fResoHigh,
		 nNumUniqueCalc,
		 ptStats->nNumRefs, ptStats->nNumRejects,
		 ptStats->nNumMults, ptStats->nNumSingles,
		 nNumUniqueFound,
		 "---",
		 min(100., 100. * (float) nNumUniqueFound / (float) nNumUniqueCalc));
	}
      else
	{
	  printf("%6.2f -%5.2f%8d%8d%7d%8d%8d%7d%7s%7s ",
		 fResoLow, fResoHigh,
		 nNumUniqueCalc,
		 ptStats->nNumRefs, ptStats->nNumRejects,
		 ptStats->nNumMults, ptStats->nNumSingles,
		 nNumUniqueFound,
		 "---", "---");
	}
      if (0.0 != ptSCumul->dDenom)
	printf("%6.1f\n", min(100.0,
			      100.* (float)(  ptSCumul->nNumMults
					+ ptSCumul->nNumSingles)
                                       / ptSCumul->dDenom));
      else
	printf("%6s\n", "---");
      fLow = fLow + *pfSlope;
    }
  nNumUniqueFound = ptSCumul->nNumMults + ptSCumul->nNumSingles;
  printf(pcLine);
  printf("%6.2f -%5.2f%8d%8d%7d%8d%8d%7d ",
	 1.0f / powf(m_a2fRangeReso[0], 0.333333f),
	 1.0f / powf(m_a2fRangeReso[1], 0.333333f),
	 (int) ptStats->dDenom,
	 ptStats->nNumRefs, ptStats->nNumRejects,
	 ptStats->nNumMults, ptStats->nNumSingles, nNumUniqueFound);
  if (0.0 != ptStats->dDenom)
    {
      printf("%6.2f ", (float) (ptStats->nNumUsed + ptStats->nNumSingles)
                            / (float) nNumUniqueFound);
      printf("%6.1f ",  min(100., 100.* (float) nNumUniqueFound / ptSCumul->dDenom));
      printf("%6.1f\n", min(100., 100.* (float) nNumUniqueFound / ptSCumul->dDenom));
    }
  else
    {
      printf("%7s%7s%7s\n", "---", "---", "---");
    }
//+jwp
  printf("\n\nRedundancy vs Resolution\n");
  printf(pcLine);
  printf("%13s%8s%47s %5s %5s\n",
	 "Resolution  ", "  Calc",
	 "Percent of reflections measured N times, N = ", "%Comp", "%Comp");
  printf("%13s%8s%5s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	 "   range    ", "unique", "0", "1", "2", "3",
	 "4", "5-8", "9-12", ">12", "shell", "cumul" );

  printf(pcLine);
  ptStats  =  m_ptStatsReso;
  ptSCumul = &m_ptStatsReso[m_nNumBinsReso];
  vInitStats(ptSCumul);                          // Initialize cumulative stats
  fLow     =  m_a2fRangeReso[0];
  pfSlope  = &m_fSlopeReso;
  int j;
  int nNumMultBinsM1 = 7;
  for (i = 0; i < m_nNumBinsReso; i++, ptStats++)
    {
      fResoLow              = 1.0f / powf(fLow, 0.333333f);
      fResoHigh             = 1.0f / powf(fLow + *pfSlope, 0.333333f);
      nNumUniqueCalc = m_poReflnlistComp->nCountUnique(fResoHigh, fResoLow,
						       m_nScaleAnom);
      nNumUniqueFound        = ptStats->nNumMults    + ptStats->nNumSingles;
      ptStats->anNumInBinMults[0]   = nNumUniqueCalc-nNumUniqueFound;
      ptSCumul->dDenom             += (float)nNumUniqueCalc;
      ptSCumul->dNumer             += (float) nNumUniqueFound;
      printf("%6.2f -%5.2f%8d", fResoLow, fResoHigh, nNumUniqueCalc);
      ptSCumul->anNumInBinMults[nNumMultBinsM1]
	           += ptStats->anNumInBinMults[nNumMultBinsM1];
      for (j = 0; j < nNumMultBinsM1 + 1; j++)
	{
//+15-Mar-1999
	  if (j < nNumMultBinsM1)
//-15-Mar-1999
	    ptSCumul->anNumInBinMults[j] += ptStats->anNumInBinMults[j];

	  if (0 == j)
	    printf("%5.1f",
		   max(0.0, min(100.0,
				(float) ptStats->anNumInBinMults[j]
				/ (float) nNumUniqueCalc * 100.0))
		   );
	  else
	    printf("%6.1f",
		   max(0.0, min(100.0,
				(float) ptStats->anNumInBinMults[j]
				/ (float) nNumUniqueCalc * 100.0))
		   );
	}
      printf("%6.1f%6.1f",
	     max(0.0, min(100., 100. * (float) nNumUniqueFound / (float) nNumUniqueCalc)),
	     max(0.0, min(100., 100.* (float)ptSCumul->dNumer / ptSCumul->dDenom))
	     );


      printf("\n");
      fLow = fLow + *pfSlope;
    }
  ptStats = ptSCumul;
  printf(pcLine);
  printf("%6.2f -%5.2f%8d%5.1f",
	 1.0f / powf(m_a2fRangeReso[0], 0.333333f),
	 1.0f / powf(m_a2fRangeReso[1], 0.333333f),
	 (int) ptStats->dDenom,
	 (float)ptStats->anNumInBinMults[0] / (float) ptStats->dDenom * 100.0);
  for (j = 1; j < nNumMultBinsM1+1; j++)
    {
      printf("%6.1f",
	     (float)ptStats->anNumInBinMults[j] / (float) ptStats->dDenom * 100.0);
    }
  printf("%6.1f%6.1f",
	 max(0.0, min(100., 100.* (float)ptSCumul->dNumer / ptSCumul->dDenom)),
	 max(0.0, min(100., 100.* (float)ptSCumul->dNumer / ptSCumul->dDenom)));
  printf("\n");

  float fDenom = 0.01f * ptStats->dDenom;
  if (0.0 < fDenom)
    {
      // Form some cumulative statistics

      printf("\n");
      printf(pcLine);
      printf("%13s          %s", "Resolution  ",
	     "Percent of reflections measured AT LEAST N times, N =\n");
      printf("%13s%8s%5s %5s %5s %5s %5s %5s %5s %5s\n",
	 "   range    ", " ", " ",  "13", "9", "5",
	 "4", "3", "2", "1" );
      printf(pcLine);
      printf("%6.2f -%5.2f%13s",
	     1.0f / powf(m_a2fRangeReso[0], 0.333333f),
	     1.0f / powf(m_a2fRangeReso[1], 0.333333f),
	     " ");
      for (j = nNumMultBinsM1; j > 0; j--)
	{
	  ptStats->anNumInBinMults[j] += ptStats->anNumInBinMults[j+1];
	  printf("%6.1f", (float)ptStats->anNumInBinMults[j] / fDenom);
	}
      printf("\n");
//      printf(pcLine);
      printf("\n");
    }
  fflush(stdout);
//-jwp
  return (0);
}

void
Cscalemerge::vApplyWeights(void)
{
  // Apply command line weights to the fSigmaI fields of m_poReflnlist

  int i;
  float fWeight;
  Crefln *poRefln;

  for (i = 0; i < m_poReflnlist->nGetNumReflns(); i++)
    {
      poRefln = m_poReflnlist->poGetRefln(i);
      fWeight = fCalcGetWeight(poRefln->fGetIntensity(),
			       poRefln->fGetSigmaI(),
			       poRefln);
      fWeight = sqrtf(fWeight);
      if (0.0 != fWeight)
	fWeight = 1.0f / fWeight;
      if (0.0 >= poRefln->fGetSigmaI())
	fWeight = -fWeight;
      poRefln->vSetSigmaI(fWeight);
    }
}


void
Cscalemerge::vSetResolution(const float fResoMin, const float fResoMax)
{
  // Set the min and max resolution in Angstroms.  Min is closest to
  // beam, max is as far away from beam, thus Min has a higher value
  // in Angstroms than Max

  m_fResolutionMin = max(fResoMin, fResoMax);
  m_fResolutionMax = min(fResoMin, fResoMax);
}

