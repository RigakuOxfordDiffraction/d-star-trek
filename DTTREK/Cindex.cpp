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
// Cindex.cc            Initial author: J.W. Pflugrath           26-May-1995
//  This file contains the member functions of class Cindex which implements
//    the autoindexing encapsulation of d*TREK.
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
#include "Cindex.h"         // Class definition and prototypes
#include "dtrekdefs.h"

#include <string.h>

#ifdef SSI_PC
#include "CrclHelper.h"
#endif  // SSI_PC

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
using std::cerr;
#endif

//+Definitions, constants, and initialization of static member variables

Cstring Cindex::ms_sfA      = "fA";
Cstring Cindex::ms_sfB      = "fB";
Cstring Cindex::ms_sfC      = "fC";
Cstring Cindex::ms_sfAlpha  = "fAlpha";
Cstring Cindex::ms_sfBeta   = "fBeta";
Cstring Cindex::ms_sfGamma  = "fGamma";
Cstring Cindex::ms_sfVolume = "fVolume";
Cstring Cindex::ms_sfRot1   = "fRot1";
Cstring Cindex::ms_sfRot2   = "fRot2";
Cstring Cindex::ms_sfRot3   = "fRot3";
Cstring Cindex::ms_sfDiffLength = "fDiffLength";
Cstring Cindex::ms_ssLattice= "sLattice";
Cstring Cindex::ms_ssLattSymm= "sLattSymm";
Cstring Cindex::ms_snLattNum= "nLattNum";
Cstring Cindex::ms_snSortOrder= "nSortOrder";
Cstring Cindex::ms_snIndexSelect = "nIndexSelect";
Cstring Cindex::ms_snTwinID = "nTwinID";
Cstring Cindex::ms_snNumIndexed= "nNumIndexed";


Cstring Cindex::ms_sDtindexOptions = D_K_DtindexOptions;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cindex::Cindex()
{
 (void) nInitValues();
}

Cindex::Cindex(Cimage_header& oHeader)
{
 (void) nInitValues(oHeader);
}

Cindex::~Cindex()
{
  if (NULL != m_poDiffVecs)
    {
      delete m_poDiffVecs;
      m_poDiffVecs = NULL;
    }
  if (NULL != m_poSolutions)
    {
      delete m_poSolutions;
      m_poSolutions = NULL;
    }

  if (m_bNewDetNames && (NULL != m_psDetectorNames) )
    {
      delete [] m_psDetectorNames;
      m_psDetectorNames = NULL;
      m_bNewDetNames    = FALSE;
    }

  if (m_bNewDetector && (NULL != m_ppoDetector) )
    {
      for (int i = 0; i < m_nNumDetectors; i++)
        {
          delete *(m_ppoDetector+i);
          *(m_ppoDetector+i) = NULL;
        }
      delete [] m_ppoDetector;
      m_ppoDetector  = NULL;
      m_bNewDetector = FALSE;
    }

  if (m_bNewSource && (NULL != m_poSource) )
    {
      delete m_poSource;
      m_poSource   = NULL;
      m_bNewSource = FALSE;
  }


  if (m_bNewCrysGonio && (NULL != m_poCrysGonio) )
    {
      delete m_poCrysGonio;
      m_poCrysGonio   = NULL;
      m_bNewCrysGonio = FALSE;
  }

  if (m_bNewRotation && (NULL != m_poRotation) )
    {
      delete m_poRotation;
      m_poRotation   = NULL;
      m_bNewRotation = FALSE;
    }

  if (NULL != m_pnIndex)
    {
      delete [] m_pnIndex;
      m_pnIndex = NULL;
    }

  if (NULL != m_pnSoln)
    {
      delete [] m_pnSoln;
      m_pnSoln = NULL;
    }

  if (NULL != m_po3DTrick)
    {
      delete m_po3DTrick;
      m_po3DTrick = NULL;
    }

  m_nNumDetectors  = 0;
}

int
Cindex::nInitValues(void)
{
  int i;  // loop counter
  m_nVerbose         = 2;
  m_nNumDetectors    = 0;
  m_nWhichDetector   = 0;

  m_psDetectorNames  = NULL;
  m_poReflnlist      = new Creflnlist();  // Create the working refln list
  m_pnIndex          = NULL;
  m_pnSoln           = NULL;

  m_bNewReflnlist    = TRUE;
  m_ppoDetector      = NULL;
  m_poSource         = NULL;
  m_poCrysGonio      = NULL;
  m_poRotation       = NULL;
  m_poDiffVecs       = NULL;
  m_poSolutions      = NULL;

  m_po3DTrick        = new C3Ddata(1,1,1,0,0,0);
  m_bNewDetNames     = FALSE;
  m_bNewDetector     = FALSE;
  m_bNewSource       = FALSE;
  m_bNewCrysGonio    = FALSE;
  m_bNewRotation     = FALSE;

  m_bReIndex         = TRUE;

  m_bKnownCell        = FALSE;
  m_fSubPeak          = 1.0;

  m_fMinAllowedResidual      = 0.00f;
  m_fMinAllowedResidualStart = 0.94f;
  m_fMinAllowedResidualEnd   = 0.50f;
  m_fMinAllowedResidualDelta = 0.01f;
  m_fResolutionMin           = 99999.0f;
  m_fResolutionMax           = 0.00001f;

  m_a6fInputErrors[0] = 3.0;
  m_a6fInputErrors[1] = 3.0;
  m_a6fInputErrors[2] = 3.0;
  m_a6fInputErrors[3] = 3.0;
  m_a6fInputErrors[4] = 3.0;
  m_a6fInputErrors[5] = 3.0;
  m_a6fInputCell[0]   = 100.0;
  m_a6fInputCell[1]   = 100.0;
  m_a6fInputCell[2]   = 100.0;
  m_a6fInputCell[3]   =  90.0;
  m_a6fInputCell[4]   =  90.0;
  m_a6fInputCell[5]   =  90.0;
  m_fCellLengthMax    = 2000.0;

  m_nSolnPick1        = -1;
  m_nSolnPick2        = -1;
  m_nLattNum          = -1;
  m_fWavelength       = 1.0;
  m_fAngleEps         = 0.4f;          // Degrees
  m_fLengthEps        = 0.005f;        // Fraction of length
  m_fFourierMax       = 10000.0;
  m_fMMSeparationMin  = 0.65f;


  for (i = 0; i < 1025; i++)
    {
      m_afCos2PI[i] = m_fFourierMax * (float)cos(2.0 * Gs_dPI * (double)i / 1024.0);
    }

  // Initialize matrices and values required for cell reduction

  // RB We need to initialize the members that used to be statics in the member functions in the old code
  m_nAtom_DPSFFT = 0;
  m_bFirstTime_IndexCos = TRUE;
  m_bFirstTime_IndexDPS = TRUE;
  m_nJ_IndexCos = 0;
  m_nDepth_DPSVecRefine = 0;

  return (0);
}

int
Cindex::nInitValues(Cimage_header& oHeader)
{
  int    i;       // Loop counter
  int    nStat, nTemp;
  Cstring sTemp;

  (void) nInitValues();

  if (!oHeader.bIsAvailable())
    {
      cout << "Cindex ERROR: image header is not valid."
           << "  Cannot construct index!\n";
      nStat = 1;
    }
  else
    {
      // Try to get required information from the header

      nStat = 0;

      m_poCrysGonio   = new Cgoniometer(oHeader, Ccrystal::ms_sCrystalPrefix);
      m_bNewCrysGonio = TRUE;

      m_poRotation    = new Crotation(oHeader);
      m_bNewRotation  = TRUE;

      // Get the number of detectors

      nTemp = oHeader.nGetValue(Cdetector::ms_sDetectorNumber, &m_nNumDetectors);

      if (0 == nTemp)
        {
          // Get the detector names

          m_psDetectorNames = new Cstring [m_nNumDetectors];
          m_bNewDetNames = TRUE;
          nTemp = oHeader.nGetValue(Cdetector::ms_sDetectorNames,
                                    m_nNumDetectors, m_psDetectorNames);
          if (0 == nTemp)
            {
              m_ppoDetector = new Cdetector* [m_nNumDetectors];
              for (i = 0; i < m_nNumDetectors; i++)
                {
                  // We do not need the non-uniformity of response information

                  m_ppoDetector[i] = new Cdetector (oHeader,
                                                    m_psDetectorNames[i],
                                                    TRUE, FALSE);
                }
              m_bNewDetector  = TRUE;
            }
          else
            nStat++;
        }
      else
        {
          nStat++;
        }

      m_poSource      = new Csource(oHeader);

      m_poSource->vCalcGetS0(m_a3fS0);
      m_fWavelength = m_poSource->fGetWavelength();

      m_bNewSource    = TRUE;
    }
  return (0);
}

int
Cindex::nList(const int nFlag)
{
  int   nStat;
  float fTempW, fRMin, fRMax;
  float a3fTempS0[3];

  fTempW = 0.0;
  cout << "Index listing:\n";

  if (NULL != m_poCrysGonio)
    {
      if (0 == nFlag)
        (void) m_poCrysGonio->nList();
    }
  else
    {
      cout << "     Crystal goniometer not defined." << endl;
    }
  if (NULL != m_poSource)
    {
      if (0 == nFlag)
        (void) m_poSource->nList();
      fTempW = m_poSource->fGetWavelength();
      m_poSource->vCalcGetS0(&a3fTempS0[0]);
    }
  else
    {
      cout << "     Source not defined." << endl;
      fTempW = 0.0;
    }
  for (int i = 0; i < m_nNumDetectors; i++)
    {
      if (NULL != m_ppoDetector[i])
        {
          (void) m_ppoDetector[i]->nList();
          if (0.0 != fTempW)
            {
              nStat = m_ppoDetector[i]->nGetResolution(a3fTempS0, &fRMin, &fRMax);
              if (0 == nStat)
                {
                  cout << "DetResolution min: " << fRMin * fTempW << endl;
                  cout << "DetResolution max: " << fRMax * fTempW << endl;
                }
            }
        }
    else
      {
        cout << "     Detector " << i << " not defined." << endl;
      }
    }
  cout << "\nIndex resol min: " << m_fResolutionMin << '\n';
  cout << "Index resol max: " << m_fResolutionMax << '\n';

  return (0);
}

int
Cindex::nExpandReflnlist(void)
{
  // Add required fields to the working reflection list if they are not
  // already there.
  // The new fields are named by the static member variables of the
  //   Creflnlist class

  m_nFI_fA       = m_poReflnlist->nExpandGetField(ms_sfA);
  m_nFI_fB       = m_poReflnlist->nExpandGetField(ms_sfB);
  m_nFI_fC       = m_poReflnlist->nExpandGetField(ms_sfC);
  m_nFI_fAlpha   = m_poReflnlist->nExpandGetField(ms_sfAlpha);
  m_nFI_fBeta    = m_poReflnlist->nExpandGetField(ms_sfBeta);
  m_nFI_fGamma   = m_poReflnlist->nExpandGetField(ms_sfGamma);
  m_nFI_fVolume  = m_poReflnlist->nExpandGetField(ms_sfVolume);
  m_nFI_fRot1    = m_poReflnlist->nExpandGetField(ms_sfRot1);
  m_nFI_fRot2    = m_poReflnlist->nExpandGetField(ms_sfRot2);
  m_nFI_fRot3    = m_poReflnlist->nExpandGetField(ms_sfRot3);
  m_nFI_sLattice = m_poReflnlist->nExpandGetField(ms_ssLattice);
  m_nFI_sLattSymm= m_poReflnlist->nExpandGetField(ms_ssLattSymm);
  m_nFI_nLattNum = m_poReflnlist->nExpandGetField(ms_snLattNum);
  m_nFI_nSortOrder = m_poReflnlist->nExpandGetField(ms_snSortOrder);
  m_nFI_nIndexSelect = m_poReflnlist->nExpandGetField(ms_snIndexSelect);
  m_nFI_nTwinID = m_poReflnlist->nExpandGetField(ms_snTwinID);
  m_nFI_nNumIndexed = m_poReflnlist->nExpandGetField(ms_snNumIndexed);


  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_snDetNum);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_snNonunfFlag);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_snOrignlReflnNum);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsPx0);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsPx1);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsXmm);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsYmm);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsRotMid);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsRotWidth);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcXmm);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcYmm);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcRotMid);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcRotWidth);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfPartiality);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfDetResolution);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfRecipCoord0);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfRecipCoord1);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfRecipCoord2);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfRecipCoordD0);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfRecipCoordD1);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfRecipCoordD2);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfFloatH);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfFloatK);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfFloatL);
  (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_snDiffFreq);


  return (0);
}


int
Cindex::nSetup(void)
{
  // Setup before index?  This may not be needed

  return (0);
}


int
Cindex::nCalcSumDiffs(bool bDoDiffVecs,const int nMaxNumRef, const float fCellLengthMaxIn,int nTwinID)
{
  // Calculate sum and difference vectors of the reflections reciprocal lattice
  // vectors.

  int    i, j, m, n;         // Loop counters
  int    nRef;               // Number of reflections in input reflection list
  float fVecI[3],fVecJ[3];
  float fVecDI[3],fVecDJ[3];
  float fVec[3],fVecD[3];
  float fLength;
  int   nSum;
  float fSum;
  int   nUseNumRef;          // Number of most intense reflections to use
  Crefln *poReflnI, *poReflnJ;  // Convenience pointers
  float fCellLengthMax;

  if ((m_nFI_nTwinID == -1) || (m_nFI_nIndexSelect == -1)) {
      printf("TwinID or nIndexSelect fields not available!\n");
      return 1;
  };

  int *pnSort;

  //  It would be wise to limit the number of reflection at this point
  //  Otherwise there could be zillions of difference vectors

  // Sort reflection list on intensities

  nRef = m_poReflnlist->nGetNumReflns();  // Number of reflections in the list

  if (nRef < 5) {
      printf("ERROR:  Not enough reflections!\n");
      return 1;
  };


  pnSort = new int[nRef];

  m_poReflnlist->vSort(eReflnField_float_type, m_poReflnlist->m_nFI_fIntensity,
                       pnSort);


  nUseNumRef = min(nMaxNumRef, nRef);
  if (bDoDiffVecs)
      nUseNumRef = min(nMaxNumRef,nRef*nRef/2);
  fCellLengthMax = fCellLengthMaxIn;
  if (0.0 >= fCellLengthMax)
      fCellLengthMax = m_fCellLengthMax;
  
  
  if (NULL != m_poDiffVecs)
      delete m_poDiffVecs;
  
  m_poDiffVecs       = new Creflnlist(*m_poReflnlist);
  m_poDiffVecs->nExpand(nUseNumRef);
  
  // First expand the fields in the m_poDiffVecs object (if required!):
  
  m_nFI_fDiffLength = m_poDiffVecs->nExpandGetField(ms_sfDiffLength);
  
  Crefln oRefln(m_poDiffVecs);

  // Don't use the top 10% most intense reflections, because they may be
  // zingers or diffraction from salt or ice  

  if (!bDoDiffVecs) {
      for (i = (nRef*9)/10-1;i>=0;i--)
      {
          
          j = pnSort[i];
          poReflnI = m_poReflnlist->poGetRefln(j);
          if ((poReflnI->nGetField(m_nFI_nTwinID)==nTwinID) &&
              (poReflnI->nGetField(m_nFI_nIndexSelect)==1)) {
              
              fVecI[0] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoord0);
              fVecI[1] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoord1);
              fVecI[2] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoord2);
              
              oRefln.vSetH(i);
              oRefln.vSetK(i);
              oRefln.vSetL(1);
              
              fLength = fLenVec3D(fVecI);  // Length is 1/ detector_resolution!?
              oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord0, fVecI[0]);
              oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord1, fVecI[1]);
              oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord2, fVecI[2]);
              
              oRefln.vSetField(m_nFI_fDiffLength, fLength);
              oRefln.vSetField(m_poDiffVecs->m_nFI_nDiffFreq, 1);
              
              fVecI[0] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoordD0);
              fVecI[1] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoordD1);
              fVecI[2] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoordD2);
              oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoordD0, fVecI[0]);
              oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoordD1, fVecI[1]);
              oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoordD2, fVecI[2]);
              
              m_poDiffVecs->nInsert(&oRefln);
          };
          if (m_poDiffVecs->nGetNumReflns()>=nUseNumRef)
              break;
      }
  } else {
      printf("...calculating difference vectors...\n");
          
      // Create the working diff vec list:         
      
      
      // Next compute the sums and difference of the input reflections.
      // Compute average length, too.
      
      Crefln oRefln(m_poDiffVecs);
      
      fSum = 0.0;
      nSum = 0;

      int nAvgSkip;

      nAvgSkip = max(1,nRef*(nRef-1)/2/nUseNumRef);
      
      for (m = 0; (m < nRef) && (m_poDiffVecs->nGetNumReflns()<nUseNumRef); m++)
      {
          i = pnSort[nRef-m-1];
          poReflnI = m_poReflnlist->poGetRefln(i);
          // h = -1 means this reflection rejected
          fVecI[0] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoord0);
          fVecI[1] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoord1);
          fVecI[2] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoord2);
          fVecDI[0] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoordD0);
          fVecDI[1] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoordD1);
          fVecDI[2] = poReflnI->fGetField(m_poReflnlist->m_nFI_fRecipCoordD2);
          
          for (n = m+1; (n < nRef) && (m_poDiffVecs->nGetNumReflns()<nUseNumRef);n+=nAvgSkip)
          {
              j = pnSort[nRef-1-n];
              poReflnJ = m_poReflnlist->poGetRefln(j);
              // h = -1 means ...
              
              fVecJ[0] = poReflnJ->fGetField(m_poReflnlist->m_nFI_fRecipCoord0);
              fVecJ[1] = poReflnJ->fGetField(m_poReflnlist->m_nFI_fRecipCoord1);
              fVecJ[2] = poReflnJ->fGetField(m_poReflnlist->m_nFI_fRecipCoord2);
              fVecDJ[0] = poReflnJ->fGetField(m_poReflnlist->m_nFI_fRecipCoordD0);
              fVecDJ[1] = poReflnJ->fGetField(m_poReflnlist->m_nFI_fRecipCoordD1);
              fVecDJ[2] = poReflnJ->fGetField(m_poReflnlist->m_nFI_fRecipCoordD2);
              
              if (fDot3D(fVecI,fVecJ)<0.0) {
                  vAddVec3DVec3D(fVecI,fVecJ,fVec);
                  vAddVec3DVec3D(fVecDI,fVecDJ,fVecD);
                  
                  
                  fLength = fLenVec3D(fVec);
                  if (fLength > 1.0 / fCellLengthMax)
                  {
                      // Do not let too short vectors into the diff set!
                      
                      oRefln.vSetH(i);
                      oRefln.vSetK(j);
                      oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord0,
                          fVec[0]);
                      oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord1,
                          fVec[1]);
                      oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord2,
                          fVec[2]);
                      oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoordD0,
                          fVecD[0]);
                      oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoordD1,
                          fVecD[1]);
                      oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoordD2,
                          fVecD[2]);
                      oRefln.vSetField(m_nFI_fDiffLength, fLength);
                      oRefln.vSetField(m_poDiffVecs->m_nFI_nDiffFreq,
                          1);
                      
                      
                      
                      m_poDiffVecs->nInsert(&oRefln);
                      fSum    = fSum + fLength;
                      nSum++;
                  } // end of fMul loop
                  
              } // end of n,j loop
          }
      }  // end of m, i loop
      
     
      if (0 == nSum)
      {
          if (NULL == m_poDiffVecs)
              m_poDiffVecs = new Creflnlist (); // Make empty list.
          return (1);                         // No Difference vectors!
      }
      
      fSum = fSum / (float) nSum;
      
  };
  
  
  delete [] pnSort;  
  
  
  
  
  
  // Do some things to fake out the rest of the procedures
  
  nSum = 1;
  (void) nCalcGetMaxCell(&fSum);
  if (0.0 != fSum) fSum = 1.0f / fSum;


  if (NULL != m_pnIndex) delete [] m_pnIndex;

  nRef      = m_poDiffVecs->nGetNumReflns();
  m_pnIndex = new int [nRef];              // Insure proper length

  m_poDiffVecs->vSort(eReflnField_float_type, m_nFI_fDiffLength, m_pnIndex);

  // Keep only the shortest 5000 reflections
  // First mark them for deletion, then delete them.
  // We have to do it this way because we have them
  // sorted on an index array.

  if (8 <= m_nVerbose )
    {
      Cstring sName;
      cout << "Enter name of PRE-CULLED difference vector file to write: ";
      cin >> sName;
      (void) m_poDiffVecs->nWrite(sName, m_pnIndex);
  }
  int nSaved = 0;

  if (5000 < nRef)
  {
      int *pnDelFlag;
      pnDelFlag = new int [nRef];

      for (i = 0; (i < nRef) && (5000 > nSaved); i++)
        {
          if (-1 != m_poDiffVecs->poGetRefln(m_pnIndex[i])->nGetH())
            {
              pnDelFlag[m_pnIndex[i]] = 0;
              nSaved++;
            }
          else
            pnDelFlag[m_pnIndex[i]] = -1;
        }
      for (; i < nRef; i++)  // i is already set above
        {
          pnDelFlag[m_pnIndex[i]] = -1;
        }
      m_poDiffVecs->nDelete(-1, pnDelFlag);
      delete [] pnDelFlag;

      // Since we deleted things, we need to re-sort

      if (NULL != m_pnIndex)
        delete [] m_pnIndex;
      nRef = m_poDiffVecs->nGetNumReflns();
      m_pnIndex = new int [nRef];         // Insure proper length
      m_poDiffVecs->vSort(eReflnField_float_type, m_nFI_fDiffLength, m_pnIndex);
    }

  // List out shortest 5 difference vector lengths, average them too.

  cout << endl;
  if (0 < nSum)
    {
      m_fDiffVecAvg = fSum / (float) nSum;
    }

  m_fDiffVecMax = m_fDiffVecAvg;

  if (8 <= m_nVerbose)
    {
      Cstring sName;
      
      cout << "Enter name of difference vector file to write: ";
      cin >> sName;
      (void) m_poDiffVecs->nWrite(sName, m_pnIndex);
    }

  // Save number of difference vectors in list, since we may add
  // to the list, but wish to distinguish between the old and newly added info

  m_nNumDiffVecs = m_poDiffVecs->nGetNumReflns();

  return (0);
}

int
Cindex::nGroupVecs(const float fTestIn)
{
  // Group the sum and difference vectors found in m_poDiffVecs
  // into clusters.  On return m_poDiffVecs holds the clusters.

  int i, j, k;
  int nRef;
  int nFreq;
  int nI;
  int nJ;
  float fVecI[3], fVecJ[3];
  float fLenI, fLenJ;
  float fTest;
  float fMult;
  int   nMult;
  bool  bMeetsCriteria;
  Crefln *poReflnI, *poReflnJ;   // Convenience pointers


  nRef = m_poDiffVecs->nGetNumReflns();

  if (nRef < 5) {
      printf("ERROR:  Not enough reflections!\n");
      return 1;
  };


  // Sort the m_poDiffVecs on frequency or length

  if (NULL != m_pnIndex)
    {
      delete [] m_pnIndex;
    }
  m_pnIndex = new int [nRef];
  m_poDiffVecs->vSort(eReflnField_float_type, m_nFI_fDiffLength, m_pnIndex);


/*
+jwp 4-Jun-1999  Re-ordering appears to be a bug, so don't do it
                 or really use the shortest vector 2 lines below

  // The re-order the list the other way because we want shortest with
  // highest sort index

  for (i = 0; i < nRef/2; i++)
    {
      j                 = m_pnIndex[i];
      m_pnIndex[i]        = m_pnIndex[nRef-i-1];
      m_pnIndex[nRef-i-1] = j;
    }
*/
  nRef = min(nRef, 5000);   // Limit grouping to 1st 5000 shortest diff vectors

  fTest = m_poDiffVecs->poGetRefln(m_pnIndex[0])->fGetField(
                                                            m_nFI_fDiffLength);
  fTest = fTest * 0.4f;
  if (.004 < fTest) fTest = 0.004f;   // .004 = 0.4 * 1./100.0

  if (fTest > fTestIn)
    fTest = fTestIn;

  if (3 <= m_nVerbose)
    {
      cout << "Grouping criteria used: " << fTest
           << "  Input: " << fTestIn << endl;
    }

  for (i = 0; i < nRef; i++)
    {
      nI = m_pnIndex[i];
      poReflnI = m_poDiffVecs->poGetRefln(nI);
      if (1 == poReflnI->nGetField(m_poDiffVecs->m_nFI_nDiffFreq))
        {
          // This reflection can be the start of a cluster

          fVecI[0] = poReflnI->fGetField(m_poDiffVecs->m_nFI_fRecipCoord0);
          fVecI[1] = poReflnI->fGetField(m_poDiffVecs->m_nFI_fRecipCoord1);
          fVecI[2] = poReflnI->fGetField(m_poDiffVecs->m_nFI_fRecipCoord2);
          fLenI    = poReflnI->fGetField(m_nFI_fDiffLength);
          nFreq    = 1;

          if (6 <= m_nVerbose)
            {
              cout << "Starting   " << fVecI[0] << ", " << fVecI[1]
                   << ", " << fVecI[2]
                   << endl;
            }
          for (j = i+1; (j < nRef); j++)
            {
              nJ = m_pnIndex[j];
              poReflnJ = m_poDiffVecs->poGetRefln(nJ);
              if (1 == poReflnJ->nGetField(m_poDiffVecs->m_nFI_nDiffFreq))
                {
                  fLenJ = poReflnJ->fGetField(m_nFI_fDiffLength);
                  fMult = fLenJ / fLenI;
                  nMult = (int) (fMult + (float) 0.5);
                  if ( fabs(((float) nMult * fLenI) - fLenJ) <= fTest)
                    {
                      // The lengths of the vectors are similar, so check if
                      // directions are the same by using a dot product

                      fVecJ[0] = poReflnJ->fGetField(
                                              m_poDiffVecs->m_nFI_fRecipCoord0);
                      fVecJ[1] = poReflnJ->fGetField(
                                              m_poDiffVecs->m_nFI_fRecipCoord1);
                      fVecJ[2] = poReflnJ->fGetField(
                                              m_poDiffVecs->m_nFI_fRecipCoord2);
                      bMeetsCriteria = TRUE;
                      for (k = 0; (k < 3) && bMeetsCriteria; k++)
                        {
                          bMeetsCriteria = fabs(fVecI[k] * (float)nMult - fVecJ[k])
                                           <= fTest;
                        }
                      // Ok, j is a member of cluster i, so average it into
                      // this cluster.
                      // First, get multiple of j's length versus i's length.
                      // We know it is an integer multiple because of the sin
                      //  test above.

                      if (bMeetsCriteria)
                        {
                          if (6 <= m_nVerbose)
                            {
                              cout << "Adding in  " << fVecJ[0] << ", "
                                   << fVecJ[1] << ", " << fVecJ[2]
                                   << " : " << fMult << " : " << nMult << endl;
                            }
                          for (k = 0; k < 3; k++)
                            {
                              // Add this vector to the group and
                              // Compute new average of the group.  By doing it
                              // this way, an outlier for the initial group
                              // member gets corrected (or an initial good
                              // value gets pulled away!)

                              fVecI[k] = ((fVecI[k] * (float) nFreq)
                                          +  fVecJ[k])
                                         / (float) (nFreq + nMult);

                              // Perhaps we should compute s.d. or spread of
                              // this group and weight the frequency by that?!
                            }
                          fLenI  = fLenVec3D(fVecI);
                          nFreq  = nFreq + nMult;

                          // Flag the jth diff vec as in a cluster
                        
                          poReflnJ->vSetField(m_poDiffVecs->m_nFI_nDiffFreq,
                                              -nI);
                        }
                    }
                }
            } // end of j loop

          poReflnI->vSetField(m_poDiffVecs->m_nFI_nDiffFreq, nFreq);
          poReflnI->vSetField(m_poDiffVecs->m_nFI_fRecipCoord0, fVecI[0]);
          poReflnI->vSetField(m_poDiffVecs->m_nFI_fRecipCoord1, fVecI[1]);
          poReflnI->vSetField(m_poDiffVecs->m_nFI_fRecipCoord2, fVecI[2]);
          poReflnI->vSetField(m_nFI_fDiffLength, fLenI);
          if (6 <= m_nVerbose)
            {
              cout << "Result     " << fVecI[0] << ", " << fVecI[1]
                   << ", " << fVecI[2] << ", " << fLenI
                   << endl;
            }
          //debug:
          poReflnI->vSetSigmaI(1.0f / fLenI * m_fWavelength);
        }
    } // end of i

  // Delete all reflections with a freq < 1 (they are in a group)
  // There is something to be said for keeping them though!

  nRef = m_poDiffVecs->nGetNumReflns();
  for (i = nRef-1; i >=0; i--)
    {
      if (1 > m_poDiffVecs->poGetRefln(i)->nGetField(
                                                m_poDiffVecs->m_nFI_nDiffFreq))
        {
          m_poDiffVecs->nDelete(i);
        }
    }

  // Since we deleted things, we need to re-sort on Frequency this time
  // Should we sort on Frequency or length?

  if (NULL != m_pnIndex) delete [] m_pnIndex;

  nRef      = m_poDiffVecs->nGetNumReflns();
  m_pnIndex = new int [nRef];      // Insure proper length
  m_poDiffVecs->vSort(eReflnField_int_type, m_poDiffVecs->m_nFI_nDiffFreq,
                      m_pnIndex);
  if (8 <= m_nVerbose)
    {
      Cstring sName;
      cout << "Enter name of grouped difference vector file to write: ";
      cin >> sName;
      (void) m_poDiffVecs->nWrite(sName, m_pnIndex);
  }

  return (0);
}


int
Cindex::nIndex(void)
{
  // Index the basis vectors!

  int   nStat;  // Local error flag
  int   i;      // Loop counter
  float fDet;

  // Make sure the basis vectors form a right-handed set of vectors.
  // The determinant will be positive.  If not, just negate the first
  // vector.  This will make things right(-handed).

  fDet = fInvMat3D(&m_a3x3fBasisVecs[0][0], &m_a3x3fInvSMat[0][0]);

  if (0.0 > fDet)
    {
      for (i = 0; i < 3; i++) m_a3x3fBasisVecs[0][i] = -m_a3x3fBasisVecs[0][i];
      fDet = fInvMat3D(&m_a3x3fBasisVecs[0][0], &m_a3x3fInvSMat[0][0]);
    }

  // If we know the cell constants or range of cell constants, calculate
  // some things here

  // Now set up the guessing loops
  // Relevant references are:
  // Busing & Levy (1967) Acta Cryst. 22, 457.
  // Sparks (1982) in Computational Crystallography (D. Sayre, ed.)
  //     Clarendon Press, Oxford.
  // Howard (1986) Proc. of EEC Cooperative Workshop on PSD Software, 89-100.

  // We have:
  //
  //1.  S       =  (s1 s2 s3), where s1, s2, s3 are the 3 basis vectors
  //    =           -- -- --         --  --  --
  //
  //2.  S       =  A * H , where H contains the indices of the basis vectors
  //    =          =   =         =                and A is the UB matrix
  //                                                  =        ==
  //     -1
  //3.  A  * S  =  H
  //    =    =     =
  //
  //     -1          -1       -1
  //4.  A       = (UB) = H * S
  //    =          ==    =   =
  //
  //     -1       a1
  //5.  A       = a2  and H = (h k l)  (h is all the h's, k is all the k's, etc)
  //    =         a3      =    - - -    -                 -
  //
  //                   -1                 -1               -1
  //6.  a1      =  h * S  ,    a2  = k * S  ,   a3  = l * S
  //    --         -   =       --    -   =      --    -   =
  //
  //7.  a1 dot s1 = integer! (from 3.), since every element of H is an integer
  //    --     --                                              =
  //
  //

  // Save number of difference vectors in list, since we are going to add
  // to the list, but wish to distinguish between the old and newly added info

  m_nNumDiffVecs = m_poDiffVecs->nGetNumReflns();

  // Sort difference vectors on frequency

  if (NULL != m_pnIndex) delete [] m_pnIndex;
  m_pnIndex = new int [m_nNumDiffVecs];  // Insure proper length
  m_poDiffVecs->vSort(eReflnField_int_type, m_poDiffVecs->m_nFI_nDiffFreq,
                      m_pnIndex);

  int   nIndexAbsMin  = 0;
  int   nIndexAbsMax  = 50;
  float fCellMin      = 2.0;
  float fCellMax      = m_fCellLengthMax;
  float fResidual     = 0.6f;
//  cout << "Enter residual cutoff: ";
//  cin >> fResidual;
  nStat = nFindIndices(fCellMin, fCellMax, nIndexAbsMin, nIndexAbsMax, fResidual);
  if (0 == nStat)
    {
      // Sort difference vectors on the integer residual found in the intensity

      delete [] m_pnIndex;
      m_pnIndex = new int [m_poDiffVecs->nGetNumReflns()];

      m_poDiffVecs->vSort(eReflnField_float_type,
                          m_poDiffVecs->m_nFI_fIntensity, m_pnIndex);

      // Here we have a set of possible aN's in
      //  m_poDiffVec[m_nNumDiffVecs] - m_poDiffVec[nNumReflns-1]
    }
  return (nStat);
}

int
Cindex::nFindIndices(const float fCellMin, const float fCellMax,
                     const int nIndexAbsMin, const int nIndexAbsMax,
                     const float fResCutoffIn)
{
  //                                            -1
  // Search for all possible solutions of  h * S   = a1,
  // where a1 dot s = integer

  int   nStat;
  int   i;                  // Loop counter
  float fI, fJ, fK, fL;     // Loop counters

  int   nSolutions;

  float fAbsMin, fAbsMax;
  bool  bIlessthanL;
  bool  bJlessthanI;
  bool  bKlessthanJ;

  float fVec[3];
  float fLenVec;
  float fRes;
  float fResCutoff;
  float fResMax;
  int   nMaxSolutions = 100;

  Crefln oRefln(m_poDiffVecs);

  nStat          = 0;
  fResMax        = -1.0;
  fResCutoff     = fResCutoffIn;
  nSolutions     = 0;

  fAbsMax = -1.0;
  for (i = 0; i < 3; i++)
    {
      fAbsMax = max(fAbsMax,fLenVec3D(&m_a3x3fBasisVecs[i][0]));
    }

  fAbsMax =  min(50,3 * (fAbsMax * fCellMax + 0.5f));

//  while ( (nSolutions <= 3) && (fResCutoff >= 0.50) ) {
    nSolutions = 0;
    fAbsMin = (float) nIndexAbsMin;
    if (0.0 >= fAbsMin) fAbsMin = 1.0;
//    fAbsMax = (float) nIndexAbsMax;

    if (fAbsMax < fAbsMin) fAbsMax = fAbsMin;

    // Arrange for an inward outward looping over fI, fJ, fK
    for (fL = fAbsMin; (fL < fAbsMax) && (nSolutions < nMaxSolutions); fL++)
      {
        for (fI = -fL; fI <= fL; fI++)
          {
            bIlessthanL = fabs(fI) < fL;
            for (fJ = -fL; fJ <= fL; fJ++)
              {
                bJlessthanI = bIlessthanL && (fabs(fJ) < fL);
                //for (fK = -fL; fK <= fL; fK++) {                    //?
                for (fK = 0; fK <= fL; fK++)
                  {
                    bKlessthanJ = bJlessthanI && (fabs(fK) < fL);
                    if (!bKlessthanJ)
                      {
                        fLenVec  = 0.0;
                        for (i = 0; i < 3; i++)
                          {
                            fVec[i] =      fI * m_a3x3fInvSMat[i][0]
                                        +  fJ * m_a3x3fInvSMat[i][1]
                                        +  fK * m_a3x3fInvSMat[i][2];
                            fLenVec = fLenVec  +  fVec[i] * fVec[i];
                          }
                        fLenVec = (float)sqrt((double)fLenVec) * m_fWavelength;
                        if ((fLenVec >= fCellMin) && (fLenVec <= fCellMax))
                          {
                            // Compute integer residual

                            if (   fIntegerResidual(fVec, m_nNumDiffVecs/7)
                                >= fResCutoff)
                              {
                                // Check the first 10% before checking
                                // the complete set

                                fRes = fIntegerResidual(fVec);
                                fResMax = max(fResMax,fRes);
                                if (fRes >= fResCutoff)
                                  {
                                    // Put in m_poDiffVecs

                                    oRefln.vSetH((int)fI);
                                    oRefln.vSetK((int)fJ);
                                    oRefln.vSetL((int)fK);
                                    oRefln.vSetField(m_nFI_fDiffLength,
                                                     fLenVec);
                                    oRefln.vSetField(
                                     m_poDiffVecs->m_nFI_fRecipCoord0, fVec[0]);
                                    oRefln.vSetField(
                                     m_poDiffVecs->m_nFI_fRecipCoord1, fVec[1]);
                                    oRefln.vSetField(
                                     m_poDiffVecs->m_nFI_fRecipCoord2, fVec[2]);
                                    oRefln.vSetIntensity(fRes);
                                    m_poDiffVecs->nInsert(&oRefln);
                                    nSolutions++;
                                  }
                              }
                          }
                      } // endif bK...
                  } // end of fK loop
              } // end of fJ loop
          } // end of fI loop
      } //end of fL loop

    fResCutoff = fResCutoff * 0.95f;

    if (5 <= m_nVerbose)
      {
        cout << "nFindIndices: fAbsMax       is: " << fAbsMax << endl;
        cout << "nFindIndices: Max residual  is: " << fResMax << endl;
        cout << "nFindIndices: Num solutions is: " << nSolutions << endl;
        cout << "nFindIndices:           fL  is: " << fL << endl;
      }
//  } // end of while

  if (3 > nSolutions)
    {
      // Delete bogus solutions added to the END of the list

      for (i = 0; i < nSolutions; i++)
        {
          m_poDiffVecs->nDelete(m_poDiffVecs->nGetNumReflns()-1);
        }
      nStat = 1;
    }
  return (nStat);
}

float
Cindex::fIntegerResidual(const float fVec[3], const int nTopRefs)
{
  //  Compute 1.0 - (SUM(cos(2PI*fVec dot every diffvec)) / nRefs)
  //  Best residuals will be close to 1.

  int   i, ni, nEnd;
  float fSum;
  float fVec2[3];
  float f1dot2;
  int nSum;

  if (0 >= m_nNumDiffVecs) return (-1.0);   // There are no vectors to use

  if (0 >= nTopRefs)
    nEnd = 0;                               // Use all available diff. vecs
  else
    {
      nEnd = m_nNumDiffVecs - min(nTopRefs, 5000); // Use the top nTopRef ones,
      nEnd = max(0, nEnd);                         //  but no more than 5000
    }                                              //  or the number available.

  fSum = 0.0;
  nSum = 0;
  for (ni = m_nNumDiffVecs - 1; ni >=  nEnd; ni--)
    {
      i = m_pnIndex[ni];
      fVec2[0] = m_poDiffVecs->poGetRefln(i)->fGetField(
                                            m_poDiffVecs->m_nFI_fRecipCoord0);
      fVec2[1] = m_poDiffVecs->poGetRefln(i)->fGetField(
                                            m_poDiffVecs->m_nFI_fRecipCoord1);
      fVec2[2] = m_poDiffVecs->poGetRefln(i)->fGetField(
                                            m_poDiffVecs->m_nFI_fRecipCoord2);

      f1dot2 = fDot3D(fVec, fVec2);

      // The following cosine function goes from -1 to +1

      fSum = fSum + (float)cos(2.0 * Gs_dPI * (double)f1dot2);
      nSum++;
    }

  // Shift residual from -1 - +1 to 0 to 1;

  fSum = 0.5f * (1.0f + fSum / (float) nSum);
  return (fSum);
}

int Cindex::nPutSolution(const int nSolnNum, Ccrystal *poCrystalInOut,double fResidual)
{
    int       nRef, nSoln;
    Crefln   *poSoln;   // Convenience pointer
    Ccrystal *poCrystal;
    
    nRef = m_poSolutions->nGetNumReflns();
    if (0 > nSolnNum)
        nSoln = nRef + nSolnNum;
    else
        nSoln = nSolnNum;
    if ( (nSoln >= nRef) || (0 > nSoln) )
    {
        cout << "Error in nLoadSolution, solution " << nRef
            << " out of bounds!\n";
        return (-1);
    }
    
    if (NULL == poCrystalInOut)
    {
        poCrystal = &m_oCrystal;
    }
    else
    {
        poCrystal = poCrystalInOut;
    }
    nRef   = m_pnSoln[nSoln];
    poSoln = m_poSolutions->poGetRefln(nRef);
    
    poSoln->vSetField(m_nFI_fA,(float) poCrystal->fGetCell(0));
    poSoln->vSetField(m_nFI_fB,(float) poCrystal->fGetCell(1));
    poSoln->vSetField(m_nFI_fC,(float) poCrystal->fGetCell(2));
    poSoln->vSetField(m_nFI_fAlpha,(float) poCrystal->fGetCell(3));
    poSoln->vSetField(m_nFI_fBeta,(float) poCrystal->fGetCell(4));
    poSoln->vSetField(m_nFI_fGamma,(float) poCrystal->fGetCell(5));
   
    poSoln->vSetField(m_nFI_fRot1,(float) poCrystal->fGetOrientAngle(0));
    poSoln->vSetField(m_nFI_fRot2,(float) poCrystal->fGetOrientAngle(1));
    poSoln->vSetField(m_nFI_fRot3,(float) poCrystal->fGetOrientAngle(2));      
    poSoln->vSetIntensity(fResidual);
    poSoln->vSetSigmaI(poCrystal->fCalcVolume());
    if (poCrystal->m_poSpacegroup->nGet()>=1)
        poSoln->vSetSigmaI(poCrystal->m_poSpacegroup->nGet());

    return 0;

};


int
Cindex::nLoadSolution(const int nSolnNum, Ccrystal *poCrystalInOut)
{
  // Load a given unit cell and unit parameters from the m_poSolutions object
  // into a Ccrystal object.

  int       nRef, nSoln;
  int       nStat;
  float     fDet;
  Crefln   *poSoln;   // Convenience pointer
  Ccrystal *poCrystal;

  nRef = m_poSolutions->nGetNumReflns();
  if (0 > nSolnNum)
    nSoln = nRef + nSolnNum;
  else
    nSoln = nSolnNum;
  if ( (nSoln >= nRef) || (0 > nSoln) )
    {
      cout << "Error in nLoadSolution, solution " << nRef
           << " out of bounds!\n";
      return (-1);
    }

  if (NULL == poCrystalInOut)
    {
      poCrystal = &m_oCrystal;
    }
  else
    {
      poCrystal = poCrystalInOut;
    }
  nRef   = m_pnSoln[nSoln];
  poSoln = m_poSolutions->poGetRefln(nRef);

  
  
  m_nLattNum                  = poSoln->nGetField(m_nFI_nLattNum);
  poCrystal->vSetCell(poSoln->fGetField(m_nFI_fA),
      poSoln->fGetField(m_nFI_fB),
      poSoln->fGetField(m_nFI_fC),
      poSoln->fGetField(m_nFI_fAlpha),
      poSoln->fGetField(m_nFI_fBeta),
      poSoln->fGetField(m_nFI_fGamma));
  poCrystal->vSetOrientAngles(poSoln->fGetField(m_nFI_fRot1),
      poSoln->fGetField(m_nFI_fRot2),
      poSoln->fGetField(m_nFI_fRot3));
  
  // Save the index into the m_poDiffVecs list that was used to create
  // the real space orientation matrix, so that these may be used to
  // re-create it elsewhere.
  
  m_nIndex0 = poSoln->nGetH();
  m_nIndex1 = poSoln->nGetK();
  m_nIndex2 = poSoln->nGetL();
  
  if (5 <= m_nVerbose)
  {
      cout << "indices in nLoad: " << m_nIndex0 << ", " << m_nIndex1 << ", " << m_nIndex2 << endl;
      cout << "a,b,c: " << poSoln->fGetField(m_nFI_fA) << ", " 
          << poSoln->fGetField(m_nFI_fB) << ", " 
          << poSoln->fGetField(m_nFI_fC) << endl;
  }
  if (!m_bKnownCell)
  {
      int nSG;
      nSG = (int)poSoln->fGetSigmaI();
      if (1 > nSG) nSG = 1;
      poCrystal->m_poSpacegroup->vSet(nSG);  // Could be bogus number!
  }
  else
  {
      poCrystal->m_poSpacegroup->vSet(1);
  }
  
  nStat = poCrystal->nCalcOrientMatrix(m_fWavelength);
  if (0 == nStat)
  {
      poCrystal->vGetOrientMatrix((float *)m_a3x3fUBMat);
      fDet = fInvMat3D((float *)m_a3x3fUBMat, (float *)m_a3x3fInvUBMat);
      if (0.0 == fDet) nStat = -1;
  }
  return (nStat);
}

int Cindex::nPruneSolutions(const int nNumToKeep) {
    int  nx;
    int  ny;
    Ccrystal oCrystal1;
    Ccrystal oCrystal2;
    CCellReduce   oReduce;

    double fVol1,fVol2,fErr1,fErr2;
    int    n1,n2;


    if (!m_poSolutions->nGetNumReflns())
        return 1;
    // Reject the top half worst residuals and/or residuals that are simply too large.Sort the solutions first based on the residualRun through the solutions, and 

    
    
    m_poSolutions->vSort(eReflnField_float_type, 0,m_pnSoln);    


    // Scan through the solutions, and prune those which have roughly the same residual,
    // but highly different volumes.  In this case, we choose the one with the lower volume.
    for (nx=0;nx<m_poSolutions->nGetNumReflns();nx++) {
        for (ny=nx+1;ny<m_poSolutions->nGetNumReflns();ny++) {
            fVol1 =  (*m_poSolutions)[m_pnSoln[nx]].fGetSigmaI();
            fErr1 = (*m_poSolutions)[m_pnSoln[nx]].fGetIntensity();
            fVol2 =  (*m_poSolutions)[m_pnSoln[ny]].fGetSigmaI();
            fErr2 = (*m_poSolutions)[m_pnSoln[ny]].fGetIntensity();       
            n1 = nx;
            n2 = ny;
            
            if ((fErr1 != 0.0) && (fErr2!=0.0)) {
                if (fVol1 < fVol2) {
                    std::swap(fVol1,fVol2);
                    std::swap(fErr1,fErr2);
                    std::swap(n1,n2);
                };
                
                if (fVol1/fVol2 > 1.8) {
                    if ((fErr1<fErr2) || (ABS(fErr1 - fErr2)/(0.5*(fErr1+fErr2))<0.2)) { 
                        // Erase the solution with the larger volume.                
                        m_poSolutions->nDeleteSorted(n1,m_pnSoln);
                        nx--;
                        break;
                    };
                };
            };
        };
    };


    // Prune solutions that have nearly the same vectors.  
    // Always favor the solution with the slightly lower residual.
    for (nx=0;nx<m_poSolutions->nGetNumReflns();nx++) {
        nLoadSolution(nx,&oCrystal1);
        for (ny=nx+1;ny<m_poSolutions->nGetNumReflns();ny++) {
            nLoadSolution(ny,&oCrystal2);
            
            // For each edge in oCrystal1, see if there is some edge nearly equal to it in oCrystal2.
            // Once that edge is found, mark it off.
            const double fMaxEdgeResid = 0.10;
            double a3x3fCrossIndex12[3][3];
            double a3x3fCrossIndex21[3][3];
            double a3x3fOrient1[3][3];
            double a3x3fOrient2[3][3];
            double a3x3fOrient2Inv[3][3];
            double a3x3fOrient1Inv[3][3];

            fVol1 =  (*m_poSolutions)[m_pnSoln[nx]].fGetSigmaI();
            fErr1 = (*m_poSolutions)[m_pnSoln[nx]].fGetIntensity();
            fVol2 =  (*m_poSolutions)[m_pnSoln[ny]].fGetSigmaI();
            fErr2 = (*m_poSolutions)[m_pnSoln[ny]].fGetIntensity();       

            oCrystal1.vGetOrientMatrix(&a3x3fOrient1[0][0]);
            oCrystal2.vGetOrientMatrix(&a3x3fOrient2[0][0]);
            fInvMat3D(&a3x3fOrient1[0][0],&a3x3fOrient1Inv[0][0]);
            fInvMat3D(&a3x3fOrient2[0][0],&a3x3fOrient2Inv[0][0]);
            vMulMat3DMat3D(a3x3fOrient2Inv,a3x3fOrient1,a3x3fCrossIndex12);
            vMulMat3DMat3D(a3x3fOrient1Inv,a3x3fOrient2,a3x3fCrossIndex21);

            for (n1 = 0; n1 < 3; n1++) {
                for (n2 = 0; n2 < 3; n2++) {
                    if (ABS(a3x3fCrossIndex12[n1][n2] - floor(0.5 + a3x3fCrossIndex12[n1][n2])) > fMaxEdgeResid)
                        break;
                    if (ABS(a3x3fCrossIndex21[n1][n2] - floor(0.5 + a3x3fCrossIndex21[n1][n2])) > fMaxEdgeResid)
                        break;
                };
                if (n2 != 3)
                    break;
            };
            if (n1 == 3) {
                // We found a duplicate cell. 
                // Since oCrystal1 has the higher residual, it should be removed.
                m_poSolutions->nDeleteSorted(nx,m_pnSoln);
                nx--;
                break;
            };
        };
    };


    ny = m_poSolutions->nGetNumReflns();
    for (nx=0;nx<ny-nNumToKeep;nx++) 
        m_poSolutions->nDeleteSorted(0,m_pnSoln);

    if ((*m_poSolutions)[0].fGetIntensity()==0.0) {
        // This is not a solution!!!!
        return 1;
    };
    return 0;
};

void
Cindex::vSetReflectionList(Creflnlist* poList) {
    if (m_poReflnlist != poList) {
        if (m_poReflnlist)
            delete m_poReflnlist;
        m_poReflnlist = poList;
    };
    
};



int Cindex::nReduceCell(float *pfCell, float *pfRot)
{

  // DON'T See Andrews & Bernstein (1988) Acta Cryst. A44, 1009-1018
    /*
    This matrix records 24 possible reductions. Do not confuse these with the 24 Buerger transforms!
    The integers 1,2,3 signify a,b, and c axes respectively.

    Col [0] =           Axis to be modified.
    Col [1] Col [3]  =  [3]* Col[1]  are added to Col [0]
    Col [2] Col [4]  =  [4]* Col[2]  are added to Col [0]

    We loop over all 24 of these operations.  When we can no longer reduce the cell, we stop.
    A cell reduction occurs when no edge get's shorter under any of the 24 operations.  Each edge
    has 8 possible ways of getting reduced.

    This cell works in non-oriented 3 dimensional space (9 quantities)

    */


   int nLCPre[][5]={
       {0,1,0,  1,0},
       {0,1,0,  -1,0},
       {1,0,0,  1,0},
       {1,0,0,  -1,0},
       {0,2,0,  1,0},
       {0,2,0,  -1,0},
       {2,0,0,  1,0},
       {2,0,0,  -1,0},
       {1,2,0,  1,0},
       {1,2,0,  -1,0},
       {2,1,0,  1,0},
       {2,1,0,  -1,0},
       {0,1,2,  1,1},
       {0,1,2,  1,-1},
       {0,1,2,  -1,1},
       {0,1,2,  -1,-1},
       {1,0,2,  1,1},
       {1,0,2,  -1,1},
       {1,0,2,  1,-1},
       {1,0,2,  -1,-1},
       {2,0,1,  1,1},
       {2,0,1,  -1,1},
       {2,0,1,  1,-1},
       {2,0,1,  -1,-1}
        };
   int nReductionLoopCt=0;
   bool bReductionmade;
   Ccrystal oCrystal;

   float a3fVecT[3];      // Temporary vector
   float a3fVecT2[3];     // Temporary vector.
   float a3x3fCellMat[3][3];
   int nx;
   float f0,f1;

   // Build the cell in a3x3fCellMat
   oCrystal.vSetCell(pfCell);
   oCrystal.nCalcGetRealMatrix(&(a3x3fCellMat[0][0]));
   vTranMat3D(a3x3fCellMat);
   
   do
     {
       bReductionmade = FALSE;
       for (nx = 0; nx < 24; nx++)
         {
           do 
             {
               f0 = fabs(
                      fDot3D(a3x3fCellMat[nLCPre[nx][0]],
                             a3x3fCellMat[nLCPre[nx][0]]));
      
               // [1]/[3]
               vCopyVec3D((float*)a3x3fCellMat[nLCPre[nx][1]], a3fVecT);
               vMulVec3DScalar(a3fVecT, nLCPre[nx][3], a3fVecT);
               vAddVec3DVec3D(a3fVecT, a3x3fCellMat[nLCPre[nx][0]], a3fVecT);
               // [2]/[4]
               vCopyVec3D((float*)a3x3fCellMat[nLCPre[nx][2]], a3fVecT2);
               vMulVec3DScalar(a3fVecT2, nLCPre[nx][4], a3fVecT2);
               vAddVec3DVec3D(a3fVecT, a3fVecT2, a3fVecT);
      
               f1 = fabs(fDot3D(a3fVecT, a3fVecT));

               // Did the edge get shorter?

               if ( f0 > f1 ) 
               { 
                   vCopyVec3D(a3fVecT, (float*)a3x3fCellMat[nLCPre[nx][0]]); 
                   bReductionmade = TRUE;
                   nReductionLoopCt++; 
#ifdef _DEBUG
                   if (500 < nReductionLoopCt)
                     cout << "f0, f1: " << f0 << ", " << f1 << endl << flush;
#endif
               }
             } while (f0 > f1);      
         }
       // Loop again if any reduction was made.
     } while ( bReductionmade && (5000 > nReductionLoopCt) );
   
   if (5000 <= nReductionLoopCt)
     {
       cout << "\n\n\nWARNING in Cindex::nReductionLoopCt: " 
           << nReductionLoopCt << endl << endl << endl << flush;
     }

   // Copy the cell over.
   pfCell[0] = fLenVec3D(a3x3fCellMat[0]);
   pfCell[1] = fLenVec3D(a3x3fCellMat[1]);
   pfCell[2] = fLenVec3D(a3x3fCellMat[2]);
   pfCell[3] = acos((double)fDot3D(a3x3fCellMat[1], a3x3fCellMat[2])
                    /(double)pfCell[1] / (double)pfCell[2])
               / Gs_dRADIANS_PER_DEGREE;
   pfCell[4] = acos((double)fDot3D(a3x3fCellMat[0], a3x3fCellMat[2])
                    / (double)pfCell[0] / (double)pfCell[2])
               / Gs_dRADIANS_PER_DEGREE;
   pfCell[5] = acos((double)fDot3D(a3x3fCellMat[0], a3x3fCellMat[1])
                    / (double)pfCell[0] / (double)pfCell[1]) 
              / Gs_dRADIANS_PER_DEGREE;  
  return (0);
}





int
Cindex::nIndexCos(const float fMinRadius, const float fMaxRadius,
                  const int nGrid, const int nDiffs, float *pfAvg)
{
  // Real space cosine Fourier transform indexing
  // Input:  m_poDiffVecs  (member) list of difference vectors sorted on
  //                              increasing length
  //         fMinRadius  (float)  Minimum radius in real cell to calculate
  //         fMaxRadius  (float)  Maximum radius in real cell (in Angstroms)
  //         nGrid       (int)    Number of grid points of output Fourier map
  //         nDiffs      (int)    Max number of difference vectors to use

  // Compute f(X) = sum over diffvecs[ cos 2pi(diffvec * X)]
  // where X is real space vector of (x,y,z)
  // Place f(X) in a 3Ddata object so it can be viewed and searched for
  // peaks.
  // TO DO:  nGrid should probably have a minimum of 50 and a maximum of ?.

  int i, j, k, m;   // Loop counters
  int nStat;
  int nVecs;
  Crefln *poRefln;

  // Save number of difference vectors in list, since we are going to add
  // to the list, but wish to distinguish between the old and newly added info

  m_nNumDiffVecs = m_poDiffVecs->nGetNumReflns();

  // Sort m_poDiffVecs on increasing difference vector length

  if (NULL != m_pnIndex)
    {
      delete [] m_pnIndex;
    }

  nVecs = min(nDiffs, m_nNumDiffVecs);

  m_pnIndex = new int [m_nNumDiffVecs];

  if (m_bFirstTime_IndexCos)
    {
//      s_bFirstTime = FALSE;
      for (m = 0; m < m_nNumDiffVecs; m++)
        {
          poRefln = m_poDiffVecs->poGetRefln(m);
          poRefln->vSetL(0);
        }

      // TODO: Choose the n shortest, the n most frequent, the n with smallest
      //       x, y, z values.

      m_poDiffVecs->vSort(eReflnField_float_type, m_nFI_fDiffLength, m_pnIndex);

      // Decrement the L value of the nVecs shortest ones

      for (m = 0; m < nVecs; m++)
        {
          j = m_pnIndex[m];
          poRefln = m_poDiffVecs->poGetRefln(j);
          poRefln->vSetL(poRefln->nGetL() - 2);
        }

      // Sort on frequency

      m_poDiffVecs->vSort(eReflnField_int_type,
                          m_poDiffVecs->m_nFI_nDiffFreq, m_pnIndex);

      // Decrement the L value of the nVecs most frequent ones

      for (m = m_nNumDiffVecs-1;
           m >= max(0, m_nNumDiffVecs - nVecs); m--)
        {
          j = m_pnIndex[m];
          poRefln = m_poDiffVecs->poGetRefln(j);
          poRefln->vSetL(poRefln->nGetL() - 1);
        }
/*
      // Sort on X value

      m_poDiffVecs->vSort(eReflnField_float_type,
                          m_poDiffVecs->m_nFI_fRecipCoord0, m_pnIndex);

      // Decrement the L value of the nVecs most shortest ones

      for (m = 0; m < nVecs; m++)
        {
          j = m_pnIndex[m];
          poRefln = m_poDiffVecs->poGetRefln(j);
          poRefln->vSetL(poRefln->nGetL() - 1);
        }

      // Sort on Y value

      m_poDiffVecs->vSort(eReflnField_float_type,
                          m_poDiffVecs->m_nFI_fRecipCoord1, m_pnIndex);

      // Decrement the L value of the nVecs most shortest ones

      for (m = 0; m < nVecs; m++)
        {
          j = m_pnIndex[m];
          poRefln = m_poDiffVecs->poGetRefln(j);
          poRefln->vSetL(poRefln->nGetL() - 1);
        }

      // Sort on Z value

      m_poDiffVecs->vSort(eReflnField_float_type,
                          m_poDiffVecs->m_nFI_fRecipCoord2, m_pnIndex);

      // Decrement the L value of the nVecs most shortest ones

      for (m = 0; m < nVecs; m++)
        {
          j = m_pnIndex[m];
          poRefln = m_poDiffVecs->poGetRefln(j);
          poRefln->vSetL(poRefln->nGetL() - 1);
        }
*/
    }

  // For speed, choose the nDiffs with smallest L values, so sort on nL

  m_poDiffVecs->vSort(eReflnField_int_type, m_poDiffVecs->m_nFI_nL,
                      m_pnIndex);


  // The code below is done in integer math for speed.  All difference
  // vectors are multiplied by scaling factor nScale to convert to integer

  int   nScale    = 8192; // * 2 * 2 * 2;    // MUST BE a power of 2
  int   nShift    = 3;    // + 1 + 1 + 1;       // (nScale / 1024)
  int   nScale3   = nScale - 1;
  float fStepSize = fMaxRadius / (float) (nGrid) / m_fWavelength

                    * (float) nScale;
//  int   anStep0[nVecs], *pnStep0; // Step size in Angstrom along 0th direction
//  int   anStep1[nVecs], *pnStep1; // Step size in Angstrom along 1st direction
//  int   anStep2[nVecs], *pnStep2; // Step size in Angstrom along 2nd direction
  int   *anStep0, *pnStep0;   // Step size in Angstrom along 0th direction
  int   *anStep1, *pnStep1;   // Step size in Angstrom along 1st direction
  int   *anStep2, *pnStep2;   // Step size in Angstrom along 2nd direction

  anStep0 = new int [nVecs];
  anStep2 = new int [nVecs];
  anStep1 = new int [nVecs];

  pnStep0 = anStep0;
  pnStep1 = anStep1;
  pnStep2 = anStep2;

  for (m = 0; m < nVecs; m++)
    {
      j = m_pnIndex[m];
      poRefln = m_poDiffVecs->poGetRefln(j);
      *pnStep0++ = (int) (poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord0)
                           * fStepSize  + 0.5);
      *pnStep1++ = (int) (poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord1)
                           * fStepSize  + 0.5);
      *pnStep2++ = (int) (poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord2)
                           * fStepSize  + 0.5);
    }

  float *pfCos2PI;
  pfCos2PI = &m_afCos2PI[0];      // Point to cosine lookup table

  int          ksq, jsq, ksqjsq, ksqjsqisq, nGridSq, nMinGridSq;
  nGridSq    = nGrid * nGrid;
  nMinGridSq = (int) (fMinRadius / fMaxRadius * (float)nGrid + 0.5);
  nMinGridSq = nMinGridSq * nMinGridSq;

  // Create a 3Ddata object to hold results have extra pixel
  // at start and end of each row,column,layer
  // so that limits do not need to be check in later peak searching.

  C3Ddata o3Ddata(nGrid+2, 2*nGrid+2, 2*nGrid+2, 0, 0, 0);
  o3Ddata.vZero();

  // Need some temporary vectors

//  int an2[nVecs],   *pn2;            // d2 . X2
//  int an21[nVecs],  *pn21;           // d2 . X2   +  d1 . X1
//  int an210[nVecs], *pn210;          // d2 . X2   +  d1 . X1  +  d0 . X0
  int *an2,   *pn2;            // d2 . X2
  int *an21,  *pn21;           // d2 . X2   +  d1 . X1
  int *an210, *pn210;          // d2 . X2   +  d1 . X1  +  d0 . X0
  int   nTemp;
  float fSum;

  an2   = new int [nVecs];
  an21  = new int [nVecs];
  an210 = new int [nVecs];

  // Initialize the 2 vector to -Step2 vector

  for (pn2 = an2, pnStep2 = anStep2; pn2 < &an2[nVecs]; pn2++)
    {
      *pn2 = -(nGrid+1) * (*pnStep2++);
    }

  for (k = -nGrid; k < nGrid; k++)
    {
      ksq = k * k;

      for (pn2 = an2, pn21 = an21, pnStep2 = anStep2, pnStep1 = anStep1;
           pn2 < &an2[nVecs];
           pn2++)
        {
          // Add Step2 vector to 2 vector
          // Initialize 21 vector to 2 vector less Step1 vector
        
          *pn2    = *pn2 + *pnStep2++;
          *pn21++ = *pn2 - (nGrid+1) * (*pnStep1++);
        }

      for (j = -nGrid; j < nGrid; j++)
        {
          jsq    = j * j;
          ksqjsq = ksq + jsq;
          if (ksqjsq > nGridSq)
            {
              if (0 > j)
                {
                  for (pn21 = an21, pnStep1 = anStep1;
                       pn21 < &an21[nVecs];
                       pn21++)
                    {
                      *pn21    = *pn21 + *pnStep1++;
                    }
                }
              else
                {
                  j = nGrid;
                }
            }
          else
            {
              // Add Step1 vector to 21 vector
              // Initialize 210 vector to 21 vector less Step0 vector

              for (pn21 = an21, pn210 = an210, pnStep1 = anStep1, pnStep0 = anStep0;
                   pn21 < &an21[nVecs];
                   pn21++)
                {
                  *pn21    = *pn21 + *pnStep1++;
                  *pn210++ = *pn21 - *pnStep0++;
                }

              for (i = 0; i < nGrid; i++)
                {
                  ksqjsqisq = ksqjsq + i * i;
                  if (ksqjsqisq > nGridSq)
                    {
                      if (0 > i)
                        {
                          for (pn210 = an210, pnStep0 = anStep0;
                               pn210 < &an210[nVecs];
                               pn210++)
                            {
                              *pn210  = *pn210 + *pnStep0++;
                            }
                        }
                      else
                        {
                          i = nGrid;
                        }
                    }
                  else if (ksqjsqisq < nMinGridSq)
                    {
                      // This real space coordinate does not meet the
                      // min radius criteria, so ...

                      // Add 0 vector to 210 vector

                      for (pn210 = an210, pnStep0 = anStep0;
                           pn210 < &an210[nVecs];
                           pn210++)
                        {
                          *pn210  = *pn210 + *pnStep0++;

                        }
                    }
                  else if (ksqjsqisq >= nMinGridSq)
                    {
                      // Ok, this real space coordinate meets the
                      // min and max radius criteria, so ...
                      // Add 0 vector to 210 vector and compute f(X)=fSum

                      fSum = 0.0;
                      for (pn210 = an210, pnStep0 = anStep0;
                           pn210 < &an210[nVecs];
                           pn210++)
                        {
                          *pn210  = *pn210 + *pnStep0++;

                          // Compute sum of cos2pi(d dot X)

                          nTemp   = *pn210;
                          if (0 > nTemp) nTemp = -nTemp;  // Symmetric around 0

                          // (nTemp modulo nScale) / (nScale/1024)
                          nTemp   = (nTemp & nScale3) >> nShift;
                          fSum    = fSum + *(pfCos2PI + nTemp);
                        }
                      fSum = fSum / nVecs;       // Normalize to m_fFourierMax
                      if (fSum < 0.0)
                        {
                          fSum = 0.0;
                        }

                      // Store in o3Ddata object with border of 1 pixel

                      o3Ddata.vSetValue(fSum, i+1, nGrid+j+1, nGrid+k+1);

                      // Possible optimizations:
                      // First, compute the sum for a portion of the
                      // difference vectors, and check it.
                    }
                } // end i
            } // end if ksqjsq
        } // end j
    } // end k

  // Set "origin" peak to m_fFourierMax+1 instead of m_fFourierMax
  // so it can be excluded in nGetPeaks

  o3Ddata.vSetValue((float)m_fFourierMax+1, 1,nGrid+1,nGrid+1);

// For debugging, write out o3Ddata as image for viewing with dtdisplay

  if (5 <= m_nVerbose)
    {
      Cimage oImage(&o3Ddata, 0, NULL);
      oImage.m_oHeader.nReplaceValue(Cimage_header::ms_sSaturatedValue,
                                   m_fFourierMax * 0.9);
      oImage.nWrite(Cstring("indexcos.img"));
    }

  fStepSize = fMaxRadius / (float) (nGrid);
  if (NULL == pfAvg)
    {
      fSum = 0.0;
      int nBeforePeakSearch = m_poDiffVecs->nGetNumReflns();
      do
        {
          nStat = nGetPeaks(&o3Ddata, fStepSize, &fSum);
        } while ( (0 != nStat) && (m_fFourierMax / 3.0 < fSum));

      if (0 != nStat)
        {
          // Get peaks with a lower min
          cout << "ERROR in nIndexCos, nGetUBMatrix2 failed!\n";
        }
    }
  else if (NULL != pfAvg)
    {
      *pfAvg = 0.0;
      *pfAvg = *pfAvg + o3Ddata.fGetValue(1,nGrid+2,nGrid+1);  // Like 0, 1, 0
      *pfAvg = *pfAvg + o3Ddata.fGetValue(2,nGrid+1,nGrid+1);  // Like 1, 0, 0
      *pfAvg = *pfAvg * 0.5f;
    }

  delete [] anStep0;
  delete [] anStep1;
  delete [] anStep2;
  delete [] an2;
  delete [] an21;
  delete [] an210;

  return (0);
}

int
Cindex::nIndexCos(const float *pfMinAngstrom, const float *pfMaxAngstrom,
                  const int nGrid, const int nDiffs, const int nFlag,
                  float *pfCOG)
{
  // Real space cosine Fourier transform indexing
  // Input:  m_poDiffVecs  (member) list of difference vectors sorted on
  //                              increasing length
  //        *pfMinAngstrom (float)Min Angstroms in real cell to calc (in 3 dirs)
  //        *pfMaxAngstrom (float)Maximum Angstroms in real cell
  //         nGrid       (int)    Number of grid points of output Fourier map
  //         nDiffs      (int)    Max number of difference vectors to use
  //         nFlag    = 0 if compute max position
  //                  = 1 if special call from nGetGrid()

  // Compute f(X) = sum over diffvecs[ cos 2pi(diffvec * X)]
  // where X is real space vector of (x,y,z)
  // Place f(X) in a 3Ddata object so it can be viewed and searched for
  // peaks.

  int   i, j, k, m;   // Loop counters
  int   nVecs;
  float fMaxAngstrom;
  float a3fCOG[3];

  fMaxAngstrom = pfMaxAngstrom[0] - pfMinAngstrom[0];
  fMaxAngstrom = max(fMaxAngstrom,(pfMaxAngstrom[1] - pfMinAngstrom[1]));
  fMaxAngstrom = max(fMaxAngstrom,(pfMaxAngstrom[2] - pfMinAngstrom[2]));

  // For speed, choose the nDiffs shortest difference vectors.
  // This code assumes that m_poDiffVecs has been sorted with index m_pnIndex.

  nVecs = min(nDiffs, m_nNumDiffVecs); // m_nNumDiffVecs comes from nIndexCos
                                   // routine and m_pnIndex had them sort there!

  // The code below is done in integer math for speed.  All difference
  // vectors are multiplied by scaling factor nScale to convert to integer

  int   nScale  = 8192; //  * 2 * 2 * 2;    // MUST BE a power of 2
  int   nShift  = 3; // + 1 + 1 + 1;       // (nScale / 1024)
  int   nScale3 = nScale - 1;

  float fStepSize = fMaxAngstrom / (float) (nGrid) / m_fWavelength
                    * (float) nScale;
//  int   anStep0[nVecs], *pnStep0;   // Step size in Angstrom along 0th direction
//  int   anStep1[nVecs], *pnStep1;   // Step size in Angstrom along 1st direction
//  int   anStep2[nVecs], *pnStep2;   // Step size in Angstrom along 2nd direction
  int   *anStep0, *pnStep0;   // Step size in Angstrom along 0th direction
  int   *anStep1, *pnStep1;   // Step size in Angstrom along 1st direction
  int   *anStep2, *pnStep2;   // Step size in Angstrom along 2nd direction

  anStep0 = new int [nVecs];
  anStep1 = new int [nVecs];
  anStep2 = new int [nVecs];

  pnStep0 = anStep0;
  pnStep1 = anStep1;
  pnStep2 = anStep2;
  for (m = 0; m < nVecs; m++)
    {
      j = m_pnIndex[m];
      *pnStep0++ = (int) (m_poDiffVecs->poGetRefln(j)->fGetField(
                                             m_poDiffVecs->m_nFI_fRecipCoord0)
                           * fStepSize  + 0.5);
      *pnStep1++ = (int) (m_poDiffVecs->poGetRefln(j)->fGetField(
                                             m_poDiffVecs->m_nFI_fRecipCoord1)
                           * fStepSize  + 0.5);
      *pnStep2++ = (int) (m_poDiffVecs->poGetRefln(j)->fGetField(
                                             m_poDiffVecs->m_nFI_fRecipCoord2)
                           * fStepSize  + 0.5);
    }

  float *pfCos2PI;
  pfCos2PI = &m_afCos2PI[0];          // Point to cosine lookup table

  // Create a 3Ddata object to hold results

  C3Ddata o3Ddata(nGrid+2, nGrid+2, nGrid+2, 0, 0, 0);
  o3Ddata.vZero(); // Zero the data array

  // Need some temporary vectors

//  int an2[nVecs],   *pn2;            // d2 . X2
//  int an21[nVecs],  *pn21;           // d2 . X2   +  d1 . X1
//  int an210[nVecs], *pn210;          // d2 . X2   +  d1 . X1  +  d0 . X0
  int *an2,   *pn2;            // d2 . X2
  int *an21,  *pn21;           // d2 . X2   +  d1 . X1
  int *an210, *pn210;          // d2 . X2   +  d1 . X1  +  d0 . X0
  int   nTemp;
  float fSum;

  an2   = new int [nVecs];
  an21  = new int [nVecs];
  an210 = new int [nVecs];

  // Calculate initial steps

  int anInitStep[3];
  for (m = 0; m < 3; m++)
    {
      // Find starting grid point

//  float fStepSize = fMaxAngstrom / (float) (nGrid) / m_fWavelength
//                    * (float) nScale;
      fSum =  pfMinAngstrom[m]  / m_fWavelength / fStepSize
                             * (float) nScale;
      if (fSum >= 0.0)
        {
          anInitStep[m] = (int) (fSum + 0.5);  // Round up
        }
      else
        {
          anInitStep[m] = (int) (fSum - 0.5);  // Round down
        }
    }

  // Initialize the 2 vector to anInitStep[2] * Step2 vector

  for (pn2 = an2, pnStep2 = anStep2; pn2 < &an2[nVecs]; pn2++)
    {
      *pn2 =  (anInitStep[2]-1) * *pnStep2++;
    }

  for (k = 0; k < nGrid; k++)
    {
      for (pn2 = an2, pn21 = an21, pnStep2 = anStep2, pnStep1 = anStep1;
           pn2 < &an2[nVecs];
           pn2++)
        {
          // Add Step2 vector to 2 vector
          // Initialize 21 vector to 2 vector + anInitStep[1] * Step1 vector
        
          *pn2    = *pn2 + *pnStep2++;
          *pn21++ = *pn2 + ((anInitStep[1]-1) * *pnStep1++);
        }

      for (j = 0; j < nGrid; j++)
        {
          // Add Step1 vector to 21 vector
          // Initialize 210 vector to 21 vector plus anInitStep[0] * Step0 vector

          for (pn21 = an21, pn210 = an210, pnStep1 = anStep1, pnStep0 = anStep0;
               pn21 < &an21[nVecs];
               pn21++)
            {
              *pn21    = *pn21 + *pnStep1++;
              *pn210++ = *pn21 + ((anInitStep[0]-1) *  *pnStep0++);
            }

          for (i = 0; i < nGrid; i++)
            {
              // Add 0 vector to 210 vector and compute f(X)=fSum

              fSum = 0.0;
              for (pn210 = an210, pnStep0 = anStep0;
                   pn210 < &an210[nVecs];
                   pn210++)
                {
                  *pn210  = *pn210 + *pnStep0++;
                
                  // Compute sum of cos2pi(d dot X)
                
                  nTemp   = *pn210;
                  if (nTemp < 0) nTemp = -nTemp;  // Symmetric around 0
                
                  // (nTemp modulo nScale) / (nScale/1024)
                  nTemp   = (nTemp & nScale3) >> nShift;
                  fSum    = fSum + *(pfCos2PI + nTemp);
                }
              fSum = fSum / nVecs;       // Normalize to 1 (or nScale)
              // Store in o3Ddata object with border of 1 pixel
              o3Ddata.vSetValue(fSum, i+1, j+1, k+1);
        
              // Possible optimizations:
              // First, compute the sum for a portion of the
              // difference vectors, and check it.
        
            } // end i
        } // end j
    } // end k

  float fInt;
  int a3nMax[3];
  if ( (0 == nFlag) && (NULL != pfCOG) )
    {
      // Get the Center of Gravity of the data in the o3Ddata object

//      o3Ddata.nCalcGetCOG(&pfCOG[0], &pfCOG[1], &pfCOG[2],
//                          &fInt, &fSigmaI);
//      cout << "nIndexCos: COG from nCalcGetCOG: "
//           << pfCOG[0] << ", "
//           << pfCOG[1] << ", "
//           << pfCOG[2] << '\n';

      o3Ddata.nCalcGetMax(&a3nMax[0], &a3nMax[1], &a3nMax[2], &fInt);

      o3Ddata.vSetValue(-fInt, a3nMax[0], a3nMax[1], a3nMax[2]); // For debugging

      // Recompute stepsize in units of Angstroms/gridpoint

      fStepSize =  fStepSize * m_fWavelength / (float) nScale;

      // Convert COG to Angstroms
      for (m = 0; m < 3; m++)
        {
          // Then return COG in Angstroms!
          // Something is wrong with this calculation!!!
/*
//  float fStepSize = fMaxAngstrom / (float) (nGrid) / m_fWavelength
//                    * (float) nScale;
//      anInitStep[m] = (int) (pfMinAngstrom[m] / fStepSize / m_fWavelength
//                             * (float) nScale + 0.5);
          afStartAngstrom[m] = (float) anInitStep[m] * fStepSize;
          a3fCOG[m] = (float) a3nMax[m];
          pfCOG[m]  = a3fCOG[m];
//          a3fCOG[m] = pfCOG[m]; // Save for later
//          pfCOG[m] = 11.0;// Temporary, put it at midpoint
          pfCOG[m] = afStartAngstrom[m] + ( (pfCOG[m]+1.0) * fStepSize) - 0.5;
*/
          a3fCOG[m] = (float) a3nMax[m] - 1.0f;
          pfCOG[m]  = pfMinAngstrom[m] +  ((float) (a3nMax[m]) - 1.0f) * fStepSize;
        }
    }
  else if ( (1 == nFlag) && (NULL != pfCOG) )
    {
      // nGrid MUST BE 4
/*      printf("nGrid is %d\n", nGrid);
      for (k = 1; k <= nGrid; k++)
      for (j = 1; j <= nGrid; j++)
      for (i = 1; i <= nGrid; i++)
        printf("i,j,k, value %d %d %d %f\n",
               i, j, k, o3Ddata.fGetValue(i, j, k));
*/
      // Set origin peak to 0 and then return max value
      o3Ddata.vSetValue((float)0.0, nGrid-1, nGrid-1, nGrid-1);
      o3Ddata.nCalcGetMax(&a3nMax[0], &a3nMax[1], &a3nMax[2], pfCOG);
    }

// For debugging, write out o3Ddata as image for viewing with dtdisplay

  if (8 <= m_nVerbose)
    {
      char   cNum[5];
      Cimage oImage(&o3Ddata, 0, NULL);
      (void) oImage.m_oHeader.nReplaceValue(Cimage_header::ms_sCog1, 3, a3fCOG);
      (void) oImage.m_oHeader.nReplaceValue(Cimage_header::ms_sFint, fInt);
      (void) oImage.m_oHeader.nReplaceValue(Cimage_header::ms_sCog2, 3, pfCOG);
      sprintf(cNum, "%04d", m_nJ_IndexCos);
      oImage.nWrite(Cstring("indexcosB") + cNum + ".img");
      m_nJ_IndexCos++;
    }

  delete [] anStep0;
  delete [] anStep1;
  delete [] anStep2;
  delete [] an2;
  delete [] an21;
  delete [] an210;

  return (0);
}

int
Cindex::nGetPeaks(C3Ddata *po3Ddata, const float fStep, float *pfPeakMin)
{
  // Find 3 decent peaks in the C3Ddata object created by nIndexCos
  // DO NOT use centroid determination, but just find them.

  // First, find lower limit.  This is easy, since origin pixel
  // has a value of m_fMaxFourier+1.
  // Then neighbors (1,1,2), (1,2,1) and (2,1,1) are down
  // a notch, so use the average of these for the peak minimum

  int i, j, k;                // Loop counters
  float fValue;
  float fLength;
  float fPeakFlag;

  int nPeaks;
  int nExt0, nExt1, nExt2;
  nExt0 = po3Ddata->nGetExt(0);
  nExt1 = po3Ddata->nGetExt(1);
  nExt2 = po3Ddata->nGetExt(2);

  int nGrid = nExt0 - 2;

  float fMin = 0.0;
  if (NULL != pfPeakMin)
    {
      // Use input peak minimum if it is greater than 0

      if (0.0 < *pfPeakMin)
        {
          fMin = *pfPeakMin;
        }
    }
  if (0.0 == fMin)
    {
      //  fMin = fMin + po3Ddata->fGetValue(1,nGrid+1,nGrid+1); //Like 0, 0, 0
      fMin = fMin + po3Ddata->fGetValue(1,nGrid+2,nGrid+1);  // Like 0, 1, 0
      fMin = fMin + po3Ddata->fGetValue(2,nGrid+1,nGrid+1);  // Like 1, 0, 0
      fMin = fMin / 2.0f;
      fMin = fMin * 0.90f;
    }

  if (fMin < (m_fFourierMax / 3.0))
    {
      cout << "WARNING, fMin is too low: " << fMin  / m_fFourierMax << endl;
      fMin = m_fFourierMax / 3.0f;
    }
  else
    {
      cout << "Start peak minimum: " << fMin / m_fFourierMax << endl;
    }

  fMin = fMin + m_fFourierMax / 20.0f;

  // Place centroids at the end of the m_poDiffVecs list
  // Be careful to note in other routines that the m_poDiffVecs has 2 kinds of
  // data:  the difference vectors from 0->(m_nNumDiffVecs-1) and
  //        the direct space peaks from m_nNumDiffVecs->last reflection
  //

  fPeakFlag = m_fFourierMax * 2.0f;

  Crefln oRefln(m_poDiffVecs);

  do {
    nPeaks     = 0;
    fMin       = fMin - (m_fFourierMax / 20.0f);
    for (k = 1; k < (nExt2-1); k++)
      {
        for (j = 1; j < (nExt1-1); j++)
          {
            for (i = 1; i < (nExt0-1); i++)
              {
                // Peak criteria, must be >=fMin and be higher than all
                // its 26 neighbors
                fValue = po3Ddata->fGetValue(i, j, k);
                if ( (fValue >= fMin) && (m_fFourierMax >= fValue) )
                  {
                    // Test neighbor criteria
                    if (   (fValue >= po3Ddata->fGetValue(i+1, j,   k-1))
                        && (fValue >= po3Ddata->fGetValue(i+1, j,   k))
                        && (fValue >= po3Ddata->fGetValue(i+1, j,   k+1))
                        && (fValue >= po3Ddata->fGetValue(i+1, j-1, k-1))
                        && (fValue >= po3Ddata->fGetValue(i+1, j-1, k))
                        && (fValue >= po3Ddata->fGetValue(i+1, j-1, k+1))
                        && (fValue >= po3Ddata->fGetValue(i+1, j+1, k-1))
                        && (fValue >= po3Ddata->fGetValue(i+1, j+1, k))
                        && (fValue >= po3Ddata->fGetValue(i+1, j+1, k+1))
                        
                        && (fValue >= po3Ddata->fGetValue(i-1, j,   k-1))
                        && (fValue >= po3Ddata->fGetValue(i-1, j,   k))
                        && (fValue >= po3Ddata->fGetValue(i-1, j,   k+1))
                        && (fValue >= po3Ddata->fGetValue(i-1, j-1, k-1))
                        && (fValue >= po3Ddata->fGetValue(i-1, j-1, k))
                        && (fValue >= po3Ddata->fGetValue(i-1, j-1, k+1))
                        && (fValue >= po3Ddata->fGetValue(i-1, j+1, k-1))
                        && (fValue >= po3Ddata->fGetValue(i-1, j+1, k))
                        && (fValue >= po3Ddata->fGetValue(i-1, j+1, k+1))
                        
                        && (fValue >  po3Ddata->fGetValue(i, j,   k-1))
                        && (fValue >  po3Ddata->fGetValue(i, j,   k+1))
                        && (fValue >= po3Ddata->fGetValue(i, j-1, k-1))
                        && (fValue >= po3Ddata->fGetValue(i, j-1, k))
                        && (fValue >= po3Ddata->fGetValue(i, j-1, k+1))
                        && (fValue >= po3Ddata->fGetValue(i, j+1, k-1))
                        && (fValue >= po3Ddata->fGetValue(i, j+1, k))
                        && (fValue >= po3Ddata->fGetValue(i, j+1, k+1)))
                      {
                        // Peak found!
                        nPeaks++;

                        // It is a new peak, so look at it more closely

                        int i1, j1, k1;
                        float f0, f1, f2;
                        i1 = i-1;
                        j1 = j-nGrid-1;
                        k1 = k-nGrid-1;
                        f0 = (float) i1 * fStep;
                        f1 = (float) j1 * fStep;
                        f2 = (float) k1 * fStep;
                        if (8 <= m_nVerbose)
                          {
                            cout << "Peak " << nPeaks << " found at i,j,k: "
                                 << i1 << ", "
                                 << j1 << ", "
                                 << k1 << endl;
                            cout << "with height " << fValue << endl;
                          }

                        // Flag as found

                        po3Ddata->vSetValue(fPeakFlag + (float)1.0, i, j, k);
                        po3Ddata->vSetValue(fPeakFlag, i, j, k-1);
                        po3Ddata->vSetValue(fPeakFlag, i, j, k+1);
                        po3Ddata->vSetValue(fPeakFlag, i-1, j, k);
                        po3Ddata->vSetValue(fPeakFlag, i+1, j, k);
                        po3Ddata->vSetValue(fPeakFlag, i, j-1, k);
                        po3Ddata->vSetValue(fPeakFlag, i, j+1, k);

                        // Now compute centroid in better way
                        // (Cube should be perhaps +-fStep)
                        // Look at an +-fStep Angstrom cube with grid of
                        // (fStep / 15) Ang
                        // This should be adjusted based on parent grid
                        // And up to 250 difference vectors

                        float fMinA[3], fMaxA[3];
                        float fHalfExtent = min(3.0f, fStep);
//                        float fGridFine   = 0.20;
                        float fGridFine   = fStep / 15.0f;
                        int   nGridFine   = (int) (0.5  +  2.0 * fHalfExtent
                                                           / fGridFine);
                        fMinA[0] = f0 - fHalfExtent;
                        fMinA[1] = f1 - fHalfExtent;
                        fMinA[2] = f2 - fHalfExtent;
                        fMaxA[0] = f0 + fHalfExtent;
                        fMaxA[1] = f1 + fHalfExtent;
                        fMaxA[2] = f2 + fHalfExtent;
                        float afCOG[3];

                        // TO DO:  Maybe we should not do the fine grid
                        //         searching on every peak.  Maybe we
                        //         should do it only in nGetUBMatrix2 on
                        //         the peaks of the reduced primitive cell.

//                        afCOG[0] = f0; afCOG[1] = f1; afCOG[2] = f2;
                        
                        nIndexCos(&fMinA[0], &fMaxA[0], nGridFine, 200,
                                  0, afCOG);

                        // Insert peak into m_poDiffVecs with special flag

                        fLength = (float) sqrt((double)(  afCOG[0] * afCOG[0]
                                         + afCOG[1] * afCOG[1]
                                         + afCOG[2] * afCOG[2]));
                        if (8 <= m_nVerbose)
                          {
                            cout << "Angstroms, len: "
                                 << f0 << ", "
                                 << f1 << ", "
                                 << f2 << ", "
                                 << (float)sqrt((double)(f0*f0 + f1*f1 + f2*f2)) <<  endl;
                                fLength =  (float)sqrt((double)(  afCOG[0] * afCOG[0]
                                             + afCOG[1] * afCOG[1]
                                             + afCOG[2] * afCOG[2]));
                            cout << "COG: "
                              << afCOG[0] << ", "
                              << afCOG[1] << ", "
                              << afCOG[2] << ", "
                              << fLength
                              << '\n'
                              << endl;
                          }
// Store the peak for later use
                        oRefln.vSetIntensity(fValue / m_fFourierMax);
                        oRefln.vSetL((int)-100);  // L is set to -100
                        oRefln.vSetField(m_poDiffVecs->m_nFI_fFloatH,
                                         fHalfExtent);
                        oRefln.vSetField(m_poDiffVecs->m_nFI_fFloatK,
                                         (float) nGridFine);
                        oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord0,
                                         afCOG[0]);
                        oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord1,
                                         afCOG[1]);
                        oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord2,
                                         afCOG[2]);
                        oRefln.vSetField(m_nFI_fDiffLength, fLength);
                        oRefln.vSetField(m_poDiffVecs->m_nFI_nDiffFreq, -100);
                        m_poDiffVecs->nInsert(&oRefln);

                        // Every peak is centrosymmetric, so add the inverse
                        // peak, too

                        oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord0,
                                         -afCOG[0]);
                        oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord1,
                                         -afCOG[1]);
                        oRefln.vSetField(m_poDiffVecs->m_nFI_fRecipCoord2,
                                         -afCOG[2]);
                        m_poDiffVecs->nInsert(&oRefln);
                      }
                  }
                else if (fValue > fPeakFlag)
                  {
                    nPeaks++;  // It was a previously found peak
                  }
              }
          }
      }
  } while ( (3 > nPeaks) && ((m_fFourierMax / 5.0) < fMin));

  cout << "Used  peak minimum: " << fMin  / m_fFourierMax << endl;
  cout << "Number of peaks found: " << nPeaks << endl;

  if (NULL != pfPeakMin)
    {
      // Return the actually used peak minimum
      *pfPeakMin = fMin;
    }

// For debugging, write out modified o3Ddata as image for viewing with dtdisplay

  if (5 <= m_nVerbose)
    {
      Cimage oImage(po3Ddata, 0, NULL);
      oImage.m_oHeader.nReplaceValue(Cimage_header::ms_sSaturatedValue,
                                   m_fFourierMax * 1.5);
      oImage.nWrite(Cstring("indexcosC.img"));
    }

  if (nPeaks < 3)
    {
      cout << "ERROR in Cindex::nGetPeaks, less than 3 peaks found!\n"
           << "  Try a longer maximum cell length.\n\n";
      return (1);   // If less than 3 peaks found, then return error
    }
  return (0);
}


int
Cindex::nCalcGetMaxCell(float *pfMaxCell)
{
  // Try to calculate a reasonable default value for the maximum cell
  // length that can/should be indexed based on the minimum spot
  // separation and the crystal to detector distance.

  float fNum;
  *pfMaxCell = 1000.0;
  if (NULL != m_ppoDetector[0])
    {
      fNum = m_ppoDetector[0]->fGetDistance();
      if (0.0 != fNum)
        {
          *pfMaxCell = min(*pfMaxCell, fabs(fNum) * m_fWavelength
                           / m_fMMSeparationMin);
        }
    }
  return (0);
}


void
Cindex::vInitSolutions(void)
{
  // Initialize the m_poSolutions object

  if (NULL != m_poSolutions)
    delete m_poSolutions;

  m_poSolutions  = new Creflnlist (*m_poReflnlist);
}



void
Cindex::vSetKnownCell(const float *pfCell)
{
  int i;
  for (i = 0; i < 6; i++)
    {
      m_a6fInputCell[i] = pfCell[i];
    }
}

void
Cindex::vSetKnownErrors(const float *pfErrors)
{
  int i;
  for (i = 0; i < 6; i++)
    {
      m_a6fInputErrors[i] = pfErrors[i];
    }
}

float
Cindex::fCalcGetGrid(void)
{
  // Try to determine the largest acceptable grid interval in Angstroms
  // Method: Compute a small portion of the indexing function near the
  //         origin.  Look at 1 grid point away from the origin to
  //                  see if peak is high enough (>90% of origin) to
  //                  be a meaningful starting minimum peak in nGetPeaks.

  int  i;
  int  nStat        = 0;
  float a10fGrid[10] = {3.0f, 2.5f, 2.0f, 1.5f, 1.0f, 0.75f, 0.50f, 0.25f, 0.20f, 0.15f };
  float fAvg        = 0.0;
  for (i = 0; (i < 10) && (fAvg < (m_fFourierMax * 0.9)); i++)
    {
      nStat = nIndexCos(0.0, 2.0f*a10fGrid[i], 2, 50, &fAvg);
      if (5 <= m_nVerbose)
        {
          cout << "Grid stepsize: " << a10fGrid[i] << "  Peak test: " << fAvg
               << endl;
        }
    }
  i--;
  return (a10fGrid[i]);
}

int
Cindex::nIndexDPS(const int nNumRefsUsed,bool bCoarse)
{
  // Index based on Steller, Bolotovsky, Rossmann algorithm.  See
  // Steller, Bolotovsky, and Rossmann (1997) J. Appl. Cryst. 30,1036-1040.

/*
  1.  For each of ~7300 directions, t(j), do
      Compute projections p,
                  p = x(i) . t(j)) for each reciprocal lattice vector x(i)
                  -      -      (or difference vector?)
  2.  Place p in a histogram array T based on |t|, so for a given distance
      from origin you will have 1, 2, or n vectors than are within that
      distance.

  3.  Invert the histogram array to discover any periodicity in it.  Save
      the directions d(j) with the largest Fourier coefficients.  Probably
      need to re-do the 30 largest answers on a finer grid.

  4.  Take the 3 non-coplanar directions and lengths to get a basis set for
      cell reduction.

  5.  Reduce cell as usual.
*/
  Cimage *poImgHist = NULL;
  Cimage *poImgFFT  = NULL;
  //Cgoniometer  *poGon;

  float a3fVecPhi[3] = {1.0, 0.0, 0.0};
  float a3x3fMatPhi[3][3];
  float a3fVecPsi[3] = {0.0, 1.0, 0.0};
  float a3x3fMatPsi[3][3];
  float a3x3fMatPsiOffset[3][3];
  float a3x3fMatRot[3][3];
  float a3fVecT[3];
  float a3fVecT2[3];

  int nMaxBins = 1024;
  if (5 <= m_nVerbose)
    {
      poImgHist = new Cimage(5 + nMaxBins, 17000);
      poImgFFT  = new Cimage(5 + nMaxBins/2, 17000);
    }

  int    i, j, m;
  float  fPsi, fPhi;  // Polar angles in degrees
  float fPsiOffset;

//  for (0 <= psi <= 90); for (0 < phi <= 360)

  float fIncPhi, fIncPsi;
  float fNominalPsiPhiStep;
  float fNumSteps;
  int   nNumSteps;
  int   nNumVecs = 0;

  float fFFTMax = 0;
  int   nFFTMax = 0;
  int   nFFTRow = 0;

  int    nRef;
  float  fCellMin      = 2.0;
  float  fCellMax      = m_fCellLengthMax;

  Crefln     *poRefln;
  Crefln     *poRefln2;
  Creflnlist *poAnswers;
  float      *pfData    = NULL;  // Scratch array for histogramming and FFT
  int         nDataSize = 0;

  poAnswers  = new Creflnlist(40000);
  m_nFI_fPhi   = poAnswers->nExpandGetField("fPhi");
  m_nFI_fPsi   = poAnswers->nExpandGetField("fPsi");
  m_nFI_fTvec0 = poAnswers->nExpandGetField("fTvec0");
  m_nFI_fTvec1 = poAnswers->nExpandGetField("fTvec1");
  m_nFI_fTvec2 = poAnswers->nExpandGetField("fTvec2");
  m_nFI_nBins  = poAnswers->nExpandGetField("nBins");
  m_nFI_fDelta = poAnswers->nExpandGetField("fDeltaBin");

  Crefln oReflnAnswer(poAnswers);

  // Prepare difference vectors
  // Save number of difference vectors in list, since we are going to add
  // to the list, but wish to distinguish between the old and newly added info

  m_nNumDiffVecs = m_poDiffVecs->nGetNumReflns();

  // Sort m_poDiffVecs on increasing difference vector length

  if (NULL != m_pnIndex)
    {
      delete [] m_pnIndex;
    }

  int nDiffs = nNumRefsUsed;
  if (0 >= nDiffs)
    nDiffs = 200;

  int nVecs = min(nDiffs, m_nNumDiffVecs);

  m_pnIndex = new int [m_nNumDiffVecs];

  if (m_bFirstTime_IndexDPS)
    {
//      s_bFirstTime = FALSE;
      for (m = 0; m < m_nNumDiffVecs; m++)
        {
          poRefln = m_poDiffVecs->poGetRefln(m);
          poRefln->vSetL(0);
        }

      // TODO: Choose the n shortest, the n most frequent, the n with smallest
      //       x, y, z values.

      m_poDiffVecs->vSort(eReflnField_float_type, m_nFI_fDiffLength, m_pnIndex);

      // Decrement the L value of the nVecs shortest ones

      for (m = 0; m < nVecs; m++)
        {
          j = m_pnIndex[m];
          poRefln = m_poDiffVecs->poGetRefln(j);
          poRefln->vSetL(poRefln->nGetL() - 2);
        }

      // Sort on frequency

      m_poDiffVecs->vSort(eReflnField_int_type,
                          m_poDiffVecs->m_nFI_nDiffFreq, m_pnIndex);

      // Decrement the L value of the nVecs most frequent ones

      for (m = m_nNumDiffVecs-1;
           m >= max(0, m_nNumDiffVecs - nVecs); m--)
        {
          j = m_pnIndex[m];
          poRefln = m_poDiffVecs->poGetRefln(j);
          poRefln->vSetL(poRefln->nGetL() - 1);
        }
    }

  // For speed, choose the nDiffs with smallest L values, so sort on nL

  m_poDiffVecs->vSort(eReflnField_int_type, m_poDiffVecs->m_nFI_nL,
                      m_pnIndex);

  nRef = m_nNumDiffVecs;

  nDiffs = min(nDiffs, nRef);
  
  if (!bCoarse) {
      cout << "Number of reflections/vectors used: " << nDiffs << endl;
      if (3 <= m_nVerbose)
          cout << "Number of reciprocal lattice vectors used: " << nDiffs << endl;
  };

  if (bCoarse)
      fNominalPsiPhiStep = 0.3f;
  else
      fNominalPsiPhiStep = 0.03f;

  fIncPsi    = fNominalPsiPhiStep / Gs_dRADIANS_PER_DEGREE;
  fPsiOffset = 0.0;
  vConvRotVec3DMat3D(fPsiOffset,
                     (float*)a3fVecPsi, &a3x3fMatPsiOffset[0][0]);
  fPsi       = -fIncPsi;
  int nOddEven = 1;
  int nStep;
  m_fSumFFTResidual = 0.0;
  while (90.0f > fPsi)
    {
      fPsi += fIncPsi;
      vConvRotVec3DMat3D(fPsi,
                         (float*)a3fVecPsi, &a3x3fMatPsi[0][0]);

      // Radius of the circle sin(fPsi), so perimeter is 2pi*r
      fNumSteps = 1.0 + (2.0 * Gs_dPI * sin(fPsi * Gs_dRADIANS_PER_DEGREE)
                         / (double) fNominalPsiPhiStep);
      nNumSteps = max(1, nint(fNumSteps));
      fIncPhi   = 360.0f / (float) nNumSteps;
      nOddEven  = 3 - nOddEven;  // 2, 1, 2, 1
      fPhi      = (fIncPhi * 0.5f) * (float)(nOddEven % 2);
      for (nStep = 0; nStep < nNumSteps; nStep++)
        {
          fPhi += fIncPhi;
          vConvRotVec3DMat3D(fPhi,
                             (float*)a3fVecPhi, &a3x3fMatPhi[0][0]);
          vMulMat3DMat3D(a3x3fMatPhi, a3x3fMatPsi, a3x3fMatRot);
          vMulMat3DMat3D(a3x3fMatPsiOffset, a3x3fMatRot, a3x3fMatPhi);
          vMulMat3DVec3D(a3x3fMatPhi, (float*)a3fVecPhi, (float*)a3fVecT);
//          vListVec3D(a3fVecT);

//          printf("ATOM %6d  CA  FRA %5d    %8.3f%8.3f%8.3f  0.00  0.00   6\n",
//                 nNumVecs, nNumVecs,
//                 25.0 * a3fVecT[0], 25.0 * a3fVecT[1], 25.0 * a3fVecT[2]);

          // OK, now have vector t

          // For all diffvecs, compute t . diffvec, and |t.diffvec|
          // save min, max length, place in histogram for FFT

          if ((0 == (nNumVecs % 200)) && (!bCoarse))
            {
              if (3 <= m_nVerbose)
                cout << "...working on dir " << nNumVecs << endl << flush;
              else
                cout << '.' << flush;
            }

          nDPSFFT(a3fVecT, nNumVecs, nDiffs, &nDataSize, &pfData,
                  &oReflnAnswer, poImgHist, poImgFFT);

          oReflnAnswer.vSetField(m_nFI_fPhi,   fPhi);
          oReflnAnswer.vSetField(m_nFI_fPsi,   fPsi);
          oReflnAnswer.vSetField(m_nFI_fTvec0, a3fVecT[0]);
          oReflnAnswer.vSetField(m_nFI_fTvec1, a3fVecT[1]);
          oReflnAnswer.vSetField(m_nFI_fTvec2, a3fVecT[2]);
          poAnswers->nInsert(&oReflnAnswer);
          nNumVecs++;

        } // end fPhi loop
    } // end fPsi loop

  cout << '\n' << flush;

  // Now sort answers on intensity (peak height of FFT)

  poAnswers->vSort(eReflnField_float_type, poAnswers->m_nFI_fIntensity,
                   NULL);

  m_fMinFFTResidual = poAnswers->poGetRefln(poAnswers->pnGetSortIndex()[0])->fGetIntensity();
  m_fSumFFTResidual = 0.0;
  for (i = 0; i < poAnswers->nGetNumReflns();i++) {
      m_fSumFFTResidual += poAnswers->poGetRefln(poAnswers->pnGetSortIndex()[i])->fGetIntensity();
  };
  if (bCoarse)
      return 0;
  
  // poAnswers->nWrite("dps.ref");
  // Loop through sorted answers and refine top nAnswers ones

  int nMaxAnswers = 30; // a number less than 30 failed in at least 1 case
  int nAnswers;

  i = 0;
  int *pnIndex;

  cout << "...refining best " << nMaxAnswers << " directions and lengths...\n"
       << flush;
  pnIndex  = poAnswers->pnGetSortIndex();
  nAnswers = nMaxAnswers;    // Count this down
  while ( (0 < nAnswers) && (i < poAnswers->nGetNumReflns()) )
    {
      // Get next T vector and try to refine it

      poRefln = poAnswers->poGetRefln(pnIndex[i]);
      poRefln->vSetK(0); // To be sure!
      a3fVecT[0] = poRefln->fGetField(m_nFI_fTvec0);
      a3fVecT[1] = poRefln->fGetField(m_nFI_fTvec1);
      a3fVecT[2] = poRefln->fGetField(m_nFI_fTvec2);

      // Is this Tvector within 5 degrees of a previously refined vector?

      bool bNoMatchFound = TRUE;
      for (j = 0; (j < i) && bNoMatchFound; j++)
        {
          poRefln2 = poAnswers->poGetRefln(pnIndex[j]);
          if (0 == poRefln2->nGetK())
            {
              // This previous refln not marked for redundancy, so check it

              a3fVecT2[0] = poRefln2->fGetField(m_nFI_fTvec0);
              a3fVecT2[1] = poRefln2->fGetField(m_nFI_fTvec1);
              a3fVecT2[2] = poRefln2->fGetField(m_nFI_fTvec2);

              // Compute the cosine of the angle between i and j with dot product
              // Note: cos(7.25deg) = 0.992

              float fCos = fabs(fDot3D(a3fVecT, a3fVecT2));
              if (0.992 < fCos)
                {
                  // The vector is within 5 degrees of previously refined vector,
                  // so mark for deletion

                  bNoMatchFound = FALSE;
                  poRefln->vSetH(j+1);
                  poRefln->vSetK(2);
                }
            }
        }
      if (bNoMatchFound)
        {


          // This reflection not within X degrees of previous result,
          // so refine it by testing several nearby directions to get
          // maximum FFT coefficient and best direction

          // fIncPhi is a function of fPsi, so recompute it

          fPsi = poRefln->fGetField(m_nFI_fPsi);
          fNumSteps = 1.0 + (2.0 * Gs_dPI
                              * sin((double)fPsi * Gs_dRADIANS_PER_DEGREE)
                              / (double)fNominalPsiPhiStep);
          nNumSteps = max(1, nint(fNumSteps));
          fIncPhi   = 360.0f / (float) nNumSteps;
          poRefln->vSetH(i);  // For debugging
          

          nDPSVecRefine(fIncPhi, fIncPsi,
              nNumVecs, nDiffs, &nDataSize, &pfData,
              poRefln, poImgHist, poImgFFT);

        

          // Check if this refined close to something previously found
        
          a3fVecT[0] = poRefln->fGetField(m_nFI_fTvec0);
          a3fVecT[1] = poRefln->fGetField(m_nFI_fTvec1);
          a3fVecT[2] = poRefln->fGetField(m_nFI_fTvec2);
          for (j = 0; (j < i) && bNoMatchFound; j++)
            {
              poRefln2 = poAnswers->poGetRefln(pnIndex[j]);
              if (0 == poRefln2->nGetK())
                {
                  // This previous refln not marked as redundant, so check it

                  a3fVecT2[0] = poRefln2->fGetField(m_nFI_fTvec0);
                  a3fVecT2[1] = poRefln2->fGetField(m_nFI_fTvec1);
                  a3fVecT2[2] = poRefln2->fGetField(m_nFI_fTvec2);

                  // Compute the cosine of the angle between i and j with dot product
                  // Note: cos(7.25deg) = 0.992

                  float fCos = fabs(fDot3D(a3fVecT, a3fVecT2));
                  if (0.992 < fCos)
                    {
                      // The vector is within 7.25 degrees of previously refined vector,
                      // so mark for deletion

//                      cout << "Match after!\n";
                      bNoMatchFound = FALSE;
                      if (poRefln->fGetIntensity() < poRefln2->fGetIntensity())
                        {
                          // This newly refined direction has better FFT coefficient,
                          // so replace earlier one
                        
                          poRefln2->vSetH(j+1);
                          poRefln2->vSetK(2);
                        }
                      else
                        {
                          poRefln->vSetH(j+1);
                          poRefln->vSetK(2);
                        }
                    }
                }
            }
          if (bNoMatchFound)
            {
              cout << '.' << flush;
              nAnswers--;              
            }
        }

      // Next reflection
      i++;
    }

  cout << "  done.\n" << flush;

  // Flag answers that were not refined

  j = poAnswers->nGetNumReflns();
  for (; i < j; i++)
    {
      poRefln = poAnswers->poGetRefln(pnIndex[i]);
      poRefln->vSetK(2);
    }

  if (3 <= m_nVerbose)
    poAnswers->nWrite("dtindexdps1.ref");

  // Delete answers that were flagged as duplicate not refined

  poAnswers->nSelect("-nK==2");
  poAnswers->nDelete((const char*)NULL, (const bool)FALSE);

  if (3 <= m_nVerbose)
    {
      poAnswers->vSort(eReflnField_float_type, poAnswers->m_nFI_fIntensity,
                       NULL);
      poAnswers->nWrite("dtindexdps2.ref");
    }

  // Copy refined answers to m_poDiffVecs for use in nGetUBMatrix2()

  j = poAnswers->nGetNumReflns();
  Crefln oReflnDF(m_poDiffVecs);

  int k;
  float fSum;
  float fDot;

  cout << "Number of vectors used for integer residual calculation: " << nDiffs << endl;
  for (i = 0; i < j; i++)
    {
      poRefln    = poAnswers->poGetRefln(i);
      
      a3fVecT[0] = poRefln->fGetField(m_nFI_fTvec0) * poRefln->fGetSigmaI();
      a3fVecT[1] = poRefln->fGetField(m_nFI_fTvec1) * poRefln->fGetSigmaI();
      a3fVecT[2] = poRefln->fGetField(m_nFI_fTvec2) * poRefln->fGetSigmaI();

      // Compute integer residual

      fSum = 0.0;

      for (k = 0; k < nDiffs; k++)
        {
          poRefln2 = m_poDiffVecs->poGetRefln(m_pnIndex[k]);
          a3fVecT2[0] = poRefln2->fGetField(m_poDiffVecs->m_nFI_fRecipCoord0);
          a3fVecT2[1] = poRefln2->fGetField(m_poDiffVecs->m_nFI_fRecipCoord1);
          a3fVecT2[2] = poRefln2->fGetField(m_poDiffVecs->m_nFI_fRecipCoord2);
          fDot        = fDot3D(a3fVecT, a3fVecT2) / m_fWavelength;
//          if (0 == i)
//            cout << "k, fDot: " << k << ", " << fDot << endl;
          fSum       += fabs(fDot - nint(fDot));
        }
      fSum = 1.0f - (fSum / (float)nDiffs);

      if (3 <= m_nVerbose)
        cout << "Soln, resid, length: " << i << ", " << fSum << ", "
             << poRefln->fGetSigmaI() << endl;

      poRefln->vSetField(m_nFI_fDelta, fSum);
      oReflnDF.vSetIntensity(fSum);

      oReflnDF.vSetL((int)-100);  // L is set to -100
      oReflnDF.vSetField(m_poDiffVecs->m_nFI_fRecipCoord0, a3fVecT[0]);
      oReflnDF.vSetField(m_poDiffVecs->m_nFI_fRecipCoord1, a3fVecT[1]);
      oReflnDF.vSetField(m_poDiffVecs->m_nFI_fRecipCoord2, a3fVecT[2]);
      oReflnDF.vSetField(m_nFI_fDiffLength, poRefln->fGetSigmaI());
      oReflnDF.vSetField(m_poDiffVecs->m_nFI_nDiffFreq, -100);
      m_poDiffVecs->nInsert(&oReflnDF);

      // Every peak is centrosymmetric, so add the inverse peak, too

      // tjn: (8/9/01)  
      // NO.  The cell reduction is good enough to avoid this.
      // oReflnDF.vSetField(m_poDiffVecs->m_nFI_fRecipCoord0, -a3fVecT[0]);
      // oReflnDF.vSetField(m_poDiffVecs->m_nFI_fRecipCoord1, -a3fVecT[1]);
      // oReflnDF.vSetField(m_poDiffVecs->m_nFI_fRecipCoord2, -a3fVecT[2]);
      // m_poDiffVecs->nInsert(&oReflnDF);
    }
  if (3 <= m_nVerbose)
    {
      poAnswers->vSort(eReflnField_float_type, m_nFI_fDelta, NULL);
      poAnswers->nWrite("dtindexdps3.ref");
      m_poDiffVecs->nWrite("dpsdiffs.ref", NULL, NULL, NULL, m_nNumDiffVecs,
                           m_nNumRowVecs-1);
    }

  delete poAnswers;

  if (0 <= m_nVerbose)
    {
      if (NULL != poImgHist)
        {
          poImgHist->m_oHeader.nReplaceValue(D_K_DtdisplayOrientation, "+X+Y");
          poImgHist->nWrite("hist.img");
        }
      if (NULL != poImgFFT)
        {
          poImgFFT->m_oHeader.nReplaceValue(D_K_DtdisplayOrientation, "+X+Y");
          poImgFFT->nWrite("FFT.img");
        }
    }
  if (NULL != poImgHist)
    delete poImgHist;
  if (NULL != poImgFFT)
    delete poImgFFT;
  if (NULL != pfData)
    delete [] pfData;

  return 0;
}

int
Cindex::nDPSVecRefine(const float fPhiInc, const float fPsiInc,
                      const int nRow, const int nNumRefsUsed, int *pnDataSize, float **ppfData,
                      Crefln *poReflnAnswer, Cimage *poImgHist, Cimage *poImgFFT)
{
  // Refine the direction a3fVecT by calculating the histogram and its FFT for
  // nearby Phi, Psi values.  Return the <nearby> direction with the largest
  // Fourier coefficient.
  // Input: fPhiInc, fPsiInc, the previously used Phi,Psi search grid increments
  // Phi, Psi are in degrees!
  // Recursively call this routine with smaller and smaller Phi,Psi increments

  int     i, j;
  int     nRow2;
  Crefln *poRefln;
  float   a3fVecT[3];
  float   a3fVecNew[3];
  float   fPhiDel, fPsiDel;
  float   fPhiOrig, fPsiOrig;

  float a3fVecPhi[3] = {1.0, 0.0, 0.0};
  float a3x3fMatPhi[3][3];
  float a3fVecPsi[3] = {0.0, 1.0, 0.0};
  float a3x3fMatPsi[3][3];
  float a3x3fMatRot[3][3];

  float   fPhiNew, fPsiNew;
  float   fPhiIncNew, fPsiIncNew;

  fPhiIncNew = fPhiInc / 2.2f;
  fPsiIncNew = fPsiInc / 2.2f;

  fPhiOrig   = poReflnAnswer->fGetField(m_nFI_fPhi);
  fPsiOrig   = poReflnAnswer->fGetField(m_nFI_fPsi);
  a3fVecT[0] = poReflnAnswer->fGetField(m_nFI_fTvec0);
  a3fVecT[1] = poReflnAnswer->fGetField(m_nFI_fTvec1);
  a3fVecT[2] = poReflnAnswer->fGetField(m_nFI_fTvec2);

  // Copy the input refln for future comparison to.

  poRefln = new Crefln(*poReflnAnswer);
  m_nDepth_DPSVecRefine++;

  // Need to create loop of ever smaller Phi,Psi increments

  nRow2 = nRow;
  if (NULL != poImgHist)
    nRow2 = min(nRow2, poImgHist->nGetDimension(1)-1);

  // Need to create loop of ever smaller Phi,Psi increments
  
  int           nIbest = 0;
  int           nJbest = 0;

  for (j = -2; j <= 2; j++)
    {
      // Compute psi rotation matrix

      fPsiDel = (float) j * fPsiIncNew;
      fPsiNew = fPsiOrig + fPsiDel;
      vConvRotVec3DMat3D(fPsiNew, (float*)a3fVecPsi, &a3x3fMatPsi[0][0]);
      for (i = -2; i <= 2; i++)
        {
          if ( (0 != i) || (0 != j) )
            {
              fPhiDel = (float) i * fPhiIncNew;
              fPhiNew = fPhiOrig + fPhiDel;

              // WARNING, there is no PsiOffset here, so in ::nIndexDPS it had better be 0!

              vConvRotVec3DMat3D(fPhiNew,(float*)a3fVecPhi, &a3x3fMatPhi[0][0]);
              vMulMat3DMat3D(a3x3fMatPhi, a3x3fMatPsi, a3x3fMatRot);
              vMulMat3DVec3D(a3x3fMatRot, (float*)a3fVecPhi, (float*)a3fVecNew);

//              cout << "test: phi, psi: " << fPhiNew << ", "  << fPsiNew << endl;

              nDPSFFT(a3fVecNew, nRow2, nNumRefsUsed, pnDataSize, ppfData,
                      poRefln, poImgHist, poImgFFT);

              // Compare with current best coefficient and save the better one
              // Perhaps put a convergence test here, too.

              if (poRefln->fGetIntensity() < poReflnAnswer->fGetIntensity())
                {
                  // Retain better fit in our answer

                  *poReflnAnswer = *poRefln;
                  poReflnAnswer->vSetField(m_nFI_fPhi,   fPhiOrig + fPhiDel);
                  poReflnAnswer->vSetField(m_nFI_fPsi,   fPsiOrig + fPsiDel);
                  poReflnAnswer->vSetField(m_nFI_fTvec0, a3fVecNew[0]);
                  poReflnAnswer->vSetField(m_nFI_fTvec1, a3fVecNew[1]);
                  poReflnAnswer->vSetField(m_nFI_fTvec2, a3fVecNew[2]);
                  nIbest = i;
                  nJbest = j;
                }
            }
        }
    }

  // Recursively refine until the angular steps are less than 0.1 degree
  // TODO: perhaps a convergence test here too.

  if ( (8 > m_nDepth_DPSVecRefine) && ( (0.02f < fPhiIncNew) || (0.02f < fPsiIncNew) ) )
    {
      nRow2++;
      if ( (2 == abs(nIbest)) || (2 == abs(nJbest) ) )
        {
          // We are on the edge of the refinement grid, so reset increments
          // back to original

          fPhiIncNew = fPhiIncNew * 2.2f;
        }
//      cout << "vecrefine depth: " << m_nDepth_DPSVecRefine << endl;
      nDPSVecRefine(fPhiIncNew, fPsiIncNew,
                    nRow2, nNumRefsUsed, pnDataSize, ppfData,
                    poReflnAnswer, poImgHist, poImgFFT);
    }

/*
  cout << "nrow, orig, delta, results: "
       << nRow << ", "
       << fPhiOrig << ", "
       << fPsiOrig << ", "
       << fPhiDel << ", "
       << fPsiDel << ", "
       << fPhiOrig + fPhiDel << ", "
       << fPsiOrig + fPsiDel << endl;
*/
  delete poRefln;
  m_nDepth_DPSVecRefine--;
  return (0);
}


int
Cindex::nDPSFFT(const float a3fVecT[3], const int nRow,
                const int nNumRefsUsed, int *pnDataSize, float **ppfData,
                Crefln *poReflnAnswer, Cimage *poImgHist, Cimage *poImgFFT)
{
  int     i, j;
  float   fP, fPmin, fPmax;
  float   a3fVecD[3];
  float   a3fVecDD[3];
  Crefln *poRefln;
  float  *pfTemp;
  float  *pfData;

  int    nNumBins;
  int    nMinBins      =   64; // Must be a power of 2
  int    nMaxBins      = 2048; // Must be a power of 2

  float  fDeltaBin;
  int    nBin;
  unsigned short int uiFreq;

  int    nNeeded;
  float  fTemp;
  float  f0;

  // Loop through the input reflnlist list and form the dot product
  // p = d* . t ,  then put p in a histogram or frequency array
  //     --   -

  // Initialize fPmin and fPmax with first dot product

  j          = m_pnIndex[0];
  poRefln    = m_poDiffVecs->poGetRefln(j);
  a3fVecD[0] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord0);
  a3fVecD[1] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord1);
  a3fVecD[2] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord2);
  a3fVecDD[0] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoordD0);
  a3fVecDD[1] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoordD1);
  a3fVecDD[2] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoordD2);
  fNormVec3D(a3fVecDD);

  f0 = (1.0 - fabs(fDot3D(a3fVecDD,a3fVecT)/fLenVec3D(a3fVecT)));
  fP = fDot3D(a3fVecT, a3fVecD);
  poRefln->vSetField(m_poDiffVecs->m_nFI_fFloatH, fP);
  poRefln->vSetField(m_poDiffVecs->m_nFI_fFloatK, f0);

  fPmin = fP;
  fPmax = fP;

  for (i = 1; i < nNumRefsUsed; i++)
    {
      j          = m_pnIndex[i];
      poRefln    = m_poDiffVecs->poGetRefln(j);
      a3fVecD[0] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord0);
      a3fVecD[1] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord1);
      a3fVecD[2] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoord2);
      a3fVecDD[0] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoordD0);
      a3fVecDD[1] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoordD1);
      a3fVecDD[2] = poRefln->fGetField(m_poDiffVecs->m_nFI_fRecipCoordD2);
      fNormVec3D(a3fVecDD);

      f0 = (1.0 - fabs(fDot3D(a3fVecDD,a3fVecT)/fLenVec3D(a3fVecT)));
      fP = fDot3D(a3fVecT, a3fVecD);

      // Store P in the fFloatH field for now

      poRefln->vSetField(m_poDiffVecs->m_nFI_fFloatH, fP);
      poRefln->vSetField(m_poDiffVecs->m_nFI_fFloatK, f0);
      if (fPmin > fP) fPmin = fP;
      if (fPmax < fP) fPmax = fP;
    }

  // Now build frequency distribution (histogram)

  float fStretch;
  int nBinsBefore;

  nNumBins  = nint((fPmax - fPmin) * 5.0 * m_fCellLengthMax);
  if (NULL != poImgHist)
    nNumBins  = min(nNumBins, poImgHist->nGetDimension(0));
  nBinsBefore = nNumBins;

  // Make nNumBins a power of 2
        
  nNumBins = nint(log((double)nNumBins) / log(2.0));
  nNumBins = (int)pow(2.0, (double)nNumBins);
  fStretch = (float) nNumBins / (float) nBinsBefore;

  // cout << "nNumBins " << nNumBins << endl;

  if (nMinBins > nNumBins)
    {
      nNumBins = nMinBins;
      if (6 <= m_nVerbose)
        cout << "WARNING, minimum frequency bins is " << nMinBins << endl;
    }
  if (nMaxBins < nNumBins)
    {
      nNumBins = nMaxBins;
      if (6 <= m_nVerbose)
        cout << "WARNING, maximum frequency bins is " << nMaxBins << endl;
    }

  fDeltaBin = (fPmax - fPmin) / (float)(nNumBins-1);

  if (0.0 >= fDeltaBin)
    {
      cout << "WARNING, minimum delta bin is 0.1" << endl;
      fDeltaBin = 0.1f;
    }
//  cout << "deltabin, nNumBins: " << fDeltaBin << ", " << nNumBins << endl;

  int nBinTestMin, nBinTestMax;
  // Minimum test bin of 1/10 the cell max?  I don't like it!
//  nBinTestMin = nint(0.1f * m_fCellLengthMax * (float)(nNumBins-1) * fDeltaBin);
  nBinTestMin = 3;
  nBinTestMin = max(1, nBinTestMin);
  nBinTestMax = nint(m_fCellLengthMax * (float)(nNumBins-1) * fDeltaBin);
  nBinTestMax = min(nNumBins, nBinTestMax);

  // Make sure enough space exists for the histogram and FFT arrays
        
  nNeeded = nNumBins * 2 + 100;

  if (nNeeded > *pnDataSize)
    {
      if (NULL != *ppfData)
        delete [] *ppfData;
      *pnDataSize = nNeeded;
      *ppfData    = new float [nNeeded];
    }

  // Zero the histogram
        
  pfTemp = *ppfData;
  for (i = 0; i < nNeeded; i++, pfTemp++)
    {
      *pfTemp = 0.0f;
    }

  // Build histogram

  pfData = *ppfData;
  for (i = 0; i < nNumRefsUsed; i++)
    {
      j       = m_pnIndex[i];
      poRefln = m_poDiffVecs->poGetRefln(j);
      fP      = poRefln->fGetField(m_poDiffVecs->m_nFI_fFloatH);
      nBin    = nint((fP - fPmin) / fDeltaBin);
      if ( (0 <= nBin) && (nNumBins > nBin) )
        {
          // Note we use 2*nBin because of FFT below expects complex numbers
          pfData[2*nBin] += poRefln->fGetField(m_poDiffVecs->m_nFI_fFloatK);
        }
    }

  // As a diagnostic tool, copy histogram to poImgHist if it is available

  if (NULL != poImgHist)
    {
      m_nAtom_DPSFFT++;
      if (m_nAtom_DPSFFT < poImgHist->nGetDimension(1))
        {
          poImgHist->nSetNextPixel(0, m_nAtom_DPSFFT);
          pfTemp = pfData;
          for (i = 0; i < nNumBins; i++, pfTemp++)
            {
              uiFreq = (unsigned short int) *pfTemp;  // Use the REAL part
              poImgHist->vSetNextPixel(uiFreq);
              pfTemp++;                               // Skip the imaginary part
            }
          for (i = nNumBins; i < poImgHist->nGetDimension(0); i++)
            poImgHist->vSetNextPixel((unsigned short int) 0);
        }
    }

  // Frequency histogram done, so FFT it
  // The FFT routine expects complex number pairs, so set the imaginary component to 0.
  // This FFT routine is from Numerical Recipes and expects pfData[1] to be the first
  // real value, so pass pfData-1 to achieve this.

  nFFT1D(nNumBins, 1, pfData);

  // Everytime we do an FFT, write the vector out for viewing with O

  if (5 <= m_nVerbose)
    {
      if (NULL == poImgHist) m_nAtom_DPSFFT++;
      printf("ATOM %6d  CA  FFT %5d    %8.3f%8.3f%8.3f  0.00  0.00   6\n",
             m_nAtom_DPSFFT, m_nAtom_DPSFFT,
             10.0 * a3fVecT[0], 10.0 * a3fVecT[1], 10.0 * a3fVecT[2]);

      cout << flush;
    }

  // Convert real and imaginary parts to a single square of magnitude,
  // save largest values. pfData[1] should be 0.0

  pfTemp    = pfData;
  for (i = 0; i < nNumBins/2; i++)
    {
      fTemp       = *pfTemp * *pfTemp++;
      fTemp      += *pfTemp * *pfTemp++;
      pfData[i]   = fTemp;
    }

  // As a diagnostic tool, copy the FFT to poImgFFT image if it is available
        
  if (NULL != poImgFFT)
    {
      // Transfer FFT results to the image

      if (m_nAtom_DPSFFT < poImgFFT->nGetDimension(1))
        {
          poImgFFT->nSetNextPixel(0, m_nAtom_DPSFFT);

          // Try to normalize from 0 to 100

          pfTemp    = pfData;
          for (i = 0; i < nNumBins/2; i++)
            {
              fTemp = (float)sqrt((double)fTemp) /  (float)nNumRefsUsed * 100.0f;
              if (65535 < fTemp) fTemp = 65535;
              uiFreq = (unsigned short int) fTemp;
              poImgFFT->vSetNextPixel(uiFreq);
            }
        }
    }

  // Find a minimum away from origin

  float fMax;

  fMax = pfData[0];
  i    = nBinTestMin;
  while ( (pfData[i] < fMax) && (i < nNumBins/2) )
    {
      fMax = pfData[i];
      i++;
    }

  // Now find maximum in rest of Fourier data

  int nIMax = i;
  while ( (i < nBinTestMax) && (i < nNumBins/2) )
    {
      if (pfData[i] > fMax)
        {
          fMax  = pfData[i];
          nIMax = i;
        }
      i++;
    }

  // This gives approximate percentage of reflns indexed

  fMax = (float)sqrt((double)fMax) /  (float)nNumRefsUsed * 100.0f;

  poReflnAnswer->vSetIntensity(-fMax); // Use negative so sort is descending
  poReflnAnswer->vSetSigmaI((float)nIMax / (float)(nNumBins-1) / fDeltaBin * m_fWavelength);

  poReflnAnswer->vSetH(m_nAtom_DPSFFT);
  poReflnAnswer->vSetL(nIMax);
  poReflnAnswer->vSetField(m_nFI_nBins, nNumBins);
  poReflnAnswer->vSetField(m_nFI_fDelta, fDeltaBin);

//  cout << "Max for row " << nAtom << " is/at " << fMax << " / " << nIMax << endl;
  return (0);
}

int
Cindex::nLoadSolutionArray(void)
{
    // Loop through all direct space peaks found in the m_poDiffVecs objects
    // Try all possible combinations of 3 vectors to yield a solution
    // Sort the solutions on ...
    
    int i, h, k, l;
    int nh, nk, nl;
    int nRef;
    float fCellOut[6];
    float a6fCellReduce[6];
    float fVec1[3], fVec2[3], fVec3[3];
    float fVecCross12[3];
    //  float fH[3][3], fInvH[3][3];
    
    float fRes;
    float fWave2;
    float fCos12, fCos23, fCos13, fCosCross;
    float fAngleDel;
    float a3x3fBMatrix[3][3];
    float a3x3fInvBMat[3][3];
    float a3x3fUMatrix[3][3];
    
    float  fVolume, fVolMin, fRect, fRectMin, fRectMul;
    float fRot1, fRot2, fRot3;
    
    bool bEnterKLoop, bEnterLLoop, bCellOK;
    int  nToCheck, nLast;
    
    // Init m_poSolution object
    
    vInitSolutions();
    
    Crefln oRefln(m_poSolutions);
    
    float fAngleMin = 14.9f;  // Some rhombohedral axes will have anglemin<60
    float fAngleMax = 165.1f;
    
    nRef        = m_poDiffVecs->nGetNumReflns();
    m_nNumRowVecs = nRef;
    nLast       = nRef - 1;
    //  fWave2      = m_fWavelength * m_fWavelength;
    fWave2      = 1.0; // m_fWavelength * m_fWavelength;
    
    fVolMin     = -1.0;
    fRectMin    = -1.0;
    fRot1       = 0.0;
    fRot2       = 0.0;
    fRot3       = 0.0;
    float fResH, fResK, fResL;

   
    
    nToCheck =  m_nNumRowVecs - m_nNumDiffVecs;
    
    
    float fMinAllowedResidual;
    float fMinAllowedResidualStart = m_fMinAllowedResidualStart;
    float fMinAllowedResidualEnd   = m_fMinAllowedResidualEnd;
    
    fMinAllowedResidual = m_fMinAllowedResidualStart;
    while (   (fMinAllowedResidualEnd <= fMinAllowedResidual)
        && (DT_NUMSOL_MIN_INIT > m_poSolutions->nGetNumReflns()))
    {
        m_poSolutions->vDeleteAll();
        // Check best ones first
        
        //      cout << "minallowed resid: " << fMinAllowedResidual << endl;
        for (nh = 0; nh < nToCheck; nh++)
        {
            //    h = m_pnIndex[nLast-nh];  I do not think they are sorted
            h = nLast-nh;
            fCellOut[0] = // m_fWavelength *
                m_poDiffVecs->poGetRefln(h)->fGetField(m_nFI_fDiffLength);
            fResH = m_poDiffVecs->poGetRefln(h)->fGetIntensity();
            bEnterKLoop = TRUE;
            if (m_bKnownCell)
            {
                // Find vector which matches our input cell
                
                bEnterKLoop =    (fCellOut[0] >= m_a6fInputCell[0]
                    - m_a6fInputErrors[0])
                    && (fCellOut[0] <= m_a6fInputCell[0]
                    + m_a6fInputErrors[0]);
            }
            bEnterKLoop = (bEnterKLoop) && (fResH >= fMinAllowedResidual);
            if (bEnterKLoop)
            {
                fVec1[0] = m_poDiffVecs->poGetRefln(h)->fGetField(
                    m_poDiffVecs->m_nFI_fRecipCoord0);
                fVec1[1] = m_poDiffVecs->poGetRefln(h)->fGetField(
                    m_poDiffVecs->m_nFI_fRecipCoord1);
                fVec1[2] = m_poDiffVecs->poGetRefln(h)->fGetField(
                    m_poDiffVecs->m_nFI_fRecipCoord2);
                
                //      fH[0][0] = (float) m_poDiffVecs->poGetRefln(h)->nGetH();
                //      fH[1][0] = (float) m_poDiffVecs->poGetRefln(h)->nGetK();
                //      fH[2][0] = (float) m_poDiffVecs->poGetRefln(h)->nGetL();
            }
            for (nk = 0; (bEnterKLoop) && (nk < nToCheck); nk++)
            {
                //      k = m_pnIndex[nLast-nk];
                k = nLast - nk;
                fCellOut[1] = // m_fWavelength *
                    m_poDiffVecs->poGetRefln(k)->fGetField(m_nFI_fDiffLength);
                fResK = m_poDiffVecs->poGetRefln(k)->fGetIntensity();
                //          cout << "k is : " << k << endl;
                //          cout << "UB2 fCellOut[1] = " << fCellOut[1]
                //               << "  Resid is: " << fResK << endl;
                
                bEnterLLoop = TRUE;
                if (m_bKnownCell)
                {
                    // Find vector which matches our input cell
                    
                    bEnterLLoop =    (fCellOut[1] >= m_a6fInputCell[1]
                        - m_a6fInputErrors[1])
                        && (fCellOut[1] <= m_a6fInputCell[1]
                        + m_a6fInputErrors[1]);
                }
                bEnterLLoop = (bEnterLLoop) && (fResK >= fMinAllowedResidual);
                if (bEnterLLoop)
                {
                    fVec2[0] = m_poDiffVecs->poGetRefln(k)->fGetField(
                        m_poDiffVecs->m_nFI_fRecipCoord0);
                    fVec2[1] = m_poDiffVecs->poGetRefln(k)->fGetField(
                        m_poDiffVecs->m_nFI_fRecipCoord1);
                    fVec2[2] = m_poDiffVecs->poGetRefln(k)->fGetField(
                        m_poDiffVecs->m_nFI_fRecipCoord2);
                    
                    //    fH[0][1] = (float) m_poDiffVecs->poGetRefln(k)->nGetH();
                    //    fH[1][1] = (float) m_poDiffVecs->poGetRefln(k)->nGetK();
                    //    fH[2][1] = (float) m_poDiffVecs->poGetRefln(k)->nGetL();
                    fCos12 = fWave2 * fDot3D(fVec1, fVec2)
                        / (fCellOut[0] * fCellOut[1]);
                    fCos12 = min(fCos12, 1.0f);
                    fCos12 = max(fCos12, -1.0f);
                    fCellOut[5] = acos((double)fCos12) / Gs_dRADIANS_PER_DEGREE;
                    vCross3D(fVec1, fVec2, fVecCross12);
                    (void) fNormVec3D(fVecCross12);  // Normalize to length 1.
                    
                    if (m_bKnownCell)
                    {
                        fAngleDel   =  fabs(fCellOut[5] - m_a6fInputCell[5]);
                        bEnterLLoop =  (fAngleDel <= m_a6fInputErrors[5]);
                        
                    }
                    else
                    {
                        bEnterLLoop =    (fCellOut[5] >= fAngleMin)
                            && (fCellOut[5] <= fAngleMax);
                    }
                }
                
                for (nl = 0; (bEnterLLoop) && (nl < nToCheck); nl++)
                {
                    //        l = m_pnIndex[nLast-nl];
                    l = nLast - nl;
                    fCellOut[2] = // m_fWavelength *
                        m_poDiffVecs->poGetRefln(l)->fGetField(m_nFI_fDiffLength);
                    fResL = m_poDiffVecs->poGetRefln(l)->fGetIntensity();
                    //              cout << "L is : " << l << endl;
                    //              cout << "UB2 fCellOut[2] = " << fCellOut[2]
                    //                   << "  Resid is: " << fResL << endl;
                    bCellOK = TRUE;
                    if (m_bKnownCell)
                    {
                        // Find vector which matches our input cell
                        
                        bCellOK =   (fCellOut[2] >= m_a6fInputCell[2]
                            - m_a6fInputErrors[2])
                            && (fCellOut[2] <= m_a6fInputCell[2]
                            + m_a6fInputErrors[2]);
                    }
                    
                    bCellOK = (bCellOK) && (fResL >= fMinAllowedResidual);
                    
                    if (bCellOK)
                    {
                        fVec3[0] = m_poDiffVecs->poGetRefln(l)->fGetField(
                            m_poDiffVecs->m_nFI_fRecipCoord0);
                        fVec3[1] = m_poDiffVecs->poGetRefln(l)->fGetField(
                            m_poDiffVecs->m_nFI_fRecipCoord1);
                        fVec3[2] = m_poDiffVecs->poGetRefln(l)->fGetField(
                            m_poDiffVecs->m_nFI_fRecipCoord2);
                        
                        //fH[0][2] = (float) m_poDiffVecs->poGetRefln(l)->nGetH();
                        //fH[1][2] = (float) m_poDiffVecs->poGetRefln(l)->nGetK();
                        //fH[2][2] = (float) m_poDiffVecs->poGetRefln(l)->nGetL();
                        fCos13 = fWave2 * fDot3D(fVec1, fVec3)
                            / (fCellOut[0] * fCellOut[2]);
                        fCos13 = min(fCos13, 1.0f);
                        fCos13 = max(fCos13, -1.0f);
                        fCellOut[4] = acos((double)fCos13) / Gs_dRADIANS_PER_DEGREE;
                        fCos23 = fWave2 * fDot3D(fVec2, fVec3)
                            / (fCellOut[1] * fCellOut[2]);
                        fCos23 = min(fCos23, 1.0f);
                        fCos23 = max(fCos23, -1.0f);
                        fCellOut[3] = acos((double)fCos23) / Gs_dRADIANS_PER_DEGREE;
                        
                        // Compute angle between fVecCross12 and fVec3
                        // We want this angle to be as close to 0 or 180 as
                        // possible.
                        
                        fCosCross =  fWave2 * fDot3D(fVec3, fVecCross12)
                            / (fCellOut[2] * 1.0f);
                        fCosCross = min(fCosCross, 1.0f);
                        fCosCross = max(fCosCross, -1.0f);
                        fCosCross = 90.0f - (float)(acos(fabs((double)fCosCross)) / Gs_dRADIANS_PER_DEGREE);
                        
                        if (m_bKnownCell)
                        {
                            fAngleDel = fabs(fCellOut[4] - m_a6fInputCell[4]);
                            
                            bCellOK   = (fAngleDel <= m_a6fInputErrors[4]);
                            
                            fAngleDel = fabs(fCellOut[3] - m_a6fInputCell[3]);
                            
                            bCellOK   = bCellOK
                                && (fAngleDel <= m_a6fInputErrors[3]);
                        }
                        else
                        {
                            bCellOK =    (fCellOut[3] >= fAngleMin)
                                && (fCellOut[4] >= fAngleMin)
                                && (fCellOut[3] <= fAngleMax)
                                && (fCellOut[4] <= fAngleMax)
                                && (fCosCross   <= fAngleMax)
                                && (fCosCross   >= fAngleMin);
                        }
                        if (bCellOK) 
                        {
                            // When you are here, you have a solution? 
                            
                            // Is it right-handed?
                            
                            for (i = 0; i < 3; i++) {
                                m_a3x3fInvUBMat[i][0] = fVec1[i];
                                m_a3x3fInvUBMat[i][1] = fVec2[i];
                                m_a3x3fInvUBMat[i][2] = fVec3[i];
                            }
                            
                            // The determinant of the Inverse UB matrix is the
                            // unit cell volume
                            
                            fVolume = fInvMat3D(&m_a3x3fInvUBMat[0][0],&m_a3x3fUBMat[0][0]);
                            
                            fRectMul = fCellOut[0] * fCellOut[1] * fCellOut[2];
                            
                            if (fVolume > 5.001) {
                                // Ok, right-handed and big enough!
                                
                                fRes = (fResH + fResK + fResL) / 3.0f;
                                if (fRes >= fMinAllowedResidual)
                                {
                                    if (5 <= m_nVerbose)
                                    {
                                        cout << "minallowed resid: " << fMinAllowedResidual << endl;
                                        cout << "h is : " << h << endl;
                                        cout << "UB2 fCellOut[0] = " << fCellOut[0]
                                            << "  Resid is: " << fResH << endl << flush;
                                        cout << "k is : " << k << endl;
                                        cout << "UB2 fCellOut[1] = " << fCellOut[1]
                                            << "  Resid is: " << fResK << endl;
                                        cout << "L is : " << l << endl;
                                        cout << "UB2 fCellOut[2] = " << fCellOut[2]
                                            << "  Resid is: " << fResL << endl;
                                        cout << "UB2 Angles: " << fCellOut[3]
                                            << ", " << fCellOut[4] << ", "
                                            << fCellOut[5] << endl;
                                        cout << "fVolume is << " << fVolume << endl;
                                    }
                                    
                                    m_oCrystal.vSetCell((float *)fCellOut);
                                    
                                    if (!m_bKnownCell)
                                    {
                                        // Reduce the cell, so fRectMin comparisons are legit
                                        
                                        vCopyVecND(6, fCellOut, a6fCellReduce);
                                        nReduceCell(a6fCellReduce);
                                        fRect    = a6fCellReduce[0] + a6fCellReduce[1] + a6fCellReduce[2];
                                    }
                                    else
                                    {
                                        fRect    = fCellOut[0] + fCellOut[1] + fCellOut[2];
                                    }
                                    
                                    // Because cell reduction works better, we can
                                    // use 1.6 instead of 1.05 here
                                    
                                    if (   (   (fVolume <= (fVolMin  * 1.6))
                                        || (fVolMin  < 0.0) )
                                        && (   (fRect   <= (fRectMin * 1.50))
                                        || (fRectMin < 0.0) ) )
                                    {
                                        fVolMin = fVolume; // Allow min vol creep up
                                        if ( (fRect < fRectMin) || (0.0f > fRectMin) )
                                            fRectMin = fRect;
                                        
                                        // Make residual large and negative, so when sorted,
                                        // best solutions come first in the list
                                        
                                        if (.999 == fRes) 
                                            fRes = .998f;
                                        oRefln.vSetIntensity(fRes*(-1000.0f));
                                        oRefln.vSetSigmaI(fRect);
                                        oRefln.vSetField(m_nFI_fA,     fCellOut[0]);
                                        oRefln.vSetField(m_nFI_fB,     fCellOut[1]);
                                        oRefln.vSetField(m_nFI_fC,     fCellOut[2]);
                                        oRefln.vSetField(m_nFI_fAlpha, fCellOut[3]);
                                        oRefln.vSetField(m_nFI_fBeta,  fCellOut[4]);
                                        oRefln.vSetField(m_nFI_fGamma, fCellOut[5]);
                                        oRefln.vSetField(m_nFI_sLattice, "?");
                                        oRefln.vSetField(m_nFI_sLattSymm, "?");
                                        oRefln.vSetField(m_nFI_nLattNum, m_nLattNum);
                                        oRefln.vSetH(h);
                                        oRefln.vSetK(k);
                                        oRefln.vSetL(l);
                                        
                                        (void) m_oCrystal.nCalcBMatrix();
                                        m_oCrystal.vGetBMatrix((float *)a3x3fBMatrix);
                                        fVolume = fInvMat3D((float *)a3x3fBMatrix,(float *)a3x3fInvBMat);
                                        vMulMat3DMat3D(m_a3x3fUBMat, a3x3fInvBMat,a3x3fUMatrix);
                                        
                                        vDeconvMat3D3XYZ(a3x3fUMatrix,&fRot1, &fRot2, &fRot3);
                                        oRefln.vSetField(m_nFI_fRot1, fRot1);
                                        oRefln.vSetField(m_nFI_fRot2, fRot2);
                                        oRefln.vSetField(m_nFI_fRot3, fRot3);
                                        m_poSolutions->nInsert(&oRefln);
                                    } // end volume and rect test
                                } // end residual test
                            } // end volume test
                        }
                    }
                }  // end l loop
            }  // end k loop
        }  // end h loop
        
        nRef = m_poSolutions->nGetNumReflns();
        fMinAllowedResidual -= m_fMinAllowedResidualDelta;
        
    } // end of while loop

    while (m_poSolutions->nGetNumReflns()>DT_NUMSOL_MAX_INIT) 
        m_poSolutions->nDelete(DT_NUMSOL_MAX_INIT);
    nRef = m_poSolutions->nGetNumReflns();

    if (1 > nRef)
    {
        return (1);  // No solutions found!
    }
    
    // Save the residual that was actually used
    
    m_fMinAllowedResidual = fMinAllowedResidual + m_fMinAllowedResidualDelta;
    
    if (4 <= m_nVerbose)
    {
        cout << "fVolMin, fRectMin: " << fVolMin << ", " << fRectMin << endl;
        m_poSolutions->nWrite("soln.UB2");
    }
    
    if (NULL != m_pnSoln) delete [] m_pnSoln;
    m_pnSoln = new int [m_poSolutions->nGetNumReflns()];
    m_poSolutions->vSort(eReflnField_float_type, 0, m_pnSoln);
    
    return 0;
}

int
Cindex::nCalcReflnResolution(Creflnlist *poReflnlist)
{
  // Calculate the resolution of the reflns in poReflnlist
  // and place in the field m_nFI_fDetResolution
  // Use the CDetector and CSource member objects

  if (NULL == m_poSource)
    {
      cout << "ERROR in Cindex::nCalcReflnResolution.  No Source object.\n"
           << endl;
      return (-1);
    }
  if (NULL == m_ppoDetector[0])
    {
      cout << "ERROR in Cindex::nCalcReflnResolution.  No Detector object.\n"
           << endl;
      return (-1);
    }
  
  m_ppoDetector[0]->vUpdateDDDN();    // Required for resolution calc below

  int    nFI0, nFI1;
  int    nx;
  int    nRefs;
  double fPx0, fPx1;
  float  fDetReso;
  Crefln *poRefln;

  nFI0 = poReflnlist->m_nFI_fObsPx0;
  if (0 > nFI0)
    nFI0 = poReflnlist->m_nFI_fCalcPx0;
  nFI1 = poReflnlist->m_nFI_fObsPx1;
  if (0 > nFI1)
    nFI1 = poReflnlist->m_nFI_fCalcPx1;

  if ( (0 > nFI0) || (0 > nFI1) )
    {
      cout << "ERROR in Cindex::nCalcReflnResolution.  No spot centroids.\n"
           << endl;
      return (-1);
    }

  if (0 > poReflnlist->m_nFI_fDetResolution)
    poReflnlist->nExpandGetField(Creflnlist::ms_sfDetResolution);

  nRefs = poReflnlist->nGetNumReflns();
  for (nx = 0; nx < nRefs; nx++)
    {
      poRefln = poReflnlist->poGetRefln(nx);
      fPx0 = (double)poRefln->fGetField(nFI0);
      fPx1 = (double)poRefln->fGetField(nFI1);
      fDetReso = m_fWavelength 
        * m_ppoDetector[0]->fCalcGetResolution(fPx0, fPx1, m_a3fS0);
      
      poRefln->vSetField(poReflnlist->m_nFI_fDetResolution, fDetReso);
    }
  return (0);
}

double Cindex::fGetKnownCellResidual(double a6fKnownCell[6],
                                     Ccrystal& oCrystal,
                                     double fPromptingResidual)
{
    // This will not use the Lattice characters
    // Only looks for cell using vectors from Cindex.
    
    double a3x3fOrientMatOut[3][3];           
    double a3x3fRealOrientMatOut[3][3];
    double a3x3fRealOrientMatIn[3][3];
    double a3x3fOrientMatIn[3][3];
    double a3x3fRealMatOut[3][3];
    double a6fUserCellFound[6];
    float  pfOrientAngles[3];
    bool  bUseFound;
    double fResultResidual = -1.0;
    CCellReduce oReduce;
    Cstring sTemp;


    oCrystal.nCalcOrientMatrix();
    oCrystal.vGetOrientMatrix(&a3x3fOrientMatIn[0][0]);
    fInvMat3D(&a3x3fOrientMatIn[0][0],&a3x3fRealOrientMatIn[0][0]);
    vTranMat3D(a3x3fRealOrientMatIn);

    // Try to find the cell, and print out the residual if found.
    oCrystal.vSetCell(a6fKnownCell);
    
    oCrystal.nCalcGetRealMatrix(&a3x3fRealMatOut[0][0]);    
    vTranMat3D(a3x3fRealMatOut);
    
    fResultResidual = oReduce.fFindOrigVecs(a3x3fRealMatOut, a3x3fRealOrientMatIn, a3x3fRealOrientMatOut);
    if( fResultResidual < 0.0 )
    {
        if( fPromptingResidual > -100.0 )
            printf("ERROR: Could not find user specified cell.\n");
        
        return -1.0; 
    }
    
    vTranMat3D(a3x3fRealOrientMatOut);
    
    // The orientation matrix is in reciprocal space.
    
    fInvMat3D(&a3x3fRealOrientMatOut[0][0], &a3x3fOrientMatOut[0][0]);
    
    oCrystal.nSetOrientMatrix(&a3x3fOrientMatOut[0][0]);
    
    oCrystal.vGetOrientAngles(pfOrientAngles);
    
    oCrystal.vGetCell(a6fUserCellFound);
    
    if( fPromptingResidual > -100.0 )
    {
        printf("\nSpec. Cell:  [%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f]\n"
                "Found Cell:  [%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f]\n"
                "Orientation: [%4.2f %4.2f %4.2f ]\n"
                "%% Error:      %4.2f\n",
                a6fKnownCell[0],a6fKnownCell[1],a6fKnownCell[2],
                a6fKnownCell[3],a6fKnownCell[4],a6fKnownCell[5],
                a6fUserCellFound[0],a6fUserCellFound[1],a6fUserCellFound[2],
                a6fUserCellFound[3],a6fUserCellFound[4],a6fUserCellFound[5],
                pfOrientAngles[0],pfOrientAngles[1],pfOrientAngles[2],
                100.0f*fabs(fResultResidual)
                );
        
        // If there is no user interaction expected _and_ the cell error is too large, return an error
        if( fabs(fResultResidual*100.0) > fabs(fPromptingResidual) && fPromptingResidual < 0.0 )
        {
            printf("ERROR: Could not find user specified cell within user specified error tolerance.\n");
            return -1.0;
        } 
        
        // It is not certain whether the user specified cell, or the
        // found cell should be selected.  If we are prompting, then
        // we can get the answer!
        
        // By default, we take the user input cell if it fit within the user specified error tol.
        // Else we use the cell we have found
        bUseFound = (fabs(fResultResidual*100.0) > ABS(fPromptingResidual)*2.0);

        while (1) 
        {
            cout << flush;
            printf("\nWrite found cell or specified cell to header?  \n"
                "             U = Write User  Cell   \n"
                "             F = Write Found Cell  \n"
                "       <Enter> = Default:    (%c)  ",
                bUseFound?('F'):'U');
            fflush(stdout);
            sTemp = "";
#ifdef SSI_PC
            if (fPromptingResidual>=0.0)
              sTemp = CCrclHelper::GetInstance()->sSendIndexUseFoundCell( a6fKnownCell,
                                                                          a6fUserCellFound, 
                                                                          pfOrientAngles, 
                                                                          100.0f * fResultResidual );
            else
              sTemp = "";
#else
            if (fPromptingResidual>=0.0)
                getline(cin, sTemp);                    
            else
                sTemp = "";
#endif
            sTemp.upcase();
            if (0 >= sTemp.length())
                break;
            else if ('U' == sTemp.GetAt(0))
            { 
                bUseFound = FALSE; 
                break;
            }
            else if ('F' == sTemp.GetAt(0))
            { 
                bUseFound = TRUE; 
                break; 
            } 
        }
    } 
    else
    {
        bUseFound = TRUE;
    }
        
    
    // We have the crystal orientation... Just change the cell.
    if( !bUseFound ) 
        oCrystal.vSetCell(a6fKnownCell);

    oCrystal.m_poSpacegroup->vSet(1);
    
    return fResultResidual;
}

                

void
Cindex::vSetResolution(const float fResoMin, const float fResoMax)
{
  // Set the min and max resolution in Angstroms.  Min is closest to
  // beam, max is as far away from beam, thus Min has a higher value
  // in Angstroms than Max

  m_fResolutionMin = max(fResoMin, fResoMax);
  m_fResolutionMax = min(fResoMin, fResoMax);

}

// Delete some reflns from the m_poReflnlist if sigma<=0.0, or
// resolution out of bounds.  Before Jan 2004, this stuff was
// in dtindex.cc.  Perhaps this routine should do the spot sharpness
// and I/sigmaI rejection also.
int Cindex::nDeleteReflns(const double fSigmaCutoff, const double fSharpnessCutoff)
{
    if( NULL == m_poReflnlist )
        return 0;

    int       nRef = m_poReflnlist->nGetNumReflns();
    if( 0 == nRef )
        return 0;

    int*      pnIndex = new int[nRef];
    Crefln*   poRefln = NULL;

    for(int nx = 0; nx < nRef; nx++)
    {
        poRefln = m_poReflnlist->poGetRefln(nx);

        if( 0.0 >= poRefln->fGetSigmaI() || 
            m_fResolutionMin < poRefln->fGetField(m_poReflnlist->m_nFI_fDetResolution) ||
            m_fResolutionMax > poRefln->fGetField(m_poReflnlist->m_nFI_fDetResolution) )
        {    
            pnIndex[nx] = 1;
        }
        else
            pnIndex[nx] = 0;
    }

    int       nDel = m_poReflnlist->nDelete(1, pnIndex);

    delete [] pnIndex;

    return nDel;
}

// RB 3/25/2004 The following 2 functions are a hack for getting the ranking information
// out to the header. 
void Cindex::vSetNumIndexedBySolution(int nSolution, int nIndexed)
{
    if( !m_poSolutions || nSolution < 0 || nSolution > m_poSolutions->nGetNumReflns() - 1 )
        return;

    (*m_poSolutions)[nSolution].vSetField(m_nFI_nNumIndexed, nIndexed);
}

void Cindex::vWriteRankingInfoToHeader(Cimage_header* poHdr, int nUsedForIndexing)
{
    if( !poHdr || !m_poSolutions || 0 == m_poSolutions->nGetNumReflns() )
        return;
    
    // Assuming that the first solution is the "chosen" solution
    int     nIndexed = (*m_poSolutions)[0].nGetField(m_nFI_nNumIndexed);
    
    int     nTemp[2] = {nUsedForIndexing, nIndexed};
    
    poHdr->nReplaceValue("DTINDEX_RANK_USED_REFLECTIONS", 2, nTemp); 
}
