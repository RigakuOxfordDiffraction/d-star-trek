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
// DTMainIndex.cpp  Initial author: RB           15-Mar-2004
//
//  This file contains the member functions of class CDTMainIndex
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

#include "DTMainIndex.h"

#include "Cindex.h"       // Includes many other include files
#include "Cimage.h"
#include "Csource.h"
#include "Cdetector.h"
#include "Ccrystal.h"
#include "Cgoniometer.h"
#include "CCellReduce.h"
#include "CEyeballIndex.h"

#if !defined(VC6)
using namespace std;
#endif

enum eIndexMethods {
  eUnknown_method,
  eFourier_method,
  eReciprocal_method
};

//+Description
//
//    Example:  % dtindex header_file ref_file [options]
//
//+ToDo
//
//   Error messages need implementing
//

int CDTMainIndex::nExecute(unsigned int argc, char* argv[])
{
  // RB temporarily
  //int*      pNULL = NULL;
  //*pNULL = 0;

  int            i, j, nx;       // Loop counters
  int            nRef;
  int            nStat;
  int            nTemp;
  float          fTemp;
  float          a3fTemp1[3];
  double         a3fTemp2[3];
  double         a3x3fTemp[3][3];
  double         f0,f1;
  Cstring        sTemp, sReflnIn, sOut;


  Cstring        sName, sNameOut;
  Cstring        sIndexOptions;
  Creflnlist    *poReflnlistIn;
  Cindex        *poIndex;
  Cspacegroup   *poSpacegroup;
  CCellReduce   oReduce;
  Ccrystal      oCrystal;

  float          fCellMax;
  float          fCellMaxB;
  bool           bCellMaxIn;

  eIndexMethods  eMethod;
  int            nNum;
  int            nChoose;
  bool           bPrompt = FALSE;
  bool           bDoList = FALSE;
  bool           bTweakObs = FALSE;
  float          fResidualMax = 3.0;
  const double   fInitRejectFraction = 0.10;
  double         fRejectFraction = fInitRejectFraction;

  double a6fCell[6], a6fErrors[6], a3fCellOrientationAngles[3];
  double a6fCellSave[6];
  double a6fCellLSQ[6];
  double a6fCellNew[6];
  double fArgument;
  double fInitMosaicity = 0;                    // Initial mosaicity from the crystal.
  Cstring sCent;
  float fGrid;
  int   nGrid;
  int   nVerbose;
  int   nSpacegroup;
  int   nNumDiffs;
  bool  bKnown;
  bool  bDPS = TRUE;
  bool  bDoBeamCheck = TRUE;
  bool  bDeleteIceRings = TRUE;
  int   nBeamSearchRadius = 10;
  int   nBeamValidRadius = 9;
  float a2fResolution[2] = {999999.00f, 0.00001f};
  double fSigma;                                // Min I/sigmaI to accept
  double fMin;
  double fSharpness;                            // Sharpness to reject.
  bool   bUseDiffVecs = FALSE;                  // Are we using difference vectors?

  // Twin variables.
  bool   bTwinCheck;                            // Take into account other crystals already in the header.
  bool   bTwinLawSearch;                        // Search for a twin law in one of the crystals already found.
  int    nTwinMultiplicity;                     // Multiplicity of solutions.  If this is !=1, then much of the default behaviour is changed.
  int    nTwinLawSearchMaxIndex;                // When searching for an alternate twin law, we construct rotation vectors.
  float  fTwinLawSearchMinFractionToIndex;      // Minimum fraction to index when searching for a twin comp.

  int    nUsedForIndexing = 0;
  
  Cstring        sMissingOption
                 = "ERROR: dtindex - missing option argument!";
  Cstring        sInvalidArgument
                 = "ERROR: dtindex - invalid argument!";
  //

  Cstring        sError("");


  vDtrekSetModuleName("dtindex");
  vPrintCopyrightInfo();
  //  cout << "\ndtindex:  Copyright (c) 2006 Rigaku\n";
  // ocout << D_K_DTREKVersion << endl;

#ifdef SSI_PC
  Cstring sCCAppVersion = (const char*) argv[argc-1];
  if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
	cout << sCCAppVersion;
	argc -= 2;
	cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
	argv[argc] = NULL;
  }
#endif

  // Copy command line to output log

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

  for (i = 1; i < argc; i++)
    {
      sIndexOptions = sIndexOptions + ' ' + (const char*) argv[i];
    }

  // Parse command line arguments

  argc--; argv++;

  nStat = 1;

  if (1 > argc)
    {
      // No header file argument and thats an error!

      DTREK_ERROR(1, "ERROR - no header filename!\n");
    }

  sTemp = (const char*) argv[0];
  if ( ("-help" == sTemp) || ("-h" == sTemp) )
    {
      DTREK_ERROR(0, "");
    }

  if (!bCreateHeader(sTemp) )
    {
      DTREK_ERROR(2, "ERROR - Header not available!\n");
    }

    m_poHeader->nClearDictionary("DTINDEX_RANK", "*");

    // Try to get a default cell from the header

  vZeroMat(6,1,a6fCell);
  vZeroMat(3,1,a3fCellOrientationAngles);
  (void) m_poHeader->nGetValue(Ccrystal::ms_sCrystalXUnitCell, 6, a6fCell);
  (void) m_poHeader->nGetValue(Ccrystal::ms_sCrystalXOrientAngles,3,a3fCellOrientationAngles);

  argc--;
  argv++;

  if (1 > argc)
    {
      // No reflection list arguments and thats an error!

      DTREK_ERROR(3, "ERROR - no reflection list filename!\n");
    }

  sReflnIn = (const char*) argv[0];

  // Read a reflection list from a file

  cout << "Reflection list: " << sReflnIn << endl << flush;
  poReflnlistIn = new Creflnlist(sReflnIn);

  if (!poReflnlistIn->bIsAvailable())
    {
      DTREK_ERROR(4, "ERROR: Reflection list not available!\n");
    }
  

  argc--;
  argv++;

  poIndex   = new Cindex(*m_poHeader);
  poIndex->vSetSubPeak(1.0f);

  // Parse rest of command line arguments


  // Parse any additional command line options

  eMethod     = eFourier_method;
  fCellMax    = 0.0;
  fCellMaxB   = 0.0;
  bCellMaxIn  = FALSE;
  fGrid       = 0.0;
  nGrid       = 0;
  nNumDiffs   = 0;
  nVerbose    = 1;
  nSpacegroup = 0;
  nChoose     = 1;
  nTwinMultiplicity = 1;
  nTwinLawSearchMaxIndex = 1;
  fTwinLawSearchMinFractionToIndex = 0.25;
  bTwinCheck = FALSE;
  bTwinLawSearch = FALSE;
  
  bKnown      = FALSE;
  sOut        = sDtrekGetPrefix() + "dtindex.head";
  fSigma            = 0.0f;
  fMin              = 0.0f;
  fSharpness        = 0.0f;

  for (i = 0; i < 6; i++)
    a6fErrors[i] = 0.0;


  for (i = 0; i < argc; i++)
    {
      sTemp = (const char*) argv[i];
      cout << "Command line string: >>" << sTemp << "<<" << endl;

      if ("-fourier" == sTemp)
        {
          bDPS    = FALSE;
          eMethod = eFourier_method;
        }
      else if ("-dps" == sTemp)
        {
          bDPS = TRUE;
          eMethod = eFourier_method;
        }
      else if ("-nodps" == sTemp)
        {
          bDPS = FALSE;
        }
      else if ("-list" == sTemp)
        {
          bDoList = TRUE;
        }
      else if ("-nofourier" == sTemp)
        {
          eMethod = eReciprocal_method;
        }
      else if ("-nodeice" == sTemp)
        {
          bDeleteIceRings = FALSE;
        }
      else if ("-known" == sTemp)
        {
          bKnown = TRUE;
        }
      else if ("-diffs" == sTemp)
        {
          bUseDiffVecs = TRUE;
        }
      else if ("-nodiffs" == sTemp)
        {
          bUseDiffVecs = FALSE;
        }
      else if ("-subpeak" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                {
                  poIndex->vSetSubPeak(fTemp);
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-minsep" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                {
                  poIndex->vSetMMSeparationMin(fTemp);
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-verbose" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                nVerbose = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-grid" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                fGrid = fTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-maxcell" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                {
                  fCellMax   = fTemp;
                  bCellMaxIn = TRUE;
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-maxresid" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                fResidualMax = fTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-spacegroup" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                nSpacegroup = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-sigma" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                fSigma = fTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
              nStat = 0;
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-sharpness" == sTemp) 
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i],"%f",&fTemp);
              if (1 == nStat)
                  fSharpness = fTemp;
              else
                  DTREK_ERROR(6,sInvalidArgument);
              nStat = 0;
          }
          else
          {
              DTREK_ERROR(5,sMissingOption);
          };
      }
      else if ("-reso" == sTemp)
        {
          i++;
          
          Cstring       sResoParameters("");
          while( i < argc && '-' != argv[i][0] )
          {
                sResoParameters += argv[i];
                sResoParameters += " ";
                i++;
          }

          Cstring       sError("");
          if( !bParseCalculateResolution(sResoParameters, a2fResolution[0], a2fResolution[1], sError) )
          {
              DTREK_ERROR(6, sError);
          }
          
          i--;
      }
      else if ("-num" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                nNumDiffs = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-choose" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                nChoose = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-out" == sTemp)
        {
          i++;
          if (i < argc)
            {
              sOut = (const char *) argv[i];
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-cell" == sTemp)
        {
          i++;
          for (j=0; (j < 6) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                {
                  a6fCell[j] = fTemp;
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (6 != j)
            {
              DTREK_ERROR(5, sMissingOption);
            }
          bKnown = TRUE;
        }
      else if ("-errors" == sTemp)
        {
          i++;
          for (j=0; (j < 6) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                {
                  a6fErrors[j] = fTemp;
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (6 != j)
            {
              DTREK_ERROR(5, sMissingOption);
            }
          bKnown = TRUE;
        }
      else if ("-prompt" == sTemp)
        {
          bPrompt = TRUE;
        }
      else if ("-notweakobs" == sTemp)
        {
          bTweakObs = FALSE;
        }
      else if ("-tweakobs" == sTemp)
        {
          bTweakObs = TRUE;
        }
      else if ("-multiplicity" == sTemp)
      {
          i++;
          if ((i < argc) && (1==sscanf(argv[i],"%f",&fTemp)) && (fTemp>=1) && (fTemp<=100))
              nTwinMultiplicity = (int) fTemp;
           else
                DTREK_ERROR(6, sInvalidArgument);
      }
      else if ("-twinindex" == sTemp)
        {
          bTwinCheck = TRUE;
        }
      else if ("-twinlawsearch" == sTemp)
      {
          bTwinLawSearch = TRUE;
          if ((i + 1 < argc) && (1==sscanf(argv[i + 1],"%d",&nTemp)) && (nTemp>=1) && (nTemp<=10)) {
              i++;
              nTwinLawSearchMaxIndex = nTemp;
              if ((i + 1 < argc) && (1==sscanf(argv[i + 1],"%f",&fTemp)) && (fTemp>=0) && (fTemp<=1.0)) {
                  fTwinLawSearchMinFractionToIndex = fTemp;
                  i++;
              };

          };
      }
      else if ("-beamcheck" == sTemp)
        {
          bDoBeamCheck = TRUE;
          i++;
          for (j=0; (j < 2) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
              {
                  if (j == 0)
                    nBeamSearchRadius = (int) fTemp;
                  else if (j == 1) 
                    nBeamValidRadius = (int) fTemp;
              }
              else
                  break;
            }
          i--;
        }
      else if ("-nobeamcheck" == sTemp) {
          bDoBeamCheck = FALSE;      
      } 
      else
        {
          DTREK_ERROR(6, "ERROR: dtindex - unknown command line option: " + sTemp);
        }
    }
#ifdef SSI_PC
    cout << flush;
#endif

    // Error checks.
    if ((bKnown) && (bTwinCheck || (nTwinMultiplicity>1))) {
        printf("ERROR:  Twin checking and known cell are incompatible options!\n");
        nStat = 1;
    } else
        nStat = 0;   

    poIndex->vSetReflectionList(poReflnlistIn);
    poIndex->nExpandReflnlist();
    poIndex->vSetResolution(a2fResolution[0], a2fResolution[1]);
    poIndex->nCalcReflnResolution(poReflnlistIn);

    // Delete reflns with sigma <= 0.0 and out of resolution bounds
    // as well as i/sigmai and sharpness cutoff

    i = poIndex->nDeleteReflns(fSigma, fSharpness);
    cout << "INFO: deleted " << i << " reflns outside of resolution bounds.\n"
         << "      This leaves " << poReflnlistIn->nGetNumReflns()
         << " reflns for indexing.\n" << flush;
    
    // Even before we start, we are going to delete potential ice rings.

    if (bDeleteIceRings) {
        double afIceRingReso[10] = { 3.917, 3.684, 3.458, 2.683, 2.261, 2.081, 1.9584,
            1.9272, 1.8927, 1.7292 };
        Csource oSource(*m_poHeader);
        int nNumStart = poReflnlistIn->nGetNumReflns();

	// Delete reflns in rings only in the resolution range required
	// (other rings will be deleted by the -reso option anyways)

        i = 0;
        while ( (i <= 9) && (a2fResolution[1] < afIceRingReso[i]) ) i++;
        poReflnlistIn->nDeleteRing(afIceRingReso, i, 1.0,
                                   oSource.m_poWavelength->fGetWavelength());
        int nNumAfter = poReflnlistIn->nGetNumReflns();
        if (0 >= nNumStart)
	  {
	    nStat     = -2;
	  }
        else if ( ((nNumStart - nNumAfter) * 100 / nNumStart) > 60)
          {
            // Deleted 60% or more of original reflections
            cout << "WARNING: too many reflns deleted by de-icing.\n"
                 << " To override: Consider changing resolution (-reso) or use the -nodeice option.\n\n";
            nStat = -5;
            delete poReflnlistIn;
            delete poIndex;
            return (nStat);
          }
    };

    if ((fSigma>0.0) && (fSigma<1.0)) {
        // Fractional I/sig value specified.
        double* pfIntensity;
        int*    pnIntensity;
        pfIntensity = new double[poReflnlistIn->nGetNumReflns()];
        pnIntensity = new int[poReflnlistIn->nGetNumReflns()];
        for (nx=0;nx<poReflnlistIn->nGetNumReflns();nx++) {
            f0 = (*poReflnlistIn)[nx].fGetIntensity();
            pfIntensity[nx] = f0;
            pnIntensity[nx] = nx;
        };
        g_pfCmpDoubles = &pfIntensity[0];
        qsort(&pnIntensity[0],poReflnlistIn->nGetNumReflns(),sizeof(int),double_cmp_rel);
        // Get the top 10 percent of the intensities.
        fMin = (*poReflnlistIn)[pnIntensity[max(0,(int) (1.0 - fSigma)*poReflnlistIn->nGetNumReflns())]].fGetIntensity();
        fSigma = 0.0;
        
        delete[] pfIntensity;
        delete[] pnIntensity;
        
    };

    // Do I/sig rejection 
    if ((fSigma>=1.0) || (fMin>=1.0)) {
        // Absolute I/sig value specified.
        for (nx=0;nx<poReflnlistIn->nGetNumReflns();nx++) {

            if ((fSigma>=1.0) && (((*poReflnlistIn)[nx].fGetSigmaI()<=0.0) ||
                (((*poReflnlistIn)[nx].fGetIntensity()/max(1.0,(*poReflnlistIn)[nx].fGetSigmaI()))<fSigma))) {
                poReflnlistIn->nDelete(nx--);
            } else if ((fMin>=1.0) && ((*poReflnlistIn)[nx].fGetIntensity()<fMin)) {
                poReflnlistIn->nDelete(nx--);
            };
        };
    } 
    // Do sharpness rejection.
    if (fSharpness>0.0) {
        // Absolute I/sig value specified.
      nStat = 0;
        for (nx=0;nx<poReflnlistIn->nGetNumReflns();nx++)
          {
            if ((poReflnlistIn->m_nFI_fObsSharpness>=0) && ((*poReflnlistIn)[nx].fGetField(poReflnlistIn->m_nFI_fObsSharpness)>fSharpness))
              {
                poReflnlistIn->nDelete(nx--);
                nStat++;
              }
          };
        printf("INFO: deleted %d reflections with sharpness > %.2lf\n", nStat, fSharpness);
        printf("INFO: there are %d reflections left for indexing.\n", poReflnlistIn->nGetNumReflns());
        nStat = 0;
    } 


    nUsedForIndexing = poReflnlistIn->nGetNumReflns();
    
    if( !nStat ) 
    {
        if( nUsedForIndexing < 8 )
        {
            printf("ERROR:  dtindex cannot index with %d reflections.\n      Terminating.\n", nUsedForIndexing);
            nStat = 1;
        } 
        else if ( nUsedForIndexing < 30 )
        {
            printf("WARNING:  dtindex might fail with only %d reflections available!\n", nUsedForIndexing);
        nStat = 0;
        }
    }
    

    CEyeballIndex oEyeIndex(*m_poHeader,*poReflnlistIn);
    

    if (!nStat) {
        if (bTwinLawSearch) {
            int nTwinID;
            for (nTwinID = 1; 1; nTwinID ++) {
                Ccrystal oCrystalIn(*m_poHeader,nTwinID);
                if (oCrystalIn.bIsAvailable()) {
                    Ccrystal oCrystalOut;
                    if (!oEyeIndex.nSearchForTwinLaw(oCrystalIn,oCrystalOut,nTwinLawSearchMaxIndex,fRejectFraction,fTwinLawSearchMinFractionToIndex)) {
                        printf("INFO: Twin law added to component %d\n",nTwinID);
                        if (fInitMosaicity)
                            oCrystalOut.vSetMosaicity(fInitMosaicity);
                        oCrystalOut.nUpdateHeader(m_poHeader);
                        m_poHeader->nWrite(sOut);
                        delete poIndex;
                        delete poReflnlistIn;
                        return 0;

                    };
                } else
                    break;
            };
        };
            
        
        if (bTweakObs)
        {
            oEyeIndex.nAddTweaks();
        } // if tweakobs
        
        // If there was a previous crystal mosaicity, get it.
        sTemp = D_K_CrystalPrefix D_K_CrystalXUnitCell;
        if (m_poHeader->bKeywordExists(sTemp)) {
            Ccrystal oCrystalIn(*m_poHeader);
            if (oCrystalIn.bIsAvailable())
                fInitMosaicity = oCrystalIn.fGetMosaicity();
        };
          
        // Delete any previous crystal solutions.
        if (!bTwinCheck) {
            sTemp = D_K_CrystalMask;
            m_poHeader->nDeleteMask(sTemp);
        };
        fArgument = -1.0;
        
        if (poReflnlistIn->nCalcRecipLatticePoints(*m_poHeader,fArgument)) {
            printf("ERROR: Something wrong with input reflection list!\n");
            nStat = 1;
        };
    };
    

    // Outer grand loop.  This allows 'restarts' of the program 
    while (!nStat) 
      {
        if (bDoList)
          poIndex->nList(1);
        if (bKnown)
          {
            // Use a known cell
            
            for (i = 0; i < 6; i++)
              a6fCellSave[i] = a6fCell[i];
            if (1000 < nSpacegroup)
              {
                                //+jwp 12-Dec-2002 No longer convert to reduced primitive cell!
                
                                // Spacegroup also given, (NOT from header, but from -spacegroup option!),
                                // so get the reduced primitive cell if it is not already primitive.
                                
                                poSpacegroup = new Cspacegroup (nSpacegroup);
                                Cstring sFlag;
                                sFlag = poSpacegroup->sNameFromNumber();
                                if ( ('P' != sFlag.GetAt(0)) && ('u' != sFlag.GetAt(0)) )
                                {
                                        // Not a primitive cell, so get the primitive cell
                                        // because it is easier to search for that.
                                        
                                        oReduce.nCentToPrim(sFlag.GetAt(0),a6fCell,false);
                                        oReduce.nReduceCell(a6fCell);
                                        printf("NON-PRIMITIVE INPUT CELL REDUCED FOR SEARCHING!\n");
                                }
                                delete poSpacegroup;
                        }
                        for (i = 0; i < 6; i++)
                        {
                if (0.0 >= a6fErrors[i])
                                        a6fErrors[i] = 0.03f * a6fCell[i];
                        }
            //poIndex->vSetKnownCell((float*)a6fCell);
            //poIndex->vSetKnownErrors((float*)a6fErrors);
            //poIndex->m_bKnownCell = TRUE;
            printf("\nSearching for a solution with a known cell:\n"
                "    a       b       c     alpha   beta  gamma\n"
                " %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f\n"
                " %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f <+-errors\n\n",
                a6fCell[0], a6fCell[1], a6fCell[2],
                a6fCell[3], a6fCell[4], a6fCell[5],
                a6fErrors[0], a6fErrors[1], a6fErrors[2],
                a6fErrors[3], a6fErrors[4], a6fErrors[5]);
                }
                
        poIndex->vSetReflectionList(poReflnlistIn);
        poIndex->nExpandReflnlist();

        for (nRef = 0; nRef<poReflnlistIn->nGetNumReflns(); nRef++) {
            (*poReflnlistIn)[nRef].vSetField(poIndex->m_nFI_nIndexSelect,1);
        };
        
        poIndex->vSetVerboseLevel(nVerbose);
        
        
        fTemp = 0.0;
        (void)poIndex->nCalcGetMaxCell(&fTemp);
        if (!bCellMaxIn)
        {
            fCellMax  = fTemp;
            fCellMaxB = fTemp;
        }
        else
            fTemp = fCellMax;
        
        cout << "Max cell length allowed for reciprocal lattice vectors: " << fTemp << endl << flush;
        
        if (0 == nNumDiffs)
        {
            // Here since may have different defaults for each method
            
            if ( bDPS)
                nNumDiffs = 1000;
            else if (eFourier_method == eMethod)
                nNumDiffs = 1000;
            else 
            {
                // Reciprocal space indexing                
                nNumDiffs = 50;
            }
        }
        
        // Next statement uses only nNumDiffs reflns!!!  We really want to use widely
        // spaced reflns in reciprocal space if possible, or the n most intense, but
        // still in 3D.
        
        nStat = poIndex->nCalcSumDiffs(bUseDiffVecs,nNumDiffs, fTemp);
        if (nStat)
            continue;                
        nStat = poIndex->nGroupVecs(0.02f);
        if (nStat)
            continue;

        
        if (bKnown)
        {
            // Check that max cell is in agreement with the known input cell.
            for (i = 0; i < 3; i++)
            {
                if (fCellMax < (a6fCell[i] + a6fErrors[i] + 1.0))
                {
                    fCellMax = a6fCell[i] + a6fErrors[i] + 1.0f;
                }
            }
        }
        poIndex->vSetCellLengthMax(fCellMax);
        
        if (eFourier_method == eMethod)
        {
            if (bDPS)
            {
                cout << "\nMethod:     " << "1D FFT with DPS algorithm";
            }
            else
            {
                cout << "\nMethod:     " << "3D direct space cosine Fourier transform";
            }
        }
        else
        {
            cout << "\nMethod:     " << "Standard reciprocal space indexing";
        }
        cout << "\nOut header: " << sOut
            << "\nMax cell:   " << fCellMax
            << "\nNum vecs:   " << nNumDiffs;
        cout << "\nSpacegroup: " << nSpacegroup
            << "\nVerbose:    " << nVerbose
            //       << "\nnChoose:    " << nChoose
            << '\n'
            << endl << flush;
        
        // Try one of 2 possible methods of indexing
        
        if (eFourier_method == eMethod)
        {
            if (bDPS)
            {
                cout << "Performing 1D FFT indexing (not cell reduction) with the DPS algorithm...see\n"
                    << "  Steller, Bolotovsky, & Rossmann (1997) J. Appl. Cryst. 30, 1036-1040.\n" 
                    << flush;
            }
            else
            {
                cout << "Performing direct space cosine Fourier transform indexing...\n"
                    << flush;
            }
            //      nNum  = min(200, nNumDiffs);
            nNum  = nNumDiffs;
            fTemp = poIndex->fCalcGetGrid();
            if (0.0 == fGrid)
            {
                fGrid = fTemp;
            }
            if (!bDPS)
                cout << "Suggested grid interval is: " << fTemp << '\n' << flush;
            if (0 == nGrid)
            {
                nGrid = 75;
                if (bDPS) nGrid = 100;
            }
            fTemp = (float) nGrid * fGrid;
            float fGridFactor;
            fGridFactor = 75.0f;
            if (bDPS) fGridFactor = 100.0f;
            if (!bCellMaxIn)
            {
                // If the shortest diff vec suggests a maxcell that can still
                // fit within 100 grid points, then set maxcell to that value.
                
                float fInvDiffMax = 200.0f;
                if (0.0 != poIndex->fGetDiffVecMax())
                    fInvDiffMax = 1.0f / poIndex->fGetDiffVecMax();
                
                //          cout << "fTemp, fGrid, DifvecMax: " << fTemp << ' '
                //               << fGrid << ' ' << poIndex->fGetDiffVecMax() << ' '
                //               << fInvDiffMax << endl;
                if (   (fGrid * fGridFactor > fInvDiffMax * 1.02)
                    && (fTemp >  fInvDiffMax * 1.02) )
                {
                    fCellMax = fTemp;
                }
                else if (   (fGrid * 100 > fInvDiffMax * 1.02)
                    && (fTemp <  fInvDiffMax * 1.02) )
                {
                    fCellMax = fInvDiffMax * 1.02f;
                }
                else if ((fGrid * 100) <= fCellMaxB)
                {
                    fCellMax = fCellMaxB;
                }
            }
            if (!bDPS)
                cout << "Suggested max cell is:      " << fTemp << flush;
            nGrid = (int) (fCellMax / fGrid + 0.5);
            if (100 < nGrid)
            {
                // WARNING, is this what we want to do?
                // Testing shows there is no reason to have more than 100 grid points!?
                
                nGrid = 100;
                fGrid = fCellMax / (float) nGrid;
            }
            if (500.0 < fCellMax)
            {
                nGrid = 125;
                fGrid = fCellMax / (float) nGrid;
            }
            if (!bDPS)
            {
                cout << "\nActual grid interval is:    " << fGrid;
                cout << "\nActual max ";
            }
            else
            {
                cout << "\nMax ";
            }
            cout << "cell is:         " << fCellMax << endl << flush;
            
            if (fCellMax < poIndex->fGetDiffVecAvg())
            {
                cout << "WARNING! The actual max cell is less than the reciprocal\n"
                    << "         of some of the shortest difference vectors.\n"
                    << "         You may wish to use max cell of "
                    << (int)poIndex->fGetDiffVecAvg() + 10 << " Angstroms." << endl << flush;
            }
            poIndex->vSetCellLengthMax(fCellMax);

            
            if (bDPS)
                nStat = poIndex->nIndexDPS(nNum,false);
            else
                nStat = poIndex->nIndexCos(0.0, fCellMax, nGrid, nNumDiffs);

            if (nStat)
                continue;
            
        }  
                oEyeIndex.m_nTwinID = 0;
                
        
        // Load the solutions.
        poIndex->nLoadSolutionArray();

        
        // Refine the solutions that were found.  Use CEyeballIndex to do this.
        {
                    
            Ccrystal oCrystalIn;
            Ccrystal oCrystalOut;
            int nWorked;
            oEyeIndex.m_nPrint = 2;
            do {
                nWorked = 0;


                for (nx=0;nx<poIndex->m_poSolutions->nGetNumReflns();nx++) 
                {
                    poIndex->nLoadSolution(nx,&oCrystalIn);
                    
                    int     nIndexed = 0;
                    f0 = oEyeIndex.fRefineSolution(fRejectFraction, oCrystalIn, oCrystalOut, &nIndexed); 
                    
                    poIndex->vSetNumIndexedBySolution(nx, nIndexed);

                    if (f0 != 0.0) {

                        if (bKnown) {
                            f0 = poIndex->fGetKnownCellResidual(a6fCell,oCrystalOut,-100.0);
                            if (f0 >= 0.0)
                                f0 = 1.0/f0;
                            else 
                                f0 = 0.0;
                        };
                        if (f0 != 0.0) {                        
                            nWorked++;                        
                            poIndex->nPutSolution(nx,&oCrystalOut,f0);
                        };
                    } 
                };


                if (nWorked<poIndex->m_poSolutions->nGetNumReflns()/4) {
                    fRejectFraction /= 2.0;
                    if (fRejectFraction<0.05)
                        fRejectFraction = 0.0;
                };
            } while ((nWorked<poIndex->m_poSolutions->nGetNumReflns()/4) && (fRejectFraction!=0.0));

            if (nWorked == 0) {
                nStat = 1;
                break;
            };

            oEyeIndex.m_nPrint = 0;
            // Sort the reflections based on residual, and only keep the top ones.
            if (poIndex->nPruneSolutions(nTwinMultiplicity)) {
                if (bTwinCheck) {
                    printf("INFO:  Failed to find a twin component.\n");
                    nStat = 0;
                    break;
                } else {
                    printf("ERROR:  No solutions are satisfactory!!!");
                    nStat = 1;
                    break;
                };
            };

            printf("\n\n");
            if( bDoBeamCheck )
            {
                printf("Executing beam refinement with\n");
                printf("beam search radius, acceptable shift radius: %d %d ...\n", nBeamSearchRadius, nBeamValidRadius);
                
                poIndex->nLoadSolution(0,&oCrystalIn);                
                
                double      dBeamPixChange = 0.0;
                
                sError = "";
                oEyeIndex.bRefineBeam(fRejectFraction,
                                      nBeamSearchRadius,
                                      nBeamValidRadius,
                                      !bTwinCheck,
                                      oCrystalIn, 
                                      oCrystalOut,
                                      &dBeamPixChange,
                                      f0,
                                      sError);

                if( f0 == 0.0 && "" != sError )
                {
                    // Catastrophic failure of the beam center refinement
                    Cstring     strTemp("ERROR:");
                    strTemp += sError;
                    strTemp += "\n";

                    printf(strTemp.string());
                    
                    nStat = 1;
                    break;
                }
                
                poIndex->nPutSolution(0, &oCrystalOut, 0.0);
                bDoBeamCheck = FALSE;
                
                if( 5.0 < dBeamPixChange ) 
                {
                    printf("\nWARNING!: Beam position moved more than 5 pixels!\n=======\n");
                    fflush(stdout);
                }
                
                if( 0.1 < dBeamPixChange )
                {
                    fArgument = -1.0;
                    poReflnlistIn->nCalcRecipLatticePoints(*m_poHeader, fArgument);
                    
                    printf("INFO: Restart with adjusted beam center.\n");
                    fflush(stdout);
                    
                    for (i = 0; i < 6; i++)
                        a6fCell[i] = a6fCellSave[i];
                    
                    continue;
                }
            }
        }
        
        poIndex->vWriteRankingInfoToHeader(m_poHeader, nUsedForIndexing);

        // Cell reduction starts here.       
        oReduce.m_fResidualMax = fResidualMax;
        
        if (NULL != poIndex->m_poSolutions)
                {
            if (4 <= poIndex->nGetVerboseLevel())
                        {
                sNameOut = "soln.A";
                (void) poIndex->m_poSolutions->nWrite(sNameOut, poIndex->m_pnSoln);
                        }
            
            if( bKnown )
            {
                nStat = poIndex->nLoadSolution(0, &oCrystal);
                
                double  fPromptingResidual = bPrompt ? fResidualMax : fResidualMax * -1.0;

                if( 0.0 > poIndex->fGetKnownCellResidual(a6fCell, 
                                                         oCrystal, 
                                                         fPromptingResidual) )
                {
                    nStat = 1;
                    break;
                }
                else
                {
                    if (nSpacegroup>0)
                        oCrystal.m_poSpacegroup->vSet(nSpacegroup);

                    poIndex->nPutSolution(0,&oCrystal,f0 = 0.0);
                }
            } // end of bKnown group
            
            if ( 0 == nStat && nTwinMultiplicity>1) 
              {
                // A multiplicity > 1 indicates that we are loading a bunch of twin solutions.
                int nSoln;
                int nBest;
                Ccrystal oCrystalIn;
                Ccrystal oCrystalOut;
                Ccrystal oCrystalTemporary;
                
                
                for (nSoln = 0; nSoln < poIndex->m_poSolutions->nGetNumReflns(); nSoln ++) {
                    oReduce.oGetSolutions().vDeleteAll();
                    poIndex->nLoadSolution(nSoln,&oCrystalIn);
                    oReduce.nAddCrystal(oCrystalIn);
                    oReduce.nCalcGetBestLattices();

                    nStat = oReduce.nListResults(&nBest);
                    if (nSpacegroup>0) 
                        nBest = oReduce.nSelectMatch(nSpacegroup) + 1;
                    if (0 == nStat) 
                        nStat = oReduce.nLoadSolution(nBest-1, &oCrystalOut);
                    if (0 == nStat) {
                        poIndex->nLoadSolution(nSoln,&oCrystalTemporary);
                        f0 = oReduce.fGetOrientInfo(oCrystalTemporary,oCrystalOut);
                        poIndex->nPutSolution(nSoln,&oCrystalOut,f0 = 0.0);
                    };                             
                    if (nStat) {
                        poIndex->m_poSolutions->nDeleteSorted(nSoln,poIndex->m_pnSoln);
                        nSoln--;
                    };
                };
            };
            
            // Perform a least-squares fit of the found cell 
            // to the 44 different lattices

            if ((!bKnown) && (nTwinMultiplicity==1))
              {
                
                int nRelist;
                int nReCalc = 1;
                int nBest = -1;
                do 
                {
                    
                    
                    if (nReCalc) {
                        Ccrystal oCrystalIn;
                        Ccrystal oCrystalOut;
                        int nPass;
                        
                        nPass = 0;
                        oEyeIndex.m_nPrint = 1;
                        do {
                            oReduce.oGetSolutions().vDeleteAll();
                            if (!nPass) {
                                for (nx = 0; nx < poIndex->m_poSolutions->nGetNumReflns(); nx++) {
                                    poIndex->nLoadSolution(nx,&oCrystalIn);
                                    oReduce.nAddCrystal(oCrystalIn);
                                };
                            } else {
                                // Add the un-doubled cell from the previous iteration.
                                oReduce.nAddCrystal(oCrystalIn);
                            };
                            oReduce.nCalcGetBestLattices();
                            
                            // Look at the triclinic cell to see if it has a doubled axis.
                            // If so, we need to reitereate.
                            oReduce.nLoadSolution(oReduce.oGetSolutions().nGetNumReflns()-1, &oCrystalOut);
                            nPass++;
                        } while (oEyeIndex.nComputeDoubledAxes(fRejectFraction,oCrystalOut,oCrystalIn));
                    };
                    
                    nRelist = 0;
                    nReCalc = 0;
                    
                    if (0 == nStat)
                    {
                        // If the user specified a spacegroup, then use that.
                        if (nSpacegroup>0) {
                            nBest = oReduce.nSelectMatch(nSpacegroup) + 1;
                            fTemp = oReduce.oGetSolutions()[nBest-1].fGetIntensity();
                            if (fTemp>fResidualMax) 
							{
								printf("WARNING:  Solution %d has residual > %5.2f\n"
                                        "          resetting maximum residual for listing solutions.\n",nBest,fResidualMax);                          
								oReduce.m_fResidualMax = fTemp;
								nStat = oReduce.nListResults(&j);
								fResidualMax = fTemp;
                            } 
							else
                                nStat = oReduce.nListResults(&j);
                        } 
						else
                            nStat = oReduce.nListResults(&nBest);
                        
                    }
                    if (0 == nStat)
                    {
                        // Select one of the Bravais lattices based on the residual or
                        // the desired spacegroup
                        
                        if (bPrompt)
                        {
                            // Prompt option on, prompt for choice
                            
                            if (100.0f > fResidualMax)
                            {
                                cout << "\nTo view least-squares fits to other lattices, "
                                    << "\nenter a new percent residual between "
                                    << max((int)fResidualMax, oReduce.nGetNumSolutions()+1) 
                                    << " and 100% "
                                    << "at the following prompt.\n";
                            }
                            printf("\n   PLEASE ANSWER THE FOLLOWING PROMPT:\n");
                            printf("\nEnter solution number (Soln num) of your choice (1-14) from the above table\n"
                                   " or a new limiting %% residual to show more (but worse) solutions (>= %d"
                                   "%%)"
                                   "\nNote: Good solutions usually have 'LeastSq resid(%%)' less than 0.5%% to 1.0%%."
                                   "\n   or enter L to get a lattice character listing.\n",
                                   oReduce.nGetNumSolutions()+1);
			    fflush(stdout);

                            printf(" <Enter> or <cr> or %d will select the '%c %s' lattice.\n",
                                   nBest, 
                                   oReduce.sGetLattice(nBest-1,false).string()[1],
                                   oReduce.sGetLattice(nBest-1,true).string());
                            
#ifdef SSI_PC
                            sTemp = CCrclHelper::GetInstance()->sSendIndexTable( oReduce );
#else
                            printf("Select> ");
                            fflush(stdout);
                            sTemp = "";
                            getline(cin, sTemp);
#endif
                            sTemp.upcase();
                            if (0 >= sTemp.length())
                                i = nBest;
                            else if ((sTemp.length()==1) && (sTemp.GetAt(0)=='L'))
                            {
                                oReduce.nListResultsL();
                                nRelist = -1;
                                i = -1;
                            } 
                            else if (oReduce.nParseSelection(sTemp,i))
                                i = atoi(sTemp.string());

                            
                            // Allow re-listing if 
                            // input > oReduce.nGetNumSolutions()
                            
                            if (oReduce.nGetNumSolutions() < i)
                            {
                                nRelist = i;
                                fResidualMax = (float) nRelist;
                            }
                        }
                        else if (0 >= nSpacegroup)
                        {
                            // Spacegroup <= 0, use program choice
                            
                            i = nBest;
                        }
                        else
                        {
                            // Spacegroup > 0, use the one that matches the spacegroup
                            
                            i = oReduce.nSelectMatch(nSpacegroup) + 1;
                        }
                        
                        if (0 == nRelist)
                        {
                                nStat = oReduce.nLoadSolution(i-1, &oCrystal);
                            
                            if (0 < nSpacegroup)
                            {
                                // If nSpacegroup is specified (i.e. not 0), and
                                // it is consistent with the symmetry of selected solution
                                
                                poSpacegroup = new Cspacegroup (nSpacegroup);
                                
                                if (  (poSpacegroup->cGetClass()
                                    == oCrystal.m_poSpacegroup->cGetClass())
                                    && (poSpacegroup->cGetCentFlag()
                                    == oCrystal.m_poSpacegroup->cGetCentFlag())
                                    )
                                {
                                    oCrystal.m_poSpacegroup->vSet(nSpacegroup);
                                    cout << "Input spacegroup " << nSpacegroup
                                        << " is consistent with selected solution, so used.\n" << flush;
                                }
                                else
                                {
                                    cout << "WARNING: Input spacegroup " << nSpacegroup
                                        << " MISMATCHES selected solution, so ignored!\n" << flush;
                                }
                                delete poSpacegroup;
                            }
                            if (0 == nStat)
                            {
                                // List unit cell and spacegroup ...
                                
                                // Now get possible orientation angles ...
                                
                                // First save spacegroup number in nGrid because the
                                // nGetUBMatrix() step uses the K field in the solutions
                                // for something other than spacegroup
                                // and the nLoadSolutions step uses the K field for the
                                // spacegroup number.
                                
                                nGrid = oCrystal.m_poSpacegroup->nGet();
                                
                                // Get cell into local variables a6fCell and a6fCellLSQ
                                
                                oCrystal.vGetCell(&a6fCellLSQ[0]);                                                   
                                oCrystal.vGetCell(&a6fCell[0]);
                                
                                cout << "...determining orientation angles...\n" << flush;
                                
                                j = 1;
                                Ccrystal        oCrystalIn;
                                Ccrystal        oCrystalOut;
                                Ccrystal        oCrystalTemporary;
                                
                                // Need to loop through all poIndex solutions to get the best one.
                                f1 = 1.0;
                                
                                for (j = 0; j < poIndex->m_poSolutions->nGetNumReflns(); j++)
                                {
                                    poIndex->nLoadSolution(j, &oCrystalTemporary);
                                
                                    oCrystal.vGetOrientMatrix(&a3x3fTemp[0][0]);
                                    
                                    oCrystalIn.nSetOrientMatrix(&a3x3fTemp[0][0]);
                                    
                                    // RB Aparently, this Thad's function calculates the orientation matrix for the reduced cell.
                                    // The original OM is in oCrystalTemporary. The resultant OM is out in oCrystalIn.
                                    f0 = oReduce.fGetOrientInfo(oCrystalTemporary, oCrystalIn);
                                    
                                    printf(".");
                                    
                                    if( f0 < f1 )
                                    {
                                        oCrystalIn.vGetOrientMatrix(&a3x3fTemp[0][0]);
                                        
                                        oCrystalOut.nSetOrientMatrix(&a3x3fTemp[0][0]);
                                    
                                        f1 = f0;
                                    }
                                }
                                
                                printf("\n");
                                
                                nStat = 0;
                                
                                if (0 == nStat)
                                {
                                    
                                    // List possible orientation angles ...
                                    // and choose just one of them...
                                    
                                    oReduce.nLoadRotationChoices(-1,&oCrystalOut);                              
                                    oReduce.nListResultsB();
                                                                                                     
                                    // Get the solution with the lowest differnece between the
                                    // input angle in the header.
                                    for (i = 0; i < oReduce.oGetSolutions().nGetNumReflns(); i++) {
                                        oReduce.nLoadSolution(i,&oCrystalIn);
                                        oCrystalIn.vGetOrientAngles(a3fTemp1);
                                        if ((i == 0) || (fabs(a3fTemp1[0] - a3fCellOrientationAngles[0]) +
                                            fabs(a3fTemp1[1] - a3fCellOrientationAngles[1]) +
                                            fabs(a3fTemp1[2] - a3fCellOrientationAngles[2]))<
                                            (fabs(a3fTemp2[0] - a3fCellOrientationAngles[0]) +
                                            fabs(a3fTemp2[1] - a3fCellOrientationAngles[1]) +
                                            fabs(a3fTemp2[2] - a3fCellOrientationAngles[2]))) {
                                            nChoose = i + 1;
                                            vCopyVec3D(a3fTemp1,a3fTemp2);
                                        };
                                    };
                                    
#ifdef SSI_PC
                                    CCrclHelper::GetInstance()->vSendIndexOrientTable( oReduce, nChoose );
#endif

                                    if (bPrompt || (0 > nChoose))
				      {
                                        // Prompt user for choice
					cout << "\n   PLEASE ANSWER THE FOLLOWING PROMPT:\n";
					cout << "\nEnter your choice 0=Abort [" << nChoose << "]: " << flush;
					
                                        sTemp = "";
                                        getline(cin, sTemp);
                                        if (0 < sTemp.length()) {
                                            nChoose = atoi(sTemp.string());
#ifdef SSI_PC
                                            if( nChoose >= 0 ) {
                                                    cout << nChoose << "\n";
                                            }
                                            else {
                                                    cout << "<<Cancelled>>\n";
                                            }
#endif
                                        }
                                    }
                                    if (0 == nChoose)
                                    {
                                        // Allow for a re-listing here here, too.
#ifdef SSI_PC
										if( fResidualMax < 1.0 ) {
											nRelist = 1;
										}
										else {
											nRelist = (int) fResidualMax;
										}
#else
                                        nRelist = (int) fResidualMax;
#endif
                                        nReCalc = 1;
                                        nStat = 0;
                                    }
#ifdef SSI_PC
                                    else if (-5 == nChoose)
                                    {
                                        // Allow for proper "canceling" of results dialog in CrystalClear (bug 425)
                                        nRelist = 0;
                                        nReCalc = 0;
                                        nStat = -5;
                                        return (nStat);
                                    }
#endif
                                    else
                                    {
                                        nStat = oReduce.nLoadSolution(nChoose-1, &oCrystal);
                                        
                                        // Restore saved least squares fit cell
                                        
                                        oCrystal.vSetCell(a6fCellLSQ);
                                        
                                        // Restore saved spacegroup number and cell parameters
                                        
                                        oCrystal.m_poSpacegroup->vSet(nGrid);
                                        cout << "Orientation angles choice " << nChoose << " selected.\n" << flush;
                                        oCrystal.nList(1);
                                        oCrystal.vGetCell(a6fCellNew);
                                        if (   ('m' == oCrystal.m_poSpacegroup->cGetClass())
                                            && (90.0 > a6fCellNew[4]) )
                                        {
                                            cout << "WARNING, monoclinic cell choice has beta angle < 90!\n" << flush;
                                        }
                                    }
                                }
                        }
                    }
                }
                
                if (0 < nRelist)
                {
                    oReduce.m_fResidualMax = fResidualMax;
                }
            } while (0 != nRelist);
            
            if (0 != nStat)
            {
                cout << "WARNING! Could not find orientation angles!\n" << flush;
                cout << "You may wish to:\n"
                    << "  check that the detector position is input correctly,\n"
                    << "  check that the beam position is input correctly,\n"
                    << "  change maximum cell with the -maxcell option,\n"
                    << "  or get different spots for indexing.\n" << flush;
                break;
            } 
            else
                poIndex->nPutSolution(0,&oCrystal,f0 = 0.0);
        }
        
        if (0 == nStat)
        {
                
            if (bTwinCheck) {
                // Try to load the default oCrystal object with the same settings as all new solutions
                // This includes mosaicity and shift vectors.
                oCrystal.nInitValues(*m_poHeader,1);
            };

            for (nx = 0; nx < poIndex->m_poSolutions->nGetNumReflns(); nx++) {
                oCrystal.vSetTwin(0);   // Informs that twin componant updated will be next available slot.
                poIndex->nLoadSolution(nx,&oCrystal);
                if (fInitMosaicity)
                    oCrystal.vSetMosaicity(fInitMosaicity);
                nStat = oCrystal.nUpdateHeader(m_poHeader);
            };
            
            m_poHeader->nReplaceValue(Cindex::ms_sDtindexOptions, sIndexOptions);
            if (0 == nStat)
            {
                nStat = m_poHeader->nWrite(sOut);
                if (0 == nStat)
                  {
                    cout << "dtindex - Wrote header file " << sOut.string() << endl;
                    cout << "\nINFO - The above indexing solution is ONLY a hypothesis."
                         << "\n       One must confirm the hypothesis by examining symmetry"
                         << "\n       of observed intensities as well as possible systematic"
                         << "\n       absences."
                         << "\n       One does this by collecting data, processing the data"
                         << "\n       and using dtcell and dtscaleaverage.\n";
                    cout << "\nINFO - The suggested next step is to run dtrefine with"
                         << "\n       " 
                         << sOut.string() 
                         << " as the input .head file: \n" 
                         << "\n       dtrefine " 
                         << sOut.string() 
                         << " " 
                         << sReflnIn.string() 
                         << " ...\n\n" 
                         << flush;
                  }
                else
                {
                    cout << "dtindex - ERROR writing header file " << sOut.string() << endl << flush;
                    break;
                }                
                
            }
            else
            {
                cout << "dtindex - ERROR updating header file!\n";
            	break;
            }
        }
    }
    else
    {
        cout << "No solutions found!\n";
        nStat = -1;
        break;
    }
    
    nStat = 0;
    break;
  };
  
  delete poIndex;
  delete poReflnlistIn;
  
  return nStat;
}

void CDTMainIndex::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
  cout << "\ndtindex - Usage:\n"
        "dtindex header_file refln_file [ options ...]\n\n"
        "Command line options: Description:\n\n"
        " header_file          Name of a file with a d*TREK header used to\n"
        "                      initialize the indexing procedure.\n\n"
        " refln_file           Name of a d*TREK reflection list file (no default)\n"
        "                      which will be used for indexing.\n\n"
        " -fourier             Use the direct space cosine Fourier transform method.\n\n"
        " -nofourier           Use a non-Fourier method.\n\n"
        " -dps                 Use the DPS indexing method.  This is the default.\n\n"
        " -diffs (-nodiffs)    Do (NOT) use difference vectors.  -nodiffs is the default.\n\n"
        " -spacegroup nSpace   Spacegroup to use.  The default of 0, lets dtindex decide.\n\n"
        " -cell a b c alp bet gam\n"
        "                      Use this as the default cell instead of the one\n"
        "                      in the input d*TREK header.  This also sets -known.\n\n"
        " -errors erra errb errc erralp errbet errgam\n"
        "                      Used as the default errors on cell parameters.\n"
        "                      Only applicable if -known is also used.\n\n"
        " -grid fGrid          The Fourier grid size in Angstroms.\n"
        "                      The default is chosen automatically.\n\n"
        " -known               Try to find a cell which matches the one in the\n"
        "                      header or specified by the -cell option.  Use\n"
        "                      the -spacegroup option to specify a spacegroup\n"
        "                      number different from the default of 1.\n\n"
        " -maxcell fCellMax    Specifies the max cell length to search for.\n"
        "                      The default is derived from the crystal to\n"
        "                      detector distance and minimum spot separation\n"
        "                      (see -minsep below).\n\n"
        " -maxresid fResidMax  Largest permissible least squares residual for a reduced\n"
        "                      cell.  All lattices with a residual higher than fResidMax\n"
        "                      are not displayed.  Default: 3%.\n\n"
        " -num nMaxNumDiffs    Maximum number of difference vectors to use in indexing.\n\n"
        " -out header_out      New header file with results.  Default: dtindex.head.\n\n"
        " -prompt              Prompt user for solution choices.  Default: do not prompt.\n\n"
        " -list                List the detector position information.\n\n"
        " -minsep fPixels      Minimum separation between spots in millimeters.\n"
        "                      This helps determine the maximum cell length to search\n"
        "                      for.  Default: 0.65 mm.\n\n"
        " -sigma fSigma        If fSigma>=1.0, then it specifies the minimum I/sigmaI\n"
        "                      required for indexing contributors.  If 0.0<fSigma<1.0\n"
        "                      then it specifies the fraction of the input list\n"
        "                      reflections to use.\n\n"
        " -nodeice             Prevents deletion of spots in resolution rings close\n"
        "                      to the resolution of the powder rings of hexagonal ice.\n"
        "                      Default: delete spots thought to be in ice rings.\n\n"
        " -reso fResoMin fResoMax\n"
        "                      Specifies the resolution range for indexing.\n"
        "                      The default is entire resolution range of the\n" 
        "                      reflection list.\n\n"
        " -reso edge           Set upper resolution to be equal to the highest\n"
        "                      *minimum* resolution found on any of the four\n"
        "                      edges (sides) of the image. The image detector\n"
        "                      information will be taken from the input header.\n\n"
        " -reso corner         Set upper resolution to be equal to the highest\n"
        "                      *maximum* resolution found on any of the four\n"
        "                      edges (sides) of the image. The image detector\n"
        "                      information will be taken from the input header.\n\n"
        " -sharpness fSharp    Filters out 'spots' that have a sharpness greater\n"
        "                      than fSharp.  This tends to improve indexing and\n"
        "                      refinement, since it excludes garbage and/or\n"
        "                      amorphous spot profiles.\n"
        "                      Default:  fSharp = 0.2\n\n"
        " -twinindex           Enables twin indexing.  Twin indexing works best\n"
        "                      when the input reflection list is generated from\n"
        "                      the -rejects command dtrefine.  This will augment the\n"
        "                      existing twin solutions with new indexing solution(s).\n\n"
        " -multiplicity nMultiplicity\n"
        "                      Generate at up to nMultiplicity different indexing\n"
        "                      solutions and store them as separate components.\n\n"
        " -nobeamcheck\n"
        " -beamcheck [nSearchRadius] [nValidRadius]\n"
        "                      After finding a pre-reduced cell solution, an algorithm is\n"
        "                      employed to verify the validity of the cell.  This is done\n"
        "                      by attempting to reindex the data with different beam\n"
        "                      centers located within nSearchRadius pixels of the input\n"
        "                      beam center.  If a local maximum is not found within\n"
        "                      nValidRadius pixels of the original center, dtindex\n"
        "                      returns with an error message. To disable this check,\n"
        "                      use the -nobeamcheck option. Default:-beamcheck 10 9\n\n"
        " -tweakobs            Enable tweaking of estimated refln observed rotation\n"
        "                      midpoint values.\n\n"
        " -verbose nLevel      Verbose level.  High values give more output.  Default: 1.\n\n"
        " -help                Print this help text.\n\n"
        "Information about scan angles and the detector is\n"
        "determined from the headers of the input image files.  If the\n"
        "headers are incorrect, they may be edited with dtheaderedit.\n\n"
        "Example:\n"
        "   dtindex dtprocess.head dtfind.ref -prompt\n\n" << flush;
}

