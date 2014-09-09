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
// DTMainIntegrate.cpp  Initial author: RB                14-Jan-2005
//
//  This file contains the member functions of class CDTMainIntegrate
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
//    Example:  % dtintegrate sHeaderfile [ options ...]
//
//    Command line options:
//                sHeaderfile
//                -sigma  fSigma
//                -window nWin0 nWin1
//                -windowmultiplier nMult
//                -pad    nPad3D
//                -nonunf sNonunfFile
//                -reso   Resolution_min Resolution_max
//                -ref    sReflectionFileName
//                -seq    nStart nEnd
//                -num    nNumImage
//                -batch  nBatchImages
//                -batchprefix sBatchPrefix
//                -mosaicitymodel fMosMul fMosAdd
//                -nooverlap
//                -profit nNumReflns nNumImages [s[equential] | c[oncurrent] | l[ater]]
//                -profit
//                -noprofit
//                -nopurge
//                -mergeheaderscan, -differentscan
//                -mergeheaderdet, -differentdet
//                -wait   fSecs
//                -help
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Cintegrate.h"
#include "Crotation.h"
#include "CinterpMesh.h"

#undef max //mrp
#include <limits> //mrp

#if defined(DTREK_TRAP_FPE) && defined(IRIX)
#include <limits.h>
#include <sigfpe.h>
#endif

#ifdef PARALLEL
#include "par_dtrek.h"
#endif

#ifdef ANL_TIMER
#include "anlTimer.h"
#endif

//+Code begin

//+Definitions, constants, and initialization of static member variables

//+Public functions

#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#endif

#include "DTMainIntegrate.h"

#include "CCommandLine.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

//+ Code begin

extern Cstring* ms_asErrorsNames;

const char*     c_pcPrefind      = "-prefind";  

int CDTMainIntegrate::nExecute(unsigned int argc,      // Command line argument count
                               char        *argv[])    // Pointers to command line args

{
  int            nStat = 0;
  Cstring        sScanFile  = sDtrekGetPrefix() + "dtrefine.head";
  Cstring        sOut;
  Cstring        sStub;
  int            nNumImages;
  int            nSeqNumStart = -1;
  int            nProfitNumReflections = -1;
  int            nProfitMaxImageRange = -1;

  Cimage        *poImage      = NULL;

  Creflnlist    *poReflnlist  = NULL;
  Cintegrate    *poIntegrate  = NULL;

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  CXprop        *poXprop      = NULL;
#endif

  int            i, j;  // Loop counters
  int            nTemp;
  int            nTemp2;
  float          fTemp;
  Cstring        sTemp;
  Cstring        sMissingOption
                    = "ERROR: dtintegrate - missing option argument!";
  Cstring        sInvalidArgument
                    = "ERROR: dtintegrate - invalid argument!";

  bool           bDeleteScratch = TRUE;
  bool           bMergeHeadersScan = FALSE;
  bool           bMergeHeadersDet  = FALSE;

  float          a2fMosaicityModel[2];
  float          a2fResolution[2];
  float          a2fExcludeResoRing[2];
  float          a2fProfileSize[2];
  float          fDelHKLIn = 5.0;
  int            nProfitFlag = DTI_DO_NOT_PROFIT;

  nNumImages      = 0;

  a2fMosaicityModel[0] = 1.0;
  a2fMosaicityModel[1] = 0.0;

  a2fResolution[0] = 0.0;
  a2fResolution[1] = 99999.99f;

  a2fProfileSize[0] = 0.5;
  a2fProfileSize[1] = 0.3f;

#if defined(DTREK_TRAP_FPE) && defined(IRIX)
  // Underflow to zero
  sigfpe_[_UNDERFL].repls = _ZERO;
  // Overflow to max
  sigfpe_[_OVERFL].repls  = _MAX;
  // Trace first 5 exceptions
   sigfpe_[_UNDERFL].trace = 5;
   sigfpe_[_OVERFL].trace = 5;
   sigfpe_[_DIVZERO].trace = 5;
   sigfpe_[_INVALID].trace = 5;
   //   sigfpe_[_INT_OVERFL].trace = 5;

  // counts at end
  sigfpe_[_UNDERFL].count = INT_MAX;
  sigfpe_[_OVERFL].count = INT_MAX;
  sigfpe_[_DIVZERO].count = INT_MAX;
  sigfpe_[_INVALID].count = INT_MAX;
  //  sigfpe_[_INT_OVERFL].count = INT_MAX;

  // abort after 100
  sigfpe_[_UNDERFL].abort = 100;
  sigfpe_[_OVERFL].abort = 100;
  sigfpe_[_INVALID].abort = 100;
  //  sigfpe_[_INT_OVERFL].abort = 100;

  // abort on first divide by 0
  sigfpe_[_DIVZERO].abort = 1;

  //  handle_sigfpes(_ON, _EN_UNDERFL | _EN_OVERFL | _EN_DIVZERO | _EN_INVALID
  handle_sigfpes(_ON, _EN_OVERFL | _EN_DIVZERO
                 // | _EN_INVALID
                 //                      | _EN_INT_OVERFL
                 , 0, _ABORT_ON_ERROR, 0);
  cout << "handle_sigfpes just called\n" << flush;
#endif

  vDtrekSetModuleName("dtintegrate");
  vPrintCopyrightInfo();
  //cout << "\ndtintegrate:  Copyright (c) 2006 Rigaku\n";
  //cout << D_K_DTREKVersion << endl;

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

  vSetCommandArguments(argc, argv);

#ifdef SSI_PC
  if( 1 == argc ) 
  {
      DTREK_ERROR(0, "");
  }
#endif

  // Parse command line arguments

#ifdef PARALLEL
  cout << "\n***EXECUTING PARALLEL dtintegrate***\n\n" << flush;
#endif

#ifdef ANL_TIMER
  cout << "\n***EXECUTING TIMED    dtintegrate***\n\n" << flush;
  anl_reset_timer(30, "Entire dtintegrate process");
  anl_start_timer(30);
#endif

  argc--; argv++;

  nStat = 1;
  
  // Look at some command line arguments that need parsing BEFORE
  // creating a Cintegrate object

  for (i = 0; i < argc; i++)
    {
      sTemp = (const char*) argv[i];
      if ( ("-help" == sTemp) || ("-h" == sTemp) )
        {
          DTREK_ERROR(0, "");
        } 
      else if (i == 0)
        {
          sScanFile = sTemp;
        } 
      else if (   ("-mergeheadersscan" == sTemp)
               || ("-differentscan" == sTemp) )
        {
          bMergeHeadersScan = TRUE;
        }
      else if (   ("-mergeheadersdet" == sTemp)
               || ("-differentdet" == sTemp) )
        {
          bMergeHeadersDet = TRUE;
        }
      else if ("-nomergeheaders" == sTemp)
        {
          bMergeHeadersScan = FALSE;
          bMergeHeadersDet  = FALSE;
        }
      else if ("-seq" == sTemp)
        {
          // This will be used only if -mergeheadersdet or -mergheadersscan
          // is set as well.  The -seq option will be parsed below as well.

          i++;
          nStat = sscanf(argv[i], "%d", &nTemp);
          if (1 == nStat)
            {
              nSeqNumStart = nTemp;
            }
          else
            DTREK_ERROR(6, sInvalidArgument);
        }
      else if ("-prefix" == sTemp)
        {
          i++;
          if (i>=argc)
              DTREK_ERROR(5,sMissingOption);
          sStub = argv[i];
          nPutEnv((Cstring) "DTREK_PREFIX",sStub);          
        }
    };

  if (argc >= 1)
    {
      argc--;
      argv++;
    };
    
  sOut = sDtrekGetPrefix() + "dtintegrate.head";
  
  // Create default scan and integrate objects;
  if( bCreateHeader(sScanFile) )
    {
      poReflnlist  = new Creflnlist();
      
      if (bMergeHeadersScan || bMergeHeadersDet)
        {
          Cimage_header oMergeImgHeader;

          Cscan *poScan = new Cscan(*m_poHeader);

          if (!poScan->bIsAvailable())
            {
              cout << "ERROR in dtintegrate: Scan not available!\n";
              DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              // Get starting image number for the command line "-seq" option.
              
              if (0 > nSeqNumStart)
                {
                  cout << "ERROR in dtintegrate: First image sequence number not valid!\n";               
                  DTREK_ERROR(6, sInvalidArgument);
                }
              poScan->vSetSeqNum(nSeqNumStart);
              Cimage_header oMergeImgHeader(poScan->sGetImageName());
              if (!oMergeImgHeader.bIsAvailable())
                {
                  nStat = 2;
                  cout << "ERROR: dtintegrate - failed to merge image header info\n      from image\n "
                       << poScan->sGetImageName() << endl;
                  DTREK_ERROR(6, sInvalidArgument);
                }
              else
                {
                  cout << "\nINFO: dtintegrate - merging image header info from image\n"
                       << "      " << poScan->sGetImageName() 
                       << " into input header.\n" << flush;
                  Cstring sCopyMask;
                  sCopyMask = "SCAN*";
                  m_poHeader->nCopyMask(oMergeImgHeader, sCopyMask);
                  sCopyMask = "ROTATION*";
                  m_poHeader->nCopyMask(oMergeImgHeader, sCopyMask);
                  sCopyMask = "CRYSTAL_GONIO_VALUES";
                  m_poHeader->nCopyMask(oMergeImgHeader, sCopyMask);
                  sCopyMask = "SOURCE_WAVELENGTH";
                  m_poHeader->nCopyMask(oMergeImgHeader, sCopyMask);

                  if (bMergeHeadersDet)
                    {
                      sCopyMask = "*_GONIO_VALUES";
                      m_poHeader->nCopyMask(oMergeImgHeader, sCopyMask);
                    }
                }
            }
          delete poScan;
          poScan = NULL;
        }

      // The proper Scan info needs to be available before a 
      // Cintegrate object is instanced.

      poIntegrate  = new Cintegrate(m_poHeader, poReflnlist);
      
      poIntegrate->m_nProfitFlag = nProfitFlag;
      
      if ( poIntegrate->m_poScan->bIsAvailable() )
            {
              nStat = 0;
            }
    }
  
  if (0 != nStat)
    {
      DTREK_ERROR(nStat, "ERROR: dtintegrate error reading files!");
    }
   
      //  {
      //  Ccrystal oCrystal(*poIntegrate->m_poHeader);
      //  poIntegrate->m_fFixedMosaicity = oCrystal.fGetMosaicity();
      //  }


  cout << "\n\n" << flush;

  // From here on we have a *poIntegrate object 
  // and a poIntegrate->m_poScan object!

  // Set poIntegrate->m_anSeqNum[0], m_anSeqNum[1],
  // poIntegrate->m_poScan->vSetImgNum(k) before continuing
  // in case no -seq option was on the command line

  poIntegrate->m_anSeqNum[0] = poIntegrate->m_poScan->nGetSeqNum(0);
  nNumImages = poIntegrate->m_poScan->nGetSeqNumImages();
  if (1 >= nNumImages) 
    nNumImages = 720;
  if ( (1000 < nNumImages) && (nNumImages > poIntegrate->m_anSeqNum[0]) )
    {
      // The default for nNumImages is probably wrong.  
      // It probably represents the last sequence number in the scan
      poIntegrate->m_anSeqNum[1] = nNumImages;
    }
  else
    poIntegrate->m_anSeqNum[1] = poIntegrate->m_poScan->nGetSeqNum(nNumImages-1);
  nNumImages = 0;  // This MUST be reset to 0, but will be set below

  //  cout << "-seq default set to: -seq " << poIntegrate->m_anSeqNum[0]
  //   << " " << poIntegrate->m_anSeqNum[1] << endl;

  // Re-parse command line options

  for (i = 0; i < argc; i++)
    {
      sTemp = (const char*) argv[i];
      cout << "Command line string: >>" << sTemp << "<<" << endl;

      if ("-pad" == sTemp)
      {
          if (i + 1 < argc)
          {
              nStat = sscanf(argv[i+1], "%d", &nTemp);
              if (1 == nStat) {
                  i++;
                  poIntegrate->m_nPad3D = nTemp;
                  if ((i + 1 < argc) && (1 == sscanf(argv[i+1],"%f",&fTemp))) {
                      i++;
                      poIntegrate->m_fPad3D = fTemp;
                  };
              } else
                  DTREK_ERROR(6,sInvalidArgument);
          } else
              DTREK_ERROR(6, sInvalidArgument);
      } 
      else if ("-ref" == sTemp)
      {
          i++;
          if (i < argc)
          {
              poIntegrate->vSetReflnFilename((Cstring)argv[i]);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if("-ref_profit" == sTemp)
      {
          i++;
          if (i < argc)
          {
              poIntegrate->vSetReflnProfitFilename((Cstring)argv[i]);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-batchprefix" == sTemp)
      {
          i++;
          if (i < argc)
          {
              sTemp = (const char *)argv[i];
              if (4 < sTemp.length())
              {
                  sTemp = sTemp.substr(0, 4);  // Truncate to 4 characters
                  cout << "Batch prefix too long truncated to : " << sTemp << '\n';
              }
              poIntegrate->vSetBatchPrefix(sTemp);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
        else if ("-mosaicitymodel" == sTemp)
      {
          i++;
          for (j=0; (j < 2) && (i < argc); j++, i++)
          {
              nTemp = sscanf(argv[i], "%f", &fTemp);
              if (1 == nTemp)
                  a2fMosaicityModel[j] = fTemp;
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          i--;
          if (2 != j)
          {
              DTREK_ERROR(5, sMissingOption);
          }
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
      else if ("-window" == sTemp)
      {
          int nWin;
          int* pnWindow;
          i++;
          if (i < argc) {
              sTemp = argv[i];
              sTemp.upcase();
              if (sTemp == "MIN") {
                  pnWindow = &poIntegrate->m_a2nWinMin[0];
                  i++;
              } else if (sTemp == "MAX") {
                  pnWindow = &poIntegrate->m_a2nWinMax[0];
                  i++;
              } else
                  pnWindow = &poIntegrate->m_a2nWinMax[0];
          };
          for (j=0; (j < 2) && (i < argc); j++, i++)
          {
              nTemp = sscanf(argv[i], "%d", &nWin);
              if (1 == nTemp)
              {
                  pnWindow[j] = nWin;
              }
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          i--;
          if (2 != j)
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-ellipsoids" == sTemp)
      {
          poIntegrate->m_bWriteEllipsoids = TRUE;
      }
      else if ("-intmask" == sTemp)
      {
          poIntegrate->m_bWriteScanBitmap = TRUE;
      }
      else if ("-spotchase" == sTemp) 
      {
          poIntegrate->m_bSpotChase = TRUE;
      }
      else if ("-nospotchase" == sTemp) 
      {
          poIntegrate->m_bSpotChase = FALSE;
      }
      else if ("-prefix" == sTemp)
      {
          i++;
      }
      else if ("-windowmultiplier" == sTemp)
      {
          i++;
          if ((i >= argc) ||  (1!=sscanf(argv[i],"%f",&poIntegrate->m_fSpotSizeMultiplier)))
              DTREK_ERROR(6,sInvalidArgument);
      }
      else if ("-ring" == sTemp)
      {
          i++;
          for (j=0; (j < 2) && (i < argc); j++, i++)
          {
              nTemp = sscanf(argv[i], "%f", &fTemp);
              if (1 == nTemp)
                  a2fExcludeResoRing[j] = fTemp;
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          i--;
          if (2 != j)
          {
              DTREK_ERROR(5, sMissingOption);
          }
          poIntegrate->vAddExcludeResoRing(a2fExcludeResoRing[0], 
              a2fExcludeResoRing[1]);
      }
      else if ("-seq" == sTemp)
      {
        i++;
        for (j=0; (j < 2) && (i < argc); j++, i++)
          {
            nStat = sscanf(argv[i], "%d", &nTemp);
            if (1 == nStat)
              {
                poIntegrate->m_anSeqNum[j] = nTemp;
              }
            else
              DTREK_ERROR(6, sInvalidArgument);
          }
        i--;
        if (2 != j)
          {
            DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-noprofit" == sTemp)
      {
          poIntegrate->m_nProfitFlag = DTI_DO_NOT_PROFIT;
          nProfitFlag = DTI_DO_NOT_PROFIT;
      }
      else if ( ("-2D" == sTemp) || ("-profit" == sTemp) )
      {
          // With no extra arguments, profile fit externally later
          
          nProfitFlag =   DTI_DO_2D_PROFIT;
          
          // Check for extra parameters and 'c' or 's'
          
          if (    (i+2 < argc) 
              && (1 == sscanf(argv[i+1],"%d",&nTemp)) 
              && (1 == sscanf(argv[i+2],"%d",&nTemp2))) 
          {
              nProfitNumReflections = nTemp;
              nProfitMaxImageRange  = nTemp2;
              poIntegrate->m_nProfitNumReflections = nTemp;
              poIntegrate->m_nProfitMaxImageRange  = nTemp2;
              i += 2;
          }
          poIntegrate->m_nProfitFlag = nProfitFlag;
      }
      else if ("-num" == sTemp)         // Only valid in -scan mode
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
              {
                  nNumImages = nTemp;
                  poIntegrate->m_anSeqNum[1]
                      = poIntegrate->m_poScan->nGetSeqNum(nNumImages-1);
              }
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-wait" == sTemp)         // Only valid in -scan mode
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
              {
                  poIntegrate->vSetWaitTime((double)fTemp);
              }
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if("-out" == sTemp)
      {
          i++;
          if(i < argc)
          {
              sOut = (const char *)argv[i];
          }
          else
          {
              DTREK_ERROR(5,sMissingOption);
          }
      }
      else if ("-batch" == sTemp)
      {
          nTemp = -1;
          nTemp2 = -1;
          j = 0;
          i++;
          if (i < argc) {
              if (1 == sscanf(argv[i], "%d", &nTemp)) {
                  j++;
                  i++;
                  if (i < argc) {
                      if (1 == sscanf(argv[i], "%d",&nTemp2)) {
                          j++;
                      } else
                          i--;
                  } else
                      i--;
              }
              else
                  DTREK_ERROR(6, sInvalidArgument);
          } else          
              DTREK_ERROR(5, sMissingOption);
          
          if (j==1) 
              nTemp2 = nTemp;
          poIntegrate->vSetBatch(nTemp,nTemp2);      
      }
      else if ("-minpeakrad" == sTemp) 
      {
          i++;
          if ((i < argc) && (1 == sscanf(argv[i],"%d",&nTemp)))
              poIntegrate->m_nMinPeakRad = nTemp;
          else
              DTREK_ERROR(6,sInvalidArgument);
      }
      else if ("-fastpredict" == sTemp)
      {
          poIntegrate->m_bDoFastPrediction = TRUE;
      }
      else if ("-weakspotflag" == sTemp)
      {
          i++;
          if ((i < argc) && (1 == sscanf(argv[i],"%d",&nTemp)))
              poIntegrate->m_nIntegrateWeakFlag = nTemp;
          else
              DTREK_ERROR(6,sInvalidArgument);
      }
      else if ("-gain" == sTemp)
      {
          i++;
          if ((i < argc) && (1 == sscanf(argv[i],"%f",&fTemp)))
              poIntegrate->vSetGain(fTemp);
          else
              DTREK_ERROR(6,sInvalidArgument);

      } else if ("-sigmaa" == sTemp) 
      {
          i++;
          if ((i < argc) && (1 == sscanf(argv[i],"%f",&fTemp)))
              poIntegrate->m_fSigmaAboveBackground = fTemp;
          else
              DTREK_ERROR(6,sInvalidArgument);
      }
      else if ("-pcell" == sTemp) {
          i++;
          float fArgs[4];
          if ((i + 3 < argc) && 
              (1 == sscanf(argv[i],"%f",&fArgs[0])) &&
              (1 == sscanf(argv[i+1],"%f",&fArgs[1])) &&
              (1 == sscanf(argv[i+2],"%f",&fArgs[2])) &&
              (1 == sscanf(argv[i+3],"%f",&fArgs[3]))) {
              poIntegrate->m_fPCellAngle = fArgs[0];
              poIntegrate->m_a3fPCellVector[0] = fArgs[1];
              poIntegrate->m_a3fPCellVector[1] = fArgs[2];
              poIntegrate->m_a3fPCellVector[2] = fArgs[3];
              i+= 3;
          } else
              DTREK_ERROR(6,sInvalidArgument);
      }
      else if ("-noprerefine" == sTemp)
      {
          poIntegrate->vSetPrerefineFlag(0);
      }
      else if ("-prerefine" == sTemp)
      {
          if ((i + 1 < argc) && (1 == sscanf(argv[i+1],"%d",&nTemp)))
          {
              i++;
              poIntegrate->vSetPrerefineFlag(nTemp);
          } else
              poIntegrate->vSetPrerefineFlag(1);
      }
      else if ("-prefind" == sTemp)
      {
          // skip this option, because it will be parsed later
          while( i + 1 < argc && !bIsDTCommandLineOptionString(argv[i+1]) )
              i++;
      }
      else if ("-verbose" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                  poIntegrate->vSetVerboseLevel(nTemp);
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-display" == sTemp)
      {
          poIntegrate->m_nDisplay = 1;
#ifdef SSI_PC
          if ((i + 1 < argc) && (1 == sscanf(argv[i+1],"%d",&nTemp)))
          {
              i++;
              poIntegrate->m_nDisplayFreq = nTemp;
          } else
              poIntegrate->m_nDisplayFreq = 1;
#endif

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
          poXprop = new CXprop("DTDISPLAY_WINDOW_ID");
          poIntegrate->m_poXprop = poXprop;
#endif
      }
      else if ("-purge" == sTemp)
      {
          bDeleteScratch = TRUE;
      }
      else if ("-nopurge" == sTemp)
      {
          bDeleteScratch = FALSE;
      }
      else if (   ("-mergeheadersscan" == sTemp)
               || ("-differentscan" == sTemp) )
        {
          // Ignored, parsed above: bMergeHeadersScan = TRUE;
        }
      else if (   ("-mergeheadersdet" == sTemp)
               || ("-differentdet" == sTemp) )
        {
          // Ignored, parsed above: bMergeHeadersDet = TRUE;
        }
      else if ("-fixmosaicity" == sTemp) 
      {
          double f0;
          if ((i + 1>=argc) || (1!=sscanf(argv[i + 1],"%lf",&f0))) {
              f0 = -1.0;
          } else {
              i++;
          };
          if (f0 == -1.0)
              poIntegrate->m_bDetermineFixedMosaicity = TRUE;
          else {
              poIntegrate->vSetFixedMosaicity(f0);
          };
      }
      else if ("-checkoverlaps" == sTemp) 
      {
          double f0;
          i++;
          if ((i>=argc) || (1!=sscanf(argv[i],"%lf",&f0)))
              DTREK_ERROR(5,sMissingOption);
          poIntegrate->m_fMaxIntersectionFraction = f0; 
      }
      else if ("-nocheckoverlaps" == sTemp)
      {
          poIntegrate->m_fMaxIntersectionFraction = 1.0;
      }
      else if ("-rogue" == sTemp) 
      {
          int nMask = 0;
          int nNegMask = 0;
          int nImageNumber = -1;
          int nMinDim = -1;
          
          bool bNotFlag;
          Cstring sTemplate;
          Cstring sReflnlist;
          Cstring sArg;
          
          printf(
              "\n\n"
              "-rogue sOutputFileTemplate [[imagenumber nImage]        \n"
              "                    [mindim nDim]               \n"
              "                    [[not] sError ...]          \n"
              "                    [all] ...]  |               \n"
              "                    [file.ref]                  \n\n"
              );
          printf("-rogue error flag options (always printed when -rogue option is used)\n");
          for (j = 0; j < nMaxErrorStatusTypes; j++) {
              if (Crefln::ms_asErrorsNames[j].string()[0])
                  printf("%s\n",Crefln::ms_asErrorsNames[j].string());
          };
          printf("Use 'all' to include/exclude all flags.\n");
          printf("\n\n");
          
          i++;
          if (i>=argc) 
              DTREK_ERROR(5,sMissingOption);
          sTemplate = argv[i];
          if (!sTemplate.contains('?')) {
              printf("Output file template must contain at least one '?' character in -rogue\n");
              DTREK_ERROR(1,"");
          };
          i++;
          while (i<argc) {
              sArg = argv[i];
              if (sArg.contains('-'))
                  break;
              if (sArg.contains(".ref")) {
                  i++;
                  sReflnlist = sArg;
                  continue;
              } else if (sArg == "imagenumber") {
                  i++;
                  if ((i>=argc) || (1!=sscanf(argv[i],"%d",&nImageNumber))) {
                      printf("'imagenumber' requires a number\n");
                      DTREK_ERROR(1,"");
                  };
                  i++;
                  continue;
              } else if (sArg=="mindim") {
                  i++;
                  if ((i>=argc) || (1!=sscanf(argv[i],"%d",&nMinDim))) {
                      printf("'mindim' requires a number\n");
                      DTREK_ERROR(1,"");
                  };
                  i++;
                  continue;
              };
              
              if (sArg == "not") {
                  bNotFlag = TRUE;
                  i++;
                  if (i<argc)
                      break;
              } else
                  bNotFlag = FALSE;
              for (j = 0; j < nMaxErrorStatusTypes; j++) {
                  if (sArg == Crefln::ms_asErrorsNames[j])
                      break;
              };
              if (sArg == "all") 
                  j = 32; 
              else if (j == nMaxErrorStatusTypes) {
                  printf("Unrecognized error status type '%s'\n",sArg.string());
                  DTREK_ERROR(1,"");
              };
              if (bNotFlag) {
                  if (j == 32)
                      nNegMask = 0x7fffffff;
                  else
                      nNegMask |= 1 << j;
              } else {
                  if (j == 32)
                      nMask = 0x7fffffff;
                  else
                      nMask |= 1 << j;
              };
              i++;
          };
          i--;
          if (poIntegrate->nSetRogueData(nMask,nNegMask,sTemplate,sReflnlist,nImageNumber,nMinDim))
              DTREK_ERROR(1,"");
      }
      else if ("-diagnostic" == sTemp) 
      {
          if ((i+1<argc) && (argv[i+1][0]!='-')) {
              poIntegrate->m_sProfitDiagnostic = argv[i+1];
              i++;
          } else
              poIntegrate->m_sProfitDiagnostic = "diagnostic";

      }
      else if ("-socketport" == sTemp)
      {
          if ((i + 1 < argc) && (1 == sscanf(argv[i+1],"%d",&nTemp)))
          {
              i++;
              if (bOpenSocket(nTemp))
              {
                  poIntegrate->m_poSocket = m_poSocket;
              }
          }
          else
          {
              printf("-socketport option must have a port number following it.\n");
              DTREK_ERROR(6, sInvalidArgument);
          }
      }
      else if ("-help" == sTemp)
      {
          DTREK_ERROR(0, "");
      }
      else
      {
          DTREK_ERROR(1, "");
      }
    }

    ////////////////////////////////////////////////////
    std::vector<double>     daArgs;
    if( bGetCommandOption(c_pcPrefind, daArgs, " ") )
    {
        switch( (int)daArgs.size() )
        {
        case 1:
            poIntegrate->vSetPrefindFlag((int)daArgs[0]);
            break;
        case 2:
            poIntegrate->vSetPrefindFlag((int)daArgs[0]);
            poIntegrate->vSetPrefindSigma(daArgs[1]);
            break;
        default:
            poIntegrate->vSetPrefindFlag(1);  // By default do prefind on 1 batch
        }
    }
   /////////////////////////////////////////////////////////////////////////////

#ifdef SSI_PC
    cout << flush;
#endif
    //+JWP 2009-02-17
    // Move validation of -seq option to here from above -seq
    // Valid sequence start and sequence end numbers
    int k, l;
//mrp    for (k = 0; (k < 20000) &&
    for (k = 0; k < std::numeric_limits<int>::max() && //mrp
           (poIntegrate->m_poScan->nGetSeqNum(k)
//mrp            != poIntegrate->m_anSeqNum[0]);
            < poIntegrate->m_anSeqNum[0]); //mrp
         k++) {}
//mrp    if (20000 == k)
    if (poIntegrate->m_poScan->nGetSeqNum(k) > poIntegrate->m_anSeqNum[0]) //mrp
      {
        DTREK_ERROR(7, "ERROR: dtintegrate - invalid starting seqnum!");
      }
    else
      {
//mrp        for (l = k; (l < 20000) &&
        for (l = k; l < std::numeric_limits<int>::max() && //mrp
               (poIntegrate->m_poScan->nGetSeqNum(l)
//mrp                != poIntegrate->m_anSeqNum[1]);
                < poIntegrate->m_anSeqNum[1]); //mrp
             l++) {}
//mrp        if (20000 == l)
        if (poIntegrate->m_poScan->nGetSeqNum(l) > poIntegrate->m_anSeqNum[1]) //mrp
          {
            DTREK_ERROR(7, "ERROR: dtintegrate - invalid ending seqnum!");
          }
        else
          {
            // We need this check in case -num option was used above
            // because nNumImages is used below

            if (0 >= nNumImages)  
              nNumImages = l - k + 1;
            poIntegrate->m_poScan->vSetImgNum(k);
          }
      }
    //cout << "nNumImages = " << nNumImages << endl;
    //cout << "-seq : " << poIntegrate->m_anSeqNum[0]
    //<< ", " << poIntegrate->m_anSeqNum[1] << endl;
    //-JWP 2009-02-17


//   // RB: THIS DEBUG CODE IS JUST TO STOP THE DTREK PROCESS WHEN STARTED FROM CC:
//   int*  pNull = NULL;
//   *pNull = 1;
    
   poIntegrate->m_sOutputHeader = sOut;

  if (0.0 < a2fResolution[0])
    {
      poIntegrate->vSetResolution(a2fResolution[0], a2fResolution[1]);
    }
  poIntegrate->vSetMosaicityModel(a2fMosaicityModel[0], a2fMosaicityModel[1]);

  poIntegrate->nList();

  //
  // Done with command line arguments so now go integrate reflections!
  //

  Crotation *poRotation = NULL;

  cout << "dtintegrate: 3D method used" << endl << flush;
#ifdef PARALLEL
  nStat = par_init();
#endif
  nStat = poIntegrate->nIntegrate(nNumImages);
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  poIntegrate->m_poXprop = poXprop;
#endif
  (void) poIntegrate->m_poHeader->nWrite(poIntegrate->m_sOutputHeader);     


//+jwp 5-Jan-2001
/*
  if (1 > poIntegrate->nGetNumWritten())
    {
      cout << "Done. There were no profiles written.\n" << flush;
    }
  else
    {
      cout << "\nDone. dtintegrate: There were "
           << poIntegrate->nGetNumWritten() << endl << flush;
    }
*/
//-jwp 5-Jan-2001

  if (NULL != poReflnlist)
    delete poReflnlist;

  if (bDeleteScratch)
    {
      // Delete a variety of scratch files created/used by dtintegrate.
      // NOTE: This is NOT the same as the dtprocess Utils/Purge command.

      // Problem with deletes, if -display is set, then we often delete
      // the reflnlist file before it has been displayed, so wait 10 seconds
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
      if (NULL != poXprop)
        nDoSystemCommand("sleep 10");
#endif

      sTemp = sDtrekGetPrefix();
      (void) nFileDelete(sTemp + "dtintref*.ref*");
      (void) nFileDelete(sTemp + "dtrefinetmp*.ref");
      (void) nFileDelete(sTemp + "dtintpred*.ref");
#ifndef WIN32
      (void) nFileDelete(sTemp + "dtintegrate_*.head");
#endif

    }
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  if (NULL != poXprop)
    delete poXprop;
#endif

#ifdef ANL_TIMER
//  cout << "TIMERS STOPPED\n" ;
  anl_stop_timer(30);
  anl_print_timer(30);
#endif

  cout << "\ndtintegrate: All done.\n" << flush;
  cout << "\nINFO - If dtintegrate was successful, the suggested next step"
       << "\n       is to run dtcell to check your spacegroup hypothesis"
       << "\n       (also be sure look at the axial reflections listed above):\n"
       << "\n         dtcell " << poIntegrate->m_sOutputHeader << " " 
       << poIntegrate->sGetReflnProfitFilename()
       << "\n    -or-"
       << "\n       run dtscaleaverage with " << poIntegrate->m_sOutputHeader 
       << "\n       as the input .head file:\n"
       << "\n         dtscaleaverage " << poIntegrate->m_sOutputHeader << " ...\n\n" << flush; 

  if (NULL != poIntegrate)
    delete poIntegrate;

  return (nStat);
}

void CDTMainIntegrate::vError(const int nErrorNum, const Cstring& sMessage)
{
  char* pcErrorMessage = 

       "\ndtintegrate - Usage:\n"
       "dtintegrate sHeaderfile [ options ...]\n\n"
       "Options:     Description:\n\n"
       "sHeaderfile  File sHeaderfile (no default) has the scan\n"
       "             definition and the crystal, detector, and goniometer\n"
       "             properties.\n\n"
       "-profit [nNumReflections nMaxImageRangePerProfile] \n"
       "             Do profile fitting.  By default, no profile fitting\n"
       "             information is generated.  The -profit option (without arguments)\n"
       "             will not generate 2D profiles (which is not recommended).\n"
       "             Good values for the parameters are 50 7 for\n"
       "             macromolecules and 20 20 for small molecules.\n\n"
       "-pcell fAngle fVec0 fVec1 fVec2\n"
       "             Assume that there is a preasure cell mounted on the\n"
       "             phi goniometer stage.  Reject reflections that are\n"
       "             masked out by the preasure cell.\n\n"
       "-pad nPad [fPadMult]\n"
       "             Instructs dtintegrate to use nPad images on either side of\n"
       "             the block of images containing a predicted spot.\n"
       "             If fPadMult is specified, the algorithm will multiply the predicted\n"
       "             reflection width by fPadMult, and use the block of images fully\n"
       "             containing the new expanded width. \n"
       "             Default:  nPad=0, fPadMult=1.0\n\n"
       "-ref sReflectionFileName\n"
       "             Specifies an alternate output reflection file\n"
       "             instead of the default integrate.ref. The file format\n"
       "             is a standard d*TREK reflection file format.\n\n"
       "-seq nStart nEnd \n"
       "             Specifies an alternate starting and ending\n"
       "             image sequence values in a scan.  The default\n"
       "             is all images in a scan.\n\n"
       "-reso fResoMin fResoMax\n"
       "             Specifies the resolution range to predict\n"
       "             reflections for.  The default range is determined\n"
       "             by the detector position and source wavelength.\n\n"
       " -reso edge  Set upper resolution to be equal to the highest\n"
       "             *minimum* resolution found on any of the four\n"
       "             edges (sides) of the image.\n\n"
       " -reso corner\n"
       "             Set upper resolution to be equal to the highest\n"
       "             *maximum* resolution found on any of the four\n"
       "             edges (sides) of the image.\n\n"
       "-ring fRingResoMin fRingResoMax\n"
       "             Reflections with predicted centroids within the specified\n"
       "             resolution range 'ring' or annulus are flagged as bad.\n"
       "             This is useful to exclude integrating reflections in ice\n"
       "             rings.  Multiple -ring options may be used.\n\n"
       "-noprefrefine\n"
       "-prerefine [0|1|2|3|4|5|nNum]\n"
       "             Select a prerefinement option where reflections are predicted,\n"
       "             centroids determined, refinement occurs BEFORE integration.\n"
       "             This is extremely useful to get the best refined properties.\n"
       "             Options: 0  no prerefinement, same as -noprerefine\n"
       "                      1,2,3... try this many refinements before starting\n"
       "                      integration.\n"
       "             Default: -prerefine 1.\n"
       "-noprofit    Do not perform profile analysis and do not create profiles.\n\n"
       "-batch nImagesPerScaleBatch nImagesPerRefineBatch \n"
       "-batch nImagesPerBatch\n"
       "             nImagesPerScaleBatch specifies the rotation width\n"
       "             per scaling batch (in images).  nImagesPerRefineBatch specifies\n"
       "             how many images should get processed between refinements.  If\n"
       "             only one number (nImagesPerBatch) is provided, dtintegrate sets,\n"
       "             nImagesPerScaleBatch = nImagesPerRefineBatch = nImagesPerBatch.\n"
       "             Default: 1 4.\n\n"
       "-prefind  [n [sigma]]\n" 
       "-prefind  -n [sigma]\n" 
       "             Before integrating a scan, perform refinement after a peak search\n"
       "             on all images in the first n batches (default: 1) of the scan.\n"
       "             If n is negative, then |n| relatively-evenly-spaced-throughout-\n"
       "             the-scan images will be searched instead of consecutive images at\n"
       "             the beginning of the scan.  Also an attempt is made to include the\n"
       "             first and last images in the peak search (if -n is used).\n"
       "             Use 'sigma' standard deviations (default 5.0) above background\n"
       "             for pixels to be considered for peaks. Then perform refinement,\n"
       "             using the obtained peak list.\n\n"
       "-batchprefix sBatchPrefix\n"
       "             Specifies a prefix to use for batch ids given to reflections\n"
       "             in the output reflection list.  The prefix should not be\n"
       "             longer than 4 characters.  Default prefix is no prefix.\n\n"
       "-display     Tell a dtdisplay process to display images and\n"
       "             reflection lists as dtintegrate proceeds.\n\n"
       "-wait fWait  Maximum number of seconds to wait for a required image\n"
       "             before stopping.  Default 0.0, that is do not wait.\n\n"
       "-nopurge     Do not purge scratch files that were created (Default).\n\n"
       "-purge       Purge scratch files that were created.\n\n"
       "-differentscan or -mergeheadersscan\n"
       "             Merge the crystal goniometer, scan, rotation, and\n"
       "             source info from the image header into input header.\n"
       "             This option uses the (refined) detector position\n"
       "             from the input header.  It is suggested to use the -prefind\n"
       "             and -prerefine options with this option.\n\n"
       "-differentdet or -mergeheadersdet\n"
       "             Merge the crystal goniometer, detector goniometer,\n"
       "             scan, rotation, and source info from the image header\n"
       "             (i.e. implies -mergeheaderscan).\n"
       "             This option uses the (unrefined) detector position\n"
       "             from the image header. It is suggested to use -prefind\n"
       "             and -prerefine options with this option.\n\n"
       "-out headerfile\n"
       "             Write a header to headerfile at the end of integration.\n\n"
       "-window [Min|Max] nWin0 nWin1\n"
       "             Maximum or Minimum box or window size around a spot in pixels.\n"
       "             The actual window size will be determined automatically and\n"
       "             may be smaller than the values given.  If 0,0 then the Find\n"
       "             spot size will be used to determine the values.  Restricting\n"
       "             the max size can make the boxes smaller and thus make the\n"
       "             integration faster.  Boxes may include portions of neighboring\n"
       "             spots, but should not extend past the centers of the neighbors.\n"
       "             The background is determined from the non-spot area of the\n"
       "             window.  Default: 0, 0 which lets the program decide by using\n"
       "             the -windowmultiplier value times value of the spot size in the\n"
       "             DTFIND_SPOTSIZE in the input header.  Note: the -pad option\n"
       "             helps set the box size in the rotation direction.\n\n"
       "-windowmultiplier fMult\n"
       "             Default is 2.8.  If the -window is 0, 0 or not specified,\n"
       "             then this value times the ellipsoid values in the\n"
       "             input header is used to set the max window size.  If the\n"
       "             ellipsoid values are not found in the input header, then\n"
       "             it has a default value of 10, so the max box size is 28,28.\n\n"
       "-mosaicitymodel fMosMul fMosAdd\n"
       "             Adjust any refined mosaicity value MosIn:\n"
       "                MosOut = (MosIn * fMosMul) + fMosAdd\n"
       "             with defaults fMosMul = 1 and fMosAdd = 0, the MosIn remains\n"
       "             unchanged.  With values fMosMul = 0, the mosaicity remains\n"
       "             fixed at the value of fMosAdd.\n\n"
       "-fixmosaicity fMosaicity\n"
       "             Disable mosaicity refinement by specifying a fixed mosaicity.\n"
       "             Note that refined values of the mosaicity will still be displayed\n"
       "             during integration refinements, but will be ignored.\n\n"
       "-checkoverlaps fMaximumFractionOverlap\n"
       "             Reject spots that have a significant portion of their integrated\n"
       "             pixels intersecting other neighboring spots.  A spot is flagged as\n"
       "             overlapped if the fraction of intensity shared by other neighboring\n"
       "             spots is greater than or equal to  fmaximumFractionOverlap.\n"
       "             Default:  fMaximumFractionOverlap=0.03\n\n"
       "-nocheckoverlaps\n"
       "             Equivalent to -checkoverlaps 1.0\n Effectively disables\n"
       "             overlap checking\n\n"
       "-minpeakrad fMinPeakRad\n"
       "             Set minimum peak radius to fMinPeakRad. Default is 3 pixels.\n\n"
       "-fastpredict\n"
       "             Enable a 'faster' prediction algorithm for use during integration.\n\n"
       "-ellipsoids\n"
       "             Enable generation of profiling ellipsoid information.  Since\n"
       "             this takes additional processing time, it is disabled by default\n\n"
       "-spotchase\n" 
       "-nospotchase\n" 
       "             Enable/Disable spot chasing of mis-predicted spots.  There is\n"
       "             always some degree of spot-chasing occuring in dtintegrate, but\n"
       "             this option, when set, increases the boldness of the algorithm.\n"
       "             In cases where alpha1/alpha2 splitting is evident, it might be\n"
       "             useful to have this option on.  In cases of twinned spots, it might\n"
       "             be helpful to have this option off, unless the twin component\n"
       "             is specified in the input header.  In cases of a large unit cell\n"
       "             with spots close together the option should be off.\n"
       "             By default this option is set to disable spotchasing, EXCEPT for\n"
       "             full Rigaku RAPID images where the default is -spotchase.\n\n"
       "-rogue sOutputFile [[imagenumber nImage]        \n"
       "                    [mindim nDim]               \n"
       "                    [[not] sError ...]          \n"
       "                    [all] ...]  |               \n"
       "                    [file.ref]                  \n"
       "             Print shoeboxes containing any one of the given                \n"
       "             error messages to file with template sOutputFile.              \n"
       "             sOutputFile should contain '?' character. '-rogue' will        \n" 
       "             print out all shoeboxes that are processed 'strong'.           \n"
       "             If instead, file.ref is specified, then only shoeboxes matching\n"
       "             some HKL in file.ref will get printed.                         \n\n"
       "-socketport  [nPort]\n"
       "             Connect to a waiting socket on nPort and send update information.\n"
       "-help        Print this help text.             \n\n"
       "Information about scan angles and non-uniformity of response is\n"
       "determined from the headers of the input image files.  If the\n"
       "headers are incorrect, they may be edited with dtheaderedit or dtprocess.\n\n"
       "Examples:\n\n"
       " dtintegrate dtrefine.head -pad 2 -profit 50 7 -batch 1 4\n"
       " dtintegrate hivrt.head -seq 2 50 -differentdet -profit 50 7 -prerefine 6\n"
       " dtintegrate smallmol.head -seq 1 600 -window 50 50 -profit 20 20 -batch 50 50\n\n";
    
    if (0 != nErrorNum)
    {
      cout << sMessage << endl;
  } else
  printf(pcErrorMessage);
}
