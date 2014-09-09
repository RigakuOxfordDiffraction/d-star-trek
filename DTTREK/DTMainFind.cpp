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
// DTMainFind.cpp       Initial author: RB          15-Mar-2004
// This file contains the member functions of class CDTMainFind

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

#include "DTMainFind.h"

#include "Cfind.h"
#include "Crotation.h"
#include "dtarray.h"
#include "Cstat.h"
#include "CBeamData.h"
#include "CBackground.h"
#include "CScanBitmap.h"
#include "CIceRingLocator.h"

#include "CCommandLine.h"

#undef max //mrp
#include <limits> //mrp

#if !defined(VC6)
using std::cout;
using std::flush;
using std::endl;
#endif

CDTMainFind::CDTMainFind()
{
    m_nNumberOfResoBins = 10;
    
    m_nMinReflnsForRankStats = 5;
}

//+Description
//
//    Example:  % dtfind sHeaderFile [ options ...]
//
//    Command line options:
//                sHeaderFile
//                -sigma  fSigma
//                -min    fMin
//                -peakmin fPeakMin
//                -circle nCen0 nCen1 nMinrad nMaxrad
//                -rect   nMinPx0 nMinPx1 nExtPx0 nExtPx1
//                -window nWin0 nWin1
//                -filter nPeakFilter
//                -bloops nLoops
//                -reso   reso1 reso2
//                -out    sOutHeader
//                -ref    sReflectionFileName
//                -seq    nStart nEnd
//                -num    nNumImage
//                -dump   nDumpRefln
//                -display
//                -help
//
//+ToDo
//
//   Error messages need implementing
//

int CDTMainFind::nExecute(unsigned int argc, char* argv[])
{
  int          nStat;
  Cstring      sScanFile        = sDtrekGetPrefix() + "dtprocess.head";
  Cstring      sReflnFile       = sDtrekGetPrefix() + "dtfind.ref";
  Cstring           sReflnFile2   = sDtrekGetPrefix() + "dtfind_2d.ref";

  Cstring      sOutHead      = "";
  Cstring      sFindOptions;
  eFind_methods eMethod, eMethod2;

  Cimage      *poImage       = NULL;
  Cfind       *poFind        = NULL;

#if (!defined(SSI_PC) && !defined(NO_X_WINDOWS))
  CXprop      *poXprop       = NULL;
  CXprop      *poXprop2D     = NULL;
#endif
  Creflnlist  *poReflnlist   = NULL;
  int          nNumImages    = 1;
  int          nEveryNthImage = 1;  // This default is overridden by sGetEnv
  int          nNumSeqOption = 0;
  int          nFirstReflnLastSeq = 0;
  int          a20nSeqInfo[20][4];
  float        a2fResolution[2];
  bool         bResolutionProvided;
  int          nMinReflections = 0;
  int          nMinImages = 0;
  int          nTotalImages = 0;
  bool         bSparseFind = FALSE;
  bool         bFindBeamCenter = FALSE;
  bool         bFindBeamCenterDiagnostics = FALSE;
  
  bool         bPurgeTmpFiles = TRUE;
  bool         bPrintTables = FALSE;

  bool         bRank = false;
  
  bool         bFindRings = false;
  float        fFindRingsMin2Theta = DT_RING_LOCATOR_TWO_THETA_START;
  float        fFindRingsWidthToExclude = DT_RING_LOCATOR_EXPECT_RING_WIDTH;
  int          nFindRingsImageSizeCompression = DT_RING_LOCATOR_IMAGE_SIZE_COMPRESSION;
  
  float        fInitBeamCenter0,fInitBeamCenter1;
  Cstat        oBeamStat0,oBeamStat1;
  int          nFirstReflnThisImage;
  int          nFirstReflnThisSeqOption;
  int          nSeqCount;
  Cstring      sSeqOptions = "";
  Cstring      sIntMask = "";


  int          i, j, k;  // Loop counters
  int          nTemp;
  float        fTemp;
  Cstring      sTemp;

  Cstring      sMissingOption
                 = "ERROR: dtfind - missing option argument!";
  Cstring      sInvalidArgument
                 = "ERROR: dtfind - invalid argument!";

  a2fResolution[0] = 999999.99f;  // Default minimum resolution
  a2fResolution[1] = 0.1f;        // Default maximum resolution
  bResolutionProvided = false;
  //+2011-07-21
  if ("" != sGetEnv("DTFIND_EVERY"))
    nEveryNthImage = (int)atoi(sGetEnv("DTFIND_EVERY").string());
  //-2011-07-21


  vDtrekSetModuleName("dtfind");
  vPrintCopyrightInfo();
  //cout << "\ndtfind:  Copyright (c) 2006 Rigaku\n";
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

  // Copy command line options to sFindOptions
  for (i = 1; i < argc; i++)
    {
      sFindOptions = sFindOptions + ' ' + (const char*) argv[i];
    }

  // Parse command line arguments

  argc--; argv++;

  nStat = 1;

  if (1 <= argc)
    {
      sTemp = (const char*) argv[0];
      if ( ("-help" == sTemp) || ("-h" == sTemp) )
        {
          DTREK_ERROR(0, "");
        }
      else
        {
          sScanFile = sTemp;
        }
      argc--;
      argv++;
    }

  // Create default scan and find objects;
  if( bCreateHeader(sScanFile) )
    {
      poReflnlist = new Creflnlist();
      poFind      = new Cfind(*m_poHeader, poReflnlist);
      if (poFind->m_poScan->bIsAvailable())
        {
          // Default method is unknown.

          poFind->m_eMethod  = eFind_unknown_method;
          eMethod            = eFind_unknown_method;
          eMethod2           = eFind_unknown_method;
          nStat = 0;
        }
      poFind->m_fMinimum = 50.0;
      poFind->m_nEveryNthImage = nEveryNthImage; // No effect as of 2011-07-21,
                                // but is needed for Cfind::nList()
    }

  if (0 != nStat)
    {
      DTREK_ERROR(nStat, "ERROR: dtfind error reading files!");
    }

  cout << "\n\n" << flush;

  // From here on we have a *poFind object and a poFind->m_poScan object!
  // Parse any additional command line options

  for (i = 0; i < argc; i++)
    {
      sTemp = (const char*) argv[i];
      cout << "Command line string: >>" << sTemp << "<<" << endl;

      if ("-sigma" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                poFind->m_fSigma = fTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-min" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                poFind->m_fMinimum = fTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-nopurge" == sTemp)
	{
	  bPurgeTmpFiles = FALSE;
	}
      else if ("-peakmin" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                poFind->m_fPeakMinimum = fTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-bloops" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                poFind->m_nBLoops = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-nosaturated" == sTemp)
      {
          poFind->m_bIncludeSat = FALSE;
      }
      else if ("-circle" == sTemp)
        {
          i++;
          for (j=0; (j < 4) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                poFind->m_a4nCircle[j] = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (4 != j)
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-rect" == sTemp)
        {
          i++;
          for (j=0; (j < 4) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                poFind->m_a4nRect[j] = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (4 != j)
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

          bResolutionProvided = true;
          i--;
        }    
      else if ("-resobins" == sTemp)
        {
          i++;
          if (i < argc)
            {
              if( 1 == sscanf(argv[i], "%d", &nTemp) )
                m_nNumberOfResoBins = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
          //////////////////////////////////////////////////////////////////////////////////////
          i++;     // the second parameter of -resobins is optional. Let's see if it is there...
          if (i < argc && !bIsDTCommandLineOptionString(argv[i]) )
            {
                  if( 1 == sscanf(argv[i], "%d", &nTemp) )
                    m_nMinReflnsForRankStats = nTemp;
                  else
                    DTREK_ERROR(6, sInvalidArgument);
            }
            else
                i--;  // The second argument is not there, so let's go back.
        }    
      else if ("-window" == sTemp)
        {
          i++;
          for (j=0; (j < 2) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                poFind->m_a2nSpotWindow[j] = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (2 != j)
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-getPMTinfo" == sTemp)
        {
          i++;
          for (j=0; (j < 2) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                poFind->m_a2nPMTInfo[j] = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (2 != j)
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-filter" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                poFind->m_nPeakFilter = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-pad" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                nTemp = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
          printf("WARNING: -pad option is now obsolete.  Please use -2D or -2Dmerge\n");
        }
      else if ("-minrefs" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%d",&nTemp);
              if (1 == nStat)
                  nMinReflections = abs(nTemp);
              else
                  DTREK_ERROR(6,sInvalidArgument);
              bSparseFind = (nTemp<0);
          } 
          else
          {
              DTREK_ERROR(6,sMissingOption);
          }
      }
      else if ("-minimages" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%d",&nTemp);
              if (1 == nStat)
                  nMinImages= nTemp;
              else
                  DTREK_ERROR(6,sInvalidArgument);
          } 
          else
          {
              DTREK_ERROR(6,sMissingOption);
          }
      }
      else if ("-dump" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                poFind->m_nDumpRefln = nTemp;
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-tables" == sTemp)
      {
          poFind->m_bIncludeBackground = TRUE;
          bPrintTables = TRUE;
      }
      else if ("-ref" == sTemp)
        {
          i++;
          if (i < argc)
            {
              sReflnFile = (const char *)argv[i];
                          sReflnFile2 = sReflnFile.before(".");
                          sReflnFile2 += "_2d.ref";
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
              sOutHead = (const char *)argv[i];
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-seq" == sTemp)
        {
          i++;
          
          for (j=0; (j < 2) && (i < argc); j++, i++)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                {
                  poFind->m_a2nSeqNum[j]        = nTemp;
                  a20nSeqInfo[nNumSeqOption][j] = nTemp;
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          i--;
          if (2 != j)
            {
              DTREK_ERROR(5, sMissingOption);
            }
          
          int       k = 0;
          int       l = 0;
          
//mrp          for (k = 0; (k < 20000) && (poFind->m_poScan->nGetSeqNum(k) != poFind->m_a2nSeqNum[0]);k++)
		  for (k = 0; k < std::numeric_limits<int>::max() && poFind->m_poScan->nGetSeqNum(k) < poFind->m_a2nSeqNum[0]; k++) //mrp
          {
          }
          
//mrp          if (20000 == k)
          if (poFind->m_poScan->nGetSeqNum(k) > poFind->m_a2nSeqNum[0]) //mrp
            {
              cout << "ERROR: dtfind: first seqnum in scan is " 
                   << poFind->m_poScan->nGetSeqNum(0) << ",\n"
                   << "               but dtfind -seq is " 
                   << poFind->m_a2nSeqNum[0]
                   << endl;
              DTREK_ERROR(7, "ERROR: dtfind - invalid starting seqnum!");
            }
          else
            {
//mrp              for (l = k; (l < 20000) && (poFind->m_poScan->nGetSeqNum(l) != poFind->m_a2nSeqNum[1]);l++)
              for (l = k; l < std::numeric_limits<int>::max() && poFind->m_poScan->nGetSeqNum(l) < poFind->m_a2nSeqNum[1]; l++) //mrp
              {
              }
              
//mrp              if (20000 == l)
              if (poFind->m_poScan->nGetSeqNum(l) > poFind->m_a2nSeqNum[1]) //mrp
                {
                  DTREK_ERROR(7, "ERROR: dtfind - invalid ending seqnum!");
                }
              else
                {
                  nNumImages = l - k + 1;
                  poFind->m_poScan->vSetImgNum(k);
                  a20nSeqInfo[nNumSeqOption][2] = k;
                  a20nSeqInfo[nNumSeqOption][3] = nNumImages;
                  sSeqOptions += " -seq ";
                  sSeqOptions += a20nSeqInfo[nNumSeqOption][0];
                  sSeqOptions += " ";
                  sSeqOptions += a20nSeqInfo[nNumSeqOption][1];
                }
            }
          if (20 > nNumSeqOption)
            nNumSeqOption++;
          
        }
      else if ("-2D" == sTemp)
        {
          eMethod2 = eMethod;
          eMethod  = eFind_2D_method;
                  
        }
      else if ("-2Dmerge" == sTemp)
        {
          eMethod2 = eMethod;
          eMethod  = eFind_2D_merge_method;
          printf("WARNING: -2Dmerge is not recommended at the present time.\n");
        }
      else if ("-3D" == sTemp)
        {
          eMethod2 = eMethod;
          //eMethod  = eFind_2D_merge_method;
          eMethod  = eFind_2D_method;
          printf("WARNING:  3D method is now obsolete.  Reset to -2D.\n");
        }
      else if ("-beamcenter" == sTemp) 
      {
          bFindBeamCenter = TRUE;
          if ((i+1 < argc) && (((Cstring) "diagnostic") == argv[i+1])) {
              bFindBeamCenterDiagnostics = TRUE;
              i++;
          };
          if ((i+2 >= argc) || 
              (1!=sscanf(argv[i+1],"%f", &fInitBeamCenter0)) ||
              (1!=sscanf(argv[i+2],"%f", &fInitBeamCenter1))) {
              fInitBeamCenter0 = -1;
              fInitBeamCenter1 = -1;
          } else
              i+=2;
      }
          else if ("-ringdetect" == sTemp)
      {
          bFindRings = TRUE;
          if ((i+1 < argc) && (1 == sscanf(argv[i+1],"%f",&fTemp))) {
                
              fFindRingsWidthToExclude = fTemp;

              if ((i+2 < argc) && (1 == sscanf(argv[i+2],"%f",&fTemp))) {
                  fFindRingsMin2Theta = fTemp;
                  if ((i+3 < argc) && (1 == sscanf(argv[i+3],"%d",&nTemp))) {
                      nFindRingsImageSizeCompression = nTemp;
                      i++;
                  };

                  i++;
              } 
              i++;
          } 
      }
      else if ("-every" == sTemp)
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                {
                  nEveryNthImage = nTemp;
		  // Next line has no effect as of 2011-07-21,
                  // but is needed for Cfind::nList()
                  poFind->m_nEveryNthImage = nEveryNthImage; 
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-num" == sTemp)
        {
          // Only valid in 3D mode

          i++;
          if (i < argc)
            {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                {
                  nNumImages = nTemp;
                  poFind->m_a2nSeqNum[1]
                    = poFind->m_poScan->nGetSeqNum(nNumImages-1);
                }
              else
                DTREK_ERROR(6, sInvalidArgument);
            }
          else
            {
              DTREK_ERROR(5, sMissingOption);
            }
        }
      else if ("-intmask" == sTemp) 
      {
          i++;
          if (i >= argc) 
              DTREK_ERROR(5,sMissingOption);
          sIntMask = argv[i];
      }
      else if ("-peakinfo" == sTemp) 
        {
          poFind->m_bPeakInfo = TRUE;
          poFind->vHeaderToDetectorAreas();

        }
      else if ("-rank" == sTemp) 
        {
          poFind->m_bRank = true;
          bRank = true;
        
          poFind->m_bIncludeBackground = true;
          bPrintTables = true;

          bFindRings = true;
      }
      else if ("-display" == sTemp)
        {
#if (!defined(SSI_PC) && !defined(NO_X_WINDOWS))
          poXprop = new CXprop("DTDISPLAY_WINDOW_ID");
#else
     poFind->m_nDisplay = 1;
#endif
        }
      else if ("-help" == sTemp)
        {
          DTREK_ERROR(0, "");
        }
      else
        {
          DTREK_ERROR(0, "");
        }
    }
#ifdef SSI_PC
    cout << flush;
#endif

  if (!poFind->m_poScan->m_poNonunf->bIsAvailable())
    {
      cout << "ERROR: problem with mask/non-uniformity info!\n" << flush;
#if (!defined(SSI_PC) && !defined(NO_X_WINDOWS))
      if (NULL != poXprop)
        delete poXprop;
#endif

      if (NULL != poFind)
        delete poFind;

      if (NULL != poReflnlist)
        delete poReflnlist;

//      if (NULL != poHeader)
//        delete poHeader;

      nStat = -3;
      return (nStat);
    }

  // Make sure that a -out is specified with -beamcenter
  if ((bFindBeamCenter) && (!sOutHead.length())) {
      printf("ERROR:  Must specifiy -out command when using -beamcenter\n");
      return 1;
  };

  //
  // Done with command line arguments so now go find reflections!
  //

  Crotation     *poRotation;
  Crefln        *poRefln;
  Cdetector     *poDetector = NULL;
  Csource       *poSource   = NULL;
  CScanBitmap    oScanBitmap;

  nStat = 0;

  // Create the source object just once

  poSource = new Csource(*m_poHeader);
  if (!poSource->bIsAvailable())
    {
      cout << "ERROR: dtfind - creating source object!\n";
      return (-1);
    }

  // Create the detector object just once and do not create
  // its non-uniformity, since we have that in the poScan object.

  poDetector = new Cdetector(*m_poHeader, "", TRUE, FALSE);

  if (!poDetector->bIsAvailable())
    {
      cout << "ERROR: dtfind - creating detector object!\n";
      return (-1);
    }

  if (sIntMask.length()) {
      // Integration mask specified.  Open the associated scan bitmap file.
      if (oScanBitmap.nOpenInput(sIntMask)) {
          cout << "ERROR:  Could not open scan bitmap file.\n";
          return (-1);
      };

  };

  if ((bFindBeamCenter) && ((fInitBeamCenter0<0) || (fInitBeamCenter1<0))) {
      poDetector->m_poSpatial->nGetBeamPosition(&fInitBeamCenter0,&fInitBeamCenter1);
  };

  float a3fS0[3];
  float fTempW;
  fTempW = poSource->fGetWavelength();
  poSource->vCalcGetS0(&a3fS0[0]);  // Required for resolution calc below
  float fTemp1, fTemp2, fResoMin, fResoMax;
//mrp  fResoMin = max(a2fResolution[0], a2fResolution[1]);
  fResoMin = std::max(a2fResolution[0], a2fResolution[1]); //mrp
  fResoMax = min(a2fResolution[0], a2fResolution[1]);

  poDetector->nGetResolution(a3fS0, &fTemp1, &fTemp2);
  fTemp1 = fTemp1 * fTempW;
  fTemp2 = fTemp2 * fTempW;
  if (fTemp1 < fResoMin)
    fResoMin = fTemp1;
  if (fTemp2 > fResoMax)
    fResoMax = fTemp2;
  poFind->vSetResolution(fResoMin, fResoMax);

  cout << "\nResolution limits of an image are " << fTemp1
       << " to " << fTemp2 << '\n';
  cout << "Resolution limits of peak search are " << fResoMin
       << " to " << fResoMax << '\n' << endl;

  // Loop through number of -seq options

  
  nSeqCount = 0;
  nTotalImages = 0;
  
  // We might have a subset of the sequence options.
  if ((eFind_2D_method == eMethod) || (eFind_2D_merge_method == eMethod))
    sSeqOptions = "";

  if (eMethod < eMethod2)
    std::swap(eMethod,eMethod2);

  if (!nNumSeqOption) 
    { 
      // No -seq option specified

      a20nSeqInfo[0][0] = poFind->m_poScan->nGetSeqNum();
      a20nSeqInfo[0][1] = poFind->m_poScan->nGetSeqNum();
      a20nSeqInfo[0][2] = 0;
      a20nSeqInfo[0][3] = 1;
      nNumSeqOption++;
    }

 //+2011-07-20 JWP
 // It is probably easiest to modify a20nSeqInfo[][] based on 
 // the expected action of nEveryNthImage here

  if ( (1 == nNumSeqOption) && (1 < abs(nEveryNthImage)) )
    {
      // See help for how this is supposed to work

      int nSeqStart, nSeqEnd;
      int nNum;
      int nTemp;
      int nMaxImages; 
      float fSeqStep;
      nSeqStart = a20nSeqInfo[0][0];
      nSeqEnd   = a20nSeqInfo[0][1];
      nNum = nSeqEnd - nSeqStart + 1; // Total num of images in the scan

      if (0 > nEveryNthImage)
	{
	  fSeqStep = float(nNum ) / float(abs(nEveryNthImage)-1);
	  nMaxImages = min(20, abs(nEveryNthImage));
	}
      else
	{
	  fSeqStep = nEveryNthImage;  // But remember only up to 20 -seq options are used
	  nMaxImages = min(20, nint((float)nNum / fSeqStep));
	}
      if (1.0 > fSeqStep) fSeqStep = 1.0;
      nMaxImages = min(nNum, nMaxImages);

      cout << "INFO: -every option modifies the -seq option ...\n      Seq:";
      int i;
      int nFirstImageInScan;
      
      nNumSeqOption = 0;  // This makes sure the first one is also overwritten and thus set properly
      nFirstImageInScan = poFind->m_poScan->nGetSeqNum(0);
      //cout << "First Image in Scan is: " << nFirstImageInScan << endl;
      int nn;
      for (nn = 0; nn < (nMaxImages-1); nn++)
	{
	  i = nSeqStart + int(fSeqStep * (float)nn);
	  a20nSeqInfo[nNumSeqOption][0] = i;
	  a20nSeqInfo[nNumSeqOption][1] = i;
	  a20nSeqInfo[nNumSeqOption][2] = i - nFirstImageInScan;
	  a20nSeqInfo[nNumSeqOption][3] = 1;
	  nNumSeqOption++;  //JWP: Probably should use sizeof() instead of 20
	  cout << "  " << i;
	}
      // And make sure last image is included, too.
      a20nSeqInfo[nNumSeqOption][0] = nSeqEnd;
      a20nSeqInfo[nNumSeqOption][1] = nSeqEnd;
      a20nSeqInfo[nNumSeqOption][2] = nSeqEnd - nFirstImageInScan;
      a20nSeqInfo[nNumSeqOption][3] = 1;
      if (20 > nNumSeqOption) nNumSeqOption++;
      cout << "  " << nSeqEnd;
      cout << "\nINFO: Total number of images to be searched: " << nNumSeqOption << "\n" << flush;
    }
 //-2011-07-20 JWP

  do
    {
      poFind->m_a2nSeqNum[0]        = a20nSeqInfo[nSeqCount][0];
      poFind->m_a2nSeqNum[1]        = a20nSeqInfo[nSeqCount][1];
      poFind->m_poScan->vSetImgNum(   a20nSeqInfo[nSeqCount][2]);        
      nNumImages                    = a20nSeqInfo[nSeqCount][3];
      nFirstReflnLastSeq            = poFind->m_poReflnlist->nGetNumReflns();
      
      // If the sequence number start and end are the same, we want to
      // force 2D search when the method is 3D.
          
      if (eFind_unknown_method == eMethod)
        //poFind->m_eMethod = eFind_2D_merge_method;
        poFind->m_eMethod = eFind_2D_method;
      else
        poFind->m_eMethod = eMethod;


      if ((eFind_2D_method == poFind->m_eMethod) || (eFind_2D_merge_method == poFind->m_eMethod))
        {
          printf("\ndtfind: 2D method %s used\n",((eFind_2D_merge_method == poFind->m_eMethod)?"with merging":""));

          // Get/read first image to search
          // There is no check that subsequent images are the same size and type.

          if (0 == nNumSeqOption)
            {
              // Read just the image specified on the command line
              
              cout << "...reading image " << sScanFile << "..." << endl << flush;
              poImage = new Cimage(sScanFile);
                          poImage->m_bReadApplyGonioMask = true;
                          poImage->m_bReadApplyEmbeddedMask = true;
              if (!poImage->bIsAvailable())
                  nStat = 1;

              sTemp   = sScanFile + " Template: " + sScanFile + " New!";
            }
          else
            {
              poImage = new Cimage();
                          poImage->m_bReadApplyGonioMask = true;
                          poImage->m_bReadApplyEmbeddedMask = true;

              (void) poFind->m_poScan->nGetImageName(&sTemp);
              cout << "...reading image " << sTemp << "..." << endl << flush;
              nStat = poFind->m_poScan->nGetImage(poImage);
              if (sIntMask.length()) {
                  if (oScanBitmap.nReadBitmapForScan(*poFind->m_poScan)) {
                      printf("ERROR:  Could not load bitmap image.\n");
                      nStat = 1;
                      break;
                  };
              };
              if (oScanBitmap.nMaskOut(*poImage))
                  nStat = 1;


              sTemp =  sTemp + " Template: "
                             + poFind->m_poScan->sGetTemplate() + " New!";
            }

          // Loop through images

#if (!defined(SSI_PC) && !defined(NO_X_WINDOWS))
          if ( (1 >= nNumImages) && (2 > nNumSeqOption) )
            {
              // If searching only a single image, do display update later.

              poXprop2D = NULL;
            }
          else
            {
              poXprop2D = poXprop;
            }
#endif
                  

          nFirstReflnThisSeqOption = poFind->m_poReflnlist->nGetNumReflns();         

          double        fLastRotStep = 0.0;
          double        fThisRotStep = 0.0;
          itr<int>      anOnEdge;
          itr<int>      anImageFirstRef;
          itr<int>      anImageLastRef;
          itr<int>      anImageSort;

          // If we are trying to 'sample' the images, then we select images in a way that
          // obtains the largest amount of reciprocal space before terminating.
          // This involves determining an image sorting order.

          -anImageSort;
          
          if( bSparseFind )
          {
              int   nBest = 0;
              int   nTest = 0;
              int   nDist = 0;
              int   nBestDist = 0;
          
              anImageSort[0] = (poFind->m_a2nSeqNum[0]);
              
              for (i = 1; (i< min(50,nNumImages)); i++) {                  
                  nBest = -1;
                  for (j = 0; (j < nNumImages); j++) {
                      nTest = (j + poFind->m_a2nSeqNum[0]);
                      nDist = 1000000;
                      for (k = 0; (k < anImageSort.size()) && (nDist); k++) 
                          nDist = min(nDist,abs(nTest - anImageSort[k]));
                      nDist = min(nDist,abs(nTest - (poFind->m_a2nSeqNum[1]+1)));

                      if ((nDist) && ((nBest==-1) || (nDist>nBestDist))) {
                          nBest = nTest;
                          nBestDist = nDist;
                      };
                  };
                  anImageSort + nBest;
              };
              // Place the rest of the images in in order.
              for (i = 0; (i< nNumImages); i++) {
                  nTest = (i + poFind->m_a2nSeqNum[0]);
                  if (-1 == anImageSort.find(nTest))
                      anImageSort + nTest;                  
              };

          } else {
              for (i = 0; (i< nNumImages); i++) 
                  anImageSort[i] = (i + poFind->m_a2nSeqNum[0]);
          };


          -anOnEdge;
          -anImageFirstRef;
          -anImageLastRef;
         
          for (i = 1; (i <= nNumImages) && (0 == nStat); i++)
            {
              nFirstReflnThisImage = poFind->m_poReflnlist->nGetNumReflns();
              anImageFirstRef + nFirstReflnThisImage;
              nTotalImages++;


              // Update any display program

              poRotation = new Crotation(poImage->m_oHeader);
              if (!poRotation->bIsAvailable())
                {
                  // Rotation info not in image header, so try to
                  // create rotation from the scan info we have

                  delete poRotation;
                  poRotation = new Crotation(poFind->m_poScan,
                                             poFind->m_poScan->nGetSeqNum());
                }

              fLastRotStep = fThisRotStep;
              fThisRotStep = poRotation->fGetMidValue();
              if (i > 1) {
                  if ((fThisRotStep - fLastRotStep)/poRotation->fGetIncrement()>1.1) {
                      anOnEdge.last() = 1;
                      anOnEdge + 1;
                  } else {
                      anOnEdge + 0;
                  };
              } else {
                  anOnEdge + 1;
              };

              
                //if( i == 1 && bFindRings ) 
                if( bFindRings ) 
                {
                    CIceRingLocator     oLocator = CIceRingLocator(m_poHeader);
                    oLocator.vLocateRings(poFind->m_poScan->nGetSeqNum(),
                                          fFindRingsMin2Theta,              
                                          fFindRingsWidthToExclude, 
                                          nFindRingsImageSizeCompression);
                } 
              
              
              if ((i == 1) && (bFindBeamCenter))
              {
                  printf("\nFinding Beam center using symmetric R-merge method");
                  printf("\n--------------------------------------------------\n");
                  fflush(stdout);

                  Cbeamdata oBeamData(*poImage);
                  int a2nBeamCenter[2];
                  
                  // Run this to get the lowest R-merge beam center.
                  oBeamData.m_oBackground.nBuildPeakPixels(30);
                  oBeamData.m_oBackground.nFatten(2);
                  oBeamData.m_oBackground.nExcludeMinimum(poFind->m_fMinimum);

                  if (!oBeamStat0.nSize()) {
                    a2nBeamCenter[0]= (int) fInitBeamCenter0;
                    a2nBeamCenter[1]= (int) fInitBeamCenter1;
                  } else {
                      a2nBeamCenter[0] = (int) oBeamStat0.fAverage();
                      a2nBeamCenter[1] = (int) oBeamStat1.fAverage();
                  };
                  oBeamData.nFindLowestRmerge(a2nBeamCenter[0],a2nBeamCenter[1]);    
                  oBeamStat0.vAdd((double) a2nBeamCenter[0]);
                  oBeamStat1.vAdd((double) a2nBeamCenter[1]);
                  if (bFindBeamCenterDiagnostics) {
                      sTemp = "beamcenter_peaks";
                      sTemp += oBeamStat0.nSize();
                      sTemp += ".img";
                      oBeamData.m_oBackground.nPrintPeaks(sTemp);
                      sTemp = "beamcenter_background";
                      sTemp += oBeamStat0.nSize();
                      sTemp += ".img";
                      oBeamData.m_oBackground.nPrintBackground(sTemp);
                  };
              } 
              
              {
                  
                  
#if (!defined(SSI_PC) && !defined(NO_X_WINDOWS))
                  poFind->m_poXprop = poXprop2D;
                  poFind->m_psXpropUpdate = &sTemp;
#endif

                  poFind->nMarkGonioNumbers();
                  nStat      = poFind->nFind2D(*poImage, poRotation->fGetMidValue(),
                                               poDetector, poSource,
                                               (eFind_2D_merge_method == poFind->m_eMethod));

                  // We might need to add goniometer information.
                  
                  if (poFind->nAddGonioNumbers(*poImage)) {
                      break;
                      nStat = 1;
                  };
                  
                  // OK, found spots, now find centroid of each one
                  DTREK_WORD        wCtrl = bRank ? DTREK_DTFIND_FC2D_GET_PEAK_AREA_INFO : 0U; 

                  for (j = nFirstReflnThisImage; j < poFind->m_poReflnlist->nGetNumReflns(); j++)
                  {
                      poRefln = poFind->m_poReflnlist->poGetRefln(j);
                      
                      
                      nStat = poFind->nFindCentroid2D(*poImage, poRefln, poDetector, wCtrl);
                      
                      if (0 != nStat)
                      {
                          // Reset local status and mark refln for deletion
                          
                          nStat = 0;
                          poRefln->vSetIntensity(-999.0);
                      }
                  }    
                  
              }

              delete poRotation;

              anImageLastRef + (poFind->m_poReflnlist->nGetNumReflns() -1);              
              if ((nMinReflections>0) && (poFind->m_poReflnlist->nGetNumReflns()-nFirstReflnLastSeq>=nMinReflections) && (nTotalImages>=nMinImages))
                  break;

              if (i < nNumImages)
                {
                  // Get/read next image
                  // Set seqnum to next image in scan

                  poFind->m_poScan->vSetImgNum(anImageSort[i]+ (a20nSeqInfo[nSeqCount][2] -a20nSeqInfo[nSeqCount][0]));
                  (void) poFind->m_poScan->nGetImageName(&sTemp);
                  cout << "...reading image " << sTemp << "..." << endl << flush;
                  nStat = poFind->m_poScan->nGetImage(poImage);

                  if (sIntMask.length()) {
                      if (oScanBitmap.nReadBitmapForScan(*poFind->m_poScan)) {
                          printf("ERROR:  Could not load bitmap image.\n");
                          nStat = 1;
                          break;
                      };
                  };
                  if (!nStat) 
                      nStat= oScanBitmap.nMaskOut(*poImage);


                  sTemp =  sTemp + " Template: "
                                 + poFind->m_poScan->sGetTemplate()
                                 + " Seq: "
                                 + Cstring(poFind->m_poScan->nGetSeqNum());
                }
      
              /*
              // Used to view masked out portions of image.
              {
                   static int nCount = 0;
                   Cstring sName;
                   sName = "thad";
                   sName += nCount++;
                   sName += ".img";
                   poImage->nWrite(sName);
               };
               */
                       

            } // end of i loop
          delete poImage;

          poFind->m_nTotalImages = nTotalImages;

          // Consolidate all -seq options.

          if (!nStat) {
              for (j = 0,k = 0; (j < min(i,nNumImages)); j++) {
                  if ((j+1>=min(i,nNumImages) || (anImageSort[j+1]!=anImageSort[j]+1))) {
                      sSeqOptions += " -seq ";
                      sSeqOptions += anImageSort[k];
                      sSeqOptions += " ";
                      sSeqOptions += anImageSort[j];
                      k = j+1;
                  };
              };
          };


          // If all of the images are on edge, then reset.
          if (0) {
              if (-1 != anOnEdge.find(i=0)) {
                  int nFI_nMark = poFind->m_poReflnlist->nGetFieldIndex(poFind->ms_sOnFirstLastImage);
                  for (j = 0; j < anImageFirstRef.size(); j++) {
                      for (k = anImageFirstRef[j]; k<= anImageLastRef[j]; k++) {
                          (*poFind->m_poReflnlist)[k].vSetField(nFI_nMark,(int) anOnEdge[j]);
                      };
                  };
              };
          };

          // Delete any reflection which has a bad intensity
          // (i.e. Centroid determination failed)

          for (j = (poFind->m_poReflnlist->nGetNumReflns()-1);
               j >= nFirstReflnThisSeqOption; j--)
            {
              poRefln = poFind->m_poReflnlist->poGetRefln(j);
              if (0.0 > poRefln->fGetIntensity())
                {
                  poFind->m_poReflnlist->nDelete(j);
                }
            }         
        }      
    } while (nSeqCount++,nSeqCount < nNumSeqOption);


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Remove reflections, if any, based on the ice rings that were found, if any.
    DTREK_WORD      wCtrl = bRank ? DTREK_REFLNLIST_DL_RANK | DTREK_REFLNLIST_DL_FROM_FIND : DTREK_REFLNLIST_DL_FROM_FIND; 

    if( bFindRings )
        poFind->m_poReflnlist->vRemoveReflnsInResolutionRingAreas(m_poHeader, (double)fFindRingsWidthToExclude, wCtrl);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if (bPrintTables)
    {
        poFind->m_poReflnlist->nResoTable(bResolutionProvided?a2fResolution[0]:-1.0,
                                          bResolutionProvided?a2fResolution[1]:-1.0,
                                          m_nNumberOfResoBins,
                                          10, // 10 I/Sigma bins
                                          eTableVsResolution | eTableVsIoverSig,
                                          -1,
                                          m_poHeader);        
    }


    if( bRank )
    {
        poFind->m_poReflnlist->vRank(m_poHeader,
                                     bResolutionProvided?a2fResolution[0]:-1.0,
                                     bResolutionProvided?a2fResolution[1]:-1.0,
                                     m_nNumberOfResoBins,
                                     m_nMinReflnsForRankStats);        
    }


    poFind->vSetOptionString(sFindOptions);
    m_poHeader->nReplaceValue((Cstring) D_K_DtfindSeqOptions, sSeqOptions);
    poFind->nUpdateHeader(m_poHeader);
    poFind->m_poReflnlist->nPutHeaderInfo(*m_poHeader,(Cstring) D_K_HeaderDtfindInfo,true);

    
    if (eFind_2D_merge_method == poFind->m_eMethod) {
        if (eMethod2 == eFind_2D_method) {
            // User requested both -2D and -2Dmerge output.
            nTemp = poFind->nWrite(*poFind->m_poReflnlist,sReflnFile2);
            cout << "\ndtfind: There were " << poFind->m_poReflnlist->nGetNumReflns()
                << " spots found." << endl << flush;
           
            if (0 == nTemp)
            {
                cout << "dtfind: Spots written to " << sReflnFile2 << endl << flush;
            };
        };


                        

        printf("Merging 2D reflections...\n");
        Creflnlist* poNewList           = NULL;

        poNewList = new Creflnlist;
        poFind->nMerge2DReflections(poNewList);
        delete poFind->m_poReflnlist;
        poReflnlist = poNewList;
        poFind->m_poReflnlist = poNewList;
        printf("\n");
    };
    
  // Write out the reflnlist with found spots, even if there were errors finding
  //   spots

  if (1 > poFind->m_poReflnlist->nGetNumReflns())
    {
      cout << "There were no spots found.\n" << flush;
      nStat = 1;
    }
  else
    {
      cout << "\ndtfind: There were " << poFind->m_poReflnlist->nGetNumReflns()
           << " spots found." << endl << flush;
      nTemp = poFind->nWrite(*poFind->m_poReflnlist,sReflnFile);

      if (0 == nTemp)
        {
          cout << "dtfind: Spots written to " << sReflnFile << endl << flush;
#if (!defined(SSI_PC) && !defined(NO_X_WINDOWS))
          if (NULL != poXprop)
            {
              // Be sure sTemp is not changed since above

              if (sReflnFile == sFileGetBasename(sReflnFile))
                  sReflnFile = sGetCWD() + sReflnFile;

              poXprop->hSetProperty("DTDISPLAY_IMAGE_UPDATE",
                  sTemp + " Reflnlist: "
                  + sReflnFile + " Header: " + sScanFile);

            }
#else
#ifndef NO_X_WINDOWS
  if(poFind->m_nDisplay)
  {
      CCrclHelper::GetInstance()->vSendFindUpdateDisplay( "All", sReflnFile );
      /*Cstring Result;
      CCoreTclParam TclCom("FindUpdateDisplay");

      TclCom.SetParamString("RefList", sReflnFile);
      TclCom.SetParamString("Server", "Dtrek");
      TclCom.SetParamString("Image", "All");
      // Tell the calling process what to display from the ref list.
      TclCom.SetParamInt("ShowObserved", 1);
      TclCom.SetParamInt("ShowCalculated", 0);
      TclCom.SetParamInt("ShowDifference", 0);

      // Send the Tcl command to be processed.
      gSocket.SendString(TclCom.ConstructCommand());
      // Block and wait for the display to be updated so that we
      // don't get ahead of the display.
      gSocket.ReceiveString(Result);*/
  }
#endif
#endif
        
      }
      else
        {
          cout << "ERROR: dtfind - writing file: " << sReflnFile << '\n' << flush;
          if (0 == nStat)
            {
              nStat = nTemp;  // Don't overwrite prev nStat if bad
            }
        }
    }

  if ("" != sOutHead)
    {

      if (bFindBeamCenter) {
          // Reject outliers.
          oBeamStat0.nMarkStat(FALSE,oBeamStat0.S_ABOVEBELOW,3.0);
          oBeamStat1.nMarkStat(FALSE,oBeamStat0.S_ABOVEBELOW,3.0);
          poDetector->m_poSpatial->nSetBeamPosition(oBeamStat0.fAverage(),oBeamStat1.fAverage());
          poDetector->m_poSpatial->nUpdateHeader(m_poHeader);
      };
     
      if (sOutHead == sReflnFile)
        {
          cout << "ERROR: dtfind - reflection file and header file are same!\n"
               << "                Output header file NOT written!\n" << flush;

        }
      else
        {
          nTemp = m_poHeader->nWrite(sOutHead);
          if (0 == nTemp)
	    {
	      cout << "dtfind - Wrote header file " << sOutHead << endl << flush;
             if (19 < poFind->m_poReflnlist->nGetNumReflns())
              cout << "\nINFO - The suggested next step is to run dtindex with"
                   << "\n       " <<  sOutHead << " as the input .head file and" 
                   << "\n       " << sReflnFile << " as the input .ref file:\n" 
                   << "\n       dtindex " << sOutHead << " " << sReflnFile << " ...\n\n" << flush;
	     else
              cout << "\nINFO - Not very many spots were found. The suggested"
                   << "\n       next step is to re-run dtfind with different"
                   << "\n       arguments to help find more spots.\n\n" << flush;
            }

          else
            {
              cout << "dtfind - ERROR writing header file " << sOutHead << endl << flush;
              if (0 == nStat)
                nStat = nTemp;  // Don't overwrite prev nStat if bad                
            }
        }
    }
#if (!defined(SSI_PC) && !defined(NO_X_WINDOWS))
  if (NULL != poXprop)
    delete poXprop;
#endif
  if (NULL != poFind)
    delete poFind;

  if (NULL != poReflnlist)
    delete poReflnlist;

  if (NULL != poDetector)
    delete poDetector;

  if (NULL != poSource)
    delete poSource;

  //if (NULL != poHeader)
  //  delete poHeader;

  // Attempt to delete temporary reflnlist files 
  if (bPurgeTmpFiles)
    {
#if !defined(SSI_PC) && !defined(WIN32)
      (void) nDoSystemCommand((Cstring)"sleep 2; rm -f *dtfind*tmp*.ref &");
#endif
    }

  return (nStat);
}

void CDTMainFind::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
  cout << "\ndtfind - Usage:\n"
       "dtfind sHeaderfile [ options ...]\n\n"
       "Command line options: Description:\n\n"
       "sHeaderfile           Name of a file with a header that has the scan\n"
       "                      definition. The search method will gbe 3D \n"
       "                      unless overridden by the -2D option.\n\n"
       " -reso fResoMin fResoMax\n"
       "                      Specifies the resolution range to search for\n"
       "                      reflections.  The default range of the entire image\n"
       "                      is determined by the detector position and source\n"
       "                      wavelength.\n\n"
       " -reso edge           Set upper resolution to be equal to the highest\n"
       "                      *minimum* resolution found on any of the four\n"
       "                      edges (sides) of the image.\n\n"
       " -reso corner         Set upper resolution to be equal to the highest\n"
       "                      *maximum* resolution found on any of the four\n"
       "                      edges (sides) of the image.\n\n"
       " -resobins nBins [nMinReflns]\n"
       "                      nBins - Number of resolution bins to be used for\n"
       "                      statistics. nMinReflns - minimum number of reflns \n"
       "                      in a bin to be used for statistics.\n\n"
       " -sigma  fSigma       Number of standard deviations (default 5.0) above\n" 
       "                      background for pixels to be considered for peaks.\n"  
       "                      This is tested AFTER any non-uniformity correction\n"
       "                      is applied.\n\n"
       " -min    fMin         Minimum pixel value to be considered for peaks.\n"
       "                      This is tested AFTER any non-uniformity correc-\n"
       "                      tion is applied.  The default is 50.\n\n"
       " -filter nPeakFilter  In a 3x3 box, nPeakFilter pixels must be above\n"
       "                      the threshold which is the higher of the value\n"
       "                      specified by the -fSigma or -fMinimum option.\n\n"
       " -peakmin fPeakMin    Reject a found peak if the difference between the peak\n"
       "                      maximum pixel value and the peak mean background pixel\n"
       "                      value is less than fPeakMin. This is tested _after_ any\n"
       "                      non-uniformity correction is applied. The default is 0.\n\n"
       " -seq    nStart nEnd  Specifies the starting and ending image sequence\n"
       "                      values in a scan to search for spots, but the images\n"
       "                      searched may also affected by the -every option.  Up\n"
       "                      to 20 -seq options may be used together.  The default\n"
       "                      is the first image of the scan and 9 more\n"
       "                      images. \n\n"
       " -every nEveryImage   Note: Only applies if one -seq option is given.  If\n"
       "                      multiple -seq options are given, then -every is ignored.\n\n"
       "                      Specifies whether every image in a -seq option is\n"
       "                      searched (nEveryImage = 1) or if only some images\n"
       "                      are searched.  If nEveryImage > 1, then only every\n"
       "                      nth nEveryImage is search (e.g. if nEveryImage = 10,\n"
       "                      then every 10th image is searched).  If nEveryImage < 1\n"
       "                      then about |nEveryImage| number of images in total are\n"
       "                      searched.  The images chosen are usually the first and\n"
       "                      last image as specified by the -seq option and some\n"
       "                      relatively evenly spaced images in-between these two\n"
       "                      images.\n"
       "                      NOTE: At most only 20 images will be searched if the\n"
       "                      -every option is in effect. Default is 1.\n\n"
       " -nosaturated         Do not include saturated reflections.  By default.\n"
       "                      saturated reflections are included.\n\n"
       " -circle nCen0 nCen1 nMinrad nMaxrad\n"
       "                      Only pixels within a circle of nMaxrad and\n"
       "                      outside a circle of nMinrad pixels will be\n"
       "                      considered for peaks.  The circles are center\n"
       "                      at nCen0, nCen1.  The default settings are\n"
       "                      the image center and radii to search the entire\n"
       "                      image.\n\n"
       " -rect   nMinPx0 nMinPx1 nMaxPx0 nMaxPx1\n"
       "                      Only pixels within a rectangle with origin at\n"
       "                      nMinPx0, nMinPx1 and opposite corner nMaxPx0, nMaxPx1\n"
       "                      will be considered for peaks.  The defaults\n"
       "                      settings are the entire image less a border\n"
       "                      of 1%.\n\n"
       " -window nWin0 nWin1  Size in pixels for a window of pixels around\n"
       "                      a peak.  For a 2D search, the default values\n"
       "                      of 0, 0 make dtfind dynamically a suitable\n"
       "                      size.  For a 3D search, the default is 21,21.\n\n"
//       " -nonunf sNonunfFile  Specifies an alternate file with mask/non-uniformity\n"
//       "                      information to use.  The default info comes from\n"
//       "                      the -scan file or the first -image file.  This\n"
//       "                      option overrides the information in these files.\n\n"
       " -ref    sReflectionFileName\n"
       "                      Specifies an alternate output reflection file\n"
       "                      instead of the default dtfind.ref. The file format\n"
       "                      is a standard d*TREK reflection file format.\n\n"
       " -num    nNumImage    Specifies an alternate number of images to search\n"
       "                      in a scan.\n\n"
       " -dump   nthRefln     Specifies when to dump or write reflection\n"
       "                      shoeboxes to disk.  The non-uniformity cor-\n"
       "                      rected data of every nthRefln will be written to\n"
       "                      disk as an image file with name \n"
       "                      dtfindref????.img, if nthRefln > 0.  This\n"
       "                      applies only if a 3D search method is used.\n\n"
       " -2D                  Uses a 2D instead of a 3D search method when -seq\n"
       "                      is specified.  This is the default.\n\n"
       " -out    sOutHead     Add the dtfind command line options to the first\n"
       "                      header read and write out a new header file with\n"
       "                      the dtfind command line options to the file sOutHead.\n"
       "                      Default: Do not write new header file.\n\n"
       " -minrefs nMinRefs    Minimum number of 2D spots to find before terminating\n"
       "                      search.  If the argument nMinRefs is <0 then dtfind\n"
       "                      will not restrict itself to use contiguous images\n"
       "                      to meet it's quota of spots.  Instead, it will use\n"
       "                      images that are spaced widely apart until it has\n"
       "                      at least ABS(nMinRefs) spots.\n"
       "                      Default: Search all images specified in -seq\n\n"
       " -minimages nImages   Minimum number of images to search for spots on.  This\n"
       "                      should only be used if you fear that the -minrefs option\n"
       "                      would cause dtfind to not use enough images.\n\n"
       " -display             Tell a dtdisplay process to display images and\n"
       "                      reflection lists as dtfind proceeds.\n\n"
       " -intmask sFile       Specify an integration mask.\n\n"
       " -help                Print this help text.\n\n"
       " -beamcenter [diagnostic] [nInit0 nInit1]\n"
       "                      Finds beam center for images in the search range using\n"
       "                      initial beamcenter at [nInit0 nInit1] (otherwise\n"
       "                      defaults to beamcenter found in header).  Output is\n"
       "                      placed in header file (see -out).  If keyword\n"
       "                      'diagnostic' is specified, debug info is written.\n\n"
       " -ringdetect [f2ThetaToExclude] [fMin2ThetaToSearch]\n"
       "                      Finds powder rings (for example, water rings) on the\n"
       "                      first image. If a ring is found at 2Theta = X, then the\n"
       "                      excluded 2theta range is from (X - f2ThetaToExclude/2) to\n"
       "                      X + f2ThetaToExclude/2).  Default: 2.0 20.0\n\n"
       " -peakinfo            Use the peak ellipsoid information in sHeaderfile.\n"
       "                      Report strong peak information in the log file.\n"
       "                      Default: no peak ellipsoid info is used, no strong peak\n"
       "                      info is reported.\n\n"
       " -nopurge             Do not delete any temporary *dtfind*tmp*.ref files.\n"
       "                      On Windows systems, files are not purged regardless.\n\n"
       " -tables              Output statistics of found peaks in terms of resolution\n"
       "                      and I/Sigma.\n\n"
       "Information about scan angles and mask/non-uniformity of response is\n"
       "determined from the headers of the input image files.  If the\n"
       "headers are incorrect, they may be edited with dtheaderedit.\n\n"
       "Examples:\n\n"
       "   dtfind dtprocess.head -seq 1 1\n"
       "   dtfind dtprocess.head -seq 1 1 -seq 90 90\n"
    "   dtfind dtprocess.head -seq 1 1 -sigma 0 -min 1000 -window 30 30\n\n" << flush;
}


