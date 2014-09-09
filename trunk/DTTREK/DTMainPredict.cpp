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
// DTMainPredict.cpp       Initial author: RB       13-Oct-2004
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
 
//+Description
//
//    Example:  % dtpredict header_file [ options ...]
//
//    Command line options:
//                header_file sFilename
//                -seq                  nStart nEnd
//                -ref                  sReflectionFileName
//                -rot                  fRotStart fRotEnd
//                -reso                 Resolution_min Resolution_max
//                -list                 List predict info before predicting
//                -nonunf               Do not predict reflections in bad nonunf areas
//                -display              Tell dtdisplay to display the refln list
//                -lscale               Scale factor to apply to cell lengths
//                -mosaicity            Override crystal mosaicity in header via
//                -mergeheaderscan, -differentscan
//                -mergeheaderdet, -differentdet
//                -out sOutheaderfilename      Output new header
//                -overlap fO1 fO2 fO3
//                -help
//

#include "Dtrek.h"
#include "dtreksys.h"
#include "Cpredict.h"
#include "Cscan.h"
#include "Crotation.h"

#ifndef NO_X_WINDOWS
    #ifndef SSI_PC
      #include "CXprop.h"
    #endif
#endif

#include "DTMainPredict.h"
#include "CExperimentalStrategy.h"

#include "CCommandLine.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

const char*     c_pcPredictCheckOverlap            = "-overlap";

int CDTMainPredict::nExecute(unsigned int argc, char* argv[])
{
  int            nStat;
  Cstring        sReflnFile   = sDtrekGetPrefix() + "dtpredict.ref";
  Cstring        sOutFile     = sDtrekGetPrefix() + "dtpredict.head";

  Creflnlist    *poReflnlist  = NULL;
  Cpredict      *poPredict    = NULL;
  Cscan         *poScan       = NULL;
#if  ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  CXprop        *poXprop      = NULL;
#endif
  int          i, j;  // Loop counters
  int          nTemp;
  int          nTemp2;
  float        fTemp,fTemp2;
  Cstring      sTemp;
  bool         bDoList;
  bool         bWriteNewHeader   = FALSE;
  bool         bMergeHeadersScan = FALSE;
  bool         bMergeHeadersDet  = FALSE;

  bool         bDoOverlap = FALSE;
  float        fScaleLength;
  float        a2fResolution[2];

  itr<double>  afRotStart;      // Rotation ranges 
  itr<double>  afRotEnd;        // Rotation ranges 
  itr<int>     anSeqNumStart;   // Sequence ranges
  itr<int>     anSeqNumEnd;     // Sequence ranges.

  int          a2nSeq[2];
  float        fMaxLorentz;

  fMaxLorentz      = 50.0;
  if ("" != sGetEnv("DTREK_MAXLORENTZ"))
    fMaxLorentz = (float)atof(sGetEnv("DTREK_MAXLORENTZ").string());
  fScaleLength     = 1.0f;        // Default length scale factor
  a2fResolution[0] = 999999.99f;  // Default minimum resolution
  a2fResolution[1] = 0.5f;        // Default maximum resolution
  -afRotStart;
  -afRotEnd;
  -anSeqNumStart;
  -anSeqNumEnd;


  bDoList          = FALSE;

  Cstring        sMissingOption
                 = "ERROR: dtpredict - missing option argument!";
  Cstring        sInvalidArgument
                 = "ERROR: dtpredict - invalid argument!";

  vDtrekSetModuleName("dtpredict");
  vPrintCopyrightInfo();
  //cout << "\ndtpredict:  Copyright (c) 2006 Rigaku\n";
  //cout << D_K_DTREKVersion << endl << flush;

#ifdef SSI_PC
  Cstring sCCAppVersion = (const char*) argv[argc-1];
  if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
	cout << sCCAppVersion;
	argc -= 2;
	cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
	argv[argc] = NULL;
  }
#endif

  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;
  
  vSetCommandArguments(argc, argv);

  // Get header file to start with

  argc--;
  argv++;

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

  if( !bCreateHeader(sTemp) )
    {
      cout << "ERROR in dtpredict: Header not available!\n";
      nStat    = 1;
    }
  else
    {
      nStat  = 0;
      poScan = new Cscan(*m_poHeader);
      if (!poScan->bIsAvailable())
	{
	  cout << "ERROR in dtpredict: Scan not available!\n";
	  nStat = 1;
	}
    }
  if (0 != nStat)
    {
      DTREK_ERROR(3, "");
      return 0;
    }


  poPredict  = new Cpredict(*m_poHeader, NULL);
  if (!poPredict->bIsAvailable()) {
      DTREK_ERROR(1, "ERROR in dtpredict: Something wrong with input image header.\n"
	     "                    Perhaps the crystal unit cell info is missing.\n");
      return 0;
  };



  // Parse any additional command line options

  for (i = 1; i < argc; i++)
    {
      sTemp = (const char*) argv[i];

      if ("-ref" == sTemp)
	{
	  i++;
	  if (i < argc)
	    {
	      sReflnFile = (const char *)argv[i];
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
	      bWriteNewHeader = TRUE;
	      sOutFile = (const char *)argv[i];
	    }
	  else
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
      
      poPredict->vSetResolution(a2fResolution[0], a2fResolution[1]);
	}
    else if ("-overlap" == sTemp)
    //	{
    //	  bDoOverlap = TRUE;
    //	}
    {
          // skip this option, because it will be parsed later
          while( i + 1 < argc && !bIsDTCommandLineOptionString(argv[i+1]) )
              i++;
    }
    else if ("-rot" == sTemp)
	{
	  i++;
	  for (j=0; (j < 2) && (i < argc); j++, i++)
	    {
	      nTemp = sscanf(argv[i], "%f", &fTemp);
          if (nTemp == 1) {
              if (j == 0)
                  -afRotStart + fTemp;
              else
                  -afRotEnd + fTemp;
          } else
		DTREK_ERROR(6, sInvalidArgument);
	    }
	  i--;
	  if (j != 2)
	    {
	      DTREK_ERROR(5, sMissingOption);
	    }
	}
      else if ("-seq" == sTemp)
	{
	  i++;
	  for (j=0; (j < 2) && (i < argc); j++, i++)
	    {
	      nTemp = sscanf(argv[i], "%d", &nTemp2);
	      if (nTemp == 1)
		a2nSeq[j] = nTemp2;
	      else
		DTREK_ERROR(6, sInvalidArgument);
	    }
	  i--;
	  if (j != 2)
	    {
	      DTREK_ERROR(5, sMissingOption);
	    }
      anSeqNumStart + a2nSeq[0];
      anSeqNumEnd + a2nSeq[1];

	}
      else if ("-list" == sTemp)
	{
	  bDoList = TRUE;
	}
      else if (   ("-mergeheadersscan" == sTemp)
	       || ("-differentscan" == sTemp) )
	{
	  bWriteNewHeader = TRUE;
	  bMergeHeadersScan = TRUE;
	}
      else if (   ("-mergeheadersdet" == sTemp)
	       || ("-differentdet" == sTemp) )
	{
	  bWriteNewHeader = TRUE;
	  bMergeHeadersDet = TRUE;
	}
      else if ("-nonunf" == sTemp)
	{
          poPredict->vSetNonunfFlag(1);
	}
      else if ("-fastpredict" == sTemp)
	{
          if ((i + 2 < argc) && (1 == sscanf(argv[i+1],"%f",&fTemp)) && (1 == sscanf(argv[i+2],"%f",&fTemp2))) {
            poPredict->vSetFastPredict(fTemp,fTemp2);
            i += 2;
          } else
              DTREK_ERROR(6,sInvalidArgument);
      }
      else if ("-display" == sTemp)
      {

#if  ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
	  poXprop = new CXprop("DTDISPLAY_WINDOW_ID");
#else
     poPredict->m_nDisplay = 1;
#endif
	}
      else if ("-lscale" == sTemp)
	{
	  i++;
	  if (i < argc)
	    {
	      nTemp = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nTemp)
		{
		  fScaleLength = fTemp;
		  poPredict->vSetScaleLength(fScaleLength);
		}
	      else
		{
		  DTREK_ERROR(6, sInvalidArgument);
		}
	    }
	  else
	    {
	      DTREK_ERROR(5, sMissingOption);
	    }
	}
      else if ("-maxlorentz" == sTemp)
	{
	  i++;
	  if (i < argc)
            {
	      nStat = sscanf(argv[i], "%f", & fTemp);
	      if (1 == nStat) 
		{
		  fMaxLorentz = fTemp;
		  nStat = 0;
		} 
	      else
		  DTREK_ERROR(6, sInvalidArgument);
            }
	  else
            {
	      DTREK_ERROR(5, sMissingOption);
            }
	}
      else if ("-mosaicity" == sTemp)
	{
	  i++;
	  if (i < argc)
	    {
	      nTemp = sscanf(argv[i], "%f", &fTemp);
	      if (1 == nTemp)
		{		  
		  bWriteNewHeader = TRUE;
		  poPredict->vSetCrysMosaicity(fTemp);
		}
	      else
		{
		  DTREK_ERROR(6, sInvalidArgument);
		}
	    }
	  else
	    {
	      DTREK_ERROR(5, sMissingOption);
	    }
	}
      else if ("-twin" == sTemp)
      {
          int nx,ny;
          for (nx=0;nx<poPredict->m_nNumCrystals;nx++)
              poPredict->m_pbActive[nx] = FALSE;
          ny = 0;
          while (1) {
              i++;
              if ((i>=argc) || (1!=sscanf(argv[i],"%d",&nx)))
                  break;
              if ((nx<=poPredict->m_nNumCrystals) && (nx>=1))
                  poPredict->m_pbActive[nx-1] = TRUE;
              else 
                  DTREK_ERROR(5,sTemp = "Component does not exists in header.");
              ny++;
          };
          if (ny == 0)
              DTREK_ERROR(5,sTemp = "At least one valid component must be specified with -twin");
          i--;
      }
      else
      {
	    DTREK_ERROR(0, "");
      }
    }

    /////////////////////////////////////////////////////
    // Parse overlap option
    Cstring     sError("");
    bDoOverlap = bParseReflnOverlapCheckParams(sError);
    /////////////////////////////////////////////////////

    // Done with command line arguments so now go predict reflections!

    nStat = 0;    

    // Create default reflection list:
    
    poReflnlist = new Creflnlist();
    
    // Check that a detector object does exist:
    
    if (poPredict->m_nNumDetectors > 0) 
        nStat = 0;  // Yes, we can go on.
    
    // All is still well, predict the reflections
    
    poPredict->vSetMaxLorentz(fMaxLorentz);
    
    if (bDoList) 
        poPredict->nList();


    sTemp = poScan->sGetImageName();
    sTemp.upcase();
    if (sFileGetBasename(sTemp).contains("SCREEN") && (!afRotStart.size())) {
        itr<int> anNewSeqNumStart;
        itr<int> anNewSeqNumEnd;
        int nSeq,nSeqSet;

        // Print message regarding screen images.
        printf("Image template contains 'screen'.\nSequence numbers used will be:  ");
        for (nSeqSet = 0; (!nStat) && (nSeqSet < anSeqNumStart.size()); nSeqSet++) {
            for (nSeq = anSeqNumStart[nSeqSet]; nSeq <= anSeqNumEnd[nSeqSet]; nSeq++) {
                anNewSeqNumStart + nSeq;
                anNewSeqNumEnd + nSeq;
                printf("-seq %d %d ",nSeq,nSeq);
            };
        };
        printf("\n");
        anSeqNumStart = anNewSeqNumStart;
        anSeqNumEnd = anNewSeqNumEnd;
    };


    // Do all predictions that involve afRotStart and afRotEnd.
    for (i = 0; (i < afRotStart.size()) && (!nStat); i++) {
        nStat = poPredict->m_poRotation->nSetRotRange(afRotStart[i],
            afRotEnd[i]);
        nStat += poPredict->nPredict(NULL, poReflnlist);
    };
      
    Cimage_header oMergeImgHeader;

    for (i = 0; (i < anSeqNumStart.size()) && (!nStat); i++)
    {
      Cgoniometer* poGoniometer;
      
      poScan->vSetSeqNum(anSeqNumStart[i]);
      Cimage_header oHeader(poScan->sGetImageName());
      Cimage_header oHeaderToUse = (oHeader.bIsAvailable())?oHeader:(*m_poHeader);
      if (!oHeader.bIsAvailable()) {
          printf("WARNING:  Could not load header for '%s'\n",poScan->sGetImageName().string());
          printf("WARNING:  Prediction might be incorrect.\n");
      };
      
      if ( (0 == i) && (bMergeHeadersScan || bMergeHeadersDet))
	{
	  // If this is the first sequence statement AND
	  // we are going to need to merge headers below, then save
	  // a copy of this header for future use.

	  oMergeImgHeader = oHeaderToUse;
	}

      poGoniometer = new Cgoniometer(oHeaderToUse,Ccrystal::ms_sCrystalPrefix);
      if (poGoniometer->nReconcileHeaderWithImageGoniometer(poPredict->m_poCrysGonio)) {
          nStat = 1;
          break;
      };
      delete poPredict->m_poCrysGonio;
      poPredict->m_poCrysGonio = poGoniometer;
      
      

      Cscan oScan(oHeaderToUse);
      Crotation oRotation(oHeaderToUse);
      double fRotationStart;
      double fRotationEnd;

      if (!poPredict->m_poCrysGonio->bIsAvailable() || ! oScan.bIsAvailable() || ! oRotation.bIsAvailable()) {
          printf("ERROR:  Could not initialize goniometer\n");
          nStat = 1;
          break;
      };
            
      if (anSeqNumStart[i] == anSeqNumEnd[i]) {
          // If we only have a single image, trust the ROTATION keyword of the image
          fRotationStart = oRotation.fGetRotStart();
          fRotationEnd = oRotation.fGetRotEnd();
      } else {
          // Otherwise, we need to look at the SCAN_ROTATION keyword.
          fRotationStart = oScan.fCalcGetRotStart(anSeqNumStart[i]);
          fRotationEnd = oScan.fCalcGetRotEnd(anSeqNumEnd[i]);
      };

      
      nStat = poPredict->m_poRotation->nSetRotRange(
      fRotationStart,
      fRotationEnd);
      nStat += poPredict->nPredict(NULL, poReflnlist);

      
  }

  // Write out the reflnlist with found spots, even if there was error
  //   predicting spots
  
  if (1 > poReflnlist->nGetNumReflns())
  {
      cout << "There were no reflections predicted.\n" << flush;
  } 
  else if (!nStat)
  {
      cout << "\ndtpredict: There were " << poReflnlist->nGetNumReflns()
          << " reflections predicted." << endl << flush;
      
      if( bDoOverlap )
      {
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Figure out the spot size
        bool                bDefaultSpotSize = false;
        bool                bCoefficientAppliedToSpotSize = false;
        const double        cfSpotSizeCoefficient = 1.2;  // 1.2 is just a coefficient to provide some 
                                                          // buffer between spots, so that spots don't 
                                                          // "touch" each other.

        if( m_stOverlapCheckParams.m_dDistTest < 0.0 )  // command line says: try and get it from the header
        {
             SPOT_SIZE_INFO              stSpotSizeInfo;   
             CExperimentalStrategy       oExpStrategy(*m_poHeader);

             if( !oExpStrategy.bGetSpotSizeFromHeader(stSpotSizeInfo) )
             {
                printf("WARNING:  Could not find ellipsoid shape information for overlap checking.\n");
                m_stOverlapCheckParams.m_dDistTest = 20; // default - 20 pixels
                bDefaultSpotSize = true;
             }
             else 
             {
                m_stOverlapCheckParams.m_dDistTest = stSpotSizeInfo.fAverageSpotSize * cfSpotSizeCoefficient;
                bCoefficientAppliedToSpotSize = true;
             }
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        nTemp = poReflnlist->nOverlapCheck(m_stOverlapCheckParams.m_dDistTest, 
                                           poPredict->m_poRotation->fGetIncrement(),
                                           poPredict->m_poRotation->fGetRotStart(),
                                           true); // true means: do display overlapped reflections
        
        printf("\nPredicting reflection overlap...\n");
        
        printf("\nTesting coordinate differences between neighboring reflections:\n");
        
        if( bDefaultSpotSize )
            printf("%d pix (default)\n", (int)m_stOverlapCheckParams.m_dDistTest);
        else if( bCoefficientAppliedToSpotSize )
            printf("%d pix = %.1f * \"average spot size\"\n", (int)m_stOverlapCheckParams.m_dDistTest, cfSpotSizeCoefficient); 
        else
            printf("%d pix (from command line)\n", (int)m_stOverlapCheckParams.m_dDistTest);
        printf("\n");
          //// Flag for overlapped reflns before writing out
          //
          //int   nResultField;
          //int   nFirstField;
          //int   nNextField;
          //int   nx;
          //C3DdataDetectorProfile oAreaProfiles;
          //C3DdataDetectorArea oArea;
          //
          //if (oAreaProfiles.nInitValues(*m_poHeader)) {
          //    printf("ERROR:  Could not find ellipsoid shape information for overlap checking.\n");
          //} else {
          //    (void) poReflnlist->nExpandGetField(poReflnlist->ms_sfEllipsoidA00);
          //    (void) poReflnlist->nExpandGetField(poReflnlist->ms_sfEllipsoidA01);
          //    (void) poReflnlist->nExpandGetField(poReflnlist->ms_sfEllipsoidA11);
          //    (void) poReflnlist->nExpandGetField(poReflnlist->ms_sfEllipsoidb0);
          //    (void) poReflnlist->nExpandGetField(poReflnlist->ms_sfEllipsoidb1);
          //    (void) poReflnlist->nExpandGetField(poReflnlist->ms_sfEllipsoidc);
          //
          //    // Set the ellipsoid fields for each.
          //    for (nx = 0; nx<poReflnlist->nGetNumReflns(); nx++) {
          //        oAreaProfiles.nGetPeakArea((int) (*poReflnlist)[nx].fGetField(poReflnlist->m_nFI_fCalcPx0),(int) (*poReflnlist)[nx].fGetField(poReflnlist->m_nFI_fCalcPx1),oArea);
          //        (*poReflnlist)[nx].nPutGetEllipse(true,&oArea.a2x2fEllipsoidA[0][0],&oArea.a2fEllipsoidb[0],&oArea.fEllipsoidc);
          //    };
          //
          //    
          //    nResultField = poReflnlist->nGetFieldIndex(poReflnlist->ms_snNonunfFlag);
          //    nFirstField = poReflnlist->nExpandGetField("nMergeHead");
          //    nNextField = poReflnlist->nExpandGetField("nMergeNext");
          //    
          //    nTemp = poReflnlist->nOverlapCheck(nFirstField,nNextField,poPredict->m_poRotation->fGetIncrement());
          //    
          //    // Check to see that we have ellipsoid information.
          //    
          //    
          //    if (nTemp == -1) {
          //        printf("ERROR:  Could not initialize overlap checking!\n");
          //    } else {
          //        // Go back and set the nonunf field.
          //        for (nx = 0; nx<poReflnlist->nGetNumReflns(); nx++) {
          //            if (((*poReflnlist)[nx].nGetField(nNextField) != -1) || ((*poReflnlist)[nx].nGetField(nFirstField)!=nx)) {
          //                (*poReflnlist)[nx].vSetField(nResultField,0x3);
          //            };
          //        };
          //        
        if( nTemp < 0 )
        {
            printf("\nFailed to predict reflection overlap\n\n");
        }
        else
        {
            printf("There were %d predicted overlaps out of %d predicted reflections\n", 
                 nTemp,
                 poReflnlist->nGetNumReflns());
          
            if( 0 != poReflnlist->nGetNumReflns() )
            {
                printf("Predicted overlap percentage is: %.1f%%\n", (double)nTemp/(double)poReflnlist->nGetNumReflns()*100.0);
            }
            printf("\n");
        } 
          //    };
          //};
      }
      
      nTemp = poReflnlist->nWrite(sReflnFile);
      
      if (0 == nTemp)
      {
          cout << "dtpredict: Reflections written to " << sReflnFile << endl << flush;
#ifndef NO_X_WINDOWS
#ifndef SSI_PC
          if (NULL != poXprop)
          {
              poXprop->hSetProperty("DTDISPLAY_REFLN_UPDATE",
                  sGetCWD() + sReflnFile);
          }
#else
          if(poPredict->m_nDisplay) {
              if(anSeqNumStart.size() > 0) {
                  poScan->vSetSeqNum(anSeqNumStart[0]);
                  CCrclHelper::GetInstance()->vSendPredictUpdateDisplay( poScan->sGetImageName(), sReflnFile );
              }
              else {
                  CCrclHelper::GetInstance()->vSendPredictUpdateDisplay( sReflnFile );
              }
              /*Cstring Result;
              CCoreTclParam TclCom("PredictUpdateDisplay");
              TclCom.SetParamString("RefList", sReflnFile);
              TclCom.SetParamString("Server", "Dtrek");
              TclCom.SetParamString("Step", "Predict");
              // Tell the calling process what to display from the ref list.
              TclCom.SetParamInt("ShowObserved", 0);
              TclCom.SetParamInt("ShowCalculated", 1);
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
          cout << "ERROR: dtpredict - writing file: " << sReflnFile << endl;
          if (0 == nStat)
          {
              nStat = nTemp;  // Don't overwrite prev nStat if bad
          }
      }
  }
  if (bWriteNewHeader)
    {
      int nx;
      // Be sure to update any updated or overridden crystal mosaicity
      for (nx=0;nx<poPredict->m_nNumCrystals;nx++)
        poPredict->m_ppoCrystal[nx]->nUpdateHeader(m_poHeader);
      
      // Check if the output header should have info merged from the
      // image header and the input header

      if (bMergeHeadersScan || bMergeHeadersDet)
	{
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
          cout << "INFO: dtpredict - merging image header info into input header.\n";
	}
      (void) m_poHeader->nWrite(sOutFile);
      cout << "INFO: dtpredict - writing header file: " << sOutFile << endl << flush;
      cout << "\n\nINFO - Examine the predictions with dtdisplay."
           << "\n       The suggested next step is to run dtmultistrategy"
           << "\n       or dtintegrate with " << sOutFile << " as the input .head file:\n"
           << "\n         dtmultistrategy " << sOutFile << " ..." 
           << "\n         dtintegrate " << sOutFile << " ...\n\n" << flush;

    }
#if  ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  if (NULL != poXprop)
    delete poXprop;
#endif
  if (NULL != poPredict)
    delete poPredict;
  if (NULL != poScan)
    delete poScan;

  delete poReflnlist;

  return (nStat);
}

void CDTMainPredict::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << '\n';
    }
  cout << "\ndtpredict - Usage:\n"
       << "dtpredict sHeaderfile [ options ...]\n\n"
       << "Command line options: Description:\n\n"
       << "sHeaderfile         Name of a file with an input header that has the scan\n"
       << "                    definition.\n\n"
       << "-mosaicity fMosaicity\n"
       << "                    Specify mosaicity constant\n\n"
       << "-differentscan  or\n"
       << "-mergeheadersscan   Merge the crystal goniometer, scan, rotation, and\n"
       << "                    source info from the image header into input header and\n"
       << "                    write out the merged header (implies -out option).\n"
       << "                    This option copies the (refined) detector position\n"
       << "                    from the input header to the output header.\n"
       << "                    The resultant header should be suitable for input to\n"
       << "                    dtintegrate to integrate the scan of images.\n\n"
       << "-differentdet  or\n"
       << "-mergeheadersdet    Merge the crystal goniometer, detector goniometer,\n"
       << "                    scan, rotation, and source info from the image header into \n"
       << "                    input header (i.e. implies -mergeheaderscan and -out).\n"
       << "                    This option copies the (unrefined) detector position\n"
       << "                    from the image header to the output header.\n"
       << "                    The resultant header should be suitable for input to\n"
       << "                    dtintegrate to integrate the scan of images.\n\n"
       << "-out  sOutputHeader Write a new header file, possibly with new mosaicity\n"
       << "                    and merged info from the action of -merge* options.\n\n"
       << "-ref  sReflectionFileName\n"
       << "                    Specifies an alternate output reflection file\n"
       << "                    instead of the default *dtpredict.ref. The file format\n"
       << "                    is a standard d*TREK reflection file format.\n\n"
       << "-reso fResoMin fResoMax\n"
       << "                    Specifies the resolution range to predict\n"
       << "                    reflections for.  The default range is determined\n"
       << "                    by the detector position and source wavelength.\n\n"
       << "-reso edge          Set upper resolution to be equal to the highest\n"
       << "                    *minimum* resolution found on any of the four\n"
       << "                    edges (sides) of the image.\n\n"
       << "-reso corner        Set upper resolution to be equal to the highest\n"
       << "                    *maximum* resolution found on any of the four\n"
       << "                    edges (sides) of the image.\n\n"
       << "-rot  fStart fEnd   Specifies a different rotation angle start\n"
       << "                    and end in degrees that overrides the values\n"
       << "                    found in the input image header or a previous\n"
       << "                    -seq option.\n"
       << "-seq  nStart nEnd   Specifies a different rotation angle start\n"
       << "                    and end in degrees that overrides the values\n"
       << "                    found in the input image header or a previous\n"
       << "                    -rot option.  Default: all images in a scan.\n\n"
       << "-nonunf             Do NOT keep reflections predicted in bad non-\n"
       << "                    uniformity regions.  Default: keep them, but flagged.\n\n"
       
       //<< "-overlap\n"
       //<< "                    Turns on overlap checking.  Default is no checking.\n"
       //<< "                    Overlapped reflections have a non-zero nNonunf_flag field.\n\n"
       
       << "-overlap [nDetCoordTest]\n"
       << "                    For each scan solution calculate the total percentage of\n"
       << "                    overlapped HKL's. Two HKL's are considered overlapped\n"
       << "                    if (1) they are predicted on the same image and (2) on that\n"
       << "                    image the absolute difference between both corresponding\n"
       << "                    detector coordinates is less than nDetCoordTest (pix). If\n"
       << "                    nDetCoordTest is not supplied or supplied as -1, the test\n"
       << "                    value will be set to \"1.2 * average_spot_size\" where spot\n"
       << "                    size comes from input header. If the header does not have\n"
       << "                    that information, the test value will be set to 20 pix.\n\n"
       
       << "-display            Tell a dtdisplay process to display the reflection\n"
       << "                    list.  Default: do not tell dtdisplay.\n\n"
       << "-twin nComp1 [ nComp2 ...]\n"
       << "                    Specify the twin component(s) to predict. 1 will predict the\n"
       << "                    first component.  By default, all components are predicted.\n\n"
       << "-maxlorentz fML     Set maximum allowed Lorentz factor.  The value will\n"
       << "                    determine which reflections will be flagged as too close\n"
       << "                    to the scan rotation axis.  Default: 50.\n\n"
       << "-list               List properties used in the prediction.\n\n"
       << "-lscale fScale      Factor to multiple cell lengths by.  Default: 1.\n\n"
       << "-help               Print this help text.\n\n"
       << "Information about scan angles and non-uniformity of response is\n"
       << "determined from the headers of the input image files.  If the\n"
       << "headers are incorrect, they may be edited with dtheaderedit or dtprocess.\n\n"
       << "Examples:\n\n"
       << "   dtpredict dtrefine.head -seq 1 1\n"
       << "   dtpredict dtrefine.head -seq 201 201 -differentscan\n"
       << "   dtpredict dtrefine.head -seq 1 1 -mosaicity 0.5 -differentdet\n"
       << "   dtpredict dtintegrate.head -rot -5 20.2\n"
       << "   dtpredict hivrt.head -ref hivrt.ref\n" << flush;
}
//////////////////////////////////////////////////////////////////////////
bool CDTMainPredict::bParseReflnOverlapCheckParams(Cstring& sError)
{
    sError = ""; // initialize
    
    Cstring                 sArgs("");
    if( !bGetCommandOption(c_pcPredictCheckOverlap, sArgs) )
        return false; // option -overlap not given    

    if( !m_stOverlapCheckParams.bParseInput(sArgs) )
    {
        sError = "Option ";
        sError += Cstring(c_pcPredictCheckOverlap);
        sError += " has wrong number of parameters.";
        return false;
    }
	
    // For now we cannot handle "zero" case
    if( 0.0 == m_stOverlapCheckParams.m_dDistTest )
    {
        sError = "Option ";
        sError += Cstring(c_pcPredictCheckOverlap);
        sError += " has zero value for nDetCoordTest";
        return false;
    }

    return true;
}
////////////////////////////////////////////////////////////////////////'


