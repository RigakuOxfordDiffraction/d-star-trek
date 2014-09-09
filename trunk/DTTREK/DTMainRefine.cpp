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
// DTMainRefine.cpp  Initial author: RB           15-Mar-2004
//
//  This file contains the member functions of class CDTMainRefine
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
//    Example:  % dtrefine input_header input_reflnlist [options]
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "DTMainRefine.h"

#include "Crefine.h"       // Includes many other include files
#include "Creflnlist.h"
#include "Cpredict.h"
#include "Cscan.h"
#include "CIceRingLocator.h"
#include "dtsvd.h"

#include "CCommandLine.h"

#if !defined(VC6)
using std::cout;
using std::cin;
using std::endl;
using std::flush;
using std::cerr;
#endif

static char*     c_pcSigma       = "-sigma"; 
static char*     c_pcETE         = "-rankete";

CDTMainRefine::CDTMainRefine()
{
    m_nNumberOfResoBins = 10;

    m_nMinReflnsForRankStats = 5;
}

int CDTMainRefine::nExecute(unsigned int argc, char* argv[])
{

  int            i = 0;
  float          fTemp = 0.0f;
  int            nStat = 0;
  Cstring        sReflnIn("");
  Cstring        sTemp("");
  Cstring        sRefineFiles("");
  Cstring        sRefineOptions("");
  Cstring        sPredReflnlist("");
  Creflnlist    *poReflnlistIn          = NULL;
  Crefine       *poRefine;
  bool           bPredict   = FALSE;
  bool           bNoPredict = FALSE;

  float          fFindRingsWidthToExclude = 2.0;
  bool           bFindRings = false;
  bool           bRankEvaluateExposureTime = false;
  double         fRankTargetIOverSigmaAtTargetReso = 10.0;
  double         fIOverSigmaTimePower  = 0.4;   // I/Sigma as a function of time changes proportionally to time ^ (Power)   

   // RB: THIS DEBUG CODE IS JUST TO STOP THE DTREK PROCESS WHEN STARTED FROM CC:
   //int*  pNull = NULL;
   //*pNull = 1;

  vDtrekSetModuleName("dtrefine");
  vPrintCopyrightInfo();
  //cout << "\ndtrefine:  Copyright (c) 2006 Rigaku\n";
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

#ifdef SSI_PC
  if( 1 == argc )
  {
      DTREK_ERROR(0, "");
  }
#endif

  // Copy command line files to sRefineFiles
  for (i = 1; (i < argc) && (i < 3); i++)
    {
      sRefineFiles = sRefineFiles + ' ' + (const char*) argv[i];
    }

  // Get header file to start with

  argc--;
  argv++;

  nStat = 0;

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

    // Read the name of a file and read the header only from that image file
  if ( !bCreateHeader(sTemp) )
    {
      DTREK_ERROR(2, "Header not available!\n");
    }

  // Get reflection list file to start with

    argc--;
    argv++;

    if (1 > argc)
    {
      // No reflection list argument and thats an error!

      DTREK_ERROR(3, "ERROR - no reflection list filename!\n");
    }

    //+jwp 24-Aug-1999

    int     nSeq = 0;
    int     nSeqTo = 0;
    int     nNumSeqOptions = 0;
    int     *pnSeqNums = NULL;

    // Look for reflnlist or -seq options after the header filename
    sTemp = (const char*) argv[0];
    sReflnIn = sTemp;

    while( "-seq" == sTemp )
    {
        // If reflnlist file has name "-seq", then ...
        // ... get the next argument as an image number,
        // ... predict spots on that image,
        // ... get centroids of predicted spots

        // Read the sequence number that follows -seq
        argc--;
        argv++;

        nSeq = atoi(argv[0]);
        argc--;
        argv++;

        if( argc > 0 && 1 == sscanf(argv[0],"%d",&nSeqTo) )
        {
            argc--;
            argv++;
        } 
        else 
            nSeqTo = nSeq;

        if( 0 == nNumSeqOptions )
        {
            // Create a reflnlist for later use
            poReflnlistIn = new Creflnlist();
        }

        int*        pnTemp = new int[nNumSeqOptions+abs(nSeqTo-nSeq+1)];

        if( NULL != pnSeqNums )
        {
            for (i = 0; i < nNumSeqOptions; i++)
                pnTemp[i] = pnSeqNums[i];
            
            delete [] pnSeqNums;
        }
        
        pnSeqNums = pnTemp;
        
        for (;nSeq<=nSeqTo;nSeq++)
        {
            pnSeqNums[nNumSeqOptions] = nSeq;
            nNumSeqOptions++;
        }

        // Look ahead to see if next argument is also a -seq option
        sTemp = (const char*) argv[0];
    }


  if (0 >= nNumSeqOptions)
    {
      argc--;
      argv++;

      // There were no -seq options, so MUST be a reflnlist filename

      // Read a reflection list from a file

      poReflnlistIn = new Creflnlist(sTemp);

      if (!poReflnlistIn->bIsAvailable())
	{
	  DTREK_ERROR(4, "ERROR: Reflection list not available!\n");
	}
  } 

  // Create refine object
  poRefine   = new Crefine(*m_poHeader, poReflnlistIn);

  float     fResoMin = 0.0f;
  float     fResoMax = 0.0f;
  bool      bResolutionProvided = false;

  double     fSig = 0.0;
  if( bGetCommandOption(c_pcSigma, fSig, false) )
      poRefine->m_fSigma = (float)fSig;
  
  ////////////////////////////////////////////////////////////
  // Determine parameters for exposure time evaluation
  std::vector<double>       afETEparams;
  if( bGetCommandOption(c_pcETE, afETEparams, " ", false) )
  {
    bRankEvaluateExposureTime = true;
    if( afETEparams.size() >= 1 )
        fRankTargetIOverSigmaAtTargetReso = afETEparams[0];
    if( afETEparams.size() >= 2 )
        fIOverSigmaTimePower  = afETEparams[1];
  }
  ///////////////////////////////////////////////////////////

  // Copy command line options to sRefineOptions
  
  for (i = 0; i < argc; i++)
  {
      if (0 == strcmp("-predict",(const char*)argv[i]))
      {
          bPredict = TRUE;
      }
      else if (0 == strcmp("-nopredict",(const char*)argv[i]))
      {
          bNoPredict = TRUE;
      }
      else if (0 == strcmp("-out",(const char*) argv[i]))
      {
          i++;
          if (i < argc) 
              poRefine->m_sOutHeader = argv[i];
      } else if (0 == strcmp("-intmask",(const char*) argv[i]))
      {
          i++;
          if (i < argc)
              poRefine->m_sMapFile = argv[i];
      } 
      else if (0 == strcmp("-pad",(const char*) argv[i]))
      {
          i++;
          if ((i < argc) && (1==sscanf(argv[i],"%f",&fTemp))) 
              poRefine->m_fFindPad = fTemp;

      }
      else if (0 == strcmp("-resobins",(const char*) argv[i]))
      {
          i++;
          int       nTemp = 0;
          if( i < argc &&  1 == sscanf(argv[i],"%d", &nTemp) ) 
          {
              m_nNumberOfResoBins = nTemp;

              poRefine->vSetNumberOfResoBins(nTemp);
          }
          ///////////////////////////////////////////////////////////
          i++;     // the second parameter of -resobins is optional. Let's see if it is there...
          if (i < argc && !bIsDTCommandLineOptionString(argv[i]) )
          {
                if( 1 == sscanf(argv[i], "%d", &nTemp) )
                  m_nMinReflnsForRankStats = nTemp;
                else
                  DTREK_ERROR(5, "Error processing -resobins option.\n");
          }
          else
              i--;
      }
      else if (0 == strcmp("-tables",(const char*) argv[i]))
      {
          poRefine->m_bPrintTables = TRUE;
      }
      else if (0 == strcmp("-useoverlapped",(const char*) argv[i]))
      {
          poRefine->m_bUseOverlapped = true;
      }
      else if (0 == strcmp("-rejects",(const char*) argv[i]))
	  {
		  i++;
		  if (i < argc) 
			  poRefine->nSetRejectReflectionFile((Cstring) argv[i]);
		  else
			  DTREK_ERROR(4, "ERROR: Please specify rejection file!\n");
	  }
      else if (0 == strcmp("-rank",(const char*) argv[i]))
	  {
          bFindRings = true;  
          poRefine->vSetRanking(true);
	  }
      else if (0 == strcmp(c_pcETE,(const char*) argv[i]))
	  {
          // skip this option, because it has been parsed earlier
          i++;
          while( i < argc && !bIsDTCommandLineOptionString(argv[i]) )
            i++;
          i--;
	  }
      else if ( 0 == strcmp("-ringdetect",(const char*) argv[i]) )
      {
          bFindRings = true;
          i++;
          double        dTemp = 0;
          if( i < argc &&  1 == sscanf(argv[i],"%f", &dTemp) ) 
          {
              fFindRingsWidthToExclude = dTemp;
          }
      }
      // RB 6/21/05  In the old code the whole -reso option was automatically appended to sRefineOptions.
      // However, Crefine object does not know how to handle new -reso parameters like "edge" etc.  So we need to set the resolution
      // limits explicitly.
      else if( 0 == strcmp("-reso",(const char*) argv[i]) )
      {
            Cstring       sResoParameters("");
            i++;
            
            while( i < argc && ('-' != argv[i][0] && '+' != argv[i][0]) )
            {
                sResoParameters += argv[i];
                sResoParameters += " ";
                i++;
            }
            
            i--;

            Cstring       sError("");
            if( !bParseCalculateResolution(sResoParameters, fResoMin, fResoMax, sError) )
            {
                DTREK_ERROR(6, sError);
            }
            else
            {
                bResolutionProvided = true;
	            
                char        szTemp[128];

                sprintf(szTemp, "-reso %.2f %.2f", fResoMin, fResoMax);
                
                sRefineOptions = sRefineOptions + ' ' + (const char*)szTemp;
            }
      }
      else
          sRefineOptions = sRefineOptions + ' ' + (const char*) argv[i];
  }

  if (0 < argc)
  {
      // Parse additional options
    if (0 < nNumSeqOptions)
	{
	  if (sRefineOptions.contains("-testmosaicity"))
	    {
	      // If -testmosaicity found among the refine options
		  
	      float fVal1 = poRefine->m_fTestMosaicityRangeMax;
	      float fTemp;
	      for (i = 0; i < argc; i++)
		{
		  if (0 == strcmp("-testmosaicity", (const char*)argv[i]) )
		    {
		      break;
		    }
		}
	      if (argc > i+1)
		{
		  nStat = sscanf(argv[i+1],"%f",&fTemp);
		  if (1 == nStat)
		    {
		      argv[i+1] = "-testmosaicity";
		      fVal1 = fTemp;
		    }
		}
	      poRefine->m_fTestMosaicityRangeMax = fVal1;
	      printf("\nTesting mosaicity values of %6.2f degrees and less, see"
		     "\n  Rossmann (1979) J. Appl. Cryst. 12, 225-238.\n",
		     poRefine->m_fTestMosaicityRangeMax);
	      fflush(stdout);
		  
	      nStat = poRefine->nGetMosaicity(nNumSeqOptions,
					      pnSeqNums, m_poHeader, fResoMin,
					      fResoMax,0,NULL,NULL,NULL);
	    } // end of -testmosaicity test
	  
	  if ( (0 < argc) && (0 == nStat) && sRefineOptions.contains("-go") )
	    {
	      // Now that the poRefine object is available, go search images
	      // for refln centroids (if requested)
		  
	      nStat = poRefine->nGetReflnsFromImages(nNumSeqOptions,
						     pnSeqNums,
						     m_poHeader, 
						     fResoMin, 
						     fResoMax);
	      if (0 != nStat)
		  {
		    DTREK_ERROR(4, "ERROR: Reflection list not available!\n");
		  }
	    }
	} 
    else
    {
	  nStat = poRefine->nGetReflnsFromFind(*poReflnlistIn,fResoMin,fResoMax);
	}
  }

  // Do refinement if sRefineOptions has options and contains the word "-go"

  if ( (0 < argc) && (0 == nStat) && sRefineOptions.contains("-go") )
    nStat = poRefine->nDoRefinement(sRefineOptions, poReflnlistIn);

  if (sRefineOptions.contains("-twinreject")) {
      printf("\n\n");
      printf("Listing of all twin rejection information.\n");
      printf("=======================================================================\n");
      printf(poRefine->m_sLogTwinReject.string());
      printf("=======================================================================\n");
      printf("\n\n");
      
  };
  poRefine->vPrintResidualLog(0);


  if (sRefineOptions.contains("-prompt"))
    {
      cout << "dtrefine - PROMPT Write new output header? (y,n) [y]: " << flush;
      getline(cin, sTemp);
      if ("" == sTemp) sTemp = 'y';
      if (!cin)
	nStat = 0;
      else if ( ('n' == sTemp.GetAt(0)) || ('N' == sTemp.GetAt(0)) )
	nStat = 1;
      else
	nStat = 0;
      cout << endl << flush;
    }

    /////////////////////////////////////////////////////
    // RB do some ranking stuff
    if( poRefine->bIsRanking() )
        m_poHeader->nClearDictionary("DTREFINE_RANK", "*");
    /////////////////////////////////////////////////////      

    if (0 == nStat)
    {
      // Write results to output

      // No longer should do this here. It is done inside nDoRefinement(), which has a pointer to the header.
      // poRefine->nUpdateHeader(poHeader);

      m_poHeader->nReplaceValue(Crefine::ms_sDtrefineFiles, sRefineFiles);
      if (sRefineOptions.contains("-go") )
	{
	  // Update refinement options string/scheme only if -go in options

	  m_poHeader->nReplaceValue(Crefine::ms_sDtrefineOptions,
				sRefineOptions);
//+jwp
	  // If using -seq options (-seq n) AND the -predict option,
	  // then also predict the spots for the
	  // first -seq option image at this time

	  if( (1 < nNumSeqOptions && bPredict) || (1 == nNumSeqOptions && !bNoPredict) )
	    {
          float a2fRotRange[2];
	      Cpredict *poPredict;
	      Cscan    *poScan;
	      Creflnlist *poReflnlistPred = NULL;
	      poScan = new Cscan(*m_poHeader);
	      poPredict = new Cpredict(*m_poHeader, NULL);

	      // Create default reflection list:

	      poReflnlistPred = new Creflnlist();

	      if ( (0.0 != fResoMin) && (0.0 != fResoMax) )
		poPredict->vSetResolution(fResoMin, fResoMax);
	      a2fRotRange[0] = poScan->fCalcGetRotStart(pnSeqNums[0]);
	      a2fRotRange[1] = poScan->fCalcGetRotEnd(pnSeqNums[0]);
	      nStat = poPredict->m_poRotation->nSetRotRange(a2fRotRange[0],
							    a2fRotRange[1]);
	      cout << "...predicting reflections for image sequence number "
                   << pnSeqNums[0] << endl;

	      nStat = poPredict->nPredict(NULL, poReflnlistPred);
	      if (0 == nStat)
		{
		  sPredReflnlist = sDtrekGetPrefix() + "dtrefpredict.ref";
		  // TODO:  Copy dtrefinetmp.ref info into this list.

		  nStat = poReflnlistPred->nWrite(sPredReflnlist);
		}
	      delete poPredict;
	      delete poScan;
	      delete poReflnlistPred;
	    }
//-jwp
         /////////////////////////////////////////////////////////////////////////
         // RB Find Rings
         if( bFindRings && NULL != pnSeqNums && nNumSeqOptions >= 1 )
         {
            // Try and locate ice rings (on the first image only )
            CIceRingLocator         oLocator = CIceRingLocator(m_poHeader);
            
            for(int ii=0; ii < nNumSeqOptions; ii++)
                oLocator.vLocateRings(pnSeqNums[ii]);
         }
         /////////////////////////////////////////////////////////////////////////
	    
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( poRefine->bIsRanking() )
        {
            // Use the existent CRefine mechanism to predict spots, and (through CFind) quick-integrate them.
            // The mosaicity value has just been refined, so we will fix it at the refined value and make sure 
            // that Crefine communicates it to Cpredict (DTREK_REFINE_GRFI_PREDICT_MOSAICITY).
            
            cout << endl << "Predicting reflections using refined values and performing quick 2D integration..." << endl;
            
            poRefine->m_fFindPad = poRefine->fGetMosaicity();
            
            nStat = poRefine->nGetReflnsFromImages(nNumSeqOptions,
			                                       pnSeqNums,
			                                       m_poHeader, 
                                                   fResoMin, 
                                                   fResoMax,
                                                   DTREK_REFINE_GRFI_FOR_RANKING | DTREK_REFINE_GRFI_PREDICT_MOSAICITY);

	        if( 0 == nStat )
		    {
                ////////////////////////////////////////////////////////////////////////////////////////////////
                poRefine->m_poReflnlistIn->vRemoveIntensityOrSigmaMarkedReflections();
                /////////////////////////////////////////////////////////////////////////////////////////////////
                
                ////////////////////////////////////////////////////////////////////////////////////////////////
                Cscan*      poScan = new Cscan(*m_poHeader);
                poRefine->m_poReflnlistIn->vRemoveOverlapReflns(poScan->m_poRotation->fGetIncrement());
                delete poScan;
                ////////////////////////////////////////////////////////////////////////////////////////////////

                if( poRefine->bIsDisplay() )
                {
                    poRefine->m_poReflnlistIn->vSetColorsForDisplay();

                    // Since ranking has been requested, we must have searched for reso rings already
                    poRefine->m_poReflnlistIn->vRemoveReflnsInResolutionRingAreas(m_poHeader, 
                                                                                  (double)fFindRingsWidthToExclude,
                                                                                  DTREK_REFLNLIST_DL_FROM_REFINE | 
                                                                                  DTREK_REFLNLIST_DL_MARK_ONLY |
                                                                                  DTREK_REFLNLIST_DL_RANK);
                    
                    poRefine->m_poReflnlistIn->nWrite(sDtrekGetPrefix() + "dtranker_refine_pred.ref");  
                }

                m_poHeader->nReplaceValueDictionary("DTREFINE_RANK", "PREDICTED_REFLN_COUNT", poRefine->m_poReflnlistIn->nGetNumReflns());
                
                // Since ranking has been requested, we must have searched for reso rings already
                poRefine->m_poReflnlistIn->vRemoveReflnsInResolutionRingAreas(m_poHeader, 
                                                                              (double)fFindRingsWidthToExclude,
                                                                              DTREK_REFLNLIST_DL_FROM_REFINE |
                                                                              DTREK_REFLNLIST_DL_RANK);
                
                poRefine->m_poReflnlistIn->nResoTable(bResolutionProvided ? fResoMin : -1, 
                                                      bResolutionProvided ? fResoMax : -1,
                                                      m_nNumberOfResoBins,
                                                      10, // 10 I/Sigma bins
                                                      eTableVsResolution | eTableVsIoverSig,
                                                      -1,
                                                      m_poHeader,
                                                      DTREK_REFLNLIST_RT_RANK);        
                
                bool    bETE = bRankEvaluateExposureTime && bResolutionProvided;

                poRefine->m_poReflnlistIn->vRank(m_poHeader,
                                                 bResolutionProvided ? fResoMin : -1,
                                                 bResolutionProvided ? fResoMax : -1,
                                                 m_nNumberOfResoBins,
                                                 m_nMinReflnsForRankStats,
                                                 bETE);

                if( bETE )
                {
                    double      fTemp = fResoMax;                  
                    poRefine->m_poReflnlistIn->bEvaluateExpTime(m_poHeader, 
                                                                fTemp, 
                                                                fRankTargetIOverSigmaAtTargetReso,
                                                                fIOverSigmaTimePower);
                }
            }
        }
    }
    
      nStat = m_poHeader->nWrite(poRefine->sGetOutHeader());
      if (0 == nStat)
	{
	  (void) poRefine->nNotifyDisplay(sPredReflnlist, sGetCWD()
                                          + poRefine->sGetOutHeader());
          cout << "dtrefine - INFO wrote header file: " << poRefine->sGetOutHeader() << flush;
	  if ("-seq" != sReflnIn )
	    {
              cout << "\n\nINFO - The suggested next step is to run dtrefine with"
                   << "\n       " << poRefine->sGetOutHeader() << " as the input"
                   << " and the -seq option:\n"
                   << "\n       dtrefine " << poRefine->sGetOutHeader() << " -seq ...\n\n" << flush;
            }
          else
            {
              cout << "\n\nINFO - The suggested next step is to run dtrefine again:"
                   << "\n         dtrefine " << poRefine->sGetOutHeader() << " -seq ..."
                   << "\n      -or-"
                   << "\n       Run dtpredict to check that the results predict the"
                   << "\n       first few images in each scan correctly: "
                   << "\n         dtpredict " << poRefine->sGetOutHeader() << " -seq ..." 
                   << "\n      -or-"
                   << "\n       Run dtmultistrategy or dtintegrate with " 
                   << "\n       " << poRefine->sGetOutHeader() << " as the input .head file:\n"
                   << "\n        dtmultistrategy " << poRefine->sGetOutHeader() << " ..." 
                   << "\n        dtintegrate " << poRefine->sGetOutHeader() << " ...\n\n" << flush;
            }
	}
      else
	{
	  cout << "dtrefine - ERROR writing header file: ";
          cout << poRefine->sGetOutHeader() << endl << flush;
	}
    }
  else
    {
      cout << "dtrefine - WARNING did not write new header file.\n" << flush;
    }

  delete poRefine;
  delete poReflnlistIn;
  if (NULL != pnSeqNums)
    {
      delete [] pnSeqNums;
      pnSeqNums = NULL;
    }
      return (nStat);
}

void CDTMainRefine::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
    {
      cout << sMessage << endl;
    }
  else {
      cout << "\ndtrefine - Usage:\n"
          "dtrefine header_file refln_file|-seq nSeq [[ options ...] -go] -go \n\n"
          "Command line options: Description:\n\n"
          " header_file          Name of a file with a d*TREK header used to\n"
          "                      initialize the refinement procedure.\n\n"
          " refln_file           Name of a d*TREK reflection list file (no default)\n"
          "                      which will be used for refinement.\n"
          "           -OR-\n"
          " -seq nSeq            Instead of a reflection list file, specify the sequence\n"
          "                      numbers of one or more images to predict reflections\n"
          "                      for and get centroids from. Example: -seq 1 -seq 100\n\n"
          " -converge            Refine until all shifts are less than 0.00003 or\n"
          "                      -cycles nCycles is reached.  Use -noconverge to always\n"
          "                      refine for nCycles.\n\n"          
          " -cycles nCycles      Number of refinement cycles.  Default: 30.\n"
          "                      Note: fewer cycles may be performed if converged.\n"
          "                      It is better to have more -go's than more cycles since\n"
          "                      each -go recalculates derivatives.\n\n"
          " -display             Tell a dtdisplay process to display the reflection\n"
          "                      list.\n\n"
          " -go                  Do refinement with available info.\n\n"
          " -out sOutputHeader   Output header file name.  Default: dtrefine.head.\n\n"
          " -prompt              Prompt for whether to write new header file or not.\n\n"
          " -rej fRejX fRejY fRejRot\n"
          "                      Spatial rejection criteria.  If calculated position\n"
          "                      of a refln differs by more by fRejX, fRejY or fRejRot\n"
          "                      in X mm, Y mm or Rot degrees from the observed\n"
          "                      position, then that refln is not used.  Defaults: 2.5,\n"
          "                      2.5, 2.0.\n\n"
          " -reso fResoMin fResoMax\n"
          "                      Resolution limits of input reflections to use.\n"
          "                      Default: 99999 0.\n\n"
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
          " -sigma fSigma        I/sigmaI minimum for input reflections.  Default: 0.\n\n"
          " -sharpness fSharp    Filters out 'spots' that have a sharpness greater\n"
          "                      than fSharp.  This tends to improve indexing and\n"
          "                      refinement, since it excludes garbage and/or\n"
          "                      amorphous spot profiles.\n"
          "                      Default:  fSharp = 0.2\n\n"
          " -verbose nVerbose    Verbosity level of output.  Higher values give\n"
          " -v nVerbose          more output.  Default: 1.\n\n"
          " -detmos              Print out detector dependent mosaicity results.\n"
          "                      Default: Disabled\n\n"
          " -wideslice\n"
          " -nowideslice         Enable or disable wide slice algorithms.\n"
          "                      By default, we use them if and only if the image\n"
          "                      widths are 3 degrees or greater.\n\n"
          " -resid fResid        Stop when below this residual.\n"
          "                      Default: 1e-20.\n\n"
          " -maxlorentz fML      Set Maximum lorentz factor of reflections used in refinement.\n"
          "                      Default: 35\n\n"
		  " -rejects sRejectFile Output rejects to a file.  This can be used to locate\n"
		  "						 and index twin components.\n\n"
          " -intmask sMapImage   Write map file containing found ellipsoids.\n\n"
          " -pad nPad            Number of images to pad by during -seq refinements.\n\n"
          " -testmosaicity [ fMosaicityMax ]\n"
          "                      Derive mosaicity from images by finding mosaicity between\n"
          "                      below fMosaicityMax that predicts at least 98-99 percent\n"
          "                      of the spots found in images specified by -seq options.\n"
          "                      Default: fMosaicityMax = 3.0\n"
          "                      Note: no parameter refinement will be performed.\n\n"
          " -weight sWeightID    Weight data according to different data properties.\n"
          "                      Options are:\n"
          "                      RotSigma   Weight with observed rotation sigma (default).\n"
          "                      Lorentz    Weight with inverse Lorentz factor.\n"
          "                      Intensity  Weight with reflection intensity.\n"
          "                      IoverSig   Weight with I/sig value.\n"
          "                      RotWidth   Weight with rotation width.\n"
          "                      None       Use unit weights (no weighting).\n\n"
          " -reflimit nMaxRefs   Maximum reflections used during refinement. Do not use\n"
          "                      this option when refining off a reflection list obtained\n"
          "                      from a 2D find, or in conjunction with either the -reso\n"
          "                      or -seq options.\n\n"
          " -twinreject [fMaxUnaccountedIntensity] [fMaxRMS-MM]\n"
          "                       Add solutions until at most a small fraction of the\n"
          "                       data is accounted for. Any solution\n"
          "                       having an RMS-MM residual greater than fMaxRMS-MM is\n"
          "                       rejected.\n\n"
          " -checkdoubled [fMinFractionInOdd]\n"
          "                       Checks all twin components for doubled axes.  Each\n"
          "                       component is examined to see how many reflections\n"
          "                       fall in odd/even parity classes.  If no more than\n"
          "                       fMinFractionInOdd lie in the odd parity class, the\n"
          "                       algorithm concludes that there is a doubled axis, and\n"
          "                       acts to remove the problem.  A similar check is made\n"
          "                       for tripled and quadrupled axes.\n\n"
          " -useoverlapped        When reflection centroids are obtained from images,\n"
          "                       do keep the reflections that are overlapped. If this\n"
          "                       option is not specified, such reflections are rejected.\n\n"
          "Refinement flags which may be preceded with a - (fix, do not refine)\n"
          "or a + (refine, do not fix):\n\n"
          " CrysAstar, CrysBstar, CrysCstar, CrysAlps, CrysBets, CrysGams\n"
          "                       Unit cell parameters.\n"
          " CrysCell, CrysLengths, CrysAngles\n"
          "                       Convenience flags for all or some cell parameters\n\n"
          " CrysRot1, CrysRot2, CrysRot3, CrysRot\n"
          "                       Crystal orientation angles, individually or collectively\n\n"
          " CrysMosaicity         Crystal effective mosaicity\n\n"
          " CrysAll               Same as CrysCell, CrysRot, CrysMosaicity\n\n"
          "                       For twinned samples using '-twin 0', each -go expands to:\n"
          "                       +CrysXXXXX -twin 1 -go -twin 2 -go [...]\n\n"
          " CrysConstraints       Use crystal lattice constraints as required by the\n"
          "                       lattice system. (For example, alp=bet=gam in orthorhombic,\n"
          "                       or a==b in tetragonal).  By default this is ENABLED.\n\n"
          " DetTrans1, DetTrans2, DetTrans3, DetTrans\n"
          "                       Detector translations, individually or collectively\n"
          "                       (crystal to detector distance is DetTrans3)\n\n"
          " DetRot1, DetRot2, DetRot3, DetRot\n"
          "                       Detector rotations, individually or collectively\n"
          "                       (detector swing angle is DetRot2)\n\n"
          " DetAll                Same as DetTrans, DetRot\n\n"
          "                       For twinned samples using '-twin 0', each -go yields:\n"
          "                       +DetXXXX -twin -1 -go\n\n"
          " GonioMiss1, GonioMiss2, GonioMiss\n"
          "                       Refine goniometer miss-setting matrices.  These matrices\n"
          "                       will adjust the value of the goniometer rotation vector\n"
          "                       attached to the physical crystal.\n\n"
          " SourceWave, SourceRot1, SourceRot2, SourceRot\n\n"
          " All                   All parameters.\n\n";
          "                       For twinned samples using '-twin 0', each -go expands to:\n"
          "                       +CrysAll -twin 1 -go -twin 2 -go [...] +DetAll\n"
          "                       -twin -1 -go -DetAll -twin -2 -go [...]\n\n";
      cout <<"Examples:\n"
          " dtrefine dtindex.head dtfind.ref +All -verbose 0 -go -go -go -go -verbose 1 -go\n"
          " dtrefine dtrefine.head -seq 1 -seq 2 +All -verbose 0 -go -go -go -verbose 1 -go\n"
          " dtrefine dtintegrate.head dtintegrate.ref -sigma 15 +All -go -go -go -go\n\n"
          << flush;
  }
}
