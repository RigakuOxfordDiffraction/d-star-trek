//
// Copyright (c) 2005 Rigaku/MSC, Inc.
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainStrategy.cpp       Initial author: RB          21-Mar-2005
// This file contains the member functions of class CDTMainStrategy

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

#include "DTMainStrategy.h"

#include <stdio.h>        // for sprintf, until we figure out sstream

#include "Dtrek.h"

#include "dtrekvec.h"
#include "Cstring.h"
#include "Creflnlist.h"
#include "Cimage_header.h"
#include "Cstrategy.h"
#include "Cscan.h"
#include "Cpredict.h"
#include "dtarray.h"
#include "CExperimentalStrategy.h"

#include "ResoBins.h"


#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::ios;
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CDTMainStrategy::CDTMainStrategy():
m_poStrategy(NULL), 
m_poPredictReflnList(NULL), 
m_bExternalPredictReflnListPtr(false),
m_poResoBins(NULL),
m_bExternalResoBinsPtr(false),
m_poExternScanInfo(NULL)
{
    vInit();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CDTMainStrategy::~CDTMainStrategy()
{
    vDeleteStrategyObject();
    
    vDeletePredictReflnList();
    
    vDeleteResoBins();
}
/////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vInit()
{
     vDeleteStrategyObject();
     
     vDeletePredictReflnList();
     
     vDeleteResoBins();

     m_fLengthScaleFactor = -1.0;

     m_bAnom = false;

     m_wVerbose = 0U;

     m_fUserInputMosaicity = -1.0;
}
//////////////////////////////////////////////////////////////////////////////////////////
double  CDTMainStrategy::dGetTotalRotWidth()
{
    if( m_poStrategy )
        return m_poStrategy->dGetTotalRotWidth();
    else
        return 0.0;
}
/////////////////////////////////////////////////////////////////////////////////////////
int CDTMainStrategy::nExecute(unsigned int argc, char* argv[])
{
  //int*      pNULL = NULL;
  //*pNULL = 0;
#ifdef WIN32
  ios::sync_with_stdio(); 
#endif

  Cstring     sError("");

  int            i, j, k;
  int            nStat;
  int            argc_start = 0;
  Cstring        sTemp,sTemp2;
  double         f0;

  Cstring        sReflnlistPredicted;
  Cstring        sReflnlistPrevious;
  Cstring        sHeaderOut;

  Cpredict      *poPredict        = NULL;
  double         fDistance        = cdStratDetUnrealisticValue;
  double         f2ThetaSwing     = cdStratDetUnrealisticValue;
          
  int            nMaxScan;
  int            nMinScan;
  double         fDebyeRingCoverage;

  float          fTemp;
  int            nTemp;
  bool           bDoList = TRUE;
  
  Cstring        sReflnFile;
  Cstring        sMissingOption
                    = "ERROR: dtstrategy - missing option argument!";
  Cstring        sInvalidArgument
                    = "ERROR: dtstrategy - invalid argument!";

  int            nNumResoBins = -1;

  DTREK_WORD        wListComplMsgCtrl = 0U;

  fDebyeRingCoverage = 20.0;
    
  sHeaderOut       = sDtrekGetPrefix() + "dtstrategy.head";

  vDtrekSetModuleName("dtstrategy");
  vPrintCopyrightInfo();
  //printf("\ndtstrategy:  Copyright (c) 2006 Rigaku\n");
  //printf("%s\n",D_K_DTREKVersion);

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
  cout << "Command line:" << endl << sGetCommandLine(argc, argv, 71).string() << endl << endl << flush; 

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
  
    if( !m_poHeader ) // This will not be NULL, if someone has set that member before calling this function
    {
        //// Read the name of a file and read the header only from that image file
        bCreateHeader(sTemp);
    }
    else
        printf("\nUsing external header object.\n\n");


    argc_start = 1;

    // Check for a second header file.  The second header file can (optionally) contain crystal information.
    if (argc >= 2)
    {
      sTemp2 = argv[1];
    
      if ((sTemp2.after('.')) == (sTemp.after('.')))
      {
          Cimage_header oHeader2(sTemp2);
      
          if (oHeader2.bIsAvailable())
          {
              // We have two header files!
              // Only copy crystal information over from the other header file.
              
              argc_start = 2;
              
              m_poHeader->nCopyMask(oHeader2, D_K_CrystalMask);
          }
      }
    }
  

    if (!m_poHeader->bIsAvailable())
    {
        DTREK_ERROR(2, "ERROR - Header not available!\n");
    }

    /////////////////////////////////////////////////////////////////////
    
    vCreateStrategyObject();
    
    /////////////////////////////////////////////////////////////////////
    // Determine resolution range
    double                        a2fResolution[2] = {-1.0, -1.0};

    if( !m_bExternalResoBinsPtr )
    {
        CExperimentalStrategy         oResoStrat(*m_poHeader);
        if( !oResoStrat.bIsAvailable() )
        {
          printf("ERROR:  Could not initialize CExperimentalStrategy object from header!\n");

          nStat = 1;
          goto end;
        }
        
        if( oResoStrat.nParseResoCommandLine(argc - argc_start, 
                                             argv + argc_start, 
                                             a2fResolution[0],
                                             a2fResolution[1], 
                                             f2ThetaSwing, 
                                             fDistance) )
        {
            nStat = 1;
            goto end;
        }
        
        if( a2fResolution[0] > 0.0 ) // If resolution is supplied on the command line, create a reso object and pass it to strategy
                                     // otherwise strategy will have to figure out resolution on its own
        {
            m_poResoBins = new CResoStrategyBins(a2fResolution[0], a2fResolution[1], nNumResoBins);

            m_poStrategy->vSetExternalResoBinsPtr(m_poResoBins);
        }
    }
    else
    {
        ///////////////////////////
        //consistency check
        if( !m_poResoBins )
        {
            nStat = 1;
            goto end;
        }
        ///////////////////////////

        a2fResolution[0] = m_poResoBins->dGetMinReso();
        a2fResolution[1] = m_poResoBins->dGetMaxReso();
        
        m_poStrategy->vSetExternalResoBinsPtr(m_poResoBins);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////
    // Discover how many scans there are in the header.  
    // This is important, since we want to adjust values in the options below.
    nMaxScan = -1;
    nMinScan = 1;

    for (i = 0; 1; i++) 
    {
        Crotation       oRotation(*m_poHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, i));
        
        oRotation.ms_nVerbose = 0;
        if( oRotation.bIsAvailable() )
        {
            nMaxScan = i;
            nMinScan = min(0,i);
            if (i > 1000)
            {
              DTREK_ERROR(2,"ERROR:  Too many scans!\n");
            }
        } 
        else if (i>=1)
          break;
    }
  
    Crotation::ms_nVerbose = 1;
    if( nMaxScan == -1 )
    {
        DTREK_ERROR(2,"ERROR:  Could not find any scans!\n");
    }
    /////////////////////////////////////////////////////////////////////////////

    // Parse any additional command line options
  
  for (i = argc_start; i < argc; i++)
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
      else if ("-maxtime" == sTemp)
      {
          if ((i + 1 < argc) && (1==sscanf(argv[i+1],"%d",&nTemp))) {
              m_poStrategy->vSetMaxTime(nTemp);
              i++;
          } else
              DTREK_ERROR(5, sMissingOption);
      }
      else if ("-out" == sTemp)
      {
          i++;
          if (i < argc)
          {
              sHeaderOut = (const char *)argv[i];
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-cmin" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
              {
                  m_poStrategy->vSetCompleteness(fTemp);
              } 
              else
                  DTREK_ERROR(6, sInvalidArgument);

              
              if( i + 1 < argc && 1 == sscanf(argv[i+1],"%f",&fTemp))
              {
                  m_poStrategy->vSetCompletenessTolerance(fTemp);
                  i++;
              }
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-rmin" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
                  m_poStrategy->vSetRedundancy(fTemp);
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if( "-zeroomegawithphioffset" == sTemp )
      {
          m_poStrategy->m_bZeroOmegaWithPhiOffset = TRUE;
          
          if ((i+1 < argc) && (1==sscanf(argv[i+1],"%f",&fTemp)))
          {
              m_poStrategy->m_fZeroOmega = fTemp;
              
              i++;
          }
      }
      else if("-rangemax" == sTemp)
      {
          i++;
          
          if((i+1 < argc) && (1==sscanf(argv[i+1], "%d",&j)) && (1==sscanf(argv[i],"%f",&fTemp))) 
          {
              i++;
          } 
          else if ((i < argc) && (1==sscanf(argv[i], "%f",&fTemp))) 
          {
              j = -1;
          } 
          else
              DTREK_ERROR(5, sMissingOption);
          
          // Go through the scans, and set the rotation.

          for( k = nMinScan; k<=nMaxScan; k++) 
          {
              if( j == -1 || j == k ) // if j scan has been inputted by the user, set rotation max range only to scan nScan. Else, set it to ALL scans.
              {
                  Crotation     oRotation(*m_poHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, k));
                  
                  oRotation.vSetRotMinMaxRange(oRotation.fGetRotMin(),oRotation.fGetRotMax(),fTemp);
                  
                  oRotation.nUpdateHeader(m_poHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, k));
              }
          }
      }
      else if( "-totalrotmax" == sTemp )
      {
          i++;
          
          if( i < argc && 1 == sscanf(argv[i], "%f", &fTemp) && fTemp > 0.0 )
              m_poStrategy->vSetMaximumTotalRotation(fTemp);
          else
          {
              DTREK_ERROR(5,sMissingOption);
          }
      }
      else if( "-totalscansmax" == sTemp )
      {
          i++;
          
          if( i < argc && 1 == sscanf(argv[i], "%d", &nTemp) && nTemp > 0 ) // The maximum number of scans must be at least one
              m_poStrategy->vSetMaximumTotalScans(nTemp);
          else
              DTREK_ERROR(6, sInvalidArgument);
      }
      else if (("-reso" == sTemp) || ("-distance" == sTemp) || ("-swing2theta" == sTemp))
      {
          
          // RB 2/16/05 Apparently Thad is skipping here -reso, -distance and -swing2theta options,
          // because they have been parsed elsewhere. If there is a dash, it might be the next option, unless it's a negative number. 
          while ((i+1 < argc) && ((argv[i+1][0]!='-') || (strspn(argv[i+1],"-+0123456789.e")==strlen(argv[i+1]))))
              i++;

      }
      else if ("-rot" == sTemp)
      {
          double        fRotRangeMax = 0;
          double        fRotRangeMin = 0;
          int           nScan = -1;
          
          i++;
          
          for (j=0; (j < 3) && (i < argc); j++, i++)
          {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (nStat == 1) {
                  if (j==0)
                      fRotRangeMin = fTemp;
                  else if (j==1)
                      fRotRangeMax = fTemp;
                  else if (j==2)
                      nScan = ((int) fTemp);
              } else if (j<2)
                  DTREK_ERROR(6, sInvalidArgument)
              else
                  i--;
          }
          if (j != 3)
          {
              DTREK_ERROR(5, sMissingOption);
          }
          
          // Go through all scans, and set the rotation.

          for(k = nMinScan; k <= nMaxScan; k++)
          {
              if( nScan == -1  ||  nScan == k )  // if nScan has been inputted by the user, set rotation only to scan nScan. Else, set it to ALL scans.
              {
                  Crotation     oRotation(*m_poHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, k));
                  
                  oRotation.vSetRotMinMaxRange(fRotRangeMin, fRotRangeMax, oRotation.fGetRotRange());

                  oRotation.nUpdateHeader(m_poHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, k));
              }
          }
          
          i--;
      }
      else if ("-verbose" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
                  m_poStrategy->vSetVerboseLevel(nTemp);
              else
                  DTREK_ERROR(6, sInvalidArgument);
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-scale" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%f", &fTemp);
              if (1 == nStat)
              {
                  m_fLengthScaleFactor = fTemp;
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
      else if ("-bins" == sTemp)
      {
          i++;
          if (i < argc)
          {
              nStat = sscanf(argv[i], "%d", &nTemp);
              if (1 == nStat)
              {
                  //m_poStrategy->vSetNumResoBins(nTemp);
                  nNumResoBins = nTemp;
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
      else if ("-choose" == sTemp) 
      {
          int nFirstSecondFixed;
          int a2nScans[2];
          i++;
          if ((i < argc) && (argv[i] == (Cstring) "first")) {
              nFirstSecondFixed = 0;
          } else if ((i < argc) && (argv[i] == (Cstring) "second")) {
              nFirstSecondFixed = 1;
          } else if ((i < argc) && (argv[i] == (Cstring) "fixed")) {
              nFirstSecondFixed = 2;
          } else {
              DTREK_ERROR(6,sInvalidArgument);
          };
          i++;
          if ((i<argc) && (1 == sscanf(argv[i],"%d",&a2nScans[0]))) {
              if ((i+1<argc) && (1 == sscanf(argv[i+1],"%d",&a2nScans[1]))) 
                  i++;
              else
                  a2nScans[1] = a2nScans[0];
          } else
              DTREK_ERROR(6,sInvalidArgument);
          if (nFirstSecondFixed == 0) 
              m_poStrategy->vSetFirstScan(a2nScans[0],a2nScans[1]);
          else if (nFirstSecondFixed == 1)
              m_poStrategy->vSetSecondScan(a2nScans[0],a2nScans[1]);
          else 
              m_poStrategy->vAddFixedScan(a2nScans[0],a2nScans[1]);
      }
      else if ("-predicted" == sTemp)
      {
          i++;
          if (i < argc)
          {
              sReflnlistPredicted = (const char *)argv[i];
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-pad" == sTemp)
      {
          i++;
          if ((i < argc) && (1 == sscanf(argv[i], "%f", &fTemp)) && (fTemp>0.0))
              m_poStrategy->vSetPad(fTemp);
          else
              DTREK_ERROR(5,sMissingOption);
      }
      else if ("-prev" == sTemp)
      {
          i++;
          if (i < argc)
          {
              sReflnlistPrevious = (const char *)argv[i];
          }
          else
          {
              DTREK_ERROR(5, sMissingOption);
          }
      }
      else if ("-anom" == sTemp)
      {
          m_poStrategy->vSetAnom(TRUE);
          m_bAnom = true;
      }
      else if ( "-fixrotstart" == sTemp )
      {
          m_poStrategy->vSetFixRotStart(true);
      }
      else
      {
          printf("Unrecognized option '%s'\n",sTemp.string());
          DTREK_ERROR(DTMAIN_ERROR_UNRECOGNIZED_COMMAND_OPTION, "");
      }

    }

    /////////////////////////////////////////////////////////
    // double-check that the strategy object is set correctly
    if( !m_poStrategy->bIsValid(sError) )
    {
        DTREK_ERROR(6, sError);        
    }
    /////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // If the resolution object already exists AND the caller passed a number of bins AND 
  // it is different from what we have, re-generate reso bins.
  if( m_poResoBins && nNumResoBins > 0 && nNumResoBins != m_poResoBins->nGetNumberOfBins() )
  {
       m_poResoBins->vSet(m_poResoBins->dGetMinReso(), m_poResoBins->dGetMaxReso(), nNumResoBins); 
       m_poStrategy->vSetExternalResoBinsPtr(m_poResoBins);
  }
  else if( nNumResoBins > 0 )
      m_poStrategy->vSetNumResoBins(nNumResoBins);
  /////////////////////////////////////////////////////////////////////////////////////////////////


  //MAKE REFLN LIST FOR PREDICTIONS
  vCreatePredictReflnList();
  vExpandPredictReflnList();
  
  if( "" == sReflnlistPredicted )
  {
      // We need to do a prediction, so ...
      nStat = 1;
      
      ///////////////////////////////////////////////
      int           nMaxRefsInPredict = 0;
      if (m_fLengthScaleFactor > 100)
      {
          nMaxRefsInPredict = (int) m_fLengthScaleFactor;
          m_fLengthScaleFactor = -1.0;
      }
      else
          nMaxRefsInPredict = 5000;
      ///////////////////////////////////////////////

      m_poStrategy->vSetAxisInfo(m_poHeader);
      int      nNumAxes = m_poStrategy->nGetNumAxes();
      
      ///////////////////////////////////////////////////////////////////////////////
      // Assume each scan has the same detector position and save it
      m_poStrategy->vSaveDetectorPosition(m_poHeader, f2ThetaSwing, fDistance);
      ///////////////////////////////////////////////////////////////////////////////
      
      // Now loop over all scans and do a predict.
      int           nScan = 0;

      // Keep track of the min and max resolution from predict, so we could set it to strategy later
      double        dMinResoFromPredictObject = -1.0;
      double        dMaxResoFromPredictObject = -1.0;

      for (nScan = 0; nScan <= nMaxScan; nScan++)
      {
          m_poStrategy->vAddStrategyScan(m_poHeader, nScan, nNumAxes);

          printf("\nPredicting reflections for scan %sS%d\n\n", D_K_StrategyTestPrefix, nScan);
          
          // If the scale length is < 0.0, then we need to automatically deterine the multiplier.
          // This is done by observing how many spots get predicted per rotation width,
          // and then restricting accordingly.
          if( m_fLengthScaleFactor < 0.0 )
          {
              // RB 9/14/05 Originally Thad used the main Cpredict object to figure out the 
              // automatic scale factor. The problem with that approach is that when a nPredict function is
              // called on that object, it changes the object state, specifically it applies a scale factor to
              // the cell lengths. To fix the problem, I now create a temporary Cpredict object. 
              
              double        a2fResolutionTemp[2] = {a2fResolution[0], a2fResolution[1]};   
              Cpredict*     poPredictTemp = pSetupPredictObjectForScan(nScan,
                                                                       f2ThetaSwing,
                                                                       fDistance,
                                                                       a2fResolutionTemp);

              // Must figure out how many reflections per degree we will get.
              // Predict 5 degrees worth of data.
              double        fTestRotRange = 5.0;
              double        fScale = 0.0;

              m_fLengthScaleFactor = poPredictTemp->m_ppoCrystal[0]->fCalcVolume();
              m_fLengthScaleFactor = min(1.0, pow(30 * 30 * 30 / m_fLengthScaleFactor, 0.333333));

              m_poPredictReflnList->nDelete();
              
              Crotation&        oRotationTemp = *(poPredictTemp->m_poRotation);
              oRotationTemp.nSetRotRange(0.0, fTestRotRange);

              poPredictTemp->vSetScaleLength(m_fLengthScaleFactor);
              
              nStat = poPredictTemp->nPredict(NULL, m_poPredictReflnList);

              if( nStat ) 
              {
                  printf("ERROR:  Failure during initial prediction.\n");
                  goto end;
              }

              f0 = 360.0 * m_poPredictReflnList->nGetNumReflns() / fTestRotRange;
              
              f0 = pow(f0 / nMaxRefsInPredict, 0.333333);
              
              m_poPredictReflnList->vDeleteAll();             
              
              if( f0 > 1.0)
                  m_fLengthScaleFactor *= 1.0 / f0;
              else
                  m_fLengthScaleFactor *= 1.0;
              
              m_fLengthScaleFactor = dRoundOff(m_fLengthScaleFactor, 2);

              printf("\n\nAutomatically setting cell multiplier (scale factor) to %.2f\n\n", m_fLengthScaleFactor);
              
              delete poPredictTemp;
              poPredictTemp = NULL;
          }
          
          ////////////////////////////////////////////////////////////////////////////////////
          poPredict = pSetupPredictObjectForScan(nScan,
                                                 f2ThetaSwing,
                                                 fDistance,
                                                 a2fResolution);
         
          ////////////////////////////////////////////////////////////////////////////////////////
          // Keep track of the reso limits from predict, so we could set it to strategy
          if( nScan == 0 )   // initialize min and max resolution from the very first scan
          {
                dMinResoFromPredictObject = poPredict->dGetResolution(0);
                dMaxResoFromPredictObject = poPredict->dGetResolution(1);
          }
          else
          {
                if( dMinResoFromPredictObject < poPredict->dGetResolution(0) )
                    dMinResoFromPredictObject = poPredict->dGetResolution(0);

                if( dMaxResoFromPredictObject > poPredict->dGetResolution(1) )
                    dMaxResoFromPredictObject = poPredict->dGetResolution(1);
          }
          ///////////////////////////////////////////////////////////////////////////////////////

          m_poStrategy->vSetCellScale(m_fLengthScaleFactor);
          
          // RB this is taken care of in pSetupPredictObjectForScan()
          //poPredict->vSetScaleLength(m_fLengthScaleFactor);

          Crotation&        oRotation = *(poPredict->m_poRotation);
          if( m_poStrategy->bScanUsed(nScan) ) 
              oRotation.nSetRotRange(oRotation.fGetRotMin(),oRotation.fGetRotMax());
          
          if( m_wVerbose & DTREK_DTMAINSTRAT_VERBOSE_PREDICT )
          {
              poPredict->m_poCrysGonio->nList();
              poPredict->m_poRotation->nList(10,poPredict->m_poCrysGonio);
              poPredict->m_ppoDetector[0]->m_poGoniometer->nList();
          }
          ///////////////////////////////////////////////////////////////////////////
          
          ////////////////////////////////////////////////////////////////////////////////
          // NOW DO THE PREDICTION
          i = m_poPredictReflnList->nGetNumReflns();  // memorize the current number of reflections
          
          nStat = poPredict->nPredict(NULL, m_poPredictReflnList);
          
          if (nStat)
          {
              printf("ERROR:  Failure during prediction of scan #%d.\n",nScan);
              goto end;
          }
          
          // Set the scan field.
          for (;i<m_poPredictReflnList->nGetNumReflns();i++)
          {
              (*m_poPredictReflnList)[i].vSetIntensity((float) nScan);
          }
          
          m_poStrategy->vRememberScanGonioSettings(poPredict->m_poCrysGonio,poPredict->m_poRotation,nScan);
          delete poPredict;
          poPredict = NULL;
      }
      
      printf(" ... done.  ");
      printf("%d reflections predicted.\n",m_poPredictReflnList->nGetNumReflns());
  
      ///////////////////////////////////////////////////////////////////////////////////////////////////////
      //RB: Reset strategy resolution to match resolution, found by predict
      if( !m_bExternalResoBinsPtr )
      {
          printf("\nSetting strategy resolution from predicted reflection list: %.2f %.2f A\n", 
                    dMinResoFromPredictObject, dMaxResoFromPredictObject);
          
          if( m_poResoBins )
            vDeleteResoBins();

          m_poResoBins = new CResoStrategyBins(dMinResoFromPredictObject, dMaxResoFromPredictObject, nNumResoBins);
          m_poStrategy->vSetExternalResoBinsPtr(m_poResoBins);
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////
  }
  else
  {
      m_poPredictReflnList->nRead(sReflnlistPredicted);
      
      // Set the scan field.
      for(i = 0; i < m_poPredictReflnList->nGetNumReflns(); i++)
      {
          (*m_poPredictReflnList)[i].vSetIntensity((float)0);
      }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if( "" != sReflnlistPrevious )
  {
      // Read in the previous reflection list.
      // Make its reflections into a new scan.
      Creflnlist*       poReflnlistPrev = new Creflnlist();
      
      m_poStrategy->nExpand(*poReflnlistPrev);
      
      poReflnlistPrev->nRead(sReflnlistPrevious);
      
      i = m_poPredictReflnList->nGetNumReflns();
      
      j = (int) (*m_poPredictReflnList)[i - 1].fGetIntensity();  // RB: Apparently to get the last scan number?
      
      j++;
      
      m_poPredictReflnList->nInsertListFrom(*poReflnlistPrev);
      
      for( ; i < m_poPredictReflnList->nGetNumReflns(); i++) 
      {
          (*m_poPredictReflnList)[i].vSetIntensity((float) j);
      }
      
      delete poReflnlistPrev;
      
      // Add this as a fixed scan.
      m_poStrategy->vAddFixedScan(j, j);
  }

  // Must get any resolution changes from poPredict BEFORE calling ::nSetup!

  nStat = m_poStrategy->nSetup(*m_poHeader, *m_poPredictReflnList);
  
  if( !nStat && bDoList )
      m_poStrategy->nList(DTREK_CSTRATEGY_LIST_TITLE |
                          DTREK_CSTRATEGY_LIST_RESO  |
                          DTREK_CSTRATEGY_LIST_REFLN_LIST);

  // Seed the random number generator so that we can have duplicatable results.
  srand(0);


  m_poStrategy->nUnPackRedundancy();

  
  printf("\n\n-----------------------------------------------------------------------\n\n");
  
  if (!nStat)
      nStat = m_poStrategy->nCalcScans();


  if( nStat == 1 )
  {
      printf("WARNING:  Could not find a scan that would meet requested conditions!\n");
  }
  printf("\n\n-----------------------------------------------------------------------\n\n");
  
  if( 0== nStat || 2 == nStat )
  {
      if( !m_poStrategy->nCalcCompletenessSensitivity(m_poStrategy->m_oBestRange, true) )
      {
      }
  }

  if( 0 == nStat || 2 == nStat ) // 2 == nStat is a special case that could happen when CrystalClear is using dtstrategy, but not dtmultistrategy.
  {
    char*         pcLine = "--------------------------------------------------------------------------------\n";
    printf("\nStrategy results\n");
    (m_poStrategy->m_oBestRange).nPrint(*m_poStrategy);
  }
  
    ///////////////////////////////////////////////////////////////////////////
    //if( 0 == nStat )  // RB: trying to still get the listing, even though strategy failed to achieve objective
    if( !m_bExternalResoBinsPtr )
        wListComplMsgCtrl = DTREK_STRATEGY_LC_MSG_LOW_COMPLETENESS; // If the reso bins object IS external, we do not want to warn 
                                                                    // users about low completeness.  The reason: in dtmultistrategy
                                                                    // the TOTAL reso range may be larger than a reso range for a SINGLE
                                                                    // detector position, so the completeness, based on the total
                                                                    // reso range for a single position may be lower than the
                                                                    // target completness because of that. But that is ok.

    if( 0 == nStat || 1 == nStat ) 
    {
      m_poStrategy->nListCompleteness(&(m_poStrategy->m_oBestRange), wListComplMsgCtrl);
    }
    else if( 2 == nStat )  // 2 == nStat is a special case that could happen when CrystalClear is using dtstrategy, but not dtmultistrategy
    {
      m_poStrategy->nListCompleteness(&(m_poStrategy->m_oBestRange), wListComplMsgCtrl | DTREK_STRATEGY_LC_MSG_USER_ACCEPTED);
    }
    ///////////////////////////////////////////////////////////////////////////

    if( 0 == nStat  && "" != sReflnFile )
    {
        m_poStrategy->bLeaveBestRangeReflectionsOnly();
        nStat = m_poStrategy->nWrite(sReflnFile);
    }
    // If reflection list ptr was set from outside, must assume that fixed scans are not needed, 
    // So pass false to bLeaveBestRangeReflectionsOnly().
    // For a more detailed explanation why we need to pass false please see a comment above a call 
    // to CDTMainMultiStrategy::bAddPreviousReflections(...) in file DTMainMultiStrategy.cpp
    else if( m_bExternalPredictReflnListPtr ) 
        m_poStrategy->bLeaveBestRangeReflectionsOnly(false);

    // Print SINGLE scan solution here, but only if it's NOT the actual solution, i.e. the actual solution is MULTIPLE scan.
    // The reason: if SINGLE is the actual solution, we may need to adjust the phi and omega first.
    if( m_poStrategy->bIsMultipleScanSolution() )
        m_poStrategy->vPrintSummary(true);
    
    //////////////////////////////////////////////////////////////////////////////
    // SEE IF WE NEED TO ADJUST PHI AND OMEGA
    if( m_poStrategy->m_bZeroOmegaWithPhiOffset )
        m_poStrategy->nAdjustPhiToZero(*m_poHeader, m_poStrategy->m_oBestRange);
    //////////////////////////////////////////////////////////////////////////////

    // Print SINGLE scan solution here, but only if it IS the actual solution, i.e. not a MULTIPLE scan.
    if( !m_poStrategy->bIsMultipleScanSolution() )
        m_poStrategy->vPrintSummary(true);
    else  // print a multiple scan solution
        m_poStrategy->vPrintSummary(false);

    if( 0 == nStat || 2 == nStat ) // 2 == nStat is a special case that could happen when CrystalClear is using dtstrategy, but not dtmultistrategy.
    {
      nStat = m_poStrategy->m_oBestRange.nUpdateHeader(*m_poHeader, *m_poStrategy);
      if( 0 == nStat && !m_bExternalHeader )  // Don't write out the header if dtmainstrategy is called from elsewhere. The caller will write out the header.
        nStat = m_poHeader->nWrite(sHeaderOut);
    }    
end:

  if (NULL != poPredict)
    delete poPredict;
  
  return nStat;
}


void CDTMainStrategy::vError(const int nErrorNum, const Cstring& sMessage)
{
  if( m_bExternalHeader )
      return; // since dtstrategy was invoked by another module, that module is in charge of the error info
    
  if (0 != nErrorNum)
    {
      printf("%s\n",sMessage.string());
    }
  printf("\ndtstrategy calculates the starting rotation angle and range for a\n"
         "diffraction experiment that gives a desired completeness.\n"
         "It writes *dtstrategy.head.\n");
  
  printf("\ndtstrategy - Usage:\n"
         "dtstrategy header_file [options...]\n\n"
         "Command line options: Description:\n\n");
        
       vPrintOptionHelp();

  printf("Examples:\n"
         " dtstrategy dtrefine.head\n"
         " dtstrategy dtrefine.head -rot 0 180 -reso 80 2\n\n");

  fflush(stdout);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vPrintOptionHelp()
{
    printf(" header_file        Name of a file with a d*TREK header used to\n"
           "                    initialize the strategy procedure.\n\n"
           " -pad fPad\n"
           "                    Specify a padding term.  The padding determines a small\n"
           "                    no-mans-land of fPad degrees at the start and end of the\n"
           "                    rotation range. Reflections in the no-mans-land are not\n"
           "                    counted when computing the completeness or redundancy.\n\n"
           " -rot  fStart fEnd  [nScan]\n"
           "                    Specifies the min and max rotation angle start and\n"
           "                    end in degrees to test for completeness. This will apply to\n"
           "                    all scans found, unless nScan is set.\n"
           "                    Default: 0 360 for all scans.\n\n"
           " -fixrotstart       Fix rotation angle start at the input value.\n\n"
           "                    Default: not fixed.\n\n"
           " -rangemax fRangeMax [nScan]\n"
           "                    Maximum rotation range allowed per scan. This will apply to\n"
           "                    all scans found, unless nScan is set.\n"
           "                    Default: obtained from -rot spec. for each scan.\n\n"
           " -totalrotmax fMaxRot\n"
           "                    Total maximum rotation width in a multiple scan solution.\n"
           "                    Default: allow any width.\n\n"
           " -totalscansmax nMaxScans\n"
           "                    Total maximum number of scans in a multiple scan solution.\n"
           "                    NOTE: This value should be at least 2. Default: Allow any\n"
           "                    number of scans in a multiple scan solution.\n\n"
           //" -cmin fCompletenessMin [fMinIncCompPerDegree]\n"
           //"                    Minimum percentage completeness for rotation range\n"
           //"                    solution(s). If fMinIncCompPerDegree is specified,\n"
           //"                    then dtstrategy will not increase the total rotation width\n"
           //"                    beyond the point where each added degree adds less than\n"
           //"                    fMinIncCompPerDegree to the total completeness.\n"
           //"                    Default: fCompletenessMin = 99 fMinIncCompPerDegree = 0\n\n"
           " -cmin fCompletenessMin [fCompletenessTolerance]\n"
           "                    Minimum percentage completeness for rotation range\n"
           "                    solution(s). If fCompletenessTolerance is specified,\n"
           "                    then dtstrategy will explore solutions within a range\n"
           "                    of completeness: fCompletenessMin +/- fCompletenessTolerance.\n"
           "                    Default: fCompletenessMin = 99; fCompletenessTolerance = 0.\n\n"
           " -rmin fRedundancyMin\n"
           "                    Minimum allowed redundancy for rotation range\n"
           "                    solution(s). If fRedundancyMin==0.0, then we\n"
           "                    assume the user wants to collect an entire\n"
           "                    hemisphere.  Default: Not used.\n\n"
           " -choose [first | second | fixed] nFirstScan [nLastScan]\n"
           "                    Specify which scan (or range of scans) to choose from\n"
           "                    for strategy.  By default, we choose from all available\n"
           "                    scans.\n\n"
           " -prev sPrevReflnlistFilename\n"
           "                    Name of a reflnlist file with previously available\n"
           "                    reflections. Default: no previous list.\n\n"
           " -predicted sPredictReflnlistFilename\n"
           "                    Name of a reflnlist file with previously\n"
           "                    predicted reflections. If this is specified no\n"
           "                    prediction is done.  Default: predict reflections.\n\n"
           " -ref sReflnlist    Write the selected reflections to the file sReflnlist.\n\n"
           " -scale fScale      Factor to multiple cell lengths by, OR if >100, maximum\n"
           "                    number of reflections to use per scan.\n\n"
           " -anom              Assume I+ != I-\n\n"
           "                    Default: Automatically determined.\n\n"
           " -maxtime fTime     Maximum number of seconds to take in strategy calculation\n"
           "                    This does not include overhead time which might\n"
           "                    be non-negligible\n\n"
           " -out sHeaderOut    Name of output file with resultant scan rotation info.\n"
           "                    Default: *dtstrategy.head.\n\n"
           " -help              Print this help text.\n\n");

           if (NULL != CExperimentalStrategy::ms_pcHelp)
           {
                printf("%s", CExperimentalStrategy::ms_pcHelp); // Prints out additional help.
           }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vCreateStrategyObject()
{
    vDeleteStrategyObject();

    if( !m_poHeader )
    {
        printf("\nERROR: header object has NULL pointer.\n");
        return;
    }

    m_poStrategy = new Cstrategy(m_poHeader);

    m_poStrategy->nList(DTREK_CSTRATEGY_LIST_CRYSTAL);

    m_poStrategy->vSetScanInfoObj(m_poExternScanInfo);
}
////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vDeleteStrategyObject()
{
    if( m_poStrategy )
    {
        delete m_poStrategy;
        m_poStrategy = NULL;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vCreatePredictReflnList()
{
    if( NULL == m_poPredictReflnList )
        m_poPredictReflnList = new Creflnlist();
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vExpandPredictReflnList()
{
    m_poStrategy->nExpand(*m_poPredictReflnList);
}
///////////////////////////////////////////////////////////////////////////////////////////////// 
void CDTMainStrategy::vDeletePredictReflnList()
{
    if( !m_bExternalPredictReflnListPtr && m_poPredictReflnList )
    {
        delete m_poPredictReflnList;
        m_poPredictReflnList = NULL;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vDeleteResoBins()
{
    if( !m_bExternalResoBinsPtr && m_poResoBins )
    {
        delete m_poResoBins;
        m_poResoBins = NULL;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vSetExternalPredictReflnListPtr(Creflnlist* pList)
{
    vDeletePredictReflnList();

    if( pList )
    {
        m_poPredictReflnList = pList;
        m_bExternalPredictReflnListPtr = true;
    }
    else
    {
        m_poPredictReflnList = NULL;
        m_bExternalPredictReflnListPtr = false;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vSetExternalResoBinsPtr(CResoStrategyBins* pBins)
{
    vDeleteResoBins();
    
    if( pBins )
    {
        m_poResoBins = pBins;
        m_bExternalResoBinsPtr = true;
    }
    else
    {
        m_poResoBins = NULL;
        m_bExternalResoBinsPtr = false;
    }

    if( m_poStrategy )
        m_poStrategy->vSetExternalResoBinsPtr(pBins);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vSetExternalScanInfoPtr(CstrategyScans* pScans)
{
    m_poExternScanInfo = pScans;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vDoExternalResoBinsStatistics(CResoStrategyBins* poResoBins, Creflnlist* poReflnList)
{
    if( !poResoBins || !poReflnList )
        return; // nothing to do

    if( m_fLengthScaleFactor < 0 )
    {
        printf("\nERROR: cannot do resolution bin statistics.\n");
        return; 
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // DO THEORETICAL REFLECTION COUNT 
    Ccrystal    oCrystal(*m_poHeader);
    if( !oCrystal.bIsAvailable() )
    {
        printf("\nERROR: cannot do resolution bin statistics. Crystal information is not available.\n");
        return;
    }

    double      dVolumeScale = 0.0;
    oCrystal.bLinearScaleToVolumeScale(m_fLengthScaleFactor, dVolumeScale);
    
    // Generate theoretical reso bin count

    // Now, calculate theoretical count
    poResoBins->vCountTheoreticalUniqueReflectionsInResoShells(&oCrystal, m_bAnom);

    // If the crystal unit cell is scaled, correct the theoretical reflection count. 
    if( 1.0 != dVolumeScale )
        poResoBins->vScaleReflnCount(dVolumeScale, DTREK_RESOBIN_SCALE_THEORETICAL_COUNT);
    //////////////////////////////////////////////////////////////////////////////////////////////

    if( m_poStrategy )
    {
        m_poStrategy->vSetExternalResoBinsPtr(poResoBins);
        m_poStrategy->vSetExternalWorkingReflnListPtr(poReflnList);

        m_poStrategy->vFillStats(NULL, false); // We pass NULL, because all reflections in the list should be considered.
                                               // We pass false, because the list may actually contain reflections outside 
                                               // of the resolution range of poResoBins.
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vDoExternalReflnListCalcSetDStarCubed(Creflnlist* poReflnList, 
                                                            double& dStarCubedMin, 
                                                            double& dStarCubedMax)
{
    if( !poReflnList  || !m_poStrategy )
        return; // nothing to do

    m_poStrategy->vSetExternalWorkingReflnListPtr(poReflnList);
    
    m_poStrategy->vReflnListCalcSetDStarCubed(dStarCubedMin, dStarCubedMax);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainStrategy::vDoExternalReflnListSortReduce(Creflnlist* poReflnList)
{
    if( !poReflnList  || !m_poStrategy )
        return; // nothing to do

    m_poStrategy->vSetExternalWorkingReflnListPtr(poReflnList);
    
    m_poStrategy->vReflnListSortReduce();
}
///////////////////////////////////////////////////////////////////////////////////////////////////
Cpredict* CDTMainStrategy::pSetupPredictObjectForScan(int nScan,
                                                      double f2ThetaSwing,
                                                      double fDistance,
                                                      double* a2fResolution)
{
    double        a6fCell[6] = {0.0};
    double        a6fCell_save[6] = {0.0};
    bool          bScale = false;
    
    if( m_fLengthScaleFactor < 1.0 && m_fLengthScaleFactor > 0.0 )  // sanity check
    {
        // Temporarily modify all cell lengths in the header object by the scale factor
        bScale = true;
        m_poHeader->nGetValue(Cstring("CRYSTAL_UNIT_CELL"), 6, a6fCell);
        for(int ii=0; ii < 6; ii++)
        {
            a6fCell_save[ii] = a6fCell[ii];
            
            if( ii < 3 )
                a6fCell[ii] *= m_fLengthScaleFactor;
        }

        m_poHeader->nReplaceValue(Cstring("CRYSTAL_UNIT_CELL"), 6, a6fCell);
    }
    
    
    Cpredict*     poPredict = new Cpredict(*m_poHeader, NULL, nScan, D_K_StrategyTestPrefix);

    if( m_wVerbose & DTREK_DTMAINSTRAT_VERBOSE_PREDICT )
        poPredict->vSetVerbose(DTREK_PREDICT_VERBOSE_MOSAICITY);
    else
        poPredict->vSetVerbose(0U);

    // Restore the cell lenghts in the header object
    if( bScale)
    {
        m_poHeader->nReplaceValue(Cstring("CRYSTAL_UNIT_CELL"), 6, a6fCell_save);
    }
    
    // Let poPredict calculate some resolutions available on the detector

    poPredict->nSetupSourceCrystal();

    if( f2ThetaSwing > cdStratDetUnrealisticValue )
      poPredict->m_ppoDetector[0]->m_poGoniometer->nSetSwing(f2ThetaSwing);

    if( fDistance > cdStratDetUnrealisticValue )
      poPredict->m_ppoDetector[0]->m_poGoniometer->nSetDistance(fDistance);

    if( f2ThetaSwing > cdStratDetUnrealisticValue || fDistance > cdStratDetUnrealisticValue )
      poPredict->m_ppoDetector[0]->nUpdateHeader(m_poHeader,"");

    poPredict->nSetupDetector();

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////          
    // RB: The following piece of code, using a2fResolution[] will have to be revised later,
    // because we want to use the resolution object exclusively
    if( 0.0 > a2fResolution[1] )
    {
      // This means go to the edge of the detector, not the corner
      poPredict->vUseEdgeReso();
  
      poPredict->vSetResolution(-1.0, -1.0);
    }

    if( 0.0 != a2fResolution[0]  && 0.0 != a2fResolution[1] )
    {
      // Restrict the Cpredict set resolution by the 
      // detector maximum resolution calculated in Cpredict()
  
      a2fResolution[1] = max(poPredict->dGetResolution(1), a2fResolution[1]);

      poPredict->vSetResolution(a2fResolution[0], a2fResolution[1]);
    }
    else if( 0.0 == a2fResolution[0] && 0.0 == a2fResolution[1] )
    {
      poPredict->vSetResolution(a2fResolution[0], a2fResolution[1]);
    }
    //////////////////////////////////////////////////////////////////////////////
    // Now that the predict resolution has been possibly set,
    // retrieve those values to use for m_poStrategy

    poPredict->nSetupDetector();
    
    if( m_fUserInputMosaicity >= 0.0 )
    {
        poPredict->vSetCrysMosaicity(m_fUserInputMosaicity);
        
        printf("Changing mosaicity for prediction to a user-specified value of %.2f deg\n", m_fUserInputMosaicity);
    }
    
    return poPredict;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
double CDTMainStrategy::fGetBestCompleteness()
{
    if( !m_poStrategy )
        return 0.0;

    return m_poStrategy->m_oBestRange.m_fCompleteness;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
double CDTMainStrategy::fGetCompletenessTolerance()
{
    if( !m_poStrategy )
        return 0.0;

    return m_poStrategy->fGetCompletenessTolerance();
}
////////////////////////////////////////////////////////////////////////////////////////////////////

