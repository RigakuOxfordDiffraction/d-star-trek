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
// DTMainMultiStrategy.cpp       Initial author: RB          21-Mar-2005
// This file contains the member functions of class CDTMainMultiStrategy

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

#include "DTMainMultiStrategy.h"
#include "dtrekvec.h"
#include "Cdetector.h"
#include "Creflnlist.h"
#include "Cpredict.h"

#include "ResoBins.h"

#include "Segment.h"
#include "Crotation.h"

#include "ImageExpTimeEvaluator.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////
const char*     c_pcScale                   = "-scale";  
const char*     c_pcAnom                    = "-anom";  
const char*     c_pcAuto2Theta              = "-auto2theta";  
const char*     c_pcSet2Theta               = "-set2theta";  
const char*     c_pcAutoDistance            = "-autodist";  
const char*     c_pc2ThetaMax               = "-2thetamax";  
const char*     c_pcSwing2Theta             = "-swing2theta";  
const char*     c_pcFirstAuto2Theta         = "-firstauto2theta";  
const char*     c_pcSetDistance             = "-distance";  
const char*     c_pcDistanceLimits          = "-distlimits";
const char*     c_pcRotationLimits          = "-rot";
const char*     c_pcReso                    = "-reso";  
const char*     c_pcPrev                    = "-prev";
const char*     c_pcHeaderDistSearchKey     = "ist";  
const char*     c_pcHeader2ThetaSearchKey   = "2Th";  
const char*     c_pcCrystalGonio            = "-crysgonio";
const char*     c_pcMosaicity               = "-mosaicity";
const char*     c_pcCheckOverlap            = "-overlap";
const char*     c_pcOverlapList             = "-overlaplist";
const char*     c_pcCompleteness            = "-cmin";
const char*     c_pcOldHeaderScanFormat     = "-oldheaderscanformat";
const char*     c_pcAutoExposure            = "-autoexposure";
const char*     c_pcImageRotWidth           = "-rotimagewidth";
////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct _tagOVERLAP_RESULTS
{
    int     nIndex;                               // scan index
    double  dResoMin;
    double  dResoMax;
    std::vector<double>     afOverlapsPercentage; // to store overlap percentages for each rotation width tested for each scan
}OVERLAP_RESULTS;

const int c_nNumberOfResoBins = 10; // 2DO: consolidate the bin count. BTW, it could also be passed on the command line! 10 must be the default.
////////////////////////////////////////////////////////////////////////////////////////////////////
CDTMainMultiStrategy::CDTMainMultiStrategy()
{
    m_sError = "";

    m_dUserResoMin   = -1.0;
    m_dUserResoMax   = -1.0;
    
    m_segUser2ThetaLimits.vSet(0.0, 0.0);
    
    m_n2ThetaDirectionSign = 1;

    m_dUser2ThetaMust = cdStratDetUnrealisticValue;
    
    m_dFractionDetectorOverlap = 0.0;

    m_poPredictReflnList = new Creflnlist();
    vExpandReflnList(m_poPredictReflnList);

	m_poPreviousReflnList = NULL;

    m_poResoBins = NULL;

    m_sHeaderOut = sDtrekGetPrefix() + "dtmultistrategy.head";
    m_sReflnListOut = "";

    m_dDistance = -1.0;

    m_d2ThetaFromHeader = cdStratDetUnrealisticValue;
    
    m_segUserDistLimits.vSet(0, 1000); // it is important to initialize it to something wide enough

    m_nUserSpotSepPix = 2;
    
    m_nActive2ThetaPositionIndex = -1;

    m_dAllDetectorPositionsAllScansTotalRotWidth = 0.0;

    m_bScaleSetFromCommandLine = false;

    m_bResolutionSet = false;

    m_pstOverlapCheckParams = NULL;

    m_stSpotSizeInfo.vInit();

    m_fUserInputMosaicity = -1.0;

    m_fUserTargetCompleteness = 99;
    m_fUserCompletnessTolerance = 0;
    
    m_bWriteOverlapReflnList = false;
    
    m_paoSaveTestScans = NULL;
    m_bUseOldHeaderScanFormat = false;
    m_bTestScansGeneratedFromCommandLine = false;
    m_bTestScansGeneratedFromOldHeaderFormat = false;

    m_pstPreviousReflnsInfo = NULL;

    m_fUserImageRotationWidth = -1.0;

    m_fImageAutoExposureTime = -1.0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
CDTMainMultiStrategy::~CDTMainMultiStrategy()
{
    if( m_poPredictReflnList )
    {
        delete m_poPredictReflnList;
        m_poPredictReflnList = NULL;
    }
    
    if( m_poPreviousReflnList )
    {
        delete m_poPreviousReflnList;
        m_poPreviousReflnList = NULL;
    }
    
    if( m_poResoBins )
    {
        delete m_poResoBins;
        m_poResoBins = NULL;
    }

    for(int ii=0; ii < (int)m_aoScanSolutionSet.size(); ii++)
    {
        if( m_aoScanSolutionSet[ii] )
        {
            delete m_aoScanSolutionSet[ii];
            m_aoScanSolutionSet[ii] = NULL;
        }
        
        m_aoScanSolutionSet.clear();
    }
    
    if( m_pstOverlapCheckParams )
    {
        delete m_pstOverlapCheckParams;
        m_pstOverlapCheckParams = NULL;
    }

    if( m_paoSaveTestScans )
    {
        delete m_paoSaveTestScans;
        m_paoSaveTestScans = NULL;
    }

    if( m_pstPreviousReflnsInfo )
    {
        delete m_pstPreviousReflnsInfo;
        m_pstPreviousReflnsInfo = NULL;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
int CDTMainMultiStrategy::nExecute(unsigned int argc, char* argv[])
{
  vDtrekSetModuleName("dtmultistrategy");
  vPrintCopyrightInfo();
  //cout << "\ndtmultistrategy:  Copyright (c) 2006 Rigaku\n";
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
    cout << "\nCommand line:\n" << sGetCommandLine(argc, argv, 71) << endl << endl << flush;


    /////////////////////////////////////////////////////////////////////////////////////////////
    
    vSetCommandArguments(argc, argv);

    if( bIsHelpRequest() )
    {
        DTREK_ERROR(0, "");  // display general help
    }

    if( argc < 2 )
    {
        DTREK_ERROR(1, "No header filename specified\n");
    }

    Cstring     sError("");
    
    // Create a header object
    Cstring     sHeader(argv[1]);  // assuming that the first argument is a header
    if( !bCreateHeader(sHeader) )
    {
        sError = "Invalid header file ";
        sError += Cstring(argv[1]);
        
        DTREK_ERROR(1, sError);
    }
    /////////////////////////////////////////////////////////////////////////////////
    
    vClearMultiStrategyOutputFromHeader();

    ///////////////////////////////////////////////////////////
    sError = "";
    if( !bInsureHeaderRotLimitsConsistency(sError) )
    {
        if( "" != sError )
            DTREK_ERROR(1, sError);
    }

    m_bUseOldHeaderScanFormat = bGetCommandOption(c_pcOldHeaderScanFormat);

    sError = "";
    if( !bSetTestScansFromCommandLineToHeader(sError) )
    {
        if( "" != sError )
            DTREK_ERROR(1, sError);
    
        if( !bInsureTestScansInHeaderAreInNewFormat(sError) )
            DTREK_ERROR(1, sError);  // bInsureTestScansInHeaderAreInNewFormat() must not fail, because there's got to be at least one scan
                                     // which we can then convert to the new format
    }
    
    if( !bApplyRestrictedRotationLimitsToScansInHeader(sError) )
    {
        if( "" != sError )
            DTREK_ERROR(1, sError);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Save test scans, because when we call DTMainStrategy it may update some info, e.g. rotation ranges
    vSaveTestScansFromHeader();

    vSetOutputHeaderPath();
    vSetOutputReflnListPath();

    sError = "";
    
    if( !bParseMosaicityInput(sError) )
    {
        if( "" != sError )
        {
            DTREK_ERROR(1, sError);
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////

    m_oMainStrategy.vSetExternalHeader(m_poHeader);
    // Since we are going to call CDTMainStrategy::nExecute(), we need to change the exe name.
    // This change is needed for the log file only, 
    vSetExecutableNameDTSTRATEGY();

    ///////////////////////////////////////////////////////////////////////////////////////////
    // DETERMINE THE DISTANCE
    bool        bAutoDistance = bSetAutoDistance(sError);
    if( !bAutoDistance || m_dDistance < 0.0 ) 
    {
        if( sError != "" )  // we have to distinguish 2 cases: (a) the user didn't request autodistance and (b) they did, but made an error
        {
            DTREK_ERROR(1, sError);
        }
        else if( !bSetUserInputDistance(sError) ) // see if the user set distance explicitly
        {
            if( sError != "" )
            {
                DTREK_ERROR(1, sError);
            }
            else 
                vGetDetectorPositionFromHeader(DTREK_MULTISTRATEGY_DET_DISTANCE);
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // Now we can try and determine required resolution and 2theta
    //////////////////////////////////////////////////////////////////////////////////////////
    // DETERMINE 2THETA
    bool        bAuto2Theta = bSetAuto2Theta(sError);
    if( !bAuto2Theta || 0 == m_ad2ThetaPositions.size() ) 
    {
        if( sError != "" )
        {
            DTREK_ERROR(1, sError);
        }
        else if( !bSetUserInput2Theta(sError) ) // see if the user set 2Theta(s) explicitly
        {
            if( sError != "" )
            {
                DTREK_ERROR(1, sError);
            }
            else
            {
                vGetDetectorPositionFromHeader(DTREK_MULTISTRATEGY_DET_2THETA); // this should set m_d2ThetaFromHeader
                m_ad2ThetaPositions.clear();
                m_ad2ThetaPositions.push_back(m_d2ThetaFromHeader);
            }
        }
    }
    
    // The resolution could have been already set inside bSetAuto2Theta(), otherwise we need to determine it now.
    if( !m_bResolutionSet ) 
    {
        Cstring     sResoOption("");
        
        if( !bSetResolutionLimitsFromCommandLine(sResoOption, sError) )
        {
            if( "" != sError && !(sResoOption == "corner"  ||
                                  sResoOption == "corners" || 
                                  sResoOption == "edge"    ||
                                  sResoOption == "edges"   ||
                                  sResoOption == "upperEqualsHighest")
                                 )
            {
                DTREK_ERROR(1, sError);
            }

            printf("\nTrying to determine the resolution from the header and/or %s option.\n", c_pcSet2Theta);
            if( !bSetResolutionLimitsFrom2ThetaPositions(sResoOption, sError) )
            {
                if( sError != "" )
                {
                    DTREK_ERROR(1, sError);
                }
            }
            else
                m_bResolutionSet = true;
        }
        else
            m_bResolutionSet = true;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    if( !m_bResolutionSet )
    {
        sError = "Failed to determine the required resolution";
        DTREK_ERROR(1, sError);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    if( !bParseReflnOverlapCheckParams(sError) )
    {
        if( "" != sError )
        {
            DTREK_ERROR(1, sError);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Calculate exposure time per image based on the exposure time analysis results recorded 
    // in the input header by dtranker.
    // One of the assumptions is that every image will have a constant rotation width passed
    // either through the header or on the command line. In the future we may have to 
    // calculate the image rotation width based on the reflection overlap check results. 
    // So the "calculate auto exposure" call may eventually be moved to after the overlap checking.
   
    // Establish an image rotation width common to all scans. In the future we may work
    // on an automatic (variable?) image rotation width
    if( !bGetImageRotationWidth() )
    {
        DTREK_ERROR(1, m_sError);
    }

    if( !bCalculateAutomaticExposureTime() && "" != m_sError )
        printf("\n\nERROR: %s\n\n", m_sError.string()); // for now let's NOT stop the whole thing
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    // Memorize target completeness and tolerance
    vSaveTargetCompletenessTolerance();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // We need to memorize whether the scale has been set. If it hasn't, we will need to ask the strategy object for an
    // automatically-generated scale. The scale is needed to scale reflection count in the resobins object. 
    Cstring     sArgs("");
    if( bGetCommandOption(c_pcScale, sArgs, false) && sArgs.string() && sArgs.GetAt(0) != '\0' )
        m_bScaleSetFromCommandLine = true;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vSetDetectorPositionToHeader(DTREK_MULTISTRATEGY_DET_DISTANCE);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    CstrategyScans*     pCurrentScans = NULL;
    Creflnlist*         poCurrentPredictList = NULL;
    double              dScaleFactor = -1.0;
    double              dResoMin = -1.0;
    double              dResoMax = -1.0;
    int                 nRet = 0;
    bool                bCompletenessAchieved = false;
    int                 nTrialsToAchieveCompleteness = 0;
    double              fRealisticCompleteness = -1.0;
    bool                bRestoreSavedInputScans = false;
    
    for(int ii=0; ii < (int)m_ad2ThetaPositions.size(); ii++)
    {
        // Since the completeness tolerance on the command line might have changed because of the
        // automatic adjustment, we must restore it to its original value.
        if( ii != 0)
            vUpdateTargetCompletenessTolerance(m_fUserTargetCompleteness, m_fUserCompletnessTolerance);

        bCompletenessAchieved = false;
        nTrialsToAchieveCompleteness = 0;
        
        while( !bCompletenessAchieved && nTrialsToAchieveCompleteness < 2 ) 
            // A check nTrialsToAchieveCompleteness < 2 is used just for safety: 
            // we don't want to end up with an infinite loop. 
        {                                                                   
            /////////////////////////////////////////////////////
            // Reset everything...
            m_oMainStrategy.vInit();
            m_oMainStrategy.vSetUserInputMosaicity(m_fUserInputMosaicity);

            if( 0 == ii && 0 == nTrialsToAchieveCompleteness )
                m_oMainStrategy.vSetVerbose(DTREK_DTMAINSTRAT_VERBOSE_PREDICT);
            else
                m_oMainStrategy.vSetVerbose(0U);

            //...but keep the old scale factor.
            // check that the scale factor was, in fact, generated automatically by dtstrategy,
            // when ii was equal to 0, and set the same scale factor for the remaining values of ii
            if( !m_bScaleSetFromCommandLine && ii > 0 && dScaleFactor > 0.0 ) 
            {
                m_oMainStrategy.vSetLengthScaleFactor(dScaleFactor);
            }

            //////////////////////////////////////////////////////////////////////////
            // Take care of the current reflection list
            if( poCurrentPredictList )
                delete poCurrentPredictList;
            m_oMainStrategy.vSetExternalPredictReflnListPtr(NULL);
            //////////////////////////////////////////////////////////////////////////
            
            //////////////////////////////////////////////////////////////////////
            // Make m_oMainStrategy use external memory to store predicted reflections
            poCurrentPredictList = new Creflnlist();
            m_oMainStrategy.vSetExternalPredictReflnListPtr(poCurrentPredictList);
            //////////////////////////////////////////////////////////////////////
        
            //////////////////////////////////////////////////////////////////////
            // Make m_oMainStrategy use external memory to store scans
            // Note, we are not deleting an existing pCurrentScans, because it will be kept and
            // eventually deleted as part of the container object m_aoScanSolutionSet.
            pCurrentScans = new CstrategyScans();
            m_oMainStrategy.vSetExternalScanInfoPtr(pCurrentScans);

            if( 0 == nTrialsToAchieveCompleteness )
            {
                //// Set 2Theta to header
                m_nActive2ThetaPositionIndex = ii;
                vSetDetectorPositionToHeader(DTREK_MULTISTRATEGY_DET_2THETA);

                //// Get resolution range for the current detector 2Theta position
                if( !bGetResoRangeForCurrentDetectorPosition(dResoMin, dResoMax, sError) )
                {
                    printf("\nCould not set resolution for the current detector position.\n");
                    if( "" != sError )
                    {
                        DTREK_ERROR(1, sError);
                    }
                }
            }

            // Save resolution into the current scans object
            pCurrentScans->vSetResoRange(dResoMin, dResoMax);
            
            /// Generate reso bins and set current required resolution to strategy
            if( !bGenerateResoBins(dResoMin, dResoMax, c_nNumberOfResoBins) )
            {
                sError = "Crystal information in header is invalid.";
            
                DTREK_ERROR(3, sError); 
            }
            m_oMainStrategy.vSetExternalResoBinsPtr(m_poResoBins);
            ////////////////////////////////////////////////////////////////////////////////////////////////

            
            //////////////////////////////////////////////////////////////////////////////////////////
            // THIS IS THE CENTRAL CALL OF THIS FUNCTION. IT PASSES THE CONTROL TO DTSTRATEGY.
            nRet =  m_oMainStrategy.nExecute(m_nDTREKCommandArgumentCount, m_ppcDTREKCommandArguments);
            //////////////////////////////////////////////////////////////////////////////////////////
            
            if ( DTMAIN_ERROR_UNRECOGNIZED_COMMAND_OPTION == nRet )
            {
                sError = "Unrecognized command option";
                DTREK_ERROR(DTMAIN_ERROR_UNRECOGNIZED_COMMAND_OPTION, sError);
            }

            // If dtstrategy failed to achieve the desired "completeness minus tolerance",
            // we will try to adjust the tolerance.
            if( 0 == nTrialsToAchieveCompleteness++ &&
                0 != nRet &&  // 0 != nRet could be either due to a low completeness or to a low redundancy.
                              // We are only interested in the case when it is due to a low completeness.
                m_oMainStrategy.fGetBestCompleteness() > 0.0 &&  // If it's zero - don't bother with tolerance adjustments - something is totally wrong
                m_fUserTargetCompleteness - m_oMainStrategy.fGetCompletenessTolerance() > m_oMainStrategy.fGetBestCompleteness()
                )
            {
                // Automatic tolerance adjustment.
                // Determine a realistic tolerance.
                if( m_oMainStrategy.fGetBestCompleteness() - floor(m_oMainStrategy.fGetBestCompleteness()) < 0.05 ) // If the best result is too close to its floor,
                    fRealisticCompleteness = floor(m_oMainStrategy.fGetBestCompleteness()) - 0.5;                 // chances are that the floor isn't a realistic target either.
                else
                    fRealisticCompleteness = floor(m_oMainStrategy.fGetBestCompleteness());
                
                printf("\nResetting completeness tolerance to %.1f%%\n", m_fUserTargetCompleteness - fRealisticCompleteness);

                vUpdateTargetCompletenessTolerance(m_fUserTargetCompleteness,
                                                   m_fUserTargetCompleteness - fRealisticCompleteness);
                
                //////////////////////////////////////////////////////////////////////////////////////////////////////
                //Since we are going to call 'continue', pCurrentScans will not be included into m_aoScanSolutionSet
                //(and deleted eventually). So we need to delete pCurrentScans right now.
                if( pCurrentScans )
                {
                    delete pCurrentScans;
                    pCurrentScans = NULL;
                }
                ///////////////////////////////////////////////////////////////////////////////////////////////////////

                continue; // re-run strategy with the larger tolerance
            }
            else
                bCompletenessAchieved = true;
        }

        // Save resultant refln list.
        if( poCurrentPredictList->nGetNumReflns() > 0 )
            m_poPredictReflnList->nInsertListFrom(*poCurrentPredictList);

        // Save resultant scan solution
        
        m_aoScanSolutionSet.push_back(pCurrentScans);
    
        // Get the automatically-generated scale factor, if necessary.
        if( !m_bScaleSetFromCommandLine && ii == 0 && ii != m_ad2ThetaPositions.size() - 1 )
        {
            dScaleFactor = m_oMainStrategy.fGetLengthScaleFactor();
            printf("\nSetting scale factor to %.2f\n", dScaleFactor);
        }
    } // end of main loop
    /////////////////////////////////////////////////////////////////////////////////////////////

    if( poCurrentPredictList )
    {
        delete poCurrentPredictList;
        poCurrentPredictList = NULL;
    }
    
    ///////////////////////////////
    // Load "previous" reflections
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // For the final statistics we will need to add "previous reflections", a.k.a "fixed scan",
    // (if they were specified in the input, of course) to our total ref list. 
    // NOTE: every time we call m_oMainStrategy.nExecute(), it takes the "previous reflections"
    // into consideration. But when m_oMainStrategy.nExecute() at the end prepares a reflection list 
    // for dtmultistrategy, it does not leave those "previous reflections" there. This makes sense, 
    // because if it did, we would end up here with as many "previous refection" lists as the number 
    // of m_oMainStrategy.nExecute() calls we made.

    // NOTE: It is important that the "previous reflections" are loaded here, but not actually added
    // to the main reflection list before the latter is (optionally) output.

    if( bLoadPreviousReflections(sError) && sError != "" )
    {
        DTREK_ERROR(1, sError);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Write out the total refln list if necessary
    if( m_sReflnListOut != "" )
    {
        printf("\nWriting %s\n", m_sReflnListOut.string());
        m_poPredictReflnList->nWrite(m_sReflnListOut);
    }
    
    ////////////////////////////////////////////////////////////////
    // Add the "previous list" to the main list.
    if( m_poPreviousReflnList )
        m_poPredictReflnList->nInsertListFrom(*m_poPreviousReflnList);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Now do statistics on the reflection list combined from ALL calls to m_oMainStrategy.nExecute() plus "previous" reflections.
    bDoExternalResoBinsStatistics(c_nNumberOfResoBins, m_poPredictReflnList, sError);
 
    printf("\n\ndtmultistrategy SCAN SOLUTION:\n");
    m_poResoBins->vListStats();
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ATTENTION: This call must be made prior to updating the header, because
    // while vPrintSummaryScanTable() prepares the summary output, it also updates m_dAllDetectorPositionsAllScansTotalRotWidth.
    vPrintSummaryScanTable();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Write out the header
    vUpdateHeader();
    
    m_poHeader->nWrite(m_sHeaderOut);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // WARNING: bCheckReflectionOverlap() call can only be made after the header has been updated (call vUpdateHeader() above)
    // WARNING: bCheckReflectionOverlap() uses m_poPredictReflnList, so we had to write that refln list out first
    // (see call m_poPredictReflnList->nWrite() above).
    // The reason we have to re-generate the list is because the list may have been generated with a scalefactor less than 1.0, 
    // so it cannot be used for overlap testing.
    if( !bCheckReflectionOverlap(sError) )
    {
        if( "" != sError )
        {
            DTREK_ERROR(1, sError);
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("\ndtmultistrategy: Done.\n");

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiStrategy::vError(const int nErrorNum, const Cstring& sMessage)
{
    if( 0 != nErrorNum ) 
    {
        Cstring     sErrorMsg = "\nERROR: ";
        sErrorMsg += sMessage;
        sErrorMsg += "\n";

        printf("%s\n",sErrorMsg.string());
    }

    printf("\ndtmultistrategy generates multiple scan information for\n"
         "a diffraction experiment that gives desired completeness and\n"
         "redundancy. It writes *dtmultistrategy.head.\n");

    printf("\ndtmultistrategy - Usage:\n"
         "dtmultistrategy header_file [options...]\n\n"
         "Command line options: Description:\n\n");

    printf(" header_file        Name of a file with a d*TREK header used to\n"
           "                    initialize the strategy procedure.\n\n"
           " -2thetamax max     Maximum allowed 2Theta value (deg) for automatic 2Theta\n"
           "                    calculation. Note: this value will be ignored, if it is\n"
           "                    not consistent with the hardware limits specified in the \n"
           "                    header\n\n"
           " -anom              Assume I+ != I-\n"
           "                    Default: Automatically determined.\n\n"
           " -auto2theta [ovlp] Generate detector 2Theta position(s) automatically,\n"
           "                    ovlp is the minimum fractional overlap between successive\n"
           "                    detector positions. Default ovlp=0.25. Note: this option\n"
           "                    is incompatible with -set2theta and requires -reso option\n" 
           "                    to be specified.\n\n"
           " -autodist [nSep]   Generate detector distance automatically. nSep is\n"
           "                    separation (pix) between borders of neighboring\n"
           "                    reflections. Default nSep=2. Note: this option is\n"
           "                    incompatible with -distance.\n\n"
           " -autoexposure [[fSigma [nShells]]]\n"
           "                    Generate image exposure time automatically. fSigma - desired\n"
           "                    I/Sigma ratio in the highest resolution shell after scaling\n"
           "                    and averaging data; nShells - number of resolution shells.\n"
           "                    Default: fSigma = 3.0, nShells = 10.\n"
           "                    Note: this option will only work if the image header already\n"
           "                    has exposure time evaluation information from a dtranker run.\n\n"
           " -cmin fCompletenessMin [fCompletenessTolerance]\n"
           "                    Minimum percentage completeness for rotation range\n"
           "                    solution(s). If fCompletenessTolerance is specified,\n"
           "                    then for each detector position dtstrategy will explore\n"
           "                    solutions within a range of completeness: fCompletenessMin\n"
           "                    +/- fCompletenessTolerance.\n"
           "                    Default: fCompletenessMin = 99; fCompletenessTolerance = 0.\n\n"
           //" -cmin fCompletenessMin [fMinIncCompPerDegree]\n"
           //"                    Minimum percentage completeness for rotation range\n"
           //"                    solution(s). If fMinIncCompPerDegree is specified,\n"
           //"                    then dtstrategy will not increase the total rotation width\n"
           //"                    beyond the point where each added degree adds less than\n"
           //"                    fMinIncCompPerDegree to the total completeness.\n"
           //"                    Default: fCompletenessMin = 100 fMinIncCompPerDegree = 0\n\n"
           " -crysgonio axis v1, v2, ... vN\n"
           "                    Setting crystal goniometer angles for scans to be checked.\n"
           "                    \"axis\" is a valid crystal goniometer axis name,\n"
           "                    v1, v2, ... vN - are axis values.\n"
           "                    NOTE: \"axis\" must match a crystal goniometer axis name\n"
           "                    in the input header file.\n\n"
           " -distance D        Set detector distance to D (mm). Note: this option is\n"
           "                    incompatible with -autodist.\n\n"
           " -distlimits d1 d2  Limits (mm) for automatic distance calculation. Note:\n"
           "                    this value will be ignored, if it is not consistent with\n"
           "                    the hardware limits specified in the header\n\n"
           " -firstauto2theta value\n"
           "                    Defines the first detector position (deg) to be generated\n"
           "                    by -auto2theta. Note: this option will be ignored, if the\n"
           "                    value is not consistent with the high 2Theta limit and/or\n"
           "                    required resolution limits.\n\n"
           " -fixrotstart       Fix rotation angle start at the input value.\n"
           "                    Default: not fixed.\n\n"
           " -help              Print this help text.\n\n"
		   " -mosaicity fMosaicity\n"
           "                    Crystal mosaicity (deg). If this option is not supplied,\n"
           "                    the mosaicity value from the header will be used.\n\n"
           " -out sHeaderOut    Name of output file with resultant scan rotation info.\n"
           "                    Default: *dtmultistrategy.head.\n\n"
           //" -checkoverlap nIntegBoxWidth fOsc | [fOsc1 fOscStep fOsc2]\n"
           //"                    For each scan solution calculate the percentage of\n"
           //"                    reflection overlap. Two reflections are considered \n"
           //"                    overlapped if their 2D-integration boxes with a width of\n"
           //"                    nIntegBoxWidth (pix) overlap or if their centers are less\n"
           //"                    than fOsc (degrees) of rotation apart.\n"
           //"                    To test several rotation widths enter fOsc1 - minimum width\n"
           //"                    fOscStep - width step and fOsc2 - maximum width.\n"
           //"                    The defaults are: nIntegBoxWidth = 20 pix, fOsc = 0.5 deg.\n\n"
           " -overlap [nDetCoordTest [fRotWidthTest [fRotWidthStepTest fRotWidthMaxTest]]]\n"
           "                    For each scan solution calculate the total percentage of\n"
           "                    overlapped HKL's. Two HKL's are considered overlapped\n"
           "                    if (1) they are predicted on the same image and (2) on that\n"
           "                    image the absolute difference between corresponding\n"
           "                    detector coordinates is less than nDetCoordTest (pix). If\n"
           "                    nDetCoordTest is not supplied or supplied as -1, the test\n"
           "                    value will be set to \"1.2 * average_spot_size\" where the\n"
           "                    average spot size comes from the input header. If the\n"
           "                    header does not have that information, the test value will\n"
           "                    be 20 pix. fRotWidthTest is the image rotation width (deg)\n"
           "                    to be tested.  To test several rotation widths supply also\n"
           "                    fRotWidthStepTest and fRotWidthMaxTest - step and maximum\n"
           "                    test rotation width. The default is to test only the image\n"
           "                    rotation width set by option -rotimagewidth or, if that\n"
           "                    is not present, the image rotation width in the input\n"
           "                    header. If the image rotation width is not available from\n"
           "                    the command line or header, the test value will be 0.5 deg.\n\n"
           " -overlaplist       Write reflection list for each overlap check.\n\n"
           " -pad fPad          Specifies a padding term.  The padding determines a small\n"
           "                    no-mans-land of fPad degrees at the start and end of the\n"
           "                    rotation range. Reflections in the no-mans-land are not\n"
           "                    counted when computing the completeness or redundancy.\n\n"
           " -predicted sPredictReflnlistFilename\n"
           "                    Name of a reflnlist file with previously\n"
           "                    predicted reflections. If this is specified no\n"
           "                    prediction is done.  Default: predict reflections.\n\n"
           " -prev sPrevReflnlistFilename\n"
           "                    Name of a reflnlist file with previously available\n"
           "                    reflections. Default: no previous list.\n\n"
           " -rangemax fRangeMax [nScan]\n"
           "                    Maximum rotation range (deg) allowed per scan. This will\n"
           "                    apply to all scans found, unless nScan is set.\n\n"
           " -ref sReflnlist    Write the selected reflections to the file sReflnlist.\n\n"
           " -reso fResoMin fResoMax\n"
           "                    Specifies the resolution range (angstrom) for completeness\n"
           "                    calculations. Note: this option must be present if\n"
           "                    option -auto2theta is present.\n\n"
           " -rmin fRedundancyMin\n"
           "                    Minimum allowed redundancy for rotation range\n"
           "                    solution(s). If fRedundancyMin==0.0, then we\n"
           "                    assume the user wants to collect an entire\n"
           "                    hemisphere.  Default: Not used.\n\n"
           " -rot  fStart fEnd  [nScan]\n"
           "                    Specifies the min and max rotation angle start and\n"
           "                    end (deg) to test for completeness. This will apply to\n"
           "                    all scans found, unless nScan is set.\n"
           "                    NOTE: option -rot can only set the rotation angle start and\n"
           "                    end limits within the rotation hardware limits specified\n"
           "                    in the input header file.\n"
           "                    \n\n"
           " -rotimagewidth  fWidth\n"
           "                    Image rotation width (deg). This value will be set\n"
           "                    to all strategy-generated scans.\n"
           "                    Note: this value overrides the image rotation width in the\n"
           "                    input header.\n\n"
           " -scale fScale      Factor to multiple cell lengths by, OR if >100, maximum\n"
           "                    number of reflections to use per scan.\n\n"
           " -set2theta T1 T2 ... TN\n"
           "                    Use these detector 2Theta positions (deg). Note: this\n"
           "                    option is incompatible with -auto2theta and\n"
           "                    -firstauto2theta.\n\n"
           " -totalrotmax fMaxRot\n"
           "                    Total maximum rotation width in a multiple scan solution\n"
           "                    per one detector position. Default: allow any width.\n\n"
           " -totalscansmax nMaxScans\n"
           "                    Total maximum number of scans in a multiple scan solution\n"
           "                    per one detector position. This value should be at least 2.\n"
           "                    Default: Allow any number of scans in a multiple scan\n"
           "                    solution.\n\n");

    printf("Examples:\n"
           " dtmultistrategy dtrefine.head -auto2theta -2thetamax -90\n"
           " -firstauto2theta 0 -crysgonio chi 0, 50 -crysgonio phi 0, 90 -reso 100 0.9\n"
           " -cmin 100 2 -overlap -1 0.5 0.5 5\n\n"
           " dtmultistrategy dtrefine.head -set2theta 46 80 -autodist 2 -distlimits 30 70\n"
           " -cmin 95 rmin 2.0 -scale 0.3\n\n");

    fflush(stdout);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bGenerate2ThetaPositions(Cstring& sError)
{
    sError = "";
    m_ad2ThetaPositions.clear(); // initialize
    
    if( !m_poHeader )
        return false;

    if( !bReadSetHardwareLimits(sError, DTREK_MULTISTRATEGY_DET_2THETA) )  // This call should set 2Theta max limit. If it fails, that's it.
    {
        if( "" == sError )
            sError = "Cannot generate automatic 2Theta. Limits are not specified.";
        return false;
    }

    Cstring              strDetectorPrefix(""); 
    if ( 0 != m_poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    // GET DETECTOR WIDTH.  
    float       afDetSizes[2] = {0.0, 0.0}; 
    if( 0 != m_poHeader->nGetValue(strDetectorPrefix+D_K_DetectorSize, 2, afDetSizes) )
        return false;
    
    // FIGURE OUT WHICH SIZE IS IN THE DIRECTION PARALLEL TO Y, i.e. perpendicular
    // to 2THETA ROT AXIS
    float       afDetVectors[6]={0.0f};
    if( 0 != m_poHeader->nGetValue(strDetectorPrefix+D_K_DetectorVectors, 6, afDetVectors) )
        return false;

    double      dYsize = afDetVectors[1] != 0.0 ? afDetSizes[0] : afDetSizes[1];
    //////////////////////////////////////////////////////////////////////////////////////////

    // CHECK DISTANCE. IT MUST BE SET IN ADVANCE
    if( 0.0 >= m_dDistance )
        return false;

    //GET 2-THETA INCREMENT BASED ON REQUIRED OVERLAP
    double      d2ThetaIncrement = atan2((double)dYsize*(0.5 - m_dFractionDetectorOverlap), m_dDistance) + 
                                   atan2((double)dYsize*0.5, m_dDistance);
    
    d2ThetaIncrement /= Gs_dRADIANS_PER_DEGREE;  // from radians to degrees
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // Figure out the maximum and minimum 2Theta based on the required resolution
    CExperimentalStrategy       oExpStrategy(*m_poHeader);
    double      dBestHighResoPossible = -1.0;
    double      dBestLowResoPossible  = -1.0;

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Starting with the SMALLEST allowed 2Theta and going UP, what is the 2theta needed to barely achieve the MAXIMUM resolution?
	oExpStrategy.vSetResoMethod(oExpStrategy.eResoOuter, 90.0); //mrp
    double      dMax2Theta = oExpStrategy.fCalcBest2ThetaSwingForTargetHighReso(m_segUser2ThetaLimits.dGetBeg(), 
                                                                                m_segUser2ThetaLimits.dGetEnd(), 
                                                                                m_dUserResoMax, 
                                                                                dBestHighResoPossible);
    
    if( cdStratDetUnrealisticValue == dMax2Theta )
    {
        printf("WARNING: Could not find a detector position to satisfy %f A resolution.\n", m_dUserResoMax);
        printf("         Will use the highest allowed 2Theta: %f deg\n",                     m_segUser2ThetaLimits.dGetEnd()*m_n2ThetaDirectionSign); 
        
        if( dBestHighResoPossible > 0.0 )  // Check if this has been calculated...
            printf("         This will give the highest resolution of %f A\n",              dBestHighResoPossible); 
        
        dMax2Theta = m_segUser2ThetaLimits.dGetEnd();

        // Reset the resolution
        m_dUserResoMax = dBestHighResoPossible;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Starting with the LARGEST allowed 2Theta and going DOWN, what is the 2theta needed to barely achieve the MINIMUM resolution?
    double      dMin2Theta = oExpStrategy.fCalcBest2ThetaSwingForTargetLowReso(m_segUser2ThetaLimits.dGetBeg(), 
                                                                               m_segUser2ThetaLimits.dGetEnd(), 
                                                                               m_dUserResoMin,
                                                                               dBestLowResoPossible);

    if( cdStratDetUnrealisticValue == dMin2Theta ) // if failed, use the limit
    {
        printf("WARNING: Could not find a detector position to satisfy %f A resolution.\n", m_dUserResoMin);
        printf("         Will use the lowest allowed 2Theta: %f deg\n",                     m_segUser2ThetaLimits.dGetBeg()*m_n2ThetaDirectionSign); 
        
        if( dBestLowResoPossible > 0.0 )  // Check if this has been calculated...
            printf("         This will give the lowest resolution of %f A\n",              dBestLowResoPossible); 

        dMin2Theta = m_segUser2ThetaLimits.dGetBeg();
        
        // Reset the resolution
        m_dUserResoMin = dBestLowResoPossible;
    }
    
    // If dMin2Theta > dMax2Theta, it means that the detector is so wide, that dMax2Theta position is enough for the whole required reso range!
    if( dMin2Theta > dMax2Theta )
       dMin2Theta = dMax2Theta;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // See if the user mandates a certain first 2Theta. If it is valid, make it the minimum 2Theta.
    if( m_segUser2ThetaLimits.bIsIn(fabs(m_dUser2ThetaMust)) ) // if m_dUser2ThetaMust is a valid value...
    {
        dMin2Theta = min(dMin2Theta, fabs(m_dUser2ThetaMust)); // if m_dUser2ThetaMust is less than the reso defined minimum, use it!  
    }

    // Take care of a special (weird) case when d2ThetaIncrement==0 
    if( d2ThetaIncrement == 0.0 )
    {
        m_ad2ThetaPositions.push_back(dMin2Theta);
        return true;
    }
    
    ///////////////////////////////////
    //Generate all necessary 2Thetas
    bool        bDone = false;
    double      dCur2Theta = dMin2Theta;
    
    while( true )
    {
        if( dCur2Theta >= dMax2Theta )
        {
            dCur2Theta = dMax2Theta;
            
            bDone = true;
        }
        
        m_ad2ThetaPositions.push_back(dCur2Theta);
        
        if( bDone )
            break;
        else
            dCur2Theta += d2ThetaIncrement;
    }

    vRedistributeDetectorOverlaps();

    // Now apply the direction sign
    int     ii=0;
    for(ii=0; ii < (int)m_ad2ThetaPositions.size(); ii++)
    {
        m_ad2ThetaPositions[ii] *= m_n2ThetaDirectionSign;
    }
    
    ////////////////////////////////////////////////
    // Round them off, keeping only the integer part
    for(ii=0; ii < (int)m_ad2ThetaPositions.size(); ii++ )
    {
        m_ad2ThetaPositions[ii] = dRoundOff(m_ad2ThetaPositions[ii], 0);
    }
    
    
    //////////////////////////////////////////////////////////////////////////
    // Now report the generated 2Theta values
    Cstring     strMsg("\nAutomatically generated 2Theta positions (deg): ");
    for(ii=0; ii < (int)m_ad2ThetaPositions.size(); ii++)
    {
        strMsg += Cstring(m_ad2ThetaPositions[ii]);
        strMsg += ' ';
    }
    strMsg += '\n';

    printf("%s", strMsg.string());
    //////////////////////////////////////////////////////////////////////////
    
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bGetResoRangeForCurrentDetectorPosition(double& dResoMin, double& dResoMax, Cstring& sError)
{
    // By the time this function is called the overall resolution (m_dUserResoMin,m_dUserResoMax) must already be determined
    if( !m_bResolutionSet )
        return false;

    // Find resolution for the current detector 2Theta position
    Cstring       sEdge = "edge";
    if( !bParseCalculateResolution(sEdge,    // base class function
                                   dResoMin, 
                                   dResoMax,
                                   sError) )    
        return false;
    
    // Make sure dResoMin is indeed the minimum resolution
    if( dResoMin < dResoMax )
        std::swap(dResoMin, dResoMax);

    // If the current 2Theta position is neither first nor last, just pass the found resolution through
    if( m_nActive2ThetaPositionIndex > 0 && m_nActive2ThetaPositionIndex < (int)m_ad2ThetaPositions.size()-1 )
        return true;
    
    if( 0 == m_nActive2ThetaPositionIndex )
    {
        if( dResoMax >= m_dUserResoMin ) // just a safety check
        {
            sError = "A problem is encountered in the algorithm, while trying to determine the current resolution.";
            return false;
        }
        
        dResoMin = m_dUserResoMin;
    }
    
    if( m_ad2ThetaPositions.size()-1 == m_nActive2ThetaPositionIndex )
    {
        if( dResoMin <= m_dUserResoMax ) // just a safety check
        {
            sError = "A problem is encountered in the algorithm, while trying to determine the current resolution.";
            return false;
        }
        
        dResoMax = m_dUserResoMax;
    }

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiStrategy::vExpandReflnList(Creflnlist* poReflnList)
{
  if( !poReflnList )
    return;
  
  poReflnList->nExpandGetField(Creflnlist::ms_snPackedHKL);

  poReflnList->nExpandGetField(Creflnlist::ms_snNonunfFlag);

  poReflnList->nExpandGetField(Cstrategy::ms_sf2STLcu);

  poReflnList->nExpandGetField(Creflnlist::ms_sfCalcRotStart);
  poReflnList->nExpandGetField(Creflnlist::ms_sfCalcRotEnd);

  poReflnList->nExpandGetField(Creflnlist::ms_sfCalcPx0);
  poReflnList->nExpandGetField(Creflnlist::ms_sfCalcPx1);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bGenerateResoBins(double dResoMin, double dResoMax, int nNumberOfBins)
{
    if( !m_poHeader )
        return false; // safety 

    printf("\nGenerating resolution bins. Resolution range: %f %f A.\n", dResoMin, dResoMax);

    if( m_poResoBins )
        delete m_poResoBins;

    m_poResoBins = new CResoStrategyBins(dResoMin, dResoMax, nNumberOfBins);  
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function detects the auto2Theta option on the command line, and generates automatic two-theta positions.
// NOTE: if the function returns false, the only case when sError can be empty is if auto2Theta option is not supplied at all.
bool CDTMainMultiStrategy::bSetAuto2Theta(Cstring& sError)
{
    char    acTemp[DTREK_MAX_LINE_LENGTH];
    sError = "";
    m_ad2ThetaPositions.clear(); // initialize
    
    Cstring                 sArgs("");     

    if( !bGetCommandOption(c_pcAuto2Theta, sArgs) )
    {
        return false; // auto2theta option is not supplied    
    }

    printf("\nSetting automatic 2Theta positions...\n");

    /////////////////////////////////////////////////////////////////////////
    // Check if -set2theta is supplied.
    Cstring     sArgs1("");
    if( bGetCommandOption(c_pcSet2Theta, sArgs1, false) )
    {
        sprintf(acTemp, "Options %s and %s are incompatible.", c_pcAuto2Theta, c_pcSet2Theta);
        sError = acTemp;

        return false;
    }
    /////////////////////////////////////////////////////////////////////////

    // Figure out parameters passed with the auto two-theta option
    if( "" != sArgs )
    {
        if( !sArgs.bIsNumeric() )
        {
            sprintf(acTemp, "Option %s does not have a numeric argument.", c_pcAuto2Theta);
            sError = acTemp;

            return false;
        }
        else
        {
            m_dFractionDetectorOverlap = atof(sArgs.string());
            if( m_dFractionDetectorOverlap < 0.0 )  // sanity check
                m_dFractionDetectorOverlap *= -1; // could return an error condition, of course, but let's help the user...
            
            if( m_dFractionDetectorOverlap < 0.0 || m_dFractionDetectorOverlap >= 1.0 )
            {
                sprintf(acTemp, "Option %s has invalid detector overlap %f. The value must be between 0 and 1.", 
                                c_pcAuto2Theta, m_dFractionDetectorOverlap);
                sError = acTemp;

                return false;
            }
        }    
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Now check if the user supplied resolution explicitly.
    Cstring     strTemp("");
    if( !bSetResolutionLimitsFromCommandLine(strTemp, sError) )
    {
        if( sError == "" )
            sError = "Cannot generate automatic 2theta positions, because explicit resolution limits are unknown. Check -reso option.";
        
        m_bResolutionSet = false;
        
        return false;
    }
    else
        m_bResolutionSet = true;


    // Check if we have  -firstauto2theta option
    if( bGetCommandOption(c_pcFirstAuto2Theta, sArgs) )
    {
        if( sArgs.bIsNumeric() )
            m_dUser2ThetaMust = atof(sArgs.string());
        //else if( sArgs == "center" || sArgs == "edge")
        //{
        //
        //}
        else
        {
            sprintf(acTemp, "Option %s has invalid parameter.", c_pcFirstAuto2Theta);
            sError = acTemp;
            
            return false;
        }
    }


    // Now generate the 2Theta positions    
    if( !bGenerate2ThetaPositions(sError) )
    {
        if( sError == "" )
            sError = "Cannot generate automatic 2theta positions.";
        return false;
    }
    
    return true;
}    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function detects the set2theta option on the command line, and sets up two-theta values
bool CDTMainMultiStrategy::bSetUserInput2Theta(Cstring& sError)
{
    char    acTemp[DTREK_MAX_LINE_LENGTH];

    sError = "";

    if( !bGetCommandOption(c_pcSet2Theta, m_ad2ThetaPositions, " ,") )
    {
        return false; // set2theta option is not supplied    
    }

    if( 0 == m_ad2ThetaPositions.size() )
    {
        printf("\nWARNING: option %s does not have arguments and will be ignored\n", c_pcSet2Theta);
        return false;
    }

    /////////////////////////////////////////////////////////////////////////
    // Just a consistency check....
    Cstring     sArgs("");
    if( bGetCommandOption(c_pcAuto2Theta, sArgs, false) )
    {
        sprintf(acTemp, "Options %s and %s are incompatible.", c_pcAuto2Theta, c_pcSet2Theta);
        sError = acTemp;

        return false;
    }
    /////////////////////////////////////////////////////////////////////////

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function detects -reso option on the command line and parses it.  If -reso option is not supplied, the
// function returns false with sError empty. The function expects at least high
// resolution to be specified on the command line, otherwise it will return false with an error.
// If the resolution option is one non-numeric token, the function will return false, but will set
// sResoOption equal to that token for possible further parsing. 
bool CDTMainMultiStrategy::bSetResolutionLimitsFromCommandLine(Cstring& sResoOption, Cstring& sError)
{
    const double    cdMaxLowReso = 999.99;
    
    sError = "";

    m_dUserResoMin = -1.0;
    m_dUserResoMax = -1.0;  // initialize


    Cstring                 sArgs("");     
    std::vector<Cstring>    asParams;
    int                     nParamsCount=0;

    // See if the user specified resolution
    if( !bGetCommandOption(c_pcReso, sArgs) )
        return false;

    sArgs.nListToVector(asParams, " ");
    nParamsCount = asParams.size();

    if( 0 == nParamsCount )
    {
        sError = "-reso option has no arguments";
        return false;
    }

    //Now we know for sure that asParams has at least one element...
    if( !asParams[0].bIsNumeric() )
    {
        sError = "-reso option does not specify resolution explicitly.";
        
        sResoOption = asParams[0];

        return false;
    }

    // Now let's try and get the resolution from the argument string
    double      dReso1 = -1.0;
    double      dReso2 = -1.0;

    if( asParams[0].bIsNumeric() )
        dReso1 = atof(asParams[0].string());

    if( nParamsCount > 1 )
        dReso2 = atof(asParams[1].string());

    // If a value is not supplied, treat it as if the low resolution is not supplied.
     if( dReso1 <= 0.0 )
        dReso1 = cdMaxLowReso; // a low resolution default

    if( dReso2 <= 0.0 )
        dReso2 = cdMaxLowReso; // a low resolution default

    // Just a safety check for things like "-reso -1"
    if( cdMaxLowReso == dReso1 && cdMaxLowReso == dReso2 )
    {
        sError = "-reso option does not specify resolution explicitly.";
        return false;
    }

    //OK, now where is the low and where is the high? Let's finally set the class members.
    m_dUserResoMin = max(dReso1, dReso2);
    m_dUserResoMax = min(dReso1, dReso2);

    // Cap the low resolution
    m_dUserResoMin = min(cdMaxLowReso, m_dUserResoMin);

    printf("\nResolution limits from command line (A): %f %f\n", m_dUserResoMin, m_dUserResoMax);

    return true;
}    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The purpose of this function is as follows. Suppose, the user specified 2Theta positions, but didn't provide resolution.
// Then we should try and figure out the required resolution, using the detector information and the user-specified
// 2Theta values. Default: edge resolution.
bool CDTMainMultiStrategy::bSetResolutionLimitsFrom2ThetaPositions(Cstring& sResoOption, Cstring& sError)
{
    if( "" == sResoOption )
        sResoOption = "edge";   // default
    
    m_dUserResoMin = -1.0;
    m_dUserResoMax = -1.0;  // initialize
    
    if( 0 == m_ad2ThetaPositions.size() )
        return false; // 2theta positions must be known!
    
    if( !m_poHeader )
        return false;  // safety check

    Cstring              strDetectorPrefix(""); 
    if ( 0 != m_poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    Cgoniometer         oDetGonio(*m_poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
        return false;
    
    float      f2Theta_sav = oDetGonio.fGetSwing(); // the value that is currently in the header.

    double      dReso1 = -1.0;
    double      dReso2 = -1.0;

    for( int ii=0; ii < (int)m_ad2ThetaPositions.size(); ii++ )
    {
        oDetGonio.nSetSwing((float)m_ad2ThetaPositions[ii]);
        
        oDetGonio.nUpdateHeader(m_poHeader, strDetectorPrefix);
        
        
        if( !bParseCalculateResolution(sResoOption,    // base class function
                                       dReso1, 
                                       dReso2,
                                       sError) )
            return false;
        
        // Order variables dReso1 and dReso2, so that dReso1 is the maximum resolution
        if( dReso1 > dReso2 )
            std::swap(dReso1, dReso2);

        // Update our max and min resolutions
        if( ii == 0 || m_dUserResoMax > dReso1 )
            m_dUserResoMax = dReso1;
        
        if( ii == 0 || m_dUserResoMin < dReso2 )
            m_dUserResoMin = dReso2;
    }
    
    // Restore the original 2Theta value in the header
    oDetGonio.nSetSwing(f2Theta_sav);
    oDetGonio.nUpdateHeader(m_poHeader, strDetectorPrefix);
    
    printf("\nCalculated total resolution range for detector position(s): %f %f A\n", m_dUserResoMin, m_dUserResoMax);
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function takes all strategy solution sets generated by repeated calls to DTMainStrategy::Execute() and
// writes them out into the header object.
void CDTMainMultiStrategy::vUpdateHeader()
{   
    if( !m_poHeader )
        return; // safety

    printf("\nUpdating multistrategy header\n");

    vRemoveMultiScanInfoFromHeader(); // clean-up
    
    // Since the original input scans in the header may have changed due to 
    // some activity in DTMainStrategy and Cstrategy - restore them.
    // but do it only if the user put the original test scans into the header
    ////////////////////////////////////////////////////////////////////////////////
    if( !m_bUseOldHeaderScanFormat )
    {
        
        if( !m_bTestScansGeneratedFromCommandLine &&
            !m_bTestScansGeneratedFromOldHeaderFormat &&
            m_paoSaveTestScans )
            
            m_paoSaveTestScans->nUpdateHeader(m_poHeader, 0, D_K_StrategyTestPrefix, false);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    int                 nCurScanOffset = 0;
    CstrategyScans*     pScans = NULL;
    int                 nScansUsed = 0;
    for(int nSet = 0; nSet < (int)m_aoScanSolutionSet.size(); nSet++)
    {
        pScans = m_aoScanSolutionSet[nSet];
        
        // It is important that we are setting the new exposure time AFTER the saved original test scans 
        // have been written into the header, because we want to preserve the original input exposure time.
        if( m_fImageAutoExposureTime > 0.0 )
            pScans->vSetImageExposureTime(m_fImageAutoExposureTime);

        if( !m_bUseOldHeaderScanFormat )
            nScansUsed = pScans->nUpdateHeader(m_poHeader, nCurScanOffset, D_K_StrategyResultPrefix);
        else
            nScansUsed = pScans->nUpdateHeader(m_poHeader, nCurScanOffset, "");
        
        nCurScanOffset += nScansUsed;
    }

    if( 0 == nCurScanOffset )
        return;    // this means no scans have been outputted
    
    ////////////////////////////////////////////////////////////////////////////////
    // The list of used scans is trivial: it's simply ALL scans we have outputted.
    int*    pnScansUsed = new int[nCurScanOffset];
    for(int ii=0; ii < nCurScanOffset; ii++)
        pnScansUsed[ii]=ii;   

    m_poHeader->nReplaceValue(D_K_ScanSeqSelected, nCurScanOffset, pnScansUsed);
    /////////////////////////////////////////////////////////////////////////////////

    m_poHeader->nReplaceValue(D_K_DtstrategyPercentComplete, m_poResoBins->dGetCompleteness());
    m_poHeader->nReplaceValue(D_K_DtstrategyRedundancy,      m_poResoBins->dGetMultiplicity());
    m_poHeader->nReplaceValue(D_K_DtstrategyRedundancyDev,   m_poResoBins->dGetMultiplicityDev());

    delete [] pnScansUsed;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Write ranking information
    m_poHeader->nReplaceValueDictionary("DTMULTISTRATEGY_RANK", "TOTAL_ROTATION_WIDTH", m_dAllDetectorPositionsAllScansTotalRotWidth);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bParseMosaicityInput(Cstring& sError)
{
    sError = ""; // initialize
    
    double                 fValue = -1.0;
    if( !bGetCommandOption(c_pcMosaicity, fValue) )
    {
        return false; // mosaicity not supplied    
    }
    
    if( 0.0 > fValue )
    {
        sError = "Input mosaicity cannot be negative";
        return false;
    }

    m_fUserInputMosaicity = fValue;

    return true;
}
/////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bParseReflnOverlapCheckParams(Cstring& sError)
{
    sError = ""; // initialize
    
    Cstring                 sArgs("");
    if( !bGetCommandOption(c_pcCheckOverlap, sArgs) )
        return false; // option -overlap not given    

    m_pstOverlapCheckParams = new STRATEGY_REFLN_OVERLAP_CHECK();
    
    if( !m_pstOverlapCheckParams->bParseInput(sArgs) )
    {
        sError = "Option ";
        sError += Cstring(c_pcCheckOverlap);
        sError += " has wrong number of parameters.";
        delete m_pstOverlapCheckParams;
        m_pstOverlapCheckParams = NULL;
        return false;
    }
    
    // For now we cannot handle "zero" case
    if( 0.0 == m_pstOverlapCheckParams->m_dDistTest )
    {
        sError = "Option ";
        sError += Cstring(c_pcCheckOverlap);
        sError += " has zero value for nDetCoordTest";
        delete m_pstOverlapCheckParams;
        m_pstOverlapCheckParams = NULL;
        return false;
    }


    m_bWriteOverlapReflnList = bGetCommandOption(c_pcOverlapList);

	return true;
}
////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bCheckReflectionOverlap(Cstring& sError)
{
    if( NULL == m_pstOverlapCheckParams )
        return false; // check not required
    
    if( 0 == m_aoScanSolutionSet.size() )
        return false; // we have no solutions
    
    // See if we need to set rotation width
    if( m_pstOverlapCheckParams->m_dOscTest_1 == 0.5 &&  
        m_pstOverlapCheckParams->m_dOscStep   == 0.0 &&  
        m_pstOverlapCheckParams->m_dOscTest_2 == m_pstOverlapCheckParams->m_dOscTest_1 ) // default setting: test just 0.5 deg
    {
        if( m_fUserImageRotationWidth >= 0.0 )
        {
            // Change the default setting: check just the rot width from the header 
            m_pstOverlapCheckParams->m_dOscTest_1 = m_fUserImageRotationWidth;  
            m_pstOverlapCheckParams->m_dOscStep   = 0.0;  
            m_pstOverlapCheckParams->m_dOscTest_2 = m_pstOverlapCheckParams->m_dOscTest_1;
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    
    printf("\n\nCalculating reflection overlap...\n\n");
    char  acTemp[DTREK_MAX_LINE_LENGTH];

    int   nNumberOfPredictedReflns = -1;
    int   nNumberOfOverlapReflns   = -1;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Figure out the spot size
    bool                bDefaultSpotSize = false;
    bool                bDistanceToTestFromCommandLine = true;
    const double        cfSpotSizeCoefficient = 1.2;  // 1.2 is just a coefficient to provide some 
                                                      // buffer between spots, so that spots don't 
                                                      // "touch" each other.
    
    if( m_pstOverlapCheckParams->m_dDistTest < 0.0 )  // command line says: try and get it from the header
    {
        bDistanceToTestFromCommandLine = false;
        
        CExperimentalStrategy       oExpStrategy(*m_poHeader);

        if( !oExpStrategy.bGetSpotSizeFromHeader(m_stSpotSizeInfo) )
        {
            m_pstOverlapCheckParams->m_dDistTest = CExperimentalStrategy::ms_fDefaultSpotSize[1]; // default for average spot size
            bDefaultSpotSize = true;
        }
        else 
        {
            m_pstOverlapCheckParams->m_dDistTest = m_stSpotSizeInfo.fAverageSpotSize;
        }
        
        m_pstOverlapCheckParams->m_dDistTest *= cfSpotSizeCoefficient;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    m_pstOverlapCheckParams->vPrint();

    Cpredict*   poPredict = NULL;
    
    double      dTwoTheta = cdStratDetUnrealisticValue;
    double      dDistance = cdStratDetUnrealisticValue;
    
    std::vector<OVERLAP_RESULTS>     vecOverlapResults; // to store all overlap percentages for each scan
    
    double              dTestOscWidth = 0.0;

    // For now we do not support more than 10 possible rotation widths. It's just a log file print-out issue.
    const int   cnMaxNumberOfOscWidths = 10;

    Creflnlist*         poCurrentPredictList = new Creflnlist();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // WARNING: nScanReferenceNumber must have the same meaning as the one we use for each scan 
    // when we save it into the header in CstrategyScans::nUpdateHeader()
    // This is critical because below we are constructing CPredict based on the scan in the header,
    // whereas the resolution limits come from m_aoScanSolutionSet array, not from the header! 
    // to-do: set crystal goniometer in the predict object from m_aoScanSolutionSet, so we don't need to read the header.
    int     nScanReferenceNumber = -1;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Cstring     sOverlapReflnListBaseName = sDtrekGetPrefix() + "dtmultistrategy_overlap_";
    Cstring     sOverlapReflnListName("");
    Cstring     sTestWidth("");

    double      fMosaicity = -1.0;

    Cstring     sPrefix = m_bUseOldHeaderScanFormat ? "" : D_K_StrategyResultPrefix;

    for(int ii=0; ii < (int)m_aoScanSolutionSet.size(); ii++)
    {
        for(int jj=0; jj < m_aoScanSolutionSet[ii]->nGetNumScans(); jj++)
        {
            if( !m_aoScanSolutionSet[ii]->bIsScanUsed(jj) )
                continue;

            nScanReferenceNumber++;
            
            poPredict = new Cpredict(*m_poHeader, NULL, nScanReferenceNumber, sPrefix);
            
            if( !poPredict->bIsAvailable() )
            {
                sprintf(acTemp, "Could not create predict object for scan #%d.", nScanReferenceNumber);
                sError = acTemp;
        
                delete  poPredict;
                poPredict = NULL;

                delete poCurrentPredictList;
                poCurrentPredictList = NULL;
        
                return false;
            }
            
            if( 0 != poPredict->nSetupSourceCrystal() )
            {
                sprintf(acTemp, "Could not set-up the source and crystal info in the predict object for scan #%d.", nScanReferenceNumber);
                sError = acTemp;
        
                delete  poPredict;
                poPredict = NULL;
                
                delete poCurrentPredictList;
                poCurrentPredictList = NULL;
        
                return false;
            }

            /////////////////////////////////////////////////////////
            if( m_fUserInputMosaicity < 0.0 )  // i.e. it was not supplied on the command line
                fMosaicity = poPredict->m_ppoCrystal[0]->fGetMosaicity(); // for now assume just one crystal
            else
            {
                fMosaicity = m_fUserInputMosaicity;
                poPredict->vSetCrysMosaicity((float)m_fUserInputMosaicity);
            }
            /////////////////////////////////////////////////////////

            OVERLAP_RESULTS     stCurrentScanOverlapResults;
            stCurrentScanOverlapResults.nIndex = nScanReferenceNumber;

            // SET DETECTOR in PREDICT
            
            dTwoTheta = m_aoScanSolutionSet[ii]->dGetTwoTheta();
            dDistance = m_aoScanSolutionSet[ii]->dGetDistance();

            poPredict->m_ppoDetector[0]->m_poGoniometer->nSetSwing((float)dTwoTheta);
            poPredict->m_ppoDetector[0]->m_poGoniometer->nSetDistance((float)dDistance);

            //poPredict->nSetupDetector();

            // Restrict the Cpredict-set resolution by the user-specified resolution. 
            poPredict->vSetResolution(m_aoScanSolutionSet[ii]->dGetResoMin(),  
	                                  m_aoScanSolutionSet[ii]->dGetResoMax());

            stCurrentScanOverlapResults.dResoMin = m_aoScanSolutionSet[ii]->dGetResoMin();
            stCurrentScanOverlapResults.dResoMax = m_aoScanSolutionSet[ii]->dGetResoMax();


            if( 0 != poPredict->nSetupDetector() )
            {
                sprintf(acTemp, "Could not set-up the detector info in the predict object for scan #%d.", nScanReferenceNumber);
                sError = acTemp;
        
                delete  poPredict;
                poPredict = NULL;
        
                delete poCurrentPredictList;
                poCurrentPredictList = NULL;
                
                return false;
            }

            //Crotation&        oRotation = *(poPredict->m_poRotation);
            //oRotation.nSetRotRange(m_aoScanSolutionSet[ii]->dGetScanRotBegin(jj), m_aoScanSolutionSet[ii]->dGetScanRotEnd(jj));
            
            //poPredict->m_poCrysGonio->nList();
            //poPredict->m_poRotation->nList(11, poPredict->m_poCrysGonio);
            //poPredict->m_ppoDetector[0]->m_poGoniometer->nList();
            ///////////////////////////////////////////////////////////////////////////
        
            //if( 0 != poPredict->nPredict(NULL, m_poPredictReflnList) ) 
            poCurrentPredictList->vDeleteAll();
            if( 0 != poPredict->nPredict(NULL, poCurrentPredictList) ) 
            {
                sprintf(acTemp, "Failure during prediction for scan #%d.\n", nScanReferenceNumber);
                sError = acTemp;
            
                delete poPredict;
                poPredict = NULL;
            
                delete poCurrentPredictList;
                poCurrentPredictList = NULL;
            
                return false;
            }

            nNumberOfPredictedReflns = poCurrentPredictList->nGetNumReflns();

            if( 0 == nNumberOfPredictedReflns )
            {
                printf("\nWARNING: Zero reflections predicted for scan #%d.\n", nScanReferenceNumber);
            
                continue;
            }
            
            printf("\nNumber of reflections predicted for scan %d: %d\n", nScanReferenceNumber, nNumberOfPredictedReflns);
            printf("\n\nCalculating reflection overlap for scan %d...\n\n", nScanReferenceNumber);
            
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Calculate overlap
            dTestOscWidth = m_pstOverlapCheckParams->m_dOscTest_1;
            
            while( dTestOscWidth <= m_pstOverlapCheckParams->m_dOscTest_2 && 
                   (stCurrentScanOverlapResults.afOverlapsPercentage).size() < cnMaxNumberOfOscWidths )
            {
                nNumberOfOverlapReflns   = poCurrentPredictList->nOverlapCheck(m_pstOverlapCheckParams->m_dDistTest, 
                                                                               dTestOscWidth,
                                                                               poPredict->m_poRotation->fGetRotStart(),
                                                                               true); // true means: do display overlapped reflections
                if( nNumberOfOverlapReflns < 0 ) // error returned!
                {
                    nNumberOfOverlapReflns = 0;
                }
                
                if( m_bWriteOverlapReflnList )
                {
                    sTestWidth = Cstring(dTestOscWidth);
                    for(int iPos=0; iPos < (int)(sTestWidth.GetLength()); iPos++)
                    {
                        if('.' == sTestWidth.GetAt(iPos) )
                            sTestWidth.SetAt(iPos, 'p'); // replace dot with letter 'p'
                    }

                    sOverlapReflnListName = sOverlapReflnListBaseName + Cstring(nScanReferenceNumber);
                    sOverlapReflnListName += "_";
                    sOverlapReflnListName += sTestWidth;
                    sOverlapReflnListName += ".ref";
                    poCurrentPredictList->nWrite(sOverlapReflnListName);
                }
            
                stCurrentScanOverlapResults.afOverlapsPercentage.push_back((double)nNumberOfOverlapReflns / (double)nNumberOfPredictedReflns * 100.0);

                if( 0.0 == m_pstOverlapCheckParams->m_dOscStep )
                    break;
            
                dTestOscWidth += m_pstOverlapCheckParams->m_dOscStep;
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Save overlap test results
            vecOverlapResults.push_back(stCurrentScanOverlapResults);
            
            delete  poPredict;
            poPredict = NULL;
        }
    }
    
    delete poCurrentPredictList;
    poCurrentPredictList = NULL;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    printf("\nReflection overlap has been tested with the following assumptions:\n");
    
    if( !bDistanceToTestFromCommandLine )
    {
        printf("- The average spot size is: %d pix",(int)(m_pstOverlapCheckParams->m_dDistTest/cfSpotSizeCoefficient));   
        if( bDefaultSpotSize )
            printf(" (default spot size)");
        else
            printf(" (from input header)");
        printf("\n");
    }

    printf("- The minimum non-overlap coordinate difference is:\n");
    printf("   %d pix ", (int)m_pstOverlapCheckParams->m_dDistTest);
    if( !bDistanceToTestFromCommandLine )
        printf("(%.1f * \"average spot size\")\n",cfSpotSizeCoefficient);   
    else
        printf("(from command line)\n");
    

    printf("- The crystal mosaicity is: %.2f deg\n", fMosaicity);
    printf("\nReflection overlap (%%) predicted as a function of image rotation width (deg)\n");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    dTestOscWidth = m_pstOverlapCheckParams->m_dOscTest_1;

    Cstring     sOut1 = " #  Resolution			 Image rotation width (deg)";
    Cstring     sOut2 = "    range (A)     ";
    
    int     nOscWidths = -1;
    int     iTestNumberCorrespToOscWidthFromHeader = -1;

    while( dTestOscWidth <= m_pstOverlapCheckParams->m_dOscTest_2 && nOscWidths < cnMaxNumberOfOscWidths - 1 )
    {
        nOscWidths++;
        
        if( m_fUserImageRotationWidth == dTestOscWidth )
            iTestNumberCorrespToOscWidthFromHeader = nOscWidths;

        sprintf(acTemp, "%5.2f ", dTestOscWidth);
        sOut2 += acTemp;
        
        if( 0.0 == m_pstOverlapCheckParams->m_dOscStep )
            break;
        
        dTestOscWidth += m_pstOverlapCheckParams->m_dOscStep;
    }
    //////////////////////////////////////////////////////////////////////////////////////
    Cstring     sTableDelineator("-----------------------------------------------------------------------------\n");
    printf(sTableDelineator.string());
    printf("%s\n", sOut1.string());
    printf("%s\n", sOut2.string());
    printf(sTableDelineator.string());
    
    int     iScan = 0;
    int     iTest = 0;
    for(iScan = 0; iScan < (int)vecOverlapResults.size(); iScan++)
    {
        sprintf(acTemp, "%2d %7.2f - %5.2f", iScan,
                                             vecOverlapResults[iScan].dResoMin,
                                             vecOverlapResults[iScan].dResoMax);
        sOut1 = acTemp;
        
        for(iTest=0; iTest < (int)(vecOverlapResults[iScan].afOverlapsPercentage.size()); iTest++)
        {
            sprintf(acTemp, "%5.1f ", vecOverlapResults[iScan].afOverlapsPercentage[iTest]);
            sOut1 += acTemp;
        }
        
        printf("%s", sOut1.string());
        printf("\n");
    }
    printf(sTableDelineator.string());
    
    // Do warnings
    const   int c_nBadOverlapPercentage = 5; 
    if( iTestNumberCorrespToOscWidthFromHeader > -1 )
    {
        for(iScan = 0; iScan < (int)vecOverlapResults.size(); iScan++)
        {
            for(iTest=0; iTest < (int)(vecOverlapResults[iScan].afOverlapsPercentage.size()); iTest++)
            {
                if( iTest == iTestNumberCorrespToOscWidthFromHeader ) 
                {
                    if( vecOverlapResults[iScan].afOverlapsPercentage[iTest] > c_nBadOverlapPercentage )
                        printf("\nWARNING: for scan #%d the image rotation width in the header (%.2f deg)\n"
                               "         shows more than %d%% overlaps\n",
                                 iScan, m_fUserImageRotationWidth, c_nBadOverlapPercentage);
                    else if( 0.0 == vecOverlapResults[iScan].afOverlapsPercentage[iTest] )
                        printf("\nINFO: for scan #%d the image rotation width in the header (%.2f deg)\n      shows no overlaps\n",
                                 iScan, m_fUserImageRotationWidth);
                    else
                        printf("\nINFO: for scan #%d the image rotation width in the header (%.2f deg)\n      shows %.1f%% overlaps\n",
                                 iScan, m_fUserImageRotationWidth, vecOverlapResults[iScan].afOverlapsPercentage[iTest]);
                    break; 
                } 
            }
        }
    }

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function takes all strategy solution sets generated by repeated calls to DTMainStrategy::Execute() and
// writes them out to the summary table. 
void CDTMainMultiStrategy::vPrintSummaryScanTable()
{   
    if( 0 == m_aoScanSolutionSet.size() )
        return; // nothing to do
    
    int                 nCurScanOffset = 0;
    CstrategyScans*     pScans = NULL;
    int                 nScansUsed = 0;
    
    bool                bPrintHeader = false;
    double              dTotalRotWidth = 0.0;
    m_dAllDetectorPositionsAllScansTotalRotWidth = 0.0;
    std::vector<Cstring> padInfoVec;
    for(int nSet = 0; nSet < (int)m_aoScanSolutionSet.size(); nSet++)
    {
        bPrintHeader = nSet == 0 ? true : false;
        
        pScans = m_aoScanSolutionSet[nSet];
        
        if( 0 == nSet )
            nScansUsed = pScans->nPrintResultScanTable(nCurScanOffset, 
                                                       bPrintHeader, 
                                                       dTotalRotWidth, 
													   padInfoVec,
                                                       0U, 
                                                       m_pstPreviousReflnsInfo);
        else
            nScansUsed = pScans->nPrintResultScanTable(nCurScanOffset, bPrintHeader, dTotalRotWidth);

        m_dAllDetectorPositionsAllScansTotalRotWidth += dTotalRotWidth;
        nCurScanOffset += nScansUsed;
    }

    Cstring     sHorizLine("----------------------------------------------------------------------------\n");
    printf(sHorizLine.string());
    
    printf("All scans                                         %8.1f %6.1f%6.2f%5.2f\n", 
                                                              m_dAllDetectorPositionsAllScansTotalRotWidth, 
                                                              m_poResoBins->dGetCompleteness(), 
                                                              m_poResoBins->dGetMultiplicity(),
                                                              m_poResoBins->dGetMultiplicityDev());
    
    printf("\n* Scan rotation axis\n");
    printf("+/- is a standard deviation of redundancy among resolution bins\n");
		for(int i=0; i<padInfoVec.size(); i++)
	{
		printf(padInfoVec.at(i));
		printf("\n");
	}
	padInfoVec.clear();
	printf("\n");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiStrategy::vSetExecutableNameDTSTRATEGY()
{
    Cstring     sCommand("dtstrategy ");
    for(int ii=1; ii < m_nDTREKCommandArgumentCount; ii++)
    {
		const bool addQuotes = !strchr(m_ppcDTREKCommandArguments[ii], '"') &&
			strchr(m_ppcDTREKCommandArguments[ii], ' ');
		if(addQuotes) sCommand += "\"";
        sCommand += m_ppcDTREKCommandArguments[ii];
		if(addQuotes) sCommand += "\"";
        sCommand += ' ';
    }

    vCommandLineToCommandArgumentArray(sCommand.string());
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// While keeping the min and max intact, redistribute the overlaps evenly
void CDTMainMultiStrategy::vRedistributeDetectorOverlaps()
{
    int     nTotalPositions = m_ad2ThetaPositions.size();
    if(  nTotalPositions < 3 )
        return; // nothing to do

    //////////////////////////////////////////////////////////////////////////
    // Calculate the mean overlap value
    int ii = 0;
    double      dOverlap = 0.0;
    for(ii=0; ii < nTotalPositions - 1; ii++)
    {
        dOverlap += fabs(m_ad2ThetaPositions[ii] - m_ad2ThetaPositions[ii+1]);    
    }
        
    dOverlap /= (nTotalPositions - 1); 
    //////////////////////////////////////////////////////////////////////////

    // Now set the mean overlap to the affected positions
    for(ii=1; ii < nTotalPositions - 1; ii++)
    {
        m_ad2ThetaPositions[ii] = m_ad2ThetaPositions[0] + dOverlap * ii;    
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// If the user supplied a previously generated reflection list, 
// load it here and calculate the resolution as well as reduced HKL's for each reflection.
bool CDTMainMultiStrategy::bLoadPreviousReflections(Cstring& sError)
{
    sError = "";

    Cstring                 sPrevList("");     
    if( !bGetCommandOption(c_pcPrev, sPrevList, false) )
        return false;
    
    if( "" == sPrevList )
        return false; // nothing to do
    
    char    acTemp[DTREK_MAX_LINE_LENGTH];

    //////////////////////////////////////////////////////
    // Read in the previous reflection list.
    if( m_poPreviousReflnList )
        delete m_poPreviousReflnList;
        
    m_poPreviousReflnList = new Creflnlist();

    if( 0 != m_poPreviousReflnList->nRead(sPrevList) )
    {
        sprintf(acTemp, "Error reading reflection file: %s.", sPrevList.string());
        sError = acTemp;
        
        delete m_poPreviousReflnList;
        m_poPreviousReflnList = NULL;

        return false;
    }
    
    if( 0 == m_poPreviousReflnList->nGetNumReflns() )
    {
        printf("/nWARNING: Reflection list in file %s is empty./n", sPrevList.string());

        return false; // nothing to do
    }
    //////////////////////////////////////////////////////

    vExpandReflnList(m_poPreviousReflnList);
    
    double  fTemp1=0.0, fTemp2=0.0;
    m_oMainStrategy.vDoExternalReflnListCalcSetDStarCubed(m_poPreviousReflnList, fTemp1, fTemp2);
    
    m_oMainStrategy.vDoExternalReflnListSortReduce(m_poPreviousReflnList);

    /////////////////////////////////////////////////////////////////////////////
    // Do statistics on the "previous list".
    bDoExternalResoBinsStatistics(c_nNumberOfResoBins, m_poPreviousReflnList, sError);

    if( m_pstPreviousReflnsInfo )
        delete m_pstPreviousReflnsInfo;
    
    m_pstPreviousReflnsInfo = new REFLN_FILE_STRAT_INFO(sPrevList);
    m_pstPreviousReflnsInfo->m_fCompleteness = m_poResoBins->dGetCompleteness();
    m_pstPreviousReflnsInfo->m_fRedundancy   = m_poResoBins->dGetMultiplicity();
    m_pstPreviousReflnsInfo->m_fRedundancyDev = m_poResoBins->dGetMultiplicityDev();
    /////////////////////////////////////////////////////////////////////////////////

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function detects the autodist option on the command line, and generates automatic distance.
// NOTE: if the function returns false, the only case when sError can be empty is if autodist option is not supplied at all.
bool CDTMainMultiStrategy::bSetAutoDistance(Cstring& sError)
{
    char    acTemp[DTREK_MAX_LINE_LENGTH];

    sError = "";
    m_dDistance = -1.0; // initialize
    
    Cstring                 sArgs("");     

    if( !bGetCommandOption(c_pcAutoDistance, sArgs) )
    {
        return false; // autodist option is not supplied    
    }

    printf("\nSetting automatic detector distance...\n");


    /////////////////////////////////////////////////////////////////////////
    // Check if -distance is supplied.
    Cstring     sArgs1("");
    if( bGetCommandOption(c_pcSetDistance, sArgs1, false) )
    {
        sprintf(acTemp, "Options %s and %s are incompatible.", c_pcAutoDistance, c_pcSetDistance);
        sError = acTemp;
        
        return false;
    }
    /////////////////////////////////////////////////////////////////////////

    // Figure out parameters passed with the auto distance option
    if( "" != sArgs )
    {
        if( !sArgs.bIsNumeric() )
        {
            sprintf(acTemp, "Option %s has a non-numeric argument.", c_pcAutoDistance);
            sError = acTemp;

            return false;
        }
        else
        {
            m_nUserSpotSepPix = atoi(sArgs.string());
            if( m_nUserSpotSepPix < 0 )  // sanity check
                m_nUserSpotSepPix *= -1; // could return an error condition, of course, but let's help the user...
        }    
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Now generate automatic distance    
    if( !bGenerateValidDistance(sError) )
    {
        if( sError == "" )
            sError = "Cannot generate automatic distance.";
        return false;
    }
    
    return true;
}    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bGenerateValidDistance(Cstring& sError)
{
    if( !m_poHeader )
        return false;   // safety
    
    sError = "";
    
    if( !bReadSetHardwareLimits(sError, DTREK_MULTISTRATEGY_DET_DISTANCE) )  // This call should set distance limits. If it fails, that's it.
    {
        if( "" == sError )
            sError = "Cannot generate automatic distance. Limits are not specified.";
        return false;
    }

    CExperimentalStrategy       oExpStrategy(*m_poHeader);
    m_dDistance = oExpStrategy.fCalcBestDistanceSpotSeparation(m_nUserSpotSepPix,
                                                               CExperimentalStrategy::eComputeDistanceSpotSeparationDefault,
                                                               CExperimentalStrategy::eComputeDistanceSpotSeparationDefault,
                                                               &m_stSpotSizeInfo);
    if( 0.0 > m_dDistance )
    {
        sError = "Failed to generate crystal-to-detector distance automatically.\n Crystal information in the header might be invalid.\n";
        return false;
    }

    // Get the nearest multiple of 2 mm
    m_dDistance = dGetNearestMultipleInteger(m_dDistance, 2);
    
    // Clamp distance
    m_segUserDistLimits.vClamp(m_dDistance);

    // UPDATE HEADER TO INCLUDE NEW DISTANCE
    
    // Get Prefix
    Cstring              strDetectorPrefix(""); 
    if ( 0 != m_poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    // Get Gonio
    Cgoniometer         oDetGonio(*m_poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
        return false;
    
    // Set Distance
    oDetGonio.nSetDistance((float)m_dDistance);
    
    // Update Header
    oDetGonio.nUpdateHeader(m_poHeader, strDetectorPrefix);
   
    //////////////////////////////////////////////////////////////////////////
    // Now report the generated distance
    printf("\nAutomatically generated distance: %.1f mm.\n", m_dDistance);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bSetUserInputDistance(Cstring& sError)
{
    char    acTemp[DTREK_MAX_LINE_LENGTH];

    sError = "";
    m_dDistance = -1.0; // initialize
    
    Cstring                 sArgs("");     
    std::vector<double>     aParams;
    int                     nParamsCount=0;

    if( !bGetCommandOption(c_pcSetDistance, sArgs) )
    {
        return false; // distance option is not supplied    
    }

    /////////////////////////////////////////////////////////////////////////
    // Just a consistency check....
    Cstring     sArgs1("");
    if( bGetCommandOption(c_pcAutoDistance, sArgs1, false) )
    {
        sprintf(acTemp, "Options %s and %s are incompatible.", c_pcAutoDistance, c_pcSetDistance);
        sError = acTemp;

        return false;
    }
    /////////////////////////////////////////////////////////////////////////

    // Figure out parameters passed with the distance option
    sArgs.nListToVector(aParams, " ");
    nParamsCount = aParams.size();
    
    if( nParamsCount != 1 )
    {
        sprintf(acTemp, "Option %s has no parameters.", c_pcSetDistance);
        sError = acTemp;

        return false;
    }
    else if( aParams[0] < 0.0 )
    {
        sprintf(acTemp, "Option %s has a negative distance.", c_pcSetDistance);
        sError = acTemp;

        return false;
    }
    
    // Set the distance to class member
    m_dDistance = aParams[0];
    
    // UPDATE HEADER TO INCLUDE NEW DISTANCE
    Cstring              strDetectorPrefix(""); 
    if ( 0 != m_poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    // GET DISTANCE
    Cgoniometer         oDetGonio(*m_poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
        return false;
    
    oDetGonio.nSetDistance((float)m_dDistance);
    oDetGonio.nUpdateHeader(m_poHeader, strDetectorPrefix);

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bReadSetHardwareLimits(Cstring& sError, DTREK_WORD wCtrl)
{
    if( !m_poHeader )
        return false;  //safety
    
    char    acTemp[DTREK_MAX_LINE_LENGTH];
    sError = "";

    // First try and get the limits from the header
    
    CSegment        segHeaderLimits(cdStratDetUnrealisticValue, cdStratDetUnrealisticValue);   // initialize to unrealistic numbers

    Cstring              strDetectorPrefix(""); 
    if ( 0 != m_poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
    {
        sprintf(acTemp, "Cannot get detector name prefix from the header.");
        sError = acTemp;
        
        return false;
    }

    Cgoniometer         oDetGonio(*m_poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
    {    
        sprintf(acTemp, "Cannot create gonio object from the header.");
        sError = acTemp;
        
        return false;
    }

    // OK, we're ready to get the limits!
    float     fMin = 0.0f;
    float     fMax = 0.0f;
    
    bool      bHeaderHasLimits = false;
    if( wCtrl == DTREK_MULTISTRATEGY_DET_DISTANCE )
        bHeaderHasLimits = oDetGonio.bGetHardwareLimits(c_pcHeaderDistSearchKey, fMin, fMax);
    else if(  wCtrl == DTREK_MULTISTRATEGY_DET_2THETA )
    {
        bHeaderHasLimits = oDetGonio.bGetHardwareLimits(c_pcHeader2ThetaSearchKey, fMin, fMax);
        
        if( bHeaderHasLimits )
        {
            // We want the 2Theta values to be positive. If the actual 2Theta limits are negative, we will still make them positive
            // and memorize the fact that they are really negative. Later, when the sign matters, we will set the correct 2Theta sign.
            // The complication is that we may have ranges like (-5,110) or (-110, 5), but we really want to treat them as (0, 110) and (-110, 0)
            bool        bMaxValueHasTheLargerAbsoluteValue = fabs(fMax) > fabs(fMin);
        
            if( fMin < 0.0 && !bMaxValueHasTheLargerAbsoluteValue )  // a case like (-110, 5) or (-110, 0) or even (-45, 45)
            {
                std::swap(fMin, fMax);
                fMin = -fMin;
                fMax = -fMax;
                m_n2ThetaDirectionSign = -1;
            }                                                        // else we may have a case like (-5, 110) or (0, 110)
        
            fMin = 0.0f; // regardless of what the header says, we always set min 2Theta to 0
        }
    }

    if( bHeaderHasLimits )
        segHeaderLimits.vSet((double)fMin, (double)fMax);

    // Now, let's see if the user supplied the limits on the command line.
    bool      bCommandLineHasLimits = false;

    CSegment  segCommandLineLimits(cdStratDetUnrealisticValue, cdStratDetUnrealisticValue);   // initialize to unrealistic numbers
    
    Cstring                 sArgs("");     
    std::vector<double>     aParams;
    int                     nParamsCount=0;
             
    Cstring     strCommandLineLimitsOption("");
    if( wCtrl == DTREK_MULTISTRATEGY_DET_DISTANCE )
        strCommandLineLimitsOption = c_pcDistanceLimits;
    else if(  wCtrl == DTREK_MULTISTRATEGY_DET_2THETA )
        strCommandLineLimitsOption = c_pc2ThetaMax;
    
    if( bGetCommandOption(strCommandLineLimitsOption.string(), sArgs) )
    {
        sArgs.nListToVector(aParams, ' ');
        nParamsCount = aParams.size();

        int     nExpectedNumberOfArguments = (wCtrl == DTREK_MULTISTRATEGY_DET_2THETA) ? 1 : 2; // For 2Theta the mimimum is assumed to be 0.0, so we only expect 1 argument
        
        
        if( nParamsCount != nExpectedNumberOfArguments )
        {
            sprintf(acTemp, "Option %s must have %d arguments.", strCommandLineLimitsOption.string(), nExpectedNumberOfArguments);
            sError = acTemp;
        
            return false;
        }
        
        if( wCtrl == DTREK_MULTISTRATEGY_DET_DISTANCE )
        {
            
            if( aParams[0] < 0.0 || aParams[1] < 0.0 )
            {
                sprintf(acTemp, "Option %s has one or two negative arguments: %f %f", strCommandLineLimitsOption.string(),
                                                                                      aParams[0], aParams[1]);
                sError = acTemp;
                return false;
            }
            segCommandLineLimits.vSet(aParams[0], aParams[1]);
        }
        else if( wCtrl == DTREK_MULTISTRATEGY_DET_2THETA )
        {
            if( aParams[0] < 0.0 )
            {
                aParams[0] = -aParams[0];
                m_n2ThetaDirectionSign = -1;
            }
            segCommandLineLimits.vSet(0.0, aParams[0]);
        }
        
        bCommandLineHasLimits = true;
    }

    // Combine command line limits with header limits, giving preference to header limits 
    if( bHeaderHasLimits )
    {
        if( wCtrl == DTREK_MULTISTRATEGY_DET_DISTANCE )
        {
            if( bCommandLineHasLimits )
            {
                CSegment    segTemp = segHeaderLimits & segCommandLineLimits;
                m_segUserDistLimits = segTemp.bIsEmpty() ? segHeaderLimits : segTemp;
            }
            else
            {
                m_segUserDistLimits   = segHeaderLimits;
            }
        }
        else if( wCtrl == DTREK_MULTISTRATEGY_DET_2THETA )
        {
            if( bCommandLineHasLimits )
            {
                CSegment    segTemp = segHeaderLimits & segCommandLineLimits;
                m_segUser2ThetaLimits = segTemp.bIsEmpty() ? segHeaderLimits : segTemp;
            }
            else
            {
                m_segUser2ThetaLimits = segHeaderLimits;
            }
        }
    }
    else if( bCommandLineHasLimits )
    {
        if( wCtrl == DTREK_MULTISTRATEGY_DET_DISTANCE )
            m_segUserDistLimits = segCommandLineLimits;
        else if( wCtrl == DTREK_MULTISTRATEGY_DET_2THETA )
            m_segUser2ThetaLimits = segCommandLineLimits;
    }
    else
        return false;  // Neither header nor command line provide limits
    
    // We need to make sure that the hardware limits are only accurate to the integer digits
    if( wCtrl == DTREK_MULTISTRATEGY_DET_DISTANCE )
        m_segUserDistLimits.bRoundOffTowardCenter(0);   // 0 - number of decimal digits to keep
    else if( wCtrl == DTREK_MULTISTRATEGY_DET_2THETA )
        m_segUser2ThetaLimits.bRoundOffTowardCenter(0); // 0 - number of decimal digits to keep

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiStrategy::vSetDetectorPositionToHeader(DTREK_WORD wCtrl)
{
    if( !m_poHeader )
        return; // safety
    
    Cstring              strDetectorPrefix(""); 
    if ( 0 != m_poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return;

    Cgoniometer         oDetGonio(*m_poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
        return;
    
    if( wCtrl & DTREK_MULTISTRATEGY_DET_DISTANCE && m_dDistance > 0.0 )
    {
        printf("\nSetting detector distance = %.1f mm to header.\n", m_dDistance);
        oDetGonio.nSetDistance((float)m_dDistance);
    }
    
    if( wCtrl & DTREK_MULTISTRATEGY_DET_2THETA && m_nActive2ThetaPositionIndex >= 0 && 
        m_nActive2ThetaPositionIndex < (int)m_ad2ThetaPositions.size() )
    {
        printf("\nSetting detector 2Theta = %.1f deg to header.\n\n", m_ad2ThetaPositions[m_nActive2ThetaPositionIndex]);
        oDetGonio.nSetSwing((float)m_ad2ThetaPositions[m_nActive2ThetaPositionIndex]);
    }

    oDetGonio.nUpdateHeader(m_poHeader, strDetectorPrefix);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function gets detector position for the first input header scan (scan with index 0)
void CDTMainMultiStrategy::vGetDetectorPositionFromHeader(DTREK_WORD wCtrl)
{
    if( !m_poHeader )
        return; // safety
    
    Cstring              strDetectorPrefix(""); 
    if ( 0 != m_poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return;

    Cgoniometer         oDetGonio(*m_poHeader, strDetectorPrefix);
    if( !oDetGonio.bIsAvailable() )
        return;
    
    double      dDetRelZero[3] = {0.0, cdStratDetUnrealisticValue, cdStratDetUnrealisticValue};
    Cstring     sScanDetRelZero(D_K_ScanDetDatum);
    Cstring     sScanPrefix(D_K_ScanPrefix);
    m_poHeader->nGetValue(sScanDetRelZero, 3, dDetRelZero); // the first value should be equal to 2, then two theta and distance

    // If the hardware ("relative zero") value is available, use it. Else use the "detector gonio" value. 
    // Of course, the detector gonio value will most likely be a refined value, but
    // the assumtion is that the difference is negligible. Otherwise, an offset must already be entered into the configuration file.
    double      dValue = 0.0;
    if( wCtrl & DTREK_MULTISTRATEGY_DET_2THETA )
    {
        if( cdStratDetUnrealisticValue != dDetRelZero[1] )
            dValue = dDetRelZero[1];
        else
        {
            dValue = (double)oDetGonio.fGetSwing();
            dValue = dRoundOff(dValue, 0); // for 2theta it's enough to just keep the integer part
        }
    
        m_d2ThetaFromHeader = dValue;
    }
    
    if( wCtrl & DTREK_MULTISTRATEGY_DET_DISTANCE )
    {
        if( cdStratDetUnrealisticValue != dDetRelZero[2] )
            dValue = dDetRelZero[2];
        else
        {
            dValue = (double)oDetGonio.fGetDistance();
            dValue = dRoundOff(dValue, 1); // let's keep one decimal digit
        }
        
        m_dDistance = dValue;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function removes all strategy scans except for the original scan with index 0. 
// The latter must stay in the header, because it carries information about the rotation axis, original detector and crystal position etc.
void CDTMainMultiStrategy::vRemoveMultiScanInfoFromHeader()
{
    ///////////////////////////////////////////////////////////////////////////////////////////
	// Clear multiscan information (i.e. those entries prefixed with "S*", where * is a number)
	
    ///////////////////////////////////////////////////////////////////
    // Old style
    m_poHeader->nDeleteMask("S*_SCAN_DET_RELZERO");
	m_poHeader->nDeleteMask("S*_SCAN_ROTATION");
	m_poHeader->nDeleteMask("S*_SCAN_ROTATION_AXIS_NAME");
	m_poHeader->nDeleteMask("S*_SCAN_ROTATION_LIMITS");
	m_poHeader->nDeleteMask("S*_SCAN_ROTATION_VECTOR");
	m_poHeader->nDeleteMask("S*_CRYSTAL_GONIO_VALUES");
	///////////////////////////////////////////////////////////////////
	
    /////////////////////////////////////////////////////////////////////////////
    // New style 
    m_poHeader->nDeleteMask(D_K_StrategyTestPrefix"S*_SCAN_DET_RELZERO");
	m_poHeader->nDeleteMask(D_K_StrategyTestPrefix"S*_SCAN_ROTATION");
	m_poHeader->nDeleteMask(D_K_StrategyTestPrefix"S*_SCAN_ROTATION_AXIS_NAME");
	m_poHeader->nDeleteMask(D_K_StrategyTestPrefix"S*_SCAN_ROTATION_LIMITS");
	m_poHeader->nDeleteMask(D_K_StrategyTestPrefix"S*_SCAN_ROTATION_VECTOR");
	m_poHeader->nDeleteMask(D_K_StrategyTestPrefix"S*_CRYSTAL_GONIO_VALUES");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Clear multstrategy output that might be left over from a previous run
void CDTMainMultiStrategy::vClearMultiStrategyOutputFromHeader()
{
    m_poHeader->nDeleteMask(D_K_StrategyResultPrefix"S*_SCAN_DET_RELZERO");
	m_poHeader->nDeleteMask(D_K_StrategyResultPrefix"S*_SCAN_ROTATION");
	m_poHeader->nDeleteMask(D_K_StrategyResultPrefix"S*_SCAN_ROTATION_AXIS_NAME");
	m_poHeader->nDeleteMask(D_K_StrategyResultPrefix"S*_SCAN_ROTATION_LIMITS");
	m_poHeader->nDeleteMask(D_K_StrategyResultPrefix"S*_SCAN_ROTATION_VECTOR");
	m_poHeader->nDeleteMask(D_K_StrategyResultPrefix"S*_CRYSTAL_GONIO_VALUES");
    
    m_poHeader->nDeleteMask(D_K_ScanSeqSelected);
    m_poHeader->nDeleteMask(D_K_DtstrategyPercentComplete);
    m_poHeader->nDeleteMask(D_K_DtstrategyRedundancy);
    m_poHeader->nDeleteMask(D_K_DtstrategyRedundancyDev);

    // Clear ranking info in the header
    m_poHeader->nClearDictionary("DTMULTISTRATEGY_RANK", "*");
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bSetTestScansFromCommandLineToHeader(Cstring& sError)
{
    sError = "";            // initialize
    
    m_bTestScansGeneratedFromCommandLine = false;
    
    char    acTemp[DTREK_MAX_LINE_LENGTH];

    // Read command line and fill out an array of pairs: "Axis Name" & "Axis Values"
    Cstring                     sParams("");
    std::vector<Cstring>        asParams;
    
    std::vector< std::pair< Cstring, std::vector<double> > >           vecCrysGonioNamesAndValues;

    std::vector<double>         adValues;

    int     ii = 0;
    int     jj = 0;

    while( bGetCommandOption(c_pcCrystalGonio, sParams) )
    {
        sParams.nListToVector(asParams, " ,");
        
        if( asParams.size() < 2 )
        {
            sprintf(acTemp, "Option %s must have an axis name and at least one value", c_pcCrystalGonio);
            sError = acTemp;
            
            return false;
        }
        
        adValues.clear();
        for(ii=1; ii < (int)asParams.size(); ii++)  // the first element of asParams must be an angle name; the rest must be that angle's values
        {
            adValues.push_back( atof( asParams[ii].string() ) );
        }
        
        vecCrysGonioNamesAndValues.push_back(std::pair<  Cstring, std::vector<double> > (asParams[0], adValues) );
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( 0 == vecCrysGonioNamesAndValues.size() )
        return false;  // Nothing to do: no scan information on the command line.
    else
        m_bTestScansGeneratedFromCommandLine = true;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // For now let's support only 2 or 1 non-rotation axes ("chi(kappa) and phi" or just phi)
    if( vecCrysGonioNamesAndValues.size() != 2 && vecCrysGonioNamesAndValues.size() != 1 )
    {
        sError = "Number of crystal gonio axes on the command line must be 1 or 2";
        return false;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Create an output scans object to hold the scans we need to generate.
    // Initialize oOutputScans from the header. Actually, all we need for the initialization is a scan with index 0 in the header.
    CstrategyScans          oOutputScans(m_poHeader,"");  
    if( 0 == oOutputScans.nGetNumScans() )
    {
        sError = "Cannot set scans. The input header must contain at least one scan.";
        return false;
    }
    
    int     nNumAxes = oOutputScans.nGetNumAxes();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Now let's clear all scans in OutputScans object. We will be adding new scans based on the command-line input _and_ scan with index 0 in the header.
    oOutputScans.vClearScans();  
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int     nCurScanIndex = -1;
    if( vecCrysGonioNamesAndValues.size() == 2 ) // chi(kappa) AND phi 
    {
        for(ii=0; ii < (int)((vecCrysGonioNamesAndValues[0].second).size()); ii++)
        {
            for(jj=0; jj < (int)((vecCrysGonioNamesAndValues[1].second).size()); jj++)
            {
                oOutputScans.vAddScan(m_poHeader, 0, nNumAxes);   // 0, because we base each new scan on index-zero scan, present in header

                nCurScanIndex = oOutputScans.nGetNumScans() - 1;  // the index of the scan we have just added

                oOutputScans.vSetScanUsed(nCurScanIndex, CstrategyScan::enMultipleScan, true); // make sure it is set as used, so that it can be written to the header

                if( !oOutputScans.bSetScanAxisValue(nCurScanIndex, (vecCrysGonioNamesAndValues[0].first).string(), 
                                                                   (vecCrysGonioNamesAndValues[0].second)[ii]) )
                {
                    sprintf(acTemp, "Cannot set axis %s value to header. Please check axis name.", 
                                    (vecCrysGonioNamesAndValues[0].first).string());
                    sError = acTemp;
                    return false;
                }

                if( !oOutputScans.bSetScanAxisValue(nCurScanIndex, (vecCrysGonioNamesAndValues[1].first).string(), 
                                                                   (vecCrysGonioNamesAndValues[1].second)[jj]) )
                {
                    sprintf(acTemp, "Cannot set axis %s value to header. Please check axis name.", 
                                    (vecCrysGonioNamesAndValues[1].first).string());
                    sError = acTemp;
                    return false;
                }
            }
        }
    }
    else if( vecCrysGonioNamesAndValues.size() == 1 ) //  chi(kappa) OR phi
    {
        for(jj=0; jj < (int)((vecCrysGonioNamesAndValues[0].second).size()); jj++)
        {
            oOutputScans.vAddScan(m_poHeader, 0, nNumAxes); // 0, because we base each new scan on index-zero scan, present in header

            nCurScanIndex = oOutputScans.nGetNumScans() - 1;  // the index of the scan we have just added

            oOutputScans.vSetScanUsed(nCurScanIndex, CstrategyScan::enMultipleScan, true); // make sure it is set as used, so that it can be written to the header

            if( !oOutputScans.bSetScanAxisValue(nCurScanIndex, (vecCrysGonioNamesAndValues[0].first).string(), 
                                                               (vecCrysGonioNamesAndValues[0].second)[jj]) )
            {
                sprintf(acTemp, "Cannot set axis %s value to header. Please check axis name.", 
                                (vecCrysGonioNamesAndValues[0].first).string());
                sError = acTemp;
                return false;
            }
        }
    }

    // Remove all multi scan info from the header
    vRemoveMultiScanInfoFromHeader();

    // Put the new scans into the header
    oOutputScans.vSetMultipleScanSolution(true);  // so that _all_ scans are taken into account
    
    oOutputScans.nUpdateHeader(m_poHeader, 0, D_K_StrategyTestPrefix, false);

    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The purpose of this function is to make sure that the image rotation and 
// scan rotation values in the header 
// (a) are consistent with the corresponding hardware limits in the header
// (b) cover a reasonably-narrow range 
bool CDTMainMultiStrategy::bInsureHeaderRotLimitsConsistency(Cstring& sError)
{
    sError = "";
    if( !m_poHeader )
        return false;  //safety check

    Crotation       oRotation(*m_poHeader);
    if( !oRotation.bIsAvailable() )
        return false;

    
    Cstring     strRotAxisName("");
    strRotAxisName = oRotation.sGetName();
 
    float       fRotAxisHardwareMin = -1.0f;
    float       fRotAxisHardwareMax = -1.0f;
    if( !bGetReasonableCrystalGonioHardwareLimits(strRotAxisName, 
                                                  fRotAxisHardwareMin, 
                                                  fRotAxisHardwareMax,
                                                  true,
                                                  sError) )
    {
        printf("\n\nWARNING: The image header does not contain information about the rotation axis hardware limits.\n\n");
        //to-do: handle errors, use sError
        return false;
    }
   
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Update the image rotation limits
    oRotation.vSetRotMinMaxRange(fRotAxisHardwareMin, fRotAxisHardwareMax, fRotAxisHardwareMax-fRotAxisHardwareMin);
    
    if( 0 != oRotation.nUpdateHeader(m_poHeader) )
        return false;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Update the scan rotation limits
    Cstring         sScanPrefix(D_K_ScanPrefix);
    Crotation       oScanRotation(*m_poHeader, sScanPrefix);
    if( !oScanRotation.bIsAvailable() )
        return false;
    
    oScanRotation.vSetRotMinMaxRange(fRotAxisHardwareMin, fRotAxisHardwareMax, fRotAxisHardwareMax-fRotAxisHardwareMin);

    if( 0 != oScanRotation.nUpdateHeader(m_poHeader, sScanPrefix) )
        return false;
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiStrategy::vUpdateTargetCompletenessTolerance(double fCompleteness, double fTolerance)
{
    Cstring     sTemp(fCompleteness);
    sTemp += ' ';
    sTemp += fTolerance;

    bUpdateCommandOption(c_pcCompleteness, sTemp);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiStrategy::vSaveTargetCompletenessTolerance()
{
    std::vector<double>     afTemp;
    if( bGetCommandOption(c_pcCompleteness, afTemp," ,", false) )
    {
        if( afTemp.size() > 0 ) 
            m_fUserTargetCompleteness = fabs(afTemp[0]);  // fabs is just for safety
        
        if( afTemp.size() > 1 )
            m_fUserCompletnessTolerance = fabs(afTemp[1]);
    }
    
    if( m_fUserCompletnessTolerance > m_fUserTargetCompleteness ) // sanity check
       m_fUserCompletnessTolerance = m_fUserTargetCompleteness;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is called if the user didn't specify scans to check on the command line,
// It could be that the user has set up some scans with the correct keywords in the header -
// then we have nothing to do. If they didn't or set them up with the old keywords, we have
// to insure the scans are using the new keywords.
bool CDTMainMultiStrategy::bInsureTestScansInHeaderAreInNewFormat(Cstring& sError)
{
    m_bTestScansGeneratedFromOldHeaderFormat = false;

    sError = "";

    Cstring     sTemp("");
    
    Cstring     sFirstScanRotation = sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, 0);
    sFirstScanRotation += D_K_Rotation;

    if( 0 == m_poHeader->nGetValue(sFirstScanRotation, &sTemp) )
        return true; // the header does have at least one scan entry with the new keyword used

    // Else - read old keyword(s) and generate new keyword entries
    
    CstrategyScans          oOldFormatScans(m_poHeader, "");
    
    if( 0 == oOldFormatScans.nGetNumScans() )
    {
        sError = "Cannot set test scans. The input header must contain at least one scan.";

        return false;
    }

    // Remove all multi scan info from the header
    vRemoveMultiScanInfoFromHeader();

    // Put the scans in the new format into the header
    
    // Make sure _all_ scans are taken into account
    oOldFormatScans.vSetMultipleScanSolution(true);  
    
    // Make sure all scans are used
    for(int ii=0; ii < oOldFormatScans.nGetNumScans(); ii++)
        oOldFormatScans.vSetScanUsed(ii, CstrategyScan::enMultipleScan, true);

    // Now we are ready to update the header
    oOldFormatScans.nUpdateHeader(m_poHeader, 0, D_K_StrategyTestPrefix, false);
    
    m_bTestScansGeneratedFromOldHeaderFormat = true;

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiStrategy::vSaveTestScansFromHeader()
{
    m_paoSaveTestScans = new CstrategyScans(m_poHeader, D_K_StrategyTestPrefix);
    
    m_paoSaveTestScans->vSetMultipleScanSolution(true);
    
    for(int ii=0; ii < m_paoSaveTestScans->nGetNumScans(); ii++)
        m_paoSaveTestScans->vSetScanUsed(ii, CstrategyScan::enMultipleScan, true);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Get rotation limits from command line, validate them, set them to the header
bool CDTMainMultiStrategy::bApplyRestrictedRotationLimitsToScansInHeader(Cstring& sError)
{
    sError = "";
   
    std::vector<double>     afArgs;
    CSegment                segRotLimitsFromCommandLine;
    int                     nChosenScan = -1;
    CSegment                segRotLimitsToBeUsed;
    CSegment                segHardwareRotLimits;
    bool                    bHardwareLimitsAvailable = false;
    bool                    bAlreadyTriedToGetHardwareLimits = false;
    
    float                   fAxisHardwareMin = -1.0f;
    float                   fAxisHardwareMax = -1.0f;
    
    while( bGetCommandOption(c_pcRotationLimits, afArgs, " ") )
    {
        if( afArgs.size() < 2 )
        {
            char    acTemp[DTREK_MAX_LINE_LENGTH];
            sprintf(acTemp, "Option %s must have at least 2 arguments.", c_pcRotationLimits);

            sError = acTemp;

            return false;
        }
    
        segRotLimitsFromCommandLine.vSet(afArgs[0], afArgs[1]);
    
        nChosenScan = -1;
        if( afArgs.size() > 2 )
            nChosenScan = (int)afArgs[2];
    
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( !bAlreadyTriedToGetHardwareLimits )
        {
            if( !m_poHeader )
                return false;  //safety check

            //////////////////////////////////////////////////////////////////////////////
            // Try and get the hardware limits and figure out which limits to use:
            // from the hardware or from command line.
    
            // Build the rotation object to get the rotation axis name
            Crotation   oRotation(*m_poHeader);
            if( !oRotation.bIsAvailable() )
            {
                sError = "Failed to construct rotation object from header.";
                return false;
            }
            Cstring     strRotAxisName = oRotation.sGetName();
    
            

            bHardwareLimitsAvailable = bGetReasonableCrystalGonioHardwareLimits(strRotAxisName,
                                                                                fAxisHardwareMin,
                                                                                fAxisHardwareMax,
                                                                                false,
                                                                                sError);
            segHardwareRotLimits.vSet((double)fAxisHardwareMin, (double)fAxisHardwareMax);

            bAlreadyTriedToGetHardwareLimits = true;
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        // If hardware limits are available, combine them with the limits from the command line.
        // Else use the limits from the command line.
        if( bHardwareLimitsAvailable )
        {
            CSegment        segTemp = segHardwareRotLimits & segRotLimitsFromCommandLine;
            segRotLimitsToBeUsed = segTemp.bIsEmpty() ? segHardwareRotLimits : segTemp;
        }
        else
            segRotLimitsToBeUsed = segRotLimitsFromCommandLine;
        /////////////////////////////////////////////////////////////////////////////
    
        /////////////////////////////////////////////////////////////////////////////
        // Go through the list of scans and set the rotation limits
        CstrategyScans          oScans(m_poHeader, D_K_StrategyTestPrefix);
        if( 0 == oScans.nGetNumScans() )
            return false;   // no scans - nothing to do, although this should never happen

        //////////////////////////////////////////////////////////////////////////////////
        oScans.vSetMultipleScanSolution(true);  
        for(int ii=0; ii < oScans.nGetNumScans(); ii++)
        {
            oScans.vSetScanUsed(ii, CstrategyScan::enMultipleScan, true);      // Make sure _all_ scans are taken into account
        
            if( nChosenScan == ii || nChosenScan == -1 )
                oScans.vSetScanRotLimits(ii, segRotLimitsToBeUsed); 
        }
        ////////////////////////////////////////////////////////////////////

        // Now we are ready to update the header
        oScans.nUpdateHeader(m_poHeader, 0, D_K_StrategyTestPrefix, false);
    }

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bDoExternalResoBinsStatistics(int nNumberOfBins, Creflnlist* pReflnList, Cstring& sError)
{
    sError = "";

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // First, sort on packed HKL
    int*        pnSortIndex = new int[pReflnList->nGetNumReflns() + 1]; //RB: not sure yet why we need one extra element in the sort array, but Creflnlist class does that, when it creates the sort array on its own.
    int         nPackedHKLIndex = pReflnList->nGetFieldIndex(Creflnlist::ms_snPackedHKL);
    pReflnList->vSort2(eReflnField_int_type, nPackedHKLIndex, eReflnField_float_type, 0, pnSortIndex);
    pReflnList->vSetExternalSortIndex(pnSortIndex);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate bins for overall resolution range
    if( !bGenerateResoBins(m_dUserResoMin, m_dUserResoMax, nNumberOfBins) )
    {
        sError = "Failed to generate resolution bins";

        return false; 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Now do the actual statistics
    m_oMainStrategy.vDoExternalResoBinsStatistics(m_poResoBins, pReflnList);

    m_poResoBins->vCalculateTotalStats();

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This function gets hardware limits for a crystal gonio axis 
//and if limits are too wide, tries to limit the range (otherwise strategy calculations will take forever).
bool CDTMainMultiStrategy::bGetReasonableCrystalGonioHardwareLimits(Cstring& sAxisName,
                                                                    float&   fAxisHardwareMin,
                                                                    float&   fAxisHardwareMax,
                                                                    bool     bOptimize,
                                                                    Cstring& sError)
{
    sError = "";
    Cgoniometer         oGonio(*m_poHeader, Ccrystal::ms_sCrystalPrefix);
    
    if( !oGonio.bIsAvailable() )
        return false;

    fAxisHardwareMin = -1.0f;
    fAxisHardwareMax = -1.0f;

    if( !oGonio.bGetHardwareLimits(sAxisName.string(), fAxisHardwareMin, fAxisHardwareMax) )
        return false;
    
    if( !bOptimize )
        return true; // done

    const   float   c_fMaxReasonableRotationWidth = 400.0f; // This is set a bit larger than 360 degrees on purpose
    if( fAxisHardwareMax - fAxisHardwareMin > 2 * 360.0f )  // so if it is wider than (-360, 360), we need to trim it
    {
        fAxisHardwareMin = max(fAxisHardwareMin, -c_fMaxReasonableRotationWidth);
        fAxisHardwareMax = min(fAxisHardwareMax, c_fMaxReasonableRotationWidth);
        fAxisHardwareMin = max(fAxisHardwareMin, fAxisHardwareMax - c_fMaxReasonableRotationWidth);
    }
    
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainMultiStrategy::bCalculateAutomaticExposureTime()
{
#if 0
    // RB 10/04/07 This is for testing plus we need to eventually switch to using the evaluator class
    CImageExposureTimeCalculator    oCalc;
    double                          dTemp =0.0;
    oCalc.bCalculate(*m_poHeader,
                     //m_dUserResoMin, 
                     0,0,//m_dUserResoMax, 
                     10,
                     1.0,
                     5.0,
                     dTemp,
                     DTREK_IETC_FULL_MERGED_REFLNS);
    Cstring     ss = oCalc.sGetLastError();
#endif

    m_sError = "";

    std::vector<double>     afOptions;
    if( !bGetCommandOption(c_pcAutoExposure, afOptions, " ") )
        return true; // option simply not found on the command-line
    
    if( !m_poHeader )
    {
        m_sError = "Failed to calculate auto-exposure time. Header object is not available";
        return false; // safety
    }

    if( !m_bResolutionSet )
    {
        m_sError = "Failed to calculate auto-exposure time. Target resolution is not available";
        return false;  // Target resolution not known
    }

    // Parse auto-exposure command-line option
    double      fTargetIOverSigmaInTheLastResoShell = 3.0;
    if( afOptions.size() > 0 && afOptions[0] > 0.0 )
        fTargetIOverSigmaInTheLastResoShell = afOptions[0];
    
    double  fSaveTargetIOverSigmaInTheLastResoShell = fTargetIOverSigmaInTheLastResoShell;
    fTargetIOverSigmaInTheLastResoShell *= 4.0; // Fudge factor. We know that I/Sigma after screening is ca 4 times larger
                                                // than the I/Sigma after dtscaleaverage.  We need to work more on this.

    int         nNumberOfResolutionShells = c_nNumberOfResoBins;
    if( afOptions.size() > 1 && afOptions[1] > 0.0 )
        nNumberOfResolutionShells = (int)afOptions[1];
    
    double      fImageRotationWidth = m_fUserImageRotationWidth;
    if( fImageRotationWidth <= 0.0 )
    {
        m_sError = "Failed to calculate auto-exposure time. Image rotation width is not known";
        return false; 
    }

    Cstring              sExposureTimeEvaluationInfo(""); 
    Cstring              sExpTimeHeaderKeyword("DTREFINE_RANK_");
    sExpTimeHeaderKeyword += Cstring(D_K_ExposureTimeEvaluationInfo);
    if ( 0 != m_poHeader->nGetValue(sExpTimeHeaderKeyword, &sExposureTimeEvaluationInfo) )
    {
        m_sError = "Failed to calculate auto-exposure time.\nThe header file does not have necessary information from dtranker (1)";
        return false; 
    }                    

    afOptions.clear();
    sExposureTimeEvaluationInfo.nListToVector(afOptions, " ");
    if( afOptions.size() < 6 )
    {
        m_sError = "Failed to calculate auto-exposure time. The header file does not have necessary information from dtranker (2)";
        return false;
    }

    // We will skip sigmaA and sigmaB for now.
    double      fParam_A = afOptions[0];
    double      fParam_B = afOptions[2];
    double      fParam_ImageRotationWidth  = afOptions[4];
    double      fParam_ImageExposureTime   = afOptions[5];
    double      fParam_IOverSigmaTimePower = afOptions[6];
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // We need to figure out the target resolution from the user resolution range and the number of resolution bins
    // Generate bins for overall resolution range
    if( !bGenerateResoBins(m_dUserResoMin, m_dUserResoMax, nNumberOfResolutionShells) )
    {
        m_sError = "Failed to calculate auto-exposure time. Failed to generate resolution bins from min and max resolution";
        return false; 
    }  

    CResoBin*   poBin = m_poResoBins->pGetBin(nNumberOfResolutionShells-1);
    if( !poBin )
    {
        m_sError = "Failed to calculate auto-exposure time. Failed to get the last resolution bin object";
        return false;
    }

    double      fTargetReso = poBin->fGetMiddleReso();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double      fExtrapolatedIOverSigmaAtTargetResolution = pow(10.0, (fParam_A +
                                                            fParam_B * (0.25/(float)fTargetReso/(float)fTargetReso) ) );

    m_fImageAutoExposureTime = pow(fParam_ImageExposureTime, fParam_IOverSigmaTimePower) *
                                fTargetIOverSigmaInTheLastResoShell/fExtrapolatedIOverSigmaAtTargetResolution;
    m_fImageAutoExposureTime = pow(m_fImageAutoExposureTime, 1.0/fParam_IOverSigmaTimePower);
    m_fImageAutoExposureTime *= fImageRotationWidth/fParam_ImageRotationWidth; // assuming that exposure time 
                                                                               // should grow proportionally 
                                                                               // to the image rotation width
    
    std::vector<int>        anExposureTimeCheckValues;
    anExposureTimeCheckValues.push_back(1);
    anExposureTimeCheckValues.push_back(2);
    anExposureTimeCheckValues.push_back(5);
    anExposureTimeCheckValues.push_back(10);
    anExposureTimeCheckValues.push_back(15);
    anExposureTimeCheckValues.push_back(20);
    anExposureTimeCheckValues.push_back(30);
    anExposureTimeCheckValues.push_back(40);
    anExposureTimeCheckValues.push_back(60);
    anExposureTimeCheckValues.push_back(90);

    m_fImageAutoExposureTime = nGetNearestIntFromVector(m_fImageAutoExposureTime, 
                                                        anExposureTimeCheckValues, 
                                                        DTREK_VEC_ROUNDOFF_CEILING);

    printf("\nAutomatic exposure time evaluation:\n");
    printf("Number of resolution shells: %d\n", nNumberOfResolutionShells);
    printf("Target resolution (mean resolution of the highest resolution shell): %.3f A\n", fTargetReso);
    printf("Target I/Sigma in the highest resolution shell: %.1f\n", fSaveTargetIOverSigmaInTheLastResoShell);
    printf("Image rotation width: %.2f deg\n", fImageRotationWidth);
    printf("Automatically generated exposure time: %.1f sec\n\n", m_fImageAutoExposureTime);
    
    // Save the exposure time evaluation results in the header
    char    cTemp[256];
    sprintf(cTemp, "%.0f %.2f %.2f %.1f", m_fImageAutoExposureTime, 
                                          fImageRotationWidth,
                                          fTargetReso,
                                          fSaveTargetIOverSigmaInTheLastResoShell);

    m_poHeader->nReplaceValueDictionary(Cstring("DTMULTISTRATEGY"), Cstring(D_K_ProposedExposureTime), cTemp);

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Create a rotation object from a header and get the image rotation width
bool CDTMainMultiStrategy::bGetImageRotationWidth()
{   
    m_sError = "";
    m_fUserImageRotationWidth = -1.0;

    double fRotWidth = -1.0;
    if( bGetCommandOption(c_pcImageRotWidth, fRotWidth) )
    {
        if( fRotWidth <= 0.0 )
        {
            m_sError = "Image rotation width on command line is negative or zero.";
            return false;
        }
        
        m_fUserImageRotationWidth = fRotWidth;
    }
    
    if( !m_poHeader )
    {    
        m_sError = "Header object not available";    
        return false; // should not happen
    
    }
    Crotation*  poRotation = new Crotation(*m_poHeader, D_K_ScanPrefix);
    if( !poRotation->bIsAvailable() )
    {
        m_sError = "Failed to construct rotation object from header.";
        delete poRotation;
        return false;
    }
    

    if( m_fUserImageRotationWidth > 0.0 )
    {
        poRotation->vSetIncrement((float)m_fUserImageRotationWidth);
        poRotation->nUpdateHeader(m_poHeader, D_K_ScanPrefix);


		//shijie yao : update the STRATEGY_INPUT_S{0-3}_SCAN_ROTATION entries as well
		//modified 03:25:2008 : retain the *_LIMITS values while only change the rotation width in _ROTATION;
		//						_LIMITS contains hardware limits and may set by user through -rot.
		Cstring sScanPrefix;
		sScanPrefix = Cstring(D_K_StrategyTestPrefix) + Cstring("S0_SCAN_");
        vUpdateStrategyEntryWithImageRotWidth(sScanPrefix, (float)m_fUserImageRotationWidth);

		sScanPrefix = Cstring(D_K_StrategyTestPrefix) + Cstring("S1_SCAN_");
        vUpdateStrategyEntryWithImageRotWidth(sScanPrefix, (float)m_fUserImageRotationWidth);

		sScanPrefix = Cstring(D_K_StrategyTestPrefix) + Cstring("S2_SCAN_");
        vUpdateStrategyEntryWithImageRotWidth(sScanPrefix, (float)m_fUserImageRotationWidth);

		sScanPrefix = Cstring(D_K_StrategyTestPrefix) + Cstring("S3_SCAN_");
        vUpdateStrategyEntryWithImageRotWidth(sScanPrefix, (float)m_fUserImageRotationWidth);
    }
    else
        m_fUserImageRotationWidth = (double)poRotation->fGetIncrement();
    
    delete  poRotation;
    
    return true;
}

void CDTMainMultiStrategy::vUpdateStrategyEntryWithImageRotWidth(Cstring EntryStr, float imageWidth)
{
	Crotation* pStrategyRot = new Crotation(*m_poHeader, EntryStr);
    if(pStrategyRot->bIsAvailable())
    {
	    pStrategyRot->vSetIncrement(imageWidth);
        pStrategyRot->nUpdateHeader(m_poHeader, EntryStr);
    }
	    
    delete pStrategyRot;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////























