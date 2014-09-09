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
// DTMainRanker.cpp  Initial author: RB           15-Mar-2004
//
//  This file contains the member functions of class CDTMainRanker
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

#include "DTMainRanker.h"

#include "DTMainExtractHeader.h"
#include "DTMainFind.h"
#include "DTMainIndex.h"
#include "DTMainRefine.h"
#include "DTMainRank.h"

#include "DTCoreTclParam.h"

#include "Cspatial.h"
#include "Csource.h"
#include "Cscan.h"

#include "CCommandLine.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::ios;
#endif

static char s_strRankingRules[] =
    " Rule  1: Spot count in resolution shells (found spots)\n"
    " Rule  2: I/Sigma in resolution shells (found spots)\n"
    " Rule  3: Spot sharpness\n"
    "  Rule  3a: Spot asymmetry\n"
    "  Rule  3b: Spot intensity per pixel\n"
    "  Rule  3c: Spot anisotropy\n"
    "  Rule  3d: Spot amorphousness\n"
    " Rule  4: Strong ice rings\n"
    " Rule  5: Diffuse ice rings\n"
    " Rule  6: Percentage of spots indexed\n"
    " Rule  7: RMS residual after refinement\n"
    " Rule  8: Mosaicity\n"
    " Rule  9: Percentage of spots refined\n"
    " Rule 10: Spot count in resolution shells (predicted and found spots)\n"
    " Rule 11: I/Sigma in resolution shells (predicted and found spots)\n"
    " Rule 12: Strategy total rotation width\n\n"
    ;

const   Cstring     c_strInputHeaderName          = "dtranker_input.head";
const   Cstring     c_strFindOutputHeaderName     = "dtranker_find.head";
const   Cstring     c_strFindOutputReflnFileName  = "dtranker_find.ref";
const   Cstring     c_strIndexOutputHeaderName    = "dtranker_index.head";
const   Cstring     c_strRefineOutputHeaderName   = "dtranker_refine.head";
const   Cstring     c_strRankOutputHeaderName     = "dtranker.head";

static char*     c_pcSeq       = "-seq";  
static char*     c_pcReso      = "-reso";  
static char*     c_pcResoBins  = "-resobins";  
static char*     c_pcRefList   = "-ref";  
static char*     c_pcOutputHdr = "-out";  
static char*     c_pcSigma     = "-sigma";

static char*     c_pcRank      = "-rank";  
static int       l_nMinNumberOfFilePathChars = 24;
static char*     c_pcResultsFile = "RankSummary.txt";
static char*     c_pcETE = "exptimeeval";
static char*     c_pcRankETE = "-rankete";


class CDTMainRankerEnvironment
{
public:
    CDTMainRankerEnvironment()
    {
        // Temporarily disable nonuniformity info, because it breaks the ice-ring finding mechanism
        m_strSaveEnvCCDNonunfInfo = sGetEnv(Cstring("CCD_NONUNF_INFO"));
        m_strSaveEnvRXNonunfInfo  = sGetEnv(Cstring("RX_NONUNF_INFO"));
        m_strSaveEnvCCDNonunfType = sGetEnv(Cstring("CCD_NONUNF_TYPE"));
        m_strSaveEnvRXNonunfType  = sGetEnv(Cstring("RX_NONUNF_TYPE"));
        
        m_strSaveEnvNonunfOKValue = sGetEnv(Cstring("DTREK_NONUNF_OKVALUE"));
    
        nPutEnv(Cstring("CCD_NONUNF_INFO"), Cstring("None") );
        nPutEnv(Cstring("RX_NONUNF_INFO"), Cstring("None") );
        nPutEnv(Cstring("CCD_NONUNF_TYPE"), Cstring("None") );
        nPutEnv(Cstring("RX_NONUNF_TYPE"), Cstring("None") );
        
        nPutEnv(Cstring("DTREK_NONUNF_OKVALUE"), Cstring(1));
    }

    ~CDTMainRankerEnvironment()
    {
        nPutEnv(Cstring("CCD_NONUNF_INFO"), m_strSaveEnvCCDNonunfInfo);
        nPutEnv(Cstring("RX_NONUNF_INFO"), m_strSaveEnvRXNonunfInfo);
        nPutEnv(Cstring("CCD_NONUNF_TYPE"), m_strSaveEnvCCDNonunfType);
        nPutEnv(Cstring("RX_NONUNF_TYPE"), m_strSaveEnvRXNonunfType);
        
        nPutEnv(Cstring("DTREK_NONUNF_OKVALUE"), m_strSaveEnvNonunfOKValue);
    }
private:
        Cstring             m_strSaveEnvCCDNonunfInfo;
        Cstring             m_strSaveEnvRXNonunfInfo; 
        Cstring             m_strSaveEnvCCDNonunfType;
        Cstring             m_strSaveEnvRXNonunfType; 
        Cstring             m_strSaveEnvNonunfOKValue;
};

CDTMainRanker::CDTMainRanker() :
    m_strImageFilePath(""),
    m_strCurrentHeaderFilePath(""),
    m_strInputHeaderFilePath(""),
    m_strResolutionLimits(""),
    m_strResolutionBins(""),
    m_nTotalNumberOfResoBins(0), 
    m_nNumberOfHighResoBins(-1),  
    m_nNumberOfLowResoBins(-1),   
    m_strTemplate(""),
    m_strBeamCenter(""),       
    m_strDetectorDistance(""), 
    m_strTwoTheta(""),         
    m_strWavelength(""),
    m_strMinPeakCount(""),
    m_strSigmaFIND(""),
    m_strSigmaREFINE(""),
    m_strAnisotropy(""),
    m_strStrategy(""),
    m_strFindOptions(""),
    m_strIndexOptions(""),
    m_strRefineOptions(""),
    m_strRankOptions(""),
    m_strResultsFileName(""),
    m_nNumberOfFilePathChars(l_nMinNumberOfFilePathChars),
    m_strExcludeRankRules(""),
    m_bUseNonunf(true),
    m_strETE(""),
    m_strSampleNameForReport(""),
    m_bRulesHelpOnly(false)
{
}

CDTMainRanker::~CDTMainRanker()
{
}

int CDTMainRanker::nExecute(unsigned int argc, char* argv[])
{
    // RB: THIS DEBUG CODE IS JUST TO STOP THE DTREK PROCESS WHEN STARTED FROM CC:
    //int*  pNull = NULL;
    //*pNull = 1;
    
#ifdef WIN32
  ios::sync_with_stdio(); 
#endif

  vDtrekSetModuleName("dtranker");
  vPrintCopyrightInfo();
  //cout << "\ndtranker:  Copyright (c) 2006 Rigaku\n";
  //cout << D_K_DTREKVersion << endl;

#ifdef SSI_PC
    Cstring sCCAppVersion = (const char*) argv[argc-1];
	if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
		cout << sCCAppVersion;
		argc -= 3; 
		cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
		argv[argc] = NULL;
	}
#endif

    // Copy command line to output log
    cout << "\nCommand line:\n" << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

    argc--; argv++;

    if( 0 == argc )
    {
        DTREK_ERROR(1, "ERROR: dtranker - input image/header file is not specified.");
    }

    if( 1 <= argc && (0==strcmp(argv[0],"-h") || 0==strcmp(argv[0],"-help")) )
    {
        if( 2 <= argc && 0 == strcmp(argv[1], "rules") ) // display rules help
        {
            CDTMainRank     oMainRank;
            Cstring         strRules("");
            oMainRank.vGetRulesDescription(strRules);
            cout << strRules.string() << flush;
            return 0;
        }
        else
        {
            DTREK_ERROR(0, "");  // display general help
        }
    }
    
    argc++; argv--;

    if( false == bCommandArgumentArrayToRankerCommandLineParameters(argc, argv) )
    {
        DTREK_ERROR(1, "ERROR: dtranker - input image/header file is not specified.");
    }

    if( m_bRulesHelpOnly )
    {
        printf("\nHelp on ranking rules has been requested.\n\n");
        
        bDoRank();
        
        return 0;
    }

    if( "" == m_strCurrentHeaderFilePath && !bCreateImageHeaderFile() )
    {
        DTREK_ERROR(2, "ERROR: dtranker - specified image/header file is not available.");
    }

    if( false == bUpdateHeaderFileAccordingToRankerCommandLine() )
    {
        DTREK_ERROR(3, "ERROR: dtranker - cannot update header file from command line options.");
    }
    
    if( false == bGetImageTemplate() )
    {
        DTREK_ERROR(4, "ERROR: dtranker - cannot process image path template information in header.");
    }
    
    if( 0 == m_vecScanImageSequenceRanges.size() )
    {
        if( false == bDeriveImageSequence() )
        {
            DTREK_ERROR(5, "ERROR: dtranker - cannot parse image sequence information in header and/or command line.");
        }
    }

    if ( false == bGenerateSampleNameForReport() )
    {
        DTREK_ERROR(6, "ERROR: dtranker - cannot generate sample name for report.");
    }

    if( "" == m_strResolutionLimits )
    {
        m_strResolutionLimits = " -reso edge";
    }
    else // put "-reso" in front of the provided resolution parameter(s)
    {
        Cstring     strResoSave(m_strResolutionLimits);
        
        m_strResolutionLimits = " -reso " + strResoSave;
    }

    enDoImageProcessing();

    if( false == bDoRank() )
    {
        DTREK_ERROR(7, "ERROR: dtranker - image ranking failed.");
    }

    if( !bOutputResults(12) )
    {
        DTREK_ERROR(8, "ERROR: dtranker - failed to write ranking results to output file.");
    }

    return 0;
}

bool CDTMainRanker::bCommandArgumentArrayToRankerCommandLineParameters(unsigned int argc, char* argv[])
{
    CCoreTclParam       tclParam(argc, argv);

    //////////////////////////////////////////////////////////////////////////
    // Check if the user requested help on the current rules:
    bool                bFound = false;
    m_bRulesHelpOnly = false;
    m_strRankOptions = tclParam.GetParamString ("rules", bFound);
    if( bFound && ( -1 != m_strRankOptions.find("-help") || -1 != m_strRankOptions.find("-HELP") ) )
    {
        m_bRulesHelpOnly = true;
    }
    /////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // We must have either an image name or a header name, unless the command line has a help request on the rules.
    bool                bFoundImageName = false;
    bool                bFoundHeaderName = false;
    m_strImageFilePath =  tclParam.GetParamString ("image", bFoundImageName);
    m_strInputHeaderFilePath =  tclParam.GetParamString ("header", bFoundHeaderName);
    m_strCurrentHeaderFilePath = bFoundHeaderName ? m_strInputHeaderFilePath : Cstring("");

    if( !m_bRulesHelpOnly )
    {
        if( !bFoundImageName && !bFoundHeaderName )
            return false;
        else if( bFoundImageName && bFoundHeaderName )
        {
            cout << "\n\nBoth image and header keywords found on command line. Image keyword will be ignored.\n\n" << flush;
            m_strImageFilePath = "";
        }
    }
    //////////////////////////////////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////////////////////
    // Options don't have to be supplied, but if they still are, let's get'em
    m_strFindOptions = tclParam.GetParamString ("find", bFound);
    m_strIndexOptions = tclParam.GetParamString ("index", bFound);
    m_strRefineOptions = tclParam.GetParamString ("refine", bFound);
    
    //m_strRankOptions = tclParam.GetParamString ("rules", bFound);

    m_strExcludeRankRules = tclParam.GetParamString ("rules_exclude", bFound);

    /////////////////////////////////////////////////////////////////////////
    m_strResolutionLimits = tclParam.GetParamString ("reso", bFound);
    
    m_strResolutionBins = tclParam.GetParamString ("resobins", bFound);
    
    if( bFound )
        vSetResoBinsParams();

    m_strBeamCenter = tclParam.GetParamString ("beampos", bFound);
    m_strDetectorDistance = tclParam.GetParamString ("distance", bFound);
    m_strTwoTheta = tclParam.GetParamString ("2Theta", bFound);
    
    m_strWavelength = tclParam.GetParamString ("wavelength", bFound);
    if( !bFound )
        m_strWavelength = tclParam.GetParamString ("lambda", bFound);  // try old name 

    m_strMinPeakCount = tclParam.GetParamString ("minpeakcount", bFound);
    m_strSigmaFIND = tclParam.GetParamString ("sigmafind", bFound);
    m_strSigmaREFINE = tclParam.GetParamString ("sigmarefine", bFound);

    m_strAnisotropy = tclParam.GetParamString ("aniso", bFound);

    m_strStrategy = tclParam.GetParamString ("strategy", bFound);

    ////////////////////////////////////////////////////////////////////
    m_strResultsFileName = tclParam.GetParamString("results", bFound);
    bParseResultsFileOption();
    ///////////////////////////////////////////////////////////////////

    Cstring         strImageSequence = tclParam.GetParamString ("sequence", bFound);
    if( bFound )
        strImageSequence.nListToVectorOfPairs(m_vecScanImageSequenceRanges, " ,", "-");

    Cstring     strYesNo = tclParam.GetParamString ("nonunf", bFound);
    if( bFound && (strYesNo == "NO" || strYesNo == "no" ) )
        m_bUseNonunf = false;

    m_strETE = tclParam.GetParamString (c_pcETE, bFound);

    return true;
}

bool CDTMainRanker::bCreateImageHeaderFile()
{
    if( "" == m_strImageFilePath )
        return false;

    m_strCurrentHeaderFilePath = (sDtrekGetPrefix() + c_strInputHeaderName);

    Cstring     strExtractHeaderCmd("\"");
    strExtractHeaderCmd += m_strImageFilePath;
    strExtractHeaderCmd += "\"";

    strExtractHeaderCmd += (" " + m_strCurrentHeaderFilePath);

    if ( !bExecuteDTREKCommand(enProcExtractHeader, strExtractHeaderCmd) )
        return false;
    /////////////////////////////////////////////////////////////////////////////

    // Make sure the non-uniformity information in the header refers to the given image.
    // For now assume there is just one detector
    Cimage_header       oHeader(sDtrekGetPrefix() + c_strInputHeaderName);
    if( !oHeader.bIsAvailable() )
        return false;

    Cstring     strKeyDetectorName("");
    if ( 0 != oHeader.nGetValue(Cstring(D_K_DetectorNames), &strKeyDetectorName) )
        return false;

    Cstring     strKeyNonunfInfo =  strKeyDetectorName + D_K_NonunfInfo;
    if( 0 != oHeader.nReplaceValue(strKeyNonunfInfo, m_strImageFilePath) )
        return false;

    if( 0 != oHeader.nWrite(sDtrekGetPrefix() + c_strInputHeaderName) )
        return false;

    return true;
}

// The return value represents the last successful processing step 
CDTMainRanker::enumProcessingStep CDTMainRanker::enDoImageProcessing()
{
    CDTMainRankerEnvironment*       poEnvironmentSetter = NULL;
    
    if( !m_bUseNonunf )   // disable nonunf temporarily
    {
        poEnvironmentSetter = new CDTMainRankerEnvironment();
    }

    if( "NO" == m_strFindOptions || "no" == m_strFindOptions || false == bDoFind() )
    {
        if( poEnvironmentSetter )
            delete poEnvironmentSetter;
        return enProcExtractHeader;
    }

    if( "NO" == m_strIndexOptions || "no" == m_strIndexOptions || false == bDoIndex() )
    {
        if( poEnvironmentSetter )
            delete poEnvironmentSetter;
        return enProcFind;
    }

    if( "NO" == m_strRefineOptions || "no" == m_strRefineOptions || false == bDoRefine() )
    {
        if( poEnvironmentSetter )
            delete poEnvironmentSetter;
        return enProcIndex;
    }

    if( poEnvironmentSetter )
        delete poEnvironmentSetter;

    return enProcRefine;
}

bool CDTMainRanker::bExecuteDTREKCommand(CDTMainRanker::enumProcessingStep enProcStep, Cstring& strOptions)
{
    CDTMain*        poDTMain = NULL;
    Cstring         strCommand("");

    switch( enProcStep )
    {
    case enProcExtractHeader:
        poDTMain = new CDTMainExtractHeader();
        strCommand += "dtextractheader";
        break;
    case enProcFind:
        poDTMain = new CDTMainFind();
        strCommand += "dfind";
        break;
    case enProcIndex:  
        poDTMain = new CDTMainIndex();
        strCommand += "dtindex";
        break;
    case enProcRefine: 
        poDTMain = new CDTMainRefine();
        strCommand += "dtrefine";
        break;
    case enProcRank: 
        poDTMain = new CDTMainRank();
        
        ((CDTMainRank*)poDTMain)->vSetSampleNameForReport(m_strSampleNameForReport);
        
        strCommand += "dtrank";
        break;
    default:
        return false;
    }
    
    Cstring         strLogFile = strCommand + ".log";
    
    strCommand += (" " + strOptions);

    vCommandLineToCommandArgumentArray(strCommand.string());

    //CCrclHelper::GetInstance()->RedirectStdout(strLogFile);
    int     nRet = poDTMain->nExecute(m_nDTREKCommandArgumentCount, m_ppcDTREKCommandArguments);
    
    delete  poDTMain;

    return nRet == 0 ? true : false;
}

bool CDTMainRanker::bDoFind()
{
    Cstring     strFindCommandLine = '"';
    strFindCommandLine += m_strCurrentHeaderFilePath;
    strFindCommandLine += "\" ";
    strFindCommandLine += m_strResolutionLimits;
    
    if( m_nTotalNumberOfResoBins > 0 ) // just verify
    {
        Cstring     strBinsNumber(m_nTotalNumberOfResoBins);
        strFindCommandLine += " -resobins ";
        strFindCommandLine += strBinsNumber;
        
        if( "" != m_strMinPeakCount )
        {
            strFindCommandLine += " ";
            strFindCommandLine += m_strMinPeakCount;
        }
    }
    
    strFindCommandLine += " ";

    ////////////////////////////////////////////////////////////////////////////////
    // Generate dtfind type sequence information from the vector of sequence ranges
    Cstring     strSeq("");
    int     nBeg = -1; 
    int     nEnd = -1;
    for(int ii=0; ii < m_vecScanImageSequenceRanges.size(); ii++)
    {
        nBeg = m_vecScanImageSequenceRanges[ii].first;
        nEnd = m_vecScanImageSequenceRanges[ii].second;
    
        strSeq += "-seq ";
        strSeq += nBeg;
        strSeq += " ";
        strSeq += nEnd;
        strSeq += " ";
    }
    
    strFindCommandLine += strSeq;
    //////////////////////////////////////////////////////////////////////////////////

    if( "" == m_strFindOptions )
    {
        ////////////////////////////////////////////////////////////////////////////////
        // Figure out -sigma to be used
        Cstring     strSigmaOption(" -sigma ");
        if( "" != m_strSigmaFIND )  // user-supllied
            strSigmaOption += m_strSigmaFIND;
        else
            strSigmaOption += "3 ";  // 3 -default I/sigma value
        
        strFindCommandLine += strSigmaOption;
        ///////////////////////////////////////////////////////////////////////////////

        strFindCommandLine += " -min 50 -2D -filter 6 -window 0 0 -rank";  // default dtfind options
    
    }
    else
    {
        /////////////////////////////////////////////////////////////////
        // Parse out the options that have to be set by dtranker
        vCommandLineToCommandArgumentArray(m_strFindOptions.string());
        Cstring     strTemp("");
        while( bGetCommandOption(c_pcSeq, strTemp) );
        while( bGetCommandOption(c_pcReso, strTemp) );
        while( bGetCommandOption(c_pcResoBins, strTemp) );
        while( bGetCommandOption(c_pcRefList, strTemp) );
        while( bGetCommandOption(c_pcOutputHdr, strTemp) );

        // Parse out -rank option, if the user did put it there. The problem is users often forget to do that, so we better put it ourselves.
        while( bGetCommandOption(c_pcRank, strTemp) );

        // Just find out if the user put -sigma on the find command line
        bool   bUserSuppliedSigmaOptionAsPartofFINDoption = bGetCommandOption(c_pcSigma, strTemp, false);
        if( !bUserSuppliedSigmaOptionAsPartofFINDoption && "" != m_strSigmaFIND )
            bUpdateCommandOption(c_pcSigma, m_strSigmaFIND);

        vGetCommandArrayString(m_strFindOptions);
        /////////////////////////////////////////////////////////////////

        strFindCommandLine += m_strFindOptions;

        strFindCommandLine += " -rank ";
    }

    strFindCommandLine += " -ref ";
    strFindCommandLine += (sDtrekGetPrefix() + c_strFindOutputReflnFileName);
    
    strFindCommandLine += " -out ";
    strFindCommandLine += (sDtrekGetPrefix() + c_strFindOutputHeaderName);
    
    if( !bExecuteDTREKCommand(enProcFind, strFindCommandLine) )
        return false;

    m_strCurrentHeaderFilePath = (sDtrekGetPrefix() + c_strFindOutputHeaderName);

    return true;
}

bool CDTMainRanker::bDoIndex()
{
    Cstring         strIndexCommandLine = m_strCurrentHeaderFilePath + " ";
    strIndexCommandLine += (sDtrekGetPrefix() + c_strFindOutputReflnFileName);
    strIndexCommandLine += " ";
    
    if( "" == m_strIndexOptions )
    {
        strIndexCommandLine += " -maxresid 3.0 -sigma 5 "; // default dtindex options
    }
    else
    {
        /////////////////////////////////////////////////////////////////
        // Parse out the options that have to be set by dtranker
        vCommandLineToCommandArgumentArray(m_strIndexOptions.string());
        Cstring     strTemp("");
        while( bGetCommandOption(c_pcOutputHdr, strTemp) );
        vGetCommandArrayString(m_strIndexOptions);
        /////////////////////////////////////////////////////////////////
        
        strIndexCommandLine += m_strIndexOptions;
    }

    strIndexCommandLine += " -out ";
    strIndexCommandLine += (sDtrekGetPrefix() + c_strIndexOutputHeaderName);
    
    if( !bExecuteDTREKCommand(enProcIndex, strIndexCommandLine) )
        return false;

    m_strCurrentHeaderFilePath = (sDtrekGetPrefix() + c_strIndexOutputHeaderName);

    return true;
}

bool CDTMainRanker::bDoRefine()
{
	const bool addQuotes = !strchr(m_strCurrentHeaderFilePath, '"') &&
		strchr(m_strCurrentHeaderFilePath, ' ');
	Cstring strRefineCommandLine;
	if(addQuotes) strRefineCommandLine += "\"";
	strRefineCommandLine += m_strCurrentHeaderFilePath;
	if(addQuotes) strRefineCommandLine += "\"";
	strRefineCommandLine += " ";

    //NOTE: since we want to refine on images, the sequence option must follow the header file name. This has to do with the way dtrefine parses the input.
    
    ///////////////////////////////////////////////////////////////////////////////////
    // Generate dtrefine type sequence information from the vector of sequence ranges
    Cstring     strSeq("");
    int     nBeg = -1; 
    int     nEnd = -1;
    for(int ii=0; ii < m_vecScanImageSequenceRanges.size(); ii++)
    {
        nBeg = m_vecScanImageSequenceRanges[ii].first;
        nEnd = m_vecScanImageSequenceRanges[ii].second;
    
        for(int jj=nBeg; jj <= nEnd; jj++)
        {
            strSeq += "-seq ";
            strSeq += jj;
            
            strSeq += " ";
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////
    
    strRefineCommandLine += strSeq;

    if( "" == m_strRefineOptions )
    {
        ////////////////////////////////////////////////////////////////////////////////
        // Figure out -sigma to be used
        Cstring     strSigmaOption(" -sigma ");
        if( "" != m_strSigmaREFINE )  // user-supllied
            strSigmaOption += m_strSigmaREFINE;
        else
            strSigmaOption += "5 ";  // 5 - default I/sigma value
        
        strRefineCommandLine += strSigmaOption;
        ///////////////////////////////////////////////////////////////////////////////

        strRefineCommandLine += " -rank +All -rej 1 1 1 -cycles 30 -verbose 0 -go -go -go -verbose 1 -go";  // default dtrefine options
    }
    else
    {
        /////////////////////////////////////////////////////////////////
        // Parse out the options that have to be set by dtranker
        vCommandLineToCommandArgumentArray(m_strRefineOptions.string());
        Cstring     strTemp("");
        while( bGetCommandOption(c_pcSeq, strTemp) );
        while( bGetCommandOption(c_pcReso, strTemp) );
        while( bGetCommandOption(c_pcResoBins, strTemp) );
        while( bGetCommandOption(c_pcOutputHdr, strTemp) );

        // Parse out -rank option, if the user did put it there. The problem is users often forget to do that, so we better put it ourselves.
        while( bGetCommandOption(c_pcRank, strTemp) );

        // Just find out if the user put -sigma on the refine command line
        bool   bUserSuppliedSigmaOptionAsPartofREFINEoption = bGetCommandOption(c_pcSigma, strTemp, false);

        strRefineCommandLine += " -rank ";

        if( m_strETE != "" )
        {
            strRefineCommandLine += " -rankete ";
            strRefineCommandLine += m_strETE;
            strRefineCommandLine += " ";
        }

        if( !bUserSuppliedSigmaOptionAsPartofREFINEoption && "" != m_strSigmaREFINE )
            bUpdateCommandOption(c_pcSigma, m_strSigmaREFINE);
        
        vGetCommandArrayString(m_strRefineOptions);
        strRefineCommandLine += m_strRefineOptions;
    }

    strRefineCommandLine += " ";
    strRefineCommandLine += m_strResolutionLimits;
    
    if( m_nTotalNumberOfResoBins > 0 ) // just verify
    {
        Cstring     strBinsNumber(m_nTotalNumberOfResoBins);
        strRefineCommandLine += " -resobins ";
        strRefineCommandLine += strBinsNumber;
        
        if( "" != m_strMinPeakCount )
        {
            strRefineCommandLine += " ";
            strRefineCommandLine += m_strMinPeakCount;
        }
    }

    strRefineCommandLine += " -out ";
    strRefineCommandLine += (sDtrekGetPrefix() + c_strRefineOutputHeaderName);

    if( !bExecuteDTREKCommand(enProcRefine, strRefineCommandLine) )
        return false;

    m_strCurrentHeaderFilePath = (sDtrekGetPrefix() + c_strRefineOutputHeaderName);

    return true;
}

bool CDTMainRanker::bDoRank()
{
    Cstring     strRankCommandLine("");
    if( "" != m_strCurrentHeaderFilePath )
        strRankCommandLine = m_strCurrentHeaderFilePath + " ";
    
    bool    bExcludeSearchRules = false;
    bool    bExcludeRefinementRules = false;

    if( "" != m_strExcludeRankRules )
    {
        std::vector<Cstring>    saExcludeRules;
        m_strExcludeRankRules.nListToVector(saExcludeRules, " ");

        for(int ii=0; ii < saExcludeRules.size(); ii++)
        {
            // Effectively disable Rules by setting the award/penalty to zero
            if( saExcludeRules[ii] == "peaksearch" )
            {
                strRankCommandLine += "-param 1 5 0 -param 1 6 0 -param 2 5 0 -param 2 6 0 ";     
                bExcludeSearchRules = true;
            }
            else if(  saExcludeRules[ii] == "icerings" )
	      {
                strRankCommandLine += "-param 4 2 0 -param 5 2 0 ";     
		cout << "RULES EXCLUDE icerings in effect\n" << flush;
	      }
            else if( saExcludeRules[ii] == "refine" ) // indexing/refinement
            {
                bExcludeRefinementRules = true;
                
                //strRankCommandLine += "-param 3 1 0 -param 3 2 0 -param 3 4 0 -param 3 5 0";  // peak shape                 
                strRankCommandLine += "-param 6 2 0 -param 6 4 0 ";      // indexing            
        
                strRankCommandLine += "-param 7 2 0 -param 8 2 0 ";   // mosaicity and r.m.s.d  
                
                strRankCommandLine += "-param 9 2 0 -param 9 4 0 ";    // refined spot count                
                
                strRankCommandLine += "-param 10 5 0 -param 10 6 0 -param 11 5 0 -param 11 6 0 ";   // predict/found  
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Modify Rules 1,2,10 and 11 to communicate the peak count and number of resolution bins
    Cstring     str1(""), str2(""), str10(""), str11(""), strNumber("");
    
    if( "" != m_strMinPeakCount )
    {
        str1  = "-param 1 4 "  + m_strMinPeakCount;
        Cstring     str2  = "-param 2 4 "  + m_strMinPeakCount;
        Cstring     str10 = "-param 10 4 " + m_strMinPeakCount;
        Cstring     str11 = "-param 11 4 " + m_strMinPeakCount;
        
        strRankCommandLine += " ";
        strRankCommandLine += str1;

        strRankCommandLine += " ";
        strRankCommandLine += str2;

        strRankCommandLine += " ";
        strRankCommandLine += str10;

        strRankCommandLine += " ";
        strRankCommandLine += str11;
    }

    ///////////////////////////////////////////////////////
    if( m_nTotalNumberOfResoBins > 0 ) 
    {
        strNumber = Cstring(m_nTotalNumberOfResoBins);
        
        str1  = "-param 1 1 "  + strNumber;
        str2  = "-param 2 1 "  + strNumber;
        str10 = "-param 10 1 " + strNumber;
        str11 = "-param 11 1 " + strNumber;
        
        strRankCommandLine += " ";
        strRankCommandLine += str1;

        strRankCommandLine += " ";
        strRankCommandLine += str2;

        strRankCommandLine += " ";
        strRankCommandLine += str10;

        strRankCommandLine += " ";
        strRankCommandLine += str11;
    }
    
    /////////////////////////////////////////////////////////
    if( m_nNumberOfHighResoBins >= 0 ) 
    {
        strNumber = Cstring(m_nNumberOfHighResoBins);
        
        str1  = "-param 1 2 "  + strNumber;
        str2  = "-param 2 2 "  + strNumber;
        str10 = "-param 10 2 " + strNumber;
        str11 = "-param 11 2 " + strNumber;
    
        strRankCommandLine += " ";
        strRankCommandLine += str1;

        strRankCommandLine += " ";
        strRankCommandLine += str2;

        strRankCommandLine += " ";
        strRankCommandLine += str10;

        strRankCommandLine += " ";
        strRankCommandLine += str11;
    }
    
    //////////////////////////////////////////////////////////
    if( m_nNumberOfLowResoBins >= 0 ) 
    {
        strNumber = Cstring(m_nNumberOfLowResoBins);
        
        str1  = "-param 1 3 "  + strNumber;
        str2  = "-param 2 3 "  + strNumber;
        str10 = "-param 10 3 " + strNumber;
        str11 = "-param 11 3 " + strNumber;
    
        strRankCommandLine += " ";
        strRankCommandLine += str1;

        strRankCommandLine += " ";
        strRankCommandLine += str2;

        strRankCommandLine += " ";
        strRankCommandLine += str10;

        strRankCommandLine += " ";
        strRankCommandLine += str11;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////
    // Modify Rules 2 and 11, if the user wants to change the rotation anisotropy penalty  
    // coefficient.
    if( "" == m_strAnisotropy )
    {
        m_strAnisotropy = nGetNumberOfImagesInSeqRanges() > 1 ? 20 : 0;   // 20 points is the default value
    }

    if( !bExcludeSearchRules )
    {
        str1 = " -param 2 6 ";
        str1 += m_strAnisotropy;
        str1 += " ";
    }

    if( !bExcludeRefinementRules )
    {
        str1 += " -param 11 6 ";
        str1 += m_strAnisotropy;
        str1 += " ";
    }

    strRankCommandLine += " ";
    strRankCommandLine += str1;
    ///////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    // Modify Rule 12, if the user wants to change the strategy penalty coefficient.
    if( "" == m_strStrategy )
    {
        m_strStrategy = 0;   // 0 points is the default value
    }

    str1 = "-param 12 1 ";
    str1 += m_strStrategy;

    strRankCommandLine += " ";
    strRankCommandLine += str1;
    ///////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    // Finally add whatever the user set explicitly as -param options
    if( "" != m_strRankOptions )
    {
        strRankCommandLine += " ";
        strRankCommandLine += m_strRankOptions;
    }

    if( !bExecuteDTREKCommand(enProcRank, strRankCommandLine) )
        return false;
    
    m_strCurrentHeaderFilePath = (sDtrekGetPrefix() + c_strRankOutputHeaderName);

    return true;
}

/////////////////////////////////////////////////////////////////////////////
bool CDTMainRanker::bGetImageTemplate()
{
    // First check if we have a valid header
    Cimage_header        oHeader(m_strCurrentHeaderFilePath);
    if( !oHeader.bIsAvailable() )
        return false;
    
    if ( 0 != oHeader.nGetValue(D_K_ScanTemplate, &m_strTemplate) )
        return false;

    return true;
}
//////////////////////////////////////////////////////////////////////////////
bool CDTMainRanker::bDeriveImageSequence()
{
    // If the user has supplied an image, but did not supply sequence information,
    // figure out the image number from the template in the extracted header and
    // the image file name. 
    int     nImageNumber = -1;
    
    if( "" != m_strImageFilePath && "" != m_strTemplate )
    {
        Cstring         strImageNumber("");
        bool            bFoundQuestionMarkInTemplate = false;
        Cstring         strTemplateFileName = m_strTemplate.sGetFileName();
        Cstring         strImageFileName = m_strImageFilePath.sGetFileName();
        
        for(int ii=0; ii < strTemplateFileName.length(); ii++)
        {
            if( '?' == strTemplateFileName.GetAt(ii) )
            {
                bFoundQuestionMarkInTemplate = true;
                strImageNumber += strImageFileName.GetAt(ii);
            }
            else if( bFoundQuestionMarkInTemplate )
                break;  // if the question marks have ended, get out
        }

        if( "" == strImageNumber )
            return false; // no question marks in scan template found

        nImageNumber = atoi(strImageNumber.string());
        
        if( nImageNumber < 0 )      // A sanity check
            return false;
    
        m_vecScanImageSequenceRanges.push_back(std::pair<int,int>(nImageNumber,nImageNumber));
    }
    // If the image was not supplied, try and get the sequence information from the find or refine options.
    // If find and refine have different sequence options, give preference to find.
    // RB 8/15/05  Actually, now we DON'T want the users to pass -seq on the dtfind/dtrefine option,
    // so the two cases below are obsolete.
    else if( "" != m_strFindOptions && ("NO" != m_strFindOptions && "no" != m_strFindOptions) )
    {
        vGetImageSequenceInformationFromDTREKCommandLine(enProcFind);
    }
    else if( "" != m_strRefineOptions && ("NO" != m_strRefineOptions && "no" != m_strRefineOptions) )
    {
        vGetImageSequenceInformationFromDTREKCommandLine(enProcRefine);
    }
    //// If header was supplied, try and get the scan info from the header itself, assuming that
    //// that seq_info means: first image, increment, last image,  
    //   TODO: go from BEGIN STEP END (in header) to BEGIN-END, BEGIN-END, BEGIN-END
    //else if( "" != m_strInputHeaderFilePath )
    //{
    //    int         a3nSeqInfo[3] = {0};
    //    if ( 0 == oHeader.nGetValue(D_K_ScanSeqInfo, 3, a3nSeqInfo) )
    //    {
    //        nImageNumber = a3nSeqInfo[0];
    //        
    //        if( nImageNumber < 0 )      // A sanity check
    //            return false;
    //
    //        m_strImageScanSequence = "-seq ";
    //        m_strImageScanSequence += Cstring(nImageNumber);
    //        m_strImageScanSequence += " ";
    //        m_strImageScanSequence += Cstring(nImageNumber);
    //    }
    //}
    
    return true;
}

void CDTMainRanker::vGetImageSequenceInformationFromDTREKCommandLine(enumProcessingStep eStep)
{
    Cstring         strCommand = ( eStep == enProcFind ) ? m_strFindOptions : m_strRefineOptions;
    
    Cstring         strTemp1("");
    Cstring         strTemp2("");

    vCommandLineToCommandArgumentArray(strCommand.string());
    
    m_vecScanImageSequenceRanges.clear();

    /////////////////////////////////////////////////////////////////////////////////////////
    // Parse the sequence information and fill up the sequence array with "begin-end" pairs
    Cstring                 strSeqArgs("");
    std::vector<int>   vecSeqs;    
    while( bGetCommandOption(c_pcSeq, strSeqArgs) )
    {
        strSeqArgs.nListToVector(vecSeqs, " ");  
    
        if( 0 == vecSeqs.size() )
            continue; // useless seq option
        else if( 1 == vecSeqs.size() )
            m_vecScanImageSequenceRanges.push_back(std::pair<int,int>(vecSeqs[0],vecSeqs[0]));
        else if( 2 == vecSeqs.size() )
            m_vecScanImageSequenceRanges.push_back(std::pair<int,int>(vecSeqs[0],vecSeqs[1]));
    }
}
///////////////////////////////////////////////////////////////////
bool CDTMainRanker::bUpdateHeaderFileAccordingToRankerCommandLine()
{
    if( "" == m_strBeamCenter        && 
        "" == m_strDetectorDistance  && 
        "" == m_strTwoTheta          && 
        "" == m_strWavelength )
        return true; // nothing to update
    
    // First check if we have a valid header
    Cimage_header        oHeader(m_strCurrentHeaderFilePath);
    if( !oHeader.bIsAvailable() )
        return false;
    
    Cstring              strDetectorPrefix(""); 
    if ( 0 != oHeader.nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    if( "" != m_strBeamCenter )
    {
        float       fBeamPos0 = 0.0, fBeamPos1 = 0.0;
        if( 0 == sscanf(m_strBeamCenter.string(), "%f %f", &fBeamPos0, &fBeamPos1) )
            return false;
        
        Cspatial             oSpatial(oHeader, strDetectorPrefix);
        if( !oSpatial.bIsAvailable() )
            return false;

        oSpatial.nSetBeamPosition(fBeamPos0, fBeamPos1);

        oSpatial.nUpdateHeader(&oHeader, strDetectorPrefix);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( "" != m_strDetectorDistance || "" != m_strTwoTheta )
    {
        Cgoniometer         oDetGonio(oHeader, strDetectorPrefix);
        if( !oDetGonio.bIsAvailable() )
            return false;
    
        if( "" != m_strDetectorDistance )
        {
            float   fDistance = 0.0;
            if( 0 == sscanf(m_strDetectorDistance.string(), "%f", &fDistance) )
                return false;
            
            oDetGonio.nSetDistance(fDistance);
        }
    
        if( "" != m_strTwoTheta )
        {
            float   fTwoTheta = 0.0;
            if( 0 == sscanf(m_strTwoTheta.string(), "%f", &fTwoTheta) )
                return false;
            
            oDetGonio.nSetSwing(fTwoTheta);
        }
        
        oDetGonio.nUpdateHeader(&oHeader, strDetectorPrefix);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( "" != m_strWavelength )
    {
        Cstring     strSourcePrefix("SOURCE_");
        float       fWavelength = 0.0;
        if( 0 == sscanf(m_strWavelength.string(), "%f", &fWavelength) )
            return false;
        
        Cwavelength             oWavelength(oHeader, strSourcePrefix, false);

        oWavelength.nSetWavelength(fWavelength);

        oWavelength.nUpdateHeader(&oHeader, strSourcePrefix);
    }

    oHeader.nWrite(m_strCurrentHeaderFilePath);

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainRanker::bOutputResults(int nNumberOfRankingRules)
{
    if( "" == m_strResultsFileName || m_nNumberOfFilePathChars < l_nMinNumberOfFilePathChars ) // sanity check
        return false;

    /////////////////////////////////////////////////////
    // Check if the file exists already
    bool            bPreExistentFile = bFileExists(m_strResultsFileName);
    /////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // Open file in append mode
    FILE*       pfOut = fopen(m_strResultsFileName,"at");
    if( !pfOut )
        return false;
    /////////////////////////////////////////////////////
    
    char        cRuleTxt[256];
    Cstring     strLine("");
    
    Cstring     strFormat = Cstring("%") + Cstring(m_nNumberOfFilePathChars);
    strFormat += "s%s";
    
    if( !bPreExistentFile )
    {
        fprintf(pfOut, "dtranker:  Copyright (c) 2006 Rigaku\n");
        fprintf(pfOut, "%s\n\n", D_K_DTREKVersion);

        fprintf(pfOut, "%s\n", s_strRankingRules);

        Cstring     strTableHeader("Sample / Rules");
        
        /////////////////////////////////////////////////////
        // Construct the header
        for(int nRule=1; nRule <= nNumberOfRankingRules; nRule++)
        {
            if( nRule != 3 ) // temporary hack
            {
                sprintf(cRuleTxt, "%5d", nRule);
                strLine += Cstring(cRuleTxt);
            }
            else
            {
                Cstring     strRule_3_str("    3   3a   3b   3c   3d");
                strLine  += strRule_3_str;
            }
        }
        strLine += "  Total";

        fprintf(pfOut, strFormat.string(), strTableHeader.string(), strLine.string());    
        fprintf(pfOut, "\n\n");
    }

    ////////////////////////////////////////////////////////////
    // Check if we have a valid header

#ifdef WIN32
    _flushall();     // otherwise the log file gets messed up, at least in debug
#endif
    
    Cimage_header        oHeader(m_strCurrentHeaderFilePath);
    if( !oHeader.bIsAvailable() )
        return false;

    int         nVal = 0;
    strLine = "";
    for(int nRule=1; nRule <= nNumberOfRankingRules; nRule++)
    {
        if( nRule != 3 ) // temporary hack
        {
            if( 0 != oHeader.nGetValueDictionary("DTRANK_RANK", "RULE_*", nVal, nRule) )
                nVal = 0;
        
            sprintf(cRuleTxt, "%5d", nVal); 
            strLine += Cstring(cRuleTxt);
        }
        else
        {
            int     anTemp[5] = {0};
            const   int cnNumberOfRule_3_Scores = 5;
            Cstring     sTemp("DTRANK_RANK_RULE_3");
            oHeader.nGetValue(sTemp, cnNumberOfRule_3_Scores, anTemp); // 1 total score for the rule + 4 subscores
            sprintf(cRuleTxt, "%5d%5d%5d%5d%5d", anTemp[0],anTemp[1],anTemp[2],anTemp[3],anTemp[4]); 
            strLine += Cstring(cRuleTxt);
        }
    }
    
    if( 0 != oHeader.nGetValue("DTRANK_RANK_CUMULATIVE", &nVal) )
        nVal = 0;
    
    sprintf(cRuleTxt, "%7d", nVal); // => cRuleTxt array size must be at least 7+1
    strLine += Cstring(cRuleTxt);

    fprintf(pfOut, strFormat.string(), m_strSampleNameForReport.string(), strLine.string());
    fprintf(pfOut, "\n");

    fclose(pfOut);

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainRanker::vSetResoBinsParams()
{
    if( "" == m_strResolutionBins )
        return;

    std::vector<int>    vecBinNumbers;
    m_strResolutionBins.nListToVector(vecBinNumbers, " ");

    int     nSize = vecBinNumbers.size();

    if( nSize > 0 )
        m_nTotalNumberOfResoBins = vecBinNumbers[0];
    
    if( nSize > 1 )
        m_nNumberOfHighResoBins  = vecBinNumbers[1];

    if( nSize > 2 )
        m_nNumberOfLowResoBins   = vecBinNumbers[2];
}
/////////////////////////////////////////////////////////////////////////////////////////////
int CDTMainRanker::nGetNumberOfImagesInSeqRanges()
{
    int     nBeg = -1;
    int     nEnd = -1;
    int     nTot = 0;
    
    for(int ii=0; ii < m_vecScanImageSequenceRanges.size(); ii++)
    {
        nBeg = m_vecScanImageSequenceRanges[ii].first;
        nEnd = m_vecScanImageSequenceRanges[ii].second;
    
        nTot += (nEnd - nBeg + 1);
    }
    
    return nTot;
}
//////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainRanker::bGenerateSampleNameForReport()
{
    if( 0 >= m_vecScanImageSequenceRanges.size() )
    {
        printf("ERROR: cannot determine the image sequence. Please provide the \"sequence\" option to dtranker.\n");
        return false;
    }

    if( "" != m_strImageFilePath && 
        1 == m_vecScanImageSequenceRanges.size() && 
        m_vecScanImageSequenceRanges[0].second == m_vecScanImageSequenceRanges[0].first
      )  // one image only
    {
        if( "" != m_strTemplate )  
        {
            Cscan       oScan(m_strTemplate, m_vecScanImageSequenceRanges[0].first, 1);    
            oScan.nGetImageName(&m_strSampleNameForReport);
        }            
        else    
            m_strSampleNameForReport = m_strImageFilePath;
    }
    else if( "" != m_strTemplate )
    {
        m_strSampleNameForReport = m_strTemplate + " ";
        
        m_strSampleNameForReport += m_vecScanImageSequenceRanges[0].first;
        
        m_strSampleNameForReport += '+';
    }    
    else
        m_strSampleNameForReport = m_strInputHeaderFilePath; 
    
    m_strSampleNameForReport.vAbbreviateFilePath(m_nNumberOfFilePathChars);
    
    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainRanker::bParseResultsFileOption()
{
    // If no user input, use default file name and default number of chars in file path.
    if( "" == m_strResultsFileName )
    {
        m_strResultsFileName = c_pcResultsFile;
        m_nNumberOfFilePathChars = l_nMinNumberOfFilePathChars;
    }
    else   
    {   
        // If one arg - it has to be the number of chars
        // If two args - the first is filename, and the second - the number of chars
        std::vector<Cstring>        asParams;
        m_strResultsFileName.nListToVector(asParams, " ");
        switch( asParams.size() )
        {
        case 1:
            if( asParams[0].bIsNumeric() )
            {
                int     nTemp = atoi(asParams[0].string());
                m_nNumberOfFilePathChars = nTemp > l_nMinNumberOfFilePathChars ? nTemp : l_nMinNumberOfFilePathChars; 
                
                m_strResultsFileName = c_pcResultsFile;
            }
            else
                m_strResultsFileName = asParams[0];
            break;
        case 2:
            m_strResultsFileName = asParams[0];
            if( asParams[1].bIsNumeric() )
            {
                int     nTemp = atoi(asParams[1].string());
                m_nNumberOfFilePathChars = nTemp > l_nMinNumberOfFilePathChars ? nTemp : l_nMinNumberOfFilePathChars; 
            }
            break;
        default:
            m_strResultsFileName = c_pcResultsFile;
            m_nNumberOfFilePathChars = l_nMinNumberOfFilePathChars;
        }
    }

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainRanker::vError(const int nErrorNum, const Cstring& sMessage)
{
    if( 0 != nErrorNum )
    {
        cout << sMessage << endl;
    }
    
    cout << "\ndtranker - Usage:\n\n"
       "dtranker image=\"sImageFile\" [header=\"sHeaderFile\"] [sequence=\"aSeqList\"]\n"
       "[reso=\"sResoParams\"] [resobins=\"[nTotalBins [nHighBins [nLowBins]]]\"]\n"
       "[beampos=\"nPix0 nPix1\"] [distance=fDist] [2Theta=f2Theta]\n"
       "[wavelength=fWavelength] [find=\"sFindOptions\"] [index=\"sIndexOptions\"]\n"
       "[refine=\"sRefineOptions\"] [rules=\"sRuleParams\"]\n"
       "[rules_exclude=\"sRulesExclude\"] [minpeakcount=nMinPeakCount]\n" 
       "[sigmafind=fSigmaFIND] [sigmarefine=fSigmaREFINE]\n" 
       "[aniso=nAniso] [strategy=nStrategy]\n"
       //"[aniso=nAniso] [strategy=nStrategy] [exptimeeval=fSigmaETE [fTimePower]]\n"
       "[results=\"[sResultsFileName] | [nImagePathChars]\"]\n"
       "\n"
       "Command line options:\n"
       "\n"
       "image=\"sImageFile\"\n"         
       "                   Image file path. If this is not specified, a header file\n"
       "                   path must be specified.\n"
       "\n"
       "header=\"sHeaderFile\"\n"
       "                   Header file path. If this is not specified, an image file\n"
       "                   path must be specified. If both header and image files are\n" 
       "                   specified, the image file specification will be ignored.\n"  
       "                   NOTE: if header is specified then image sequence must be\n"
       "                   specified as well.\n"
       "\n"
       "sequence=\"aSeqList\"\n"
       "                   List of image numbers to be substituted to the image path\n"
       "                   template to be passed to dtfind and dtrefine. The list may\n"
       "                   have one or more group(s) of images. Each group of images\n"
       "                   may contain one or more images in it.\n"
       "                   EXAMPLE: sequence=\"1-2 5 10-11\"\n" 
       "\n"
       "reso=\"sResoParams\"\n"
       "                   Resolution parameters. Either two numbers for the minimum\n"
       "                   and maximum resolution, or a keyword \"edge\" or \"corner\".\n"
       "                   If \"edge\" is specified, then the maximum resolution will be\n"
       "                   set to the highest *minimum* resolution found on any of the\n"
       "                   four edges (sides) of the image. If \"corner\" is specified,\n"
       "                   then the maximum resolution will be set to the highest\n"
       "                   *maximum* resolution found on any of the four edges (sides)\n"
       "                   of the image.\n"
       "                   DEFAULT: edge\n"
       "\n"
       "resobins=\"nTotalBins [nHighBins [nLowBins]]\"\n"
       "       nTotalBins  Number of resolution bins used when processing image.\n"
       "                   DEFAULT: 10 bins.\n"
       "       nHighBins   Number of high resolution bins used in ranking.\n"          
       "                   DEFAULT: 7 bins.\n"
       "       nLowBins    Number of low resolution bins used in ranking.\n"          
       "                   DEFAULT: 0 bins.\n"
       "                   NOTE: Increasing the total number of resolution bins\n"
       "                   will reduce the number of reflections found in a bin.\n"
       "                   As a result, Rules 1,2,10,11 may be unable to award\n"
       "                   any points for reflections found on the image.\n"
       "                   Please also see the description of \"minpeakcount\".\n"
       "\n"
       "beampos=\"nPix0 nPix1\"\n"
       "                   Direct beam position, pixels, d*TREK convention.\n"
       "                   DEFAULT: beam position from the image header file.\n"
       "\n"
       "distance=fDist     Detector distance, mm\n"
       "                   DEFAULT: distance from the image header file.\n"
       "\n"
       "2Theta=f2Theta     Detector swing angle, deg\n"
       "                   DEFAULT: swing angle from the image header file.\n"
       "\n"
       "wavelength=fWavelength\n"        
       "                   X-ray source wavelength, A\n"
       "                   DEFAULT: wavelength from the image header file.\n"
       "\n"
       "find=\"sFindOptions\"\n"
       "                   dtfind command-line options. Type \"dtfind -h\" for help.\n"
       "                   NOTE: do not specify input/output header files.\n"
       "                   NOTE: do not specify output reflection file.\n"
       "                   NOTE: do not specify -seq option.\n"
       "                   NOTE: do not specify -reso or -resobins option.\n"
       "                   NOTE: to skip peak search, indexing and refinement pass\n"
       "                   sFindOptions as NO.\n"
       "                   DEFAULT: \"-sigma 3 -min 50 -2D -filter 6 -window 0 0\n"
       "                   -rank [-reso sResoParams] [-resobins nTotalBins]\"\n"
       "\n"
       "index=\"sIndexOptions\"\n"
       "                   dtindex command-line options. Type \"dtindex -h\" for help.\n"
       "                   NOTE: do not specify input/output header files.\n"
       "                   NOTE: do not specify input reflection file.\n"
       "                   NOTE: to skip indexing and refinement, pass\n"
       "                   sIndexOptions as NO.\n"
       "                   DEFAULT: \"-maxresid 3.0 -sigma 5\"\n"
       "\n"
       "refine=\"sRefineOptions\"\n"
       "                   dtrefine command-line options. Type \"dtrefine -h\" for help.\n"
       "                   NOTE: do not specify input/output header files.\n"
       "                   NOTE: do not specify -seq option.\n"
       "                   NOTE: do not specify -reso or -resobins option.\n"
       "                   NOTE: to skip refinement, pass sRefineOptions as NO.\n"
       "                   DEFAULT: \"-rank +All -sigma 5.0 -rej 1 1 1 \n"
       "                   -cycles 30 -verbose 0 -go -go -go -verbose 1 -go [-reso sResoParams]\n"
       "                   [-resobins nTotalBins]\"\n" 
       "\n"
       "rules=\"sRuleParams\"\n"
       "                   Ranking rules parameters. Type \"dtranker -help rules\"\n"         
       "                   for help on the rules with default parameters. To see\n"
       "                   the rules with your custom parameters include string \"-help\"\n"
       "                   in sRuleParams, for example:\n"
       "                   dtranker aniso=80 rules=\"-param 2 1 15 -param 9 1 11 -help\"\n"
       "\n"
       "rules_exclude=\"sRulesExclude\"\n"
       "                   A list of rule groups to be excluded, separated by spaces.\n"         
       "                   Rule groups are:\n"
       "                   \"peaksearch\" - exclude rules, related to peak search.\n"
       "                   \"icerings\" - exclude rules, related to ice rings.\n"
       "                   \"refine\" - exclude rules, related to indexing and\n"
       "                   refinement.\n"
       "                   NOTE: excluding rules should not be confused with skipping\n"
       "                   image processing steps like \"find=NO\".\n"
       "                   DEFAULT: Exclude no rules.                                   \n"
       "\n"
       "minpeakcount=nMinPeakCount\n"
       "                   Minimum peak count per resolution shell. This will be passed\n"                                                                        
       "                   as a parameter to ranking rules 1-2 and 10-11.\n"
       "\n"
       "sigmafind=fSigmaFIND\n"
       "                   Minimum ratio Intensity/Sigma for a reflection to be\n"
       "                   selected by the peak search (dtfind).\n"
       "                   NOTE: this option is equivalent to '-sigma' in option 'find'.\n"
       "\n"
       "sigmarefine=fSigmaREFINE\n"
       "                   Minimum ratio Intensity/Sigma for a reflection to be\n"
       "                   used in the refinement (dtrefine).\n"
       "                   NOTE: this option is equivalent to '-sigma' in option\n" 
       "                   'refine'.\n"
       "\n"
       "aniso=nAniso       A non-negative integer coefficient, used for calculating\n"
       "                   a penalty for rotational anisotropy. This will be passed as\n"
       "                   a parameter to ranking rules 2 and 11.\n"
       "                   NOTE: If nAniso is 0, no penalty will be applied.\n"
       "                   DEFAULT: 20.\n"
       "\n"
       "strategy=nStrategy\n" 
       "                   A non-negative integer coefficient, used for calculating\n"
       "                   a penalty for rotation width from strategy. This will be\n"
       "                   passed as a parameter to ranking rule 12.\n"
       "                   NOTE: If nStrategy is 0, no penalty will be applied.\n"
       "                   DEFAULT: 0.\n"
       "\n"
       //"exptimeeval=fSigmaETE [fTimePower]\n" 
       //"                   Presence of this option indicates dtranker will attempt\n"
       //"                   to evaluate the minimum exposure time required to see\n"
       //"                   reflections in the highest resolution shell. fSigmaETE is\n"
       //"                   Intensity/Sigma ratio required for the highest resolution shell.\n"
       //"                   The minimum recommended value for fSigmaETE is 10.
       //"                   fTimePower = is a power parameter of function I/Sigma vs time.
       //"                   This is a detector-dependent parameter and needs to be
       //"                   determined experimentally for each detector.
       //"                   DEFAULT: no exposure time evaluation will be attempted.\n"
       //"\n"
       "results=\"[sResultsFileName] | [nImagePathChars]\"\n"
       "                   sResultsFileName is the name of an ASCII file with summary\n"                                                                        
       "                   of the ranking results. nImagePathChars is the number\n"
       "                   of characters for the image template in that file.\n" 
       "                   NOTE: If the file already exists, it will be appended.\n"
       "                   DEFAULTS: sResultsFileName = RankSummary.txt.\n"
       "                   nImagePathChars = 24.\n"
       "\n"
       "\n"
       "Examples:\n"
       "\n"
       "   dtranker image=\"C:\\Projects\\Rank\\MySample101__screen0001.img\"\n"
       "\n"
       "   dtranker image=\"C:\\Projects\\Rank\\MySample101__screen0001.img\"\\\n"
       "            sequence=\"1 90\"\n"
       "\n"
       "   dtranker image=\"C:\\Projects\\Rank\\MySample101__screen0001.img\"\\\n"
       "            beampos=\"514 515\" distance=40 sigmafind=5 results=50\n"
       "\n"
       "   dtranker header=\"C:\\Projects\\Rank\\MySample101.head\"\n"
       "            sequence=\"15 135\"\\\n"
       "            find=\"-sigma 3 -min 50 -2D -filter 6 -window 0 0\"\\\n"
       "            index=\"-maxresid 3.0 -sigma 5\"\\\n"
       "            refine=NO\\\n"
       "            results=\"Results.txt\"\n"
       "\n"
       "   dtranker header=\"C:\\Projects\\Rank\\MySample101.head\"\\\n"
       "            sequence=\"1 90\" rules_exclude=\"icerings refine\"\\\n"
       "            results=\"Results.txt 55\"\n"
       "\n\n"
       << flush;
}
///////////////////////////////////////////////////////////////////////////////////////////////
