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
// DTMainExposureTimeEvaluator.cpp       Initial author: RB          06-Jan-2006
// This file contains the member functions of class CDTMainExposureTimeEvaluator

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
 
#include "DTMainExposureTimeEvaluator.h"
#include "ImageExpTimeEvaluator.h"
#include "Cscan.h"


/////////////////////////////////////////////////////////////////////////////////////////////////
const char*     c_pcSaturatedPercentage     = "-sat";  
const char*     c_pcTargetReso              = "-reso";  
const char*     c_pcStartTests              = "-start";  
const char*     c_pcTargetIOverSigma        = "-sigma";  
const char*     c_pcImageSequence           = "-seq";  

/////////////////////////////////////////////////////////////////////////////////////////////////
CDTMainExposureTimeEvaluator::CDTMainExposureTimeEvaluator()
{
    m_poEvaluator = new CImageExposureTimeEvaluator();

    m_sHeaderOut = sDtrekGetPrefix() + "dtete.head";
    m_sReflnListOut = "";
}
/////////////////////////////////////////////////////////////////////////////////////////////////
CDTMainExposureTimeEvaluator::~CDTMainExposureTimeEvaluator()
{
    if( m_poEvaluator )
        delete m_poEvaluator;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
int CDTMainExposureTimeEvaluator::nExecute(unsigned int argc, char* argv[])
{
    vDtrekSetModuleName("dtete");
    std::cout << "\ndtete:  Copyright (c) 2006 Rigaku\n";
    std::cout << D_K_DTREKVersion << std::endl;


#ifdef SSI_PC
    Cstring sCCAppVersion = (const char*) argv[argc-1];
	if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1)
    {
		std::cout << sCCAppVersion;
		argc -= 2; 
		std::cout << CCrclHelper::GetInstance()->GetTimestamp() << std::endl;
	}
#endif

    // Copy command line to output log
    std::cout << "\nCommand line:\n" << sGetCommandLine(argc, argv, 71) << std::endl << std::endl << std::flush;

    /////////////////////////////////////////////////////////////////////////////////////////////
    
    vSetCommandArguments(argc, argv);

    if( bIsHelpRequest() )
    {
        DTREK_ERROR(0, "");  // display general help
    }

    if( argc < 2 )
    {
        DTREK_ERROR(1, "ERROR - no header filename!\n");
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
    
    m_poHeader->nClearDictionary("DTETE", "RECOMMENDED_EXP_TIME");
    
    if( bGetCommandOption(c_pcStartTests) )
    {
        m_poHeader->nReplaceValueDictionary("DTETE", "MAXIMUM_REFLN_COUNT_IN_TARGET_SHELL", 0);
        m_poEvaluator->vSetMaxObservedPeakCountInTargetShell(0); 
    }
    else
    {
        int     nMaxReflnCount = 0;
        if( 0 == m_poHeader->nGetValueDictionary("DTETE", "MAXIMUM_REFLN_COUNT_IN_TARGET_SHELL", nMaxReflnCount) )
            m_poEvaluator->vSetMaxObservedPeakCountInTargetShell(nMaxReflnCount); 
    }

    //////////////////////////////////////////////////////////////////////////////////
    vSetOutputHeaderPath();
    vSetOutputReflnListPath();
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    // Get input targets
    char        acTemp[DTREK_MAX_LINE_LENGTH];

    double    dTargetResolution = -1.0;
    bGetCommandOption(c_pcTargetReso, dTargetResolution);

    double    dTargetIOverSigma = -1.0;
    bGetCommandOption(c_pcTargetIOverSigma, dTargetIOverSigma);
    
    double    dMaxSaturatedPeak_Percent = 0.0;
    if( !bGetCommandOption(c_pcSaturatedPercentage, dMaxSaturatedPeak_Percent) )
    {
        sprintf(acTemp, "Option %s must be specified!", c_pcSaturatedPercentage);
        sError = acTemp;

        DTREK_ERROR(1, sError);
    }
    ///////////////////////////////////////////////////////////////////////////////////

    // Get image sequence option
    if( !bGetImageSequenceInformationFromCommandLine() )
    {
        sprintf(acTemp, "Option %s must be specified!", c_pcImageSequence);
        sError = acTemp;

        DTREK_ERROR(1, sError);
    }

    // Call Evaluator to process ALL images
    Cscan       oScan(*m_poHeader);
    Cstring     sName("");
    for(int ii=0; ii < m_vecScanImageSequenceRanges.size(); ii++)
    {
        for(int jj = m_vecScanImageSequenceRanges[ii].first; jj <= m_vecScanImageSequenceRanges[ii].second; jj++ )
        {
            oScan.vSetSeqNum(jj);

            oScan.nGetImageName(&sName);
            
            m_poEvaluator->ProcessImage(sName.string());
        }
    }
    
    m_poEvaluator->vStartImageTestSequence(dTargetResolution, dTargetIOverSigma, false); // false means: do not reset max number of peaks in target shell.
                                                                                         // It was set earlier by calling vSetMaxObservedPeakCountInTargetShell()
    double  dRecommendedExposureTime_sec = 0.0;    
    m_poEvaluator->bSuggestImageExposureTime("", 
                                             dMaxSaturatedPeak_Percent,
                                             dRecommendedExposureTime_sec);

    int     nMaxReflnCount = m_poEvaluator->nGetMaxObservedPeakCountInTargetShell(); 
    m_poHeader->nReplaceValueDictionary("DTETE", "MAXIMUM_REFLN_COUNT_IN_TARGET_SHELL", nMaxReflnCount);

    m_poHeader->nReplaceValueDictionary("DTETE", "RECOMMENDED_EXP_TIME", dRecommendedExposureTime_sec);

    m_poHeader->nWrite(m_sHeaderOut);
    
    printf("\n\nRecommended optimal exposure time: %.1f sec\n\n", dRecommendedExposureTime_sec);

    printf("\ndtete: Done.\n");

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainExposureTimeEvaluator::vError(const int nErrorNum, const Cstring& sMessage)
{
    if( 0 != nErrorNum ) 
    {
        Cstring     sErrorMsg = "\nERROR: ";
        sErrorMsg += sMessage;
        sErrorMsg += "\n";

        printf("%s\n",sErrorMsg.string());
    }

    printf("\ndtete suggests an optimal exposure time based on sample image(s) and user\n"
           "targets for percentage of saturated reflections and I/sigma in a target\n"
           "resolution shell. It is best to use dtete iteratively, i.e. take sample\n"
           "image(s), have dtete suggest exposure time, take sample image(s) again, now\n"
           "with the suggested exposure time, call dtete again using the output header\n"
           "from the previous run etc. until the suggested exposure time values converge.\n\n");
    printf("\ndtete - Usage:\n"
         "dtete header_file [options...]\n\n"
         "Command line options: Description:\n\n");

    printf(" header_file        Qualified path of a file with a d*TREK header\n\n"
           " -seq nStart nEnd   Starting and ending image sequence values\n"
           "                    in the scan, specified by the header\n\n"
           " -sat    value      Target maximum percentage of saturated reflections.\n\n"
           " -reso   value      Target resolution value (A).\n\n"
           " -sigma  value      Target I/Sigma value for the target resolution.\n\n"
           " -start             Start a new sample testing sequence. This option\n"
           "                    should be used when calling dtete for the first time\n"
           "                    on a given sample. It signals the program to start\n"
           "                    keeping track of the maximum number of observed peaks\n"
           "                    in the target resolution shell. Subseqent calls to dtete\n"
           "                    on the same sample should not have that option.\n\n"
           " -out    string     Name of output file with the recommended exposure time.\n"
           "                    Default: *dtete.head.\n\n"
           //" -ref sReflnlist    Write the selected reflections to the file sReflnlist.\n\n"
           " -help              Print this help text.\n\n");

    printf("Examples:\n"
           "dtete mydata.img -seq 1 2 -seq 41 42 -start -sat 0.1 -reso 2 -sigma 50\n"
           "-out dtete.head\n\n");

    fflush(stdout);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CDTMainExposureTimeEvaluator::bGetImageSequenceInformationFromCommandLine()
{
    m_vecScanImageSequenceRanges.clear();

    /////////////////////////////////////////////////////////////////////////////////////////
    // Parse the sequence information and fill up the sequence array with "begin-end" pairs
    Cstring                 strSeqArgs("");
    std::vector<int>   vecSeqs;    
    
    bool    bOptionSpecified = false;
    
    while( bGetCommandOption(c_pcImageSequence, strSeqArgs) )
    {
        bOptionSpecified = true;

        strSeqArgs.nListToVector(vecSeqs, " ");  
    
        if( 0 == vecSeqs.size() )
            continue; // useless seq option
        else if( 1 == vecSeqs.size() )
            m_vecScanImageSequenceRanges.push_back(std::pair<int,int>(vecSeqs[0],vecSeqs[0]));
        else if( 2 == vecSeqs.size() )
            m_vecScanImageSequenceRanges.push_back(std::pair<int,int>(vecSeqs[0],vecSeqs[1]));
    }

    return bOptionSpecified;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
