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
// DTMainProfit.cpp  Initial author: RB                14-Jan-2005
//
//  This file contains the member functions of class CDTMainProfit
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

#include "DTMainProfit.h"

#include "Cprofit2D.h"
#include "Cscan.h"
#include "CScanBitmap.h"

#include <string.h>

// Add the following three lines
#ifdef SSI_PC
#include "CrclIncludes.h"
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

char* g_pcUsage = 

"dtprofit Usage:                                                      \n"
"dtprofit [options]                                                   \n"
"Command line options: Description:                                   \n"
"                                                                     \n"
" -map [sHeaderFile] [nFirstImage nLastImage]                         \n"
"                                                                     \n"
"            Maps shoeboxes onto images.  The output images will have \n"
"            the same name as the integrated images (specified by the \n"
"            scan template) with _map extended to the base name.      \n"
"                                                                     \n"
" -profit [sHeaderFile] [sIntegrateFile] [sPartialsFile] [sScaleFile] \n"
"            Attempt to profile fit the data so that each profile     \n"
"            uses at least nNumReflections.                           \n"
"            Profiles will not be constructed for spots which lie more\n"
"            than nMaxImageRangePerProfile apart.                     \n"
"                                                                     \n"
" -postrefine [nMinRefsPerProfile] [nResolutionBins]                  \n"
"            Post refine the data.  We require that a minimum of      \n"
"            nMinRefsPerProfile reflections are used for every        \n"
"            mosaicity profile and that at least nResolutionBins data \n"
"            bins are used.                                           \n"
"                                                                     \n"
" -checkoverlaps fMaximumFractionOverlap                              \n"
"            Reject spots that have a significant portion of their    \n"
"            integrated pixels intersecting other neighboring spots.  \n"
"            A spot is flagged as overlapped if the fraction of       \n"
"            intensity shared by other neighboring spots is greater   \n"
"            than or equal to  fmaximumFractionOverlap.               \n"
"            Default:  fMaximumFractionOverlap=0.03                   \n"
"                                                                     \n"
" -rmergemosaicity [fFixMosaicity]                                    \n"
"            Generate mosaicity plots of Rmerge vs mosaicity for      \n"
"            resolution shells on the detector.  If fFixMosaicity is  \n"
"            specified, cut off the mosaicity at fFixMosaicity.       \n"
"                                                                     \n"
" -ellipsoids                                                         \n"
"            Generate ellipsoid file ellipsoids.ref                   \n"
" -prefix sPrefix                                                     \n"
"           Equivalent to setting the DTREK_PREFIX environment        \n"
"           variable                                                  \n"
"                                                                     \n"
" Examples:                                                           \n"
" dtprofit -prefix 1_ -profit                                         \n"
" dtprofit -profit 1_dtintegrate.head 1_dtintegrate.ref 1_dtintpartials.ref\n"
;


#define SYNTAXERROR(x) { if ((x)==0) printf("Syntax error with option '%s'\n",sArg.string()); else printf("Syntax error with option '%s' (near '%s')\n",sArg.string(),argv[nArgStart+(x)-1]); return 1; }


int CDTMainProfit::nExecute(unsigned int argc,char* argv[]) 
{

    int nArg,nArgCount,nArgStart,nArgEnd;
    Cstring sArg;
    Cstring sTemp;
    Cstring sTemp2;
    int nStat = 0;
    double fTemp;
    int nTemp;

    // User defined settings.
    bool    bDoProfit = false;
    bool    bMapImages_map = false;
    bool    bWriteEllipsoids = false;
    bool    bCalcWidthsHKL = false;
    int     nFirstImage_map;
    int     nLastImage_map;
    double  fCalcWidthsHKLFixedMosaicity = 0.0;

    Cprofit2DSlices oProfit;


    Cstring sHeaderFile_map;
    Cstring sProfitFile = "";
    Cstring sIntegrateFile = "";
    Cstring sPartialsFile = "";
    Cstring sScaleFile = "";
    Cstring sIntegrateHeader = "";
    Cstring sStub = sTransSymbol(sDtrekGetPrefix());

    vDtrekSetModuleName("dtprofit");
    vPrintCopyrightInfo();
    //cout << "\ndtprofit:  Copyright (c) 2006 Rigaku\n";
    //cout << D_K_DTREKVersion << endl;
    // Copy command line to output log

    cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;

    sTemp = (const char*) argv[1];
    if ("-help" == sTemp || "-h" == sTemp || argc<2) {
        printf(g_pcUsage);
        return 0;
    };

    argv++;
    argc--;

    // Process Options.
    for (nArg=0;nArg<argc;) {
        if (argv[nArg][0]!='-') {
            printf("- must preceed option %s\n",argv[nArg]);
            return 1;
        };
        sArg = argv[nArg];
        printf("Processing Option: %s ",argv[nArg]);
        for (nArg++,nArgStart=nArg,nArgEnd=nArgStart-1;nArg<argc;nArg++) {
            bool bIsMinus;
            bIsMinus = FALSE;
            
            if ((argv[nArg][0]=='-') && (strspn(argv[nArg],"-0123456789.+:")!=strlen(argv[nArg])))  {
                // there are certain keywords that also have a '-' in front of them.
                sTemp = argv[nArg];
                sTemp.upcase();
                bIsMinus = TRUE;
            };
            if (bIsMinus)                  
                break;
            else {
                nArgEnd=nArg;
                printf("%s ",argv[nArg]);
            };
        };
        nArgCount = nArgEnd-nArgStart+1;

        printf("\n");

        if (sArg == "-map") {
            sTemp = argv[nArgStart];
            sTemp.downcase();
            bMapImages_map = true;

            if (nArgCount == 3) {
                sHeaderFile_map = argv[nArgStart + 0];
                if ((1!= sscanf(argv[nArgStart + 1],"%d",&nFirstImage_map)) ||
                    (1!= sscanf(argv[nArgStart + 2],"%d",&nLastImage_map))) 
                    SYNTAXERROR(1);
            } else if (nArgCount == 2) {
                sHeaderFile_map = "dtintegrate.head";
                if ((1!= sscanf(argv[nArgStart + 0],"%d",&nFirstImage_map)) ||
                    (1!= sscanf(argv[nArgStart + 1],"%d",&nLastImage_map))) 
                    SYNTAXERROR(1);
            } else if (nArgCount == 1) {
                sHeaderFile_map = argv[nArgStart + 0];
                nFirstImage_map = -1;
                nLastImage_map = -1;
            } else {
                sHeaderFile_map = "dtintegrate.head";
                nFirstImage_map = -1;
                nLastImage_map = -1;
            };
        } else if (sArg == "-ellipsoids") {
            bWriteEllipsoids = true;
        } else if (sArg == "-gain") {
            if ((nArgCount != 1) || (1 !=sscanf(argv[nArgStart],"%lf",&fTemp)))
                SYNTAXERROR(0);
            oProfit.m_fGainEst = fTemp;
        } else if (sArg == "-prefix") {
            if (nArgCount != 1)
                SYNTAXERROR(0);
            sStub = argv[nArgStart];
            nPutEnv((Cstring) "DTREK_PREFIX",sStub);
        } else if (sArg == "-profit") {
            if (nArgCount >= 1)
                sIntegrateHeader = argv[nArgStart + 0];
            if (nArgCount >=2)
                sIntegrateFile = argv[nArgStart + 1];
            if (nArgCount >= 3) 
                sPartialsFile = argv[nArgStart + 2];
            if (nArgCount >= 4)
                sScaleFile = argv[nArgStart + 3];
            if (nArgCount >= 5)
                SYNTAXERROR(0);
            
           bDoProfit = true;
        } else if (sArg == "-rmergemosaicity") {
            if ((nArgCount == 1) && (1 == sscanf(argv[nArgStart],"%lf",&fTemp)))
                fCalcWidthsHKLFixedMosaicity = fTemp;
            bCalcWidthsHKL = true;
        } else if (sArg == "-postrefine") {
            if ((nArgCount >= 1) && (1 == sscanf(argv[nArgStart + 0],"%d",&nTemp))) {
                oProfit.m_nMinRefsPerProfile = nTemp;
                if ((nArgCount >= 2) && (1 == sscanf(argv[nArgStart + 1],"%d",&nTemp))) {
                    oProfit.m_nResoGroups = nTemp;
                    
                }
            }
        } else if (sArg == "-checkoverlaps") {
            if ((nArgCount >=1) && (1== sscanf(argv[nArgStart + 0],"%lf",&fTemp))) {
                oProfit.m_fIntersectMax = fTemp;
            } else
                SYNTAXERROR(0);
        } else
            SYNTAXERROR(0);

    };
    
    if (sIntegrateFile == "") {
        sTemp = sStub;
        sTemp += "dtintegrate.ref";
        sIntegrateFile = sTemp;
    };
    if (sPartialsFile == "") {
        sTemp = sStub;
        sTemp += "dtintpartials.ref";
        sPartialsFile = sTemp;
    };
    if (sIntegrateHeader == "") {
        sTemp = sStub;
        sTemp += "dtintegrate.head";
        sIntegrateHeader = sTemp;
    };
    
    if (sProfitFile == "") {
        sTemp = sStub;
        sTemp += "dtprofit.ref";
        sProfitFile = sTemp;
    };


    if (bDoProfit) {
        nStat = oProfit.nLoadHeader(sIntegrateHeader);
        

        int nFileParts = oProfit.nSplitFileProcessSize(sIntegrateHeader,g_nMaxReflnsToProfit);
        int nFirstPart,nLastPart;
        int nFilePart;

        nFirstPart = min(nFileParts,0);
        nLastPart = max(0,nFileParts) - 1;
        for (nFilePart = nFirstPart; nFilePart <= nLastPart; nFilePart++) {
            if (!nStat)
                nStat = oProfit.nLoad(sIntegrateFile,sPartialsFile,sScaleFile,nFilePart);
            oProfit.vSetIntegrateOutputName(sProfitFile);
            if (!nStat)
                nStat = oProfit.nSyncLists(false);
            if (!nStat)
                nStat = oProfit.nPrintAxial(nFilePart==nLastPart);
            if (!nStat)
                nStat = oProfit.nRefineMosaicityModel();
            if (!nStat)
                nStat = oProfit.nPrintMosaicityModel(nFilePart);
            if ((!nStat) && (bCalcWidthsHKL))
                nStat = oProfit.nCalcWidthsHKL(fCalcWidthsHKLFixedMosaicity);
            if (!nStat)
                nStat = oProfit.nFindIntersections();
            if (!nStat)
                nStat = oProfit.nCalcRemoveNoise(nFilePart == nFirstPart,nFilePart == nLastPart);
            if (!nStat)
                nStat = oProfit.nWrite(nFilePart,true,false,bWriteEllipsoids);  
        };
        oProfit.nWriteParts(true,false,bWriteEllipsoids);
    };

    
    // If user requested an image map, generate it.
    if (bMapImages_map)
    {
        CScanBitmap oScanMap;
        Cstring sOutName;
        
        if (sHeaderFile_map.length()) {
            Cimage_header oHeader(sHeaderFile_map);
            Cscan oScan(oHeader);
            if ((oHeader.bIsAvailable()) && (oScan.bIsAvailable()))
                nStat = oScanMap.nOpenInput(oHeader);
            else
                nStat = 1;

            if (!nStat) {
                if (nFirstImage_map < 0) {
                    nFirstImage_map = oScan.nGetSeqNum();
                    nLastImage_map = nFirstImage_map + oScan.nGetSeqNumImages() - 1;
                };
                
                
                oScan.vSetSeqNum(nFirstImage_map);
                nStat = 0;
                while (!nStat) {
                    nStat = oScanMap.nReadBitmapForScan(oScan);

                    if ((!nStat) && (bMapImages_map)) {
                        sTemp = oScan.sGetImageName();
                        Cimage oImageOld(sTemp);
                        if (oImageOld.bIsAvailable()) {
                            sTemp2 = sFileGetBasename(sTemp) + "_map";
                            oScanMap.nWriteMaskedImage(sTemp2,&oImageOld);
                        };
                    } 
                    if (oScan.nGetSeqNum() == nLastImage_map)
                        break;
                    oScan.vNextSeqNum();
                };
            }           
        };
    };

    if (sStub) {
        nPutEnv((Cstring) "DTREK_PREFIX",(Cstring) "");
    };

    cout << "\ndtprofit: Done.\n" << flush;
    return (nStat);
}

void CDTMainProfit::vError(const int nErrorNum, const Cstring& sMessage)
{
}

