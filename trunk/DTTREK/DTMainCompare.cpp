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
// DTMainCompare.cpp  Initial author: RB                14-Jan-2005
//
//  This file contains the member functions of class CDTMainCompare
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

#include "DTMainCompare.h"
#include "Ccompare.h"

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
#endif

char* pcUsage = 
"Usage:  [header.head] file1.ref [file2.ref] [options]              \n"
"                                                                   \n"
" Program compares files file1.ref and file2.ref by trimming        \n"
" non-overlaped data so that a fair comparison is made.             \n"
" Only one reflection file need be specified.                       \n"
"                                                                   \n"
" * * * Action Options * * *                                        \n"
"                                                                   \n"
" -trim [fPercent]     Output trimmed reflection files.  fPercent   \n"
"                      percent of the data is rejected from one or  \n"
"                      both files (see -benefit) before outputing.  \n"
"                      The files are also trimmed so that the       \n"
"                      requirements specified by the -min-redundancy\n"
"                      -equal-redundancy and -apples-??? commands   \n"
"                      are met.  New reflection lists are generated.\n"
"                                                                   \n"
" -rejectplot sName [fMaxFraction] [fFractionStep]                  \n"
"                      Generates a plot of Rmerge vs. %% Rejections \n"
"                      for both data sets and places the result in  \n"
"                      the data file 'sName'.                       \n"
"                      Default: fMaxFraction=1.0 fFractionStep=0.01 \n"
"                                                                   \n"
" -intensityplot sName [fRejectPercent] [nNumberBins]               \n" 
"                      Generates a plot of Rmerge vs. Intensity for \n"
"                      both data sets and places the result in the  \n"
"                      data file 'sName'.  An optional number       \n"
"                      'fRejectPercent' will reject data before     \n"
"                      making the plot.                             \n"
"                                                                   \n"
" * * * Settings Options * * *                                      \n"
"                                                                   \n"
" -apples-oranges      No attempt is made to find a 1-1 mapping     \n"
"                      of reflections in file1.ref to file2.ref     \n"
"                                                                   \n"
" -apples-macintoshes  Every group of equivalent reflections in     \n"
"                      file1.ref must exist in file2.ref and vica   \n"
"                      versa. If using -equal-redundancy, matching  \n"
"                      groups of reflections are trimmed to have the\n"
"                      same number of reflections.                  \n"
"                      This option is the default                   \n"
"                                                                   \n"
" -apples-apples       Every unique HKL reflection in file1.ref must\n"
"                      exist in file2.ref and vica versa.  If using \n"
"                      -equal-redundancy, redundant reflections are \n"
"                      excluded until equal redundancy is achieved. \n"
"                                                                   \n"
" -min-redundancy nN   Set the minimum redundancy of each group of  \n"
"                      reflections contributing to the R-merge.     \n"
"                      Default: nN = 2                              \n"
"                                                                   \n"
" -equal-redundancy    Determines if data sets should be trimmed so \n"
"                      that equal redundancy is used for all sets of\n"
"                      reflections contributing to the R-merge.     \n"
"                      This is dependent upon the setting of the    \n"
"                      -apples-??? commands.                        \n"
"                      Default:  Disabled                           \n"
"                                                                   \n"
" -benefit nDataSet    Use to specify which data set should be      \n"
"                      given the benefit of the doubt.  The file    \n"
"                      with the benefit of the doubt gets it's      \n"
"                      worst reflections removed.                   \n"
"                      nDataSet = 1:  Reject data from file1.ref    \n"
"                      nDataSet = 2:  Reject data from file2.ref    \n"
"                      nDataSet = 0:  Reject 50/50 from both.       \n"
"                      Default:  nDataSet = 0.                      \n"
"                                                                   \n"
" -dupsigmas           Use to set the fSigmaI field to the observed \n"
"                      statistical sigma in the output lists        \n"
"                      generated by the -trim command.              \n"
" -usesigmas nDataSet  Use the sigma values from one of the data    \n"
"                      sets.                                        \n"
"                                                                   \n"
;


int CDTMainCompare::nExecute(unsigned int argc,char* argv[]) 
{
    int nStat;
    int nPass;
    int nx;
    Cstring sTemp,sTemp2;
    Cstring sArg;

    int nArg,nArgStart,nArgEnd,nArgCount;
    int nActionCount = 0;

    Creflnlist* poReflnlist1    = NULL;
    Creflnlist* poReflnlist2    = NULL;
    Cimage_header* poHeader     = NULL;
    Ccompare* poCompare         = NULL;
    bool bEqualReflectionLists  = FALSE;
    Cstring sFileName1;
    Cstring sFileName2;
    Cstring sHeaderName;

   

    vDtrekSetModuleName("dtcompare");
    cout << "\ndtcompare: Copyright (c) 2006 Rigaku\n";
    cout << D_K_DTREKVersion << endl << flush;
    cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;

    nStat = 1;
    argc--;
    argv++;
    if (argc<2) {
        printf(pcUsage);
        goto error_place;
    };

    for (nPass=0;nPass<2;nPass++) {

        // Settings.
        eComparisonType     eType                   = APPLES_MAC;
        int                 nMinRedundancy          = 2;
        bool                bEqualRedundancy        = FALSE;
        bool                bDupSigmas              = FALSE;
        int                 nBenefit                = 0;
        int                 nUseSigmasIn            = 0;

        

        if ((nPass==1) && (nActionCount==0)) {
            printf("No action options specified!\n");
            goto error_place;
        } else
            nActionCount = 0;
        
        nArg = 0;
        sHeaderName = argv[nArg++];
        if (!sHeaderName.contains(".head")) {
            nArg--;
            sHeaderName = "";
        };
        sFileName1 = argv[nArg++];
        if ((argv[nArg]!=NULL) && (argv[nArg][0]!='-'))
            sFileName2 = argv[nArg++];
        else
            sFileName2 = sFileName1;

        // If this is the second pass, then we actually open the data items.
        if (nPass==1) {

            if (sHeaderName.length()) {
                poHeader = new Cimage_header(sHeaderName);
                if (!poHeader->bIsAvailable()) {
                    printf("Could not open header file '%s'\n",sHeaderName.string());
                    goto error_place;
                };            
            };
            
            poReflnlist1 = new Creflnlist(sFileName1);
            if (!poReflnlist1->bIsAvailable()) {
                printf("Could not open reflection file '%s'\n",sFileName1.string());
                goto error_place;
            };

            if (sFileName1 == sFileName2) {
                printf("Reflection lists are equal.\n");
                bEqualReflectionLists = TRUE;
                poReflnlist2 = poReflnlist1;
            } else {
                poReflnlist2 = new Creflnlist(sFileName2);
                if (!poReflnlist2->bIsAvailable()) {
                    printf("Could not open reflection file '%s'\n",sFileName2.string());
                    goto error_place;
                };
            };
            if (!poHeader) {
                Ccrystal* poCrystal;

                poCrystal = poReflnlist1->poGetCrystal();
                if (!poCrystal) {
                    poCrystal = poReflnlist2->poGetCrystal();
                    if (!poCrystal) {
                        printf("Could not extract embedded crystal object from file '%s' or file '%s'.\n",sFileName1.string(),sFileName2.string());
                        goto error_place;
                    } else {
                        printf("Using crystal found in file '%s'.\n",sFileName2.string());
                        poHeader = poReflnlist2->poGetHeader();
                        delete poCrystal;
                    };
                } else {
                    printf("Using crystal found in file '%s'.\n",sFileName1.string());
                    poHeader = poReflnlist1->poGetHeader();
                    delete poCrystal;
                };
            };
            poCompare = new Ccompare(*poReflnlist1,*poReflnlist2,*poHeader);
            if (!poCompare->bIsAvailable()) {
                printf("Could not initialize comparison");
                goto error_place;
            };
        } else {
            if (!bFileExists(sFileName1)) {
                printf("File '%s' does not exist.\n",sFileName1.string());
                goto error_place;
            };
            if (!bFileExists(sFileName2)) {
                printf("File '%s' does not exist.\n",sFileName2.string());
                goto error_place;
            };
        };

        
        for (;nArg<argc;) {
            if (argv[nArg][0]!='-') {
                printf("- must preceed option %s\n",argv[nArg]);
                goto error_place;
            };
            
			if( bIsSkipCommandLineArgument(argv[nArg]) )
			{
				nArg++;
				continue;
			}

			sArg = argv[nArg];
            if (nPass == 1)
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
                    if (nPass==1)
                        printf("%s ",argv[nArg]);
                };
            };
            nArgCount = nArgEnd-nArgStart+1;
            
            if (nPass==1)
                printf("\n");

            if (sArg == "-apples-oranges") {
                eType = APPLES_ORANGES;
            } else if (sArg == "-apples-macintoshes") {
                eType = APPLES_MAC; 
            } else if (sArg == "-apples-apples") {
                eType = APPLES_APPLES;
            } else if (sArg == "-min-redundancy") {
                if ((nArgCount!=1) || (1!=sscanf(argv[nArgStart],"%d",&nMinRedundancy))) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };
            } else if (sArg == "-equal-redundancy" ) {
                bEqualRedundancy = TRUE;
            } else if (sArg == "-benefit") {                
                if ((nArgCount!=1) || (1!=sscanf(argv[nArgStart],"%d",&nMinRedundancy)) ||
                    (nBenefit>2) || (nBenefit<0)) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };
            } else if (sArg == "-trim") {
                double fReject;
                
                nActionCount++;

                if (nArgCount==0) {
                    fReject = 0.0;
                } else {
                    if (1!=sscanf(argv[nArgStart],"%lf",&fReject)) {
                        printf("Could not read option '%s'\n",sArg.string());
                        goto error_place;
                    };
                };
                if (nPass == 1) {
                    poCompare->nClearRejects();
                    if (poCompare->nEquivocate(eType,bEqualRedundancy,nMinRedundancy,nUseSigmasIn)) {
                        nStat = 1;
                        goto error_place;
                    };
                    poCompare->nComputeRMerge(nMinRedundancy,bDupSigmas);

                    if (fReject!=0.0) {
                        if (nBenefit==0) {                            
                            poCompare->nReject(fReject*0.5,0);
                            poCompare->nReject(fReject*0.5,1);
                        } else if (nBenefit==1) {
                            poCompare->nReject(fReject*0.5,0);
                        } else if (nBenefit==2) {
                            poCompare->nReject(fReject*0.5,1);
                        };
                        poCompare->nEquivocate(eType,bEqualRedundancy,nMinRedundancy,nUseSigmasIn);
                        poCompare->nComputeRMerge(nMinRedundancy,bDupSigmas);
                    };

                    sTemp = sFileName1.before(".");
                    sTemp += "_comp.ref";
                    poCompare->nWriteNonRejected(0,sTemp);
                    if (!bEqualReflectionLists) {
                        sTemp = sFileName2.before(".");
                        sTemp += "_comp.ref";
                        poCompare->nWriteNonRejected(1,sTemp);
                    };                       
                };
            } else if (sArg == "-rejectplot") {
                Cstring sName;
                double fReject;
                double fMaxFraction;
                double fFractionStep;
                FILE* pFOut;
                
                nActionCount++;
                
                if (nArgCount>3) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };
                fMaxFraction = 1.0;
                fFractionStep = 0.01;
                if ((nArgCount>=2) && (1!=sscanf(argv[nArgStart+1],"%lf",&fMaxFraction))) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };
                if ((nArgCount>=3) && (1!=sscanf(argv[nArgStart+2],"%lf",&fFractionStep))) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };

                sName = argv[nArgStart];

                if (nPass == 1) {
                    poCompare->nClearRejects();
                    if (poCompare->nEquivocate(eType,bEqualRedundancy,nMinRedundancy,nUseSigmasIn)) {
                        nStat = 1;
                        goto error_place;
                    };
                    pFOut = fopen(sName.string(),"w+t");
                    if (!pFOut) {
                        printf("Could not open file '%s' for writting.\n",sName.string());
                    } else {
                        poCompare->nComputeRMerge(nMinRedundancy,bDupSigmas);
                        for (fReject=0.0;fReject<fMaxFraction;fReject+=fFractionStep) {
                            if (nBenefit==0) {
                                poCompare->nReject(fReject*0.5,0);
                                poCompare->nReject(fReject*0.5,1);
                            } else if (nBenefit==1) {
                                poCompare->nReject(fReject*0.5,0);
                            } else if (nBenefit==2) {
                                poCompare->nReject(fReject*0.5,1);
                            };
                            poCompare->nEquivocate(eType,bEqualRedundancy,nMinRedundancy,nUseSigmasIn,2);
                            poCompare->nComputeRMerge(nMinRedundancy,bDupSigmas);
                            fprintf(pFOut,"%f %f %f\n",fReject,poCompare->m_a2fRMerge[0],poCompare->m_a2fRMerge[1]);
                        };
                    };
                    if (pFOut)
                        fclose(pFOut);
                };

            } else if (sArg == "-intensityplot") {
                double fReject;
                Cstring sName;
                float fHighestIntensity;
                int nIntensityBins;

                fHighestIntensity = 0.95f;
                nIntensityBins = 20;
                fReject = 0.0;
                
                nActionCount++;

                if ((nArgCount<1) || (nArgCount>3)) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };

                if ((nArgCount>=2) && (1!=sscanf(argv[nArgStart+1],"%lf",&fReject))) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };
                if ((nArgCount>=3) && (1!=sscanf(argv[nArgStart+2],"%d",&nIntensityBins))) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };

                sName = argv[nArgStart];

                if (nPass == 1) {
                    if (poCompare->nEquivocate(eType,bEqualRedundancy,nMinRedundancy,nUseSigmasIn)) {
                        nStat = 1;
                        goto error_place;
                    };
                    poCompare->nClearRejects();
                    poCompare->nComputeRMerge(nMinRedundancy,bDupSigmas);
                    if (nBenefit==0) {
                        poCompare->nReject(fReject*0.5,0);
                        poCompare->nReject(fReject*0.5,1);
                    } else if (nBenefit==1) {
                        poCompare->nReject(fReject*0.5,0);
                    } else if (nBenefit==2) {
                        poCompare->nReject(fReject*0.5,1);
                    };                    
                    if (poCompare->nEquivocate(eType,bEqualRedundancy,nMinRedundancy,nUseSigmasIn)) {
                        nStat = 1;
                        goto error_place;
                    };
                    poCompare->nInitIntensityBins(fHighestIntensity,nIntensityBins);
                    poCompare->nComputeRMerge(nMinRedundancy,bDupSigmas);
                    poCompare->nWriteIntensityBins(sName);
                };
            } else if (sArg == "-dupsigmas") {
                bDupSigmas = TRUE;
            } else if (sArg == "-usesigmas") {
                if ((nArgCount!=1) || (1!=sscanf(argv[nArgStart],"%d",&nx))) {
                    printf("Could not read option '%s'\n",sArg.string());
                    goto error_place;
                };
                nUseSigmasIn = nx;
            } else {
                printf("Could not recognize option '%s'\n",sArg.string());
                goto error_place;
            };
        };
    };
        
  
    nStat = 0;
error_place:

    delete poReflnlist1;
    if (!bEqualReflectionLists)
        delete poReflnlist2;
    delete poHeader;
    delete poCompare;
    return nStat;   
}

void CDTMainCompare::vError(const int nErrorNum, const Cstring& sMessage)
{
}






