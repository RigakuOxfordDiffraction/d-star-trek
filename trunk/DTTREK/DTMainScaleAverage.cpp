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
// DTMainScaleAverage.cpp       Initial author: RB     20-Oct-2004
// This file contains the member functions of class CDTMainScaleAverage

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

#ifdef SSI_PC
#include "CrclIncludes.h"
#endif

#ifdef SSI_PC
  #define DTREKEXIT(x) exit(x)
#elif defined DEBUG
  #include <conio.h>
  #define DTREKEXIT(x) { getch(); return(x); }
#elif !defined(EXIT)
  #define DTREKEXIT(x) return(x)
#endif

#include "CCommandLine.h"
#include "DTMainScaleAverage.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::ios;
#endif

char* pcHelp =
" dtscaleaverage [input.head] [options] ref-file1 [ ref-file2 ...]       \n"
"                [ output.ref ]                                          \n"
"                                                                        \n"
" Note:  Output reflection file should be last argument on command line. \n"
"                                                                        \n"
" Options are:                                                           \n"
"                                                                        \n"
" -help | -h             Print this help screen                          \n"
" -help errormodel       Print help on error modeling options.           \n"
" -help reject           Print help on reflection rejection options.     \n"
"                                                                        \n"
" -spacegroup nNum       Specify spacegroup.  If header is absent,       \n"
"                        spacegroup number must be specified for         \n"
"                        meaningful output.                              \n"
"                                                                        \n"
" -cell a b c al bt gm   Specify unit cell.  If header is absent, failing\n"
"                        to specify a cell will result in meaningless    \n"
"                        resolution data.                                \n"
"                                                                        \n"
" -scaleanom             Specify that I+ != I-.  Affects absorption      \n"
"                        correction as well as batch scaling.            \n"
"                        Default: I+ == I-                               \n"
"                                                                        \n"
" -batchscale            Enables and runs batch scaling. This utilizes   \n"
"                        reflections which are within limits             \n"
"                        imposed by -reject sigma. Use -fix to fix a      \n"
"                        batch scaling multiplier.                       \n"
"                                                                        \n"
" -fix sBatch            Specify which batch to fix during batch scaling.\n"
"                        Default: Fix the batch that has the most        \n"
"                        overlaps with other batches.  Note that this    \n"
"                        option affects ALL -batchscale options.         \n"
"                                                                        \n"
" -bfactor               Enables B factor scaling. Only reflections which\n"
"                        are within limits imposed by -reject sigma are  \n"
"                        used.                                           \n"
"                                                                        \n"
" -batchrestrain [options]                                               \n"
"                        Restrains scale factors for specified batch     \n"
"                        ranges, using specified weights.                \n"
"                        If no options are provided, all batches within  \n"
"                        each scan will be restrained with default weight.\n"
"                        The default weight is 0.1.                      \n"
"                        Allowed options are: W |                        \n" 
"                        sB1s-sB1e[-W1],sB2s-sB2e[-W2],...,sBNs-sBNe[-WN]\n"
"                        where W - a common weight to restrain all scans,\n"
"                        sBNs,sBNe are start and end batch names for N-th \n"
"                        batch range; WN - restraining weight for N-th   \n"
"                        batch range. Whenever a weight is not provided, \n"
"                        the default weight will be used.                \n"
"                        Note: you cannot restrain batches across        \n"
"                        different scans.                                \n"
"                                                                        \n"
" -reso fMin fMax        Resolution limits of input reflections to use.  \n"
"                        Default: 1000 0.0                               \n"
"                                                                        \n"
" -reso edge             Set upper resolution to be equal to the highest\n"
"                        *minimum* resolution found on any of the four\n"
"                        edges (sides) of the image.\n\n"
"                                                                        \n"
" -reso corner           Set upper resolution to be equal to the highest\n"
"                        *maximum* resolution found on any of the four\n"
"                        edges (sides) of the image.\n\n"
"                                                                        \n"
" -countoverlap          Count and show overlaps among the batches. This \n"
"                        option will also write out a file named         \n"
"                        'batchoverlaps.img' that contains a graphical   \n"
"                        map of the intersection matrix.                 \n"
"                        Default: Disabled.                              \n"
"                                                                        \n"
" -sabsorb [options]     Spherical absorption correction (an alternate   \n"
"                        empirical absorption algorithm).                \n"
"                        Type '-help sabsorb' for help on options.       \n"
"                        Default: Do not run absorption correction.      \n"
"                                                                        \n"
" -errormodel [options]  Error modeling options.                         \n"
"                        Type -errormodel off' to avoid error modeling   \n"
"                        which is an adjustment of the input sigmas.     \n"
"                        Type '-help errormodel' for help on options.    \n"
"                                                                        \n"
" -chisq fDesiredChiSq   Error model (Emul, Eadd) will be adjusted so that\n"
"                        the overall reduced Chi-squared nearly equals the\n"
"                        fDesiredChiSq value. If this option is NOT used,\n"
"                        then the default target or desired reduced Chi- \n"
"                        squared is set to 1 unless -anom is used.  In the\n"
"                        the case where -anom is used, the target of 1.2 \n"
"                        is used, since then the reduced Chi-Squared_Anom\n"
"                        will be closer to 1.                          \n\n"
" -reject  [options]     Rejection options.                              \n"
"                        Determines which reflections are rejected.      \n"
"                        Type '-help reject' for help on options.        \n"
"                                                                        \n"
" -partialrefs sFile     After scaling, attempt to re-scale reflections  \n"
"                        from the list of partials.                      \n"
"                                                                        \n"
" -list sTable           List a table. The (case insensitive) options are\n"
"                           All                 Print all tables.        \n"
"                           RMergeVsBatch       Rmerge vs Batch.         \n"
"                           RMergeVsIntensity   RMerge vs Intensity.     \n"
"                           RMergeVsResolution  RMerge vs Resolution.    \n"
"                           Multiplicity        Print Multiplicities.    \n"
"                           Completeness        Data Completeness.       \n"
"                           ChiSquaredDistrib   Chi Squared Distribution.\n"
"                           InputStats          Input statistics tables. \n"
"                        By default, all tables (with the exception of   \n"
"                        'InputStats') are printed if a -list option is  \n"
"                        omitted.                                        \n"
"                                                                        \n"
" -ref sUnavgReflnFile   This options writes the scaled, but unaveraged, \n"
"                        reflection list to the file sUnavgReflnFile.    \n"
"                        Rejected reflections will have standard         \n"
"                        deviations less than 0.                         \n"
"                        Default: Do not write such a file.              \n"
"                        See also -texsan, -texsan2, and -shelx          \n"
"                                                                        \n"
" -anom                  Output anomalous information as Intensity,      \n"
"                        sigmaI,Intensity+,sigmaI+,Intensity-,sigmaI-.   \n"
"                        Use -1 for missing information.                 \n"
"                        Default: Output Intensity,sigmaI; do not        \n"
"                        output I+,sigI+,I-,sigI- on the same line.      \n"
"                                                                        \n"
" -texsan                Output only h,k,l,Intensity,sigmaI fields with  \n"
"                        no header and no rejected measurements to the   \n"
"                        scaled, unaveraged file specified by the -ref   \n"
"                        option (default f2plus.dat).                    \n"
"                                                                        \n"
" -texsan2               Output only h,k,l, unscaled Intensity, unscaled \n"
"                        SigmaI, unscaled and non-LP-corrected Intensity,\n"
"                        transmission factor, scale factor to the file   \n"
"                        specified by the -ref option (default f2plus.dat).\n"
"                                                                        \n"
" -shelx                 Output only h,k,l,Intensity,SigmaI and a batchID\n"
"                        no header and no rejected measurements in the   \n"
"                        scaled, unaveraged file specified by the -ref   \n"
"                        option (default shelxl.hkl) in Fortran (what's  \n"
"                        that?) format (3I4, 2F8.2, I4) suitable for input\n"
"                        to SHELX.  The batchID is always the integer 1. \n"
"                                                                        \n"
" -noheader              Do NOT write standard d*TREK reflnlist header on\n"
"                        the scaled and averaged output file.            \n"
"                        Default: write a header.                        \n"
"                                                                        \n"
" -nounavgheader         Do NOT write standard d*TREK reflnlist header on\n"
"                        the unaveraged output file.                     \n"
"                        Default: write a header.                        \n"
"                                                                        \n"
" -out headerfile        Write a header to the headerfile at termination \n"
"                        of dtscaleaverage                               \n"
"                                                                        \n"
" -rebatch [imageboundary] nAvgOverlaps                                  \n"
"                        Rebatches all scans that do not have a sufficent\n"
"                        number of overlaps.  Rebatching will group      \n"
"                        batches together so that an averge of           \n"
"                        nAvgOverlaps occur between batches.             \n"
"                        Using the 'imageboundary' command will create   \n"
"                        batches that lie on image boundaries.           \n"
"                        See also:  -rebatch bookends                    \n"
"                                                                        \n"
" -rebatch bookends                                                      \n"
"                        The 'bookends' version of -rebatch reassigns    \n"
"                        only those batches who's size is less than one  \n"
"                        half the average batch size                     \n"
"                                                                        \n"
" -completenessreject fPercentComplete [+ | -] sSortField1 [sSortField2] \n"
"                        Determine the fraction of data needed to obtain \n"
"                        a particular value of the completeness.         \n"
"                        The data is sorted on reflection field          \n"
"                        sSortField1 (secondarily on field sSortField2)  \n"
"                        and rejections are made from the sorted list    \n"
"                        until the desired completeness is obtained.     \n"
"                        Specifying the [-] flag reverses the sort order \n"
"                        Example:                                        \n"
"                        -completenessreject 95.0 sBatch fObs_rot_mid    \n"
"                                                                        \n"
" Note: options, header and reflection files must occur first on the     \n"
"       the command line.  If header file does not end in .head extension\n"
"       then the header file must be specified first.  Order of options  \n"
"       on the command line dictates the order in which they are applied \n"
"       to the data.                                                     \n"
"\n"
" Example: dtscaleaverage dtintegrate.head dtprofit.ref                \\\n"
"                         -reject sigma 5                              \\\n"
"                         -reject symabs                               \\\n"
"                         -reject deviation 5                          \\\n"
"                         -reject fraction .005                        \\\n"
"                         -batchscale                                  \\\n"
"                         -batchrestrain 0001-0049-0.001,10051-10360   \\\n"
"                          dtscaleaverage.ref                            \n"
"\n REFERENCE - Pflugrath, JW (1999) Acta Cryst. D55, 1718-1725.\n";

char* pcSphericalAbsorbHelp = 
" Spherical absorption correction algorithm.                               \n"
" -sabsorb [nHarmonicOrderDiffracted] [nHarmonicOrderIncident] [nMaxIter]  \n"
"                                                                          \n"
" nHarmonicOrderDiffracted            Order of harmonics for modeling      \n"
"                                     diffracted intensities. Default: 4   \n"
"                                                                          \n"
" nHarmonicOrderIncident              Order of (planar)harmonics for       \n"
"                                     modeling incident beam absorption.   \n"
"                                     Default: 1                           \n"
"                                                                          \n"
" nMaxIter                            Maximum number of iterations to      \n"
"                                     execute for each full scan of data.  \n"
"                                     Default: 2                           \n"
"                                                                          \n"
;

char* pcErrorModelHelp =
" -errormodel [mul #fMul ] [add #fAdd] [addmult #fAddMult]                 \n"
"                                                                          \n"
" For a valid error model, the variance calculated for each reflection     \n"
" should match the observed variance obtained from multiple measurements   \n"
" of the reflection.  When this holds, a reduced chi squared (|ChiSq|)     \n"
" statistic has an expected value of 1.0.  Since calculated sigmas are     \n"
" affected by detector gain, correction parameters Eadd and Emul are       \n"
" required to scale the calculated sigmas to obtain |ChiSq| == 1.0.  The   \n"
" parameters Eadd and Emul can be set manually, or determined automatically.\n"
;


char* pcRejectHelp = 
" -reject symabs                                                         \n"
" -reject list sFile                                                     \n"
" -reject deviation fRejSigma                                            \n"
" -reject fraction  fFractionReject                                      \n"
" -reject sigma fExcludeSigma                                            \n"
"                                                                        \n"
" -rejectbatch percentbad fPercent                                       \n"
" -rejectbatch rmerge     fDeviation                                     \n"
" -rejectbatch intensity  fDeviation                                     \n"
" -rejectbatch chisq      fMaxChiSq                                      \n"
" -rejectbatch id sBatch                                                 \n"
" -rejectbatch id sBatch1,sBatch2,...,sBatchN                            \n"
" -rejectbatch id sBatch1s-sBatch1e,sBatch2s-sBatch2e,...,sBatchN        \n"
" -rejectbatch                                                           \n"
"                                                                        \n"
" Options for -reject:                                                   \n"
"                                                                        \n"
" deviation fRejSigma   Sets the rejection level for reflections.        \n"
"                       A reflection with scaled intensity that differs  \n"
"                       by more than fRejSigma from the weighted average \n"
"                       intensity calculated from other symmetry-related \n"
"                       reflections is flagged as rejected.              \n"
"                                                                        \n"
" fraction fFracReject  Sets the ChiSq cutoff to reject a fraction       \n"
"                       (fFractionReject) of the input reflections,      \n"
"                       e.g. 0.005 rejects 0.5%% of reflns.              \n"
"                                                                        \n"
" symabs                Symmetry-forbidden input reflections are disregarded\n"
"                       at all stages of scaling, averaging and outputting\n"
"																		 \n"
" sigma fExcludeSigma   Input reflections with I/sigmaI less than       \n"
"                       fExcludeSigma are excluded from contributing to \n"
"                       the scale factor refinement.  However, these    \n"
"                       reflections are included in the final statistics\n"
"                       Default: 3.                                     \n"
"																		 \n"
" list sFile            Exclude any reflections found in sFile from     \n"
"                       scaling, but not necessarily averaging.  Although\n"
"                       one would typically use the rejects list output \n"
"                       from a previous dtscaleaverage job for sFile, this\n"
"                       does NOT pre-reject these reflections.  While the\n"
"                       reflections are not pre-rejected, they may be   \n"
"                       rejected in the normal course of averaging and  \n"
"                       rejection.                                      \n"
"                                                                        \n"
" Options for -rejectbatch:                                              \n"
"                                                                        \n"
"                       Batch rejection allows you to reject an entire   \n"
"                       batch of reflections as opposed to the           \n"
"                       'per-reflection' rejection described above. A    \n"
"                       batch is rejected when too many reflections have \n"
"                       dubious statistics:                              \n"
"                                                                        \n"
" percentbad fRejects   If the percentage of reflections rejected in the \n"
"                       batch is greater than fPercentRejects, the entire\n"
"                       batch is rejected.                               \n"
"                                                                        \n"
" rmerge fValue         For fValue >= 2.0                                \n"
"                       If the Rmerge value of a batch is greater than   \n"
"                       fValue statistical sigmas above <Rmerge>,        \n"
"                       the batch is rejected.   <Rmerge> is the average \n"
"                       Rmerge value over all batches.                   \n"
"                                                                        \n"
"                       For fValue < 1.0                                 \n"
"                       If the Rmerge value is greater than fValue then  \n"
"                       the batch is rejected.                           \n"
"                                                                        \n"
" chisq fValue          For fValue > 2.0                                 \n"
"                       If the Chi-Squared value of a batch is greater   \n"
"                       than fValue statistical sigmas above the average \n"
"                       Chi-squared statistic of the entire data set,    \n"
"                       the batch is rejected.                           \n"
"                                                                        \n"
"                       For fValue < 0.0                                 \n"
"                       If the Chi-Squared value for the reflections     \n"
"                       in a single batch is greater than ABS(fValue),   \n"
"                       then the batch is rejected.                      \n"
"                                                                        \n"
" intensity fDeviation  If the <I/sig> value of a batch is greater than  \n"
"                       fDeviation statistical sigmas above <<I/sig>>,   \n"
"                       the batch is rejected.   <<I/sig>> is the average\n"
"                       <I/sig> value over all batches.                  \n"
"                                                                        \n"
" id sBatchTemplate     Reject reflns in specific batch name(s).  sBatchTemplate\n"
"                       can have wildcards ('?', '*') to reject a number \n"
"                       of batches.  Multiple batch id's must be separated\n"
"                       by a comma.  Ranges are valid as well (but not   \n"
"                       wildcards in ranges!). Only the last '-rejectbatch\n"
"                       id' is accepted.                                 \n"
"                       Example:  -rejectbatch id 10001-10004,10010,1002*\n" 
"                                                                        \n"
" -rejectbatch without options will clear all batch rejection settings.  \n"
"                       No batch rejection will be done.                 \n"
"                                                                        \n"
;

char* tScaleProcessOptions::ms_apcTableNames[] = { "rmergevsbatch","rmergevsintensity","rmergevsresolution","completeness","multiplicity","chisquareddistrib","anomsignal","inputstats","all",NULL};

int CDTMainScaleAverage::nExecute(unsigned int argc, char* argv[])    
{
#ifdef WIN32
  ios::sync_with_stdio(); 
#endif
  int nx,ny;
  float f0;
  Cstring sTemp;
    

    Cstring         sError("");
    // Structures used.
    Cstring         sArg;
    Cscaleaverage   *poScale=NULL;
    Creflnlist      *poReflnlist = NULL;

    // File output format variables
    Cstring         sOutFile;       // Output reflection file.
    Cstring         sUnavgOutFile;  // Output reflection file for unaveraged reflections.
    Cstring         sOutputHeader;  // Output header.
    Cstring         sOutputRejects; // Output rejects file.
    Cstring         sScaleOutFile;  // Scale reflection file.  Used for profile fitting based on R-merge with iterative pass.
    bool            bWriteHeader;
    bool            bWriteUnavgHeader;
    bool            bTexsan;
    bool            bTexsan2;
    bool            bShelx;
    bool            bAnom;
#ifdef SSI_PC
    bool            bTexsan2andShelx;
#endif
    bool            bScatterPlot = false;
    bool            bIgnoreReadError = FALSE;
    
    // General variables
    bool            bScaleAnom;     // Assume I+ != I-
    bool            bComputeOverlap;  // Count and print overlap information.
    Ccrystal*       poCrystal = NULL; // Unit cell and such.
    float           fResoMin;       // Upper resolution bound.
    float           fResoMax;       // Lower resolution bound.
    double          fDesiredChiSq = 1.0;   // See new code below for actual default used (bAnom affects it)
    bool            bChiSqTargetDefault = TRUE;  // If no -chisq option found in arg list, TRUE
    Cstring         sFixedBatch;    // Fixed batch for batchscaling.
    
    
    
    // Use specified activities.
    
    const int       nMAXPROCESSOPTIONS = 200;
    static tScaleProcessOptions aoOptions[nMAXPROCESSOPTIONS];
    int             nOptionCount;
    
    nModifyCommandLine(argc,&argv);
    vDtrekSetModuleName("dtscaleaverage");
    vPrintCopyrightInfo();
    //cout << "\ndtscaleaverage: Copyright (c) 2000-2006 Rigaku\n";
    //cout << D_K_DTREKVersion << endl << flush;

#ifdef SSI_PC
  Cstring sCCAppVersion = (const char*) argv[argc-1];
  if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
    printf(sCCAppVersion.string());
	argc -= 2;
	printf(CCrclHelper::GetInstance()->GetTimestamp().string());
	argv[argc] = NULL;
  }
#endif

    printf("Command line:\n %s\n\n", sGetCommandLine(argc, argv, 71).string());
    
    vSetCommandArguments(argc, argv);

    // Initialize variables.
    
    sOutFile            = "";
    sUnavgOutFile       = "";
    bTexsan = FALSE;
    bTexsan2 = FALSE;
    bShelx = FALSE;
    bAnom = FALSE;
#ifdef SSI_PC
    bTexsan2andShelx = FALSE;
#endif
    bWriteHeader=TRUE;
    bWriteUnavgHeader=TRUE;
    sOutputHeader  = sDtrekGetPrefix() + "dtscaleaverage.head";
    sOutputRejects = sDtrekGetPrefix() + "dtscaleaverage_rejects.ref";
    nOptionCount = 0;
    sFixedBatch = "";
    
    bScaleAnom          = FALSE;
    bComputeOverlap       = FALSE;
    fResoMin            = 999999.000f;
    fResoMax            = 0.00001f;
    
    // Clear all of the options (we might be in CC where the static array "aoOptions" is not reconstructed.
    for (nx=0;nx<nMAXPROCESSOPTIONS;nx++)
        aoOptions[nx].vReset();
    
    // RB Since the syntax of -reject option has changed, need to set default values for some aoOptions[a0] here:
    aoOptions[0].fReject             = 50.0;
    aoOptions[0].fFractionChiRejects   = 0.0;
    bool    bExcludeSysAbsences = false;

    bool    bBatchRestrainRequested = false;

    argc--; argv++; // Skip program name
    int nArg = 0;
    int nStat= 1;
    
    if( nArg >= argc )
    {
        DTREK_ERROR(0, pcHelp);
    }

    // Read in the options.
    printf("\n--------- Command line option parsing --------\n\n");
    while (nArg<argc) {
        
        sArg = argv[nArg];
        
        if (sArg.contains(".head"))
        {
            // Read in a header, and get relevant data.
            if (!bCreateHeader(sArg) ) 
            {
                printf("Header %s not available!\n",sArg.string());
                DTREKEXIT(1);
            };
            
            // Read the old crystal object and make new one.
            if (NULL != poCrystal)
                delete poCrystal;
            poCrystal = new Ccrystal(*m_poHeader);
        } 
        else if ((sArg.GetAt(0)!='-') || (sArg.contains(".ref")))
        {
            // Perhaps this is the output reflection file.
            
            
            if ((poReflnlist) && (nArg+1>=argc)) {
                printf("Output reflection file    %s\n", argv[nArg]);
                sOutFile = argv[nArg];
            } 
            else
            {                
                Ccrystal* poCrystalIn;
                Creflnlist* poNewList;
                
                // Read in reflection file.
                
                printf("Reading reflection list %s ... \n", sArg.string());
                
		// Look ahead in argument list to see if ignorereaderror is present
		int ia;
		for (ia = nArg; ia < argc; ia++)
		  {
		    sTemp = argv[ia];
		    if ("-ignorereaderror" == sTemp)
		      {
			bIgnoreReadError = TRUE;
			printf("Ignore any errors reading the input reflnlist.\n");
			break;
		      }
		  }
 
                poNewList = new Creflnlist(sArg, -1, 0, bIgnoreReadError);

                if (!poNewList->bIsAvailable())
                {
                    delete (poNewList);
                    delete (poReflnlist);
                    poReflnlist = NULL;
                    printf("Could not read reflnlist file: %s\n", sArg.string());
                    DTREKEXIT(1);
                }

                if (!poReflnlist) 
                    poReflnlist = poNewList;
                else
                {
                    poReflnlist->nInsertListFrom(*poNewList);
                    delete (poNewList);
                }           
                
                // Make sure oReflnlist will have all the fields needed
                // for scaling before we insert the other list.
                
                Cscaleaverage::vExpandReflnlist(poReflnlist);
                if (NULL == poCrystal)
                {
                    poCrystalIn = poReflnlist->poGetCrystal();
                    if (poCrystalIn && poCrystalIn->bIsAvailable()) 
                    {
                        poCrystal = poCrystalIn;  // No previous crystal, so use Reflnlist
                    }
                }
                
                // This seems to take a long time for large reflnlists
            };
            
        } 
        else if ( (sArg=="-help") || (sArg == "-h"))
        {
            nArg++;
        
            if ((nArg<argc) && ((sTemp=argv[nArg])=="sabsorb"))
            {
                DTREK_ERROR(0, pcSphericalAbsorbHelp);
            }
            else if ((nArg<argc) && ((sTemp=argv[nArg])=="errormodel"))
            {
                DTREK_ERROR(0, pcErrorModelHelp);
            }
            else if ((nArg<argc) && ((sTemp=argv[nArg])=="reject"))
            {
                DTREK_ERROR(0, pcRejectHelp);
            }
            else
            {
                DTREK_ERROR(0, pcHelp);
            }
        } 
        else if (sArg=="-ignorereaderror") 
	  {
	    // Ignore this since it was treated before
	    bIgnoreReadError = TRUE;
            printf("Ignore any errors reading the input reflnlist.\n");
	  }
        else if (sArg=="-spacegroup") 
        {
            nArg++;
            if (nArg>=argc) {
                printf("Bad syntax for option -spacegroup\n");
                DTREKEXIT(1);
            };
            if (NULL == poCrystal)
            {
                poCrystal = new Ccrystal();
                poCrystal->vSetCell(1., 1., 1., 90., 90., 90.);
            }
            if (strlen(argv[nArg]) == strspn(argv[nArg],"0123456789"))
                sscanf(argv[nArg],"%d",&nx);
            else {
                poCrystal->m_poSpacegroup->vSet(argv[nArg]);
                nx = poCrystal->m_poSpacegroup->nGet();
            };
            if ((nx>230) || (nx<1)) {
                printf("Invalid spacegroup specified\n");
                DTREKEXIT(1);
            };
            poCrystal->m_poSpacegroup->vSet(nx);
            printf("Spacegroup  %d (%s)\n",nx,poCrystal->m_poSpacegroup->sGetName().string());
        } 
        else if (sArg=="-cell") {
            float a6fCell[6];

            for (ny=0;ny<6;ny++) {
                nArg++;
                if ((nArg>=argc) || (1!=sscanf(argv[nArg],"%f",&a6fCell[ny]))) {
                    printf("Bad syntax for option -cell\n");
                    DTREKEXIT(1);
                };
                if (NULL == poCrystal)
                {
                    poCrystal = new Ccrystal();
                    poCrystal->m_poSpacegroup->vSet(1);
                }
                poCrystal->vSetCell(a6fCell);
                printf("Unit cell will be set to:    [%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f]\n",
                    a6fCell[0],a6fCell[1],a6fCell[2],a6fCell[3],a6fCell[4],a6fCell[5]);
            };
        } else if ( (sArg=="-ref") || (sArg=="-unavg") ) {
            nArg++;
            if (nArg>=argc) {
                printf("Bad syntax for option -ref\n");
                DTREKEXIT(1);
            };
            sUnavgOutFile=argv[nArg];
            printf("Unaveraged output file    %s\n",sUnavgOutFile.string());
        } else if (sArg=="-out") {
            nArg++;
            if (nArg>=argc) {
                printf("Bad syntax for option -out\n");
                DTREKEXIT(1);
            };
            sOutputHeader=argv[nArg];
            printf("Output header file             %s\n",sOutputHeader.string());
        } else if (sArg=="-scaleanom") {
            bScaleAnom=TRUE;

            printf("====================================================================INFO\n"
                   "INFO: -scaleanom used, is there enough redundancy?\n"
                   "Scaling set to assume  I+ != I-, that is: Friedel's Law does not apply.\n"
		   //+JWP 2008-05-23        
		   // Because of changes in Cscaleaverage.cpp remove the next 2 lines
		   //                   "Thus, the number of unique reflections with be roughly doubled,\n"
		   //                   "and the average redundancy will be roughly halved.\n"
		   //-JWP 2008-05-23        
                   "====================================================================INFO\n");
        } else if (sArg=="-texsan") {
            bTexsan=TRUE;
            printf("TEXSAN output selected.\n");
        } else if (sArg=="-texsan2") {
            bTexsan2=TRUE;
            printf("TEXSAN2 output selected.\n");
        } else if (sArg =="-shelx") {
            bShelx = TRUE;
            printf("SHELX output selected.\n");
#ifdef SSI_PC
        } else if (sArg =="-texsan2andshelx") {
            bTexsan2andShelx = TRUE;
            printf("Both TEXSAN2 and SHELX output selected.\n");
#endif
        } else if (sArg=="-anom") {
            bAnom=TRUE;
            printf("ANOM output set:               TRUE\n");
        } else if (sArg=="-noheader") {
            bWriteHeader=FALSE;
            printf("No header will be used in output reflection file.\n");
        } else if (sArg=="-nounavgheader") {
            bWriteUnavgHeader=FALSE;
            printf("No header will be used in unaveraged output reflection file.\n");
        } else if (sArg=="-batchscale") {
            aoOptions[nOptionCount].bBatchScale = TRUE;
            nOptionCount++;
            printf("Batch scaling enabled.\n");
        } else if (sArg=="-fix") {
            if (nArg+1<argc) {
                nArg++;
                sFixedBatch = argv[nArg];
                printf("Batch %s is fixed.\n",sFixedBatch.string());
            } else {
                printf("Bad syntax for option -fix\n");
            };
        } else if (sArg=="-bfactor") {
            aoOptions[nOptionCount].bBFactorScale = TRUE;
            nOptionCount++;
            printf("B factor scaling enabled.\n");
        } else if (sArg=="-list") {
            aoOptions[nOptionCount].bList = TRUE;
            if (nArg+1<argc) {
                nArg++;
                sArg=argv[nArg];
                sArg.downcase();
                for (nx=0;tScaleProcessOptions::ms_apcTableNames[nx];nx++) {
                    if (sArg==tScaleProcessOptions::ms_apcTableNames[nx])
                        break;
                };
                if (!tScaleProcessOptions::ms_apcTableNames[nx]) {
                    printf("Could not find table name '%s'\n",sArg.string());
                    DTREKEXIT(1);
                };
                aoOptions[nOptionCount].eTable = (_eTable) nx;
            } else {
                printf("Bad syntax for option -list\n");
            };
            nOptionCount++;
        } 
        else if (sArg=="-reso")
        {
            Cstring       sResoParameters("");
            nArg++;
            
            while( nArg < argc && ('-' != argv[nArg][0] && '+' != argv[nArg][0]) )
            {
                sResoParameters += argv[nArg];
                sResoParameters += " ";
                nArg++;
            }

            sError = "";
            if( !bParseCalculateResolution(sResoParameters, fResoMin, fResoMax, sError) )
            {
                DTREK_ERROR(6, sError);
            }
            nArg--;
        } 
        else if( sArg == "-reject" ) 
        { 
	        if ( ++nArg > argc )
	        {
		        printf("Bad syntax for option -reject\n");
		        DTREKEXIT(1);
	        }
            
            if( 0 == strcmp(argv[nArg], "deviation") ) 
	        {
		        nArg++;
		        if( nArg < argc  && 1 == sscanf(argv[nArg], "%f", &f0) ) 
                    aoOptions[0].fReject = f0;
		        else 
		        {
			        printf("Bad syntax for option -reject deviation\n");
                    DTREKEXIT(1);
                }
	        }
            else if( 0 == strcmp(argv[nArg], "fraction") ) 
	        {
		        nArg++;
		        if( nArg < argc  && 1 == sscanf(argv[nArg], "%f", &f0) ) 
                    aoOptions[0].fFractionChiRejects = f0;
		        else 
		        {
			        printf("Bad syntax for option -reject fraction\n");
                    DTREKEXIT(1);
                }
	        }
            else if( 0 == strcmp(argv[nArg], "sigma") ) 
	        {
		        nArg++;
		        if( nArg < argc  && 1 == sscanf(argv[nArg], "%f", &f0) ) 
                    aoOptions[0].fSigma = f0;
		        else 
		        {
			        printf("Bad syntax for option -reject sigma\n");
                    DTREKEXIT(1);
                }
	        }
            else if( 0 == strcmp(argv[nArg], "symabs") ) 
		        bExcludeSysAbsences = true;
		    else if( 0 == strcmp(argv[nArg], "list") )
            {
                if( ++nArg < argc )
                    aoOptions[0].sPrevRejectFile = argv[nArg];
                else
                {
                    printf("Bad syntax for option -reject list\n");
                    DTREKEXIT(1);
                }
            }
		    else if ( 1 == sscanf(argv[nArg], "%f", &f0) ) // support for old syntax like -reject 0.0005
            {
                  aoOptions[0].fFractionChiRejects = f0;
            }
            else 
		    {
			    printf("Bad syntax for option -reject\n");
                DTREKEXIT(1);
            }
        }
        else if (sArg == "-scalerefs") {
            if (nArg+1 < argc) {
                sScaleOutFile = argv[nArg+1];
                nArg++;
            } else {
                printf("Bad syntax for option -scalerefs\n");
                DTREKEXIT(1);

            };
        } else if (sArg.contains("-rejectbatchid=")) { 
            aoOptions[0].sBatchRejectTemplate = sArg.after("-rejectbatchid=");
            aoOptions[0].bBatchRejectFlags[BATCH_BATCHNAME] = TRUE;
        } else if (sArg=="-rejectbatch") { 
            if (nArg+1<argc) {
                nArg++;
                if ((sTemp=argv[nArg])=="rmerge") {
                    nArg++;
                    if ((nArg<argc) && (1 == sscanf(argv[nArg],"%f",&f0))) {
                        aoOptions[0].fBatchRejectTypes[BATCH_RMERGE] = f0;
                        aoOptions[0].bBatchRejectFlags[BATCH_RMERGE] = TRUE;
                    } else {
                        printf("Bad syntax for option -rejectbatch\n");
                        DTREKEXIT(1);
                    };
                } else if ((sTemp=argv[nArg]) == "intensity") {
                    nArg++;
                    if ((nArg<argc) && (1 == sscanf(argv[nArg],"%f",&f0))) {
                        aoOptions[0].fBatchRejectTypes[BATCH_INTENSITY] = f0;
                        aoOptions[0].bBatchRejectFlags[BATCH_INTENSITY] = TRUE;
                    } else {
                        printf("Bad syntax for option -rejectbatch\n");
                        DTREKEXIT(1);
                    };
                } else if ((sTemp=argv[nArg]) == "percentbad") {
                    nArg++;
                    if ((nArg<argc) && (1 == sscanf(argv[nArg],"%f",&f0))) {
                        aoOptions[0].fBatchRejectTypes[BATCH_PERCENT] = f0;
                        aoOptions[0].bBatchRejectFlags[BATCH_PERCENT] = TRUE;
                    } else {
                        printf("Bad syntax for option -rejectbatch\n");
                        DTREKEXIT(1);
                    };
                } else if ((sTemp=argv[nArg]) == "chisq") {
                    nArg++;
                    if ((nArg<argc) && (1 == sscanf(argv[nArg],"%f",&f0))) {
                        aoOptions[0].fBatchRejectTypes[BATCH_CHISQ] = f0;
                        aoOptions[0].bBatchRejectFlags[BATCH_CHISQ] = TRUE;
                    } else {
                        printf("Bad syntax for option -rejectbatch\n");
                        DTREKEXIT(1);
                    };
                } else if ((sTemp = argv[nArg]) == "id") {
                    nArg++;
                    if (nArg<argc) {
                        aoOptions[0].sBatchRejectTemplate = argv[nArg];
                        aoOptions[0].bBatchRejectFlags[BATCH_BATCHNAME] = TRUE;
                    } else {
                        printf("Bad syntax for option -rejectbatch\n");
                        DTREKEXIT(1);
                    };
                } else
                    nArg--;
            };                         
        } else if (sArg=="-sigma") // support for old syntax "-sigma" 
        {
            nArg++;
            if ((nArg>=argc) || (1!=sscanf(argv[nArg],"%f",&f0))) {
                printf("Bad syntax for option -sigma\n");
                DTREKEXIT(1);
            };
            aoOptions[0].fSigma = f0;
        } else if (sArg=="-nochisq") {
            fDesiredChiSq = (double) -1.0;
        } else if (sArg== "-scatterplot") {
            bScatterPlot = true;
        } else if (sArg=="-chisq") {
            nArg++;
            if ((nArg>=argc) || (1!=sscanf(argv[nArg],"%f",&f0))) {
                printf("Bad syntax for option -chisq\n");
                DTREKEXIT(1);
            };
            fDesiredChiSq = (double) f0;
	    bChiSqTargetDefault = FALSE;
        } else if (sArg=="-countoverlap") {
            bComputeOverlap=TRUE;
            printf("Overlap counting ENABLED.\n");
        } 
        else if( sArg == "-sabsorb" )
        {            
            aoOptions[nOptionCount].bSphericalAbsorb = TRUE;
            
            if( (nArg + 1 < argc) && (1 == sscanf(argv[nArg+1],"%d",&nx)) )
            {
                aoOptions[nOptionCount].nHarmonicOrderDifrac = nx;
                nArg++;
            }
        
            if( (nArg + 1 < argc) && (1 == sscanf(argv[nArg+1],"%d",&nx)) )
            {
                aoOptions[nOptionCount].nHarmonicOrderIncident = nx;
                nArg++;
            }
            
            if( (nArg + 1 < argc) && (1 == sscanf(argv[nArg+1],"%d",&nx)) )
            {
                aoOptions[nOptionCount].nMaxIter = nx;
                nArg++;
            }
            
            nOptionCount++;
        } else if (sArg=="-errormodel") {
            while (1) {
                if ((nArg+2<argc) && (argv[nArg+1]==(Cstring) "add") && sscanf(argv[nArg+2],"%f",&f0)) {
                    if (f0>=0.0) {
                        aoOptions[0].fAdd = f0;
			            aoOptions[0].fAddMult = -1.0;
		            };
                    nArg+=2;
                } else if ((nArg+2<argc) && (argv[nArg+1]==(Cstring) "mul") && sscanf(argv[nArg+2],"%f",&f0)) {
                    if (f0>0.0) {
                        aoOptions[0].fMul = f0;
					};
                    nArg+=2;
				} else if ((nArg+2<argc) && (argv[nArg+1]==(Cstring) "addmult") && (sscanf(argv[nArg+2],"%f",&f0))) {
					if (f0>0.0) {
						aoOptions[0].fAddMult = f0;
						aoOptions[0].fAdd = -1.0;
					};
					nArg+=2;
                } else
                    break;
            };           
            if ((nArg+1<argc) && (1==sscanf(argv[nArg+1],"%f",&f0))) {
                nArg++;
                aoOptions[0].fMul = f0;
                if ((nArg+1<argc) && (1==sscanf(argv[nArg+1],"%f",&f0))) {
                    nArg++;
                    aoOptions[0].fAdd = f0;
					aoOptions[0].fAddMult = -1.0;
                };
            };
        } else if (sArg=="-rebatch") {
            nArg++;
            if ((nArg<argc) && ((sTemp = argv[nArg])== "imageboundary")) {
                nArg++;
                aoOptions[nOptionCount].bRebatchImageBoundary = true;
            } else
                aoOptions[nOptionCount].bRebatchImageBoundary = false;

            if ((nArg<argc) && ((sTemp = argv[nArg])== "bookends")) 
                aoOptions[nOptionCount].bRebatchBookends = true;
            else
                aoOptions[nOptionCount].bRebatchBookends = false;

            if (aoOptions[nOptionCount].bRebatchBookends && aoOptions[nOptionCount].bRebatchImageBoundary) {
                printf("WARNING: bookends and imageboundary are mutually exclusive!\n");
                aoOptions[nOptionCount].bRebatchImageBoundary = false;
            };
            if (!aoOptions[nOptionCount].bRebatchBookends) {
                if ((nArg>=argc) || (1!=sscanf(argv[nArg],"%d",&aoOptions[nOptionCount].nMinOverlaps))) {
                    printf("Incorrect parameters for -rebatch.\n");
                    nStat = 1;
                    DTREKEXIT(nStat);
                }
            } else
                aoOptions[nOptionCount].nMinOverlaps = 0;
            aoOptions[nOptionCount].bRebatch = TRUE;
            nOptionCount++;
        } 
        else if (sArg=="-completenessreject")
        {
            nArg++;
            if ((nArg>=argc) || (1!=sscanf(argv[nArg],"%f",&f0)))
            {
                printf("Incorrect command line for %s\n",sArg.string()); 
                DTREKEXIT(1);
            }
            nArg++;
            nx  = 1;
            if ((nArg<argc) && ((sTemp=argv[nArg])=="+"))
            {
                nx = 1;
                nArg++;
            } 
            else if ((nArg<argc) && ((sTemp=argv[nArg])=="-"))
            {
                nx = -1;
                nArg++;
            }
            aoOptions[nOptionCount].nSign = nx;
            aoOptions[nOptionCount].fCompleteness = f0;

            if ((nArg+1<argc) && (!(sTemp=argv[nArg+1]).contains('-')) && (!(sTemp=argv[nArg]).contains('-')))
            {
                // Two Fields.
                aoOptions[nOptionCount].sField1 = argv[nArg];
                aoOptions[nOptionCount].sField2 = argv[nArg+1];
                nArg+=1;
            } 
            else if ((nArg<argc) && (!(sTemp=argv[nArg]).contains('-')))
            {
                aoOptions[nOptionCount].sField1 = argv[nArg];
                aoOptions[nOptionCount].sField2 = "";
            }
            else
            {
                printf("Incorrect command line for %s\n",sArg.string()); 
                DTREKEXIT(1);
            }
            
            aoOptions[nOptionCount].bCompletenessReject = TRUE;        
        } 
        else if( sArg=="-batchrestrain" )
        {
            nArg++;
            
            Cstring         strTemp("");
            
            // If -batchrestrain is the LAST command line argument, strTemp should be left empty.
            // If there is another argument after -batchrestrain  and it is NOT an option keyword
            // grab it into strTemp
            if( nArg < argc )
            {
                if( '-' != argv[nArg][0] ) // if next token is NOT a new option...
                    strTemp = argv[nArg];
                else
                    nArg--; // no restrain option supplied, so no need to increment nArg
            }
            else if( nArg > argc )
            {
                printf("Incorrect command line for %s\n", sArg.string()); 
                DTREKEXIT(1);
            }

            bBatchRestrainRequested = true;
            aoOptions[0].sBatchRestrInfo = strTemp;
        } 
        else
        {
            printf("Unrecognized option %s\n",sArg.string());
            nStat = 1;
            DTREKEXIT(nStat);
        }        
        nArg++;
    }
    // Add default options here.
    int     nOption = 0;
    for (nOption=0;nOption<nOptionCount;nOption++)
    {
        if (aoOptions[nOption].bList)
            break;
    }
    
    if( nOption == nOptionCount )
    {
        aoOptions[nOptionCount].bList = TRUE;
        aoOptions[nOptionCount].eTable = TABLE_ALL;
        nOptionCount++;
    }
    
    printf("----- End of command line options -----\n");
    
    if ( (NULL == poCrystal) || !poCrystal->bIsAvailable() ) {
        printf("ERROR, no crystal info available!\n");
        DTREKEXIT(1);
    }
    poCrystal->nList(1);
    if (NULL==poReflnlist) {
        printf("ERROR:  No reflection list specified!\n");
        DTREKEXIT(1);
    };

    // Delete any -rejectbatch id batches from the input reflnlist

    if (aoOptions[0].bBatchRejectFlags[BATCH_BATCHNAME])
      {
        printf("\nINFO: Rejecting (deleting) reflns in batches matching template:\n      %s\n\n", aoOptions[0].sBatchRejectTemplate.string());
        poReflnlist->nDeleteBatchID(aoOptions[0].sBatchRejectTemplate);
        aoOptions[0].bBatchRejectFlags[BATCH_BATCHNAME] = FALSE;
      }
    

    Creflnlist& oReflnlist = *poReflnlist;

    // Make sure that we read in a header file, and at least one reflection file.

    if (!oReflnlist.nGetNumReflns()) {
        printf("No Reflections found.\n");
        goto main_termination;
    };
    if( !m_poHeader )
    {
        // If we don't have a header, it's okay... just build a dummy.
        vCreateHeader();
    };
    // Place the crystal object inside the header and into the input
    // reflnlist.
    // Note: if the -spacegroup option was used, then the unit cell 
    //       parameters may be inconsistent with the spacegroup.

    poCrystal->nUpdateHeader(m_poHeader);
    oReflnlist.nPutCrystal(*poCrystal);

    // Build the scaleaverage structure.

    poScale = new Cscaleaverage(oReflnlist,*m_poHeader);

    poScale->m_fResolutionMin = max(fResoMin,fResoMax);
    poScale->m_fResolutionMax = min(fResoMin,fResoMax);
    poScale->m_fReject = aoOptions[0].fReject; 
    poScale->m_fMaxDevReject = aoOptions[0].fReject; 
    poScale->m_fSigma = aoOptions[0].fSigma;
    poScale->m_fFractionReject = aoOptions[0].fFractionChiRejects;

//+2011-07-17 JWP
    if (bChiSqTargetDefault)
      {
//  If no -chisq option was used, set the default target to 1.0 unless the -anom option is used.
//  If -anom is used, then the default target is 1.2.
	fDesiredChiSq = 1.0;
	if (bAnom) fDesiredChiSq = 1.2;
      }
//-2011-07-17 JWP
    poScale->m_fTargetChiSquare = fDesiredChiSq;

    poScale->m_fBatchPercentReject = 1000.0;
    poScale->m_fBatchRmergeReject = 1000.0;
    poScale->m_fBatchIoverSigReject = 1000.0;    
    poScale->m_fBatchChiSqReject = 1000.0;
    poScale->m_sBatchRejectTemplate = "";
    poScale->m_fUserEAdd = aoOptions[0].fAdd;
    poScale->m_fUserEMul = aoOptions[0].fMul;
    poScale->m_fUserEAddMult = aoOptions[0].fAddMult;
    poScale->m_bExcludeSysAbsences = bExcludeSysAbsences;

    //////////////////////////////////////////////////////////////////
    // RB Add restrain information
    if( bBatchRestrainRequested )
        poScale->vSetScaleRestrainInfo(aoOptions[0].sBatchRestrInfo);
    //////////////////////////////////////////////////////////////////

    if (aoOptions[0].sPrevRejectFile.length())
      {
	printf("\n\nINFO: -reject list option used:\n");
        poScale->m_poPrevRejectList = new Creflnlist(aoOptions[0].sPrevRejectFile);
      }

    if (aoOptions[0].bBatchRejectFlags[BATCH_PERCENT]) {
        printf("\nRejecting batches with more than %5.3f%% rejected reflections.\n",100.0*aoOptions[0].fBatchRejectTypes[BATCH_PERCENT]);
        poScale->m_fBatchPercentReject = aoOptions[0].fBatchRejectTypes[BATCH_PERCENT];        
        poScale->m_bRejectBatches = TRUE;
    } 
    if (aoOptions[0].bBatchRejectFlags[BATCH_RMERGE]) {
        if (aoOptions[0].fBatchRejectTypes[BATCH_RMERGE]>=2.0)
            printf("\nRejecting batches at end with\n   Rmerge > %5.2f sigma above average batch Rmerge.\n",aoOptions[0].fBatchRejectTypes[BATCH_RMERGE]);
        else if (aoOptions[0].fBatchRejectTypes[BATCH_RMERGE]< 1.0)
            printf("\nRejecting batches at end with\n   Rmerge > %5.2f%%\n",100*aoOptions[0].fBatchRejectTypes[BATCH_RMERGE]);
	else
            printf("\nWARNING: Illegal value for batch Rmerge reject option: %5.2f\n",aoOptions[0].fBatchRejectTypes[BATCH_RMERGE]);
        poScale->m_fBatchRmergeReject = aoOptions[0].fBatchRejectTypes[BATCH_RMERGE];
        poScale->m_bRejectBatches = TRUE;
    } 
    if (aoOptions[0].bBatchRejectFlags[BATCH_INTENSITY]) {
        printf("\nRejecting batches with <I/sig> outside %5.2f deviations of <<I/sig>>.\n",aoOptions[0].fBatchRejectTypes[BATCH_INTENSITY]);
        poScale->m_fBatchIoverSigReject = aoOptions[0].fBatchRejectTypes[BATCH_INTENSITY];
        poScale->m_bRejectBatches = TRUE;
    } 
    if (aoOptions[0].bBatchRejectFlags[BATCH_CHISQ]) {
        poScale->m_fBatchChiSqReject = aoOptions[0].fBatchRejectTypes[BATCH_CHISQ];
        poScale->m_bRejectBatches = TRUE;
	if (0.0 > poScale->m_fBatchChiSqReject)
	  printf("\nRejecting batches at end with\n   batch |ChiSq| > %.2f \n",
		 ABS(aoOptions[0].fBatchRejectTypes[BATCH_CHISQ]));

	else if (2.0 <= poScale->m_fBatchChiSqReject)
	  printf("\nRejecting batches at end with\n   batch |ChiSq| > %5.2f sigma above overall average batch |ChiSq|.\n",aoOptions[0].fBatchRejectTypes[BATCH_CHISQ]);
	else
            printf("\nWARNING: Illegal value for batch |ChiSq| reject option: %5.2f\n",aoOptions[0].fBatchRejectTypes[BATCH_CHISQ]);
    } 
    if (aoOptions[0].bBatchRejectFlags[BATCH_BATCHNAME]) {
        printf("\nRejecting batches matching template '%s'.\n",aoOptions[0].sBatchRejectTemplate.string());
        poScale->m_sBatchRejectTemplate = aoOptions[0].sBatchRejectTemplate;

    };

    printf("Rejection sigma set to       %6.2f\n",poScale->m_fReject);
    if (poScale->m_fFractionReject>0.0) {
        printf("Maximum |ChiSq| selected so that %6.2f%% of data is rejected\n",poScale->m_fFractionReject*100.0);
    };
    printf("Reflections excluded from scale factor calculation (NOT results and statistics) with I/SigmaI less than %6.2f\n",poScale->m_fSigma);
    if ((poScale->m_fUserEMul>0.0) && (poScale->m_fUserEAdd>=0.0)) {
        printf("Using Explicit Emul = %6.2f Explicit Eadd = %6.2f.\n",poScale->m_fUserEMul,poScale->m_fUserEAdd);
    } else if (poScale->m_fUserEMul>0.0) {
        printf("Using Explicit Emul = %6.2f Auto Eadd.\n",poScale->m_fUserEMul);
    } else {
        printf("Using Auto Emul Auto Eadd.\n");
    };
        
    poScale->m_bScaleAnom=bScaleAnom;

    // nInitArrays() MUST BE CALLED BEFORE OTHER ROUTINES.  IT MAKES IMPORTANT CHANGES TO REFLECTION FILE
    // ORDERING AND DELETIONS.
    if( 0 != poScale->nInitArrays() )
        return 1;

    // Turn on overlap computation just once.
    poScale->m_bComputeOverlap = TRUE;
    poScale->fCalcChiSquared(1);
    poScale->m_bComputeOverlap = bComputeOverlap;
    poScale->nFitChiSquared();

    for (nOption=0;nOption<nOptionCount;nOption++)
    {
        tScaleProcessOptions& oOption = aoOptions[nOption];

        if (oOption.bBatchScale) { 
            printf("Calling Batch scale...\n");
            fflush(stdout);
            poScale->nExcludePriorReflections();            
            if (poScale->nBatchScale(sFixedBatch,TRUE)) {
                printf("\nWARNING:  Problem with batch scaling.\n");
            };
        };
        if (oOption.bBFactorScale) {
            printf("Calling B factor scale...\n");
            poScale->nExcludePriorReflections();
            if (poScale->nBFactorScale(sFixedBatch,TRUE)) {
                printf("\nWARNING:  Problem with B factor scaling.\n");
            };
            poScale->nFitChiSquared();
            poScale->fCalcChiSquared(1);

        };

        if (oOption.bSphericalAbsorb)
        {
            printf("Calling Spherical harmonic scale...\n");
            fflush(stdout);
            poScale->nExcludePriorReflections();            
            if (poScale->nSphericalHarmonics(oOption)) {
                printf("WARNING: Problem with spherical harmonic scaling.\n");
            };
            poScale->nFitChiSquared();
            poScale->fCalcChiSquared(1);
        };
       
/*+JWP 2007-12-15
        if (oOption.bList) {
            switch (oOption.eTable) {
            case TABLE_ANOMSIGNAL:
	      // (disabled for now) poScale->nTableAnomSignal(); 
	      break;
            case TABLE_RMERGEVSBATCH:       
                poScale->nTableRMergeVsBatch(); break;
            case TABLE_RMERGEVSINTENSITY:   
                poScale->nTableRMergeVsIntensity(); break;
            case TABLE_RMERGEVSRESOLUTION:  
                poScale->nTableRMergeVsResolution(); break;
            case TABLE_COMPLETENESS:        
                poScale->nTableCompleteness(); break;
            case TABLE_MULTIPLICITY:
                poScale->nTableMultiplicity(); break;
            case TABLE_CHIDISTRIB:
                poScale->nTableChiDistribution(); break;
            case TABLE_INPUTSTATS:
                poScale->nTableIntensityResolutionPosition(); break;
            case TABLE_ALL:
                // (disabled for now) poScale->nTableAnomSignal();
                
                // One last chi-squared so that rejects are correct.

                poScale->fCalcChiSquared(2);

                poScale->nTableMultiplicity();
                poScale->nTableChiDistribution();
                poScale->nTableRMergeVsBatch();
                poScale->nTableRMergeVsIntensity();
                poScale->nTableCompleteness();
                poScale->nTableRMergeVsResolution();
                break;
            };
        };
-JWP 12-15-2007 */
        if (oOption.bRebatch) {
            nStat = poScale->nReBatch(oOption.nMinOverlaps,oOption.bRebatchImageBoundary,oOption.bRebatchBookends,TRUE);
            if (nStat == 2) {
                printf("Error in Rebatch()!\n");
                goto main_termination;
            } else if (nStat == 1) {
                poScale->nInitArrays();  
                poScale->m_bComputeOverlap = TRUE;
                poScale->fCalcChiSquared(1);
                poScale->m_bComputeOverlap = bComputeOverlap;
                poScale->nReBatch(oOption.nMinOverlaps,oOption.bRebatchImageBoundary,oOption.bRebatchBookends,FALSE);
            };
        };

        if (oOption.bCompletenessReject) {
            poScale->nCompletenessCutoff(oOption.fCompleteness,oOption.sField1,oOption.sField2,oOption.nSign==-1);
        };

    };
    if (bScatterPlot)
        poScale->nChiPlot();

    // Write out the reflection lists.
#ifdef SSI_PC
    poScale->nWriteReflns(sOutFile,sUnavgOutFile,sOutputRejects,sScaleOutFile,bWriteHeader,bWriteUnavgHeader,bTexsan,bTexsan2,bShelx,bAnom,bTexsan2andShelx);
#else
    poScale->nWriteReflns(sOutFile,sUnavgOutFile,sOutputRejects,sScaleOutFile,bWriteHeader,bWriteUnavgHeader,bTexsan,bTexsan2,bShelx,bAnom);
#endif
    // Write out header.
        m_poHeader->nWrite(sOutputHeader);

//+jwp  Moved ALL table printing AFTER writing out the reflnlist and header
//      Note that the anomalous scattering signal analysis printout MUST OCCUR AFTER writing the reflnlist.

        if (TRUE) {
            switch (TABLE_ALL) {
            case TABLE_ANOMSIGNAL:
	      // (disabled for now) poScale->nTableAnomSignal();
	      break;
            case TABLE_RMERGEVSBATCH:       
                poScale->nTableRMergeVsBatch(); break;
            case TABLE_RMERGEVSINTENSITY:   
                poScale->nTableRMergeVsIntensity(); break;
            case TABLE_RMERGEVSRESOLUTION:  
                poScale->nTableRMergeVsResolution(); break;
            case TABLE_COMPLETENESS:        
                poScale->nTableCompleteness(); break;
            case TABLE_MULTIPLICITY:
                poScale->nTableMultiplicity(); break;
            case TABLE_CHIDISTRIB:
                poScale->nTableChiDistribution(); break;
            case TABLE_INPUTSTATS:
                poScale->nTableIntensityResolutionPosition(); break;
            case TABLE_ALL:
	      // (disabled for now) poScale->nTableAnomSignal();
                
	      //+JWP 2008-06-07 Remove per Ken Yates
                // One last chi-squared so that rejects are correct.

	      //poScale->fCalcChiSquared(2);
	      //-JWP 2008-06-07

                poScale->nTableMultiplicity();
                poScale->nTableChiDistribution();
                poScale->nTableRMergeVsBatch();
                poScale->nTableRrimVsBatch();
                //+-JWP 2008-05-23 We do not need this table: 
		//poScale->nTableRpimVsBatch();
                //-JWP 2008-05-23 
                poScale->nTableRMergeVsIntensity();
                poScale->nTableRrimVsIntensity();
		//+-JWP 2008-05-23 We do not need this table: 
                //poScale->nTableRpimVsIntensity();
                //-JWP 2008-05-23 
		        if (bAnom)
                {
                    poScale->m_oRASBins.vPrintStats();   // MUST happen AFTER poScale->nWriteReflns()
                }
                poScale->nTableCompleteness();
                poScale->nTableRMergeVsResolution(false);   //do not print summary
                poScale->nTableRrimVsResolution(false);     //do not print summary
		//+-JWP 2008-05-23 We do not need this table: 
                //poScale->nTableRpimVsResolution(false);     //do not print summary
                //-JWP 2008-05-23 
                poScale->nTable3RVsResolution();            //print summary only once here

		        poScale->nReportRejects();
                break;
            };
        };
//-jwp

    nStat = 0;
main_termination:
    
    if (poScale)
        delete poScale;
    delete poCrystal;
    
    if( poReflnlist )
    {
        delete poReflnlist;
        poReflnlist = NULL;
    }

    printf ("\ndtscaleaverage: Done.\n");

    printf ("\nINFO - The suggested next step is to use jdtplot to plot the"
            "\n       results, then perhaps reject outlier or problem batches"
            "\n       in a subsequent run of dtscaleaverage.  If successful,"
            "\n       the next steps are to use the %s file with subsequent"
            "\n       crystallographic programs.\n\n", sOutFile.string());
    printf ("\nREFERENCE - Pflugrath, JW (1999) Acta Cryst. D55, 1718-1725."
            "\n       Although that is the official publication, Thad Niemeyer"
            "\n       of Rigaku/MSC has re-written the bulk of the crystal-"
            "\n       lographic code in d*TREK.  Robert Bolotovsky of Rigaku"
            "\n       has made numerous contributions as well."
            "\n       We also want to acknowledge the works of many others.\n\n");
    fflush(stdout);

    return nStat;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainScaleAverage::vError(const int nErrorNum, const Cstring& sMessage)
{
    if( 0 != nErrorNum ) 
    {
        Cstring     sErrorMsg = "\nERROR: ";
        sErrorMsg += sMessage;
        sErrorMsg += "\n";

        printf("%s\n",sErrorMsg.string());
    }
    else
        printf("%s\n",sMessage.string());
}

