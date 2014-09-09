#ifndef DT_CSCALEAVERAGE_H
#define DT_CSCALEAVERAGE_H
//
// Copyright (c) 2000 Molecular Structure Corporation
//                    9009 New Trails Drive
//                    The Woodlands, TX, USA  77381
//
// The contents are unpublished proprietary source
// code of Molecular Structure Corporation
//
// All rights reserved.
//
// Cscaleaverage.h     Initial author: T.J. Niemeyer  Spring 2000
//             Based on Cscalemerge.h
//    Header file for Cscaleaverage.cc
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

#include <string.h>
#include "Dtrek.h"

#include "Creflnlist.h"
#include "Cscan.h"
#include "CinterpMesh.h"


const int g_nErrorModelParams = 4;
const int g_nErrorModelBins = 20; // Must be a multiple of the number of resolution bins.

namespace hac
{
class DTREK_abscor_ref_item;
class ReflnMemoryPool; 
}

struct tagStatsScaleAverage
{
  int   nNumRefs;       // Number of reflns
  int   nNumRejects;    // Number of rejected reflections.
  int   nNumUsed;       // Number of NON rejected reflections that are not singles.
  int   nNumUsed2;      // Number of NON rejected reflections that participate in the same group.
  int   nNumExcluded;   // Number of excluded reflections.
  int   nNumMults;      // Number of sym-equiv groups.  Used in Chi-Squared.
  int   nNumSingles;    // Number of singles (nNumSingles+nNumUsed = nNumRefs)

  double fNumer;         // Numerator of R-merge (Only over mult.>=2)

  //shijie yao -- 04:01:08 added members for additional merging R factors
  double fNumer_rim;     // Numerator of reduncance-independent merging R factor (Only over mult.>=2)
  double fNumer_pim;     // Numerator of precision-indicating  merging R factor (Only over mult.>=2)
  double fNumer_pcv;     // Numerator of pooled coefficient of variation (Only over mult.>=2)

  //jwp - 2008-06-10 add members for separate R factors with I+ = I- and I+ != I-
  double fNumer_rimplus;
  double fNumer_rimminus;
  double fNumer_plus;
  double fNumer_minus;

  double fDenom;         // Denominator of R-merge. (Only over mult.>=2)
  double fDenomAnom;     // Denominator of R's, but used when keeping I+ and I- separate
  double fSum;           // Same as fDenom, except includes multiplicites <2. Used to calculate the average.
  double fSumIoverSig;   // Sum of I/Sigma(I)
  double fSumIoverSig2;  // Sum of I/Sigma(I) (Averaged Reflection I/sig)
  double fSumChiSq;      // Chi-squared sum.  Must divide by (nSumChiContrib - nSumChiGroups)
  double fSumChiSqA;      // Chi-squared sum with I+ and I- kept separate.  Must divide by (nSumChiContrib - nSumChiGroups)
  double fEAddNumer;     // Eadd for each chi-squared group.  Must divide by fEaddDenom to get average Eadd)
  double fEAddDenom;
  int    nSumChiGroups;
  int    nSumChiContrib;
  int    nSumChiAGroups;
  int    nSumChiAContrib;


  int   nLastRefFlag;   // Used for computing nNumMults.  We store a cookie here to see when the
                        // structure was last modified.

  int anNumInBinMults[10];   // Multiplicities of reflections (see redundancy vs resolution table for usage)

  void vInitStats();
  tagStatsScaleAverage& operator += (tagStatsScaleAverage& oTag);
};
struct tScaleProcessOptions;
class  Cabsorbmodel;

struct tagBatch {
    Cstring m_sBatchName;            // Batch names in reflecion file.
    bool    m_bBatchReject;          // Rejected batches.  (superset of RejectsUser) Filled in when fCalcChiSquared() is called.
    bool    m_bBatchRejectUser;      // User defined flags rejecting certain batches.
    float   m_fBatchScale;           // Batch scale factors. (1.0 by default)
    float   m_fBFactorScale;         // B factor scales.     (0.0 by default)
    
    
    float   m_fBFactorConst;         // exponent argument scale Factor, refined simultaneously with B-factor (default 0, i.e. scale=e**0=1)
    
    int     m_nNumReflns;            // Number of reflections in batch.
    
    Cstring m_sScan;
    int     m_nBatch;
};

typedef struct _tagSCALE_RESTRAIN_INFO
{
    Cstring        m_sBatchBegin;
    Cstring        m_sBatchEnd;
    double         m_dWeight;
    
    _tagSCALE_RESTRAIN_INFO(Cstring& s1, Cstring& s2, double dWeight){m_sBatchBegin=s1; m_sBatchEnd=s2; m_dWeight=dWeight;}
}SCALE_RESTRAIN_INFO;


class DTREK_EXPORT  Cscaleaverage
{

public:
    Creflnlist& m_oReflnlist;
    Ccrystal m_oCrystal;
    Cscan*   m_poScan;
    float    m_a3x3fBMat[3][3];      // Var for crystal B matrix
    
    float*   m_pfAbsorbFactors;      // Used to apply/unapply a set of scale factors.
    float*   m_pfAbsorbFactorsSig;   // Used to apply/unapply a set of scale factors.

    int m_nStat;
    int m_nNumRefs;
    int m_nTotalAllRejects;         // ALL rejects (deviation rejects)
    int m_nTotalSigmaExcluded;      // ALL sigma excluded.
    int m_nTotalAdjustedChiRejects;


    // Reflection list indices.
    int m_nBatchI;
    int m_nFI_nIntensityI;
    int m_nFI_nSigmaI;
    int m_nFI_nBackgroundSigma;
    int m_nFI_fIntensityPlus;
    int m_nFI_fSigmaIPlus;
    int m_nFI_fIntensityMinus;
    int m_nFI_fSigmaIMinus;

    // Information arrays for batches.
    tagBatch* m_ptBatchInfo;            // Batch Info.
    float*  m_pfMin[5];
    float*  m_pfMax[5];                 // Min, max of 5 batch variables

    // Information arrays.

    float* m_pfMergeSigma;              // The sigma calculated for the merge of each reflection with others in it's group.
    float* m_pfChi;                     // The average intensity for each reflection.
    float* m_pfDev;
    float* m_pfDev2;
    float* m_pfOrigIntensity;           // Array of original (input) intensities.
    int*   m_pnGroup;                   // Sym group.
    int*   m_pnGroupSize;               // Size of group.

    // More Information arrays
    int* m_pnNumOverlap;                // Num-batches x Num-batches intersections.
    int* m_pnIndex;                     // Index for HKL ordering of reflections.
    int* m_pnRevIndex;                  // Sorted HKL index for each reflection.
    int* m_pnIndexBatch;                // Index for Batch ordering of reflections.

    int* m_pnIndexIntensity;            // highest to lowest intensity.
    
    float* m_pfRevIndexReso;            // Resolution rank of each reflection.
    int* m_pnIndexReso;                 // Sort reflections on increasing resolution (floating point gives bins that point is in).
    int* m_pnReject;                    // Reject flags. (WRT sorted index... i.e., Do not dereference with m_pnIndex)
    
    // Old information arrays (now used to plot chisq vs redundancy)

    int  m_a10nMult[10];                // Multiplicities for reflections.
    int  m_a10nMultLE[10];              // Multiplicities <= 
    int  m_a10nUsedLE[10];              // Number used <= 
    float m_a10fSumChi[10];             // chi-squared for reflections.

    // Arrays used to compute chi-squared distribution.    
    int    m_nChiBins;                  // Number of distribution bins.
    float* m_pfChiDistribBins;          // Chi distribution bins.  These are allocated once
    int*   m_pnChiDistribGroups;        // How many groups fall into each bin.
    float* m_pfChiDistribMultiplicity;  // Multiplicity of each bin.

    tagStatsScaleAverage* m_ptStatsReso;
    tagStatsScaleAverage* m_ptStatsInt;
    tagStatsScaleAverage* m_ptStatsBatch;

    
    int m_nFI_f2STLsq;
    int m_nFI_nBatchIdx;
    float m_a2fRangeReso[2];            // Observed maximum and minimum resolution.
    float m_a2fRangeInt[2];             // Observed maximum and minimum intensities.
    float m_fOverallPercentComplete;    // Overall completeness (loaded after calling completeness table).
    float m_fSlopeReso;
    float m_fSlopeInt;
    int m_nNumBinsReso;                 // Number of resolution bins. (Fixed)
    int m_nNumBinsInt;                  // Number of intensity bins.  (Fixed)

    float   m_fWavelength;
    int     m_nNumUniqueFoundCumul;
    int     m_nNumUniqueCalcCumul;
    int     m_nNumUniqueFoundShell;
    int     m_nNumUniqueCalcShell;

    friend int nCscaleaverage_UserFunction1();            // External functions that can access Cscaleaverage variables (used for testing).
    friend int nCscaleaverage_UserFunction2();            // External functions that can access Cscaleaverage variables (used for testing).

    double* m_pfErrorModelEAdd;
    double* m_pfErrorModelEMul;
    double  m_fErrorModelEAddMinIntensity;

public:
    // Public variables.

    int m_nVerbose;
    float m_fChiSquared;            // Store the overall chi squared here.
    int m_nGroups;                  // Number of sym-groups overall.
    int m_nNumBatches;              // Number of batches in the reflection file.
    int m_nNumUsed;                 // Store the overall number of reflections participating in some sym-related-group of size >=2

    double m_fTargetChiSquare;              // Desired chi-squared value.
    double m_fUserEAdd;                     // Set by user before we calculate the error model using nFitChiSquared().
    double m_fUserEAddMult;					// Multiplier to adjust Eadd.  If this is specified, we use the Rmerge in each resolution shell.
    double m_fUserEMul;                     // Set by user before we calculate the error model using nFitChiSquared().
    double m_fLastEAddUsed;                 // Loaded by fCalcSigma().  Gives the last Eadd used.
    double m_fLastEMulUsed;                 // Loaded by fCalcSigma().  Gives the last Emul used.

    double m_fResolutionMin;
    double m_fResolutionMax;
    double m_fFractionReject;                // Fraction of data to reject.  This sets m_fDevReject if not equal to 0.0.
    double m_fReject;                        // Equivalent reflection merge sigma reject.
    double m_fMaxDevReject;                  // Maximum m_fDevReject.  This really only affects very bad data sets.
    double m_fDevReject;                     // Maximum deviation (normalized by adjusted sigma) 
    double m_fBatchIoverSigReject;           // Reject a batch who's <I/sig> value deviates by this many sigmas from average.
    double m_fBatchRmergeReject;             // Reject a batch who's Rmerge value deviates by this many sigmas from average.
    double m_fBatchPercentReject;            // Reject a batch that has this many or more rejections.
	double m_fBatchChiSqReject;				 // Reject a batch that has a bad chi-squared value.
    double m_fSigma;                         // I/sigma exclcusion.
    double m_fMaxScaleIntensity;             // Maximum intensity to use in scaling.
    bool   m_bRejectBatches;                 // Is batch rejection enabled?
    Cstring m_sBatchRejectTemplate;          // Batches explicitly rejected.    
    Creflnlist* m_poPrevRejectList;          // Previously rejected reflections.
    Creflnlist* m_poPartialList;             // Partial list used by nTruncate()

    bool  m_bExcludeSysAbsences;
    bool  m_bScaleAnom;
    bool  m_bComputeOverlap;
    bool  m_bUseBackgroundSigma;

    static Cstring ms_sf2STLsq;
    static Cstring ms_snBatchIdx;
    tagBatch* poGetBatchInfo() { return m_ptBatchInfo; };

private:

public:
    CResoRASBins      m_oRASBins;
    // Error-scaling code.

    double fCalcChiSquared(int nRecomputeRejects);
    double fRmerge(bool bIgnore1=TRUE,bool bIgnore2=FALSE); // Quick and dirty R-merge.

    //shijie yao : 03:26:08 the other merge factors
    double fRmerge_RIM(bool bIgnore1=TRUE,bool bIgnore2=FALSE); // Redundance-independent merging R factor.
    double fRmerge_PIM(bool bIgnore1=TRUE,bool bIgnore2=FALSE); // Precision-indicating merging R factor.
    double fRmerge_PCV(bool bIgnore1=TRUE,bool bIgnore2=FALSE); // Pooled coefficient of variation factor.

    double fCalcSigma(Crefln& oRefln,int nRefNumber);       // Calculates the sigma for a given reflection using the type of sigma model.
    double fFindDevRejection(double fRejectPercent,double* pfChiSq = NULL);
    int nRejectBatches(const Cstring& sBatchTemplateStart,const Cstring& sBatchTempateEnd);            // Rejects certain batch names.  If psBatchTemplate == NULL, clears all rejections.
    int nCompletenessCutoff(double fCutoff,Cstring& sSort1,Cstring& sSort2,bool bReverse);
    int nExcludePriorReflections();

    bool bCalcAnomalousSignal(Creflnlist *poReflnlist );
    bool bDoAnomalousStatisticsOnSymEquivGroup(std::vector<Crefln*>& vecSymEquivReflns,
                                               int nResoField,
                                               int nPlusMinusFlagField,
                                               int nCentricityField);
    int nReportRejects();
    int nChiPlot();
    int nScatterPlot();
    int nPrintEaddEmul();
    void vInitErrorModel();
    int nFitChiSquared(bool bPrint = true);


    // Scaling code.

    int nBatchScale(Cstring& sFixedBatch,bool bPrint);
    int nBFactorScale(Cstring& sFixedBatch,bool bPrintResults);
    int nReBatch(int nMinBatchOverlaps,bool bRebatchImageBoundry,bool bRebatchBookends,bool bRebatch);
    int nSphericalHarmonics(tScaleProcessOptions& oOption);
    int nRememberIntensities(int nScaleType,int nMode);  // 0==remember, 1=get scale and unapply, 2= apply scale

    // Tables.
    int nTableIntensityResolutionPosition();
    int nTableMultiplicity();
    int nTableOverlap();
    int nTableRMergeVsBatch();
    int nTableRMergeVsResolution(bool bSummary=true);
    int nTableRMergeVsIntensity();

    //shijie yao : additional tables for the 3 new merging factors
    int nTableRrimVsBatch();
    int nTableRrimVsResolution(bool bSummary=true);
    int nTableRrimVsIntensity();
    int nTableRpimVsBatch();
    int nTableRpimVsResolution(bool bSummary=true);
    int nTableRpimVsIntensity();
    int nTable3RVsResolution(bool bSummary=true);     //a table of Rmerge, Rmeas, Rpim vs Resolution

    int nTableCompleteness();
    int nTableChiDistribution();
    int nTableAnomSignal();

    // Output
#ifdef SSI_PC
    int nWriteReflns(Cstring& sName,Cstring& sUnavgName,Cstring& sRejectsName,Cstring& sScaleFileName,bool bWriteHeader,bool bWriteUnavgHeader,bool bTexsan,bool bTexsan2,bool bShelx,bool bAnom,bool bTexsan2andShelx);
#else
    int nWriteReflns(Cstring& sName,Cstring& sUnavgName,Cstring& sRejectsName,Cstring& sScaleFileName,bool bWriteHeader,bool bWriteUnavgHeader,bool bTexsan,bool bTexsan2,bool bShelx,bool bAnom);
#endif

    bool bIsAvailable() { return (!m_nStat); };

    int nInitArrays();
    static void vExpandReflnlist(Creflnlist *poReflnlistIn);

    Cscaleaverage(Creflnlist& oListIn,Cimage_header& oHeader);

    ~Cscaleaverage();

private:
    std::vector<SCALE_RESTRAIN_INFO>       m_vecScaleRestrInfo;
public:
    void vSetScaleRestrainInfo(Cstring strInfo);
    void vGenerateDefaultBatchRestrainInfo(double dWeight);
    void vConvertBatchRestrainInfoToBatchIndices();
};

enum            _eBatchRejects{ BATCH_INTENSITY,BATCH_RMERGE,BATCH_PERCENT,BATCH_CHISQ,BATCH_BATCHNAME};
enum            _eTable {TABLE_RMERGEVSBATCH,TABLE_RMERGEVSINTENSITY,TABLE_RMERGEVSRESOLUTION,
                        TABLE_COMPLETENESS,TABLE_MULTIPLICITY,TABLE_CHIDISTRIB,TABLE_ANOMSIGNAL,TABLE_INPUTSTATS,TABLE_ALL};

struct DTREK_EXPORT tScaleProcessOptions
{
    static char* ms_apcTableNames[];

    bool bBatchScale;
    bool bBFactorScale;
    bool bSphericalAbsorb;
    bool bList;
    bool bRebatch;
    bool bCompletenessReject;

    // Sabsorb options.
    int             nHarmonicOrderDifrac;
    int             nHarmonicOrderIncident;
    int             nMaxIter;

    // Linear absorption specific variables.
    float           fDegreesPerS0;

    // Rebatching.
    int             nMinOverlaps;
    bool            bRebatchImageBoundary;
    bool            bRebatchBookends;

    // Truncation
    Cstring         sTruncateFile;

    // Table data.
    _eTable         eTable;

    // Sigma / Reject
    float           fSigma;
    float           fReject;
    float           fFractionChiRejects;     // This is set to 0.0 unless the 'auto' mode is used.
    float           fBatchRejectTypes[5];    // See enumerations above for indices.
    bool            bBatchRejectFlags[5];
    Cstring         sBatchRejectTemplate;    // Rejecting batches a la carte from user.
    Cstring         sPrevRejectFile;         // Previous rejections written to a file.

    Cstring         sBatchRestrInfo;         // A list of batches to be restrained while scaling 


    // Error scaling variables.
    float           fAdd;
    float           fMul;
	float			fAddMult;

    // Completeness Reject variables;
    int             nSign;
    double          fCompleteness;
    Cstring         sField1;
    Cstring         sField2;


    void vReset()
    {
        bSphericalAbsorb =      FALSE;

        bBatchScale   =         FALSE;
        bBFactorScale =         FALSE;
        bList         =         FALSE;
        bRebatch           =    FALSE;
        bCompletenessReject =   FALSE;

        
        // Initialize table data.
        eTable = TABLE_ALL;

        // Rejection
        fReject             = 200000000;
        fSigma              = 3.0;
        fFractionChiRejects  = 0.005f;
        bBatchRejectFlags[BATCH_INTENSITY] = FALSE;
        bBatchRejectFlags[BATCH_RMERGE] = FALSE;
        bBatchRejectFlags[BATCH_PERCENT] = FALSE;
        bBatchRejectFlags[BATCH_BATCHNAME] = FALSE;
        bBatchRejectFlags[BATCH_CHISQ] = FALSE;
        sBatchRejectTemplate = "";
        sBatchRestrInfo = "";
       
        fDegreesPerS0 = -1.0;
        
        // Initialize Spherical Absorption Correction
        nHarmonicOrderDifrac = 4;
        nHarmonicOrderIncident = 1;
        nMaxIter = 2;
        
        fAdd=-1.0;        // Signifies that values have not been specified.
		fAddMult = 1.0;	  // Signifies that we are determining the Eadd as 1.0*Rmerge of resolution shell.
        fMul=-1.0;         
    };

    tScaleProcessOptions() {
        vReset();
    };
};

#endif   // DT_CSCALEAVERAGE_H
