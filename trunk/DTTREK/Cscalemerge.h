#ifndef DT_CSCALEMERGE_H
#define DT_CSCALEMERGE_H
//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cscalemerge.h    Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class Cscalemerge
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

//+Include files

#include "Dtrek.h"
#include "Cstring.h"
#include "Crefln.h"
#include "Creflnlist.h"
#include "Ccrystal.h"  
#include "Cimage_header.h"

//+Definitions and constants
enum eScalemerge_states {
  eScalemerge_unknown_state,
  eScalemerge_available_state
};

typedef struct _tagStats
{
  int   nNumRefs;      // Number of reflns
  int   nNumRejects;   
  int   nNumUsed;
  int   nNumMults;
  int   anNumInBinMults[10];
  int   anNumInBin[10];
  int   nNumSingles;

  double dNumer;
  double dDenom;
  double dWtNumer;
  double dWtDenom;
  double dSumIoverSig;
  double dSumChiSq;
  double dRmerge;
} tagStats;

//+Code begin

class DTREK_EXPORT Cscalemerge {

protected:
  eScalemerge_states 
              m_eThe_State;       // The state of the reflection list

  Cstring     m_sScalemerge_key;   // Unique Scalemerge key for database lookup
  Cstring     m_sBatchFixed;
  Cstring    *m_psBatchNames;
  float       m_fScaleFixed;
  float       m_fBValueFixed;
  int         m_nBatchFixed;
  Ccrystal   *m_poCrystal;
//made public  Creflnlist *m_poReflnlist;
  Creflnlist *m_poReflnlistAve;
  Creflnlist *m_poReflnlistComp;  // Reflnlist object for completeness
  bool        m_bNewList;

  bool        m_bSigmaOverEstimate;
  int         m_nAnomFlag;
  int         m_nFixBFlag;
  int         m_nMethod;
  int         m_nVerbose;
  int         m_nCountOverlap;
  int         m_nScaleAnom;    // Must be 0 (treat I+ and I- anom separately)
                               // or 1 (treat anom as same)

  float       m_fResolutionMin;// Min resolution limit in A (99999 is min)
  float       m_fResolutionMax;// Max resolution limit in A (0.0001 A is max)
  float       m_fRejCriteria;
  int         m_nMaxCycles;
  bool        m_bLastCycle;

// Some variables used by the Kabsch algorithm 
// See Kabsch (1988) J.Appl.Cryst. 21, 916-924.
// From p. 923:

  float      *m_pfAab;
  float      *m_pfBa;
  float      *m_pfIh;
  float      *m_pfWhla;
  float      *m_pfRha;
  float      *m_pfUha;
  float      *m_pfVha;
  float      *m_pfUh;
  float      *m_pfVh;
  float      *m_pfGa;
  int        *m_pnGaFixed;
//
  float      *m_pfW1;            // Some temporary work vectors
  int        *m_pnW2;
  float      *m_pfW3;
  float      *m_pfW4;
  int        *m_pnW5;
  int        *m_pnIndex;         // An array for sorting reflns
  int         m_nMaxSymRelated;
  int         m_nNumBatches;
  int         m_nNumSubBatchesPerBatch;
  int         m_nNumParams;
  int        *m_pnNumContrib;
  int        *m_pnNumInBatch;
  int        *m_pnNumAccepted;
  int        *m_pnNumRejected;
  int        *m_pnNumDeselected;
  int        *m_pnNumOverlap;

  int        *m_pnNumMultsBatch;

  float       m_a2fRangeReso[2];
  float       m_fSlopeReso;
  float       m_a2fRangeInt[2];
  float       m_fSlopeInt;
  int         m_nNumBinsReso;
  int         m_nNumBinsInt;
  tagStats   *m_ptStatsReso;
  tagStats   *m_ptStatsInt;
  tagStats   *m_ptStatsBatch;

  float       m_fErrorMul;
  float       m_fErrorAdd;
  float       m_fEigRho;

  static Cstring ms_sf2STLsq;
  static Cstring ms_snBatchIdx;

public:
  int         m_nFI_f2STLsq;
  int         m_nFI_nBatchIdx;
  int         m_nFI_nAnomFlag;
  int         m_nFI_fIntensityPlus;
  int         m_nFI_fSigmaIPlus;
  int         m_nFI_fIntensityMinus;
  int         m_nFI_fSigmaIMinus;
//  int         m_nFI_fOtherSig;
//  int         m_nFI_fOtherInt;

  Creflnlist *m_poReflnlist;

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments
////////////////////////////////////////////////////////////////////////

  Cscalemerge(Creflnlist *poReflnlistIn, Cimage_header *poHeaderIn);
  ~Cscalemerge ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes (none are virtual and cannot be overridden!?!
////////////////////////////////////////////////////////////////////////

int  nInitValues(void);
int  nList(const int nNum = 0);
int  nScaleSetup(void);
int  nScale1(const int nMaxCycles = -1, Creflnlist *poReflnlistAve = NULL);

float fAverageSymRelated(const int nMethod, const int nNumRef, 
			 const int *pnIndex, const float fRejCriteria,
			 Crefln *poReflnAve,
			 int *pnAccepted, int *pnNewRejects, float *pfSumResid);

virtual float fCalcGetWeight(const float fIntensity, const float fSigma,
			     Crefln *poRefln);
virtual float fCalcGetScale(Crefln *poRefln);
virtual void  vCalcGetWeights(Crefln *poRefln);
virtual void  vApplyWeights(void);

void vKabsch(const int nNumRefs, const int *pnIndex);
void vKandB(const int nNumRefs, const int *pnIndex);
void vKandB2(const int nNumRefs, const int *pnIndex);
void vRoughScale(void);
void vInitStats(tagStats *ptStats);
void vListStats(void);
void vListOverlap(void);
int  nListCompleteness(void);
void vSetResolution(const float fResoMin, const float fResoMax);

inline bool bIsAvailable(void)
  {  return (m_eThe_State != eScalemerge_unknown_state); }

inline void vSetRejCrit(const float fRej) { m_fRejCriteria = fRej; }
inline void vSetScale(const float fScale) { m_fScaleFixed = fScale; }
inline void vSetTemp(const float fBValue) { m_fBValueFixed = fBValue; }
inline void vSetBatch(const Cstring& sBatch) { m_sBatchFixed = sBatch; }
inline void vSetError(const float fMul, const float fAdd)
           { m_fErrorMul = fMul; m_fErrorAdd = fAdd; }
inline void vSetMaxCycles(const int nCyc) { m_nMaxCycles = nCyc; }
inline void vSetOutputAnomFlag(const int nFlag) { m_nAnomFlag = nFlag; }
inline void vSetScaleAnomFlag(const int nFlag) { m_nScaleAnom = nFlag; }
inline void vSetFixBFlag(const int nFlag) { m_nFixBFlag = nFlag; }
inline void vSetVerbose (const int nFlag) { m_nVerbose  = nFlag; }
inline void vSetCountOverlap (const int nFlag) { m_nCountOverlap  = nFlag; }
inline void vSetSigmaOverEstimate(const bool bFlag) 
			       { m_bSigmaOverEstimate = bFlag; }
private:

void vDeleteArrays(void);

};  // end of class Cscalemerge

#endif   // DT_CSCALEMERGE_H
