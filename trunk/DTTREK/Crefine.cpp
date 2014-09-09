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
// Crefine.cc            Initial author: J.W. Pflugrath           01-May-1995
//  This file contains the member functions of class Crefine which implements
//    the refinement encapsulation of d*TREK.
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
//+ToDo
//
//   Error messages need implementing
//
//+Include files
#if (defined(SSI_PC) || (defined(NO_X_WINDOWS)))
#include <string.h>
#endif
#include "Dtrek.h"
#include "dtreksys.h"
#include "Crefine.h"         // Class definition and prototypes
#include "dtrekdefs.h"
#include "Cpredict.h"
#include "Cscan.h"
#include "Cfind.h"
#include "Cstat.h"
#include "CScanBitmap.h"
#include "Cprofit2D.h"

#ifdef SSI_PC
#include "CrclHelper.h"
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::cerr;
#endif

int Crefine::ms_nDetTrans1  =  0;
int Crefine::ms_nDetTrans2  =  1;
int Crefine::ms_nDetTrans3  =  2;
int Crefine::ms_nDetRot1    =  3;
int Crefine::ms_nDetRot2    =  4;
int Crefine::ms_nDetRot3    =  5;
int Crefine::ms_nGonioMiss1 =  6;
int Crefine::ms_nGonioMiss2 =  7;
int Crefine::ms_nSourceRot1 =  8;
int Crefine::ms_nSourceRot2 =  9;
int Crefine::ms_nSourceWave = 10;
int Crefine::ms_nCrysRot1   = 11;
int Crefine::ms_nCrysRot2   = 12;
int Crefine::ms_nCrysRot3   = 13;
int Crefine::ms_nCrysAstar  = 14;
int Crefine::ms_nCrysBstar  = 15;
int Crefine::ms_nCrysCstar  = 16;
int Crefine::ms_nCrysAlps   = 17;
int Crefine::ms_nCrysBets   = 18;
int Crefine::ms_nCrysGams   = 19;
int Crefine::ms_nCrysMosaicity = 20;
int Crefine::ms_nCrysConstraints = 21;
int Crefine::ms_nCrysRecip00 = 22;
int Crefine::ms_nCrysRecip01 = 23;
int Crefine::ms_nCrysRecip02 = 24;
int Crefine::ms_nCrysRecip10 = 25;
int Crefine::ms_nCrysRecip11 = 26;
int Crefine::ms_nCrysRecip12 = 27;
int Crefine::ms_nCrysRecip20 = 28;
int Crefine::ms_nCrysRecip21 = 29;
int Crefine::ms_nCrysRecip22 = 30;
int Crefine::ms_nCrysRecip30 = 31;
int Crefine::ms_nCrysRecip31 = 32;
int Crefine::ms_nCrysRecip32 = 33;
int Crefine::ms_nCrysTwinLaw10 = 34;
int Crefine::ms_nCrysTwinLaw11 = 35;
int Crefine::ms_nCrysTwinLaw12 = 36;
int Crefine::ms_nCrysTwinLaw20 = 37;
int Crefine::ms_nCrysTwinLaw21 = 38;
int Crefine::ms_nCrysTwinLaw22 = 39;
int Crefine::ms_nCrysTwinLaw30 = 40;
int Crefine::ms_nCrysTwinLaw31 = 41;
int Crefine::ms_nCrysTwinLaw32 = 42;
int Crefine::ms_nFirstCrysParam = ms_nCrysRot1;
int Crefine::ms_nLastCrysParam = ms_nCrysTwinLaw32;


Cstring Crefine::ms_sDetTrans1       = "DetTrans1";
Cstring Crefine::ms_sDetTrans2       = "DetTrans2";
Cstring Crefine::ms_sDetTrans3       = "DetTrans3";
Cstring Crefine::ms_sDetTrans        = "DetTrans";
Cstring Crefine::ms_sDetRot1         = "DetRot1";
Cstring Crefine::ms_sDetRot2         = "DetRot2";
Cstring Crefine::ms_sDetRot3         = "DetRot3";
Cstring Crefine::ms_sDetRot          = "DetRot";
Cstring Crefine::ms_sDetAll          = "DetAll";
Cstring Crefine::ms_sGonioMiss1      = "GonioMiss1";
Cstring Crefine::ms_sGonioMiss2      = "GonioMiss2";
Cstring Crefine::ms_sGonioMiss       = "GonioMiss";
Cstring Crefine::ms_sSourceRot1      = "SourceRot1";
Cstring Crefine::ms_sSourceRot2      = "SourceRot2";
Cstring Crefine::ms_sSourceRot       = "SourceRot";
Cstring Crefine::ms_sSourceWave      = "SourceWave";
Cstring Crefine::ms_sCrysRot1        = "CrysRot1";
Cstring Crefine::ms_sCrysRot2        = "CrysRot2";
Cstring Crefine::ms_sCrysRot3        = "CrysRot3";
Cstring Crefine::ms_sCrysRot         = "CrysRot";
Cstring Crefine::ms_sCrysAstar       = "CrysAstar";
Cstring Crefine::ms_sCrysBstar       = "CrysBstar";
Cstring Crefine::ms_sCrysCstar       = "CrysCstar";
Cstring Crefine::ms_sCrysLengths     = "CrysLengths";
Cstring Crefine::ms_sCrysAlps        = "CrysAlps";
Cstring Crefine::ms_sCrysBets        = "CrysBets";
Cstring Crefine::ms_sCrysGams        = "CrysGams";
Cstring Crefine::ms_sCrysAngles      = "CrysAngles";
Cstring Crefine::ms_sCrysCell        = "CrysCell";
Cstring Crefine::ms_sCrysMosaicity   = "CrysMosaicity";
Cstring Crefine::ms_sCrysAll         = "CrysAll";
Cstring Crefine::ms_sCrysConstraints = "CrysConstraints";
Cstring Crefine::ms_sCrysRecipShift  = "CrysShift";
Cstring Crefine::ms_sCrysRecipShift0 = "CrysShift0";
Cstring Crefine::ms_sCrysRecipShift1 = "CrysShift1";
Cstring Crefine::ms_sCrysRecipShift2 = "CrysShift2";
Cstring Crefine::ms_sCrysRecipShift3 = "CrysShift3";
Cstring Crefine::ms_sCrysTwinLaw     = "TwinLaw";
Cstring Crefine::ms_sCrysTwinLaw1    = "TwinLaw1";
Cstring Crefine::ms_sCrysTwinLaw2    = "TwinLaw2";
Cstring Crefine::ms_sCrysTwinLaw3    = "TwinLaw3";
  
Cstring Crefine::ms_sAll             = "All";


Cstring* Crefine::ms_psFlags[nMAXPARAMS] = {
        &ms_sDetTrans1,  
        &ms_sDetTrans2,  
        &ms_sDetTrans3,  
        &ms_sDetRot1,    
        &ms_sDetRot2,    
        &ms_sDetRot3,    
        &ms_sGonioMiss1,
        &ms_sGonioMiss2,
        &ms_sSourceRot1, 
        &ms_sSourceRot2, 
        &ms_sSourceWave, 
        &ms_sCrysRot1,   
        &ms_sCrysRot2,   
        &ms_sCrysRot3,   
        &ms_sCrysAstar,  
        &ms_sCrysBstar,  
        &ms_sCrysCstar,  
        &ms_sCrysAlps,   
        &ms_sCrysBets,   
        &ms_sCrysGams,   
        &ms_sCrysMosaicity,
        &ms_sCrysConstraints,
        &ms_sCrysRecipShift0,   // It's okay to have these guys point to the same entry since we don't allow
        &ms_sCrysRecipShift0,   // indivisual refinement of these parameters.
        &ms_sCrysRecipShift0,
        &ms_sCrysRecipShift1,
        &ms_sCrysRecipShift1,
        &ms_sCrysRecipShift1,
        &ms_sCrysRecipShift2,
        &ms_sCrysRecipShift2,
        &ms_sCrysRecipShift2,
        &ms_sCrysRecipShift3,
        &ms_sCrysRecipShift3,
        &ms_sCrysRecipShift3,
                &ms_sCrysTwinLaw1,
                &ms_sCrysTwinLaw1,
                &ms_sCrysTwinLaw1,
                &ms_sCrysTwinLaw2,
                &ms_sCrysTwinLaw2,
                &ms_sCrysTwinLaw2,
                &ms_sCrysTwinLaw3,
                &ms_sCrysTwinLaw3,
                &ms_sCrysTwinLaw3,
};

double Crefine::ms_pfRefineConvergence[nMAXPARAMS] = {
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003,
    0.000003
};


Cstring Crefine::ms_sDtrefineFiles   = D_K_DtrefineFiles;
Cstring Crefine::ms_sDtrefineOptions = D_K_DtrefineOptions;
Cstring Crefine::ms_sDtrefineRmsMm   = D_K_DtrefineRmsMm;
Cstring Crefine::ms_sDtrefineRmsDeg  = D_K_DtrefineRmsDeg;
Cstring Crefine::ms_sDtrefineReflectionNumbers = D_K_DtrefineReflectionNumbers;
Cstring Crefine::ms_sDtrefineReflectionNumbersPhotons = D_K_DtrefineReflectionNumbersPhotons;
Cstring Crefine::ms_sDtrefineTwinsRefined = D_K_DtrefineTwinsRefined;
Cstring Crefine::ms_sDtrefineTwinRejected = D_K_DtrefineTwinRejected;
Cstring Crefine::ms_sDtrefineLambda         = D_K_DtrefineLambda;
Cstring Crefine::ms_sDtrefineLambdaProfile  = D_K_DtrefineLambdaProfile;
Cstring Crefine::ms_sDtrefineResidProfile   = D_K_DtrefineResidProfile;
Cstring Crefine::ms_sDtrefineProfileCount   = D_K_DtrefineProfileCount;
Cstring Crefine::ms_sDtrefineResoUsed       = D_K_DtrefineResoUsed;

char* Crefine::m_pcVarNames[] = {
      "Det Trans 1",
      "Det Trans 2",
      "Det Trans 3",
      "Det Rot 1",
      "Det Rot 2",
      "Det Rot 3",
      "Gonio Miss 1",
      "Gonio Miss 2",
      "Src Rot 1",
      "Src Rot 2",
      "Src Wave",
      "Crystal Rot1",
      "Crystal Rot2",
      "Crystal Rot3",
      "Crys-Param1",
      "Crys-Param2",
      "Crys-Param3",
      "Crys-Param4",
      "Crys-Param5",
      "Crys-Param6",
      "Mosaicity",
      "Crys-Constraints",
      "CrysShift00",
      "CrysShift01",
      "CrysShift02",
      "CrysShift10",
      "CrysShift11",
      "CrysShift12",
      "CrysShift20",
      "CrysShift21",
      "CrysShift22",
      "CrysShift30",
      "CrysShift31",
      "CrysShift32",
          "TwinLaw10",
          "TwinLaw11",
          "TwinLaw12",
          "TwinLaw20",
          "TwinLaw21",
          "TwinLaw22",
          "TwinLaw30",
          "TwinLaw31",
          "TwinLaw32",
      "Bogus1",
      "Bogus2"
};
//+Code begin

//+Public functions

// Constructors, destructors and assignments

Crefine::Crefine()
{
 if (nInitValues())
     m_bIsAvailable = false;
}

Crefine::Crefine(Cimage_header& oHeader, Creflnlist *poReflnlistInput)
{
 if (nInitValues(oHeader))
     m_bIsAvailable = false;
 m_poReflnlistIn   = poReflnlistInput;
 m_poHeader = &oHeader;
}

Crefine::Crefine (Csource *poSourceIn, Cdetector *poDetectorIn,
                  Ccrystal *poCrystalIn, Cgoniometer *poCrysGonioIn,
                  Crotation *poRotationIn, Creflnlist *poReflnlistIn)
{
  (void) nInitValues();

  // There should be some error checking below to make sure that the
  // pointers actually point to real objects.

  m_poSource      = poSourceIn;
  m_bNewSource    = FALSE;
  m_ppoCrystals   = new Ccrystal*[1];
  m_poComp        = new CrefineTwinData[1];
  m_ppoCrystals[0]= poCrystalIn;
  m_bNewCrystal   = FALSE;
  m_poCrysGonio   = poCrysGonioIn;
  m_bNewCrysGonio = FALSE;
  m_poRotation    = poRotationIn;
  m_bNewRotation  = FALSE;

  m_nNumDetectors = 1;
  m_ppoDetector   = new Cdetector* [m_nNumDetectors];
  *m_ppoDetector  = poDetectorIn;
  if (NULL != poDetectorIn)
    {
      int i;
      for (i = 0; i < 6; i++)
        {
          m_a6sUserDetNames[i] = "Det"
                              + poDetectorIn->m_poGoniometer->sGetName(i);
          if (m_a6sUserDetNames[i].contains("/"))
            m_a6sUserDetNames[i] = m_a6sUserDetNames[i].before("/");
        }
    }

  m_bNewDetector  = FALSE;
  m_fLambda       = -1;
  m_poReflnlistIn = poReflnlistIn;

}

Crefine::~Crefine()
{
    if (m_bNewDetNames && (NULL != m_psDetectorNames) )
    {
        delete [] m_psDetectorNames;
        m_psDetectorNames = NULL;
        m_bNewDetNames    = FALSE;
    }
    
    if (NULL != m_poReflnlist)
    {
        delete m_poReflnlist;
        m_poReflnlist = NULL;
    }
    
    if (m_bNewDetector && (NULL != m_ppoDetector) )
    {
        for (int i = 0; i < m_nNumDetectors; i++)
        {
            delete *(m_ppoDetector+i);
            *(m_ppoDetector+i) = NULL;
        }
        m_bNewDetector = FALSE;
    }
    if (NULL != m_ppoDetector)
    {
        delete [] m_ppoDetector;
        m_ppoDetector = NULL;
    }
    
    if (m_bNewSource && (NULL != m_poSource) )
    {
        delete m_poSource;
        m_poSource   = NULL;
        m_bNewSource = FALSE;
    }
    
    if (m_bNewCrystal && (NULL != m_ppoCrystals) )
    {
        for (int i = 0; i < m_nNumDetectors; i++) {
            delete m_ppoCrystals[i];
            m_ppoCrystals[i] = NULL;
        };
        m_bNewCrystal = FALSE;
    }
    if (NULL != m_ppoCrystals)
    {
        delete [] m_ppoCrystals;
        m_ppoCrystals = NULL;
        m_ppoDetector = NULL;
    }
    if (NULL != m_poComp)
    {
        delete[] m_poComp;
        m_poComp = NULL;
    };
    
    if (m_bNewCrysGonio && (NULL != m_poCrysGonio) )
    {
        delete m_poCrysGonio;
        m_poCrysGonio   = NULL;
        m_bNewCrysGonio = FALSE;
    }
    
    if (m_bNewRotation && (NULL != m_poRotation) )
    {
        delete m_poRotation;
        m_poRotation     = NULL;
        m_bNewRotation = FALSE;
    }
    
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
    if (NULL != m_poXprop)
    {
        delete m_poXprop;
        m_poXprop = NULL;
    }
#endif
    if (m_pfResidProfile)
        delete[] m_pfResidProfile;
    if (m_pfLambdaProfile)
        delete[] m_pfLambdaProfile;
    
    m_nNumDetectors  = 0;
    if (NULL != m_psVarNames)
    {
        delete [] m_psVarNames;
        m_psVarNames = NULL;
    }
}

int Crefine::nInitValues(void)
{
  int i , j;  // loop counter
  m_bIsAvailable     = TRUE;
  m_bPrintTables     = FALSE;
  m_bUseOverlapped   = false;
  m_nDisplay         = 0;
  m_nVerbose         = 1;
  m_nAutoCorrCheck   = 1;
  m_nCycles          = 30;
  m_nAdjustRotation  = TRUE;
  m_nNumDetectors    = 0;
  m_nWhichDetector   = 0;
  m_fMaxCorrCoeff    = 0.98f;
  m_fRMSMM           = 0.0;
  m_fRMSDEG          = 0.0;
  m_fRMSANG          = 0.0;

  m_bCheckConvergence= TRUE;
  m_bUseWideSlice    = FALSE;
  m_bNotUseWideSlice = FALSE;
  m_eWeightType      = eWeightNone;

  m_psDetectorNames  = NULL;
  m_poReflnlist      = new Creflnlist();  // Create the working refln list
  m_ppoDetector      = NULL;
  m_poSource         = NULL;
  m_ppoCrystals      = NULL;
  m_poComp          = NULL;
  m_poCrysGonio      = NULL;
  m_poRotation       = NULL;
  m_poHeader         = NULL;

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  m_poXprop          = NULL;
#endif

  m_asWeightNames[eWeightRotSigma] = "Rot sigma";
  m_asWeightNames[eWeightIntensity] = "Intensity";
  m_asWeightNames[eWeightIoverSig]  = "Int/sigInt";
  m_asWeightNames[eWeightLorentz]   = "Inverse Lorentz";
  m_asWeightNames[eWeightRotWidth]  = "Rot width";
  m_asWeightNames[eWeightNone]      = "Unit weights";
  m_bNewDetNames     = FALSE;
  m_bNewDetector     = FALSE;
  m_bNewSource       = FALSE;
  m_bNewCrystal      = FALSE;
  m_bNewCrysGonio    = FALSE;
  m_bNewRotation     = FALSE;
  m_bNoReflectionMerges = FALSE;

  m_bReIndex         = TRUE;
  m_bIsRhomboAsHex   = FALSE;
  m_bPrintDetectorMosaicity = FALSE;

  m_fResidMin        = 1e-20f;
//+JWP 2010-02-23
  m_fLorentzMax      = 50.0;
  m_fLorentzMax      = 35.0; // Use a different default max for refinement
//-JWP 2010-02-23
  if ("" != sGetEnv("DTREK_MAXLORENTZ"))
    m_fLorentzMax = (float)atof(sGetEnv("DTREK_MAXLORENTZ").string());
  m_pfResidProfile   = NULL;
  m_pfLambdaProfile  = NULL;
  m_nProfileCt       = 0;

  m_a4fRejectLimits[0] = 2.5f;
  m_a4fRejectLimits[1] = 2.5f;
  m_a4fRejectLimits[2] = 2.0f;
  m_a4fRejectLimits[3] = m_a4fRejectLimits[2];
  m_fSigma           = 0.0f;
  m_fSharpness       = 0.2f;
  m_fTweakSensitivity = 1.0;
  m_fResolutionMin   = 99999.0f;
  m_fResolutionMax   = 0.00001f;
  m_nTwinID          = 1;
  m_nNumTwins        = 1;
  m_nMaxTwinLaws     = 0;
  m_fFindPad         = 0.5;
  m_sRejectFile          = "";
  m_sMapFile         = "";
  m_sLogTwinReject   = "";
  m_sLogRefine       = "";
  m_sLogRefineLast   = "";
  m_anReflnOverlapMat.clear();        
  m_afReflnOverlapIntensities.clear();
  m_afReflnMinResid.clear();
  m_anReflnBestTwin.clear();

  m_fTestMosaicityRangeMax = 6.0;
  m_fTestMosaicityRangeMin = 0.05f;

  for (i = 0; i < 6; i++)
    {
      m_a6fCrysUnitCell[i]        = 90.0;
      m_a6fCrysTurds[i]           = 0.0;
      m_a6fCrysUnitCellSigs[i]    = 0.0;
    }
  for (i = 0; i < 3; i++)
    {
      m_a3fCrysRot[i]             = 0.0;
      m_a3fCrysRotSigs[i]         = 0.0;
      m_a3fCrysRotShifts[i]       = 0.0;
    }
  for (i = 0; i < 10; i++)
      for (j = 0; j< 3; j++)
        {
          m_a10x3fDetRotSigs[i][j]   = 0.0;
          m_a10x3fDetTransSigs[i][j] = 0.0;
        };
  for (i = 0; i < 2; i++)
      m_a2fSourceRotSigs[i] = 0.0;

  for (i = 0; i < 3; i ++) {
    m_a3fCrysMosaicity[i]       = 0.0;
    m_a3fCrysMosaicitySig[i]    = 0.0;
      };
  m_a3fCrysMosaicity[0]       = 0.3f;

  for (i = 0; i <= g_nMaxTwinLaws; i++) {
          vZeroMat(3,1,&m_aa3fCrysRecipShift[i][0]);
          vZeroMat(3,1,&m_aa3fCrysRecipShiftSig[i][0]);
          vZeroMat(3,1,&m_aa3fCrysRecipShiftShift[i][0]);
          vZeroMat(3,1,&m_aa3fCrysTwinLaw[i][0]);
          vZeroMat(3,1,&m_aa3fCrysTwinLawSig[i][0]);
          vZeroMat(3,1,&m_aa3fCrysTwinLawShift[i][0]);
          m_afCrysTwinFraction[i] = 0.0;
          m_afCrysTwinFractionSig[i] = 0.0;
          m_afCrysTwinFractionShift[i] = 0.0;
  };




  m_fSourceWave                 = 1.54178f;
  m_fSourceWaveSig              = 0.0;
  m_fSourceWaveShift            = 0.0;

  m_nFI0                        = -1;
  m_nFI1                        = -1;
  m_nFI2                        = -1;
  m_nFI3                        = -1;
  m_nFI4                        = -1;
  m_nFIMM0                      = -1;
  m_nFIMM1                      = -1;
  m_nFIMM2                      = -1;
  m_fLambda                     = -1;

  m_sOutHeader                  = sDtrekGetPrefix() + "dtrefine.head";

  m_a6sUserDetNames[0] = ms_sDetRot1;
  m_a6sUserDetNames[1] = ms_sDetRot2;
  m_a6sUserDetNames[2] = ms_sDetRot3;
  m_a6sUserDetNames[3] = ms_sDetTrans1;
  m_a6sUserDetNames[4] = ms_sDetTrans2;
  m_a6sUserDetNames[5] = ms_sDetTrans3;

  m_psVarNames = new Cstring [nMAXPARAMS];
  for (i = 0; i < nMAXPARAMS; i++) {
    m_psVarNames[i] = m_pcVarNames[i];
    m_nUserRefineFlags[i] = 0;
  };
  m_nUserRefineFlags[ms_nCrysConstraints] = 1;

  // Fix all user refine flags to begin with

  (void) nSetRefineFlags(Cstring('-') + ms_sAll);

  m_bRanking = false;

  m_nNumberOfResoBins = 10;

  return 0;
}

int Crefine::nInitValues(Cimage_header& oHeader)
{
  int    i;       // Loop counter
  int    nStat, nTemp;
  Cstring sTemp;

  if (nInitValues())
      return 1;

  if (!oHeader.bIsAvailable())
    {
      cout << "Crefine ERROR: image header is not valid."
           << "  Cannot construct refine!\n";
      nStat = 1;
    }
  else
    {
      // Try to get required information from the header

      nStat = 0;

      Ccrystal oCrystal(oHeader,1);
      if (!oCrystal.bIsAvailable())
          nStat = 1;
      else {
          m_nNumTwins = oCrystal.nNumTwins();
                  m_nMaxTwinLaws = oCrystal.nMaxTwinLaws();
          m_ppoCrystals = new Ccrystal*[m_nNumTwins];
          m_poComp      = new CrefineTwinData[m_nNumTwins];
          for (i = 0; i < m_nNumTwins; i ++)
              m_ppoCrystals[i] = new Ccrystal(oHeader,i+1);
      };

      m_bNewCrystal   = TRUE;

      m_poCrysGonio   = new Cgoniometer(oHeader, Ccrystal::ms_sCrystalPrefix);
      m_bNewCrysGonio = TRUE;

      // Should m_poRotation come from SCAN_ROTATION or ROTATION keywords?????

      m_poRotation    = new Crotation(oHeader);
      m_bNewRotation  = TRUE;

      // Get the number of detectors

      nTemp = oHeader.nGetValue(Cdetector::ms_sDetectorNumber, &m_nNumDetectors);

      if ((0 == nTemp) && (nStat==0)) {
          // Get the detector names

          m_psDetectorNames = new Cstring [m_nNumDetectors];
          m_bNewDetNames = TRUE;
          nTemp = oHeader.nGetValue(Cdetector::ms_sDetectorNames, m_nNumDetectors,
                                    m_psDetectorNames);
          if (0 == nTemp)
            {
              m_ppoDetector = new Cdetector* [m_nNumDetectors];
              for (i = 0; i < m_nNumDetectors; i++)
                {
                  // We do not need the non-uniformity information

                  m_ppoDetector[i] = new Cdetector (oHeader,
                                                    m_psDetectorNames[i],
                                                    TRUE, FALSE);
                  if (0 == i)
                    {
                      int j;
                      for (j = 0; j < 6; j++)
                        {
                          m_a6sUserDetNames[j] = "Det"
                            + m_ppoDetector[0]->m_poGoniometer->sGetName(j);
                          if (m_a6sUserDetNames[j].contains("/"))
                            m_a6sUserDetNames[j] = m_a6sUserDetNames[j].before("/");
                          m_psVarNames[j] = m_a6sUserDetNames[j];
                        }
                    }
                }
              m_bNewDetector  = TRUE;
            }
          else
            {
              nStat++;
            }
        }
      else
        {
          nStat++;
        }
      m_poSource      = new Csource(oHeader);
      m_bNewSource    = TRUE;
  }
  
  return (nStat);
}

int Crefine::nList(const int nFlag)
{
  float fTempW, fRMin, fRMax, fRMaxEdge;
  float fTempS0[3];
  int i;

  cout << "Refine listing:\n";
  if (NULL != m_ppoCrystals)
  {
      
      for (i = 0; i < m_nNumTwins; i++)
          (void) m_ppoCrystals[i]->nList(nFlag);
  }
  else
  {
      cout << "     Crystal not defined.\n";
  }
  if (NULL != m_poCrysGonio)
    {
      (void) m_poCrysGonio->nList(nFlag);
    }
  else
    {
      cout << "     Crystal goniometer not defined.\n";
    }
  if (NULL != m_poSource)
    {
      (void) m_poSource->nList(nFlag);
      fTempW = m_poSource->fGetWavelength();
      m_poSource->vCalcGetS0(&fTempS0[0]);
    }
  else
    {
      cout << "     Source not defined.\n";
      fTempW = 0.0;
    }
  for (i = 0; i < m_nNumDetectors; i++)
    {
      if (NULL != m_ppoDetector[i])
        {
          (void) m_ppoDetector[i]->nList(nFlag);
          if (fTempW != 0.0)
            {
              (void) m_ppoDetector[i]->nGetResolution(fTempS0,
                                                 &fRMin, &fRMax, &fRMaxEdge);
              cout << "DetResolution min:  " << fRMin * fTempW
                   << "\nDetResolution max:  " << fRMax * fTempW
                   << "\nDetResolution edge: " << fRMaxEdge * fTempW << '\n';
            }
        }
      else
        {
          cout << "     Detector " << i << " not defined.\n";
        }
    }
  cout << "\nRefine resol min: " << m_fResolutionMin << '\n';
  cout << "Refine resol max: " << m_fResolutionMax << '\n';
  cout << "I/sigma cutoff:   " << m_fSigma << '\n';
  cout << "Sharpness cutoff: " << m_fSharpness << '\n';
  cout << "Rejection limits: " << m_a4fRejectLimits[0] << ", "
                               << m_a4fRejectLimits[1] << ", "
                               << m_a4fRejectLimits[2] << '\n';
  cout << "Weighting scheme: " << m_asWeightNames[m_eWeightType] << "\n\n" << flush;
  fflush(stdout);
  return (0);
}


int Crefine::nExpandReflnlist(Creflnlist *poReflnlistIn)
{

// Add required fields to the working reflection list if they are not
// already there..
// The new fields are named by static member variables of the Creflnlist class.
// Most fields below were not used, so they were commented out.
// Eventually they should be used.

  Creflnlist *poList;
  if (NULL != poReflnlistIn)
    {
      poList = poReflnlistIn;
    }
  else
    {
      poList = m_poReflnlist;
    }

  (void) poList->nExpandGetField(poList->ms_snDetNum);
  (void) poList->nExpandGetField(poList->ms_snNonunfFlag);
  (void) poList->nExpandGetField(poList->ms_sfObsPx0);
  (void) poList->nExpandGetField(poList->ms_sfObsPx1);
  (void) poList->nExpandGetField(poList->ms_sfObsRotMid);
  (void) poList->nExpandGetField(poList->ms_sfObsRotWidth);
  (void) poList->nExpandGetField(poList->ms_sfObsXmm);
  (void) poList->nExpandGetField(poList->ms_sfObsYmm);
  (void) poList->nExpandGetField(poList->ms_sfCalcPx0);
  (void) poList->nExpandGetField(poList->ms_sfCalcPx1);
  (void) poList->nExpandGetField(poList->ms_snPackedHKL);

  // Expand twin ID field if required.
  if ((m_ppoCrystals) && ((m_nNumTwins>1) || (m_nMaxTwinLaws>0))) {
      (void) poList->nExpandGetField(poList->ms_snTwinID);
          (void) poList->nExpandGetField(poList->ms_snOrignlReflnNum);
      (void) poList->nExpandGetField(poList->ms_snSourceRefNum);
      (void) poList->nExpandGetField(poList->ms_snCompeteRefNum);
  };

  return (0);
}

int
Crefine::nCopyReflns(Creflnlist *poReflnlistOther)
{
  int   i, nStat, nDetNum, nRef;
  float f1PX, f2PX, f1MM, f2MM, f3MM;

  // Copy from the input reflection list to the working reflection list
  // and calculate their mm coordinates from the pixel coordinates

  // Get rid of old reflections ...

  m_poReflnlist->vDeleteAll();
/*
  if (NULL != m_poReflnlist)
    {
      // Get rid of old reflections ...

      delete m_poReflnlist;
      m_poReflnlist = new Creflnlist();
    }
*/
  // Next expand the reflection list to the fields we need

  (void) nExpandReflnlist();

  // If the observed centroid fields are missing from the input list,
  // then add calculated centroid fields to the working list in the hopes
  // that calculated centroids are present in the input list and will
  // be transferred to the working list.

  Creflnlist *poList;

  if (NULL != poReflnlistOther)
    {
      poList = poReflnlistOther;
    }
  else
    {
      poList = m_poReflnlistIn;
    }

  m_nFI0 = poList->m_nFI_fObsPx0;
  m_nFI1 = poList->m_nFI_fObsPx1;
  m_nFI2 = poList->m_nFI_fObsRotMid;
  m_nFI3 = poList->m_nFI_fObsRotWidth;
  m_nFI4 = poList->m_nFI_fObsRotSigma;

  if (poList->m_nFI_fCalcMosCoeffA>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcMosCoeffA);
  if (poList->m_nFI_fCalcMosCoeffB>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcMosCoeffB);


  if (0 > m_nFI0)
    {
      (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcPx0);
    }
  if (0 > m_nFI1)
    {
      (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcPx1);
    }
  if (0 > m_nFI2)
    {
      (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcRotMid);
    }
  if (0 > m_nFI3)
  {
      (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfCalcRotWidth);
  }

  if (0 < m_nFI4) 
  {
      (void) m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsRotSigma);
      
  };

  m_nFIMM0 = m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsXmm);
  m_nFIMM1 = m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsYmm);
  m_nFIMM2 = m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfObsZmm);

  // Then copy the list to the new one.
  
  if (poList->m_nFI_fGonio1>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfGonio1);
  if (poList->m_nFI_fGonio2>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfGonio2);
  if (poList->m_nFI_fGonio3>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfGonio3);
  if (poList->m_nFI_fGonio4>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfGonio4);
  if (poList->m_nFI_fGonio5>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_sfGonio5);
  if (poList->m_nFI_nGonioRotAxis>=0)
      m_poReflnlist->nExpandGetField(m_poReflnlist->ms_snGonioRotAxis);

  m_poReflnlist->nInsertListFrom(*poList);

  

  // Now use observed centroids or calculated centroids

  if (0 > m_nFI0)
    {
      m_nFI0 = m_poReflnlist->m_nFI_fCalcPx0;
    }
  else
    {
      m_nFI0 = m_poReflnlist->m_nFI_fObsPx0;
    }
  if (0 > m_nFI1)
    {
      m_nFI1 = m_poReflnlist->m_nFI_fCalcPx1;
    }
  else
    {
      m_nFI1 = m_poReflnlist->m_nFI_fObsPx1;
    }
  if (0 > m_nFI2)
    {
      m_nFI2 = m_poReflnlist->m_nFI_fCalcRotMid;
    }
  else
    {
      m_nFI2 = m_poReflnlist->m_nFI_fObsRotMid;
    }
  if (0 > m_nFI3)
    {
      m_nFI3 = m_poReflnlist->m_nFI_fCalcRotWidth;
    }
  else
    {
      m_nFI3 = m_poReflnlist->m_nFI_fObsRotWidth;
    }

  m_nFI4 = m_poReflnlist->m_nFI_fObsRotSigma;

  
  


  if ( (0 > m_nFI0) || (0 > m_nFI1) || (0 > m_nFI2) )
    {
      cout << "ERROR, missing spot centroids from reflection list!\n";
      return (-1);
    }


  nRef  = m_poReflnlist->nGetNumReflns();
  for (i = 0; i < nRef; i++)
  {

      // In this loop over reflections, we should probably reject
      //  reflections outside resolution and I/sigma cutoffs, too.
      //  (or use dtreflnmerge to filter the list beforehand)
      // Get detector number
      
      nDetNum = m_poReflnlist->poGetRefln(i)->nGetField(m_poReflnlist->m_nFI_nDetNum);
      
      // Get observed pixel coordinates
      
      f1PX = m_poReflnlist->poGetRefln(i)->fGetField(m_nFI0);
      f2PX = m_poReflnlist->poGetRefln(i)->fGetField(m_nFI1);
      
      // Calculate observed millmeter coordinates
      
      nStat = m_ppoDetector[nDetNum]->m_poSpatial->nPixeltoMM(f1PX, 
                                                              f2PX,
                                                              &f1MM,
                                                              &f2MM,
                                                              &f3MM);
  
      if (0 != nStat)
      {
          cout << "WARNING in Crefine::nCopyReflns, invalid "
              << "reflection coordinates!\n";
          (*m_poReflnlist)[i].vSetIntensity(-9999.0);  // Effectively delete from list.
          
          // tjn:  Don't delete the reflection:  we can benefit from the assumption that
          // the m_poReflnlist and m_poReflnlistIn objects have corresponding fields.
          //m_poReflnlist->nDelete(i); 
          //nRef--;
          //i--;
      }
      else
      {
          m_poReflnlist->poGetRefln(i)->vSetField(m_poReflnlist->m_nFI_fObsXmm,
              f1MM);
          m_poReflnlist->poGetRefln(i)->vSetField(m_poReflnlist->m_nFI_fObsYmm,
              f2MM);
          m_poReflnlist->poGetRefln(i)->vSetField(m_poReflnlist->m_nFI_fObsZmm,
              f3MM);
      }
  }

   return (0);
}

int Crefine::nRoundIndices(double fMatrix[3][3],// Crystal UB matrix (CB)
                           const double a3fR0[3],     // Input float indices (starting position)
                           const double a3fR1[3],     // Input float indices (ending position)
                           double& fLowest,           // Current lowest residual.
                           double a3fLowestResid[3],  // Lowest residual vector.
                           int   *pnHKL)             // Returned nearest indices
{
    

    /*  This function was modified from a previous version of nRoundIndices() so that it now checks all HKL values
        between a3fRO[] and a3fR1[] (instead of rounding a single a3fR[] vector).  Typically, these two arrays will represent 
        the floating HKL values for the reflection  at it's start and end rotation points.  
        We want to search the line between these two points for an HKL
        point that is nearby.  Thus, for wide sliced images, the algorithm will choose the best reflection.

        The algorithm works by searching a box of reflections (chosen by rounding the extrema a3fRO and a3fR1 up and down accordingly).


    */
    //   Rounds non-integral indices 'a3fR' to integral indices 'pnHKL[0]...[2]'
    //   consistent with absences for space group number.
    //  'fMatrix' is matrix to convert from HKL to reciprocal lattice
    //   coordinates.
    
    int   a3x2nLimits[3][2];   // Min, max limits on H, K, L to test
    int   nH,nK,nL,nx;
    static int nEntry = 0;
    double a3fR01[3];
    double a3fDelta[3];
    double a3fTempVec[3];
    double a3fDiff[3];
    double fR01Length;
    double fDist,f0;
    
    // Round indices to nearest integers:

    for (nx = 0; nx<3; nx++) {
        a3x2nLimits[nx][0] = (int) floor( 0.5 + min(a3fR0[nx],a3fR1[nx]));
        a3x2nLimits[nx][1] = (int) floor( 0.5 + max(a3fR0[nx],a3fR1[nx]));
    };
    vSubVec3DVec3D(a3fR1,a3fR0,a3fR01);
    fR01Length = fNormVec3D(a3fR01);
    nEntry++;
  


    // Find indices of closest possible lattice point to original indices 'a3fR'    
    
    int nStat = -1;
    for (nH = a3x2nLimits[0][0]; nH <= a3x2nLimits[0][1]; nH++)
    {
        a3fDelta[0] = nH - a3fR0[0];
        for (nK = a3x2nLimits[1][0]; nK <= a3x2nLimits[1][1]; nK++)
        {
            a3fDelta[1] = nK - a3fR0[1];
            for (nL = a3x2nLimits[2][0]; nL <= a3x2nLimits[2][1]; nL++)
            {
                a3fDelta[2] = nL - a3fR0[2];

                
                // Test systematic absence
                
                if (   ( (0 != nH) || (0 != nK) || (0 != nL) )
                    && (!m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->bIsExtinct(nH, nK, nL)) ) {
                    
                    // Calculate reciprocal lattice distance to the line between a3fR0 and a3fR1

                    // Project a3fDelta down to the line between a3fR0 and a3fR1.
                    // (a3fR01 represents a normalized direction vector for this line).
                    vMulVec3DScalar(a3fR01,f0 = fDot3D(a3fDelta,a3fR01),a3fTempVec);

                    // Acceptable solution points should project onto the line between a3fR0 and a3fR1
                    // But for now, we are not checking this, since it might not necc. hold.

                    vSubVec3DVec3D(a3fDelta,a3fTempVec,a3fDiff);

                    // a3fDiff[] now contains the HKL difference vector for spot indexed by nH nK nL
                                        
                    if ((f0 = fLenVec3D(a3fDiff))<0.2) {
                        vMulMat3DVec3D(fMatrix, a3fDiff, a3fTempVec);
                        fDist = fLenVec3D(a3fTempVec);
                        
                        if (fDist < fLowest) {
                            fLowest = fDist;
                            vCopyVec3D(a3fTempVec,a3fLowestResid);
                            pnHKL[0] = nH;
                            pnHKL[1] = nK;
                            pnHKL[2] = nL;
                            nStat    = 0;
                        }
                    };
                }
            }  // end l loop
        }  // end k loop
    }  // end h loop

    return (nStat);
}

int Crefine::nSetup(void)
{
  // Setup before refinement loops

  int i;      // Loop counter

  //
  for (i = 0; i < (nMAXPARAMS); i++)
    {
      // Make sure all parameters are fixed

      m_anRefineFlags[i] = 0;
    }

  // Transfer user refinement flags to internal ones

  for (i = 0; i < (nMAXPARAMS); i++)
    {
      m_anRefineFlags[i] = m_nUserRefineFlags[i];
    }

  // Change some crystal refinement flags from user to their internal ones:
  // However:  Don't do so unless we are refining crystal constraints.

  if (!m_ppoCrystals[m_nTwinID-1]->bIsAvailable())
    {
      cout << "Crystal is not available!\n" << endl;
      return (-1);
    }

  if (m_anRefineFlags[ms_nCrysConstraints]) {
      
      if (  (m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass() == eSpacegroup_trigonal_class)
          ||(m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass() == eSpacegroup_hexagonal_class)
          ||(m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass() == eSpacegroup_tetragonal_class))
      {
          // a* = b*, so c* slot is moved to b* slot
          
          m_anRefineFlags[ms_nCrysAstar] = 
              m_anRefineFlags[ms_nCrysAstar] ||
              m_anRefineFlags[ms_nCrysBstar];
          m_anRefineFlags[ms_nCrysBstar] = m_anRefineFlags[ms_nCrysCstar];
          m_anRefineFlags[ms_nCrysCstar] = 0;
          m_anRefineFlags[ms_nCrysAlps] = 0;  // Do not refine angles
          m_anRefineFlags[ms_nCrysBets] = 0;
          m_anRefineFlags[ms_nCrysGams] = 0;
          
      }
      else if (eSpacegroup_orthorhombic_class
          == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass())
      {
          m_anRefineFlags[ms_nCrysAlps] = 0;  // Do not refine angles
          m_anRefineFlags[ms_nCrysBets] = 0;
          m_anRefineFlags[ms_nCrysGams] = 0;
          
      }
      else if (eSpacegroup_rhombohedral_class
          == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass())
      {
          // a* = b* = c*, so alp* slot is moved to b*
          m_anRefineFlags[ms_nCrysAstar] = 
              m_anRefineFlags[ms_nCrysAstar] ||
              m_anRefineFlags[ms_nCrysBstar] ||
              m_anRefineFlags[ms_nCrysCstar];
          m_anRefineFlags[ms_nCrysBstar] = 
              m_anRefineFlags[ms_nCrysAlps] ||
              m_anRefineFlags[ms_nCrysBets] ||
              m_anRefineFlags[ms_nCrysGams];
          m_anRefineFlags[ms_nCrysCstar] = 0;
          m_anRefineFlags[ms_nCrysAlps] = 0;
          m_anRefineFlags[ms_nCrysBets] = 0;
          m_anRefineFlags[ms_nCrysGams] = 0;
          
      }
      else if (eSpacegroup_cubic_class
          == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass())
      {
          m_anRefineFlags[ms_nCrysAstar] = 
              m_anRefineFlags[ms_nCrysAstar] ||
              m_anRefineFlags[ms_nCrysBstar] ||
              m_anRefineFlags[ms_nCrysCstar];
          m_anRefineFlags[ms_nCrysBstar] = 0;
          m_anRefineFlags[ms_nCrysCstar] = 0;
          m_anRefineFlags[ms_nCrysAlps] = 0;
          m_anRefineFlags[ms_nCrysBets] = 0;
          m_anRefineFlags[ms_nCrysGams] = 0;
      }
      else if (eSpacegroup_monoclinic_class
          == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass())
      {
          // alp* = 90, so bet* slot is moved to alp*
          
          m_anRefineFlags[ms_nCrysAlps] = 0;
          if (1 == m_anRefineFlags[ms_nCrysBets])
          {
              m_anRefineFlags[ms_nCrysAlps] = 1;
          }
          m_anRefineFlags[ms_nCrysBets] = 0;
          m_anRefineFlags[ms_nCrysGams] = 0;
          
      }
  };

  // If this is a RAPID, then don't refine the third translation axis (=='s the rapid radius)
  // Or the second translation (which is dependent upon the curvature of the plate).

  if (m_ppoDetector[0]->m_poSpatial->eSpatialType() == eSpatial_cylindrical_geometry) {
      m_anRefineFlags[ms_nDetTrans3] = 0;
//      m_anRefineFlags[ms_nDetRot2] = 0;
//      m_anRefineFlags[ms_nSourceRot2] = 0;
      m_anRefineFlags[ms_nDetTrans2] = 0;
//      m_anRefineFlags[ms_nDetTrans1] = 0;
  };

  // Disable some variables if this is not twinned.
  for (i = m_nMaxTwinLaws + 1; i <= g_nMaxTwinLaws; i++) {
      m_anRefineFlags[ms_nCrysTwinLaw10 + i*3 - 3] = 0;
      m_anRefineFlags[ms_nCrysTwinLaw11 + i*3 - 3] = 0;
      m_anRefineFlags[ms_nCrysTwinLaw12 + i*3 - 3] = 0;
      m_anRefineFlags[ms_nCrysRecip00 + i*3] = 0;
      m_anRefineFlags[ms_nCrysRecip01 + i*3] = 0;
      m_anRefineFlags[ms_nCrysRecip02 + i*3] = 0;
  };

  
  
  // These are disabled for now, until they get supported.
  m_anRefineFlags[ms_nGonioMiss1] = 0;
  m_anRefineFlags[ms_nGonioMiss2] = 0;


  return (0);
}

void Crefine::vPrintResidualLog(int nPrintCode) {
    
          
  if ( (0 == nPrintCode) && (0 >= m_sLogRefine.length()) )
    {
      printf("\nNo refinement iterations to list.\n");
      printf("Perhaps your input file has invalid CRYSTAL information?\n");
      printf("Did you forget to input the previous dtrefine or dtindex (.head) results?\n\n");
      return;
    }
  if ((nPrintCode==0) || (nPrintCode==1))
    {
      printf("\n");    
      printf("Listing for all refinement iterations\n");
      printf("-----------------------------------------------------------------\n");
      printf("Iter  Twin   Volume    Rms     Rms  Total Ignored Rejects Rejects\n");
      printf("                      (MM)   (DEG)  refln   refln    (MM)   (DEG)\n");
      printf("-----------------------------------------------------------------\n");
    }
  if (nPrintCode==0) 
    {
      printf("%s",m_sLogRefine.string());
    } 
  else
    {
      printf("%s",m_sLogRefineLast.string());
    }
    
  if ((nPrintCode==0) || (nPrintCode==3)) 
    {
      printf("-----------------------------------------------------------------\n");
      printf("\n\n");
    }
}

void Crefine::vLogResiduals(int nNumIter) {
    char pcTemp[1000];    
    char pcTwinString[20];
    sprintf(pcTwinString,"%2.2d/%2.2d",m_nTwinID,m_nNumTwins);
    sprintf(pcTemp,"%4d %4s %8.0lf %6.4lf %7.4lf %6d %7d %7d %7d\n",
        nNumIter,pcTwinString,(double) m_ppoCrystals[m_nTwinID-1]->fCalcVolume(),
        m_fRMSMM,m_fRMSDEG,        
        m_a7nReflnNumbers[eReflectionNumberAccept],
        m_a7nReflnNumbers[eReflectionNumberIgnored],
        m_a7nReflnNumbers[eReflectionNumberMMRejects],
        m_a7nReflnNumbers[eReflectionNumberDegRejects]
        );
    m_sLogRefineLast = pcTemp;
    m_sLogRefine += pcTemp;
};



int Crefine::nListResults(void)
{
  // List some results

  int i;
  int nShiftStat = 0;
  Cstring sTemp;
  char a100cTemp[100];
  char *acFixed = "   fixed ";
  a100cTemp[9] = '\0';
  static Cstring sLine = "==================================================================";
  
  if (1 < m_nVerbose) {
      int nAxis;
      int nPlot;
      double fTotal;
      double f0;
      float* pfRefs;
      char a3cHKL[3] = { 'h','k','l'};
      for (nPlot = 0; nPlot < 2; nPlot++) {
          if (nPlot == 0)
              printf("Percentage of total reflections vs HKL integer residual\n");
          else
              printf("Percentage of total intensity vs HKL integer residual\n");
          printf("============================================================================\n");           
          printf("%2s <=0.05   0.10   0.15   0.20   0.25   0.30   0.35   0.40   0.45  >0.50\n","");
          for (nAxis = 0; nAxis < 3; nAxis++) {
              if (nPlot == 0)
                  pfRefs = &m_a3x10fIndexHKLErrorRefs[nAxis][0];
              else
                  pfRefs = &m_a3x10fIndexHKLErrorIntensity[nAxis][0];
              // Sum the total.
              fTotal = 0.0;
              for (i = 0; i < 10 ;i++)
                  fTotal += pfRefs[i];
              printf("%c ",a3cHKL[nAxis],"");
              for (i = 0;i < 10; i++) {
                  f0 = pfRefs[i]/max(1.0,fTotal);
                  if (f0 == 0.0)
                      printf("%7s","---");
                  else
                      printf("%7.2lf",100.0*f0);
              };              
              printf("\n");
          };
          printf("============================================================================\n");           
      };
  };

  cout << "\nRefinement results\n" << flush;
  printf("%s", sLine.string());
  printf("\nCrystal\n                         a, b, c:  %9.4f  %9.4f  %9.4f\n",
         m_a6fCrysUnitCell[0], m_a6fCrysUnitCell[1], m_a6fCrysUnitCell[2]);

  printf("                          Sigmas:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysAstar])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysUnitCellSigs[0]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysBstar])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysUnitCellSigs[1]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysCstar])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysUnitCellSigs[2]);
  printf("%s", a100cTemp);
  printf("\n");

  //  printf("%9.4f  %9.4f  %9.4f\n",
  // m_a6fCrysUnitCellSigs[0], m_a6fCrysUnitCellSigs[1], m_a6fCrysUnitCellSigs[2]);

  printf("                          Shifts:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysAstar])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysTurds[0]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysBstar])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysTurds[1]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysCstar])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysTurds[2]);
  printf("%s", a100cTemp);
  printf("\n");

  //  printf("                          Shifts:  %9.4f  %9.4f  %9.4f\n",
  //         m_a6fCrysTurds[0], m_a6fCrysTurds[1], m_a6fCrysTurds[2]);

  printf("\n              alpha, beta, gamma:  %9.4f  %9.4f  %9.4f\n",
         m_a6fCrysUnitCell[3], m_a6fCrysUnitCell[4], m_a6fCrysUnitCell[5]);

  printf("                          Sigmas:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysAlps])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysUnitCellSigs[3]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysBets])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysUnitCellSigs[4]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysGams])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysUnitCellSigs[5]);
  printf("%s", a100cTemp);
  printf("\n");

  //  printf("                          Sigmas:  %9.4f  %9.4f  %9.4f\n",
  // m_a6fCrysUnitCellSigs[3], m_a6fCrysUnitCellSigs[4], m_a6fCrysUnitCellSigs[5]);

  printf("                          Shifts:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysAlps])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysTurds[3]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysBets])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysTurds[4]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysGams])
    sprintf(a100cTemp, "%9.4f", m_a6fCrysTurds[5]);
  printf("%s", a100cTemp);
  printf("\n");

  //  printf("                          Shifts:  %9.4f  %9.4f  %9.4f\n",
  // m_a6fCrysTurds[3], m_a6fCrysTurds[4], m_a6fCrysTurds[5]);

  printf("\n           Crys Rot1, Rot2, Rot3: %9.3f  %9.3f  %9.3f\n",
         m_a3fCrysRot[0], m_a3fCrysRot[1], m_a3fCrysRot[2]);

  printf("                          Sigmas:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysRot1])
    sprintf(a100cTemp, "%9.4f", m_a3fCrysRotSigs[0]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysRot2])
    sprintf(a100cTemp, "%9.4f", m_a3fCrysRotSigs[1]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysRot3])
    sprintf(a100cTemp, "%9.4f", m_a3fCrysRotSigs[2]);
  printf("%s", a100cTemp);
  printf("\n");

  //  printf("                          Sigmas:  %9.4f  %9.4f  %9.4f\n",
  //         m_a3fCrysRotSigs[0], m_a3fCrysRotSigs[1], m_a3fCrysRotSigs[2]);

  printf("                          Shifts:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysRot1])
    sprintf(a100cTemp, "%9.4f", m_a3fCrysRotShifts[0]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysRot2])
    sprintf(a100cTemp, "%9.4f", m_a3fCrysRotShifts[1]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nCrysRot3])
    sprintf(a100cTemp, "%9.4f", m_a3fCrysRotShifts[2]);
  printf("%s", a100cTemp);
  printf("\n");

  //  printf("                          Shifts:  %9.4f  %9.4f  %9.4f\n",
  //         m_a3fCrysRotShifts[0], m_a3fCrysRotShifts[1], m_a3fCrysRotShifts[2]);


    printf("\n                       Mosaicity: ");
  if (   (1 == m_nUserRefineFlags[ms_nCrysMosaicity])
      && (0.0 >= m_a3fCrysMosaicitySig[0]) )
    {
      printf(" %9.4f Fixed.\n", m_a3fCrysMosaicity[0]);
    }
  else
    {
      printf(" %9.4f", m_a3fCrysMosaicity[0]);
      printf("\n");
    }

  printf("                           Sigma:  ");
  strcpy(a100cTemp, acFixed);
  if (   (1 == m_nUserRefineFlags[ms_nCrysMosaicity])
      && (0.0 < m_a3fCrysMosaicitySig[0]) )
    sprintf(a100cTemp, "%9.4f ",  m_a3fCrysMosaicitySig[0]);
  printf("%s\n", a100cTemp);


  printf("                           Shift:  ");
  strcpy(a100cTemp, acFixed);
  if (   (1 == m_nUserRefineFlags[ms_nCrysMosaicity])
      && (0.0 < m_a3fCrysMosaicitySig[0]) )
    sprintf(a100cTemp, "%9.4f ", m_a3fCrysMosaicityShift[0]);
  printf("%s\n", a100cTemp);


  printf("%s", sLine.string());
  for (i = 0; i < m_nNumDetectors; i++)
    {
      printf("\nDetector:  %d\n", i);
      sTemp  = "DetTrans: ";
      sTemp += m_ppoDetector[i]->m_poGoniometer->sGetName(3);
      sTemp += ", " + m_ppoDetector[i]->m_poGoniometer->sGetName(4);
      sTemp += ", " + m_ppoDetector[i]->m_poGoniometer->sGetName(5);

      if (32 < sTemp.length())
        {
          sTemp = sTemp.substr(0, 32);
        }
      else if (32 > sTemp.length())
        {
          while (32 != sTemp.length())
            sTemp = Cstring(' ') + sTemp;
        }
      printf("%s: %9.3f  %9.3f  %9.3f\n", sTemp.string(),
             m_a10x3fDetTrans[i][0], m_a10x3fDetTrans[i][1], m_a10x3fDetTrans[i][2]);

      //printf("                          Sigmas:  %9.4f  %9.4f  %9.4f\n",
      //     m_a10x3fDetTransSigs[i][0], m_a10x3fDetTransSigs[i][1], m_a10x3fDetTransSigs[i][2]);
      printf("                          Sigmas:  ");
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetTrans1])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetTransSigs[i][0]);
      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetTrans2])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetTransSigs[i][1]);
      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetTrans3])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetTransSigs[i][2]);
      printf("%s\n", a100cTemp);

      //      printf("                Shifts:  %9.4f  %9.4f  %9.4f\n",
      //     m_a10x3fDetTransShifts[i][0], m_a10x3fDetTransShifts[i][1], m_a10x3fDetTransShifts[i][2]);
      printf("                          Shifts:  ");
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetTrans1])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetTransShifts[i][0]);
      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetTrans2])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetTransShifts[i][1]);

      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetTrans3])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetTransShifts[i][2]);
      printf("%s\n", a100cTemp);

      sTemp = "DetRots: ";
      sTemp += m_ppoDetector[i]->m_poGoniometer->sGetName(0);
      sTemp += ", " + m_ppoDetector[i]->m_poGoniometer->sGetName(1);
      sTemp += ", " + m_ppoDetector[i]->m_poGoniometer->sGetName(2);

      if (32 < sTemp.length())
        sTemp = sTemp.substr(0, 32);
      else if (32 > sTemp.length())
        {
          while (32 != sTemp.length())
            sTemp = Cstring(' ') + sTemp;
        }

      printf("\n%s: %9.3f  %9.3f  %9.3f\n", sTemp.string(),
             m_a10x3fDetRot[i][0], m_a10x3fDetRot[i][1], m_a10x3fDetRot[i][2]);
      //      printf("                Sigmas:  %9.4f  %9.4f  %9.4f\n",
      //     m_a10x3fDetRotSigs[i][0], m_a10x3fDetRotSigs[i][1], m_a10x3fDetRotSigs[i][2]);
      printf("                          Sigmas:  ");
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetRot1])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetRotSigs[i][0]);
      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetRot2])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetRotSigs[i][1]);
      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetRot3])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetRotSigs[i][2]);
      printf("%s\n", a100cTemp);

      //printf("                          Shifts:  %9.4f  %9.4f  %9.4f\n",
      //     m_a10x3fDetRotShifts[i][0], m_a10x3fDetRotShifts[i][1], m_a10x3fDetRotShifts[i][2]);
      printf("                          Shifts:  ");
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetRot1])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetRotShifts[i][0]);
      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetRot2])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetRotShifts[i][1]);
      printf("%s  ", a100cTemp);
      strcpy(a100cTemp, acFixed);
      if (1 == m_nUserRefineFlags[ms_nDetRot3])
        sprintf(a100cTemp, "%9.4f", m_a10x3fDetRotShifts[i][2]);
      printf("%s\n", a100cTemp);
    }

  printf("%s", sLine.string());
  printf("\nSource\n");
  printf("          Wavelength, Rot1, Rot2:  %11.6f %7.3f  %9.3f\n",
         m_fSourceWave, m_a2fSourceRot[0], m_a2fSourceRot[1]);
  //  printf("                Sigmas:  %9.4f  %9.4f  %9.4f\n",
  // m_fSourceWaveSig, m_a2fSourceRotSigs[0], m_a2fSourceRotSigs[1]);

  printf("                          Sigmas:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nSourceWave])
    sprintf(a100cTemp, "%11.6f", m_fSourceWaveSig);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nSourceRot1])
    sprintf(a100cTemp, "%9.4f", m_a2fSourceRotSigs[0]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nSourceRot2])
    sprintf(a100cTemp, "%9.4f", m_a2fSourceRotSigs[1]);
  printf("%s\n", a100cTemp);

  //  printf("                          Shifts:  %9.4f  %9.4f  %9.4f\n",
  // m_a2fSourceRotShifts[0], m_a2fSourceRotShifts[1], m_fSourceWaveShift);
  printf("                          Shifts:  ");
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nSourceWave])
    sprintf(a100cTemp, "%11.6f", m_fSourceWaveShift);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nSourceRot1])
    sprintf(a100cTemp, "%9.4f", m_a2fSourceRotShifts[0]);
  printf("%s  ", a100cTemp);
  strcpy(a100cTemp, acFixed);
  if (1 == m_nUserRefineFlags[ms_nSourceRot2])
    sprintf(a100cTemp, "%9.4f", m_a2fSourceRotShifts[1]);
  printf("%s\n", a100cTemp);
  printf("%s", sLine.string());
  printf("\nRefinement residuals "
         "\nrmsResid (A-1) = %10.5f"
         "\nrmsResid (mm)  = %9.4f", m_fRMSANG, m_fRMSMM);
  float fPxSize
      = m_ppoDetector[0]->m_poSpatial->fGetNominalPixelSize(0);
  if (m_fRMSMM <= 0.0)
    printf("  WARNING: TOO GOOD! Perhaps there are no reflections!");
  else if (m_fRMSMM <= (1.0 * fPxSize))
    printf("  EXCELLENT, less than or equal to 1 pixel (%.3f mm).", fPxSize);
  else if (m_fRMSMM <= (2.0 * fPxSize))
    printf("  GOOD, less than or equal to 2 pixels (%.3f mm).", 2.0 * fPxSize);
  else if (m_fRMSMM <= (3.0 * fPxSize))
    printf("  FAIR, less than or equal to 3 pixels (%.3f mm).", 3.0 * fPxSize);
  else
    printf("  NOT GOOD, more than 3 pixels (%.3f mm).", 3.0 * fPxSize);

  printf("\nrmsResid (Deg) = %9.4f", m_fRMSDEG);
  if (m_a3fCrysMosaicity[0] < (0.667 * m_poRotation->fGetIncrement()) )
    {
      if (m_fRMSDEG <= 0.0)
        printf("  WARNING: TOO GOOD! Perhaps there are no reflections!\n");
      else if (m_fRMSDEG <= (0.333 * m_poRotation->fGetIncrement()))
        printf("  EXCELLENT, less than 1/3rd the image rot width.\n");
      else if (m_fRMSDEG <= (0.667 * m_poRotation->fGetIncrement()))
        printf("  GOOD, less than 2/3rds the image rot width.\n");
      else
        printf("  NOT GOOD, more than 2/3rds the image rot width.\n");
    }
  else
    {
      if (m_fRMSDEG <= 0.0)
        printf("  WARNING: TOO GOOD! Perhaps there are no reflections!\n");
      else if (m_fRMSDEG <= (0.25 * (m_a3fCrysMosaicity[0] ) ))
        printf("  EXCELLENT, less than 1/4 the mosaicity.\n");
      else if (m_fRMSDEG <= (0.5 * (m_a3fCrysMosaicity[0]) ))
        printf("  GOOD, less than 1/2 the mosaicity.\n");
      else
        printf("  NOT GOOD, more than 1/2 the mosaicity.\n");
    }

  printf("%s\n", sLine.string());
  printf("\n");

  if ((0.0 < m_a3fCrysMosaicitySig[0]) && (m_bPrintDetectorMosaicity))
  {
    printf("================================\n");
  };

  fflush(stdout);
  return (0);
}

int Crefine::nUpdateHeader(Cimage_header *poHeader, int nTwinID,const Cstring &sPre)
{
    int nx,ny;
    float a6fSigmas[6];
    float afBuf[100];
    int   anBuf[700];
    // Place the current values of the crystal, detector and source
    // properties in the given header

    
    m_ppoCrystals[nTwinID-1]->nUpdateHeader(poHeader);
        m_ppoCrystals[nTwinID-1]->nUpdateHeaderSigmaValues(poHeader,
                &m_a6fCrysUnitCellSigs[0],
                &m_a3fCrysRotSigs[0],
                & m_a3fCrysMosaicitySig[0],
                & m_aa3fCrysRecipShiftSig[0][0],
                & m_aa3fCrysTwinLawSig[1][0],
                & m_afCrysTwinFractionSig[1]);

    
    // Which detector(s)?
    
    m_ppoDetector[0]->nUpdateHeader(poHeader, sPre);
    m_poSource->nUpdateHeader(poHeader);
    m_poRotation->nUpdateHeader(poHeader);
    m_poCrysGonio->nSaveOffsets(poHeader);
    
    for (nx=0;nx<m_nNumTwins;nx++) {
        afBuf[nx] = m_poComp[nx].m_fRMSMM;
        anBuf[nx] = 1;
    };
    poHeader->nReplaceValue(Crefine::ms_sDtrefineRmsMm,m_nNumTwins, afBuf, 3);
    poHeader->nReplaceValue(Crefine::ms_sDtrefineTwinsRefined,m_nNumTwins,anBuf);
    for (nx=0;nx<m_nNumTwins;nx++)
        afBuf[nx] = m_poComp[nx].m_fRMSDEG;
    poHeader->nReplaceValue(Crefine::ms_sDtrefineRmsDeg,m_nNumTwins, afBuf, 3);
    for (nx=0;nx<m_nNumTwins;nx++) {
        for (ny=0;ny<7;ny++) {
            anBuf[nx*7+ny] = m_poComp[nx].m_a7nReflnNumbers[ny];
            afBuf[nx*7+ny] = m_poComp[nx].m_a7fReflnNumbers[ny];
        };
    };
    poHeader->nReplaceValue(Crefine::ms_sDtrefineReflectionNumbers,7*m_nNumTwins,anBuf);
    poHeader->nReplaceValue(Crefine::ms_sDtrefineReflectionNumbersPhotons,7*m_nNumTwins,afBuf);

    // Save the profile information.
    
    poHeader->nReplaceValue(Crefine::ms_sDtrefineLambda,m_fLambda);
    poHeader->nReplaceValue(Crefine::ms_sDtrefineProfileCount,m_nProfileCt);
    if (m_nProfileCt)
    {
        poHeader->nReplaceValue(Crefine::ms_sDtrefineResidProfile,m_nProfileCt,m_pfResidProfile);
        poHeader->nReplaceValue(Crefine::ms_sDtrefineLambdaProfile,m_nProfileCt,m_pfLambdaProfile);
        
    }
    Cstring sSigmas = "_SIGMA";
    
    
    a6fSigmas[0] = m_a10x3fDetRotSigs[0][0];
    a6fSigmas[1] = m_a10x3fDetRotSigs[0][1];
    a6fSigmas[2] = m_a10x3fDetRotSigs[0][2];
    a6fSigmas[3] = m_a10x3fDetTransSigs[0][0];
    a6fSigmas[4] = m_a10x3fDetTransSigs[0][1];
    a6fSigmas[5] = m_a10x3fDetTransSigs[0][2];
    Cstring sPrefix = sPre;
    if ("" == sPrefix)
    {
        if (NULL == m_psDetectorNames)
        {
            // Get the number of detectors
            
            int nTemp;
            nTemp = poHeader->nGetValue(Cdetector::ms_sDetectorNumber,
                &m_nNumDetectors);
            if (0 == nTemp)
            {
                // Get the detector names
                
                m_psDetectorNames = new Cstring [m_nNumDetectors];
                m_bNewDetNames = TRUE;
                nTemp = poHeader->nGetValue(Cdetector::ms_sDetectorNames,
                    m_nNumDetectors,
                    m_psDetectorNames);
            }
            sPrefix = m_psDetectorNames[0];
        }
    }
    poHeader->nReplaceValue(sPrefix + Cgoniometer::ms_sGonioValues + sSigmas, 6,
        a6fSigmas, 4);
    
    poHeader->nReplaceValue(Csource::ms_sSourceValues + sSigmas, 2,
        m_a2fSourceRotSigs, 4);
    poHeader->nReplaceValue(Csource::ms_sSourcePrefix
        + Cwavelength::ms_sWavelength + sSigmas,
        m_fSourceWaveSig, 6);
    poHeader->nReplaceValue(Crefine::ms_sDtrefineResoUsed,2,&m_a2fResoUsed[0]);

  
    return (0);
}

int
Crefine::nSetRefineFlags(const Cstring& sFlags)
{
  // Look in the string sFlags and set the parameters to be refined based
  // on what is found in sFlags.  Work from left to right.
  // +DetTrans1   refine detector translation 1
  // -DetTrans1   fix    detector translation 1
  // ...

  int     i;
  int     nRefined;
  int     nTemp;
  Cstring sLower;
  Cstring sTemp;

  sLower = sFlags;

  // Strip off leading white space

  while (   (' '  == sLower.GetAt(0))
         || ('\n' == sLower.GetAt(0))
         || ('\t' == sLower.GetAt(0)))
    sLower = sLower.after(0);

  if ('-' == sLower.GetAt(0))
    nRefined = 0;
  else if ('+' == sLower.GetAt(0))
    nRefined = 1;
  else
    return (1); // Error

  // Get the first word in sTemp

  nTemp = 0;
  for (i = 1; ( (i < sLower.length()) && (0 == nTemp)); i++)
    {
      if (   (' '  == sLower.GetAt(i))
          || ('\n' == sLower.GetAt(i))
          || ('\t' == sLower.GetAt(i)))
        nTemp = i-1;
    }
  if (0 == nTemp) nTemp = sLower.length()-1;
  sTemp  = sLower.substr(1, nTemp);

  // sTemp contains a single flag


  if (  (ms_sDetTrans1 == sTemp) || (m_a6sUserDetNames[3] == sTemp) )
    {
      m_nUserRefineFlags[ms_nDetTrans1]  = nRefined;
    }
  else if ( (ms_sDetTrans2 == sTemp) || (m_a6sUserDetNames[4] == sTemp) )
    {
      m_nUserRefineFlags[ms_nDetTrans2]  = nRefined;
    }
  else if ( (ms_sDetTrans3 == sTemp) || (m_a6sUserDetNames[5] == sTemp) )
    {
      m_nUserRefineFlags[ms_nDetTrans3]  = nRefined;
    }
  else if ( (ms_sDetRot1 == sTemp) || (m_a6sUserDetNames[0] == sTemp) )
    {
      m_nUserRefineFlags[ms_nDetRot1]    = nRefined;
    }
  else if ( (ms_sDetRot2 == sTemp) || (m_a6sUserDetNames[1] == sTemp) )
    {
      m_nUserRefineFlags[ms_nDetRot2]    = nRefined;
    }
  else if ( (ms_sDetRot3 == sTemp) || (m_a6sUserDetNames[2] == sTemp) )
    {
      m_nUserRefineFlags[ms_nDetRot3]    = nRefined;
    }
  else if (ms_sGonioMiss1 == sTemp)
    {
      m_nUserRefineFlags[ms_nGonioMiss1]    = nRefined;
    }
  else if (ms_sGonioMiss2 == sTemp)
    {
      m_nUserRefineFlags[ms_nGonioMiss2]    = nRefined;
    }
  else if (ms_sGonioMiss == sTemp)
    {
      m_nUserRefineFlags[ms_nGonioMiss1]    = nRefined;
      m_nUserRefineFlags[ms_nGonioMiss2]    = nRefined;
    }
  else if (ms_sSourceRot1 == sTemp)
    {
      m_nUserRefineFlags[ms_nSourceRot1] = nRefined;
    }
  else if (ms_sSourceRot2 == sTemp)
    {
      m_nUserRefineFlags[ms_nSourceRot2] = nRefined;
    }
  else if (ms_sSourceRot == sTemp)
    {
      m_nUserRefineFlags[ms_nSourceRot1] = nRefined;
      m_nUserRefineFlags[ms_nSourceRot2] = nRefined;
    }
  else if (ms_sSourceWave == sTemp)
    {
      m_nUserRefineFlags[ms_nSourceWave] = nRefined;
    }
  else if (ms_sCrysRot1 == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysRot1]   = nRefined;
    }
  else if (ms_sCrysRot2 == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysRot2]   = nRefined;
    }
  else if (ms_sCrysRot3 == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysRot3]   = nRefined;
    }
  else if (ms_sCrysAstar == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysAstar]  = nRefined;
    }
  else if (ms_sCrysBstar == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysBstar]  = nRefined;
    }
  else if (ms_sCrysCstar == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysCstar]  = nRefined;
    }
  else if (ms_sCrysAlps == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysAlps]   = nRefined;
    }
  else if (ms_sCrysBets == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysBets]   = nRefined;
    }
  else if (ms_sCrysGams == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysGams]   = nRefined;
    }
  else if (ms_sCrysMosaicity == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysMosaicity] = nRefined;
    }
  else if ( (ms_sCrysLengths == sTemp) || ("CrysLength" == sTemp) )
    {
      m_nUserRefineFlags[ms_nCrysAstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysBstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysCstar]  = nRefined;
    }
  else if (ms_sCrysAngles == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysAlps]   = nRefined;
      m_nUserRefineFlags[ms_nCrysBets]   = nRefined;
      m_nUserRefineFlags[ms_nCrysGams]   = nRefined;
    }
  else if (ms_sCrysRot == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysRot1]   = nRefined;
      m_nUserRefineFlags[ms_nCrysRot2]   = nRefined;
      m_nUserRefineFlags[ms_nCrysRot3]   = nRefined;
    }
  else if (ms_sCrysCell == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysAstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysBstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysCstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysAlps]   = nRefined;
      m_nUserRefineFlags[ms_nCrysBets]   = nRefined;
      m_nUserRefineFlags[ms_nCrysGams]   = nRefined;
    }
  else if (ms_sCrysAll == sTemp)
    {
      m_nUserRefineFlags[ms_nCrysAstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysBstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysCstar]  = nRefined;
      m_nUserRefineFlags[ms_nCrysAlps]   = nRefined;
      m_nUserRefineFlags[ms_nCrysBets]   = nRefined;
      m_nUserRefineFlags[ms_nCrysGams]   = nRefined;
      m_nUserRefineFlags[ms_nCrysRot1]   = nRefined;
      m_nUserRefineFlags[ms_nCrysRot2]   = nRefined;
      m_nUserRefineFlags[ms_nCrysRot3]   = nRefined;
      m_nUserRefineFlags[ms_nCrysMosaicity] = nRefined;
    }
  else if (ms_sCrysConstraints == sTemp) {
      m_nUserRefineFlags[ms_nCrysConstraints] = nRefined;
    }
  else if (ms_sCrysRecipShift == sTemp)
  {
          m_nUserRefineFlags[ms_nCrysRecip00] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip01] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip02] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip10] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip11] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip12] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip20] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip21] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip22] = nRefined;
      m_nUserRefineFlags[ms_nCrysRecip30] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip31] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip32] = nRefined;

  }
  else if (ms_sCrysRecipShift0 == sTemp) 
  {
          m_nUserRefineFlags[ms_nCrysRecip00] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip01] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip02] = nRefined;
  }
  else if (ms_sCrysRecipShift1 == sTemp) 
  {
          m_nUserRefineFlags[ms_nCrysRecip10] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip11] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip12] = nRefined;
  }
  else if (ms_sCrysRecipShift2 == sTemp) 
  {
          m_nUserRefineFlags[ms_nCrysRecip20] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip21] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip22] = nRefined;
  }
  else if (ms_sCrysRecipShift3 == sTemp) 
  {
      m_nUserRefineFlags[ms_nCrysRecip30] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip31] = nRefined;
          m_nUserRefineFlags[ms_nCrysRecip32] = nRefined;
  }
  else if (ms_sCrysTwinLaw == sTemp)
  {
          m_nUserRefineFlags[ms_nCrysTwinLaw10] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw11] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw12] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw20] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw21] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw22] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw30] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw31] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw32] = nRefined;
  }
  else if (ms_sCrysTwinLaw1 == sTemp)
  {
          m_nUserRefineFlags[ms_nCrysTwinLaw10] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw11] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw12] = nRefined;

  }
  else if (ms_sCrysTwinLaw2 == sTemp)
  {
          m_nUserRefineFlags[ms_nCrysTwinLaw20] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw21] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw22] = nRefined;

  }
  else if (ms_sCrysTwinLaw3 == sTemp)
  {
          m_nUserRefineFlags[ms_nCrysTwinLaw30] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw31] = nRefined;
          m_nUserRefineFlags[ms_nCrysTwinLaw32] = nRefined;

  }
  else if (ms_sDetAll == sTemp)
    {
      m_nUserRefineFlags[ms_nDetTrans1]  = nRefined;
      m_nUserRefineFlags[ms_nDetTrans2]  = nRefined;
      m_nUserRefineFlags[ms_nDetTrans3]  = nRefined;
      m_nUserRefineFlags[ms_nDetRot1]    = nRefined;
      m_nUserRefineFlags[ms_nDetRot2]    = nRefined;
      m_nUserRefineFlags[ms_nDetRot3]    = nRefined;
    }
  else if (ms_sDetTrans == sTemp)
    {
      m_nUserRefineFlags[ms_nDetTrans1]  = nRefined;
      m_nUserRefineFlags[ms_nDetTrans2]  = nRefined;
      m_nUserRefineFlags[ms_nDetTrans3]  = nRefined;
    }
  else if (ms_sDetRot == sTemp)
    {
      m_nUserRefineFlags[ms_nDetRot1]    = nRefined;
      m_nUserRefineFlags[ms_nDetRot2]    = nRefined;
      m_nUserRefineFlags[ms_nDetRot3]    = nRefined;
    }
  else if (ms_sAll == sTemp)
  {
      if (nRefined == 0) {
          for (i = 0; i <= nMAXPARAMS; i++) {
              if (i != ms_nCrysConstraints)
                  m_nUserRefineFlags[i] = nRefined;
          };
      } else {
          for (i = 0; i <= CELL_MOSAICITY; i++) {
              if (i != ms_nCrysConstraints)
                  m_nUserRefineFlags[i] = nRefined;
          };
          m_nUserRefineFlags[ms_nSourceWave] =  0; // Wavelength refinement always OFF
          
          for (i = CELL_TWINLAW10; i <= CELL_TWINLAW32; i++)
              m_nUserRefineFlags[i] = nRefined;
      };

    }
  else 
  {
      printf("ERROR:  Unrecognized option '%s'\n",sTemp.string());
      return 1;
  }


  vBuildRefineFlags();

  // Now process rest of string if it exists

  sLower = sLower.after(sTemp);
  if (0 < sLower.length())
    return (nSetRefineFlags(sLower));   // Note: recursion!
  else
    return (0);
}

void
Crefine::vBuildRefineFlags() {
    int i;
    
    // Set the refine flag strings based on which parameters were set.
    m_sUserRefineNonCrysFlags = "";
    m_sUserRefineCrysFlags = "";
    
    for (i=0; i<nMAXPARAMS;i++) {
        if (m_nUserRefineFlags[i]) {
            if ((i<=ms_nLastCrysParam) && (i>=ms_nFirstCrysParam)) {
                m_sUserRefineCrysFlags += " +";
                m_sUserRefineCrysFlags += *(ms_psFlags[i]);
                m_sUserRefineCrysFlags += " ";
            } else {
                m_sUserRefineNonCrysFlags += " +";
                m_sUserRefineNonCrysFlags += *(ms_psFlags[i]);
                m_sUserRefineNonCrysFlags += " ";
            };
        };      
    };
};


void
Crefine::vSetResolution(const float fResoMin, const float fResoMax)
{
  // Set the min and max resolution in Angstroms.  Min is closest to
  // beam, max is as far away from beam, thus Min has a higher value
  // in Angstroms than Max

  m_fResolutionMin = max(fResoMin, fResoMax);
  m_fResolutionMax = min(fResoMin, fResoMax);

}
void
Crefine::vSetRejectLimits(const float fRejLimX, const float fRejLimY,
                          const float fRejLimR)
{
  m_a4fRejectLimits[0] = fRejLimX;
  m_a4fRejectLimits[1] = fRejLimY;
  m_a4fRejectLimits[2] = fRejLimR;
  m_a4fRejectLimits[3] = fRejLimR;
}

int
Crefine::nDoRefinement(const Cstring& sRefineOptions,
                       Creflnlist *poReflnlistOther)
{
    nDeleteHeaderEntries();
    
    // Perform the refinement with the command options in sRefineOptions
    
    int      i, j;
    int      nx,ny;
    int      nStat;
    Cstring  sTemp;
    float    fTemp;
    int      nTemp;                        
    int      anSavedUserRefineFlags[nMAXPARAMS]; // Save the previous flags.
    bool     bSubCall;

   
    Cstring *psWords;

    
    
    static Cstring ssErrorMissing
        = "ERROR: Crefine::nDoRefinement - missing option argument!";
    static Cstring ssErrorInvalid
        = "ERROR: Crefine::nDoRefinement - invalid argument> ";

    
    char     a2cSpace[2]= {' ', '\0'};
    int      argc;
    
    float    a2fReso[2];
    float    a3fRejectLimits[3];
    bool     bFirstTime = TRUE;
    int      nTwin;
    
    // Remove leading whitespace from sRefineOptions
    
    bSubCall = false;
    sTemp = sRefineOptions;
    while (   (' '  == sTemp.GetAt(0))
        || ('\n' == sTemp.GetAt(0))
        || ('\t' == sTemp.GetAt(0)))
        sTemp = sTemp.after(0);
    
    // Find number of words in sRefineOptions
    
    argc   = 0;
    i      = -1;
    do
    {
        i = sTemp.find(a2cSpace[0], i+1);
        argc++;
    } while (0 <= i);
    
    // argc will be at least 1
    // Split sTemp into words
    
    psWords = new Cstring [argc];
    split(sTemp, psWords, argc, a2cSpace);
    
    nStat = 0;


    // Save the previous settings.
    for (nx=0;nx<nMAXPARAMS;nx++)
        anSavedUserRefineFlags[nx] = m_nUserRefineFlags[nx];

    
    // Now parse each word
    
    for (i = 0; (i < argc) && (0 == nStat); i++)
    {
        sTemp = psWords[i];
        if (3 < m_nVerbose)
            cout << "Command line string: >>" << sTemp << "<<" << endl;
        if ( ("-verbose" == sTemp) || ("-v" == sTemp) )
        {
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(), "%d", &nTemp);
                if (1 == nStat)
                {
                    m_nVerbose = nTemp;
                    nStat = 0;
                }
                else
                    cout << ssErrorInvalid << psWords[i] << " < !";
            }
            else
            {
                cout << ssErrorMissing;
            }
        }
        else if ("-cycles" == sTemp)
        {
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(), "%d", &nTemp);
                if (1 == nStat)
                {
                    m_nCycles = nTemp;
                    nStat   = 0;
                }
                else
                    cout << ssErrorInvalid << psWords[i] << " < !";
            }
            else
            {
                cout << ssErrorMissing;
            }
        }
        else if ("-noautocorr" == sTemp)
        {
            i++;
            m_nAutoCorrCheck = 0;
            
        } 
        else if ("-converge" == sTemp)
        {
            m_bCheckConvergence = TRUE;
        } 
        else if ("-noconverge" == sTemp)
        {
            m_bCheckConvergence = FALSE;
        } 
        else if ("-detmos" == sTemp)
        {
            m_bPrintDetectorMosaicity = TRUE;
        }
        else if ("-resid" == sTemp)
        {
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(), "%f", &fTemp);
                if ((1 == nStat) || (fTemp<=0.0))
                {
                    m_fResidMin=fTemp;
                    nStat   = 0;
                }
                else
                    cout << ssErrorInvalid << psWords[i] << " < !";
            }
            else
            {
                cout << ssErrorMissing;
            }
        } else if ("-maxlorentz" == sTemp) {
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(), "%f", & fTemp);
                if (1==nStat) {
                    m_fLorentzMax = fTemp;
                    nStat = 0;
                } else
                    cout << ssErrorInvalid << psWords[i] << " < !";
            }
            else
            {
                cout << ssErrorMissing;
            };
        } else if ("-weight" == sTemp) {
            i++;
            nStat = 0;
            if (i < argc) {
                sTemp = psWords[i];
                sTemp.downcase();
                if (sTemp == "rotsigma")
                    m_eWeightType = eWeightRotSigma;
                else if (sTemp == "intensity")
                    m_eWeightType = eWeightIntensity;
                else if (sTemp == "ioversig")
                    m_eWeightType = eWeightIoverSig;
                else if (sTemp == "lorentz")
                    m_eWeightType = eWeightLorentz;
                else if (sTemp == "rotwidth")
                    m_eWeightType = eWeightRotWidth;
                else if (sTemp == "none")
                    m_eWeightType = eWeightNone;
                else
                    nStat = 1;
            } else
                nStat = 2;
            if (1==nStat)
                cout << ssErrorInvalid << psWords[i] << " < !";
            else if (2==nStat)
                cout << ssErrorMissing;
            
        } else if ("-corr" == sTemp) {
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(), "%f", & fTemp);
                if (1==nStat) {
                    m_fMaxCorrCoeff = fTemp;
                    nStat = 0;
                } else
                    cout << ssErrorInvalid << psWords[i] << " < !";
            }
            else
            {
                cout << ssErrorMissing;
            };
        }
        else if ("-sigma" == sTemp)
        {
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(), "%f", &fTemp);
                if (1 == nStat)
                {
                    vSetSigma(fTemp);
                    nStat   = 0;
                }
                else
                    cout << ssErrorInvalid << psWords[i].string() << " < !";
            }
            else
            {
                cout << ssErrorInvalid;
            }
        }
        else if ("-sharpness" == sTemp) 
        {
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(),"%f",&fTemp);
                if (1 == nStat) {
                    vSetSharpness(fTemp);
                    nStat = 0;
                } else
                    cout << ssErrorInvalid << psWords[i].string() << " < !";
            }
            else
            {
                cout << ssErrorInvalid;
            };
        }
        else if ("-checkdoubled" == sTemp) {

            int nTwinID;
            bool bReductionMade;
            float fMinFractionInOdd;
            if ((i+1 < argc) && (1 == sscanf(psWords[i+1],"%f",&fMinFractionInOdd))) {
                i++;
            } else
                fMinFractionInOdd = 0.1f;
            
            bReductionMade = false;
            for (nTwinID = 1; nTwinID <= m_nNumTwins;nTwinID++) {
                if (nDoubledAxes(nTwinID,fMinFractionInOdd)) {
                    bReductionMade = true;
                };
            };
            if (bReductionMade)
                nStat = nDoRefinement((Cstring) "-go",poReflnlistOther);
            

        }
        else if ("-reso" == sTemp)
        {
            i++;
            for (j=0; (j < 2) && (0 == nStat) && (i < argc); j++, i++)
            {
                nStat = sscanf(psWords[i].string(), "%f", &fTemp);
                if (1 == nStat)
                {
                    a2fReso[j] = fTemp;
                    nStat   = 0;
                }
                else
                    cout << ssErrorInvalid << psWords[i] << " < !";
            }
            i--;
            if (2 != j)
            {
                cout << ssErrorMissing;
            }
            else
            {
                vSetResolution(a2fReso[0], a2fReso[1]);
            }
        }
        else if ("-rej" == sTemp)
        {
            i++;
            for (j=0; (j < 3) && (0 == nStat) && (i < argc); j++, i++)
            {
                nStat = sscanf(psWords[i].string(), "%f", &fTemp);
                if (1 == nStat)
                {
                    a3fRejectLimits[j] = fTemp;
                    nStat   = 0;
                }
                else
                    cout << ssErrorInvalid << psWords[i] << " < !";
            }
            i--;
            if (3 != j)
            {
                cout << ssErrorMissing;
            }
            else
            {
                vSetRejectLimits(a3fRejectLimits[0],
                    a3fRejectLimits[1], a3fRejectLimits[2]);
            }
        }
        else if ("-reflimit" == sTemp) 
        {
            
            i++;
            if (i < argc)
            {
                nStat = sscanf(psWords[i].string(), "%d", &nTemp);
                if (1 == nStat)
                {
                    if ((g_sCommandLine.contains("dtrefine ") || g_sCommandLine.contains("dtrefine.exe ")) &&
                        (!g_sCommandLine.contains("-reso")) && (!g_sCommandLine.contains("-seq"))) {
                       
                        int nNumReflections;
                        // Find out how many reflections are above I/sig.
                        nNumReflections = 0;
                        for (nx=0;nx<poReflnlistOther->nGetNumReflns();nx++) {
                            if (((*poReflnlistOther)[nx].fGetSigmaI()>0.0) && ((*poReflnlistOther)[nx].fGetIntensity()/(*poReflnlistOther)[nx].fGetSigmaI()>=m_fSigma)) {
                                nNumReflections++;
                            };
                        };
                        if (nNumReflections>nTemp) {
                            for (nx=0;nx<nNumReflections-nTemp;nx++) {
                                ny = rand() % poReflnlistOther->nGetNumReflns();
                                for (;ny<poReflnlistOther->nGetNumReflns();ny++) {
                                    if (((*poReflnlistOther)[ny].fGetSigmaI()>0.0) && ((*poReflnlistOther)[ny].fGetIntensity()/(*poReflnlistOther)[ny].fGetSigmaI()>=m_fSigma)) 
                                        break;
                                };
                                if (ny<poReflnlistOther->nGetNumReflns())
                                    (*poReflnlistOther)[ny].vSetIntensity(-1.0);
                            };
                        };
                        
                    } else {
                        printf("WARNING:  Option -reflimit incompatible with current command line!  Ignored.\n");
                        
                    };
                    nStat   = 0;
                }
                else
                    cout << ssErrorInvalid << psWords[i].string() << " < !";
            }
            else
            {
                cout << ssErrorInvalid;
            }
            
        }
        else if (("-" == sTemp) || ("--" == sTemp))
        {
        }
        else if ("-wideslice" == sTemp) 
        {
            m_bUseWideSlice = TRUE;
            m_bNotUseWideSlice = FALSE;
        }
        else if ("-nowideslice" == sTemp) 
        {
            m_bNotUseWideSlice = TRUE;
            m_bUseWideSlice = FALSE;
        
        } else if ("-go" == sTemp)
        {


            if (0 >= poReflnlistOther->nGetNumReflns())
            {
                cout << "No reflections available for refinement.\n" << flush;
                delete [] psWords;
                return (-1);
            }
            
            nStat = nCopyReflns(poReflnlistOther);

            if (0 == nStat) {
                if (((1 < m_nVerbose) || (bFirstTime)))
                {
                    nList(1);
                    bFirstTime = FALSE;
                }
                
                nSetup();
                for (nTwin = 1; (nTwin <= m_nNumTwins) && (!nStat); nTwin++) {
                    m_nTwinID = nTwin;
                    if (nTwin > 1) {
                        sTemp = " -levelinfo ";
                        sTemp += " -All ";
                        sTemp += m_sUserRefineCrysFlags;
                        nStat = nDoRefinement(sTemp,poReflnlistOther);
                    };
                    nStat = nRefineLoop(m_nCycles);
                    if (m_poHeader)
                        nUpdateHeader(m_poHeader,nTwin);

                };
                
            }
    
        }
        else if ("-twinreject" == sTemp) {
            double fMaxUnaccountedIntensity;
            double fMaxRMSMM;
            int  nPrevTwins;

            if ((i+1 < argc) && (1 == sscanf(psWords[i+1],"%lf",&fMaxUnaccountedIntensity))) {
                i++;
            } else
                fMaxUnaccountedIntensity = 0.1;
            if ((i+1 < argc) && (1 == sscanf(psWords[i+1],"%lf",&fMaxRMSMM))) {
                i++;
            } else
                fMaxRMSMM = 0.30;

            nPrevTwins = m_nNumTwins;
            nStat = nTwinReject(fMaxUnaccountedIntensity,fMaxRMSMM);
            if ((!nStat) && (m_poHeader)) {
                m_poHeader->nDeleteMask((Cstring) D_K_CrystalMask);
                                for (nTwin = 1; (nTwin <= m_nNumTwins) && (!nStat); nTwin++) 
                                        nUpdateHeader(m_poHeader,nTwin);
            };
            if (nPrevTwins != m_nNumTwins) {
                nStat = nDoRefinement((Cstring) "-go",poReflnlistOther);
                nStat = nDoRefinement((Cstring) "-go",poReflnlistOther);
            };

        }
        else if ("-maxmosaicity" == sTemp)
        {
            if ((i+1 < argc) && (1 == sscanf(psWords[i+1],"%f",&fTemp))) {
                i++;
                m_fTestMosaicityRangeMax = fTemp;
            };
                
        }
        else if ("-display" == sTemp)
        {
            
#if ( !defined(SSI_PC) && !defined(NO_X_WINDOWS) )
        if (NULL == m_poXprop)
        {
            m_poXprop = new CXprop("DTDISPLAY_WINDOW_ID");
        }
#else
            m_nDisplay = 1;
#endif
        }
        else if ("-prompt" == sTemp)
        {
            // Ignore this
        }
        else if ("-testmosaicity" == sTemp)
        {
            // Ignore this
        }
        else if ("-levelinfo" == sTemp)
        {
            // Ignore this
            bSubCall = true;
        }
        else if ("" != sTemp)
        {
            // Assume all others are refinement flags
            
            nStat = nSetRefineFlags(sTemp);
            if (0 != nStat)
            {
                cout << ssErrorInvalid << psWords[i] << " < !" << flush;
            }
        }
    }
    // Restore the previous settings for the refinement flags.
    if (bSubCall) {
        for (nx=0;nx<nMAXPARAMS;nx++)
            m_nUserRefineFlags[nx] = anSavedUserRefineFlags[nx];
        
        vBuildRefineFlags();
    };
    
    

    delete [] psWords;
    return (nStat);

}

int
Crefine::nNotifyDisplay(const Cstring sReflnlistIn, const Cstring sHeaderIn)
{
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  if (NULL != m_poXprop)
    {
      // Tell dtdisplay to display the list

      Cstring sReflnlist;

      if ("" == sReflnlistIn)
        {
          sReflnlist =  sGetCWD() + sDtrekGetPrefix() + "dtrefinetmp.ref";
        }
      else
        {
          sReflnlist = sReflnlistIn;
        }
      if ("" == sHeaderIn)
        m_poXprop->hSetProperty("DTDISPLAY_REFLN_UPDATE",
                              sReflnlist);
      else
        m_poXprop->hSetProperty("DTDISPLAY_REFLN_UPDATE",
                              sReflnlist + " Header: " + sHeaderIn);
    }
  return (0);
#else
  return (1);
#endif
}

#if ((defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))

void Crefine::ReportResults( long Total, long Rejected, long Accepted, long Excluded )
{
    char *acFixed = "fixed";
    CCrclHelper* pHelper = CCrclHelper::GetInstance();
    pHelper->vUtilCmdInit( "UpdateRefinementResults" );
    char pcTemp[60];
    //CCoreTclParam TclCom("UpdateRefinementResults");

    sprintf(pcTemp, "%.4f", m_a6fCrysUnitCell[0]);
    pHelper->vUtilCmdSetParam("a", pcTemp);
    sprintf(pcTemp, "%.4f", m_a6fCrysUnitCell[1]);
    pHelper->vUtilCmdSetParam("b", pcTemp);
    sprintf(pcTemp, "%.4f", m_a6fCrysUnitCell[2]);
    pHelper->vUtilCmdSetParam("c", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysAstar]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysUnitCellSigs[0]);
    }
    pHelper->vUtilCmdSetParam("Sigmaa", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysBstar]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysUnitCellSigs[1]);
    }
    pHelper->vUtilCmdSetParam("Sigmab", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysCstar]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysUnitCellSigs[2]);
    }
    pHelper->vUtilCmdSetParam("Sigmac", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysAstar]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysTurds[0]);
    }
    pHelper->vUtilCmdSetParam("Shifta", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysBstar]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysTurds[1]);
    }
    pHelper->vUtilCmdSetParam("Shiftb", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysCstar]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysTurds[2]);
    }
    pHelper->vUtilCmdSetParam("Shiftc", pcTemp);

    sprintf(pcTemp, "%.4f", m_a6fCrysUnitCell[3]);
    pHelper->vUtilCmdSetParam("alpha", pcTemp);
    sprintf(pcTemp, "%.4f", m_a6fCrysUnitCell[4]);
    pHelper->vUtilCmdSetParam("beta", pcTemp);
    sprintf(pcTemp, "%.4f", m_a6fCrysUnitCell[5]);
    pHelper->vUtilCmdSetParam("gamma", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysAlps]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysUnitCellSigs[3]);
    }
    pHelper->vUtilCmdSetParam("Sigmaalpha", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysBets]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysUnitCellSigs[4]);
    }
    pHelper->vUtilCmdSetParam("Sigmabeta", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysGams]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysUnitCellSigs[5]);
    }
    pHelper->vUtilCmdSetParam("Sigmagamma", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysAlps]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysTurds[3]);
    }
    pHelper->vUtilCmdSetParam("Shiftalpha", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysBets]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysTurds[4]);
    }
    pHelper->vUtilCmdSetParam("Shiftbeta", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysGams]) {
        sprintf(pcTemp, "%.4f", m_a6fCrysTurds[5]);
    }
    pHelper->vUtilCmdSetParam("Shiftgamma", pcTemp);

    sprintf(pcTemp, "%.4f", m_a3fCrysRot[0]);
    pHelper->vUtilCmdSetParam("CrysRot1", pcTemp);
    sprintf(pcTemp, "%.4f", m_a3fCrysRot[1]);
    pHelper->vUtilCmdSetParam("CrysRot2", pcTemp);
    sprintf(pcTemp, "%.4f", m_a3fCrysRot[2]);
    pHelper->vUtilCmdSetParam("CrysRot3", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysRot1]) {
        sprintf(pcTemp, "%.4f", m_a3fCrysRotSigs[0]);
    }
    pHelper->vUtilCmdSetParam("SigmaCrysRot1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysRot2]) {
        sprintf(pcTemp, "%.4f", m_a3fCrysRotSigs[1]);
    }
    pHelper->vUtilCmdSetParam("SigmaCrysRot2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysRot3]) {
        sprintf(pcTemp, "%.4f", m_a3fCrysRotSigs[2]);
    }
    pHelper->vUtilCmdSetParam("SigmaCrysRot3", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysRot1]) {
        sprintf(pcTemp, "%.4f", m_a3fCrysRotShifts[0]);
    }
    pHelper->vUtilCmdSetParam("ShiftCrysRot1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysRot2]) {
        sprintf(pcTemp, "%.4f", m_a3fCrysRotShifts[1]);
    }
    pHelper->vUtilCmdSetParam("ShiftCrysRot2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nCrysRot3]) {
        sprintf(pcTemp, "%.4f", m_a3fCrysRotShifts[2]);
    }
    pHelper->vUtilCmdSetParam("ShiftCrysRot3", pcTemp);

    sprintf(pcTemp, "%.4f ", m_a3fCrysMosaicity[0]);
    pHelper->vUtilCmdSetParam("Mos", pcTemp);

    strcpy(pcTemp, acFixed);
    if (   (1 == m_nUserRefineFlags[ms_nCrysMosaicity])
        && (0.0 < m_a3fCrysMosaicitySig[0]) )
    {
        sprintf(pcTemp, "%.4f ", m_a3fCrysMosaicitySig[0]);
    }
    pHelper->vUtilCmdSetParam("SigmaMos", pcTemp);

    strcpy(pcTemp, acFixed);
    if (   (1 == m_nUserRefineFlags[ms_nCrysMosaicity])
        && (0.0 < m_a3fCrysMosaicitySig[0]) )
    {
        sprintf(pcTemp, "%.4f ", m_a3fCrysMosaicityShift[0]);
    }
    pHelper->vUtilCmdSetParam("ShiftMos", pcTemp);

        // Only handling one detector right now.
    sprintf(pcTemp, "%.4f", m_a10x3fDetTrans[0][0]);
    pHelper->vUtilCmdSetParam("Trans1", pcTemp);
    sprintf(pcTemp, "%.4f", m_a10x3fDetTrans[0][1]);
    pHelper->vUtilCmdSetParam("Trans2", pcTemp);
    sprintf(pcTemp, "%.4f", m_a10x3fDetTrans[0][2]);
    pHelper->vUtilCmdSetParam("Trans3", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetTrans1]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetTransSigs[0][0]);
    }
    pHelper->vUtilCmdSetParam("SigmaTrans1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetTrans2]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetTransSigs[0][1]);
    }
    pHelper->vUtilCmdSetParam("SigmaTrans2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetTrans3]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetTransSigs[0][2]);
    }
    pHelper->vUtilCmdSetParam("SigmaTrans3", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetTrans1]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetTransShifts[0][0]);
    }
    pHelper->vUtilCmdSetParam("ShiftTrans1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetTrans2]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetTransShifts[0][1]);
    }
    pHelper->vUtilCmdSetParam("ShiftTrans2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetTrans3]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetTransShifts[0][2]);
    }
    pHelper->vUtilCmdSetParam("ShiftTrans3", pcTemp);

        // Only handling one detector right now.
    sprintf(pcTemp, "%.4f", m_a10x3fDetRot[0][0]);
    pHelper->vUtilCmdSetParam("DetRot1", pcTemp);
    sprintf(pcTemp, "%.4f", m_a10x3fDetRot[0][1]);
    pHelper->vUtilCmdSetParam("DetRot2", pcTemp);
    sprintf(pcTemp, "%.4f", m_a10x3fDetRot[0][2]);
    pHelper->vUtilCmdSetParam("DetRot3", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetRot1]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetRotSigs[0][0]);
    }
    pHelper->vUtilCmdSetParam("SigmaDetRot1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetRot2]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetRotSigs[0][1]);
    }
    pHelper->vUtilCmdSetParam("SigmaDetRot2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetRot3]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetRotSigs[0][2]);
    }
    pHelper->vUtilCmdSetParam("SigmaDetRot3", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetRot1]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetRotShifts[0][0]);
    }
    pHelper->vUtilCmdSetParam("ShiftDetRot1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetRot2]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetRotShifts[0][1]);
    }
    pHelper->vUtilCmdSetParam("ShiftDetRot2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nDetRot3]) {
        sprintf(pcTemp, "%.4f", m_a10x3fDetRotShifts[0][2]);
    }
    pHelper->vUtilCmdSetParam("ShiftDetRot3", pcTemp);

    sprintf(pcTemp, "%.6f", m_fSourceWave);
    pHelper->vUtilCmdSetParam("Wave", pcTemp);
    sprintf(pcTemp, "%.4f", m_a2fSourceRot[0]);
    pHelper->vUtilCmdSetParam("SRot1", pcTemp);
    sprintf(pcTemp, "%.4f", m_a2fSourceRot[1]);
    pHelper->vUtilCmdSetParam("SRot2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nSourceWave]) {
        sprintf(pcTemp, "%.6f", m_fSourceWaveSig);
    }
    pHelper->vUtilCmdSetParam("SigmaWave", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nSourceRot1]) {
        sprintf(pcTemp, "%.4f", m_a2fSourceRotSigs[0]);
    }
    pHelper->vUtilCmdSetParam("SigmaSRot1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nSourceRot2]) {
        sprintf(pcTemp, "%.4f", m_a2fSourceRotSigs[1]);
    }
    pHelper->vUtilCmdSetParam("SigmaSRot2", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nSourceWave]) {
        sprintf(pcTemp, "%.6f", m_fSourceWaveShift);
    }
    pHelper->vUtilCmdSetParam("ShiftWave", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nSourceRot1]) {
        sprintf(pcTemp, "%.4f", m_a2fSourceRotShifts[0]);
    }
    pHelper->vUtilCmdSetParam("ShiftSRot1", pcTemp);

    strcpy(pcTemp, acFixed);
    if (1 == m_nUserRefineFlags[ms_nSourceRot2]) {
        sprintf(pcTemp, "%.4f", m_a2fSourceRotShifts[1]);
    }
    pHelper->vUtilCmdSetParam("ShiftSRot2", pcTemp);

    sprintf(pcTemp, "%.4f", m_fRMSMM);
    pHelper->vUtilCmdSetParam("mm", pcTemp);
    sprintf(pcTemp, "%.4f", m_fRMSDEG);
    pHelper->vUtilCmdSetParam("deg", pcTemp);
    sprintf(pcTemp, "%ld", Total);
    pHelper->vUtilCmdSetParam("Total", pcTemp);
    sprintf(pcTemp, "%ld", Accepted);
    pHelper->vUtilCmdSetParam("Acc", pcTemp);
    sprintf(pcTemp, "%ld", Rejected);
    pHelper->vUtilCmdSetParam("Rej", pcTemp);
    sprintf(pcTemp, "%ld", Excluded);
    pHelper->vUtilCmdSetParam("Exc", pcTemp);
    // Send the Tcl command to be processed.
    //gSocket.SendString(TclCom.ConstructCommand());
    pHelper->vUtilCmdSend();
}


#endif


int
Crefine::nTwinReject(double fMaxUnaccountedIntensity,double fMaxRMSMM) {
    // This algorithm works by adding components one at a time, until we
    // don't get an appreciable increase in total predicted intensity.

    itr<int> anTwinState;   // State of each component.  0 == available, 1 == taken  2 == remove
    int nTwinID;
    int nNextTwinID;
    int nRef;
    int nGoodTwins;
    int nTwinsInUse;
    double fTotalIntensity;
    char pcBuf[200];

    g_psPrintfLog = &m_sLogTwinReject;

    if (m_anReflnOverlapMat.size() == 0) {
        lprintf("Cannot use the Twin Reject command with -seq refinement!\n"         
               "Please refine off one image.\n");
        return 0;
    };

    // Count the total intensity in the find list
    fTotalIntensity = 0.0;
    for (nRef = 0; nRef < m_anReflnOverlapMat.size(); nRef++) {
        fTotalIntensity += m_afReflnOverlapIntensities[nRef];
    };

    lprintf("Checking twin components for coverage...\n\n");

    nGoodTwins = 0;
    for (nTwinID = 1; nTwinID <= m_nNumTwins; nTwinID++) {
        if ((nTwinID == m_nNumTwins) && (nGoodTwins==0))
            anTwinState[nTwinID] = 0;
        else if ((m_poComp[nTwinID-1].m_fRMSMM>fMaxRMSMM) || (m_poComp[nTwinID-1].m_a7nReflnNumbers[eReflectionNumberAccept]<10)) {
                          
            lprintf("Component %d has an RMS residual of %.4lf and %d accepted reflections.\n"
                   "This component will be rejected.\n\n",
                   nTwinID,
                   (double) m_poComp[nTwinID-1].m_fRMSMM,
                   (double) m_poComp[nTwinID-1].m_a7nReflnNumbers[eReflectionNumberAccept]);
            anTwinState[nTwinID] = 2;
        } else {
            anTwinState[nTwinID] = 0;
            nGoodTwins++;
        };
    };
    nTwinsInUse = 0;

    double fMaxUnaccounted;
    int    nMaxUnaccounted;
    double fUnaccounted;

    do {
        fMaxUnaccounted = 0.0;
        nMaxUnaccounted = 0;
        for (nTwinID = 1; nTwinID <= m_nNumTwins; nTwinID++) {
            if (anTwinState[nTwinID]==0) {
                fUnaccounted = 0.0;
                for (nRef = 0; nRef < m_anReflnOverlapMat.size(); nRef++) {
                    if ((m_anReflnOverlapMat[nRef] & (1 << (nTwinID - 1))) &&
                        (!(nTwinsInUse & m_anReflnOverlapMat[nRef])))
                        fUnaccounted += m_afReflnOverlapIntensities[nRef];
                };
                if (fUnaccounted>fMaxUnaccounted) {
                    fMaxUnaccounted = fUnaccounted;
                    nMaxUnaccounted = nTwinID;
                };
                
            };            
        };
        if (fMaxUnaccounted/fTotalIntensity > fMaxUnaccountedIntensity) {
            lprintf("Component %d accounts for %.2lf percent of the data.\n"
                   "This component will be used.\n\n",nMaxUnaccounted,100.0*fMaxUnaccounted/fTotalIntensity);
            anTwinState[nMaxUnaccounted] = 1;
            nTwinsInUse |= (1 << (nMaxUnaccounted - 1));
            if (fMaxUnaccountedIntensity == 0.0) {
                // This is a flag. It say's: only select the first component.
                fMaxUnaccountedIntensity = 1.0;
            };
        } else if (nTwinsInUse == 0) {
            lprintf("ERROR:  Could not find any twin components that indexed %.2lf percent of the data!\n",
                fMaxUnaccountedIntensity*100.0);
            lprintf("Component %d accounts for only %.2lf percent of the data.\n"
                   "This component will be used.\n\n",
                   nMaxUnaccounted,100.0*fMaxUnaccounted/fTotalIntensity);
                        anTwinState[nMaxUnaccounted] = 1;
            nMaxUnaccounted = 0;
        } else
            nMaxUnaccounted = 0;
    } while (nMaxUnaccounted);

    for (nTwinID = 1; nTwinID <= m_nNumTwins; nTwinID++) {
        if (anTwinState[nTwinID] == 0) {
            lprintf("Component %d accounted for less than %.2lf percent of the data.\n"
                    "This component will be rejected.\n\n",
                    nTwinID,100.0*fMaxUnaccountedIntensity);
        };
    };
    lprintf("End of twin rejection.\n\n");

    // Go through and remove components that were not used.
    
    for (nTwinID = 1; nTwinID <= m_nNumTwins; nTwinID++) {
        if (anTwinState[nTwinID] != 1) {
            
            for (nRef = 0; nRef < m_poReflnlistIn->nGetNumReflns(); nRef++) {
                if ((m_poReflnlistIn->m_nFI_nTwinID>=0) &&
                    (((*m_poReflnlistIn)[nRef].nGetField(m_poReflnlistIn->m_nFI_nTwinID)  % 1000)==nTwinID))
                    (*m_poReflnlistIn)[nRef].vSetField(m_poReflnlistIn->m_nFI_nTwinID,0);
            };
        };
    };


    for (nTwinID = 1,nNextTwinID = 1; nTwinID <= m_nNumTwins; nTwinID++) {
        if (anTwinState[nTwinID] == 1) {
            m_poComp[nNextTwinID-1] = m_poComp[nTwinID-1];
            if ((nNextTwinID != nTwinID) && (m_ppoCrystals[nNextTwinID-1]))
                delete m_ppoCrystals[nNextTwinID-1];
            m_ppoCrystals[nNextTwinID-1] = m_ppoCrystals[nTwinID-1];
            m_ppoCrystals[nNextTwinID-1]->vSetTwin(nNextTwinID);
            if (nNextTwinID != nTwinID)
                m_ppoCrystals[nTwinID-1] = NULL;
            for (nRef = 0; nRef < m_poReflnlistIn->nGetNumReflns(); nRef++) {
                if ((m_poReflnlistIn->m_nFI_nTwinID>=0) &&
                    (((*m_poReflnlistIn)[nRef].nGetField(m_poReflnlistIn->m_nFI_nTwinID) % 1000)==nTwinID))
                    (*m_poReflnlistIn)[nRef].vSetField(m_poReflnlistIn->m_nFI_nTwinID,nNextTwinID);
            };
            nNextTwinID++;
        };
    };

    for (nTwinID = nNextTwinID; nTwinID <= m_nNumTwins; nTwinID++) {
        delete m_ppoCrystals[nTwinID - 1];
        m_ppoCrystals[nTwinID - 1] = NULL;
    };

    if (m_nNumTwins != nNextTwinID - 1) {
        for (nTwinID = 1; nTwinID <= m_nNumTwins; nTwinID++) {
            if (anTwinState[nTwinID] != 1) {
                sprintf(pcBuf,"INFO %2.2d/%2.2d  Removed\n",nTwinID,m_nNumTwins);
                m_sLogRefine += (Cstring) pcBuf;
            };
        };
    };

    m_nNumTwins = nNextTwinID - 1;
    m_nTwinID = 1;

    return 0;
};

int Crefine::nDoubledAxes(int nTwinID,double fMinFractionInOdd)
{
    int         nx,ny;
    int         nRef;
    int         nStat;
    double      fIntensity;
    itr<int>    anHKL;
    itr<int>    anIntensity;
    itr<int>    anRef;
    int*        pnHKL;
    int*        pnIntensity;
    double      fIntensityTotal,fIntensityFails,fIntensityWorks;
    int         nCountWorks,nCountFails;
    char        pcBuf[100];

    -anHKL;
    -anIntensity;
    -anRef;

    g_psPrintfLog = &m_sLogTwinReject;

    printf("Calling doubled axis detection algorithm for component %d...\n",nTwinID);

    if (m_ppoCrystals[nTwinID-1]->m_poSpacegroup->cGetCentFlag() != 'P') {
        lprintf("WARNING:  Will not check for doubled axes in component %d\n"
                "          This cell is already centered.\n",nTwinID);
        return 0;
    };


    // Get the entries that we will be using.
    fIntensityTotal = 0.0;
    for (nRef = 0; nRef < m_poReflnlist->nGetNumReflns(); nRef++) {
        if ((m_poReflnlist->m_nFI_nTwinID<0) ||
            (nTwinID == ((*m_poReflnlist)[nRef].nGetField(m_poReflnlist->m_nFI_nTwinID) % 1000))) {
           
            pnHKL = (*m_poReflnlist)[nRef].pnGetHKL();
            fIntensity = (*m_poReflnlist)[nRef].fGetIntensity();
            for (nx=0;nx<3;nx++)
                anHKL + pnHKL[nx];
            anIntensity + (int) min((double) 0xfffffff, (double) fIntensity);
            anRef + nRef;
            fIntensityTotal += anIntensity.last();
        };
    };
    pnIntensity = & anIntensity[0];
    pnHKL = & anHKL[0];


    const int nMaxIndex1 = 1;
    const int nMaxIndex2 = 4;
    const int nMaxVol = 6;
    int a3x3nTrans[3][3];
    int a3x3nInvTrans[3][3];
    int a3x3nMaxIndex[3][3];
    int a3x3nMinIndex[3][3];
    int nDet;
    int nMatCount;
    int a3nTransHKL[3];
    bool bFound = false;
    
    int a3x3nBestTrans[3][3];
    double a3x3fBestTrans[3][3];
    double a3x3fOrientMatIn[3][3];
    double a3x3fOrientMatOut[3][3];
    
    int         nBestDet = 0;
    
    
    double fVolStart,fVolEnd;


    double fMinFail = 1.0;

    for (nx=0;nx<3;nx++) {
        for (ny=0;ny<3;ny++) {
            a3x3nMaxIndex[nx][ny] = nMaxIndex1;
            a3x3nMinIndex[nx][ny] = -nMaxIndex1;
        };
    };
    for (nx=0;nx<3;nx++) {
        a3x3nMaxIndex[nx][nx] = nMaxIndex2;
        a3x3nMinIndex[nx][nx] = 0;
    };

    // Initialize the transformation matrix.
    for (nx=0;nx < 3; nx++) {
        for (ny=0; ny < 3; ny++) {
            a3x3nTrans[nx][ny] = a3x3nMinIndex[nx][ny];
        };
    };
    a3x3nTrans[0][0]--;
    
    nMatCount = 0;
    do {
        // Increment to the next transformation matrix.
        nStat = 0;
        for (nx=0;(nx < 3) && (!nStat); nx++) {
            for (ny=0; (ny < 3) && (!nStat); ny++) {
                if (a3x3nTrans[nx][ny] == a3x3nMaxIndex[nx][ny])
                    a3x3nTrans[nx][ny] = a3x3nMinIndex[nx][ny];
                else {
                    a3x3nTrans[nx][ny]++;
                    nStat = 1;
                };                
            };
        };
        if (!nStat)
            continue;

        // Invert the reindexing matrix.
        nDet = nInvMat3D(&a3x3nTrans[0][0],&a3x3nInvTrans[0][0]);
        if ((nDet <=1) || (nDet>nMaxVol))
            continue;

        fIntensityWorks = 0.0;
        nCountWorks = 0;
        for (nRef = anIntensity.size() - 1; nRef >= 0; nRef--) {
            vMulMat3DVec3D(&a3x3nInvTrans[0][0],pnHKL + nRef*3,a3nTransHKL);
            // Does the reflection transform into integer entries?
            if (bDivides(3,a3nTransHKL,nDet)) {
                fIntensityWorks += pnIntensity[nRef];
                nCountWorks++;
            };
        };
        fIntensityFails = fIntensityTotal - fIntensityWorks;
        nCountFails = anIntensity.size() - nCountWorks;
        
        fMinFail = min(fMinFail,fIntensityFails/fIntensityTotal);
        if ((fIntensityFails/fIntensityTotal < fMinFractionInOdd) && (nCountWorks>=20)) {
            if ((!bFound) || (nDet>=nBestDet)) {
                vCopyMat3D(&a3x3nTrans[0][0],&a3x3nBestTrans[0][0]);
                nBestDet = nDet;
                bFound = true;
            };
            
        };

        nMatCount++;

        // Loop through each HKL and count the ones that transform correctly.
    } while (nStat);

    if (bFound) {
        // Transform the crystal
        lprintf("\nTransformation matrix found for component %d that\n"
            "preserved %.2lf percent of the input intensity.\n\n",
            nTwinID,100.0*fMinFractionInOdd);
        
        for (nx=0;nx<3;nx++) {
            lprintf("[%c'] = [ ",'h' + nx);
            for (ny=0;ny<3;ny++) {
                lprintf("%4d ",a3x3nBestTrans[nx][ny]);
            };
            lprintf("]%c[%c]\n",(nx==1)?'*':' ','h' + nx);
        };
        lprintf("\n");
        nDet = nInvMat3D(&a3x3nBestTrans[0][0],&a3x3nInvTrans[0][0]);
        lprintf("Matrix determinant is %d\n",nDet);
        lprintf("Former crystal vol:  %.1lf\n",fVolStart = m_ppoCrystals[nTwinID-1]->fCalcVolume());
        lprintf("Former crystal cell: [%6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf]\n",
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(0),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(1),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(2),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(3),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(4),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(5));
        m_ppoCrystals[nTwinID-1]->nCalcOrientMatrix();
        m_ppoCrystals[nTwinID-1]->vGetOrientMatrix(&a3x3fOrientMatIn[0][0]);
        vCopyMat3D(&a3x3nBestTrans[0][0],&a3x3fBestTrans[0][0]);
        vMulMat3DMat3D(a3x3fOrientMatIn,a3x3fBestTrans,a3x3fOrientMatOut);
        m_ppoCrystals[nTwinID-1]->nSetOrientMatrix(&a3x3fOrientMatOut[0][0]);
        lprintf("Transformed vol:     %.1lf\n",fVolEnd = m_ppoCrystals[nTwinID-1]->fCalcVolume());
        lprintf("Transformed cell:    [%6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf]\n",
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(0),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(1),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(2),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(3),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(4),
            (double) m_ppoCrystals[nTwinID-1]->fGetCell(5));
        

        for (nRef = anIntensity.size() - 1; nRef >= 0; nRef--) {
            vMulMat3DVec3D(&a3x3nInvTrans[0][0],pnHKL + nRef*3,a3nTransHKL);
            // Does the reflection transform into integer entries?
            if (!bDivides(3,a3nTransHKL,nDet)) {
                // No?  Then delete the original reflection in the input list.
                // (effectively exclude it by setting it's intensity very low)
                (*m_poReflnlistIn)[anRef[nRef]].vSetIntensity(-99999.0);
            };
        };
        sprintf(pcBuf,"INFO %2.2d/%2.2d %8.0lf was %lf\n",nTwinID,m_nNumTwins,
            (double) fVolEnd,(double) fVolStart);
        m_sLogRefine += pcBuf;
       

        return 1;
    };          

    return 0;
};

int Crefine::nGetReflnsFromFind(Creflnlist& oListFind,const float fResoMin, const float fResoMax)
{
    Cpredict* poPredict;
    Creflnlist* poNewList;
    Creflnlist* poPredList;

    float fPad;
    int nStat;
    double f0;
    bool bIsIntersecting;

    int nRefPred;
    int nRefFind;
    int nRefFindSearch;
    int nPass;
    int* pnSortFind;
    double fCalc0,fCalc1;
    double fCalcRotStart,fCalcRotEnd,fCalcRotWidth,fCalcRotMid;
    double fFind0,fFind1;
    double a2x2fFindEllipseA[2][2];
    double a2fFindEllipseb[2];
    double fFindEllipsec;

    double fFindRotStart,fFindRotEnd,fFindRotWidth,fFindRotMid;
    const double fMaxPixSearch = 40.0;

    itr<double> afObsRotMid;
    itr<int>    anRefln;
    int*        pnOverlapMat;
    int         nObsRef;

    if ((!m_poReflnlistIn->nGetNumReflns()) || (!(*m_poReflnlistIn)[0].bEllipsoidAvailable()))
        return 0;

    poNewList = new Creflnlist(oListFind);
    poPredList = new Creflnlist(oListFind);
    poPredict = new Cpredict (*m_poHeader, poPredList);

    // Create a place to put reflections ...        
    poPredict->vSetMosaicityLinearCoeffs();
    poPredict->vSetMaxLorentz(m_fLorentzMax);
    poPredict->vSetNonunfFlag(1);
    (void) poPredict->nExpandReflnlist(poNewList);
    (void) nExpandReflnlist(poNewList);
    poNewList->nExpandGetField(poNewList->ms_snSourceRefNum);
    
    Crefln oNewRefln(poNewList);

    // Initialize the reflection overlap matrix.
    // This stores a bits for every reflection indicating which components account for that reflection.
    m_afReflnOverlapIntensities.clear();
    m_anReflnOverlapMat.clear();
    m_anReflnOverlapMat.setsize(oListFind.nGetNumReflns());
    pnOverlapMat = & m_anReflnOverlapMat[0];
    memset(pnOverlapMat,0,sizeof(int)*oListFind.nGetNumReflns());

    // Don't worry about duplicates.  We only grab one spot when predicting anyway.
    for (nRefFind = 0;nRefFind < oListFind.nGetNumReflns(); nRefFind++)
        m_afReflnOverlapIntensities[nRefFind] = oListFind[nRefFind].fGetIntensity();

    fPad = max(1.0,m_fFindPad);
        

        if ((0.0f != fResoMin) && (0.0f != fResoMax)
                && (fResoMin != fResoMax) )
                poPredict->vSetResolution(fResoMin, fResoMax);

        // Find the average spot width.
        f0 = 0.0;
        for (nRefFind = 0; nRefFind < oListFind.nGetNumReflns(); nRefFind++) {
                f0 += oListFind[nRefFind].fGetField(oListFind.m_nFI_fObsRotWidth);
        };
        f0 /= nRefFind;
        fPad = max(f0/2.0,fPad);

        
        nStat = poPredict->nPredict(oListFind,fPad);
        if (nStat)
                goto end;

        // Sort the incoming reflection list based on observed pixel coordinates.
        oListFind.vSort(eReflnField_float_type,oListFind.m_nFI_fObsPx0,NULL);
        pnSortFind = oListFind.pnGetSortIndex();

        for (nRefPred = 0; nRefPred < poPredList->nGetNumReflns(); nRefPred++) {
                fCalc0 = (*poPredList)[nRefPred].fGetField(poPredList->m_nFI_fCalcPx0);
                fCalc1 = (*poPredList)[nRefPred].fGetField(poPredList->m_nFI_fCalcPx1);
                fCalcRotStart = (*poPredList)[nRefPred].fGetField(poPredList->m_nFI_fCalcRotStart);
                fCalcRotEnd = (*poPredList)[nRefPred].fGetField(poPredList->m_nFI_fCalcRotEnd);
                fCalcRotWidth = fCalcRotEnd - fCalcRotStart;
                fCalcRotMid = (fCalcRotStart + fCalcRotEnd)*0.5;
                nRefFind = 
                        oListFind.nFindFirst(oListFind.m_nFI_fObsPx0,
                        fCalc0,0.0,NULL,true);
        
        // The loops below will select reflections from the lists.
        // Since we are looping through *predicted* reflections in the outer loop,
        // we must take into account that more than one dtfind.ref reflection might have contribute
        // to the reflection.  We build a list of unique ones here.
        -anRefln;
        -afObsRotMid;

                for (nPass = 0; nPass < 2; nPass++) {
                        nRefFindSearch = nRefFind + ((nPass==0)?0:(-1));
                        while ((nRefFindSearch>=0) && (nRefFindSearch< oListFind.nGetNumReflns())) {
                                fFind0 = oListFind[pnSortFind[nRefFindSearch]].fGetField(oListFind.m_nFI_fObsPx0);
                                fFind1 = oListFind[pnSortFind[nRefFindSearch]].fGetField(oListFind.m_nFI_fObsPx1);
                if (oListFind[pnSortFind[nRefFindSearch]].nPutGetEllipse(false,&a2x2fFindEllipseA[0][0],&a2fFindEllipseb[0],&fFindEllipsec))
                    nBuildEllipse(10.0,&a2x2fFindEllipseA[0][0],&a2fFindEllipseb[0],&fFindEllipsec);

                fFindRotMid = oListFind[pnSortFind[nRefFindSearch]].fGetField(oListFind.m_nFI_fObsRotMid);
                                fFindRotWidth = oListFind[pnSortFind[nRefFindSearch]].fGetField(oListFind.m_nFI_fObsRotWidth);                  
                                fFindRotStart = fFindRotMid - 0.5*fFindRotWidth - fPad;
                                fFindRotEnd = fFindRotMid + 0.5*fFindRotWidth + fPad;
                
                bIsIntersecting = ((ABS(fFindRotEnd - fFindRotStart) + ABS(fCalcRotEnd - fCalcRotStart))>(max(fFindRotEnd,fCalcRotEnd) - min(fFindRotStart,fCalcRotStart)));

                                if ((fEvalEllipse(fCalc0 - fFind0,fCalc1 - fFind1,&a2x2fFindEllipseA[0][0],&a2fFindEllipseb[0],&fFindEllipsec)<=1.0) &&
                    bIsIntersecting)
                                        {
                        anRefln + pnSortFind[nRefFindSearch];
                        afObsRotMid + fFindRotMid;
                    
                    
                                };
                                if (ABS(fFind0-fCalc0)>fMaxPixSearch)
                                        break;
                                nRefFindSearch += ((nPass==0)?1:(-1));
                        };
            for (nObsRef = 0; nObsRef < anRefln.size(); nObsRef++) {
                // Add the spot(s) to the reflection list.
                oNewRefln = (*poPredList)[nRefPred];
                oNewRefln = oListFind[anRefln[nObsRef]];
                                oNewRefln.vSetH((*poPredList)[nRefPred].nGetH());
                                oNewRefln.vSetK((*poPredList)[nRefPred].nGetK());
                                oNewRefln.vSetL((*poPredList)[nRefPred].nGetL());
                                if (poPredList->m_nFI_nTwinID>=0)
                                        oNewRefln.vSetField(poNewList->m_nFI_nTwinID,(*poPredList)[nRefPred].nGetField(poPredList->m_nFI_nTwinID));
                oNewRefln.vSetField(poNewList->m_nFI_nSourceRefNum,anRefln[nObsRef]);
                poNewList->nInsert(&oNewRefln);
                pnOverlapMat[anRefln[nObsRef]] |= 
                    (1 << ((poPredList->m_nFI_nTwinID>=0)?(((*poPredList)[nRefPred].nGetField(poPredList->m_nFI_nTwinID) % 1000) - 1):0));
            };
                };
        };

        // Create 'rejects' file from the reflections that didn't get included.
        // This can be passed on to indexing for another pass if we want.
        if (m_sRejectFile.length()) {
                Creflnlist oListFindReject(oListFind);
                for (nObsRef = 0; nObsRef < oListFind.nGetNumReflns(); nObsRef++) {
                        if (!(pnOverlapMat[nObsRef])) 
                                oListFindReject.nInsert(&oListFind[nObsRef]);
                };
                oListFindReject.nWrite(m_sRejectFile);
        };



    m_poReflnlistIn->vDeleteAll();
    poPredict->nExpandReflnlist(m_poReflnlistIn);

    // Merge the reflections here (but include twin merging)
    {
    Cfind oFind(*m_poHeader,poNewList);
    oFind.nMerge2DReflectionsWithTwins(m_poReflnlistIn);
    };


end: 
        delete poPredict;
        delete poNewList;
        delete poPredList;

        return nStat;
};

int Crefine::nGetReflnsFromImages(const int nNumImages, 
                                  const int *pnSeqNum,
                                  Cimage_header *poHeader,
                                  const float fResoMin, 
                                  const float fResoMax,
                                  DTREK_WORD wCtrl)
{
    int nStat;
    
    int nNumSeqOptions = 0;
    int nNumOldRefs    = 0;
    Cpredict *poPredict         = NULL;
    Cscan    *poScan            = NULL;
    Crotation*poRotation        = NULL;
    Cfind    *poFind            = NULL;
    Cimage   *poImage           = NULL;
    Crefln   *poRefln           = NULL;
    Cgoniometer* poGoniometer   = NULL;
    int      i, j;
    float    fRotStart, fRotEnd, fRotMid,fRotWidth;
    int      nSeq;
    int      nNumImagesUsed;
    
    if ( (0 >= nNumImages) || (NULL == pnSeqNum) )
        return (-1);
    
    float fTest = 0.98f;
    float fDelta      = 0.06f;
    CScanBitmap oScanBitmap;
      

    nStat = 0;
    nNumImagesUsed = nNumImages;
    
    if( 0.0f >= m_fSigma )
        m_fSigma = 5.0f;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Figure out rotation pad.
    float fPad = 0.0f;
    // When we predict and search for reflections for ranking purposes, we should not pad the rotation width of an image.
    // Otherwise, if predicted reflections share the same real spot on the image, we will end up with the
    // overestimated count of the real reflections on the image.
    if( !(wCtrl & DTREK_REFINE_GRFI_FOR_RANKING) )
    {
        // Thad's (?) comment: "For mosaicity predictions, we set an extremely large mosaicity."
        if( wCtrl & DTREK_REFINE_GRFI_PREDICT_MOSAICITY )
            fPad = m_fFindPad;
        else
            fPad = max(1.0, m_fFindPad);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Clean up the list
    m_poReflnlistIn->vDeleteAll();
    
    // For each specified image ...
    int nCountNonRejected = 0;
    for (i = 0; (i < nNumImagesUsed) && (0 == nStat); i++)
    {
        // ... predict spots on that image,
        // ... get centroids of predicted spots
        
        nSeq = pnSeqNum[i];
        
        nNumSeqOptions++;
        
        if (0 == i)
        {
            // First -seq option, so do some initializations
            
            cout << "\nReflection centroids will come from a search of image(s),"
                << "\nand not from a reflnlist file.\n" << flush;
            poPredict = new Cpredict (*poHeader, NULL);

            poPredict->vSetMosaicityLinearCoeffs();
            
            // Do not predict refs in bad nonunf regions
            
            poPredict->vSetMaxLorentz(m_fLorentzMax);
            poPredict->vSetNonunfFlag(1);
            
            // Create a place to put reflections ...
            
            (void) poPredict->nExpandReflnlist(m_poReflnlistIn);
            
            poFind = new Cfind(*poHeader, m_poReflnlistIn, 
                               poPredict->poGetNonunf());
            
            poFind->m_bIncludeBackground = wCtrl & DTREK_REFINE_GRFI_FOR_RANKING ? true : false;
            
            poFind->nExpandReflnlist(m_poReflnlistIn,false,true);


            poScan  = new Cscan(*poHeader);
            poImage = new Cimage ();
                        poImage->m_bReadApplyGonioMask = true;
                        poImage->m_bReadApplyEmbeddedMask = true;

        }
        
        
        poScan->vSetSeqNum(nSeq);
        nStat = poScan->nGetImage(poImage);
        if (0 != nStat)
        {
            cout << "ERROR getting image!\n" << flush;
            break;
        }
        if (m_sMapFile.length()) {
            if (i==0) {
                nStat = oScanBitmap.nOpenOutput(m_sMapFile,poImage->nGetDimension(0),poImage->nGetDimension(0),pnSeqNum[0],pnSeqNum[nNumImagesUsed-1],3,0);
                if (nStat) {
                    printf("ERROR:  Opening map file for write.\n");
                    break;
                };
            } 
            oScanBitmap.vClearLocalBitmap();
        };

                if (C3Ddata::m_poProfit) {
                        poImage->nGetDimensions(&C3Ddata::m_poProfit->m_nDim0,&C3Ddata::m_poProfit->m_nDim1);
                        C3Ddata::m_poProfit->vSetCurrentImage(nSeq);
                };
        
        // Reload the image header from the image... NOT from the scan.
        // This solves problems with screening images.
        if (poRotation)
            delete poRotation;
        poRotation = new Crotation(poImage->m_oHeader);

        // Update the crystal goniometer as well.
        poGoniometer = new Cgoniometer(poImage->m_oHeader,Ccrystal::ms_sCrystalPrefix);
        if (poGoniometer->nReconcileHeaderWithImageGoniometer(poPredict->m_poCrysGonio)) {
            nStat = 1;
            break;
        };
        delete poPredict->m_poCrysGonio;
        poPredict->m_poCrysGonio = poGoniometer;

        if(m_ppoDetector != NULL)
                m_ppoDetector[0]->nCalcPixelShiftArray(*m_poSource,*poRotation,*poGoniometer);
           
                      
        cout << "\n... predicting reflections and getting centroids for image "
            << nSeq << " ...\n" << flush;

        if ("" == sGetEnv("DTREK_DONT_TRUST_IMAGE_SCAN_INFO")) {
            fRotStart = poRotation->fGetRotStart();
            fRotEnd   = poRotation->fGetRotEnd();
        } else {      
            fRotStart = poScan->fCalcGetRotStart(nSeq);
            fRotEnd   = poScan->fCalcGetRotEnd(nSeq);
        };

        fRotMid   = (fRotStart + fRotEnd) * 0.5f;
        fRotWidth = (fRotEnd - fRotStart);
        nStat     = poPredict->m_poRotation->nSetRotRange(fRotStart - fPad, fRotEnd + fPad);
        
        nNumOldRefs = m_poReflnlistIn->nGetNumReflns();
        
        if (   (0.0f != fResoMin) && (0.0f != fResoMax)
            && (fResoMin != fResoMax) )
            poPredict->vSetResolution(fResoMin, fResoMax);
        
        poFind->nMarkGonioNumbers();
        nStat = poPredict->nPredict(NULL, m_poReflnlistIn);                      
    
        poFind->m_fSigma = m_fSigma;

        poFind->nAddGonioNumbers(*poImage);
        
        for (j = nNumOldRefs; j < poFind->m_poReflnlist->nGetNumReflns(); j++)
        {
            poRefln = poFind->m_poReflnlist->poGetRefln(j);

            /// FOR DEBUG PURPOSES
            //float fJunk = poRefln->fGetField(poFind->m_poReflnlist->m_nFI_fResolution);
            //
            //int   nH = poRefln->nGetH();
            //int   nK = poRefln->nGetK();
            //int   nL = poRefln->nGetL();
            //
            //if(nH==-11&&nK==-9&&nL==36)
            //{
            //    int     nTrap = 0;
            //}
            ////////////////////////////////

            // Copy calc posn to observed posn, except for calcrotmid
            poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx0,
                               poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx0));
            
            poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsPx1,
                poRefln->fGetField(poRefln->m_poReflnlist->m_nFI_fCalcPx1));
            
            poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotMid, fRotMid);
            poRefln->vSetField(poRefln->m_poReflnlist->m_nFI_fObsRotWidth, fRotWidth);
                       
            DTREK_WORD  wFC2DCtrl = wCtrl & DTREK_REFINE_GRFI_FOR_RANKING ? DTREK_DTFIND_FC2D_REPORT_WEAK_SPOTS : 0U; 
            
            if( m_bRanking )
                wFC2DCtrl |= DTREK_DTFIND_FC2D_GET_PEAK_AREA_INFO;

            nStat = poFind->nFindCentroid2D(*poImage, 
                                            poRefln,
                                            NULL,                   // ptr to Detector object not necessary
                                            wFC2DCtrl);
            
            if( 0 != nStat )
            {
                // Reset local status and mark refln for deletion                

                nStat = 0;
                poRefln->vSetIntensity(-999.0);
                nCountNonRejected++;
            } 
            else if (m_sMapFile.length())
                oScanBitmap.vSetLocalBitmap(poFind->m_oShoebox,0);

        }
        if (m_sMapFile.length())
            oScanBitmap.nWriteNextBitmap();

    } // end of image loop i
        
    
    // Merge 2D centroids (but only if we are not predicting mosaicity)
    if( !(wCtrl & DTREK_REFINE_GRFI_PREDICT_MOSAICITY) )
    {
        Creflnlist* poNewList = new Creflnlist(*m_poReflnlistIn);            

        // When we do a refinement, we don't want to use spots that have more than one predicted HKL, 
        // because that makes the HKL assignment of the spot ambiguous, so we reject all such spots.
        bool    bRejectOverlappedDifferentHKL = m_bUseOverlapped ? false : true;
        poFind->nMerge2DReflectionsWithTwins(poNewList, bRejectOverlappedDifferentHKL);
        
        m_poReflnlistIn->vDeleteAll();
        m_poReflnlistIn->nInsertListFrom(*poNewList);
        
        delete poNewList;
    }
    
    if( !(wCtrl & DTREK_REFINE_GRFI_PREDICT_MOSAICITY) )
        printf("\nTotal number of predicted merged reflections: %d\n", m_poReflnlistIn->nGetNumReflns());
    else
        printf("\nTotal number of predicted (non-merged) reflections: %d\n", m_poReflnlistIn->nGetNumReflns());

    if (m_sMapFile.length()) {
        oScanBitmap.nCloseOutput();
        oScanBitmap.nUpdateHeader(*m_poHeader);
    };
    
    delete poImage;
    delete poFind;
    delete poScan;
    delete poRotation;
    delete poPredict;
    if (0 >= m_poReflnlistIn->nGetNumReflns())
        nStat = -1;
    return (nStat);
}


int Crefine::nGetMosaicity(const int nNumImages, const int *pnSeqNum,
                           Cimage_header *poHeader,
                           const float fResoMin,
                           const float fResoMax,
                           const int   nGonioSet,
                           const float* pfCalcRotMid,
                           float* pfRotOffset,
                           int* pnRejectFlag,
                           int* pnGonioMatNum
                           )
{

  const double fMosaicityTestStep=0.02;            // Testp mosaicity step.
        double fMaxMosaicity;                       // Maximum mosaicity to test.
        double fMinMosaicity;
        
        double fRotOffsetRange = 0.4;             // Maximum rot offset to test.
  const double fRotOffsetStep = 0.05;             // Minimum rot offset to test.
        double fRotOffsetOffset = 0.0;            // Used if a large rot offset is needed.  Outer most loop changes this.
  const int    nRotBins = (int) 500;
  const bool   bUseWeightedReflns = true;         // Should we use weighted reflections?
  const bool   bUseZeroMosaicityReflns = false;   // Should we use reflections which are predicted even when mosaicity is 0.0?
  const double fMosaicityFraction = 0.99;

        double fMosaicity = 0.0;                  // Computed from nMosaicity.
        int    nRotOffset = 0;                    // Rot counter (used to index into anRot arrays.
        double fRotOffset = 0.0;                  // Rotation offset.
        double fTotalReflns = 0.0;
        int    nTotalReflns = 0;

        double          fRotationWidth = 0.0;     // Scan information.
        double          fRotationStart = 0.0;     // Scan information.
        Creflnlist*     poReflnlist = NULL;

        int             anRotReflns[nRotBins];    // These three arrays store rot-offset dependent info.
        double          afRotValues[nRotBins];
        double          afRotMosaicity[nRotBins];

        bool            bUseThisMosaicity = false;
        int             nBestReflns = 0;         // A copy of nReflns for the best mosaicity.
        double          fBestMosaicity = 0.0;        // A copy of fMosaicity for the best mosaicity.
        double          fBestOffset = 0.0;           // A copy of fRotOffset for the best mosaicity.
        double          fBogusOffset = -100.0;
        bool            bLastLoop = false;
        double          fReflns = 0.0;
        int             nReflns = 0;

    // Temporary variables.

    int         nx,ny;
    double      f0,f1,f2;
    int         nStat;
    int         nRef;

    // Set the maximum mosaicity that we would allow.

    fMaxMosaicity = m_fTestMosaicityRangeMax;
    fMinMosaicity = 0.05;

    if (pfCalcRotMid)
      {
        poReflnlist = m_poReflnlist;
        // Use the reflection list already generated for refinement.
      } 
    else 
      {
        // Predict reflections.
        nStat = nGetReflnsFromImages(nNumImages,
                                     pnSeqNum,
                                     poHeader,
                                     fResoMin,
                                     fResoMax,
                                     DTREK_REFINE_GRFI_PREDICT_MOSAICITY);
        if (nStat)
          return nStat;
        poReflnlist = m_poReflnlistIn;
      };

    if ((poReflnlist->m_nFI_fCalcMosCoeffA<0) || (poReflnlist->m_nFI_fCalcMosCoeffB<0))
      return 1;

    
    fRotationWidth = m_poRotation->fGetIncrement();
    fRotationStart = m_poRotation->fGetRotStart();
    bLastLoop = false;

    if (fRotOffsetRange < fRotationWidth*0.05)
      fRotOffsetRange = fRotationWidth*0.05;    // 1/20th of image rotation width

    // Count the total # of reflections in use.
    fTotalReflns = 0.0;
    nTotalReflns = 0;
    for (nRef = 0; nRef<poReflnlist->nGetNumReflns(); nRef++)
      {
        if (   (poReflnlist->m_nFI_nTwinID<0) 
            || (m_nTwinID == ((*poReflnlist)[nRef].nGetField(poReflnlist->m_nFI_nTwinID) % 1000))) 
          {
            if (   ((!pnRejectFlag) || (!pnRejectFlag[nRef]))
                && ((!pnGonioMatNum) || (pnGonioMatNum[nRef] == nGonioSet)))
              {
                double fObsRotStart;
                double fObsRotEnd;
                double fObsRotOffset;
                                
                // Fill in basic fields.
                fObsRotStart = (*poReflnlist)[nRef].fGetField(poReflnlist->m_nFI_fObsRotMid)
                                - (f2 = (*poReflnlist)[nRef].fGetField(poReflnlist->m_nFI_fObsRotWidth)/2);
                fObsRotEnd = fObsRotStart +  f2 * 2.0;
                                
                f0 = (floor(720.95 + (fObsRotStart - fRotationStart)/fRotationWidth) - 720)*fRotationWidth + fRotationStart;
                f1 = (floor(720.95 + (fObsRotEnd - fRotationStart)/fRotationWidth) - 720)*fRotationWidth + fRotationStart;
                if (f0 == f1) 
                  {
                    if (ABS(fObsRotStart - f0) < ABS(fObsRotEnd - f0))
                      f1 += fRotationWidth;
                    else
                      f0 -= fRotationWidth;
                  };
                                
                fObsRotStart = f0;
                fObsRotEnd = f1;
                                
                for (fObsRotOffset = 0.0; 
                     ABS(fObsRotOffset - (fObsRotEnd - fObsRotStart))/fRotationWidth>0.01; 
                     fObsRotOffset += fRotationWidth) 
                  {
                    nTotalReflns++;
                    fTotalReflns += (*poReflnlist)[nRef].fGetIntensity();
                  };
              };
          };
      }; //     for (nRef = 0; ...

//+2010-01-27 JWP
    // Prevent too many of the next do loops
    int nLoopNum = 0;
//-2010-01-27 JWP
    
    do {
      double fStartMosaicity = fMinMosaicity;
      double fEndMosaicity   = fMaxMosaicity;

      if (bLastLoop)
        fStartMosaicity = fEndMosaicity = fBestMosaicity;
      else 
        {
          fBestMosaicity = fMaxMosaicity;
          fBestOffset = fBogusOffset;
        };
      //cout << "DBG: after A1: fMosaicity, fEndMosaicity, fBestMosaicity: "
      //<< fMosaicity << ", " << fEndMosaicity << ", " << fBestMosaicity << endl;


      nBestReflns = 0;
      nStat = 0;
      // Load the mosaicity values.
      for (fRotOffset = -fRotOffsetRange + fRotOffsetOffset,nRotOffset = 0; 
           fRotOffset <= fRotOffsetRange + fRotOffsetOffset; 
           fRotOffset += fRotOffsetStep,nRotOffset++)
        {
          int nBisectionStat = 0;
          int nBisectionState = 0;
          double fBisectionTotalReflns = 0.0;

          while (!(nBisectionStat 
                   = nBisection(nBisectionState,fStartMosaicity,fEndMosaicity,fMosaicity,fBisectionTotalReflns,(bUseWeightedReflns)?(fMosaicityFraction*fTotalReflns):(fMosaicityFraction*nTotalReflns),fMosaicityTestStep)))
            {
              fReflns = 0.0;
              nReflns = 0;

              for (nRef = 0; nRef<poReflnlist->nGetNumReflns(); nRef++)
                {
                  if (   (poReflnlist->m_nFI_nTwinID<0)
                      || (m_nTwinID == ((*poReflnlist)[nRef].nGetField(poReflnlist->m_nFI_nTwinID) % 1000)))
                    {
                      double fCalcRotStart;
                      double fCalcRotEnd;
                      double fCalcRotMid;
                      double fObsRotStart;
                      double fObsRotEnd;
                      double fObsRotStartOffset;
                      double fObsRotEndOffset;
                      double fObsRotOffset;
                      double fCoeffA;
                      double fCoeffB;
                        
                      if (   ((!pnRejectFlag) || (!pnRejectFlag[nRef]))
                          && ((!pnGonioMatNum) || (pnGonioMatNum[nRef] == nGonioSet)))
                        {
                          // Fill in basic fields.
                          fObsRotStart = (*poReflnlist)[nRef].fGetField(poReflnlist->m_nFI_fObsRotMid) 
                                       - (f2 = (*poReflnlist)[nRef].fGetField(poReflnlist->m_nFI_fObsRotWidth)/2);
                          fObsRotEnd = fObsRotStart + f2*2;
                            
                          f0 = (floor(720.95 + (fObsRotStart - fRotationStart)/fRotationWidth) - 720)*fRotationWidth + fRotationStart;
                          f1 = (floor(720.95 + (fObsRotEnd - fRotationStart)/fRotationWidth) - 720)*fRotationWidth + fRotationStart;
                          if (f0 == f1) {
                            if (ABS(fObsRotStart - f0) < ABS(fObsRotEnd - f0))
                              f1 += fRotationWidth;
                            else
                              f0 -= fRotationWidth;
                          };
                            
                          fObsRotStart = f0;
                          fObsRotEnd = f1;
                            
                          if ((f2>0.0) && (fObsRotStart != fObsRotEnd))
                            {
                              fCoeffA = (*poReflnlist)[nRef].fGetField(poReflnlist->m_nFI_fCalcMosCoeffA)/Gs_dRADIANS_PER_DEGREE;
                              fCoeffB = (*poReflnlist)[nRef].fGetField(poReflnlist->m_nFI_fCalcMosCoeffB)/Gs_dRADIANS_PER_DEGREE;
                              if (pfCalcRotMid)
                                fCalcRotMid = pfCalcRotMid[nRef];
                              else
                                fCalcRotMid = (*poReflnlist)[nRef].fGetField(poReflnlist->m_nFI_fCalcRotMid);
                                
                              f0 = fCoeffA*fMosaicity + fCoeffB;
                              fCalcRotStart = fRotOffset + fCalcRotMid -  f0 / 2.0;
                              fCalcRotEnd   = fRotOffset + fCalcRotMid +  f0 / 2.0;
                                
                              // The reflection might be the merge of several reflections. Step over each image.

                              for (fObsRotOffset = 0.0; 
                                   ABS(fObsRotOffset - (fObsRotEnd - fObsRotStart))/fRotationWidth>0.01; 
                                   fObsRotOffset += fRotationWidth)
                                {
                                  fObsRotStartOffset = fObsRotStart + fObsRotOffset;
                                  fObsRotEndOffset = fObsRotStartOffset + fRotationWidth;
                                    
                                  // See if there is an overlap.
                                  if ((max(fObsRotEndOffset,fCalcRotEnd) 
                                           - min(fObsRotStartOffset,fCalcRotStart) 
                                      < fObsRotEndOffset - fObsRotStartOffset + fCalcRotEnd - fCalcRotStart))
                                    {
                                      // This last line makes sure to only use reflections which are partial
                                      // from a different slice.
                                      // if ((0.5*(fCalcRotStart+fCalcRotEnd)>fObsRotEndOffset) 
                                      //   || (0.5*(fCalcRotStart+fCalcRotEnd)<fObsRotStartOffset) || bUseZeroMosaicityReflns) {
                                      nReflns ++;
                                      fReflns += (*poReflnlist)[nRef].fGetIntensity();
                                    } 
                                };
                            };
                        };
                    };
                }; //         for (nRef = 0; nRef<poReflnlist->nGetNumReflns(); nRef++)
              //cout << "DBG: after D ...\n";   
                                
              if (bUseWeightedReflns)
                {
                  fBisectionTotalReflns = fReflns;
                  if ((fMosaicity == fStartMosaicity) && (fBisectionTotalReflns >= fMosaicityFraction*fTotalReflns))
                    break;
                } 
              else 
                {
                  fBisectionTotalReflns = nReflns;
                  if ((fMosaicity == fStartMosaicity) && (fBisectionTotalReflns >= fMosaicityFraction*nTotalReflns))
                    break;
                };
            };  // while (!(nBisectionStat ...
          //cout << "DBG: after C ...\n";   
          if (nBisectionStat == 2)
            fMosaicity = fEndMosaicity;
          else 
            fMosaicity = fMosaicity;
          
          //cout << "DBG: after C1: fMosaicity, fEndMosaicity, fBestMosaicity: "
          // << fMosaicity << ", " << fEndMosaicity << ", " << fBestMosaicity << endl;

          if (bLastLoop)
            {
              // We are only interested in modifying the number of reflections using this mosaicity.
              anRotReflns[nRotOffset] = nReflns;                      
            } 
          else 
            {
              if (!m_nAdjustRotation)
                {
                  if (fMosaicity < fBestMosaicity)
                    // if (ABS(fRotOffset) < ABS(fBestOffset)) 
                    bUseThisMosaicity = true;
                  else
                    bUseThisMosaicity = false;
                } 
              else if (   (fMosaicity < fBestMosaicity) 
                       || (    ((fBestMosaicity == fMosaicity) 
                            || (fBestMosaicity == fMaxMosaicity)) && (nReflns >= nBestReflns)))
                {
                                        
                  if ((nReflns != nBestReflns) || (ABS(fRotOffset) < ABS(fBestOffset)))
                    {
                      bUseThisMosaicity = true;
                    }
                  else
                    bUseThisMosaicity = false;
                } 
              else
                bUseThisMosaicity = false;

              if (bUseThisMosaicity)
                {
                  fBestMosaicity = fMosaicity;
                  nBestReflns = nReflns;
                  if (m_nAdjustRotation) 
                    fBestOffset = fRotOffset;
                  else 
                    fBestOffset = 0.0;
                };
                
              anRotReflns[nRotOffset] = nReflns;
              afRotValues[nRotOffset] = fRotOffset;
              afRotMosaicity[nRotOffset] = fMosaicity;

            }; // ... else not bLastLoop

        }; // for (fRotOffset = -fRotOffsetRange + fRotOffsetOffset,nRotOffset = 0; 
      //cout << "DBG: after B ...\n";   
      if ((bLastLoop) || (fRotOffsetRange==0.0))
        break;
      else if (fBestOffset - fRotOffsetOffset > 0.6*fRotOffsetRange)
        {
          fRotOffsetOffset += fRotOffsetRange;
        } 
      else if (fRotOffsetOffset - fBestOffset > 0.6*fRotOffsetRange)
        {
          fRotOffsetOffset -= fRotOffsetRange;
        } 
      else 
        bLastLoop = true;

//+2010-01-27 JWP
      // Continue looping until we find the local min or 10 loops.
    } while (1 && (10 > nLoopNum++));    // do { ...
    //cout << "DBG: after A ...\n";
//+2010-01-27 JWP

    // Print the table for the best rot. offset.

    //+JWP 2010-02-01    
    //  Do not change offset if it was never really calculated
    //  Do not change mosaicity if there was a problem

    bool bNoErrorSoPrint = true;
    if (fBestOffset == fBogusOffset)
      {
        fBestOffset = 0.0;
        bNoErrorSoPrint = false;
      }
    if (fBestMosaicity >= m_fTestMosaicityRangeMax)
      {
        fBestOffset = 0.0;
        bNoErrorSoPrint = false;
      }
    
    //-JWP 2010-02-01
    
    if ( (0 < m_nVerbose) && (bNoErrorSoPrint) )
      {
        printf("Number of found reflections and best mosaicity vs Trial rotation offset\n"
               "------------------------------------------------------------------------------\n"
               "                             %% of Max number of observed reflections*\n");
        printf(" Offset Mosaicity^ reflns 0...10...20...30...40...50...60...70...80...90...100\n"
               "------------------------------------------------------------------------------\n");
        int nMaxRotOffsetReflns = anRotReflns[0];
        int nMinRotOffsetReflns = anRotReflns[0];
        for (nx=1;nx<nRotOffset;nx++)
          {
            nMaxRotOffsetReflns = max(nMaxRotOffsetReflns,anRotReflns[nx]);
            nMinRotOffsetReflns = min(nMinRotOffsetReflns,anRotReflns[nx]);
          };
        for (nx=0;nx<nRotOffset;nx++)
          {
            printf(" %6.2lf %9.2lf %6d ",
                   afRotValues[nx],
                   afRotMosaicity[nx],
                   anRotReflns[nx]);
            f0 = 0.2+ 0.8*max(0.0,min(1.0,((double) (anRotReflns[nx] - nMinRotOffsetReflns))/max(1,(nMaxRotOffsetReflns - nMinRotOffsetReflns))));
            for (ny = 0;ny < f0*52;ny++) 
              printf("*");
            printf("\n");
          };
        printf("------------------------------------------------------------------------------\n");
        printf(" ^Mosaicity value needed to get maximum number of reflns with this offset.\n");
        printf(" *When mosaicity is set to %4.2f degrees.\n", 
               (float)fBestMosaicity);
        printf("\n\n");
    
        vSetMosaicity(fBestMosaicity);
      }
    else if (0 < m_nVerbose)
      {
        printf("\nWARNING determining Mosaicity value or offset.  Left unchanged.\n\n");
        // Note: no call to ::vSetMosaicity()
      }

    if (0 != nStat)
      {
        printf("\nWARNING:  Could not estimate mosaicity from image(s).\n"
               "          Higher mosaicity always gives more reflections!\n"
               " -> Consider making max mosaicity (now at %4.2f) larger.\n\n", fMaxMosaicity);
      }
    else if (0 < m_nVerbose)
      {
        if (  (1 == m_nUserRefineFlags[ms_nCrysMosaicity])  && bNoErrorSoPrint)
          {
            printf("INFO: Mosaicity set to %4.2f degrees.\n", (float)fBestMosaicity);

            // Only change the offset if mosaicity is refined.

            if (pfRotOffset)
              {
                *pfRotOffset = (float)(ABS(fBestOffset)<fRotOffsetStep/2)?0.0:fBestOffset;
                printf("INFO: Rotation Offset of %4.2f degrees applied to obtain lowest mosaicity.\n\n",
                       *pfRotOffset);
              }
          }
        else if (!bNoErrorSoPrint)
          {
            printf("INFO: Mosaicity kept fixed at %4.2f deg, though estimated to be %4.2f deg.\n",
                   m_a3fCrysMosaicity[0], (float)fBestMosaicity);
            if (pfRotOffset)
              printf("      Rotation Offset not changed: %4.2f degrees\n\n", *pfRotOffset);
          }
        else
          {
            printf("INFO: Mosaicity kept fixed at %4.2f deg, because no new estimate available.\n",
                   m_a3fCrysMosaicity[0]);
            if (pfRotOffset)
              printf("      Rotation Offset not changed: %4.2f degrees\n\n", *pfRotOffset);

          }
        fflush(stdout);
      }
    return 0;
};

void Crefine::vSetMosaicity(double fMosaicity) {
    // Update mosaicity values. (This is usefull if the routine is called from some other program).
    m_ppoCrystals[m_nTwinID-1]->vSetMosaicity(fMosaicity);
    return;
};

int Crefine::nForceCell(double a6fCell[6])
{
  eSpacegroup_classes eClass;
  float f0;

  eClass = m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->eGetClass();

  if (eClass == eSpacegroup_triclinic_class)
    {
    }
  else if (eClass == eSpacegroup_monoclinic_class)
    {
      a6fCell[3] = 90.0;
      a6fCell[5] = 90.0;
    }
  else if (eClass == eSpacegroup_orthorhombic_class)
    {
      a6fCell[3] = 90.0;
      a6fCell[4] = 90.0;
      a6fCell[5] = 90.0;
    }
  else if (eClass == eSpacegroup_tetragonal_class)
    {
      a6fCell[3] = 90.0;
      a6fCell[4] = 90.0;
      a6fCell[5] = 90.0;
      f0 = 0.5 * (a6fCell[0] + a6fCell[1]);
      a6fCell[0] = a6fCell[1] = f0;
    }
  else if (   (eClass == eSpacegroup_trigonal_class)
           || (eClass == eSpacegroup_hexagonal_class) )
    {
      a6fCell[3] = 90.;
      a6fCell[4] = 90.;
      a6fCell[5] = 120.;
      f0 = 0.5 * (a6fCell[0] + a6fCell[1]);
      a6fCell[0] = a6fCell[1] = f0;
    }
  else if (eClass == eSpacegroup_cubic_class)
    {
      a6fCell[3] = 90.0;
      a6fCell[4] = 90.0;
      a6fCell[5] = 90.0;
      f0 = (a6fCell[0] + a6fCell[1] + a6fCell[2]) / 3.0;
      a6fCell[0] = a6fCell[1] = a6fCell[2] = f0;
    }
  else if (eClass == eSpacegroup_rhombohedral_class)
    {
      f0         = (a6fCell[0] + a6fCell[1] + a6fCell[2]) / 3.0;
      a6fCell[0] = a6fCell[1] = a6fCell[2] = f0;
      f0         = (a6fCell[3] + a6fCell[4] + a6fCell[5]) / 3.0;
      a6fCell[3] = a6fCell[4] = a6fCell[5] = f0;
    }
  return (0);
}


/*  Notes on fields read in from the header.


CCD_DETECTOR_VECTORS=1.0000 0.0000 0.0000 0.0000 1.0000 0.0000;

    Should always be fixed to these values.

CCD_GONIO_NAMES=RotZ RotX/2Theta RotY TransX TransY TransZ/Distance;
CCD_GONIO_NUM_VALUES=6;
CCD_GONIO_UNITS=deg deg deg mm mm mm;
CCD_GONIO_VALUES=0.0000 0.1207 4.0000 -0.0666 -0.0014 26.2035;
CCD_GONIO_VECTORS=0.000 0.000 1.000 1.000 0.000 0.000 0.000 1.000 0.000 1.000 0.000 0.000 0.000 1.000 0.000 0.000 0.000 -1.000;

    These apply to the Cgoniometer object attached to the Cdetector object.
    Rotz applied first, then 2theta swing then rotation about y.
    Rotations applied to translation vector.


CCD_SPATIAL_BEAM_POSITION=269.9915 267.0191;
CCD_SPATIAL_DISTORTION_INFO=269.9915 267.0191 0.1367 0.1367;

    Distortion info overrides any settings in CCD_SPATIAL_BEAM_POSITION.
    This is read in by the Cspatial class and is used to define the pixel values
    which map to the tip of the translation vector.

CCD_SPATIAL_DISTORTION_TYPE=Simple_spatial;

    Indicates that no distortion infor will be used.

CCD_SPATIAL_DISTORTION_VECTORS=0.0000 -1.0000 1.0000 0.0000;

    2x2 PIXEL transformation.
    These are used to get any one of 8 possible raster strategies to match the
    CCD_DETECTOR_VECTORS above.  These vary from detector to detector.


CRYSTAL_DESCRIPTION=unknown;
CRYSTAL_GONIO_DESCRIPTION=AFC8: Eulerian 3-circle;
CRYSTAL_GONIO_NAMES= Omega Chi Phi;
CRYSTAL_GONIO_NUM_VALUES=3;
CRYSTAL_GONIO_UNITS= deg deg deg;
CRYSTAL_GONIO_VALUES=0.0000 45.0000 0.0000;
CRYSTAL_GONIO_VECTORS=1.0000 0.0000 0.0000 0.0000 -1.0000 0.0000 1.0000 0.0000 0.0000;

      Depending upon the type of detector, these might differ.  The first vector is
      almost always an omega rotation about x.  The second might be a chi rotation
      about z if the circle is perpendicular to the beam vector, or (as in this case),
      a chi rotation about the y axis.
 
CRYSTAL_MOSAICITY=2.8924;
CRYSTAL_MOSAICITY_SIGMA=0.1446;
CRYSTAL_ORIENT_ANGLES=-39.8366 -14.5992 107.8116;
CRYSTAL_ORIENT_ANGLES_SIGMA=0.2660 0.1478 0.1447;
CRYSTAL_SPACEGROUP=5;
CRYSTAL_UNIT_CELL=10.9724 20.1750 12.0593 90.0000 95.0470 90.0000;
CRYSTAL_UNIT_CELL_SIGMA=0.03614 0.13007 0.06044 0.00000 0.13442 0.00000;

     Crystal rotation angles are relative to the standard basis.

ROTATION=-90.0000 90.0000 0.5000 15.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
ROTATION_AXIS_NAME=Omega;
ROTATION_VECTOR=1.0000 0.0000 0.0000;

     The rotation vector might change for a given experiment.  Per-image rotation information
     is contained in the ROTATION vector.  We have START,STOP and STEP values along with the
     number of seconds scanned.  It should be noted that the values for START and STOP should
     be on the order of 3 degrees.


SCAN_ROTATION=-90.0000 90.0000 0.5000 15.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;
SCAN_ROTATION_AXIS_NAME=Omega;
SCAN_ROTATION_VECTOR=1.0000 0.0000 0.0000;

     Same as ROTATION, but contains information for the entire scan.  Used by predict.
     The first 3 elements are Start Finish and Step.  Given an image sequence number, uses
     the information to obtain the rotation.

*/

#include <string.h>

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::cerr;
#endif


int Crefine::nUpdateDisplay()
{
#if ( !defined(SSI_PC) && !defined(NO_X_WINDOWS) )
    if( NULL != m_poXprop )
    {
      // Tell dtdisplay to display the list
      m_poReflnlist->vSetOverwrite(TRUE);
      m_poReflnlist->nWrite(sDtrekGetPrefix() + "dtrefinetmp.ref");
      // Calling program will do this by calling nNotifyDisplay.
      //      m_poXprop->hSetProperty("DTDISPLAY_REFLN_UPDATE",
      //                              sGetCWD() + sDtrekGetPrefix() + "dtrefinetmp.ref");
    }
#else
    if( 0 != m_nDisplay )
    {
      // Do not append a version number when writing the following file
      m_poReflnlist->vSetOverwrite(TRUE);
      m_poReflnlist->nWrite(sDtrekGetPrefix() + "dtrefinetmp.ref");

#ifdef SSI_PC
      CCrclHelper::GetInstance()->vSendRefineUpdateDisplay("dtrefinetmp.ref");
#endif
    }
#endif
    
    return 0;
}

Crefine::eTransTypes Crefine::m_aeTransWeighting[] = {
      SOURCE_WAVE,
      SOURCE_ROT1,
      SOURCE_ROT2,
      DET_ROT1,
      DET_ROT2,
      DET_ROT3,
      DET_TRANS1,
      DET_TRANS2,
      DET_TRANS3,
      ARRAYFLAG
};



int Crefine::nComputeAbMats(double* ppfDeriv[nMAXPARAMS],double* pfVal,float* pfSigs,int anRefineFlags[nMAXPARAMS],int* pnRejectFlag,int nNumParams,int nNumRefs,double* afA,double* afB,int* pnNumParamsInc) {
    // We have all derivatives so that we can compute the normal equations.
    // Although afA and afB will accept up to nMAXPARAMS, we are not refining all
    // parameters.

    int nx,ny,nz;
    int nParamCt0,nParamCt1;

    // Some variables have been excluded.
    // nNumParamsInc will count each valid parameter.

    *pnNumParamsInc = 0;
    for (nx = 0; nx < nNumParams; nx++) {
        if (anRefineFlags[nx]==1)
            (*pnNumParamsInc)++;
    };

    // Build fA Matrix and afB Vector
    for (nx = 0, nParamCt0 = 0; nx < nNumParams; nx++) {
        if (anRefineFlags[nx]==1)
        {
            afB[nParamCt0]= 0.0;
            for (nz= 0; nz<nNumRefs; nz++) {
                if (!pnRejectFlag[nz])
                {
                    float fD;

                    fD = ppfDeriv[nx][nz];
                    afB[nParamCt0]+=(-pfVal[nz]*fD)  /(max(m_fLowestSig,pfSigs[nz])*max(m_fLowestSig,pfSigs[nz])) ;
                };
            };

            // Compute fA.
            for (ny= 0,nParamCt1= 0; ny<=nx; ny++) if (anRefineFlags[ny]==1) {
                afA[nParamCt0+ *pnNumParamsInc *nParamCt1]= 0;
                if (1) {
                    for (nz= 0; nz<nNumRefs; nz++) if (!pnRejectFlag[nz]) {
                        float fD[2];
                        
                        fD[0]=ppfDeriv[nx][nz];
                        fD[1]=ppfDeriv[ny][nz];
                        afA[nParamCt0+*pnNumParamsInc*nParamCt1]+=
                            fD[0]*fD[1]/((max(m_fLowestSig,pfSigs[nz])*max(m_fLowestSig,pfSigs[nz])) );
                    };
                };
                nParamCt1++;
            };
            nParamCt0++;
        };
    };
    // Fill in symmetric elements of matrix.

    for (nx=1; nx<*pnNumParamsInc; nx++) {
        for (ny= 0; ny<nx; ny++)
            afA[ny+*pnNumParamsInc*nx]=
            afA[nx+*pnNumParamsInc*ny];
    };
    
    return 0;
};


int
Crefine::nOpLevenbergMarquardt(bool bComputeSigmas,
                               int nNumParams,
                               int nNumRefs,
                               double *ppfDeriv[nMAXPARAMS],
                               double *pfVal,
                               float *pfSigs,
                               int   *pnRejectFlag,
                               float afDelta[nMAXPARAMS],
                               float afSigma[nMAXPARAMS]
                               )
{
    int     nx,ny,nz;
    float   f0;
    double  fScale;

    double afScale[nMAXPARAMS];
    double afA[nMAXPARAMS*nMAXPARAMS];
    double afAInv[nMAXPARAMS*nMAXPARAMS];
    double afB[nMAXPARAMS];
    double afX[nMAXPARAMS];
    double fSumWeightsDiffSq;
    int    nDiffCount;


    double afDeriv[nMAXPARAMS];          // Computation of the derivative w.r.t. each parameter.
    int   nNumParamsInc;                 // Number of included parameters.  This parameter is returned by nComputeAbMats()



    static Cstat oStat;


    // Compute the scale vactor for each variable.
    // In this loop: nx Counts the parameter, ny counts the reflection.

    // Calculate derivatives (for printing mainly)
    // Also, calculate the scale factor.
    for (nx = 0; nx < nNumParams; nx++) {
        if (m_anRefineFlags[nx]==1)
        {
            afDeriv[nx] = 0.0;
            fScale = 0.0;
            for (ny=0;ny<nNumRefs;ny++) if (!pnRejectFlag[ny])
            {
                afDeriv[nx] += ppfDeriv[nx][ny];
                fScale += fabs(ppfDeriv[nx][ny]);
            }
            fScale/=max(1,nNumRefs);
            afScale[nx] = fScale;
         }
    };

    // Apply scale values to the data.  These must be
    // taken into account when we get sigma and delta values.
    for (nx = 0; nx < nNumParams; nx++) {
        if ((m_anRefineFlags[nx]==1) && (afScale[nx]))
        {
            for (ny=0;ny<nNumRefs;ny++) if (!pnRejectFlag[ny])
            {
                ppfDeriv[nx][ny]/=afScale[nx];
            }
        }
    };

    // Compute sum of deviations in the input.
    fSumWeightsDiffSq = 0.0;
    nDiffCount = 0;
    oStat.vClear();
    for (nx=0,ny=0;nx<nNumRefs;nx++) {
        if (!pnRejectFlag[nx])
        {
            oStat.vAdd(pfSigs[nx]);
            fSumWeightsDiffSq += pfVal[nx]*pfVal[nx]/(pfSigs[nx]*pfSigs[nx]);
            nDiffCount++;
        }
    };
    fSumWeightsDiffSq/=nDiffCount;

    // Compute a standard deviation.
    // We do this so that reflections with VERY low sigmas are not weighted
    // disproportionally high.  Reflections with high sigmas will just
    // get ignored which is okay.


    
    m_fLowestSig = max(0.0,oStat.fAverage()- oStat.fStandardDeviation()*3.0);
    /*
    if (m_fLowestSig == 0.0)
        m_fLowestSig = oStat.fAverage();
    */



    // Print the variable derivatives.
    if (m_nVerbose > 2)
    {
        printf("\n---------------Derivatives--------------\n");
        for (nx = 0; nx < nNumParams; nx++)
            if (m_anRefineFlags[nx]==1)
            {
                printf("%-20s %6.4f\n",m_pcVarNames[nx], afDeriv[nx]);
            }
            printf("------------------------------------------\n");
    }



    // Compute the sigmas.  This must be done before distorting the matrix
    // with the lambda factors.

    for (nx=0;nx<nMAXPARAMS;nx++)
      m_anSigmaFlags[nx] = m_anRefineFlags[nx];

    if (bComputeSigmas) {
        int nRejectCt = 0;
        bool bVariablesRemoved;

        // Initially, we use all variables for computing the sigma values.
        for (nx=0;nx<nNumParams;nx++)
            m_anSigmaFlags[nx] = m_anRefineFlags[nx];

        bVariablesRemoved = TRUE;
        while ((bVariablesRemoved) && (m_nAutoCorrCheck)) {

            // Compute the A and b matrix and vector.
            nComputeAbMats(ppfDeriv,pfVal,pfSigs,m_anSigmaFlags,pnRejectFlag,nNumParams,nNumRefs,afA,afB,&nNumParamsInc);


            if (0 == nInvMatND_svd(nNumParamsInc, afA, afAInv)) {

                // On the first pass, we reject variables that were highly correlated.  These are found by
                // investigating the covariance matrix.
                // We are looking for entries that have a variance nearly equal to the covariance with some other
                // variable: i.e., we want to find A and B such that Cov(A,B) = sqrt(Var(A))*sqrt(Var(B))
                // Another way of looking at this, is we want a high correlation coeff. between the two variables.
                // Once the highly correlated pairs are found, we reject the variable with the lowest weighting.


                bVariablesRemoved = FALSE;
                for (nx=0;(bVariablesRemoved==FALSE) && (nx<nNumParamsInc);nx++) {
                    for (ny=nx+1;(bVariablesRemoved==FALSE) && (ny<nNumParamsInc);ny++) {
                        f0 = afAInv[nx+nNumParamsInc*ny]
                    / (sqrt(fabs( afAInv[nx+nNumParamsInc*nx] * afAInv[ny+nNumParamsInc*ny])));
                        if (fabs(f0)>m_fMaxCorrCoeff) {
                            int nRealx,nRealy;
                            int nWeightx,nWeighty;
                            int nReject;
                            int nNotReject;
                            // Get the "real" variable numbers.
                            vUsed2RealVarNum(nx,&nRealx);
                            vUsed2RealVarNum(ny,&nRealy);
                            // Look up both of the variables in the table.  Assign a weight (100 means we did not find the variable).
                            for (nWeightx=0;(m_aeTransWeighting[nWeightx]!=ARRAYFLAG) && ((int)m_aeTransWeighting[nWeightx]!=nRealx);nWeightx++);
                            if (m_aeTransWeighting[nWeightx]==ARRAYFLAG)
                                nWeightx = 100;
                            for (nWeighty=0;(m_aeTransWeighting[nWeighty]!=ARRAYFLAG) && ((int)m_aeTransWeighting[nWeighty]!=nRealy);nWeighty++);
                            if (m_aeTransWeighting[nWeighty]==ARRAYFLAG)
                                nWeighty = 100;
                            // If one of the variables was found in the weighting table AND both variables are turned ON, then we must remove
                            // the variable with the higher weighting. (that is, further up in the table).
                            if (((nWeightx!=100) || (nWeighty!=100)) && (m_anSigmaFlags[nRealx]==1) && (m_anSigmaFlags[nRealy]==1)) {
                                if (nWeightx>nWeighty) {
                                    nReject = nRealy;
                                    nNotReject = nRealx;
                                } else {
                                    nReject = nRealx;
                                    nNotReject = nRealy;
                                };
                                if (2 < m_nVerbose) {
                                  if (1.0 < fabs(f0))
                                    printf("WARNING: |correlation| > 1 : %5.3f\n", f0);
                                };
                                if (-1.0 > f0)
                                  f0 = -1.0;
                                else if (1.0 < f0)
                                  f0 = 1.0;
                                if (1 < m_nVerbose)
                                  {
                                    printf("INFO: '%s' sigma estimated, "
                                           "correlation of %5.3f with '%s'\n",
                                           m_psVarNames[nReject].string(), f0, m_psVarNames[nNotReject].string());
                                    fflush(stdout);
                                  }
                                m_anSigmaFlags[nReject] = 0;
                                nRejectCt++;
                                bVariablesRemoved = TRUE;
                              }
                          }
                      };
                  };
            } else {
                for (nx = 0; nx<nNumParams; nx++)
                    m_anSigmaFlags[nx] = 0;
                bVariablesRemoved = FALSE;
            };
        }

        int nNumParamsIncBefore;
        int nNumParamsIncAfter;

        nNumParamsIncAfter = nNumParamsInc;
                nComputeAbMats(ppfDeriv,pfVal,pfSigs,m_anRefineFlags,pnRejectFlag,nNumParams,nNumRefs,afA,afB,&nNumParamsIncBefore);

        for (nx = 0,ny= 0, nz = 0; nx<nNumParams; nx++) {
            if (m_anSigmaFlags[nx]==1) {
                afSigma[nx]=sqrt(fabs(afAInv[ny*nNumParamsIncAfter+ny])*fSumWeightsDiffSq);
                ny++;
            } else
                afSigma[nx]=sqrt(1.0/fabs(afA[nz*nNumParamsIncBefore+nz])*fSumWeightsDiffSq);
            if (m_anRefineFlags[nx]==1)
                nz++;
        };
    }

    nComputeAbMats(ppfDeriv,pfVal,pfSigs,m_anRefineFlags,pnRejectFlag,
                   nNumParams,nNumRefs,afA,afB,&nNumParamsInc);



    // Check for m_fLambda==-1.

    if (m_fLambda < 0)
      {
        m_fLambda = 10;
      }

  // Add in the lambda factor.

  for (nx = 0; nx < nNumParamsInc; nx++)
    afA[nx +  nNumParamsInc * nx]*=(1 + m_fLambda);


    if (m_nVerbose>3)
      {
        printf("-------------------- Matrix Elements ------------------------\n");
        for (nx = 0; nx < nNumParamsInc; nx++)
          {
            printf("[ ");
            for (ny = 0; ny < nNumParamsInc; ny++)
              printf("%f ", afA[ny+nNumParamsInc*nx]);
            printf(" ]   [%f]\n", afB[nx]);
          }
        printf("-------------------------------------------------------------\n");
      }

    // Solve the system.
    
    if (0 == nSolveAxB_svd(nNumParamsInc, afA, afX, afB))
      {
        if (3 < m_nVerbose)
          printf("Matrix inverted successfully.\n");
      }
    else
      {
        if (3 < m_nVerbose)
          printf("WARNING: Matrix inversion FAILED!\n");
      }


    // Reinstate variables in the afDelta array.
    // nx counts the parameter.
    // ny counts the refined parameter.

    for (nx = 0, ny = 0; nx < nNumParams; nx++)
      if (m_anRefineFlags[nx]==1)
        {
          afDelta[nx]=afX[ny++];
        }
      else
        afDelta[nx]= 0.0;

  // Print the variable deltas.
      
     if (m_nVerbose>2)
     {
          printf("\n-----------------Deltas------------------\n");
          for (nx = 0; nx<nNumParams; nx++)
              if (m_anRefineFlags[nx]==1)
              {
                  printf("%-20s %6.4f\n",m_pcVarNames[nx], afDelta[nx]);
              }
              printf("-------------------------------------------\n");
      }

      
      // Unapply the scale factors.
      for (nx = 0; nx<nNumParams; nx++) {
          if (m_anRefineFlags[nx]==1)
          {
              afDelta[nx]/=afScale[nx];
              afSigma[nx]/=afScale[nx];
          };
      };
    return 0;
}

/*  Equation to minimize:

    ResidVec = (R*s0 + (DD*Pos + Rot*Gonio*RecipShift)/Length(DD*Pos + Rot*Gonio*RecipShift) - Rot*Gonio*T*U*B*HKL)
    Min(Trans(ResidVec)*Trans(P)*P*ResidVec)

    R     3x3 = Source Rotation 
    s0    3   = Source Vector (unrotated)
    DD    3x3 = Detector matrix (with rotation matrix applied)
    Pos   3   = MM pixel coordinates.
    Rot   3x3 = Scan rotation matrix
    Gonio 3x3 = Crystal Goniometer matrix
        T     3x3 = Crystal Twin Law.
    U     3x3 = Crystal rotation matrix (missetting)
    B     3x3 = Crystal B matrix
    HKL   3   = HKL coordinates of reflection.
    P     3x3 = Projection matrix.

*/


int Crefine::nRefineLoop(const int nLoopIn)
{

  int   nStat;                  // General error status
  int   nx,ny;
  double fDet;                   // General purpose determinant variable.
  double f0, f1;
  double a3x3fTempMat1[3][3];
  double a3x3fTempMat2[3][3];
  double a3fTempVec1[3];
  double a3fTempVec2[3];
  double a3fTempVec3[3];
  double a3fTempVec4[3];
  float a3x3fFloatBuf[3][3];
  float a3fFloatBuf[3];
  double a6fTemp[6];
  double a3fTemp[3];

  // General variables
  int   nVar;
  int   nNumReflns;             // Total number of reflections in the list
  int   nRef;                   // Counter for reflection loop.
  int   nNumParams;             // Total Number of parameters to refine.
  Crefln *poRefln;              // Used to point to a reflection.

  // Reflection rejection varaibles.

  int  *pnRejectFlag=NULL;      // Array of rejection reasons.
  float fMaxResStar;            // Maximum resolution.
  float fMinResStar;            // Minimum resolution.
  int   nNumRejects;            // The number of reject reflections (in this twin comp)
  int   nNumIgnored;            // The number of ignored reflections (in this twin comp)
  int   nNumMMRejects;          // The number of rejected because MM position was off.
  int   nNumDegRejects;         // The number of rejected because PHI centroid was off.
  int   nNumOtherComp;          // The number of reflections in another twin component.
  int   nNumMosaic;             // Number of reflections with width information.
  bool  bHKLChanged;            // Did one of the reflections change rounded HKL value?
  bool   bOneGoodRefln;                 // Used to detect the presence of a 'predicted' input list.

  // Source variables.

  double a3fS0[3];               // Normalized source vector after rotations.
  float fSpectralDispersion;    // delta(lambda)/lambda.
  double a2fOrigS0Rot[2];        // Original Source Vector rotations (before tweaking)
  double a2x3x3fS0Deriv[2][3][3];// Derivatives.

  // Crystal variables.

  double a3x3fCrysRot[3][3];                                                    // U
  double a3x3fCrysRotInv[3][3];                                                 // Inverse of U
  double a3x3fB[3][3];                                                                  // B
  double a3x3fBInv[3][3];                                                               // Inverse of B
  double a3fRecipVec[3];                                                                // Reciprocal shift vector rotated (actually modifies the detector portion of refinement equ.)
  double aa3x3fCrysTwinLaws[g_nMaxTwinLaws + 1][3][3];  // Twin laws (matrix form T) for this crystal.
  double aa3fOrigCrysTwinLaws[g_nMaxTwinLaws + 1][3];   // Twin laws (three rots form) for this cyrstal.
  double aa3x3fCrysOrient[g_nMaxTwinLaws + 1][3][3];     // T*U*B
  double aa3x3fCrysOrientInv[g_nMaxTwinLaws + 1][3][3];  // Inverse of T*U*B
  double aa3fOrigRecipShift[g_nMaxTwinLaws + 1][3];
  double afOrigTwinFraction[g_nMaxTwinLaws + 1];
  double a6fOrigCell[6];                                                                // Original crystal cell (before tweaking).
  double a7x3x3fCrysDeriv[7][3][3];                                             // Crystal derivatives.
  double a3x3x3fCrysRotDeriv[3][3][3];                                  // Crystal orientation derivatives.
  double a3x3x3fCrysTwinLawDeriv[g_nMaxTwinLaws + 1][3][3][3];  // Crystal twin law derivatives.
  double a3fOrigCrysRot[3];                                                             // Original crystal rotation values.
  int   nNumCrystalParams;                                                              // The number of refinement parameters for the cyrstal.  
  double fOrigSourceWave;                                                               // Original source wave length.


  // Detector Information

  double a3x3fDD[3][3];              // Only collumns [0] and [1] are used. [2] is set to the rotated translation vector (a3fTrans)
  double a3x3fDDInv[3][3];           // Inverse of Detector DD matrix
  double a3fDN[3];                   // Normal vector
  double a3fOrigRot[3];              // Rotations.
  double a3fOrigTrans[3];            // Translation
  double a3fOrigTransStandard[3];    // Translation (standard parameters)
  double a6x3x3fDDDeriv[6][3][3];    // Detector Derivatives.


  // Rotation information.

  int       nNumRotations;                          // Number of goniometer rotations.
  float     aa5fGonio[nMaxGonioRotationOffsets][5]; // Goniometer settings.
  int       anGonio[nMaxGonioRotationOffsets];      // Goniometer rotation axis.
  double    aa3fE3[nMaxGonioRotationOffsets][3];    // Rotation vector.
  double    a3fE3CrossS0[3];                        // Find inverse of the distance of fXRO onto this plane for the Lorentz factor.
  double    a3fE3[3];                               // Rotation vector (inside the reflection processing loop).
  double    fRotationWidth;                         // Used in mosaicity refinement.

  // Goniometer Information.

  double aa3x3fGonMatrix[nMaxGonioRotationOffsets][3][3];       // Crystal goniometer rotation matrix G
  double aa3x3fGonMatrixInv[nMaxGonioRotationOffsets][3][3];    // Inverse of crystal goniometer matrix G

  // Mosaicity information.

  double a3fOrigMosaicity[3];           // The original input mosaicity.
  bool   bRefineMosaicity;              // Combination of user flag and loop counter.
  bool   bMosaicityWasChanged;          // If mosaicity is not refined, then it's
                                        // sigma is set to zero at termination.
  const float fMinPercentAcceptedInMosaicity = 0.1f; // Minimum number of reflections with mosaicity information required for a mosaicity refinement.


  // Control Variables.


  double fInitLambda;             // Set lambda back to it's original value after running through loops.
  int nLoop, nLoopCt;
  int nTwinLaw;                                   // Which twin law for a given reflection?  (Nominally zero)
  float afDelta[nMAXPARAMS];      // Delta Values loaded in by the optimization routine.
  float afSigma[nMAXPARAMS];      // Sigma values for each parameter.
  double fTotalError;             // Total Error after each pass.
  double fLastTotal;              // Total for last loop. (used to detect a 'faulty' optimization pass).
  int    nLastAccept;             // Total number of reflections accepted in last loop.
  double a3fDeriv[3];             // Used when calculating derivatives.
  bool bAllParametersConverged;   // Flag indicating that all parameters have converged.
  Cstat oStatOutlierROT;          // Used to reject outliers.
  Cstat oStatOutlierMM;
  double fAverageMMResid;
  double fDevMMResid;
  double fAverageROTResid;
  double fDevROTResid;



  // Per-reflection memory arrays.

  double *ppfDeriv[nMAXPARAMS];             // Value for each reflection with a positive tweaking delta.
  float  *pfSigs;                           // Sigmas for each value.
  double  *pfVal;                           // Values.
  double *aa3x3fPhiRotMat;                  // Rotation matricies for each reflection about a3fE3.
  double *aa3fVDCO;                         // Virtual detector coordinates.
  double *aa3fResidVec;                     // Residual S vector 
  double *aa3fXProjVec;                     // projection vector used to discount effect of wide images.
  double *afResidVecLen;                    // Length of residual vector 
  float  *afObsRotMid;                  
  float  *afObsRotSigma;                
  float  *afObsRotWidth;
  float  *afObsRotStart;
  float  *afCalcRotMid;
  float  *afObsXmm;
  float  *afObsYmm;
  float  *afObsZmm;
  float  *afObsIntensity;
  float  *afObsIntensitySig;
  float  *afObsSharpness;
  int    *anFlopCount;
  int    *anGonioMatNum;
  int    *anTwinID;                                                     // Twin ID (with 1000*TwinLaw).
  int    *anCompete;
  int    *anLostCompetition;



  // Internal computational loop variables.
  
  double a3fVNTS[3];                // Virtual normalized T (i.e. S) vector.
  double a3fVTS[3];                 // Virtual T (unnormalized) vector .
  double a3x3fPhiGon[3][3];         // Combination of goniometer, and phi rotation.
  double a3x3fPhiGonInv[3][3];      // Inverse of fPhiGon
  double a3x3fCrysPhiGonInv[3][3];  // This takes X into HKL vector but WARNING:  It does not take into account reciprocal shifts!
  double a3fHKL[3];                 // Fractional HKL entries.
  int    a3nHKL[3];                 // HKL entries.
  int    a3nPrevHKL[3];             // Used to detect if any HKL's changed.
  double a3fXRO[3];                 // Observed reciprocal x-vector.
  double a3fXROR[3];                            // Rotated (out of diffracting condition) reciprocal x-vector
  double a3fXRC[3];                 // Calculated reciprocal x-vector.
  double a3fXRCR[3];                // Calculated reciprocal x-vector rotated into diffracting condition.
  double a3fXProjVector[3];

  double a3fXROXRCR[3];              // Difference of XRO and XRCR
  double a3fXROXRC[3];               // Difference of XRO and XRC
  double fPHO;                      // Observed rotation about rotation vector.
  double fRotO;                     // Observed rotation width (radians)
  
  // At end of each loop, these variables will hold the crucial data.
  
  double a3fVDCO[3];          // Observed virtual detector coordinates.
  double a3fVDCC[3];          // Calculated virtual detector coordinates.
  double a2fVDErr[2];         // Observed-Calculated error for vitual detector coordinates.
  double fRotErr;             // Rotation error.  Does not use fPHO.
  double fROTError,fMMError;
  double fLorentz;

  nNumReflns  = m_poReflnlist->nGetNumReflns();

  nLoop       = nLoopIn;
  nLoopCt     = 0;
  fInitLambda = m_fLambda;

  // Get the degrees of freedom for the crystal, and a basis set spanning the space
  // (up to 6 dimensional for triclinic)

  if ((m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet()<= 0) || (!m_anRefineFlags[ms_nCrysConstraints])) {
      
      if (m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet()<= 0)
          printf("WARNING: Spacegroup <= 0... Assuming triclinic.\n");
      nx = m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet();
      m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->vSet(1);
      m_ppoCrystals[m_nTwinID-1]->nGetDerivVecs(nNumCrystalParams, m_a6x6fCrystalVecs);
      m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->vSet(nx);

  } else
    m_ppoCrystals[m_nTwinID-1]->nGetDerivVecs(nNumCrystalParams, m_a6x6fCrystalVecs);

  nNumParams =  nMAXPARAMS;
  // Change certain parameters to have a flag of 2, indicating that they are refined seperately or not at all.
  if (m_anRefineFlags[ms_nCrysMosaicity])
          m_anRefineFlags[ms_nCrysMosaicity] = 2;
  if (m_anRefineFlags[ms_nCrysConstraints])
      m_anRefineFlags[ms_nCrysConstraints] = 2;

  // This is only used when rounding indices.

  m_bIsRhomboAsHex =  (
             (146 == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet())
          || (148 == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet())
          || (155 == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet())
          || (160 == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet())
          || (161 == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet())
          || (166 == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet())
          || (167 == m_ppoCrystals[m_nTwinID-1]->m_poSpacegroup->nGet()));
  fRotationWidth = m_poRotation->fGetIncrement();
  if (fRotationWidth<=0.0)
      fRotationWidth = 0.5;

  // Initialize the delta arrays and some other arrays.
  
  for (nx = 0; nx < nNumParams; nx++)
  {
      if (m_anRefineFlags[nx]==1) {
          // Only allocate for memory if we are actually refining this variable.
          ppfDeriv[nx] = new double[nNumReflns];
          for (ny = 0; ny < nNumReflns; ny++)
              ppfDeriv[nx][ny] = 0.0f;
      } else
          ppfDeriv[nx] = NULL;
      afDelta[nx]        = 0.0f;
      afSigma[nx]        = 0.0f;
  }

  fLastTotal = -1;
  nLastAccept = 0;

  m_nProfileCt = 0;

  bRefineMosaicity = ( 0 != m_anRefineFlags[ms_nCrysMosaicity]);
  bMosaicityWasChanged = FALSE;


  // Set up other per-refleciton arrays.
  pnRejectFlag =        new int   [nNumReflns];
  anFlopCount  =        new int   [nNumReflns];
  anGonioMatNum=        new int   [nNumReflns];
  anTwinID     =                new int   [nNumReflns];
  anCompete    =        new int   [nNumReflns];
  anLostCompetition =   new int   [nNumReflns];
  pfSigs       =        new float [nNumReflns];
  pfVal        =        new double[nNumReflns];
  aa3x3fPhiRotMat =     new double [nNumReflns*9]; 
  aa3fVDCO        =     new double [nNumReflns*3];
  aa3fResidVec    =     new double [nNumReflns*3];
  aa3fXProjVec    =     new double [nNumReflns*3];
  afResidVecLen   =     new double [nNumReflns];
  afObsRotMid   =       new float[nNumReflns];
  afObsRotSigma =       new float[nNumReflns];
  afObsRotWidth =       new float[nNumReflns];
  afObsRotStart =       new float[nNumReflns];
  afCalcRotMid  =       new float[nNumReflns];  
  afObsXmm      =       new float[nNumReflns];
  afObsYmm      =       new float[nNumReflns];
  afObsZmm      =       new float[nNumReflns];
  afObsIntensity=       new float[nNumReflns];
  afObsIntensitySig =   new float[nNumReflns];
  afObsSharpness    =   new float[nNumReflns];
  

  // Figure out the number of distinct goniometer settings we have in out reflection list.
  nNumRotations = 0;
  nStat = 0;
  for (nRef=0;nRef<nNumReflns;nRef++)
  {
      double a5fGonioSet[5];
      int    nGonioValues;
      int    nRotAxis;
      
      nGonioValues = 0;
      nRotAxis = -1;
      if (m_poReflnlist->m_nFI_fGonio1>=0)
          a5fGonioSet[nGonioValues++] = (*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio1);
      if (m_poReflnlist->m_nFI_fGonio2>=0)
          a5fGonioSet[nGonioValues++] = (*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio2);
      if (m_poReflnlist->m_nFI_fGonio3>=0)
          a5fGonioSet[nGonioValues++] = (*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio3);
      if (m_poReflnlist->m_nFI_fGonio4>=0)
          a5fGonioSet[nGonioValues++] = (*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio4);
      if (m_poReflnlist->m_nFI_fGonio5>=0)
          a5fGonioSet[nGonioValues++] = (*m_poReflnlist)[nRef].fGetField(m_poReflnlist->m_nFI_fGonio5);
      if (m_poReflnlist->m_nFI_nGonioRotAxis>=0)
          nRotAxis = (*m_poReflnlist)[nRef].nGetField(m_poReflnlist->m_nFI_nGonioRotAxis) - 1;
      
      if ((!nGonioValues) || (nRotAxis<0)) {
          // Use the default goniometer values.
          for (nx=0;nx<m_poCrysGonio->m_nNumRotValues;nx++)
              a5fGonioSet[nx] = m_poCrysGonio->m_pfDatumValue[nx];
          nRotAxis = m_poCrysGonio->nGetNum(m_poRotation->sGetName());
      } 
      // Check and see if this goniomter constant is in use.
      for (ny = 0; ny < nNumRotations; ny++) {
          for (nx=0;nx<m_poCrysGonio->m_nNumRotValues;nx++) {
              if (ABS(a5fGonioSet[nx] - aa5fGonio[ny][nx])>0.01)
                  break;
          };
          // Did we find a match?  If so, then break out.
          if (nx == m_poCrysGonio->m_nNumRotValues)
              break;
      };
      if (ny == nNumRotations) {
          // We did not find a match.  Create a new entry.
          for (nx=0;nx<m_poCrysGonio->m_nNumRotValues;nx++) {
              m_poCrysGonio->m_pfDatumValue[nx] = a5fGonioSet[nx];
              // Always save nominal values.
              aa5fGonio[ny][nx] = a5fGonioSet[nx];
          };
          anGonio[ny] = nRotAxis;

          m_poCrysGonio->vCalcGetRotMatrix(& a3x3fFloatBuf[0][0],nRotAxis); 
          vCopyMat3D(&a3x3fFloatBuf[0][0],&aa3x3fGonMatrix[ny][0][0]);
          fDet  = fInvMat3D(&aa3x3fGonMatrix[ny][0][0], &aa3x3fGonMatrixInv[ny][0][0]);
          if (m_poCrysGonio->nGetRotVector(nRotAxis,&a3fFloatBuf[0])) {
              nStat =1;
              break;
          };
          vCopyVec3D(&a3fFloatBuf[0],&aa3fE3[ny][0]);
          
          nNumRotations++;
      };
      anGonioMatNum[nRef] = ny;
          if (m_poReflnlist->m_nFI_nTwinID>=0) 
                  anTwinID[nRef] = (*m_poReflnlist)[nRef].nGetField(m_poReflnlist->m_nFI_nTwinID);     
          else
                  anTwinID[nRef] = 1;

  };
  
  // Save the observed rotation midpoint, and observed rotation sigma (if available).
  // We will use the afXXXX[] arrays instead of changing values in the reflection list itself.
  bOneGoodRefln = false;
  for (nx=0;nx<nNumReflns;nx++) {
      anFlopCount[nx] = 0;
      afObsRotMid[nx] = (*m_poReflnlist)[nx].fGetField(m_nFI2);
      afObsRotWidth[nx] = (*m_poReflnlist)[nx].fGetField(m_nFI3);
      if (afObsRotWidth[nx]<0.0)
          afObsRotWidth[nx] = fRotationWidth;
      afObsRotStart[nx] = -999;  
      afCalcRotMid[nx] = -999;
      if (m_nFI4>=0) 
          afObsRotSigma[nx] = (*m_poReflnlist)[nx].fGetField(m_nFI4);
      else
          afObsRotSigma[nx] = -999;
      
      afObsXmm[nx] = (*m_poReflnlist)[nx].fGetField(m_nFIMM0);
      afObsYmm[nx] = (*m_poReflnlist)[nx].fGetField(m_nFIMM1);
      afObsZmm[nx] = (*m_poReflnlist)[nx].fGetField(m_nFIMM2);
      afObsIntensity[nx] = (*m_poReflnlist)[nx].fGetIntensity();
      afObsIntensitySig[nx] = (*m_poReflnlist)[nx].fGetSigmaI();
          if ((afObsIntensity[nx]>0.0) && (afObsIntensitySig[nx]>0.0))
                  bOneGoodRefln = true;
      afObsSharpness[nx] = (m_poReflnlist->m_nFI_fObsSharpness>=0)?((*m_poReflnlist)[nx].fGetField(m_poReflnlist->m_nFI_fObsSharpness)):0.0;
      
      if (m_poReflnlist->m_nFI_nCompeteRefNum>=0)
          anCompete[nx] = (*m_poReflnlist)[nx].nGetField(m_poReflnlist->m_nFI_nCompeteRefNum);
      else
          anCompete[nx] = -1;
      anLostCompetition[nx] = 0;
      afResidVecLen[nx] = 0.0;
  };

  if (!bOneGoodRefln) {
          // We have a predicted list.
          // Go ahead and update the reflection fields to 'fool' the refinement into using everything.
          for (nx=0;nx<nNumReflns;nx++) {
                  afObsIntensity[nx] = 10.0 + 1.0;
                  afObsIntensitySig[nx] = 10.0/m_fSigma;
          };
  };

  // Load the original values of most variables.

  /*1*/
  a2fOrigS0Rot[0] = m_poSource->fGetRotation(0);
  a2fOrigS0Rot[1] = m_poSource->fGetRotation(1);
  fOrigSourceWave = m_poSource->fGetWavelength();
  /*2*/
  m_ppoCrystals[m_nTwinID-1]->vGetCell(a6fTemp);
  if (m_anRefineFlags[ms_nCrysConstraints])
      nForceCell(a6fTemp);              // Force the cell to match the spacegroup.
  m_ppoCrystals[m_nTwinID-1]->vSetCell(a6fTemp);
  m_ppoCrystals[m_nTwinID-1]->vGetRecipCell(a6fOrigCell);  // We are refining reciprocal space parameters.
  /*3*/
  m_ppoCrystals[m_nTwinID-1]->nCalcRotMatrix();
  m_ppoCrystals[m_nTwinID-1]->vGetRotMatrix(&a3x3fTempMat1[0][0]);
  a3fOrigCrysRot[0] = m_ppoCrystals[m_nTwinID-1]->fGetOrientAngle(0);
  a3fOrigCrysRot[1] = m_ppoCrystals[m_nTwinID-1]->fGetOrientAngle(1);
  a3fOrigCrysRot[2] = m_ppoCrystals[m_nTwinID-1]->fGetOrientAngle(2);
  vDeconvMat3D3XYZ(a3x3fTempMat1, &a3fOrigCrysRot[0], &a3fOrigCrysRot[1], &a3fOrigCrysRot[2]);
  for (nTwinLaw = 0; nTwinLaw <= m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws;nTwinLaw++) {
          m_ppoCrystals[m_nTwinID-1]->nGetTwinLaw(nTwinLaw,&aa3fOrigCrysTwinLaws[nTwinLaw][0],false);
          m_ppoCrystals[m_nTwinID-1]->nGetRecipShift(nTwinLaw,&aa3fOrigRecipShift[nTwinLaw][0]);
          m_ppoCrystals[m_nTwinID-1]->nGetTwinFraction(nTwinLaw,&afOrigTwinFraction[nTwinLaw]);
  };  
  /*4*/
  m_ppoDetector[0]->nGetModifiedDetectorParams(a6fTemp);
  vCopyVec3D(a6fTemp,a3fOrigRot);
  vCopyVec3D(&a6fTemp[3],a3fOrigTrans);
  m_ppoDetector[0]->nGetStandardDetectorParams(a6fTemp);
  vCopyVec3D(&a6fTemp[3],a3fOrigTransStandard);
  /*5*/
  m_ppoCrystals[m_nTwinID-1]->vGetMosaicity(a3fOrigMosaicity);


  // The OrigXXXXX variables in this routine are only as "original" as the last loop iteration
  // When we get better values, they are updated.  Thus, we must save some values in the
  // Shift variables (located in Crefine).  At the bottom when we are finished, we will update
  // the shift variables to reflect our changes.

  /*1*/
  vCopyVecND(2, &a2fOrigS0Rot[0], & m_a2fSourceRotShifts[0]);
  m_fSourceWaveShift     = fOrigSourceWave;
  /*2*/
  m_ppoCrystals[m_nTwinID-1]->vGetCell(& m_a6fCrysTurds[0]);             // Make sure to save the real space parameters here.
  /*3*/
  m_a3fCrysRotShifts[0] =  a3fOrigCrysRot[0];
  m_a3fCrysRotShifts[1] =  a3fOrigCrysRot[1];
  m_a3fCrysRotShifts[2] =  a3fOrigCrysRot[2];
  vCopyVecND(3*(1+m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws),&aa3fOrigCrysTwinLaws[0][0],&m_aa3fCrysTwinLawShift[0][0]);
  vCopyVecND(3*(1+m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws),&aa3fOrigRecipShift[0][0],&m_aa3fCrysRecipShiftShift[0][0]);
  vCopyVecND(1+m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws,&afOrigTwinFraction[0],&m_afCrysTwinFraction[0]);
  /*4*/
  m_ppoDetector[0]->nGetStandardDetectorParams(a6fTemp);  // Make sure to save the standard parameters, here.
  vCopyVec3D(&a6fTemp[3], m_a10x3fDetTransShifts[0]);
  vCopyVec3D(&a6fTemp[0], m_a10x3fDetRotShifts[0]);
  /*5*/
  vCopyVec3D(a3fOrigMosaicity,m_a3fCrysMosaicityShift);


  // Resolution information.
  if (0.0  ==  m_fResolutionMin) {
      fMinResStar = 9999999.0;
  } else {
      fMinResStar = m_poSource->fGetWavelength() / m_fResolutionMin;
  }
  if (0.0 == m_fResolutionMax)
  {
      fMaxResStar = 9999999.0;
  } else {
      fMaxResStar = m_poSource->fGetWavelength() / m_fResolutionMax;
  }

  // Clear all rejection flags.
  for (nx = 0; nx < nNumReflns; nx++)
      pnRejectFlag[nx] = 0;
  bHKLChanged = FALSE;
  fAverageROTResid = 0.0;
  fDevROTResid = 0.0;
  fAverageMMResid = 0.0;
  fDevMMResid = 0.0;

  // Main loop, Each pass is one optimization loop.
  do
  {
      nLoopCt++;

      if (3 < m_nVerbose)
      {
          printf("LOOP num: %d\n", nLoopCt);
      }
      else if (nLoopCt < nLoop) 
      {
          if (bHKLChanged)
              printf("*");
          else
              printf(".");
          if (0 == (nLoopCt % 40))
              printf("\n");
          fflush(stdout);
      }
      else
          printf("\n");


      int nPass;

      
      for (nPass=0;nPass<2;nPass++) {

          // This loop will possibly go around twice.  The second call will be necessary if
          // the variables fail to give better values than last time.
          // It should only get called once on the first pass.

          // Source information.

          m_poSource->vCalcGetS0(&a3fS0[0],&a2x3x3fS0Deriv[0][0][0]);            // Get normalized source dir. vector
          fSpectralDispersion = m_poSource->fGetSpectralDispersion();

          // Crystal information.

          nStat = m_ppoCrystals[m_nTwinID-1]->nCalcOrientMatrix(m_poSource->fGetWavelength());    // Crystal orientation matrix with miss-settings.
                  for (nTwinLaw = 0; nTwinLaw <= m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws;nTwinLaw++) {
                          m_ppoCrystals[m_nTwinID-1]->vGetOrientMatrix(& a3x3fTempMat2[0][0]);
                          m_ppoCrystals[m_nTwinID-1]->nGetTwinLaw(nTwinLaw,&aa3x3fCrysTwinLaws[nTwinLaw][0][0],true);
                          vMulMat3DMat3D(aa3x3fCrysTwinLaws[nTwinLaw],a3x3fTempMat2,aa3x3fCrysOrient[nTwinLaw]);
                          fDet = fInvMat3D(&aa3x3fCrysOrient[nTwinLaw][0][0], &aa3x3fCrysOrientInv[nTwinLaw][0][0]);
                          if (0.0f == fDet)
                          {
                                  printf("ERROR in Crefine::nRefineLoop: Could not invert crystal matrix!\n");
                                  return (1);
                          }
                  };
                  for (nTwinLaw = 1; nTwinLaw <= m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws;nTwinLaw++) 
                          m_ppoCrystals[m_nTwinID-1]->nCalcTwinLawDeriv(nTwinLaw,&a3x3x3fCrysTwinLawDeriv[nTwinLaw][0]);

                  
          m_ppoCrystals[m_nTwinID-1]->vGetRotMatrix(&a3x3fCrysRot[0][0]);      // Loaded when we called nCalcOrientMatrix
          fDet = fInvMat3D(&a3x3fCrysRot[0][0], &a3x3fCrysRotInv[0][0]);
          m_ppoCrystals[m_nTwinID-1]->vGetBMatrix(&a3x3fB[0][0]);              // Loaded when we called nCalcOrientMatrix
          fDet = fInvMat3D(&a3x3fB[0][0],&a3x3fBInv[0][0]);
                  
          m_ppoCrystals[m_nTwinID-1]->nCalcGetBDeriv(a7x3x3fCrysDeriv,&nx,fOrigSourceWave);
          m_ppoCrystals[m_nTwinID-1]->nCalcGetRotDeriv(a3x3x3fCrysRotDeriv);

          // Detector information.
          (void) m_ppoDetector[0]->nCalcGetDDDN(&a3x3fDD[0][0], a3fDN);   // Get detector information and normal vector.
          fDet  = fInvMat3D(&a3x3fDD[0][0], &a3x3fDDInv[0][0]);
          m_ppoDetector[0]->nCalcGetDDDNDerivModified(a6x3x3fDDDeriv);
          m_ppoDetector[0]->vUpdateDDDN();

          // Set the error variables.

          fTotalError= 0.0;
          m_fRMSDEG = 0.0;
          m_fRMSMM = 0.0;
          m_fRMSANG = 0.0;
          m_a7nReflnNumbers[0] = 0;
          m_a7nReflnNumbers[1] = 0;
          m_a7nReflnNumbers[2] = 0;
          m_a7nReflnNumbers[3] = 0;
          m_a7nReflnNumbers[4] = 0;
          m_a7nReflnNumbers[5] = 0;
          m_a7nReflnNumbers[6] = 0;
          m_a7fReflnNumbers[0] = 0.0;
          m_a7fReflnNumbers[1] = 0.0;
          m_a7fReflnNumbers[2] = 0.0;
          m_a7fReflnNumbers[3] = 0.0;
          m_a7fReflnNumbers[4] = 0.0;
          m_a7fReflnNumbers[5] = 0.0;
          m_a7fReflnNumbers[6] = 0.0;
          
          // Initialize mosaicity variables.
          bHKLChanged = FALSE;

          /*
          Steps in refinement loop:

            1) Calculate matrix derivatives.

              A) DET_XXXXX:   Calculate derivatives from the Cdetector::nCalcGetDDDNDerivXXXXX() routines.
              B) SOURCE_ROTX: Calculate Source rotations for Csource::vCalcGetS0()
              C) CRYS_ROTX:   Calculate rotation matrix derivatives from Ccrystal::nCalcRotDeriv()
              D) CELL_XXXX:   Calculate cell matrix derivatives from Ccrystal::nCalcGetBDeriv()
              E) GONIO_MISSX: Calculate rotation matrix derivatives from vConvRotVec3DMat3D()


                  2)  Do preliminary work.

                    A)  Calculate Rejections.  These will be flagged as rejected, and will not participate in the refinement.
                    B)  Calculate reciprocal space positions (Svec) of reflections using detector plate positions and rotations.
                    This information is used extensively for derivative calculations, and should only get computed once per loop.
                    C)  Calculate Rounded HKL entries for each non-rejected reflection.  (This is stored in the reflection list).
                    D)  Calculate the phi rotation matrix for each reflection.
                    E)  Calculate the rotational and detector plate errors.  Of course, we need this for the rejections as well.
                    F)  Calculate the residual vector.  This is defined as the X vector (from the observed position) minus the calculated X vector (from the rounded HKL).
                    G)  Calculate residual vector length.

                      Thus, the reflection dependent variables are pre-computed.  These are:
                      Rejections                Storage: pnRejectFlag
                      Phi Rotation matrix       Storage: aa3x3fPhiRotMat
                      Phi Rotation matrix Deriv Storage: aa3x3fPhiRotMatD[]
                      Rounded HKL value.        Storage: m_poReflnlist
                      NonNormalize S == DD*VD   Storage: aa3fVDCO
                      Residual Vector           Storage: aa3fResidVec
                      Residual Vector Length    Storage: afResidVecLen
                                          Twin ID                                       Storage  anTwinID

                        3) Calculate reflection derivatives.

                          A) All derivatives can be computed quickly once this information is known.

          */


          nNumRejects = 0;
          nNumIgnored = 0;
          nNumMMRejects = 0;
          nNumDegRejects = 0;
          nNumOtherComp = 0;
          nNumMosaic = 0;
          // Clear all rejection flags (except redundant flags).
          for (nx = 0; nx < nNumReflns; nx++) {
              if (pnRejectFlag[nx] != 3) 
                pnRejectFlag[nx] = 0;
          };
          m_a2fResoUsed[0] = 0.0;
          m_a2fResoUsed[1] = 999999999.0f;
          for (nx= 0 ; nx < 10; nx++) {
              for (ny = 0; ny < 3; ny++) {
                m_a3x10fIndexHKLErrorRefs[ny][nx] = 0.0;
                m_a3x10fIndexHKLErrorIntensity[ny][nx] = 0.0;
              };
          };


          oStatOutlierMM.vClear();
          oStatOutlierROT.vClear();
          for (nRef = 0; nRef < nNumReflns; nRef++)
          {
             poRefln = m_poReflnlist->poGetRefln(nRef);

             if ((anTwinID[nRef] % 1000)!= m_nTwinID)
             {
                 pnRejectFlag[nRef] = 1;
                 nNumOtherComp++;
                 continue;
             }
             
             nTwinLaw = anTwinID[nRef]/1000;

             // Get the rotation information.

             // Rotation information.
             vCopyVec3D(&aa3fE3[anGonioMatNum[nRef]][0],&a3fE3[0]); // The crystal rotation vector
             vCross3D(a3fE3, a3fS0, a3fE3CrossS0);

             if( anFlopCount[nRef] >= 4)
             {
                 pnRejectFlag[nRef] = 1;
                 nNumRejects++;
                 continue;
             }

              // Get the observed virtual detector coordinates for this reflection
              // Millimeter coordinates were calculated in nCopyReflns

              a3fVDCO[0] = afObsXmm[nRef];
              a3fVDCO[1] = afObsYmm[nRef];
              a3fVDCO[2] = afObsZmm[nRef];

              // Calculate the vector TS = matrix DD * vector VDCO
              // Normalize it.  Since the radius of the Ewald sphere is 1 in this
              // routine, then VNTS is the diffracting vector S with unit length.

              vCopyVec3D(a3fVDCO,aa3fVDCO+3*nRef);
              vMulMat3DVec3D(a3x3fDD, &a3fVDCO[0], &a3fVNTS[0]);            
                          // Add in (and compute) the precessing reciprocal shift 
                                                  
              /////////////////////////////////////////////////////////////////////////////////
              // Add in the reciprocal shift 
              if( afCalcRotMid[nRef] > -999.0 ) 
              {
                 vCopyMat3D(aa3x3fPhiRotMat + 9 * nRef, &a3x3fTempMat1[0][0]);
                 vMulMat3DVec3D(a3x3fTempMat1, aa3fOrigRecipShift[nTwinLaw], a3fRecipVec);
              }
              else
                 vZeroMat(3, 1, &a3fRecipVec[0]);

              vAddVec3DVec3D(a3fVNTS,a3fRecipVec,a3fVNTS);
              (void) fNormVec3D(a3fVNTS);
  
              vAddVec3DVec3D(&a3fVNTS[0], &a3fS0[0], &a3fXRO[0]);           
              /////////////////////////////////////////////////////////////////////////////////
              
              fLorentz = 1.0f / max(fabs(fDot3D(a3fE3CrossS0, a3fXRO)), 1e-10);


              // Get the observed rotation midpoint.

              fPHO  = afObsRotMid[nRef];
              fRotO = afObsRotWidth[nRef]* Gs_dRADIANS_PER_DEGREE;


              // Fill in sigma values.

              switch (m_eWeightType) 
              {
              case eWeightRotSigma:
                  if (afObsRotSigma[nRef]>0.0)
                      pfSigs[nRef] = afObsRotSigma[nRef];
                  else
                      pfSigs[nRef] = 1.0;
                  break;
              case eWeightIntensity:
                  f0 = fabs(afObsIntensity[nRef]);
                  if (f0 > 0.0)
                      pfSigs[nRef] = 1.0/f0;
                  else
                      pfSigs[nRef] = 1.0;
                  break;
              case eWeightIoverSig:
                  f0 = afObsIntensitySig[nRef];
                  if (f0> 0.0)
                      pfSigs[nRef] = 1.0/max(0.1,fabs(afObsIntensity[nRef]/f0));
                  else
                      pfSigs[nRef] = 1.0;
                  break;
              case eWeightLorentz:
                  pfSigs[nRef] = fLorentz;
                  break;
              case eWeightRotWidth:
                  pfSigs[nRef] = fRotO;
                  break;
              case eWeightNone:
                  pfSigs[nRef] = 1.0;
                  break;
              }
                                    
                               
              double afObsRots[20];
              const double fObsRotStep = 5.0;
              double a3fLastObsHKL[3];
              double a3fThisObsHKL[3];
              double fBestObsRotResid;
              double a3fBestObsRotResid[3];
              int nObsRot;
              int nMaxObsRot;
              vZeroMat(3,1,&a3fThisObsHKL[0]);
              
              
              afObsRots[0] = afObsRotMid[nRef] - afObsRotWidth[nRef]/2;
              for (nObsRot=1,fPHO = fObsRotStep; 
              fPHO < afObsRotWidth[nRef]; fPHO+= fObsRotStep) 
                  afObsRots[nObsRot++] = afObsRotMid[nRef] - afObsRotWidth[nRef]/2 + fPHO;
              afObsRots[nObsRot++] = afObsRotMid[nRef] + afObsRotWidth[nRef]/2;
              afObsRots[nObsRot++] = afObsRotMid[nRef];
              
              if( 20 <= nObsRot )
              {
                  pnRejectFlag[nRef] = 1;
                  nNumRejects++;
                  nNumDegRejects++;
                  continue;
              }
              nMaxObsRot = nObsRot;
              
              fBestObsRotResid = 100.0;
              nStat = 1;
              for (nObsRot = 0; nObsRot< nMaxObsRot; nObsRot++)
              {
                  fPHO = afObsRots[nObsRot];
                  
                  // Compute an initial 3x3fPhiGon
                  // Later, we will recompute this to account for the actual calculated
                  // rotation rather than the observed rotation.
                  vConvRotVec3DMat3D(fPHO, a3fE3, & a3x3fTempMat1[0][0]);
                  vMulMat3DMat3D(a3x3fTempMat1, aa3x3fGonMatrix[anGonioMatNum[nRef]], a3x3fPhiGon);
                  
                  // Compute the inverse fPhiGon matrix.
                  
                  vCopyMat3D(& a3x3fPhiGon[0][0], & a3x3fPhiGonInv[0][0]);
                  vTranMat3D(a3x3fPhiGonInv);
                  
                  // Calculate the fractional hkl entries.
                  //              -1        -1
                  //       R = (TUB)  * ( PHR  * (VNTS + S0))
                  
                  vMulMat3DVec3D(a3x3fPhiGonInv,a3fXRO,a3fXROR);
                  vMulMat3DVec3D( aa3x3fCrysOrientInv[nTwinLaw], a3fXROR, a3fHKL);

                  // This is computed and used later for shift vector refinement.
                  vMulMat3DMat3D(aa3x3fCrysOrientInv[nTwinLaw],a3x3fPhiGonInv,a3x3fCrysPhiGonInv);
                                    
                  vCopyVec3D(a3fThisObsHKL,a3fLastObsHKL);
                  vCopyVec3D(a3fHKL,a3fThisObsHKL);
                  if ((nObsRot==0) || (nObsRot+1==nMaxObsRot))
                      continue;
                  
                  // Get the proper rounded indices.
                  
                  if (!nRoundIndices(aa3x3fCrysOrient[nTwinLaw], a3fThisObsHKL,a3fLastObsHKL,fBestObsRotResid,a3fBestObsRotResid, a3nHKL))
                      nStat = 0;
              }

              if( nStat )
              {
                  // Could not find rounded indices
                  pnRejectFlag[nRef] = 1;
                  nNumRejects++;
                  nNumDegRejects++;
                  continue;
              }

              // Save the previous HKL values.

              a3nPrevHKL[0] = poRefln->nGetH();
              a3nPrevHKL[1] = poRefln->nGetK();
              a3nPrevHKL[2] = poRefln->nGetL();

              // Copy them over to the reflection.
              poRefln->vSetH(a3nHKL[0]);
              poRefln->vSetK(a3nHKL[1]);
              poRefln->vSetL(a3nHKL[2]);

              // Move them over to floating point variables.

              a3fHKL[0] = (float)a3nHKL[0];
              a3fHKL[1] = (float)a3nHKL[1];
              a3fHKL[2] = (float)a3nHKL[2];


              // REJECTION TEST:  Sigma/I

              // Do not want to look at I/sig for predicted list, since these are sometimes not provided.
              if (   (0.0 >= afObsIntensitySig[nRef])
                  || (m_fSigma > (afObsIntensity[nRef]
                  / afObsIntensitySig[nRef])) ||
                  (m_fSharpness < (afObsSharpness[nRef])))
              {
                  pnRejectFlag[nRef] = 2;
                  nNumIgnored++;
                  continue;
              }

              // REJECTION TEST:  Lorentz too high
              if (fLorentz > m_fLorentzMax)
              {
                  pnRejectFlag[nRef] = 2;
                  nNumIgnored++;
                  continue;
              }


              // REJECTION TEST: [0,0,0] reflection.

              if ((0 == a3nHKL[0]) && (0 == a3nHKL[1]) && (0 == a3nHKL[2]))
              {
                  pnRejectFlag[nRef] = 1;
                  nNumRejects++;
                  continue;
              }

              // Compute XRC
              // XRC is the reciprocal lattice vector X rotated into diffracting
              // position in diffractometer coordinate system.

              vMulMat3DVec3D(aa3x3fCrysOrient[nTwinLaw], a3fHKL, a3fTempVec1);
              vMulMat3DVec3D(a3x3fPhiGon, a3fTempVec1, a3fXRC);
              vSubVec3DVec3D(a3fXRO, a3fXRC, a3fXROXRC);

              // REJECTION TEST:  Out of resolution bounds.

              f0 = fLenVec3D(a3fXRC);
              if ((f0 < fMinResStar) || (f0 > fMaxResStar) )
              {
                  pnRejectFlag[nRef] = 2;
                  nNumIgnored++;
                  continue;
              }
              m_a2fResoUsed[0] = max(m_a2fResoUsed[0],fOrigSourceWave/f0);
              m_a2fResoUsed[1] = min(m_a2fResoUsed[1],fOrigSourceWave/f0);

             
              // We want to discover the amount of rotation necc. to bring XRC
              // back onto the Ewald sphere.
              // a3fXRO contains the 'ideal' Ewald sphere position (calculated from (x,y) points on plate)
              // a3fXRC contains the 'calculated' Ewald sphere position (probably NOT on the sphere).

              double a3fDStarV1[3];
              double a3fDStarV2[3];
              double a3fDStarV3[3];
              double fk1,fk2,fb;
              double a2fSolutions[2];
              
              // Find the three vectors v1,v2 and v3.
              // Project the offset d* vector of the spot down into the rotation plane to get v1
              vMulVec3DScalar(a3fE3,fDot3D(a3fE3,a3fXRC),a3fTempVec1);
              vSubVec3DVec3D(a3fXRC,a3fTempVec1,a3fDStarV1);
              // Take cross product of m_a3fRotVector and v1 to get v2.
              // This gives v1 cross v2 == m_a3fRotVector
              vCross3D(a3fE3,a3fDStarV1,a3fDStarV2);
              // Calculate v3
              vSubVec3DVec3D(a3fXRC,a3fDStarV1,a3fDStarV3);
              
              // Calculate constants a (int a3fTemp1), k1,k2 and b
              
              vSubVec3DVec3D(a3fDStarV3,a3fS0,a3fTempVec1);
              fb = fDot3D(a3fDStarV1,a3fDStarV1) + fDot3D(a3fTempVec1,a3fTempVec1) - 1.0;
              fk1 = 2.0* fDot3D(a3fDStarV1,a3fTempVec1);
              fk2 = 2.0* fDot3D(a3fDStarV2,a3fTempVec1);
              
              // Solve quadratic system.  Choose solution with the lowest absolute value.
              nx = nSolveQuadratic(fk1*fk1 + fk2*fk2, fk2*fb,fb*fb - fk1*fk1,a2fSolutions);
              if (nx) {
                  f0 = 1000.0;
                  for (ny=0;ny<nx;ny++) {
                      if ((a2fSolutions[ny]>1.0) || (a2fSolutions[ny]<-1.0)) {
                          break;
                      };
                      f1 = asin(a2fSolutions[ny]);
                      if (ABS(f1)<ABS(f0))
                          f0 = f1;
                  };
              };
              if ((nx==0) || (ny !=nx) || (f0==1000.0)) {
                  pnRejectFlag[nRef] = 2;
                  nNumIgnored++;
                  continue;
              };
              fRotErr = f0/fRADIANS_PER_DEGREE;

              // Set the new calculated midpoint here.  This is used in mosaicity refinements.              
              afCalcRotMid[nRef] = afObsRotMid[nRef] + fRotErr;

              // Calculate the rotated X-calc
              vConvRotVec3DMat3D(fRotErr,&a3fE3[0],&a3x3fTempMat1[0][0]);
              vMulMat3DVec3D(a3x3fTempMat1,a3fXRC,a3fXRCR);

              // Calculate the rotated difference.
              vSubVec3DVec3D(a3fXRO,a3fXRCR,a3fXROXRCR);

              // Calculate the X projection vector.
              vMulVec3DScalar(a3fE3,fDot3D(a3fXRCR,a3fE3),a3fTempVec1);
              vSubVec3DVec3D(a3fXRC,a3fTempVec1,a3fTempVec1);
              vCross3D(a3fE3,a3fTempVec1,a3fXProjVector);
              fNormVec3D(a3fXProjVector);             


              // We split the reflections into two categories.
              // In the first category are reflections who's XRC-XRO error residual
              // is due mainly to an imperfectly known rotation centroid.
              // Reflections in the other category have residuals due mainly
              // to imprefectly known MM residuals.
              // We treat the reflections differently based on the category they are in:
              // Category 1)
              // Rotate the reflection so that it is in diffracting condition (i.e. use XRCR-XRO error residual).
              // Also, use a projection vector in the derivative calculations.
              // Category 2)
              // Use XRC-XRO residual.  Don't use a projection vector in the derivative calculations.

              // Calculate ROT error and MM error
              fROTError = fDot3D(a3fXProjVector,a3fXROXRC);
                          fMMError = fLenVec3D(a3fXROXRCR);
              fROTError = ABS(fROTError);


              // Add the residuals to the stat arrays.
              oStatOutlierMM.vAdd((double) fMMError);
              oStatOutlierROT.vAdd((double) fROTError);
              if (((fDevROTResid!=0.0) && (fROTError > fAverageROTResid + 3.0*fDevROTResid)) ||
                  ((fDevMMResid!=0.0) && (fMMError > fAverageMMResid + 3.0*fDevMMResid))) {
                  pnRejectFlag[nRef] = 1;
                  nNumRejects++;
                  nNumMMRejects++;
                  continue;
              };

                                          
              if (((fRotationWidth>=3.0) || (m_bUseWideSlice)) && (!m_bNotUseWideSlice)) {
                  // Category 1
                                  
                                  // Add a little bit of contribution from phi rotation.
                                  vMulVec3DScalar(a3fXProjVector,fAverageMMResid*0.1,a3fTempVec1);
                                  if (fDot3D(a3fXProjVector,a3fXROXRC)<0.0)
                                          vMulVec3DScalar(a3fTempVec1,-1.0,a3fTempVec1);
                                  vAddVec3DVec3D(a3fXROXRCR,a3fTempVec1,a3fTempVec1);


                  vCopyVec3D(a3fTempVec1,aa3fResidVec + 3*nRef);
                  afResidVecLen[nRef] = pfVal[nRef] = fLenVec3D(a3fTempVec1);
                  vConvRotVec3DMat3D(afCalcRotMid[nRef], a3fE3, & a3x3fTempMat1[0][0]);
                  vMulMat3DMat3D(a3x3fTempMat1, aa3x3fGonMatrix[anGonioMatNum[nRef]], a3x3fPhiGon);                  
                  vCopyMat3D(&a3x3fPhiGon[0][0],aa3x3fPhiRotMat+9*nRef);
                                  vZeroMat(3,1,aa3fXProjVec +nRef*3);


              } else {
                  // Cateogry 2
                  vCopyVec3D(a3fXROXRC,aa3fResidVec + 3*nRef);
                  afResidVecLen[nRef] = pfVal[nRef] = fLenVec3D(a3fXROXRC);
                  vConvRotVec3DMat3D(afObsRotMid[nRef], a3fE3, & a3x3fTempMat1[0][0]);
                  vMulMat3DMat3D(a3x3fTempMat1, aa3x3fGonMatrix[anGonioMatNum[nRef]], a3x3fPhiGon);                  
                  vCopyMat3D(&a3x3fPhiGon[0][0],aa3x3fPhiRotMat+9*nRef);
                  vZeroMat(3,1,aa3fXProjVec +nRef*3);
              }
                          

              // Calculate the HKL difference between the observed and calculated down the projected portion.
              // Reflections for which this is too large should be rejected.
              vMulMat3DVec3D(a3x3fCrysPhiGonInv,a3fXROXRCR,a3fTempVec1);
              for (nx=0;nx<3;nx++) {
                  if (fabs(a3fTempVec1[nx])>0.3) {
                      pnRejectFlag[nRef] = 1;
                      nNumRejects++;
                      nNumDegRejects++;
                      break;
                  };
              };
              if (nx!=3)
                  continue;

              // This code transforms X-calc (rotated into difracting condition) into pixel coordinates
              // by first: adding in the S0 vector
              vSubVec3DVec3D(a3fXRCR, a3fS0, a3fTempVec1);
              
              // next:  accounting for the precessing reciprocal shift. (we recalculate a3fRecipVec)
              vCopyMat3D(aa3x3fPhiRotMat+9*nRef,&a3x3fTempMat1[0][0]);

              vMulMat3DVec3D(a3x3fTempMat1,aa3fOrigRecipShift[nTwinLaw],a3fRecipVec);
              
              
              // RB: I was concerned that Thad didn't check here the return values of nScaleSOnPlateItr() and nMMtoPixel().
              // So I am adding those safety checks now, so that the "bad" reflections are ignored.
              bool bScatteredVecIntersectsDetector = ( 0 == m_ppoDetector[0]->nScaleSOnPlateItr(a3fTempVec1,a3fVTS,a3fRecipVec) );

              if( bScatteredVecIntersectsDetector )
              {
                  vMulMat3DVec3D(a3x3fDDInv, a3fVTS, a3fVDCC);
                  
                  float     fPx0 = 0.0f;
                  float     fPx1 = 0.0f;
                  if( 0 == m_ppoDetector[0]->m_poSpatial->nMMtoPixel(a3fVDCC[0],
                                                                     a3fVDCC[1],
                                                                     a3fVDCC[2],
                                                                    &fPx0, 
                                                                    &fPx1) )
                  {
                        poRefln->vSetField(m_poReflnlist->m_nFI_fCalcPx0, fPx0);
                        poRefln->vSetField(m_poReflnlist->m_nFI_fCalcPx1, fPx1);
                  }
                  else
                  {
                      pnRejectFlag[nRef] = 2;
                      nNumIgnored++;
                      continue;
                  }
              }
              else
              {
                  pnRejectFlag[nRef] = 2;
                  nNumIgnored++;
                  continue;
              }
              
              // Calculate the difference in length for the virtual detector coordinates

              vSubVec3DVec3D(a3fVDCO, a3fVDCC, a3fTempVec1);
              a2fVDErr[0] = a3fTempVec1[0];
              a2fVDErr[1] = a3fTempVec1[1];

              // REJECTION TEST:  Are the displacements of x,y,theta within bounds?
              if ( (ABS(a3fVDCO[0] - a3fVDCC[0]) > m_a4fRejectLimits[0]) ||
                   (ABS(a3fVDCO[1] - a3fVDCC[1]) > m_a4fRejectLimits[1]) )
              {
                  pnRejectFlag[nRef] = 1;
                  nNumRejects++;
                  nNumMMRejects++;
                  continue;
              }
              
              if (ABS(fRotErr) > max( max(fRotationWidth,m_a4fRejectLimits[2]),afObsRotWidth[nRef]))
              {
                  pnRejectFlag[nRef] = 1;
                  nNumRejects++;
                  nNumDegRejects++;
                  continue;
              }

              // Check that the residual is smaller than all other residuals.
              if (anCompete[nRef]>-1) {
                  int nRef2;
                  anLostCompetition[nRef] = 0;
                  for (nRef2 = anCompete[nRef]; nRef2 != nRef; nRef2 = anCompete[nRef2]) {
                      if ((pnRejectFlag[nRef2] == 0) && (afResidVecLen[nRef2] > 0.0) && (afResidVecLen[nRef]>afResidVecLen[nRef2])) {
                          // There is a better reflection.
                          anLostCompetition[nRef] = 1;
                          pnRejectFlag[nRef] = 2;
                          nNumIgnored++;                          
                          break;
                      };
                  };
                  if (pnRejectFlag[nRef])
                      continue;
              } else
                  anLostCompetition[nRef] = 0;
              
              for (nx=0;nx<3;nx++) {
                m_a3x10fIndexHKLErrorRefs[nx][min(10-1,(int) floor(0.5 + ABS(a3fBestObsRotResid[nx])/0.05))]++;
                m_a3x10fIndexHKLErrorIntensity[nx][min(10-1,(int) floor(0.5 + ABS(a3fBestObsRotResid[nx])/0.05))]+=afObsIntensity[nRef];
              };

              m_fRMSDEG += fRotErr * fRotErr;
              m_fRMSMM  += a2fVDErr[0] * a2fVDErr[0]  +  a2fVDErr[1] * a2fVDErr[1];
              if (pfVal[nRef]!=0.0)
                m_fRMSANG += (pfVal[nRef]*pfVal[nRef]);

              fTotalError += pfVal[nRef];

              // See if any of the HKL values have changed from their settings.
              if ((a3nPrevHKL[0]!=a3nHKL[0]) || (a3nPrevHKL[1]!=a3nHKL[1]) || (a3nPrevHKL[2]!=a3nHKL[2])) {
                  bHKLChanged = TRUE;
                  anFlopCount[nRef]++;
              };

              // REJECTION TEST:  The residual is too small.
              if (afResidVecLen[nRef]<=0.0) {
                  pnRejectFlag[nRef] = 2;
                  // DO NOT COUNT THESE AS IGNORED OR REJECTED.
                  // We just don't want to refine on them because they are "too good"
                  // Really though, they cause a divide by zero when calculating derivatives.
                  continue;
              };

      } // end of reflection loop.

          oStatOutlierMM.nMarkStat(false,oStatOutlierMM.S_ABOVEBELOW_ITER,3.0);
          oStatOutlierROT.nMarkStat(false,oStatOutlierROT.S_ABOVEBELOW_ITER,3.0);
      fAverageMMResid = oStatOutlierMM.fAverage();
      fDevMMResid = oStatOutlierMM.fStandardDeviation();
      fAverageROTResid = oStatOutlierROT.fAverage();
      fDevROTResid = oStatOutlierROT.fStandardDeviation();

      // Compute RMS values.
      if (0 < (nNumReflns - nNumRejects -nNumIgnored -nNumOtherComp) )
      {
          m_fRMSDEG = sqrt((double)((m_fRMSDEG) / (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp)));
          m_fRMSMM  = sqrt((double)((m_fRMSMM)  / (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp)));
          m_fRMSANG = sqrt((double)((m_fRMSANG) / (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp)));                

          m_a7nReflnNumbers[eReflectionNumberTotal] = nNumReflns;
          m_a7nReflnNumbers[eReflectionNumberAccept] = (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp);
          m_a7nReflnNumbers[eReflectionNumberRejects] = nNumRejects;
          m_a7nReflnNumbers[eReflectionNumberIgnored] = nNumIgnored;
          m_a7nReflnNumbers[eReflectionNumberDegRejects] = nNumDegRejects;
          m_a7nReflnNumbers[eReflectionNumberMMRejects] = nNumMMRejects;
          m_a7nReflnNumbers[eReflectionNumberOtherComp] = nNumOtherComp;
          
          fTotalError *= 1000.0 / (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp);
      }
      else
      {
          m_fRMSDEG   = 0.0;
          m_fRMSMM    = 0.0;
          m_fRMSANG   = 0.0;
          fTotalError = 0.0;
      }


      if (fLastTotal == -1) {
          // Change nothing.
      } else if ((nPass==0) && (fTotalError<fLastTotal) && (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp > 0.9*nLastAccept)) {
          // We were able to decrease the residual error.
          // Set the new "original" values to the current state of the objects.
          /*1*/
          a2fOrigS0Rot[0] = m_poSource->fGetRotation(0);
          a2fOrigS0Rot[1] = m_poSource->fGetRotation(1);
          fOrigSourceWave = m_poSource->fGetWavelength();
          /*2*/
          m_ppoCrystals[m_nTwinID-1]->vGetRecipCell(a6fOrigCell);
          /*3*/
          m_ppoCrystals[m_nTwinID-1]->nCalcRotMatrix();
          m_ppoCrystals[m_nTwinID-1]->vGetRotMatrix(&a3x3fTempMat1[0][0]);
              a3fOrigCrysRot[0] = m_ppoCrystals[m_nTwinID-1]->fGetOrientAngle(0);
              a3fOrigCrysRot[1] = m_ppoCrystals[m_nTwinID-1]->fGetOrientAngle(1);
              a3fOrigCrysRot[2] = m_ppoCrystals[m_nTwinID-1]->fGetOrientAngle(2);

          vDeconvMat3D3XYZ(a3x3fTempMat1, &a3fOrigCrysRot[0], &a3fOrigCrysRot[1], &a3fOrigCrysRot[2]);
                  for (nTwinLaw = 0; nTwinLaw <= m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws;nTwinLaw++) {
                          m_ppoCrystals[m_nTwinID-1]->nGetTwinLaw(nTwinLaw,&aa3fOrigCrysTwinLaws[nTwinLaw][0],false);
                          m_ppoCrystals[m_nTwinID-1]->nGetRecipShift(nTwinLaw,&aa3fOrigRecipShift[nTwinLaw][0]);
                          m_ppoCrystals[m_nTwinID-1]->nGetTwinFraction(nTwinLaw,&afOrigTwinFraction[nTwinLaw]);
                  };  
          /*4*/
          m_ppoDetector[0]->nGetModifiedDetectorParams(a6fTemp);
          vCopyVec3D(a6fTemp,a3fOrigRot);
          vCopyVec3D(&a6fTemp[3],a3fOrigTrans);
          /*5*/

          m_fLambda /= 10;

          // We can break out, since the calculated values are relevant.
          if (!bHKLChanged)
              break;          
            
      } else if (nPass==0) {

          // We were not able to decrease the residual error.
          // Put the "original" values back into the objects.

          /*1*/
          m_poSource->vSetRotation(a2fOrigS0Rot[0],a2fOrigS0Rot[1]);
          m_poSource->m_poWavelength->nSetWavelength(fOrigSourceWave);
          /*2*/
          m_ppoCrystals[m_nTwinID-1]->vSetRecipCell(a6fOrigCell);
          /*3*/
          m_ppoCrystals[m_nTwinID-1]->vSetOrientAngles(a3fOrigCrysRot[0],a3fOrigCrysRot[1],a3fOrigCrysRot[2]);

                  for (nTwinLaw = 0; nTwinLaw <= m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws;nTwinLaw++) {
                          m_ppoCrystals[m_nTwinID-1]->vSetTwinLaw(nTwinLaw,&aa3fOrigCrysTwinLaws[nTwinLaw][0],false);
                          m_ppoCrystals[m_nTwinID-1]->vSetRecipShift(nTwinLaw,&aa3fOrigRecipShift[nTwinLaw][0]);
                          m_ppoCrystals[m_nTwinID-1]->vSetTwinFraction(nTwinLaw,afOrigTwinFraction[nTwinLaw]);
                  };  

          /*4*/
          vCopyVec3D(a3fOrigTrans,&a6fTemp[3]);
          vCopyVec3D(a3fOrigRot,&a6fTemp[0]);
          m_ppoDetector[0]->nSetModifiedDetectorParams(a6fTemp);
          /*5*/

          m_fLambda *= 10;
          // Must go around a second time.
      };

    }; // End of double check loop.

    if (nLoopCt > nLoop)
        goto breakout_place;

    fLastTotal = fTotalError;
    nLastAccept = (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp);


    // Print current state if user is debugging.

    if (2 < m_nVerbose)
    {
        printf("\nReflections in list:  %d", nNumReflns);
        printf("\nRefinement residuals ");
        printf("\nrmsResid (1000 * A0^-1) = %f\n",fLastTotal);
        printf("\nrmsResid (Angstrom-1)   = %f\n",m_fRMSANG);
        printf("\nrmsResid (mm)           = %f",m_fRMSMM);
        printf("\nrmsResid (Deg)          = %f\n",m_fRMSDEG);

        printf("\nReflections accepted      :  %d", nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp);
        printf("\nReflections rejected (mm) :  %d", nNumRejects-nNumDegRejects);
        printf("\nReflections ignored       :  %d", nNumIgnored);
        printf(" (Outside I/sigI or Resolution limits)\n");
        printf("\nReflections rejected (phi):  %d", nNumDegRejects);
        if (nNumOtherComp)
            printf("\nReflections twinned   :  %d", nNumOtherComp);
        printf(" (Part of another twin component)\n");
    }

    
    {
        // Calculate derivatives...
        
        // Detector Derivatives.
        
        for (nVar = DET_TRANS1; nVar <= DET_ROT3; nVar++) {
            double a3fVD[3];        // DD*Vpos
            double a3fVDP[3];       // (DD*Vpos)' = DD'*Vpos
            double a3fVDNP[3];      // S' = Normalized(DD*Vpos)'
            
            if (m_anRefineFlags[nVar]) {
                for (nRef = 0; nRef < nNumReflns;nRef++) {
                    if (!pnRejectFlag[nRef])  {
                                                nTwinLaw = anTwinID[nRef]/1000;
                        // Compute the derivative of S = Normalized(DD*VPos) = (DD*Vpos)/Length(DD*Vpos)
                        // = (DD*Vpos)'/Length(DD*Vpos + RecipShift) - (DD*Vpos + RecipShift) * (dot((DD*Vpos + RecipShift),(DD*Vpos)')/(Length(DD*Vpos + RecipShift)^3))
                        vMulMat3DVec3D(a3x3fDD,aa3fVDCO + nRef*3,a3fVD);
                                                // Add in the reciprocal shift vector.
                                                vCopyVecND(9,(double*) (aa3x3fPhiRotMat + nRef*9),&a3x3fTempMat1[0][0]);
                        vMulMat3DVec3D(a3x3fTempMat1,aa3fOrigRecipShift[nTwinLaw],a3fTempVec1);
                                                vAddVec3DVec3D(a3fTempVec1,a3fVD,a3fVD);

                        if (nVar<=DET_TRANS3)
                            vMulMat3DVec3D(a6x3x3fDDDeriv[nVar - DET_TRANS1+3],aa3fVDCO + nRef*3,a3fVDP);
                        else
                            vMulMat3DVec3D(a6x3x3fDDDeriv[nVar - DET_ROT1],aa3fVDCO + nRef*3,a3fVDP);
                        f0 = fLenVec3D(a3fVD);
                        vMulVec3DScalar(a3fVDP,1.0/f0,a3fTempVec1);
                        vMulVec3DScalar(a3fVD,-fDot3D(a3fVD,a3fVDP)/(f0*f0*f0),a3fTempVec2);
                        vAddVec3DVec3D(a3fTempVec1,a3fTempVec2,a3fVDNP);                  
                        
                        vCopyVec3D(a3fVDNP,a3fDeriv);
                        ppfDeriv[nVar][nRef] = fDot3D(a3fDeriv,aa3fResidVec + nRef*3)/afResidVecLen[nRef];
                    };
                };
            };
        };

                // Crystal shift vector refinement (These actually  modify the 'detector' half of the refinement equation)
                for (nVar = CELL_RECIP00; nVar <= CELL_RECIP_MAX;nVar++) {
            double a3fVD[3];        // DD*Vpos
            double a3fVDP[3];       // (DD*Vpos)' = DD'*Vpos
            double a3fVDNP[3];      // S' = Normalized(DD*Vpos)'

                        if (m_anRefineFlags[nVar]) {
                                for (nRef = 0; nRef < nNumReflns;nRef++) {
                    if (!pnRejectFlag[nRef])  {
                                                nTwinLaw = anTwinID[nRef]/1000;
                                                if (nTwinLaw == (nVar - CELL_RECIP00)/3) {
                                                        vCopyVecND(9,(double*) (aa3x3fPhiRotMat + nRef*9),&a3x3fTempMat1[0][0]);
                                                        a3fTempVec3[0] = 0.0;
                                                        a3fTempVec3[1] = 0.0;
                                                        a3fTempVec3[2] = 0.0;
                                                        a3fTempVec3[(nVar - CELL_RECIP00) % 3] = 1.0;
                                                        vMulMat3DVec3D(a3x3fTempMat1,a3fTempVec3,a3fVDP);
                                                        

                                                        // Compute the derivative of S = Normalized(DD*VPos) = (DD*Vpos)/Length(DD*Vpos)
                                                        // = (Recip)'/Length(DD*Vpos + RecipShift) - (DD*Vpos + RecipShift) * (dot((DD*Vpos + RecipShift),(Recip)')/(Length(DD*Vpos + RecipShift)^3))
                                                        vMulMat3DVec3D(a3x3fDD,aa3fVDCO + nRef*3,a3fVD);
                                                        // Add in the reciprocal shift vector.
                                                        
                                                        vMulMat3DVec3D(a3x3fTempMat1,aa3fOrigRecipShift[nTwinLaw],a3fTempVec1);
                                                        vAddVec3DVec3D(a3fTempVec1,a3fVD,a3fVD);
                                                        
                                                        f0 = fLenVec3D(a3fVD);
                                                        vMulVec3DScalar(a3fVDP,1.0/f0,a3fTempVec1);
                                                        vMulVec3DScalar(a3fVD,-fDot3D(a3fVD,a3fVDP)/(f0*f0*f0),a3fTempVec2);
                                                        vAddVec3DVec3D(a3fTempVec1,a3fTempVec2,a3fVDNP);                  
                                                        
                                                        vCopyVec3D(a3fVDNP,a3fDeriv);
                                                        ppfDeriv[nVar][nRef] = fDot3D(a3fDeriv,aa3fResidVec + nRef*3)/afResidVecLen[nRef];
                                                } else
                                                        ppfDeriv[nVar][nRef] = 0.0;
                                        };
                                };
                        };
                };

        
        // Source Parameters.
        for (nVar = SOURCE_ROT1; nVar<= SOURCE_ROT2; nVar++) {
            if (m_anRefineFlags[nVar]) {
                double a3fS0P[3];
                // This one is easy.  We compute the derivative outside the reflection loop,
                // and just dot it in to each residual.
                m_poSource->vCalcGetUnrotatedS0(a3fTempVec1);
                vMulMat3DVec3D(a2x3x3fS0Deriv[nVar - SOURCE_ROT1],a3fTempVec1,a3fS0P);
                for (nRef = 0; nRef < nNumReflns;nRef++) {
                    if (!pnRejectFlag[nRef])  {
                        vCopyVec3D(a3fS0P,a3fDeriv);
                        ppfDeriv[nVar][nRef] = fDot3D(a3fDeriv,aa3fResidVec + nRef*3)/afResidVecLen[nRef];
                    };
                };
            };
        };
        
        
        // Crystal rotations.
        for (nVar = CRYS_ROT1; nVar<= CRYS_ROT3; nVar++) {
            if (m_anRefineFlags[nVar])  {
                for (nRef = 0; nRef < nNumReflns;nRef++) {
                    if (!pnRejectFlag[nRef])   {
                        int* pnHKL = (*m_poReflnlist)[nRef].pnGetHKL();
                                                nTwinLaw = anTwinID[nRef]/1000;
                        a3fTempVec1[0] = pnHKL[0];
                        a3fTempVec1[1] = pnHKL[1];
                        a3fTempVec1[2] = pnHKL[2];
                        vMulMat3DVec3D(a3x3fB,a3fTempVec1,a3fTempVec2);         
                        vMulMat3DVec3D(a3x3x3fCrysRotDeriv[nVar - CRYS_ROT1],a3fTempVec2,a3fTempVec4);
                                                vMulMat3DVec3D(aa3x3fCrysTwinLaws[nTwinLaw],a3fTempVec4,a3fTempVec3);
                        // Project this vector (remove the portion parrallel to the projection vector).
                                                vMulVec3DScalar(aa3fXProjVec + nRef*3,fDot3D(a3fTempVec3,aa3fXProjVec + nRef*3),a3fTempVec2);
                        vSubVec3DVec3D(a3fTempVec3,a3fTempVec2,a3fTempVec3);                  

                        // Complete the rest of derivative.
                        vCopyVecND(9,(double*) (aa3x3fPhiRotMat + nRef*9),&a3x3fTempMat1[0][0]);
                        vMulMat3DVec3D(a3x3fTempMat1,a3fTempVec3,a3fDeriv);
                        
                        // Remember:  This is on the *negative* side of the residual difference.
                        ppfDeriv[nVar][nRef] = - fDot3D(a3fDeriv,aa3fResidVec + nRef*3)/afResidVecLen[nRef];
                        
                        // Unfortunatly, with the advent of projective transformations, it is undesirable to
                        // refine all three crys-rots at the same time, since a singularity MUST occur.
                        // To solve this problem, we disable one of the CRYS_ROT? variables in each cycle (determined by nLoopCt).
                        // This is done by setting it's derivatives to be something stupid.
                        //if ((nLoopCt % 3)== nVar - CRYS_ROT1)
                        //    ppfDeriv[nVar][nRef] = 1.0;
                        
                    };
                };
            };
        };

                // Twin Law refinement.
                for (nVar = CELL_TWINLAW10; nVar <= CELL_TWINLAW_MAX; nVar++) {
                        if (m_anRefineFlags[nVar]) {
                                for (nRef = 0; nRef < nNumReflns; nRef++) {
                                        if (!pnRejectFlag[nRef]) {
                        int* pnHKL = (*m_poReflnlist)[nRef].pnGetHKL();
                                                nTwinLaw = anTwinID[nRef]/1000;
                                                if (nTwinLaw == 1 + (nVar - CELL_TWINLAW10)/3) {
                                                        a3fTempVec1[0] = pnHKL[0];
                                                        a3fTempVec1[1] = pnHKL[1];
                                                        a3fTempVec1[2] = pnHKL[2];
                                                        vMulMat3DVec3D(a3x3fB,a3fTempVec1,a3fTempVec2);         
                                                        vMulMat3DVec3D(a3x3fCrysRot,a3fTempVec2,a3fTempVec4);
                                                        vMulMat3DVec3D(a3x3x3fCrysTwinLawDeriv[nTwinLaw][(nVar - CELL_TWINLAW10) % 3],a3fTempVec4,a3fTempVec3);
                                                        // Project this vector (remove the portion parrallel to the projection vector).
                                                        vMulVec3DScalar(aa3fXProjVec + nRef*3,fDot3D(a3fTempVec3,aa3fXProjVec + nRef*3),a3fTempVec2);
                                                        vSubVec3DVec3D(a3fTempVec3,a3fTempVec2,a3fTempVec3);                  
                                                        
                                                        // Complete the rest of derivative.
                                                        vCopyVecND(9,(double*) (aa3x3fPhiRotMat + nRef*9),&a3x3fTempMat1[0][0]);
                                                        vMulMat3DVec3D(a3x3fTempMat1,a3fTempVec3,a3fDeriv);
                                                        
                                                        // Remember:  This is on the *negative* side of the residual difference.
                                                        ppfDeriv[nVar][nRef] = - fDot3D(a3fDeriv,aa3fResidVec + nRef*3)/afResidVecLen[nRef];
                                                } else
                                                        ppfDeriv[nVar][nRef] = 0.0;
                                        };
                                };
                        };
                };

                        
        // Cell parameters.
        for (nVar = CELL_A; nVar<= CELL_GAM; nVar++) {
            if (m_anRefineFlags[nVar]) {
                for (nRef = 0; nRef < nNumReflns;nRef++) {
                    if (!pnRejectFlag[nRef])  {
                        int* pnHKL = (*m_poReflnlist)[nRef].pnGetHKL();
                                                nTwinLaw = anTwinID[nRef]/1000;
                        a3fTempVec1[0] = pnHKL[0];
                        a3fTempVec1[1] = pnHKL[1];
                        a3fTempVec1[2] = pnHKL[2];
                        vMulMat3DVec3D(a7x3x3fCrysDeriv[nVar - CELL_A +1],a3fTempVec1,a3fTempVec2);
                        vMulMat3DVec3D(a3x3fCrysRot,a3fTempVec2,a3fTempVec4);
                                                vMulMat3DVec3D(aa3x3fCrysTwinLaws[nTwinLaw],a3fTempVec4,a3fTempVec3);
                        
                        // Project this vector (remove the portion parrallel to the projection vector).
                                                vMulVec3DScalar(aa3fXProjVec + nRef*3,fDot3D(a3fTempVec3,aa3fXProjVec + nRef*3),a3fTempVec2);
                        vSubVec3DVec3D(a3fTempVec3,a3fTempVec2,a3fTempVec3);                  
                       
                        // Compute the rest of the derivative.
                        vCopyVecND(9,(double*) (aa3x3fPhiRotMat + nRef*9),&a3x3fTempMat1[0][0]);
                        vMulMat3DVec3D(a3x3fTempMat1,a3fTempVec3,a3fDeriv);
                        ppfDeriv[nVar][nRef] = - fDot3D(a3fDeriv,aa3fResidVec + nRef*3)/afResidVecLen[nRef];
                    };
                };
            };
        };
        
        // Wavelength refinement.
        if (m_anRefineFlags[SOURCE_WAVE]) {
            nVar = SOURCE_WAVE;
            // The only matrix that is affected by the wavelength is the B matrix,
            // so we simply divide the B matrix by the wavelength to get the derivative w.r.t. wavelength.
            
            vMulVecNDScalar(9,&a3x3fB[0][0],1.0/fOrigSourceWave,&a3x3fTempMat1[0][0]);
            for (nRef = 0; nRef < nNumReflns;nRef++) {
                if (!pnRejectFlag[nRef])  {
                                        nTwinLaw = anTwinID[nRef]/1000;
                    int* pnHKL = (*m_poReflnlist)[nRef].pnGetHKL();
                    a3fTempVec1[0] = pnHKL[0];
                    a3fTempVec1[1] = pnHKL[1];
                    a3fTempVec1[2] = pnHKL[2];
                    vMulMat3DVec3D(a3x3fTempMat1,a3fTempVec1,a3fTempVec2);      
                    vMulMat3DVec3D(a3x3fCrysRot,a3fTempVec2,a3fTempVec3);
                                        vMulMat3DVec3D(aa3x3fCrysTwinLaws[nTwinLaw],a3fTempVec3,a3fTempVec4);
                    vCopyVecND(9,(double*) (aa3x3fPhiRotMat + nRef*9),&a3x3fTempMat2[0][0]);
                    vMulMat3DVec3D(a3x3fTempMat2,a3fTempVec4,a3fDeriv);
                    // Remember:  This is on the *negative* side of the residual difference.
                    ppfDeriv[nVar][nRef] = - fDot3D(a3fDeriv,aa3fResidVec + nRef*3)/afResidVecLen[nRef];
                };
            };
        };
        
        
        // Solve for shift values.
        
        if (0 >= nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp)
        {
            printf("\nERROR: NO reflections accepted for refinement!"
                "\n       Check I/sigI, resolution and rejection limits.\n");
            goto breakout_place;
        }
        
        // We must have enough reflections here.
        if (m_a7nReflnNumbers[eReflectionNumberAccept]<nNumParams*1.5) {
            vZeroMat(1,nNumParams,&afDelta[0]);
            vZeroMat(1,nNumParams,&afSigma[0]);
        } else {
            nOpLevenbergMarquardt((nLoopCt==nLoop), nNumParams, nNumReflns, ppfDeriv,
                pfVal, pfSigs, pnRejectFlag,
                afDelta, afSigma);
        };
    };



    // Next, we change the data values.  When we run the next loop, we will
    // try these to see if they give better results.  If they do, they are kept.
    // If not, we always have the "Orig" variables to fall back on.

    // DET_TRANS1
    // DET_TRANS2
    // DET_TRANS3
    // DET_ROT1
    // DET_ROT2
    // DET_ROT3

    for (nx = DET_TRANS1; nx <= DET_TRANS3; nx++)
        a6fTemp[nx - DET_TRANS1 + 3] = a3fOrigTrans[nx - DET_TRANS1] + afDelta[nx];
    for (nx = DET_ROT1; nx <= DET_ROT3; nx++)
        a6fTemp[nx - DET_ROT1] = a3fOrigRot[nx - DET_ROT1] + afDelta[nx];
    m_ppoDetector[0]->nSetModifiedDetectorParams(a6fTemp);

    // Sometimes we have problems with these when using modified detector parameters.
    // A tiny residual will show up when transforming from the modified parameters.
    if (m_anRefineFlags[DET_TRANS1]==0) 
        m_ppoDetector[0]->nSetStandardDetectorParams(&a3fOrigTransStandard[0]-3,3 + 0,1);
    if (m_anRefineFlags[DET_TRANS2]==0) 
        m_ppoDetector[0]->nSetStandardDetectorParams(&a3fOrigTransStandard[0]-3,3 + 1,1);
    if (m_anRefineFlags[DET_TRANS3]==0) 
        m_ppoDetector[0]->nSetStandardDetectorParams(&a3fOrigTransStandard[0]-3,3 + 2,1);

        // CELL_RECIP00
        // CELL_RECIP01
        // CELL_RECIP02
        // CELL_RECIP10
        // CELL_RECIP11

        
        for (nx = CELL_RECIP00; nx < min(CELL_RECIP00 + 3*(1+m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws),CELL_RECIP12+1);nx++) {
        
                a3fTemp[(nx - CELL_RECIP00) % 3] = afDelta[nx] + aa3fOrigRecipShift[(nx - CELL_RECIP00)/3][(nx - CELL_RECIP00) % 3];
                if (2 == ((nx - CELL_RECIP00) % 3)) {
                        m_ppoCrystals[m_nTwinID-1]->vSetRecipShift((nx - CELL_RECIP00)/3,&a3fTemp[0]);
                };
        };

    // SOURCE_ROT1
    // SOURCE_ROT2

    a3fTempVec1[0] = a2fOrigS0Rot[0]+afDelta[SOURCE_ROT1];
    a3fTempVec1[1] = a2fOrigS0Rot[1]+afDelta[SOURCE_ROT2];
    m_poSource->vSetRotation(a3fTempVec1[0], a3fTempVec1[1]);

    // SOURCE_WAVE

    m_poSource->m_poWavelength->nSetWavelength(fOrigSourceWave + afDelta[SOURCE_WAVE]);

    // CRYS_ROT1
    // CRYS_ROT2
    // CRYS_ROT3

    for (nx = CRYS_ROT1; nx <= CRYS_ROT3; nx++) {
        a3fTemp[nx - CRYS_ROT1] = a3fOrigCrysRot[nx - CRYS_ROT1] + afDelta[nx];
    };
    m_ppoCrystals[m_nTwinID-1]->vSetOrientAngles(a3fTemp[0],a3fTemp[1],a3fTemp[2]);
    
    // CELL_A
    // CELL_A+1
    // CELL_A+2
    // CELL_A+3
    // CELL_A+4
    // CELL_A+5
    
    vCopyVecND(6, a6fOrigCell, a6fTemp);
    for (nx = CELL_A; nx < CELL_A + nNumCrystalParams; nx++)
    {
        for (ny = 0; ny < 6; ny++)
        {
            a6fTemp[ny] += m_a6x6fCrystalVecs[nx-CELL_A][ny] * afDelta[nx];
        }
    }
    m_ppoCrystals[m_nTwinID-1]->vSetRecipCell(a6fTemp);
    m_ppoCrystals[m_nTwinID-1]->vGetCell(a6fTemp);
    if (m_anRefineFlags[ms_nCrysConstraints])
        nForceCell(a6fTemp);
    m_ppoCrystals[m_nTwinID-1]->vSetCell(a6fTemp);   
    
    // CELL_TWINLAW10
    // CELL_TWINLAW11
    // CELL_TWINLAW12
    for (nTwinLaw = 1; nTwinLaw <= m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws;nTwinLaw++) {
        m_ppoCrystals[m_nTwinID-1]->nGetTwinLaw(nTwinLaw,a3fTempVec1,false);
        for (nVar = CELL_TWINLAW10 + (nTwinLaw - 1)*3,nx=0; nVar < CELL_TWINLAW10 + nTwinLaw*3; nVar++,nx++) 
            a3fTempVec1[nx] += afDelta[nVar];
        m_ppoCrystals[m_nTwinID-1]->vSetTwinLaw(nTwinLaw,a3fTempVec1,false);
    };

    // Check all of the delta values.  
    // If they are all below the minimum shift, 
    // then we can halt (after the next iteration)
    // since we are within convergence.
    // This only applies if the -converge option is turned on.

      if (m_bCheckConvergence)
      {
          bAllParametersConverged = TRUE;
          for (nx=0;nx<nNumParams;nx++) 
          {
              if (fabs(afDelta[nx])>ms_pfRefineConvergence[nx]) 
              {
                  bAllParametersConverged = FALSE;
                  if (3 < m_nVerbose)
                  {
                      printf("Loop: %d, not yet converged, a delta > %f\n",
                          nLoopCt, ms_pfRefineConvergence[nx]);
                  }
              }
          };
          if (bAllParametersConverged) 
          {
              // We halt the refinement by setting nLoopCt equal to nLoop-1
              
              if (nLoopCt != nLoop)
              {
                  printf("INFO: converged after %d cycles.\n", nLoopCt);
                  nLoop   = nLoopCt + 1;
              }
          }
          else if (nLoopCt == nLoopIn)
          {
               printf("INFO: no convergence after %d cycles.\n", nLoopCt);
          }
      } else
          bAllParametersConverged = (nLoopCt+1 == nLoop);



  } while (1);    // BIG loop.

breakout_place:
  
  // Compute the mosaicity if we are refining that.
  {
      float fRotOffset;
      int   nRotSet;
      int   nMaxRotOffset;
      double a3x3fOrientAdjustMat[3][3];
      double fMaxMosaicity;
      double fMosaicity;
    
      vZeroMat(3,1,a3fOrigMosaicity);
      vZeroMat(3,1,m_a3fCrysMosaicitySig);

      nStat = 0;
      fMaxMosaicity = 0.0;
      nMaxRotOffset = -1;
      // Need to loop over each scan that we are refining off.
      for (nRotSet = 0; nRotSet < nNumRotations; nRotSet++)
        {
          fRotOffset = 0.0;  // Make sure this offset starts off as 0.
          nStat += nGetMosaicity(0,NULL,NULL,0.0,0.0,nRotSet,afCalcRotMid,&fRotOffset,pnRejectFlag,&anGonioMatNum[0]);

          fMosaicity = m_ppoCrystals[m_nTwinID-1]->fGetMosaicity();
          fMaxMosaicity = max(fMaxMosaicity,fMosaicity);
      } 

      if ((!nStat) && (nNumRotations==1)) {
          // Apply the Phi-centroid correction.  This is effected by modifying the orientation angles.
          // Use the nDeconvMat3DMat3D() function to get small changes in the angles.          
        
//+JWP2010
// Probably should only do this if the offset != 0.
//      cout << "fRotOffset:: "  << fRotOffset
//             << endl;
//-JWP2010
          vConvRotVec3DMat3D(-fRotOffset,&a3fE3[0],&a3x3fOrientAdjustMat[0][0]); 
          vMulMat3DMat3D(a3x3fOrientAdjustMat,aa3x3fGonMatrix[0],a3x3fTempMat1);
          vMulMat3DMat3D(aa3x3fGonMatrixInv[0],a3x3fTempMat1,a3x3fOrientAdjustMat);

          vCopyVec3D(a3fOrigCrysRot,a3fTempVec1);

              nDeconvMat3DMat3D(a3fTempVec1,m_ppoCrystals[m_nTwinID-1]->m_fOrientVectors,a3x3fOrientAdjustMat,true);
              vCopyVec3D(a3fTempVec1,a3fOrigCrysRot);

          m_ppoCrystals[m_nTwinID-1]->vSetOrientAngles(a3fOrigCrysRot[0],a3fOrigCrysRot[1],a3fOrigCrysRot[2]);
      };

      
      if ((!nStat) && (bRefineMosaicity)) {
          a3fOrigMosaicity[0] = fMaxMosaicity;
          a3fOrigMosaicity[1] = 0.0;
          a3fOrigMosaicity[2] = 0.0;
          m_a3fCrysMosaicitySig[0] = 0.1f;
          m_a3fCrysMosaicitySig[1] = 0.0;
          m_a3fCrysMosaicitySig[2] = 0.0;
      } else {
          vCopyVec3D(m_a3fCrysMosaicityShift,a3fOrigMosaicity);
          vZeroMat(3,1,m_a3fCrysMosaicitySig);
      };
      // Might need to override the mosaicity set.
      vSetMosaicity(a3fOrigMosaicity[0]);
  };

  // Compute the photons, and number of accepted reflections in this twin.
  for (nRef=0;nRef<nNumReflns;nRef++) {
      f0 = afObsIntensity[nRef];
      m_a7fReflnNumbers[eReflectionNumberTotal] += f0;
      if (pnRejectFlag[nRef] == 0) 
          m_a7fReflnNumbers[eReflectionNumberAccept] += f0;
      else if (pnRejectFlag[nRef] == 1) 
          m_a7fReflnNumbers[eReflectionNumberRejects] += f0;
      else if (pnRejectFlag[nRef] == 2)
          m_a7fReflnNumbers[eReflectionNumberIgnored] += f0;
  };



  // Reset the lambda parameter.
  m_fLambda = fInitLambda;

  // Calculate the shifts induced.   This is used when printing.

  /*1*/
  vCopyVecND(2, a2fOrigS0Rot, m_a2fSourceRot);
  vSubVecNDVecND(2,m_a2fSourceRot, m_a2fSourceRotShifts, m_a2fSourceRotShifts);
  m_fSourceWave = fOrigSourceWave;
  m_fSourceWaveShift = m_fSourceWave - m_fSourceWaveShift;
  /*2*/
  m_ppoCrystals[m_nTwinID-1]->vGetCell(m_a6fCrysUnitCell);
  vSubVecNDVecND(6, m_a6fCrysUnitCell, m_a6fCrysTurds, m_a6fCrysTurds);
  /*3*/
  vCopyVec3D(a3fOrigCrysRot, m_a3fCrysRot);
  vSubVec3DVec3D(m_a3fCrysRot, m_a3fCrysRotShifts, m_a3fCrysRotShifts);
  /*4*/
  m_ppoDetector[0]->nGetStandardDetectorParams(a6fTemp);    // Make sure to reload the standard parameters, since that is what the user wants to see.
  vCopyVec3D(&a6fTemp[0],m_a10x3fDetRot[0]);
  vSubVec3DVec3D(m_a10x3fDetRot[0], m_a10x3fDetRotShifts[0], m_a10x3fDetRotShifts[0]);
  vCopyVec3D(&a6fTemp[3],m_a10x3fDetTrans[0]);
  vSubVec3DVec3D(m_a10x3fDetTrans[0], m_a10x3fDetTransShifts[0], m_a10x3fDetTransShifts[0]);
  for (nTwinLaw = 0; nTwinLaw <= m_ppoCrystals[m_nTwinID-1]->m_nTwinLaws;nTwinLaw++) {
          m_ppoCrystals[m_nTwinID-1]->nGetTwinLaw(nTwinLaw,&aa3fOrigCrysTwinLaws[nTwinLaw][0],false);
          m_ppoCrystals[m_nTwinID-1]->nGetRecipShift(nTwinLaw,&aa3fOrigRecipShift[nTwinLaw][0]);
          m_ppoCrystals[m_nTwinID-1]->nGetTwinFraction(nTwinLaw,&afOrigTwinFraction[nTwinLaw]);
  };  
  /*5*/
  vCopyVec3D(a3fOrigMosaicity,m_a3fCrysMosaicity);
  vSubVec3DVec3D(m_a3fCrysMosaicity,m_a3fCrysMosaicityShift,m_a3fCrysMosaicityShift);

  // Now do the sigmas
  // Note:  If there is a non-zero value for a sigma, then it should not be made zero,
  //         since multiple -go statements call nRefineLoop() multiple times.  The sigmas
  //         from previous loops should not be overwritten with zero.
  //


  // Cell sigmas are a little more tricky because they are in reciprocal angstroms rather than angstroms.
  // To solve this problem, we calculate a cell with sigma values added, and see how much it differs
  // from the refined cell.  The real mathematical way to do this is to use the equations transforming Crystal_Param* <==> Crystal_Param
  // (i.e. a* = b*c*sin(Alp)/Volume etc...).
  // This would allow the calculation of sigmas by propagation of errors.
  // Presently, I am too lazy to do this.

  Ccrystal oModifiedCrystal;
  m_ppoCrystals[m_nTwinID-1]->vGetRecipCell(a6fTemp);
  for (ny = 0; ny < nNumCrystalParams; ny++)
    {
      for (nx = 0; nx < 6; nx++)
        {
          if ((afSigma[CELL_A+ny] != 0.0) && (m_a6x6fCrystalVecs[ny][nx] != 0.0))
            a6fTemp[nx] += m_a6x6fCrystalVecs[ny][nx]*afSigma[CELL_A+ny];
        }
    }
  oModifiedCrystal.vSetRecipCell(a6fTemp);
  oModifiedCrystal.vGetCell(a6fTemp);
  vZeroMat(6,1,m_a6fCrysUnitCellSigs);
  for (ny = 0; ny < nNumCrystalParams; ny++)
    {
      for (nx = 0; nx < 6; nx++)
        {
          if ((m_a6x6fCrystalVecs[ny][nx] != 0.0))
            m_a6fCrysUnitCellSigs[nx] = fabs(a6fTemp[nx] - m_a6fCrysUnitCell[nx]);
        }
  }

  for (nx = 0; nx < 3; nx++)
    {
      if (afSigma[nx+CRYS_ROT1])
        m_a3fCrysRotSigs[nx] = afSigma[nx+CRYS_ROT1];
    }

  // These sigmas are difficult, since we might be refining in Modified cell parameters
  // I have not yet made these sigma values correct.
  for (nx = 0; nx < 3; nx++)
    {
      if (afSigma[nx+DET_TRANS1])
        m_a10x3fDetTransSigs[0][nx] = afSigma[nx+DET_TRANS1];
    }

  for (nx = 0; nx < 3; nx++)
    {
      if (afSigma[nx+DET_ROT1])
        m_a10x3fDetRotSigs[0][nx] = afSigma[nx+DET_ROT1];
    }

  for (nx = 0; nx < 2; nx++)
    {
      if (afSigma[nx+SOURCE_ROT1])
        m_a2fSourceRotSigs[nx] = afSigma[nx+SOURCE_ROT1];
    }

  if (afSigma[SOURCE_WAVE])
    m_fSourceWaveSig = afSigma[SOURCE_WAVE];


  // Store rejection information fields for display.
  for (nRef = 0; nRef<nNumReflns; nRef++)
    {
      (*m_poReflnlist)[nRef].vSetField(m_poReflnlist->m_nFI_nNonunfFlag,
                                       pnRejectFlag[nRef]);
    }

  // Save the residuals in arrays.
  m_poComp[m_nTwinID-1].m_fRMSDEG = m_fRMSDEG;
  m_poComp[m_nTwinID-1].m_fRMSMM  = m_fRMSMM;
  m_poComp[m_nTwinID-1].m_fRMSANG = m_fRMSANG;
  for (nx=0;nx<7;nx++) {
      m_poComp[m_nTwinID-1].m_a7nReflnNumbers[nx] = m_a7nReflnNumbers[nx];
      m_poComp[m_nTwinID-1].m_a7fReflnNumbers[nx] = m_a7fReflnNumbers[nx];
  };

  // Display stuff.
  nUpdateDisplay();

  if (0 < m_nVerbose)
  {
      nx = nListResults();
      if (0 != nx)
      {
          cout << "ERROR in Crefine::nRefineLoop: problem with listing results\n" << flush;
      }
      cout << "Reflections in list:         " << nNumReflns;
      cout << "\nReflections accepted [blue]: " << (nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp)
           << " (Resolution limits of these: " << m_a2fResoUsed[0] << " to " 
           << m_a2fResoUsed[1] << " A)"
           << "\nReflections rejected [ red]: " << (nNumRejects - nNumDegRejects)
          << "   (|Xobs-Xcalc| >= " << m_a4fRejectLimits[0] << " mm), or"
          << "\n                                 (|Yobs-Ycalc| >= " << m_a4fRejectLimits[1] << " mm)"
          << "\nReflections rejected [ red]: " << (nNumDegRejects)
          << "   (|Rotobs-Rotcalc| >= "
          << m_a4fRejectLimits[3] << " deg)"
          << "\nReflections ignored [green]: " << nNumIgnored
          << "   (I/sigI < " << m_fSigma 
          << "\n                                 or outside resolution"
          << " of " << m_fResolutionMin << " to " << m_fResolutionMax << ","
          << "\n                                 or Lorentz factor too large,"
          << "\n                                 or sharpness > " << m_fSharpness << ", etc)\n"
          << endl << flush;
      if (nNumOtherComp) 
      {
          cout << "\nReflections in other twins:  " << nNumOtherComp;
          cout << endl << flush;
      }
      vLogResiduals(nLoopCt);
  } else {
      if (!m_sLogRefineLast.length())
          vPrintResidualLog(1);      
      vLogResiduals(nLoopCt);
      vPrintResidualLog(2);
  };

#if ((defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  ReportResults(nNumReflns,nNumRejects,nNumReflns - nNumRejects - nNumIgnored - nNumOtherComp,nNumIgnored);
#endif
  if (m_bPrintTables)
    nPrintTables(pnRejectFlag);

    for (nx = 0; nx < nNumParams; nx++) {
      delete[] ppfDeriv[nx];
    }

    delete[] pnRejectFlag;
    delete[] anFlopCount;
    delete[] anGonioMatNum;
        delete[] anTwinID;
    delete[] anLostCompetition;
    delete[] anCompete;
    delete[] pfSigs;
    delete[] pfVal;
    delete[] aa3x3fPhiRotMat;
    delete[] aa3fVDCO;
    delete[] aa3fResidVec;
    delete[] aa3fXProjVec;
    delete[] afResidVecLen;
    delete[] afObsRotMid;
    delete[] afObsRotSigma;
    delete[] afObsRotWidth;
    delete[] afObsRotStart;
    delete[] afCalcRotMid;
    delete[] afObsXmm;
    delete[] afObsYmm;
    delete[] afObsZmm;
    delete[] afObsIntensity;
    delete[] afObsIntensitySig;
    delete[] afObsSharpness;

  return (0);
};


void
Crefine::vUsed2RealVarNum(int x, int *py)
{
 int _nx;
 int y;
 for (_nx=0,y=0; ((_nx < x) || (!m_anSigmaFlags[y])); y++)
   {
     if (m_anSigmaFlags[y])
       {
         _nx++;
       }
   }
 *py = y;
}


void
Crefine::vGetRMS(float *pfRMSDeg, float *pfRMSMM, float *pfRMSInvAngstrom)
{
  if (NULL != pfRMSDeg)
    *pfRMSDeg = m_fRMSDEG;
  if (NULL != pfRMSMM)
    *pfRMSMM = m_fRMSMM;
  if (NULL != pfRMSInvAngstrom)
    *pfRMSInvAngstrom = m_fRMSANG; 
  return;
} 





int Crefine::nPrintTables(int* pnReject)
{
    for(int nx = 0; nx < m_poReflnlistIn->nGetNumReflns(); nx++)
    {
        if (pnReject[nx] == 1)
            (*m_poReflnlistIn)[nx].m_nSelect = 1;
        else
            (*m_poReflnlistIn)[nx].m_nSelect = 0;
    }
    
    m_poReflnlistIn->m_pcSelect[m_poReflnlistIn->m_nTotalFieldsPlus1-1] = 0;
    m_poReflnlistIn->nResoTable(m_a2fResoUsed[0],
                                m_a2fResoUsed[1],
                                m_nNumberOfResoBins,
                                10,   // 10 I/Sigma bins
                                eTableVsResolution | eTableVsIoverSig,
                                -1,
                                NULL);

    return 0;
}

/// RB: it is possible that the input header already contains some "refinement results" entries.
// These entries need to be removed, so that if the refinement fails, they are not present in the header.
int Crefine::nDeleteHeaderEntries()
{
    if( !m_poHeader )
        return 1;

    m_poHeader->nDelete(Crefine::ms_sDtrefineRmsMm);
    m_poHeader->nDelete(Crefine::ms_sDtrefineRmsDeg);
    m_poHeader->nDelete(Crefine::ms_sDtrefineReflectionNumbers);
    m_poHeader->nDelete(Crefine::ms_sDtrefineReflectionNumbersPhotons);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Crefine::bIsDisplay()
{
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
    return (NULL != m_poXprop);
#else
    return (0 != m_nDisplay);
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


