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
// DTMainRank.cpp  Initial author: RB           15-Mar-2004
//
//  This file contains the member functions of class CDTMainRank
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

#include "DTMainRank.h"

#include "CCommandLine.h"
#include "Ccrystal.h"
#include "Crefine.h"


#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#endif

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


#define         SIZE_OF_PARAMS_ARRAY            100

static const char* s_pcDescriptionFormat = 
        "Ranking procedure awards/penalizes points based on a set of rules.\n"
        "A description of each rule that is used is given below.\n"
        "\n"
        "Rule  1: Spot count in resolution shells (found spots)\n"
        "Divide the total resolution range into [%.0lf] shells of equal volume.\n"
        "Consider [%.0lf] highest and [%.0lf] lowest shells. If there are [%.0lf]\n"
        "or more reflections found in a shell, a [%.0lf] point bonus is awarded.\n"
        "However, if there are less than %.0lf reflections found in the highest\n"
        "resolution shell, [%.0lf] points are subtracted from the award.\n"
        "\n"
        "Rule  2: I/Sigma in resolution shells (found spots)\n"
        "Divide the total resolution range into [%.0lf] shells of equal volume.\n"
        "Consider [%.0lf] highest and [%.0lf] lowest shells. If there are [%.0lf]\n"
        "or more reflections found in a shell,\n"
        "([%.0lf] * <I/sig>  - [%.0lf] * dev) / [%.2lf]\n"
        "points are awarded, where <I/sig> is the average I/sig value in that shell,\n"
        "and dev is a standard deviation of I/sig with respect to the rotation\n"
        "angle.\n"
        "\n"
        "Rule  3: Spot sharpness\n"
        "   Rule  3a: Spot asymmetry\n"  
        "   A penalty of [%.0lf] times the average spot profile asymmetry is made.\n"
        "   Rule  3b: Spot intensity per pixel\n"
        "   An award of [%.0lf] times <Ipix>/[%.0lf] is made, where Ipix is an average\n"
        "   peak intensity per pixel\n"
        "   Rule  3c: Spot anisotropy\n"
        "   A penalty of [%.0lf] times the average spot anisotropy is made.\n"
        "   Rule  3d: Spot amorphousness\n"
        "   A penalty of [%.0lf] times the average spot amorphousness is made.\n"
        "\n"
        "Rule  4: Strong ice rings\n"
        "For each strong (peak to background ratio > [%.2lf] percent) ice ring, [%.0lf]\n" 
        "points are deducted.\n"
        "\n"
        "Rule  5: Diffuse ice rings\n"
        "For each diffuse (peak to background ratio > [%.2lf] percent) ice ring, [%.0lf]\n"
        "points are deducted.\n"
        "\n"
        "Rule  6: Percentage of spots indexed\n"
        "If more than [%.0lf] reflections are indexed,\n"
        "[%.0lf] * <Indexed reflections / Used Reflections> points are awarded.\n"
        "If less than %.0lf or less than [%.0lf] percent reflections are indexed,\n"
        "[%.0lf] points are subtracted.\n"
        "\n"
        "Rule  7: RMS residual after refinement\n"
        "If a sample fails to refine, or has an RMS residual > [%.2lf] mm,\n"
        "%.0lf points are deducted. Otherwise, a penalty of\n"
        "[%.0lf] * <RMS-MM-Residual> is made.\n"
        "\n"
        "Rule  8: Mosaicity\n"
        "If a sample fails to refine, or has a mosaicity > [%.2lf] deg,\n"
        "%.0lf points are deducted. Otherwise, a penalty of\n"
        "[%.0lf] * Mosaicity is made.\n"
        "\n"
        "Rule  9: Percentage of spots refined\n"
        "If more than [%.0lf] reflections are refined,\n"
        "[%.0lf] * <Refined reflections / Predicted reflections> points are awarded.\n"
        "If less than %.0lf or less than [%.0lf] percent reflections are refined,\n"
        "[%.0lf] points are subtracted.\n"
        "\n"
        "Rule 10: Spot count in resolution shells (predicted and found spots)\n"
        "After refinement, divide the total resolution range into [%.0lf] shells of\n"
        "equal volume. Consider [%.0lf] highest and [%.0lf] lowest shells. If there are\n"
        "[%.0lf] or more predicted reflections found in a shell, a [%.0lf] point\n"
        "bonus is awarded. However, if there are less than %.0lf predicted reflections\n"
        "found in the highest resolution shell, [%.0lf] points are subtracted from\n"
		"the award.\n"
        "\n"
        "Rule 11: I/Sigma in resolution shells (predicted and found spots)\n"
        "After refinement, divide the total resolution range into [%.0lf] shells of\n"
        "equal volume. Consider [%.0lf] highest and [%.0lf] lowest shells. If there are\n"
        "[%.0lf] or more reflections found in a shell,\n" 
        "([%.0lf] * <I/sig>  - [%.0lf] * dev) / [%.2lf]\n"
        "points are awarded, where <I/sig> is the average I/sig value in that shell,\n"
        "and dev is a standard deviation of I/sig with respect to the rotation\n"
        "angle.\n"
        "\n"
        "Rule 12: Strategy total rotation width\n"
        "A penalty of [%.0f] times RotWidth/[%.2lf] is applied, where RotWidth is the\n"
        "total rotation width suggested by strategy.\n"
        "\n"
        "Each number in braces is a modifiable parameter which can \n"
        "be altered by using the -param command.  If for example\n"
        "you wish to modify the first parameter in Rule (3) to be 100, \n"
        "you specify -param 3 1 100.\n\0"
        ;

// RULE 1
double fMinSpotsInResoRing_Total_Bins = 10;
double fMinSpotsInResoRing_High_Bins = 7;
double fMinSpotsInResoRing_Low_Bins = 0;
double fMinSpotsInResoRing_Value = 5.0;
double fMinSpotsInResoRing_Points = 10;
double fMinSpotsInResoRing_Penalty_Points = 0;

// RULE 2
double fMinIoverSigInResoRing_Total_Bins = 10;
double fMinIoverSigInResoRing_High_Bins = 7;
double fMinIoverSigInResoRing_Low_Bins = 0;
double fMinIoverSigInResoRing_SpotsCount_Value = 5.0;
double fMinIoverSigInResoRing_Points = 20;
double fMinIoverSigRotDevInResoRing_Points = 20;
double fMinIoverSigInResoRing_Value = 60.0;

// RULE 3
double fSharpness_Asymmetry_Points     = 10;     // peak asymmetry     
double fSharpness_Density_Points      = 10;      // peak density
double fSharpness_Density_Value       = 300;     // peak density value
double fSharpness_Anisotropy_Points   = 10;      // spot anisotropy
double fSharpness_Irregularity_Points = 100;    // spot irregularity

// RULE 4
double fIceRingStrong_Value = 4.0;
double fIceRingStrong_Points = 10;

// RULE 5
double fIceRingMinimum_Value = 1.0;
double fIceRingDiffuse_Points = 5;

// RULE 6
double fMinIndexed_Spots_Count_Value = 20.0;
double fIndexed_Spots_Percentage_Points = 100.0;
double fIndexed_Spots_Minimum_Percentage = 0.0;
double fIndexed_Spots_Penalty_Points = 0.0;

// RULE 7
double fRMSResidualMM_Bad_Value = 0.30;
double fRMSResidualMM_Bad_Points = 30;

// RULE 8
double fMosaicity_Bad_Value = 1.0;
double fMosaicity_Bad_Points = 70;

// RULE 9
double fCount_Spots_Min_Value = 20.0;
double fPercentage_Spots_Min_Points = 100;
double fRefined_Spots_Minimum_Percentage = 0.0;
double fRefined_Spots_Penalty_Points = 0.0;
                                                
// RULE 10
double fMinSpotsInResoRing_AfterRefinement_Total_Bins = 10;
double fMinSpotsInResoRing_AfterRefinement_High_Bins  = 7;
double fMinSpotsInResoRing_AfterRefinement_Low_Bins   = 0;
double fMinSpotsInResoRing_AfterRefinement_Value      = 5.0;
double fMinSpotsInResoRing_AfterRefinement_Points     = 10;
double fMinSpotsInResoRing_AfterRefinement_Penalty_Points = 0;

// RULE 11
double fMinIoverSigInResoRing_AfterRefinement_Total_Bins = 10;
double fMinIoverSigInResoRing_AfterRefinement_High_Bins  = 7;
double fMinIoverSigInResoRing_AfterRefinement_Low_Bins   = 0;
double fMinIoverSigInResoRing_SpotsCount_AfterRefinement_Value = 5.0;
double fMinIoverSigInResoRing_AfterRefinement_Points = 20;
double fMinIoverSigRotDevInResoRing_AfterRefinement_Points = 20;
double fMinIoverSigInResoRing_AfterRefinement_Value = 60.0;

// RULE 12
double fStrategy_Points = 10;
double fStrategy_Value  = 360;

static void vInitializeRankingParameters()
{
    fMinSpotsInResoRing_Total_Bins = 10;
    fMinSpotsInResoRing_High_Bins = 7;
    fMinSpotsInResoRing_Low_Bins = 0;
    fMinSpotsInResoRing_Value = 5;                              
    fMinSpotsInResoRing_Points = 10;
    fMinSpotsInResoRing_Penalty_Points = 0;
                                                            
    fMinIoverSigInResoRing_Total_Bins = 10;
    fMinIoverSigInResoRing_High_Bins = 7;
    fMinIoverSigInResoRing_Low_Bins = 0;
    fMinIoverSigInResoRing_SpotsCount_Value = 5;                
    fMinIoverSigInResoRing_Value = 60.0;                        
    fMinIoverSigInResoRing_Points = 20;                         
    fMinIoverSigRotDevInResoRing_Points = 20;                         

    fSharpness_Asymmetry_Points     = 10;    // peak asymmetry     
    fSharpness_Density_Points      = 10;    // peak density
    fSharpness_Density_Value       = 300.0; // peak density value
    fSharpness_Anisotropy_Points   = 10;    // spot anisotropy
    fSharpness_Irregularity_Points = 100;    // spot irregularity
                                                            
    fIceRingStrong_Value = 4.0;                                 
    fIceRingStrong_Points = 10;                                 
                                                            
    fIceRingMinimum_Value = 1.0;                                
    fIceRingDiffuse_Points = 5;                                 
                                                            
    fMinIndexed_Spots_Count_Value = 20.0;                       
    fIndexed_Spots_Percentage_Points = 100.0;                   
    fIndexed_Spots_Minimum_Percentage = 0.0;
    fIndexed_Spots_Penalty_Points = 0.0;

    fRMSResidualMM_Bad_Value = 0.30;                            
    fRMSResidualMM_Bad_Points = 30;                             
                                                            
    fMosaicity_Bad_Value = 1.0;                                 
    fMosaicity_Bad_Points = 70;                                 
                                                            
    fCount_Spots_Min_Value = 20.0;                              
    fPercentage_Spots_Min_Points = 100;                         
    fRefined_Spots_Minimum_Percentage = 0.0;
    fRefined_Spots_Penalty_Points = 0.0;
                                                            
    fMinSpotsInResoRing_AfterRefinement_Total_Bins = 10;
    fMinSpotsInResoRing_AfterRefinement_High_Bins  = 7;
    fMinSpotsInResoRing_AfterRefinement_Low_Bins   = 0;
    fMinSpotsInResoRing_AfterRefinement_Value = 5;              
    fMinSpotsInResoRing_AfterRefinement_Points = 10;            
    fMinSpotsInResoRing_AfterRefinement_Penalty_Points = 0;

    fMinIoverSigInResoRing_AfterRefinement_Total_Bins = 10;
    fMinIoverSigInResoRing_AfterRefinement_High_Bins  = 7;
    fMinIoverSigInResoRing_AfterRefinement_Low_Bins   = 0;
    fMinIoverSigInResoRing_SpotsCount_AfterRefinement_Value = 5;
    fMinIoverSigInResoRing_AfterRefinement_Value = 60.0;        
    fMinIoverSigInResoRing_AfterRefinement_Points = 20;
    fMinIoverSigRotDevInResoRing_AfterRefinement_Points = 20;
    
    fStrategy_Points = 10;
    fStrategy_Value  = 360;
}

const int apnParamRules[100] =
{
    1,1,1,1,1,1,
    
    2,2,2,2,2,2,2,
    
    3,3,3,3,3,
    
    4,4,
    
    5,5,
    
    6,6,6,6,   // indexing
    
    7,7,
    
    8,8,
    
    9,9,9,9,  // refinement
    
    10,10,10,10,10,10,
    
    11,11,11,11,11,11,11,

    12,12,
    
    0
};

enum eRankingRules
{
    enSpotCountFound=1,
    enIOverSigmaFound,
    enSpotShape,
    enStrongIceRings,
    enDiffuseIceRings,
    enSpotsIndexed,
    enRMS,  
    enMosaicity,    
    enSpotsRefined,
    enSpotCountPredictedFound,
    enIOverSigmaPredictedFound,
    enStrategy
};

// Attention! Changing array apfParamValues may require changing 
// "-param" commands, constructed in the CDTMainRanker class.

double* apfParamValues[SIZE_OF_PARAMS_ARRAY] =
{
    &fMinSpotsInResoRing_Total_Bins,                          // 1
    &fMinSpotsInResoRing_High_Bins,                           // 2
    &fMinSpotsInResoRing_Low_Bins,                            // 3
    &fMinSpotsInResoRing_Value,                               // 4
    &fMinSpotsInResoRing_Points,                              // 5
    &fMinSpotsInResoRing_Penalty_Points,                      // 6

    &fMinIoverSigInResoRing_Total_Bins,                       // 1
    &fMinIoverSigInResoRing_High_Bins,                        // 2
    &fMinIoverSigInResoRing_Low_Bins,                         // 3
    &fMinIoverSigInResoRing_SpotsCount_Value,                 // 4
    &fMinIoverSigInResoRing_Points,                           // 5
    &fMinIoverSigRotDevInResoRing_Points,                     // 6
    &fMinIoverSigInResoRing_Value,                            // 7

    &fSharpness_Asymmetry_Points,         
    &fSharpness_Density_Points,
    &fSharpness_Density_Value,
    &fSharpness_Anisotropy_Points,    
    &fSharpness_Irregularity_Points,
    
    &fIceRingStrong_Value,
    &fIceRingStrong_Points,
    
    &fIceRingMinimum_Value,
    &fIceRingDiffuse_Points,
    
    &fMinIndexed_Spots_Count_Value,
    &fIndexed_Spots_Percentage_Points,
    &fIndexed_Spots_Minimum_Percentage,
    &fIndexed_Spots_Penalty_Points,
    
    &fRMSResidualMM_Bad_Value,                               
    &fRMSResidualMM_Bad_Points,                              
    
    &fMosaicity_Bad_Value,                                   
    &fMosaicity_Bad_Points,                                  
    
    &fCount_Spots_Min_Value,
    &fPercentage_Spots_Min_Points,
    &fRefined_Spots_Minimum_Percentage,
    &fRefined_Spots_Penalty_Points,

    &fMinSpotsInResoRing_AfterRefinement_Total_Bins,
    &fMinSpotsInResoRing_AfterRefinement_High_Bins,
    &fMinSpotsInResoRing_AfterRefinement_Low_Bins,
    &fMinSpotsInResoRing_AfterRefinement_Value,
    &fMinSpotsInResoRing_AfterRefinement_Points,
    &fMinSpotsInResoRing_AfterRefinement_Penalty_Points,

    &fMinIoverSigInResoRing_AfterRefinement_Total_Bins,
    &fMinIoverSigInResoRing_AfterRefinement_High_Bins,
    &fMinIoverSigInResoRing_AfterRefinement_Low_Bins,
    &fMinIoverSigInResoRing_SpotsCount_AfterRefinement_Value,
    &fMinIoverSigInResoRing_AfterRefinement_Points,
    &fMinIoverSigRotDevInResoRing_AfterRefinement_Points,
    &fMinIoverSigInResoRing_AfterRefinement_Value,

    &fStrategy_Points,
    &fStrategy_Value
};
////////////////////////////////////////////////////////////////////////

CDTMainRank::CDTMainRank()
{    
    m_nCumulativePointCount = 0;
    m_nActiveRuleNumber = 0;
    m_nActiveRulePoints = 0;

    m_bRefineRankingInfoMustBeInHeader = false;

    m_strActiveRuleSubpoints = "";
    m_strSampleNameForReport = "";
}

int CDTMainRank::nExecute(unsigned int argc, char* argv[])
{
  vInitializeRankingParameters();
  nFreeCommandLine();

  bool bHelp                      = false;
  
  bool bRank                      = false;
  bool bSort                      = false;
  bool bSortReverse               = false;
  int nArg                        = -1;
  int nStat;
  int nx,ny,nz;
  Cstring sRankLogFile;
  Cstring sHeaderName;
  itr<double> afParamMods;
        
  const char* apcOptions[] = 
  {
    "-START","s[sFileSpec]","r1",
    "-START","n",
    "-sort","n.reverse","s[sFileSpec]",
    "-sort","s[sFileSpec]",
    "-param","d[nRuleNumber]","d[nParamNumber]","f[fValue]","r3",
    "-help", "n",
    NULL
  };

  vDtrekSetModuleName("dtrank");
  vPrintCopyrightInfo();
  //cout << "\ndtrank:  Copyright (c) 2006 Rigaku\n";
  //cout << D_K_DTREKVersion << endl;

  // Copy command line to output log
  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << flush;

  // Generate command options for the header
  for (int ii = 1; ii < argc; ii++)
  {
     m_strCommandLineOptions = m_strCommandLineOptions + ' ' + (const char*) argv[ii];
  }

  while (0 == (nStat = nReadCommandLine(nArg,argc,argv,apcOptions)))
    {
      if (g_sOptionArg == "-START")
        {
          if (g_asOptionBuf.size())
            {
              for (nx = 0; nx < g_asOptionBuf.size(); nx++)
                {
                  
                  Cimage_header oHeaderNew(g_asOptionBuf[nx]);
                  
                  if (!oHeaderNew.bIsAvailable())
                    {
                      printf("WARNING: loading '%s'\n",g_asOptionBuf[nx].string());
                    } 
                  else
                  {
                      if( !m_poHeader )
                          vCreateHeader();
                   
                      m_poHeader->nCopyMask(oHeaderNew,(Cstring) "*");
                  }

                  sHeaderName = g_asOptionBuf[nx];
                }
              bRank = true;
            } 
          else
            bRank = FALSE;
        } 
      else if (g_sOptionArg == "-param")
        {
          for (nx = 0; nx < g_afOptionBuf.size(); nx++)
            {
              afParamMods + ((double) g_anOptionBuf[nx*2 + 0]) 
                + ((double) g_anOptionBuf[nx*2 + 1]) + g_afOptionBuf[nx];
            }
        } 
      else if (g_sOptionArg == "-sort-reverse")
        {
          bSortReverse = true;
        } 
      else if (g_sOptionArg == "-sort")
          {
            sRankLogFile = g_asOptionBuf[0];
            bSort = true;
          }
      else if (g_sOptionArg == "-help")
          {
            bHelp = true;
          }
    }


    // Make modifications to the parameters.

  if (afParamMods.size())
    {
      for (nx = 0; nx + 2 < afParamMods.size(); nx += 3)
        {
            for (ny = 0;     (apnParamRules[ny] != 0) 
                          && (apnParamRules[ny] != afParamMods[nx]);ny++);
            if (apnParamRules[ny] != 0)
              {
                for (nz = 1;(apnParamRules[ny] != 0) && (nz != afParamMods[nx+1]);ny++,nz++);
              }
            if (apnParamRules[ny] == 0)
              {
                printf("WARNING:  Could not modify parameter %d in rule %d.\n",
                       (int) afParamMods[nx + 0],(int) afParamMods[nx + 1]);
              } 
            else 
              {
                *apfParamValues[ny] = afParamMods[nx +2];
              }            
        }
    }

    //////////////////////////////////////////////////////////////////
    if( bHelp )
    {
        vPrintDtRankDescription();
        nFreeCommandLine();
        return 0;
    }
    /////////////////////////////////////////////////////////////////
    
    //RB adding this safety feature in case the header name is wrong
    if( !m_poHeader )
        return 1;

    if (bRank)
    {
        Cstring     sRankKeywordsFromDtRefine("DTREFINE_RANK_NAMES");
        Cstring     sTemp("");
        m_bRefineRankingInfoMustBeInHeader = 0 == m_poHeader->nGetValue(sRankKeywordsFromDtRefine, &sTemp);
        
        ////////////////////////////////////////////
        vPrintDtRankDescription();
        //cout << flush;
        vRankToHeader(DTREK_RANK_START);
        vRankToOutput(-1, NULL, 0, 0, DTREK_RANK_START);
        ////////////////////////////////////////////

        Ccrystal        oCrystal(*m_poHeader);
        int             nNewPoints = 0;
        double          f0 = 0.0;
        double          f2 = 0.0;
        char            pcTitle[100];           
        
	    
        Cstring             sFind = "dtfind";
        vRankResolutionShellInfo(sFind,
                                 (int)fMinSpotsInResoRing_Total_Bins,                      
                                 (int)fMinSpotsInResoRing_High_Bins,   
                                 (int)fMinSpotsInResoRing_Low_Bins,

                                 (int)fMinSpotsInResoRing_Value,
                                 (int)fMinSpotsInResoRing_Points,
                                 (int)fMinIoverSigInResoRing_SpotsCount_Value,
                                 (int)fMinIoverSigInResoRing_Points,
                                 (int)fMinIoverSigRotDevInResoRing_Points,
                                 fMinIoverSigInResoRing_Value,
                                 (int)fMinSpotsInResoRing_Penalty_Points);
      
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        vRankPeakSharpness();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Rank ice ring info. dtrefine results have preference
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;

        Cstring       strDTREKModulePrefix = m_bRefineRankingInfoMustBeInHeader ? "DTREFINE_RANK" : "DTFIND_RANK";

        nx = -1;
        while ( 0 == m_poHeader->nGetValueDictionary(strDTREKModulePrefix, "RESO_RING_MID_RESO_*", f0, ++nx) )
        {
            if ( 0 == m_poHeader->nGetValueDictionary(strDTREKModulePrefix,"RESO_RING_PERCENT_*",f2,nx) )
            {
                if( f2 > fIceRingMinimum_Value )
                {
                    if( f2 < fIceRingStrong_Value ) // diffuse ring - skip
                        continue;

                    sprintf(pcTitle,"Penalty for %s ring (%.2lf%%) near resln. %.3lf", "strong",f2,f0);

                    nNewPoints = -(int)fIceRingStrong_Points;
                    vRankToOutput((int)enStrongIceRings, pcTitle, nNewPoints, -(int)fIceRingStrong_Points);
                    m_nActiveRulePoints += nNewPoints;
                }
            }
        }
        
        vRankToHeader();
      
        ///////////////////////////////////////
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;

        nx = -1;
        while ( 0 == m_poHeader->nGetValueDictionary(strDTREKModulePrefix, "RESO_RING_MID_RESO_*", f0, ++nx) )
        {
            if ( 0 == m_poHeader->nGetValueDictionary(strDTREKModulePrefix,"RESO_RING_PERCENT_*",f2,nx) )
            {
                if( f2 > fIceRingMinimum_Value )
                {
                    if( f2 > fIceRingStrong_Value ) // strong ring - skip
                        continue;

                    sprintf(pcTitle,"Penalty for %s ring (%.2lf%%) near resln. %.3lf","diffuse",f2,f0);

                    nNewPoints = -(int)fIceRingDiffuse_Points;
                    vRankToOutput((int)enDiffuseIceRings, pcTitle, nNewPoints, -(int)fIceRingDiffuse_Points);
                    m_nActiveRulePoints += nNewPoints;
                }
            }
       }
        
        vRankToHeader();
        /////////////////////////////////////

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Rank indexing results
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;
        
        bool    bApplyPenaltyForBadIndexingResults = (fMinIndexed_Spots_Count_Value     > 0.0 &&
                                                      fIndexed_Spots_Minimum_Percentage > 0.0 &&
                                                      fIndexed_Spots_Penalty_Points     > 0.0);

        double  afTemp[2] = {0.0,0.0};
        if ( 0 == m_poHeader->nGetValue(Cstring("DTINDEX_RANK_USED_REFLECTIONS"), 2 , afTemp) )
        {
            if( afTemp[0] > 0.0 && afTemp[1] >= fMinIndexed_Spots_Count_Value )
            {
                f0 = afTemp[1] / afTemp[0];  
                
                nNewPoints = (int)(f0 * fIndexed_Spots_Percentage_Points);
        
                sprintf(pcTitle, "Indexed %.0f%% of all reflections used in indexing", 
                                 f0 * 100.0);
                vRankToOutput((int)enSpotsIndexed, pcTitle, nNewPoints, (int)fIndexed_Spots_Percentage_Points);
                m_nActiveRulePoints += nNewPoints;
            }
            else if( bApplyPenaltyForBadIndexingResults )
            {
                nNewPoints = -(int)fIndexed_Spots_Penalty_Points;
                sprintf(pcTitle, "Unacceptable number and/or percentage of indexed reflections");
                vRankToOutput((int)enSpotsIndexed, pcTitle, nNewPoints, -(int)fIndexed_Spots_Penalty_Points);
                m_nActiveRulePoints += nNewPoints;
             }
        }
        else if( bApplyPenaltyForBadIndexingResults ) // indexing results not found
        {
            nNewPoints = -(int)fIndexed_Spots_Penalty_Points;
            sprintf(pcTitle, "Unacceptable number and/or percentage of indexed reflections");
            vRankToOutput((int)enSpotsIndexed, pcTitle, nNewPoints, -(int)fIndexed_Spots_Penalty_Points);
            m_nActiveRulePoints += nNewPoints;
        }
        
        vRankToHeader();
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Rank RMS residual results
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;

        if ( 0 == m_poHeader->nGetValue(D_K_DtrefineRmsMm, 1, &f0) )
        {
            if (f0 == 0.0)
                f0 = 100.0;
        } 
        else 
        {
            f0 = 100.0;
        }

        if( f0 > fRMSResidualMM_Bad_Value )
        {
            nNewPoints = -(int)(fRMSResidualMM_Bad_Points * fRMSResidualMM_Bad_Value);
            sprintf(pcTitle,"Penalty for RMS residual > %.2lf (or no refinement)", fRMSResidualMM_Bad_Value);
        } 
        else  
        {
            nNewPoints = -(int)(fRMSResidualMM_Bad_Points * f0);
            sprintf(pcTitle,"Penalty for RMS residual value of %.3lf",f0);
        }

        vRankToOutput((int)enRMS, pcTitle, nNewPoints, -(int)(fRMSResidualMM_Bad_Points * fRMSResidualMM_Bad_Value));
        m_nActiveRulePoints += nNewPoints;
        vRankToHeader();
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ///////////////////////////////////////////
        // Rank mosaicity results
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;

        if (oCrystal.bIsAvailable())
        {
            f0 = oCrystal.fGetMosaicity();
        } 
        else
            f0 = 100.0;

        if( f0 > fMosaicity_Bad_Value )
        {
            nNewPoints = -(int)(fMosaicity_Bad_Points * fMosaicity_Bad_Value);
            sprintf(pcTitle,"Penalty for Mosaicity > %.2f (or no refinement)", fMosaicity_Bad_Value);
        } 
        else  
        {
            nNewPoints = -(int)(fMosaicity_Bad_Points * f0);
            sprintf(pcTitle,"Penalty for Mosaicity value of %.2f", f0);
        }
        
        vRankToOutput((int)enMosaicity, pcTitle, nNewPoints, -(int)(fMosaicity_Bad_Points * fMosaicity_Bad_Value));
        m_nActiveRulePoints += nNewPoints;
        vRankToHeader();
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;

        bool    bApplyPenaltyForBadRefinementResults = (fCount_Spots_Min_Value            > 0.0 &&
                                                        fRefined_Spots_Minimum_Percentage > 0.0 &&
                                                        fRefined_Spots_Penalty_Points     > 0.0);

        double      afRefineReflectionNumbers[7] = {0.0};
        if( 0 == m_poHeader->nGetValue(D_K_DtrefineReflectionNumbers, 7, afRefineReflectionNumbers) )
        {
            if( afRefineReflectionNumbers[0] > 0.0 && afRefineReflectionNumbers[1] >= fCount_Spots_Min_Value )
            {
                f0 = afRefineReflectionNumbers[1] / afRefineReflectionNumbers[0];  
                
                nNewPoints = (int)(f0 * fPercentage_Spots_Min_Points);
            
                sprintf(pcTitle, "Refined %.0f%% of all reflections used in refinement", 
                                  f0 * 100.0);
                
                vRankToOutput((int)enSpotsRefined, pcTitle, nNewPoints, (int)fPercentage_Spots_Min_Points);
                m_nActiveRulePoints += nNewPoints;
            }
            else if( bApplyPenaltyForBadIndexingResults )
            {
                nNewPoints = -(int)fRefined_Spots_Penalty_Points;
                sprintf(pcTitle, "Unacceptable number and/or percentage of refined reflections");
                vRankToOutput((int)enSpotsRefined, pcTitle, nNewPoints, -(int)fRefined_Spots_Penalty_Points);
                m_nActiveRulePoints += nNewPoints;
             }
        }
        else if( bApplyPenaltyForBadIndexingResults ) // refinement results not found
        {
            nNewPoints = -(int)fRefined_Spots_Penalty_Points;
            
            sprintf(pcTitle, "Unacceptable number and/or percentage of refined reflections");
            
            vRankToOutput((int)enSpotsRefined, pcTitle, nNewPoints, -(int)fRefined_Spots_Penalty_Points);
            m_nActiveRulePoints += nNewPoints;
        }

        vRankToHeader();

	    Cstring         sRefine = "dtrefine";
        vRankResolutionShellInfo(sRefine,
                                 (int)fMinSpotsInResoRing_AfterRefinement_Total_Bins,                      
                                 (int)fMinSpotsInResoRing_AfterRefinement_High_Bins,   
                                 (int)fMinSpotsInResoRing_AfterRefinement_Low_Bins,

                                 (int)fMinSpotsInResoRing_AfterRefinement_Value,
                                 (int)fMinSpotsInResoRing_AfterRefinement_Points,
                                 (int)fMinIoverSigInResoRing_SpotsCount_AfterRefinement_Value,
                                 (int)fMinIoverSigInResoRing_AfterRefinement_Points,
                                 (int)fMinIoverSigRotDevInResoRing_AfterRefinement_Points,
                                 fMinIoverSigInResoRing_AfterRefinement_Value,
                                 (int)fMinSpotsInResoRing_AfterRefinement_Penalty_Points);

        
        vRankStrategyResults();
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        vRankToOutput(-1, NULL, 0, 0, DTREK_RANK_TALLY);
        vRankToHeader(DTREK_RANK_TALLY);
    }
        

  // Sort the ranking output.
  if (bSort)
    {
      FILE*                 pFIn = NULL;
      FILE*                 pFOut = NULL;
      int                   nx = 0;
      int                   ny = 0;
      const int             nBufLength = 200;
      char                  pcBuf[nBufLength+1];
      std::vector<Cstring>  asTables;
      itr<int>              anRank;
      itr<int>              anSort;
      Cstring               sCumulTable;
      Cstring               sTemp;
      Stringreg             oReg;

      pFIn = fopen(sRankLogFile.string(),"rt");

      if (pFIn)
        {
          while (!feof(pFIn))
            {
              fgets(pcBuf,nBufLength,pFIn);
              sTemp = pcBuf;
              if ((sTemp.contains('\n')) && (sTemp.contains("Rank for"))) {
                sCumulTable = "";
              }
              sCumulTable += sTemp;
              if (!oReg.nParse(sTemp,".Cumulative.#."))
                {
                  if (!oReg.nLoadInt(0,ny))
                    {
                      sCumulTable += "\n\n";
                      asTables.push_back(sCumulTable);
                      anRank + ny;
                      sCumulTable = "";
                    }
                }

            }
          fclose(pFIn);

                        // Sort the output.
          if (asTables.size())
            {
              for (nx = 0; nx < asTables.size(); nx++)
                anSort + nx;
              g_pnCmpInts = &anRank[0];
              qsort((void*) &anSort[0],asTables.size(),sizeof(int),int_cmp_rel);
              if (!bSortReverse)
                {
                  for (nx = 0; nx < anSort.size()/2;nx++)
                    std::swap(anSort[nx],anSort[anSort.size() - nx - 1]);
                }
              pFOut = fopen(sRankLogFile.string(),"wt");
              if (pFOut)
                {
                  vPrintDtRankDescription(pFOut);

                  for (nx =0; nx < asTables.size(); nx++)
                    {
                      fprintf(pFOut,"%s",asTables[anSort[nx]].string());
                    }
                  fclose(pFOut);
                } 
              else 
                {
                  printf("ERROR:  Could not write to file '%s'\n",sRankLogFile.string());
                }
            } 
          else 
            {
              printf("ERROR:  The file '%s' did not appear to be a dtrank log file.\n",sRankLogFile.string());
            }
          
        } 
      else 
        {
          printf("ERROR:  Could not open file '%s'\n",sRankLogFile.string());
        }
    }
  
  nFreeCommandLine();
  
  return 0;
}

void CDTMainRank::vError(const int nErrorNum, const Cstring& sMessage)
{
}

void CDTMainRank::vGetRulesDescription(Cstring& strDesc)
{
    char      acTemp[4096];
    sprintf(acTemp, s_pcDescriptionFormat,
            (!apfParamValues[0])?((double) 0.0):(*apfParamValues[0]),  // number of bins
            (!apfParamValues[1])?((double) 0.0):(*apfParamValues[1]),  // number of high bins
            (!apfParamValues[2])?((double) 0.0):(*apfParamValues[2]),  // number of low bins
            (!apfParamValues[3])?((double) 0.0):(*apfParamValues[3]),  // minimum number of spots 
            (!apfParamValues[4])?((double) 0.0):(*apfParamValues[4]),  // bonus
            (!apfParamValues[3])?((double) 0.0):(*apfParamValues[3]),  // minimum number of spots 
            (!apfParamValues[5])?((double) 0.0):(*apfParamValues[5]),  // penalty

            (!apfParamValues[6])?((double) 0.0):(*apfParamValues[6]),    // number of bins                 
            (!apfParamValues[7])?((double) 0.0):(*apfParamValues[7]),    // number of high bins
            (!apfParamValues[8])?((double) 0.0):(*apfParamValues[8]),    // number of low bins
            (!apfParamValues[9])?((double) 0.0):(*apfParamValues[9]),    // minimum number of spots
            (!apfParamValues[10])?((double) 0.0):(*apfParamValues[10]),  // <I/sig> coefficient
            (!apfParamValues[11])?((double) 0.0):(*apfParamValues[11]),  // <I/sig> dev coefficient
            (!apfParamValues[12])?((double) 0.0):(*apfParamValues[12]),  // denominator (scale) coefficient
            
            (!apfParamValues[13])?((double) 0.0):(*apfParamValues[13]), // asymmetry
            (!apfParamValues[14])?((double) 0.0):(*apfParamValues[14]), // Intensity per pixel coeff1
            (!apfParamValues[15])?((double) 0.0):(*apfParamValues[15]), // Intensity per pixel coeff 2
            (!apfParamValues[16])?((double) 0.0):(*apfParamValues[16]), // Anisotropy
            (!apfParamValues[17])?((double) 0.0):(*apfParamValues[17]), // amorphousness
            
            
            (!apfParamValues[18])?((double) 0.0):(*apfParamValues[18]),  // strong ring I/sigma
            (!apfParamValues[19])?((double) 0.0):(*apfParamValues[19]),  // strong ring points
            
            (!apfParamValues[20])?((double) 0.0):(*apfParamValues[20]),  // weak ring I/sigma
            (!apfParamValues[21])?((double) 0.0):(*apfParamValues[21]),  // weak ring points
            
            (!apfParamValues[22])?((double) 0.0):(*apfParamValues[22]),  // minimum points for indexing
            (!apfParamValues[23])?((double) 0.0):(*apfParamValues[23]),  // coefficient for indexing
            (!apfParamValues[22])?((double) 0.0):(*apfParamValues[22]),  // minimum points for indexing
            (!apfParamValues[24])?((double) 0.0):(*apfParamValues[24]),  // minimum percentage for indexing
            (!apfParamValues[25])?((double) 0.0):(*apfParamValues[25]),  // penalty for indexing

            // rmsd
            (!apfParamValues[26])?((double) 0.0):(*apfParamValues[26]),  // threshold residual
            (!apfParamValues[26] || !apfParamValues[27] )?((double) 0.0):(*apfParamValues[26] * *apfParamValues[27]),
            (!apfParamValues[27])?((double) 0.0):(*apfParamValues[27]),

            // mosaicity
            (!apfParamValues[28])?((double) 0.0):(*apfParamValues[28]),
            (!apfParamValues[28] || !apfParamValues[29] )?((double) 0.0):(*apfParamValues[28] * *apfParamValues[29]),
            (!apfParamValues[29])?((double) 0.0):(*apfParamValues[29]),

            (!apfParamValues[30])?((double) 0.0):(*apfParamValues[30]),  // minimum points for refinement
            (!apfParamValues[31])?((double) 0.0):(*apfParamValues[31]),  // coefficient for refinement
            (!apfParamValues[30])?((double) 0.0):(*apfParamValues[30]),  // minimum points for refinement
            (!apfParamValues[32])?((double) 0.0):(*apfParamValues[32]),  // minimum percentage for refinement
            (!apfParamValues[33])?((double) 0.0):(*apfParamValues[33]),  // penalty for refinement

            (!apfParamValues[34])?((double) 0.0):(*apfParamValues[34]),  // number of bins
            (!apfParamValues[35])?((double) 0.0):(*apfParamValues[35]),  // number of high bins
            (!apfParamValues[36])?((double) 0.0):(*apfParamValues[36]),  // number of low bins 
            (!apfParamValues[37])?((double) 0.0):(*apfParamValues[37]),  // minimum number of spots  
            (!apfParamValues[38])?((double) 0.0):(*apfParamValues[38]),  // bonus 
            (!apfParamValues[37])?((double) 0.0):(*apfParamValues[37]),  // minimum number of spots  
            (!apfParamValues[39])?((double) 0.0):(*apfParamValues[39]),  // penalty 
            
            (!apfParamValues[40])?((double) 0.0):(*apfParamValues[40]),   // number of bins                 
            (!apfParamValues[41])?((double) 0.0):(*apfParamValues[41]),   // number of high bins
            (!apfParamValues[42])?((double) 0.0):(*apfParamValues[42]),   // number of low bins
            (!apfParamValues[43])?((double) 0.0):(*apfParamValues[43]),   // minimum number of spots
            (!apfParamValues[44])?((double) 0.0):(*apfParamValues[44]),   // <I/sig> coefficient
            (!apfParamValues[45])?((double) 0.0):(*apfParamValues[45]),   // <I/sig> dev coefficient
            (!apfParamValues[46])?((double) 0.0):(*apfParamValues[46]),   // denominator (scale) coefficient
            
            (!apfParamValues[47])?((double) 0.0):(*apfParamValues[47]),    // strategy coeff
            (!apfParamValues[48])?((double) 0.0):(*apfParamValues[48])    // strategy scal coeff
          );

    strDesc = acTemp; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  CDTMainRank::vPrintDtRankDescription(FILE* pFile)
{
    char      acTemp[4096];
    char      c_szLine[] = "\n-----------------------------------------------------------------------\n";

    strcpy(acTemp, c_szLine);  

    Cstring     strDesc("");
    vGetRulesDescription(strDesc);

    strcat(acTemp, strDesc.string());
    strcat(acTemp, c_szLine);
    
    if( pFile )  
    {
        fprintf(pFile, acTemp);
        fflush(pFile);
    }
    else
    {
        printf(acTemp);
        cout << flush;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainRank::vRankResolutionShellInfo(Cstring& strDTREKModuleName,
                                           int    nTotalBins,   
                                           int    nHighBins,
                                           int    nLowBins, 

                                           int    nMinAcceptableReflnCountPerShell,
                                           int    nPointsAwardForReflnCountPerShell,
                                           int    nMinAcceptableReflnCountPerShellForIOverSigmaCalc,
                                           int    nPointsAwardForIOverSigPerShell,
                                           int    nPointsAwardForIOverSigRotDevPerShell,
                                           double dIOverSigScaleCoefficient,
                                           int    nPointsPenaltyMaxReso)
{
    if( nHighBins > nTotalBins || nLowBins > nTotalBins || 
        nTotalBins < 1 || nLowBins < 0 || nHighBins < 0 )
        return;  // safety
    
    Cstring                     strTemp("");
    int                         nx = 0;
    int                         nNewPoints = 0;
    char                        pcTitle[256];
    int                         nRuleNumberSpotCount = -1;
    int                         nRuleNumberIOverSigma = -1;
    
    /////////////////////////////////////////////////////////
    Cstring                     sDTREKModuleReference(""); 
    if( strDTREKModuleName=="dtfind" )
    {
        sDTREKModuleReference = "Found:";
        nRuleNumberSpotCount = (int)enSpotCountFound;
        nRuleNumberIOverSigma = (int)enIOverSigmaFound;
    }
    else if ( strDTREKModuleName=="dtrefine" )
    {
        sDTREKModuleReference = "Predicted&Found:";
        nRuleNumberSpotCount = (int)enSpotCountPredictedFound;
        nRuleNumberIOverSigma = (int)enIOverSigmaPredictedFound;
    }
    /////////////////////////////////////////////////////////

    strDTREKModuleName.upcase();
    strDTREKModuleName += "_RANK";

    // Read shell info from header and fill in an array with information 
    RESOLUTION_SHELL_INFO*      astResoShellInfo = new RESOLUTION_SHELL_INFO[nTotalBins];
    
    bool        bHeaderHasShellInfo = true;

    for(nx = 0; nx < nTotalBins; nx++)
    {
        astResoShellInfo[nx].vInitialize();
        if( 0 != m_poHeader->nGetValueDictionary(strDTREKModuleName,"RESO_SHELL_*", strTemp, nx) )
        {
            if( 0 == nx )
            {
                bHeaderHasShellInfo = false;
                break; // no information available in header!
            }
            else
                continue;
        }    
        /////////////////////////////////////////
        // Parse header value string
        if( !astResoShellInfo[nx].bSet(strTemp) )
            continue;
    }
    
    ///////////////////////////////////////////////////////////
    if( !bHeaderHasShellInfo )  // just ouput zeros to header
    {
        /// RULE peak count
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;
        vRankToHeader();
        
        /// RULE I/Sigma
        m_nActiveRulePoints = 0;
        m_nActiveRuleNumber++;
        vRankToHeader();

        delete [] astResoShellInfo;
        
        return;
    }
    ////////////////////////////////////////////////////////////

    m_nActiveRulePoints = 0;
    m_nActiveRuleNumber++;

    //////////////////////
    // Do Peak Count Rule

    // Do High Reso Shells
    for(nx = nTotalBins-1; nx >= nTotalBins - nHighBins; nx--)
    {
        if( astResoShellInfo[nx].m_nReflnCount >= nMinAcceptableReflnCountPerShell )
        {
           sprintf(pcTitle,"%s >=%d reflns in HRS #%d (%.2f-%.2f)Å",
                           sDTREKModuleReference.string(),
                           nMinAcceptableReflnCountPerShell,
                           nx+1, 
                           astResoShellInfo[nx].m_dResoStart_A, 
                           astResoShellInfo[nx].m_dResoEnd_A);
           
           vRankToOutput(nRuleNumberSpotCount, pcTitle, nPointsAwardForReflnCountPerShell, nPointsAwardForReflnCountPerShell);
           m_nActiveRulePoints += nPointsAwardForReflnCountPerShell;
        }
    }
    
    // Do the highest shell penalty
    if( nPointsPenaltyMaxReso > 0 && astResoShellInfo[nTotalBins-1].m_nReflnCount < nMinAcceptableReflnCountPerShell )
    {
        sprintf(pcTitle,"%s <%d reflns in HRS #%d (%.2f-%.2f)Å",
                       sDTREKModuleReference.string(),
                       nMinAcceptableReflnCountPerShell,
                       nTotalBins, 
                       astResoShellInfo[nTotalBins-1].m_dResoStart_A, 
                       astResoShellInfo[nTotalBins-1].m_dResoEnd_A);

        vRankToOutput(nRuleNumberSpotCount, pcTitle, -nPointsPenaltyMaxReso, -nPointsPenaltyMaxReso);
        m_nActiveRulePoints -= nPointsPenaltyMaxReso;
    }

    // Do Low Reso Shells
    for(nx = 0; nx < nLowBins; nx++)
    {
        if( astResoShellInfo[nx].m_nReflnCount >= nMinAcceptableReflnCountPerShell )
        {
           sprintf(pcTitle,"%s >=%d reflns in LRS #%d (%.2f-%.2f)Å",
                           sDTREKModuleReference.string(),
                           nMinAcceptableReflnCountPerShell,
                           nx+1, 
                           astResoShellInfo[nx].m_dResoStart_A, 
                           astResoShellInfo[nx].m_dResoEnd_A);
           
           vRankToOutput(nRuleNumberSpotCount, pcTitle, nPointsAwardForReflnCountPerShell, nPointsAwardForReflnCountPerShell);
           m_nActiveRulePoints += nPointsAwardForReflnCountPerShell;
        }
    }

    vRankToHeader();

    //////////////////////
    // Do I/Sigma Rule
    double      dTemp = 0.0;

    m_nActiveRulePoints = 0;
    m_nActiveRuleNumber++;
    
    bool        bRankAnisotropy = nPointsAwardForIOverSigRotDevPerShell > 0;


    // Do High Reso Shells
    for(nx = nTotalBins - 1; nx >= nTotalBins - nHighBins; nx--)
    {
        if( astResoShellInfo[nx].m_nReflnCount < nMinAcceptableReflnCountPerShellForIOverSigmaCalc )
            continue;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // If I/Sigma e.s.d. is negative, it means we could not evaluate it. Therefore, no award will be given at all.
        // Note: not being able to evaluate the e.s.d. here means not having enough reflections in all of the "rotation bins".
        // For example, we are ranking 2 images 45 degrees apart (so we have 2 rotation bins), and our minimum statistically-accepted 
        // number of reflections is 5. Suppose, we have 5 reflections in a resolution shell, 2 on the first image and 3 - on the second. 
        // Even though "5" is enough to award the resolution shell for the presence of reflections there (Peak Count Rule), 
        // still we do not have enough reflections (2 and 3) to evaluate the mean I/sigma for either of the rotation bins. 
        // For details refer to CResoRankingBin::vCalcIOverSigmaMeanAndRotationDev().
        if( bRankAnisotropy && astResoShellInfo[nx].m_dIOverSigDev < 0.0 )
            continue;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        dTemp = (double)nPointsAwardForIOverSigPerShell * astResoShellInfo[nx].m_dIOverSig -
                (double)nPointsAwardForIOverSigRotDevPerShell * astResoShellInfo[nx].m_dIOverSigDev;

        nNewPoints = min( nPointsAwardForIOverSigPerShell, (int)(dTemp / dIOverSigScaleCoefficient) ); 

        if( nNewPoints < 0 )  // The award for I/Sigma cannot be negative
            nNewPoints = 0;

        if( bRankAnisotropy )  // print I/Sigma's with e.s.d's in brackets
            sprintf(pcTitle,"%s I/S=%.1f(%.1f) in HRS #%d (%.2f-%.2f)Å",
                            sDTREKModuleReference.string(),
                            astResoShellInfo[nx].m_dIOverSig,
                            astResoShellInfo[nx].m_dIOverSigDev,
                            nx+1,
                            astResoShellInfo[nx].m_dResoStart_A, 
                            astResoShellInfo[nx].m_dResoEnd_A);
        else
            sprintf(pcTitle,"%s I/S=%.1f in HRS #%d (%.2f-%.2f)Å",
                            sDTREKModuleReference.string(),
                            astResoShellInfo[nx].m_dIOverSig,
                            nx+1,
                            astResoShellInfo[nx].m_dResoStart_A, 
                            astResoShellInfo[nx].m_dResoEnd_A);

        vRankToOutput(nRuleNumberIOverSigma, pcTitle, nNewPoints, nPointsAwardForIOverSigPerShell);
        m_nActiveRulePoints += nNewPoints;
    }
    
    // Do Low Reso Shells
    for(nx = 0; nx < nLowBins; nx++)
    {
        if( astResoShellInfo[nx].m_nReflnCount < nMinAcceptableReflnCountPerShellForIOverSigmaCalc )
            continue;
        
        // If I/Sigma e.s.d. is negative, it means we could not evaluate it. Therefore, no award will be given at all.
        if( bRankAnisotropy && astResoShellInfo[nx].m_dIOverSigDev < 0.0 )
            continue;
        
        dTemp = (double)nPointsAwardForIOverSigPerShell * astResoShellInfo[nx].m_dIOverSig -
                (double)nPointsAwardForIOverSigRotDevPerShell * astResoShellInfo[nx].m_dIOverSigDev;

        nNewPoints = min( nPointsAwardForIOverSigPerShell, (int)(dTemp / dIOverSigScaleCoefficient) ); 

        if( nNewPoints < 0 )  // The award for I/Sigma cannot be negative
            nNewPoints = 0;

        if( bRankAnisotropy )  // print I/Sigma's with e.s.d's in brackets
            sprintf(pcTitle,"%s I/S=%.1f(%.1f) in LRS #%d (%.2f-%.2f)Å",
                            sDTREKModuleReference.string(),
                            astResoShellInfo[nx].m_dIOverSig,
                            astResoShellInfo[nx].m_dIOverSigDev,
                            nx+1,
                            astResoShellInfo[nx].m_dResoStart_A, 
                            astResoShellInfo[nx].m_dResoEnd_A);
        else
            sprintf(pcTitle,"%s I/S=%.1f in LRS #%d (%.2f-%.2f)Å",
                            sDTREKModuleReference.string(),
                            astResoShellInfo[nx].m_dIOverSig,
                            nx+1,
                            astResoShellInfo[nx].m_dResoStart_A, 
                            astResoShellInfo[nx].m_dResoEnd_A);

        vRankToOutput(nRuleNumberIOverSigma, pcTitle, nNewPoints, nPointsAwardForIOverSigPerShell);
        m_nActiveRulePoints += nNewPoints;
    }

    vRankToHeader();
    
    delete [] astResoShellInfo;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Rank sharpness info. dtrefine results have preference
void CDTMainRank::vRankPeakSharpness()
{
    if( !m_poHeader )
        return; //safety
    
    int                         nNewPoints = 0;
    char                        pcTitle[256];

    m_nActiveRulePoints = 0;
    m_nActiveRuleNumber++;

    m_strActiveRuleSubpoints = "";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( fSharpness_Asymmetry_Points > 0.0 )
    {
        double      dProfileAsymmetry = 0.0;
        bool        bProfileAsymmetryAvailable = false;
    
        if( 0 != m_poHeader->nGetValueDictionary("DTREFINE_RANK","PEAK_PROFILE_ASYMMETRY", dProfileAsymmetry) )
        {
            if( !m_bRefineRankingInfoMustBeInHeader && 
                0 == m_poHeader->nGetValueDictionary("DTFIND_RANK",  "PEAK_PROFILE_ASYMMETRY", dProfileAsymmetry) )
            bProfileAsymmetryAvailable = true;
        }
        else
           bProfileAsymmetryAvailable = true;

        if( bProfileAsymmetryAvailable )
        {
            nNewPoints = max(-(int) fSharpness_Asymmetry_Points,-(int) floor(0.5 + (fSharpness_Asymmetry_Points * dProfileAsymmetry)));
            sprintf(pcTitle,"Penalty for peak profile asymmetry of %.2f", dProfileAsymmetry);       
        } 
        else
        {
            nNewPoints = -(int)fSharpness_Asymmetry_Points;
            sprintf(pcTitle,"Penalty for absent profile asymmetry info");
        }
        
        vRankToOutput((int)enSpotShape, pcTitle, nNewPoints, -(int)fSharpness_Asymmetry_Points);
        m_nActiveRulePoints += nNewPoints;                                                     
    
        m_strActiveRuleSubpoints += Cstring(nNewPoints);
    }
    else
        m_strActiveRuleSubpoints += Cstring(0);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( fSharpness_Density_Points > 0.0 )
    {
        double      dIntensityDensity = 0.0;
        bool        bIntensityDensityAvailable  = false;
        
        if( 0 != m_poHeader->nGetValueDictionary("DTREFINE_RANK","PEAK_INTENSITY_DENSITY", dIntensityDensity) )
        {
            if( !m_bRefineRankingInfoMustBeInHeader && 
                0 == m_poHeader->nGetValueDictionary("DTFIND_RANK",  "PEAK_INTENSITY_DENSITY", dIntensityDensity) )
            bIntensityDensityAvailable = true;
        }
        else
           bIntensityDensityAvailable = true;
        
        if( bIntensityDensityAvailable )
        {
            nNewPoints = min((int)fSharpness_Density_Points,
                             (int)(fSharpness_Density_Points * dIntensityDensity / fSharpness_Density_Value) );
            sprintf(pcTitle,"Award for peak intensity per pixel of %.2f", dIntensityDensity);       
        } 
        else
        {
            nNewPoints = -(int)fSharpness_Density_Points;
            sprintf(pcTitle, "Penalty for absent peak intensity density info");
        }
        
        vRankToOutput((int)enSpotShape,pcTitle, nNewPoints, (int)fSharpness_Density_Points);
        m_nActiveRulePoints += nNewPoints;                                                     
        
        m_strActiveRuleSubpoints += Cstring(" ");
        m_strActiveRuleSubpoints += Cstring(nNewPoints);
    }
    else
    {
        m_strActiveRuleSubpoints += Cstring(" 0");
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( fSharpness_Anisotropy_Points > 0.0 )
    {
        double      dShapeAnisotropy = 0.0;
        bool        bShapeAnisotropyAvailable = false;
        
        if( 0 != m_poHeader->nGetValueDictionary("DTREFINE_RANK", "PEAK_SHAPE_ANISOTROPY", dShapeAnisotropy) )
        {
            if( !m_bRefineRankingInfoMustBeInHeader && 
                0 == m_poHeader->nGetValueDictionary("DTFIND_RANK", "PEAK_SHAPE_ANISOTROPY", dShapeAnisotropy) )
            bShapeAnisotropyAvailable = true;
        }
        else
           bShapeAnisotropyAvailable = true;

        if( bShapeAnisotropyAvailable )
        {
            nNewPoints = max(-(int) fSharpness_Anisotropy_Points, -(int)(fSharpness_Anisotropy_Points * dShapeAnisotropy) );
            sprintf(pcTitle,"Penalty for spot shape anisotropy of %.2lf", dShapeAnisotropy);       
        } 
        else
        {
            nNewPoints = -(int)fSharpness_Anisotropy_Points;
            sprintf(pcTitle,"Penalty for absent spot shape anisotropy info");
        }
        
        vRankToOutput((int)enSpotShape, pcTitle, nNewPoints, -(int)fSharpness_Anisotropy_Points);
        m_nActiveRulePoints += nNewPoints;                                                     
        
        m_strActiveRuleSubpoints += Cstring(" ");
        m_strActiveRuleSubpoints += Cstring(nNewPoints);
    }
    else
    {
        m_strActiveRuleSubpoints += Cstring(" 0");
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( fSharpness_Irregularity_Points > 0.0 )
    {
        double      dShapeIrregularity = 0.0;
        bool        bShapeIrregularityAvailable = false;

        if( 0 != m_poHeader->nGetValueDictionary("DTREFINE_RANK", "PEAK_SHAPE_IRREGULARITY", dShapeIrregularity) )
        {
            if( !m_bRefineRankingInfoMustBeInHeader && 
                0 == m_poHeader->nGetValueDictionary("DTFIND_RANK", "PEAK_SHAPE_IRREGULARITY", dShapeIrregularity) )
            bShapeIrregularityAvailable = true;
        }
        else
           bShapeIrregularityAvailable = true;

        if( bShapeIrregularityAvailable )
        {
            nNewPoints = max(-(int) fSharpness_Irregularity_Points, -(int)(fSharpness_Irregularity_Points * dShapeIrregularity) );

            sprintf(pcTitle,"Penalty for spot shape amorphousness of %.3f", dShapeIrregularity);       
        } 
        else
        {
            nNewPoints = -(int)fSharpness_Irregularity_Points;
            sprintf(pcTitle,"Penalty for absent spot shape amorphousness info");
        }
        
        vRankToOutput((int)enSpotShape, pcTitle, nNewPoints, -(int)fSharpness_Irregularity_Points);
        m_nActiveRulePoints += nNewPoints;                                                     
        
        m_strActiveRuleSubpoints += Cstring(" ");
        m_strActiveRuleSubpoints += Cstring(nNewPoints);
    }
    else
    {
        m_strActiveRuleSubpoints += Cstring(" 0");
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vRankToHeader(DTREK_RANK_PRINT_SUBSCORES);                                                                   
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainRank::vRankStrategyResults()
{
    if( !m_poHeader )
        return; //safety
    
    int                         nNewPoints = 0;
    char                        pcTitle[256];

    m_nActiveRulePoints = 0;
    m_nActiveRuleNumber++;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( fStrategy_Points > 0.0 )
    {
        double      dRotWidth = 0.0;
        if( 0 == m_poHeader->nGetValueDictionary("DTMULTISTRATEGY_RANK","TOTAL_ROTATION_WIDTH", dRotWidth) )
        {
            nNewPoints = max(-(int)fStrategy_Points, -(int) floor(0.5 + (fStrategy_Points * dRotWidth/fStrategy_Value)));
            sprintf(pcTitle,"Penalty for rotation width from strategy of %.0f deg", dRotWidth);       
        } 
        else
        {
            nNewPoints = -(int)fStrategy_Points;
            sprintf(pcTitle,"Penalty for absent strategy info");
        }
        
        vRankToOutput((int)enStrategy, pcTitle, nNewPoints, -(int)fStrategy_Points);
        m_nActiveRulePoints += nNewPoints;                                                     
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    vRankToHeader();                                                                   
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainRank::vRankToOutput(int nRuleNumber,
                                char* pcDescription,
                                int nRulePoints, 
                                int nRuleMaxPoints,
                                DTREK_WORD wCtrl)
{   
    if( wCtrl == DTREK_RANK_START )
    {
        printf("--------------------------------------------------------------------------------\n");
        printf("Rule %-57s Pts Cumul    Max\n","Category");
        printf("--------------------------------------------------------------------------------\n");

        return;
    }
    else if( wCtrl == DTREK_RANK_TALLY )
    {
        Cstring     strTemp("Cumulative for ");
        strTemp += m_strSampleNameForReport;
        
        printf("--------------------------------------------------------------------------------\n");
        printf("     %-57s %3s %5d\n",strTemp.string(),"",m_nCumulativePointCount);
        printf("\n\n\n");
        
        return;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    m_nCumulativePointCount += nRulePoints;
    printf(" %2d  %-57s%4d %5d %6d\n", nRuleNumber, pcDescription, nRulePoints, m_nCumulativePointCount, nRuleMaxPoints);
    /////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainRank::vRankToHeader(DTREK_WORD wCtrl)
{
    if( wCtrl == DTREK_RANK_START )
    {
        m_poHeader->nClearDictionary("DTRANK_RANK", "*");
        m_poHeader->nReplaceValue("DTRANK_OPTIONS", m_strCommandLineOptions);

        return;
    }
    else if ( wCtrl == DTREK_RANK_TALLY )
    {
        m_poHeader->nReplaceValueDictionary("DTRANK_RANK", "CUMULATIVE", m_nCumulativePointCount);
        m_poHeader->nWrite(sDtrekGetPrefix() + "dtranker.head"); // Save the header
        
        return;
    }
    
    if( !(wCtrl & DTREK_RANK_PRINT_SUBSCORES) )
        m_poHeader->nReplaceValueDictionary("DTRANK_RANK", "RULE_*", m_nActiveRulePoints, m_nActiveRuleNumber);
    else
    {
        Cstring     strOutput(m_nActiveRulePoints);
        strOutput += " ";
        strOutput += m_strActiveRuleSubpoints;
        m_poHeader->nReplaceValueDictionary("DTRANK_RANK", "RULE_*", strOutput, m_nActiveRuleNumber);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
