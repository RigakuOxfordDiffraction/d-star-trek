//
// Copyright (c) 1997 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cstrategy.cc        Initial author: J.W. Pflugrath           10-May-1997
//  This file contains the member functions of class Cstrategy
//    which implements the strategy encapsulation of d*TREK.
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
#include <algorithm>

#include "Dtrek.h"
#include "Cstrategy.h"         // Class definition and prototypes
#include <time.h>
#include "Cstat.h"
#include "Crotation.h"
#include "Cscaleaverage.h"

#ifdef SSI_PC
#include "CrclHelper.h"
#endif


#define     MAX_SEARCH_STEP_SINGLE_SCAN     5.0
static bool bRangesCompareConsiderTargets(CstrategyRange& oFirstRange, 
                                          CstrategyRange& oSecondRange, 
                                          double dTargetCompl, 
                                          double dTargetRed,
                                          bool& bIsBetterCompleteness);
static void vSortRangeVector(std::vector<CstrategyRange>& vRanges, bool (*pfbCompare)(CstrategyRange&, CstrategyRange&));
static bool bRangesCompareFindShortest(CstrategyRange& first, CstrategyRange& second);
static bool bRangesCompareFindMostComplete(CstrategyRange& first, CstrategyRange& second);


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

int CstrategyRange::nPrint(Cstrategy& oStrategy)
{
    int nx,ny,nz;
    double fRotCumul;
    printf("Scan Solution\n"
           "--------------------------------------\n"
           "Completeness      %-8.2lf\n"
           "Redundancy        %-8.2lf +/- %8.2lf\n"
           "--------------------------------------\n",
           m_fCompleteness,
           m_fRedundancy,
           m_fRedundancyDev);


    fRotCumul = 0.0;
    
    //RB 05/24/04 By Thad's design, the array of scan numbers, m_anScanNumbers, starts with the fixed scans             
    int     nScanNumber = -1;
    for (nx = oStrategy.m_anFixedScans.size(); nx < m_anScanNumbers.size(); nx++)
    {
        nScanNumber = m_anScanNumbers[nx];
        
        printf("Scan              %-8.2d   (%sS%d_SCAN_*)\n"
               "Rotation Start    %-8.2lf\n"
               "Rotation End      %-8.2lf\n"
               "Rotation Width    %-8.2lf\n",
               nScanNumber,
               D_K_StrategyTestPrefix,
               nScanNumber,
               m_afRotStart[nx],
               m_afRotEnd[nx],
               m_afRotEnd[nx]-m_afRotStart[nx]);

        fRotCumul += m_afRotEnd[nx] - m_afRotStart[nx];

        if (nx - oStrategy.m_anFixedScans.size() < oStrategy.m_asScanGonioAxis.size())
        {
            printf("Rotation Axis:    %-8s\n",oStrategy.m_asScanGonioAxis[nScanNumber].string());
            ny = oStrategy.m_afScanGonioAngles.size()/oStrategy.m_afScanRotMax.size();
            printf("Crys Goniometer    ");
            for (nz = 0; nz < ny; nz++)
                printf("%c%.2lf",nz?' ':'(',oStrategy.m_afScanGonioAngles[nScanNumber*ny + nz]);
            printf(")\n");
        };
        printf("--------------------------------------\n");
    };
    
    printf("Cumulative Width: %-8.2lf\n",fRotCumul);
    printf("--------------------------------------\n");
    return 0;
};

bool CstrategyRange::operator < (CstrategyRange& oOther)
{
    if( m_bIsAvailable &&  !oOther.m_bIsAvailable )
        return true;
    
    if( !m_bIsAvailable && oOther.m_bIsAvailable )
        return false;
    
    if( !m_bIsAvailable && !oOther.m_bIsAvailable )
    {
        printf("ERROR:  Two uninitialized strategy ranges compared!\n");
        return false;
    }
    
    double      dThisTotalRot = dGetTotalRotWidth();
    double      dOtherTotalRot = oOther.dGetTotalRotWidth();
    
    if( dThisTotalRot == dOtherTotalRot )
    {
        if( m_fCompleteness == oOther.m_fCompleteness )
            return (m_fRedundancy > oOther.m_fRedundancy);
        else
            return (m_fCompleteness > oOther.m_fCompleteness);
    }
    else
        return (dThisTotalRot < dOtherTotalRot);
}

int Cstrategy::nAdjustPhiToZero(Cimage_header& oHeader, CstrategyRange& oTestRange)
{
    if( !m_bZeroOmegaWithPhiOffset )
        return 0;  // do nothing
    
    int         nScan = -1;

    // Figure out the crystal gonio axes
    Cgoniometer         oGonio(oHeader, Ccrystal::ms_sCrystalPrefix);
    if( !oGonio.bIsAvailable() )
        return 1;

    const char*     cpcOmega = "OMEGA";
    const char*     cpcPhi = "PHI";

    int         nOmegaIndex = oGonio.nGetNum(Cstring(cpcOmega));
    int         nPhiIndex   = oGonio.nGetNum(Cstring(cpcPhi));
    
    if( -1 == nPhiIndex || -1 == nOmegaIndex )
        return 1; // must have both omega and phi!
    
    // Figure out the chi(kappa) index, assuming the the max number of axes is 3
    int         nChiorKappaIndex = -1;
    for(int ii=0; ii < 3; ii++)
    {
        if( ii != nOmegaIndex && ii != nPhiIndex )
        {
            nChiorKappaIndex = ii;
            break;
        }
    }
    
    ////////////////////////////////
    double          a3fValues[3];
    double          a3x3fVecs[3][3];
    
    ////////////////////////////////

    Cstring         sRotAxisName("");
    //RB 05/24/04 By Thad's design, the array of scan numbers, m_anScanNumbers, starts with the fixed scans             
    for (int nScani = m_anFixedScans.size(); nScani < oTestRange.m_anScanNumbers.size(); nScani++)
    {
        nScan = oTestRange.m_anScanNumbers[nScani];
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Crotation           oRotation(oHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, nScan));
        if( !oRotation.bIsAvailable() )
        {
            printf("\nERROR: Cannot adjust phi and omega for scan with index: %d, because the rotation object cannot be created.", nScan);
            continue;
        }
        sRotAxisName = oRotation.sGetName();
        sRotAxisName.upcase();
        if( sRotAxisName != "OMEGA" )
        {
            printf("\nERROR: Cannot adjust phi and omega for scan with index: %d, because the rotation axis is not omega.", nScan);
            continue;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Cgoniometer         oGonio(oHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, Ccrystal::ms_sCrystalPrefix, nScan));
        if( !oGonio.bIsAvailable() )
        {
            printf("\nERROR: Cannot adjust phi and omega for scan with index: %d, because the crystal gonio object cannot be created.", nScan);
            continue;
        }
        
        oGonio.nGetDatum(3, &a3fValues[0]);

        if( nChiorKappaIndex != -1 && a3fValues[nChiorKappaIndex] != 0.0 )
        {
            printf("\nERROR: Cannot adjust phi and omega for scan with index: %d, because these axes are not parallel for that scan.", nScan);
            continue;
        }

        oGonio.nGetVectors(3, &a3x3fVecs[0][0]);
        
        // Check the rotation directions.
        // and adjust phi value accordingly
        if( ABS(fDot3D(a3x3fVecs[0], a3x3fVecs[2]) - 1.0) < 1e-4 ) // Turning in the same direction.
        {
            a3fValues[nPhiIndex] += oTestRange.m_afRotStart[nScani] - m_fZeroOmega;
        } 
        else if( ABS(fDot3D(a3x3fVecs[0], a3x3fVecs[2]) + 1.0) < 1e-4 )  // Turning in opposite directions.
        {
            a3fValues[nPhiIndex] += m_fZeroOmega - oTestRange.m_afRotStart[nScani];
        } 
        else 
        {
            printf("\nERROR: Cannot adjust phi and omega for scan with index: %d, because these axes are not parallel for that scan.", nScan);
            continue;
        }
        
        // Now adjust begin and end of the rotation range
        double      dOmegaShift = m_fZeroOmega - oTestRange.m_afRotStart[nScani];
        
        vShiftReflnRotValueInScan(nScan, dOmegaShift);
        
        oTestRange.m_afRotStart[nScani] = oTestRange.m_afRotStart[nScani] + dOmegaShift;                                       
        oTestRange.m_afRotEnd[nScani]   = oTestRange.m_afRotEnd[nScani]   + dOmegaShift;
        
        // Save adjusted phi value in the scan info object
        m_poScanInfo->bSetScanAxisValue(nScan, cpcPhi, a3fValues[nPhiIndex]);
        
        // Save the adjusted phi in the gonio object and update header
        oGonio.nSetDatum(3, &a3fValues[0]);
        oGonio.nUpdateHeader(&oHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, Ccrystal::ms_sCrystalPrefix, nScan));

        printf("\nINFO:  Adjusting phi to %.2lf so that omega range is %.2lf to %.2lf for scan %d\n",
               a3fValues[nPhiIndex], oTestRange.m_afRotStart[nScani], oTestRange.m_afRotEnd[nScani], nScan);
    }
    
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CstrategyRange::nUpdateHeader(Cimage_header& oHeader, Cstrategy& oStrategy)
{   
    //RB 05/24/04 By Thad's design, the array of scan numbers, m_anScanNumbers, starts with the fixed scans             
    int     nTotalScansForUpdate = m_anScanNumbers.size() - oStrategy.m_anFixedScans.size();

    if( 0 == nTotalScansForUpdate )
        return 0; // nothing to do
    
    int     nFirstScanForUpdate = oStrategy.m_anFixedScans.size(); 

    int     nStat = 0;
    int     nScan = 0;
    
    for(int nScani = nFirstScanForUpdate; nScani < m_anScanNumbers.size(); nScani++)
    {
        nScan = m_anScanNumbers[nScani];

        Crotation       oRotation(oHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, nScan));
        
        if( !oRotation.bIsAvailable() )
        {
            printf("WARNING:  Not writing scan number %.2d to header.\n", nScan);
            nStat = 1;
            break;
        }

        oRotation.vSetRotStart(m_afRotStart[nScani]);
        oRotation.vSetRotEnd(m_afRotEnd[nScani]);
        
        nStat = oRotation.nUpdateHeader(&oHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, nScan));
        if( 0 != nStat )
            break;
    }
    
    if( 0 == nStat )
    {
        oHeader.nReplaceValue(D_K_ScanSeqSelected, nTotalScansForUpdate, &m_anScanNumbers[nFirstScanForUpdate]);
        oHeader.nReplaceValue(D_K_DtstrategyPercentComplete, m_fCompleteness);
        oHeader.nReplaceValue(D_K_DtstrategyRedundancy, m_fRedundancy);
        oHeader.nReplaceValue(D_K_DtstrategyRedundancyDev, m_fRedundancyDev);
    }

    return nStat;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




//+Definitions, constants, and initialization of static member variables

Cstring Cstrategy::ms_sf2STLcu    = "f2STLcu";

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cstrategy::Cstrategy(Cimage_header *poHeader)
{
    m_oBestRange.vSetStrategyPtr(this);
    
    nInitValues(*poHeader);
}

Cstrategy::~Cstrategy()
{

  if (NULL != m_poCrystal)
    {
      delete m_poCrystal;
      m_poCrystal = NULL;
    }

  if (NULL != m_pnIndex)
    {
      delete [] m_pnIndex;
      m_pnIndex = NULL;
    }
  if (NULL != m_pnIndexA)
    {
      delete[] m_pnIndexA;
      m_pnIndexA = NULL;
    }
  if (NULL != m_pnIndex2)
    {
      delete [] m_pnIndex2;
      m_pnIndex2 = NULL;
    }

    if( !m_bExternalResoBinsPtr && m_poResoBins )
    {
        delete m_poResoBins;
        m_poResoBins = NULL;
    }

    if( !m_bExternalScanInfoPtr && m_poScanInfo )
    {
        delete m_poScanInfo;
        m_poScanInfo = NULL;
    }
} 

int Cstrategy::nInitValues(void)
{
  m_bFixRotStart = false;

  m_poResoBins = NULL;
  m_bExternalResoBinsPtr = false;

  m_poScanInfo = NULL;
  m_bExternalScanInfoPtr = false;
  m_dTotalRotWidthAllowed = 100000.0;          // just a large number
  m_nTotalNumberOfScansAllowed = 100000;       // just a large number

  m_dTotalRotWidth = 0.0;
    
  m_poCrystal           = NULL;
  m_poReflnlist         = NULL;
  m_pnIndex             = NULL;
  m_pnIndexA            = NULL;
  m_pnIndex2            = NULL;
  m_nVerbose            = 3;
  m_bAnom               = FALSE;
  m_bZeroOmegaWithPhiOffset = FALSE;
  m_fZeroOmega          = 0.0;
  
  m_nNumResoBins        = -1;
  
  m_fCompletenessMin    = 99.0; 
  
  m_fCompletenessTol    = 0.0;       // default completeness tolerance
  m_fRotationTol        = 10.0;
  m_dMinCompletenessIncreaseOnExpansion = 0.5; 

  m_fMinCompPerDegree   = 0.0;
  
  m_fRedundancyMin      =  1.0;
  
  m_fCompletenessTolOnFailure = 0.5;
  m_fRedundancyTolOnFailure = 0.3;

  m_a2nFirstScan[0]     = 0;
  m_a2nFirstScan[1]     = 100000;
  m_a2nSecondScan[0]    = -1;
  m_a2nSecondScan[1]    = -1;
  
  -m_anFixedScans;
  -m_anExcludeScans;

  m_fCellLengthFactor   =   1.0;
  m_fCellVolumeFactor   =   1.0;
  m_nFI_fRotStart       = -1;
  m_nFI_fRotEnd         = -1;
  m_nMaxTime            = 30;
  m_fPad                = 0.0;

  return 0;
}

bool Cstrategy::bIsValid(Cstring& sError)
{
    if( !(m_fCompletenessTol >= 0.0 && m_fCompletenessTol < m_fCompletenessMin) )
    {
        sError = "\n\nERROR: completeness tolerance must be greater or equal to 0, and less than your target completeness.\n\n";
        return false;
    }
    
    return true;
}


int Cstrategy::nInitValues(Cimage_header& oHeader)
{
  // TODO:  fill in results from image header

  int nStat;
  nStat = nInitValues();

 m_poCrystal     = new Ccrystal(oHeader);

 if (!m_poCrystal->bIsAvailable())
   {
     printf("ERROR in Cstrategy::nInitValues(), invalid crystal info!\n");
     nStat = -1;
#ifdef SSI_PC
     return (nStat);
#else
     exit (nStat);
#endif
   }
  return (nStat);
}

int Cstrategy::nList(DTREK_WORD wCtrl)
{
    if( wCtrl & DTREK_CSTRATEGY_LIST_TITLE )
        printf("Cstrategy:: listing:\n");

    if( wCtrl & DTREK_CSTRATEGY_LIST_CRYSTAL && NULL != m_poCrystal )
    {
        m_poCrystal->nList(1);
    }

    if( wCtrl & DTREK_CSTRATEGY_LIST_RESO && NULL != m_poResoBins )
    {
        printf("\nResolution min:   %.3f"  
                "\nResolution max:   %.3f"  
                "\nResolution bins:  %d"  
                "\nRequired completeness min: %.2f"
                "\nAllowed completeness tolerance:  %.2f"
                "\nAllowed rotation expansion: %.1f%%"
                "\nRequired minimum completeness increase on rotation expansion: %.2f%%"
                "\nRequired minimum redundancy: %.2f",
                m_poResoBins->dGetMinReso(), 
                m_poResoBins->dGetMaxReso(), 
                m_poResoBins->nGetNumberOfBins(),
                m_fCompletenessMin,
                m_fCompletenessTol,
                m_fRotationTol,
                m_dMinCompletenessIncreaseOnExpansion,
                m_fRedundancyMin);
    }

    if( wCtrl & DTREK_CSTRATEGY_LIST_REFLN_LIST && NULL != m_poReflnlist)
    {
        printf("\nNum reflns:      %d\n",m_poReflnlist->nGetNumReflns());
    }

    return 0;
}

//////////////////////////////////////////////////////////////////////////
// Go through the reduced reflection list and fill up the statistics table.
void Cstrategy::vFillStats(CstrategyRange* poTestRange, bool bMustBeResolutionConsistent)
{
    if( !m_poResoBins )
        return;  // safety
    
    m_poResoBins->vInitBins(DTREK_RESOBIN_SCALE_EXPERIMENTAL_COUNT);

    int                 nScan = -1;

    const   int         cn_ArraySize = 200;

    bool                abInUse[cn_ArraySize];
    double              afRotStart[cn_ArraySize];
    double              afRotEnd[cn_ArraySize];

    if( poTestRange )
    {
        // Compute InUse,RotStart and RotEnd fields for all reflections.
        if( m_anScanStart.size() > cn_ArraySize )
        {
            printf("ERROR:  Serious.  Fixed internal array dimensions exceeded.\n");
                return;
        }

        int         nScani = -1;
        double      fRotStart = 0.0;
        double      fRotEnd   = 0.0;
        for(nScan=0; nScan < m_anScanStart.size(); nScan++) 
        {
            if( bIsFixedScan(nScan) )
                continue;  // Fixed scans are always used 

            nScani = poTestRange->m_anScanNumbers.find(nScan);
            if (nScani != -1) 
            {
                // The user specified this scan. Use it.
                fRotStart = poTestRange->m_afRotStart[nScani];
                fRotEnd   = poTestRange->m_afRotEnd[nScani];
            } 
            else 
            {
                // The user did not specify this scan.
                abInUse[nScan] = false;
                continue;
            }

            abInUse[nScan] = true;
            afRotStart[nScan] = fRotStart;
            afRotEnd[nScan] = fRotEnd;
        }
    }

    CResoStrategyBin*       pBin = NULL;
    CResoStrategyBin*       pBinPrev = NULL;
    
    int       nPackPrev = -1;
    int       nMults    = 0;

    bool      bInRange = false;
    int       nNumSelect = 0;
    int       nNumUnselect = 0;

    const double c_dSTL_tolerance = 0.0001;
    Crefln*             poRefln = NULL;
    double              fSTLcu = 0.0;
    int                 nPack = 0;
    int                 nReject = 0;

    ////RB debug
    //int* pna = new int[m_poReflnlist->nGetNumReflns()];
    //
    ////m_poReflnlist->vSort(eReflnField_float_type, m_poReflnlist->m_nFI_fIntensity, 
    //m_poReflnlist->vSort2(eReflnField_float_type, m_poReflnlist->m_nFI_fIntensity, 
    //                      eReflnField_int_type, m_poReflnlist->m_nFI_nPackedHKL, pna); 
    //m_poReflnlist->nWrite("C:\\Temp\\test.ref", pna);
    int*        pReflnListSortIndex = m_poReflnlist->pnGetSortIndex();
    int         nSortedIndex = -1;
    int         nFI_f2STLcu = m_poReflnlist->nGetFieldIndex(ms_sf2STLcu); 
    for(int ii = 0; ii < m_poReflnlist->nGetNumReflns(); ii++)
    {
        nSortedIndex = pReflnListSortIndex ? pReflnListSortIndex[ii] : m_pnIndexA[ii];
        poRefln  = m_poReflnlist->poGetRefln(nSortedIndex);

        if( poTestRange )
        {
            bInRange = false;
            
            nScan = (int)poRefln->fGetIntensity();

            if( bIsFixedScan(nScan) )
            {
                bInRange = true;  // All fixed scan reflections must pass this test  
            }
            else
            {
                if( !abInUse[nScan] )
                    continue;

                bInRange  =   poRefln->fGetField(m_nFI_fRotStart) >= afRotStart[nScan] + m_fPad &&
                              poRefln->fGetField(m_nFI_fRotEnd)   <= afRotEnd[nScan]   - m_fPad;

                bInRange |=   poRefln->fGetField(m_nFI_fRotStart) >= afRotStart[nScan]  + m_fPad - 360.0 &&
                              poRefln->fGetField(m_nFI_fRotEnd)   <= afRotEnd[nScan]    - m_fPad - 360.0;

                bInRange |=   poRefln->fGetField(m_nFI_fRotStart) >= afRotStart[nScan]  + m_fPad + 360.0 &&
                              poRefln->fGetField(m_nFI_fRotEnd)   <= afRotEnd[nScan]    - m_fPad + 360.0;
            }
        }
        else
            bInRange = true;

        if( bInRange )
        {
            nNumSelect++;

            fSTLcu = poRefln->fGetField(nFI_f2STLcu);

            pBin = m_poResoBins->pGetBin(fSTLcu);

            if( NULL == pBin ) // We need to skip that reflection, but first: let's check for possible algorithm errors:
            {
                if( fabs(fSTLcu - m_poResoBins->dGetDStarCubedMax()) < c_dSTL_tolerance )
                {     // Cpredict may have predicted reflections with resolutions 
                      // equal or slightly higher than the maximum allowed value -
                      // it's not an error.
                }
                else if( bMustBeResolutionConsistent )
                {
                    printf("WARNING: Internal problem with predicted reflection resolution!\n"
                           "Resolution, Resolution Limits, Tolerance (A**(-1)): %.2f %.2f %.2f %2f\n",
                       fSTLcu, m_poResoBins->dGetDStarCubedMin(), m_poResoBins->dGetDStarCubedMax(), c_dSTL_tolerance);
                }

                continue; 
            }

            pBin->vAddRefl();

            nPack   = poRefln->nGetField(m_poReflnlist->m_nFI_nPackedHKL);
            nReject = poRefln->nGetField(m_poReflnlist->m_nFI_nNonunfFlag);

            if( 0 != nReject )
            {
                // Rejected for bad nonunf or overlap or other reason
                pBin->vAddRej();
                continue;
            }
            else
            {
                pBin->vAddUsed();
            }

            if( nPack == nPackPrev )
            {
                // Count redundant HKL
                nMults++;
                // RB: Notice that since the same reduced HKL must have the same resolution, Thad doesn't have to set nBinPrev to nBinReso
                // because that must hold true anyway.
            }
            else  // The reduced HKL has changed
            {
              //// New HKL, first deal with any old set of HKLs
              if( pBinPrev )
              {
                  if( nMults > 0 )
                  {
                      pBinPrev->vAddMult();
                  }
                  else
                  {
                      pBinPrev->vAddSingle();
                  }
              }

              // Reset flags and counters for this set of unique reflections
              nPackPrev = nPack;
              
              
              pBinPrev = pBin;

              nMults    = 0;
            }
        }
        else
        {
            nNumUnselect++;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( pBinPrev )
    {
        if( nMults > 0 )
        {
          pBinPrev->vAddMult();
        }
        else
        {
          pBinPrev->vAddSingle();
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    if( 1.0 != m_fCellVolumeFactor )
    {
        m_poResoBins->vScaleReflnCount(m_fCellVolumeFactor, DTREK_RESOBIN_SCALE_EXPERIMENTAL_COUNT);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
int Cstrategy::nCalcCompletenessSensitivity(CstrategyRange& oRange, bool bPrint)
{
    if( oRange.m_anScanNumbers.size() <= m_anFixedScans.size() )
        return 1;  // RB safety feature to make sure range is not empty
    
    int nx,ny;


    CstrategyRange oTestRange(this);
    const double fDegSensitivityStep = 4.0;
    const double fDegSensitivityMicroStep = 0.5;
    const double fMinPercentCompletePerDegree = 0.025;
    int nDegSensitivityStep = (int) (fDegSensitivityStep / fDegSensitivityMicroStep);

    double f0;
    double fMinWidth;
    double fOrigWidth;
    double fWidth;
    double a2fCompleteness[2];
    double fLastCompleteness = 0.0;
    itr<double> afRangeStart;
    itr<double> afRangeEnd;
    itr<double> afPercentComplete;
    itr<double> afRedundancy;

    if ( oRange.m_anScanNumbers.size() - m_anFixedScans.size() > 1 || 
         m_anScanStart.size() - m_anFixedScans.size() != 1 )
        return 1;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    int  nBest = 0;

   
    -afRangeStart;
    -afRangeEnd;
    -afPercentComplete;
    -afRedundancy;
    oTestRange = oRange;

    // Restrict range if requested.
    if (oTestRange.m_afRotEnd.last() - oTestRange.m_afRotStart.last()> m_afScanRotRange[0])
        oTestRange.m_afRotEnd.last() = oTestRange.m_afRotStart.last() + m_afScanRotRange[0];

    // Decrease until we are 10 degrees or less.

    fOrigWidth = oTestRange.m_afRotEnd.last() - oTestRange.m_afRotStart.last();
    fMinWidth = max(10.0,fOrigWidth - 90.0);
    while ((fWidth = oTestRange.m_afRotEnd.last() - oTestRange.m_afRotStart.last())>fMinWidth) {
        nCompleteness(oTestRange);
        afRangeStart + oTestRange.m_afRotStart.last();
        afRangeEnd + oTestRange.m_afRotEnd.last();
        afPercentComplete + oTestRange.m_fCompleteness;
        afRedundancy + oTestRange.m_fRedundancy;
        
        oTestRange.m_afRotEnd.last() -= fDegSensitivityMicroStep;
        nCompleteness(oTestRange);
        a2fCompleteness[0] = oTestRange.m_fCompleteness;
        oTestRange.m_afRotStart.last() += fDegSensitivityMicroStep;
        oTestRange.m_afRotEnd.last() += fDegSensitivityMicroStep;
        nCompleteness(oTestRange);
        oTestRange.m_afRotStart.last() -= fDegSensitivityMicroStep;
        a2fCompleteness[1] = oTestRange.m_fCompleteness;
        f0 = fDegSensitivityMicroStep;
        if ((fOrigWidth - fWidth)/(fOrigWidth - fMinWidth) > 0.666)
            f0*=4;
        else if ((fOrigWidth - fWidth)/(fOrigWidth - fMinWidth) > 0.25)
            f0*=2;
        if (a2fCompleteness[0]>a2fCompleteness[1])
            oTestRange.m_afRotStart.last() += f0;
        else
            oTestRange.m_afRotEnd.last() -= f0;
    } 

    if (!afRangeStart.size())
        return 0;

       
    // Only display this table if we have a singe scan.
    if (bPrint) {
        printf("Completeness and Redundancy vs Rotation range\n");
        printf("---------------------------------------------\n");
        printf("    Rotation     Cumul  %%Comp   Incr%%C  Redu-\n");
        printf("      range      range   cumul   /deg   nd'y\n");
        printf("---------------------------------------------\n");
    };

    double fLastPrintedCompleteness;
    double fLastPrintedWidth;
    bool   bSolutionFound;
    bool   bPrintThisLine;

    fLastPrintedCompleteness = 0.0;
    fLastPrintedWidth = 0.0;
    bSolutionFound = false;
    for (nx=afRangeStart.size()-1,ny=nDegSensitivityStep;nx>=0;nx--,ny++) {


        bPrintThisLine =(ny >= nDegSensitivityStep);
            
        if ((afPercentComplete[nx]>=m_fCompletenessMin) && 
            (afRedundancy[nx]>=m_fRedundancyMin) && 
            (!bSolutionFound))
        {
            bPrintThisLine = true;
            bSolutionFound = true;

            nBest = nx;            

        } 
        else if (afPercentComplete[nx]-fLastPrintedCompleteness>=fMinPercentCompletePerDegree*fDegSensitivityStep)\
        {
            nBest = nx;            
        }
        else
            bPrintThisLine = false;

        if( bPrintThisLine && bPrint )
        {
            ny = 0;

            printf(" %5.1lf - %5.1lf %7.1lf %6.1lf  %6.2lf %5.1lf\n",
                afRangeStart[nx],
                afRangeEnd[nx],
                afRangeEnd[nx] - afRangeStart[nx],
                afPercentComplete[nx],
                (afPercentComplete[nx] - fLastPrintedCompleteness)/(afRangeEnd[nx] - afRangeStart[nx] - fLastPrintedWidth),
                afRedundancy[nx]
                );       

            fLastPrintedWidth = afRangeEnd[nx] - afRangeStart[nx];
            fLastPrintedCompleteness = afPercentComplete[nx];

        };
    };
    nx = max(0,nBest);
    if (bPrint) {
        printf("---------------------------------------------\n");
        printf(" %5.1lf - %5.1lf %7.1lf %6.1lf  %6s %5.1lf\n\n",
            afRangeStart[nx],
            afRangeEnd[nx],
            afRangeEnd[nx] - afRangeStart[nx],
            afPercentComplete[nx],
            "----",
            afRedundancy[nx]
            );
        printf("(Incr%%C/deg is the incremental %%completeness"
               " per additional degree of rotation.)\n\n");
    };

    if (bPrint) {
        printf("Completeness vs Rotation range (selected scan)\n");
        printf("---------------------------------\n");
        printf("    Rotation     Cumul  %%Comp\n");
        printf("      range      range   cumul\n");
        printf("---------------------------------\n");

        for (nx=0;nx<10;nx++) {
            oTestRange = oRange;
            oTestRange.m_afRotEnd.last() = oTestRange.m_afRotStart.last() + ((nx+1)*(oTestRange.m_afRotEnd.last() - oTestRange.m_afRotStart.last()))/10;
            nCompleteness(oTestRange);
            printf(" %5.1lf - %5.1lf %7.1lf %6.1lf\n",
            oTestRange.m_afRotStart.last(),
            oTestRange.m_afRotEnd.last(),
            oTestRange.m_afRotEnd.last() - oTestRange.m_afRotStart.last(),
            oTestRange.m_fCompleteness);            
        };
        printf("---------------------------------\n");
    };


    return 0;
}

int Cstrategy::nListCompleteness(CstrategyRange* poRange, DTREK_WORD wCtrl)
{
   char*         pcLine = "--------------------------------------------------------------------------------\n";

    // The solution should already be loaded. Just recompute statistics here to be sure.
    vFillStats(poRange);

    // List completeness statistics
    printf("\n\nExpected Completeness vs Resolution\n");
    printf(pcLine);
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
    "Resolution  ", "  Calc", "Num", " Num", " Num", "   Num", "Num",
    " Avg", "%Comp", "%Comp");
    printf("%14s %7s %7s %5s %7s %7s %6s %6s %6s %6s\n",
    "   range    ", "unique", "reflns", "rejs", "mults", "single",
    "unique", "mult", "shell", "cumul");

    printf(pcLine);

    CResoStrategyBin*       pBin = NULL;
    CResoStrategyBin        oCumulBin;

    for (int ii = 0; ii < m_poResoBins->nGetNumberOfBins(); ii++)
    {
        pBin = m_poResoBins->pGetBin(ii);

        if( !pBin )
          continue;

        oCumulBin.vAddStats(*(pBin->pGetStats()));

        if( 0 != pBin->nGetNumUniqueCalc()   && 
            0 != pBin->nGetNumUsed()         &&
            0 != pBin->nGetNumUniqueFound() )
        {
         printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d %6.2f %6.1f ",
                 min(pBin->dGetLowReso(), 999.99), 
                 pBin->dGetHighReso(),
      
                 pBin->nGetNumUniqueCalc(),
                 pBin->nGetNumRefls(),
                 pBin->nGetNumRejs(),
                 pBin->nGetNumMults(),
                 pBin->nGetNumSingles(),
                 pBin->nGetNumUniqueFound(),
                 (double) (pBin->nGetNumUsed()) / (double) pBin->nGetNumUniqueFound(),
                 min(100., 100. * (double) pBin->nGetNumUniqueFound() / (double) pBin->nGetNumUniqueCalc()));
        }
        else if( 0 < pBin->nGetNumUniqueCalc() )
        {
          printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d %6s %6.1f ",
                  min(pBin->dGetLowReso(), 999.99), 
                  pBin->dGetHighReso(),
                  pBin->nGetNumUniqueCalc(),
                  pBin->nGetNumRefls(), 
                  pBin->nGetNumRejs(),
                  pBin->nGetNumMults(), 
                  pBin->nGetNumSingles(),
                  pBin->nGetNumUniqueFound(),
                  "---",
                  min(100., 100. * (double) pBin->nGetNumUniqueFound() / (double) pBin->nGetNumUniqueCalc()));
        }
        else
        {
          printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d %6s %6s ",
                  min(pBin->dGetLowReso(), 999.99), 
                  pBin->dGetHighReso(),
                  pBin->nGetNumUniqueCalc(),
                  pBin->nGetNumRefls(), 
                  pBin->nGetNumRejs(),
                  pBin->nGetNumMults(),
                  pBin->nGetNumSingles(),
                  pBin->nGetNumUniqueFound(),
                  "---", 
                  "---");
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if( 0 != oCumulBin.nGetNumUniqueCalc() )
          printf("%6.1f\n", 
            min(100., 100. * (double) oCumulBin.nGetNumUniqueFound() / (double)oCumulBin.nGetNumUniqueCalc()));
        else
          printf("%6s\n", "---");

        //fLow = fLow + *pfSlope;
    }

    printf(pcLine);

    printf("%6.2f - %5.2f %7d %7d %5d %7d %7d %6d ",
           min(m_poResoBins->dGetMinReso(), 999.99), m_poResoBins->dGetMaxReso(),
           oCumulBin.nGetNumUniqueCalc(),
               oCumulBin.nGetNumRefls(), 
           oCumulBin.nGetNumRejs(),
               oCumulBin.nGetNumMults(), 
           oCumulBin.nGetNumSingles(),
           oCumulBin.nGetNumUniqueFound());
  
    if( 0.0 != oCumulBin.nGetNumUniqueCalc() )
    {
        printf("%6.2f ", (double) oCumulBin.nGetNumUsed() / (double) max(1, oCumulBin.nGetNumUniqueFound() ) );
        
        printf("%6.1f ",  min(100., 100.* (double) oCumulBin.nGetNumUniqueFound() / (double) oCumulBin.nGetNumUniqueCalc()));
        printf("%6.1f\n", min(100., 100.* (double) oCumulBin.nGetNumUniqueFound() / (double) oCumulBin.nGetNumUniqueCalc()));
    }
    else
    {
        printf("%6s %6s %6s\n", "---", "---", "---");
    }

    if (1.0 != m_fCellVolumeFactor)
    {
        printf("\n*** Note: cell length factor was %4.2f, so the results above "
        "are approximate. ***\n", m_fCellLengthFactor);
    }

    double          fLow =  min(100.0f, 100.0f* (double) oCumulBin.nGetNumUniqueFound() / oCumulBin.nGetNumUniqueCalc());
    
    if( 100.0 < m_fCompletenessMin )
    {
        printf("\n*** It is impossible to get more than 100.0%% completeness! ***");
    }
    
    if( m_fCompletenessMin > fLow + m_fCompletenessTol )
    {
        if( wCtrl & DTREK_STRATEGY_LC_MSG_LOW_COMPLETENESS )
        {
            // Did not achieve requested completeness, so inform user
            printf("\n*** The requested completeness of %.3lf%% was not achieved! ***\n", m_fCompletenessMin);
        }

        if( wCtrl & DTREK_STRATEGY_LC_MSG_USER_ACCEPTED )
        {
            printf("\n*** Instead, the user accepted solution with best possible completeness of %.3lf%% and redundancy of %.3lf! ***\n",
            m_oBestRange.m_fCompleteness, m_oBestRange.m_fRedundancy);
        }
    
        if( !(wCtrl & DTREK_STRATEGY_LC_MSG_USER_ACCEPTED) && (wCtrl & DTREK_STRATEGY_LC_MSG_LOW_COMPLETENESS) )
        {
            fLow = (double) ( (int) (fLow - 0.5));
            printf("*** You may wish to try the following:\n"
            "    * Decrease the requested completeness to %.1lf%%.\n"
            "    * Adjust the Resolution range.\n",
            fLow);
        }
    }
  
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cstrategy::nClearRange(CstrategyRange& oRange)
{
    -oRange.m_anScanNumbers;
    -oRange.m_afRotStart;
    -oRange.m_afRotEnd;
    
    for(int nx=0; nx < m_anFixedScans.size(); nx++)
    {
        oRange.m_anScanNumbers + m_anFixedScans[nx];
        oRange.m_afRotStart    + m_afScanRotMin[m_anFixedScans[nx]];
        oRange.m_afRotEnd      + m_afScanRotMax[m_anFixedScans[nx]];
    }
    
    return 0;
}

// This function will first call nCompleteness for each scan in the range and save the results in m_poScanInfo.
// Then it will call nCompleteness for the whole range.
int Cstrategy::nCompletenessMultiScan(CstrategyRange& oTestRange)
{
    int     nNumberOfScans = oTestRange.m_anScanNumbers.size();
    
    if( nNumberOfScans <= 0 )
        return 1; // safety
    
    if( nNumberOfScans == 1 )
    {
        int     nRet = nCompleteness(oTestRange);
        nFromRangeToStrategyScansArray(oTestRange, DTREK_DTSTRATEGY_FRTSSA_UPDATE_ALL_RESULTS); 
        return nRet;
    }

    CstrategyRange      oSingleScanRange(this);
    
    nClearRange(oSingleScanRange);

    int     nFixedScanOffset = m_anFixedScans.size();

    int     nScan = -1;
    for(int ii=0; ii < nNumberOfScans; ii++)
    {
        nScan = oTestRange.m_anScanNumbers[ii];
        
        if( -1 != m_anFixedScans.find(nScan) )
           continue;

        oSingleScanRange.m_anScanNumbers[nFixedScanOffset] = oTestRange.m_anScanNumbers[ii]; 
        oSingleScanRange.m_afRotStart[nFixedScanOffset]    = oTestRange.m_afRotStart[ii];   
        oSingleScanRange.m_afRotEnd[nFixedScanOffset]      = oTestRange.m_afRotEnd[ii]; 
        
        nCompleteness(oSingleScanRange);

        nFromRangeToStrategyScansArray(oSingleScanRange, DTREK_DTSTRATEGY_FRTSSA_UPDATE_ALL_RESULTS); 
    }
    
    return nCompleteness(oTestRange);
}

int Cstrategy::nCompleteness(CstrategyRange& oTestRange)
{
  //// Calculate completeness
  Cstat                     oStat;
  
  vFillStats(&oTestRange);

  oStat.vClearRaw();
  
  CResoStrategyBin*     pBin = NULL;
  CResoStrategyBin      oCumulBin;

  for (int ii = 0; ii < m_poResoBins->nGetNumberOfBins(); ii++)
  {
      pBin = m_poResoBins->pGetBin(ii);

      oCumulBin.vAddStats(*(pBin->pGetStats()));
      
      oStat.vAddRaw( (double)pBin->nGetNumUsed() / max(1.0, (double)pBin->nGetNumUniqueFound()) );
  }     

  if( 0 == oCumulBin.nGetNumUniqueFound() )
  {
    oTestRange.m_fRedundancy = 0.0;
    oTestRange.m_fRedundancyDev = 0.0;
  }
  else
  {
      oTestRange.m_fRedundancy =   (double) (oCumulBin.nGetNumUsed()) / (double) oCumulBin.nGetNumUniqueFound();
      oTestRange.m_fRedundancyDev = oStat.fStandardDeviationRaw();
  }

  //////////////////////////////////////////////////////////////////////////
  // RB: not sure why Thad had this. Let's just set it to zero in this case.
  //if( oTestRange.m_fRedundancyDev < 0 )
  //  oTestRange.m_fRedundancyDev = oStat.fStandardDeviationRaw();
  if( oTestRange.m_fRedundancyDev < 0 )
    oTestRange.m_fRedundancyDev = 0.0;
  ///////////////////////////////////////////////////////////////////////////

  if( 0 == oCumulBin.nGetNumUniqueCalc() )
      oTestRange.m_fCompleteness = 0.0;
  else
  {
      oTestRange.m_fCompleteness = min(100.0, 100.0 * ((double)oCumulBin.nGetNumUniqueFound() / (double)oCumulBin.nGetNumUniqueCalc()));
  }
  
  return 0;
}

int Cstrategy::nExpand(Creflnlist& oList)
{
    (void) oList.nExpandGetField(oList.ms_ssBatch);

    m_nFI_nPackedHKL  = oList.nExpandGetField(oList.ms_snPackedHKL);
    
    m_nFI_nPackedHKL2 = oList.nExpandGetField("nPackedHKL2");
    
    m_nFI_f2STLcu = oList.nExpandGetField(ms_sf2STLcu);
    
    m_nFI_fRotStart = oList.nExpandGetField(oList.ms_sfCalcRotStart); 
    
    m_nFI_fRotEnd = oList.nExpandGetField(oList.ms_sfCalcRotEnd); 
    
    return 0;
}


int Cstrategy::nSetup(Cimage_header& oHeader, Creflnlist& oReflnlist)
{
    int       i,j;
    itr<int>  anTemp;
    int       nStat;

    // Check if m_poReflnlist and m_poCrystal are valid
    if (nExpand(oReflnlist))
      return 1;
    m_poReflnlist = & oReflnlist;

    nStat = 0;
    if (   (NULL == m_poReflnlist)
      || (NULL == m_poCrystal) )
    {
      nStat = -1;
    }
    else if (   !m_poReflnlist->bIsAvailable()
       || !m_poCrystal->bIsAvailable())
    {
      nStat = -1;
    }
    else if (   (0 >= m_poReflnlist->m_nFI_fCalcRotEnd)
       || (0 >= m_poReflnlist->m_nFI_fCalcRotStart)
       || (0 >= m_poReflnlist->m_nFI_fResolution) )
    {
      nStat = -2;
    }

    if (0 != nStat)
    {
      printf("ERROR in Cstrategy::nSetup: invalid crystal or reflnlist!\n");
      return (nStat);
    }

    //nList(DTREK_CSTRATEGY_LIST_CRYSTAL);

    if( !m_poCrystal->bLinearScaleToVolumeScale(m_fCellLengthFactor, m_fCellVolumeFactor) )
    {
        double      fVolume = m_poCrystal->fCalcVolume();

        if ( (50.0 * 50.0 * 25.0) < fVolume )
        {
            printf("NOTE: Unit cell volume is %.1lf\n",fVolume);
            printf( "      you may wish to use a cell length factor of about %.3lf next time.\n",
                    pow((50.0 * 50.0 * 25.0) / fVolume, 0.33333));
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Figure out the sin(theta)/lambda limits for the reflection list and SAVE it for every refln
    double        dStarCubedMin     = 0.0;
    double        dStarCubedMax     = 0.0;
    vReflnListCalcSetDStarCubed(dStarCubedMin, dStarCubedMax);
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Create (if not yet created) and set a reso bins object. 
    if( !bSetResoBinsObject(dStarCubedMin, dStarCubedMax) )
    {
        printf("ERROR: Cannot set resolution object in Cstrategy\n");
        return 1;
    }

    nStat = m_poReflnlist->nSelect(m_nFI_f2STLcu, "<", (float) m_poResoBins->dGetDStarCubedMin() , FALSE);
    nStat = m_poReflnlist->nSelect(m_nFI_f2STLcu, ">", (float) m_poResoBins->dGetDStarCubedMax()   , FALSE);

    nStat = m_poReflnlist->nSelect((Cstring)'+' + m_poReflnlist->ms_snDetNum + "<0");

    // Delete the deselected reflections
    m_poReflnlist->nDelete((const char*)NULL, (const bool)FALSE);

    if( !m_poReflnlist->nGetNumReflns() )
    {
      printf("ERROR:  No reflections predicted!\n");
      return 1;
    }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Make sure that every scan has a representative.  If not, then we have to add some dummy reflections.
  j = 0;
  
  for (i =0; i < m_poReflnlist->nGetNumReflns(); i++)
  {
      j = max(j , (int)(*m_poReflnlist)[i].fGetIntensity());
  }
  
  -anTemp;
  
  for (i = 0; i <= j; i++)
      anTemp[i] = 0;
  
  for (i =0; i< m_poReflnlist->nGetNumReflns(); i++)
  {
      anTemp[(int) (*m_poReflnlist)[i].fGetIntensity()] = 1;
  }

  for (i = 0; i <= j; i++)
  {
      if (anTemp[i]!=1)
      {
          // Add a dummy reflection.
          Crefln oRefln((*m_poReflnlist)[0]);
          oRefln.vSetIntensity((float) i);
          m_poReflnlist->nInsert(&oRefln);
      }
  }

  // Sort and reduce reflnlist to asymmetric unit

  printf("Reducing reflnlist to unique HKL and sorting ...\n");
  
  /////////////////////////////////////////////////////////////
  Crefln*       poRefln = NULL;

  for (int nx=0;nx<m_poReflnlist->nGetNumReflns();nx++)
  {
      poRefln = m_poReflnlist->poGetRefln(nx);
      poRefln->vSetField(m_nFI_nPackedHKL2,poRefln->nPackHKL(poRefln->pnGetHKL()));
  }
  
  if (m_pnIndex2) 
      delete[] m_pnIndex2;
  
  m_pnIndex2 = new int[m_poReflnlist->nGetNumReflns()];
  
  m_poReflnlist->vSort2(eReflnField_float_type, 0, eReflnField_int_type, m_nFI_nPackedHKL2, m_pnIndex2);

  printf("Reducing reflnlist to asymmetric unit and sorting ...\n");
  printf("Pass 1...\n");
  nStat = m_poReflnlist->nReduce(*m_poCrystal, m_bAnom);

  if (m_pnIndex)
      delete[] m_pnIndex;

  m_pnIndex = new int[m_poReflnlist->nGetNumReflns()];

  // The intensity field contains the scan number.  Sort on packed HKL within each scan.
  m_poReflnlist->vSort2(eReflnField_float_type,0,eReflnField_int_type, m_poReflnlist->m_nFI_nPackedHKL, m_pnIndex);
  
  printf("Pass 2...\n");
  if (m_pnIndexA)
      delete[] m_pnIndexA;

  m_pnIndexA = new int[m_poReflnlist->nGetNumReflns()];
  m_poReflnlist->vSort2(eReflnField_int_type, m_nFI_nPackedHKL, eReflnField_float_type, 0, m_pnIndexA);

  // For each scan, load the scan start and the scan end.
  -m_anScanStart;
  -m_anScanEnd;

  m_anScanStart + 0;
  
  for (i =0; i< m_poReflnlist->nGetNumReflns(); i++) 
  {
      bool  bIsScanEndReached = i == m_poReflnlist->nGetNumReflns() - 1 || 
                                (*m_poReflnlist)[m_pnIndex[i+1]].fGetIntensity() != (*m_poReflnlist)[m_pnIndex[i]].fGetIntensity();
      
      if( bIsScanEndReached )
      {
          m_anScanEnd + i;

          if( i != m_poReflnlist->nGetNumReflns() - 1 )
              m_anScanStart + (i+1);
      }
  }
  
  
  // Next, load the range/min/max arrays.  This information is obtained from the Crotation objects in the header, which have the appropriate fields.
  -m_afScanRotMin;
  -m_afScanRotMax;
  -m_afScanRotRange;
  
  for (i = 0; i < m_anScanStart.size(); i++)
    {
      if( -1 == m_anFixedScans.find(i) ) 
      {
          Crotation         oRotation(oHeader, sBuildStrategyScanPrefix(D_K_StrategyTestPrefix, D_K_ScanPrefix, i));
          
          if (!oRotation.bIsAvailable())
          {
              printf("ERROR:  Not all predicted scans were available!\n");
              return 1;
          }
      
          m_afScanRotMin + ((double) oRotation.fGetRotMin());
          m_afScanRotMax + ((double) oRotation.fGetRotMax());
          m_afScanRotRange + ((double) oRotation.fGetRotRange());
      } 
      else 
      {
          // We should use this scan.  Find the max and min ranges of it to fill in
          // these arrays.
          double    fRotMax = -1000;
          double    fRotMin = 1000.0;
          
          for(j = m_anScanStart[i]; j <= m_anScanEnd[i]; j++)
          {
              // Yes, this looks incorrect, but actually it's the right way to do it.  
              // dtpredict predicts reflections which are overlap the scan, so we can really only trust
              // the lower calculated range as an estimate for the max scan range, and
              // the upper calcllated range as an estimate for the min scan range.
              fRotMax = max(fRotMax, (double) (*m_poReflnlist)[j].fGetField(m_poReflnlist->m_nFI_fCalcRotStart));
              fRotMin = min(fRotMin, (double) (*m_poReflnlist)[j].fGetField(m_poReflnlist->m_nFI_fCalcRotEnd));
          }

          // Round the values.  They are probably a bit off.
          fRotMin = floor(0.5 + fRotMin);
          fRotMax = floor(0.5 + fRotMax);
          m_afScanRotMin + fRotMin;
          m_afScanRotMax +  fRotMax;
          m_afScanRotRange + (fRotMax - fRotMin);
      }
  }
  
  if (0 != nStat)
    {
      printf("ERROR in Cstrategy::nSetup: failure to reduce dataset!\n");
      nStat = -1;
#ifdef SSI_PC
      return (nStat);
#else
      exit (nStat);
#endif
    }
  printf(" ... done.\n");

  return 0;
}

void Cstrategy::vRememberScanGonioSettings(Cgoniometer* poCrysGonio,Crotation* poRotation,int nScan)
{
    int nx;
    // Save goniometer settings.
    m_asScanGonioAxis.push_back(poRotation->sGetName());
    for (nx = 0; nx < poCrysGonio->m_nNumRotValues;nx++)
        m_afScanGonioAngles + ((double) poCrysGonio->m_pfDatumValue[nx]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Delete predicted reflections that did not make it into the best scan range. 
bool Cstrategy::bLeaveBestRangeReflectionsOnly(bool bLeaveFixedScans)
{
    int               nScanNumber = -1;
    float             fRotStart   =  0.0f;
    float             fRotEnd     =  0.0f;

    for(int nx=0; nx < m_poReflnlist->nGetNumReflns(); nx++)
    {
        (*m_poReflnlist)[nx].vSetSigmaI(0.0);

        nScanNumber = (int)((*m_poReflnlist)[nx].fGetIntensity());

        if( !bLeaveFixedScans && -1 != m_anFixedScans.find(nScanNumber) )
            continue;

        fRotStart = (*m_poReflnlist)[nx].fGetField(m_poReflnlist->m_nFI_fCalcRotStart);
        fRotEnd   = (*m_poReflnlist)[nx].fGetField(m_poReflnlist->m_nFI_fCalcRotEnd);

        for(int ny = 0; ny < m_oBestRange.m_anScanNumbers.size(); ny++)
        {          
            if( nScanNumber == m_oBestRange.m_anScanNumbers[ny] )
            {
                if ( fRotStart >= m_oBestRange.m_afRotStart[ny] && fRotEnd <= m_oBestRange.m_afRotEnd[ny] )
                {
                    // This reflection is good.  Keep it.
                    (*m_poReflnlist)[nx].vSetSigmaI(1.0);

                    break;
                }
            }
        }
    }

    m_poReflnlist->nSelect(m_poReflnlist->m_nFI_fSigmaI, "==", (const float) 0.0, FALSE);

    // Delete the deselected reflections
    m_poReflnlist->nDelete((const char *)NULL, (const bool)FALSE);
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The purpose of this function is as follows. Originally, the refln list is generated with the original phi value,
// However, when we change phi, we need to adjust rotation values as well.
void Cstrategy::vShiftReflnRotValueInScan(int iScan, double dRotShift)
{
    int               nScanNumber = -1;
    float             fRotStart   =  0.0f;
    float             fRotEnd     =  0.0f;

    for(int nx=0; nx < m_poReflnlist->nGetNumReflns(); nx++)
    {
        nScanNumber = (int)((*m_poReflnlist)[nx].fGetIntensity());

        if( -1 != m_anFixedScans.find(nScanNumber) )
            continue;  // we don't want to touch the fixed scans

        fRotStart = (*m_poReflnlist)[nx].fGetField(m_poReflnlist->m_nFI_fCalcRotStart);
        fRotEnd   = (*m_poReflnlist)[nx].fGetField(m_poReflnlist->m_nFI_fCalcRotEnd);

        fRotStart += (float)dRotShift;
        fRotEnd   += (float)dRotShift;

        (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_fCalcRotStart, fRotStart);
        (*m_poReflnlist)[nx].vSetField(m_poReflnlist->m_nFI_fCalcRotEnd  , fRotEnd);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Cstrategy::nWrite(const Cstring& sName)
{
    Cstring     sTemp(sName);
    
    if("" == sTemp)
        sTemp = sDtrekGetPrefix() + "dtstrategy.ref";

    return m_poReflnlist->nWrite(sTemp);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
////////////// Hardcore strategy routines begin here //////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



bool Cstrategy::bScanUsed(int nScan, int nOrder)
{
    if( -1 != m_anExcludeScans.find(nScan) )
        return false;
    
    if( -1 != m_anFixedScans.find(nScan) )
        return false;

    bool        bRet = false;
    
    // RB 7/20/05  If I understand Thad correctly, ANY scan number will fall between m_a2nFirstScan[0] && m_a2nFirstScan[1],
    // and NO scan number will fall between m_a2nSecondScan[0] && m_a2nSecondScan[1],
    // unless the user used the -choose option
    if( nOrder == 0 )
    {
        bRet =  nScan >= m_a2nFirstScan[0]  && nScan <= m_a2nFirstScan[1] ||
                nScan >= m_a2nSecondScan[0] && nScan <= m_a2nSecondScan[1];
    }
    else if( nOrder == 1)
    {
        bRet =  nScan >= m_a2nFirstScan[0]  && nScan <= m_a2nFirstScan[1];
    }
    else if( nOrder == 2)
    {
        bRet =  nScan >= m_a2nSecondScan[0] && nScan <= m_a2nSecondScan[1];
    }

    return bRet;
}
////////////////////////////////////////////////////////////////////////////
bool Cstrategy::bIsFixedScan(int nScan)
{
    bool    bFound = false;
    
    for(int ii=0; ii < m_anFixedScans.size(); ii++)
    {
        if( nScan == m_anFixedScans[ii] )
        {
            bFound = true;
            break;
        }
    }
    
    return bFound;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RB If I understand Thad's code correctly, the purpose of this function is to get the best rotation
// range for oTestRange. Technically, the range object may start with fixed scans (i.e. scan numbers, scan begins and scan ends arrays
// will start from the fixed scans). That is why we need parameter nIndexToChange. It tells us where in the arrays the "variable"
// scan parameters start. We are not supposed to touch the fixed scans at all, we are not supposed to change the scan numbers either, but we
// can change the rotation start and end for the scan referenced in the arrays by nIndexToChange.
int Cstrategy::nCalcScanSingle(CstrategyRange& oTestRange, double fMaxSearchStep, int nIndexToChange)
{
    int         nUsedScan = oTestRange.m_anScanNumbers[nIndexToChange]; // this is the scan for which we are going to find the best rotation start and end.
    
    double      fAbsRotStart = oTestRange.m_afRotStart[nIndexToChange];
    double      fAbsRotEnd   = oTestRange.m_afRotEnd[nIndexToChange];
    double      fAbsRotWidth = min(360.0, m_afScanRotRange[nUsedScan]);
    
    nCompleteness(oTestRange);
    
    bool        bScanSolvable = oTestRange.m_fCompleteness >= m_fCompletenessMin - m_fCompletenessTol &&
                                oTestRange.m_fRedundancy   >= m_fRedundancyMin;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // If we have a range of 0 to 360, then we assume we can do a 'wraparound'
    bool        bWrapAround =  fAbsRotStart == 0.0 && fAbsRotEnd == 360.0;
    
    bool        bSolutionFound = false;
    double      fWidthStep = 0.0;
    double      fTestRotStart = 0.0;
    bool        bSaturated = false;

    double      a2fLastSolution[2];

    const  double       fWidthMultiplierTol = 0.001;              // How close does the width multiplier need to be to 1.0
    const  double       fInitWidthMultiplier = 0.8;               // Initial width multiplier (implies we are decreasing the width).

    
    double      fWidthMultiplier = fInitWidthMultiplier;
    double      fWidth = min(fAbsRotWidth, ceil((fAbsRotEnd - fAbsRotStart) * fWidthMultiplier));
    
    int         nTotalPasses = 0;
    
    if( !bScanSolvable )
        goto GET_BEST_WHEN_NO_SOLUTIONS;
    
    do {
        // Test to make sure that we arn't maxing out.
        if( fWidth > fAbsRotWidth ) 
        {
            fWidth = fAbsRotWidth;
            bSaturated = TRUE;
        } 
        else
            bSaturated = FALSE;
        
        fWidthStep = max(fMaxSearchStep,floor(fabs(1.0 - fWidthMultiplier)*fWidth));
        
        bSolutionFound = FALSE;
        
        for(fTestRotStart = fAbsRotStart; (!bSolutionFound) && ((fTestRotStart + fWidth<fAbsRotEnd) || 
                                          ((bWrapAround) && (fTestRotStart<fAbsRotEnd))); fTestRotStart+= fWidthStep)
        {
            ////////////////////////////////////////////////////////////////////////////////////////
            // RB: This is a "hack" into Thad's code, so we don't have to change a lot
            if(m_bFixRotStart)
               fTestRotStart = fAbsRotStart;
            //////////////////////////////////////

            oTestRange.m_afRotStart[nIndexToChange] = fTestRotStart;
            oTestRange.m_afRotEnd[nIndexToChange] = min((fTestRotStart + min(fAbsRotWidth,fWidth)),fAbsRotEnd);
            ////////////////////////////////////////////////////////////////////////////////////////
            
            nCompleteness(oTestRange);

            // Have we found a solution?
            if( oTestRange.m_fCompleteness >= m_fCompletenessMin - m_fCompletenessTol && 
                oTestRange.m_fRedundancy   >= m_fRedundancyMin )
            {
                bSolutionFound = TRUE;
                
                // Save this solution
                a2fLastSolution[0] = oTestRange.m_afRotStart[nIndexToChange];
                a2fLastSolution[1] = oTestRange.m_afRotEnd[nIndexToChange];
            }
        
            if( m_bFixRotStart )
                break;    // we don't need to loop in this case, because  the start cannot be moved
        }
        
        if (bSolutionFound)
        {
            // We found a solution.  
            if( fWidthMultiplier > 1.0 )
            {
                // We were increasing the step.  We can now decrease.
                fWidthMultiplier = 1.0/fWidthMultiplier;
                // But now the width multiplier should be closer to 1.0
                fWidthMultiplier = (1.0 + fWidthMultiplier) / 2.0;
            } 
            else 
            {
                // We were decreasing the step.  Continue decreasing at the same rate.
            }
        }
        else 
        {
            if( fWidthMultiplier < 1.0 ) 
            {
                // We were decreasing the step.  We should not increase it.
                fWidthMultiplier = 1.0 / fWidthMultiplier;
                // But now the width multiplier should be closer to 1.0
                fWidthMultiplier = (1.0 + fWidthMultiplier) / 2.0;
            } 
            else 
            {
                // We were increasing the step.  Continue to do so at the same rate.
            }
        }
        
        if (fWidth != ceil(fWidth*fWidthMultiplier))
            fWidth = ceil(fWidth*fWidthMultiplier);
        else
            fWidthMultiplier = 1.0;
        
        nTotalPasses++;
    
    } while ((fabs(fWidthMultiplier - 1.0) > fWidthMultiplierTol) && (!bSaturated));
    
    // If we don't need to adjust the completeness, just break out.


    // RB 2/16/05 We need to find the best possible solution even if it does not satisfy the user's desire for completeness/redundancy.
    // Therefore, instead of "&&" we have to use "||"
    //if( !bSolutionFound && !bScanSolvable ) 
GET_BEST_WHEN_NO_SOLUTIONS:;
    if( !bSolutionFound || !bScanSolvable ) 
    {
        CstrategyRange      oBestTestRange(this);
        nClearRange(oBestTestRange);
        
        oBestTestRange.m_anScanNumbers + nUsedScan;
        
        oBestTestRange.m_afRotStart[nIndexToChange] = fAbsRotStart;
        oBestTestRange.m_afRotEnd[nIndexToChange]   = fAbsRotStart + fAbsRotWidth < fAbsRotEnd ? fAbsRotStart + fAbsRotWidth : fAbsRotEnd; // let's start with the most left adjusted range
        
        oBestTestRange.m_bIsAvailable = true;

        fWidthStep = MAX_SEARCH_STEP_SINGLE_SCAN;
        fWidth = fAbsRotWidth;
        
        bool    bIsBetterBecauseOfCompleteness = false;  // a flag to tell us if a scan is better because of its completeness
        
        for( fTestRotStart = fAbsRotStart; fTestRotStart + fWidth < fAbsRotEnd || (bWrapAround && fTestRotStart < fAbsRotEnd); fTestRotStart += fWidthStep)
        {
            ////////////////////////////////////////////////////////////////////////////////////////
            // RB: This is a "hack" into Thad's code, so we don't have to change a lot
            if( m_bFixRotStart )
               fTestRotStart = fAbsRotStart;
            //////////////////////////////////////

            oTestRange.m_afRotStart[nIndexToChange] = fTestRotStart;
            oTestRange.m_afRotEnd[nIndexToChange] = min((fTestRotStart + min(fAbsRotWidth,fWidth)),fAbsRotEnd);
            
            nCompleteness(oTestRange);
            
            //// First consider completeness. If completeness is the same, go for the higher redundancy.
            if( !bRangesCompareConsiderTargets(oBestTestRange,          // false means oBesttestRange is worse than oTestrange
                                               oTestRange, 
                                               m_fCompletenessMin - m_fCompletenessTol, 
                                               m_fRedundancyMin,
                                               bIsBetterBecauseOfCompleteness) ) // is the "better" scan better because of completeness?
            {
                oBestTestRange = oTestRange;
            }
            // If there is no solution and a test range loses to the currently best range, 
            // because the the test range's completeness is lower,
            // we can stop looping right there.
            else if( bIsBetterBecauseOfCompleteness  && !bScanSolvable )
                break;

            if( m_bFixRotStart )
                break;    // we don't need to loop in this case, because  the start cannot be moved
        }

        oTestRange = oBestTestRange;
    }
    
    nCompleteness(oTestRange);
    
    if( oTestRange.m_fCompleteness >= m_fCompletenessMin - m_fCompletenessTol && oTestRange.m_fRedundancy >= m_fRedundancyMin )
        return 0;
    else
        return 1;
}

int Cstrategy::nCalcScanMultiple(CstrategyRange& oUserRange) 
{
    int         nx;
    itr<int>    anTemp;

    time_t      oTimeStart;
    time_t      oTimeNow;
    
    const int nMaxSaveRanges = 200;
    int         nSaveRanges;
    
    CstrategyRange  oBestRangesSave[nMaxSaveRanges];
    for(int ii=0; ii < nMaxSaveRanges; ii++)
    {
        oBestRangesSave[ii].vSetStrategyPtr(this);    
    }
    
    CstrategyRange  oTestRange(this);
    CstrategyRange  oBestRangeInner(this);
    CstrategyRange  oBestRangeOuter(this);
    
    itr<int>        anLastWorkingWidth;

    int nWidth;
    int nVariable;
    int nVariables;
    int nInnerLoop;
    int nOuterLoop;
    int nLastOuterLoopWorked;

    bool bConverged;
    bool bSolutionFound;
    bool bTimeOut;
    double fWidth;
    double fWidthStep;

    double fTestRotStart;
    double fAbsRotStart;
    double fAbsRotEnd;
    double fAbsRotWidth;
    bool   bWrapAround;

    int anWidth[]       = { 30,15,8,3,-1};
    int anWidthStep[]   = { 10,10,5,5,-1};
  
    nVariables = oUserRange.m_anScanNumbers.size();
    if( nVariables <= 0 )
        return 1; //safety
    
    oBestRangeOuter = oUserRange;
    nOuterLoop = 0;

    time(&oTimeStart);

    do {

        for (nx=0;nx<nVariables;nx++)
            anLastWorkingWidth[nx] = 0;

        bConverged = false;
        oBestRangeInner = oUserRange;
        nInnerLoop = 0;

        while (!bConverged) 
        {
            
            nSaveRanges = 0;
            
            // Choose a variable to tweak pseudo-randomly.
            do 
            {
                nVariable = rand() % nVariables;
            } while (anWidth[anLastWorkingWidth[nVariable]]==-1);
            
            for (nWidth=anLastWorkingWidth[nVariable];anWidth[nWidth]!=-1;nWidth++) 
            {
                oTestRange = oBestRangeInner;
                bSolutionFound = false;
                fAbsRotStart = oUserRange.m_afRotStart[nVariable];
                fAbsRotEnd = oUserRange.m_afRotEnd[nVariable];
                fAbsRotWidth= m_afScanRotRange[oTestRange.m_anScanNumbers[nVariable]];
                fWidth = (double) min(fAbsRotWidth,oTestRange.m_afRotEnd[nVariable] - oTestRange.m_afRotStart[nVariable] - anWidth[nVariable]);
                if (fWidth>20.0)
                {
                    fWidthStep = max((double) anWidthStep[nVariable],fWidth/4.0);
                    
                    // If we have a range of 0 to 360, then we assume we can do a 'wraparound'
                    bWrapAround =  ((fAbsRotStart == 0.0) && (fAbsRotEnd == 360.0));
                    
                    
                    for (fTestRotStart = fAbsRotStart; ((fTestRotStart + fWidth<fAbsRotEnd) || 
                                                       ((bWrapAround) && (fTestRotStart<fAbsRotEnd)));
                                                        fTestRotStart+= fWidthStep)
                    {
                        ////////////////////////////////////////////////////////////////////////////////////////
                        // RB: This is a "hack" into Thad's code, so we don't have to change a lot
                         if(m_bFixRotStart)
                            fTestRotStart = fAbsRotStart;
                        ////////////////////////////////////
                        
                        oTestRange.m_afRotStart[nVariable] = fTestRotStart;
                        oTestRange.m_afRotEnd[nVariable] = (fTestRotStart + fWidth);
                        
                        nCompleteness(oTestRange);
                        
                        // Have we found a solution?
                        if( oTestRange.m_fCompleteness >= m_fCompletenessMin  - m_fCompletenessTol && 
                            oTestRange.m_fRedundancy   >= m_fRedundancyMin ) 
                        {
                            // Save this solution, and
                            bSolutionFound = true;
                            if (nSaveRanges>=nMaxSaveRanges)
                                break;
                            
                            oBestRangesSave[nSaveRanges] = oTestRange;

                            nSaveRanges++;
                        }
                        
                        if(m_bFixRotStart)
                            break;    // we don't need to loop in this case, because  the start cannot be moved
                    }
                }
                
                // If we found a solution, don't look at finer intervals just yet.
                if (bSolutionFound)
                    break;
                
                // Otherwise, don't check this coursness again for this variable.
                anLastWorkingWidth[nVariable] = nWidth + 1;
            }
            
            // Choose randomly which solution to use.
            if (bSolutionFound)
            {
                oBestRangeInner = oBestRangesSave[rand() % nSaveRanges];
                oBestRangeInner.m_bIsAvailable = true;
            }

            // Check for convergence.
            bConverged = true;
            for (nVariable = 0; nVariable < nVariables; nVariable++) 
            {
                if (anWidth[anLastWorkingWidth[nVariable]]!=-1)
                    bConverged = false;
            }
            
            nInnerLoop++;
        }

        if( oBestRangeInner < oBestRangeOuter )
        {

            oBestRangeOuter = oBestRangeInner;
            
            printf("Better solution found...\n");

            //double      dTemp = 0.0;
            //m_poScanInfo->nPrintResultScanTable(0, true, dTemp, 0U);

            oBestRangeOuter.nPrint(*this);
            
            nLastOuterLoopWorked = nOuterLoop;
        }

        time(&oTimeNow);
        bTimeOut =  (((int) (oTimeNow - oTimeStart)>m_nMaxTime) || (nOuterLoop - nLastOuterLoopWorked>5));
        nOuterLoop++;
           
    } while (!bTimeOut);

    oUserRange = oBestRangeOuter;

    return 0;
}

int Cstrategy::nComputeBestScanToAdd(itr<int>& anScansInUse) 
{
    int nx;
    int nScan,nPackedHKL,nNewPackedHKL;
    int nNewScan;
    itr<int> anScanIndex;
    itr<int> anScanEnd;
    int  a2nTestNew[2];         // [0] == spacegroup HKL, [1] == triclinic HKL
    int  nBestScan;             // Best Scan Found.
    int  a2nBestNew[2];         // Best a2nTestNew values for above.
    int nPass;
    int nFound;
    int nRedundancyToCheck;

    if (!anScansInUse.size())
        return 1;
    nRedundancyToCheck = 0;    // We count only reflections that have a redundancy <= this number

    do {
        nRedundancyToCheck++;
        
        nBestScan = -1;
        // Loop through every scan that is not already in use.
        for (nNewScan = 0; nNewScan < m_anScanStart.size(); nNewScan++)
        {
            if( !bScanUsed(nNewScan, 0) )
                continue;
            
            // See if this new scan is in use.
            if( -1 != anScansInUse.find(nNewScan) )
                continue;

            for (nPass = 0; nPass < 2; nPass++)
            {
                int* pnIndex;
                int  nFI_nPackedHKL;
                
                a2nTestNew[nPass] = 0;
                
                if (nPass == 0) 
                {
                    // Use the spacegroup's HKL values.
                    pnIndex = m_pnIndex;
                    nFI_nPackedHKL = m_nFI_nPackedHKL;
                } 
                else 
                {
                    // Use the triclinic/anom HKL values.
                    pnIndex = m_pnIndex2;
                    nFI_nPackedHKL = m_nFI_nPackedHKL2;
                }
                
                // Initialize the index pointers
                -anScanIndex + m_anScanStart[nNewScan];
                -anScanEnd + m_anScanEnd[nNewScan];
                
                for (nx=0;nx<anScansInUse.size();nx++)
                {
                    anScanIndex + m_anScanStart[anScansInUse[nx]];
                    anScanEnd + m_anScanEnd[anScansInUse[nx]];
                }
                
                while (anScanIndex[0]<=anScanEnd[0])
                {
                    
                    // Make all packed HKL values of existing scans <= the packed HKL value of the new scan.
                    nPackedHKL = (*m_poReflnlist)[pnIndex[anScanIndex[0]]].nGetField(nFI_nPackedHKL);
                    
                    nFound = 0;
                    
                    for (nScan = 1; nScan<anScanIndex.size(); nScan++)
                    {
                        while ((anScanIndex[nScan]+1<=anScanEnd[nScan]) && 
                            ((nNewPackedHKL = (*m_poReflnlist)[pnIndex[anScanIndex[nScan]+1]].nGetField(nFI_nPackedHKL))
                            <=nPackedHKL)) 
                        {
                            anScanIndex[nScan]++;
                        
                            if (nPackedHKL == nNewPackedHKL)
                            {
                                // If the HKL is equal, then we can just drop out of this loop now.
                                nFound++;
                            }
                        }
                    }
                    
                    if (nFound<=nRedundancyToCheck) 
                        a2nTestNew[nPass]++;
                    
                    // Increment the scan index until we get a different packed HKL.
                    while( (anScanIndex[0] + 1 <= anScanEnd[0]) && 
                           ((*m_poReflnlist)[pnIndex[anScanIndex[0]+1]].nGetField(nFI_nPackedHKL)==nPackedHKL)) 
                    {
                        anScanIndex[0]++;
                    }

                    anScanIndex[0]++;
                }
            }
            
            // Now see if we have a better solution.
            
            if( nBestScan == -1 )
            {
                // We have a new scan.
                a2nBestNew[0] = a2nTestNew[0];
                a2nBestNew[1] = a2nTestNew[1];
                nBestScan = nNewScan;
            } 
            else if ( a2nTestNew[0] > a2nBestNew[0] ) 
            {
                // This covers more spacegroup HKL values.  Use it.
                a2nBestNew[0] = a2nTestNew[0];
                a2nBestNew[1] = a2nTestNew[1];
                nBestScan = nNewScan;
            } 
            else if( a2nBestNew[0] == 0 && a2nTestNew[1] > a2nBestNew[1] ) 
            {
                // This covergs more triclinic HKL values.  Use it.
                a2nBestNew[0] = a2nTestNew[0];
                a2nBestNew[1] = a2nTestNew[1];
                nBestScan = nNewScan;
            }
        }
        
        // It is possible that we already are covering every single unique HKL value.  In that case, we really want
        // to count up to the next higher redundancy.
        if( nBestScan == -1 ) 
            return -1;

    } while( a2nBestNew[0] == 0 && a2nBestNew[1] == 0 );

    return nBestScan;
}

int Cstrategy::nUnPackRedundancy()
{   
    int     nSphereRedundancy = 0;   
    
    // A redundancy of 0 indicates the user wanted to collect an entire hemisphere.

    nSphereRedundancy = m_poCrystal->m_poSpacegroup->nGetNumLaue() * 2;
    
    if( m_bAnom ) 
        nSphereRedundancy /=2;

    if( m_fRedundancyMin == 0.0 )
    {
        m_fRedundancyMin = nSphereRedundancy;
    } 
    else if (m_fRedundancyMin < 0.0)
    {
        m_fRedundancyMin = nSphereRedundancy*ABS(m_fRedundancyMin);
    }

    m_fRedundancyMin -= 0.01;

    return 0;
}

int Cstrategy::nCalcScans()
{
    CstrategyRange      oTestRange(this);      
    CstrategyRange      oBestRange(this);

    CstrategyRange      oBestSingleScanRange(this);  // we need that to see if the best multiple scan range any better tahn the single scan range

    int         nScan = 0;
    int         ni = 0;
    int         nx = 0;
    int         nStat = 0;

    itr<int>    anUse;

    /// RB 02/05/07 I am commenting out Thad's outside loop that handles redundancy, as it does not make much sense.
    // All redundancy can be handled in the inside loop. Even if I left Thad's code in place, there would be a bug in that
    // the algorithm does not check the real cumulative redundancy after each pass. For example, if the user requires 
    // a redundancy of 3 for a triclinic crystal, Thad would break it down into 2 (sphere) plus 1, but would not check after
    // the first pass whether redundancy already achieved is exactly 2 or higher.
    
    //int         nSphereRedundancy = m_poCrystal->m_poSpacegroup->nGetNumLaue()*2;
    //if (m_bAnom)
    //nSphereRedundancy /= 2;
    //double      fInputRedundancy = m_fRedundancyMin;
    //double      fCumulRedundancy = 0.0;

    -m_oBestRange;
    
    int     nPass = -1;
    
    double  dRotationAboveLimit = 0.0;
    
    std::vector<CstrategyRange>     vCandidateRanges;

    bool    bTemp = false;
    
    //// In this outer most loop, we must handle redundancy.//
    //do
    //{
    //    // We need to handle this in pieces.
    //    m_fRedundancyMin = min(nSphereRedundancy - 0.01,fInputRedundancy - fCumulRedundancy);
    //    
    //    if( m_fRedundancyMin <= 0.2 && m_fRedundancyMin > -1.0 )
    //    {
    //      //break;
    //    }
        
        nPass++;

        //-anUse;
        vCandidateRanges.clear();
    
        oTestRange.m_bIsAvailable = TRUE;
        oBestRange.m_bIsAvailable = FALSE;

        // Determine if we have to do anything.
        if( nPass == 0 && m_anFixedScans.size() )
        {
            nClearRange(oTestRange);
            nCompleteness(oTestRange);

            //if( oTestRange.m_fCompleteness >= m_fCompletenessMin  - m_fCompletenessTol && oTestRange.m_fRedundancy >= fInputRedundancy )
            if( oTestRange.m_fCompleteness >= m_fCompletenessMin  - m_fCompletenessTol && oTestRange.m_fRedundancy >= m_fRedundancyMin )
            {
                printf("INFO:  Input previous reflection list is sufficient to meet conditions.\n       No new scans will be generated.\n");
        
                nStat = 0;
        
                m_oBestRange += oTestRange;
                nCompleteness(m_oBestRange);
        
                //break;
                return 0;
            }
        }

        // First, we determine if one scan (from the first set of scans) is enough to satisfy the criteron.    
        printf("Determining if one scan will solve system...\n");
        m_poScanInfo->vSetMultipleScanSolution(false);

        for(nScan = 0; nScan < m_anScanStart.size(); nScan++)
        {
            if( !bScanUsed(nScan, 1) )
                continue;
            
            nClearRange(oTestRange);
        
            oTestRange.m_anScanNumbers + nScan;   // RB: notice that the scan numbers for non-fixed scans start from 0
            oTestRange.m_afRotStart + m_afScanRotMin[nScan];
            oTestRange.m_afRotEnd + m_afScanRotMax[nScan];
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // RB: Although a scan may fail the completeness/redundancy checks, we still need to save the best range for each scan to report
            nCalcScanSingle(oTestRange, MAX_SEARCH_STEP_SINGLE_SCAN, oTestRange.m_anScanNumbers.size()-1); // m_anScanNumbers.size()-1==m_anFixedScans.size(), because this range object has a single scan after the fixed scans
            nFromRangeToStrategyScansArray(oTestRange, DTREK_DTSTRATEGY_FRTSSA_UPDATE_ALL_RESULTS); // Update strategy scans array.
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if( oTestRange.m_fCompleteness >= m_fCompletenessMin - m_fCompletenessTol && 
                oTestRange.m_fRedundancy   >= m_fRedundancyMin )
            {
                //anUse + nScan;
                vCandidateRanges.push_back(oTestRange);
            }
            
            // If we don't find a solution, we will use this as the first solution.
            //if( !oBestRange.m_bIsAvailable || oTestRange < oBestRange )
            if( !oBestRange.m_bIsAvailable || bRangesCompareConsiderTargets(oTestRange, 
                                                                            oBestRange, 
                                                                            m_fCompletenessMin - m_fCompletenessTol, 
                                                                            m_fRedundancyMin,
                                                                            bTemp) )
            {
                oBestRange = oTestRange;
                oBestRange.m_bIsAvailable = TRUE;
            }
        }

        oBestSingleScanRange = oBestRange; // memorize the best single scan range 

        if (!oBestRange.m_bIsAvailable)
        {
            printf("ERROR:  No scans found, or not enough scans to meet requested comleteness and/or redundancy.\n");
            nStat = 1;
        } 
        //else if( anUse.size() > 0 ) 
        else if( vCandidateRanges.size() > 0 ) 
        {
            //printf("%d scans were found which may meet completeness and redundancy \n      in a single scan.\n", anUse.size());
            printf("%d scan(s) were found which meet target completeness (with tolerance)\n"
                   "and redundancy in a single scan.\n", 
                   vCandidateRanges.size());
            
            if( m_fCompletenessTol > 0.0 )
            {
                double      dTemp=0.0;
                m_poScanInfo->nPrintResultScanTable(0, true, dTemp, DTREK_STRATSCANS_PRST_ALL);
                
                vSortFindBestExpandableRange(vCandidateRanges);
                oBestRange = vCandidateRanges[0];
                //printf("INFO:  Best scan solution:\n");
                //oBestRange.nPrint(*this);
            }
        }
        else if(m_nTotalNumberOfScansAllowed > 1 && m_poScanInfo && m_poScanInfo->nGetNumScans() > 1) // check if we could do multiple scans
        {
            // We require multiple scans.  
            // Compute the intersections matrix.  This will encode information about how much
            // different scans intersect.
            
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // RB: This is a fix to the maximum rotation range problem. The purpose of this call is to save what the StrategyScans array 
            // knows about the "best restricted" ranges for all scans. For example, the allowed rotation range is -90 to 90, but the maximum
            // user-requested range is only 50. The StrategyScans array would have the best 50 degree ranges for all scans, and we need 
            // to store those now, so when the multiple scan strategy is invoked, we are not using (-90 - 90) for any of the scans, as Thad's
            // algorithm would do. 
            nFromStrategyScansArrayToScanRotationArrays();
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            printf("Cannot satisfy conditions with a single scan.\n");
            
            // Now that we have called nFromStrategyScansArrayToScanRotationArrays, we are resetting the state:
            m_poScanInfo->vSetMultipleScanSolution(true);
            //m_bMultipleScanSolution = true;
            
            printf("Computing intersections of available scans...\n");
            
            // Add the longest scan for starters.
            -anUse; 
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // RB: Well, there could be a problem here if m_dTotalRotWidthAllowed is shorter than the rotation width of the "longest" Thad's scan 
            // The users should be educated to have m_dTotalRotWidthAllowed large enough.
            // For now, let's just fix that longest scan rotation.
            nx = oBestRange.m_anScanNumbers.last();
            
            if( m_afScanRotMax[nx] - m_afScanRotMin[nx] > m_dTotalRotWidthAllowed )
            {
                printf("\n\nWARNING: Scan #%d rotation width exceeds the total maximum rotation width allowed for a multiple scan solution.\n", nx);
                printf("Please consider increasing the total maximum rotation width allowance.\n\n");
                m_afScanRotMax[nx] = m_afScanRotMin[nx] + m_dTotalRotWidthAllowed;
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            anUse + oBestRange.m_anScanNumbers.last();
            printf("Using Scan  %5d\n",anUse.last());

            // We need to continue choosing extra scans until we meet the requirements.
            nStat = 0;
            do
            {
                if( oTestRange.m_anScanNumbers.size() >= m_nTotalNumberOfScansAllowed ||
                    oTestRange.dGetTotalRotWidth()    >= m_dTotalRotWidthAllowed      ||
                    -1 == (nScan = nComputeBestScanToAdd(anUse))
                  )
                {
                    printf("INFO:  Could not satisfy completeness and/or redundancy.\n");

                    printf("INFO:   Best completeness was %.2lf\n",oTestRange.m_fCompleteness);
                    printf("INFO:   Best redundancy was %.2lf\n",oTestRange.m_fRedundancy);
                    
                    nStat = 1;
                    break;
                }
                
                anUse + nScan;
                
                printf("Adding Scan %5d\n",anUse.last());
                
                nClearRange(oTestRange);
                oTestRange.m_anScanNumbers + anUse;
                
                // RB: We want to put all scans "in use" to oTestRange, but we don't want to exceed the total allowed rotation range
                for(nx=0; nx < anUse.size(); nx++)
                {
                    // Figure out whether we will exceed the limit if we add this scan
                    dRotationAboveLimit = ( oTestRange.dGetTotalRotWidth() + m_afScanRotMax[anUse[nx]] - m_afScanRotMin[anUse[nx]] ) - 
                                          m_dTotalRotWidthAllowed;
                    
                    if( dRotationAboveLimit <= 0.0 )
                        dRotationAboveLimit = 0.0;    // we still have room for more rotation
                    else
                    {
                        // This should not really happen, because it would mean that oTestRange has exceeded the total rotation
                        // limit even before we tried to add this scan... But just for safety:
                        if( m_afScanRotMax[anUse[nx]] - dRotationAboveLimit < m_afScanRotMin[anUse[nx]] ) 
                        {
                            dRotationAboveLimit = m_afScanRotMax[anUse[nx]] - m_afScanRotMin[anUse[nx]];  // force this scan to have zero width, when added to oTestRange 
                        }
                    }

                    // RB: OK, now add the scan begin and end to oTestRange
                    oTestRange.m_afRotStart + m_afScanRotMin[anUse[nx]];
                    // RB: We will shorten the added scan by dRotationAboveLimit...
                    // Todo: test what is better: chop that piece off the end or off the begin of the scan
                    oTestRange.m_afRotEnd + (m_afScanRotMax[anUse[nx]] - dRotationAboveLimit); 
                }
                
                nCompleteness(oTestRange);
            
            }while( oTestRange.m_fCompleteness < m_fCompletenessMin - m_fCompletenessTol || 
                    oTestRange.m_fRedundancy < m_fRedundancyMin );

            if( nStat != 1 ) // 0 means OK by the algorithm; 2 means OK by the user
            {
                printf("Calculating restricted ranges for scans...\n");

                if( !nCalcScanMultiple(oTestRange) ) 
                {
                    if( m_fCompletenessTol > 0.0 )
                        bTryExpandRange(oTestRange, oTestRange.dGetTotalRotWidth() * m_fRotationTol / 100.0 );
                    
                    oBestRange = oTestRange;
                }
                else 
                {
                    printf("WARNING:  Could not compute scan!\n");
                    nStat = 1;
                }
            }
            else  
            {    
                // Since we are going to exit with an error, so we need to make sure 
                // that our multi-scan solution is better than the single scan.
                // Also, we need to save the best range anyway, just to be able 
                // to show the best possible solution to the user.
                nCompletenessMultiScan(oTestRange);
                bool    bTemp = false;
                if( bRangesCompareConsiderTargets(oBestSingleScanRange,
                                                  oTestRange, 
                                                  m_fCompletenessMin - m_fCompletenessTol, 
                                                  m_fRedundancyMin,
                                                  bTemp) )
                {
                    printf("\nThe best multiple scan solution is no better than the best single scan solution.\n"
                           "Keeping the best single scan solution.\n");
                    m_poScanInfo->vSetMultipleScanSolution(false);
                    m_oBestRange += oBestSingleScanRange;
                }
                else
                    m_oBestRange += oTestRange;
                
                nCompletenessMultiScan(m_oBestRange);
            }
        } 
        else //RB 06/19/2006 Adding this case to handle a situation when the total allowed number of scans is 1.
        {
            nStat = 1;
            m_oBestRange += oBestRange;
            nCompletenessMultiScan(m_oBestRange);
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if( nPass > 0 ) // if we've had more than one pass, then we have a multiple scan situation, even if during the last pass we had a single scan situation
            m_poScanInfo->vSetMultipleScanSolution(true);

        if( nStat == 1 )
        {
            return 1;
        }

        m_oBestRange += oBestRange;
        m_anExcludeScans += oBestRange.m_anScanNumbers;  // so that we don't use them during the next pass
        
        ////////////////////////////////////////////////////////////////////////////
        // Take care of the maximum allowed number of scans and total rotation range
        m_nTotalNumberOfScansAllowed -= oBestRange.nGetNumberOfScans();
        if( m_nTotalNumberOfScansAllowed < 0 )
            m_nTotalNumberOfScansAllowed = 0; // safety

        m_dTotalRotWidthAllowed -= oBestRange.dGetTotalRotWidth();
        if( m_dTotalRotWidthAllowed < 0.0 )
            m_dTotalRotWidthAllowed = 0.0;  // safety
        /////////////////////////////////////////////////////////////////////////////

        if( m_poScanInfo->bIsMultipleScanSolution() )
            nCompletenessMultiScan(m_oBestRange);
        else
            nCompleteness(m_oBestRange);// If it is a single scan solution, then the statistics for the whole range 
                                        // are the same as for the individual scan
                                        // Therefore, we do not need to do statistics for the scan in this case.


        //fCumulRedundancy += oBestRange.m_fRedundancy;
    
        if( oBestRange.m_fRedundancy <= 0.0 )
        {
            printf("ERROR:  Problem with data redundancy???\n");
            return 1;
        }

        if (nStat == 2)
        {
            return 2;
        }
    //}while (fCumulRedundancy < fInputRedundancy );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////
// Read the header, find out how many axes we have and save that information
void Cstrategy::vSetAxisInfo(Cimage_header* poHeader)
{
    if( !poHeader )
        return; // safety

    if( !m_poScanInfo )
        return; //safety

    m_poScanInfo->vClearAxes(); // initialize


    // How many crystal gonio values are there?
    Cstring         strCrystalPrefix(D_K_CrystalPrefix);

    int             nAxes = 0;
    if( 0 != poHeader->nGetValue(strCrystalPrefix + D_K_GonioNumValues, &nAxes) )
        return;

    if( nAxes <= 0 )
        return;

    // Get axes names
    Cstring*        psGonioNames = new Cstring[nAxes];
    if( 0 != poHeader->nGetValue(strCrystalPrefix + D_K_GonioNames, nAxes, psGonioNames) )
    {
        delete [] psGonioNames;
        return;
    }

    // Get rotation axis name
    Cstring         strScanPrefix(D_K_ScanPrefix);

    Cstring         strRotAxis("");
    if( 0 != poHeader->nGetValue(strScanPrefix + D_K_RotAxisName, &strRotAxis) )
        return;

    int     ii = 0;
    // Get axes collision offsets
    double*         pdGonioCollisionOffsets = new double[nAxes];
    for(ii=0; ii < nAxes; ii++)
        pdGonioCollisionOffsets[ii] = 0.0;
    poHeader->nGetValue(strCrystalPrefix + D_K_GonioCollisionOffsets, nAxes, pdGonioCollisionOffsets);
    
    
    // Now add crystal axes objects 
    for(ii=0; ii < nAxes; ii++)
    {
        m_poScanInfo->vAddAxis(psGonioNames[ii], strRotAxis==psGonioNames[ii], pdGonioCollisionOffsets[ii]);
    }

    delete [] psGonioNames;
    delete [] pdGonioCollisionOffsets;
}
/////////////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vSaveDetectorPosition(Cimage_header* poHeader, double dTwoTheta, double dDistance)
{
     Cdetector      oDetector(*poHeader);
     
     double         dTwoThetaOrig = oDetector.fGetSwing();
     double         dDistOrig = oDetector.fGetDistance();
    
     m_poScanInfo->vSetTwoTheta(dTwoTheta > cdStratDetUnrealisticValue ? dTwoTheta : dTwoThetaOrig);
     m_poScanInfo->vSetDistance(dDistance > cdStratDetUnrealisticValue ? dDistance : dDistOrig);
}
/////////////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vAddStrategyScan(Cimage_header* poHeader, int nScanIndex, int nNumAxes)
{
    if( !m_poScanInfo )
        return;

    m_poScanInfo->vAddScan(poHeader, nScanIndex, nNumAxes, D_K_StrategyTestPrefix);
    
    int     nJustAddedScanIndex = m_poScanInfo->nGetNumScans() - 1;

    // We need to initialize the rotation range
    m_poScanInfo->vSetScanRotation(nJustAddedScanIndex, CstrategyScan::enSingleScan, 0.0, 0.0);
    m_poScanInfo->vSetScanRotation(nJustAddedScanIndex, CstrategyScan::enMultipleScan, 0.0, 0.0);
}
//////////////////////////////////////////////////////////////////////////////////////////////////
// If the best range is a single scan solution, this function will print the single scan summary.
// If, however, the best range is a multiple scan solution, this function will print both
// the single scan and multiple scan summaries.
void Cstrategy::vPrintSummary(bool bPrintSingleScanSolution)
{
    if( !bPrintSingleScanSolution && !m_poScanInfo->bIsMultipleScanSolution() )
        return; // nothing to do
    
    bool        bPrintTheActualSolution = ( bPrintSingleScanSolution && !m_poScanInfo->bIsMultipleScanSolution() ) ||
                                          ( !bPrintSingleScanSolution && m_poScanInfo->bIsMultipleScanSolution() );

    
    DTREK_WORD      wCtrl = 0U;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( bPrintTheActualSolution )
        wCtrl |= DTREK_DTSTRATEGY_FRTSSA_SET_USED;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( m_bZeroOmegaWithPhiOffset && bPrintTheActualSolution )
        wCtrl |= DTREK_DTSTRATEGY_FRTSSA_UPDATE_ROT_RANGE; 
    
    //////////////////////////////////////////////////////////////////////////////////////
    if( wCtrl != 0U )
        nFromRangeToStrategyScansArray(m_oBestRange, wCtrl);
    
    vPrintSummaryScanTable(bPrintSingleScanSolution);        
}
////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vPrintSummaryScanTable(bool bSingleScanSummary)
{
    if( !m_poScanInfo )
        return; // safety
    
    if( !bSingleScanSummary && !m_poScanInfo->bIsMultipleScanSolution() )
        return; // inconsistency: we don't have a multiple scan solution, so it should not be printed
    
    if( 0 == m_poScanInfo->nGetNumScans() )
        return; // nothing to do

    if( 1 > nGetNumAxes() )
        return; // Should not happen. There should be at least one rotation axis. 
    //////////////////////////////////////////////////////////////////////////////
    
    CstrategyScan::enStrategyScanType      eScanType = bSingleScanSummary ? CstrategyScan::enSingleScan : CstrategyScan::enMultipleScan;
    if( eScanType == CstrategyScan::enMultipleScan )
        m_poScanInfo->vSortScans(); // set priorities as to which scan need to be collected first

    Cstring     sTableTitle = bSingleScanSummary ? "\n\nSingle Scan Summary\n" : "\n\nBest Multiple Scan Solution\n";
    printf(sTableTitle.string());

    printf("Detector position: 2Theta=%7.2f deg Distance=%7.2f mm\n", m_poScanInfo->dGetTwoTheta(), m_poScanInfo->dGetDistance());
    
    Cstring     sHorizLine("-------------------------------------------------------------------------------\n");
    printf(sHorizLine.string());
    
    // We need to figure out the axes names for the table header.
    // Assume every scan has the same axes names.
    
    Cstring     sHeader ("Scan     Gonio axes names        Rot     Rot    Scan   %%Comp  Redun   +/-  Used\n");

    Cstring     sHeader2(" num");
    
    Cstring     sTemp("");
    Cstring     sTemp2("");

    int                 nAxisNameLength = 0;
    int                 nNumberOfPaddingSpaces = 0;

    const int           c_nMaxAxisColumnWidth = 8;// assuming that an axis name would not be longer than 7 characters + a space for separation
    
    int     ii=0, jj=0, kk=0;

    for(jj=0; jj < nGetNumAxes(); jj++)
    {
        vGetAxisName(jj, sTemp);
        
        nAxisNameLength = sTemp.length();
        nNumberOfPaddingSpaces = c_nMaxAxisColumnWidth - nAxisNameLength;
        
        sTemp2 = "";
        for(int kk=0; kk < nNumberOfPaddingSpaces; kk++)
        {
            sTemp2 += " ";
        }
        sTemp2 += sTemp;
        
        sHeader2 += sTemp2; 
    }
    
    sHeader2 += "   Start     End   Width          dancy        scan?\n";

    printf(sHeader.string());
        
    printf(sHeader2.string());
    
    printf(sHorizLine.string());
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    bool        bUsed = false;
    
    double      dRotBegin = 0.0;
    double      dRotEnd = 0.0;
    double      dWidth = 0.0;
	
	//shijie 2008.03.17 : var for padding info
	std::vector<Cstring> sPaddingMsg;
	char msgStr[2048];
	msgStr[0] = 0;
	int scanNum;
	CstrategyScan::enStrategyScanRangePaddingType enPadType;
	CstrategyScan::enStrategyScanRangePaddingType enShiftType;
	Cstring sFormat;

    m_dTotalRotWidth = 0.0;
    
    Cstring     strYesNo("");

    bool        bAtLeastOneRowExist = false;

    for(ii=0; ii < m_poScanInfo->nGetNumScans(); ii++)                                 
    {
        sFormat = "%7.2f %7.2f %7.2f";
        bUsed = m_poScanInfo->bIsScanUsed(ii, eScanType);

        // Multiple scan solution should only list used scans
        if( !bSingleScanSummary && !bUsed )
            continue;
        else
            bAtLeastOneRowExist = true;
        
		scanNum = m_poScanInfo->nGetOriginalScanIndex(ii);
        printf("%4d ", scanNum); // scan number
        for(jj=0; jj < nGetNumAxes(); jj++)
        {
            printf("%7.2f", m_poScanInfo->dGetScanAxisValue(ii, jj) );

            ////////////////////////////////////////
            // Specify whether it is a rotation axis
            if( bIsRotationAxis(jj) )
                printf("*");
            else
                printf(" ");
        }
        ////////////////////////////////////////////////////////////////////////////
        
        dRotBegin = m_poScanInfo->dGetScanRotBegin(ii, eScanType);
        dRotEnd   = m_poScanInfo->dGetScanRotEnd  (ii, eScanType);     
        dWidth    = fabs(dRotBegin - dRotEnd);
        
		//shijie 2008.03.17 : for padding info
		enPadType = m_poScanInfo->enGetPaddingType(ii, eScanType);
		enShiftType = m_poScanInfo->enGetShiftingType(ii, eScanType);

		if(enShiftType != CstrategyScan::enNone)	//shifted
		{
			sprintf(msgStr, "%d : range shifted to start at whole angle value with integral number of images;",
							scanNum);
			sFormat = "%7.2f~%7.2f~%7.2f";
		}
		else if(enPadType != CstrategyScan::enNone)
		{
			if(enPadType == CstrategyScan::enFrontExtended)
			{
				sprintf(msgStr, "%d : range starting angle extended to contain integral number of images;",
								scanNum);
				sFormat = "%7.2f^%7.2f %7.2f";
			}
			else if(enPadType == CstrategyScan::enEndExtended)
			{
				sprintf(msgStr, "%d : range ending angle extended to contain integral number of images;",
								scanNum);
				sFormat = "%7.2f %7.2f^%7.2f";
			}
			else if(enPadType == CstrategyScan::enBothEndExtended)
			{
				sprintf(msgStr, "%d : range angles extended to contain integral number of images;",
								scanNum);
				sFormat = "%7.2f^%7.2f^%7.2f";
			}
			else	//both ends shrinked 
			{
				sprintf(msgStr, "%d : range angles shrinked to contain integral number of images;",
								scanNum);
				sFormat = "%7.2f#%7.2f#%7.2f";
			}
		}
        else
            msgStr[0] = 0;
		
		Cstring tmpStr = msgStr;
		if(tmpStr.length() > 0)
			sPaddingMsg.push_back(tmpStr);
		
        if( bUsed )
            m_dTotalRotWidth += dWidth;
        
        ///////////////////////////////////////////////////////////////////////////
        // Format Rotation Begin, End and Width
        if( bUsed || dWidth > 0.0 )  // either "used" or "tested, but not accepted". For multiple scan solutions "bUsed" should be true, because we skip non-used scans
        {
			printf(sFormat.string(), dRotBegin, dRotEnd, dWidth);
        }
        else
        {
            printf("%7s %7s %7s",  "    ---", 
                                   "    ---",
                                   "    ---"); 
        }
        ///////////////////////////////////////////////////////////////////////////
        
        if( bUsed || dWidth > 0.0 ) // a single scan that is either "used" or "was tested, but not accepted"
        {
            printf(" %7.2f %6.2f %5.2f", m_poScanInfo->dGetScanCompleteness (ii, eScanType), 
                                         m_poScanInfo->dGetScanRedundancy   (ii, eScanType), 
                                         m_poScanInfo->dGetScanRedundancyDev(ii, eScanType));
        }
        else
        {
            printf(" %7s %6s %5s", "    ---", 
                                   "   ---",
                                   "  ---");
        }
        ////////////////////////////////////////////////////////////////////////////////////

        strYesNo = bUsed ? " Yes" : "   ";
        printf(" %4s\n", strYesNo.string());
    }
    
    printf(sHorizLine.string());
    ////////////////////////////////////////////////////////////////////////////////////////
    
    if( !bSingleScanSummary )
    {
        printf("All used scans                               %7.2f %7.2f %6.2f %5.2f\n", 
                                                          m_dTotalRotWidth, 
                                                          m_oBestRange.m_fCompleteness, 
                                                          m_oBestRange.m_fRedundancy,
                                                          m_oBestRange.m_fRedundancyDev);
    }

    printf("\n* Scan rotation axis\n");
    printf("+/- is a standard deviation of redundancy among resolution bins\n");
	for(int i=0; i<sPaddingMsg.size(); i++)
	{
		printf(sPaddingMsg.at(i));
		printf("\n");
	}
	sPaddingMsg.clear();
	printf("\n");
   

    if( (bSingleScanSummary && !m_poScanInfo->bIsMultipleScanSolution()) || (!bSingleScanSummary && m_poScanInfo->bIsMultipleScanSolution()) )
    {
        printf("Totals for used scans:\n");
        printf("Cumulative width: %8.2f\n", m_dTotalRotWidth);
        printf("Completeness:     %8.2f\n", m_oBestRange.m_fCompleteness);
        printf("Redundancy:       %8.2f +/- %8.2f\n", m_oBestRange.m_fRedundancy, m_oBestRange.m_fRedundancyDev);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IMPORTANT: If there are fixed scans, then it is assumed that the range, passed to this function must have them. 
// This is actually achieved in the code by calling CStrategy::nClearRange() on that range at some point before calling this function.
int Cstrategy::nFromRangeToStrategyScansArray(CstrategyRange& oRange, DTREK_WORD wCtrl)
{
    if( !m_poScanInfo )
        return 0; // safety check;
    
    CstrategyScan::enStrategyScanType      eScanType = !m_poScanInfo->bIsMultipleScanSolution() ? 
                                                       CstrategyScan::enSingleScan : CstrategyScan::enMultipleScan;
    int     ii = 0;
    int     nNumberOfScansModified = 0;
    
    int     nFixedScansNumberOffset = m_anFixedScans.size();

    int     nUsedScansNumber = oRange.m_anScanNumbers.size() - nFixedScansNumberOffset; // fixed scans are supposed to be at the beginning of the m_anScanNumbers array

    if( 0 >= nUsedScansNumber || nUsedScansNumber > m_poScanInfo->nGetNumScans() )  // nUsedScansNumber could be negative, if oRange is empty AND there exist a fixed scan
        return 0; // just a safety feature

    //if( wCtrl & DTREK_DTSTRATEGY_FRTSSA_RESET_USED )
    //{
    //    for(ii=0; ii < m_poScanInfo->nGetNumScans(); ii++)
    //    {
    //        m_poScanInfo->vSetScanUsed(ii, eScanType, false);
    //    }
    //}

    for(ii = 0; ii < nUsedScansNumber; ii++)
    {
        int     nUsedScanIndex = oRange.m_anScanNumbers[ii+nFixedScansNumberOffset];  
        
        if( nUsedScanIndex < 0 || nUsedScanIndex >= m_poScanInfo->nGetNumScans() )
            continue;   // safety check, plus if there are fixed scans they are not in m_poScanInfo

        if( oRange.m_afRotStart[ii+nFixedScansNumberOffset] == oRange.m_afRotEnd[ii+nFixedScansNumberOffset] )
            continue;  // cannot use a scan with rotation width equal to zero

        nNumberOfScansModified++;

        m_poScanInfo->vSetScanUsed(nUsedScanIndex, eScanType, wCtrl & DTREK_DTSTRATEGY_FRTSSA_SET_USED);
        
        if( wCtrl & DTREK_DTSTRATEGY_FRTSSA_UPDATE_ROT_RANGE )
        {
            m_poScanInfo->vSetScanRotation(nUsedScanIndex, eScanType,
                                           oRange.m_afRotStart[ii+nFixedScansNumberOffset], 
                                           oRange.m_afRotEnd  [ii+nFixedScansNumberOffset]);
        }

        if( wCtrl & DTREK_DTSTRATEGY_FRTSSA_UPDATE_STATS )
        {
            m_poScanInfo->vSetScanCompleteness  (nUsedScanIndex, eScanType, (double)oRange.m_fCompleteness);
            m_poScanInfo->vSetScanRedundancy    (nUsedScanIndex, eScanType, (double)oRange.m_fRedundancy);
            m_poScanInfo->vSetScanRedundancyDev (nUsedScanIndex, eScanType, (double)oRange.m_fRedundancyDev);
        }
    }

    return nNumberOfScansModified;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cstrategy::nFromStrategyScansArrayToScanRotationArrays(bool bCheckRotWidthOnly)
{
    if( !m_poScanInfo )
        return false;  // safety
    
    CstrategyScan::enStrategyScanType      eScanType = !m_poScanInfo->bIsMultipleScanSolution() ? 
                                                       CstrategyScan::enSingleScan : CstrategyScan::enMultipleScan;
    
    int     nFixedScansNumberOffset = m_anFixedScans.size();

    int     nTotalNumberOfScansInStrategyScanArray = m_poScanInfo->nGetNumScans();
    
    int     nTotalNumberOfScansInRotationStartArray   = m_afScanRotMin.size()     - nFixedScansNumberOffset;
    int     nTotalNumberOfScansInRotationEndArray     = m_afScanRotMax.size()     - nFixedScansNumberOffset;
    int     nTotalNumberOfScansInRotationWidthArray   = m_afScanRotRange.size()   - nFixedScansNumberOffset;

    if( nTotalNumberOfScansInRotationStartArray != nTotalNumberOfScansInRotationEndArray  )
        return false;  // safety

    if( nTotalNumberOfScansInRotationEndArray   != nTotalNumberOfScansInStrategyScanArray )
        return false; // safety

    if( nTotalNumberOfScansInRotationWidthArray != nTotalNumberOfScansInStrategyScanArray )
        return false; // safety
    
    double      dInputRotStart   = 0.0;
    double      dInputRotEnd     = 0.0;
    double      dInputRotWidth   = 0.0;
    
    for(int iScan = 0; iScan < nTotalNumberOfScansInStrategyScanArray; iScan++)
    {
        dInputRotStart = m_afScanRotMin[iScan];
        dInputRotEnd   = m_afScanRotMax[iScan];
        dInputRotWidth = m_afScanRotRange[iScan];
        
        // If the input rotation width between start and end does satisfy 
        // the input maximum allowed rotation width, there is no need to modify this scan.
        if( bCheckRotWidthOnly && fabs(dInputRotStart-dInputRotEnd) <= dInputRotWidth )
            continue;

        m_afScanRotMin[iScan] = m_poScanInfo->dGetScanRotBegin(iScan, eScanType);
        m_afScanRotMax[iScan] = m_poScanInfo->dGetScanRotEnd(iScan, eScanType);
    }    
    
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vAddFixedScan(int nFixedScanStart, int nFixedScanEnd)
{ 
    for(int nx = nFixedScanStart; nx <= nFixedScanEnd; nx++) 
        m_anFixedScans + nx; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vGetAxisName(int iIndex, Cstring& strAxisName)const
{
    strAxisName = "";
    
    if( !m_poScanInfo )
        return;  // safety
    
    if( iIndex < 0 || iIndex > m_poScanInfo->nGetNumAxes() - 1 )
        return;

    m_poScanInfo->vGetAxisName(iIndex, strAxisName);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cstrategy::bIsRotationAxis(int iIndex)
{
    if( !m_poScanInfo )
        return false;  // safety
    
    if( iIndex < 0 || iIndex > m_poScanInfo->nGetNumAxes() - 1 )
        return false;

    return m_poScanInfo->bIsRotationAxis(iIndex);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vSetExternalResoBinsPtr(CResoStrategyBins* pBins)
{
    if( !m_bExternalResoBinsPtr && m_poResoBins )
    {
        delete m_poResoBins;
        m_poResoBins = NULL;
    }

    if( pBins )
    {
        m_poResoBins = pBins;
        m_bExternalResoBinsPtr = true;
        
        printf("\nSetting strategy resolution range: %f %f A\n", pBins->dGetMinReso(), pBins->dGetMaxReso());
    }
    else
    {
        m_poResoBins = NULL;
        m_bExternalResoBinsPtr = false;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vSetScanInfoObj(CstrategyScans* pScans)
{
    if( m_poScanInfo && !m_bExternalScanInfoPtr )
        delete m_poScanInfo;
    
    if( pScans )
    {
        m_poScanInfo = pScans;
        m_bExternalScanInfoPtr = true;
    }
    else
    {
        m_poScanInfo = new CstrategyScans();
        m_bExternalScanInfoPtr = false;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cstrategy::bSetResoBinsObject(double dStarCubedMin, double dStarCubedMax)
{
    if( m_bExternalResoBinsPtr && !m_poResoBins )
        return false; // inconsistency
    
    if( !m_poResoBins )
        m_poResoBins = new CResoStrategyBins(1.0 / pow(dStarCubedMin, 1.0/3.0), 1.0 / pow(dStarCubedMax, 1.0/3.0), m_nNumResoBins);

    if( !m_poCrystal )
        return false;  // Crystal object must exist

    m_poResoBins->vCountTheoreticalUniqueReflectionsInResoShells(m_poCrystal, m_bAnom);
    
    m_poResoBins->vScaleReflnCount(m_fCellVolumeFactor, DTREK_RESOBIN_SCALE_THEORETICAL_COUNT);

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function will go over the reflection list, calculate "d*cubed" for every reflection and save it in the list.
// It will also find the limits of the d*Cubed within the reflection list.
void Cstrategy::vReflnListCalcSetDStarCubed(double& dStarCubedMin, double& dStarCubedMax)
{
    if( !m_poCrystal || !m_poReflnlist )
        return; // should not happen
    
    int     n2STLcu = m_poReflnlist->nGetFieldIndex(ms_sf2STLcu);
    if( 0 > n2STLcu )
        return;

    double    a3x3fBMat[3][3];      // Local var for crystal B matrix
    double    a3fX[3];              // Local var for recip lattice vector Bh = x
    double    a3fH[3];              // hkl of refln

    /////////////////////////////////////////////////////////////////////
    // Get the crystal B matrix into local variable
    m_poCrystal->nCalcBMatrix();
    m_poCrystal->vGetBMatrix(&a3x3fBMat[0][0]);

  
    dStarCubedMin = 999999.0;  //initialize
    dStarCubedMax = -1.0;      //initialize

    double        dStarCubed = 0.0;  // a temp
    
    Crefln*       poRefln = NULL;

    int         nNumReflns   = m_poReflnlist->nGetNumReflns();
    for(int ii = 0; ii < nNumReflns; ii++)
    {
        poRefln  = m_poReflnlist->poGetRefln(ii);

        // Compute reciprocal lattice vector

        a3fH[0] = (double)poRefln->nGetH();
        a3fH[1] = (double)poRefln->nGetK();
        a3fH[2] = (double)poRefln->nGetL();

        vMulMat3DVec3D(a3x3fBMat, a3fH, a3fX);

        dStarCubed = pow(fDot3D(a3fX, a3fX), 3.0/2.0);   // the dot product gives us d*squared, so we need to raise it to the power of 3/2

        poRefln->vSetField(n2STLcu, (float)dStarCubed);

        if( dStarCubed < dStarCubedMin )
            dStarCubedMin = dStarCubed;

        if( dStarCubed > dStarCubedMax )
            dStarCubedMax = dStarCubed;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reduce reflnlist to asymmetric unit and sort on reduced HKL
void Cstrategy::vReflnListSortReduce()
{
    nExpand(*m_poReflnlist); // just to make sure it has all necessary fields
    m_poReflnlist->nReduce(*m_poCrystal, m_bAnom);  // this will set reduced packed HKLs into CReflnlist::m_nFI_nPackedHKL
                                                    // and sort on that
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Cstrategy::vSetExternalWorkingReflnListPtr(Creflnlist* pList)
{
    if( !pList )
        return; //safety

    m_poReflnlist = pList;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int  Cstrategy::nGetNumAxes()const
{
    if( m_poScanInfo )
        return m_poScanInfo->nGetNumAxes();

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// It is assumed that all ranges in the vector satisfy the minimum completeness requirement as well as the minimum redundancy 
// requirement. We will try to expand every range within the allowed rotation tolerance and see if that improves completeness.
void Cstrategy::vSortFindBestExpandableRange(std::vector<CstrategyRange>& vRanges)
{
    if( 1 > vRanges.size() )
        return; // nothing to do

    int     ii = -1;

    // Sort: key 1: total rotation witdh
    //       key 2: completeness
    //       key 3: redundancy
    // The *shortest* range will have index 0. 
    //std::sort(vRanges.begin(), vRanges.end(), bRangesCompareFindShortest);

    vSortRangeVector(vRanges, bRangesCompareFindShortest);

    // Find the tolerance rotation increment with the accuracy of 1 degree.
    double      dShortestWidth = vRanges[0].dGetTotalRotWidth();
    double      dTolRotIncr = dShortestWidth * m_fRotationTol / 100.0;
    dTolRotIncr = dRoundOff(dTolRotIncr, 0);
    
    Cstring     strScans("");
    Cstring     strScansTemp("");
    double      dRotWidth = 0.0;
    bool        bFirstScanNumberWritten = false;  // just a trick to have commas between numbers, but not after the last number
    for(ii=0; ii < vRanges.size(); ii++ )
    {
        dRotWidth = vRanges[ii].dGetTotalRotWidth();

        if( dRotWidth > dShortestWidth ) // obviously, for index ii=0 that won't be true, because dShortestWidth = vRanges[0].dGetTotalRotWidth();
            break;

        vRanges[ii].vGetScansList(strScansTemp);
        
        if( bFirstScanNumberWritten )
            strScans += ',';

        strScans += strScansTemp;
        
        bFirstScanNumberWritten = true; // reset the flag, so that once we got a scan number, 
                                        // the rest of the scan numbers will be preceded by a comma
    }
    
    printf("\nThe shortest solution width of %.1f degrees was achieved on scan(s) %s.\n", dShortestWidth, strScans.string());
    printf("\nSolutions longer than %.1f degrees will not be considered.\n", dShortestWidth + dTolRotIncr);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Get rid of the ranges that are longer than "the shortest width plus tolerance"
    std::vector<CstrategyRange>::iterator       oIt;
    
    for(oIt = vRanges.begin() + 1; oIt != vRanges.end(); oIt++)
    {
        if( (*oIt).dGetTotalRotWidth() > dShortestWidth + dTolRotIncr )
            break;
    }

    if( oIt != vRanges.end() )
    {
        vRanges.erase(oIt, vRanges.end());
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    // Try and expand all remaining ranges within the tolerance rotation increment
    double      dCurWidth = 0.0;
    for(ii=0; ii < vRanges.size(); ii++)
    {
        dCurWidth = vRanges[ii].dGetTotalRotWidth();
        
        bTryExpandRange(vRanges[ii],  dShortestWidth + dTolRotIncr - dCurWidth);
    }

    // Sort: key 1: completeness
    //       key 2: total rotation witdh
    //       key 3: redundancy
    // The *most complete* range will have index 0. This will be our *best* range. 
    //std::sort(vRanges.begin(), vRanges.end(), bRangesCompareFindMostComplete);
    vSortRangeVector(vRanges, bRangesCompareFindMostComplete);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Cstrategy::bTryExpandRange(CstrategyRange& oInputRange, double dExtraRotValue)
{
    if( m_fCompletenessTol <= 0.0 )
        return false; // nothing to do, because there is no completeness tolerance!

    dExtraRotValue = dRoundOff(dExtraRotValue, 0); // we don't need accuracy higher than 1 degree of rotation
    
    // The *minimum* value by which we will try to increase the total rotation in the input range
    const double    cdRotWidthBit = 1.0; // the *granularity* of the extra rotation      
    if( dExtraRotValue < cdRotWidthBit )
        return false;  // nothing to do if the extra rotation is less than the minimum bit we could add
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Log what we are doing...
    Cstring     strScans("");
    oInputRange.vGetScansList(strScans);
    printf("\nTesting if expanding scan(s) [%s] by up to %.0f"
           " degrees could increase completeness...\n", strScans.string(), dExtraRotValue);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( m_dTotalRotWidthAllowed < oInputRange.dGetTotalRotWidth() )
        return false; // Inconsistency! The total width of the best range should not be larger than the allowed value!

    double      dAllowedIncrease = m_dTotalRotWidthAllowed - oInputRange.dGetTotalRotWidth();  
    dExtraRotValue = min(dExtraRotValue, dAllowedIncrease);
    
    // Now we need to try and distribute the extra rotation width, bit by bit,
    // among the scans. The priority will be given to scans that give most extra total completeness.
    double      dCurLargestCompletenessIncrement = 0.0;
    
    CstrategyRange      oTestRange(this);                       // a temporary range object for testing purposes 
    CstrategyRange      oInputRangeCopy(oInputRange);           // a working copy of the input range object 

    double          dTestRot = cdRotWidthBit;
    double          dRotWidthBank = dExtraRotValue - cdRotWidthBit;

    int             nScanOriginalIndex = -1;
    
    double          dRotBeg = 0.0;
    double          dRotEnd = 0.0;

    struct
    {
        int     nScan;   // to remember the scan number
        bool    bBeg;    // to remember whether that was the beginning or end of the scan
    }stScanWinner;

    CSegment        segScanLimits(-999.0, 999.0);

    // Do it until the "bank" runs out or we achieve the highest desired completeness - whichever comes first.
    bool            bSuccess = false;
    while( dRotWidthBank >= 0.0 && oInputRangeCopy.m_fCompleteness < min(100.0, m_fCompletenessMin + m_fCompletenessTol) )
    {                                                            
        stScanWinner.nScan = -1;
        dCurLargestCompletenessIncrement = 0.0;

        oTestRange = oInputRangeCopy;

        for(int ii=0; ii < oTestRange.nGetNumberOfScans(); ii++)
        {
            nScanOriginalIndex = oTestRange.nGetScanOriginalIndex(ii);

            if( -1 != m_anFixedScans.find(nScanOriginalIndex) )
                continue;

            dRotBeg = oInputRangeCopy.dGetScanRotBegin(ii);
            dRotEnd = oInputRangeCopy.dGetScanRotEnd(ii);

            if( dRotEnd - dRotBeg + dTestRot > m_poScanInfo->dGetScanMaxRotRange(nScanOriginalIndex) )
                continue; // this scan width is already maxed out

            m_poScanInfo->vGetScanRotLimits(nScanOriginalIndex, segScanLimits);

            // If we are ok with the limits, we can try and attach that extra bit and see what happens to the overall completeness. 
            
            // Try the beginning of the scan
            if( !m_bFixRotStart )
            {
                if( dRotBeg - dTestRot >= segScanLimits.dGetBeg() )
                {
                    oTestRange.vSetScanRotBeg(ii, dRotBeg - dTestRot);
                    nCompleteness(oTestRange);
                
                    if( oTestRange.m_fCompleteness - oInputRangeCopy.m_fCompleteness > dCurLargestCompletenessIncrement )   // what if this is negative??? should not be!!!
                    {
                        dCurLargestCompletenessIncrement = oTestRange.m_fCompleteness - oInputRangeCopy.m_fCompleteness;    
                
                        // Remember the winner
                        stScanWinner.nScan = ii;
                        stScanWinner.bBeg = true;
                    }    
                
                    // Restore the test range back to whatever the input range had
                    oTestRange.vSetScanRotBeg(ii, dRotBeg);
                }
            }
            
            // Try the end of the scan
            if( dRotEnd + dTestRot <= segScanLimits.dGetEnd() )
            {
                oTestRange.vSetScanRotEnd(ii, dRotEnd + dTestRot);
                nCompleteness(oTestRange);
                
                if( oTestRange.m_fCompleteness - oInputRangeCopy.m_fCompleteness > dCurLargestCompletenessIncrement )   // what if this is negative??? should not be!!!
                {
                    dCurLargestCompletenessIncrement = oTestRange.m_fCompleteness - oInputRangeCopy.m_fCompleteness;    
                
                    // Remember the winner
                    stScanWinner.nScan = ii;
                    stScanWinner.bBeg = false;
                    
                }
                // Restore the test range back to whatever the input range had
                oTestRange.vSetScanRotEnd(ii, dRotEnd);
            }
        }

        if( stScanWinner.nScan == -1 )
        {
            dTestRot += cdRotWidthBit; // since there are no takers, increase the offer
        }
        else
        {
            // We have a "winner" scan! Change that scan's rot range accordingly.
            bSuccess = true;
            if( stScanWinner.bBeg )
            {
                dRotBeg = oInputRangeCopy.dGetScanRotBegin(stScanWinner.nScan);
                oInputRangeCopy.vSetScanRotBeg(stScanWinner.nScan, dRotBeg - dTestRot);
            }
            else
            {
                dRotEnd = oInputRangeCopy.dGetScanRotEnd(stScanWinner.nScan);
                oInputRangeCopy.vSetScanRotEnd(stScanWinner.nScan, dRotEnd + dTestRot);
            }

            // We need to update the best range completeness, We do not need 
            // to calculate the completeness directly, since we know the previous value
            // and the increment.
            oInputRangeCopy.m_fCompleteness += dCurLargestCompletenessIncrement;

            dTestRot = cdRotWidthBit;
        }
        
        dRotWidthBank -= cdRotWidthBit;
    }
    
    double      dCompletenessGain =  oInputRangeCopy.m_fCompleteness - oInputRange.m_fCompleteness;
    printf("Completeness increased by: %.2f%%\n", bSuccess ? dCompletenessGain : 0.0);
    
    if( bSuccess && dCompletenessGain >= m_dMinCompletenessIncreaseOnExpansion ) 
    {
        printf("Updating scan solution...\n\n");
        oInputRange = oInputRangeCopy;
        
        if( bIsMultipleScanSolution() )
        {
            nCompletenessMultiScan(oInputRange); // NOTE: This call automatically updates Strategy Scans Array, 
                                                  // so we do not need to call nFromRangeToStrategyScansArray() here.
        }
        else
        {
            nCompleteness(oInputRange);
            nFromRangeToStrategyScansArray(oInputRange, DTREK_DTSTRATEGY_FRTSSA_UPDATE_ALL_RESULTS);    
        }
    }
    else
    {
        bSuccess = false;
        printf("Expanding scan(s) was not effective.\n\n");
    }

    return bSuccess;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function returns true if the first range is better.
static bool bRangesCompareConsiderTargets(CstrategyRange& oFirstRange, 
                                          CstrategyRange& oSecondRange, 
                                          double dTargetCompl, 
                                          double dTargetRed,
                                          bool& bIsBetterCompleteness)
{
    // Try to find the better range based on the completeness.
    // If both are the same or both are complete ABOVE target, no distinction can be made.
    bIsBetterCompleteness = false;
    if( oFirstRange.m_fCompleteness != oSecondRange.m_fCompleteness &&
        !(oFirstRange.m_fCompleteness >= dTargetCompl && oSecondRange.m_fCompleteness >= dTargetCompl) )
    {
        bIsBetterCompleteness = true;
        return (oFirstRange.m_fCompleteness > oSecondRange.m_fCompleteness);
    }

    // Try to find the better range based on the redundancy
    // If both are the same or both are redundant ABOVE target, no distinction can be made.
    if( oFirstRange.m_fRedundancy != oSecondRange.m_fRedundancy &&
        !(oFirstRange.m_fRedundancy >= dTargetRed && oSecondRange.m_fRedundancy >= dTargetRed) )
    {
        return (oFirstRange.m_fRedundancy > oSecondRange.m_fRedundancy);
    }

    // The better range is the one that is shorter
    return (oFirstRange.dGetTotalRotWidth() < oSecondRange.dGetTotalRotWidth());
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static bool bRangesCompareFindShortest(CstrategyRange& oFirstRange, CstrategyRange& oSecondRange)
{
    double      dRotWidth_1 = oFirstRange.dGetTotalRotWidth();
    double      dRotWidth_2 = oSecondRange.dGetTotalRotWidth();

    if( dRotWidth_1 == dRotWidth_2 )
    {
        if( oFirstRange.m_fCompleteness == oSecondRange.m_fCompleteness )
            return (oFirstRange.m_fRedundancy > oSecondRange.m_fRedundancy);    
        else 
            return (oFirstRange.m_fCompleteness > oSecondRange.m_fCompleteness);
    }

    return (dRotWidth_1 < dRotWidth_2);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static bool bRangesCompareFindMostComplete(CstrategyRange& oFirstRange, CstrategyRange& oSecondRange)
{
    if( oFirstRange.m_fCompleteness == oSecondRange.m_fCompleteness )
    {
        double      dRotWidth_1 = oFirstRange.dGetTotalRotWidth();
        double      dRotWidth_2 = oSecondRange.dGetTotalRotWidth();
        
        if( dRotWidth_1 == dRotWidth_2 )    
            return (oFirstRange.m_fRedundancy > oSecondRange.m_fRedundancy);    
        else
            return dRotWidth_1 < dRotWidth_2;
    }
    else 
        return (oFirstRange.m_fCompleteness > oSecondRange.m_fCompleteness);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Bubble sort of ranges
static void vSortRangeVector(std::vector<CstrategyRange>& vRanges, bool (*pfbCompare)(CstrategyRange&, CstrategyRange&))
{
    CstrategyRange      oTempRange(NULL);
    
    for(int ii = 0; ii < vRanges.size() - 1; ii++)
    {
        for(int jj = 0; jj < vRanges.size() - 1 - ii; jj++)
        {
            if( pfbCompare(vRanges[jj+1], vRanges[jj]) ) // compare the two neighbors
            {  
                // swap vRanges[jj] and vRanges[jj+1]
                oTempRange = vRanges[jj];         
                vRanges[jj] = vRanges[jj+1];
                vRanges[jj+1] = oTempRange;
            }
        }
    }    
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CstrategyRange::dGetTotalRotWidth()
{ 
    double    dd = 0.0;
    int       nScanNumber = -1;
    for(int ii=0; ii < m_afRotStart.size(); ii++)
    {
        nScanNumber = m_anScanNumbers[ii];

        if( m_poStrategy->bIsFixedScan(nScanNumber) )
            continue;

        dd += fabs(m_afRotStart[ii] - m_afRotEnd[ii]); 
    }

    return dd;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CstrategyRange::vGetScansList(Cstring& strScans, bool bIncludeFixedScans)
{
    strScans = "";
    Cstring     strIndex("");
    int         iScan = -1;
    int         nScan = -1;

    bool        bFirstScanNumberWritten = false;  // just a trick to have commas between numbers, but not after the last number

    while( (nScan=nGetScanOriginalIndex(++iScan)) != -1 )
    {
        if( !bIncludeFixedScans && m_poStrategy->bIsFixedScan(nScan) )
            continue;

        strIndex = nScan;

        if( bFirstScanNumberWritten )
            strScans += ',';
        
        strScans += strIndex;

        bFirstScanNumberWritten = true; // reset the flag, so that once we got a scan number, 
                                        // the rest of the scan numbers will be preceded by a comma
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////



    

    
