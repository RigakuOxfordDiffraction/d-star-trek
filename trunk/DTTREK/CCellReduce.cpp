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
// CCellReduce.cc            Initial author: Thaddeus J Niemeyer
//  This file contains the member functions of class CCellReduce.
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

#include "CCellReduce.h"
#ifdef SSI_PC
#include "CrclHelper.h"
#endif
#include <string.h>

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
using std::swap;
using std::setw;
#endif

//#if defined(VC6)
//#define 'std::setw' 'setw'
//#endif

CCellReduce::CCellReduce() {
    m_bKnownCell        = FALSE;
    m_nVerbose          = 2;
    m_poSolutions       = NULL;
    m_fResidualMax      = 20.0;
    m_fWavelength       = 1.0;
    m_nLattNum          = -1;
    vInitSolutions();
    vInitLattChar();
};

CCellReduce::~CCellReduce() {
    if (m_poSolutions)
        delete m_poSolutions;
    m_poSolutions = NULL;
};

Cstring CCellReduce::ms_sfA      = "fA";
Cstring CCellReduce::ms_sfB      = "fB";
Cstring CCellReduce::ms_sfC      = "fC";
Cstring CCellReduce::ms_sfAlpha  = "fAlpha";
Cstring CCellReduce::ms_sfBeta   = "fBeta";
Cstring CCellReduce::ms_sfGamma  = "fGamma";
Cstring CCellReduce::ms_sfVolume = "fVolume";
Cstring CCellReduce::ms_sfRot1   = "fRot1";
Cstring CCellReduce::ms_sfRot2   = "fRot2";
Cstring CCellReduce::ms_sfRot3   = "fRot3";
Cstring CCellReduce::ms_sfDiffLength = "fDiffLength";
Cstring CCellReduce::ms_ssLattice= "sLattice";
Cstring CCellReduce::ms_ssLattSymm= "sLattSymm";
Cstring CCellReduce::ms_snLattNum= "nLattNum";
Cstring CCellReduce::ms_snSortOrder= "nSortOrder";



void
CCellReduce::vInitSolutions(void)
{
    // Initialize the m_poSolutions object
    
    if (NULL == m_poSolutions) {
        m_poSolutions  = new Creflnlist ();
    };
    m_poSolutions->vDeleteAll();
    
    m_nFI_fA       = m_poSolutions->nExpandGetField(ms_sfA);
    m_nFI_fB       = m_poSolutions->nExpandGetField(ms_sfB);
    m_nFI_fC       = m_poSolutions->nExpandGetField(ms_sfC);
    m_nFI_fAlpha   = m_poSolutions->nExpandGetField(ms_sfAlpha);
    m_nFI_fBeta    = m_poSolutions->nExpandGetField(ms_sfBeta);
    m_nFI_fGamma   = m_poSolutions->nExpandGetField(ms_sfGamma);
    m_nFI_fVolume  = m_poSolutions->nExpandGetField(ms_sfVolume);
    m_nFI_fRot1    = m_poSolutions->nExpandGetField(ms_sfRot1);
    m_nFI_fRot2    = m_poSolutions->nExpandGetField(ms_sfRot2);
    m_nFI_fRot3    = m_poSolutions->nExpandGetField(ms_sfRot3);
    m_nFI_sLattice = m_poSolutions->nExpandGetField(ms_ssLattice);
    m_nFI_sLattSymm= m_poSolutions->nExpandGetField(ms_ssLattSymm);
    m_nFI_nLattNum = m_poSolutions->nExpandGetField(ms_snLattNum);
    m_nFI_nSortOrder = m_poSolutions->nExpandGetField(ms_snSortOrder);
}


void CCellReduce::vCalcGetG6(const double *pfCell, double *pfG6)
{
  // Calculate and get (return) the 6 dimensional
  // vector G6:
  //
  //  a dot a
  //  -     -
  //
  //  b dot b
  //  -     -
  //
  //  c dot c
  //  -     -
  //
  //  b dot c * 2
  //  -     -
  //
  //  a dot c * 2
  //  -     -
  //
  //  a dot b * 2
  //  -     -
  //
  int i, j, k;
  double fCos;
  for (i = 0; i < 3; i++)
    {
      pfG6[i]   = pfCell[i] * pfCell[i];

      j = (i+1) % 3;             // j = 1, 2, 0
      k = (i+2) % 3;             // k = 2, 0, 1

      fCos = cos(pfCell[i+3] * Gs_dRADIANS_PER_DEGREE);
      if (fabs(fCos) <= 0.000001) fCos = 0.0;
      pfG6[i+3] = 2.0f * pfCell[j] * pfCell[k] * fCos;
    }
}

void
CCellReduce::vCalcGetG6(const double a3x3fOrientMat[3][3], double *pfG6)
{
  // Calculate and get (return) the 6 dimensional
  // vector G6:
  //
  //  a dot a
  //  -     -
  //
  //  b dot b
  //  -     -
  //
  //  c dot c
  //  -     -
  //
  //  b dot c * 2
  //  -     -
  //
  //  a dot c * 2
  //  -     -
  //
  //  a dot b * 2
  //  -     -
  //

  pfG6[0] = fDot3D(a3x3fOrientMat[0], a3x3fOrientMat[0]);
  pfG6[1] = fDot3D(a3x3fOrientMat[1], a3x3fOrientMat[1]);
  pfG6[2] = fDot3D(a3x3fOrientMat[2], a3x3fOrientMat[2]);
  pfG6[3] = 2.0f * fDot3D(a3x3fOrientMat[1], a3x3fOrientMat[2]);
  pfG6[4] = 2.0f * fDot3D(a3x3fOrientMat[0], a3x3fOrientMat[2]);
  pfG6[5] = 2.0f * fDot3D(a3x3fOrientMat[0], a3x3fOrientMat[1]);
}
////////////////////////////////////////////////////////////////

int CCellReduce::nStandardPres(double *pfCell)
{

  double fG[6], fGprime[6];
  double fMat6[6][6];
  bool  bChanged;
  int   i;
  
  vCalcGetG6(&pfCell[0], fG);    

  // Reduce a crystallographic unit cell to a standard
  // reduced cell
    do
    {
        // Convert to standard presentation
        
        bChanged = FALSE;
        vZeroMat(6, 6, &fMat6[0][0]);
        if (    (fG[0] - fG[1] > 0.0001 )
            || ((fG[0] == fG[1]) && (fabs(fG[3]) > fabs(fG[4])) ) )
        {
            // SP1: Interchange a and b axes
            
            fMat6[0][1] = 1.0;
            fMat6[1][0] = 1.0;
            fMat6[2][2] = 1.0;
            fMat6[3][4] = 1.0;
            fMat6[4][3] = 1.0;
            fMat6[5][5] = 1.0;
            //              cout << "reduce: switch a and b\n";
            bChanged    = TRUE;
        }
        else if (    (fG[1] - fG[2] > 0.0001)
            || ((fG[1] == fG[2]) && (fabs(fG[4]) > fabs(fG[5])) ) )
        {
            // SP2: Interchange b and c axes
            
            fMat6[0][0] = 1.0;
            fMat6[1][2] = 1.0;
            fMat6[2][1] = 1.0;
            fMat6[3][3] = 1.0;
            fMat6[4][5] = 1.0;
            fMat6[5][4] = 1.0;
            //              cout << "reduce: switch b and c\n";
            bChanged    = TRUE;
        }
        else if (   (fG[3] * fG[4] * fG[5] > 0.0)
            && ( (fG[3] < 0.0) || (fG[4] < 0.0) || (fG[5] < 0.0) ) )
        {
            // SP3: Set them all positive (make angles <= 90)
            
            fMat6[0][0] = 1.0;
            fMat6[1][1] = 1.0;
            fMat6[2][2] = 1.0;
            fMat6[3][3] = fSign(fG[3]);
            fMat6[4][4] = fSign(fG[4]);
            fMat6[5][5] = fSign(fG[5]);
            //              cout << "reduce: make all angles <= 90\n";   // m = -m?
            bChanged = TRUE;
        }
        else if (     (fG[3] * fG[4] * fG[5] <= 0.0)
            && ( (fG[3] > 0.0) || (fG[4] > 0.0) || (fG[5] > 0.0) ) )
        {
            // SP4: Set them all negative
            
            fMat6[0][0] =  1.0;
            fMat6[1][1] =  1.0;
            fMat6[2][2] =  1.0;
            fMat6[3][3] = -fSign(fG[3]);
            fMat6[4][4] = -fSign(fG[4]);
            fMat6[5][5] = -fSign(fG[5]);
            //              cout << "reduce: make all angles > 90\n";   // n = -n?
            bChanged = TRUE;
        }
        if (bChanged)
        {
            vMulMatNDVecND(6,  &fMat6[0][0], &fG[0], &fGprime[0]);
            for (i = 0; i < 6; i++) fG[i] = fGprime[i];
        }
    } while (bChanged);
    
    // End of conversion to standard presentation
    
    nCalcGetCellFromG6(&fG[0], &pfCell[0]);
    
    return 0;
}


int CCellReduce::nReduceCell(double *pfCell)
{

  // DON'T See Andrews & Bernstein (1988) Acta Cryst. A44, 1009-1018

    /*
    This matrix records 24 possible reductions. Do not confuse these with the 24 Buerger transforms!
    The integers 1,2,3 signify a,b, and c axes respectively.

    Col [0] =           Axis to be modified.
    Col [1] Col [3]  =  [3]* Col[1]  are added to Col [0]
    Col [2] Col [4]  =  [4]* Col[2]  are added to Col [0]

    We loop over all 24 of these operations.  When we can no longer reduce the cell, we stop.
    A cell reduction occurs when no edge get's shorter under any of the 24 operations.  Each edge
    has 8 possible ways of getting reduced.

    This cell works in non-oriented 3 dimensional space (9 quantities)

    */


   int nLCPre[][5]={
       {0,1,0,  1,0},
       {0,1,0,  -1,0},
       {1,0,0,  1,0},
       {1,0,0,  -1,0},
       {0,2,0,  1,0},
       {0,2,0,  -1,0},
       {2,0,0,  1,0},
       {2,0,0,  -1,0},
       {1,2,0,  1,0},
       {1,2,0,  -1,0},
       {2,1,0,  1,0},
       {2,1,0,  -1,0},
       {0,1,2,  1,1},
       {0,1,2,  1,-1},
       {0,1,2,  -1,1},
       {0,1,2,  -1,-1},
       {1,0,2,  1,1},
       {1,0,2,  -1,1},
       {1,0,2,  1,-1},
       {1,0,2,  -1,-1},
       {2,0,1,  1,1},
       {2,0,1,  -1,1},
       {2,0,1,  1,-1},
       {2,0,1,  -1,-1}
        };
   int nReductionLoopCt=0;
   bool bReductionmade;
   Ccrystal oCrystal;

   double a3fVecT[3];      // Temporary vector
   double a3fVecT2[3];     // Temporary vector.
   double a3x3fCellMat[3][3];
   int nx;
   double f0,f1;

   // Build the cell in a3x3fCellMat
   oCrystal.vSetCell(pfCell);
   oCrystal.nCalcGetRealMatrix(&(a3x3fCellMat[0][0]));
   vTranMat3D(a3x3fCellMat);
   
   do
     {
       bReductionmade = FALSE;
       for (nx = 0; nx < 24; nx++)
         {
           do 
             {
               f0 = fabs(
                      fDot3D(a3x3fCellMat[nLCPre[nx][0]],
                             a3x3fCellMat[nLCPre[nx][0]]));
      
               // [1]/[3]
               vCopyVec3D(&a3x3fCellMat[nLCPre[nx][1]][0], a3fVecT);
               vMulVec3DScalar(a3fVecT, nLCPre[nx][3], a3fVecT);
               vAddVec3DVec3D(a3fVecT, a3x3fCellMat[nLCPre[nx][0]], a3fVecT);
               // [2]/[4]
               vCopyVec3D(&a3x3fCellMat[nLCPre[nx][2]][0], a3fVecT2);
               vMulVec3DScalar(a3fVecT2, nLCPre[nx][4], a3fVecT2);
               vAddVec3DVec3D(a3fVecT, a3fVecT2, a3fVecT);
      
               f1 = fabs(fDot3D(a3fVecT, a3fVecT));

               // Did the edge get shorter?

               if ( f0 > f1 ) 
                 { 
                   vCopyVec3D(a3fVecT, &a3x3fCellMat[nLCPre[nx][0]][0]); 
                   bReductionmade = TRUE;
                   nReductionLoopCt++; 
                   if (500 < nReductionLoopCt)
                     cout << "f0, f1: " << f0 << ", " << f1 << endl << flush;
                 }
             } while (f0 > f1);      
         }
       // Loop again if any reduction was made.
     } while ( bReductionmade && (5000 > nReductionLoopCt) );
   
   if (5000 <= nReductionLoopCt)
     {
       cout << "WARNING in CCellReduce::nReductionLoopCt: " 
            << nReductionLoopCt << endl << flush;
     }

   // Copy the cell over.
   pfCell[0] = fLenVec3D(a3x3fCellMat[0]);
   pfCell[1] = fLenVec3D(a3x3fCellMat[1]);
   pfCell[2] = fLenVec3D(a3x3fCellMat[2]);
   pfCell[3] = acos(fDot3D(a3x3fCellMat[1], a3x3fCellMat[2])
                    /pfCell[1] / pfCell[2])
               / Gs_dRADIANS_PER_DEGREE;
   pfCell[4] = acos(fDot3D(a3x3fCellMat[0], a3x3fCellMat[2])
                    / pfCell[0] / pfCell[2])
               / Gs_dRADIANS_PER_DEGREE;
   pfCell[5] = acos(fDot3D(a3x3fCellMat[0], a3x3fCellMat[1])
                    / pfCell[0] / pfCell[1]) 
              / Gs_dRADIANS_PER_DEGREE;  
  return (0);
}

int CCellReduce::nCalcGetCellFromG6(const double *pfG, double *pfCell)
{
  // Calculate the reduced primitive cell from the G6 vector.
  // Also calculate any orthogonal rotation applied to the cell

  int   i, j, k;
  double fTemp;
  double fTemp2;

  // Convert *pfG back to *pfCell

  if ( (pfG[0] <= 0.0f) || (pfG[1] <= 0.0f) || (pfG[2] <= 0.0f) )
    {
      // This is sort of impossible: a cell with non-positive cell lengths

      pfCell[0] = 0.1f;     // Return some wacko cell that could not reflect
      pfCell[1] = 0.1f;     // reality.
      pfCell[2] = 0.1f;
      pfCell[3] = 75.0;
      pfCell[4] = 80.0;
      pfCell[5] = 85.0;
      return (-1);
    }

  for (i = 0; i < 3; i++) pfCell[i] = sqrt(pfG[i]);

  for (i = 0; i < 3; i++)
    {
      j = (i+1) % 3;             // j = 1, 2, 0
      k = (i+2) % 3;             // k = 2, 0, 1
      fTemp2    = pfCell[j] * pfCell[k];
      if (0.0f != fTemp2)
        {
          fTemp     = pfG[i+3] * 0.5f / fTemp2;
          if (-1.0 > fTemp) fTemp = -1.0;
          if ( 1.0 < fTemp) fTemp =  1.0;
          pfCell[i+3] = acos(fTemp) / Gs_dRADIANS_PER_DEGREE;
        }
      else
        {
          pfCell[0] = 0.1f;     // Return some wacko cell that could not reflect
          pfCell[1] = 0.1f;     // reality.
          pfCell[2] = 0.1f;
          pfCell[3] = 75.0;
          pfCell[4] = 80.0;
          pfCell[5] = 85.0;
          return (-1);
        }
    }
  return (0);
}

int
CCellReduce::nCalcGetCellFromG6(const double *pfG, double a3x3fOrientMat[3][3])
{
  // Calculate the reduced primitive cell from the G6 vector.
  // Also calculate any orthogonal rotation applied to the cell

  //int   i, j, k;
  //double fTemp;

  // Convert *pfG back to orientation matrix

  if ( (pfG[0] <= 0.0) || (pfG[1] <= 0.0) || (pfG[2] <= 0.0) )
    {
      // This is sort of impossible: a cell with non-positive cell lengths

      vIdentMat3D(a3x3fOrientMat);
      return (-1);
    }

  //

  return (0);
}

int CCellReduce::ms_anOrder[] = {44,  3,  5,  1,                                            // cubic
                                                                                            
                                  2,  4,  9, 24,                                            // rhombohedral
                                                                                            
                                  12, 22,                                                   // hexagonal
                                                                                            
                                  6,  7, 15, 18, 11, 21,                                    // tetragonal
                                  
                                  8, 19, 42, 13, 23, 36, 38, 40, 16, 26, 32,                // orthorombic
                                 
                                 10, 14, 17, 20, 25, 27, 28, 29, 30, 37, 39, 41, 43,        // monoclinic mC
                                                             10, 14, 17, 20, 25, 27, 28, 29, 30, 37, 39, 41, 43,        // monoclinic mC
                                 
                                 33, 34, 35,                                                // monoclinic mP
                                                             33, 34, 35,                                                // monoclinic mP
                                                             
                                 31, 44, -1};                                               // triclinic

int CCellReduce::ms_anOrderSlot[] = {0, 0, 0, 0, 
                                     
                                     0, 0, 0, 0,
                                     
                                     0, 0,

                                     0, 0, 0, 0, 0, 0,
                                     
                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     
                                     0, 0, 0,
                                                                 1, 1, 1,

                                                                 0, 0, 0};


int CCellReduce::nCalcGetBestLattices()
{
    // Try to compute a residual of a chosen cell when fit to a
    // Bravais lattice.
    
    // This is based on
    // Anderson & Bernstein (1988) Acta Cryst. (1988) A44, 1009-1018 and
    // Paciorek & Bonin (1992) J. Appl. Cryst. 25, 632-637.
    // International Tables Vol A, Table 9.3.*, 746-747.
    
    // Many matrices are initialized in ::vInitLattChar()

    int     ii = 0;
    int     nx = 0;
    int     nStat;
    int     nNumStart,nEnd;
    int     nLoop;
    int     nBravisLattice;
    int         nLatt;
        int             nSlot;

    double   f0;
    double   a6fG[6], a6fGprime[6], a6fGdiff[6], fError;
    double   a6fGorig[6];
    double   a6fCell[6];
    double   fPrevVolume;
    double   fLenG;
    double   fLowest;          // Set to the lowest fError  in the nLatt loop.  
    
    double   fLowestInLattice = 1e20;
    
    Cstring sLatticeType;
    Cstring sLattSymm;
    Crefln  oRefln(m_poSolutions);
    double   fRot1, fRot2, fRot3;
    char*   pcLattString;
    
    
    double*  pfReducedCells=NULL;   // We reduce each cell OUTSIDE the main loop, and store it here.  Thus nReduceCell() is called only once per cell.
    double*  pfCellIntensity=NULL;
    
    Ccrystal oCrystal;


    
    // (1.0-fErrorLimit) gives % decrease in fError 
    // that excludes a cell from edge length, 
    // and sum of angles critereons.
    // In other words, if a cell has a (1.0-fErrorLimit)% 
    //  decrease in error, then we consider this 
    // alone to be enough justification
    // for keeping the cell.

    double fErrorLimit = 0.20;  
    double fErrorLimitTol = 0.001;    
    double fErrorLimitRot= 0.02;     
    double fErrorLimitSumSin=0.02;
    double fMonoclinicDegreeLimit = 1.0;
    double a2fVolumeFractionLimits[2] = {0.8,1.2};

    
    // The solution loaded should already be a non-centered primitive cell.
    // It need not be a reduced primitive cell.
    // It usually comes from nGetUBMatrix or nGetUBMatrix2.
    
    
    int nSol, anSol[DT_NUMSOL];
    
    nNumStart = m_poSolutions->nGetNumReflns();   // Save starting size of m_poSol
    if (0 >= nNumStart)
    {
        cout << "ERROR, no triple basis vectors for consideration in cell reduction!\n"
            << flush;
        return (1);
    }
    
    cout << "Least square fit to lattice characters...see"
        << "\n  Andrews & Bernstein (1988) Acta Cryst. A44, 1009-1018 and"
        << "\n  Paciorek & Bonin (1992) J. Appl. Cryst. 25, 632-637.\n" << flush;
    
    // Save indices of best DT_NUMSOL basis vectors found in m_poSolutions
    
    nSol = nNumStart;
    int    jj = min(DT_NUMSOL, nSol);
    
    for(ii = 0; ii < DT_NUMSOL; ii++)
    {
        // At first, set all to -1, indicating no solution available
        
        anSol[ii]               = -1;
        m_a15x3nBIndices[ii][0] = -1;
    }
    
    for (ii = 1; ii <= 44; ii++)
        m_a44fLattResid[ii]=1e20;
    
    pfReducedCells  = new double[6 * jj];
    pfCellIntensity = new double[jj];
    
    for (ii = 0; ii < jj; ii++)
    {
        anSol[ii] = ii;
        nStat    = nLoadSolution(anSol[ii]);
        if (0 != nStat)
        {
            // A flag indicating the cell did not load properly.
            
            anSol[ii]       = -1;
            
            // Basis cell did not load correctly, so continue from this loop.
            continue;
        }
        
        m_a15x3nBIndices[ii][0] = m_poSolutions->poGetRefln(anSol[ii])->nGetH();
        m_a15x3nBIndices[ii][1] = m_poSolutions->poGetRefln(anSol[ii])->nGetK();
        m_a15x3nBIndices[ii][2] = m_poSolutions->poGetRefln(anSol[ii])->nGetL();
        pfCellIntensity[ii] = m_poSolutions->poGetRefln(anSol[ii])->fGetIntensity();      
        m_oCrystal.vGetCell(&a6fCell[0]);
        
        if (m_bKnownCell && ('P' == m_oCrystal.m_poSpacegroup->cGetCentFlag()) )
        {
            // cout << "DID NOT call nReduceCell in nCalcGetBestLattices!\n" << flush;
        }
        else
        {
            // cout << "DID CALL nReduceCell in nCalcGetBestLattices!\n" << flush;             
            nStat = nReduceCell(a6fCell);
        }
        
        //      cout << "indices in nBest: " << m_nIndex0 << ", " << m_nIndex1 << ", " << m_nIndex2 << endl;
        if (5 <= m_nVerbose)
        {
            // List starting cell
            cout << "\nStarting cell in CalcBestLattices is: \n";
            (void) m_oCrystal.nList(1);
        }
        
        // Load the cell into the memory buffer.
        
        for (nx = 0; nx < 6; nx++) 
            pfReducedCells[ii * 6 + nx] = a6fCell [nx];
    }
    
    // For each of the 44 lattice characters
    //   for each of the selected basis vectors
    //     for each of the nearly reduced Buerger cells of the basis vectors
    //        get the least squares fit to that lattice character
    
    // cP  cI cF hR hP tI tP oI oC oF oP mC mP aP //
    
    // TJN: 1/13/03:
    // To solve some of the problems we have been having with cell reduction, the algorithm
    // has been modified to loop through twice.  The first loop behaves like normal:  Each lattice
    // is given it's wack at the input crystal, and saves a solution, resulting in 14 output solutions.
    // In the second loop however, each lattice is given the *best* cell for that lattice computed in the
    // first pass, and each transformation is again allowed to convert the crystal.  All correct solutions will of course
    // have an (almost) 0.0 residual.  It is hoped that this second loop will remove the interplay between projection error
    // and beta angles that has plagued cell reduction.  In the second loop, projection error is no longer an issue,
    // allowing each lattice to choose the most orthodox cell.  (In monoclinic for example, the beta angle closest to 90.0)
    
    bool    bSpecLatType = false;
    
    for(nLoop = 0; nLoop < 2; nLoop++)
    {
        oRefln.vSetIntensity(1000.0f);  // Set high for first of this lattchar
        
        fLowest = 1000.0;
        
        nBravisLattice = -1;
        
        pcLattString = "";
        
        for(ii = 1; ms_anOrder[ii] != -1; ii++)
        {
                        // what we want to do is create new array to correspond to ms_anOrder[].
                        // i will therefore go further than 44.
                        // for mP and mC cells, we will be storing "best" and "second best" entries.  The "second best" must be greater than the
                        // error tollerance, so it's much more restrictive.
            
            nLatt = ms_anOrder[ii];
            nSlot = ms_anOrderSlot[ii];
            
            if( m_a44tLattChar[nLatt].sLatticeType != pcLattString || ii == 1 || ms_anOrderSlot[ii] != ms_anOrderSlot[ii-1] )
            {
                // This is a new type of Bravis lattice.  Reset the errors.
                oRefln.vSetIntensity(1000.0f);  // Set high for first of this lattchar
                fLowest = 1000.0;
                nBravisLattice++;

                // Compute the lowest residual in the lattice ahead of time.
                // This is important during the second pass for "slot2" lattices.  The first slot is reserved for the
                // best lattice fit, whereas, the second slot is reserved for all others.
                                
                fLowestInLattice = 1e20;
                for( jj = 1; jj <= 44; jj++) 
                {
                    if( m_a44tLattChar[jj].sLatticeType == m_a44tLattChar[nLatt].sLatticeType )
                        fLowestInLattice = min(fLowestInLattice, m_a44fLattResid[jj]);
                }  
            }
            
            pcLattString = m_a44tLattChar[nLatt].sLatticeType.string();
            
            bool bKeepCurrent = FALSE;
            
            printf(".");
            
            fError = 100.0;
            
            // Keep looping over different basis vectors only until a low residual is found
            for(nSol = 0; nSol < DT_NUMSOL  && anSol[nSol] >= 0; nSol++)
            {
                // Check top DT_NUMSOL solutions in reverse order going from worst to best,
                // so that the best get to override the worst
                if (0 > anSol[nSol]) 
                {
                    cout << "Doing a break in ::nCalcGetBestLattices().  Please contact author.\n";
                    break;
                    
                }
                
                if (nLoop == 0) 
                {
                    // Load the reduced cell into a6fCell            
                    for (nx = 0; nx < 6; nx++ )
                    {
                        a6fCell[nx] = pfReducedCells[ nSol*6 + nx ];
                    }
                    fPrevVolume = 0.0;
                } 
                else 
                {
                    // Load cell from previous
                    a6fCell[0] = (*m_poSolutions)[nNumStart + nBravisLattice].fGetField(m_nFI_fA);
                    a6fCell[1] = (*m_poSolutions)[nNumStart + nBravisLattice].fGetField(m_nFI_fB);
                    a6fCell[2] = (*m_poSolutions)[nNumStart + nBravisLattice].fGetField(m_nFI_fC);
                    a6fCell[3] = (*m_poSolutions)[nNumStart + nBravisLattice].fGetField(m_nFI_fAlpha);
                    a6fCell[4] = (*m_poSolutions)[nNumStart + nBravisLattice].fGetField(m_nFI_fBeta);
                    a6fCell[5] = (*m_poSolutions)[nNumStart + nBravisLattice].fGetField(m_nFI_fGamma);
                    
                    nCentToPrim(pcLattString[1], a6fCell, false);
                    
                    nReduceCell(a6fCell);
                    
                    fPrevVolume = (*m_poSolutions)[nNumStart + nBravisLattice].fGetField(m_nFI_fVolume);
                }
                
                // Even though unit cell is reduced, we need to loop over the
                // nearly reduced Buerger cells, too
                
                double      a6fCellorig[6];              
                vCopyVecND(6, a6fCell, a6fCellorig);             

                for(jj = 0; jj < 25; jj++)
                {    
                    /*  We want to include all possible cell orientations.  
                    I am using the fact that vCalcGetG6() will properly compute six
                    dimensional vectors EVEN WHEN SOME OF THE CELL LENGTHS ARE NEGATIVE!
                
                     Since the cell is obtained by calling nCalcGetCellFromG6(), the 
                     cell in a6fCell should not matter.
                     */
                     for(int nPermCt = 0; nPermCt < 9; nPermCt++)
                     {
                         vCopyVecND(6, a6fCellorig, a6fCell);
                         switch (nPermCt)
                         {
                         case 0:
                             break;
                         case 1: 
                             std::swap(a6fCell[0],a6fCell[1]); 
                             std::swap(a6fCell[3],a6fCell[4]); 
                             std::swap(a6fCell[0],a6fCell[2]); 
                             std::swap(a6fCell[3],a6fCell[5]); 
                             break;
                         case 2: 
                             std::swap(a6fCell[0],a6fCell[2]); 
                             std::swap(a6fCell[3],a6fCell[5]); 
                             std::swap(a6fCell[0],a6fCell[1]); 
                             std::swap(a6fCell[3],a6fCell[4]); 
                             break;
                         case 3: 
                             a6fCell[0] *= -1; 
                             break;
                         case 4: 
                             a6fCell[1] *= -1; 
                             break;
                         case 5: 
                             a6fCell[2] *= -1; 
                             break;
                         case 6: 
                             a6fCell[0] *= -1; 
                             a6fCell[1] *= -1; 
                             break;
                         case 7: 
                             a6fCell[0] *= -1; 
                             a6fCell[2] *= -1; 
                             break;
                         case 8: 
                             a6fCell[1]*=-1; 
                             a6fCell[2]*=-1; 
                             break;
                         } // end case
                     
                     
                         vCalcGetG6(a6fCell, a6fGorig);
                         fLenG  = fLenVecND(6, a6fGorig);
                     
                         // Convert G6orig to the 24 nearly Buerger-reduced cells.
                         // j = 0, j=25, are the identity, while j = 1...24 are the operations
                         // in Table 3 of Andrews and Berstein.
                         // Multiply original G6 to get G6 to test
                     
                         vMulMatNDVecND(6,  &m_a26x6x6fBMat[jj][0][0], &a6fGorig[0], &a6fG[0]);
                     
                         /*  Beware:  For One of the elements a.b, b.c or a.c elements close to zero, ignore
                         the type I and type II constraints.
                         */
                     
                         bool bPositive = (a6fG[3] * a6fG[4] * a6fG[5] > 0.0f);
                     
                         fRot1 =   sqrt(fabs(a6fG[0]))
                             + sqrt(fabs(a6fG[1])) 
                             + sqrt(fabs(a6fG[2]));
                     
                     
                         if ((m_a44tLattChar[nLatt].nTypeDependent) && 
                             (ABS(a6fCell[3] - 90.0)>2.0) &&
                             (ABS(a6fCell[4] - 90.0)>2.0) &&
                             (ABS(a6fCell[5] - 90.0)>2.0))
                         {
                         
                         
                             if( m_a44tLattChar[nLatt].nType == 1 && !bPositive ) 
                                 continue;
                             if( m_a44tLattChar[nLatt].nType == 2 &&  bPositive ) 
                                 continue;
                         }
                     
                         // The least squares fit is a matrix multiplication
                     
                         vMulMatNDVecND(6, &m_a44tLattChar[nLatt].a6x6fLSMat[0][0], &a6fG[0], &a6fGprime[0]);
                         vSubVecNDVecND(6, a6fGprime, a6fG, a6fGdiff);
                         fError = 100.0f * fLenVecND(6, a6fGdiff) / fLenG; // Error as %                   
                     
                         if (8 <= m_nVerbose)
                         {
                             cout << ii << ": " << sLatticeType << ": "
				  << setw(9) << a6fGprime[0] << ", "
				  << setw(9) << a6fGprime[1] << ", "
				  << setw(9) << a6fGprime[2] << ", "
				  << setw(9) << a6fGprime[3] << ", "
				  << setw(9) << a6fGprime[4] << ", "
				  << setw(9) << a6fGprime[5] << endl;
                         }
                     
                         // Convert cell from reduced primitive to actual cell 
                         // (possibly a centered cell), so we can skip cells with
                         // strange/impossible angles and for monoclinic cells those
                         // with beta angle less than about 90 deg.
                     
                         // Calculate |a| + |b| + |c| of the reduced primitive cell
                         // and save in fRot1 for now
                     
                         fRot1 =   sqrt(fabs(a6fGprime[0]))
                             + sqrt(fabs(a6fGprime[1])) + sqrt(fabs(a6fGprime[2]));
                     
                     
                         bKeepCurrent = false;
                         
                         // Next statement returns non-zero nStat if cell is bogus
                         // Note that i is the lattice character
                         nStat = nReducedPToCent(nLatt, a6fGprime);
                         
                         if (0 != nStat)
                            continue;

                         nStat = nCalcGetCellFromG6(a6fGprime, &a6fCell[0]);
                         oCrystal.vSetCell(a6fCell);

                         if( 0 != nStat )
                            continue;

                         // We are saving the best lattice character, regardless of other information.
                         if( fError < m_a44fLattResid[nLatt] && nLoop == 0 && nSlot==0 )
                         {
                             m_a44fLattResid[nLatt] = fError;
                         
                             oCrystal.vGetCell(&m_a44x6fLattCells[nLatt][0]);
                         }
                         //////////////////////////////////////////////////////////////////////////////

                         fRot2 = (double)nSol;
                     
                         fRot3 = pfCellIntensity[nSol];
                                              
                         // Similar to fErrorLimit, except that it applies 
                         // to the sum of length param.  (currently stored in fRot1)
                     
                         // Determine whether this least squares fit 
                         // is better than a previous one
                         // Place the better fit into oRefln.
                         
                         // RB 12/20/06 Trying to rationalize what Thad is doing... Some monoclinic and triclinic solutions do not have to pass
                         // the "lowest error in group" test. That's because without refinement it is hard to say which solution
                         // is the right one based on the error only. So if a solution is close enough
                         // in error to the lowest error, but has other advantages such as a lower beta angle, we
                         // will prefer that solution.
                         bSpecLatType = "mC" == m_a44tLattChar[nLatt].sLatticeType ||                   
                                        "mP" == m_a44tLattChar[nLatt].sLatticeType || 
                                        "aP" == m_a44tLattChar[nLatt].sLatticeType;

                         if ((nSlot==0) && (nLoop == 1) && (fError > fErrorLimitTol +  fLowestInLattice))
                            bKeepCurrent = FALSE;
                         else if ((nSlot==0) && (nLoop == 0) && (fError < oRefln.fGetIntensity()))
                            bKeepCurrent = TRUE;
                         else  if ((nLoop == 1) && ((oCrystal.fCalcVolume()> a2fVolumeFractionLimits[1]*fPrevVolume) || (oCrystal.fCalcVolume() <a2fVolumeFractionLimits[0]*fPrevVolume)))
                            bKeepCurrent = FALSE;                             
                        else if( (fErrorLimit * min(fError, fLowest) >= fabs(fLowest - fError) ||   
                                   (fError < m_fResidualMax && bSpecLatType))
                                   && fabs(fLowest - fError) < 10.0 )
                         {
                             // Difference between previous & current is <2% of min(errCurr,errPrev)
                             // or less than 2%.  This means they are close.
                         
                             // Next compare CurrSumCell with PrevSumCell, pick smaller of the two.
                             // unless cell is monoclinic and beta is out of the allowed region

                             if( ("mP" == m_a44tLattChar[nLatt].sLatticeType || "aP" == m_a44tLattChar[nLatt].sLatticeType) &&
                                 (90.0f > oCrystal.fGetCell(4) ||  120.5f < oCrystal.fGetCell(4)) )
                             {
                                 bKeepCurrent = FALSE; 
                             }
                             else if( "mC" == m_a44tLattChar[nLatt].sLatticeType && 90.0f > oCrystal.fGetCell(4) )
                             {
                                 bKeepCurrent = FALSE; 
                             }
                             else if( "mC" == m_a44tLattChar[nLatt].sLatticeType &&
                                      oRefln.fGetField(m_nFI_fBeta) - fMonoclinicDegreeLimit > oCrystal.fGetCell(4) )
                             {
                                 // Keep the one with the lowest beta angle!
                                 bKeepCurrent=TRUE;
                             }
                             else if( "mC" == m_a44tLattChar[nLatt].sLatticeType )
                             {
                                 // Need this one too so does not use the ELSE part of this big IF-ELSE...
                             
                                 // Keep the one with the lowest beta angle!
                             
                                 bKeepCurrent = FALSE;

                             }
                             else if(  fErrorLimitRot * min(fRot1, oRefln.fGetField(m_nFI_fRot1)) >=
                                       fabs(fRot1 - oRefln.fGetField(m_nFI_fRot1)) )
                             {
                                     // Well, sum cell lengths are close, try to pick based on angles
                                     double     fSumSin0 = sin(oCrystal.fGetCell(3) * Gs_dRADIANS_PER_DEGREE)
                                                           + sin(oCrystal.fGetCell(4) * Gs_dRADIANS_PER_DEGREE)
                                                           + sin(oCrystal.fGetCell(5) * Gs_dRADIANS_PER_DEGREE);
                                 
                                     double     fSumSin1 = sin(oRefln.fGetField(m_nFI_fAlpha)* Gs_dRADIANS_PER_DEGREE)
                                                           + sin(oRefln.fGetField(m_nFI_fBeta)* Gs_dRADIANS_PER_DEGREE)
                                                           + sin(oRefln.fGetField(m_nFI_fGamma)* Gs_dRADIANS_PER_DEGREE);
                                 
                                     // If the errors for the sum of sin's are within tolerance, then we will resort
                                     // back to the best cell error.
                                 
                                     if (fabs((fSumSin0 - fSumSin1) / (fSumSin0 + fSumSin1))
                                         < fErrorLimitSumSin)
                                     
                                         bKeepCurrent = (fError < oRefln.fGetIntensity());
                                     else
                                         bKeepCurrent = ( fSumSin0 > fSumSin1 );
                                 
                              }
                              else if (fRot1 < oRefln.fGetField(m_nFI_fRot1))
                              {
                                // Current has lower sum(cell lengths), so is better
                                 
                                bKeepCurrent = TRUE;
                              }
                                 // TODO:  Compare on fIntensity and pick better one!
                         }
                         else if (fError < oRefln.fGetIntensity())
                         {
                             // Residuals not close, CURRENT is better than PREVIOUS
                         
                             bKeepCurrent = TRUE;
                         }

                         if (43 == nLatt) 
                            bKeepCurrent = FALSE;
                         //////////////////////////////////////////////////////////////////////////////////////////////////////

                        if( bKeepCurrent )
                        {
                            // Only keep monoclinic cells if beta > 89 deg.
                    
                            if ( ("mP" == m_a44tLattChar[nLatt].sLatticeType ||  "aP" == m_a44tLattChar[nLatt].sLatticeType) &&
                                 (90.0f > oCrystal.fGetCell(4) || 120.0f < oCrystal.fGetCell(4)) )
                            {
                                bKeepCurrent = FALSE; 
                            }
                            else if( "mC" == m_a44tLattChar[nLatt].sLatticeType &&
                                     90.0f > oCrystal.fGetCell(4) )
                            {
                                bKeepCurrent = FALSE; 
                            }
                        }
                
                
                        if( bKeepCurrent &&  oRefln.fGetIntensity() <= m_fResidualMax &&  fError > m_fResidualMax )
                        {
                            bKeepCurrent = FALSE;
                        }

                        if( !bKeepCurrent )
                            continue;

                        // Now, convert to standard pres on the crystal cell.  Since we have simple lattice types, this amounts to
                        // checking the order of edges for orthorhombic, (increasing in size),
                        oCrystal.vGetCell(a6fCell);
                
                        // cP  cI cF hR hP tI tP oI oC oF oP mC mP aP //
                        // 0   1  2  3  4  5  6  7  8  9  10 11 12 13
                
                        //                  if (nLatt == 43) 
                        //                    continue;                  
                        switch (pcLattString[0]) 
                        {
                        case 'o':
                            if( pcLattString[1] == 'C' )
                            {
                                if( a6fCell[0] > a6fCell[1] ) 
                                    std::swap(a6fCell[0], a6fCell[1]); 
                            } 
                            else
                            {
                                if (a6fCell[0]>a6fCell[2]) 
                                    std::swap(a6fCell[0],a6fCell[2]);
                                if (a6fCell[1]>a6fCell[2]) 
                                    std::swap(a6fCell[1],a6fCell[2]);
                                if (a6fCell[0]>a6fCell[1]) 
                                    std::swap(a6fCell[0],a6fCell[1]);
                            }
                            break;
                        case 'm':
                            if (pcLattString[1] == 'P')
                            {
                                if (a6fCell[0]>a6fCell[2]) 
                                    std::swap(a6fCell[0],a6fCell[2]);
                            } 
                            break;
                        case 'a':
                            nStandardPres(a6fCell);
                        } // end switch
                
                        oCrystal.vSetCell(a6fCell);
                
                        // Keep the curent solution which is better than the previous one
                        
                        fLowest = min(fLowest, fError);

                        if (nLoop == 0) 
                        {                        
                            oRefln.vSetH(m_a15x3nBIndices[nSol][0]);
                            oRefln.vSetK(m_a15x3nBIndices[nSol][1]);
                            oRefln.vSetL(m_a15x3nBIndices[nSol][2]);                        
                        }
                        
                        oRefln.vSetIntensity((float)fError);

                        oRefln.vSetSigmaI((float)m_a44tLattChar[nLatt].nSpacegroupMin);
                        oRefln.vSetField(m_nFI_fA,     (float)oCrystal.fGetCell(0));
                        oRefln.vSetField(m_nFI_fB,     (float)oCrystal.fGetCell(1));
                        oRefln.vSetField(m_nFI_fC,     (float)oCrystal.fGetCell(2));
                        oRefln.vSetField(m_nFI_fAlpha, (float)oCrystal.fGetCell(3));
                        oRefln.vSetField(m_nFI_fBeta,  (float)oCrystal.fGetCell(4));
                        oRefln.vSetField(m_nFI_fGamma, (float)oCrystal.fGetCell(5));
                        oRefln.vSetField(m_nFI_fVolume, (float)oCrystal.fCalcVolume());
                        oRefln.vSetField(m_nFI_fRot1,  (float)fRot1);
                        oRefln.vSetField(m_nFI_fRot2,  (float)fRot2);
                        oRefln.vSetField(m_nFI_fRot3,  (float)fRot3);
                        oRefln.vSetField(m_nFI_sLattice, 
                            m_a44tLattChar[nLatt].sLatticeType);
                        oRefln.vSetField(m_nFI_sLattSymm, 
                            m_a44tLattChar[nLatt].sLattSymm);
                        oRefln.vSetField(m_nFI_nSortOrder, 
                            m_a44tLattChar[nLatt].nSortOrder);
                        oRefln.vSetField(m_nFI_nLattNum, nLatt);
                
                        // Resave the lattice again.
                        if (nLoop == 0)
                        {
                            m_a44fLattResid[nLatt]=fError;
                            oCrystal.vGetCell(&m_a44x6fLattCells[nLatt][0]);                    
                        }
                    } // end permutations loop
                } // end j loop           
            } // end nSol loop
        
            // Insert the ONE best solution for each lattice type
            // There should be only 14 insertions
        
            if (   (ms_anOrder[ii+1]==-1)
                || (   m_a44tLattChar[nLatt].sLatticeType
                != m_a44tLattChar[ms_anOrder[ii+1]].sLatticeType) ||
                            (ms_anOrderSlot[ii+1]!=ms_anOrderSlot[ii]))
            {
                if (nLoop == 0)
                {
                    f0 = oRefln.fGetIntensity();
                    m_poSolutions->nInsert(&oRefln);
                } 
                else 
                {
                    f0 = (*m_poSolutions)[nNumStart + nBravisLattice].fGetIntensity();
                    // TJN: 7/16/03
                    // Discovered that sometimes no crystal shows up on the second pass?!?
                    // I do not know why, but this extra 'if' should gaurd against that contingency.
                    if (oRefln.fGetIntensity() < 1000.0)
                    {
                        (*m_poSolutions)[nNumStart + nBravisLattice] = oRefln;
                        
                        // RB 12/20/06 Making sure that a larger error from the second loop does not get wiped out
                        // by a smaller error from the first loop.
                        if( f0 > oRefln.fGetIntensity() )
                            (*m_poSolutions)[nNumStart + nBravisLattice].vSetIntensity((float) f0);
                    }
                }
            
                oRefln.vSetIntensity(1000.0f);  
            
                // Make residual high so bKeepCurrent will be TRUE
            }
        } // end of i (latt char) loop

        printf("\n");
        
        // RB Debug
        //{
        //    Cstring     sLoop(nLoop);
        //
        //    Cstring     sFileName("solutions");
        //    sFileName += '_';
        //    sFileName += sLoop;
        //    sFileName += ".ref";
        //    
        //    nEnd   = m_poSolutions->nGetNumReflns();
        //    m_poSolutions->nWrite(sFileName, NULL, NULL, NULL, nNumStart, nEnd-1);
        //}
    } // end of nLoop trial loop (goes twice)


    cout << "done.\n" << flush;
    nEnd   = m_poSolutions->nGetNumReflns();

    
    if (4 <= m_nVerbose)
    {
        m_poSolutions->nWrite("soln44a.ref", NULL, NULL, NULL,nNumStart, nEnd-1);
    }
    
    // Delete unwanted reflections.
    int*        pnDelFlag = new int[nEnd];
    

    // Delete solutions not inserted at this time
    for (ii = 0; ii < nEnd; ii++)
    {
        pnDelFlag[ii] = 0;
        if (ii < nNumStart) 
            pnDelFlag[ii] = -1;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Remove redundant solutions.  This is important, because we have two
    // solutions for mP and mC.  We might only need one solution for each.
    
    // mC
    int             nSolution_A = nNumStart + 11;
    int             nSolution_B = nNumStart + 12;
    
    if( bEquivalentCellSolutions((*m_poSolutions)[nSolution_A], (*m_poSolutions)[nSolution_B]) )
    {
        if( (*m_poSolutions)[nSolution_A].fGetIntensity() < (*m_poSolutions)[nSolution_B].fGetIntensity() )
            pnDelFlag[nSolution_B] = -1;
        else
            pnDelFlag[nSolution_A] = -1;
    }
    
    // mP
    nSolution_A = nNumStart + 13;
    nSolution_B = nNumStart + 14;
    
    if( bEquivalentCellSolutions((*m_poSolutions)[nSolution_A], (*m_poSolutions)[nSolution_B]) )
    {
        if( (*m_poSolutions)[nSolution_A].fGetIntensity() < (*m_poSolutions)[nSolution_B].fGetIntensity() )
            pnDelFlag[nSolution_B] = -1;
        else
            pnDelFlag[nSolution_A] = -1;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    m_poSolutions->nDelete(-1, pnDelFlag);
    delete [] pnDelFlag;
    
    
    if (4 <= m_nVerbose)
    {
        m_poSolutions->nWrite("soln44d.ref");
    }
    
    delete [] pfReducedCells;
    delete [] pfCellIntensity;
    
    return 0;
}

int CCellReduce::nCentToPrim(char cCentering, double *pfIO,bool bIsG6)
{
    // Convert a G6 (or cell) representation from a centered one to
    // a primitive one.  The primitive one is not reduced!
    // See Andrews & Bernstein (1988) Table 1.
    
    double fMat6[6][6];
    double fGprime[6];
    double fMatFactor;
    double pfG6[6];
    
    int a6x2nCellAxis[6][2] = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};
    double a3x3fRealSpacePToC[3][3];
    double a3x3fRealSpaceCToP[3][3];

    double fDet;
    int nx,ny;
    int nRow,nCol;


    {
        
        vZeroMat(6, 6, &fMat6[0][0]);
        
        if ('I' == cCentering)
        {
            double a3x3fRealSpacePToC_[3][3] = {{1,0,-1},{0,1,-1},{0,0,2}};
            vCopyMat3D(&a3x3fRealSpacePToC_[0][0],&a3x3fRealSpacePToC[0][0]);
            
            fMatFactor  = 1.0 / 4.0;
            fMat6[0][0] = 4.0;
            fMat6[0][2] = 1.0;
            fMat6[0][4] = 4.0;
            fMat6[1][1] = 4.0;
            fMat6[1][2] = 1.0;
            fMat6[1][3] = 4.0;
            fMat6[2][2] = 1.0;
            fMat6[3][2] = 1.0;
            fMat6[3][3] = 2.0;
            fMat6[4][2] = 1.0;
            fMat6[4][4] = 2.0;
            fMat6[5][2] = 1.0;
            fMat6[5][3] = 2.0;
            fMat6[5][4] = 2.0;
            fMat6[5][5] = 4.0;
        }
        else if ('A' == cCentering)
        {
            double a3x3fRealSpacePToC_[3][3] = {{1,0,0},{0,1,-1},{0,0,2}};
            vCopyMat3D(&a3x3fRealSpacePToC_[0][0],&a3x3fRealSpacePToC[0][0]);

            fMatFactor  = 1.0 / 4.0;
            fMat6[0][0] = 4.0;
            fMat6[1][1] = 4.0;
            fMat6[1][2] = 1.0;
            fMat6[1][3] = 4.0;
            fMat6[2][2] = 1.0;
            fMat6[3][2] = 1.0;
            fMat6[3][3] = 2.0;
            fMat6[4][4] = 2.0;
            fMat6[5][4] = 2.0;
            fMat6[5][5] = 4.0;
        }
        else if ('C' == cCentering)
        {
            double a3x3fRealSpacePToC_[3][3] = {{1,-1,0},{0,2,0},{0,0,1}};
            vCopyMat3D(&a3x3fRealSpacePToC_[0][0],&a3x3fRealSpacePToC[0][0]);
            
            fMatFactor  = 1.0 / 4.0;
            fMat6[0][0] = 4.0;
            fMat6[0][1] = 1.0;
            fMat6[0][5] = 4.0;
            fMat6[1][1] = 1.0;
            fMat6[2][2] = 4.0;
            fMat6[3][3] = 2.0;
            fMat6[4][3] = 2.0;
            fMat6[4][4] = 4.0;
            fMat6[5][1] = 1.0;
            fMat6[5][5] = 2.0;
        }
        else if ('B' == cCentering)
        {
            double a3x3fRealSpacePToC_[3][3] = {{2,0,0},{0,1,0},{-1,0,1}};
            vCopyMat3D(&a3x3fRealSpacePToC_[0][0],&a3x3fRealSpacePToC[0][0]);

            fMatFactor  = 1.0 / 4.0;
            fMat6[0][0] = 4.0;
            fMat6[0][2] = 1.0;
            fMat6[0][4] = 4.0;
            fMat6[1][1] = 4.0;
            fMat6[2][2] = 1.0;
            fMat6[3][3] = 2.0;
            fMat6[4][2] = 1.0;
            fMat6[4][4] = 2.0;
            fMat6[5][3] = 2.0;
            fMat6[5][5] = 4.0;
        }
        else if ('F' == cCentering)
        {
            double a3x3fRealSpacePToC_[3][3] = {{1,1,-1},{-1,1,1},{1,-1,1}};
            // The matrix below will yield an equivalent 6x6 transformation as that specified.  But it's got a negative determinant!
            // double a3x3fRealSpacePToC_[3][3] = {{1,1,-1},{1,-1,1},{-1,1,1}};
            vCopyMat3D(&a3x3fRealSpacePToC_[0][0],&a3x3fRealSpacePToC[0][0]);

            fMatFactor  = 1.0 / 4.0;
            fMat6[0][0] = 1.0;
            fMat6[0][1] = 1.0;
            fMat6[0][5] = 2.0;
            fMat6[1][0] = 1.0;
            fMat6[1][2] = 1.0;
            fMat6[1][4] = 2.0;
            fMat6[2][1] = 1.0;
            fMat6[2][2] = 1.0;
            fMat6[2][3] = 2.0;
            fMat6[3][2] = 1.0;
            fMat6[3][3] = 1.0;
            fMat6[3][4] = 1.0;
            fMat6[3][5] = 1.0;
            fMat6[4][1] = 1.0;
            fMat6[4][3] = 1.0;
            fMat6[4][4] = 1.0;
            fMat6[4][5] = 1.0;
            fMat6[5][0] = 1.0;
            fMat6[5][3] = 1.0;
            fMat6[5][4] = 1.0;
            fMat6[5][5] = 1.0;
        }
        else if ('R' == cCentering)
        {
            
            double a3x3fRealSpacePToC_[3][3] = {{1,0,1},{-1,1,1},{0,-1,1}};
            vCopyMat3D(&a3x3fRealSpacePToC_[0][0],&a3x3fRealSpacePToC[0][0]);
            
            fMatFactor  = 1.0f / 9.0f;
            fMat6[0][0] = 4.0;
            fMat6[0][1] = 1.0;
            fMat6[0][2] = 1.0;
            fMat6[0][3] = -2.0;
            fMat6[0][4] = -4.0;
            fMat6[0][5] = -4.0;
            
            fMat6[1][0] = 4.0;
            fMat6[1][1] = 1.0;
            fMat6[1][2] = 4.0;
            fMat6[1][3] = -4.0;
            fMat6[1][4] = -8.0;
            fMat6[1][5] = 4.0;
            
            fMat6[2][0] = 4.0;
            fMat6[2][1] = 1.0;
            fMat6[2][2] = 1.0;
            fMat6[2][3] = 2.0;
            fMat6[2][4] = 2.0;
            fMat6[2][5] = 2.0;
            
            fMat6[3][0] = 2.0;
            fMat6[3][1] = 1.0;
            fMat6[3][2] = -2.0;
            fMat6[3][3] = -1.0;
            fMat6[3][4] = 0.0;
            fMat6[3][5] = 3.0;
            
            fMat6[4][0] = 2.0;
            fMat6[4][1] = -1.0;
            fMat6[4][2] = -1.0;
            fMat6[4][3] = 0.0;
            fMat6[4][4] = 1.0;
            fMat6[4][5] = 1.0;
            
            fMat6[5][0] = 4.0;
            fMat6[5][1] = -1.0;
            fMat6[5][2] = -2.0;
            fMat6[5][3] = 1.0;
            fMat6[5][4] = -6.0;
            fMat6[5][5] = 0.0;
            
        }
        else
        {
            double a3x3fRealSpacePToC_[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
            vCopyMat3D(&a3x3fRealSpacePToC_[0][0],&a3x3fRealSpacePToC[0][0]);

            fMatFactor  = 1.0;
            fMat6[0][0] = 1.0;
            fMat6[1][1] = 1.0;
            fMat6[2][2] = 1.0;
            fMat6[3][3] = 1.0;
            fMat6[4][4] = 1.0;
            fMat6[5][5] = 1.0;
        }
        
        
        // Routine to generate matrix from original reindex matrix
        {
            
            fDet = fInvMat3D(&a3x3fRealSpacePToC[0][0],&a3x3fRealSpaceCToP[0][0]);
            
            for (nx=0;nx<3;nx++) {
                for (ny=0;ny<3;ny++) {
                    a3x3fRealSpaceCToP[nx][ny]*=fDet;
                };
            };
            for (nRow = 0; nRow < 6; nRow++) {
                for (nCol = 0; nCol < 6; nCol++) {
                    if (nCol >=3) {
                        fMat6[nCol][nRow] =                  
                            ((nRow>=3)?(2):(1))*((nCol>=3)?(0.5):(1))*(
                            a3x3fRealSpaceCToP[a6x2nCellAxis[nCol][0]][a6x2nCellAxis[nRow][0]]*
                            a3x3fRealSpaceCToP[a6x2nCellAxis[nCol][1]][a6x2nCellAxis[nRow][1]] +
                            a3x3fRealSpaceCToP[a6x2nCellAxis[nCol][1]][a6x2nCellAxis[nRow][0]]*
                            a3x3fRealSpaceCToP[a6x2nCellAxis[nCol][0]][a6x2nCellAxis[nRow][1]]);
                    } else {
                        fMat6[nCol][nRow] =                  
                            ((nRow>=3)?(2):(1))*((nCol>=3)?(0.5):(1))*(
                            a3x3fRealSpaceCToP[a6x2nCellAxis[nCol][0]][a6x2nCellAxis[nRow][0]]*
                            a3x3fRealSpaceCToP[a6x2nCellAxis[nCol][1]][a6x2nCellAxis[nRow][1]]);
                    };
                };
            };
            fMatFactor = 1.0/(fDet*fDet);
        };
        
        
        if (bIsG6)
            vCopyVecND(6,pfIO,&pfG6[0]);
        else
            vCalcGetG6(pfIO,&pfG6[0]);
        
        
        vMulVecNDScalar(6*6, &fMat6[0][0], fMatFactor, &fMat6[0][0]);
        vMulMatNDVecND(6,  &fMat6[0][0], &pfG6[0], &fGprime[0]);
        vCopyVecND(6, fGprime, pfG6);
        
        if (bIsG6)
            vCopyVecND(6,&pfG6[0],pfIO);
        else
            nCalcGetCellFromG6(&pfG6[0],pfIO);

        
    }
    
    
    return (0);
}

int 
CCellReduce::nReducedPToCent(const int nLattNumIn, double *pfG6)
{
  // Convert a cell represented as a G6 vector from the reduced primitive one
  // back to the centered one.

  double   a3x3fTMat[3][3];
  double   a3x3fTempMat[3][3];
  int     nLattNum;

  nLattNum = nLattNumIn;
  if (1 > nLattNum)
    nLattNum = m_nLattNum;  // Use the member variable lattice number

  // Get the transformation matrix needed

  vCopyMat3D(&m_a44tLattChar[nLattNum].a3x3fTMat[0][0], &a3x3fTMat[0][0]);

  int   j;
  double a6fCell[6];
  double a3fVec1[3], a3fVec2[3], a3fVec3[3];
  double fCos;
  double a3x3fRealMat[3][3];
  Ccrystal oCrystal;

  (void) nCalcGetCellFromG6(pfG6, a6fCell);
  oCrystal.vSetCell(a6fCell); 

  if (0 > oCrystal.nCalcGetRealMatrix(&a3x3fRealMat[0][0]))
    {
      // It is a bogus cell
//      cout << "WARNING bogus cell in ::nReducedPToCent!\n" << flush;
      pfG6[0] = 1.0f;
      pfG6[1] = 1.0f;
      pfG6[2] = 1.0f;
      pfG6[3] = 0.0f;
      pfG6[4] = 0.0f;
      pfG6[5] = 0.0f;
      return (-1);
    }

  // Multiply T by real matrix to get transformed direct cell matrix
  // for cell determination

  vMulMat3DMat3D(a3x3fTMat, a3x3fRealMat, a3x3fTempMat);

  // Get direct space cell constants (what about wavelength?)

  for (j = 0; j < 3; j++)
    {
      a3fVec1[j] = a3x3fTempMat[j][0];
      a3fVec2[j] = a3x3fTempMat[j][1];
      a3fVec3[j] = a3x3fTempMat[j][2];
    }

  a6fCell[0] = fLenVec3D(a3fVec1);
  a6fCell[1] = fLenVec3D(a3fVec2);
  a6fCell[2] = fLenVec3D(a3fVec3);

  fCos = fDot3D(a3fVec2, a3fVec3)
          / (a6fCell[1] * a6fCell[2]);
  fCos = min(fCos, 1.0f);
  fCos = max(fCos, -1.0f);
  a6fCell[3] = acos(fCos) / Gs_dRADIANS_PER_DEGREE;

  fCos = fDot3D(a3fVec1, a3fVec3)
                 / (a6fCell[0] * a6fCell[2]);
  fCos = min(fCos, 1.0f);
  fCos = max(fCos, -1.0f);
  a6fCell[4] = acos(fCos) / Gs_dRADIANS_PER_DEGREE;

  fCos = fDot3D(a3fVec1, a3fVec2)
                 / (a6fCell[0] * a6fCell[1]);
  fCos = min(fCos, 1.0f);
  fCos = max(fCos, -1.0f);
  a6fCell[5] = acos(fCos) / Gs_dRADIANS_PER_DEGREE;

  vCalcGetG6(a6fCell, pfG6);  

  return (0);
}



int
CCellReduce::nReducedPToCent(const Cstring& sFlag, double *pfG6, int *pnLattNum)
{
  // Convert a reduced primitive G6 vector back to the centered
  // G6 vector according to the centering in sFLag
  // These transformations have been compiled from
  // Roof, R.B. (1967) Los Alamos Report 4038, UC-4, Chemistry, TID-4500
  // Lawton, S.L. & Jacobson, R.A (1965) IS-1141, UC-4, Chemistry, TID-4500
  // and probably come originally from some Niggli publication.
  // If pnLattNum is not NULL, then also return the lattice number according
  // to Table 9.3.1 in IntTable Vol A, p. 746.

  // Probably need to determine whether cell is TYPE I or TYPE II and use this
  // below.

  Cstring sTemp;
  sTemp = sFlag;
  double fTemp;
  int   nType;
  int   nLattNum = 0;
  // Be sure to add divide by 0.0 check below
  if (pfG6[0] <= 0.0) return (1);

  // Compute reduced cell Type.  See pp. 741-743 Int Tables.
  nType = 1;
  if (0.0 >= (pfG6[3] * pfG6[4] * pfG6[5]) )
    nType = 2;

  if (sTemp == "cF")
    {
      nLattNum = 1;
      pfG6[0] = 2.0f * (pfG6[0] + pfG6[1] + pfG6[2]) / 3.0f;  // Average lengths!
      pfG6[1] = pfG6[0];
      pfG6[2] = pfG6[0];
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "cI")
    {
      nLattNum = 5;
      pfG6[0] = (4.0f / 3.0f) * (pfG6[0] + pfG6[1] + pfG6[2]) / 3.0f;
      pfG6[1] = pfG6[0];
      pfG6[2] = pfG6[0];
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "cP")
    {
      nLattNum = 3;
      pfG6[0] = (pfG6[0] + pfG6[1] + pfG6[2]) / 3.0f; // Average them
      pfG6[1] = pfG6[0];
      pfG6[2] = pfG6[0];
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "tI")
    {
      fTemp = fabs((pfG6[0] - pfG6[1]) / pfG6[0]);  // if g0=g1
      if (fTemp < 0.002)
        {
          fTemp = fabs((pfG6[1] - pfG6[2]) / pfG6[1]);  // if g1=g2
          if (fTemp < 0.002)
            {
              // Here g0 = g1 = g2
              int nA, nB, nC;
              fTemp = fabs((pfG6[3] - pfG6[4]) / pfG6[3]);  // if g3=g4
              if (fTemp < 0.1)
                {
                  nLattNum = 6;
                  nA = 3;
                  nB = 4;
                  nC = 5;
                }
              else if (fabs((pfG6[4] - pfG6[5]) / pfG6[4]) < 0.002)
                {  // if g4=g5
                  nLattNum = 7;
                  nA = 5;
                  nB = 4;
                  nC = 3;
                }
              else
                {   // if g3 = g5
                  nLattNum = -1;  // Impossible??
                  nA = 3;
                  nB = 5;
                  nC = 4;
                }
              pfG6[2] = 2.0f * pfG6[2] + pfG6[nC];
              pfG6[0] = (pfG6[0] + pfG6[1]) + 0.5f * (pfG6[nA] + pfG6[nB]);
              pfG6[1] = pfG6[0];
            }
          else
            {
              // Here g0 = g1 != g2
              nLattNum = 15;
              pfG6[0] = (pfG6[0] + pfG6[1]) / 2.0f;
              pfG6[1] = pfG6[0];
              pfG6[2] = 4.0f * pfG6[2] - 2.0f * pfG6[0];
            }
        }
      else
        {
          // Assume g0 != g1 = g2 here
          nLattNum = 18;
          fTemp   = pfG6[0];
          pfG6[0] = pfG6[1] + pfG6[2] - 0.5f * pfG6[0];
          pfG6[1] = pfG6[0];
          pfG6[2] = fTemp;
        }
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "oF")
    {
      nLattNum = 26;
      fTemp = fabs((pfG6[0] - pfG6[1]) / pfG6[0]);  // Could also test if g0=g1
      if (fTemp < 0.002)
        {
          nLattNum = 16;
          pfG6[0] = -2.0f * pfG6[3];
        }
      pfG6[1] = 4.0f * pfG6[1] - pfG6[0];
      pfG6[2] = 4.0f * pfG6[2] - pfG6[0];
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "oI")
    {
      //
      fTemp = fabs((pfG6[0] - pfG6[1]) / pfG6[0]);  // if g0=g1
      if (fTemp < 0.002)
        {
          // A = B
          fTemp = fabs((pfG6[1] - pfG6[2]) / pfG6[1]);  // if g1=g2
          if (fTemp < 0.002)
            {
              // Here g0 = g1 = g2 (A = B = C)
              nLattNum = 8;
              double fMax = max(pfG6[3], max(pfG6[4], pfG6[5]));
              double fMin = min(pfG6[3], min(pfG6[4], pfG6[5]));
              double fMid = pfG6[3] + pfG6[4] + pfG6[5] - fMin - fMax;
              pfG6[0] = 2.0f * pfG6[0] + fMin;
              pfG6[1] = 2.0f * pfG6[1] + fMid;
              pfG6[2] = 2.0f * pfG6[2] + fMax;
            }
          else if (fabs(pfG6[5]/pfG6[1]) < 0.002)
            {
              nLattNum = 42;
              pfG6[2] = 4.0f * pfG6[2] - pfG6[0] - pfG6[1];
            }
          else
            {
              // Here g0 = g1 != g2; make sure g5 != 0  (probably impossible?!)

              nLattNum = 19;
              pfG6[1] = 2.0f * pfG6[2] - pfG6[5];
              pfG6[2] = 2.0f * pfG6[2] + pfG6[5] - pfG6[0];
            }
        }
      else if ( fabs((pfG6[1] - pfG6[2]) / pfG6[1]) < 0.002)
        {
          // g1 = g2 (A != (B = C))
          nLattNum = 19;
          pfG6[1] = 2.0f * pfG6[2] - pfG6[3];
          pfG6[2] = 2.0f * pfG6[2] + pfG6[3] - pfG6[0];
        }
      else
        {
          // Here g0 != g1 != g2
          nLattNum = 42;
          pfG6[2] = 4.0f * pfG6[2] - pfG6[0] - pfG6[1];
        }
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "oC")
    {
      // Find which of g3, g4, g5 has abs() farthest from 0.
      int   nFar = 3;
      double fFar = fabs(pfG6[3]);
      if (fFar < fabs(pfG6[4]))
        {
          nFar = 4;
          fFar = fabs(pfG6[4]);
        }
      if (fFar < fabs(pfG6[5]))
        {
          nFar = 5;
          fFar = fabs(pfG6[5]);
        }

      fTemp = fabs((pfG6[0] - pfG6[1]) / pfG6[0]);
      if (fTemp < 0.002)
        {
          // if g0=g1
          nLattNum = 13;
          fTemp   = pfG6[0];
          pfG6[0] = pfG6[0] + pfG6[1] - fFar;  // Average g0 and g1
          pfG6[1] =   fTemp + pfG6[1] + fFar;
        }
      else if (fabs((pfG6[1] - pfG6[2]) / pfG6[1]) < 0.002)
        {
          // g0 != g1 = g2
          nLattNum = 23;
          fTemp   = pfG6[0];
          pfG6[0] = pfG6[1] + pfG6[2] - fFar;  // Average g1 and g2
          pfG6[1] = pfG6[1] + pfG6[2] + fFar;
          pfG6[2] = fTemp;
        }
      else
        {
          // Here g0 != g1 != g2
          if (5 == nFar)
            {
              nLattNum = 38;
              pfG6[1] = 4.0f * pfG6[1] - pfG6[0];
            }
          else if (4 == nFar)
            {
              nLattNum = 36;
              fTemp   = pfG6[1];
              pfG6[1] = 4.0f * pfG6[2] - pfG6[0];
              pfG6[2] = fTemp;
            }
          else
            {  // nFar == 3
              nLattNum = 40;
              fTemp   = pfG6[1];
              pfG6[1] = 4.0f * pfG6[2] - pfG6[1];
              pfG6[2] = pfG6[0];
              pfG6[0] = fTemp;
            }
        }
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "tP")
    {
      fTemp = fabs((pfG6[0] - pfG6[1]) / pfG6[0]);
      if (fTemp < 0.002)
        {
          // A = B
          nLattNum = 11;
          // if g0 = g1
          pfG6[1] = (pfG6[0] + pfG6[1]) * 0.5f;  // Average them
        }
      else if (fabs((pfG6[1] - pfG6[2]) / pfG6[1]) < 0.002)
        {
          // B = C
          nLattNum = 21;
          pfG6[1] = (pfG6[2] + pfG6[1]) * 0.5f;  // Average them
          pfG6[2] = pfG6[0];
        }
      pfG6[0] = pfG6[1];
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = 0.0;
    }
  else if (sTemp == "hP")
    {
      if (fabs(pfG6[5]) > fabs(pfG6[3]))
        {
          // if g0 = g1
          nLattNum = 12;
          pfG6[1] = (pfG6[0] + pfG6[1]) * 0.5f;  // Average them
        }
      else
        {
          nLattNum = 22;
          pfG6[1] = (pfG6[2] + pfG6[1]) * 0.5f;  // Average them
          pfG6[2] = pfG6[0];
        }
      pfG6[0] = pfG6[1];
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = -1.0f * pfG6[1];
    }
  else if (sTemp == "hR")
    {   // Convert primitive rhombohedral to hexagonal
      fTemp = fabs((pfG6[0] - pfG6[1]) / pfG6[0]);
      if (fTemp < 0.002)
        {
          // if g0 = g1
          if (fabs((pfG6[1] - pfG6[2]) / pfG6[1]) < 0.002)
            {
              // g0 = g1 = g2
              nLattNum = 2;  // Could be 4 too, see TYPE column in I.T.
              if (2 == nType)
                nLattNum = 4;
              pfG6[0] = (pfG6[0] + pfG6[1] + pfG6[2]) / 3.0f; // Average them
              pfG6[3] = (pfG6[3] + pfG6[4] + pfG6[5]) / 3.0f; // Average them
              pfG6[2] = 3.0f * (pfG6[0] + pfG6[3]);
              pfG6[0] = 2.0f * pfG6[0] - pfG6[3];
            }
          else
            {
              // g0 = g1 != g2
              nLattNum = 9;
              pfG6[0] = (pfG6[0] + pfG6[1]) * 0.5f;
              pfG6[2] = 3.0f * (3.0f * pfG6[2] - pfG6[0]);
            }
        }
      else if (fabs((pfG6[1] - pfG6[2]) / pfG6[1]) < 0.002)
        {
          // g0 != g1 = g2
          nLattNum = 24;
          fTemp   = (pfG6[1] + pfG6[2]) * 0.5f; // Average them
          pfG6[2] = pfG6[0];
          pfG6[0] = 3.0f * fTemp - pfG6[0] / 3.0f;
        }
      pfG6[1] = pfG6[0];
      pfG6[3] = 0.0;
      pfG6[4] = 0.0;
      pfG6[5] = -1.0f * pfG6[1];
    }
  else if (sTemp == "mP")
    {
      // Find which of g3, g4, g5 has abs() farthest from 0.
      int   nFar = 3;
      double fFar = fabs(pfG6[3]);
      if (fFar < fabs(pfG6[4]))
        {
          nFar = 4;
          fFar = fabs(pfG6[4]);
        }
      if (fFar < fabs(pfG6[5]))
        {
          nFar = 5;
          fFar = fabs(pfG6[5]);
        }
      if (3 == nFar)
        {
          nLattNum = 35;
          fTemp   = pfG6[1];
          pfG6[1] = pfG6[0];
          pfG6[0] = fTemp;
          pfG6[4] = pfG6[3];
        }
      else if (5 == nFar)
        {
          nLattNum = 34;
          fTemp   = pfG6[1];
          pfG6[1] = pfG6[2];
          pfG6[2] = fTemp;
          pfG6[4] = pfG6[5];
        }
      else
        {
          nLattNum = 33;
        }
      pfG6[3] = 0.0;
      pfG6[5] = 0.0;
      if (pfG6[4] > 0.0) pfG6[4] = -pfG6[4];  // Make sure beta >90!
    }
  else if (sTemp == "mC")
    {
      // Find which of g3, g4, g5 has abs() farthest from 0.
      double fMax = max(fabs(pfG6[3]), max(fabs(pfG6[4]), fabs(pfG6[5])));
      double fMin = min(fabs(pfG6[3]), min(fabs(pfG6[4]), fabs(pfG6[5])));
      double fMid = fabs(pfG6[3]) + fabs(pfG6[4]) + fabs(pfG6[5]) - fMin - fMax;
      int   nFar = 3;
      double fFar = fabs(pfG6[3]);
      if (fFar < fabs(pfG6[4]))
        {
          nFar = 4;
          fFar = fabs(pfG6[4]);
        }
      if (fFar < fabs(pfG6[5]))
        {
          nFar = 5;
          fFar = fabs(pfG6[5]);
        }
      // See if g3, g4, or g5 = 0
      if (fabs(pfG6[3] / pfG6[0]) < 0.002)
        {
          nLattNum = -1;
          if (8 <= m_nVerbose)
            cout << "g3 is 0: impossible" << endl;
        }
      else if (fabs((pfG6[0] - pfG6[1]) / pfG6[0]) < 0.002)
        {
          //
          if (8 <= m_nVerbose)
            cout << "g0 = g1" << endl;
          if (fabs((pfG6[3] - pfG6[4]) / pfG6[3]) < 0.002)
            {
              // g3 = g4
              if (8 <= m_nVerbose)
                cout << "and g3 = g4; p. 36A" << endl;
              nLattNum = 14;  // COULD ALSO BE 10!
              if (1 == nType)
                nLattNum = 10;
              fTemp   = pfG6[0] + pfG6[1];  // Average them
              pfG6[0] = fTemp + pfG6[5];
              pfG6[1] = fTemp - pfG6[5];
              pfG6[4] = 2.0f * pfG6[4];
            }
          else
            {
              nLattNum = 17;   // TODO TO DO
              if (8 <= m_nVerbose)
                cout << "but g3 != g4; p. 38B" << endl;
            }
        }
      else if (fabs((pfG6[1] - pfG6[2]) / pfG6[1]) < 0.002)
        {
          //
          if (8 <= m_nVerbose)
            cout << "g1 = g2" << endl;
          if (fabs((pfG6[4] - pfG6[5]) / pfG6[5]) < 0.002)
            {
              // g4 = g5
              if (8 <= m_nVerbose)
                cout << "and g4 = g5; p. 36B" << endl;
              nLattNum = 20;
              fTemp   = pfG6[1] + pfG6[2];   // Average them
              pfG6[2] = pfG6[0];
              pfG6[0] = fTemp + pfG6[3];
              pfG6[1] = fTemp - pfG6[3];
              pfG6[4] = 2.0f * pfG6[4];
            }
          else
            {
              if (8 <= m_nVerbose)
                cout << "but g3 != g4; so impossible?" << endl;
              nLattNum = 25; // TODO TO DO
            }
        }
      else if (fabs(pfG6[4] / pfG6[0]) < 0.002)
        {
          if (8 <= m_nVerbose)
            cout << "g4 is 0; p.35A" << endl;
          nLattNum = 39;
          fTemp   = pfG6[1];
          pfG6[1] = pfG6[0];
          pfG6[0] = 4.0f * fTemp - pfG6[0];
          pfG6[4] = -2.0f * pfG6[3];  // Should this be fMid instead of g3?
        }
      else if (fabs(pfG6[5] / pfG6[0]) < 0.002)
        {
          if (fabs(pfG6[0] + pfG6[4]) < 0.002)
//          if (fabs(pfG6[0] - pfG6[4]) < 0.002)
            {
              if (8 <= m_nVerbose)
                cout << "g5 is 0; p. 35C" << endl;
              nLattNum = 37;
              fTemp   = pfG6[0];
              pfG6[0] = 4.0f * pfG6[2] - pfG6[0];
              pfG6[2] = pfG6[1];
              pfG6[1] = fTemp;
            }
          else
            {
              if (8 <= m_nVerbose)
                cout << "g5 is 0; p. 35B" << endl;
              nLattNum = 41;
              fTemp   = pfG6[0];
              pfG6[0] = 4.0f * pfG6[2] - pfG6[1];        
              pfG6[2] = fTemp;
            }
          pfG6[4] = -2.0f * fMid;
        }
      else
        {
          // g0 != g1 != g2
          if (8 <= m_nVerbose)
            cout << "g0 != g1 != g2" << endl;
          if (fabs((pfG6[0] -  pfG6[4]) / pfG6[0]) < 0.002)
            {
              if (8 <= m_nVerbose)
                cout << "g0 = g4; p. 37A" << endl;
              if (fabs((pfG6[4] -  pfG6[5]) / pfG6[0]) < 0.002)
                {
                  nLattNum = 27;
                }
              else
                {
                  nLattNum = 28;
                }
              fTemp   = pfG6[1];
              pfG6[1] = 4.0f *  pfG6[2] - pfG6[0];
              pfG6[2] = fTemp;
              pfG6[4] = -pfG6[5];
            }
          else if (fabs((pfG6[1]- pfG6[3]) / pfG6[1]) < 0.005)
            {
              if (8 <= m_nVerbose)
                cout << "g1 = g3; p. 37B" << endl;
              nLattNum = 30;
              fTemp   = pfG6[0];
              pfG6[0] = pfG6[1];
              pfG6[1] = 4.0f * pfG6[2] - pfG6[0];
              pfG6[2] = fTemp;
              pfG6[4] = -pfG6[5];
            }
          else if (fabs((pfG6[0] - pfG6[5]) / pfG6[0]) < 0.005)
            {
              if (8 <= m_nVerbose)
                cout << "g0 = g5; p. 37C" << endl;  // 10 20 30 90 92 90
              nLattNum = 29;
              pfG6[1] = 4.0f * pfG6[1] - pfG6[0];
            }
          else if (fabs((pfG6[4] - pfG6[5]) / pfG6[4]) < 0.005)
            {
              if (8 <= m_nVerbose)
                cout << "g4 = g5; p. 38C" << endl;    //
            }
          else
            {
              if (8 <= m_nVerbose)
                cout << "p.38A" << endl;
              nLattNum = 43;
              pfG6[1] = 4.0f * pfG6[1] - pfG6[0] - pfG6[2] - pfG6[4];
              pfG6[2] = pfG6[0] + pfG6[2] + pfG6[4];
              pfG6[4] = - pfG6[0] - pfG6[4];
            }
        }
      pfG6[3] = 0.0;
      pfG6[5] = 0.0;
      if (pfG6[4] > 0.0) pfG6[4] = -pfG6[4];  // Make sure beta >90!
    }
  else if (sTemp == "aP")
    {
      nLattNum = 44;  //TYPE II
      if (1 == nType)
        nLattNum = 31;  //TYPE I
    }
  else if (sTemp == "oP")
    {
      nLattNum = 32;
    }


  if (NULL != pnLattNum)
    *pnLattNum = nLattNum;
  return (0);
}

int CCellReduce::nParseSelection(Cstring& sInput,int& nEntry) {

        Cstring sTemp;
        bool bHasA = false;
        bool bHasB = false;
    int  nMonoC[2];
    int  nMonoP[2];
    int  nTriclinic;
    int  nx;
    Cstring sCent;

    nTriclinic = -1;
    nMonoC[0] = nMonoC[1] = -1;
    nMonoP[0] = nMonoP[1] = -1;
    for (nx=0;nx<m_poSolutions->nGetNumReflns();nx++) {
        sCent = (*m_poSolutions)[nx].sGetField(m_nFI_sLattice);
        if ((sCent == "mP") && (nMonoP[0] == -1))
            nMonoP[0] = nx + 1;
        else if ((sCent == "mP") && (nMonoP[1] == -1))
            nMonoP[1] =  nx + 1;
        else if ((sCent == "mC") && (nMonoC[0] == -1))
            nMonoC[0] = nx + 1;
        else if ((sCent == "mC") && (nMonoC[1] == -1))
            nMonoC[1] =  nx + 1;
        else if ((sCent == "aP") && (nTriclinic == -1))
            nTriclinic = nx + 1;
    };
    if (nMonoP[1] == -1)
        nMonoP[1] = nMonoP[0];
    if (nMonoC[1] == -1)
        nMonoC[1] = nMonoC[0];
    
    sTemp = sInput;
    if (sTemp.contains('A')) {
        sTemp = sTemp.before('A');
        bHasA = true;
    } else if (sTemp.contains('B')) {
        sTemp = sTemp.before('B');
        bHasB = true;
    };
    if (strspn(sTemp.string(),"0123456789+-")==sTemp.length())
        nEntry = atoi(sTemp.string());
    else
        return 1;
    
    if ((nEntry==12) && (bHasB))
        nEntry = nMonoC[1];
    else if (nEntry==12)
        nEntry = nMonoC[0];
    else if ((nEntry==13) && (bHasB))
        nEntry = nMonoP[1];
    else if (nEntry==13)
        nEntry = nMonoP[0];
    else if (nEntry == 14)
        nEntry = nTriclinic;
    
    return 0;
};



int
CCellReduce::nListResults(int *pnBest)
{
  int     i;
  int     nNumSoln;
  int     nLattCount,nLattCount2;
  int     nBest;
  Crefln *poSoln;
  bool    bSolnMarked = FALSE;

  // Check if solutions actually are available

  if (NULL == m_poSolutions)
    {
      cout << "There are no solutions present!" << endl;
      return (-1);
    }

  nNumSoln = m_poSolutions->nGetNumReflns();
  if (1 > nNumSoln)
    {
      cout << "There are no solutions present!" << endl;
      return (-1);
    }

  Cstring sLine
    = "=======================================================================";


  Cstring sFormat = "%5s %8s %5s %4s %15s %9s %9s %9s\n";
  Cstring sCent;
  Cstring sBravais;
  Cstring sMarkSolution = ' ';
  Cspacegroup oSpacegroup;

  printf("\nLeast-squares fit of reduced primitive cell to 44 lattice characters"
         "\nsorted on decreasing (highest to lowest) symmetry."
         "\nBest possible least-squares residual is 0.000%%, worst residual is 100.0%%."
         "\nOnly solutions with %% residuals <= %5.1lf%% are listed.\n",
         m_fResidualMax);
  printf ("%s\n", (const char *)sLine);
  printf((const char *)sFormat,
         "Soln",  "LeastSq", "Spgrp",  "Cent", "Bravais type",
         "a", "b", "c");

  printf((const char *)sFormat,
         "num", "resid(%)", "num*", "type", "Cell volume",
         "alpha", "beta", "gamma");

  printf ("%s\n", (const char *)sLine);
  
  nBest = -1;
  nLattCount = -1;
  nLattCount2 = 0;
  for (i = 0; i < nNumSoln; i++)
  {
      poSoln = m_poSolutions->poGetRefln(i);


      if ((i==0) || ((*m_poSolutions)[i].sGetField(m_nFI_sLattice) != (*m_poSolutions)[i-1].sGetField(m_nFI_sLattice))) {
          nLattCount++;
          nLattCount2 = 0;
      } else 
          nLattCount2++;

      if (   ("?" != poSoln->sGetField(m_nFI_sLattice))
          && (m_fResidualMax >= poSoln->fGetIntensity()))
      {
          // Legitimate solution
          
          if (!bSolnMarked
              && (5.0 >= poSoln->fGetIntensity()))
          {
              //              printf("%s\n\n", (const char *)sLine);
              //              sMarkSolution = '+';
              //              m_nSolnPick1 = i;
              if (-1 == nBest) nBest = i+1;
              bSolnMarked = TRUE;
          }
          sCent = poSoln->sGetField(m_nFI_sLattice);
          sCent = sCent.GetAt(1);
          oSpacegroup.vSet((int)poSoln->fGetSigmaI());
          sBravais = oSpacegroup.sGetClass();
                       
          printf(" %s%2d%c %8.3lf %5d %4s %15s %9.3lf %9.3lf %9.3lf\n",
              sMarkSolution.string(),
              nLattCount+1,
              (nLattCount2==1)?'b':' ',
              (double) poSoln->fGetIntensity(),
              (int)poSoln->fGetSigmaI(),        // Holds possible spacegroup number
              (const char *)sCent,
              (const char *)sBravais,
              (double) poSoln->fGetField(m_nFI_fA),
              (double) poSoln->fGetField(m_nFI_fB),
              (double) poSoln->fGetField(m_nFI_fC)
              );
          printf( "                                %9.0lf %9.3lf %9.3lf %9.3lf\n",
              //          printf( "                      %3d %3d  %9.0lf %9.3lf %9.3lf %9.3lf\n",
              //                 poSoln->nGetField(m_nFI_nSortOrder),
              //                 poSoln->nGetField(m_nFI_nLattNum),
              (double) poSoln->fGetField(m_nFI_fVolume),
              (double) poSoln->fGetField(m_nFI_fAlpha),
              (double) poSoln->fGetField(m_nFI_fBeta),
              (double) poSoln->fGetField(m_nFI_fGamma));
          
              /*          printf( "                                         %9.3lf %9.3lf %9.3lf\n",
              (double) poSoln->fGetField(m_nFI_fRot1),
              (double) poSoln->fGetField(m_nFI_fRot2),
              (double) poSoln->fGetField(m_nFI_fRot3));
          */
          printf("\n");
          if (bSolnMarked) sMarkSolution = ' ';
      }
  }
  printf ("%s\n", (const char *)sLine);
  printf ("*Suggested spacegroup number until systematic "
      "absences are examined.\n");
  if (bSolnMarked)
  {
      //      printf ("+Suggested solution number until symmetry "
      //              "equivalences are examined.\n");
  }
  
  fflush(stdout);
  // Return the first solution listed 
  
  if (NULL != pnBest)
      *pnBest = nBest;
  return (0);
}

int CCellReduce::nListResultsB()
{
  int     i=0;
  int     nNumSoln;
  double   fResidual;
  Crefln *poSoln;

  // Check if solutions actually are available

  if (NULL == m_poSolutions)
    {
      cout << "There are no solutions present!" << endl;
      return (-1);
    }

  nNumSoln = m_poSolutions->nGetNumReflns();
  if (1 > nNumSoln)
    {
      cout << "There are no solutions present!" << endl;
      return (-1);
    }

  Cstring sLine
    = "======================================================================";

  // Remove duplicate/similar solutions

  /*
  for (i = nNumSoln-1; i > 0; i--)
    {
      poSoln = m_poSolutions->poGetRefln(i);
      if ("?" == poSoln->sGetField(m_nFI_sLattice))
        {
          // Reflection not marked for deletion yet, so check it ...

          for (j = i-1; j >= 0; j--)
            {
              poSoln2 = m_poSolutions->poGetRefln(j);
              if (   (fCloseEnough >= fabs(poSoln->fGetField(m_nFI_fRot1)
                                          - poSoln2->fGetField(m_nFI_fRot1)))
                  && (fCloseEnough >= fabs(poSoln->fGetField(m_nFI_fRot2)
                                          - poSoln2->fGetField(m_nFI_fRot2)))
                  && (fCloseEnough >= fabs(poSoln->fGetField(m_nFI_fRot3)
                                          - poSoln2->fGetField(m_nFI_fRot3)))
                  )
                {
                  // Rotation angles are about the same so mark one for deletion

                  poSoln2->vSetField(m_nFI_sLattice, "D");
                }
            }
        }
    }
    */

  m_poSolutions->nSelect("-sLattice==D");
  m_poSolutions->nDelete((const char*)NULL, (const bool)FALSE);

  nNumSoln = m_poSolutions->nGetNumReflns();

  printf("\nUnit cell parameters and orientation angles\n");
  printf ("%s\n", (const char *)sLine);
  printf("%3s %9s %9s %9s %9s\n",
         " ", "Integer",
         "a", "b", "c" );
  printf("%3s %9s %9s %9s %9s %8s %8s %8s\n",
         "Num", "residual",
         "alpha", "beta", "gamma", "Rot1", "Rot2", "Rot3");

  printf ("%s\n", (const char *)sLine);

  for (i = 0; i < nNumSoln; i++)
    {
      poSoln    = m_poSolutions->poGetRefln(i);
      fResidual = poSoln->fGetIntensity();
      if (0.0 > fResidual)
        {
          // Convert Fourier residual so range is 0 -> 1.0
          fResidual = 1.0f + fResidual / 1000.0f;
        }
      else
        {
          // Otherwise convert residual from 1 -> 0 to 0 -> 1.0
          fResidual = 1.0f - fResidual;
        }
      printf("%3d %9.3lf %9.3lf %9.3lf %9.3lf %8.3lf %8.3lf %8.3lf\n",
             i+1,
             (double) fResidual,
             (double) poSoln->fGetField(m_nFI_fA),
             (double) poSoln->fGetField(m_nFI_fB),
             (double) poSoln->fGetField(m_nFI_fC),
             (double) poSoln->fGetField(m_nFI_fRot1),
             (double) poSoln->fGetField(m_nFI_fRot2),
             (double) poSoln->fGetField(m_nFI_fRot3));

      printf("%3s %9s %9.3lf %9.3lf %9.3lf\n",
             " ", " ",
             (double) poSoln->fGetField(m_nFI_fAlpha),
             (double) poSoln->fGetField(m_nFI_fBeta),
             (double) poSoln->fGetField(m_nFI_fGamma));
      printf("\n");
    }
  printf ("%s\n", (const char *)sLine);
  printf ("The above table shows symmetry EQUIVALENT crystal orientation angles\n"
          "for the indexing orientation.  All the solutions are equivalent for\n"
          "the selected Bravais lattice.  The default selection usually has the\n"
          "values closest to crystal orientation found in the input .head file or\n"
          "the one where (|Rot1| + |Rot2| + |Rot3|) is a minimum.\n");
  fflush(stdout);
  return (0);
};

int
CCellReduce::nListResultsL() {
    int nx,ny;
    Ccrystal oCrystal;
        itr<int> anUsed;

        for (nx=0;nx<=44;nx++)
                anUsed[nx] = 0;
    // List the 44 lattice characters.                
    printf("========================================================================\n");
    printf("Latt  LeastSq  Type  Cent    Bravais type         a         b         c\n");
    printf("num  residual        type     Cell volume     alpha      beta     gamma\n");
    printf("========================================================================\n");
    for (nx=1; ms_anOrder[nx]!=-1;nx++)  {
        ny = ms_anOrder[nx];
                if (!anUsed[ny]) {
                        anUsed[ny] = 1;
                        if (m_a44fLattResid[ny]<100.0) {
                                oCrystal.vSetCell(m_a44x6fLattCells[ny]);
                                printf("%4d %8.2lf  %4s  %4c %15s  %8.3lf  %8.3lf  %8.3lf\n"
                                        "                               %10d  %8.3lf  %8.3lf  %8.3lf\n",
                                        ny,
                                        m_a44fLattResid[ny],
                                        (m_a44tLattChar[ny].nType==1)?("I"):("II"),
                                        m_a44tLattChar[ny].sLatticeType.string()[1],
                                        m_a44tLattChar[ny].sLattSymm.string(),
                                        m_a44x6fLattCells[ny][0],m_a44x6fLattCells[ny][1],m_a44x6fLattCells[ny][2],
                                        (int) oCrystal.fCalcVolume(),
                                        m_a44x6fLattCells[ny][3],m_a44x6fLattCells[ny][4],m_a44x6fLattCells[ny][5]
                                        );
                        };
                };
    };
    printf("========================================================================\n\n");
    return 0;
};

int
CCellReduce::nSelectMatch(const int nSpacegroupNum)
{
  // Select the solution with a lattice that matches the input space group number

  int     i;
  int     nNumSoln;
  Crefln *poSoln;
  Cstring sLattice;
  Cstring sTemp;

  if (   (NULL == m_poSolutions)
      || (0 > nSpacegroupNum)
      || (230 < nSpacegroupNum) )
    return (-1);

  Cspacegroup *poSpacegroup;
  poSpacegroup = new Cspacegroup (nSpacegroupNum);
  sTemp += poSpacegroup->cGetClass();
  sTemp += poSpacegroup->cGetCentFlag();
  delete poSpacegroup;

// cP  cI cF hR hP tI tP oI oC oF oP mC mP aP //

  nNumSoln = m_poSolutions->nGetNumReflns();
  for (i = 0; i < nNumSoln; i++)
    {
      poSoln   = m_poSolutions->poGetRefln(i);
      sLattice = poSoln->sGetField(m_nFI_sLattice);
      if (sTemp == sLattice)
        {
          // We have a match
          return (i);
        }
    }
  return (-2);
}



int
CCellReduce::nLoadRotationChoices(const int nLattNumIn,Ccrystal* poCrystal)
{
  // Convert all the solutions in m_poSolutions from the reduced primitive ones
  // back to the conventional ones.  Convert the orientation angles too.
  // The solutions should all have the same unit cell dimensions, but
  // potentially different orientation angles.

  int     i;
  int     nLattNum;
  int     nNumSoln;
  double   a6fCell[6];
  float    fRot1,fRot2,fRot3;


  nLattNum = nLattNumIn;
  if (1 > nLattNum)
    nLattNum = m_nLattNum;  // Use the member variable lattice number

  if (   (NULL == m_poSolutions)
      || (1 > nLattNum)
      || (44 < nLattNum) )
    return (-1);

  double fMatIn[3][3];
  double fMatOut[24][3][3];
  double a3x3fTempMat[3][3];
  double a3x3fTempMat2[3][3];

  poCrystal->nCalcOrientMatrix();
  poCrystal->vGetOrientMatrix(&a3x3fTempMat[0][0]);
  fInvMat3D(&a3x3fTempMat[0][0],&a3x3fTempMat2[0][0]);
  vTranMat3D(a3x3fTempMat2);
 
  // First, load and obtain solutions for this orientation.
  m_poSolutions->vDeleteAll();
  vCopyMat3D(&a3x3fTempMat2[0][0],&fMatIn[0][0]);
  nNumSoln = nEquivRotations(fMatIn,fMatOut,m_a44tLattChar[nLattNum].sLatticeType); 

  Ccrystal oCrystal;

  // Invert transformation matrix

  for (i = 0; i < nNumSoln; i++)
    {
      // Load solution values into a crystal object
      vCopyMat3D(&fMatOut[i][0][0], &a3x3fTempMat[0][0]);
      vTranMat3D(a3x3fTempMat);
      fInvMat3D(&a3x3fTempMat[0][0],&a3x3fTempMat2[0][0]);
      oCrystal.nSetOrientMatrix(&a3x3fTempMat2[0][0]);
      oCrystal.vGetCell(a6fCell);
      oCrystal.vGetOrientAngles(&fRot1,&fRot2,&fRot3);

      // ... and put them back into m_poSolutions

      Crefln oSoln(m_poSolutions);
      oSoln.vSetField(m_nFI_fA,     (float)a6fCell[0]);
      oSoln.vSetField(m_nFI_fB,     (float)a6fCell[1]);
      oSoln.vSetField(m_nFI_fC,     (float)a6fCell[2]);
      oSoln.vSetField(m_nFI_fAlpha, (float)a6fCell[3]);
      oSoln.vSetField(m_nFI_fBeta,  (float)a6fCell[4]);
      oSoln.vSetField(m_nFI_fGamma, (float)a6fCell[5]);
      oSoln.vSetField(m_nFI_fRot1,  (float)fRot1);
      oSoln.vSetField(m_nFI_fRot2,  (float)fRot2);
      oSoln.vSetField(m_nFI_fRot3,  (float)fRot3);
      m_poSolutions->nInsert(&oSoln);
    }
  return (0);
}


void CCellReduce::vInitLattChar(void)
{
  // Initialize the lattice character structures from
  // Table 9.3.1 in IntTable Vol A, p. 746.
  //
  // In deference to Fortran and humans, lattchar[0] is not used, though
  // it is set to triclinic.

//  cout << "vInitLattChar called!\n";

  m_a44tLattChar[ 0].nLattChar      = 44;
  m_a44tLattChar[ 0].nSortOrder     = 45;
  m_a44tLattChar[ 0].nSpacegroupMin = 1;
  m_a44tLattChar[ 0].nType          = 2;
  m_a44tLattChar[ 0].nTypeDependent = 0;
  m_a44tLattChar[ 0].fMatFactor     = 1.0f;
  m_a44tLattChar[ 0].sLattSymm      = "triclinic";
  m_a44tLattChar[ 0].sLatticeType   = "aP";

  m_a44tLattChar[ 1].nLattChar      = 1;
  m_a44tLattChar[ 1].nSortOrder     = 3;  // was 1;  
  m_a44tLattChar[ 1].nSpacegroupMin = 196;
  m_a44tLattChar[ 1].nType          = 1;
  m_a44tLattChar[ 1].nTypeDependent = 1;
  m_a44tLattChar[ 1].fMatFactor     = 1.0f;
  m_a44tLattChar[ 1].sLattSymm      = "cubic";
  m_a44tLattChar[ 1].sLatticeType   = "cF";

  m_a44tLattChar[ 2].nLattChar      = 2;
  m_a44tLattChar[ 2].nSortOrder     = 4;
  m_a44tLattChar[ 2].nSpacegroupMin = 146;
  m_a44tLattChar[ 2].nType          = 1;
  m_a44tLattChar[ 2].nTypeDependent = 1;
  m_a44tLattChar[ 2].fMatFactor     = 1.0f;
  m_a44tLattChar[ 2].sLattSymm      = "rhomb/hexagonal";
  m_a44tLattChar[ 2].sLatticeType   = "hR";

  m_a44tLattChar[ 3].nLattChar      = 3;
  m_a44tLattChar[ 3].nSortOrder     = 1; // was 3;
  m_a44tLattChar[ 3].nSpacegroupMin = 195;
  m_a44tLattChar[ 3].nType          = 2;
  m_a44tLattChar[ 3].nTypeDependent = 0;
  m_a44tLattChar[ 3].fMatFactor     = 1.0f;
  m_a44tLattChar[ 3].sLattSymm      = "cubic";
  m_a44tLattChar[ 3].sLatticeType   = "cP";

  m_a44tLattChar[ 4].nLattChar      = 4;
  m_a44tLattChar[ 4].nSortOrder     = 5;
  m_a44tLattChar[ 4].nSpacegroupMin = 146;
  m_a44tLattChar[ 4].nType          = 2;
  m_a44tLattChar[ 4].nTypeDependent = 1;
  m_a44tLattChar[ 4].fMatFactor     = 1.0f;
  m_a44tLattChar[ 4].sLattSymm      = "rhomb/hexagonal";
  m_a44tLattChar[ 4].sLatticeType   = "hR";

  m_a44tLattChar[ 5].nLattChar      = 5;
  m_a44tLattChar[ 5].nSortOrder     = 2;
  m_a44tLattChar[ 5].nSpacegroupMin = 197;
  m_a44tLattChar[ 5].nType          = 2;
  m_a44tLattChar[ 5].nTypeDependent = 1;
  m_a44tLattChar[ 5].fMatFactor     = 1.0f;
  m_a44tLattChar[ 5].sLattSymm      = "cubic";
  m_a44tLattChar[ 5].sLatticeType   = "cI";

  m_a44tLattChar[ 6].nLattChar      = 6;
  m_a44tLattChar[ 6].nSortOrder     = 12; // was 6;
  m_a44tLattChar[ 6].nSpacegroupMin = 79;
  m_a44tLattChar[ 6].nType          = 2;
  m_a44tLattChar[ 6].nTypeDependent = 1;
  m_a44tLattChar[ 6].fMatFactor     = 1.0f;
  m_a44tLattChar[ 6].sLattSymm      = "tetragonal";
  m_a44tLattChar[ 6].sLatticeType   = "tI";

  m_a44tLattChar[ 7].nLattChar      = 7;
  m_a44tLattChar[ 7].nSortOrder     = 13; // was 7;
  m_a44tLattChar[ 7].nSpacegroupMin = 79;
  m_a44tLattChar[ 7].nType          = 2;
  m_a44tLattChar[ 7].nTypeDependent = 1;
  m_a44tLattChar[ 7].fMatFactor     = 1.0f;
  m_a44tLattChar[ 7].sLattSymm      = "tetragonal";
  m_a44tLattChar[ 7].sLatticeType   = "tI";

  m_a44tLattChar[ 8].nLattChar      = 8;
  m_a44tLattChar[ 8].nSortOrder     = 22; // 8;
  m_a44tLattChar[ 8].nSpacegroupMin = 23;
  m_a44tLattChar[ 8].nType          = 2;
  m_a44tLattChar[ 8].nTypeDependent = 1;
  m_a44tLattChar[ 8].fMatFactor     = 1.0f;
  m_a44tLattChar[ 8].sLattSymm      = "orthorhombic";
  m_a44tLattChar[ 8].sLatticeType   = "oI";

  m_a44tLattChar[ 9].nLattChar      = 9;
  m_a44tLattChar[ 9].nSortOrder     = 6; // 9;
  m_a44tLattChar[ 9].nSpacegroupMin = 146;
  m_a44tLattChar[ 9].nType          = 1;
  m_a44tLattChar[ 9].nTypeDependent = 1;
  m_a44tLattChar[ 9].fMatFactor     = 1.0f;
  m_a44tLattChar[ 9].sLattSymm      = "rhomb/hexagonal";
  m_a44tLattChar[ 9].sLatticeType   = "hR";

  m_a44tLattChar[10].nLattChar      = 10;
  m_a44tLattChar[10].nSortOrder     = 30; // 10;
  m_a44tLattChar[10].nSpacegroupMin = 5;
  m_a44tLattChar[10].nType          = 1;
  m_a44tLattChar[10].nTypeDependent = 1;
  m_a44tLattChar[10].fMatFactor     = 1.0f;
  m_a44tLattChar[10].sLattSymm      = "monoclinic";
  m_a44tLattChar[10].sLatticeType   = "mC";

  m_a44tLattChar[11].nLattChar      = 11;
  m_a44tLattChar[11].nSortOrder     = 10; // 11;
  m_a44tLattChar[11].nSpacegroupMin = 75;
  m_a44tLattChar[11].nType          = 2;
  m_a44tLattChar[11].nTypeDependent = 0;
  m_a44tLattChar[11].fMatFactor     = 1.0f;
  m_a44tLattChar[11].sLattSymm      = "tetragonal";
  m_a44tLattChar[11].sLatticeType   = "tP";

  m_a44tLattChar[12].nLattChar      = 12;
  m_a44tLattChar[12].nSortOrder     = 8; // 12;
  m_a44tLattChar[12].nSpacegroupMin = 143;
  m_a44tLattChar[12].nType          = 2;
  m_a44tLattChar[12].nTypeDependent = 1;
  m_a44tLattChar[12].fMatFactor     = 1.0f;
  m_a44tLattChar[12].sLattSymm      = "trig/hexagonal";
  m_a44tLattChar[12].sLatticeType   = "hP";

  m_a44tLattChar[13].nLattChar      = 13;
  m_a44tLattChar[13].nSortOrder     = 17; // 13;
  m_a44tLattChar[13].nSpacegroupMin = 21;
  m_a44tLattChar[13].nType          = 2;
  m_a44tLattChar[13].nTypeDependent = 1;
  m_a44tLattChar[13].fMatFactor     = 1.0f;
  m_a44tLattChar[13].sLattSymm      = "orthorhombic";
  m_a44tLattChar[13].sLatticeType   = "oC";

  m_a44tLattChar[14].nLattChar      = 14;
  m_a44tLattChar[14].nSortOrder     = 31; // 16;
  m_a44tLattChar[14].nSpacegroupMin = 5;
  m_a44tLattChar[14].nType          = 2;
  m_a44tLattChar[14].nTypeDependent = 1;
  m_a44tLattChar[14].fMatFactor     = 1.0f;
  m_a44tLattChar[14].sLattSymm      = "monoclinic";
  m_a44tLattChar[14].sLatticeType   = "mC";

  m_a44tLattChar[15].nLattChar      = 15;
  m_a44tLattChar[15].nSortOrder     = 14;
  m_a44tLattChar[15].nSpacegroupMin = 79;
  m_a44tLattChar[15].nType          = 2;
  m_a44tLattChar[15].nTypeDependent = 1;
  m_a44tLattChar[15].fMatFactor     = 1.0f;
  m_a44tLattChar[15].sLattSymm      = "tetragonal";
  m_a44tLattChar[15].sLatticeType   = "tI";

  m_a44tLattChar[16].nLattChar      = 16;
  m_a44tLattChar[16].nSortOrder     = 25; // 15;
  m_a44tLattChar[16].nSpacegroupMin = 22;
  m_a44tLattChar[16].nType          = 2;
  m_a44tLattChar[16].nTypeDependent = 1;
  m_a44tLattChar[16].fMatFactor     = 1.0f;
  m_a44tLattChar[16].sLattSymm      = "orthorhombic";
  m_a44tLattChar[16].sLatticeType   = "oF";

  m_a44tLattChar[17].nLattChar      = 17;
  m_a44tLattChar[17].nSortOrder     = 32; //17;
  m_a44tLattChar[17].nSpacegroupMin = 5;
  m_a44tLattChar[17].nType          = 2;
  m_a44tLattChar[17].nTypeDependent = 1;
  m_a44tLattChar[17].fMatFactor     = 1.0f;
  m_a44tLattChar[17].nType          = 2;
  m_a44tLattChar[17].sLattSymm      = "monoclinic";
  m_a44tLattChar[17].sLatticeType   = "mC";

  m_a44tLattChar[18].nLattChar      = 18;
  m_a44tLattChar[18].nSortOrder     = 15;
  m_a44tLattChar[18].nSpacegroupMin = 79;
  m_a44tLattChar[18].nType          = 1;
  m_a44tLattChar[18].nTypeDependent = 1;
  m_a44tLattChar[18].fMatFactor     = 1.0f;
  m_a44tLattChar[18].sLattSymm      = "tetragonal";
  m_a44tLattChar[18].sLatticeType   = "tI";

  m_a44tLattChar[19].nLattChar      = 19;
  m_a44tLattChar[19].nSortOrder     = 23;
  m_a44tLattChar[19].nSpacegroupMin = 23;
  m_a44tLattChar[19].nType          = 1;
  m_a44tLattChar[19].nTypeDependent = 1;
  m_a44tLattChar[19].fMatFactor     = 1.0f;
  m_a44tLattChar[19].sLattSymm      = "orthorhombic";
  m_a44tLattChar[19].sLatticeType   = "oI";

  m_a44tLattChar[20].nLattChar      = 20;
  m_a44tLattChar[20].nSortOrder     = 33;
  m_a44tLattChar[20].nSpacegroupMin = 5;
  m_a44tLattChar[20].nType          = 1;
  m_a44tLattChar[20].nTypeDependent = 1;
  m_a44tLattChar[20].fMatFactor     = 1.0f;
  m_a44tLattChar[20].sLattSymm      = "monoclinic";
  m_a44tLattChar[20].sLatticeType   = "mC";

  m_a44tLattChar[21].nLattChar      = 21;
  m_a44tLattChar[21].nSortOrder     = 11;
  m_a44tLattChar[21].nSpacegroupMin = 75;
  m_a44tLattChar[21].nType          = 2;
  m_a44tLattChar[21].nTypeDependent = 0;
  m_a44tLattChar[21].fMatFactor     = 1.0f;
  m_a44tLattChar[21].sLattSymm      = "tetragonal";
  m_a44tLattChar[21].sLatticeType   = "tP";

  m_a44tLattChar[22].nLattChar      = 22;
  m_a44tLattChar[22].nSortOrder     = 9;
  m_a44tLattChar[22].nSpacegroupMin = 143;
  m_a44tLattChar[22].nType          = 2;
  m_a44tLattChar[22].nTypeDependent = 1;
  m_a44tLattChar[22].fMatFactor     = 1.0f;
  m_a44tLattChar[22].sLattSymm      = "trig/hexagonal";
  m_a44tLattChar[22].sLatticeType   = "hP";

  m_a44tLattChar[23].nLattChar      = 23;
  m_a44tLattChar[23].nSortOrder     = 18;
  m_a44tLattChar[23].nSpacegroupMin = 21;
  m_a44tLattChar[23].nType          = 2;
  m_a44tLattChar[23].nTypeDependent = 1;
  m_a44tLattChar[23].fMatFactor     = 1.0f;
  m_a44tLattChar[23].sLattSymm      = "orthorhombic";
  m_a44tLattChar[23].sLatticeType   = "oC";

  m_a44tLattChar[24].nLattChar      = 24;
  m_a44tLattChar[24].nSortOrder     = 7;
  m_a44tLattChar[24].nSpacegroupMin = 146;
  m_a44tLattChar[24].nType          = 2;
  m_a44tLattChar[24].nTypeDependent = 1;
  m_a44tLattChar[24].fMatFactor     = 1.0f;
  m_a44tLattChar[24].sLattSymm      = "rhomb/hexagonal";
  m_a44tLattChar[24].sLatticeType   = "hR";

  m_a44tLattChar[25].nLattChar      = 25;
  m_a44tLattChar[25].nSortOrder     = 34;
  m_a44tLattChar[25].nSpacegroupMin = 5;
  m_a44tLattChar[25].nType          = 2;
  m_a44tLattChar[25].nTypeDependent = 1;
  m_a44tLattChar[25].fMatFactor     = 1.0f;
  m_a44tLattChar[25].sLattSymm      = "monoclinic";
  m_a44tLattChar[25].sLatticeType   = "mC";

  m_a44tLattChar[26].nLattChar      = 26;
  m_a44tLattChar[26].nSortOrder     = 26;
  m_a44tLattChar[26].nSpacegroupMin = 22;
  m_a44tLattChar[26].nType          = 1;
  m_a44tLattChar[26].nTypeDependent = 1;
  m_a44tLattChar[26].fMatFactor     = 1.0f;
  m_a44tLattChar[26].sLattSymm      = "orthorhombic";
  m_a44tLattChar[26].sLatticeType   = "oF";

  m_a44tLattChar[27].nLattChar      = 27;
  m_a44tLattChar[27].nSortOrder     = 35;
  m_a44tLattChar[27].nSpacegroupMin = 5;
  m_a44tLattChar[27].nType          = 1;
  m_a44tLattChar[27].nTypeDependent = 1;
  m_a44tLattChar[27].fMatFactor     = 1.0f;
  m_a44tLattChar[27].sLattSymm      = "monoclinic";
  m_a44tLattChar[27].sLatticeType   = "mC";

  m_a44tLattChar[28].nLattChar      = 28;
  m_a44tLattChar[28].nSortOrder     = 36;
  m_a44tLattChar[28].nSpacegroupMin = 5;
  m_a44tLattChar[28].nType          = 1;
  m_a44tLattChar[28].nTypeDependent = 1;
  m_a44tLattChar[28].fMatFactor     = 1.0f;
  m_a44tLattChar[28].sLattSymm      = "monoclinic";
  m_a44tLattChar[28].sLatticeType   = "mC";

  m_a44tLattChar[29].nLattChar      = 29;
  m_a44tLattChar[29].nSortOrder     = 37;
  m_a44tLattChar[29].nSpacegroupMin = 5;
  m_a44tLattChar[29].nType          = 1;
  m_a44tLattChar[29].nTypeDependent = 1;
  m_a44tLattChar[29].fMatFactor     = 1.0f;
  m_a44tLattChar[29].sLattSymm      = "monoclinic";
  m_a44tLattChar[29].sLatticeType   = "mC";

  m_a44tLattChar[30].nLattChar      = 30;
  m_a44tLattChar[30].nSortOrder     = 38;
  m_a44tLattChar[30].nSpacegroupMin = 5;
  m_a44tLattChar[30].nType          = 1;
  m_a44tLattChar[30].nTypeDependent = 1;
  m_a44tLattChar[30].fMatFactor     = 1.0f;
  m_a44tLattChar[30].sLattSymm      = "monoclinic";
  m_a44tLattChar[30].sLatticeType   = "mC";

  m_a44tLattChar[31].nLattChar      = 31;
  m_a44tLattChar[31].nSortOrder     = 43;
  m_a44tLattChar[31].nSpacegroupMin = 1;
  m_a44tLattChar[31].nType          = 1;
  m_a44tLattChar[31].nTypeDependent = 1;
  m_a44tLattChar[31].fMatFactor     = 1.0f;
  m_a44tLattChar[31].sLattSymm      = "triclinic";
  m_a44tLattChar[31].sLatticeType   = "aP";

  m_a44tLattChar[32].nLattChar      = 32;
  m_a44tLattChar[32].nSortOrder     = 16;
  m_a44tLattChar[32].nSpacegroupMin = 16;
  m_a44tLattChar[32].nType          = 2;
  m_a44tLattChar[32].nTypeDependent = 0;
  m_a44tLattChar[32].fMatFactor     = 1.0f;
  m_a44tLattChar[32].sLattSymm      = "orthorhombic";
  m_a44tLattChar[32].sLatticeType   = "oP";

  m_a44tLattChar[33].nLattChar      = 33;
  m_a44tLattChar[33].nSortOrder     = 28;
  m_a44tLattChar[33].nSpacegroupMin = 3;
  m_a44tLattChar[33].nType          = 2;
  m_a44tLattChar[33].nTypeDependent = 1;
  m_a44tLattChar[33].fMatFactor     = 1.0f;
  m_a44tLattChar[33].sLattSymm      = "monoclinic";
  m_a44tLattChar[33].sLatticeType   = "mP";

  m_a44tLattChar[34].nLattChar      = 34;
  m_a44tLattChar[34].nSortOrder     = 29;
  m_a44tLattChar[34].nSpacegroupMin = 3;
  m_a44tLattChar[34].nType          = 2;
  m_a44tLattChar[34].nTypeDependent = 1;
  m_a44tLattChar[34].fMatFactor     = 1.0f;
  m_a44tLattChar[34].sLattSymm      = "monoclinic";
  m_a44tLattChar[34].sLatticeType   = "mP";

  m_a44tLattChar[35].nLattChar      = 35;
  m_a44tLattChar[35].nSortOrder     = 27;
  m_a44tLattChar[35].nSpacegroupMin = 3;
  m_a44tLattChar[35].nType          = 2;
  m_a44tLattChar[35].nTypeDependent = 1;
  m_a44tLattChar[35].fMatFactor     = 1.0f;
  m_a44tLattChar[35].sLattSymm      = "monoclinic";
  m_a44tLattChar[35].sLatticeType   = "mP";

  m_a44tLattChar[36].nLattChar      = 36;
  m_a44tLattChar[36].nSortOrder     = 20;
  m_a44tLattChar[36].nSpacegroupMin = 21;
  m_a44tLattChar[36].nType          = 2;
  m_a44tLattChar[36].nTypeDependent = 1;
  m_a44tLattChar[36].fMatFactor     = 1.0f;
  m_a44tLattChar[36].sLattSymm      = "orthorhombic";
  m_a44tLattChar[36].sLatticeType   = "oC";

  m_a44tLattChar[37].nLattChar      = 37;
  m_a44tLattChar[37].nSortOrder     = 40;
  m_a44tLattChar[37].nSpacegroupMin = 5;
  m_a44tLattChar[37].nType          = 2;
  m_a44tLattChar[37].nTypeDependent = 1;
  m_a44tLattChar[37].fMatFactor     = 1.0f;
  m_a44tLattChar[37].sLattSymm      = "monoclinic";
  m_a44tLattChar[37].sLatticeType   = "mC";

  m_a44tLattChar[38].nLattChar      = 38;
  m_a44tLattChar[38].nSortOrder     = 21;
  m_a44tLattChar[38].nSpacegroupMin = 21;
  m_a44tLattChar[38].nType          = 2;
  m_a44tLattChar[38].nTypeDependent = 1;
  m_a44tLattChar[38].fMatFactor     = 1.0f;
  m_a44tLattChar[38].sLattSymm      = "orthorhombic";
  m_a44tLattChar[38].sLatticeType   = "oC";

  m_a44tLattChar[39].nLattChar      = 39;
  m_a44tLattChar[39].nSortOrder     = 41;
  m_a44tLattChar[39].nSpacegroupMin = 5;
  m_a44tLattChar[39].nType          = 2;
  m_a44tLattChar[39].nTypeDependent = 1;
  m_a44tLattChar[39].fMatFactor     = 1.0f;
  m_a44tLattChar[39].sLattSymm      = "monoclinic";
  m_a44tLattChar[39].sLatticeType   = "mC";

  m_a44tLattChar[40].nLattChar      = 40;
  m_a44tLattChar[40].nSortOrder     = 19;
  m_a44tLattChar[40].nSpacegroupMin = 21;
  m_a44tLattChar[40].nType          = 2;
  m_a44tLattChar[40].nTypeDependent = 1;
  m_a44tLattChar[40].fMatFactor     = 1.0f;
  m_a44tLattChar[40].sLattSymm      = "orthorhombic";
  m_a44tLattChar[40].sLatticeType   = "oC";

  m_a44tLattChar[41].nLattChar      = 41;
  m_a44tLattChar[41].nSortOrder     = 39;
  m_a44tLattChar[41].nSpacegroupMin = 5;
  m_a44tLattChar[41].nType          = 2;
  m_a44tLattChar[41].nTypeDependent = 1;
  m_a44tLattChar[41].fMatFactor     = 1.0f;
  m_a44tLattChar[41].sLattSymm      = "monoclinic";
  m_a44tLattChar[41].sLatticeType   = "mC";

  m_a44tLattChar[42].nLattChar      = 42;
  m_a44tLattChar[42].nSortOrder     = 24;
  m_a44tLattChar[42].nSpacegroupMin = 23;
  m_a44tLattChar[42].nType          = 2;
  m_a44tLattChar[42].nTypeDependent = 1;
  m_a44tLattChar[42].fMatFactor     = 1.0f;
  m_a44tLattChar[42].nType          = 2;
  m_a44tLattChar[42].sLattSymm      = "orthorhombic";
  m_a44tLattChar[42].sLatticeType   = "oI";

  m_a44tLattChar[43].nLattChar      = 43;
  m_a44tLattChar[43].nSortOrder     = 42;
  m_a44tLattChar[43].nSpacegroupMin = 5;
  m_a44tLattChar[43].nType          = 2;
  m_a44tLattChar[43].nTypeDependent = 1;
  m_a44tLattChar[43].fMatFactor     = 1.0f;
  m_a44tLattChar[43].sLattSymm      = "monoclinic";
  m_a44tLattChar[43].sLatticeType   = "mC";

  m_a44tLattChar[44].nLattChar      = 44;
  m_a44tLattChar[44].nSortOrder     = 44;
  m_a44tLattChar[44].nSpacegroupMin = 1;
  m_a44tLattChar[44].nType          = 2;
  m_a44tLattChar[44].nTypeDependent = 1;
  m_a44tLattChar[44].fMatFactor     = 1.0f;
  m_a44tLattChar[44].sLattSymm      = "triclinic";
  m_a44tLattChar[44].sLatticeType   = "aP";

  // Build the transformation matrices needed.  There may be errors
  // in the International Tables.

  int i;
  double a3x3fTMat[3][3];
  double a6x6fMat[6][6];

  for (i = 0; i < 45; i++)
    {
      switch (m_a44tLattChar[i].nLattChar)
        {
        case 1:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  1.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] = -1.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  1.0;
          break;

        case 2:
        case 4:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  1.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] = -1.0;
          a3x3fTMat[2][2] = -1.0;
          break;

        case 5:
        case 7:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  1.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  1.0;
          break;

        case 6:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  1.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  1.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 8:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] = -1.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] = -1.0;
          a3x3fTMat[2][2] = -1.0;
          break;

        case 9:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] = -1.0;
          a3x3fTMat[2][2] =  3.0;
          break;
          
        case 10:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] = -1.0;
          break;
          
        case 13:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  1.0;
          break;
          
        case 14:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  1.0;
          break;
          
        case 15:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  2.0;
          break;
          
        case 16:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  2.0;
          break;

        case 17:                  // Changed from IT, see Kabsch (1993) J Appl Cryst 26, 798
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] = -1.0;
          break;

        case 18:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  1.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] = -1.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;
          
        case 19:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  1.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  1.0;
          break;

        case 20:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  1.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] = -1.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 21:
        case 22:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  1.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;
          
    case 23:
        case 25:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  1.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  1.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;
          
        case 24:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  2.0;
          a3x3fTMat[2][0] =  1.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  1.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;
          
        case 26:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  2.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  2.0;
          break;
          
        case 27:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] =  2.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] = -1.0;
          a3x3fTMat[2][2] =  1.0;
          break;

        case 28:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  2.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 29:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] = -2.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] = -1.0;
          break;

        case 30:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] =  1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] = -2.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 34:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] = -1.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] = -1.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 35:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] = -1.0;
          break;

        case 36:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] = -2.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  0.0;
          break;
          
        case 37:
          a3x3fTMat[0][0] =  1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  2.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 38:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  1.0;
          a3x3fTMat[1][1] =  2.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  1.0;
          break;

        case 39:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] = -2.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] =  0.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] = -1.0;
          break;

        case 40:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] =  1.0;
          a3x3fTMat[2][1] =  2.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 41:
          a3x3fTMat[0][0] =  0.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] = -2.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] = -1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
          break;

        case 42:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] =  0.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] =  0.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  1.0;
          a3x3fTMat[2][2] =  2.0;
          break;

        case 43:
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] =  0.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] = -2.0;
          a3x3fTMat[0][2] =  0.0;
          a3x3fTMat[1][2] = -1.0;
          a3x3fTMat[2][2] =  0.0;
/*
          a3x3fTMat[0][0] = -1.0;
          a3x3fTMat[1][0] = -1.0;
          a3x3fTMat[2][0] =  0.0;
          a3x3fTMat[0][1] = -1.0;
          a3x3fTMat[1][1] = -1.0;
          a3x3fTMat[2][1] = -2.0;
          a3x3fTMat[0][2] =  1.0;
          a3x3fTMat[1][2] =  0.0;
          a3x3fTMat[2][2] =  0.0;
*/
//Note:
//       Christer Svenson suggests the matrix for 43 should be:
//      -1 -1 0 /  -1 -1 -2 / 1 0 0 (in Fortran)

          break;

        case 0:
        case 3:
        case 11:
        case 12:
        case 31:
        case 32:
        case 33:
        case 44:
        default:
          // Identity
          vIdentMat3D(a3x3fTMat);
          break;

        } // end of switch(nLattNum)

      // Copy from temporarory matrix to the member variable

      vCopyMat3D(&a3x3fTMat[0][0], &m_a44tLattChar[ i].a3x3fTMat[0][0]);

      // Initialize the Least Squares 6x6 matrices of Paciorek & Bonin (1992)
      
      /*

        The following lines give the projection spaces for each of the 44 lattices.  The projection matrix
        is computed from these.


      */
        char* cpId = NULL;

      vZeroMat(6, 6, &a6x6fMat[0][0]);
      switch (m_a44tLattChar[i].nLattChar) {
      case 3:
        cpId = "r 1 , r 1 , r 1 , , , ,";
        break;
      case 5:
        cpId = "r 3 , r 3 , r 3 , r -2  , r -2 , r -2 ,";
        break;
      case 1:
        cpId = "r 1 , r 1 , r 1 , r 1 , r 1 , r 1 ,";
        break;
      case 11:
        cpId = "r 1 , r 1 , s 1 , , , , ";
        break;
      case 21:
        cpId = "r 1 , s 1 , s 1 , , , , ";
        break;
      case 15:
        cpId = "r 1 , r 1 , s 1 , r -1 , r -1 , , ";
        break;
      case 6:
        cpId = "r 1 , r 1 , r 1 , r -1 s -1 , r -1 s -1 , s 2 , ";
        break;
      case 7:
        cpId = "r 1 , r 1 , r 1 , s 2 , r -1 s -1 , r -1  s -1 , ";
        break;
      case 18:
        cpId = "r 2 , s 2 , s 2 , r 1 , r 2 , r 2 , ";
        break;
      case 12:
        cpId = "r 1 , r 1 , s 1 , , , r -1 , ";
        break;
      case 22:
        cpId = "r 1 , s 1 , s 1 , s -1 , , , ";
        break;
      case 9:
        cpId = "r 1 , r 1 , s 1 , r 1 , r 1 , r 1 , ";
        break;
      case 2:
        cpId = "r 1 , r 1 , r 1 , s 1 , s 1 , s 1 , ";
        break;
      case 4:
        cpId = "r 1 , r 1 , r 1 , s 1 , s 1 , s 1 , ";
        break;
      case 24:
        cpId = "r 3 , s 1 , s 1 , s -1 r 1 , r -2 , r -2 , ";
        break;
      case 32:
        cpId = "r 1 , s 1 , t 1 , , , , ";
        break;
      case 36:
        cpId = "r 1 , s 1 , t 1 , , r -1 , , ";
        break;
      case 38:
        cpId = "r 1 , s 1 , t 1 , , , r -1 , ";
        break;
      case 13:
        cpId = "r 1 , r 1 , s 1 , , , t 1 , ";
        break;
      case 23:
        cpId = "r 1 , s 1 , s 1 , t 1 , , , ";
        break;
      case 40:
        cpId = "r 1 , s 1 , t 1 , s -1 , , , ";
        break;
      case 16:
        cpId = "r 1 , r 1 , s 1 , t 1 , t 1  , r -2 t -2 , ";
        break;
      case 26:
        cpId = "r 2 , s 2 , t 2 , r 1 , r 2 , r 2 , ";
        break;
      case 8:
        cpId = "r 1 , r 1 , r 1 , s 1 , t 1 , r -2 s -1 t -1 , ";
        break;
      case 19:
        cpId = "r 1 , s 1 , s 1 , t 1 , r 1 , r 1 , ";
        break;
      case 42:
        cpId = "r 1 , s 1 , t 1 , s -1 , r -1 , , ";
        break;
      case 33:
        cpId = "r 1 , s 1 , t 1 , , u 1 , , ";
        break;
      case 35:
        cpId = "r 1 , s 1 , t 1 , u 1 , , , ";
        break;
      case 34:
        cpId = "r 1 , s 1 , t 1 , , , u 1 , ";
        break;
      case 39:
        cpId = "r 1 , s 1 , t 1 , u 1 , , r -1 , ";
        break;
      case 41:
        cpId = "r 1 , s 1 , t 1 , s -1 , u 1 , , ";
        break;
      case 37:
        cpId = "r 1 , s 1 , t 1 , u 1 , r -1 , , ";
        break;
      case 10:
        cpId = "r 1 , r 1 , s 1 , t 1 , t 1 , u 1 , ";
        break;
      case 14:
        cpId = "r 1 , r 1 , s 1 , t 1 , t 1 , u 1 , ";
        break;
      case 20:
        cpId = "r 1 , s 1 , s 1 , t 1 , u 1 , u 1 , ";
        break;
      case 25:
        cpId = "r 1 , s 1 , s 1 , t 1 , u 1 , u 1 , ";
        break;
      case 28:
        cpId = "r 1 , s 1 , t 1 , u 1 , r 1 , u 2 , ";
        break;
      case 30:
        cpId = "r 1 , s 1 , t 1 , s 1 , u 1 , u 2 , ";
        break;
      case 29:
        cpId = "r 1 , s 1 , t 1 , u 1 , u 2 , r 1 , ";
        break;
      case 43:
        cpId = "r 1 , s 1 , t 1 , s -1 u -1 , r -1 u -1 , u 2 , ";
        break;
      case 17:
        cpId = "r 1 , r 1 , s 1 , t 1 , u 1  , r -2 t -1 u -1 , ";
        break;
      case 27:
        cpId = "r 1 , s 1 , t 1 , u 1 , r 1 , r 1 , ";
        break;
      case 31:
        cpId = "r 1 , s 1 , t 1 , u 1 , v 1 , w 1 , ";
        break;
      case 0:
      case 44:
        cpId = "r 1 , s 1 , t 1 , u 1 , v 1 , w 1 , ";
        break;
        };

        int nx,ny;
        double a6x6fProjVecs[6][6];
        double a6x6fProjVecsTran[6][6];
        double a6fVec[6];
        char* cpChars="rstuvw";
        char cpBuf[60];
        char* cp;
        double f0,f1;

        vZeroMat(6, 6, &a6x6fProjVecs[0][0]);

        /*  In this loop, nx loops over each variable (r,s,t,u,v, or w)
            We build the projection vector in f6Vec, and then place it in the correct col. of a6x6ProjVecs.
            after proper normalization.
            Then, Compute P^T*P matrix.
        */
      
        for (nx=0; nx<6; nx++) {
            
            for (ny=0; ny<6; ny++) 
                a6fVec[ny]=0.0;
            
            // Copy over the string because strtok puts '\x0' in some places during parse.
            strcpy(cpBuf,cpId);

            for (ny = 0, f1 = 0,cp = strtok(cpBuf," \t"); cp ; cp=strtok(NULL," \t")) {
                if (cp[0]==',') {
                    ny++; 
                } else
                if (cp[0]==cpChars[nx]) {
                    if (1 != sscanf( (cp = strtok(NULL, " \t")),"%lf", &f0)) {
                        cout << "Logic Error";
                        exit(0);
                    };
                    a6fVec[ny]=f0;
                    f1+=f0*f0;
                } else 
                    cp=strtok(NULL," \t");
            };
            if (ny != 6) {
                cout << "Logic Error";
                exit(0);
            };
            // We have a (non-normalized) vector representing the 6 dimensional space spanned by the letter. (r,s,t,u,v, or w)
            if (f1 != 0.0) {
                // Grahm Schmidt with previous vectors.  This gives orthogonal basis.
                for (ny = 0; ny < nx; ny++) {
                    double a6fTemp[6];
                    vCopyVecND(6,&(a6x6fProjVecs[ny][0]),&(a6fTemp[0]));
                    f0=fDotND(6,&(a6fVec[0]),&(a6fTemp[0]));
                    vMulVecNDScalar(6,&(a6fTemp[0]),f0,&(a6fTemp[0]));
                    vSubVecNDVecND(6,&(a6fVec[0]),&(a6fTemp[0]),&(a6fVec[0]));
                };
                // Normalize the vector.
                vMulVecNDScalar(6,a6fVec,1.0/fLenVecND(6,a6fVec),a6fVec);
            };
            // Copy Vector Over.
            vCopyVecND(6, &(a6fVec[0]), &(a6x6fProjVecs[nx][0]));
        };
        // Compute projection matrix.
        vCopyVecND(6*6, &a6x6fProjVecs[0][0], & a6x6fProjVecsTran[0][0]);
        vTranMatND(6,&(a6x6fProjVecsTran[0][0]));
        vMulMatNDMatND(6,&(a6x6fProjVecs[0][0]),&(a6x6fProjVecsTran[0][0]),&(m_a44tLattChar[ i].a6x6fLSMat[0][0]));
        m_a44tLattChar[i].fMatFactor =1.0;

      

    }  // i loop end

  // Initialize matrices to convert a reduced primitive cell to the 24 other
  // possibly nearly Buerger-reduced cells.
  // i = 0, i=25, are the identity, while i = 1...24 are the operations
  // in Table 3 of Andrews and Berstein.

  for (i = 0; i < 26; i++)
    {
      vZeroMat(6, 6, &a6x6fMat[0][0]);
              
      switch (i)
        {
        case 1:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][2]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;
          
        case 2:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][3]  = -1.0;
          a6x6fMat[4][2]  = -1.0;
          a6x6fMat[4][4]  = -1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;
                  
        case 3:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = 2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;
                  
        case 4:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = 2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][2]  = -1.0;
          a6x6fMat[3][3]  = -1.0;
          a6x6fMat[4][4]  = -1.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;
                  
        case 5:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = 2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][2]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;
                  
        case 6:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;
                  
        case 7:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[3][3]  = -1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][2]  = -1.0;
          a6x6fMat[5][4]  = -1.0;
          break;
                  
        case 8:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = 2.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;

        case 9:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = -2.0;
          a6x6fMat[3][2]  = -1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][4]  = -1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;
                  
        case 10:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = 2.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][2]  = 1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;

        case 11:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][5]  = 2.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = 2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[5][1]  = 1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;

        case 12:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][3]  = 2.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[0][5]  = 2.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = 2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][2]  = 1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[5][1]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][3]  = 2.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;

        case 13:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][5]  = 2.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = -2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][2]  = -1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][3]  = -1.0;
          a6x6fMat[4][4]  = -1.0;
          a6x6fMat[5][1]  = -1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = -1.0;
          break;

        case 14:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][3]  = 2.0;
          a6x6fMat[0][4]  = 2.0;
          a6x6fMat[0][5]  = 2.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = 2.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][1]  = 1.0;
          a6x6fMat[4][2]  = 1.0;
          a6x6fMat[4][3]  = 2.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;
                  
        case 15:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][5]  = 2.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = 2.0;
          a6x6fMat[3][2]  = -1.0;
          a6x6fMat[3][3]  = -1.0;
          a6x6fMat[4][1]  = -1.0;
          a6x6fMat[4][3]  = -1.0;
          a6x6fMat[4][4]  = -1.0;
          a6x6fMat[4][5]  = -1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;
                  
        case 16:
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[1][0]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][4]  = 2.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[3][5]  = 1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;
          
        case 17:
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[1][0]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][4]  = -2.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[3][3]  = -1.0;
          a6x6fMat[3][5]  = 1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[5][2]  = -1.0;
          a6x6fMat[5][4]  = 1.0;
          break;

        case 18:
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[1][0]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = -2.0;
          a6x6fMat[3][4]  = -1.0;
          a6x6fMat[3][5]  = 1.0;
          a6x6fMat[4][2]  = -1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[5][4]  = 1.0;
          break;

        case 19:
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][3]  = 2.0;
          a6x6fMat[1][0]  = 1.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = 2.0;
          a6x6fMat[1][4]  = 2.0;
          a6x6fMat[1][5]  = 2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[3][4]  = 1.0;
          a6x6fMat[4][2]  = 1.0;
          a6x6fMat[4][3]  = 1.0;
          a6x6fMat[5][1]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][3]  = 2.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;

        case 20:
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][3]  = 2.0;
          a6x6fMat[1][0]  = 1.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][5]  = -2.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[3][4]  = -1.0;
          a6x6fMat[4][2]  = -1.0;
          a6x6fMat[4][3]  = -1.0;
          a6x6fMat[5][1]  = -1.0;
          a6x6fMat[5][3]  = -1.0;
          a6x6fMat[5][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;

        case 21:
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[1][0]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][4]  = 2.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = 2.0;
          a6x6fMat[3][2]  = -1.0;
          a6x6fMat[3][3]  = -1.0;
          a6x6fMat[3][4]  = -1.0;
          a6x6fMat[3][5]  = -1.0;
          a6x6fMat[4][1]  = -1.0;
          a6x6fMat[4][3]  = -1.0;
          a6x6fMat[5][3]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;

        case 22:
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][3]  = 2.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[2][0]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][4]  = -2.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[3][5]  = -1.0;
          a6x6fMat[4][2]  = -1.0;
          a6x6fMat[4][3]  = -1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][1]  = -1.0;
          a6x6fMat[5][3]  = -1.0;
          break;

        case 23:
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[1][3]  = -2.0;
          a6x6fMat[2][0]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][4]  = -2.0;
          a6x6fMat[3][2]  = -1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[3][4]  = 1.0;
          a6x6fMat[3][5]  = -1.0;
          a6x6fMat[4][3]  = -1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][1]  = -1.0;
          a6x6fMat[5][3]  = 1.0;
          break;

        case 24:
          a6x6fMat[0][1]  = 1.0;
          a6x6fMat[0][2]  = 1.0;
          a6x6fMat[0][3]  = 2.0;
          a6x6fMat[1][2]  = 1.0;
          a6x6fMat[2][0]  = 1.0;
          a6x6fMat[2][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[2][3]  = 2.0;
          a6x6fMat[2][4]  = 2.0;
          a6x6fMat[2][5]  = 2.0;
          a6x6fMat[3][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[3][4]  = 1.0;
          a6x6fMat[4][1]  = 1.0;
          a6x6fMat[4][2]  = 1.0;
          a6x6fMat[4][3]  = 2.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[4][5]  = 1.0;
          a6x6fMat[5][2]  = 1.0;
          a6x6fMat[5][3]  = 1.0;
          break;
          
        case 0:
        case 25:
        default:
          a6x6fMat[0][0]  = 1.0;
          a6x6fMat[1][1]  = 1.0;
          a6x6fMat[2][2]  = 1.0;
          a6x6fMat[3][3]  = 1.0;
          a6x6fMat[4][4]  = 1.0;
          a6x6fMat[5][5]  = 1.0;
          break;
        } // end of switch(i)

      vCopyVecND(6*6, &a6x6fMat[0][0], &m_a26x6x6fBMat[i][0][0]);
    }
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

/* nEquivRotations:  Given a cell, it loads all equivalent cells based on the Bravis lattice type.

    5 --- 7
   /|    /|
  / |   / |     a  
 6 --- 8  |     | 
 |  |  |  |     |
 |  1 -|- 3     +- c 
 | /   | /     /
 |/    |/     b
 2 --- 4

   Interpretation of table:
   1,2,3 symbolize the a,b,and c axes respectively.  a negative sign implies inversion of a direction.
   Thus, the direction {-2,-3,1} is the Right handed coordianate system defined by (-b,-c,a).  
   Note that all triples should be right handed.

   My reasoning for these equivalent coordinate systems are as follows:

  uP,aP:  No other representation.
  mP,mC:  We can negate a and c to get from (1) to (7)
  oP,oC,oF,oI: (1)=(4)=(6)=(7).  Negating every pair of 2 directions for 4 total.
  hR (rhombohedral Basis):  Opposite vertices (1)=(8), with three right handed selections at both (1) and (8).
                            We must flip a pair of vertices for eacho of the three (8) systems, since (-1,-2,-3) is not right handed.
  hR (hexagonal Basis) See hP. 
  hP:  (1)=(6) and (3)=(8), but we must flip a and b with both (3) and (8) since these are not right handed.
  tP,tI:  (1)=(4)=(6)=(7) and (2)=(3)=(5)=(8) but we need to flip a and b with (2),(3),(5),(8)
  cP,cI,cF:  All possible right handed coordinate systems.  24 total.

  */


#define C_TYPES 14
/* For (Hexagonal) */     
int nRotEquivsEquiv_Ct[C_TYPES+1]={    24,24,24,4,4,8,8,4,4,4,4,2,2,1,1};
// /* For (Rhombohedral) */  int nEquiv_Ct[C_TYPES]={ 24,24,24,6,4,8,8,4,4,4,4,2,2,1,1};
int nRotEquivs[C_TYPES+1][24][3]={
   /* cP */   {{1,2,3},{3,1,2},{2,3,1},   {-1,-2,3},{3,-1,-2},{-2,3,-1},   {1,-2,-3},{-3,1,-2},{-2,-3,1},   {-1,2,-3},{-3,-1,2},{2,-3,-1},    {2,-1,3},{3,2,-1},{-1,3,2},   {-2,1,3},{3,-2,1},{1,3,-2},  {2,1,-3},{-3,2,1},{1,-3,2},   {-2,-1,-3},{-3,-2,-1},{-1,-3,-2}},
   /* cI */   {{1,2,3},{3,1,2},{2,3,1},   {-1,-2,3},{3,-1,-2},{-2,3,-1},   {1,-2,-3},{-3,1,-2},{-2,-3,1},   {-1,2,-3},{-3,-1,2},{2,-3,-1},    {2,-1,3},{3,2,-1},{-1,3,2},   {-2,1,3},{3,-2,1},{1,3,-2},  {2,1,-3},{-3,2,1},{1,-3,2},   {-2,-1,-3},{-3,-2,-1},{-1,-3,-2}},
   /* cF */   {{1,2,3},{3,1,2},{2,3,1},   {-1,-2,3},{3,-1,-2},{-2,3,-1},   {1,-2,-3},{-3,1,-2},{-2,-3,1},   {-1,2,-3},{-3,-1,2},{2,-3,-1},    {2,-1,3},{3,2,-1},{-1,3,2},   {-2,1,3},{3,-2,1},{1,3,-2},  {2,1,-3},{-3,2,1},{1,-3,2},   {-2,-1,-3},{-3,-2,-1},{-1,-3,-2}},

    // /* hR (Rhombohedral) */   {{1,2,3},{3,1,2},{2,3,1}, {-1,-3,-2},{-2,-1,-3},{-3,-2,-1}},
    
   /* hR (Hexagonal)*/  {{1,2,3},{-1,-2,3},{2,1,-3},{-2,-1,-3}},
   
   /* hP */   {{1,2,3},{-1,-2,3},{2,1,-3},{-2,-1,-3}},
   /* tP */   {{1,2,3},{-1,-2,3},{1,-2,-3},{-1,2,-3},  {2,-1,3},{-2,1,3},{2,1,-3},{-2,-1,-3}},
   /* tI */   {{1,2,3},{-1,-2,3},{1,-2,-3},{-1,2,-3},  {2,-1,3},{-2,1,3},{2,1,-3},{-2,-1,-3}},
   /* oP */   {{1,2,3},{-1,-2,3},{1,-2,-3},{-1,2,-3}},
   /* oC */   {{1,2,3},{-1,-2,3},{1,-2,-3},{-1,2,-3}},
   /* oI */   {{1,2,3},{-1,-2,3},{1,-2,-3},{-1,2,-3}},
   /* oF */   {{1,2,3},{-1,-2,3},{1,-2,-3},{-1,2,-3}},
   /* mP */   {{1,2,3},{-1,2,-3}},
   /* mC */   {{1,2,3},{-1,2,-3}},
   /* aP */   {{1,2,3}},
   /* uP */   {{1,2,3}}
};
char* CCellReduce::ms_cpRotEquivsStrings [C_TYPES+1] = {
    "cP",
    "cI",
    "cF",
    "hR",
    "hP",
    "tP",
    "tI",
    "oP",
    "oC",
    "oI",
    "oF",
    "mP",
    "mC",
    "aP",
    "uP"
    };

int CCellReduce::ms_pnLatticeChars[C_TYPES+1] =  { 
        44,  
        33,  
        10,  
        32,  
        13,  
        16,  
        8,   
        4,   
        12,  
        11,  
        6,   
        3,   
        5,   
        1,
        1
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CCellReduce::nEquivRotations(double fMatIn[3][3], double fMatOut[24][3][3], Cstring sLattID)
{
    // RB 01/05/2005 Since the approach, used here by Thad, just does not work for rhombohedral cell in hexagonal settings,
    // we need to do the hR case separately. 
    if( sLattID == "hR" )
        return nEquivRotationsRhombohedralHexagonal(fMatIn, fMatOut);
    
    int      nx=0;
    int      ny=0;
    int      nz=0;

    int      nType = 0;

    for(nType = 0; nType < C_TYPES; nType++) 
    {
        if( sLattID == ms_cpRotEquivsStrings[nType] )
            break;
    }

    if( nType == C_TYPES ) // lattice ID not found in ms_cpRotEquivsStrings 
        return 0;

    for(nx = 0; nx < nRotEquivsEquiv_Ct[nType]; nx++) 
    {
        for(ny = 0; ny < 3; ny++)
        {
            nz = nRotEquivs[nType][nx][ny];

            vCopyVec3D(fMatIn[abs(nz)-1], fMatOut[nx][ny]);
            
            vMulVec3DScalar(fMatOut[nx][ny], nz < 0 ? -1 : 1, fMatOut[nx][ny]);        
        }
    }

    return nx;  // at this point nx == nRotEquivsEquiv_Ct[nType]
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CCellReduce::nEquivRotationsRhombohedralHexagonal(double fMatIn[3][3], double fMatOut[24][3][3])
{
    const int   nNumberOfEquivalentCells = 6;
    const double a3x3f_Rhombohedral_Hexagonal_Equivalent_Cell_Transforms[nNumberOfEquivalentCells][3][3] =
                                                                                   {1.0,   0.0,   0.0,     // a, b, c
                                                                                    0.0,   1.0,   0.0,
                                                                                    0.0,   0.0,   1.0,

                                                                                    0.0,   1.0,   0.0,     // b, -a-b, c
                                                                                   -1.0,  -1.0,   0.0,
                                                                                    0.0,   0.0,   1.0,

                                                                                   -1.0,  -1.0,   0.0,     // -a-b, a, c
                                                                                    1.0,   0.0,   0.0,
                                                                                    0.0,   0.0,   1.0,

                                                                                    0.0,   1.0,   0.0,     // b, a, -c
                                                                                    1.0,   0.0,   0.0,
                                                                                    0.0,   0.0,  -1.0,

                                                                                   -1.0,  -1.0,   0.0,     // -a-b, b, -c
                                                                                    0.0,   1.0,   0.0,
                                                                                    0.0,   0.0,  -1.0,

                                                                                    1.0,   0.0,   0.0,     // a, -a-b, -c
                                                                                   -1.0,  -1.0,   0.0,
                                                                                    0.0,   0.0,  -1.0
                                                                                    };
    double      a3x3f_Temp_1[3][3] = {0.0};
    double      a3x3f_Temp_2[3][3] = {0.0};
                                                                                    
    for(int ii=0; ii < nNumberOfEquivalentCells; ii++)
    {
        vCopyMat3D(&a3x3f_Rhombohedral_Hexagonal_Equivalent_Cell_Transforms[ii][0][0], &a3x3f_Temp_1[0][0]); 
        
        vMulMat3DMat3D(fMatIn, a3x3f_Temp_1, a3x3f_Temp_2);
        
        vCopyMat3D(&a3x3f_Temp_2[0][0], &fMatOut[ii][0][0]); 
    }

    return nNumberOfEquivalentCells;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RB 01/05/2005 When the solution, loaded into CCellReduce is hR, we have a problem with obverse vs reverse setting 
// when CCellReduce::fFindOrigVecs() is called. The function does not know it is dealing with a rhombohedral lattice, 
// so there is no guarantee the result will be in obverse setting, as required by IT.
//
// The fix is as follows.
// 
// If the loaded solution is rhombohedral:
// 1. Calculate the primitive cell parameters.  
// 2. Calculate a3x3fInCentered based on that primitive cell.
// 3. Calculate the real OM by calling fFindOrigVecs() as the code normally does.
// 4. Transform the resultant OM to an obverse hexagonal setting.

double CCellReduce::fFindOrigVecsA(double a3x3fInCentered[3][3], 
			                       double a3x3fInPreReduced[3][3],
                                   double a3x3fOutCoords[3][3],
                                   char*  pcLatticeType)
{
    double      a3x3f_Temp[3][3] = {0.0};
    vCopyMat3D(&a3x3fInCentered[0][0], &a3x3f_Temp[0][0]); 
    
    ///////////////////////////////////////////////////////////////////////////////////
    // Determine if it's "hR" lattice type
    
    // The lattice type is either passed in or embedded in the CCellReduce object
    Cstring     sLatticeType("");
    
    if( pcLatticeType )
        sLatticeType = pcLatticeType;
    else if( m_nLattNum > 0 && m_nLattNum <= 44 )
        sLatticeType = m_a44tLattChar[m_nLattNum].sLatticeType;

    bool    bIsRhombohedralSolution = sLatticeType == "hR";

    if( bIsRhombohedralSolution )
    {
        double          f_A_Hexagonal = a3x3fInCentered[0][0] * a3x3fInCentered[0][0] +
                                        a3x3fInCentered[0][1] * a3x3fInCentered[0][1] +
                                        a3x3fInCentered[0][2] * a3x3fInCentered[0][2];
        f_A_Hexagonal = sqrt(f_A_Hexagonal);
        
        
        ////////////////////////////////////////////////////////////////////////////////
        // RB For debugging calculate B as well. B is assumed to be equal A.
        double          f_B_Hexagonal = a3x3fInCentered[1][0] * a3x3fInCentered[1][0] +
                                        a3x3fInCentered[1][1] * a3x3fInCentered[1][1] +
                                        a3x3fInCentered[1][2] * a3x3fInCentered[1][2];
        f_B_Hexagonal = sqrt(f_B_Hexagonal);
        ////////////////////////////////////////////////////////////////////////////////


        double          f_C_Hexagonal = a3x3fInCentered[2][0] * a3x3fInCentered[2][0] +
                                        a3x3fInCentered[2][1] * a3x3fInCentered[2][1] +
                                        a3x3fInCentered[2][2] * a3x3fInCentered[2][2];
        f_C_Hexagonal = sqrt(f_C_Hexagonal);
        
        // Calculate rhombohedral cell parameters in primitive setting
        double          f_A_Hex_squared = f_A_Hexagonal * f_A_Hexagonal;
        double          f_C_Hex_squared = f_C_Hexagonal * f_C_Hexagonal;
        double          fTemp1 = 3.0 * f_A_Hex_squared + f_C_Hex_squared;

        double          f_A_Primitive = sqrt(fTemp1) / 3.0;

        double          f_Alpha_Primitive = (f_C_Hex_squared - 1.5 * f_A_Hex_squared) / fTemp1;
        f_Alpha_Primitive = acos(f_Alpha_Primitive) / acos(-1.0) * 180.0;

        Ccrystal        oRhombohedralPrimitiveCrystal;
        oRhombohedralPrimitiveCrystal.vSetCell(f_A_Primitive, 
                                               f_A_Primitive,
                                               f_A_Primitive,
                                               f_Alpha_Primitive,
                                               f_Alpha_Primitive,
                                               f_Alpha_Primitive);
    
        oRhombohedralPrimitiveCrystal.vSetOrientAngles(0.0, 0.0, 0.0); // these don't matter by Thad's design

            oRhombohedralPrimitiveCrystal.nCalcGetRealMatrix(&a3x3f_Temp[0][0]);
        
        vTranMat3D(a3x3f_Temp);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    double      dRet = fFindOrigVecs(a3x3f_Temp, a3x3fInPreReduced, a3x3fOutCoords);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( bIsRhombohedralSolution )
    {
        double      a3x3f_Primitive_To_Hexagonal[3][3] = { 1.0, -1.0,  0.0,    // R1 obverse
                                                           0.0,  1.0, -1.0,
                                                           1.0,  1.0,  1.0 };
    
        vMulMat3DMat3D(a3x3fOutCoords, a3x3f_Primitive_To_Hexagonal, a3x3f_Temp);
        vCopyMat3D(&a3x3f_Temp[0][0], &a3x3fOutCoords[0][0]); 
    }
    
    return dRet;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
double CCellReduce::fFindOrigVecs(double a3x3fInCentered[3][3], 
                                              double a3x3fInPreReduced[3][3],
                                  double a3x3fOutCoords[3][3])
{
  int   nx, ny;
  double f0, f1;
  double a3fVec[3];
  double a3fVecT[3];

  int  nMaxScan = 4;
  enum {eBestCt = 20};
  int  nBestCt = (int) eBestCt;
  double afBestVec[3][eBestCt][3];    // [nVec][nBestCt][vec]
  double afBest[3][eBestCt];          // [nVec][nBestCt]

  int   a3nTest[3];
  int   nVec;
  int   nWorst;
  double fVecLen;
  double fDetInputSign;
  double a3x3fInPreReducedOrig[3][3];



  // Calculate the sign of the input determinant.  
  fDetInputSign=fDetMat3D(& a3x3fInCentered[0][0]);
  if (fDetInputSign<0)
      fDetInputSign=-1;
  else
      fDetInputSign=1;

  // Do a cell reduction on the pre-reduced vectors to make sure that they
  // are as small as possible.

  int
    anLCPre[][3]={{0,1,1},{0,1,-1},{1,0,1},{1,0,-1},
                  {0,2,1},{0,2,-1},{2,0,1},{2,0,-1},
                  {1,2,1},{1,2,-1},{2,1,1},{2,1,-1}};
  int nReductionLoopCt = 0;
  bool bReductionmade;

  
  // Save the original pre-reduced cell.
  vCopyMat3D(&a3x3fInPreReduced[0][0],&a3x3fInPreReducedOrig[0][0]);

  
  do 
  {
      bReductionmade = FALSE;
      for (nx = 0; nx<12; nx++)
          {
                  do 
                  {
                          f0 = fabs(fDot3D(a3x3fInPreReduced[anLCPre[nx][1]], 
                                  a3x3fInPreReduced[anLCPre[nx][1]]));
                          vCopyVec3D(a3x3fInPreReduced[anLCPre[nx][0]], a3fVecT); 
                          vMulVec3DScalar(a3fVecT, anLCPre[nx][2], a3fVecT);
                          vAddVec3DVec3D(a3fVecT, a3x3fInPreReduced[anLCPre[nx][1]], 
                                  a3fVecT);
                          f1 = fabs(fDot3D(a3fVecT, a3fVecT));
                          if (f0 > f1) 
                          { 
                                  vCopyVec3D(a3fVecT, a3x3fInPreReduced[anLCPre[nx][1]]);
                                  bReductionmade = TRUE;
                                  nReductionLoopCt++; 
                          }
                  } while (f0 > f1);      
          }
  } while (bReductionmade);
  

  // In this loop, a3nTest gives the linear combination numbers.
  // Do loop for all three centered basis vectors.

  for (nVec = 0; nVec < 3; nVec++)
    {
      // Initialize variables for this centered basis vector.

      fVecLen = fLenVec3D(a3x3fInCentered[nVec]);
      for (nx = 0; nx < nBestCt; nx++) 
        afBest[nVec][nx]=1e20;
      nWorst = 0;
      for (a3nTest[0] = -nMaxScan; a3nTest[0] <= nMaxScan; a3nTest[0]++)
        {
          for (a3nTest[1] = -nMaxScan; a3nTest[1] <= nMaxScan; a3nTest[1]++)
            {
              for (a3nTest[2] = -nMaxScan; a3nTest[2] <=nMaxScan; a3nTest[2]++)
                {
                  if ((a3nTest[0]) || (a3nTest[1]) || (a3nTest[2]))
                    {
                      for (nx = 0; nx < 3; nx++)
                        a3fVec[nx] = 0.0f;
                      for (nx = 0; nx < 3; nx++)
                        { 
                          vCopyVec3D(a3x3fInPreReduced[nx], a3fVecT); 
                          vMulVec3DScalar(a3fVecT, (double) a3nTest[nx], a3fVecT);
                          vAddVec3DVec3D(a3fVec, a3fVecT, a3fVec);
                        }

                      // See if this is a better LC than the 'worst' one. If not, then continue

                      f0 = fabs(fLenVec3D(a3fVec) - fVecLen);
                      if (f0 > afBest[nVec][nWorst]) 
                        continue;

                      afBest[nVec][nWorst] = f0;
                      vCopyVec3D(a3fVec, afBestVec[nVec][nWorst]);

                      // Find the best vector that is the 'worst' for the next pass.

                      for (nx = 0, nWorst = 0; nx < nBestCt; nx++)
                        {
                          if (afBest[nVec][nx] > afBest[nVec][nWorst])
                            nWorst = nx;
                        }
                    }
                }
            }
        }
    }

  // In this loop, a3nTest gives the coordinates of the solution w.r.t. nBest
  // and afBestVec
  // Next, we try linear combinations until we get one that fits the cell in
  // the best possible way.

  double    fBestErr = 1e20;
  double    fBestDiffFromIdent = 0.0;
  double    fErrTol = 0.001;
  double    a3x3fTMat[3][3];
  double    a3x3fTMatInv[3][3];
  double    a6fTempCell[6];
  double    fTempMat[3][3];
  double    fTempMatInv[3][3];
  double    a6fCenteredCell[6];
  double    fLenCenteredCell;        // Used in % error calculation.
  double    fError = 0.0;
  double    fDiffFromIdent;
  int      nSign;
  bool     bUse;

  // Compute the cell parameters for the centered cell.

  a6fCenteredCell[0] = fLenVec3D(a3x3fInCentered[0]);
  a6fCenteredCell[1] = fLenVec3D(a3x3fInCentered[1]);
  a6fCenteredCell[2] = fLenVec3D(a3x3fInCentered[2]);

  a6fCenteredCell[3] = acos(min(1.0,max(-1.0, fDot3D(a3x3fInCentered[1], a3x3fInCentered[2])
                                      / a6fCenteredCell[1] / a6fCenteredCell[2])))
                     / Gs_dRADIANS_PER_DEGREE;

  a6fCenteredCell[4] = acos(min(1.0,max(-1.0, fDot3D(a3x3fInCentered[0], a3x3fInCentered[2])
                                      / a6fCenteredCell[0] / a6fCenteredCell[2])))
                     / Gs_dRADIANS_PER_DEGREE;

  a6fCenteredCell[5] = acos(min(1.0,max(-1.0, fDot3D(a3x3fInCentered[0], a3x3fInCentered[1])
                                      / a6fCenteredCell[0] / a6fCenteredCell[1])))
                     / Gs_dRADIANS_PER_DEGREE;

  fLenCenteredCell= fLenVecND(6,& a6fCenteredCell[0]);

  for (a3nTest[0] = 0; a3nTest[0] < nBestCt; a3nTest[0]++) 
    {
      for (a3nTest[1] = 0; a3nTest[1] < nBestCt; a3nTest[1]++) 
        {
          for (a3nTest[2] = 0; a3nTest[2] < nBestCt; a3nTest[2]++) 
            {
              if (   (afBest[0][a3nTest[0]] < 1e20) 
                  && (afBest[1][a3nTest[1]] < 1e20)
                  && (afBest[2][a3nTest[2]] < 1e20)) 
                {
                  for (nSign = 0; nSign < 8; nSign++) 
                    {
                      vCopyVec3D(afBestVec[0][a3nTest[0]], fTempMat[0]);
                      vCopyVec3D(afBestVec[1][a3nTest[1]], fTempMat[1]);
                      vCopyVec3D(afBestVec[2][a3nTest[2]], fTempMat[2]);

                      // Multiply by Sign factor

                      if (nSign % 2)
                        vMulVec3DScalar(fTempMat[0], -1.0f, fTempMat[0]);
                      if ((nSign >> 1) % 2) 
                        vMulVec3DScalar(fTempMat[1], -1.0f, fTempMat[1]);
                      if ((nSign >> 2) % 2) 
                        vMulVec3DScalar(fTempMat[2], -1.0f, fTempMat[2]);

                      // Check for a non-zero determinant.

                      if (fabs(fDetMat3D(&(fTempMat[0][0]))) < 0.001f) 
                        continue; // Determinant is close to 0, so skip 
                      
                      // Evaluate cell parameters.

                      a6fTempCell[0] = fLenVec3D(fTempMat[0]);
                      a6fTempCell[1] = fLenVec3D(fTempMat[1]);
                      a6fTempCell[2] = fLenVec3D(fTempMat[2]);

                      a6fTempCell[3] = acos(min(1.0,max(-1.0,fDot3D(fTempMat[1],fTempMat[2])
                                                    /(a6fTempCell[1]*a6fTempCell[2])))) / Gs_dRADIANS_PER_DEGREE;

                      a6fTempCell[4] = acos(min(1.0,max(-1.0,fDot3D(fTempMat[0],fTempMat[2])
                                                    /(a6fTempCell[0]*a6fTempCell[2])))) / Gs_dRADIANS_PER_DEGREE;

                      a6fTempCell[5] = acos(min(1.0,max(-1.0,fDot3D(fTempMat[0],fTempMat[1])
                                                    /(a6fTempCell[0]*a6fTempCell[1])))) / Gs_dRADIANS_PER_DEGREE;

                      // Look for error between these cell parameters, and the ones that
                      // are already calculated for the centered cell.

                      for (nx = 0, fError = 0.0; nx < 6; nx++)
                        {
                          fError += fabs(a6fCenteredCell[nx] - a6fTempCell[nx])*fabs(a6fCenteredCell[nx] - a6fTempCell[nx]);
                        }
                      fError = sqrt(fError) / fLenCenteredCell;

                      
                      // Figure out the relative transformation matrix
                                          // [T] = inv([New])*[Old]  so [Old] * [T] = [New]  
                      fInvMat3D(&fTempMat[0][0],&fTempMatInv[0][0]);
                      vMulMat3DMat3D(fTempMatInv,a3x3fInPreReducedOrig,a3x3fTMat);
                                          fInvMat3D(&a3x3fTMat[0][0],&a3x3fTMatInv[0][0]);
                      fDiffFromIdent = 0.0;
                      for (nx=0;nx<3;nx++) {
                          for (ny=0;ny<3;ny++) {
                              fDiffFromIdent += ABS(a3x3fTMat[nx][ny] - (nx==ny));
                          };
                      };

                      
                                          if ((fabs(fError - fBestErr) < fErrTol) || (fError < fBestErr))
                      {
                          // Make sure that the determinant is okay.
                          f1=fDetMat3D(& fTempMat[0][0]);
                          if (0 > f1)
                              f1 = -1;
                          else
                              f1 = 1;
                          if (f1 == fDetInputSign) {
                              
                              if (ABS(fError - fBestErr) < fErrTol)
                                  bUse = (fDiffFromIdent < fBestDiffFromIdent);
                              else
                                  bUse = true;
                              if (bUse) {
                                  
                                  vCopyMat3D(&(fTempMat[0][0]), &(a3x3fOutCoords[0][0]));
                                  fBestErr = fError;
                                  fBestDiffFromIdent = fDiffFromIdent;
                              };
                          };
                      }                        
                  }
              }
            }
        }
    }
  if (fBestErr>=1e10)
    return -1.0;
  else
    return (fBestErr);
}


double
CCellReduce::fGetSolutionResidual(const int nSolnNum) {
    return (*m_poSolutions)[nSolnNum].fGetIntensity();
};

// Load a given unit cell and unit parameters from the m_poSolutions object into a Ccrystal object.
int CCellReduce::nLoadSolution(const int nSolnNum, Ccrystal *poCrystalIn)
{
    /////////////////////////////////////////////////////////////////////////////////
    int     nRef = m_poSolutions->nGetNumReflns();
    
    int     nSoln = 0;
    
    if( 0 > nSolnNum )
        nSoln = nRef + nSolnNum;
    else
        nSoln = nSolnNum;

    if( nSoln >= nRef || 0 > nSoln )
    {
        cout << "Error in nLoadSolution, solution " << nRef << " out of bounds!\n";
        return (-1);
    }
    /////////////////////////////////////////////////////////////////////////////////

    Ccrystal*   poCrystal = NULL;
    if( NULL == poCrystalIn )
    {
        poCrystal = &m_oCrystal;
    }
    else
    {
        poCrystal = poCrystalIn;
    }
    
    nRef = nSoln;
    
    Crefln*     poSoln = m_poSolutions->poGetRefln(nRef);
    
    m_nLattNum = poSoln->nGetField(m_nFI_nLattNum);
    
    poCrystal->vSetCell(poSoln->fGetField(m_nFI_fA),
                        poSoln->fGetField(m_nFI_fB),
                        poSoln->fGetField(m_nFI_fC),
                        poSoln->fGetField(m_nFI_fAlpha),
                        poSoln->fGetField(m_nFI_fBeta),
                        poSoln->fGetField(m_nFI_fGamma));
    
    poCrystal->vSetOrientAngles(poSoln->fGetField(m_nFI_fRot1),
                                poSoln->fGetField(m_nFI_fRot2),
                                poSoln->fGetField(m_nFI_fRot3));
   
    if( !m_bKnownCell )
    {
        int     nSG = (int)poSoln->fGetSigmaI();
        
        if( 1 > nSG ) 
            nSG = 1;
        
        poCrystal->m_poSpacegroup->vSet(nSG);  // Could be bogus number!
    }
    else
    {
        poCrystal->m_poSpacegroup->vSet(m_oKnown.m_poSpacegroup->nGet());
    }
    
    int     nStat = poCrystal->nCalcOrientMatrix(m_fWavelength);
    if( 0 == nStat )
    {
        poCrystal->vGetOrientMatrix(&m_a3x3fUBMat[0][0]);
        double      fDet = fInvMat3D(&m_a3x3fUBMat[0][0], &m_a3x3fInvUBMat[0][0]);
        if( 0.0 == fDet )
            nStat = -1;
    }
    
    //  cout << "nLoadSolns: crystal is:\n";
    //  (void)poCrystal->nList(1);
    return nStat;
}

int
CCellReduce::nAddSolutionL(const int nLattice) {
    if (m_a44fLattResid[nLattice]==100.0)
        return -1;
    Crefln oRefln(m_poSolutions);
    
    
    // nLattice contains the lattice number
    oRefln.vSetField(m_nFI_nLattNum,nLattice);
    oRefln.vSetField(m_nFI_fA,(float) m_a44x6fLattCells[nLattice][0]);
    oRefln.vSetField(m_nFI_fB,(float) m_a44x6fLattCells[nLattice][1]);
    oRefln.vSetField(m_nFI_fC,(float) m_a44x6fLattCells[nLattice][2]);
    oRefln.vSetField(m_nFI_fAlpha,(float) m_a44x6fLattCells[nLattice][3]);
    oRefln.vSetField(m_nFI_fBeta,(float) m_a44x6fLattCells[nLattice][4]);
    oRefln.vSetField(m_nFI_fGamma,(float) m_a44x6fLattCells[nLattice][5]);
    oRefln.vSetField(m_nFI_fRot1,0.0f);
    oRefln.vSetField(m_nFI_fRot2,0.0f);
    oRefln.vSetField(m_nFI_fRot3,0.0f);
    oRefln.vSetSigmaI((float)m_a44tLattChar[nLattice].nSpacegroupMin);  // The SigmaI value contains the spacegroup number.
    m_poSolutions->nInsert(&oRefln);
    return m_poSolutions->nGetNumReflns();    
};


int
CCellReduce::nAddCrystal(Ccrystal& oCrystalIn) {
    Crefln oRefln(m_poSolutions);

    oRefln.vSetField(m_nFI_fA,     (float)oCrystalIn.fGetCell(0));
    oRefln.vSetField(m_nFI_fB,     (float)oCrystalIn.fGetCell(1));
    oRefln.vSetField(m_nFI_fC,     (float)oCrystalIn.fGetCell(2));
    oRefln.vSetField(m_nFI_fAlpha, (float)oCrystalIn.fGetCell(3));
    oRefln.vSetField(m_nFI_fBeta,  (float)oCrystalIn.fGetCell(4));
    oRefln.vSetField(m_nFI_fGamma, (float)oCrystalIn.fGetCell(5));
    oRefln.vSetField(m_nFI_fRot1,0.0f);
    oRefln.vSetField(m_nFI_fRot2,0.0f);
    oRefln.vSetField(m_nFI_fRot3,0.0f);
    m_poSolutions->nInsert(&oRefln);
    return 0;
};

int
CCellReduce::nAddCrystal(double* a3x3fOrientMatIn) {
    Ccrystal     oCrystal;

    oCrystal.nSetOrientMatrix(a3x3fOrientMatIn,m_fWavelength);
    return nAddCrystal(oCrystal);
};

int CCellReduce::nGetOrientInfo(Ccrystal& oCrystalIn,Ccrystal& oCrystalOut) {
    fGetOrientInfo(oCrystalIn,oCrystalOut);
    return 0;
};

double CCellReduce::fGetOrientInfo(Ccrystal& oCrystalIn,Ccrystal& oCrystalOut,Ccrystal* poCrystalPrev)
{
    double        a3x3fOrientMatIn[3][3];
    double        a3x3fRealOrientMatIn[3][3];
    double        a3x3fRealMatOut[3][3];
    double        a3x3fRealOrientMatOut[3][3];
    double        a3x3fOrientMatOut[3][3];
    double        a6fDesiredCellParams[6];
    double        f0;


    oCrystalOut.vGetCell(&a6fDesiredCellParams[0]);
    oCrystalIn.nCalcOrientMatrix();
    oCrystalIn.vGetOrientMatrix(&a3x3fOrientMatIn[0][0]);
    fInvMat3D(&a3x3fOrientMatIn[0][0],&a3x3fRealOrientMatIn[0][0]);
    vTranMat3D(a3x3fRealOrientMatIn);

    // Get the real orientation matrix of the new cell.
    oCrystalOut.nCalcGetRealMatrix(&a3x3fRealMatOut[0][0]);
    vTranMat3D(a3x3fRealMatOut);
    
    // Find the real orientation matrix.
    f0 = fFindOrigVecsA(a3x3fRealMatOut, a3x3fRealOrientMatIn, a3x3fRealOrientMatOut);

        if (poCrystalPrev) {
                // We know the lattice, so we will try to find the 'best' answer.
                // The 'best' answer will be the one that has a reindexing matrix close to identity (see below)

        Cstring       sLattice;
                double            a3x3fEquivRealOrientMatOut[24][3][3];         // These are returned by the nEquivRotations routine.
                int                       nNumEquivOrientations;
                int                       nEquivCt;
                double            a3x3fEquivOrientMatOut[3][3];
        double        a3x3fPrevOrientMat[3][3];
        double        fSum;
        double        fBestSum;
        char          a2cLattice[3];
        int           nx,ny;



        poCrystalPrev->nCalcOrientMatrix();
        poCrystalPrev->vGetOrientMatrix(&a3x3fPrevOrientMat[0][0]);
        a2cLattice[0]=oCrystalOut.m_poSpacegroup->cGetClass();
        a2cLattice[1]=g_cpSpaceGroupNames[oCrystalOut.m_poSpacegroup->nGet()-1][0];
        a2cLattice[2] = 0;
        sLattice = a2cLattice;

        nNumEquivOrientations = nEquivRotations(a3x3fRealOrientMatOut,a3x3fEquivRealOrientMatOut,sLattice);
        fBestSum = 1000.0;

                for (nEquivCt=0; nEquivCt <nNumEquivOrientations; nEquivCt ++) {
            vTranMat3D(a3x3fEquivRealOrientMatOut[nEquivCt]);
            
            // The orientation matrix is in reciprocal space.
            fInvMat3D(&a3x3fEquivRealOrientMatOut[nEquivCt][0][0],&a3x3fEquivOrientMatOut[0][0]);
            fSum = 0.0;
            for (nx=0;nx<3;nx++) {
                for (ny = 0; ny < 3; ny++)
                    fSum += ABS(a3x3fPrevOrientMat[nx][ny] - a3x3fEquivOrientMatOut[nx][ny]);
            };
            if ((nEquivCt == 0) || (fSum < fBestSum)) {
                vCopyMat3D(&a3x3fEquivOrientMatOut[0][0],&a3x3fOrientMatOut[0][0]);
                fBestSum = fSum;
            };
                };
    } else {
        vTranMat3D(a3x3fRealOrientMatOut);            
        // The orientation matrix is in reciprocal space.
        fInvMat3D(&a3x3fRealOrientMatOut[0][0],&a3x3fOrientMatOut[0][0]);
    };
    
    oCrystalOut.nSetOrientMatrix(&a3x3fOrientMatOut[0][0]);
    oCrystalOut.vSetCell(&a6fDesiredCellParams[0]);
    oCrystalOut.nCalcOrientMatrix();

    return f0;
};

Cstring CCellReduce::sGetLattice(const int nSolnNum,bool bLongName) { 
    if (bLongName) {
        Cspacegroup oSpacegroup;
        oSpacegroup.vSet((int)(*m_poSolutions)[nSolnNum].fGetSigmaI());
        return oSpacegroup.sGetClass();
    } else
        return (*m_poSolutions)[nSolnNum].sGetField(m_nFI_sLattice); 
};


void CCellReduce::vConvertToPrimitive(Ccrystal& oCrystalIn,char* pcLattice) {
        
        int   pnLatticeChars[14] =  { 44,  33,  10,  32,  13,  16,  8,   4,   12,  11,  6,   3,   5,   1};
        char* pcValidLattices[14] = { "aP","mP","mC","oP","oC","oF","oI","hR","hP","tP","tI","cP","cI","cF"};
        double a3x3fOrientCenteredMatIn[3][3];
        double a3x3fTempMat1[3][3];
        double a3x3fTempMat2[3][3];
        int nx;
        char a2cLattice[3];
        
        if ((pcLattice) && (pcLattice[0])) {
                a2cLattice[0] = pcLattice[0];
                a2cLattice[1] = pcLattice[1];
                a2cLattice[2] = 0;
        } else {
                a2cLattice[0]=oCrystalIn.m_poSpacegroup->cGetClass();
                a2cLattice[1]=g_cpSpaceGroupNames[oCrystalIn.m_poSpacegroup->nGet()-1][0];
                a2cLattice[2] = 0;
        };
        
        
        // Keep a centered version of the input UB matrix.         
        oCrystalIn.nCalcOrientMatrix();
        oCrystalIn.vGetOrientMatrix(&a3x3fOrientCenteredMatIn[0][0]);
        
        // Reduce to primative if centered cell is input.       
        for (nx=0; nx<14; nx++) {
                if (!strcmp(a2cLattice,pcValidLattices[nx]))
                        break;
        };
        
        if (nx != 14) {
                // The transformation matrices in the international tables move us from primative to centered for
                // a real cell.  
        // That is, [T] left multiplies the real cell matrix (stored in rows) [T] * [RealCell].
        // Thus, tran([T]) right multiplies the reciprocal space matrix (stored in collumns)
        // That is, tran(inv([Uc][Bc])) = tran(inv([Up][Bp])) * tran([T])
                // So,  inv([Uc][Bc]) = [T] * inv([Up][Bp])
                // So,  [Uc][Bc] = [Up][Bp] * inv([T])
                // So,  [Uc][Bc] * [T] = [Up][Bp]
                        
                vCopyMat3D(&m_a44tLattChar[pnLatticeChars[nx]].a3x3fTMat[0][0],&a3x3fTempMat1[0][0]);
                
        //TJN:  don't transpose!  vTranMat3D(a3x3fTempMat1);
                vMulMat3DMat3D(a3x3fOrientCenteredMatIn,a3x3fTempMat1,a3x3fTempMat2);
                oCrystalIn.nSetOrientMatrix(&a3x3fTempMat2[0][0]);
        }       
};


int CCellReduce::nAddCreateTwinLaw(Ccrystal& oCrystalInOut,Cstring& sRotAbout,double fRot) {
        char* pcAxes[3] = {"a*","b*","c*"};
        double a3x3fTwinMat[3][3];
        double a3x3fOrient[3][3];
        double a3fRotVec[3];
        int nAxis;
        int nTwinLawIndex;
        
        oCrystalInOut.nCalcOrientMatrix();
        oCrystalInOut.vGetOrientMatrix(&a3x3fOrient[0][0]);
        vZeroMat(3,1,&a3fRotVec[0]);
        for (nAxis = 0; nAxis < 3; nAxis ++) {
                if (sRotAbout.contains(pcAxes[nAxis])) 
                        vAddVec3DVec3D(a3fRotVec,a3x3fOrient[nAxis],a3fRotVec);
        };
        
        if (fNormVec3D(a3fRotVec)) {
                nTwinLawIndex = 1;
                vConvRotVec3DMat3D(fRot,&a3fRotVec[0],&a3x3fTwinMat[0][0]);
                oCrystalInOut.vSetTwinLaw(nTwinLawIndex,&a3x3fTwinMat[0][0],true);
        } else {
                while (oCrystalInOut.m_nTwinLaws)
                        oCrystalInOut.nClearTwinLaw(1);
        };

        return 0;
};



#ifdef SSI_PC
void CCellReduce::vBuildTableCmd(Cstring &csCommand)
{
    int iNumUsed;
    Crefln *poSoln;
    Cstring sCent;
    Cstring sLastCent;
    int iNumSolns;
    Cstring sSelectString;
    int nSelectIndex;
    int k;

    iNumUsed = 0;
    iNumSolns = m_poSolutions->nGetNumReflns();

    // Build the Tcl Command parameter list.
    for (k = 0,nSelectIndex = 1; k < iNumSolns; k++) {
        poSoln = m_poSolutions->poGetRefln(k);
        sCent = poSoln->sGetField(m_nFI_sLattice);
        
        if ((sCent == "mP") && (sLastCent == "mP")) {
            nSelectIndex--;
            sSelectString = nSelectIndex;
            sSelectString += "b";
        } else if ((sCent == "mC") && (sLastCent == "mC")) {
            nSelectIndex--;
            sSelectString = nSelectIndex;
            sSelectString += "b";
        } else 
            sSelectString = nSelectIndex;

        nSelectIndex++;
        sLastCent = sCent;

        if(("?" != poSoln->sGetField(m_nFI_sLattice)) &&
            (m_fResidualMax >= poSoln->fGetIntensity()))
        {
            iNumUsed++;
            
            csCommand += CCrclHelper::GetInstance()->sBuildIndexTableRow( poSoln, iNumUsed, sSelectString, m_nFI_sLattice, 
                m_nFI_sLattSymm, m_nFI_fA, m_nFI_fB, m_nFI_fC, m_nFI_fVolume, m_nFI_fAlpha, m_nFI_fBeta, m_nFI_fGamma );
        }
        
    }
    csCommand += "NumSolutions = ";
    csCommand += iNumUsed;
}

Cstring CCellReduce::sBuildOrientTableCmd(int nChoose)
{
    int     i;
    int     nNumSoln;
    double   fResidual;
    Crefln *poSoln;
    Cstring csCmd, csSoln, csDefaultSoln;

    // Check if solutions actually are available

    if (NULL == m_poSolutions) {
        cout << "There are no solutions present!" << endl;
        return "";
    }

    nNumSoln = m_poSolutions->nGetNumReflns();
    if (1 > nNumSoln) {
        cout << "There are no solutions present!" << endl;
        return "";
    }

    m_poSolutions->nSelect("-sLattice==D");
    m_poSolutions->nDelete((const char*)NULL, (const bool)FALSE);

    nNumSoln = m_poSolutions->nGetNumReflns();

    csCmd.Format( "IndexOrientTable NumSolns = %d", nNumSoln );
    csDefaultSoln.Format( " DefaultSoln = %d", nChoose );
        csCmd += csDefaultSoln;

    for (i = 0; i < nNumSoln; i++) {
        poSoln    = m_poSolutions->poGetRefln(i);
        fResidual = poSoln->fGetIntensity();
        if (0.0 > fResidual) {
            // Convert Fourier residual so range is 0 -> 1.0
            fResidual = 1.0f + fResidual / 1000.0f;
        }
        else {
            // Otherwise convert residual from 1 -> 0 to 0 -> 1.0
            fResidual = 1.0f - fResidual;
        }

        csSoln.Format(/*" SolnID%d %3d*/" Resid%d = %9.3lf Rot1%d = %8.3lf Rot2%d = %8.3lf Rot3%d = %8.3lf ",
            //i+1, i,
            i, (double) fResidual,
            i, (double) poSoln->fGetField(m_nFI_fRot1),
            i, (double) poSoln->fGetField(m_nFI_fRot2),
            i, (double) poSoln->fGetField(m_nFI_fRot3));

        /*printf("%3d %9.3lf %9.3lf %9.3lf %9.3lf %8.3lf %8.3lf %8.3lf\n",
            i+1,
            (double) fResidual,
            (double) poSoln->fGetField(m_nFI_fA),
            (double) poSoln->fGetField(m_nFI_fB),
            (double) poSoln->fGetField(m_nFI_fC),
            (double) poSoln->fGetField(m_nFI_fRot1),
            (double) poSoln->fGetField(m_nFI_fRot2),
            (double) poSoln->fGetField(m_nFI_fRot3));

        printf("%3s %9s %9.3lf %9.3lf %9.3lf\n",
            " ", " ",
            (double) poSoln->fGetField(m_nFI_fAlpha),
            (double) poSoln->fGetField(m_nFI_fBeta),
            (double) poSoln->fGetField(m_nFI_fGamma));
            printf("\n");*/
        csCmd += csSoln;
    }
    //printf ("%s\n", (const char *)sLine);
    //fflush(stdout);
    return csCmd;
};
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Given cell parameters and a lattice type, determine 
// the transformation matrix, necessary to transform the lattice to a standard presentation
void CCellReduce::vGetTransformToStandardPresentation(double a6fCellIn[6], double a3x3fTransformMat[3][3], Cstring& sLattice)
{
    // cP  cI cF hR hP tI tP oI oC oF oP mC mP aP //
    // 0   1  2  3  4  5  6  7  8  9  10 11 12 13
    if( sLattice.length() != 2 )
        return; // safety feature

    double      a6fCell[6];
    for(int ii=0; ii < 6; ii++)
        a6fCell[ii] = a6fCellIn[ii];


    // First determine the new cell:
    switch( sLattice.GetAt(0) ) 
    {
    case 'o':
        if( sLattice.GetAt(1) == 'C' )
        {
            if( a6fCell[0] > a6fCell[1] ) 
                std::swap(a6fCell[0], a6fCell[1]); 
        } 
        else
        {
            if( a6fCell[0] > a6fCell[2] ) 
                std::swap(a6fCell[0], a6fCell[2]);

            if( a6fCell[1] > a6fCell[2] ) 
                std::swap(a6fCell[1], a6fCell[2]);

            if (a6fCell[0] > a6fCell[1]) 
                std::swap(a6fCell[0], a6fCell[1]);
        }
        break;
    case 'm':
        if( sLattice.GetAt(1) == 'P' )
        {
            if (a6fCell[0]>a6fCell[2]) 
                std::swap(a6fCell[0],a6fCell[2]);
        } 
        break;
    case 'a':
        nStandardPres(a6fCell);
    } // end switch

    ////////////////////////////////////////////////////////////////////////
    // Now find the transformation matrix between the new and old cells

    /////////////////////////////////////////////////////////////////////////
    // Get a real orientation matrix with the old cell.
    double   a3x3fRealOrientMatIn[3][3];   
    
    Ccrystal        oCrystalIn;
    oCrystalIn.vSetCell(a6fCellIn);
    oCrystalIn.vSetOrientAngles(0.0, 0.0, 0.0);
    oCrystalIn.nCalcGetRealMatrix(&a3x3fRealOrientMatIn[0][0]);
    vTranMat3D(a3x3fRealOrientMatIn);

    /////////////////////////////////////////////////////////////////////////
    // Get a real orientation matrix with the new cell.
    double   a3x3fRealOrientMatOut[3][3];   
    
    Ccrystal        oCrystalOut;
    oCrystalOut.vSetCell(a6fCell);
    oCrystalOut.vSetOrientAngles(0.0, 0.0, 0.0);
    oCrystalOut.nCalcGetRealMatrix(&a3x3fRealOrientMatOut[0][0]);
    vTranMat3D(a3x3fRealOrientMatOut);
    ///////////////////////////////////////////////////////////////////////

    fFindOrigVecsA(a3x3fRealOrientMatIn, a3x3fRealOrientMatOut, a3x3fTransformMat, sLattice.string());
}

//This function looks at two cell solutions and decides whether they are equivalent or not
bool CCellReduce::bEquivalentCellSolutions(Crefln& oSoln_A, Crefln& oSoln_B)
{
    // RB These two factors may have to be fine-tuned as we will be testing this routine
    const double      c_fCellEquivalenceFactorLength_percent = 1.0;
    const double      c_fCellEquivalenceFactorCosAngle = 0.0001;
    
    int     ii = 0; // counter

    double      afcell_A[6];
    double      afG6_A[6];
    
    afcell_A[0] = oSoln_A.fGetField(m_nFI_fA);
    afcell_A[1] = oSoln_A.fGetField(m_nFI_fB);
    afcell_A[2] = oSoln_A.fGetField(m_nFI_fC);
    afcell_A[3] = oSoln_A.fGetField(m_nFI_fAlpha);
    afcell_A[4] = oSoln_A.fGetField(m_nFI_fBeta);
    afcell_A[5] = oSoln_A.fGetField(m_nFI_fGamma);
    
    vCalcGetG6(afcell_A, afG6_A);
    
    double     fLength_A = 0.0; 
    for(ii=0; ii < 6; ii++)
    {
        fLength_A += afG6_A[ii] * afG6_A[ii];
    }
    
    fLength_A = fLength_A > 0.0 ? sqrt(fLength_A) : 0.0;


    double      afcell_B[6];
    double      afG6_B[6];

    afcell_B[0] = oSoln_B.fGetField(m_nFI_fA);
    afcell_B[1] = oSoln_B.fGetField(m_nFI_fB);
    afcell_B[2] = oSoln_B.fGetField(m_nFI_fC);
    afcell_B[3] = oSoln_B.fGetField(m_nFI_fAlpha);
    afcell_B[4] = oSoln_B.fGetField(m_nFI_fBeta);
    afcell_B[5] = oSoln_B.fGetField(m_nFI_fGamma);
    
    vCalcGetG6(afcell_B, afG6_B);

    double     fLength_B = 0.0; 
    for(ii=0; ii < 6; ii++)
    {
        fLength_B += afG6_B[ii] * afG6_B[ii];
    }
    
    fLength_B = fLength_B > 0.0 ? sqrt(fLength_B) : 0.0;

    double      fMeanLength = (fLength_A + fLength_B)/2.0;

    if( 0.0 == fMeanLength )
        return true;

    if( fabs(fLength_A - fMeanLength) / fMeanLength * 100.0 > c_fCellEquivalenceFactorLength_percent )
        return false;
    
    double      fScalarProduct = 0.0;
    for(ii=0; ii < 6; ii++)
    {
        fScalarProduct += afG6_A[ii] * afG6_B[ii];
    }
    
    double      fCosAngle = fScalarProduct / fLength_A / fLength_B;
    
    if( fCosAngle < 0.0 || fabs(1.0 - fCosAngle) > c_fCellEquivalenceFactorCosAngle )
        return false;

    return true;
}


/*  Basic structure of CCellReduce Calls:

nAddCrystal()
nAddCrystal()
...
nCalcGetBestLattices();
nListResults()
nListResultsB()
nGetOrientInfo()


*/





