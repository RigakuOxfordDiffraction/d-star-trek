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
// CCellReduce.h             Initial author: Thaddeus J Niemeyer
//  This file contains the declarations of class CCellReduce.
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

#ifndef DT_CCELLREDUCE_H
#define DT_CCELLREDUCE_H

#include "Dtrek.h"
#include "Cimage_header.h"
#include "Cdetector.h"
#include "Csource.h"
#include "Ccrystal.h"
#include "Crotation.h"
#include "Creflnlist.h"
#include "C3Ddata.h"
#include "Cimage.h"
#include "dtrekvec.h"

const int DT_NUMSOL = 10;

typedef struct DTREK_EXPORT _tagLattChar
{
  int     nLattChar;
  int     nSortOrder;
  int     nSpacegroupMin;
  int     nType;
  int     nTypeDependent;
  double   fMatFactor;
  double   a3x3fTMat[3][3];
  double   a6x6fLSMat[6][6];
  Cstring sLattSymm;
  Cstring sLatticeType;
} tagLattChar;


class DTREK_EXPORT CCellReduce {
    int    m_nFI_fA;
    int    m_nFI_fB;
    int    m_nFI_fC;
    int    m_nFI_fAlpha;
    int    m_nFI_fBeta;
    int    m_nFI_fGamma;
    int    m_nFI_fVolume;
    int    m_nFI_fRot1;
    int    m_nFI_fRot2;
    int    m_nFI_fRot3;
    int    m_nFI_sLattice;
    int    m_nFI_sLattSymm;
    int    m_nFI_nLattNum;
    int    m_nFI_nSortOrder;


    static Cstring ms_sfA;
    static Cstring ms_sfB;
    static Cstring ms_sfC;
    static Cstring ms_sfAlpha;
    static Cstring ms_sfBeta;
    static Cstring ms_sfGamma;
    static Cstring ms_sfVolume;
    static Cstring ms_sfRot1;
    static Cstring ms_sfRot2;
    static Cstring ms_sfRot3;
    static Cstring ms_sfDiffLength;
    static Cstring ms_ssLattice;
    static Cstring ms_ssLattSymm;
    static Cstring ms_snLattNum;
    static Cstring ms_snSortOrder;
	static char*	ms_cpRotEquivsStrings [14+1];
	static int		ms_pnLatticeChars[14+1];

    Creflnlist* m_poSolutions;      // A place for solutions;
    Ccrystal    m_oCrystal;         // A place for the solution
    Ccrystal    m_oKnown;           // The known crystal.

    double       m_a26x6x6fBMat[26][6][6];
    double       m_a3x3fUBMat[3][3];
    double       m_a3x3fInvUBMat[3][3];


    void vInitSolutions(void);
    void vInitLattChar(void);

    static int  ms_anOrder[];				 // Ordering of 44 lattice chars. w.r.t. Int. Tables.    
	static int  ms_anOrderSlot[];			 // We reseve two slots for some of the lattices (mP, mC).  This array gives the slot #.
    double      m_a44x6fLattCells[45][6];    // These hold the best lattice cells for all 44 lattice characters 
    double      m_a44fLattResid[45];         // These hold the lattice residuals.
    int         m_a15x3nBIndices[15][3];


public:
    int nGetCurrentLatticeNumber()const{return m_nLattNum;}
    
    tagLattChar m_a44tLattChar[45];

    int nCalcGetCellFromG6(const double *pfG6, double *pfCell);
    int nCalcGetCellFromG6(const double *pfG6, double a3x3fOrientMat[3][3]);
    void vCalcGetG6(const double *pfCell, double *pfG6);
    void vCalcGetG6(const double a3x3fOrientMat[3][3], double *pfG6);
    int  nReduceCell(double *pfCell);
    
    int  nStandardPres(double *pfCell);
    void vGetTransformToStandardPresentation(double a6fCell[6], double a3x3fOrientMat[3][3], Cstring& sLattice);

    int nCentToPrim(char cCentering, double *pfIO,bool bIsG6);
    int nReducedPToCent(const Cstring& sFlag, double *pfG6, int *pnLattNum=NULL);
    int nReducedPToCent(const int nLattNumIn, double *pfG6);

#ifdef SSI_PC
    void vBuildTableCmd(Cstring &csCommand);
    Cstring sBuildOrientTableCmd(int nChoose);
#endif

    // Loads the m_poSolutions with equivalent rotations.
    int nLoadRotationChoices(const int nLattNum = -1, Ccrystal *poCrystalIn = NULL);

    int nEquivRotations(double fMatIn[3][3], double fMatOut[24][3][3], Cstring sLattId); // Returns the number of equivalent cells.  Each loaded into the rotation matrix.
    int nEquivRotationsRhombohedralHexagonal(double fMatIn[3][3], double fMatOut[24][3][3]);

    double fFindOrigVecs(double fInCentered[3][3], double fInPreReduced[3][3],double fOutCoords[3][3]);// Returns the best fit residual (% number)
                      
    double fFindOrigVecsA(double a3x3fInCentered[3][3], double a3x3fInPreReduced[3][3], double a3x3fOutCoords[3][3], char*  pcLatticeType=NULL); // RB 1/6/2005 A wrapper around fFindOrigVecs

    bool        m_bKnownCell;
    int         m_nVerbose;
    int         m_nLattNum;
    double      m_fResidualMax;
    double      m_fWavelength;
    
    int nLoadSolution(const int nSolnNum, Ccrystal *poCrystalIn = NULL);
    int nAddSolutionL(const int nLattice);
    Cstring sGetLattice(const int nSolnNum,bool bLongName = false);
    double fGetSolutionResidual(const int nSolnNum);
    int nGetNumSolutions() { return m_poSolutions->nGetNumReflns(); };
    int nParseSelection(Cstring& sInput,int& nEntry);
    Creflnlist& oGetSolutions() { return *m_poSolutions; };

    int nListResults(int *pnBest=NULL);
    int nListResultsB();
    int nListResultsL();
    int nSelectMatch(const int nSpacegroupNum);
    int nAddCrystal(Ccrystal& oCrystalIn);
    int nAddCrystal(double* a3x3fOrientMatIn);
    int nGetOrientInfo(Ccrystal& oCrystalIn,Ccrystal& oCrystalOut);
	int nAddCreateTwinLaw(Ccrystal& oCrystalInOut,Cstring& sRotAbout,double fRot);
    double fGetOrientInfo(Ccrystal& oCrystalIn,Ccrystal& oCrystalOut,Ccrystal* poCrystalPrev = NULL); // If you specify the lattice, it can use symmetry to get the best answer
    int nCalcGetBestLattices();
    void vConvertToPrimitive(Ccrystal& oCrystal,char* pcLattice = NULL);


    Ccrystal& oGetKnownCell() { return m_oKnown; };
    
    CCellReduce();
    ~CCellReduce();
    private:
        bool bEquivalentCellSolutions(Crefln& oSoln_A, Crefln& oSoln_B);
};


#endif
