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
// Cindex.cc            Initial author: Thaddeus J Niemeyer     11-Nov-1999
//  This file contains the member functions of class Cindex which implements
//    the autoindexing encapsulation of d*TREK.
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

#include "Cstring.h"
#include "Dtrek.h"
#include "Cimage.h"
#include "dtrekvec.h"
#include "dtsvd.h"
#include "Crefln.h"
#include "Creflnlist.h"
#include "dtcell.h"
#include "CCellReduce.h"

#include <string.h>
#include <stdio.h>
#include <ctype.h>

#ifdef SSI_PC
#include "CrclHelper.h"
#endif

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
#endif


char* g_pcTexsanSysAbsString =
" Odd/Even parity class                     \n"
" --------------------------------          \n"
" Class  Total  Observed   <I/sig>          \n"
" --------------------------------          \n"
"   eee  %  5d    %  5d      %5.1f%-4c      \n"
"   eeo  %  5d    %  5d      %5.1f%-4c      \n"
"   eoe  %  5d    %  5d      %5.1f%-4c      \n"
"   eoo  %  5d    %  5d      %5.1f%-4c      \n"
"   oee  %  5d    %  5d      %5.1f%-4c      \n"
"   oeo  %  5d    %  5d      %5.1f%-4c      \n"
"   ooe  %  5d    %  5d      %5.1f%-4c      \n"
"   ooo  %  5d    %  5d      %5.1f%-4c      \n"
" ---------------------------------         \n"
"                                           \n"
" --------------------------------   --------------------------------       \n"
"   hk0  Total  Observed   <I/sig>     h0l  Total  Observed   <I/sig>       \n"
" --------------------------------   --------------------------------       \n"
"    ee  %  5d    %  5d      %5.1f%-4c e e  %  5d    %  5d      %5.1f%-4c   \n"
"    eo  %  5d    %  5d      %5.1f%-4c e o  %  5d    %  5d      %5.1f%-4c   \n"
"    oe  %  5d    %  5d      %5.1f%-4c o e  %  5d    %  5d      %5.1f%-4c   \n"
"    oo  %  5d    %  5d      %5.1f%-4c o o  %  5d    %  5d      %5.1f%-4c   \n"
" --------------------------------   --------------------------------       \n"
"                                                                           \n"
" --------------------------------   --------------------------------       \n"
"   0kl  Total  Observed   <I/sig>     hhl  Total  Observed   <I/sig>       \n"
" --------------------------------   --------------------------------       \n"
"    ee  %  5d    %  5d      %5.1f%-4c e e  %  5d    %  5d      %5.1f%-4c   \n"
"    eo  %  5d    %  5d      %5.1f%-4c e o  %  5d    %  5d      %5.1f%-4c   \n"
"    oe  %  5d    %  5d      %5.1f%-4c o e  %  5d    %  5d      %5.1f%-4c   \n"
"    oo  %  5d    %  5d      %5.1f%-4c o o  %  5d    %  5d      %5.1f%-4c   \n"
" --------------------------------   --------------------------------       \n"
"                                                                           \n"
" --------------------------------                                          \n"
"   h-hl Total  Observed   <I/sig>                                          \n"
" --------------------------------                                          \n"
"    ee  %  5d    %  5d      %5.1f%-4c                                      \n"
"    eo  %  5d    %  5d      %5.1f%-4c                                      \n"
"    oe  %  5d    %  5d      %5.1f%-4c                                      \n"
"    oo  %  5d    %  5d      %5.1f%-4c                                      \n"
" --------------------------------                                          \n"
"                                                                           \n"
" --------------------------------   --------------------------------       \n"
"   hhh  Total  Observed   <I/sig>     hh0  Total  Observed   <I/sig>       \n"
" --------------------------------   --------------------------------       \n"
"     e  %  5d    %  5d      %5.1f%-4c   e  %  5d    %  5d      %5.1f%-4c   \n"
"     o  %  5d    %  5d      %5.1f%-4c   o  %  5d    %  5d      %5.1f%-4c   \n"
" --------------------------------   --------------------------------       \n"
"                                                                           \n"
" --------------------------------   --------------------------------       \n"
"   h00  Total  Observed   <I/sig>     0k0  Total  Observed   <I/sig>       \n"
" --------------------------------   --------------------------------       \n"
"     e  %  5d    %  5d      %5.1f%-4c   e  %  5d    %  5d      %5.1f%-4c   \n"
"     o  %  5d    %  5d      %5.1f%-4c   o  %  5d    %  5d      %5.1f%-4c   \n"
"   o/e     --       --      %5.1f     o/e     --       --      %5.1f       \n"
" --------------------------------   --------------------------------       \n"
"                                                                           \n"
" --------------------------------                                          \n"
"   00l  Total  Observed   <I/sig>                                          \n"
" --------------------------------                                          \n"
"     e  %  5d    %  5d      %5.1f%-4c                                      \n"
"     o  %  5d    %  5d      %5.1f%-4c                                      \n"
"   o/e     --       --      %5.1f                                          \n"
" --------------------------------                                          \n"
"                                                                           \n"
" ---------------------------------------------------------------           \n"
"  Zone           a+b=4n                  a+b not equal 4n                  \n"
"        Total  Observed    <I/sig>    Total  Observed    <I/sig>           \n"
" ---------------------------------------------------------------           \n"
"   0kl  %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
"   h0l  %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
"   hk0  %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
"                                                                           \n"
"                                                                           \n"
" ---------------------------------------------------------------           \n"
"  Zone            a=4n                    a not equal 4n                   \n"
"        Total  Observed    <I/sig>    Total  Observed    <I/sig>           \n"
" ---------------------------------------------------------------           \n"
"   0k0  %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
"   00l  %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
"   h00  %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
" ---------------------------------------------------------------           \n"
"                                                                           \n"
" ---------------------------------------------------------------           \n"
"  Zone           2h+l=4n                2h+l not equal 4n                  \n"
"        Total  Observed    <I/sig>    Total  Observed    <I/sig>           \n"
" ---------------------------------------------------------------           \n"
"   hhl  %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
" ---------------------------------------------------------------           \n"
"                                                                           \n"
" --------------------------------------------------------------------------\n"
"           h+l=3n;l odd             h+l=3n             h+l not equal 3n    \n"
"        Total Obsvd <I/sig>    Total Obsvd <I/sig>    Total Obsvd <I/sig>  \n"
" --------------------------------------------------------------------------\n"  
" h-h0l  %  5d %  5d  %5.1f%-4c %  5d %  5d  %5.1f%-4c %  5d %  5d  %5.1f%-4c\n"
" -------------------------------------------------------------------------- \n"
"          -h+l=3n;l even                 -h+l=3n           -h+l not equal 3n\n"
" -------------------------------------------------------------------------- \n"     
" h-h0l  %  5d %  5d  %5.1f%-4c %  5d %  5d  %5.1f%-4c %  5d %  5d  %5.1f%-4c\n" 
" -------------------------------------------------------------------------- \n"
"                                                                           \n"
"                                                                           \n"
" ---------------------------------------------------------------           \n"
"                  l=3n                    l not equal 3n                   \n"
"        Total  Observed    <I/sig>    Total  Observed    <I/sig>           \n"
" ---------------------------------------------------------------           \n"
" 000l   %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
" ---------------------------------------------------------------           \n"
"                  l=6n                    l not equal 6n                   \n"
" ---------------------------------------------------------------           \n"
" 000l   %  5d    %  5d      %5.1f%-4c %  5d   %  5d       %5.1f%-4c        \n"
" ---------------------------------------------------------------           \n"
;


#define ERROR {  printf("Bad Logic Line %d of file '%s'\n", __LINE__ , __FILE__ ); GETCH; exit(0); };
//template <class T> void swap(T& x,T& y) { T z; z=x; x=y; y=z; };
#define STR(x) #x


int nLaueTableCompare(const void* pvT1,const void* pvT2)
{
	tagLaueTable* poT1 = (tagLaueTable*) pvT1;
	tagLaueTable* poT2 = (tagLaueTable*) pvT2;
	int nLatt1,nLatt2;

	
	for (nLatt1 = 0;nLatt1 < 14; nLatt1++) {
		if (g_pcValidLattices[nLatt1] == (Cstring) poT1->acLattice)
			break;
	};
	for (nLatt2 = 0;nLatt2 < 14; nLatt2++) {
		if (g_pcValidLattices[nLatt2] == (Cstring) poT2->acLattice)
			break;
	};
	

    
	if (nLatt1 < nLatt2)
		return -1;
	else if (nLatt1 > nLatt2)
		return 1;
	else if (poT1->nLaueGroup < poT2->nLaueGroup)
		return -1;
	else if (poT1->nLaueGroup > poT2->nLaueGroup)
		return 1;
	else if (poT1->cUniqueAxis < poT2->cUniqueAxis)
		return -1;
	else if (poT1->cUniqueAxis > poT2->cUniqueAxis)
		return 1;
	else
		return 0;
};


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

/* This is used to map the Laue class to a space group.
   Values of 400 indicate the user defined space groups provided below
*/


int Claue::ms_nlaue_group[]={
	2,
	10,
	47,
	83,
	123,
	147,
	164,
	162,
	175,
	191,
	200,
	221
};

Claue::Claue() {
    int nx;

    // Set to a cubic cell, since we don't have any other cell available.
    fCellParams[0]=10;
    fCellParams[1]=10;
    fCellParams[2]=10;
    fCellParams[3]=90;
    fCellParams[4]=90;
    fCellParams[5]=90;
    // Set the sigmas.
    for (nx=0; nx< 6; nx++)
        fCellSigmas[nx]=0.0;

    fOrientAngles[0]=0.0;
    fOrientAngles[1]=0.0;
    fOrientAngles[2]=0.0;
    a2cLattice[0]=0;
    a2cLattice[1]=0;
    a2cLattice[2]=0;
    fMosaicity = 0.0;
    nDup();
}

Claue::~Claue()
{
}


int Claue::nDup()
{
    fRMergeGreat=0.08;
    fRMergeAwful=0.20;
    fMaxRMergeRatio = 1.5;
    fAbsentSigma=4;
    fLaueSigma=3;
    bAnom = FALSE;
    bScaleForLaue = FALSE;
	nScaleForLaueMinRefs = 50;
    bCollectRejects = FALSE;
    nMaxCollectRejects = 0;
    eMaxLaueCheck = LAUE_UNKNOWN;
    nMinConclusiveAbsences = 3;
    nMinConclusiveExistences = 5;
    nSetChirality(FALSE);
    nTwinID = 1;
    fResidMax = 5.0;

    m_bIsChiral = FALSE;
    m_bIsCentric = FALSE;
    m_eLaue=LAUE_UNKNOWN; // The default is to NOT update these parameters.
    m_nSpaceGroup=-1;
    m_strInconclusiveConditions = "";
    return 0;
};


int Claue::nSetChirality(bool bUseCrystalVolumeForDefault) {   
    if ((sGetEnv("DTREK_DTCELL_CHIRAL").length()!=0)) {
        m_bIsChiral = TRUE;
    } else {
        if (bUseCrystalVolumeForDefault) {
            Ccrystal oCrystal;
            oCrystal.vSetCell(fCellParams);
            if (oCrystal.fCalcVolume()>=25000.0) {
                printf("Assuming chiral cell since cell volume is greater than 25000.0 cubic angstroms.\n");
                printf("Please use -chiral or -chiralunknown to override this default.\n");
                m_bIsChiral = TRUE;
            } else {
                printf("No assumptions on chirality have been made since cell volume is less than 25000.0 cubic angstroms.\n");
                printf("Please use -chiral or -chiralunknown to override this default.\n");
                m_bIsChiral = FALSE;
            };
        } else 
            m_bIsChiral = FALSE;
    };
    return 0;
};

eLaueType Claue::eLauefromSpacegroup(int nSpacegroupIn) {
    eLaueType eType;
    
    if ((nSpacegroupIn>=1) && (nSpacegroupIn<=2))
        eType=LAUE_1_BAR;
    else if ((nSpacegroupIn>=3) && (nSpacegroupIn<=15))
        eType=LAUE_2_OVER_M;
    else if ((nSpacegroupIn>=16) && (nSpacegroupIn<=74))
        eType=LAUE_MMM;
    else if ((nSpacegroupIn>=75) && (nSpacegroupIn<=88))
        eType=LAUE_4_OVER_M;
    else if ((nSpacegroupIn>=89) && (nSpacegroupIn<=142))
        eType=LAUE_4_OVER_MMM;
    else if ((nSpacegroupIn>=143) && (nSpacegroupIn<=148))
        eType=LAUE_3_BAR;
    else if ((nSpacegroupIn>=149) && (nSpacegroupIn<=167)) {
        eType=LAUE_3_BAR_M_1;
        switch (nSpacegroupIn) {
        case 149: 
        case 151:
        case 153:
        case 157:
        case 159:
        case 162:
        case 163: eType= LAUE_3_BAR_1_M;
        };
    } else if ((nSpacegroupIn>=168) && (nSpacegroupIn<=176))
        eType=LAUE_6_OVER_M;
    else if ((nSpacegroupIn>=177) && (nSpacegroupIn<=194))
        eType=LAUE_6_OVER_MMM;
    else if ((nSpacegroupIn>=195) && (nSpacegroupIn<=206))
        eType=LAUE_M_3_BAR;
    else if ((nSpacegroupIn>=207) && (nSpacegroupIn<=230))
        eType=LAUE_M_3_BAR_M;

    return eType;
};

int Claue::nSpacegroupFromLaue(eLaueType eLaueIn){
    if (!a2cLattice[0])
        a2cLattice[1] = 0;

    switch (eLaueIn) {
    case LAUE_UNKNOWN:
        return -1;
        break;
    case LAUE_1_BAR:
        return 1;
        break;
    case LAUE_2_OVER_M:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 3;
        else
            return 5;
        break;
    case LAUE_MMM:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 16;
        else if (a2cLattice[1] == 'C')
            return 21;
        else if (a2cLattice[1] == 'F')
            return 22;
        else if (a2cLattice[1] == 'I')
            return 23;
        break;
    case LAUE_4_OVER_M:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 75;
        else if (a2cLattice[1] == 'I')
            return 79;
        break;
    case LAUE_4_OVER_MMM:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 89;
        else if (a2cLattice[1] == 'I')
            return 97;
        break;
    case LAUE_3_BAR: 

        /*  Notes on Rhombohedral and hexagonal spacegroups:
            1)  Rhombohedral and trigonal spacegroups are distinct, but they are both
                listed in the "trigonal" Laue table because they share the same Laue classes.
            2)  The Rhombohedral ones can be indexed on both Rhombodedral (primative) or hexagonal axes (triplely redundant)
            3)  In the spacegroup tables for example, H3, R3 and R3(h) are all names for "R3".
            4)  The entire set of Rhombohedrals are: { 146 148 155 160 161 166 167 }
        */

        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 143;
        else if (a2cLattice[1] == 'R')
            return 146;
        break;

    case LAUE_3_BAR_M_1:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 150;
        else if (a2cLattice[1] == 'R')
            return 155;
        break;
    case LAUE_3_BAR_1_M:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 149;
        break;
    case LAUE_6_OVER_M:
        return  168;
        break;
    case LAUE_6_OVER_MMM:
        return 177;
        break;
    case LAUE_M_3_BAR:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 195;
        else if (a2cLattice[1] == 'F')
            return 196;
        else if (a2cLattice[1] == 'I')
            return 197;
        break;
    case LAUE_M_3_BAR_M:
        if ((a2cLattice[1] == 'P') || (a2cLattice[1] == 0))
            return 207;
        else if (a2cLattice[1] == 'F')
            return 209;
        else if (a2cLattice[1] == 'I')
            return 211;
        break;
    }
    return -1;
};

int Claue::nSetCrystal(Cimage_header& oHeader) {
    int nx;    
    Cstring sTemp;

    Ccrystal oCrystal(oHeader,nTwinID);
    if (oCrystal.bIsAvailable()) {
        
        fMosaicity = oCrystal.fGetMosaicity();
        oCrystal.vGetCell(fCellParams);
        oCrystal.vGetOrientAngles(fOrientAngles);        
        // Also want to load the sigmas if they are available.
        oHeader.nGetValue(sTemp="CRYSTAL_UNIT_CELL_SIGMA",6,fCellSigmas);
        
        // If a 'spacegroup' number was already in the data, we want to set the Bravis lattice.
        // to be compatible with that spacegroup.
        if (oCrystal.m_poSpacegroup) {
            nx = oCrystal.m_poSpacegroup->nGet();
            if ((nx<=0) || (nx>230)) {
                printf("Space group number %d in header is invalid.\n",nx);
                return 1;
            };
            a2cLattice[0]=oCrystal.m_poSpacegroup->cGetClass();
            a2cLattice[1]=g_cpSpaceGroupNames[nx-1][0];
            printf("Lattice %s assumed for space group %d.\n",a2cLattice,nx);
            m_nSpaceGroup = nx;
        };
        return 0;
    } else
        return 1;    
};



int Claue::nPrint(Creflnlist& oRef) {
	int nx,ny,nz;
	int nlastpackedHKL=0L;
	int npackedHKL;

	for (nx=0,nz=0;nx<oRef.nGetNumReflns();nx++) {
		ny=oRef.pnGetSortIndex()[nx];
		npackedHKL=oRef[ny].nGetField(oRef.m_nFI_nPackedHKL);
		if ((nx==0) || (npackedHKL!=nlastpackedHKL)) {
			printf("--------------------------[%d]\n",nz++);		
		};
		printf("%d %d %d\n",oRef[ny].nGetH(),oRef[ny].nGetK(),oRef[ny].nGetL());
		nlastpackedHKL=npackedHKL;
	};
	return 0;
};


int Claue::nWriteHeader(Cstring* psName,Cimage_header& oHeader) {
   Cstring sTemp;
   Cstring sLaue;
   Cstring sPrefix;
   Ccrystal oCrystalOut;
   int nx;
   
   if (m_eLaue!=LAUE_UNKNOWN) {
        sLaue=g_cpLaueName[m_eLaue];
        oHeader.nReplaceValue(sTemp="CRYSTAL_LAUE",sLaue);
   };

   oCrystalOut.vSetOrientAngles(&fOrientAngles[0]);
   oCrystalOut.vSetCell(fCellParams);
   oCrystalOut.vSetTwin(nTwinID);
   oCrystalOut.nLoadKeywords(nTwinID);
   if (fMosaicity)
    oCrystalOut.vSetMosaicity(fMosaicity);
   if (m_nSpaceGroup>=1)
       oCrystalOut.m_poSpacegroup->vSet(m_nSpaceGroup);
   else if ((nx = nSpacegroupFromLaue(m_eLaue))!=LAUE_UNKNOWN)
       oCrystalOut.m_poSpacegroup->vSet(nx);
   else {
       printf("WARNING:  Do not know what the spacegroup should be!\n");
       printf("          Assuming triclinic.\n");
       oCrystalOut.m_poSpacegroup->vSet(1);
   };
   oCrystalOut.nUpdateHeader(&oHeader);

   oHeader.nReplaceValue(oCrystalOut.ms_sCrystalXUnitCell + "_SIGMA", 6, fCellSigmas, 4);
   
   if( m_strInconclusiveConditions != "" )
        oHeader.nReplaceValueDictionary("DTCELL", "INCONCLUSIVE_HKL_CONDITIONS", m_strInconclusiveConditions);

   if (psName) {
       oHeader.nWrite(*psName);
       printf("Header \"%s\" written.\n",psName->string());
   };
   return 0;
  
}

// We want to get the reindexing transformation needed when moving from one orientation (UB) matrix to another (UB) matrix.

int Claue::nCalcReindex(double a3x3fOrientMatIn[3][3], double a3x3fOrientMatOut[3][3], double a3x3fReindexMat[3][3]) {
    
    // Figure out the transformation that is needed to take the old UB matrix to the new UB matrix.
    
    double a3x3fTempMat[3][3];
    double a3x3fTrans[3][3];
    double a3x3fInvTrans[3][3];
    bool  bCanTrans=TRUE;
    int nx,ny;

    // [UB]old * h = [UB]old * [T] * [T]inv * h  = [UB]new * [T]inv * h
    // With [UB]new = [UB]old * [T]
    // Thus, [T] = [UB]old-inv * [UB]new and 
    // the reindexing matrix is [T]inv
    
    fInvMat3D(&a3x3fOrientMatIn[0][0],&a3x3fTempMat[0][0]);
    vMulMat3DMat3D(a3x3fTempMat,a3x3fOrientMatOut,a3x3fTrans);
    fInvMat3D(&a3x3fTrans[0][0],&a3x3fInvTrans[0][0]);
    

    for (nx=0;nx<3;nx++)
        for (ny=0;ny<3;ny++) {

            // the code below has been commented out because non-integral transformation matrices ARE valid.
            /*
            // We must have integral values in this matrix.  Otherwise it is invalid.
            f1=a3x3fInvTrans[nx][ny];
            f0=floor(a3x3fInvTrans[nx][ny]);
            if (a3x3fInvTrans[nx][ny]-f0>0.5)
                a3x3fInvTrans[nx][ny]=f0+1;
            else
            
                a3x3fInvTrans[nx][ny]=f0;

            if (fabs(f1-a3x3fInvTrans[nx][ny])>0.01) 
                bCanTrans=FALSE;

            */
            if (fabs(a3x3fInvTrans[nx][ny])<0.0001)
                a3x3fInvTrans[nx][ny]=0.0;
            
            a3x3fReindexMat[nx][ny]=a3x3fInvTrans[nx][ny];
        };    

    return 1;
        
};

int Claue::nReindexError(double a3x3fReindexMat[3][3],Creflnlist& oReflnlist,bool bPrompt)
{
    int nx;
    double f0,f1;
    
    int nRef;
    double fSumOffend;
    int nNumOffend;
    double fIoverSigOffend;

    double a3fHKLOld[3];
    double a3fHKLNew[3];

    int* pnDeleteArray;
    int nRetValue;
    bool bWarnAlways = TRUE;

    Cstring sTemp;

    nRetValue = 0;

    fSumOffend=0.0;
    nNumOffend=0;
    fIoverSigOffend=0;

    pnDeleteArray =new int[oReflnlist.nGetNumReflns()];

    // We need to go through all of the reflections, and count the number of offending reflections.

    for (nRef=0;nRef<oReflnlist.nGetNumReflns();nRef++) {
        a3fHKLOld[0]=oReflnlist[nRef].nGetH();
        a3fHKLOld[1]=oReflnlist[nRef].nGetK();
        a3fHKLOld[2]=oReflnlist[nRef].nGetL();
        vMulMat3DVec3D(a3x3fReindexMat,&a3fHKLOld[0],&a3fHKLNew[0]);
        // See if the HKL has non HKL entries.
        for (f0=0,nx=0;nx<3;nx++) {
            f1 = floor(a3fHKLNew[nx]);
            f0 += min(fabs(f1-a3fHKLNew[nx]),fabs(f1+1-a3fHKLNew[nx]));
        };
        if (f0>0.001) {
            f0 = oReflnlist[nRef].fGetIntensity();
            f1 = oReflnlist[nRef].fGetSigmaI();
            if ((f1==0.0) || (f0<0.0) || (f0/f1>fAbsentSigma)) {
                nNumOffend++;
                fSumOffend+=f0;
                if ((f0!=0.0) && (f1!=0.0))
                    fIoverSigOffend+=f0/f1;
                pnDeleteArray[nx] = 1;
            } else 
                pnDeleteArray[nx] = 0;

        } else
            pnDeleteArray[nx] = 0;
    };
    
    if( fDetMat3D(&a3x3fReindexMat[0][0]) <= 0.0 )
    {

#ifdef SSI_PC
        CCrclHelper::GetInstance()->vBuildCmdCellLeftHandedTransform(true);
#endif
        printf("\nWARNING:  Left Handed transformation detected.\nAre you sure you want to reindex? Y/[N]: ");
        
        getline(cin,sTemp);

#ifdef SSI_PC
        printf("%s\n", sTemp.string());
#endif
        sTemp.upcase();						 
    
        if (sTemp.length() <=0)
        {
            goto exit_out;
        }
        else if (sTemp.GetAt(0) != 'Y')
        {
            goto exit_out;
        }
    }
    else
    {
#ifdef SSI_PC
        CCrclHelper::GetInstance()->vBuildCmdCellLeftHandedTransform(false);
        CCrclHelper::GetInstance()->SendTclCmd();
#endif
    }

    if ( (bWarnAlways) || (nNumOffend!=0) )
      {
       if (nNumOffend>0) 
	 {
	  fSumOffend      /= nNumOffend;
	  fIoverSigOffend /= nNumOffend;
        }
       else 
	 {
	   fSumOffend     = 0.0;
	   fIoverSigOffend = 0.0;
	 }

       if (nNumOffend > 0)
	 {
	   printf("\nWARNING: %d out of %d reflections had\n"
		  "           I/sig(I)>= %4f AND failed to reindex:\n\n"
		  " Number of offending reflections:           %d\n"
		  " Average value of offending reflections:    %6.4f\n"
		  " Average I/sig(I) value of offenders:       %6.4f\n\n",
		  nNumOffend, oReflnlist.nGetNumReflns(), fAbsentSigma, 
		  nNumOffend,fSumOffend,fIoverSigOffend);
	 }
       else
	 {
	   printf("\nINFO: No reflections out of %d had\n"
		  "           I/sig(I)>= %4f and failed to reindex.\n\n",
		  oReflnlist.nGetNumReflns(), fAbsentSigma);
	 }
       
       if( bPrompt )
       {
#ifdef SSI_PC
           CCrclHelper::GetInstance()->vBuildCmdCellReindex(oReflnlist.nGetNumReflns(),
                                                            nNumOffend,
                                                            (float)fSumOffend,
                                                            (float)fIoverSigOffend);
#endif
           printf(" Do you wish to reindex? [Y]/N: ");
           getline(cin,sTemp);
#ifdef SSI_PC
		   printf("%s\n", sTemp.string());
#endif
           sTemp.downcase();
           if (!((sTemp.length() == 0) || (sTemp.GetAt(0) == 'y'))) {
               goto exit_out;
           };
       } else if (nNumOffend!=0)
           goto exit_out;
       if (nNumOffend > 0)
       {
           printf("Deleting all offending reflections...\n");
           oReflnlist.nDelete(1, pnDeleteArray);
       }
    };
    nRetValue = 1;
    
exit_out:
    delete[] pnDeleteArray;
    if (nRetValue == 0) {
        printf("No reindexing done.\n");
    };
    return nRetValue;
};

bool Claue::bCellChanges(double fPercentError,double* a6fCell,double a3x3fTransform[3][3],char* pcOutputString) {
    Ccrystal oCrystal;
    double a3x3fBMatrix[3][3];
    double a3x3fInvTransform[3][3];
    double a3x3fTempMat[3][3];

    double a6fNewCell[6];
    int nNumOutOfBounds;

    // We multiply the B matrix from the right by the inverse of a3x3fTransform.

    oCrystal.vSetCell(a6fCell);
    oCrystal.nCalcBMatrix();
    oCrystal.vGetBMatrix(&a3x3fBMatrix[0][0]);

    fInvMat3D(&a3x3fTransform[0][0],&a3x3fInvTransform[0][0]);
    vMulMat3DMat3D(a3x3fBMatrix,a3x3fInvTransform,a3x3fTempMat);
    oCrystal.nSetOrientMatrix(&a3x3fTempMat[0][0]);
    oCrystal.vGetCell(a6fNewCell);

    nNumOutOfBounds = 0;
    sprintf(pcOutputString,"WARNING:  Reindexing causes cell parameter(s) ( %s %s %s %s %s %s) to change by more than %4.2f percent.",
        (fabs(a6fNewCell[0]-a6fCell[0])/max(a6fNewCell[0],a6fCell[0])>fPercentError)?((nNumOutOfBounds++),"a"):"",
        (fabs(a6fNewCell[1]-a6fCell[1])/max(a6fNewCell[1],a6fCell[1])>fPercentError)?((nNumOutOfBounds++),"b"):"",
        (fabs(a6fNewCell[2]-a6fCell[2])/max(a6fNewCell[2],a6fCell[2])>fPercentError)?((nNumOutOfBounds++),"c"):"",
        (fabs(a6fNewCell[3]-a6fCell[3])/max(a6fNewCell[3],a6fCell[3])>fPercentError)?((nNumOutOfBounds++),"alpha"):"",
        (fabs(a6fNewCell[4]-a6fCell[4])/max(a6fNewCell[4],a6fCell[4])>fPercentError)?((nNumOutOfBounds++),"beta"):"",
        (fabs(a6fNewCell[5]-a6fCell[5])/max(a6fNewCell[5],a6fCell[5])>fPercentError)?((nNumOutOfBounds++),"gamma"):"",
        fPercentError*100.0
        );
    if (nNumOutOfBounds>0)
        return TRUE;
    else {
        pcOutputString[0]=0;
        return FALSE;
    };
};


int Claue::nCalcCell(double fReindexMat[3][3]) {
    double fTransMat[3][3];
    double fOrientMat[3][3];
    double fNewOrientMat[3][3];
    double fOrigCellParams[6];         // Saved for the sigma transformations.
    Ccrystal oCrystal;
    int nx,ny;
    double f0;
	
	for (nx=0;nx<6;nx++)
		fOrigCellParams[nx]=fCellParams[nx];


    // Get the transform matrix.  We have U B inv(T) T h
    fInvMat3D(&fReindexMat[0][0],&fTransMat[0][0]);

    // set a cell.
    oCrystal.vSetOrientAngles(&fOrientAngles[0]);
    oCrystal.vSetCell(fOrigCellParams);
    oCrystal.nCalcOrientMatrix();
    oCrystal.vGetOrientMatrix(&fOrientMat[0][0]);
    // Multiply the orient matrix from the right.
    vMulMat3DMat3D(fOrientMat,fTransMat,fNewOrientMat);

    // Set the new orientation matrix.
    oCrystal.nSetOrientMatrix(&fNewOrientMat[0][0]);
    // Read the new cell out.
	oCrystal.vGetCell(fCellParams);
	oCrystal.vGetOrientAngles(fOrientAngles);


    // Transform the sigmas.

    {
        double fAvgSigmaPerAngle;
        double fNewCellSigmas[3];

        // Find the average value for the sigma/angle
        for (nx=0, fAvgSigmaPerAngle = 0.0 ; nx<3; nx++) 
            fAvgSigmaPerAngle+=fCellSigmas[nx+3]/fOrigCellParams[nx+3];
        fAvgSigmaPerAngle/=3.0;

        for (nx=0;nx<3;nx++)
            fCellSigmas[nx+3]=fCellParams[nx+3]*fAvgSigmaPerAngle;
        
        // We are using a "quadrature" transformation.  That is,
        // we add errors in quadrature.

        for (nx=0;nx<3;nx++) {
            f0=0.0;
            for (ny=0;ny<3;ny++) 
                f0+=(fReindexMat[ny][nx]*fCellSigmas[ny])*(fReindexMat[ny][nx]*fCellSigmas[ny]);
            fNewCellSigmas[nx]=sqrt(f0);
        };
        for (nx=0;nx<3;nx++)
            fCellSigmas[nx]=fNewCellSigmas[nx];

    };


    printf("New Crystal Parameters:  [%4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %4.2f]\n",
            fCellParams[0],fCellParams[1],fCellParams[2],
            fCellParams[3],fCellParams[4],fCellParams[5]
            );
    printf("New Cell Sigmas:         [%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ]\n",
        fCellSigmas[0],fCellSigmas[1],fCellSigmas[2],
        fCellSigmas[3],fCellSigmas[4],fCellSigmas[5]
        );

    printf("New Crystal Orientation: [%4.2f %4.2f %4.2f]\n",
        fOrientAngles[0],fOrientAngles[1],fOrientAngles[2]); 
        
    return 0;
};

int Claue::nLaueCheck(Creflnlist& oRefIn,bool& bReIndexingDone,bool bPrompt)
{
    Cstring sTemp;
	int nx,ny;
	int nStat;
	int nLaueGroup;
	int nSel;
    int nBestSel;     // The best row in the table.
    int nPass;

    int     nIdealReflnsPerGroup[LAUE_TOTAL] = { 2,4,8,8,16,6,12,12,12,24,24,48 };    

	double f0,f1;
	double fScale;
	double fRmergeRaw,fRmergeScale;
    bool bGetUnique;
	int nCell;
    int a3nHKLOffend[3];					// Used to compute reflections that do not index properly.
    itr<tagLaueTable>	aoLaueTable;		// Laue table entries.
	tagLaueTable		oTableEntry;		// A temporary.
    itr<int>			aa3nHKLSave[3];

	Creflnlist oRef;
	Cspacegroup oSpace;
	CCellReduce  oReduce;

	oRef.nExpandRefln(oRefIn);
	oRef.nExpandGetField(oRef.ms_snPackedHKL);
	oRef.nInsertListFrom(oRefIn);	
	oReduce.m_fResidualMax = fResidMax;
	Ccrystal oCrystalIn;
	Ccrystal oCrystalOut;
	double a3x3fOrientMatIn[3][3];
	double a3x3fOrientMatOut[3][3];
	double a3x3fReindexMat[3][3];


    // Variables used in Rmerge calculation.
    int nGroups;
    int nGroup;
    int nLastHKL,nThisHKL;
    int nStartIndex;
    int nEndIndex;
    int nRefSort,nRef;
    double fAvgGroupSize;
    double fGroupAvg;
    double fNumer;
    double fDenom;
    bool bLocalAnom;				// Possibly reset in the loop to be false if we don't have enough reflections.
	bool bComputeScaleForLaue;		// Set to false on a second pass.
    bool bOutliersReject;           // Are we rejecting outliers yet?
    

    // Variables for scaling.
    itr<double> afSumNumer;
    itr<double> afSumDenom;
    itr<double> afSumVarianceNumer;
    itr<double> afSumVarianceDenom;
    itr<double> afCorrMat;
    itr<double> afCorrMatVariance;
    itr<double> afRMerge;
    itr<double> afRMergePrev;
    itr<int>    anRMergeSort;
	itr<int>	anReject;
    itr<int>    anBatches;
    itr<int>    anBatchStart;
    itr<int>    anBatchEnd;
    itr<double> afScale;
    itr<int>    anBatchIndex;
	itr<int>	anBatchSort;
    itr<Cstring> asBatchName;
    int         nNumBatches;
    int         nBatchi1,nBatchi2;
    int         nBatch,nBatch1,nBatch2;
    int         nIndex;
	int			nBatchCompression;
    double      fMaxRMerge;
    double      fRejectPercent = 0.01;

    
	// Preliminaries if we are doing batch scaling.
    if ((bScaleForLaue) && (oRef.m_nFI_sBatch>=0)) {
		printf("Sorting data on batch ID.\n");
		anBatchSort.setsize(oRef.nGetNumReflns());
		oRef.vSort(eReflnField_Cstring_type,oRef.m_nFI_sBatch,&anBatchSort[0]);
		printf("...done.\n");

        for (nRefSort = 0; nRefSort<oRef.nGetNumReflns();nRefSort++) {
			nRef = anBatchSort[nRefSort];
            if ((!asBatchName.size()) || (asBatchName.last() != oRef[nRef].sGetField(oRef.m_nFI_sBatch))) {
                nx = asBatchName.size();
                asBatchName[nx] = oRef[nRef].sGetField(oRef.m_nFI_sBatch);
            } else {
				nx = asBatchName.size() - 1;
			};
            anBatchIndex[nRef] = nx;
        };
        nNumBatches = asBatchName.size();
		nBatchCompression  = min(nNumBatches,max(1,(int) ceil(nScaleForLaueMinRefs * nNumBatches/ ((double) oRef.nGetNumReflns()))));
		nNumBatches = nNumBatches/nBatchCompression + ((nNumBatches % nBatchCompression)!=0);
		
            
        afSumNumer.setsize(nNumBatches*nNumBatches);
        afSumDenom.setsize(nNumBatches*nNumBatches);
        afSumVarianceNumer.setsize(nNumBatches*nNumBatches);
        afSumVarianceDenom.setsize(nNumBatches*nNumBatches);
        afCorrMat.setsize(nNumBatches*nNumBatches);
        afCorrMatVariance.setsize(nNumBatches*nNumBatches);
        afScale.setsize(nNumBatches);
    } else {
		bScaleForLaue = false;		
	};
    
    
	
	// Do a cell reduction.
	oCrystalIn.vSetCell(fCellParams[0],fCellParams[1],fCellParams[2],fCellParams[3],fCellParams[4],fCellParams[5]);
	oCrystalIn.vSetOrientAngles(fOrientAngles);
	oCrystalIn.nCalcOrientMatrix();
	oCrystalIn.vGetOrientMatrix(&a3x3fOrientMatIn[0][0]);
    Ccrystal oCrystalTemp(oCrystalIn);
	
    oReduce.vConvertToPrimitive(oCrystalIn, a2cLattice);

    oReduce.nAddCrystal(oCrystalIn);
	
    oReduce.nCalcGetBestLattices();
	
    oReduce.nListResults();
    
    ///////////////////////////////////////////////////////////////////////////////////////
    /// RB 01/08/06 This is an attempt to fix an old bug. When an orthorombic solution is
    // present among the solutions, Thad checks the data set for 2/m symmetry along all three axes.
    // The problem arises when axis 'b' is shorter than axis 'a', so the order (a,b,c) is different 
    // for monoclinic and orthorombic cells, because of the standard presentation rule.
    // Since the tests for 'a' and 'c' are performed in the "orthorombic" block of code, while 
    // the test for 'b' is performed in the "monoclinic" block of code, this causes confusion 
    // as to which axis is which. So the fix is to detect if there is an orthorombic solution
    // ahead of time and then perform the 'a' and 'c' tests for 2/m in the "monoclinic" block of
    // code.
    
    bool    b_oP_PassedLatticeResidualTest = false;
    for(nCell = 0; nCell < oReduce.nGetNumSolutions(); nCell++)
    {
		if( oReduce.oGetSolutions()[nCell].fGetIntensity() <= oReduce.m_fResidualMax )
        {
            if( (sTemp = oReduce.sGetLattice(nCell)) == "oP" )
            {
                b_oP_PassedLatticeResidualTest = true;
                break;
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////

    // Starting from the bottom of the list, we extract all cells that appear to have reasonably close lattice residuals.
	for(nCell = oReduce.nGetNumSolutions() - 1; nCell >= 0; nCell--)
    {
		if( oReduce.oGetSolutions()[nCell].fGetIntensity() <= oReduce.m_fResidualMax )
        {
			oReduce.nLoadSolution(nCell, &oCrystalOut);

			strcpy(oTableEntry.acLattice, (sTemp = oReduce.sGetLattice(nCell)).string());
			
            oReduce.fGetOrientInfo(oCrystalIn, oCrystalOut, &oCrystalTemp);
			
            oCrystalOut.nCalcOrientMatrix();
			oCrystalOut.vGetOrientMatrix(&a3x3fOrientMatOut[0][0]);
			
            nCalcReindex(a3x3fOrientMatIn,a3x3fOrientMatOut,a3x3fReindexMat);
			
            oCrystalOut.vGetCell(&oTableEntry.a6fCell[0]);
			
            vCopyVecND(3*3,&a3x3fReindexMat[0][0],&oTableEntry.a3x3fReindexMat[0][0]);
            
            fInvMat3D(&oTableEntry.a3x3fReindexMat[0][0],&oTableEntry.a3x3fReindexMatInv[0][0]);
			
            oTableEntry.pfOrientReindexMat = NULL;
			oTableEntry.pfOrientReindexMatInv = NULL;
							
			if (sTemp == "aP")
            {
				oTableEntry.nLaueGroup = LAUE_1_BAR;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;
			} 
            else if ((sTemp == "mC") || (sTemp == "mP"))
            {
				oTableEntry.nLaueGroup = LAUE_2_OVER_M;
				oTableEntry.cUniqueAxis = 'b';
				aoLaueTable + oTableEntry;
			
                if( sTemp == "mP" && b_oP_PassedLatticeResidualTest )
                {
				    oTableEntry.nLaueGroup = LAUE_2_OVER_M;
				    oTableEntry.cUniqueAxis = 'a';
				    strcpy(oTableEntry.acLattice,"mP");
				    oTableEntry.pfOrientReindexMat = &g_fMatReindexAB[0][0];
				    oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				    aoLaueTable + oTableEntry;
    				
                    oTableEntry.nLaueGroup = LAUE_2_OVER_M;
				    oTableEntry.cUniqueAxis = 'c';
				    strcpy(oTableEntry.acLattice,"mP");
				    oTableEntry.pfOrientReindexMat = &g_fMatReindexBC[0][0];
				    oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				    aoLaueTable + oTableEntry;
                }
           }
            else if (sTemp == "oP")
            {
				/// RB 01/08/06 Commenting these tests out. Please see explanation above.
				//oTableEntry.nLaueGroup = LAUE_2_OVER_M;
				//oTableEntry.cUniqueAxis = 'a';
				//strcpy(oTableEntry.acLattice,"mP");
				//oTableEntry.pfOrientReindexMat = &g_fMatReindexAB[0][0];
				//oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				//aoLaueTable + oTableEntry;
				//
                //oTableEntry.nLaueGroup = LAUE_2_OVER_M;
				//oTableEntry.cUniqueAxis = 'c';
				//strcpy(oTableEntry.acLattice,"mP");
				//oTableEntry.pfOrientReindexMat = &g_fMatReindexBC[0][0];
				//oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				//aoLaueTable + oTableEntry;
				
                oTableEntry.nLaueGroup = LAUE_MMM;
				oTableEntry.cUniqueAxis = '-';
				strcpy(oTableEntry.acLattice,"oP");
				oTableEntry.pfOrientReindexMat = NULL;
                oTableEntry.pfOrientReindexMatInv = NULL;
				aoLaueTable + oTableEntry;						
			} 
            else if ((sTemp == "oC") || (sTemp == "oF") || (sTemp == "oI"))
            {
				oTableEntry.nLaueGroup = LAUE_MMM;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;				
			} 
            else if ((sTemp == "tP") || (sTemp == "tI"))
            {
				oTableEntry.nLaueGroup = LAUE_4_OVER_M;
				oTableEntry.cUniqueAxis = 'c';
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_4_OVER_MMM;
				oTableEntry.cUniqueAxis = 'c';
				aoLaueTable + oTableEntry;				
			} 
            else if (sTemp == "cP")
            {
				oTableEntry.nLaueGroup = LAUE_4_OVER_M;
				oTableEntry.cUniqueAxis = 'a';
				strcpy(oTableEntry.acLattice,"tP");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexAC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_4_OVER_M;
				oTableEntry.cUniqueAxis = 'b';
				strcpy(oTableEntry.acLattice,"tP");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexBC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_4_OVER_MMM;
				oTableEntry.cUniqueAxis = 'a';
				strcpy(oTableEntry.acLattice,"tP");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexAC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_4_OVER_MMM;
				oTableEntry.cUniqueAxis = 'b';
				strcpy(oTableEntry.acLattice,"tP");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexBC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_M_3_BAR;
				oTableEntry.cUniqueAxis = '-';
				oTableEntry.pfOrientReindexMat = NULL;
                oTableEntry.pfOrientReindexMatInv = NULL;
				strcpy(oTableEntry.acLattice,"cP");
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_M_3_BAR_M;
				oTableEntry.cUniqueAxis = '-';
				strcpy(oTableEntry.acLattice,"cP");
				oTableEntry.pfOrientReindexMat = NULL;
                oTableEntry.pfOrientReindexMatInv = NULL;
				aoLaueTable + oTableEntry;		
			}
            else if (sTemp == "cI")
            {
				oTableEntry.nLaueGroup = LAUE_4_OVER_M;
				oTableEntry.cUniqueAxis = 'a';
				strcpy(oTableEntry.acLattice,"tI");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexAC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_4_OVER_M;
				oTableEntry.cUniqueAxis = 'b';
				strcpy(oTableEntry.acLattice,"tI");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexBC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_4_OVER_MMM;
				oTableEntry.cUniqueAxis = 'a';
				strcpy(oTableEntry.acLattice,"tI");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexAC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;				
				oTableEntry.nLaueGroup = LAUE_4_OVER_MMM;
				oTableEntry.cUniqueAxis = 'b';
				strcpy(oTableEntry.acLattice,"tI");
				oTableEntry.pfOrientReindexMat = &g_fMatReindexBC[0][0];
				oTableEntry.pfOrientReindexMatInv = oTableEntry.pfOrientReindexMat;
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_M_3_BAR;
				oTableEntry.cUniqueAxis = '-';
				strcpy(oTableEntry.acLattice,"cI");
				oTableEntry.pfOrientReindexMat = NULL;
                oTableEntry.pfOrientReindexMatInv = NULL;
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_M_3_BAR_M;
				oTableEntry.cUniqueAxis = '-';
				strcpy(oTableEntry.acLattice,"cI");
				oTableEntry.pfOrientReindexMat = NULL;
                oTableEntry.pfOrientReindexMatInv = NULL;
				aoLaueTable + oTableEntry;		
			} 
            else if (sTemp == "cF")
            {
				oTableEntry.nLaueGroup = LAUE_M_3_BAR;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_M_3_BAR_M;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;		
			}
            else if ((sTemp == "hP") || (sTemp == "hR"))
            {
				// Trigonal and Hexagonal (Assume hexagonal basis for trigonal)
				oTableEntry.nLaueGroup = LAUE_3_BAR;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_3_BAR_M_1;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_3_BAR_1_M;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_6_OVER_M;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;		
				oTableEntry.nLaueGroup = LAUE_6_OVER_MMM;
				oTableEntry.cUniqueAxis = '-';
				aoLaueTable + oTableEntry;		

				// The 'hP' lattice will give rise to 3 other oC lattices.
				if (sTemp == "hP")
                {
					for (nx=0;nx< aoLaueTable.size();nx++)
                    {
						if (((Cstring) "oC") == aoLaueTable[nx].acLattice)
							break;
					}
					
                    if (nx < aoLaueTable.size())
						aoLaueTable.remove(nx,1);
					
                    // Add three new operators.
					oTableEntry.nLaueGroup = LAUE_MMM;
					
					strcpy(oTableEntry.acLattice,"oC");
					oTableEntry.pfOrientReindexMat = &g_fMatReindexhPoC_1[0][0];
					oTableEntry.pfOrientReindexMatInv = &g_fMatReindexoChP_1[0][0];
                    oTableEntry.cUniqueAxis = '1';
					aoLaueTable + oTableEntry;				
					oTableEntry.pfOrientReindexMat = &g_fMatReindexhPoC_2[0][0];
					oTableEntry.pfOrientReindexMatInv = &g_fMatReindexoChP_2[0][0];
                    oTableEntry.cUniqueAxis = '2';
					aoLaueTable + oTableEntry;				
					oTableEntry.pfOrientReindexMat = &g_fMatReindexhPoC_3[0][0];
					oTableEntry.pfOrientReindexMatInv = &g_fMatReindexoChP_3[0][0];
                    oTableEntry.cUniqueAxis = '3';
					aoLaueTable + oTableEntry;				
				}
			}
		}
	}
    
    if (eMaxLaueCheck != LAUE_UNKNOWN)
    {
        for (nx=0;nx<aoLaueTable.size();nx++)
        {
            if (aoLaueTable[nx].nLaueGroup>eMaxLaueCheck)
            {
                aoLaueTable.remove(nx,1);
                nx--;
            }
        }
    }
	
    // Sort the Laue table entries.
	if( aoLaueTable.size() > 0 ) // safety
        qsort(&aoLaueTable[0],aoLaueTable.size(),sizeof(tagLaueTable),nLaueTableCompare);

    printf("\nLaue Check");
    printf("\n------------------------------------------------------------------------------");
    printf("\n%s %s %s  %s  %s  %s  %s %s %s %s %s",  "   Laue", "Unique","Latt"," Refln","  Non","Calc","  Obs","*Rmerge","**Rmerge","Incres"," Pass?");
    printf("\n%s %s %s  %s  %s  %s  %s %s %s %s %s",  "  class", "  axis","used","groups","index","mult"," mult","    raw","  scaled","  frac","      ");
    printf("\n------------------------------------------------------------------------------\n");

#ifdef SSI_PC
	CCrclHelper::GetInstance()->vResetLaueTable();
#endif

    // Save the HKL values.
    for(nx = 0; nx < oRef.nGetNumReflns(); nx++)
    {
        aa3nHKLSave[0] + oRef[nx].nGetH();
        aa3nHKLSave[1] + oRef[nx].nGetK();
        aa3nHKLSave[2] + oRef[nx].nGetL();
    }


	for (nBestSel = -1,nSel = 0; nSel < aoLaueTable.size(); nSel++)
    {
		tagLaueTable& oLaueEntry = aoLaueTable[nSel];
		
		nLaueGroup = oLaueEntry.nLaueGroup;
		
		strcpy(oLaueEntry.acLaueClass,g_cpLaueNameSmall[nLaueGroup]);
		
		oSpace.m_nNumber=ms_nlaue_group[nLaueGroup];
		oSpace.m_sName=oSpace.sNameFromNumber();
		nStat=oSpace.nReadEquivPos();

		if (nStat)
        { 
			cout << "Error Reading Spacegroup #" << oSpace.m_nNumber << "\n"; 
            m_eLaue = LAUE_UNKNOWN;
			return 1; 
		}

		// Apply the reindexing matrix.
        oRef.nReindex(&oLaueEntry.a3x3fReindexMat[0][0]);
		
        if (oLaueEntry.pfOrientReindexMat)
			oRef.nReindex(oLaueEntry.pfOrientReindexMat);
		
		// Sort the reflection list on the space group.  Do not worry about annomolous reflections.
		bLocalAnom = bAnom;
		bComputeScaleForLaue = bScaleForLaue;
		fRmergeRaw = -1.0;
		fRmergeScale = -1.0;
        bOutliersReject = false;
        fMaxRMerge = 2.0;

		oRef.nReduce(oSpace,bLocalAnom,0);
        afRMergePrev.clear();
		
		if (bScaleForLaue)
        {
			memset(&afSumNumer[0],0,sizeof(double)*nNumBatches*nNumBatches);
			memset(&afSumDenom[0],0,sizeof(double)*nNumBatches*nNumBatches);
			memset(&afSumVarianceNumer[0],0,sizeof(double)*nNumBatches*nNumBatches);
			memset(&afSumVarianceDenom[0],0,sizeof(double)*nNumBatches*nNumBatches);
			memset(&afCorrMat[0],0,sizeof(double)*nNumBatches*nNumBatches);
			memset(&afCorrMatVariance[0],0,sizeof(double)*nNumBatches*nNumBatches);
			memset(&afScale[0],0,sizeof(double)*nNumBatches);
		}

        nPass = 0;
		do
        {
			nGroups=0;
			nLastHKL = oRef[oRef.pnGetSortIndex()[0]].nGetField(oRef.m_nFI_nPackedHKL);
			nStartIndex = 0;
			fNumer = 0.0;
			fDenom = 0.0;
			fAvgGroupSize=0;
            
            if (bComputeScaleForLaue)
            {
                for (nBatch = 0; nBatch < nNumBatches; nBatch++)
                    afScale[nBatch] = 1.0;
            }

            
            if (!bOutliersReject)
            {
                afRMerge.setsize(oRef.nGetNumReflns());                
                memset(&afRMerge[0],0,sizeof(double)*oRef.nGetNumReflns());
            }
            
			for (nRefSort=0;nRefSort<oRef.nGetNumReflns()+1;nRefSort++)
            {
				bool bHKLEqual;
				nRef=oRef.pnGetSortIndex()[nRefSort];               
				if (nRefSort<oRef.nGetNumReflns())
					nThisHKL = oRef[nRef].nGetField(oRef.m_nFI_nPackedHKL);
				bHKLEqual = ((((unsigned int) nThisHKL ))==(((unsigned int) nLastHKL )));
				
				if ((nRefSort!=oRef.nGetNumReflns()) && (bHKLEqual))
					continue;
				nEndIndex=nRefSort-1;

                
				// Find the average of the reflections in the group.
				fGroupAvg = 0.0;
				nGroup =0;
				for(nx=nStartIndex; nx <= nEndIndex; nx++)
                {
					fScale = (bScaleForLaue)?afScale[anBatchIndex[oRef.pnGetSortIndex()[nx]]/nBatchCompression]:1.0;
					f0 = oRef[oRef.pnGetSortIndex()[nx]].fGetIntensity()*fScale;
					f1 = oRef[oRef.pnGetSortIndex()[nx]].fGetSigmaI()*fScale;
					if ((f1<=0.0) || (f0/(f1)<fLaueSigma) || ((afRMergePrev.size()>0) && (afRMergePrev[nx]>fMaxRMerge)))
						continue;
					fGroupAvg += f0;
					nGroup++;
				}

                if( bScaleForLaue )
                {
                    // Calculate on-the-fly scale correlation coeffs.
                    -anBatches;
                    -anBatchStart;
                    -anBatchEnd;
                    for (nIndex = nStartIndex; nIndex <= nEndIndex; nIndex++)
                    {
                        nx = anBatchIndex[oRef.pnGetSortIndex()[nIndex]]/nBatchCompression;

                        if ((anBatches.size()) && (nx < anBatches.last()))
                        {
                            printf("ERROR:  In batch scaling.\n");
                            return 1;
                        }
                        
                        if ((!anBatches.size()) || (nx > anBatches.last()))
                        {
                            anBatches + nx;
                            anBatchStart + nIndex;
                            anBatchEnd + nIndex;
                        }
                        else
                            anBatchEnd.last() = nIndex;
                    }
                    
                    for (nBatchi1 = 0; nBatchi1 < anBatches.size(); nBatchi1++)
                    {                
                        for (nBatchi2 = 0; nBatchi2 < anBatches.size(); nBatchi2++)
                        {
                            if (nBatchi1 != nBatchi2)
                            {
                                double fBatch1Intensity = 0.0;
                                double fBatch1Variance  = 0.0;
                                int    nCount1          = 0;
                                double fBatch2Intensity = 0.0;
                                double fBatch2Variance  = 0.0;
                                int    nCount2          = 0;
                                double fIntensity,fSigma;
                                
                                nBatch1 = anBatches[nBatchi1];
                                nBatch2 = anBatches[nBatchi2];
                                
                                // Average and get sigma of all elements in each batch.
                                for (nIndex = anBatchStart[nBatchi1]; nIndex <= anBatchEnd[nBatchi1];nIndex++)
                                {
                                   
                                    fIntensity = oRef[oRef.pnGetSortIndex()[nIndex]].fGetIntensity()*afScale[nBatch1];
                                    fSigma = oRef[oRef.pnGetSortIndex()[nIndex]].fGetSigmaI()*afScale[nBatch1];
									if ((fSigma>0.0) && (fIntensity/fSigma>=fLaueSigma) && ((afRMergePrev.size()==0) || (afRMergePrev[nIndex]<=fMaxRMerge))) {
										fBatch1Intensity += fIntensity/max(1.0,(fSigma*fSigma));
										fBatch1Variance += 1.0/max(1.0,(fSigma*fSigma));
										nCount1 ++;
									}
                                }
                                if (!nCount1)
                                    continue;
                                
                                fBatch1Intensity/=fBatch1Variance;
                                fBatch1Variance = 1.0/max(1e-20,fBatch1Variance);
                                
                                for (nIndex = anBatchStart[nBatchi2]; nIndex <= anBatchEnd[nBatchi2];nIndex++)
                                {
                                    fIntensity = oRef[oRef.pnGetSortIndex()[nIndex]].fGetIntensity()*afScale[nBatch2];
                                    fSigma = oRef[oRef.pnGetSortIndex()[nIndex]].fGetSigmaI()*afScale[nBatch2];
                                    fBatch2Intensity += fIntensity/max(1.0,(fSigma*fSigma));
									if ((fSigma>0.0) && (fIntensity/fSigma>=fLaueSigma) && ((afRMergePrev.size()==0) || (afRMergePrev[nIndex]<=fMaxRMerge))) {
										fBatch2Variance += 1.0/max(1.0,(fSigma*fSigma));
										nCount2 ++;
									}
                                }
                                
                                if (!nCount2)
                                    continue;
                                
                                fBatch2Intensity/=max(1e-20,fBatch2Variance);
                                fBatch2Variance = 1.0/max(1e-20,fBatch2Variance);
                                
                                f0 = 1.0/(fBatch1Variance+fBatch2Variance);
                                if (min(fBatch1Intensity,fBatch2Intensity)/sqrt(max(fBatch1Variance,fBatch2Variance)) > fLaueSigma) {
                                    afSumNumer[nBatch1*nNumBatches + nBatch2] += fBatch1Intensity*fBatch2Intensity*f0;
                                    afSumDenom[nBatch1*nNumBatches + nBatch2] += fBatch1Intensity*fBatch1Intensity*f0;
                                    afSumVarianceNumer[nBatch1*nNumBatches + nBatch2] += (fBatch1Variance*fBatch2Intensity*fBatch2Intensity+fBatch2Variance*fBatch1Intensity*fBatch1Intensity)*f0*f0;
                                    afSumVarianceDenom[nBatch1*nNumBatches + nBatch2] += 2.0*fBatch1Variance*fBatch1Intensity*fBatch1Intensity*f0*f0;
                                }
                            }
                        }
                    }
                }
				
				if (nGroup>1)
                {
					fGroupAvg/=nGroup;
					for (nx=nStartIndex;nx<=nEndIndex;nx++)
                    {
						fScale = (bScaleForLaue)?afScale[anBatchIndex[oRef.pnGetSortIndex()[nx]]/nBatchCompression]:1.0;
						f0 = oRef[oRef.pnGetSortIndex()[nx]].fGetIntensity()*fScale;
						f1 = oRef[oRef.pnGetSortIndex()[nx]].fGetSigmaI()*fScale;
						if ((f1 <= 0.0) || (f0/f1<fLaueSigma) || ((afRMergePrev.size()>0) && (afRMergePrev[nx]>fMaxRMerge)))
							continue;
						fNumer+=fabs(f0 - fGroupAvg);
						fDenom+=fabs(f0);
                        afRMerge[nx] = fNumer/max(0.000001,fDenom);
					}
					
                    nGroups++;
					fAvgGroupSize += nGroup;
				}
				
				nStartIndex=nRefSort;
				if (nRefSort<oRef.nGetNumReflns())
					nLastHKL =  oRef[nRef].nGetField(oRef.m_nFI_nPackedHKL);               
			}
			
            fAvgGroupSize/=max(1,nGroups);


            if (bScaleForLaue)
            {
                // Now calculate correlation constants K[i,j] where i < j and K[i,j]*Batch(i) = Batch(j) 
                // For Matrix X we assume X[i,j] = X[i*nNumBatches + j]
                
                double fNumer,fDenom,fVarNumer,fVarDenom;
                for (nBatch1 = 0; nBatch1 < nNumBatches; nBatch1++)
                {
                    for (nBatch2 = 0; nBatch2 < nNumBatches; nBatch2++)
                    {
                        if (afSumDenom[nBatch1*nNumBatches + nBatch2])
                        {
                            afCorrMat[nBatch1*nNumBatches + nBatch2] = (fNumer = afSumNumer[nBatch1*nNumBatches + nBatch2])/(fDenom = max(1e-20,(afSumDenom[nBatch1*nNumBatches + nBatch2])));
                            fVarNumer = afSumVarianceNumer[nBatch1*nNumBatches + nBatch2];
                            fVarDenom = afSumVarianceDenom[nBatch1*nNumBatches + nBatch2];
                            afCorrMatVariance[nBatch1*nNumBatches + nBatch2] = fVarNumer/max(1e-20,(fDenom*fDenom)) + fVarDenom*fNumer*fNumer/max(1e-20,(fDenom*fDenom*fDenom*fDenom));
                        }
                    }
                }
                
                // Apply scale factors.
                for (nBatch1 = 0; nBatch1 < nNumBatches; nBatch1++)
                {
                    double fScaleNumer = 0.0;
                    double fScaleDenom = 0.0;
                    double fScale;
                    for (nBatch2 = 0; nBatch2 < nNumBatches; nBatch2++) {
                        if ((nBatch2 != nBatch1) && (afCorrMat[nBatch2*nNumBatches + nBatch1])) {
                            fScaleNumer += afCorrMat[nBatch2*nNumBatches + nBatch1]/max(1e-20,afCorrMatVariance[nBatch2*nNumBatches + nBatch1]);
                            fScaleDenom += 1.0/max(1e-20,afCorrMatVariance[nBatch2*nNumBatches + nBatch1]);
                        }; 
                    };
                    fScale = fScaleNumer/max(1e-20,fScaleDenom);
					if (fScale>0.0)
						afScale[nBatch1] *= 1.0/fScale;
                }
            }

           
            if ((bLocalAnom) && (nGroups<50) && (nLaueGroup==0)) {
				bLocalAnom = false;
				oRef.nReduce(oSpace,bLocalAnom,0);
				continue;
            } else if ((!bOutliersReject) && (bScaleForLaue)) {
                // Sort the R-merge array.
                g_pfCmpDoubles = &afRMerge[0];
                anRMergeSort.setsize(afRMerge.size());
                for (nx = 0; nx < anRMergeSort.size(); nx++)
                    anRMergeSort[nx] = nx;
                
                if( anRMergeSort.size() > 0 ) // safety check
                    qsort(&anRMergeSort[0],anRMergeSort.size(),sizeof(int),double_cmp_rel);
                
                for (nx = 0; nx < anRMergeSort.size(); nx++) {
                    if (afRMerge[anRMergeSort[nx]]>0.0)
                        break;
                };
                
                ny = (int)((double)(anRMergeSort.size() - nx) * fRejectPercent);
                
                fMaxRMerge = afRMerge[anRMergeSort[anRMergeSort.size() - 1 - ny]];
                
                afRMergePrev = afRMerge;
                
                bOutliersReject = true;
                
                continue;
			} else if (bComputeScaleForLaue) {
				fRmergeRaw = fNumer/max(1.0,fDenom);
				bComputeScaleForLaue = false;
                continue;
			} else {
				fRmergeScale = fNumer/max(1.0,fDenom);
				break;
			};
			
		} while (++nPass);


		// Apply the reindexing matrix.
		if (oLaueEntry.pfOrientReindexMatInv)
			oRef.nReindex(oLaueEntry.pfOrientReindexMatInv);
        
        oRef.nReindex(&oLaueEntry.a3x3fReindexMatInv[0][0]);

		// Check for indexed/non-indexed reflections.
        a3nHKLOffend[0] = 0;
        a3nHKLOffend[1] = 0;
        a3nHKLOffend[2] = 0;
        for (nRef = 0; nRef < oRef.nGetNumReflns();nRef++) {
            if (oRef[nRef].nGetH() != aa3nHKLSave[0][nRef])
                a3nHKLOffend[0]++;
            
            oRef[nRef].vSetH(aa3nHKLSave[0][nRef]);
            if (oRef[nRef].nGetK() != aa3nHKLSave[1][nRef])
                a3nHKLOffend[1]++;
            
            oRef[nRef].vSetK(aa3nHKLSave[1][nRef]);
            if (oRef[nRef].nGetL() != aa3nHKLSave[2][nRef])
                a3nHKLOffend[2]++;
            
            oRef[nRef].vSetL(aa3nHKLSave[2][nRef]);
        };
        oLaueEntry.nNonIndexingRefs = a3nHKLOffend[0] + a3nHKLOffend[1] + a3nHKLOffend[2];

		oLaueEntry.fRmerge = fRmergeScale;
		oLaueEntry.fRmergeRaw = (fRmergeRaw==-1.0)?fRmergeScale:fRmergeRaw;
		oLaueEntry.fRmergeMin = min(oLaueEntry.fRmergeRaw,oLaueEntry.fRmerge);
		oLaueEntry.nGroups = nGroups;
		oLaueEntry.fCalcMult = nIdealReflnsPerGroup[nLaueGroup]/((bLocalAnom)?2:1);
		oLaueEntry.fObserveMult = fAvgGroupSize;
		oLaueEntry.fRatioTo1Bar = 0.0;
		
		// Find the last valid entry.
		// Fill in the ratio to that entry if we have it.
		for (nx=nSel - 1;(nx>=0);nx--) {
			if ((!(aoLaueTable[nx] == oLaueEntry)) && (aoLaueTable[nx].nPass == oLaueEntry.ePass)) {
				oLaueEntry.fRatioTo1Bar = oLaueEntry.fRmergeMin/aoLaueTable[nx].fRmergeMin;
				break;
			};
		};
		
        Cstring     sPassQualifier("");
		
        if( oLaueEntry.nGroups < 10 )
        {
            oLaueEntry.nPass = oLaueEntry.eInvalid;
            sPassQualifier = "[BAD]";
        }
        else if (nLaueGroup == 0)
        {
            oLaueEntry.nPass = oLaueEntry.ePass;  // Triclinic always passes...
            sPassQualifier = "N/A";               // so it does not make sense to say it has passed
        }
        else if( oLaueEntry.fRmergeMin <= fRMergeGreat )
        {
            oLaueEntry.nPass = oLaueEntry.ePass;
            sPassQualifier = "[PASS]";
        }
        else if( oLaueEntry.fRmergeMin > fRMergeAwful )
        {
            if( 0.0 < oLaueEntry.fRatioTo1Bar && oLaueEntry.fRatioTo1Bar < fMaxRMergeRatio )
            {
                oLaueEntry.nPass = oLaueEntry.ePass;  // we let it pass, because it's close enough to the lower symmetry
                sPassQualifier = "[POOR]";  // poor, because Rmerge is poor
            }
            else
            {
                oLaueEntry.nPass = oLaueEntry.eFail;
                sPassQualifier = "------";
            }
        }
        else if( oLaueEntry.fRatioTo1Bar == 0.0 )
        {
            oLaueEntry.nPass = oLaueEntry.eFail;
            sPassQualifier = "------";
        }
        else if (oLaueEntry.fRatioTo1Bar < fMaxRMergeRatio)
        {
            oLaueEntry.nPass = oLaueEntry.ePass;
            sPassQualifier = "[PASS]";
        }
        else
        {
            oLaueEntry.nPass = oLaueEntry.eFail;
            sPassQualifier = "------";
        }
		
        strcpy(oLaueEntry.acPassFail, sPassQualifier.string()); 
        
        aoLaueTable[nSel] = oLaueEntry;
		
		
#ifdef SSI_PC
		CCrclHelper::GetInstance()->vAddLaueTableRow( oLaueEntry );
#endif

        Cstring sRmergeRaw;
        Cstring sRmergeScale;
        char pcBuf[20];
        
        sRmergeRaw = " ";
        sRmergeScale = " ";
        if (bScaleForLaue) {
            if (oLaueEntry.fRmergeRaw <= oLaueEntry.fRmerge)  {
                sRmergeRaw = "[";
            } else
                sRmergeScale = "[";
        };
        sprintf(pcBuf,"%.2lf",oLaueEntry.fRmergeRaw);
        sRmergeRaw += pcBuf;
        sprintf(pcBuf,"%.2lf",oLaueEntry.fRmerge);
        sRmergeScale += pcBuf;
        if (bScaleForLaue) {
            if (oLaueEntry.fRmergeRaw <= oLaueEntry.fRmerge) {
                sRmergeRaw += "]";
                sRmergeScale += " ";
            } else {
                sRmergeScale += "]";
                sRmergeRaw += " ";
            };
        };

		
		printf("%7s      %c   %s  ",oLaueEntry.acLaueClass,oLaueEntry.cUniqueAxis,oLaueEntry.acLattice);               
		printf("%6d  %5d  %4d  %5.2f %7s %8s ",
			oLaueEntry.nGroups,
            oLaueEntry.nNonIndexingRefs,
			(int) oLaueEntry.fCalcMult,
			oLaueEntry.fObserveMult,
            sRmergeRaw.string(),
            sRmergeScale.string()		
			);
		if (oLaueEntry.fRatioTo1Bar == 0.0)
			printf("%6s","N/A");
		else
			printf("%6.2lf",oLaueEntry.fRatioTo1Bar);
		printf(" %6s\n",oLaueEntry.acPassFail);
	    
        //printf(" a=%.3f b=%.3f c=%.3f alpha=%.3f beta=%.3f gamma=%.3f\n", 
                             //oLaueEntry.a6fCell[0], 
                             //oLaueEntry.a6fCell[1], 
                             //oLaueEntry.a6fCell[2], 
                             //oLaueEntry.a6fCell[3], 
                             //oLaueEntry.a6fCell[4], 
                             //oLaueEntry.a6fCell[5]);

        if (oLaueEntry.nPass == oLaueEntry.ePass){ 
			// If nLaueGroup is the same as the "best" nLaueGroup then we need to compare Rmerges.
			if ((nBestSel == -1) || (aoLaueTable[nSel].nLaueGroup != aoLaueTable[nBestSel].nLaueGroup) || (aoLaueTable[nSel].fRmergeMin < aoLaueTable[nBestSel].fRmergeMin)) {
				if (nBestSel == -1) {
					nBestSel = nSel;
				} else if ((aoLaueTable[nBestSel].nLaueGroup == LAUE_2_OVER_M) && (aoLaueTable[nSel].nLaueGroup == LAUE_2_OVER_M)) {
					if ((aoLaueTable[nSel].cUniqueAxis == 'b') || (aoLaueTable[nSel].fRmergeMin<aoLaueTable[nBestSel].fRmergeMin/fMaxRMergeRatio)) {
						nBestSel = nSel;
					};
					
				} else if ((aoLaueTable[nBestSel].nLaueGroup == LAUE_4_OVER_M) && (aoLaueTable[nSel].nLaueGroup == LAUE_4_OVER_M)) {
					if ((aoLaueTable[nSel].cUniqueAxis == 'c') || (aoLaueTable[nSel].fRmergeMin<aoLaueTable[nBestSel].fRmergeMin/fMaxRMergeRatio)) {
						nBestSel = nSel;
					};
				} else if ((aoLaueTable[nBestSel].nLaueGroup == LAUE_4_OVER_MMM) && (aoLaueTable[nSel].nLaueGroup == LAUE_4_OVER_MMM)) {
					if ((aoLaueTable[nSel].cUniqueAxis == 'c') || (aoLaueTable[nSel].fRmergeMin<aoLaueTable[nBestSel].fRmergeMin/fMaxRMergeRatio)) {
						nBestSel = nSel;
					};
				} else {
					nBestSel = nSel;
				};
			};
		};
		if (((oLaueEntry.nLaueGroup==LAUE_2_OVER_M) && (oLaueEntry.cUniqueAxis != 'b'))  ||
			((oLaueEntry.nLaueGroup==LAUE_4_OVER_M) && (oLaueEntry.cUniqueAxis != 'c')) ||
			((oLaueEntry.nLaueGroup==LAUE_4_OVER_MMM) && (oLaueEntry.cUniqueAxis != 'c')))                
			oLaueEntry.bUserSelectable = FALSE;
		else
			oLaueEntry.bUserSelectable = TRUE;
	
        

    };
    printf("\n------------------------------------------------------------------------------\n");
    
    if (nBestSel != -1) {
        f0=aoLaueTable[nBestSel].fRmergeMin;
        printf("%7s selected with Rmerge of .......................... %5.2f\n\n",g_cpLaueNameSmall[aoLaueTable[nBestSel].nLaueGroup],f0);
        printf(" *   Rmerge is from reflections with I/sig >= %.2lf.\n     %.2lf percent of data rejected.\n",fLaueSigma,fRejectPercent*100.0);
        if (bScaleForLaue)
            printf(" **  Rmerge-scaled uses batch scaling.\n     Lower of Rmerge-scaled or Rmerge-raw used (enclosed in brackets)\n");
        else
            printf(" **  Rmerge-scaled value will equal Rmerge-raw value\n     because batch scaling has been disabled.\n");
    } else {
        printf("Laue test inconclusive.\n");
    };
    printf("\n");

    if( bPrompt ) 
    {
		char        cUniqueAxis = '\0';
#ifdef SSI_PC
        if( nBestSel != -1 )
        {
            CCrclHelper::GetInstance()->vBuildCmdLaueTable( g_cpLaueNameSmall[aoLaueTable[nBestSel].nLaueGroup], 
															aoLaueTable[nBestSel].acLattice,
															(float)aoLaueTable[nBestSel].fRmergeMin );
        } 
        else
        {
            CCrclHelper::GetInstance()->vBuildCmdLaueTable( "none", "none", -1.0f );
        }

        printf("Enter your choice [cancel]:\n");
#else
        if( nBestSel != -1 )
        {
            printf("Enter your choice [%s]: ", g_cpLaueNameSmall[aoLaueTable[nBestSel].nLaueGroup]);
        }		
        else
        {
            printf("Enter your choice [cancel]: ");
        }
        
        fflush(stdout);
#endif

        getline(cin, sTemp);

#ifdef SSI_PC
        printf("%s\n\n", sTemp.string());
#endif 

        if( sTemp.length() > 0 ) 
        {
            if( sTemp == "cancel" )
			{
				printf("Laue Check canceled by user request.\n");
				m_eLaue = LAUE_UNKNOWN;
				return 1;
			}
			
			// Else parse sTemp
			std::vector<MSCString>          saArray;
			sTemp.nListToVector(saArray, " "); // 0 - Laue group; 1 - Unique Axis; 2 - Lattice used
        
			int     nNumberOfInputParams = (int)saArray.size();

			if( nNumberOfInputParams == 0 )
			{
				printf("Invalid input. Laue Check canceled.\n");
				m_eLaue = LAUE_UNKNOWN;
				return 1;
			}

			Cstring     sLaueGroup = saArray[0];  
			cUniqueAxis = nNumberOfInputParams > 1 && saArray[1].length() > 0 ? saArray[1].GetAt(0) : '-';
			Cstring     sLatticeUsed = nNumberOfInputParams > 2 ? saArray[2] : Cstring("");

			for(nx = 0; nx < aoLaueTable.size(); nx++)
			{
				if( sLatticeUsed != "" )
				{
					if( sLaueGroup   == g_cpLaueNameSmall[aoLaueTable[nx].nLaueGroup] && 
						cUniqueAxis  == aoLaueTable[nx].cUniqueAxis &&
						sLatticeUsed == aoLaueTable[nx].acLattice )
						break;
				}
				else
				{
					if( sLaueGroup  == g_cpLaueNameSmall[aoLaueTable[nx].nLaueGroup] && 
						cUniqueAxis == aoLaueTable[nx].cUniqueAxis )
						break;
				}
			}
			
			if( nx < aoLaueTable.size() )
				nBestSel = nx;
		}
        else if( nBestSel == -1 )
        {
		    m_eLaue = LAUE_UNKNOWN;
		    return 1;
        }
    }


	double a3x3fTotalReindex[3][3];
	double a6fOldCell[6];
	double a6fOldCellSigma[6];
	double a3fOldOrientAngles[3];

	if (nBestSel != -1) {
		vCopyVecND(6,fCellParams,a6fOldCell);
		vCopyVecND(6,fCellSigmas,a6fOldCellSigma);
		vCopyVecND(3,fOrientAngles,a3fOldOrientAngles);
		
        if( aoLaueTable[nBestSel].pfOrientReindexMat )
        {
			vMulMatNDMatND(3, 
                           aoLaueTable[nBestSel].pfOrientReindexMat,
                           &aoLaueTable[nBestSel].a3x3fReindexMat[0][0],
                           &a3x3fTotalReindex[0][0]);
			
            nCalcCell(a3x3fTotalReindex);
		} 
        else 
        {
			vCopyMat3D(&aoLaueTable[nBestSel].a3x3fReindexMat[0][0], &a3x3fTotalReindex[0][0]);
			nCalcCell(a3x3fTotalReindex);
		}
		
        // Check to see if we have an identity matrix.
		bReIndexingDone=FALSE;
		for (nx = 0; nx < 3; nx++) 
        {
			for (ny = 0; ny < 3; ny++) 
            {
				if (ABS(a3x3fTotalReindex[nx][ny] - (nx == ny))>0.01)
					bReIndexingDone=TRUE;
			}
		}

		
		if( bPrompt && bReIndexingDone ) 
        {
			// This print statement required by StructureStudio.
			printf("Reindexing is recommended.\n");
			printf("Do you wish to Reindex? [Y]/N: ");
#ifdef SSI_PC
			CCrclHelper::GetInstance()->vBuildCmdCellReindexForAxis(' ');
#else
			fflush(stdout);
#endif
			getline(cin,sTemp);
#ifdef SSI_PC
			printf("%s\n", sTemp.string());
#endif
			if (((sTemp.length()==0) || (sTemp=="Y") || (sTemp=="y"))) 
				bGetUnique = TRUE;
			else 
            {
				bGetUnique = FALSE;
				bReIndexingDone=FALSE;
			}
			
		}
        else
        {
#ifdef SSI_PC
            if( !bReIndexingDone )
            {
		        // Let CrystalClear know that no reindexing is required
                CCrclHelper::GetInstance()->vBuildCmdCellReindexForAxis('N');
                CCrclHelper::GetInstance()->SendTclCmd();
            }
#endif
            bGetUnique = TRUE;
        }

		if( bGetUnique && bReIndexingDone )
        {
			oRefIn.nReindex(&a3x3fTotalReindex[0][0]);
		}

        for(nx=0; nx < 6; nx++)
        {
			if(2.0*ABS(a6fOldCell[nx] - fCellParams[nx])/ABS(a6fOldCell[nx] + fCellParams[nx]) > 0.001)
			    break;
		}
		
        if (nx != 6)
        {
			int     nResponse = 1; // 1 if want to use new cell, 0 if otherwise.
			
			printf("\n\n");
			
            printf("WARNING: The Laue class '%s' is incompatible with the present cell.\n\n",
                   g_cpLaueNameSmall[aoLaueTable[nBestSel].nLaueGroup]);
			
            printf("Present Cell:          [ %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f]\n",
				a6fOldCell[0],a6fOldCell[1],a6fOldCell[2],a6fOldCell[3],a6fOldCell[4],a6fOldCell[5]);
			
            printf("Suggested Cell:        [ %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f]\n",
				fCellParams[0],fCellParams[1],fCellParams[2],fCellParams[3],fCellParams[4],fCellParams[5]);
			
			if (bPrompt) {
				printf("Should I change the cell? [Y]/N: ");
#ifdef SSI_PC
				CCrclHelper::GetInstance()->vBuildCmdCellNewCellFromLaue(a6fOldCell, 
                                                                         fCellParams, 
                                                                         a3fOldOrientAngles, 
                                                                         fOrientAngles,
                                                                         aoLaueTable[nBestSel].nLaueGroup);
#else
				fflush(stdout);
#endif            
				getline(cin,sTemp);

#ifdef SSI_PC
				printf("%s\n", sTemp.string());
#else
				printf("\n");
#endif
				sTemp.upcase();
				if(sTemp.length() > 0) {
					nResponse = (sTemp.GetAt(0) == 'Y');
				}
			} else
				nResponse = 1;
			
			if (nResponse != 1) {
				vCopyVecND(6,a6fOldCell,fCellParams);
				vCopyVecND(6,a6fOldCellSigma,fCellSigmas);
				vCopyVecND(3,a3fOldOrientAngles,fOrientAngles);

				printf("\nCrystal parameters reset to original input parameters:\n");
				printf("Crystal Parameters:  [%4.2f  %4.2f  %4.2f  %4.2f  %4.2f  %4.2f]\n",
						fCellParams[0],fCellParams[1],fCellParams[2],
						fCellParams[3],fCellParams[4],fCellParams[5]
						);
				printf("Cell Sigmas:         [%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ]\n",
					fCellSigmas[0],fCellSigmas[1],fCellSigmas[2],
					fCellSigmas[3],fCellSigmas[4],fCellSigmas[5]
					);

				printf("Crystal Orientation: [%4.2f %4.2f %4.2f]\n",
					fOrientAngles[0],fOrientAngles[1],fOrientAngles[2]);

			};

		}
        else
        {
#ifdef SSI_PC
            CCrclHelper::GetInstance()->vBuildCmdCellNewCellFromLaue(a6fOldCell, 
                                                                     a6fOldCell, 
                                                                     a3fOldOrientAngles, 
                                                                     a3fOldOrientAngles,
                                                                     -1);
            CCrclHelper::GetInstance()->SendTclCmd();
#endif 
        }
    } 
    else
    {
        m_eLaue = LAUE_UNKNOWN;
        return 1;
    }

    m_eLaue = (eLaueType)aoLaueTable[nBestSel].nLaueGroup;
    strncpy(&a2cLattice[0],&aoLaueTable[nBestSel].acLattice[0],3);
    // Allow the Laue check spacegroup to over-ride the space group already present.
    m_nSpaceGroup = -1;
    
    return 0;
};

#if 0 
// This centricity check is experimental.
// Please refer to tjn for details.
// Do not use it in it's present form.
int Claue::nCentricTest(Creflnlist& oRefList) {
    int nx,ny;
    double f0,f1;
    
    double fAverage; 
    double fAverageDeviation;
    double fAverageAbsDeviation;
    double a2x10fTheoreticalPercent[2][10]={
        {0.248,0.345,0.419,0.479,0.520,0.561,0.597,0.629,0.657,0.683},
        {0.095,0.181,0.259,0.330,0.394,0.451,0.503,0.551,0.593,0.632}
    };
    double a10fMeasuredPercent[10];
    double a10fMeasuredDeviation[10];
    int nNumReflnsUsed,nNumReflnsUsedLoop;
    bool bRet;
    
    char* cpPrintString =
        "N(Z) test: percent of intensities less than Z x <I>\n"
        "----------------------------------------------------------------------------\n"
        "       Z=          0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0\n"
        "  centric         24.8  34.5  41.9  47.9  52.0  56.1  59.7  62.9  65.7  68.3\n"
        " acentric          9.5  18.1  25.9  33.0  39.4  45.1  50.3  55.1  59.3  63.2\n"
        "----------------------------------------------------------------------------\n"
        " fit-to-centric  %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n"
        " deviations*     %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n"
        "----------------------------------------------------------------------------\n"
        " fit-to-acentric %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n"
        " deviations**    %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n"
        "----------------------------------------------------------------------------\n"
        " * = Ignored weakest %.1lf%% of data and shifted <I> by factor of %.1lf\n"
        " **= Ignored weakest %.1lf%% of data and shifted <I> by factor of %.1lf\n\n"
        "The N(Z) test is most informative for intensities below 0.3 x <I>\n";
        
    
    const double a3fAvgMultRange[3] = { 1.0,1.0,1.2 };
    const double a3fIoverSigTruncRange[3] = { 0.0,0.5,0.025};
    
    //const double a3fAvgMultRange[3] = { 1.0,1.0,2.0 };
    //const double a3fIoverSigTruncRange[3] = { 0.00001,0.00001,2.0};

    struct tagCentricParams {
        double a2fParams[2];
        double fAverageDeviation;
        double fAverageAbsDeviation;
        bool   bInUse;
        double a10fMeasuredDeviation[10];
        double a10fMeasuredPercent[10];

    };
    tagCentricParams aoBestFits[2];
    double a2fParams[2];
    itr<double> afIntensity;
    itr<double> afSigma;
    itr<int>    anRank;
    int nAverage;
    double* pfIntensity;
    double* pfSigma;
    int* pnRank;

    nNumReflnsUsed = 0;
    for (nx=0; nx < oRefList.nGetNumReflns(); nx++) {
        if ((oRefList[nx].nGetH()!=0) && (oRefList[nx].nGetK()!=0) && (oRefList[nx].nGetL()!=0)) {
            afIntensity + (double) oRefList[nx].fGetIntensity();
            afSigma + (double) oRefList[nx].fGetSigmaI();
            nNumReflnsUsed++;
        };
    };
    pfIntensity = & afIntensity[0];
    pfSigma = & afSigma[0];

    anRank.setsize(nNumReflnsUsed);
    pnRank = &anRank[0];
    g_pfCmpDoubles = pfIntensity;
    for (nx=0;nx < nNumReflnsUsed;nx++)
        pnRank[nx] = nx;
    qsort(pnRank,nNumReflnsUsed,sizeof(int),double_cmp_rel);
    for (nx= 0; nx < nNumReflnsUsed/2; nx++)
        std::swappnRank[nx],pnRank[nNumReflnsUsed - 1 - nx]);

    aoBestFits[0].bInUse = false;
    aoBestFits[1].bInUse = false;
      
    for (a2fParams[0] = a3fAvgMultRange[0]; a2fParams[0] <= a3fAvgMultRange[1]; a2fParams[0]*=a3fAvgMultRange[2]) {
        for (a2fParams[1] = a3fIoverSigTruncRange[0]; a2fParams[1] <= a3fIoverSigTruncRange[1]; a2fParams[1] +=a3fIoverSigTruncRange[2]) {
            
            double fPercentReject = a2fParams[1];
            double fAdjustedAverage;

            fAverage = 0.0;
            nAverage = 0;
            nNumReflnsUsedLoop = (1.0 - fPercentReject)*nNumReflnsUsed;
            for (nx=0; nx < nNumReflnsUsedLoop; nx++) {
                    fAverage += pfIntensity[pnRank[nx]];
                    nAverage++;
            };
            fAverage /= max(1,nAverage);
            fAdjustedAverage = a2fParams[0]*fAverage;
            
            for (nx=0; nx < 10; nx++) 
                a10fMeasuredPercent[nx]=0;
            
            // Figure out percent of reflections in each category.
            for (nx=0; nx < nNumReflnsUsedLoop; nx++) {
                
                f0 = pfIntensity[pnRank[nx]];
                f1 = pfSigma[pnRank[nx]];
                
                for (ny=0; ny<10; ny++) {
                    if (f0<=(0.1*(ny+1)*fAdjustedAverage))
                        a10fMeasuredPercent[ny]+=1;
                };
                
            };

            for (ny=0;ny<10; ny++) {
                a10fMeasuredPercent[ny]/=max(1,nNumReflnsUsed);        
            };

            for (nx = 0; nx < 2; nx++) {
                fAverageDeviation = 0.0;
                fAverageAbsDeviation = 0.0;
                for (ny = 0; ny < 10; ny++) {
                    a10fMeasuredDeviation[ny]=a2x10fTheoreticalPercent[nx][ny]-a10fMeasuredPercent[ny];
                    fAverageDeviation+=a10fMeasuredDeviation[ny];
                    fAverageAbsDeviation += ABS(a10fMeasuredDeviation[ny]);
                };
                fAverageDeviation /= 10.0;
                fAverageAbsDeviation /= 10.0;
                if ((!aoBestFits[nx].bInUse) || (aoBestFits[nx].fAverageAbsDeviation > fAverageAbsDeviation)) {
                    aoBestFits[nx].bInUse = true;
                    aoBestFits[nx].a2fParams[0] = a2fParams[0];
                    aoBestFits[nx].a2fParams[1] = a2fParams[1];
                    aoBestFits[nx].fAverageAbsDeviation = fAverageAbsDeviation;
                    aoBestFits[nx].fAverageDeviation = fAverageDeviation;
                    memcpy(&aoBestFits[nx].a10fMeasuredDeviation[0],&a10fMeasuredDeviation[0],10*sizeof(double));
                    memcpy(&aoBestFits[nx].a10fMeasuredPercent[0],&a10fMeasuredPercent[0],10*sizeof(double));
                };
            };
        };
    };

    
    
    printf(cpPrintString,
        100.0*aoBestFits[0].a10fMeasuredPercent[0],100.0*aoBestFits[0].a10fMeasuredPercent[1],100.0*aoBestFits[0].a10fMeasuredPercent[2],100.0*aoBestFits[0].a10fMeasuredPercent[3],100.0*aoBestFits[0].a10fMeasuredPercent[4],
        100.0*aoBestFits[0].a10fMeasuredPercent[5],100.0*aoBestFits[0].a10fMeasuredPercent[6],100.0*aoBestFits[0].a10fMeasuredPercent[7],100.0*aoBestFits[0].a10fMeasuredPercent[8],100.0*aoBestFits[0].a10fMeasuredPercent[9],
        100.0*aoBestFits[0].a10fMeasuredDeviation[0],100.0*aoBestFits[0].a10fMeasuredDeviation[1],100.0*aoBestFits[0].a10fMeasuredDeviation[2],100.0*aoBestFits[0].a10fMeasuredDeviation[3],100.0*aoBestFits[0].a10fMeasuredDeviation[4],
        100.0*aoBestFits[0].a10fMeasuredDeviation[5],100.0*aoBestFits[0].a10fMeasuredDeviation[6],100.0*aoBestFits[0].a10fMeasuredDeviation[7],100.0*aoBestFits[0].a10fMeasuredDeviation[8],100.0*aoBestFits[0].a10fMeasuredDeviation[9],
        100.0*aoBestFits[1].a10fMeasuredPercent[0],100.0*aoBestFits[1].a10fMeasuredPercent[1],100.0*aoBestFits[1].a10fMeasuredPercent[2],100.0*aoBestFits[1].a10fMeasuredPercent[3],100.0*aoBestFits[1].a10fMeasuredPercent[4],
        100.0*aoBestFits[1].a10fMeasuredPercent[5],100.0*aoBestFits[1].a10fMeasuredPercent[6],100.0*aoBestFits[1].a10fMeasuredPercent[7],100.0*aoBestFits[1].a10fMeasuredPercent[8],100.0*aoBestFits[1].a10fMeasuredPercent[9],
        100.0*aoBestFits[1].a10fMeasuredDeviation[0],100.0*aoBestFits[1].a10fMeasuredDeviation[1],100.0*aoBestFits[1].a10fMeasuredDeviation[2],100.0*aoBestFits[1].a10fMeasuredDeviation[3],100.0*aoBestFits[1].a10fMeasuredDeviation[4],
        100.0*aoBestFits[1].a10fMeasuredDeviation[5],100.0*aoBestFits[1].a10fMeasuredDeviation[6],100.0*aoBestFits[1].a10fMeasuredDeviation[7],100.0*aoBestFits[1].a10fMeasuredDeviation[8],100.0*aoBestFits[1].a10fMeasuredDeviation[9],
        aoBestFits[0].a2fParams[1]*100.0,aoBestFits[0].a2fParams[0],
        aoBestFits[1].a2fParams[1]*100.0,aoBestFits[1].a2fParams[0]
        );
        
    if (aoBestFits[1].fAverageAbsDeviation < aoBestFits[0].fAverageAbsDeviation)
        bRet=FALSE;
    else
        bRet=TRUE;


    if (bRet==FALSE) 
        printf("\nAcentric Distribution Selected.\n");
    else if (bIsChiral)
      { 
        printf("\nAlthough statistics might suggest a centric distribution,\nAcentric Distribution Selected because CHIRAL has precedence.\n");
        bRet = FALSE;
      }
    else
        printf("\nCentric Distribution Selected.\n");

#ifdef SSI_PC
	CCrclHelper::GetInstance()->vBuildCmdCellCentricityTable( a10fMeasuredPercent, a10fMeasuredDeviation, fAverageDeviation, bRet );
	CCrclHelper::GetInstance()->SendTclCmd();
#endif

    bIsCentric = bRet;
    return 0;
};
#endif
 
int Claue::nCentricTest(Creflnlist& oRefList) {
    int nx,ny;
    double f0;

    double fAverage; 
    double fAverageDeviation;
    double af10TheoreticalPercent[10]={0.248,0.345,0.419,0.479,0.520,0.561,0.597,0.629,0.657,0.683};
    double af10MeasuredPercent[10];
    double af10MeasuredDeviation[10];
    int nNumReflns;
    int nNumReflnsUsed;
    bool bRet;

char* cpPrintString =
    "N(Z) test: fraction of intensities less than Z x <I>\n"
    "-------------------------------------------------------------------------------\n"
    "Z=           0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1.0\n"
    "centric     .248   .345   .419   .479   .520   .561   .597   .629   .657   .683\n"
    "acentric    .095   .181   .259   .330   .394   .451   .503   .551   .593   .632\n"
    "deviation  -.153  -.164  -.160  -.149  -.126  -.110  -.094  -.078  -.064  -.051\n"
    "                             theoretical average deviation ==>   -.115\n\n"
    "measured  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
    "deviation %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n"
    "                                measured average deviation ==>   %6.3f\n"
    "-------------------------------------------------------------------------------\n";

    nNumReflns=oRefList.nGetNumReflns();
    // Get average value for intensities.
    fAverage=0.0;
    nNumReflnsUsed = 0;
    for (nx=0;nx<nNumReflns;nx++) {
        if ((oRefList[nx].nGetH()!=0) && (oRefList[nx].nGetK()!=0) && (oRefList[nx].nGetL()!=0)) {
            nNumReflnsUsed++;
            fAverage+=oRefList[nx].fGetIntensity();
        };
    };
    fAverage/=max(1,nNumReflnsUsed);
    for (nx=0; nx < 10; nx++) 
        af10MeasuredPercent[nx]=0;

    // Figure out percent of reflections in each category.
    for (nx=0; nx < nNumReflns; nx++) {
        if ((oRefList[nx].nGetH()!=0) && (oRefList[nx].nGetK()!=0) && (oRefList[nx].nGetL()!=0)) {
            f0=oRefList[nx].fGetIntensity();
            for (ny=0; ny<10; ny++)
                if (f0<=(0.1*(ny+1)*fAverage))
                    af10MeasuredPercent[ny]+=1;
        };
    };
    for (ny=0;ny<10; ny++) {
        af10MeasuredPercent[ny]/=max(1,nNumReflnsUsed);
        af10MeasuredDeviation[ny]=af10TheoreticalPercent[ny]-af10MeasuredPercent[ny];
    };
    fAverageDeviation=0.0;
    for (ny=0; ny < 10; ny++) 
        fAverageDeviation+=af10MeasuredDeviation[ny];
    fAverageDeviation/=10.0;
    
    printf(cpPrintString,
        af10MeasuredPercent[0],af10MeasuredPercent[1],af10MeasuredPercent[2],af10MeasuredPercent[3],af10MeasuredPercent[4],
        af10MeasuredPercent[5],af10MeasuredPercent[6],af10MeasuredPercent[7],af10MeasuredPercent[8],af10MeasuredPercent[9],
        af10MeasuredDeviation[0],af10MeasuredDeviation[1],af10MeasuredDeviation[2],af10MeasuredDeviation[3],af10MeasuredDeviation[4],
        af10MeasuredDeviation[5],af10MeasuredDeviation[6],af10MeasuredDeviation[7],af10MeasuredDeviation[8],af10MeasuredDeviation[9],
        fAverageDeviation);

    
    if (fAverageDeviation>=0.05)
        bRet=FALSE;
    else
        bRet=TRUE;


    if (bRet==FALSE) 
        printf("\nAcentric Distribution Selected.\n");
    else if (m_bIsChiral)
      { 
        printf("\nAlthough statistics might suggest a centric distribution,\nAcentric Distribution Selected because CHIRAL has precedence.\n");
        bRet = FALSE;
      }
    else
        printf("\nCentric Distribution Selected.\n");

#ifdef SSI_PC
	CCrclHelper::GetInstance()->vBuildCmdCellCentricityTable( af10MeasuredPercent, af10MeasuredDeviation, fAverageDeviation, bRet );
	CCrclHelper::GetInstance()->SendTclCmd();
#endif

    m_bIsCentric = bRet;
    return 0;
};



#define TP(x)     m_nTexsanTotal[x],m_nTexsanObsvd[x],m_fTexsanIoverSig[x],((m_fTexsanIoverSig[x]<=fAbsentSigma)?'*':' ')


int Claue::nPrintTexsan() {
    printf("BEGIN Texsan Output\n"
           "-----------------------------------------------------------------------------\n");

    printf(g_pcTexsanSysAbsString,
        TP(TEX_EEE),TP(TEX_EE0),TP(TEX_EOE),TP(TEX_EOO),TP(TEX_OEE),TP(TEX_OEO),TP(TEX_OOE),TP(TEX_OOO),
        TP(TEX_HK0_EE),TP(TEX_H0L_EE),TP(TEX_HK0_EO),TP(TEX_H0L_EO),TP(TEX_HK0_OE),TP(TEX_H0L_OE),TP(TEX_HK0_OO),TP(TEX_H0L_OO),
        TP(TEX_0KL_EE),TP(TEX_HHL_EE),TP(TEX_0KL_EO),TP(TEX_HHL_EO),TP(TEX_0KL_OE),TP(TEX_HHL_OE),TP(TEX_0KL_OO),TP(TEX_HHL_OO),
        TP(TEX_H_HL_EE),TP(TEX_H_HL_EO),TP(TEX_H_HL_OE),TP(TEX_H_HL_OO),
        TP(TEX_HHH_E),TP(TEX_HH0_E),TP(TEX_HHH_O),TP(TEX_HH0_O),
        TP(TEX_H00_E),TP(TEX_0K0_E),TP(TEX_H00_O),TP(TEX_0K0_O), m_fTexsanIoverSig[TEX_H00_O]/max(1.0,m_fTexsanIoverSig[TEX_H00_E]),m_fTexsanIoverSig[TEX_0K0_O]/max(1.0,m_fTexsanIoverSig[TEX_0K0_E]),
        TP(TEX_00L_E),TP(TEX_00L_O),m_fTexsanIoverSig[TEX_00L_O]/max(1.0,m_fTexsanIoverSig[TEX_00L_E]),
        TP(TEX_0KL_A),TP(TEX_0KL_NA),TP(TEX_H0L_A),TP(TEX_H0L_NA),TP(TEX_HK0_A),TP(TEX_HK0_NA),
        TP(TEX_0K0_A),TP(TEX_0K0_NA),TP(TEX_00L_A),TP(TEX_00L_NA),TP(TEX_H00_A),TP(TEX_H00_NA),
        TP(TEX_HHL_A),TP(TEX_HHL_NA),
        TP(TEX_H_H0L_A1),TP(TEX_H_H0L_A2),TP(TEX_H_H0L_A3),TP(TEX_H_H0L_B1),TP(TEX_H_H0L_B2),TP(TEX_H_H0L_B3),
        TP(TEX_000L_A),TP(TEX_000L_NA),TP(TEX_000L_B),TP(TEX_000L_NB)        
        );
    printf("\n-----------------------------------------------------------------------------\n"
           "END Texsan Output\n\n\n"
        );
    return 0;
};


int Claue::nReindexRelative(Creflnlist& CReflnFrom,Creflnlist& CReflnTo,double a3x3fOutputReindex[3][3],double& fScaleFactor) {
    int nx,ny,nz,nw;
    double f0;
    Cstring sTemp;
    
    Creflnlist CReflnTemp(CReflnTo);
    
    double a3x3fReindexMat[3][3];
    double a3x3fTempMat1[3][3];
    double a3x3fTempMat2[3][3];
    double a3x3fTempMat3[3][3];
    double fLowest = 1e20;
    Ccrystal oCrystal;
    eIndexTransform  eClass1Lowest;
    eIndexTransform  eClass2Lowest;
    const int nNumClass1 = 4;
    const int nNumClass2 = 12;
    const double fSelectTol = 0.01;
    int nClass1,nClass2;
    int nBatchIndex;
    double fScale;
    eIndexTransform eClass1[nNumClass1] = {TRANS_ABC,TRANS_F12,TRANS_F23,TRANS_F13};
    eIndexTransform eClass2[nNumClass2] = {TRANS_ABC,TRANS_CAB,TRANS_BCA,TRANS_BAC,TRANS_ACB,TRANS_CBA,TRANS_C21,TRANS_C31,TRANS_C12,TRANS_C32,TRANS_C13,TRANS_C23};

    // These arrays contain the selected matrices.  Note that at the termination of the checking loop, 
    // not all of the "selected" matrices will be valid, since the ones lower in the list might give appreciably
    // better Rmerge.
    eIndexTransform eSelected1[nNumClass1*nNumClass2];
    eIndexTransform eSelected2[nNumClass1*nNumClass2];
    double           fSelectedScale[nNumClass1*nNumClass2];
    double           fSelectedR[nNumClass1*nNumClass2];
    double           a3x3fSelected[nNumClass1*nNumClass2][3][3];
    int nNumSelected;

    // These arrays contain a list of all UNIQUE
    // When communicating with CrystalClear, we must provided a count of the TOTAL number of matrices.  Since
    // This requires counting them beforehand.
    eIndexTransform eClass1Use[nNumClass1*nNumClass2];
    eIndexTransform eClass2Use[nNumClass1*nNumClass2];
    int nNumUse;
    int nUseCount;


       
    
    nNumSelected = 0;
    // Initialize a crystal object.
    oCrystal.m_poSpacegroup->vSet(m_nSpaceGroup);               
    // Add a batch identifier field.
    nBatchIndex = CReflnTemp.nExpandGetField("nTemp");
    
    // Do a count of the number of matrices that we will be checking.
    nNumUse = 0;
    for (nClass1=0;nClass1<nNumClass1;nClass1++) {
        for (nClass2=0;nClass2<nNumClass2;nClass2++) {
            vCopyMat3D(g_fReindexTransMats[eClass1[nClass1]],&a3x3fTempMat1[0][0]);
            vCopyMat3D(g_fReindexTransMats[eClass2[nClass2]],&a3x3fTempMat2[0][0]);
            vMulMat3DMat3D(a3x3fTempMat1,a3x3fTempMat2,a3x3fReindexMat);
            // Check all previous matrices to see if we already did this transformation.
            for (nx=0;nx<nNumClass1;nx++) {
                for (ny=0;ny<nNumClass2;ny++) {
                    if ((nx<nClass1) || ((nx==nClass1) && (ny<nClass2))) {
                        vCopyMat3D(g_fReindexTransMats[eClass1[nx]],&a3x3fTempMat1[0][0]);
                        vCopyMat3D(g_fReindexTransMats[eClass2[ny]],&a3x3fTempMat2[0][0]);
                        vMulMat3DMat3D(a3x3fTempMat1,a3x3fTempMat2,a3x3fTempMat3);
                        for (nw=0;nw<3;nw++) {
                            for (nz=0;nz<3;nz++) {
                                if (a3x3fTempMat3[nw][nz]!=a3x3fReindexMat[nw][nz])
                                    break;
                            };
                            if (nz!=3)
                                break;
                        };
                        if (nw==3)
                            break;
                    };
                };
                if (ny<nNumClass2)
                    break;
            };
            if (nx<nNumClass1)
                continue;
            eClass1Use[nNumUse] = eClass1[nClass1];
            eClass2Use[nNumUse] = eClass2[nClass2];
            nNumUse++;
        };
    };


#ifdef SSI_PC
	CCrclHelper::GetInstance()->vBuildCmdCellNumReindexRelativeChecks( nNumUse );
	CCrclHelper::GetInstance()->SendTclCmd();
#endif
    
    // We want to reindex sFileFrom to sFileTo
    for (nUseCount = 0,nClass1=0;nClass1<nNumClass1;nClass1++) {
        for (nClass2=0;nClass2<nNumClass2;nClass2++) {

            if ((nUseCount>=nNumUse) || (eClass1Use[nUseCount]!=eClass1[nClass1]) || (eClass2Use[nUseCount]!=eClass2[nClass2]))
                continue;
            nUseCount++;

            vCopyMat3D(g_fReindexTransMats[eClass1[nClass1]],&a3x3fTempMat1[0][0]);
            vCopyMat3D(g_fReindexTransMats[eClass2[nClass2]],&a3x3fTempMat2[0][0]);
            vMulMat3DMat3D(a3x3fTempMat1,a3x3fTempMat2,a3x3fReindexMat);
            
            printf("Checking h'=[A][B]h  A = %-8s B = %-8s ",
                g_cpTransforms[eClass1[nClass1]],
                g_cpTransforms[eClass2[nClass2]]);
            
            CReflnTemp.vDeleteAll();
            CReflnTemp.nInsertListFrom(CReflnFrom);
            for (nx=0;nx<CReflnTemp.nGetNumReflns();nx++)
                CReflnTemp[nx].vSetField(nBatchIndex,0);
            CReflnTemp.nReindex(&a3x3fReindexMat[0][0]);
            CReflnTemp.nInsertListFrom(CReflnTo);
            for (;nx<CReflnTemp.nGetNumReflns();nx++)
                CReflnTemp[nx].vSetField(nBatchIndex,1);
            CReflnTemp.nReduce(oCrystal, 0);
            
            // Compute an Rmerge.
            
            int nGroups;
            int nGroup;
            int nLastHKL,nThisHKL;
            int nStartIndex;
            int nEndIndex;
            double fNumer;
            double fDenom;
            double fGroupAvg;
            double fAvgGroupSize;
            double fSumTo,fSumFrom;      // Two sums... are used in pass #1 to compute an optimum scale factor.
            int nRefSort,nRef;
            int nPass;
            
            // We make two passes.  On the first, we calculate the scale factor necc. to apply
            // to the first data set that yeilds the smallest Rmerge.
            // On the second, we compute the Rmerge.
            
            for (nPass=0;nPass<2;nPass++) {
                fSumTo = 0.0;
                fSumFrom = 0.0;
                
                nGroups=0;
                nLastHKL = CReflnTemp[CReflnTemp.pnGetSortIndex()[0]].nGetField(CReflnTemp.m_nFI_nPackedHKL);
                nStartIndex = 0;
                fNumer = 0.0;
                fDenom = 0.0;
                fAvgGroupSize=0;
                
                
                for (nRefSort=0;nRefSort<CReflnTemp.nGetNumReflns()+1;nRefSort++)
                {
                    bool bHKLEqual;
                    
                    bool bAnomalous = true;
                    
                    nRef=CReflnTemp.pnGetSortIndex()[nRefSort];               
                    
                    if (nRefSort<CReflnTemp.nGetNumReflns())
                        nThisHKL = CReflnTemp[nRef].nGetField(CReflnTemp.m_nFI_nPackedHKL);
                    
                    bHKLEqual = ( ((bAnomalous == FALSE) && ((((unsigned int) nThisHKL )>> 1)==(((unsigned int) nLastHKL ) >> 1))) ||
                        ((bAnomalous == TRUE) && ((((unsigned int) nThisHKL ))==(((unsigned int) nLastHKL )))) );
                    
                    if ((nRefSort!=CReflnTemp.nGetNumReflns()) && (bHKLEqual))
                        continue;
                    nEndIndex=nRefSort-1;                      
                    
                    
                    // Find the average of the reflections in the group.
                    fGroupAvg = 0.0;
                    nGroup =0;
                    for (nx=nStartIndex;nx<=nEndIndex;nx++) {
                        f0 = CReflnTemp[CReflnTemp.pnGetSortIndex()[nx]].fGetIntensity();
                        fGroupAvg += f0;
                        nGroup++;
                        
                        if (CReflnTemp[CReflnTemp.pnGetSortIndex()[nx]].nGetField(nBatchIndex)==0) {
                            for (ny=nStartIndex;ny<=nEndIndex;ny++) {
                                if (CReflnTemp[CReflnTemp.pnGetSortIndex()[ny]].nGetField(nBatchIndex)==1) {
                                    fSumFrom += CReflnTemp[CReflnTemp.pnGetSortIndex()[nx]].fGetIntensity();
                                    fSumTo += CReflnTemp[CReflnTemp.pnGetSortIndex()[ny]].fGetIntensity();
                                };
                            };
                        };
                    };
                    
                    if (nGroup>1) {
                        
                        fGroupAvg/=nGroup;
                        for (nx=nStartIndex;nx<=nEndIndex;nx++) {
                            f0 = CReflnTemp[CReflnTemp.pnGetSortIndex()[nx]].fGetIntensity();
                            fNumer+=fabs(f0 - fGroupAvg);
                            fDenom+=fabs(fGroupAvg);
                        };
                        nGroups++;
                        fAvgGroupSize+=nGroup;
                    };
                    
                    nStartIndex=nRefSort;
                    if (nRefSort<CReflnTemp.nGetNumReflns())
                        nLastHKL =  CReflnTemp[nRef].nGetField(CReflnTemp.m_nFI_nPackedHKL);               
                }

                fAvgGroupSize/=nGroups;
                
                if (nPass==0) {
                    fScale = fSumTo/fSumFrom;
                    for (nx=0;nx<CReflnTemp.nGetNumReflns();nx++) {
                        if (CReflnTemp[nx].nGetField(nBatchIndex)==0)
                            CReflnTemp[nx].vSetIntensity(CReflnTemp[nx].fGetIntensity()*(float)fScale);
                        CReflnTemp[nx].vSetSigmaI(CReflnTemp[nx].fGetSigmaI()*(float)fScale);
                    };
                };
                
                f0 = fNumer/fDenom;
                
            };
            
            printf("Scale = %5.2f R = %5.2f ",fScale,f0);
            if ((f0-fLowest)/(max(f0,fLowest))<fSelectTol) { 
                
                if (f0<fLowest) {
                    fLowest = f0;
                    printf("*\nLower residual found R=%4.2f\n",fLowest);
                    eClass1Lowest=eClass1[nClass1];
                    eClass2Lowest=eClass2[nClass2];                    
                } else
                    printf("*\n");
                eSelected1[nNumSelected]=eClass1[nClass1];
                eSelected2[nNumSelected]=eClass2[nClass2];
                fSelectedScale[nNumSelected]=fScale;
                fSelectedR[nNumSelected]=f0;
                vCopyMat3D(&a3x3fReindexMat[0][0],&a3x3fSelected[nNumSelected][0][0]);
                nNumSelected++;
            } else
                printf("\n");
            
#ifdef SSI_PC
			CCrclHelper::GetInstance()->vBuildCmdCellReindexRelativeStep();
			CCrclHelper::GetInstance()->SendTclCmd();
#endif
            };

        };
        
        // For those transformations that were selected, we will print the HKL transform matrix.        
        for (nx=0,ny=0;nx<nNumSelected;nx++) {
            if ((fSelectedR[nx]-fLowest)/(max(fSelectedR[nx],fLowest))<fSelectTol) {
                printf("#%-3d  Scale = %4.2f  h' = [%s] * [%s] * h \n\n"
                "   [h']   [%5.2f %5.2f %5.2f]   [h]\n"
                "   [k'] = [%5.2f %5.2f %5.2f] * [k]\n"
                "   [l']   [%5.2f %5.2f %5.2f]   [l]\n\n",++ny,fSelectedScale[nx],
                g_cpTransforms[eSelected1[nx]],g_cpTransforms[eSelected2[nx]],
                a3x3fSelected[nx][0][0], a3x3fSelected[nx][1][0], a3x3fSelected[nx][2][0],
                a3x3fSelected[nx][0][1], a3x3fSelected[nx][1][1], a3x3fSelected[nx][2][1],
                a3x3fSelected[nx][0][2], a3x3fSelected[nx][1][2], a3x3fSelected[nx][2][2]);

                //JOHNE  fSelectedR[nx] contains the Rmerge error.
#ifdef SSI_PC
				CCrclHelper::GetInstance()->vAddReindexRelativeMat(a3x3fSelected[nx], 
                                                                  (float)fSelectedScale[nx], 
                                                                  (float)fSelectedR[nx]);
#endif

            };
        };
        do {
            printf("Please enter selection, or type 'q' to quit: ");
#ifdef SSI_PC
			CCrclHelper::GetInstance()->vBuildCmdCellReindexRelativeMats();
#endif
            getline(cin,sTemp);

#ifdef SSI_PC
			printf("%s\n", sTemp.string());
#endif

			sTemp.upcase();
            if (sTemp.length() <=0)
                continue;
            else if ( sTemp.GetAt(0) == 'Q')
                return 1;
            else
                nx = atoi(sTemp.string());
        } while ((nx>ny) || (nx<1));

        // We only listed a subset of the total number of matrices.  Thus, the number the user typed in
        // should be searched.
        for (nz=0,ny=0;ny<nx;nz++) {
            if ((fSelectedR[nz]-fLowest)/(max(fSelectedR[nz],fLowest))<fSelectTol)
                ny++;
        };
        nx=nz-1;

        vCopyMat3D(&a3x3fSelected[nx][0][0],&a3x3fOutputReindex[0][0]);

        printf("\n\nDo you wish to apply the scale factor of %4.2f to the data? Y/[N]: ",fSelectedScale[nx]);
        getline(cin,sTemp);
        if ((sTemp=="Y") || (sTemp=="y")) {
            fScaleFactor = fSelectedScale[nx];
        } else
            fScaleFactor = 1.0;
        return 0;
};

#ifdef SSI_PC
void Claue::vBuildCellInputTableCmd(Cstring &csCommand)
{    
	Ccrystal oCrystal;
	oCrystal.vSetCell(fCellParams);	
    csCommand += CCrclHelper::GetInstance()->sBuildReduceCellInputRow(
					fCellParams[0], fCellParams[1], fCellParams[2], 
					fCellParams[3], fCellParams[4], fCellParams[5],
					fOrientAngles[0], fOrientAngles[1], fOrientAngles[2], oCrystal.fCalcVolume() );
}
#endif
