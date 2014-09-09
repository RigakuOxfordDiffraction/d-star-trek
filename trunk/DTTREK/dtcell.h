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

#ifndef DT_CELL_H_INCLUDE
#define DT_CELL_H_INCLUDE

#ifdef CONIO
#include "conio.h"
#define GETCH getch()
#else
#define GETCH
#endif

#ifdef SSI_PC
#define getline(x,y) y = CCrclHelper::GetInstance()->SendTclCmd(x)
//#define getline(x,y) CCrclHelper::GetInstance()->SendTclCmd(x,y)
#endif

#define NUM_REINDEX_TRANSFORMS 25


enum eLaueType { 
	LAUE_UNKNOWN=-1,
	LAUE_1_BAR=0,
	LAUE_2_OVER_M,
	LAUE_MMM,
	LAUE_4_OVER_M,
	LAUE_4_OVER_MMM,
	LAUE_3_BAR,
	LAUE_3_BAR_M_1,
	LAUE_3_BAR_1_M,
	LAUE_6_OVER_M,
	LAUE_6_OVER_MMM,
	LAUE_M_3_BAR,
	LAUE_M_3_BAR_M,
	LAUE_TOTAL
};
#ifdef CCELL_DEFINE_DATA
double g_fMatReindexI[3][3]       ={{1,0,0},{0,1,0},{0,0,1}};
double g_fMatReindexP1[3][3]      ={{0,1,0},{0,0,1},{1,0,0}};
double g_fMatReindexP2[3][3]      ={{0,0,1},{1,0,0},{0,1,0}};
double g_fMatReindexAB[3][3]      ={{0,1,0},{1,0,0},{0,0,-1}};
double g_fMatReindexBC[3][3]      ={{-1,0,0},{0,0,1},{0,1,0}};
double g_fMatReindexAC[3][3]      ={{0,0,1},{0,-1,0},{1,0,0}};
double g_fMatReindex21[3][3]      ={{0,0,-1},{0,1,0},{1,0,-1}};
double g_fMatReindex31[3][3]      ={{-1,0,1},{0,1,0},{-1,0,0}};
// Non standard reductions.
double g_fMatReindex12[3][3]      ={{-1,0,1},{0,1,0},{-1,0,0}};
double g_fMatReindex32[3][3]      ={{0,0,-1},{0,1,0},{1,0,-1}};
double g_fMatReindex13[3][3]      ={{0,0,-1},{0,1,0},{1,0,-1}};
double g_fMatReindex23[3][3]      ={{-1,0,1},{0,1,0},{-1,0,0}};

double g_fMatReindexmImC[3][3]    ={{-1,0,1},{0,1,0},{-1,0,0}};
double g_fMatReindexmCmI[3][3]    ={{0,0,-1},{0,1,0},{1,0,-1}};
double g_fMatReindexRP[3][3]      ={{1,0,1},{-1,1,1},{0,-1,1}};

double g_fMatReindexPR[3][3]      ={{2.0/3.0,-1.0/3.0,-1.0/3.0},{1.0/3.0,1.0/3.0,-2.0/3.0},{1.0/3.0,1.0/3.0,1.0/3.0}};

double g_fMatReindexhPoC_1[3][3]    ={{1,1,0},{0,2,0},{0,0,1}};
double g_fMatReindexhPoC_2[3][3]    ={{1,-1,0},{1,1,0},{0,0,1}};
double g_fMatReindexhPoC_3[3][3]    ={{0,-2,0},{1,-1,0},{0,0,1}};

double g_fMatReindexoChP_1[3][3]    ={{1.0,-0.5,0},{0,0.5,0},{0,0,1}};
double g_fMatReindexoChP_2[3][3]    ={{0.5,0.5,0},{-0.5,0.5,0},{0,0,1}};
double g_fMatReindexoChP_3[3][3]    ={{-0.5,1,0},{-0.5,0,0},{0,0,1}};

double g_fMatReindexF12[3][3]     ={{-1,0,0},{0,-1,0},{0,0,1}};
double g_fMatReindexF23[3][3]     ={{1,0,0},{0,-1,0},{0,0,-1}};
double g_fMatReindexF13[3][3]     ={{-1,0,0},{0,1,0},{0,0,-1}};


double* g_fReindexTransMats[NUM_REINDEX_TRANSFORMS] =
                                    {&(g_fMatReindexI[0][0]),
                                     &(g_fMatReindexP1[0][0]),
                                     &(g_fMatReindexP2[0][0]),
                                     &(g_fMatReindexAB[0][0]),
                                     &(g_fMatReindexBC[0][0]),
                                     &(g_fMatReindexAC[0][0]),
                                     &(g_fMatReindex21[0][0]),
                                     &(g_fMatReindex31[0][0]),
                                     &(g_fMatReindex12[0][0]),
                                     &(g_fMatReindex32[0][0]),
                                     &(g_fMatReindex13[0][0]),
                                     &(g_fMatReindex23[0][0]),
                                     &(g_fMatReindexF12[0][0]),
                                     &(g_fMatReindexF23[0][0]),
                                     &(g_fMatReindexF13[0][0]),
                                     &(g_fMatReindexmImC[0][0]),
                                     &(g_fMatReindexmCmI[0][0]),
                                     &(g_fMatReindexRP[0][0]),
                                     &(g_fMatReindexPR[0][0]),
                                     &(g_fMatReindexhPoC_1[0][0]),
                                     &(g_fMatReindexhPoC_2[0][0]),
                                     &(g_fMatReindexhPoC_3[0][0]),
                                     &(g_fMatReindexoChP_1[0][0]),
                                     &(g_fMatReindexoChP_2[0][0]),
                                     &(g_fMatReindexoChP_3[0][0])
};
char* g_cpTransforms[]          = { "I (Identity)", "cab", "bca", "bac", "acb", "cba" , // These transformations are only used if we want to transform 
                                    "C2->C1", "C3->C1", "C1->C2","C3->C2","C1->C3","C2->C3",
                                   "Invert-AB","Invert-BC","Invert-AC",
                                   "mI->mC","mC->mI","hR->hP","hP->hR",
                                   "hP->oC (C1)","hP->oC (C2)","hP->oC (C3)",
                                   "oC->hP (C1)", "oC->hP (C2)", "oC->hP (C3)" };

bool g_bCrystalClearUses[]      = { TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
                                    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
                                    FALSE,FALSE,FALSE,
                                    TRUE,TRUE,TRUE,TRUE,
                                    TRUE,TRUE,TRUE,
                                    TRUE,TRUE,TRUE};

char* g_cpTransformStrings[]  = { "I","ABC",NULL,"P1","CAB",NULL,"P2","BCA",NULL,"AB","BAC",NULL,"BC","ACB",NULL,"AC","CBA",NULL,
                                   "C2-C1",NULL, "C3-C1",NULL, "C1-C2",NULL, "C3-C2",NULL, "C1-C3",NULL, "C2-C3",NULL,
                                    "F12",NULL,"F23",NULL,"F13",NULL,
                                    "MI-MC",NULL,"MC-MI",NULL,"HR-HP",NULL,"HP-HR",NULL,
                                    "HP-OC1",NULL,"HP-OC2",NULL,"HP-OC3",NULL,
                                    "OC1-HP",NULL,"OC2-HP",NULL,"OC3-HP",NULL
};

char* g_cpLaueName[] = {
"LAUE_1_BAR",
"LAUE_2_OVER_M",
"LAUE_MMM",
"LAUE_4_OVER_M",
"LAUE_4_OVER_MMM",
"LAUE_3_BAR",
"LAUE_3_BAR_M_1",
"LAUE_3_BAR_1_M",
"LAUE_6_OVER_M",
"LAUE_6_OVER_MMM",
"LAUE_M_3_BAR",
"LAUE_M_3_BAR_M"
};

char* g_cpLaueNameSmall[] = {
"-1","2/m","mmm","4/m","4/mmm","-3","-3m1","-31m","6/m","6/mmm","m-3","m-3m"
};

// These are used for determining the lattice character to apply when passing from centered to
// primitive before doing a cell reduction.
char* g_pcValidLattices[14] = { "aP","mP","mC","oP","oC","oF","oI","hR","hP","tP","tI","cP","cI","cF"};



#else
extern double g_fMatReindexI[3][3];
extern double g_fMatReindexP1[3][3];
extern double g_fMatReindexP2[3][3];
extern double g_fMatReindexAB[3][3];
extern double g_fMatReindexBC[3][3];
extern double g_fMatReindexAC[3][3];
extern double g_fMatReindex21[3][3];
extern double g_fMatReindex31[3][3];
extern double g_fMatReindexF12[3][3];
extern double g_fMatReindexF23[3][3];
extern double g_fMatReindexF13[3][3];

extern double g_fMatReindexmImC[3][3];
extern double g_fMatReindexmCmI[3][3];
extern double g_fMatReindexRP[3][3];
extern double g_fMatReindexPR[3][3];
extern double g_fMatReindexhPoC_1[3][3];
extern double g_fMatReindexhPoC_2[3][3];
extern double g_fMatReindexhPoC_3[3][3];
extern double g_fMatReindexoChP_1[3][3];
extern double g_fMatReindexoChP_2[3][3];
extern double g_fMatReindexoChP_3[3][3];

extern double* g_fReindexTransMats[NUM_REINDEX_TRANSFORMS];
extern char* g_cpTransforms[];
extern bool g_bCrystalClearUses[];
extern char* g_cpTransformStrings[];
extern char* g_cpLaueName[];
extern char* g_cpLaueNameSmall[];
extern char* g_pcValidLattices[14];
#endif

extern char* g_pcTexsanSysAbsString;
extern int g_anCentricSpacegroups[];

// Do NOT change the order of the restrictions or absences.  This is important for some table operations.

enum eIndexTransform { TRANS_ABC,TRANS_CAB,TRANS_BCA,TRANS_BAC,TRANS_ACB,TRANS_CBA,
                       TRANS_C21,TRANS_C31,TRANS_C12,TRANS_C32,TRANS_C13,TRANS_C23,
                       TRANS_F12,TRANS_F23,TRANS_F13 };

#define TEXSANCLASSES 80

class DTREK_EXPORT Claue
{
private:
    static int   ms_nlaue_group[];   // Data for the above.
    int m_nTexsanTotal[TEXSANCLASSES];
    int m_nTexsanObsvd[TEXSANCLASSES];
    double m_fTexsanIoverSig[TEXSANCLASSES];

public:
	// User Provided Settings
    char        a2cLattice[3];      // Lattice type.  Used in converting to a centered cell.
	double		fCellParams[6];		// Cell parameters.
    double      fCellSigmas[6];     // Sigma values for above.
    double      fOrientAngles[3];   // Orientation angles for crystal.
    double      fMosaicity;         // Mosaicity of the incoming crystal.
	
    int         m_nSpaceGroup;      // These are NOT written by any routine found below.  They are read by nWriteHeader()
    
    bool        m_bIsCentric;


    double      fResidMax;          // Residual used in cell reduction.
	double		fRMergeGreat;		// "Great" RMerge.
    double      fRMergeAwful;       // Maximum allowable RMerge.
    double      fMaxRMergeRatio;    // Maximum allowable ratio of RMerge of a laue class to the RMerge of a lower Laue class.
	double		fAbsentSigma;		// Maximum # of I/sigma(I) quantities above zero for an absent reflection.
    double      fLaueSigma;         // I/Sigma(I) for Laue check.
    bool        bAnom;              // If set, assumes I+ != I- in Laue check
    bool        bScaleForLaue;      // Do we do a batch scale before the Laue check.
	int			nScaleForLaueMinRefs; // Minimum # of reflections per scale batch.
    bool        bCollectRejects;    // Boolean flag.  Set if we collect rejects.
    int         nMaxCollectRejects; // Maximum number of rejects to collect per table entry.
    eLaueType   eMaxLaueCheck;      // Maximum laue group to check.
    int         nMinConclusiveAbsences;   // Minimum number of absences required before we decide conclusively whether we have systematic absences in a zone.
    int         nMinConclusiveExistences; // Minimum number of non-absences required before we decide conclusively whether we have systematic absences in a zone.
    
    bool        m_bIsChiral;          // Use only spacegroups that are compatible with chiral spacegroups.
    
    Cstring     m_strInconclusiveConditions;  // A list of reflection conditions which have too few reflections
    
    int         nTwinID;            // These are NOT written by any routine found below.  They are read by nWriteHeader()
    eLaueType   m_eLaue;              // These are NOT written by any routine found below.  They are read by nWriteHeader()

	int             nLaueCheck(Creflnlist&,bool& bReindexingDone,bool bPrompt);
    eLaueType       eLauefromSpacegroup(int nSpacegroup);
    int             nSpacegroupFromLaue(eLaueType eLaueIn);
    int             nSetCrystal(Cimage_header& oHeader);
    int             nSetChirality(bool bUseCrystalVolumeForDefault);
    int             nChooseSpaceGroup(Creflnlist& oRef, bool& bIndexingDone, bool bPrompt);
    int             nCentricTest(Creflnlist& oRefList);  // Test identical to the one found in Texsan
    bool            bCellChanges(double fPercentError,double* a6fCell,double a3x3fTransform[3][3],char* pcOutputString);               // Retruns TRUE if the given transformation changes the cell more than fPercentError.  Used by CrystalClear dialog box
    
    int             nReindexRelative(Creflnlist& oRefFrom,Creflnlist& oRefTo,double a3x3fOutputReindex[3][3],double& fScaleFactor);   // Reindexes one file relative to another.
    int             nCalcCell(double fReindexMat[3][3]);                                                   // Updates the cell induced by a reindexing matrix.
    int             nCalcReindex(double UBfrom[3][3],double UBto[3][3],double fReindexMat[3][3]);                                     // Calculates reindexing matrix induced.  Returns 0 if fails.
    int             nReindexError(double a3x3fReindexMat[3][3],Creflnlist& oReflnlist,bool bPrompt);                                 // Checks to see that reflections are not lost during reindexing.
    int             nWriteHeader(Cstring* psName,Cimage_header& poPrvHead);                                                          // Updates a header with information.
    int             nPrintTexsan();

#ifdef SSI_PC
	void			vBuildCellInputTableCmd(Cstring &csCommand);
#endif

	int             nPrint(Creflnlist&);
    int             nDup();
    Claue();
    ~Claue();
};

struct tagLaueTable {
    
	// These are filled in from the cell reduction.
	double a6fCell[6];
	double a3x3fReindexMat[3][3];
    double a3x3fReindexMatInv[3][3];
	double* pfOrientReindexMat;			// Matrix to get unqiue axis and/or hP->oC.
	double* pfOrientReindexMatInv;		// Inverse of above.
	char   acLattice[3];
	int	   nLaueGroup;					// Laue group to check.

	// These are filled in by eLaueCheck().
	bool bUserSelectable;
    char acLaueClass[20];
    char cUniqueAxis;
    int  nGroups;
    int  nPass;
    double fCalcMult;
    double fObserveMult;
    double fRmerge;
	double fRmergeRaw;
	double fRmergeMin;
    double fRatioTo1Bar;
	enum { ePass,eInvalid,eFail };
    char  acPassFail[20];
    int nNonIndexingRefs;

	friend int nLaueTableCompare(const void* pvT1,const void* pvT2);

    int nClear() { memset(this,0,sizeof(*this)); return 0;};
    int operator = (int nData) { if (nData==0) nClear(); return 0;};
    int operator == (tagLaueTable& oOther) { return !strcmp(acLaueClass,oOther.acLaueClass); };    
    tagLaueTable() { nClear(); };
};



enum eTexsanTableClasses {
TEX_EEE,
TEX_EE0,
TEX_EOE,
TEX_EOO,
TEX_OEE,        
TEX_OEO,
TEX_OOE,
TEX_OOO,
TEX_HK0_EE,
TEX_HK0_EO,
TEX_HK0_OE,
TEX_HK0_OO,
TEX_H0L_EE,
TEX_H0L_EO,
TEX_H0L_OE,
TEX_H0L_OO,
TEX_0KL_EE,
TEX_0KL_EO,
TEX_0KL_OE,
TEX_0KL_OO,
TEX_HHL_EE,
TEX_HHL_EO,
TEX_HHL_OE,
TEX_HHL_OO,
TEX_H_HL_EE,
TEX_H_HL_EO,
TEX_H_HL_OE,
TEX_H_HL_OO,
TEX_HHH_E,
TEX_HHH_O,
TEX_HH0_E,
TEX_HH0_O,
TEX_H00_E,
TEX_H00_O,
TEX_0K0_E,
TEX_0K0_O,
TEX_00L_E,
TEX_00L_O,
TEX_0KL_A,
TEX_0KL_NA,
TEX_H0L_A,
TEX_H0L_NA,
TEX_HK0_A,
TEX_HK0_NA,
TEX_0K0_A,
TEX_0K0_NA,
TEX_00L_A,
TEX_00L_NA,
TEX_H00_A,
TEX_H00_NA,
TEX_HHL_A,
TEX_HHL_NA,
TEX_H_H0L_A1,
TEX_H_H0L_A2,
TEX_H_H0L_A3,
TEX_H_H0L_B1,
TEX_H_H0L_B2,
TEX_H_H0L_B3,
TEX_000L_A,
TEX_000L_NA,
TEX_000L_B,
TEX_000L_NB
};
#endif // !DT_CELL_H_INCLUDE
